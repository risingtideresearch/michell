//! Hydrostatic equilibrium: solve the platform's sinkage and pitch so the
//! fleet displaces a given weight with its centre of buoyancy under a given
//! centre of gravity.
//!
//! Newton iteration on (sinkage s, pitch τ) with residuals
//! `(∇ − W/ρ, M_x − lcg·∇)`; the Jacobian is analytic from waterplane
//! properties (area, first and second longitudinal moments — the classical
//! tons-per-cm / moment-to-trim relations). Iterations run on a coarsened
//! sampling of the hulls, then the solution is polished on the requested
//! resolution. A hull lifting fully out of the water during iteration is
//! handled (its contribution drops to zero), and being dry at equilibrium is
//! reported, not an error.
//!
//! All of this is *hydrostatic*: no speed-dependent (dynamic) sinkage/trim.
//!
//! Transverse stability rides on [`solve_equilibrium_heeled`]: it holds the
//! displacement at a prescribed heel and reads the righting arm off the true
//! inclined-waterplane cut ([`crate::inclined`]).

use crate::body::{Body, BodyOptions};
use crate::error::{Error, Result};
use crate::hull::Hull;
use crate::iges::{HullPose, ImportOptions, Platform, SourceFleet};
use crate::inclined::{fleet_inclined, InclinedGrid};
use crate::michell::Placement;

/// What the platform must carry.
#[derive(Debug, Clone, Copy)]
pub struct LoadCase {
    /// Total mass [kg].
    pub mass: f64,
    /// Longitudinal centre of gravity [m], in the fleet's x coordinates.
    /// `None` locks the platform pitch at zero and balances weight only.
    pub lcg: Option<f64>,
}

/// The wetted fleet at some state: what the solver (and resistance) consume.
#[derive(Debug)]
pub struct FleetState {
    pub members: Vec<(Hull, Placement)>,
    /// Number of source hulls entirely above the water.
    pub dry: usize,
    /// Total wetted samples that fell above the members' modelled bands
    /// (bodies only): non-zero means missing topside geometry and
    /// under-counted buoyancy at this state.
    pub band_exceeded: usize,
}

/// A solved floating condition.
#[derive(Debug)]
pub struct Equilibrium {
    /// Solved additional immersion of the platform [m] (relative to the base
    /// waterline; negative = riding higher).
    pub sinkage: f64,
    /// Solved platform pitch [rad]; positive raises the +x end.
    pub trim: f64,
    /// The fleet situated at the solution, at full requested resolution.
    pub fleet: FleetState,
    /// Achieved displaced volume [m³].
    pub volume: f64,
    /// Achieved longitudinal centre of buoyancy [m].
    pub lcb: f64,
    /// Total waterplane area at the solution [m²].
    pub waterplane_area: f64,
    /// Newton iterations used (both phases).
    pub iterations: usize,
    /// |∇ − target| / target at the solution.
    pub volume_residual: f64,
    /// |LCB − lcg| [m] at the solution (0 when lcg is None).
    pub lcb_residual: f64,
}

struct Totals {
    volume: f64,
    moment_x: f64,
    wp_area: f64,
    wp_moment: f64,
    wp_second: f64,
    draft: f64,
}

fn totals(fleet: &FleetState) -> Totals {
    let mut t = Totals {
        volume: 0.0,
        moment_x: 0.0,
        wp_area: 0.0,
        wp_moment: 0.0,
        wp_second: 0.0,
        draft: 0.0,
    };
    for (hull, place) in &fleet.members {
        let v = hull.displaced_volume();
        t.volume += v;
        t.moment_x += (hull.lcb_x() + place.x) * v;
        t.wp_area += hull.waterplane_area();
        t.wp_moment += hull.waterplane_moment() + place.x * hull.waterplane_area();
        t.wp_second += hull.waterplane_second_moment()
            + 2.0 * place.x * hull.waterplane_moment()
            + place.x * place.x * hull.waterplane_area();
        t.draft = t.draft.max(hull.draft());
    }
    t
}

/// Generic equilibrium core. `situate(sinkage, trim, coarse)` produces the
/// fleet at a platform state (coarse = reduced sampling for iterations).
pub fn solve_equilibrium_with(
    mut situate: impl FnMut(f64, f64, bool) -> Result<FleetState>,
    load: &LoadCase,
    density: f64,
) -> Result<Equilibrium> {
    if !(load.mass.is_finite() && load.mass > 0.0) {
        return Err(Error::InvalidConditions(format!(
            "load mass must be finite and positive, got {}",
            load.mass
        )));
    }
    if !(density.is_finite() && density > 0.0) {
        return Err(Error::InvalidConditions("density must be positive".into()));
    }
    if let Some(l) = load.lcg {
        if !l.is_finite() {
            return Err(Error::InvalidConditions("lcg must be finite".into()));
        }
    }
    let v_target = load.mass / density;
    let z_guess = v_target.cbrt().max(1e-3);
    let pivot_x = load.lcg.unwrap_or(0.0);

    let mut s = 0.0f64;
    let mut tau = 0.0f64;
    let mut iterations = 0usize;

    for (coarse, max_iters, tol_v) in [(true, 30usize, 1e-3f64), (false, 20, 2e-4)] {
        let mut converged = false;
        // Adaptive relaxation: the waterplane-property Jacobian can
        // underestimate the true sensitivity (e.g. flare or structure
        // entering the water), which turns plain Newton into a limit cycle.
        // Halve the step whenever the volume residual flips sign, recover
        // gently while it doesn't.
        let mut relax = 1.0f64;
        let mut last_sign = 0.0f64;
        // Best state seen this phase, by tolerance-normalised residual.
        // The situate → loft model carries small-scale roughness (~0.1% of
        // volume), so Newton can stall dithering across a tolerance edge; a
        // near-miss is accepted with its residuals reported rather than
        // failing the whole point.
        let mut best = (f64::INFINITY, s, tau);
        for _ in 0..max_iters {
            iterations += 1;
            let fleet = situate(s, tau, coarse)?;
            if fleet.members.is_empty() {
                // Everything dry: sink until something gets wet.
                s += 0.5 * z_guess;
                continue;
            }
            let t = totals(&fleet);
            let z_scale = t.draft.max(z_guess);
            let l_scale = fleet
                .members
                .iter()
                .map(|(h, _)| h.length())
                .fold(0.0f64, f64::max)
                .max(1e-6);
            let r1 = t.volume - v_target;
            let lcb = t.moment_x / t.volume;
            if t.wp_area <= 1e-12 {
                return Err(Error::InvalidGeometry(
                    "waterplane area vanished during the equilibrium solve \
                     (hull fully submerged or degenerate)"
                        .into(),
                ));
            }
            let (mut ds, mut dtau) = match load.lcg {
                None => (r1 / t.wp_area, 0.0),
                Some(lcg) => {
                    let r2 = t.moment_x - lcg * t.volume;
                    // d/ds and d/dτ of (V, M); rotation about (lcg, waterline).
                    let dv_ds = t.wp_area;
                    let dv_dt = -(t.wp_moment - pivot_x * t.wp_area);
                    let dm_ds = t.wp_moment;
                    let dm_dt = -(t.wp_second - pivot_x * t.wp_moment);
                    let dr2_ds = dm_ds - lcg * dv_ds;
                    let dr2_dt = dm_dt - lcg * dv_dt;
                    let det = dv_ds * dr2_dt - dv_dt * dr2_ds;
                    if det.abs() < 1e-12 * (dv_ds.abs() * dr2_dt.abs()).max(1e-30) {
                        // Degenerate trim sensitivity: fall back to sinkage only.
                        (r1 / t.wp_area, 0.0)
                    } else {
                        (
                            (r1 * dr2_dt - dv_dt * r2) / det,
                            (dv_ds * r2 - dr2_ds * r1) / det,
                        )
                    }
                }
            };
            // Oscillation-adaptive relaxation, then absolute damping caps.
            let sign = r1.signum();
            if last_sign != 0.0 && sign != last_sign {
                relax = (relax * 0.5).max(0.1);
            } else {
                relax = (relax * 1.25).min(1.0);
            }
            last_sign = sign;
            ds = (ds * relax).clamp(-0.3 * z_scale, 0.3 * z_scale);
            dtau = (dtau * relax).clamp(-0.05, 0.05);
            let done_v = r1.abs() <= tol_v * v_target;
            let done_m = match load.lcg {
                None => true,
                Some(lcg) => (lcb - lcg).abs() <= 1e-4 * l_scale,
            };
            let metric = (r1.abs() / (tol_v * v_target)).max(match load.lcg {
                None => 0.0,
                Some(lcg) => (lcb - lcg).abs() / (1e-4 * l_scale),
            });
            if metric < best.0 {
                best = (metric, s, tau);
            }
            if std::env::var("MICHELL_DEBUG_FLOAT").is_ok() {
                eprintln!(
                    "DBG iter={iterations} s={s:.6} tau={tau:.6} V={:.6} lcb={lcb:.6} \
                     Aw={:.4} Mw={:.4} Iw={:.4} ds={ds:.6} dtau={dtau:.6}",
                    t.volume, t.wp_area, t.wp_moment, t.wp_second
                );
            }
            if done_v && done_m {
                converged = true;
                break;
            }
            s -= ds;
            tau -= dtau;
            // Beyond ~20 degrees of pitch the load case is unreachable (and
            // thin-ship theory meaningless): fail with a diagnosis instead of
            // wandering.
            if tau.abs() > 0.35 {
                return Err(Error::InvalidConditions(format!(
                    "equilibrium trim exceeded 20 degrees; the requested lcg \
                     {:?} appears unreachable for this geometry (LCB range is \
                     limited by the hull; current LCB {lcb:.3} m)",
                    load.lcg
                )));
            }
        }
        if !converged {
            // Accept a stalled near-miss (within 10x tolerance — for the
            // fine phase, volume within 0.2% and LCB within 1e-3 of the
            // length); the achieved residuals are reported in the result.
            if best.0 <= 10.0 {
                (_, s, tau) = best;
            } else {
                return Err(Error::InvalidConditions(format!(
                    "equilibrium did not converge after {iterations} iterations \
                     (mass {} kg, lcg {:?}, best residual {:.1}x tolerance); the \
                     load may be outside what the geometry can float",
                    load.mass, load.lcg, best.0
                )));
            }
        }
    }

    // Final state at full resolution.
    let fleet = situate(s, tau, false)?;
    let t = totals(&fleet);
    let lcb = if t.volume > 0.0 {
        t.moment_x / t.volume
    } else {
        0.0
    };
    Ok(Equilibrium {
        sinkage: s,
        trim: tau,
        volume: t.volume,
        lcb,
        waterplane_area: t.wp_area,
        iterations,
        volume_residual: (t.volume - v_target).abs() / v_target,
        lcb_residual: load.lcg.map_or(0.0, |l| (lcb - l).abs()),
        fleet,
    })
}

/// Rigid platform heel applied to a set of body poses — the transverse
/// reposition of each demihull, used to build the wetted **geometry** at heel
/// (e.g. for resistance). Heel is a rotation about the platform's longitudinal
/// axis through the centerline at the design floatplane; **positive heel puts
/// the +y side down**. Each hull's transverse station `y = centerplane + dy`
/// moves to `y·cos φ` and gains `y·sin φ` of immersion. A half-breadth surface
/// cannot rotate about its own x axis, so this reposition does not tilt the
/// sections. For **hydrostatics and the righting arm** — which need the true
/// tilted cut — use [`solve_equilibrium_heeled`] instead (it applies heel as an
/// exact inclined-waterplane rotation, not this reposition).
pub fn heel_poses(bodies: &[&Body], poses: &[HullPose], heel: f64) -> Result<Vec<HullPose>> {
    if bodies.len() != poses.len() {
        return Err(Error::InvalidInput(format!(
            "{} poses supplied for {} bodies",
            poses.len(),
            bodies.len()
        )));
    }
    if !heel.is_finite() || heel.abs() >= std::f64::consts::FRAC_PI_2 {
        return Err(Error::InvalidConditions(format!(
            "heel angle must be finite and within ±90 degrees, got {} rad",
            heel
        )));
    }
    let (sin, cos) = heel.sin_cos();
    Ok(bodies
        .iter()
        .zip(poses)
        .map(|(body, pose)| {
            let y = body.centerplane() + pose.dy;
            HullPose {
                dy: y * cos - body.centerplane(),
                dz: pose.dz + y * sin,
                ..*pose
            }
        })
        .collect())
}

/// Equilibrium of an IGES source fleet at design poses (see
/// [`SourceFleet::situate`]). `waterline_z` is the base waterline the sinkage
/// is measured from.
pub fn solve_equilibrium(
    src: &SourceFleet,
    waterline_z: f64,
    poses: &[HullPose],
    load: &LoadCase,
    density: f64,
    opts: &ImportOptions,
) -> Result<Equilibrium> {
    let pivot_x = load.lcg.unwrap_or(0.0);
    let mut coarse_opts = *opts;
    coarse_opts.stations = (opts.stations / 2).clamp(31, opts.stations.max(31));
    coarse_opts.waterlines = (opts.waterlines / 2).clamp(11, opts.waterlines.max(11));
    coarse_opts.fit.n_ctrl_x = opts.fit.n_ctrl_x.min(10).max(opts.fit.degree_x + 1);
    coarse_opts.fit.n_ctrl_z = opts.fit.n_ctrl_z.min(7).max(opts.fit.degree_z + 1);
    solve_equilibrium_with(
        |s, tau, coarse| {
            let fl = src.situate(
                waterline_z,
                poses,
                &Platform {
                    sinkage: s,
                    trim: tau,
                    pivot_x,
                },
                if coarse { &coarse_opts } else { opts },
            )?;
            Ok(FleetState {
                dry: fl.dry.len(),
                band_exceeded: 0, // IGES situates carry the full geometry
                members: fl
                    .members
                    .into_iter()
                    .map(|m| (m.hull, m.placement))
                    .collect(),
            })
        },
        load,
        density,
    )
}

/// Equilibrium of an assembly of full-band bodies at design poses.
/// `water_offset` is the base water position below the design floatplane
/// (normally 0) that the solved sinkage is measured from.
pub fn solve_equilibrium_bodies(
    bodies: &[&Body],
    water_offset: f64,
    poses: &[HullPose],
    load: &LoadCase,
    density: f64,
    opts: &BodyOptions,
) -> Result<Equilibrium> {
    if bodies.len() != poses.len() {
        return Err(Error::InvalidInput(format!(
            "{} poses supplied for {} bodies",
            poses.len(),
            bodies.len()
        )));
    }
    let pivot_x = load.lcg.unwrap_or(0.0);
    let mut coarse_opts = *opts;
    coarse_opts.stations = (opts.stations / 2).clamp(31, opts.stations.max(31));
    coarse_opts.waterlines = (opts.waterlines / 2).clamp(11, opts.waterlines.max(11));
    coarse_opts.fit.n_ctrl_x = opts.fit.n_ctrl_x.min(10).max(opts.fit.degree_x + 1);
    coarse_opts.fit.n_ctrl_z = opts.fit.n_ctrl_z.min(7).max(opts.fit.degree_z + 1);
    solve_equilibrium_with(
        |s, tau, coarse| {
            let platform = Platform {
                sinkage: s,
                trim: tau,
                pivot_x,
            };
            let o = if coarse { &coarse_opts } else { opts };
            let mut members = Vec::new();
            let mut dry = 0usize;
            let mut band_exceeded = 0usize;
            for (body, pose) in bodies.iter().zip(poses) {
                match body.situate(water_offset, pose, &platform, o)? {
                    Some(sb) => {
                        band_exceeded += sb.band_exceeded;
                        members.push((sb.hull, sb.placement));
                    }
                    None => dry += 1,
                }
            }
            Ok(FleetState {
                members,
                dry,
                band_exceeded,
            })
        },
        load,
        density,
    )
}

/// A solved heeled floating condition with a true inclined-waterplane righting
/// arm (see [`crate::inclined`]).
#[derive(Debug)]
pub struct HeeledEquilibrium {
    /// Solved additional immersion below the design floatplane [m].
    pub sinkage: f64,
    /// Solved platform pitch [rad]; positive raises the +x end.
    pub trim: f64,
    /// Heel the fleet was solved at [rad].
    pub heel: f64,
    /// Displaced volume at the solution (inclined cut) [m³].
    pub volume: f64,
    /// Longitudinal centre of buoyancy at the solution [m].
    pub lcb: f64,
    /// Righting arm `GZ` [m] from the inclined cut — form stability included,
    /// no metacentric approximation. Positive rights the platform.
    pub gz: f64,
    /// Wetted fleet for the resistance path, situated at the solved attitude.
    /// (Heel enters the geometry as the rigid reposition of [`heel_poses`];
    /// the tilt of each hull is carried by the wave kernel, not re-clipped —
    /// the hydrostatics above use the exact inclined cut instead.)
    pub fleet: FleetState,
    /// Newton iterations used.
    pub iterations: usize,
    /// |∇ − target| / target at the solution.
    pub volume_residual: f64,
    /// |LCB − lcg| [m] at the solution (0 when lcg is None).
    pub lcb_residual: f64,
    /// Band-top-immersed samples in the inclined cut (deck under water).
    pub band_exceeded: usize,
}

/// Heel-aware equilibrium of an assembly of full-band bodies, holding the
/// load's displacement (and LCG, if given) at a **prescribed** heel angle, with
/// the hydrostatics — volume, trim balance, and the righting arm `gz` for the
/// centre of gravity `vcg` — taken from the true inclined-waterplane cut
/// ([`crate::inclined`]) rather than the metacentric approximation of
/// [`heel_poses`] + [`righting_arm`].
///
/// `poses` are the design (un-heeled) poses. The (sinkage, pitch) solve runs on
/// the proven equilibrium core over the rigid heel reposition ([`heel_poses`])
/// — robust for any geometry — and the reported volume, LCB, and `gz` are then
/// evaluated from the exact inclined cut at that attitude, so the righting arm
/// carries the full (nonlinear) form stability. `grid` sets the section
/// integration resolution. `heel = 0` reproduces [`solve_equilibrium_bodies`]
/// with `gz = 0`.
#[allow(clippy::too_many_arguments)]
pub fn solve_equilibrium_heeled(
    bodies: &[&Body],
    water_offset: f64,
    poses: &[HullPose],
    load: &LoadCase,
    density: f64,
    heel: f64,
    vcg: f64,
    opts: &BodyOptions,
    grid: InclinedGrid,
) -> Result<HeeledEquilibrium> {
    if bodies.len() != poses.len() {
        return Err(Error::InvalidInput(format!(
            "{} poses supplied for {} bodies",
            poses.len(),
            bodies.len()
        )));
    }
    if !(load.mass.is_finite() && load.mass > 0.0) {
        return Err(Error::InvalidConditions(format!(
            "load mass must be finite and positive, got {}",
            load.mass
        )));
    }
    if !(density.is_finite() && density > 0.0) {
        return Err(Error::InvalidConditions("density must be positive".into()));
    }
    if !heel.is_finite() || heel.abs() >= std::f64::consts::FRAC_PI_2 {
        return Err(Error::InvalidConditions(format!(
            "heel angle must be finite and within ±90 degrees, got {heel} rad"
        )));
    }
    if let Some(l) = load.lcg {
        if !l.is_finite() {
            return Err(Error::InvalidConditions("lcg must be finite".into()));
        }
    }

    let pivot_x = load.lcg.unwrap_or(0.0);

    // Solve (sinkage, pitch) with the proven equilibrium core on the
    // horizontal-cut reposition (`heel_poses`). This is robust across
    // geometries — including hulls with no freeboard, where a naive inclined
    // Newton stalls (sinking a fully immersed hull adds no volume and the
    // waterplane Jacobian vanishes). The reported hydrostatics below then come
    // from the exact inclined cut at the solved attitude, so the righting arm
    // carries the full (nonlinear) form stability.
    let heeled = heel_poses(bodies, poses, heel)?;
    let mut coarse_opts = *opts;
    coarse_opts.stations = (opts.stations / 2).clamp(31, opts.stations.max(31));
    coarse_opts.waterlines = (opts.waterlines / 2).clamp(11, opts.waterlines.max(11));
    coarse_opts.fit.n_ctrl_x = opts.fit.n_ctrl_x.min(10).max(opts.fit.degree_x + 1);
    coarse_opts.fit.n_ctrl_z = opts.fit.n_ctrl_z.min(7).max(opts.fit.degree_z + 1);
    let eq = solve_equilibrium_with(
        |s, tau, coarse| {
            let platform = Platform {
                sinkage: s,
                trim: tau,
                pivot_x,
            };
            let o = if coarse { &coarse_opts } else { opts };
            let mut members = Vec::new();
            let mut dry = 0usize;
            let mut band_exceeded = 0usize;
            for (body, pose) in bodies.iter().zip(&heeled) {
                match body.situate(water_offset, pose, &platform, o)? {
                    Some(sb) => {
                        band_exceeded += sb.band_exceeded;
                        members.push((sb.hull, sb.placement));
                    }
                    None => dry += 1,
                }
            }
            Ok(FleetState {
                members,
                dry,
                band_exceeded,
            })
        },
        load,
        density,
    )?;

    // The righting arm — and only it — uses the exact inclined cut at the
    // solved attitude. Displacement, LCB, sinkage, and trim stay on the same
    // (lofted) basis the solve balanced, so those columns match the upright
    // solver exactly; the inclined cut just replaces the metacentric GZ.
    let platform = Platform {
        sinkage: eq.sinkage,
        trim: eq.trim,
        pivot_x,
    };
    let f = fleet_inclined(bodies, water_offset, poses, &platform, heel, grid);
    let gz = if f.volume > 0.0 {
        f.moment_y / f.volume - vcg * heel.sin()
    } else {
        0.0
    };

    Ok(HeeledEquilibrium {
        sinkage: eq.sinkage,
        trim: eq.trim,
        heel,
        volume: eq.volume,
        lcb: eq.lcb,
        gz,
        fleet: eq.fleet,
        iterations: eq.iterations,
        volume_residual: eq.volume_residual,
        lcb_residual: eq.lcb_residual,
        band_exceeded: f.band_exceeded,
    })
}

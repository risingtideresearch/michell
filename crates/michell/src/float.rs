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
//! The balance is *hydrostatic* by default. At speed the hull also feels a
//! hydrodynamic vertical force and pitch moment (thin-ship dynamic sinkage
//! and trim); [`solve_equilibrium_dynamic_with`] folds a caller-supplied
//! [`DynamicLoad`] into the same Newton core as a forcing term.
//!
//! Transverse stability rides on [`solve_equilibrium_heeled`]: it holds the
//! displacement at a prescribed heel and reads the righting arm off the true
//! inclined-waterplane cut ([`crate::inclined`]).

use crate::body::{Body, BodyOptions};
use crate::conditions::STANDARD_GRAVITY;
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

/// Speed-dependent hydrodynamic load on the fleet at a platform state, as
/// handed to [`solve_equilibrium_dynamic_with`] by the caller's closure. Both
/// entries are in the fleet frame: `force_up` is the net vertical
/// hydrodynamic force [N, positive **up**]; `moment_bow_up` is the pitch
/// moment [N·m, positive **bow (+x) up**] about the pivot station
/// `load.lcg.unwrap_or(0.0)` at the waterline. A hull sucked down at speed
/// reports a negative `force_up`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct DynamicLoad {
    pub force_up: f64,
    pub moment_bow_up: f64,
}

/// A solved floating condition at speed: the hydrostatics of [`Equilibrium`]
/// balanced against a [`DynamicLoad`] (see [`solve_equilibrium_dynamic_with`]).
#[derive(Debug)]
pub struct DynamicEquilibrium {
    /// Solved additional immersion of the platform [m] (relative to the base
    /// waterline; negative = riding higher).
    pub sinkage: f64,
    /// Solved platform pitch [rad]; positive raises the +x end.
    pub trim: f64,
    /// The fleet situated at the solution, at full requested resolution.
    pub fleet: FleetState,
    /// Achieved displaced volume [m³] — the hydrostatic displacement alone,
    /// which differs from `mass / density` by the dynamic force's share.
    pub volume: f64,
    /// Achieved longitudinal centre of buoyancy [m].
    pub lcb: f64,
    /// The dynamic load at the solution, as evaluated on `fleet`.
    pub dynamic: DynamicLoad,
    /// `force_up` as a fraction of the weight (negative = suction / sinkage).
    /// Past a few tenths the quasi-static balance is being asked a lot of;
    /// that is reported here, never refused.
    pub lift_fraction: f64,
    /// Newton iterations used (both phases).
    pub iterations: usize,
    /// |ρg(∇ − target) + F_z| / (ρg·target) at the solution: the vertical
    /// imbalance as a fraction of the weight.
    pub volume_residual: f64,
    /// |LCB − lcg + M/(ρg∇)| [m] at the solution: the moment imbalance as an
    /// LCB shift (0 when lcg is None).
    pub lcb_residual: f64,
}

/// What the Newton core hands back, before it is packaged as an
/// [`Equilibrium`] or a [`DynamicEquilibrium`].
struct Solved {
    sinkage: f64,
    trim: f64,
    fleet: FleetState,
    totals: Totals,
    lcb: f64,
    /// Zero when the solve carried no dynamic closure.
    dynamic: DynamicLoad,
    iterations: usize,
    volume_residual: f64,
    lcb_residual: f64,
}

/// `x` plus the dynamic forcing, when there is one. A branch rather than
/// `x + f.unwrap_or(0.0)` so the hydrostatic path's arithmetic is untouched
/// by the dynamic variant (bit-for-bit, signed zeros included).
fn forced(x: f64, f: Option<f64>) -> f64 {
    match f {
        None => x,
        Some(f) => x + f,
    }
}

/// The caller's dynamic-load closure, as the Newton core borrows it.
type DynamicFn<'a> = &'a mut dyn FnMut(&FleetState) -> Result<DynamicLoad>;

/// The Newton core shared by [`solve_equilibrium_with`] and
/// [`solve_equilibrium_dynamic_with`]. With `dynamic`, every iteration adds
/// the closure's load — scaled to volume-metre units by `1/(ρg)` — to the
/// hydrostatic residuals: `(∇ − W/ρ + F_z/ρg, M_x − lcg·∇ + M/ρg)`. The
/// Jacobian stays the analytic waterplane one, i.e. the dynamic load is
/// treated as a slowly varying forcing (quasi-Newton). The closure is
/// expensive (a resistance evaluation), so it runs exactly once per
/// iteration, on the fleet already situated for that iteration, and never on
/// a dry fleet. `warm_start` seeds `(sinkage, trim)` and skips the coarse
/// phase.
fn equilibrium_core(
    mut situate: impl FnMut(f64, f64, bool) -> Result<FleetState>,
    mut dynamic: Option<DynamicFn<'_>>,
    load: &LoadCase,
    density: f64,
    gravity: f64,
    warm_start: Option<(f64, f64)>,
) -> Result<Solved> {
    if !(load.mass.is_finite() && load.mass > 0.0) {
        return Err(Error::InvalidConditions(format!(
            "load mass must be finite and positive, got {}",
            load.mass
        )));
    }
    if !(density.is_finite() && density > 0.0) {
        return Err(Error::InvalidConditions("density must be positive".into()));
    }
    if !(gravity.is_finite() && gravity > 0.0) {
        return Err(Error::InvalidConditions("gravity must be positive".into()));
    }
    if let Some(l) = load.lcg {
        if !l.is_finite() {
            return Err(Error::InvalidConditions("lcg must be finite".into()));
        }
    }
    if let Some((s0, tau0)) = warm_start {
        if !(s0.is_finite() && tau0.is_finite()) {
            return Err(Error::InvalidConditions(format!(
                "warm start (sinkage {s0}, trim {tau0}) must be finite"
            )));
        }
    }
    let v_target = load.mass / density;
    let z_guess = v_target.cbrt().max(1e-3);
    let pivot_x = load.lcg.unwrap_or(0.0);
    // The dynamic load in the residuals' units: volume and volume-metres.
    let per_rho_g = 1.0 / (density * gravity);
    let to_forcing = |d: DynamicLoad| (d.force_up * per_rho_g, d.moment_bow_up * per_rho_g);

    let (mut s, mut tau) = warm_start.unwrap_or((0.0, 0.0));
    let mut iterations = 0usize;
    // The latest dynamic evaluation with the (state, phase) it was made at,
    // so the final report can reuse it rather than pay for another.
    let mut last_dynamic: Option<((f64, f64, bool), DynamicLoad)> = None;

    // (coarse sampling, iteration cap, volume tolerance) per phase. A warm
    // start is presumed close (the previous speed's solution) and goes
    // straight to the fine phase.
    const PHASES: [(bool, usize, f64); 2] = [(true, 30, 1e-3), (false, 20, 2e-4)];
    let phases = if warm_start.is_some() {
        &PHASES[1..]
    } else {
        &PHASES[..]
    };
    for (phase_idx, &(coarse, max_iters, tol_v)) in phases.iter().enumerate() {
        let mut converged = false;
        // Adaptive relaxation: the waterplane-property Jacobian can
        // underestimate the true sensitivity (e.g. flare or structure
        // entering the water), which turns plain Newton into a limit cycle.
        // Halve the step whenever the volume residual flips sign, recover
        // gently while it doesn't.
        let mut relax = 1.0f64;
        let mut last_sign = 0.0f64;
        // Every phase but the very first, from-scratch one starts from a
        // state that was converged *somewhere else* — the coarse phase here
        // (a coarsened loft cannot resolve fine stern detail, a transom or a
        // chine, the way the full-resolution one does, so `situate` can hand
        // back a visibly different fleet at the very same `(s, tau)`), or a
        // warm start from a *different* speed's solution, whose dynamic
        // force can be a different scale entirely. Either way the state
        // hasn't actually been validated against what this phase will now
        // evaluate, so its first Newton step is speculative; damped below
        // once (only when the dynamic load turns out genuinely nonzero) so
        // it can't overshoot correcting for what is really a model or
        // operating-point shift rather than a residual to chase. Ordinary
        // within-phase oscillation detection is untouched and recovers full
        // speed within a couple more iterations regardless.
        let is_handoff_phase = phase_idx > 0 || warm_start.is_some();
        let mut first_iter_of_phase = true;
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
            let forcing = match dynamic.as_mut() {
                None => None,
                Some(d) => {
                    let dl = d(&fleet)?;
                    if !(dl.force_up.is_finite() && dl.moment_bow_up.is_finite()) {
                        return Err(Error::InvalidConditions(format!(
                            "dynamic load must be finite, got {dl:?} at sinkage {s}, trim {tau}"
                        )));
                    }
                    last_dynamic = Some(((s, tau, coarse), dl));
                    Some(to_forcing(dl))
                }
            };
            let t = totals(&fleet);
            let z_scale = t.draft.max(z_guess);
            let l_scale = fleet
                .members
                .iter()
                .map(|(h, _)| h.length())
                .fold(0.0f64, f64::max)
                .max(1e-6);
            let r1 = forced(t.volume - v_target, forcing.map(|f| f.0));
            let lcb = t.moment_x / t.volume;
            // The moment imbalance as an LCB shift: what the trim converges on.
            let lcb_err = |lcg: f64| forced(lcb - lcg, forcing.map(|f| f.1 / t.volume));
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
                    let r2 = forced(t.moment_x - lcg * t.volume, forcing.map(|f| f.1));
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
            // The handoff damping described above: only this phase's first
            // iteration, only when there is a genuinely nonzero dynamic load
            // right now. Bit-for-bit unaffected when `dynamic` is `None`, or
            // returns exactly zero (e.g. a caller probing the hydrostatic
            // path through the dynamic API) — `handoff_damp` is then always
            // 1.0, on every phase, warm-started or not.
            let handoff_damp = if first_iter_of_phase
                && is_handoff_phase
                && forcing.is_some_and(|(f, m)| f != 0.0 || m != 0.0)
            {
                0.3
            } else {
                1.0
            };
            first_iter_of_phase = false;
            ds = (ds * relax * handoff_damp).clamp(-0.3 * z_scale, 0.3 * z_scale);
            dtau = (dtau * relax * handoff_damp).clamp(-0.05, 0.05);
            let done_v = r1.abs() <= tol_v * v_target;
            let done_m = match load.lcg {
                None => true,
                Some(lcg) => lcb_err(lcg).abs() <= 1e-4 * l_scale,
            };
            let metric = (r1.abs() / (tol_v * v_target)).max(match load.lcg {
                None => 0.0,
                Some(lcg) => lcb_err(lcg).abs() / (1e-4 * l_scale),
            });
            if metric < best.0 {
                best = (metric, s, tau);
            }
            if std::env::var("MICHELL_DEBUG_FLOAT").is_ok() {
                eprintln!(
                    "DBG iter={iterations} s={s:.6} tau={tau:.6} V={:.6} lcb={lcb:.6} \
                     Aw={:.4} Mw={:.4} Iw={:.4} ds={ds:.6} dtau={dtau:.6} dyn={forcing:?}",
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
    let dynamic_load = match dynamic.as_mut() {
        None => DynamicLoad::default(),
        Some(_) if fleet.members.is_empty() => DynamicLoad::default(),
        Some(d) => match last_dynamic {
            // Converged: the last iteration evaluated the closure on this
            // very state at full resolution, so the fleet — and the load —
            // are the same. Only a near-miss fallback moves the state and
            // has to pay for one more evaluation.
            Some(((ls, lt, false), dl)) if ls == s && lt == tau => dl,
            _ => d(&fleet)?,
        },
    };
    let forcing = dynamic.is_some().then(|| to_forcing(dynamic_load));
    let volume_residual = forced(t.volume - v_target, forcing.map(|f| f.0)).abs() / v_target;
    let lcb_residual = load.lcg.map_or(0.0, |l| {
        let shift = match forcing {
            Some((_, fm)) if t.volume > 0.0 => Some(fm / t.volume),
            _ => None,
        };
        forced(lcb - l, shift).abs()
    });
    Ok(Solved {
        sinkage: s,
        trim: tau,
        fleet,
        totals: t,
        lcb,
        dynamic: dynamic_load,
        iterations,
        volume_residual,
        lcb_residual,
    })
}

/// Generic equilibrium core. `situate(sinkage, trim, coarse)` produces the
/// fleet at a platform state (coarse = reduced sampling for iterations).
pub fn solve_equilibrium_with(
    situate: impl FnMut(f64, f64, bool) -> Result<FleetState>,
    load: &LoadCase,
    density: f64,
) -> Result<Equilibrium> {
    // Gravity cancels from a purely hydrostatic balance; any value serves.
    let sol = equilibrium_core(situate, None, load, density, STANDARD_GRAVITY, None)?;
    Ok(Equilibrium {
        sinkage: sol.sinkage,
        trim: sol.trim,
        volume: sol.totals.volume,
        lcb: sol.lcb,
        waterplane_area: sol.totals.wp_area,
        iterations: sol.iterations,
        volume_residual: sol.volume_residual,
        lcb_residual: sol.lcb_residual,
        fleet: sol.fleet,
    })
}

/// Equilibrium at speed: [`solve_equilibrium_with`] with a hydrodynamic
/// vertical force and pitch moment in the balance — thin-ship dynamic
/// sinkage and trim. `dynamic(&fleet)` returns the [`DynamicLoad`] on the
/// fleet as situated at the state under test, and the solve satisfies
///
/// ```text
/// ρg(∇ − W/ρ)     + F_z = 0
/// ρg(M_x − lcg·∇) + M   = 0        M_x = ∫ x d∇
/// ```
///
/// so a suction (`force_up < 0`) sinks the hull below its hydrostatic
/// displacement and a bow-up moment trims it bow-up. The closure runs once
/// per iteration (never more), on the coarsened fleets of the first phase as
/// well as the fine ones, and the reported `dynamic` is its value on the
/// returned fleet. `warm_start` seeds `(sinkage, trim)` — the previous
/// speed's solution, when sweeping — and skips the coarse phase. Tolerances
/// are the hydrostatic solver's; the dynamic terms are carried as a forcing
/// against the waterplane Jacobian, so a load that varies strongly with
/// attitude simply takes more iterations. A force past ~30 % of the weight
/// is solved like any other and shows up in `lift_fraction`.
pub fn solve_equilibrium_dynamic_with(
    situate: impl FnMut(f64, f64, bool) -> Result<FleetState>,
    mut dynamic: impl FnMut(&FleetState) -> Result<DynamicLoad>,
    load: &LoadCase,
    density: f64,
    gravity: f64,
    warm_start: Option<(f64, f64)>,
) -> Result<DynamicEquilibrium> {
    let sol = equilibrium_core(
        situate,
        Some(&mut dynamic),
        load,
        density,
        gravity,
        warm_start,
    )?;
    Ok(DynamicEquilibrium {
        sinkage: sol.sinkage,
        trim: sol.trim,
        volume: sol.totals.volume,
        lcb: sol.lcb,
        dynamic: sol.dynamic,
        lift_fraction: sol.dynamic.force_up / (load.mass * gravity),
        iterations: sol.iterations,
        volume_residual: sol.volume_residual,
        lcb_residual: sol.lcb_residual,
        fleet: sol.fleet,
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

/// One pose per body, or the mismatch as an error.
fn check_poses(bodies: &[&Body], poses: &[HullPose]) -> Result<()> {
    if bodies.len() != poses.len() {
        return Err(Error::InvalidInput(format!(
            "{} poses supplied for {} bodies",
            poses.len(),
            bodies.len()
        )));
    }
    Ok(())
}

/// The reduced sampling the Newton iterations run on: about half the
/// stations and waterlines (never below 31 × 11) over a capped control net.
fn coarse_body_options(opts: &BodyOptions) -> BodyOptions {
    let mut coarse_opts = *opts;
    coarse_opts.stations = (opts.stations / 2).clamp(31, opts.stations.max(31));
    coarse_opts.waterlines = (opts.waterlines / 2).clamp(11, opts.waterlines.max(11));
    coarse_opts.fit.n_ctrl_x = opts.fit.n_ctrl_x.min(10).max(opts.fit.degree_x + 1);
    coarse_opts.fit.n_ctrl_z = opts.fit.n_ctrl_z.min(7).max(opts.fit.degree_z + 1);
    coarse_opts
}

/// The `situate` closure over a body fleet that the equilibrium solvers run
/// on: every body at its pose under the platform state, the dry ones counted
/// rather than lofted; coarsened options for the iterations, the caller's
/// for the polish.
fn body_situator<'a>(
    bodies: &'a [&'a Body],
    water_offset: f64,
    poses: &'a [HullPose],
    pivot_x: f64,
    opts: &'a BodyOptions,
) -> impl FnMut(f64, f64, bool) -> Result<FleetState> + 'a {
    let coarse_opts = coarse_body_options(opts);
    move |s, tau, coarse| {
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
    }
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
    check_poses(bodies, poses)?;
    let pivot_x = load.lcg.unwrap_or(0.0);
    solve_equilibrium_with(
        body_situator(bodies, water_offset, poses, pivot_x, opts),
        load,
        density,
    )
}

/// [`solve_equilibrium_bodies`] at speed: the same body fleet balanced
/// against the caller's [`DynamicLoad`] (see
/// [`solve_equilibrium_dynamic_with`] for the balance, the closure contract,
/// and `warm_start`).
#[allow(clippy::too_many_arguments)]
pub fn solve_equilibrium_bodies_dynamic(
    bodies: &[&Body],
    water_offset: f64,
    poses: &[HullPose],
    load: &LoadCase,
    density: f64,
    gravity: f64,
    opts: &BodyOptions,
    dynamic: impl FnMut(&FleetState) -> Result<DynamicLoad>,
    warm_start: Option<(f64, f64)>,
) -> Result<DynamicEquilibrium> {
    check_poses(bodies, poses)?;
    let pivot_x = load.lcg.unwrap_or(0.0);
    solve_equilibrium_dynamic_with(
        body_situator(bodies, water_offset, poses, pivot_x, opts),
        dynamic,
        load,
        density,
        gravity,
        warm_start,
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
    check_poses(bodies, poses)?;
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
    let eq = solve_equilibrium_with(
        body_situator(bodies, water_offset, &heeled, pivot_x, opts),
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

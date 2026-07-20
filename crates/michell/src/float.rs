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

use crate::error::{Error, Result};
use crate::iges::{HullPose, ImportOptions, Platform, SituatedFleet, SourceFleet};

/// What the platform must carry.
#[derive(Debug, Clone, Copy)]
pub struct LoadCase {
    /// Total mass [kg].
    pub mass: f64,
    /// Longitudinal centre of gravity [m], in the file's x coordinates.
    /// `None` locks the platform pitch at zero and balances weight only.
    pub lcg: Option<f64>,
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
    pub fleet: SituatedFleet,
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

fn totals(fleet: &SituatedFleet) -> Totals {
    let mut t = Totals {
        volume: 0.0,
        moment_x: 0.0,
        wp_area: 0.0,
        wp_moment: 0.0,
        wp_second: 0.0,
        draft: 0.0,
    };
    for m in &fleet.members {
        let v = m.hull.displaced_volume();
        t.volume += v;
        t.moment_x += m.hull.lcb_x() * v;
        t.wp_area += m.hull.waterplane_area();
        t.wp_moment += m.hull.waterplane_moment();
        t.wp_second += m.hull.waterplane_second_moment();
        t.draft = t.draft.max(m.hull.draft());
    }
    t
}

/// Solve for the platform sinkage (and pitch, when `lcg` is given) that
/// floats `load` at the given design poses. `waterline_z` is the base
/// waterline the sinkage is measured from.
pub fn solve_equilibrium(
    src: &SourceFleet,
    waterline_z: f64,
    poses: &[HullPose],
    load: &LoadCase,
    density: f64,
    opts: &ImportOptions,
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

    // Coarse options for the iterations; the final answer is re-situated at
    // the caller's resolution.
    let mut coarse = *opts;
    coarse.stations = (opts.stations / 2).clamp(31, opts.stations.max(31));
    coarse.waterlines = (opts.waterlines / 2).clamp(11, opts.waterlines.max(11));
    coarse.fit.n_ctrl_x = opts.fit.n_ctrl_x.min(10).max(opts.fit.degree_x + 1);
    coarse.fit.n_ctrl_z = opts.fit.n_ctrl_z.min(7).max(opts.fit.degree_z + 1);

    let mut s = 0.0f64;
    let mut tau = 0.0f64;
    let mut iterations = 0usize;

    for (phase_opts, max_iters, tol_v) in [(&coarse, 25usize, 1e-3f64), (opts, 8, 2e-4)] {
        let mut converged = false;
        for _ in 0..max_iters {
            iterations += 1;
            let fleet = src.situate(
                waterline_z,
                poses,
                &Platform {
                    sinkage: s,
                    trim: tau,
                    pivot_x,
                },
                phase_opts,
            )?;
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
                .map(|m| m.hull.length())
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
            // Damping.
            ds = ds.clamp(-0.3 * z_scale, 0.3 * z_scale);
            dtau = dtau.clamp(-0.05, 0.05);
            let done_v = r1.abs() <= tol_v * v_target;
            let done_m = match load.lcg {
                None => true,
                Some(lcg) => (lcb - lcg).abs() <= 1e-4 * l_scale,
            };
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
            return Err(Error::InvalidConditions(format!(
                "equilibrium did not converge after {iterations} iterations \
                 (mass {} kg, lcg {:?}); the load may be outside what the \
                 geometry can float",
                load.mass, load.lcg
            )));
        }
    }

    // Final state at full resolution.
    let fleet = src.situate(
        waterline_z,
        poses,
        &Platform {
            sinkage: s,
            trim: tau,
            pivot_x,
        },
        opts,
    )?;
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

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

use crate::conditions::STANDARD_GRAVITY;
use crate::error::{Error, Result};
use crate::iges::{HullPose, Platform, SectionalOptions, SectionalState};
use crate::michell::Placement;
use crate::sectional::SectionalHull;
use crate::source::SourceHull;

/// What the platform must carry.
#[derive(Debug, Clone, Copy)]
pub struct LoadCase {
    /// Total mass [kg].
    pub mass: f64,
    /// Longitudinal centre of gravity [m], in the fleet's x coordinates.
    /// `None` locks the platform pitch at zero and balances weight only.
    pub lcg: Option<f64>,
}

/// An extra point mass mounted on a hull, positioned as offsets from the hull's
/// **centerpoint** (midship station, design floatplane): `dx` forward, `dz`
/// **down** (deeper) — the same axis conventions as [`HullPose`]. It rides with
/// the hull through its pose just like the hull's own CG, and adds to the
/// mass-weighted fleet CG ([`fleet_cg`]). Its height matters only through
/// design trim, which swings a raised mass fore or aft.
#[derive(Debug, Clone, Copy, Default)]
pub struct PointLoad {
    /// Point mass [kg].
    pub mass: f64,
    /// Longitudinal offset [m] from the hull midship (+forward).
    pub dx: f64,
    /// Vertical offset [m] from the design floatplane, positive **down**.
    pub dz: f64,
}

/// A single hull's contribution to the load, with its centre of gravity given
/// in the hull's own **local** design frame — the same frame [`HullPose`] maps
/// into the platform. `lcg` is the local longitudinal CG station; `vcg` is
/// metres **above the design floatplane**. `points` are extra discrete masses
/// mounted on the hull (batteries, crew, ballast). Because every contribution
/// is expressed pre-pose, mounting the hull carries its weight with it:
/// `dx`/`dz` translate each CG, design `trim` rotates it. The fleet CG is then
/// always the mass-weighted sum over all hull loads and point loads
/// ([`fleet_cg`]) — never specified directly — so it moves realistically as
/// the hulls are repositioned.
#[derive(Debug, Clone, Default)]
pub struct HullLoad {
    /// Structural mass carried on this hull [kg], at `(lcg, vcg)`.
    pub mass: f64,
    /// Longitudinal CG [m] in the hull's local x (before the pose's dx/trim).
    pub lcg: f64,
    /// Vertical CG [m] above the design floatplane (before the pose's dz).
    pub vcg: f64,
    /// Extra point masses mounted on the hull.
    pub points: Vec<PointLoad>,
}

/// The platform centre of gravity, mass-weighted over the per-hull loads at
/// their posed positions. `lcg` is in the platform (fleet) frame; `vcg` is
/// metres above the design floatplane. Zero throughout for a massless fleet.
#[derive(Debug, Clone, Copy, Default)]
pub struct FleetCg {
    /// Total mass [kg] = Σ hull mass.
    pub mass: f64,
    /// Longitudinal CG [m] in the fleet frame.
    pub lcg: f64,
    /// Vertical CG [m] above the design floatplane.
    pub vcg: f64,
}

/// Derive the fleet CG by summing every per-hull load — the hull's structural
/// CG and each mounted [`PointLoad`] — carried through the hull's pose. Each
/// contribution is a local point (longitudinal `x`, height `vcg` above the
/// floatplane) mapped through the design pose exactly as the geometry is —
/// design `trim` about the pose pivot, then `+dx`, and deepened by `+dz` —
/// then mass-weighted. Massless contributions (and, if the whole fleet is
/// massless, the fleet) add nothing. `midships` (each hull's x mid, the
/// default trim pivot and point-load origin), `loads`, and `poses` must share
/// their length and order.
pub fn fleet_cg(midships: &[f64], loads: &[HullLoad], poses: &[HullPose]) -> FleetCg {
    let mut mass = 0.0;
    let mut mx = 0.0;
    let mut mz = 0.0; // Σ m · (height above floatplane), up-positive.
    for ((&midship, load), pose) in midships.iter().zip(loads).zip(poses) {
        let pivot = pose.pivot_x.unwrap_or(midship);
        // Every contribution as (mass, local x, height above the floatplane).
        // The structural load sits at (lcg, vcg); a point load at
        // (midship+dx, −dz) — dz is +down.
        let structural = std::iter::once((load.mass, load.lcg, load.vcg));
        let points = load.points.iter().map(|p| (p.mass, midship + p.dx, -p.dz));
        for (m, lx, vcg) in structural.chain(points) {
            if !(m.is_finite() && m > 0.0) {
                continue;
            }
            // The design-pose map (z positive DOWN, as the pose's dz): the
            // point sits at down-coord z = −vcg. Rotate by design trim about
            // (pivot, 0), then shift +dx / +dz.
            let mut x = lx;
            let mut zd = -vcg;
            if pose.trim != 0.0 {
                let (s, c) = pose.trim.sin_cos();
                let (rx, rz) = (x - pivot, zd);
                x = pivot + rx * c + rz * s;
                zd = rz * c - rx * s;
            }
            x += pose.dx;
            zd += pose.dz;
            mass += m;
            mx += m * x;
            mz += m * (-zd); // back to up-positive height above floatplane.
        }
    }
    if mass > 0.0 {
        FleetCg {
            mass,
            lcg: mx / mass,
            vcg: mz / mass,
        }
    } else {
        FleetCg::default()
    }
}

/// What the equilibrium solver reads from a wetted hull: its buoyancy and
/// its waterplane (in the hull's own x coordinates), as a [`SectionalHull`]
/// provides it.
pub trait Floating {
    fn displaced_volume(&self) -> f64;
    fn lcb_x(&self) -> f64;
    fn waterplane_area(&self) -> f64;
    fn waterplane_moment(&self) -> f64;
    fn waterplane_second_moment(&self) -> f64;
    fn draft(&self) -> f64;
    fn length(&self) -> f64;
}

macro_rules! floating {
    ($t:ty) => {
        impl Floating for $t {
            fn displaced_volume(&self) -> f64 {
                <$t>::displaced_volume(self)
            }
            fn lcb_x(&self) -> f64 {
                <$t>::lcb_x(self)
            }
            fn waterplane_area(&self) -> f64 {
                <$t>::waterplane_area(self)
            }
            fn waterplane_moment(&self) -> f64 {
                <$t>::waterplane_moment(self)
            }
            fn waterplane_second_moment(&self) -> f64 {
                <$t>::waterplane_second_moment(self)
            }
            fn draft(&self) -> f64 {
                <$t>::draft(self)
            }
            fn length(&self) -> f64 {
                <$t>::length(self)
            }
        }
    };
}
floating!(SectionalHull);

/// The wetted fleet at some state: what the solver (and resistance) consume.
#[derive(Debug)]
pub struct FleetState<M = SectionalHull> {
    pub members: Vec<(M, Placement)>,
    /// Number of source hulls entirely above the water.
    pub dry: usize,
}

/// A solved floating condition.
#[derive(Debug)]
pub struct Equilibrium<M = SectionalHull> {
    /// Solved additional immersion of the platform [m] (relative to the base
    /// waterline; negative = riding higher).
    pub sinkage: f64,
    /// Solved platform pitch [rad]; positive raises the +x end.
    pub trim: f64,
    /// The fleet situated at the solution, at full requested resolution.
    pub fleet: FleetState<M>,
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

fn totals<M: Floating>(fleet: &FleetState<M>) -> Totals {
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
pub struct DynamicEquilibrium<M = SectionalHull> {
    /// Solved additional immersion of the platform [m] (relative to the base
    /// waterline; negative = riding higher).
    pub sinkage: f64,
    /// Solved platform pitch [rad]; positive raises the +x end.
    pub trim: f64,
    /// The fleet situated at the solution, at full requested resolution.
    pub fleet: FleetState<M>,
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
struct Solved<M> {
    sinkage: f64,
    trim: f64,
    fleet: FleetState<M>,
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
type DynamicFn<'a, M> = &'a mut dyn FnMut(&FleetState<M>) -> Result<DynamicLoad>;

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
fn equilibrium_core<M: Floating>(
    mut situate: impl FnMut(f64, f64, bool) -> Result<FleetState<M>>,
    mut dynamic: Option<DynamicFn<'_, M>>,
    load: &LoadCase,
    density: f64,
    gravity: f64,
    warm_start: Option<(f64, f64)>,
) -> Result<Solved<M>> {
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
        // (a coarsened cut cannot resolve fine stern detail, a transom or a
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
        // The situate → cut model carries small-scale roughness (~0.1% of
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
            let dl_base = match dynamic.as_mut() {
                None => None,
                Some(d) => {
                    let dl = d(&fleet)?;
                    if !(dl.force_up.is_finite() && dl.moment_bow_up.is_finite()) {
                        return Err(Error::InvalidConditions(format!(
                            "dynamic load must be finite, got {dl:?} at sinkage {s}, trim {tau}"
                        )));
                    }
                    last_dynamic = Some(((s, tau, coarse), dl));
                    Some(dl)
                }
            };
            let forcing = dl_base.map(to_forcing);
            let t = totals(&fleet);
            let z_scale = t.draft.max(z_guess);
            let l_scale = fleet
                .members
                .iter()
                .map(|(h, _)| h.length())
                .fold(0.0f64, f64::max)
                .max(1e-6);
            // The dynamic load's own sensitivity to (s, τ), by one-sided finite
            // difference against `dl_base`, falling back to the other side if
            // the perturbed fleet goes dry (only relevant right at the edge of
            // floating). The analytic Jacobian below (`t.wp_area` etc.) treats
            // the dynamic load as a *constant* added to the residual — exact
            // for the hydrostatic terms, but only zeroth-order for a force that
            // genuinely curves with attitude (a large transom's immersion
            // changing character, say). That mismatch is what turns plain
            // quasi-Newton into a limit cycle no matter how precisely the force
            // itself is resolved — confirmed on the motivating case by finding
            // loose- and tight-quadrature evaluations agreeing to <1% across
            // the very range the solver was oscillating in, which rules out
            // quadrature noise as the cause. `h_s`/`h_tau` sit comfortably
            // above that noise floor while staying well inside the existing
            // step-size clamps below.
            // Scoped to the same phases the handoff damping above covers —
            // never the initial, cold, far-from-solution coarse phase. That
            // phase already converges reliably on the pure hydrostatic
            // Jacobian alone (it always has; the mismatch this section
            // exists for only shows up once precision matters, in the fine
            // phase), and empirically, probing a finite difference from a
            // wild starting guess is actively harmful: on the motivating
            // case, applying it unconditionally sent the cold coarse phase
            // to a multi-metre "sinkage" and a trim past the 20° abort
            // limit, diverging outright where the plain analytic Jacobian
            // had always converged in under ten iterations.
            let dyn_jac = if is_handoff_phase { dl_base } else { None }.map(|dl| {
                let h_s = 0.02 * z_scale;
                let h_tau = 0.01f64;
                let mut probe = |ds: f64, dtau: f64| -> Result<Option<DynamicLoad>> {
                    let f = situate(s + ds, tau + dtau, coarse)?;
                    if f.members.is_empty() {
                        Ok(None)
                    } else {
                        Ok(Some(dynamic.as_mut().expect("dl_base implies a closure")(
                            &f,
                        )?))
                    }
                };
                let mut one_sided = |h: f64, along_s: bool| -> Result<Option<(f64, f64)>> {
                    let (fwd_ds, fwd_dt) = if along_s { (h, 0.0) } else { (0.0, h) };
                    let sample = match probe(fwd_ds, fwd_dt)? {
                        Some(v) => Some((h, v)),
                        None => probe(-fwd_ds, -fwd_dt)?.map(|v| (-h, v)),
                    };
                    Ok(sample.map(|(signed_h, v)| {
                        (
                            (v.force_up - dl.force_up) / signed_h,
                            (v.moment_bow_up - dl.moment_bow_up) / signed_h,
                        )
                    }))
                };
                let ds_slope = one_sided(h_s, true)?;
                // Trim isn't being solved without an lcg, so skip those evals.
                let dt_slope = match load.lcg {
                    Some(_) => one_sided(h_tau, false)?,
                    None => None,
                };
                Ok::<_, Error>((ds_slope, dt_slope))
            });
            let dyn_jac = dyn_jac.transpose()?;
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
            // Fold the dynamic load's own sensitivity into the hydrostatic
            // Jacobian, in the same volume/volume-metre units `forced` already
            // uses (i.e. scaled by `per_rho_g`). `None` (no closure, a zero
            // load, or a dry perturbed fleet on that axis) leaves the pure
            // hydrostatic term untouched, so a hydrostatic-only solve — and a
            // dynamic one whose load has no local sensitivity to reach for —
            // is completely unaffected.
            let (dfz_ds, dm_ds_dyn) = dyn_jac
                .and_then(|(ds_slope, _)| ds_slope)
                .map_or((0.0, 0.0), |(dfz, dm)| (per_rho_g * dfz, per_rho_g * dm));
            let (dfz_dt, dm_dt_dyn) = dyn_jac
                .and_then(|(_, dt_slope)| dt_slope)
                .map_or((0.0, 0.0), |(dfz, dm)| (per_rho_g * dfz, per_rho_g * dm));
            // d/ds and d/dτ of (V, M); rotation about (lcg, waterline), plus
            // the dynamic terms above.
            let dv_ds = t.wp_area + dfz_ds;
            let dv_dt = -(t.wp_moment - pivot_x * t.wp_area) + dfz_dt;
            let (mut ds, mut dtau) = match load.lcg {
                None => {
                    if dv_ds > 1e-12 * t.wp_area {
                        (r1 / dv_ds, 0.0)
                    } else {
                        // The dynamic sensitivity overwhelmed the waterplane
                        // area (a very steep local force gradient): fall back
                        // to the pure hydrostatic step rather than divide by a
                        // near-zero or negative denominator.
                        (r1 / t.wp_area, 0.0)
                    }
                }
                Some(lcg) => {
                    let r2 = forced(t.moment_x - lcg * t.volume, forcing.map(|f| f.1));
                    let dm_ds = t.wp_moment + dm_ds_dyn;
                    let dm_dt = -(t.wp_second - pivot_x * t.wp_moment) + dm_dt_dyn;
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
pub fn solve_equilibrium_with<M: Floating>(
    situate: impl FnMut(f64, f64, bool) -> Result<FleetState<M>>,
    load: &LoadCase,
    density: f64,
) -> Result<Equilibrium<M>> {
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
pub fn solve_equilibrium_dynamic_with<M: Floating>(
    situate: impl FnMut(f64, f64, bool) -> Result<FleetState<M>>,
    mut dynamic: impl FnMut(&FleetState<M>) -> Result<DynamicLoad>,
    load: &LoadCase,
    density: f64,
    gravity: f64,
    warm_start: Option<(f64, f64)>,
) -> Result<DynamicEquilibrium<M>> {
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

/// The `situate` closure for hulls cut into sections from their source
/// geometry: each of `hulls` re-cut at the platform state, warm-started from
/// its last cut. The coarse phase of the solve runs on half the stations and
/// rays.
pub fn sectional_situator<'a>(
    hulls: &'a [SourceHull<'a>],
    pivot_x: f64,
    opts: &'a SectionalOptions,
) -> impl FnMut(f64, f64, bool) -> Result<FleetState<SectionalHull>> + 'a {
    let coarse_opts = SectionalOptions {
        stations: (opts.stations / 2).max(31),
        rays: (opts.rays / 2).max(17),
        ..*opts
    };
    let mut states = vec![(SectionalState::default(), SectionalState::default()); hulls.len()];
    move |s, tau, coarse| {
        let platform = Platform {
            sinkage: s,
            trim: tau,
            pivot_x,
        };
        let o = if coarse { &coarse_opts } else { opts };
        let mut members = Vec::new();
        let mut dry = 0usize;
        for (h, st) in hulls.iter().zip(states.iter_mut()) {
            let state = if coarse { &mut st.0 } else { &mut st.1 };
            match h.source.situate_sectional_warm(
                h.index,
                h.waterline_z,
                &h.pose,
                &platform,
                o,
                state,
            )? {
                Some(h) => members.push((h.hull, h.placement)),
                None => dry += 1,
            }
        }
        Ok(FleetState { members, dry })
    }
}

/// Equilibrium of hulls cut into sections (see [`sectional_situator`]); the
/// platform trims about `load.lcg` (0 without one).
pub fn solve_equilibrium_sectional(
    hulls: &[SourceHull],
    load: &LoadCase,
    density: f64,
    opts: &SectionalOptions,
) -> Result<Equilibrium<SectionalHull>> {
    let pivot_x = load.lcg.unwrap_or(0.0);
    solve_equilibrium_with(sectional_situator(hulls, pivot_x, opts), load, density)
}

/// [`solve_equilibrium_sectional`] at speed, against the caller's
/// [`DynamicLoad`] (for thin-ship sinkage and trim,
/// [`crate::sectional::dynamic_load_closure`]); see
/// [`solve_equilibrium_dynamic_with`] for the balance and `warm_start`.
#[allow(clippy::too_many_arguments)]
pub fn solve_equilibrium_sectional_dynamic(
    hulls: &[SourceHull],
    load: &LoadCase,
    density: f64,
    gravity: f64,
    opts: &SectionalOptions,
    dynamic: impl FnMut(&FleetState<SectionalHull>) -> Result<DynamicLoad>,
    warm_start: Option<(f64, f64)>,
) -> Result<DynamicEquilibrium<SectionalHull>> {
    let pivot_x = load.lcg.unwrap_or(0.0);
    solve_equilibrium_dynamic_with(
        sectional_situator(hulls, pivot_x, opts),
        dynamic,
        load,
        density,
        gravity,
        warm_start,
    )
}

#[cfg(test)]
mod sectional_tests {
    use super::*;
    use crate::iges::SourceFleet;

    fn e12() -> Option<(SourceFleet, f64, usize)> {
        let text = crate::cad_fixture("e12.igs")?;
        let wl = -0.95;
        let src = crate::iges::source_fleet(&text, wl).unwrap();
        let idx = (0..src.len()).next().unwrap();
        Some((src, wl, idx))
    }

    /// Loaded to its own displacement and LCB at the design waterline, a
    /// sectional hull floats there; 5% heavier, it sinks by about the extra
    /// volume over the waterplane area.
    #[test]
    fn sectional_equilibrium_recovers_the_design_waterline() {
        let Some((src, wl, idx)) = e12() else {
            return;
        };
        let opts = SectionalOptions {
            waterline_z: wl,
            ..Default::default()
        };
        let design = src
            .situate_sectional(idx, wl, &HullPose::default(), &Platform::default(), &opts)
            .unwrap()
            .unwrap()
            .hull;
        let rho = 1025.0;
        let hulls = [SourceHull {
            source: &src,
            index: idx,
            waterline_z: wl,
            pose: HullPose::default(),
        }];
        let load = LoadCase {
            mass: rho * design.displaced_volume(),
            lcg: Some(design.lcb_x()),
        };
        let eq = solve_equilibrium_sectional(&hulls, &load, rho, &opts).unwrap();
        assert!(eq.sinkage.abs() < 1e-5, "sinkage {}", eq.sinkage);
        assert!(eq.trim.abs() < 1e-5, "trim {}", eq.trim);
        let heavy = LoadCase {
            mass: 1.05 * load.mass,
            ..load
        };
        let eq = solve_equilibrium_sectional(&hulls, &heavy, rho, &opts).unwrap();
        let expect = 0.05 * design.displaced_volume() / design.waterplane_area();
        assert!(
            (eq.sinkage - expect).abs() < 0.1 * expect,
            "sinkage {} vs ~{expect}",
            eq.sinkage
        );
        assert!(eq.volume_residual < 1e-6, "{}", eq.volume_residual);
    }
}

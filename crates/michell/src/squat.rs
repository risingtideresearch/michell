//! Thin-ship **dynamic sinkage and trim**: the vertical force and pitching
//! moment the steady near-field pressure exerts on a hull (or fleet) held at
//! its hydrostatic attitude, from the same exact per-span transforms the wave
//! resistance is built on.
//!
//! ## Where the force lives in the spectrum
//!
//! Write the Kelvin source in 2-D Fourier form. With the stream toward −x and
//! `z` down, the free-surface condition fixes the image amplitude
//!
//! ```text
//! A(k_x, k) = (k_x² + νk)/(k_x² − νk),      ν = g/U²,
//! ```
//!
//! which is −1 (rigid wall) for long longitudinal modes and +1 (free surface,
//! φ = 0) for short ones, with the steady-wave dispersion curve `k = k_x²/ν`
//! as the pole between them. The linearised pressure `p = ρUφ_x` integrated
//! over the hull — pressure on the sloped bottom, or equivalently the volume
//! form plus the **waterplane face** `∬_WP p dA` (the "defect of vertical
//! pressure", Havelock 1939) — reduces, after integration by parts along the
//! hull, to wavenumber integrals of the near-field transforms
//! [`SquatTransforms`]:
//!
//! ```text
//! F_up  = −(ρU²/2π²) PV∬ d²k [ A·Re(q̄q) − ((A−1)/k)·Re(w̄q) ]
//! M_pv  = −(ρU²/2π²) PV∬ d²k [ A·Re(r̄q) − ((A−1)/k)·Re(r̄_wl q) ]
//! M_res = −(4ρU²ν/π) ∫₁^∞ dλ λ/√(λ²−1) · Im[ νλ² r̄q − r̄_wl q ]   on k = k_x²/ν
//! ```
//!
//! with `r = p + q1 + (x_c − x_ref) q` the transform of `∂_x[(x − x_ref) f]`
//! and the `_wl` quantities their waterline (z = 0) counterparts. Two facts
//! shape the split. The radiation condition puts a half-residue
//! `−iπ sgn(k_x) δ(k_x² − νk)` alongside each principal value; it is odd in
//! `k_x`, so it drops out of the **force** (whose integrand is even) and
//! survives only in the **moment** — sinkage is a local-field effect, trim is
//! the fore-aft asymmetry of the wave pattern. And the unbounded-fluid Rankine
//! part of the source (`−1/r`) contributes no force at all (its kernel is odd
//! in `z − ζ` against a symmetric pair of weights: d'Alembert), but does
//! contribute a moment for a fore-aft asymmetric hull — the linearised Munk
//! moment — which needs ordered z-moments and is added separately.
//!
//! For a fore-aft symmetric hull `M_pv` vanishes by parity and the trim is
//! pure wave effect; the Wigley trim's sign change near Fn 0.35 is exactly this
//! term and is one of the tests below.
//!
//! ## What the numbers mean
//!
//! Sinkage is `−F_up/(ρ g A_w)` at leading order, the moment trims against
//! the waterplane's longitudinal inertia; both are handed to the equilibrium
//! solver as speed-dependent load terms so the platform floats at its dynamic
//! attitude. Thin-ship theory overstates both by 20–40% against surface-panel
//! (Neumann–Michell) linear theory at Fn 0.3–0.4 for B/L ≈ 0.1 — the same
//! overshoot it shows on R_w — and the whole linearisation expires once
//! dynamic lift carries a real share of the weight (volumetric Froude ≳ 3);
//! [`DynamicForce::lift_fraction`] is reported so that boundary is visible.

use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::hull::Hull;
use crate::michell::{InnerIntegral, Placement, SquatTransforms, WaveOptions, ZContracted};
use crate::moments::C64;
use crate::quadrature::gauss_legendre;
use std::f64::consts::{FRAC_PI_2, PI};

/// Near-field vertical force and pitching moment on a fleet.
#[derive(Debug, Clone, Copy)]
pub struct DynamicForce {
    /// Vertical hydrodynamic force [N], positive **upward** (negative is
    /// suction: the hull sinks).
    pub force_up: f64,
    /// Pitching moment [N·m] about `x_ref` at the waterline, positive
    /// **bow-up**.
    pub moment_bow_up: f64,
    /// The moment's local (principal-value) share, for diagnostics.
    pub moment_local: f64,
    /// The moment's wave (residue) share.
    pub moment_wave: f64,
    /// `force_up / (ρ g ∇)`: the dynamic force as a fraction of the buoyancy
    /// at this attitude. Negative means the hull is being pulled down.
    pub lift_fraction: f64,
    /// Estimated relative quadrature error on the force (last two refinement
    /// passes).
    pub est_rel_error: f64,
    /// Transform evaluations spent.
    pub evaluations: usize,
}

/// Quadrature controls for the near-field integrals.
#[derive(Debug, Clone, Copy)]
pub struct SquatOptions {
    /// Relative tolerance on the force between refinement passes.
    pub rel_tol: f64,
    /// Maximum panel-doubling passes.
    pub max_refinements: usize,
    /// Wave options — the transom closure is taken from here so the body the
    /// force acts on is the same closed composite the resistance sees.
    pub wave: WaveOptions,
}

impl Default for SquatOptions {
    fn default() -> Self {
        SquatOptions {
            rel_tol: 1e-4,
            max_refinements: 3,
            wave: WaveOptions::default(),
        }
    }
}

/// A closure adapting [`multihull_dynamic_force`] to the shape
/// `FnMut(&FleetState) -> Result<DynamicLoad>` the dynamic equilibrium solver
/// wants (`crate::float::solve_equilibrium_dynamic_with` /
/// `solve_equilibrium_bodies_dynamic`): the fleet's situated hulls are
/// borrowed fresh from `FleetState` on every call, so the closure has no
/// stale state and works unchanged across Newton iterations and warm starts.
///
/// A dry fleet (`fleet.members` empty — everything lifted clear of the water)
/// reports zero load rather than erroring: the hydrostatic side of the solver
/// already handles that case by sinking until something gets wet.
pub fn dynamic_load_closure<'a>(
    cond: &'a Conditions,
    x_ref: f64,
    opts: &'a SquatOptions,
) -> impl FnMut(&crate::float::FleetState) -> Result<crate::float::DynamicLoad> + 'a {
    move |fleet: &crate::float::FleetState| {
        if fleet.members.is_empty() {
            return Ok(crate::float::DynamicLoad {
                force_up: 0.0,
                moment_bow_up: 0.0,
            });
        }
        let members: Vec<(&Hull, Placement)> = fleet.members.iter().map(|(h, p)| (h, *p)).collect();
        let d = multihull_dynamic_force(&members, cond, x_ref, opts)?;
        Ok(crate::float::DynamicLoad {
            force_up: d.force_up,
            moment_bow_up: d.moment_bow_up,
        })
    }
}

/// Near-field force and moment of a single hull about `x_ref`.
/// Near-field force and moment of a single hull about `x_ref`.
pub fn dynamic_force(
    hull: &Hull,
    cond: &Conditions,
    x_ref: f64,
    opts: &SquatOptions,
) -> Result<DynamicForce> {
    multihull_dynamic_force(&[(hull, Placement::default())], cond, x_ref, opts)
}

/// Near-field force and moment of a fleet — every member's pressure includes
/// what every other member's source sheet induces on it, through the same
/// placement phases the wave superposition uses (`e^{ik_x Δx} cos(k_y Δy)`),
/// so demihull interaction in sinkage and trim is carried exactly within the
/// theory. `x_ref` is the pitch pivot in fleet coordinates.
pub fn multihull_dynamic_force(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    x_ref: f64,
    opts: &SquatOptions,
) -> Result<DynamicForce> {
    cond.validate()?;
    if members.is_empty() {
        return Err(Error::InvalidGeometry("empty fleet".into()));
    }
    if !(opts.rel_tol.is_finite() && opts.rel_tol > 0.0) {
        return Err(Error::InvalidConditions("rel_tol must be positive".into()));
    }
    let u = cond.speed;
    let g = cond.gravity;
    let nu = g / (u * u);
    let rho = cond.fluid.density;

    let mut fleet = Fleet::new(members, nu, x_ref, opts.wave);

    // Length scales for the quadrature layout.
    let l_max = members
        .iter()
        .map(|(h, _)| h.length())
        .fold(0.0f64, f64::max);
    let t_max = members
        .iter()
        .map(|(h, _)| h.draft())
        .fold(0.0f64, f64::max)
        .max(1e-6);

    // Local (principal-value) part: F and M_pv together, refined until the
    // force settles.
    let mut evals = 0usize;
    let mut level = 1usize;
    let (mut f_int, mut m_int) = fleet.local_integral(level, l_max, t_max, &mut evals);
    let mut est_rel = f64::INFINITY;
    for _ in 0..opts.max_refinements {
        level *= 2;
        let (f2, m2) = fleet.local_integral(level, l_max, t_max, &mut evals);
        est_rel = (f2 - f_int).abs() / f2.abs().max(f64::MIN_POSITIVE);
        f_int = f2;
        m_int = m2;
        if est_rel <= opts.rel_tol {
            break;
        }
    }
    // Wave (residue) part of the moment: a 1-D λ integral on the dispersion
    // curve, refined the same way.
    let mut level = 1usize;
    let mut w_int = fleet.wave_moment_integral(level, l_max, t_max, &mut evals);
    for _ in 0..opts.max_refinements {
        level *= 2;
        let w2 = fleet.wave_moment_integral(level, l_max, t_max, &mut evals);
        let rel = (w2 - w_int).abs() / w2.abs().max(f64::MIN_POSITIVE);
        w_int = w2;
        if rel <= opts.rel_tol {
            break;
        }
    }

    let pref = -(rho * u * u) / (2.0 * PI * PI);
    let force_up = pref * f_int;
    let moment_local = pref * m_int;
    let moment_wave = -(4.0 * rho * u * u * nu / PI) * w_int;
    let volume: f64 = members.iter().map(|(h, _)| h.displaced_volume()).sum();
    Ok(DynamicForce {
        force_up,
        moment_bow_up: moment_local + moment_wave,
        moment_local,
        moment_wave,
        lift_fraction: if volume > 0.0 {
            force_up / (rho * g * volume)
        } else {
            0.0
        },
        est_rel_error: est_rel,
        evaluations: evals,
    })
}

/// One member's evaluator plus its fleet-frame placement.
#[derive(Clone)]
struct Member<'h> {
    inner: InnerIntegral<'h>,
    /// Fleet-frame x of the hull's own phase centre.
    cx: f64,
    y: f64,
}

/// `Clone` so each worker thread of a θ fan-out gets its own scratch (the
/// hulls themselves are shared by reference).
#[derive(Clone)]
struct Fleet<'h> {
    members: Vec<Member<'h>>,
    nu: f64,
    x_ref: f64,
    scratch: Vec<SquatTransforms>,
    /// One contraction per member for points off the shared k-grid.
    zc_tmp: Vec<ZContracted>,
}

/// The pair-summed bilinear forms at one wavenumber, real parts taken over
/// the quadrant `k_x, k_y > 0` (the full plane is four copies by the
/// symmetries of real coefficient nets).
#[derive(Clone, Copy)]
struct Forms {
    /// Re Σ_ij q̄_j q_i E_ij
    qq: f64,
    /// Re Σ_ij w̄_j q_i E_ij
    wq: f64,
    /// Σ_ij r̄_j q_i E_ij (complex: the residue needs its imaginary part)
    rq: C64,
    /// Σ_ij r̄_wl,j q_i E_ij
    rwq: C64,
}

impl Default for Forms {
    fn default() -> Self {
        Forms {
            qq: 0.0,
            wq: 0.0,
            rq: C64::ZERO,
            rwq: C64::ZERO,
        }
    }
}

impl<'h> Fleet<'h> {
    fn new(members: &[(&'h Hull, Placement)], nu: f64, x_ref: f64, wave: WaveOptions) -> Self {
        Fleet {
            members: members
                .iter()
                .map(|(h, p)| Member {
                    inner: InnerIntegral::new(h, nu, wave.transom),
                    cx: h.x_center() + p.x,
                    y: p.y,
                })
                .collect(),
            nu,
            x_ref,
            scratch: vec![SquatTransforms::default(); members.len()],
            zc_tmp: vec![ZContracted::default(); members.len()],
        }
    }

    /// Transforms of every member at `(k_x, κ)`, then the pair sums with the
    /// placement phase `e^{ik_x(c_j − c_i)} cos(k_y(y_j − y_i))`.
    fn forms(&mut self, kx: f64, ky: f64, kappa: f64) -> Forms {
        for (m, zc) in self.members.iter_mut().zip(self.zc_tmp.iter_mut()) {
            m.inner.contract_z(kappa, zc);
        }
        for ((m, zc), t) in self
            .members
            .iter_mut()
            .zip(self.zc_tmp.iter())
            .zip(self.scratch.iter_mut())
        {
            *t = m.inner.transforms_at(zc, kx);
        }
        self.pair_sums(kx, ky)
    }

    /// The same from contractions already made at this `κ` — one per member.
    fn forms_cached(&mut self, zcs: &[ZContracted], kx: f64, ky: f64) -> Forms {
        for ((m, zc), t) in self
            .members
            .iter_mut()
            .zip(zcs)
            .zip(self.scratch.iter_mut())
        {
            *t = m.inner.transforms_at(zc, kx);
        }
        self.pair_sums(kx, ky)
    }

    fn pair_sums(&self, kx: f64, ky: f64) -> Forms {
        let mut out = Forms::default();
        for (j, mj) in self.members.iter().enumerate() {
            let tj = self.scratch[j];
            let shift_j = mj.cx - self.x_ref;
            // r = p + q1 + (c_j − x_ref) q, transform of ∂_x[(x − x_ref) f].
            let rj = tj.p + tj.q1 + tj.q.scale(shift_j);
            let rwj = tj.p_wl + tj.q1_wl + tj.w.scale(shift_j);
            for (i, mi) in self.members.iter().enumerate() {
                let qi = self.scratch[i].q;
                let e = C64::cis(kx * (mj.cx - mi.cx)).scale((ky * (mj.y - mi.y)).cos());
                let qi_e = qi * e;
                out.qq += (conj(tj.q) * qi_e).re;
                out.wq += (conj(tj.w) * qi_e).re;
                out.rq = out.rq + conj(rj) * qi_e;
                out.rwq = out.rwq + conj(rwj) * qi_e;
            }
        }
        out
    }

    /// `∫₀^{π/2} dθ PV∫₀^∞ dk [ k·A·Re(q̄q) − (A−1)·Re(w̄q) ]` and the same with
    /// the moment forms — the quadrant integrals with the `k dk dθ` Jacobian
    /// folded in — times 4 for the full plane.
    ///
    /// **Cost structure.** The transforms' z-contraction depends on `κ = k`
    /// only, so the k-nodes sit on one log-spaced grid shared by every θ and
    /// each is contracted once; a (θ, k) point then costs one `osc_moments`
    /// per x-span. That is what makes the 2-D integral affordable on a hull
    /// with hundreds of spans.
    ///
    /// **Pole.** With `A − 1 = (2ν sec²θ)/(k − k₀)`, `k₀ = ν sec²θ`, the
    /// singular part is `h(k)/(k − k₀)`; `∫₀^K [h(k) − h(k₀)]/(k − k₀)` is
    /// smooth and the remainder `h(k₀)·ln((K − k₀)/k₀)` is exact. `h(k₀)`
    /// needs one off-grid evaluation per θ. Where the pole lies beyond the
    /// shared grid (θ near π/2) the range is extended for that θ alone.
    ///
    /// **Range.** The x-transforms of a hull whose `∂f/∂x` jumps at the ends
    /// decay only as `1/k_x²`, so the grid reaches `k_x ≈ 400/L` at every θ;
    /// since `k_x = k cos θ`, θ is cut at `cos θ_c = Q/K`. As θ → π/2 the
    /// integrand per θ tends to a finite constant (the `∫|w(q)|²/q dq` of the
    /// rigid-wall limit), so the remaining sliver is added as
    /// `(π/2 − θ_c)·I(θ_c)`, accurate to `O(cos²θ_c)`.
    fn local_integral(&mut self, level: usize, l: f64, t: f64, evals: &mut usize) -> (f64, f64) {
        let nu = self.nu;
        let (gx, gw) = gauss_legendre(8);
        // Shared k-grid: a linear panel at the origin, then log panels.
        let k_lo = 0.05 / l;
        let k_max = (1.0e4 / l).max(60.0 / t);
        let n_log = 12 * level;
        let mut knodes: Vec<(f64, f64)> = Vec::with_capacity(8 * (n_log + 1));
        let push_panel = |ka: f64, kb: f64, out: &mut Vec<(f64, f64)>| {
            for (gj, &xk) in gx.iter().enumerate() {
                out.push((
                    0.5 * (kb - ka) * xk + 0.5 * (ka + kb),
                    0.5 * (kb - ka) * gw[gj],
                ));
            }
        };
        push_panel(0.0, k_lo, &mut knodes);
        let ratio = (k_max / k_lo).powf(1.0 / n_log as f64);
        let mut ka = k_lo;
        for _ in 0..n_log {
            let kb = ka * ratio;
            push_panel(ka, kb, &mut knodes);
            ka = kb;
        }
        // Contract every member once per k-node.
        let nm = self.members.len();
        let mut zcs: Vec<ZContracted> = vec![ZContracted::default(); knodes.len() * nm];
        for (i, &(k, _)) in knodes.iter().enumerate() {
            for (j, m) in self.members.iter_mut().enumerate() {
                m.inner.contract_z(k, &mut zcs[i * nm + j]);
            }
        }
        *evals += knodes.len();

        // θ range and grid.
        let q_cap = 400.0 / l;
        let theta_c = (q_cap / k_max).min(1.0).acos().min(FRAC_PI_2 - 1e-7);
        let n_theta = 24 * level;

        // The full k-integral for one θ; returns (force, moment) parts.
        let at_theta = |this: &mut Self, theta: f64, evals: &mut usize| -> (f64, f64) {
            let c = theta.cos();
            let sn = theta.sin();
            let k0 = nu / (c * c);
            let two_nu_sec2 = 2.0 * nu / (c * c);
            // Pole beyond the shared grid: extend the range for this θ.
            let k_hi = if k0 >= k_max { 4.0 * k0 } else { k_max };
            // h(k₀) — one off-grid point.
            *evals += 1;
            let fm0 = this.forms(k0 * c, k0 * sn, k0);
            let h0_f = two_nu_sec2 * (k0 * fm0.qq - fm0.wq);
            let h0_m = two_nu_sec2 * (k0 * fm0.rq.re - fm0.rwq.re);
            let mut pf = 0.0;
            let mut pm = 0.0;
            let acc = |k: f64, wk: f64, fm: Forms, pf: &mut f64, pm: &mut f64| {
                let d = k - k0;
                let h_f = two_nu_sec2 * (k * fm.qq - fm.wq);
                let h_m = two_nu_sec2 * (k * fm.rq.re - fm.rwq.re);
                *pf += wk * (k * fm.qq + (h_f - h0_f) / d);
                *pm += wk * (k * fm.rq.re + (h_m - h0_m) / d);
            };
            for (i, &(k, wk)) in knodes.iter().enumerate() {
                *evals += 1;
                let fm = this.forms_cached(&zcs[i * nm..(i + 1) * nm], k * c, k * sn);
                acc(k, wk, fm, &mut pf, &mut pm);
            }
            if k_hi > k_max {
                let n_ext = 4 * level;
                let r = (k_hi / k_max).powf(1.0 / n_ext as f64);
                let mut ka = k_max;
                for _ in 0..n_ext {
                    let kb = ka * r;
                    for (gj, &xk) in gx.iter().enumerate() {
                        let k = 0.5 * (kb - ka) * xk + 0.5 * (ka + kb);
                        let wk = 0.5 * (kb - ka) * gw[gj];
                        *evals += 1;
                        let fm = this.forms(k * c, k * sn, k);
                        acc(k, wk, fm, &mut pf, &mut pm);
                    }
                    ka = kb;
                }
            }
            let log_term = ((k_hi - k0) / k0).ln();
            (pf + h0_f * log_term, pm + h0_m * log_term)
        };

        // θ nodes and weights; the last entry is the sliver θ ∈ (θ_c, π/2),
        // where the integrand is near its θ → π/2 limit.
        let mut nodes: Vec<(f64, f64)> = Vec::with_capacity(n_theta * gx.len() + 1);
        for it in 0..n_theta {
            let (a, b) = (
                theta_c * it as f64 / n_theta as f64,
                theta_c * (it + 1) as f64 / n_theta as f64,
            );
            for (gi, &x) in gx.iter().enumerate() {
                nodes.push((0.5 * (b - a) * x + 0.5 * (a + b), 0.5 * (b - a) * gw[gi]));
            }
        }
        nodes.push((theta_c, FRAC_PI_2 - theta_c));

        // Every node's k-integral is independent given the shared
        // contractions; each worker runs on its own clone of the fleet.
        let this: &Self = self;
        let at_theta = &at_theta;
        let per_node: Vec<(f64, f64, usize)> = crate::parallel::map_indexed(
            nodes.len(),
            || this.clone(),
            |fleet, i| {
                let mut ev = 0usize;
                let (pf, pm) = at_theta(fleet, nodes[i].0, &mut ev);
                (pf, pm, ev)
            },
        );
        let mut f_sum = 0.0;
        let mut m_sum = 0.0;
        for (&(_, wth), &(pf, pm, ev)) in nodes.iter().zip(&per_node) {
            f_sum += wth * pf;
            m_sum += wth * pm;
            *evals += ev;
        }
        (4.0 * f_sum, 4.0 * m_sum)
    }

    /// `∫₁^∞ dλ λ/√(λ²−1) · Im[νλ² r̄q − r̄_wl q]` on the dispersion curve,
    /// via λ = sec θ so the end-point singularity becomes `sec²θ dθ`.
    fn wave_moment_integral(&mut self, level: usize, l: f64, t: f64, evals: &mut usize) -> f64 {
        let nu = self.nu;
        let _ = t;
        // On the curve k_x = ν sec θ; stop where the x-transforms have died.
        let kx_cap = (400.0 / l).max(2.0 * nu);
        let theta_max = (nu / kx_cap).min(1.0).acos().min(FRAC_PI_2 - 1e-7);
        let n_theta = 48 * level;
        let (gx, gw) = gauss_legendre(8);
        let mut nodes: Vec<(f64, f64)> = Vec::with_capacity(n_theta * gx.len());
        for it in 0..n_theta {
            let (a, b) = (
                theta_max * it as f64 / n_theta as f64,
                theta_max * (it + 1) as f64 / n_theta as f64,
            );
            for (gi, &x) in gx.iter().enumerate() {
                nodes.push((0.5 * (b - a) * x + 0.5 * (a + b), 0.5 * (b - a) * gw[gi]));
            }
        }
        let this: &Self = self;
        let vals: Vec<f64> = crate::parallel::map_indexed(
            nodes.len(),
            || this.clone(),
            |fleet, i| {
                let theta = nodes[i].0;
                let lam = 1.0 / theta.cos();
                let kx = nu * lam;
                let k = nu * lam * lam;
                let ky = (k * k - kx * kx).max(0.0).sqrt();
                let fm = fleet.forms(kx, ky, k);
                (fm.rq.scale(k) - fm.rwq).im
            },
        );
        *evals += nodes.len();
        // Same expression, same association, as the serial loop had.
        let mut sum = 0.0;
        for (&(theta, w), &val) in nodes.iter().zip(&vals) {
            let lam = 1.0 / theta.cos();
            sum += w * val * lam * lam;
        }
        sum
    }
}

fn conj(v: C64) -> C64 {
    C64::new(v.re, -v.im)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hulls::wigley;

    /// Sinkage of the Wigley hull: downward, ∝ Fn² at low speed with the
    /// rigid-wall coefficient s/L ≈ 0.021·Fn² (Havelock's 1939 half-ellipsoid
    /// interpolated to these proportions gives ≈ 0.019), and growing with
    /// the free-surface effect to ≈ 1.5e-3 at Fn 0.25 (linear Neumann–Michell:
    /// 1.2e-3; experiment 1.3–1.4e-3 — thin-ship overshoots, as it does on
    /// R_w). Values pinned against an independent numpy implementation of the
    /// same integrals.
    #[test]
    fn wigley_sinks_with_the_expected_coefficient() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let hull = wigley(l, b, t).unwrap();
        let g = 9.80665;
        let a_w = 2.0 * b * l / 3.0;
        let opts = SquatOptions::default();
        for (fn_, want, tol) in [
            (0.10, 0.0210, 0.02),
            (0.25, 0.0236, 0.02),
            (0.40, 0.0322, 0.03),
        ] {
            let u = fn_ * (g * l).sqrt();
            let cond = Conditions::seawater(u);
            let d = dynamic_force(&hull, &cond, 0.0, &opts).unwrap();
            assert!(
                d.force_up < 0.0,
                "Fn {fn_}: force should pull the hull down, got {}",
                d.force_up
            );
            let s_over_l = -d.force_up / (cond.fluid.density * g * a_w) / l;
            let coeff = s_over_l / (fn_ * fn_);
            assert!(
                (coeff - want).abs() < tol * want,
                "Fn {fn_}: s/L/Fn² = {coeff:.4} (want {want:.4}); est err {:.1e}",
                d.est_rel_error
            );
        }
    }

    /// Wigley trim: the local part vanishes by fore-aft symmetry; the wave
    /// part gives the measured pattern — a small bow-down dip near Fn 0.32,
    /// a zero crossing between Fn 0.34 and 0.36 (experiment 0.355), then
    /// bow-up. Magnitude pinned to the numpy reference (+0.73° at Fn 0.40).
    #[test]
    fn wigley_trim_changes_sign_near_fn_0_35() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let hull = wigley(l, b, t).unwrap();
        let g = 9.80665;
        let a = l / 2.0;
        let i_l = 4.0 * a * a * a * b / 15.0;
        let opts = SquatOptions::default();
        let trim_deg = |fn_: f64| {
            let u = fn_ * (g * l).sqrt();
            let cond = Conditions::seawater(u);
            let d = dynamic_force(&hull, &cond, 0.0, &opts).unwrap();
            assert!(
                d.moment_local.abs() < 1e-6 * d.moment_wave.abs().max(1e-3),
                "Fn {fn_}: local moment {} should vanish by symmetry",
                d.moment_local
            );
            (d.moment_bow_up / (cond.fluid.density * g * i_l)).to_degrees()
        };
        assert!(trim_deg(0.32) < 0.0, "expected bow-down near Fn 0.32");
        assert!(trim_deg(0.34) < 0.02, "still ~zero at Fn 0.34");
        assert!(trim_deg(0.36) > 0.05, "bow-up by Fn 0.36");
        let t40 = trim_deg(0.40);
        assert!(
            (t40 - 0.73).abs() < 0.05,
            "Fn 0.40: {t40:.3}° (reference +0.73°)"
        );
    }

    /// A fleet shifted rigidly in x feels the same force and moment.
    #[test]
    fn fleet_force_is_translation_invariant() {
        let hull = wigley(8.0, 0.8, 0.5).unwrap();
        let cond = Conditions::seawater(2.5);
        let opts = SquatOptions::default();
        let at = |dx: f64| {
            let m = [
                (&hull, Placement { x: dx, y: 1.2 }),
                (&hull, Placement { x: dx, y: -1.2 }),
            ];
            let d = multihull_dynamic_force(&m, &cond, dx, &opts).unwrap();
            (d.force_up, d.moment_bow_up)
        };
        let (f0, m0) = at(0.0);
        let (f1, m1) = at(17.0);
        assert!(
            (f0 - f1).abs() < 1e-8 * f0.abs(),
            "force moved: {f0} vs {f1}"
        );
        assert!(
            (m0 - m1).abs() < 1e-8 * m0.abs().max(1e-6),
            "moment moved: {m0} vs {m1}"
        );
    }

    /// Demihull interaction weakens with spacing: the pair's force approaches
    /// twice the solo force as the hulls separate, monotonically over a range
    /// the transverse phase `cos(k_y Δy)` is still resolved on.
    #[test]
    fn demihull_interaction_weakens_with_spacing() {
        let hull = wigley(8.0, 0.8, 0.5).unwrap();
        let cond = Conditions::seawater(2.5);
        let opts = SquatOptions::default();
        let solo = dynamic_force(&hull, &cond, 0.0, &opts).unwrap().force_up;
        let excess = |y: f64| {
            let m = [
                (&hull, Placement { x: 0.0, y }),
                (&hull, Placement { x: 0.0, y: -y }),
            ];
            let pair = multihull_dynamic_force(&m, &cond, 0.0, &opts)
                .unwrap()
                .force_up;
            (pair - 2.0 * solo).abs() / solo.abs()
        };
        let (e1, e2, e4) = (excess(0.6), excess(1.2), excess(2.4));
        assert!(
            e1 > 0.02,
            "hulls 1.2 m apart should interact noticeably, got {e1:.3}"
        );
        assert!(
            e1 > e2 && e2 > e4,
            "interaction should weaken: {e1:.3} > {e2:.3} > {e4:.3}"
        );
    }
}

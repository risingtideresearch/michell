//! Michell's thin-ship wave resistance integral.
//!
//! With ν = g/U², half-beam f(x, z), z downward from the waterline:
//!
//! ```text
//! R_w = (4 ρ g²)/(π U²) ∫_1^∞ (I² + J²) λ²/√(λ² − 1) dλ
//! I(λ) + i J(λ) = ∬ (∂f/∂x) exp(−ν λ² z) exp(i ν λ x) dx dz
//! ```
//!
//! (Tuck 1989; Dambrine, Pierre & Rousseaux 2016.) Substituting λ = sec θ
//! turns the outer integral into `∫_0^{π/2} (I² + J²) sec³θ dθ`, removing the
//! integrable singularity at λ = 1.
//!
//! Because the hull is piecewise polynomial, I and J are evaluated **exactly**
//! per knot span via the closed-form moments in [`crate::moments`]; the only
//! numerical error lives in the smooth outer θ-integral, which is integrated
//! with Gauss–Legendre panels sized to the local oscillation rate and then
//! refined until the requested tolerance is met.

use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::hull::Hull;
use crate::moments::{exp_moments, osc_moments, C64};
use crate::quadrature::gauss_legendre;
use std::f64::consts::{FRAC_PI_2, PI};

/// Options for the outer-integral quadrature.
#[derive(Debug, Clone, Copy)]
pub struct WaveOptions {
    /// Target relative tolerance on the resistance.
    pub rel_tol: f64,
    /// Maximum number of panel-halving refinement passes.
    pub max_refinements: usize,
}

impl Default for WaveOptions {
    fn default() -> Self {
        WaveOptions {
            rel_tol: 1e-5,
            max_refinements: 4,
        }
    }
}

/// Wave resistance result with quadrature diagnostics.
#[derive(Debug, Clone, Copy)]
pub struct WaveResistance {
    /// Wave resistance R_w [N].
    pub resistance: f64,
    /// Estimated relative quadrature error (difference between the last two
    /// refinement passes).
    pub est_rel_error: f64,
    /// Total number of inner-integral evaluations performed.
    pub inner_evaluations: usize,
    /// Largest λ = sec θ reached before truncation.
    pub max_lambda: f64,
}

/// Position of one hull of a multihull, in the global (fleet) frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Placement {
    /// Longitudinal shift **added** to the hull's own x coordinates [m]
    /// (0 keeps the coordinates from the hull's file).
    pub x: f64,
    /// Transverse position of the hull's centerplane [m].
    pub y: f64,
}

/// Compute Michell wave resistance with default options.
pub fn wave_resistance(hull: &Hull, cond: &Conditions) -> Result<WaveResistance> {
    wave_resistance_with(hull, cond, &WaveOptions::default())
}

/// Compute Michell wave resistance with explicit quadrature options.
pub fn wave_resistance_with(
    hull: &Hull,
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    multihull_wave_resistance_with(&[(hull, Placement::default())], cond, opts)
}

/// Combined wave resistance of several thin hulls (multihull), with default
/// options. See [`multihull_wave_resistance_with`].
pub fn multihull_wave_resistance(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
) -> Result<WaveResistance> {
    multihull_wave_resistance_with(members, cond, &WaveOptions::default())
}

/// Combined wave resistance of several thin hulls.
///
/// In thin-ship theory the far-field free-wave amplitudes superpose: hull `j`
/// at longitudinal offset `Δx_j` and transverse position `y_j` contributes
/// `F_j(λ) · exp(i ν (λ Δx_j + λ √(λ²−1) y_j))`, and the resistance integrand
/// uses `|Σ_j …|²`. For two identical hulls separated by `s` this reduces to
/// the classical catamaran interference factor `4 cos²(½ ν s λ √(λ²−1))`.
///
/// Each member hull must itself be symmetric about its own centerplane (the
/// crate's geometry contract); asymmetric demihulls are not modelled.
pub fn multihull_wave_resistance_with(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    cond.validate()?;
    if !(opts.rel_tol.is_finite() && opts.rel_tol > 0.0) {
        return Err(Error::InvalidConditions(
            "rel_tol must be finite and positive".into(),
        ));
    }
    if members.is_empty() {
        return Err(Error::InvalidConditions(
            "at least one hull is required".into(),
        ));
    }
    if members
        .iter()
        .any(|(_, p)| !(p.x.is_finite() && p.y.is_finite()))
    {
        return Err(Error::InvalidConditions(
            "hull placements must be finite".into(),
        ));
    }
    let u = cond.speed;
    let g = cond.gravity;
    let nu = g / (u * u);
    let rho = cond.fluid.density;

    // Phase references (constant overall phase is irrelevant; keeping the
    // oscillatory arguments centred keeps them small).
    let n = members.len() as f64;
    let cx_ref = members
        .iter()
        .map(|(h, p)| h.x_center() + p.x)
        .sum::<f64>()
        / n;
    let y_ref = members.iter().map(|(_, p)| p.y).sum::<f64>() / n;
    let params = OuterParams {
        nu,
        x_half: members
            .iter()
            .map(|(h, p)| (h.x_center() + p.x - cx_ref).abs() + h.x_half_extent())
            .fold(0.0, f64::max),
        y_half: members
            .iter()
            .map(|(_, p)| (p.y - y_ref).abs())
            .fold(0.0, f64::max),
        t_max: members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max),
    };

    let mut inners: Vec<(InnerIntegral, f64, f64)> = members
        .iter()
        .map(|(h, p)| (InnerIntegral::new(h, nu), h.x_center() + p.x - cx_ref, p.y - y_ref))
        .collect();
    let mut amp = |lambda: f64| -> C64 {
        let kx = nu * lambda;
        let ky = nu * lambda * (lambda * lambda - 1.0).max(0.0).sqrt();
        let mut a = C64::ZERO;
        for (inner, dx, dy) in inners.iter_mut() {
            a = a + inner.eval(lambda) * C64::cis(kx * *dx + ky * *dy);
        }
        a
    };

    let mut evals_total = 0usize;
    let mut frac = 1.0;
    let mut evals = 0usize;
    let (mut integral, mut max_lambda) = integrate_outer(&params, frac, &mut amp, &mut evals);
    evals_total += evals;
    let mut est_rel = f64::INFINITY;
    for _ in 0..opts.max_refinements {
        frac *= 0.5;
        let mut evals = 0usize;
        let (refined, ml) = integrate_outer(&params, frac, &mut amp, &mut evals);
        evals_total += evals;
        let scale = refined.abs().max(f64::MIN_POSITIVE);
        est_rel = (refined - integral).abs() / scale;
        integral = refined;
        max_lambda = ml;
        if est_rel <= opts.rel_tol {
            break;
        }
    }

    let coeff = 4.0 * rho * g * g / (PI * u * u);
    Ok(WaveResistance {
        resistance: coeff * integral,
        est_rel_error: est_rel,
        inner_evaluations: evals_total,
        max_lambda,
    })
}

/// The Michell inner integrals `(I(λ), J(λ))` — the free-wave amplitude
/// functions — evaluated exactly (per-span closed forms) for `λ >= 1`.
///
/// Phases are taken relative to the hull's x-midpoint, so I and J individually
/// depend on that (physically irrelevant) choice of origin while `I² + J²`
/// does not.
pub fn inner_integrals(hull: &Hull, cond: &Conditions, lambda: f64) -> Result<(f64, f64)> {
    cond.validate()?;
    if !(lambda.is_finite() && lambda >= 1.0) {
        return Err(Error::InvalidConditions(format!(
            "lambda must be >= 1, got {lambda}"
        )));
    }
    let nu = cond.gravity / (cond.speed * cond.speed);
    let mut inner = InnerIntegral::new(hull, nu);
    let f = inner.eval(lambda);
    Ok((f.re, f.im))
}

/// Geometry-derived phase-rate parameters for the outer quadrature.
struct OuterParams {
    nu: f64,
    /// Half-extent of the whole fleet about its longitudinal phase centre.
    x_half: f64,
    /// Largest transverse offset from the fleet's phase centre.
    y_half: f64,
    /// Deepest draft in the fleet.
    t_max: f64,
}

/// Marching-panel Gauss–Legendre integration of
/// `∫_0^{π/2} |A(sec θ)|² sec³θ dθ` for a combined amplitude `A`.
fn integrate_outer(
    params: &OuterParams,
    frac: f64,
    amp: &mut impl FnMut(f64) -> C64,
    evals: &mut usize,
) -> (f64, f64) {
    const GL_N: usize = 16;
    /// Truncate once a full quiet window contributes below this fraction.
    const STOP_REL: f64 = 1e-9;
    /// Width of the quiet window in accumulated phase — several full periods
    /// of the oscillating integrand, so a trough of cos²-type oscillation
    /// (width < π) can never trigger truncation on its own. Phase-based, so
    /// the criterion is independent of the panel-refinement level.
    const STOP_WINDOW_PHASE: f64 = 8.0 * PI;
    const LAMBDA_HARD_CAP: f64 = 1e4;
    const MAX_EVALS_PER_PASS: usize = 4_000_000;

    let (gx, gw) = gauss_legendre(GL_N);
    let nu = params.nu;
    let (x_half, y_half, t_max) = (params.x_half, params.y_half, params.t_max);

    // Local phase rate of |A|² in θ: the x-oscillation contributes
    // 2 ν x_half d(sec θ)/dθ, the z-decay envelope 2 ν T d(sec²θ)/dθ, and the
    // transverse separation phase ν y λ√(λ²−1) = ν y sec θ tan θ contributes
    // 2 ν y_half d(sec θ tan θ)/dθ = 2 ν y_half sec θ (sec²θ + tan²θ).
    let rate = |theta: f64| -> f64 {
        let sec = 1.0 / theta.cos();
        let tan = theta.tan();
        2.0 * nu * sec * tan * (x_half + t_max * sec)
            + 2.0 * nu * y_half * sec * (sec * sec + tan * tan)
            + 4.0
    };
    // Near θ = 0 the longitudinal phase grows like ν x_half θ², so cap the
    // first panels at one period of that quadratic phase (the transverse
    // phase is linear near 0 and already covered by rate(0)).
    let cap = (2.0 * PI / (2.0 * nu * x_half).sqrt().max(1.0)).min(0.12);

    let mut theta = 0.0f64;
    let mut total = 0.0f64;
    let mut window_sum = 0.0f64;
    let mut window_phase = 0.0f64;
    let mut pass_evals = 0usize;
    while theta < FRAC_PI_2 - 1e-12 {
        let local_rate = rate(theta);
        let dt = (frac * 2.0 * PI / local_rate)
            .min(frac * cap)
            .min(FRAC_PI_2 - theta)
            .max(1e-15);
        let half = dt / 2.0;
        let mid = theta + half;
        let mut panel = 0.0;
        for (i, &xi) in gx.iter().enumerate() {
            let th = mid + half * xi;
            let sec = 1.0 / th.cos();
            let a = amp(sec);
            panel += gw[i] * a.abs_sq() * sec * sec * sec;
        }
        panel *= half;
        total += panel;
        pass_evals += GL_N;
        theta += dt;

        // Truncation: only past λ = 2, and only when an entire window of
        // accumulated oscillation phase contributed negligibly.
        if 1.0 / theta.cos() > 2.0 {
            window_sum += panel;
            window_phase += local_rate * dt;
            if window_phase >= STOP_WINDOW_PHASE {
                if window_sum.abs() <= STOP_REL * total.abs() + f64::MIN_POSITIVE {
                    break;
                }
                window_sum = 0.0;
                window_phase = 0.0;
            }
        }
        if 1.0 / theta.cos() > LAMBDA_HARD_CAP || pass_evals > MAX_EVALS_PER_PASS {
            break;
        }
    }
    *evals += pass_evals;
    (total, 1.0 / theta.cos().max(1e-300))
}

/// Exact (per-span closed-form) evaluation of I + iJ at λ.
struct InnerIntegral<'h> {
    hull: &'h Hull,
    nu: f64,
    /// Scratch: z-moments including the e^{−κ z0} shift, [n_spans_z][q+1].
    zm: Vec<f64>,
    /// Scratch: raw x-moments for one span.
    xm: Vec<C64>,
    /// Scratch: raw z-moments for one span.
    zm_raw: Vec<f64>,
    p: usize,
    q: usize,
}

impl<'h> InnerIntegral<'h> {
    fn new(hull: &'h Hull, nu: f64) -> Self {
        let p = hull.surface().degree_x();
        let q = hull.surface().degree_z();
        InnerIntegral {
            hull,
            nu,
            zm: vec![0.0; hull.spans_z().len() * (q + 1)],
            xm: Vec::with_capacity(p),
            zm_raw: Vec::with_capacity(q + 1),
            p,
            q,
        }
    }

    /// I + iJ at λ = sec θ. Phases use x relative to the hull midpoint (a pure
    /// phase factor on I + iJ that leaves |I + iJ|² unchanged, but keeps the
    /// oscillatory arguments as small as possible).
    fn eval(&mut self, lambda: f64) -> C64 {
        let nu = self.nu;
        let kx = nu * lambda;
        let kappa = nu * lambda * lambda;
        let (p, q) = (self.p, self.q);
        let spans_z = self.hull.spans_z();
        let spans_x = self.hull.spans_x();
        let nsz = spans_z.len();

        // z-moments per span, shifted by the decay to the span start.
        let mut any = false;
        for (t, sz) in spans_z.iter().enumerate() {
            let decay = (-kappa * sz.start).exp();
            if decay == 0.0 {
                for b in 0..=q {
                    self.zm[t * (q + 1) + b] = 0.0;
                }
                continue;
            }
            any = true;
            exp_moments(kappa, sz.len, q, &mut self.zm_raw);
            for b in 0..=q {
                self.zm[t * (q + 1) + b] = decay * self.zm_raw[b];
            }
        }
        if !any {
            return C64::ZERO;
        }

        let coeff = self.hull.fx_coeff();
        let x_center = self.hull.x_center();
        let mut f = C64::ZERO;
        for (s, sx) in spans_x.iter().enumerate() {
            osc_moments(kx, sx.len, p - 1, &mut self.xm);
            let phase = C64::cis(kx * (sx.start - x_center));
            let mut span_sum = C64::ZERO;
            for (a, &xma) in self.xm.iter().enumerate() {
                // g_a = Σ_t Σ_b c[s,t,a,b] · Zm[t][b]
                let mut g_a = 0.0;
                for t in 0..nsz {
                    let base = ((s * nsz + t) * p + a) * (q + 1);
                    let zrow = &self.zm[t * (q + 1)..(t + 1) * (q + 1)];
                    let crow = &coeff[base..base + q + 1];
                    for b in 0..=q {
                        g_a += crow[b] * zrow[b];
                    }
                }
                span_sum = span_sum + xma.scale(g_a);
            }
            f = f + phase * span_sum;
        }
        f
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hulls::wigley;

    #[test]
    fn inner_integral_is_translation_invariant_in_modulus() {
        let l = 10.0;
        let hull = wigley(l, 1.0, 0.625).unwrap();
        // Same hull, shifted +37 m in x.
        let s = hull.surface();
        let shifted = crate::bspline::BSplineSurface::new(
            s.degree_x(),
            s.degree_z(),
            s.knots_x().iter().map(|k| k + 37.0).collect(),
            s.knots_z().to_vec(),
            s.control().to_vec(),
        )
        .unwrap();
        let hull2 = Hull::new(shifted).unwrap();
        let nu = 9.80665 / (3.0f64 * 3.0);
        let mut i1 = InnerIntegral::new(&hull, nu);
        let mut i2 = InnerIntegral::new(&hull2, nu);
        for lambda in [1.0, 1.4, 2.5, 6.0] {
            let a = i1.eval(lambda).abs_sq();
            let b = i2.eval(lambda).abs_sq();
            assert!(
                (a - b).abs() <= 1e-9 * a.abs().max(1e-300),
                "lambda={lambda}: {a} vs {b}"
            );
        }
    }
}

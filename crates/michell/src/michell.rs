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
//!
//! ## Asymmetric hulls
//!
//! A hull built with [`Hull::new_asymmetric`] is split into a symmetric
//! thickness part `f_sym = (f₊ + f₋)/2` and an antisymmetric camber part
//! `f_a = (f₊ − f₋)/2`. The thickness part is the classical **source** system
//! above; the camber part adds a centreplane **y-dipole** system whose
//! free-wave amplitude is the same inner integral over `∂f_a/∂x` weighted by
//! the transverse wavenumber (see [`dipole_weight`]). Because the source
//! amplitude is even in θ and the dipole amplitude is odd, their cross term
//! integrates to zero over the Kelvin fan, so
//!
//! ```text
//! R_w = R_source(f_sym) + R_dipole(f_a)
//! ```
//!
//! with the symmetric hull (`f_a ≡ 0`) recovering classical Michell exactly.
//! The dipole *magnitude* uses an approximate closure and should be read
//! qualitatively — see [`DIPOLE_WEIGHT_C`] for the caveat.

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
/// `F_j(λ) · exp(i ν (λ Δx_j ± λ √(λ²−1) y_j))` to the wave system
/// propagating at ±θ off the track, and the resistance integrand is the mean
/// of the two systems, `½ (|Σ_j …₊|² + |Σ_j …₋|²)` (the θ < 0 half of the
/// free-wave spectrum; the two differ only for fleets that are not
/// mirror-symmetric about their mean centerplane, e.g. staggered pairs). For
/// two identical hulls separated by `s` this reduces to the classical
/// catamaran interference factor `4 cos²(½ ν s λ √(λ²−1))`.
///
/// A member built with [`Hull::new_asymmetric`] additionally contributes a
/// centreplane y-dipole wave system from its antisymmetric half-beam `f_a`
/// (weighted by `dipole_weight`); this superposes on the source system exactly
/// like any other free-wave amplitude, so demihull asymmetry and demihull
/// interference are handled together. The dipole *magnitude* rests on an
/// approximate strip closure — see the `DIPOLE_WEIGHT_C` constant.
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
    let mut amp_sq = |lambda: f64| -> f64 {
        let kx = nu * lambda;
        let ky = nu * lambda * (lambda * lambda - 1.0).max(0.0).sqrt();
        // Dipole spectral weight: a y-normal centreplane doublet radiates the
        // source integral scaled by the transverse wavenumber (see
        // [`dipole_weight`]). It is *odd* in θ — the +θ wave system carries
        // −w·G and the −θ system +w·G — so the ½(|A₊|² + |A₋|²) average cancels
        // the source–dipole cross term, leaving R_w = R_source + R_dipole
        // additively (the θ-parity argument confirmed in the module docs).
        let w = dipole_weight(lambda);
        let mut plus = C64::ZERO;
        let mut minus = C64::ZERO;
        for (inner, dx, dy) in inners.iter_mut() {
            let (f, g) = inner.eval_pair(lambda);
            let wg = g.map_or(C64::ZERO, |g| g.scale(w));
            if f == C64::ZERO && wg == C64::ZERO {
                continue;
            }
            let a_plus = f - wg;
            let a_minus = f + wg;
            plus = plus + a_plus * C64::cis(kx * *dx + ky * *dy);
            minus = minus + a_minus * C64::cis(kx * *dx - ky * *dy);
        }
        0.5 * (plus.abs_sq() + minus.abs_sq())
    };

    let mut evals_total = 0usize;
    let mut frac = 1.0;
    let mut evals = 0usize;
    let (mut integral, mut max_lambda) = integrate_outer(&params, frac, &mut amp_sq, &mut evals);
    evals_total += evals;
    let mut est_rel = f64::INFINITY;
    for _ in 0..opts.max_refinements {
        frac *= 0.5;
        let mut evals = 0usize;
        let (refined, ml) = integrate_outer(&params, frac, &mut amp_sq, &mut evals);
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

/// Closure constant for the antisymmetric (dipole) free-wave amplitude — see
/// [`dipole_weight`]. `1.0` is the value self-consistent with this crate's
/// source normalisation *under the prescribed strip closure* `μ = 2U f_a` (the
/// naive antisymmetric analogue of Michell's `σ = 2U ∂f_sym/∂x`).
///
/// **This magnitude is approximate.** Unlike the source strength — which the
/// boundary condition fixes *pointwise*, because a source sheet's normal-
/// velocity *jump* equals its local strength — the dipole density is fixed by
/// the *mean* of the two-sided normal velocities, and a doublet sheet's mean
/// normal velocity is a hypersingular (finite-part) integral of `μ`. So `μ` is
/// genuinely **non-local** in `∂f_a/∂x`; no exact pointwise closure exists, and
/// the rigorous density solves a hypersingular Fredholm equation of the first
/// kind (Kaklis & Papanikolaou; 21st Symp. Naval Hydro., Appendix A). `μ = 2U
/// f_a` is the crudest strip estimate of that solve — the *structure* below
/// (weight, θ-parity, additive separation) is exact, but the overall constant
/// should be validated against a reference before the dipole magnitude is
/// trusted quantitatively.
const DIPOLE_WEIGHT_C: f64 = 1.0;

/// Real spectral weight relating the antisymmetric dipole amplitude to the
/// companion source integral `G(λ) = ∬ (∂f_a/∂x) e^{−κz} e^{iνλx} dx dz`.
///
/// A y-normal centreplane doublet of density `μ(x,z)` radiates a free wave
/// carrying one extra factor of the **transverse wavenumber**
/// `k_y = νλ√(λ²−1)` relative to a source:
/// `A_d(θ) = i k_y ∬ μ e^{−κz} e^{iνλx} dx dz`. With the prescribed closure
/// `μ = 2U f_a`, `f_a` closing at bow and stern, integration by parts in `x`
/// rewrites `∬ f_a e^{iνλx} = (i/νλ) G`, and the two factors of `i` collapse to
/// a **real** weight on `G` in this crate's normalisation:
///
/// ```text
/// A_d(λ) = i · k_y · (2U) · (i/νλ) G(λ) / (2U)   (÷2U → crate amplitude units)
///        = −√(λ²−1) · DIPOLE_WEIGHT_C · G(λ).
/// ```
///
/// The sign is carried in [`multihull_wave_resistance_with`] (the +θ system
/// gets `−w·G`); this function returns the magnitude `w = √(λ²−1)·C`. It
/// vanishes as `λ → 1` (θ → 0): a wave running dead ahead carries no transverse
/// wavenumber and is blind to side-to-side asymmetry, so `R_dipole` is fed
/// entirely by the diverging (large-λ) part of the spectrum.
#[inline]
fn dipole_weight(lambda: f64) -> f64 {
    DIPOLE_WEIGHT_C * (lambda * lambda - 1.0).max(0.0).sqrt()
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
/// `∫_0^{π/2} |A(sec θ)|² sec³θ dθ`, where `amp_sq` supplies the combined
/// `|A|²` of the fleet's two (±θ) wave systems.
fn integrate_outer(
    params: &OuterParams,
    frac: f64,
    amp_sq: &mut impl FnMut(f64) -> f64,
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
            panel += gw[i] * amp_sq(sec) * sec * sec * sec;
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
pub(crate) struct InnerIntegral<'h> {
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
    pub(crate) fn new(hull: &'h Hull, nu: f64) -> Self {
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

    /// Source free-wave amplitude `I + iJ` at λ = sec θ. Phases use x relative
    /// to the hull midpoint (a pure phase factor on I + iJ that leaves
    /// |I + iJ|² unchanged, but keeps the oscillatory arguments as small as
    /// possible).
    pub(crate) fn eval(&mut self, lambda: f64) -> C64 {
        self.eval_pair(lambda).0
    }

    /// The source amplitude `F(λ)` and, for an asymmetric hull, the companion
    /// **dipole** amplitude `G(λ) = ∬ (∂f_a/∂x) e^{−κz} e^{iνλx} dx dz` from the
    /// antisymmetric half-beam. `G` is `None` for a symmetric hull. Both share
    /// the single z-moment pass; only the x-span accumulation is repeated with
    /// the second coefficient array, so an asymmetric evaluation costs little
    /// more than a symmetric one.
    pub(crate) fn eval_pair(&mut self, lambda: f64) -> (C64, Option<C64>) {
        let hull = self.hull;
        let kx = self.nu * lambda;
        let kappa = self.nu * lambda * lambda;
        if !self.fill_zm(kappa) {
            return (C64::ZERO, hull.fx_a_coeff().map(|_| C64::ZERO));
        }
        let f = self.accumulate(kx, hull.fx_coeff());
        let g = hull.fx_a_coeff().map(|c| self.accumulate(kx, c));
        (f, g)
    }

    /// Fill the per-span z-moments (shifted by the decay to each span start) at
    /// decay rate `kappa = ν λ²`. Returns `false` if every span has underflowed
    /// (the whole hull is below the exponential's support), in which case the
    /// amplitudes are zero.
    fn fill_zm(&mut self, kappa: f64) -> bool {
        let hull = self.hull;
        let q = self.q;
        let mut any = false;
        for (t, sz) in hull.spans_z().iter().enumerate() {
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
        any
    }

    /// Accumulate `∬ (Σ coeff·XᵃZᵇ) e^{−κz} e^{iνλx} dx dz` over all span pairs
    /// using the already-filled z-moments and the per-span x-oscillation
    /// moments. `coeff` is either the source `∂f/∂x` net or the antisymmetric
    /// `∂f_a/∂x` net (identical layout), so source and dipole amplitudes reuse
    /// this one kernel.
    fn accumulate(&mut self, kx: f64, coeff: &[f64]) -> C64 {
        let hull = self.hull;
        let (p, q) = (self.p, self.q);
        let nsz = hull.spans_z().len();
        let x_center = hull.x_center();
        let mut f = C64::ZERO;
        for (s, sx) in hull.spans_x().iter().enumerate() {
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
    fn staggered_fleet_resistance_is_mirror_symmetric() {
        // A staggered pair and its mirror image about y = 0 are the same
        // physical system; the resistance must not change. (Regression: the
        // integrand once used only the +θ wave system, which broke this.)
        let hull = wigley(8.0, 0.8, 0.5).unwrap();
        let cond = crate::Conditions::seawater(2.5);
        let a = [
            (&hull, Placement { x: 0.0, y: 1.4 }),
            (&hull, Placement { x: 1.7, y: -1.4 }),
        ];
        let b = [
            (&hull, Placement { x: 0.0, y: -1.4 }),
            (&hull, Placement { x: 1.7, y: 1.4 }),
        ];
        let ra = multihull_wave_resistance(&a, &cond).unwrap().resistance;
        let rb = multihull_wave_resistance(&b, &cond).unwrap().resistance;
        assert!(
            (ra - rb).abs() <= 1e-9 * ra,
            "staggered {ra} vs mirrored {rb}"
        );
    }

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

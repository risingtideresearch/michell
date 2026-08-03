//! Low-Froude reduction of Michell's integral to waterline endpoint waves.
//!
//! Repeated integration by parts rewrites the exact transform of every
//! polynomial span as a finite sum of terms
//!
//! ```text
//! c exp(i nu lambda x_e) exp(-nu lambda^2 z_e) / lambda^n.
//! ```
//!
//! At low Froude number the terms with `z_e > 0` are exponentially small.
//! Keeping the waterline terms turns the outer integral into pairwise Bickley
//! kernels.  Their oscillatory phase is integrated on its exact steepest-
//! descent contour, so the work does not grow with `nu |x_i-x_j|`.

use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::hull::Hull;
use crate::moments::C64;
use crate::quadrature::gauss_legendre;
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_1_SQRT_2, PI};

/// Michell resistance from the low-Froude waterline-endpoint reduction.
#[derive(Debug, Clone)]
pub struct LowFroudeResistance {
    /// Approximate wave resistance, in newtons.
    pub resistance: f64,
    /// Estimated relative error: coefficient construction, endpoint-map
    /// accumulation, omitted submerged endpoints, and contour quadrature.
    pub est_rel_error: f64,
    /// Resistance-level relative error floor for the shared span coefficients.
    pub coefficient_rel_error_bound: f64,
    /// Absolute resistance bound for roundoff while combining endpoint terms.
    pub endpoint_summation_abs_error_bound: f64,
    /// Analytical absolute resistance bound, in newtons, for all pair terms
    /// involving at least one omitted endpoint below the waterline.
    pub omitted_abs_error_bound: f64,
    /// Absolute resistance error estimate, in newtons, from coarse/fine
    /// quadrature on the steepest-descent contours.
    pub quadrature_abs_error_estimate: f64,
    /// Number of nonzero terms in the exact endpoint representation.
    pub endpoint_terms: usize,
    /// Number of those terms located at the waterline and retained.
    pub waterline_terms: usize,
    /// Total scalar quadrature nodes used by the coarse and fine contour
    /// passes. Zero-frequency kernels are analytic and use no nodes.
    pub kernel_evaluations: usize,
    /// Signed resistance attribution for every retained waterline term pair.
    /// These entries sum to [`Self::resistance`].
    pub endpoint_pairs: Vec<EndpointPairContribution>,
}

/// Longitudinal location represented by a low-Froude endpoint wave.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndpointKind {
    /// Low-x end; the stern because the ship advances toward +x.
    Stern,
    /// A non-empty-span boundary inside the hull domain.
    InteriorKnot,
    /// High-x end; the bow because the ship advances toward +x.
    Bow,
}

/// One retained waterline term in the endpoint expansion of the exact inner
/// amplitude.
#[derive(Debug, Clone, Copy)]
pub struct EndpointWave {
    pub kind: EndpointKind,
    /// Hull-coordinate endpoint location \[m].
    pub x: f64,
    /// Endpoint depth \[m]. Retained terms are currently all at `z = 0`.
    pub z: f64,
    /// Power `n` in `coefficient · exp(i ν λ x) / λⁿ`.
    pub lambda_power: usize,
    /// Complex coefficient in the endpoint representation.
    pub coefficient: C64,
}

/// Signed contribution of one upper-triangular endpoint-wave pair to the
/// retained low-Froude resistance.
#[derive(Debug, Clone, Copy)]
pub struct EndpointPairContribution {
    pub left: EndpointWave,
    pub right: EndpointWave,
    /// Pair contribution to resistance \[N], including the factor of two for
    /// an off-diagonal pair. Interference contributions may be negative.
    pub resistance: f64,
    /// Signed fraction `resistance / total resistance`.
    pub resistance_fraction: f64,
    /// Absolute coarse/fine contour difference attributed to this pair \[N].
    pub quadrature_abs_error_estimate: f64,
    /// Scalar quadrature nodes used for this kernel (zero for equal-x pairs).
    pub kernel_evaluations: usize,
}

#[derive(Debug, Clone, Copy)]
struct EndpointTerm {
    x: f64,
    z: f64,
    lambda_power: usize,
    coeff: C64,
    coeff_abs_error_bound: f64,
}

#[derive(Debug, Clone, Copy)]
struct ComplexAccumulator {
    sum: C64,
    correction: C64,
    absolute_sum: f64,
    additions: usize,
}

impl Default for ComplexAccumulator {
    fn default() -> Self {
        Self {
            sum: C64::ZERO,
            correction: C64::ZERO,
            absolute_sum: 0.0,
            additions: 0,
        }
    }
}

impl ComplexAccumulator {
    fn add(&mut self, value: C64) {
        let adjusted = value - self.correction;
        let next = self.sum + adjusted;
        self.correction = (next - self.sum) - adjusted;
        self.sum = next;
        self.absolute_sum += value.abs();
        self.additions += 1;
    }

    fn error_bound(self) -> f64 {
        let operations = 2.0 * self.additions as f64 + 2.0;
        let gamma = operations * f64::EPSILON;
        if gamma < 1.0 {
            gamma / (1.0 - gamma) * self.absolute_sum
        } else {
            f64::INFINITY
        }
    }
}

/// Evaluate the low-Froude waterline-endpoint reduction of Michell's integral.
///
/// This is an asymptotic companion to [`crate::wave_resistance`], not a blind
/// replacement for it. The exact piecewise-polynomial inner amplitude is first
/// decomposed into endpoint waves. Terms based below `z = 0` are omitted; the
/// returned [`LowFroudeResistance::omitted_abs_error_bound`] bounds their total
/// possible contribution (including cross terms) before oscillatory
/// cancellation. Callers should accept the result only when
/// [`LowFroudeResistance::est_rel_error`] meets their accuracy requirement.
///
/// The retained waterline pair kernels are integrated after
/// `lambda = 1 + t^2` and the exact contour rotation
/// `t = exp(i pi/4) y / sqrt(omega)`. The resulting Gaussian-decaying integral
/// has a cost independent of the longitudinal frequency `omega`, which is the
/// source of the speed advantage as Froude number decreases.
///
/// Currently restricted to upright symmetric hulls whose spline domain begins
/// exactly at `z = 0`.
pub fn low_froude_wave_resistance(hull: &Hull, cond: &Conditions) -> Result<LowFroudeResistance> {
    cond.validate()?;
    if hull.is_asymmetric() {
        return Err(Error::Unsupported(
            "low-Froude endpoint reduction currently requires a symmetric hull".into(),
        ));
    }
    if hull.surface().z_domain().0 != 0.0 {
        return Err(Error::Unsupported(
            "low-Froude endpoint reduction requires the spline domain to start exactly at z = 0"
                .into(),
        ));
    }
    if hull.surface().degree_x() == 0 {
        return Err(Error::Unsupported(
            "low-Froude endpoint reduction requires degree_x >= 1".into(),
        ));
    }

    let nu = cond.gravity / (cond.speed * cond.speed);
    let terms = endpoint_terms(hull, nu);
    if terms
        .iter()
        .any(|term| !(term.x.is_finite() && term.z.is_finite() && term.coeff.is_finite()))
    {
        return Err(Error::Unsupported(
            "non-finite endpoint coefficient encountered in low-Froude resistance".into(),
        ));
    }
    let waterline: Vec<_> = terms.iter().copied().filter(|term| term.z == 0.0).collect();
    if waterline.is_empty() {
        return Err(Error::Unsupported(
            "low-Froude endpoint reduction found no nonzero waterline terms".into(),
        ));
    }
    const MIN_STEEPEST_DESCENT_FREQUENCY: f64 = 25.0;
    for (i, left) in waterline.iter().enumerate() {
        for right in &waterline[i + 1..] {
            let frequency = nu * (left.x - right.x).abs();
            if frequency > 0.0 && frequency < MIN_STEEPEST_DESCENT_FREQUENCY {
                return Err(Error::Unsupported(format!(
                    "low-Froude endpoint reduction needs every nonzero endpoint phase frequency \
                     to be at least {MIN_STEEPEST_DESCENT_FREQUENCY}; found {frequency:.6}"
                )));
            }
        }
    }

    let (coarse_x, coarse_w) = gauss_legendre(24);
    let (fine_x, fine_w) = gauss_legendre(48);
    let physical_coeff =
        4.0 * cond.fluid.density * cond.gravity * cond.gravity / (PI * cond.speed.powi(2));
    let (x0, x1) = hull.surface().x_domain();
    let mut integral = 0.0;
    let mut quadrature_error = 0.0;
    let mut kernel_evaluations = 0usize;
    let mut endpoint_pairs = Vec::with_capacity(waterline.len() * (waterline.len() + 1) / 2);
    for (i, left) in waterline.iter().enumerate() {
        for right in &waterline[i..] {
            let product = left.coeff * conjugate(right.coeff);
            let s = left.lambda_power + right.lambda_power - 2;
            let omega = nu * (left.x - right.x);
            let (kernel, kernel_error, evaluations) =
                oscillatory_kernel(s, omega, &coarse_x, &coarse_w, &fine_x, &fine_w);
            let multiplicity = if std::ptr::eq(left, right) { 1.0 } else { 2.0 };
            let pair_integral = multiplicity * real_product(product, kernel);
            let pair_error = multiplicity * product.abs() * kernel_error;
            integral += pair_integral;
            quadrature_error += pair_error;
            kernel_evaluations += evaluations;
            endpoint_pairs.push(EndpointPairContribution {
                left: endpoint_wave(*left, x0, x1),
                right: endpoint_wave(*right, x0, x1),
                resistance: physical_coeff * pair_integral,
                resistance_fraction: 0.0,
                quadrature_abs_error_estimate: physical_coeff * pair_error,
                kernel_evaluations: evaluations,
            });
        }
    }

    let mut omitted_bound = 0.0;
    let mut endpoint_summation_bound = 0.0;
    for left in &terms {
        for right in &terms {
            let s = left.lambda_power + right.lambda_power - 2;
            let decay = (-nu * (left.z + right.z)).exp();
            let kernel_bound = decay * zero_frequency_kernel(s);
            if !(left.z == 0.0 && right.z == 0.0) {
                omitted_bound += left.coeff.abs() * right.coeff.abs() * kernel_bound;
            }
            endpoint_summation_bound += (left.coeff_abs_error_bound * right.coeff.abs()
                + left.coeff.abs() * right.coeff_abs_error_bound
                + left.coeff_abs_error_bound * right.coeff_abs_error_bound)
                * kernel_bound;
        }
    }

    let resistance = physical_coeff * integral;
    let omitted_abs_error_bound = physical_coeff * omitted_bound;
    let quadrature_abs_error_estimate = physical_coeff * quadrature_error;
    let endpoint_summation_abs_error_bound = physical_coeff * endpoint_summation_bound;
    let coefficient_rel_error_bound = hull.fx_coeff_rel_error_bound();
    let abs_error = omitted_abs_error_bound
        + quadrature_abs_error_estimate
        + endpoint_summation_abs_error_bound;
    let est_rel_error =
        coefficient_rel_error_bound + abs_error / resistance.abs().max(f64::MIN_POSITIVE);
    for pair in &mut endpoint_pairs {
        pair.resistance_fraction = pair.resistance / resistance;
    }

    if !(resistance.is_finite()
        && est_rel_error.is_finite()
        && omitted_abs_error_bound.is_finite()
        && quadrature_abs_error_estimate.is_finite()
        && endpoint_summation_abs_error_bound.is_finite()
        && endpoint_pairs.iter().all(|pair| {
            pair.resistance.is_finite()
                && pair.resistance_fraction.is_finite()
                && pair.quadrature_abs_error_estimate.is_finite()
        }))
    {
        return Err(Error::Unsupported(
            "non-finite value encountered in low-Froude resistance evaluation".into(),
        ));
    }

    Ok(LowFroudeResistance {
        resistance,
        est_rel_error,
        coefficient_rel_error_bound,
        endpoint_summation_abs_error_bound,
        omitted_abs_error_bound,
        quadrature_abs_error_estimate,
        endpoint_terms: terms.len(),
        waterline_terms: waterline.len(),
        kernel_evaluations,
        endpoint_pairs,
    })
}

fn endpoint_wave(term: EndpointTerm, x0: f64, x1: f64) -> EndpointWave {
    let kind = if term.x == x0 {
        EndpointKind::Stern
    } else if term.x == x1 {
        EndpointKind::Bow
    } else {
        EndpointKind::InteriorKnot
    };
    EndpointWave {
        kind,
        x: term.x,
        z: term.z,
        lambda_power: term.lambda_power,
        coefficient: term.coeff,
    }
}

fn endpoint_terms(hull: &Hull, nu: f64) -> Vec<EndpointTerm> {
    let p = hull.surface().degree_x();
    let q = hull.surface().degree_z();
    let nsz = hull.spans_z().len();
    let mut combined: BTreeMap<(usize, usize, usize), ComplexAccumulator> = BTreeMap::new();

    for (sx_index, sx) in hull.spans_x().iter().enumerate() {
        for (sz_index, sz) in hull.spans_z().iter().enumerate() {
            for a in 0..p {
                for b in 0..=q {
                    let index = ((sx_index * nsz + sz_index) * p + a) * (q + 1) + b;
                    let polynomial_coeff = hull.fx_coeff()[index];
                    if polynomial_coeff == 0.0 {
                        continue;
                    }
                    for r in 0..=a {
                        let derivative_x = falling_factorial(a, r);
                        let right_x =
                            (-1.0f64).powi(r as i32) * derivative_x * sx.len.powi((a - r) as i32);
                        let mut x_endpoints = [(sx_index + 1, right_x), (sx_index, 0.0)];
                        if r == a {
                            x_endpoints[1].1 = -(-1.0f64).powi(r as i32) * derivative_x;
                        }
                        for u in 0..=b {
                            let derivative_z = falling_factorial(b, u);
                            let right_z = -derivative_z * sz.len.powi((b - u) as i32);
                            let mut z_endpoints = [(sz_index + 1, right_z), (sz_index, 0.0)];
                            if u == b {
                                z_endpoints[1].1 = derivative_z;
                            }
                            let lambda_power = r + 2 * u + 3;
                            let scale = polynomial_coeff / nu.powi((r + u + 2) as i32);
                            let complex_scale = inverse_i_power(r + 1).scale(scale);
                            for &(ix, x_factor) in &x_endpoints {
                                if x_factor == 0.0 {
                                    continue;
                                }
                                for &(iz, z_factor) in &z_endpoints {
                                    if z_factor == 0.0 {
                                        continue;
                                    }
                                    let key = (ix, iz, lambda_power);
                                    let contribution = complex_scale.scale(x_factor * z_factor);
                                    combined.entry(key).or_default().add(contribution);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let mut x_endpoints: Vec<f64> = hull.spans_x().iter().map(|span| span.start).collect();
    let last_x = hull.spans_x().last().expect("validated hull has x spans");
    x_endpoints.push(last_x.start + last_x.len);
    let mut z_endpoints: Vec<f64> = hull.spans_z().iter().map(|span| span.start).collect();
    let last_z = hull.spans_z().last().expect("validated hull has z spans");
    z_endpoints.push(last_z.start + last_z.len);

    combined
        .into_iter()
        .filter_map(|((ix, iz, lambda_power), accumulator)| {
            let error_bound = accumulator.error_bound();
            (accumulator.sum != C64::ZERO || error_bound > 0.0).then_some(EndpointTerm {
                x: x_endpoints[ix],
                z: z_endpoints[iz],
                lambda_power,
                coeff: accumulator.sum,
                coeff_abs_error_bound: error_bound,
            })
        })
        .collect()
}

#[inline]
fn falling_factorial(n: usize, r: usize) -> f64 {
    ((n - r + 1)..=n).fold(1.0, |product, value| product * value as f64)
}

#[inline]
fn inverse_i_power(power: usize) -> C64 {
    match power % 4 {
        0 => C64::new(1.0, 0.0),
        1 => C64::new(0.0, -1.0),
        2 => C64::new(-1.0, 0.0),
        _ => C64::new(0.0, 1.0),
    }
}

#[inline]
fn conjugate(value: C64) -> C64 {
    C64::new(value.re, -value.im)
}

#[inline]
fn real_product(left: C64, right: C64) -> f64 {
    left.re * right.re - left.im * right.im
}

/// `integral_1^inf lambda^-s / sqrt(lambda^2-1) d lambda`
/// = `integral_0^(pi/2) cos(theta)^(s-1) d theta`.
fn zero_frequency_kernel(s: usize) -> f64 {
    let exponent = s - 1;
    let even = exponent.is_multiple_of(2);
    let mut value = if even { PI / 2.0 } else { 1.0 };
    let start = if even { 2 } else { 3 };
    for power in (start..=exponent).step_by(2) {
        value *= (power - 1) as f64 / power as f64;
    }
    value
}

fn oscillatory_kernel(
    s: usize,
    omega: f64,
    coarse_x: &[f64],
    coarse_w: &[f64],
    fine_x: &[f64],
    fine_w: &[f64],
) -> (C64, f64, usize) {
    if omega == 0.0 {
        return (C64::new(zero_frequency_kernel(s), 0.0), 0.0, 0);
    }
    let sign = omega.signum();
    let frequency = omega.abs();
    let coarse = steepest_descent_kernel(s, frequency, coarse_x, coarse_w);
    let fine = steepest_descent_kernel(s, frequency, fine_x, fine_w);
    let value = if sign > 0.0 { fine } else { conjugate(fine) };
    (value, (fine - coarse).abs(), coarse_x.len() + fine_x.len())
}

/// Complex Bickley kernel on the exact steepest-descent contour.
fn steepest_descent_kernel(s: usize, omega: f64, gx: &[f64], gw: &[f64]) -> C64 {
    const Y_MAX: f64 = 8.0;
    let half = Y_MAX / 2.0;
    let rotation = C64::new(FRAC_1_SQRT_2, FRAC_1_SQRT_2).scale(omega.sqrt().recip());
    let mut sum = C64::ZERO;
    for (&node, &weight) in gx.iter().zip(gw) {
        let y = half * (node + 1.0);
        let t_squared = C64::new(0.0, y * y / omega);
        let one_plus = C64::new(1.0, 0.0) + t_squared;
        let two_plus = C64::new(2.0, 0.0) + t_squared;
        let envelope = complex_pow(one_plus, s).recip() * complex_sqrt(two_plus).recip() * rotation;
        sum = sum + envelope.scale(weight * half * (-y * y).exp() * 2.0);
    }
    C64::cis(omega) * sum
}

fn complex_pow(mut base: C64, mut exponent: usize) -> C64 {
    let mut result = C64::new(1.0, 0.0);
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = result * base;
        }
        base = base * base;
        exponent >>= 1;
    }
    result
}

fn complex_sqrt(value: C64) -> C64 {
    let magnitude = value.abs();
    let re = ((magnitude + value.re) * 0.5).sqrt();
    let im = value.im.signum() * ((magnitude - value.re) * 0.5).sqrt();
    C64::new(re, im)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hulls;
    use crate::michell::InnerIntegral;
    use crate::BSplineSurface;
    use crate::STANDARD_GRAVITY;

    fn evaluate_terms(terms: &[EndpointTerm], nu: f64, lambda: f64, x_center: f64) -> C64 {
        terms.iter().fold(C64::ZERO, |sum, term| {
            let phase = C64::cis(nu * lambda * (term.x - x_center));
            let decay = (-nu * lambda * lambda * term.z).exp();
            sum + phase
                * term
                    .coeff
                    .scale(decay / lambda.powi(term.lambda_power as i32))
        })
    }

    #[test]
    fn endpoint_decomposition_reproduces_exact_inner_amplitude() {
        let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let fn_ = 0.05;
        let speed = fn_ * (STANDARD_GRAVITY * hull.length()).sqrt();
        let nu = STANDARD_GRAVITY / speed.powi(2);
        let terms = endpoint_terms(&hull, nu);
        let mut direct = InnerIntegral::new(&hull, nu);
        for lambda in [1.0, 1.2, 2.0, 5.0, 20.0] {
            let got = evaluate_terms(&terms, nu, lambda, hull.x_center());
            let want = direct.eval(lambda);
            let scale = want.abs().max(1e-14);
            assert!(
                (got - want).abs() <= 2e-10 * scale,
                "lambda={lambda}: endpoint {got:?}, direct {want:?}"
            );
        }
    }

    #[test]
    fn endpoint_decomposition_handles_multiple_spans_and_a_chine() {
        let surface = BSplineSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 0.8, 0.8, 2.0, 2.0, 2.0],
            vec![0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0],
            vec![
                0.0, 0.0, 0.0, 0.5, 0.4, 0.1, 0.7, 0.5, 0.1, 0.6, 0.45, 0.05, 0.3, 0.2, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        let hull = Hull::new(surface).unwrap();
        let nu = 30.0;
        let terms = endpoint_terms(&hull, nu);
        let mut direct = InnerIntegral::new(&hull, nu);
        for lambda in [1.0, 1.1, 1.7, 4.0, 12.0] {
            let got = evaluate_terms(&terms, nu, lambda, hull.x_center());
            let want = direct.eval(lambda);
            let scale = want.abs().max(1e-14);
            assert!(
                (got - want).abs() <= 5e-10 * scale,
                "lambda={lambda}: endpoint {got:?}, direct {want:?}"
            );
        }

        let speed = (STANDARD_GRAVITY / 40.0).sqrt();
        let reduced = low_froude_wave_resistance(&hull, &Conditions::freshwater(speed)).unwrap();
        assert!(reduced.endpoint_pairs.iter().any(|pair| {
            pair.left.kind == EndpointKind::InteriorKnot
                || pair.right.kind == EndpointKind::InteriorKnot
        }));
    }

    #[test]
    fn reported_error_charges_shared_coefficients_and_endpoint_accumulation() {
        let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let speed = 0.05 * (STANDARD_GRAVITY * hull.length()).sqrt();
        let reduced = low_froude_wave_resistance(&hull, &Conditions::freshwater(speed)).unwrap();
        assert!(reduced.coefficient_rel_error_bound > 0.0);
        assert!(reduced.endpoint_summation_abs_error_bound > 0.0);
        let accounted = reduced.coefficient_rel_error_bound
            + (reduced.omitted_abs_error_bound
                + reduced.quadrature_abs_error_estimate
                + reduced.endpoint_summation_abs_error_bound)
                / reduced.resistance.abs();
        assert_eq!(reduced.est_rel_error, accounted);
    }

    #[test]
    fn steepest_descent_kernel_matches_resolved_real_axis_reference() {
        let (gx, gw) = gauss_legendre(64);
        let (reference_x, reference_w) = gauss_legendre(16);
        for s in [4, 7, 10] {
            for omega in [25.0, 100.0, 400.0] {
                let got = steepest_descent_kernel(s, omega, &gx, &gw);
                // Independent real-axis quadrature. Panel boundaries are
                // uniform in t^2, so the phase advances by at most pi on each
                // panel. lambda_max=1000 leaves an absolute tail below
                // 3e-13 even for the slowest-decaying s=4 kernel.
                let lambda_max = 1_000.0;
                let du = PI / omega;
                let panels = ((lambda_max - 1.0) / du).ceil() as usize;
                let mut reference = C64::ZERO;
                for panel in 0..panels {
                    let u0 = panel as f64 * du;
                    let u1 = ((panel + 1) as f64 * du).min(lambda_max - 1.0);
                    let t0 = u0.sqrt();
                    let t1 = u1.sqrt();
                    let half = 0.5 * (t1 - t0);
                    let mid = 0.5 * (t0 + t1);
                    for (&node, &weight) in reference_x.iter().zip(&reference_w) {
                        let t = mid + half * node;
                        let lambda = 1.0 + t * t;
                        let envelope = 2.0 / lambda.powi(s as i32) / (2.0 + t * t).sqrt();
                        reference =
                            reference + C64::cis(omega * lambda).scale(weight * half * envelope);
                    }
                }
                let error = (got - reference).abs();
                assert!(
                    error <= 1e-9 * got.abs().max(1e-14),
                    "s={s} omega={omega}: got {got:?}, reference {reference:?}, error {error}"
                );
            }
        }
    }
}

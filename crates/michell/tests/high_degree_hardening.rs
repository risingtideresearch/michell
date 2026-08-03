//! Adversarial high-degree regressions from the coordinated R1/R2/R3 review.
//!
//! The reference paths deliberately share no `Hull::fx_coeff` code with the
//! production solvers. The original reviewer programs are preserved under
//! `probes/adversarial/` and documented there.

use michell::{
    wave_resistance_with, BSplineSurface, Conditions, Error, Hull, WaveOptions, WaveOutcome,
    STANDARD_GRAVITY,
};
use std::f64::consts::PI;

const DEFAULT_REL_TOL: f64 = 1.0e-5;
const R1_P2_INDEPENDENT_REFERENCE: f64 = 2.565_891_321_290e-3;
const R2_RIGOROUS_UPPER_BOUND_N: f64 = 4_989.990_8;
const R3_INDEPENDENT_REFERENCE_N: f64 = 6.382_528_496_934_648_5e-3;

fn clamped_knots(degree: usize, lower: f64, upper: f64) -> Vec<f64> {
    std::iter::repeat_n(lower, degree + 1)
        .chain(std::iter::repeat_n(upper, degree + 1))
        .collect()
}

fn r1_hull(degree_x: usize) -> michell::Result<Hull> {
    let degree_z = 2;
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let mut control = Vec::with_capacity((degree_x + 1) * (degree_z + 1));
    for i in 0..=degree_x {
        let u = i as f64 / degree_x as f64;
        let bx = 2.0 * beam * u * (1.0 - u);
        for j in 0..=degree_z {
            let v = j as f64 / degree_z as f64;
            control.push(bx * (1.0 - v * v));
        }
    }
    Hull::new(BSplineSurface::new(
        degree_x,
        degree_z,
        clamped_knots(degree_x, -0.5 * length, 0.5 * length),
        clamped_knots(degree_z, 0.0, draft),
        control,
    )?)
}

fn r2_hull(degree_x: usize) -> michell::Result<Hull> {
    let degree_z = 2;
    let mut control = vec![0.0; (degree_x + 1) * (degree_z + 1)];
    for i in 1..degree_x {
        let x = i as f64 / degree_x as f64;
        for j in 0..degree_z {
            let z = j as f64 / degree_z as f64;
            control[i * (degree_z + 1) + j] = 4.0 * x * (1.0 - x) * (1.0 - z * z);
        }
    }
    Hull::new(BSplineSurface::new(
        degree_x,
        degree_z,
        clamped_knots(degree_x, -5.0, 5.0),
        clamped_knots(degree_z, 0.0, 5.0),
        control,
    )?)
}

fn r3_hull(degree_x: usize, degree_z: usize) -> michell::Result<Hull> {
    let mut control = Vec::with_capacity((degree_x + 1) * (degree_z + 1));
    for i in 0..=degree_x {
        let x = (PI * i as f64 / degree_x as f64).sin();
        for j in 0..=degree_z {
            let z = 1.0 - j as f64 / degree_z as f64;
            control.push(0.5 * x * z);
        }
    }
    Hull::new(BSplineSurface::new(
        degree_x,
        degree_z,
        clamped_knots(degree_x, -5.0, 5.0),
        clamped_knots(degree_z, 0.0, 2.0),
        control,
    )?)
}

fn speed_at_froude(hull: &Hull, froude: f64) -> f64 {
    froude * (STANDARD_GRAVITY * hull.length()).sqrt()
}

fn accepted_or_explicitly_unsupported(
    hull: michell::Result<Hull>,
    froude: f64,
) -> Option<michell::WaveResistance> {
    let hull = match hull {
        Ok(hull) => hull,
        Err(Error::Unsupported(_)) => return None,
        Err(error) => panic!("unexpected construction error: {error}"),
    };
    let conditions = Conditions::freshwater(speed_at_froude(&hull, froude));
    match wave_resistance_with(&hull, &conditions, &WaveOptions::default()) {
        Ok(result) => Some(result),
        Err(Error::Unsupported(_)) => None,
        Err(error) => panic!("unexpected resistance error: {error}"),
    }
}

fn gauss_legendre(order: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; order];
    let mut weights = vec![0.0; order];
    for i in 0..order.div_ceil(2) {
        let mut root = (PI * (i as f64 + 0.75) / (order as f64 + 0.5)).cos();
        loop {
            let mut p0 = 1.0;
            let mut p1 = root;
            for degree in 2..=order {
                let p2 = ((2 * degree - 1) as f64 * root * p1 - (degree - 1) as f64 * p0)
                    / degree as f64;
                p0 = p1;
                p1 = p2;
            }
            let derivative = order as f64 * (root * p1 - p0) / (root * root - 1.0);
            let next = root - p1 / derivative;
            if (next - root).abs() < 2.0e-16 {
                root = next;
                let weight = 2.0 / ((1.0 - root * root) * derivative * derivative);
                nodes[i] = -root;
                nodes[order - 1 - i] = root;
                weights[i] = weight;
                weights[order - 1 - i] = weight;
                break;
            }
            root = next;
        }
    }
    (nodes, weights)
}

/// Closed-form transform for the R1 single-span Bernstein surface. This is
/// the independent amplitude derived in the R1 `analytic_ref.rs` probe, not a
/// reconstruction from production coefficients.
fn r1_exact_inner_squared(lambda: f64, degree_x: usize, nu: f64) -> f64 {
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let half_length = 0.5 * length;
    let k = nu * lambda;
    let a = nu * lambda * lambda * draft;
    let decay = (-a).exp();
    let i0 = (1.0 - decay) / a;
    let i1 = (1.0 - decay * (1.0 + a)) / (a * a);
    let i2 = (2.0 - decay * (a * a + 2.0 * a + 2.0)) / (a * a * a);
    let z_transform = draft * (i0 - 0.5 * i1 - 0.5 * i2);
    let x_transform =
        2.0 * ((k * half_length).sin() / (k * k) - half_length * (k * half_length).cos() / k);
    let amplitude = (-4.0 * beam * (1.0 - 1.0 / degree_x as f64) / (length * length))
        * x_transform
        * z_transform;
    amplitude * amplitude
}

/// Independent, phase-resolved Gauss-Legendre quadrature from R1. A cutoff
/// of 1000 is already stable to far below the requested 1e-5 relative error.
fn r1_independent_reference(degree_x: usize, froude: f64) -> f64 {
    let length = 10.0;
    let speed = froude * (STANDARD_GRAVITY * length).sqrt();
    let nu = STANDARD_GRAVITY / (speed * speed);
    let phase_panel = PI / (nu * length);
    let lambda_max = 1_000.0;
    let panel_count = ((lambda_max - 1.0) / phase_panel).ceil() as usize;
    let (nodes, weights) = gauss_legendre(32);
    let mut integral = 0.0;
    let mut correction = 0.0;
    for panel in 0..panel_count {
        let u0 = panel as f64 * phase_panel;
        let u1 = ((panel + 1) as f64 * phase_panel).min(lambda_max - 1.0);
        let t0 = u0.sqrt();
        let t1 = u1.sqrt();
        let half = 0.5 * (t1 - t0);
        let mid = 0.5 * (t1 + t0);
        for (&node, &weight) in nodes.iter().zip(&weights) {
            let t = mid + half * node;
            let lambda = 1.0 + t * t;
            let term = weight * half * 2.0 * lambda * lambda / (2.0 + t * t).sqrt()
                * r1_exact_inner_squared(lambda, degree_x, nu);
            let adjusted = term - correction;
            let next = integral + adjusted;
            correction = (next - integral) - adjusted;
            integral = next;
        }
    }
    4.0 * 999.1 * STANDARD_GRAVITY * STANDARD_GRAVITY / (PI * speed * speed) * integral
}

#[test]
fn r1_degree_48_matches_independent_closed_form_quadrature_or_refuses() {
    let Some(result) = accepted_or_explicitly_unsupported(r1_hull(48), 0.05) else {
        return;
    };
    let reference = r1_independent_reference(48, 0.05);
    let relative_error = (result.resistance - reference).abs() / reference;
    assert!(result.resistance.is_finite() && result.est_rel_error.is_finite());
    assert!(
        relative_error <= DEFAULT_REL_TOL,
        "accepted p=48 result has relative error {relative_error:e} against independent R1 reference {reference:e}; reported {:e}",
        result.est_rel_error
    );
}

#[test]
fn r1_degree_sweep_has_no_silently_accepted_tolerance_violation() {
    let mut accepted = 0;
    for degree_x in 2..=48 {
        let Some(result) = accepted_or_explicitly_unsupported(r1_hull(degree_x), 0.05) else {
            continue;
        };
        accepted += 1;
        assert!(result.resistance.is_finite() && result.est_rel_error.is_finite());
        if result.outcome == WaveOutcome::Converged {
            // All members of this R1 family have the same exact quadratic
            // geometry up to this analytic Bernstein-control scaling.
            let scale = 2.0 * (1.0 - 1.0 / degree_x as f64);
            let reference = R1_P2_INDEPENDENT_REFERENCE * scale * scale;
            let relative_error = (result.resistance - reference).abs() / reference;
            assert!(
                relative_error <= DEFAULT_REL_TOL,
                "first accepted tolerance violation at longitudinal degree {degree_x}: actual={relative_error:e}, reported={:e}",
                result.est_rel_error
            );
        }
    }
    assert!(
        accepted > 0,
        "the implementation must retain a validated degree envelope"
    );
}

#[test]
fn r2_degree_48_respects_rigorous_variation_bound_or_refuses() {
    let Some(result) = accepted_or_explicitly_unsupported(r2_hull(48), 0.10) else {
        return;
    };
    assert!(result.resistance.is_finite() && result.est_rel_error.is_finite());
    assert!(
        result.resistance <= R2_RIGOROUS_UPPER_BOUND_N,
        "accepted R2 p=48 result {} N exceeds the independent variation bound {} N",
        result.resistance,
        R2_RIGOROUS_UPPER_BOUND_N
    );
}

#[test]
fn r2_degree_192_never_returns_ok_nan() {
    let Some(result) = accepted_or_explicitly_unsupported(r2_hull(192), 0.10) else {
        return;
    };
    assert!(
        result.resistance.is_finite() && result.est_rel_error.is_finite(),
        "accepted degree-192 result contains a non-finite value: {result:?}"
    );
}

#[test]
fn r3_degree_24x16_meets_requested_error_or_refuses() {
    let Some(result) = accepted_or_explicitly_unsupported(r3_hull(24, 16), 0.05) else {
        return;
    };
    assert!(result.resistance.is_finite() && result.est_rel_error.is_finite());
    if result.outcome == WaveOutcome::Converged {
        let relative_error =
            (result.resistance - R3_INDEPENDENT_REFERENCE_N).abs() / R3_INDEPENDENT_REFERENCE_N;
        assert!(
            relative_error <= DEFAULT_REL_TOL,
            "accepted R3 degree-24x16 result has relative error {relative_error:e}; reported {:e}",
            result.est_rel_error
        );
    }
}

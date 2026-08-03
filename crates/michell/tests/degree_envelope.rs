//! Independent validation of the declared spline-degree envelope.
//!
//! Bézier degree elevation is an exact geometry operation. These tests build
//! each elevated control net directly, without using the production
//! coefficient extraction, so resistance changes expose representation error.

use michell::{
    wave_resistance_with, BSplineSurface, Conditions, Error, Hull, WaveOptions, WaveOutcome,
    MAX_SUPPORTED_SPLINE_DEGREE, STANDARD_GRAVITY,
};

const REQUESTED_REL_TOL: f64 = 1.0e-5;
const R3_DEGREE_16X12_REFERENCE_N: f64 = 6.337_513_610_463_49e-3;

#[derive(Clone)]
struct BezierSurface {
    degree_x: usize,
    degree_z: usize,
    control: Vec<f64>,
}

impl BezierSurface {
    fn wigley() -> Self {
        let gx = [0.0, 1.0, 0.0];
        let gz = [1.0, 1.0, 0.0];
        Self {
            degree_x: 2,
            degree_z: 2,
            control: gx
                .into_iter()
                .flat_map(|x| gz.into_iter().map(move |z| x * z))
                .collect(),
        }
    }

    fn elevate_x(&self) -> Self {
        let old_nx = self.degree_x + 1;
        let nz = self.degree_z + 1;
        let new_degree = self.degree_x + 1;
        let mut control = vec![0.0; (new_degree + 1) * nz];
        for j in 0..nz {
            control[j] = self.control[j];
            control[new_degree * nz + j] = self.control[(old_nx - 1) * nz + j];
            for i in 1..new_degree {
                let alpha = i as f64 / new_degree as f64;
                control[i * nz + j] = alpha * self.control[(i - 1) * nz + j]
                    + (1.0 - alpha) * self.control[i * nz + j];
            }
        }
        Self {
            degree_x: new_degree,
            degree_z: self.degree_z,
            control,
        }
    }

    fn elevate_z(&self) -> Self {
        let nx = self.degree_x + 1;
        let old_nz = self.degree_z + 1;
        let new_degree = self.degree_z + 1;
        let new_nz = new_degree + 1;
        let mut control = vec![0.0; nx * new_nz];
        for i in 0..nx {
            control[i * new_nz] = self.control[i * old_nz];
            control[i * new_nz + new_degree] = self.control[i * old_nz + old_nz - 1];
            for j in 1..new_degree {
                let alpha = j as f64 / new_degree as f64;
                control[i * new_nz + j] = alpha * self.control[i * old_nz + j - 1]
                    + (1.0 - alpha) * self.control[i * old_nz + j];
            }
        }
        Self {
            degree_x: self.degree_x,
            degree_z: new_degree,
            control,
        }
    }

    fn hull(&self) -> michell::Result<Hull> {
        Hull::new(BSplineSurface::new(
            self.degree_x,
            self.degree_z,
            clamped(self.degree_x, -5.0, 5.0),
            clamped(self.degree_z, 0.0, 0.625),
            self.control.clone(),
        )?)
    }
}

fn clamped(degree: usize, lower: f64, upper: f64) -> Vec<f64> {
    std::iter::repeat_n(lower, degree + 1)
        .chain(std::iter::repeat_n(upper, degree + 1))
        .collect()
}

fn evaluate(hull: &Hull, froude: f64) -> michell::WaveResistance {
    let speed = froude * (STANDARD_GRAVITY * hull.length()).sqrt();
    wave_resistance_with(
        hull,
        &Conditions::freshwater(speed),
        &WaveOptions::default(),
    )
    .unwrap()
}

fn assert_estimate_covers(label: &str, result: michell::WaveResistance, reference: f64) {
    assert_eq!(result.outcome, WaveOutcome::Converged, "{label}");
    let actual_rel_error = (result.resistance - reference).abs() / reference;
    eprintln!(
        "{label}: resistance={:.15e} actual_rel={actual_rel_error:.6e} reported={:.6e} method={:?}",
        result.resistance, result.est_rel_error, result.method
    );
    assert!(
        actual_rel_error <= result.est_rel_error,
        "{label}: actual relative error {actual_rel_error:e} exceeds reported {:e}",
        result.est_rel_error
    );
    assert!(
        actual_rel_error <= REQUESTED_REL_TOL,
        "{label}: accepted relative error {actual_rel_error:e} exceeds requested {REQUESTED_REL_TOL:e}"
    );
}

#[test]
fn exact_wigley_degree_elevation_is_covered_through_the_supported_envelope() {
    let base = BezierSurface::wigley();
    let reference = evaluate(&base.hull().unwrap(), 0.05).resistance;

    let mut elevated_x = base.clone();
    for degree_x in 2..=MAX_SUPPORTED_SPLINE_DEGREE {
        assert_eq!(elevated_x.degree_x, degree_x);
        let result = evaluate(&elevated_x.hull().unwrap(), 0.05);
        assert_estimate_covers(&format!("degree ({degree_x}, 2)"), result, reference);
        elevated_x = elevated_x.elevate_x();
    }

    let mut elevated_z = base.clone();
    for degree_z in 2..=MAX_SUPPORTED_SPLINE_DEGREE {
        assert_eq!(elevated_z.degree_z, degree_z);
        let result = evaluate(&elevated_z.hull().unwrap(), 0.05);
        assert_estimate_covers(&format!("degree (2, {degree_z})"), result, reference);
        elevated_z = elevated_z.elevate_z();
    }

    // Exercise the tensor corner too: this reaches the maximum endpoint-kernel
    // order generated by the supported degree pair (16, 16).
    let mut elevated_both = base;
    for degree in 2..=MAX_SUPPORTED_SPLINE_DEGREE {
        assert_eq!(
            (elevated_both.degree_x, elevated_both.degree_z),
            (degree, degree)
        );
        let result = evaluate(&elevated_both.hull().unwrap(), 0.05);
        assert_estimate_covers(&format!("degree ({degree}, {degree})"), result, reference);
        elevated_both = elevated_both.elevate_x().elevate_z();
    }
}

#[test]
fn general_marcher_degree_elevation_is_covered_at_the_tensor_corner() {
    let mut elevated = BezierSurface::wigley();
    let reference = evaluate(&elevated.hull().unwrap(), 0.35).resistance;
    for degree in 2..=MAX_SUPPORTED_SPLINE_DEGREE {
        let result = evaluate(&elevated.hull().unwrap(), 0.35);
        assert_eq!(result.method, michell::WaveMethod::GeneralMarcher);
        assert_estimate_covers(
            &format!("general degree ({degree}, {degree})"),
            result,
            reference,
        );
        elevated = elevated.elevate_x().elevate_z();
    }
}

#[test]
fn degree_above_the_validated_envelope_is_explicitly_unsupported() {
    let mut elevated_x = BezierSurface::wigley();
    while elevated_x.degree_x <= MAX_SUPPORTED_SPLINE_DEGREE {
        elevated_x = elevated_x.elevate_x();
    }
    assert!(matches!(elevated_x.hull(), Err(Error::Unsupported(_))));

    let mut elevated_z = BezierSurface::wigley();
    while elevated_z.degree_z <= MAX_SUPPORTED_SPLINE_DEGREE {
        elevated_z = elevated_z.elevate_z();
    }
    assert!(matches!(elevated_z.hull(), Err(Error::Unsupported(_))));
}

#[test]
fn independent_degree_16x12_reference_is_covered_by_the_reported_error() {
    let degree_x = 16;
    let degree_z = 12;
    let mut control = Vec::with_capacity((degree_x + 1) * (degree_z + 1));
    for i in 0..=degree_x {
        let x = (std::f64::consts::PI * i as f64 / degree_x as f64).sin();
        for j in 0..=degree_z {
            let z = 1.0 - j as f64 / degree_z as f64;
            control.push(0.5 * x * z);
        }
    }
    let hull = Hull::new(
        BSplineSurface::new(
            degree_x,
            degree_z,
            clamped(degree_x, -5.0, 5.0),
            clamped(degree_z, 0.0, 2.0),
            control,
        )
        .unwrap(),
    )
    .unwrap();
    let result = evaluate(&hull, 0.05);
    assert_estimate_covers(
        "independent R3 degree (16, 12)",
        result,
        R3_DEGREE_16X12_REFERENCE_N,
    );
}

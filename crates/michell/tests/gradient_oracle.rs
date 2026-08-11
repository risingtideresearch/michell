//! Independent polynomial-basis oracle for wave resistance and its control
//! derivatives. This deliberately avoids the production B-spline span
//! extraction, moment recurrences, outer quadrature, and reverse pass.

use michell::{
    hulls, wave_resistance_gradient_with, BSplineSurface, Conditions, Hull, WaveOptions,
    STANDARD_GRAVITY,
};

const LENGTH: f64 = 10.0;
const BEAM: f64 = 1.0;
const DRAFT: f64 = 0.625;
const LAMBDA_MAX: f64 = 200.0;
const FINE_INTERVALS: usize = 2_000_000;
const TENT_SPLIT: f64 = -1.5;
const WIGLEY_CONTROLS: [f64; 9] = [0.0, 0.0, 0.0, BEAM, BEAM, 0.0, 0.0, 0.0, 0.0];
const TENT_CONTROLS: [f64; 9] = [0.0, 0.0, 0.0, BEAM, 0.9, 0.0, 0.0, 0.0, 0.0];

fn quadratic_x_basis(k: f64) -> [(f64, f64); 3] {
    let a = LENGTH / 2.0;
    let ka = k * a;
    // M0 = integral exp(ikx) dx and M1 = integral x exp(ikx) dx = i*m1
    // on [-L/2, L/2].
    let m0 = 2.0 * ka.sin() / k;
    let m1 = 2.0 * (ka.sin() / (k * k) - a * ka.cos() / k);
    let l2 = LENGTH * LENGTH;
    [
        (-2.0 * a * m0 / l2, 2.0 * m1 / l2),
        (0.0, -4.0 * m1 / l2),
        (2.0 * a * m0 / l2, 2.0 * m1 / l2),
    ]
}

fn exp_interval(k: f64, start: f64, end: f64) -> (f64, f64) {
    (
        ((k * end).sin() - (k * start).sin()) / k,
        ((k * start).cos() - (k * end).cos()) / k,
    )
}

fn tent_x_basis(k: f64) -> [(f64, f64); 3] {
    let a = LENGTH / 2.0;
    let left = exp_interval(k, -a, TENT_SPLIT);
    let right = exp_interval(k, TENT_SPLIT, a);
    let left_len = TENT_SPLIT + a;
    let right_len = a - TENT_SPLIT;
    [
        (-left.0 / left_len, -left.1 / left_len),
        (
            left.0 / left_len - right.0 / right_len,
            left.1 / left_len - right.1 / right_len,
        ),
        (right.0 / right_len, right.1 / right_len),
    ]
}

/// Resistance integrand followed by its nine control derivatives for the
/// supplied three-function longitudinal basis and quadratic vertical basis.
fn basis_integrands(
    theta: f64,
    nu: f64,
    controls: &[f64; 9],
    longitudinal_basis: fn(f64) -> [(f64, f64); 3],
) -> [f64; 10] {
    let sec = 1.0 / theta.cos();
    let k = nu * sec;
    let kappa = nu * sec * sec;
    let x_basis = longitudinal_basis(k);

    // E_m = integral z^m exp(-kappa*z) dz on [0, T], m = 0, 1, 2.
    let q = kappa * DRAFT;
    let decay = (-q).exp();
    let e0 = (1.0 - decay) / kappa;
    let e1 = (1.0 - decay * (1.0 + q)) / (kappa * kappa);
    let e2 = (2.0 - decay * (2.0 + 2.0 * q + q * q)) / kappa.powi(3);
    let z_basis = [
        e0 - 2.0 * e1 / DRAFT + e2 / (DRAFT * DRAFT),
        2.0 * e1 / DRAFT - 2.0 * e2 / (DRAFT * DRAFT),
        e2 / (DRAFT * DRAFT),
    ];

    let mut basis = [(0.0, 0.0); 9];
    for i in 0..3 {
        for (j, &z) in z_basis.iter().enumerate() {
            basis[3 * i + j] = (x_basis[i].0 * z, x_basis[i].1 * z);
        }
    }

    let amplitude = controls
        .iter()
        .zip(&basis)
        .fold((0.0, 0.0), |sum, (&control, &(re, im))| {
            (sum.0 + control * re, sum.1 + control * im)
        });
    let weight = sec * sec * sec;
    let mut values = [0.0; 10];
    values[0] = (amplitude.0 * amplitude.0 + amplitude.1 * amplitude.1) * weight;
    for (index, &(re, im)) in basis.iter().enumerate() {
        values[index + 1] = 2.0 * (amplitude.0 * re + amplitude.1 * im) * weight;
    }
    values
}

fn simpson_oracle(
    intervals: usize,
    lambda_max: f64,
    conditions: &Conditions,
    controls: &[f64; 9],
    longitudinal_basis: fn(f64) -> [(f64, f64); 3],
) -> [f64; 10] {
    assert_eq!(intervals % 2, 0);
    let nu = conditions.gravity / conditions.speed.powi(2);
    let theta_max = (1.0 / lambda_max).acos();
    let step = theta_max / intervals as f64;
    let mut sums = basis_integrands(0.0, nu, controls, longitudinal_basis);
    let endpoint = basis_integrands(theta_max, nu, controls, longitudinal_basis);
    for i in 0..10 {
        sums[i] += endpoint[i];
    }
    for sample in 1..intervals {
        let values = basis_integrands(sample as f64 * step, nu, controls, longitudinal_basis);
        let weight = if sample % 2 == 0 { 2.0 } else { 4.0 };
        for i in 0..10 {
            sums[i] += weight * values[i];
        }
    }
    let prefactor = 4.0 * conditions.fluid.density * conditions.gravity.powi(2)
        / (std::f64::consts::PI * conditions.speed.powi(2));
    for value in &mut sums {
        *value *= prefactor * step / 3.0;
    }
    sums
}

fn tent_hull() -> Hull {
    let a = LENGTH / 2.0;
    Hull::new(
        BSplineSurface::new(
            1,
            2,
            vec![-a, -a, TENT_SPLIT, a, a],
            vec![0.0, 0.0, 0.0, DRAFT, DRAFT, DRAFT],
            TENT_CONTROLS.to_vec(),
        )
        .unwrap(),
    )
    .unwrap()
}

fn check_oracle(
    case: &str,
    hull: &Hull,
    controls: &[f64; 9],
    conditions: &Conditions,
    longitudinal_basis: fn(f64) -> [(f64, f64); 3],
    extend_reference_tail: bool,
) {
    let production = wave_resistance_gradient_with(
        hull,
        conditions,
        &WaveOptions {
            rel_tol: 1e-9,
            max_refinements: 7,
        },
    )
    .unwrap();
    // The Wigley reference extends beyond production's truncation point to
    // preserve the existing continuous-integral check. The asymmetric case
    // uses the production domain to isolate amplitude and control-map errors.
    let lambda_max = if extend_reference_tail {
        LAMBDA_MAX
    } else {
        production.wave.max_lambda
    };
    if extend_reference_tail {
        assert!(
            lambda_max > 2.0 * production.wave.max_lambda,
            "{case} reference must extend well beyond production: {lambda_max} vs {}",
            production.wave.max_lambda
        );
    }
    let coarse = simpson_oracle(
        FINE_INTERVALS / 2,
        lambda_max,
        conditions,
        controls,
        longitudinal_basis,
    );
    let fine = simpson_oracle(
        FINE_INTERVALS,
        lambda_max,
        conditions,
        controls,
        longitudinal_basis,
    );

    for index in 0..10 {
        let convergence = (fine[index] - coarse[index]).abs() / fine[index].abs().max(1.0);
        assert!(
            convergence < 1e-8,
            "{case} oracle component {index} did not self-converge: {convergence:.3e}"
        );
    }

    let mut production_values = [0.0; 10];
    production_values[0] = production.wave.resistance;
    production_values[1..].copy_from_slice(&production.control_gradient);
    for index in 0..10 {
        let relative = (production_values[index] - fine[index]).abs() / fine[index].abs().max(1.0);
        assert!(
            relative < 1e-6,
            "{case} component {index}: production={}, oracle={}, scaled error={relative:.3e}",
            production_values[index],
            fine[index]
        );
    }
}

fn check_quadratic_x_basis(nu: f64) {
    const INTERVALS: usize = 100_000;
    let k = 1.37 * nu;
    let a = LENGTH / 2.0;
    let step = LENGTH / INTERVALS as f64;
    let l2 = LENGTH * LENGTH;
    let mut reference = [(0.0, 0.0); 3];
    for sample in 0..=INTERVALS {
        let x = -a + sample as f64 * step;
        let weight = if sample == 0 || sample == INTERVALS {
            1.0
        } else if sample % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let derivatives = [2.0 * (x - a) / l2, -4.0 * x / l2, 2.0 * (x + a) / l2];
        let (sin, cos) = (k * x).sin_cos();
        for i in 0..3 {
            reference[i].0 += weight * derivatives[i] * cos;
            reference[i].1 += weight * derivatives[i] * sin;
        }
    }
    for value in &mut reference {
        value.0 *= step / 3.0;
        value.1 *= step / 3.0;
    }
    for (index, (got, want)) in quadratic_x_basis(k).iter().zip(reference).enumerate() {
        let error = (got.0 - want.0).hypot(got.1 - want.1) / want.0.hypot(want.1).max(1.0);
        assert!(
            error < 1e-10,
            "quadratic x basis {index}: error={error:.3e}"
        );
    }
}

#[test]
fn gradient_matches_independent_basis_oracle() {
    let speed = 0.35 * (STANDARD_GRAVITY * LENGTH).sqrt();
    let conditions = Conditions::freshwater(speed);
    check_quadratic_x_basis(conditions.gravity / conditions.speed.powi(2));
    let wigley = hulls::wigley(LENGTH, BEAM, DRAFT).unwrap();
    check_oracle(
        "Wigley",
        &wigley,
        &WIGLEY_CONTROLS,
        &conditions,
        quadratic_x_basis,
        true,
    );

    // The off-centre internal knot breaks fore-aft symmetry while the zero end
    // rows and keel column keep the hull closed at both ends and the bottom.
    let tent = tent_hull();
    let a = LENGTH / 2.0;
    for z in [0.0, DRAFT / 2.0, DRAFT] {
        assert_eq!(tent.surface().eval(-a, z), 0.0);
        assert_eq!(tent.surface().eval(a, z), 0.0);
    }
    for x in [-a, TENT_SPLIT, 0.0, a] {
        assert_eq!(tent.surface().eval(x, DRAFT), 0.0);
    }
    let peak = tent.surface().eval(TENT_SPLIT, 0.0);
    let mirror = tent.surface().eval(-TENT_SPLIT, 0.0);
    assert!(peak > 0.5 * BEAM);
    assert!((peak - mirror).abs() > 0.1 * BEAM);
    check_oracle(
        "closed off-centre tent",
        &tent,
        &TENT_CONTROLS,
        &conditions,
        tent_x_basis,
        false,
    );
}

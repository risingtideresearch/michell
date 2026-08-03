use michell::{
    hulls, low_froude_wave_resistance, multihull_wave_resistance_with, wave_resistance_with,
    Conditions, EndpointKind, Placement, WaveMethod, WaveOptions, WaveOutcome, WaveResistance,
    STANDARD_GRAVITY,
};

fn marching_wave_resistance(
    hull: &michell::Hull,
    conditions: &Conditions,
    options: &WaveOptions,
) -> WaveResistance {
    multihull_wave_resistance_with(&[(hull, Placement::default())], conditions, options).unwrap()
}

fn wigley_j(length: f64, beam: f64, draft: f64, nu: f64, lambda: f64) -> f64 {
    let half_length = length / 2.0;
    let k = nu * lambda;
    let kappa = nu * lambda * lambda;
    let x = 2.0 * ((k * half_length).sin() / k.powi(2) - half_length * (k * half_length).cos() / k);
    let decay = (-kappa * draft).exp();
    let z = (1.0 - decay) / kappa
        - (2.0 - decay * (kappa.powi(2) * draft.powi(2) + 2.0 * kappa * draft + 2.0))
            / (kappa.powi(3) * draft.powi(2));
    -(4.0 * beam / length.powi(2)) * x * z
}

/// Independent positive-real-axis reference: analytic Wigley inner amplitude,
/// 16-point Gauss panels whose boundaries advance the bow/stern phase by at
/// most pi/2, and lambda_max=500 (leading-tail fraction below 2e-11).
fn resolved_wigley_reference(
    length: f64,
    beam: f64,
    draft: f64,
    conditions: &Conditions,
) -> (f64, usize) {
    let nu = conditions.gravity / conditions.speed.powi(2);
    let half_length = length / 2.0;
    let lambda_max = 500.0;
    let du = std::f64::consts::FRAC_PI_2 / (nu * half_length);
    let panels = ((lambda_max - 1.0) / du).ceil() as usize;
    const NODES: [f64; 8] = [
        0.095_012_509_837_637_44,
        0.281_603_550_779_258_9,
        0.458_016_777_657_227_4,
        0.617_876_244_402_643_8,
        0.755_404_408_355_003,
        0.865_631_202_387_831_8,
        0.944_575_023_073_232_6,
        0.989_400_934_991_649_9,
    ];
    const WEIGHTS: [f64; 8] = [
        0.189_450_610_455_068_5,
        0.182_603_415_044_923_6,
        0.169_156_519_395_002_54,
        0.149_595_988_816_576_73,
        0.124_628_971_255_533_87,
        0.095_158_511_682_492_78,
        0.062_253_523_938_647_89,
        0.027_152_459_411_754_095,
    ];
    let mut integral = 0.0;
    for panel in 0..panels {
        let u0 = panel as f64 * du;
        let u1 = ((panel + 1) as f64 * du).min(lambda_max - 1.0);
        let t0 = u0.sqrt();
        let t1 = u1.sqrt();
        let half = 0.5 * (t1 - t0);
        let mid = 0.5 * (t0 + t1);
        for (&node, &weight) in NODES.iter().zip(&WEIGHTS) {
            for sign in [-1.0, 1.0] {
                let t = mid + sign * half * node;
                let lambda = 1.0 + t * t;
                let j = wigley_j(length, beam, draft, nu, lambda);
                integral += weight * half * j * j * 2.0 * lambda.powi(2) / (2.0 + t * t).sqrt();
            }
        }
    }
    let coeff = 4.0 * conditions.fluid.density * conditions.gravity.powi(2)
        / (std::f64::consts::PI * conditions.speed.powi(2));
    (coeff * integral, panels * 16)
}

#[test]
fn endpoint_reduction_converges_as_froude_number_falls() {
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(length, beam, draft).unwrap();
    let reference_options = WaveOptions {
        rel_tol: 1e-8,
        max_refinements: 6,
    };
    for fn_ in [0.08, 0.05, 0.03, 0.02] {
        let speed = fn_ * (STANDARD_GRAVITY * hull.length()).sqrt();
        let conditions = Conditions::freshwater(speed);
        let direct = marching_wave_resistance(&hull, &conditions, &reference_options);
        let (reference, reference_evaluations) =
            resolved_wigley_reference(length, beam, draft, &conditions);
        let reduced = low_froude_wave_resistance(&hull, &conditions).unwrap();
        let actual_relative_error = (reduced.resistance - reference).abs() / reference;
        let direct_relative_error = (direct.resistance - reference).abs() / reference;
        eprintln!(
            "Fn={fn_:.2}: reference={reference:.12e}, reduced={:.12e}, reduced_error={actual_relative_error:.3e}, reduced_est={:.3e}, marching_error={direct_relative_error:.3e}, marching_est={:.3e}, marching_evals={}, kernel_evals={}, reference_evals={reference_evaluations}, terms={}/{}",
            reduced.resistance,
            reduced.est_rel_error,
            direct.est_rel_error,
            direct.inner_evaluations,
            reduced.kernel_evaluations,
            reduced.waterline_terms,
            reduced.endpoint_terms,
        );
        assert!(actual_relative_error <= 1.05 * reduced.est_rel_error.max(2e-10));
        if fn_ <= 0.05 {
            assert!(reduced.kernel_evaluations * 100 < direct.inner_evaluations);
            assert!(actual_relative_error < 2e-10);
            assert!(actual_relative_error * 100.0 < direct_relative_error);
        }
    }
}

#[test]
fn low_froude_reported_error_covers_actual_error() {
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(length, beam, draft).unwrap();
    let fn_ = 0.05;
    let speed = fn_ * (STANDARD_GRAVITY * length).sqrt();
    let conditions = Conditions::freshwater(speed);
    let result = wave_resistance_with(&hull, &conditions, &WaveOptions::default()).unwrap();
    let (reference, _) = resolved_wigley_reference(length, beam, draft, &conditions);
    let actual_relative_error = (result.resistance - reference).abs() / reference;
    assert!(
        actual_relative_error <= 10.0 * result.est_rel_error.max(1e-12),
        "reported {:.3e}, actual {actual_relative_error:.3e}, evaluations={}, max_lambda={:.3}",
        result.est_rel_error,
        result.inner_evaluations,
        result.max_lambda,
    );
}

#[test]
fn default_solver_accepts_only_a_reduction_within_tolerance() {
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(length, beam, draft).unwrap();

    for (fn_, should_use_reduction) in [(0.08, false), (0.05, true)] {
        let speed = fn_ * (STANDARD_GRAVITY * length).sqrt();
        let conditions = Conditions::freshwater(speed);
        let result = wave_resistance_with(&hull, &conditions, &WaveOptions::default()).unwrap();
        assert_eq!(
            result.method == WaveMethod::EndpointReduction,
            should_use_reduction,
            "unexpected default route at Fn={fn_}: estimated error {:.3e}",
            result.est_rel_error,
        );
        assert_eq!(result.outcome, WaveOutcome::Converged);
    }
}

#[test]
fn endpoint_pair_contributions_sum_to_reduced_resistance() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let speed = 0.05 * (STANDARD_GRAVITY * hull.length()).sqrt();
    let result = low_froude_wave_resistance(&hull, &Conditions::freshwater(speed)).unwrap();

    assert_eq!(
        result.endpoint_pairs.len(),
        result.waterline_terms * (result.waterline_terms + 1) / 2,
    );
    let resistance_sum: f64 = result
        .endpoint_pairs
        .iter()
        .map(|pair| pair.resistance)
        .sum();
    let fraction_sum: f64 = result
        .endpoint_pairs
        .iter()
        .map(|pair| pair.resistance_fraction)
        .sum();
    let quadrature_error_sum: f64 = result
        .endpoint_pairs
        .iter()
        .map(|pair| pair.quadrature_abs_error_estimate)
        .sum();
    assert!((resistance_sum - result.resistance).abs() <= 2e-13 * result.resistance);
    assert!((fraction_sum - 1.0).abs() <= 2e-13);
    assert!(
        (quadrature_error_sum - result.quadrature_abs_error_estimate).abs()
            <= 2e-13 * result.quadrature_abs_error_estimate.max(1e-300)
    );
    assert!(result.endpoint_pairs.iter().any(|pair| {
        pair.left.kind == EndpointKind::Bow || pair.right.kind == EndpointKind::Bow
    }));
    assert!(result.endpoint_pairs.iter().any(|pair| {
        pair.left.kind == EndpointKind::Stern || pair.right.kind == EndpointKind::Stern
    }));
}

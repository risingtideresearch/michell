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

fn gauss_legendre(order: usize) -> Vec<(f64, f64)> {
    let mut rule = Vec::with_capacity(order);
    for index in 0..order.div_ceil(2) {
        let mut root = (std::f64::consts::PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
        let (root, derivative) = loop {
            let mut previous = 1.0;
            let mut current = root;
            for degree in 2..=order {
                let next = ((2 * degree - 1) as f64 * root * current
                    - (degree - 1) as f64 * previous)
                    / degree as f64;
                previous = current;
                current = next;
            }
            let derivative = order as f64 * (root * current - previous) / (root * root - 1.0);
            let next_root = root - current / derivative;
            if (next_root - root).abs() <= 4.0 * f64::EPSILON {
                break (next_root, derivative);
            }
            root = next_root;
        };
        let weight = 2.0 / ((1.0 - root * root) * derivative * derivative);
        rule.push((-root, weight));
        if rule.len() < order {
            rule.push((root, weight));
        }
    }
    rule
}

fn add_kahan(sum: &mut f64, correction: &mut f64, value: f64) {
    let corrected = value - *correction;
    let next = *sum + corrected;
    *correction = (next - *sum) - corrected;
    *sum = next;
}

fn phase_resolved_panel_count(length: f64, nu: f64, lambda_max: f64) -> (f64, usize) {
    let delta_lambda = std::f64::consts::FRAC_PI_2 / (nu * length / 2.0);
    let panels = ((lambda_max - 1.0) / delta_lambda).ceil() as usize;
    (delta_lambda, panels)
}

fn physical_coefficient(conditions: &Conditions) -> f64 {
    4.0 * conditions.fluid.density * conditions.gravity.powi(2)
        / (std::f64::consts::PI * conditions.speed.powi(2))
}

/// Independent positive-real-axis reference using `lambda = 1 + t^2`.
/// Panel boundaries advance the bow/stern phase by at most pi/2; each panel
/// uses a locally generated Gauss rule and the total uses Kahan accumulation.
fn quadratic_endpoint_reference(
    length: f64,
    beam: f64,
    draft: f64,
    conditions: &Conditions,
    order: usize,
    lambda_max: f64,
) -> (f64, usize) {
    let nu = conditions.gravity / conditions.speed.powi(2);
    let (delta_lambda, panels) = phase_resolved_panel_count(length, nu, lambda_max);
    let rule = gauss_legendre(order);
    let mut integral = 0.0;
    let mut correction = 0.0;
    for panel in 0..panels {
        let lambda0 = 1.0 + panel as f64 * delta_lambda;
        let lambda1 = (lambda0 + delta_lambda).min(lambda_max);
        let t0 = (lambda0 - 1.0).sqrt();
        let t1 = (lambda1 - 1.0).sqrt();
        let half = 0.5 * (t1 - t0);
        let mid = 0.5 * (t0 + t1);
        for &(node, weight) in &rule {
            let t = mid + half * node;
            let lambda = 1.0 + t * t;
            let j = wigley_j(length, beam, draft, nu, lambda);
            let value = weight * half * j * j * 2.0 * lambda.powi(2) / (2.0 + t * t).sqrt();
            add_kahan(&mut integral, &mut correction, value);
        }
    }
    (physical_coefficient(conditions) * integral, panels * order)
}

/// A separately transformed endpoint-aware reference using `lambda=sec(theta)`
/// as recommended by Tuck. Phase-resolved lambda panels map independently to
/// theta panels, and the transformed weight is `sec(theta)^3`.
fn secant_endpoint_reference(
    length: f64,
    beam: f64,
    draft: f64,
    conditions: &Conditions,
    order: usize,
    lambda_max: f64,
) -> (f64, usize) {
    let nu = conditions.gravity / conditions.speed.powi(2);
    let (delta_lambda, panels) = phase_resolved_panel_count(length, nu, lambda_max);
    let rule = gauss_legendre(order);
    let mut integral = 0.0;
    let mut correction = 0.0;
    for panel in 0..panels {
        let lambda0 = 1.0 + panel as f64 * delta_lambda;
        let lambda1 = (lambda0 + delta_lambda).min(lambda_max);
        let theta0 = (1.0 / lambda0).acos();
        let theta1 = (1.0 / lambda1).acos();
        let half = 0.5 * (theta1 - theta0);
        let mid = 0.5 * (theta0 + theta1);
        for &(node, weight) in &rule {
            let theta = mid + half * node;
            let lambda = 1.0 / theta.cos();
            let j = wigley_j(length, beam, draft, nu, lambda);
            add_kahan(
                &mut integral,
                &mut correction,
                weight * half * j * j * lambda.powi(3),
            );
        }
    }
    (physical_coefficient(conditions) * integral, panels * order)
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
            quadratic_endpoint_reference(length, beam, draft, &conditions, 16, 4_000.0);
        let (quadratic_order_8, _) =
            quadratic_endpoint_reference(length, beam, draft, &conditions, 8, 4_000.0);
        let (quadratic_order_24, _) =
            quadratic_endpoint_reference(length, beam, draft, &conditions, 24, 4_000.0);
        let (secant_reference, secant_evaluations) =
            secant_endpoint_reference(length, beam, draft, &conditions, 16, 4_000.0);
        let (secant_order_8, _) =
            secant_endpoint_reference(length, beam, draft, &conditions, 8, 4_000.0);
        let (secant_order_24, _) =
            secant_endpoint_reference(length, beam, draft, &conditions, 24, 4_000.0);
        let endpoint_map_difference = (reference - secant_reference).abs() / reference;
        let quadratic_order_8_16 = (reference - quadratic_order_8).abs() / reference;
        let quadratic_order_16_24 = (reference - quadratic_order_24).abs() / quadratic_order_24;
        let secant_order_8_16 = (secant_reference - secant_order_8).abs() / secant_reference;
        let secant_order_16_24 = (secant_reference - secant_order_24).abs() / secant_order_24;
        eprintln!(
            "Fn={fn_:.2} reference convergence: endpoint_map(t2_16-vs-sec_16)={endpoint_map_difference:.3e}, panel_order_t2(8-vs-16)={quadratic_order_8_16:.3e}, panel_order_t2(16-vs-24)={quadratic_order_16_24:.3e}, panel_order_sec(8-vs-16)={secant_order_8_16:.3e}, panel_order_sec(16-vs-24)={secant_order_16_24:.3e}, lambda_max=4000, evaluations_t2={reference_evaluations}, evaluations_sec={secant_evaluations}"
        );
        assert!(endpoint_map_difference < 5e-12);
        assert!(quadratic_order_16_24 < 5e-13);
        assert!(secant_order_16_24 < 5e-13);
        if fn_ == 0.02 {
            let (doubled_reference, _) =
                quadratic_endpoint_reference(length, beam, draft, &conditions, 16, 8_000.0);
            let (doubled_secant, _) =
                secant_endpoint_reference(length, beam, draft, &conditions, 16, 8_000.0);
            eprintln!(
                "Fn=0.02 cutoff convergence: t2(4000-vs-8000)={:.3e}, sec(4000-vs-8000)={:.3e}",
                (reference - doubled_reference).abs() / doubled_reference,
                (secant_reference - doubled_secant).abs() / doubled_secant,
            );
        }
        let reduced = low_froude_wave_resistance(&hull, &conditions).unwrap();
        let actual_relative_error = (reduced.resistance - reference).abs() / reference;
        let direct_relative_error = (direct.resistance - reference).abs() / reference;
        let abs_error = reduced.omitted_abs_error_bound
            + reduced.quadrature_abs_error_estimate
            + reduced.endpoint_summation_abs_error_bound;
        let denominator = reduced.resistance.abs() - abs_error;
        eprintln!(
            "Fn={fn_:.2}: reference={reference:.12e}, reduced={:.12e}, reduced_error={actual_relative_error:.3e}, reduced_est={:.3e}, omission_rel={:.3e}, contour_rel={:.3e}, rounding_rel={:.3e}, safety_rel={:.3e}, omitted_physical_rel={:.3e}, marching_error={direct_relative_error:.3e}, marching_est={:.3e}, marching_status={}, marching_evals={}, kernel_evals={}, reference_evals={reference_evaluations}, terms={}/{}",
            reduced.resistance,
            reduced.est_rel_error,
            reduced.omitted_abs_error_bound / denominator,
            reduced.quadrature_abs_error_estimate / denominator,
            reduced.endpoint_summation_abs_error_bound / denominator,
            reduced.coefficient_rel_error_bound,
            reduced.omitted_abs_error_bound / reference,
            direct.est_rel_error,
            direct.outcome.as_str(),
            direct.inner_evaluations,
            reduced.kernel_evaluations,
            reduced.waterline_terms,
            reduced.endpoint_terms,
        );
        assert!(actual_relative_error <= 1.05 * reduced.est_rel_error.max(2e-10));
        if fn_ <= 0.05 {
            assert!(reduced.kernel_evaluations * 100 < direct.inner_evaluations);
            assert!(actual_relative_error < 2e-12);
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
    let (reference, _) =
        quadratic_endpoint_reference(length, beam, draft, &conditions, 16, 4_000.0);
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
    let fn_: f64 = 0.05;
    let speed = fn_ * (STANDARD_GRAVITY * hull.length()).sqrt();
    let result = low_froude_wave_resistance(&hull, &Conditions::freshwater(speed)).unwrap();

    assert_eq!(hull.surface().degree_x(), 2);
    assert_eq!(hull.surface().degree_z(), 2);
    let x_spans = hull
        .surface()
        .knots_x()
        .windows(2)
        .filter(|pair| pair[1] > pair[0])
        .count();
    let z_spans = hull
        .surface()
        .knots_z()
        .windows(2)
        .filter(|pair| pair[1] > pair[0])
        .count();
    assert_eq!(x_spans, 1);
    assert_eq!(z_spans, 1);
    assert_eq!(result.endpoint_terms, 16);
    assert_eq!(result.waterline_terms, 8);
    let nonzero_frequency_pairs = result
        .endpoint_pairs
        .iter()
        .filter(|pair| pair.kernel_evaluations > 0)
        .count();
    assert_eq!(nonzero_frequency_pairs, 16);
    assert!(result
        .endpoint_pairs
        .iter()
        .filter(|pair| pair.kernel_evaluations > 0)
        .all(|pair| pair.kernel_evaluations == 72));
    assert_eq!(result.kernel_evaluations, 16 * 72);
    let nu = STANDARD_GRAVITY / speed.powi(2);
    assert!((nu * hull.length() - 1.0 / fn_.powi(2)).abs() < 1e-12);
    eprintln!(
        "Wigley arithmetic: degrees=({},{}), spans=({},{}), endpoint_terms={}, waterline_terms={}, nonzero_frequency_pairs={}, nodes_per_pair=72, total_nodes={}, omega_bow_stern={:.1}",
        hull.surface().degree_x(),
        hull.surface().degree_z(),
        x_spans,
        z_spans,
        result.endpoint_terms,
        result.waterline_terms,
        nonzero_frequency_pairs,
        result.kernel_evaluations,
        nu * hull.length(),
    );

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

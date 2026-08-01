use michell::{
    hulls, wave_resistance_gradient_with, wave_resistance_with, BSplineSurface, Conditions, Hull,
    WaveOptions, STANDARD_GRAVITY,
};

fn rebuild(surface: &BSplineSurface, control: Vec<f64>) -> Hull {
    Hull::new(
        BSplineSurface::new(
            surface.degree_x(),
            surface.degree_z(),
            surface.knots_x().to_vec(),
            surface.knots_z().to_vec(),
            control,
        )
        .unwrap(),
    )
    .unwrap()
}

fn positive_wigley() -> Hull {
    let wigley = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    rebuild(
        wigley.surface(),
        wigley
            .surface()
            .control()
            .iter()
            .map(|value| value + 0.2)
            .collect(),
    )
}

#[test]
fn control_gradient_matches_centered_finite_differences() {
    let hull = positive_wigley();
    let speed = 0.35 * (STANDARD_GRAVITY * hull.length()).sqrt();
    let conditions = Conditions::freshwater(speed);
    let options = WaveOptions {
        rel_tol: 1e-9,
        max_refinements: 6,
    };
    let analytic = wave_resistance_gradient_with(&hull, &conditions, &options).unwrap();
    let direct = wave_resistance_with(&hull, &conditions, &options).unwrap();
    assert_eq!(analytic.wave.resistance, direct.resistance);

    let step = 1e-5;
    for index in 0..hull.surface().control().len() {
        let mut plus = hull.surface().control().to_vec();
        let mut minus = plus.clone();
        plus[index] += step;
        minus[index] -= step;
        let r_plus = wave_resistance_with(&rebuild(hull.surface(), plus), &conditions, &options)
            .unwrap()
            .resistance;
        let r_minus = wave_resistance_with(&rebuild(hull.surface(), minus), &conditions, &options)
            .unwrap()
            .resistance;
        let finite_difference = (r_plus - r_minus) / (2.0 * step);
        let scale = finite_difference.abs().max(1.0);
        let relative = (analytic.control_gradient[index] - finite_difference).abs() / scale;
        assert!(
            relative < 2e-6,
            "control {index}: analytic={}, finite difference={finite_difference}, scaled error={relative:.3e}",
            analytic.control_gradient[index]
        );
    }
}

#[test]
fn constant_half_breadth_shift_is_a_null_direction() {
    let hull = positive_wigley();
    let speed = 0.30 * (STANDARD_GRAVITY * hull.length()).sqrt();
    let gradient = wave_resistance_gradient_with(
        &hull,
        &Conditions::freshwater(speed),
        &WaveOptions::default(),
    )
    .unwrap();
    let directional: f64 = gradient.control_gradient.iter().sum();
    let scale: f64 = gradient
        .control_gradient
        .iter()
        .map(|value| value.abs())
        .sum();
    assert!(directional.abs() < 1e-11 * scale.max(1.0));
}

#[test]
fn reverse_pass_cost_is_independent_of_control_count() {
    let hull = positive_wigley();
    let speed = 0.35 * (STANDARD_GRAVITY * hull.length()).sqrt();
    let gradient = wave_resistance_gradient_with(
        &hull,
        &Conditions::freshwater(speed),
        &WaveOptions::default(),
    )
    .unwrap();
    let finite_difference_evaluations =
        2 * hull.surface().control().len() * gradient.wave.inner_evaluations;
    let analytic_evaluations = gradient.wave.inner_evaluations + gradient.gradient_evaluations;
    assert!(analytic_evaluations < finite_difference_evaluations / 5);
}

#[test]
fn asymmetric_control_net_requires_an_explicit_two_side_gradient() {
    let base = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let clone_surface = || {
        BSplineSurface::new(
            base.surface().degree_x(),
            base.surface().degree_z(),
            base.surface().knots_x().to_vec(),
            base.surface().knots_z().to_vec(),
            base.surface().control().to_vec(),
        )
        .unwrap()
    };
    let asymmetric = Hull::new_asymmetric(clone_surface(), clone_surface()).unwrap();
    let conditions = Conditions::freshwater(3.0);
    assert!(
        wave_resistance_gradient_with(&asymmetric, &conditions, &WaveOptions::default()).is_err()
    );
}

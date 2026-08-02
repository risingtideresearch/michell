use michell::{
    hulls, multihull_wave_resistance_gradient_with, multihull_wave_resistance_with,
    wave_resistance_gradient_with, wave_resistance_with, BSplineSurface, Conditions,
    ControlNetGradient, Hull, Placement, WaveMethod, WaveOptions, STANDARD_GRAVITY,
};

fn surface_with_control(surface: &BSplineSurface, control: Vec<f64>) -> BSplineSurface {
    BSplineSurface::new(
        surface.degree_x(),
        surface.degree_z(),
        surface.knots_x().to_vec(),
        surface.knots_z().to_vec(),
        control,
    )
    .unwrap()
}

fn rebuild(surface: &BSplineSurface, control: Vec<f64>) -> Hull {
    Hull::new(surface_with_control(surface, control)).unwrap()
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
    let scale: f64 = gradient.control_gradient.iter().map(|value| value.abs()).sum();
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
    assert!(wave_resistance_gradient_with(&asymmetric, &conditions, &WaveOptions::default()).is_err());
}

#[test]
fn multihull_placement_gradient_matches_centered_finite_differences() {
    let hull = positive_wigley();
    let placements = [
        Placement { x: 0.31, y: 0.74 },
        Placement { x: -0.17, y: -0.91 },
    ];
    let conditions = Conditions::freshwater(3.0);
    let options = WaveOptions {
        rel_tol: 1e-9,
        max_refinements: 6,
    };
    let analytic = multihull_wave_resistance_gradient_with(
        &[(&hull, placements[0]), (&hull, placements[1])],
        &conditions,
        &options,
    )
    .unwrap();
    let step = 5e-3;

    for member in 0..placements.len() {
        for transverse in [false, true] {
            let mut plus = placements;
            let mut minus = placements;
            if transverse {
                plus[member].y += step;
                minus[member].y -= step;
            } else {
                plus[member].x += step;
                minus[member].x -= step;
            }
            let r_plus = multihull_wave_resistance_with(
                &[(&hull, plus[0]), (&hull, plus[1])],
                &conditions,
                &options,
            )
            .unwrap()
            .resistance;
            let r_minus = multihull_wave_resistance_with(
                &[(&hull, minus[0]), (&hull, minus[1])],
                &conditions,
                &options,
            )
            .unwrap()
            .resistance;
            let finite_difference = (r_plus - r_minus) / (2.0 * step);
            let derivative = if transverse {
                analytic.members[member].placement.transverse
            } else {
                analytic.members[member].placement.longitudinal
            };
            let scaled_error =
                (derivative - finite_difference).abs() / finite_difference.abs().max(1.0);
            assert!(
                scaled_error < 1e-4,
                "member {member} {}: analytic={derivative}, finite difference={finite_difference}, scaled error={scaled_error:.3e}",
                if transverse { "y" } else { "x" },
            );
        }
    }
}

#[test]
fn multihull_member_control_gradients_include_interference() {
    let hulls = [positive_wigley(), positive_wigley()];
    let placements = [
        Placement { x: 0.31, y: 0.74 },
        Placement { x: -0.17, y: -0.91 },
    ];
    let conditions = Conditions::freshwater(3.0);
    let options = WaveOptions {
        rel_tol: 1e-9,
        max_refinements: 6,
    };
    let analytic = multihull_wave_resistance_gradient_with(
        &[(&hulls[0], placements[0]), (&hulls[1], placements[1])],
        &conditions,
        &options,
    )
    .unwrap();
    let step = 1e-4;
    let count = hulls[0].surface().control().len();

    for member in 0..hulls.len() {
        let ControlNetGradient::Symmetric(expected) = &analytic.members[member].control else {
            panic!("symmetric member returned an asymmetric control gradient");
        };
        for index in [0, count / 2, count - 1] {
            let mut plus_control = hulls[member].surface().control().to_vec();
            let mut minus_control = plus_control.clone();
            plus_control[index] += step;
            minus_control[index] -= step;
            let plus = rebuild(hulls[member].surface(), plus_control);
            let minus = rebuild(hulls[member].surface(), minus_control);
            let other = 1 - member;
            let plus_members = if member == 0 {
                [(&plus, placements[0]), (&hulls[other], placements[1])]
            } else {
                [(&hulls[other], placements[0]), (&plus, placements[1])]
            };
            let minus_members = if member == 0 {
                [(&minus, placements[0]), (&hulls[other], placements[1])]
            } else {
                [(&hulls[other], placements[0]), (&minus, placements[1])]
            };
            let r_plus =
                multihull_wave_resistance_with(&plus_members, &conditions, &options)
                    .unwrap()
                    .resistance;
            let r_minus =
                multihull_wave_resistance_with(&minus_members, &conditions, &options)
                    .unwrap()
                    .resistance;
            let finite_difference = (r_plus - r_minus) / (2.0 * step);
            let scaled_error =
                (expected[index] - finite_difference).abs() / finite_difference.abs().max(1.0);
            assert!(
                scaled_error < 1e-4,
                "member {member} control {index}: analytic={}, finite difference={finite_difference}, scaled error={scaled_error:.3e}",
                expected[index],
            );
        }
    }
}

#[test]
fn asymmetric_side_gradients_match_centered_finite_differences() {
    let base = positive_wigley();
    let port_control: Vec<f64> = base
        .surface()
        .control()
        .iter()
        .map(|value| 0.9 * value)
        .collect();
    let starboard_control: Vec<f64> = base
        .surface()
        .control()
        .iter()
        .map(|value| 1.1 * value)
        .collect();
    let make_hull = |port: Vec<f64>, starboard: Vec<f64>| {
        Hull::new_asymmetric(
            surface_with_control(base.surface(), port),
            surface_with_control(base.surface(), starboard),
        )
        .unwrap()
    };
    let hull = make_hull(port_control.clone(), starboard_control.clone());
    let conditions = Conditions::freshwater(3.0);
    let options = WaveOptions {
        rel_tol: 1e-9,
        max_refinements: 6,
    };
    let analytic = multihull_wave_resistance_gradient_with(
        &[(&hull, Placement::default())],
        &conditions,
        &options,
    )
    .unwrap();
    let ControlNetGradient::Asymmetric { port, starboard } =
        &analytic.members[0].control
    else {
        panic!("asymmetric hull returned a symmetric control gradient");
    };
    let step = 1e-5;

    for side_is_starboard in [false, true] {
        let expected = if side_is_starboard { starboard } else { port };
        for index in 0..expected.len() {
            let mut port_plus = port_control.clone();
            let mut port_minus = port_control.clone();
            let mut starboard_plus = starboard_control.clone();
            let mut starboard_minus = starboard_control.clone();
            if side_is_starboard {
                starboard_plus[index] += step;
                starboard_minus[index] -= step;
            } else {
                port_plus[index] += step;
                port_minus[index] -= step;
            }
            let plus = make_hull(port_plus, starboard_plus);
            let minus = make_hull(port_minus, starboard_minus);
            let r_plus = multihull_wave_resistance_with(
                &[(&plus, Placement::default())],
                &conditions,
                &options,
            )
            .unwrap()
            .resistance;
            let r_minus = multihull_wave_resistance_with(
                &[(&minus, Placement::default())],
                &conditions,
                &options,
            )
            .unwrap()
            .resistance;
            let finite_difference = (r_plus - r_minus) / (2.0 * step);
            let scaled_error = (expected[index] - finite_difference).abs()
                / finite_difference.abs().max(1.0);
            assert!(
                scaled_error < 3e-6,
                "{} control {index}: analytic={}, finite difference={finite_difference}, scaled error={scaled_error:.3e}",
                if side_is_starboard { "starboard" } else { "port" },
                expected[index],
            );
        }
    }
}

#[test]
fn gradient_api_uses_general_marcher_when_primal_uses_endpoint_reduction() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let speed = 0.05 * (STANDARD_GRAVITY * hull.length()).sqrt();
    let conditions = Conditions::freshwater(speed);
    let options = WaveOptions::default();
    let primal = wave_resistance_with(&hull, &conditions, &options).unwrap();
    let gradient = wave_resistance_gradient_with(&hull, &conditions, &options).unwrap();

    assert_eq!(primal.method, WaveMethod::EndpointReduction);
    assert_eq!(gradient.wave.method, WaveMethod::GeneralMarcher);
}

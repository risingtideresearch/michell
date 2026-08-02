use michell::{hulls, BSplineSurface, ConstraintGradient, Hull, HullConstraintGradients};

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

fn assert_constraint_component(
    label: &str,
    index: usize,
    analytic: &ConstraintGradient,
    plus: &Hull,
    minus: &Hull,
    step: f64,
) {
    let checks = [
        (
            "volume",
            analytic.displaced_volume[index],
            (plus.displaced_volume() - minus.displaced_volume()) / (2.0 * step),
            2e-8,
        ),
        (
            "LCB",
            analytic.lcb_x[index],
            (plus.lcb_x() - minus.lcb_x()) / (2.0 * step),
            2e-7,
        ),
        (
            "wetted surface",
            analytic.wetted_surface[index],
            (plus.wetted_surface() - minus.wetted_surface()) / (2.0 * step),
            2e-7,
        ),
    ];
    for (metric, derivative, finite_difference, tolerance) in checks {
        let scaled_error =
            (derivative - finite_difference).abs() / finite_difference.abs().max(1.0);
        assert!(
            scaled_error < tolerance,
            "{label} control {index} {metric}: analytic={derivative}, finite difference={finite_difference}, scaled error={scaled_error:.3e}",
        );
    }
}

#[test]
fn symmetric_constraint_gradients_match_centered_finite_differences() {
    let hull = positive_wigley();
    let HullConstraintGradients::Symmetric(analytic) = hull.constraint_gradients().unwrap() else {
        panic!("symmetric hull returned asymmetric constraint gradients");
    };
    let step = 1e-6;

    for index in 0..hull.surface().control().len() {
        let mut plus_control = hull.surface().control().to_vec();
        let mut minus_control = plus_control.clone();
        plus_control[index] += step;
        minus_control[index] -= step;
        let plus = rebuild(hull.surface(), plus_control);
        let minus = rebuild(hull.surface(), minus_control);
        assert_constraint_component("symmetric", index, &analytic, &plus, &minus, step);
    }
}

#[test]
fn asymmetric_constraint_gradients_match_centered_finite_differences() {
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
    let HullConstraintGradients::Asymmetric { port, starboard } =
        hull.constraint_gradients().unwrap()
    else {
        panic!("asymmetric hull returned symmetric constraint gradients");
    };
    let step = 1e-6;

    for starboard_side in [false, true] {
        let analytic = if starboard_side { &starboard } else { &port };
        for index in 0..analytic.displaced_volume.len() {
            let mut port_plus = port_control.clone();
            let mut port_minus = port_control.clone();
            let mut starboard_plus = starboard_control.clone();
            let mut starboard_minus = starboard_control.clone();
            if starboard_side {
                starboard_plus[index] += step;
                starboard_minus[index] -= step;
            } else {
                port_plus[index] += step;
                port_minus[index] -= step;
            }
            let plus = make_hull(port_plus, starboard_plus);
            let minus = make_hull(port_minus, starboard_minus);
            assert_constraint_component(
                if starboard_side { "starboard" } else { "port" },
                index,
                analytic,
                &plus,
                &minus,
                step,
            );
        }
    }
}

#[test]
fn zero_volume_reports_undefined_lcb_gradient() {
    let base = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let hull = rebuild(base.surface(), vec![0.0; base.surface().control().len()]);
    assert!(hull.constraint_gradients().is_err());
}

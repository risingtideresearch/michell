use michell::{
    low_froude_wave_resistance, wave_resistance_with, BSplineSurface, Conditions, Hull,
    WaveOptions, STANDARD_GRAVITY,
};

fn clamped(degree: usize, a: f64, interior: &[f64], b: f64) -> Vec<f64> {
    std::iter::repeat_n(a, degree + 1)
        .chain(interior.iter().copied())
        .chain(std::iter::repeat_n(b, degree + 1))
        .collect()
}

fn hull_with_knots(interior: &[f64], draft: f64) -> Hull {
    let (p, q) = (2, 2);
    let x_knots = clamped(p, -5.0, interior, 5.0);
    let z_knots = clamped(q, 0.0, &[], draft);
    let nx = x_knots.len() - p - 1;
    let mut control = Vec::new();
    for i in 0..nx {
        let x = (x_knots[i + 1] + x_knots[i + 2]) * 0.5;
        let bx = 0.5 * (1.0 - (x / 5.0).powi(2)).max(0.0);
        for j in 0..=q {
            let v = j as f64 / q as f64;
            control.push(bx * (1.0 - v * v));
        }
    }
    Hull::new(BSplineSurface::new(p, q, x_knots, z_knots, control).unwrap()).unwrap()
}

fn report(label: &str, hull: &Hull, fn_: f64) {
    let speed = fn_ * (STANDARD_GRAVITY * hull.length()).sqrt();
    let cond = Conditions::freshwater(speed);
    let low = low_froude_wave_resistance(hull, &cond);
    let default = wave_resistance_with(hull, &cond, &WaveOptions::default()).unwrap();
    match low {
        Ok(low) => println!(
            "{label} Fn={fn_:.4} low=OK low_est={:.3e} default={:?} default_est={:.3e} terms={}/{} nodes={}",
            low.est_rel_error,
            default.method,
            default.est_rel_error,
            low.waterline_terms,
            low.endpoint_terms,
            low.kernel_evaluations,
        ),
        Err(error) => println!(
            "{label} Fn={fn_:.4} low=REJECT error={error} default={:?} default_est={:.3e}",
            default.method, default.est_rel_error,
        ),
    }
}

fn main() {
    let deep = hull_with_knots(&[], 5.0);
    for fn_ in [0.201, 0.200_000_000_1, 0.2, 0.199_999_999_9, 0.199, 0.15, 0.10] {
        report("deep-single-span", &deep, fn_);
    }

    let smooth_close = hull_with_knots(&[-1.0, -0.9, 1.0], 0.625);
    for fn_ in [0.08, 0.05, 0.03, 0.02] {
        report("smooth-close-knots", &smooth_close, fn_);
    }

    let chine_close = hull_with_knots(&[-1.0, -1.0, -0.9, -0.9, 1.0, 1.0], 0.625);
    for fn_ in [0.05, 0.03, 0.02] {
        report("full-multiplicity-chines", &chine_close, fn_);
    }
}

use michell::{wave_resistance_with, BSplineSurface, Conditions, Hull, WaveOptions, STANDARD_GRAVITY};

fn hull(p: usize) -> Hull {
    let q = 2;
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let mut control = Vec::new();
    for i in 0..=p {
        let u = i as f64 / p as f64;
        let bx = 2.0 * beam * u * (1.0 - u);
        for j in 0..=q {
            let v = j as f64 / q as f64;
            control.push(bx * (1.0 - v * v));
        }
    }
    let knots_x = std::iter::repeat_n(-0.5 * length, p + 1)
        .chain(std::iter::repeat_n(0.5 * length, p + 1))
        .collect();
    let knots_z = std::iter::repeat_n(0.0, q + 1)
        .chain(std::iter::repeat_n(draft, q + 1))
        .collect();
    Hull::new(BSplineSurface::new(p, q, knots_x, knots_z, control).unwrap()).unwrap()
}

fn main() {
    const P2_REFERENCE: f64 = 2.565_891_321_290e-3;
    for p in 20..=48 {
        let hull = hull(p);
        let speed = 0.05 * (STANDARD_GRAVITY * hull.length()).sqrt();
        let result = wave_resistance_with(
            &hull,
            &Conditions::freshwater(speed),
            &WaveOptions::default(),
        )
        .unwrap();
        let scale = 2.0 * (1.0 - 1.0 / p as f64);
        let exact = P2_REFERENCE * scale * scale;
        let actual_rel = (result.resistance - exact).abs() / exact;
        println!(
            "p={p} production={:.12e} exact={exact:.12e} reported={:.3e} actual_rel={actual_rel:.3e} method={:?}",
            result.resistance, result.est_rel_error, result.method
        );
    }
}

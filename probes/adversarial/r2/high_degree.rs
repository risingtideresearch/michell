use michell::{low_froude_wave_resistance, multihull_wave_resistance_with, wave_resistance_with, BSplineSurface, Conditions, Hull, Placement, WaveOptions, STANDARD_GRAVITY};

fn make_hull(p: usize) -> Hull {
    let q = 2usize;
    let mut kx = vec![-5.0; p + 1];
    kx.extend(std::iter::repeat_n(5.0, p + 1));
    let mut kz = vec![0.0; q + 1];
    kz.extend(std::iter::repeat_n(5.0, q + 1));
    let mut control = vec![0.0; (p + 1) * (q + 1)];
    for i in 1..p {
        let x = i as f64 / p as f64;
        for j in 0..q {
            let z = j as f64 / q as f64;
            control[i * (q + 1) + j] = 4.0 * x * (1.0 - x) * (1.0 - z * z);
        }
    }
    Hull::new(BSplineSurface::new(p, q, kx, kz, control).unwrap()).unwrap()
}

fn main() {
    for p in [
        30usize, 32, 34, 36, 38, 40, 42, 44, 46, 48, 64, 96, 128, 160, 192,
    ] {
        let hull = make_hull(p);
        let speed = 0.1 * (STANDARD_GRAVITY * hull.length()).sqrt();
        let cond = Conditions::freshwater(speed);
        let low = low_froude_wave_resistance(&hull, &cond);
        let hybrid = wave_resistance_with(&hull, &cond, &WaveOptions::default());
        let general = if p == 48 {
            Some(multihull_wave_resistance_with(
                &[(&hull, Placement::default())],
                &cond,
                &WaveOptions { rel_tol: 1e-8, max_refinements: 10 },
            ))
        } else {
            None
        };
        let low_summary = match low {
            Ok(x) => format!(
                "ok finite={} R={:.9e} est={:.3e} terms={}",
                x.resistance.is_finite() && x.est_rel_error.is_finite(),
                x.resistance,
                x.est_rel_error,
                x.endpoint_terms
            ),
            Err(e) => format!("err={e}"),
        };
        let hybrid_summary = match hybrid {
            Ok(x) => format!(
                "ok finite={} R={:.9e} est={:.3e} method={:?} outcome={:?}",
                x.resistance.is_finite() && x.est_rel_error.is_finite(),
                x.resistance,
                x.est_rel_error,
                x.method,
                x.outcome
            ),
            Err(e) => format!("err={e}"),
        };
        println!("p={p} low=[{low_summary}] hybrid=[{hybrid_summary}]");
        if let Some(general) = general {
            match general {
                Ok(x) => println!(
                    "p={p} forced-general R={:.15e} est={:.3e} outcome={:?} evals={}",
                    x.resistance, x.est_rel_error, x.outcome, x.inner_evaluations
                ),
                Err(e) => println!("p={p} forced-general err={e}"),
            }
        }
    }
}

use michell::{hulls, resistance, wave_resistance, Conditions, STANDARD_GRAVITY};
use std::hint::black_box;
use std::time::{Duration, Instant};

fn sample_count() -> usize {
    std::env::var("MICHELL_BENCH_SAMPLES")
        .ok()
        .and_then(|value| value.parse().ok())
        .filter(|&value| value > 0)
        .unwrap_or(30)
}

fn measure(samples: usize, mut run: impl FnMut() -> f64) -> (Duration, Duration, f64) {
    let mut times = Vec::with_capacity(samples);
    let mut checksum = 0.0;
    for _ in 0..samples {
        let start = Instant::now();
        checksum += black_box(run());
        times.push(start.elapsed());
    }
    times.sort_unstable();
    let median = times[times.len() / 2];
    let best = times[0];
    (median, best, checksum)
}

fn millis(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1e3
}

fn main() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).expect("valid Wigley hull");
    let samples = sample_count();

    // Warm the instruction/data caches before collecting samples.
    let warm_speed = 0.30 * (STANDARD_GRAVITY * hull.length()).sqrt();
    black_box(resistance(&hull, &Conditions::freshwater(warm_speed)).unwrap());

    let (sweep_median, sweep_best, sweep_checksum) = measure(samples, || {
        let mut total = 0.0;
        for index in 0..=20 {
            let fn_ = 0.10 + 0.02 * index as f64;
            let speed = fn_ * (STANDARD_GRAVITY * hull.length()).sqrt();
            total += resistance(&hull, &Conditions::freshwater(speed))
                .unwrap()
                .total;
        }
        total
    });

    let low_fn = 0.05;
    let low_speed = low_fn * (STANDARD_GRAVITY * hull.length()).sqrt();
    let low_cond = Conditions::freshwater(low_speed);
    let low_diagnostics = wave_resistance(&hull, &low_cond).unwrap();
    let (low_median, low_best, low_checksum) = measure(samples, || {
        wave_resistance(&hull, &low_cond).unwrap().resistance
    });

    println!("michell Wigley benchmark ({samples} samples, default rel_tol=1e-5)");
    println!("case                         median_ms      best_ms");
    println!(
        "21 speeds Fn=0.10..0.50     {:>10.3}   {:>10.3}",
        millis(sweep_median),
        millis(sweep_best)
    );
    println!(
        "low Froude Fn={low_fn:.2}           {:>10.3}   {:>10.3}",
        millis(low_median),
        millis(low_best)
    );
    println!(
        "low-Fn diagnostics: evaluations={}, max_lambda={:.3}, est_rel_error={:.3e}",
        low_diagnostics.inner_evaluations,
        low_diagnostics.max_lambda,
        low_diagnostics.est_rel_error
    );
    println!(
        "checksums: sweep={:.12e}, low={:.12e}",
        sweep_checksum, low_checksum
    );
}

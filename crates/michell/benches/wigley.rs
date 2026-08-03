use michell::{
    hulls, low_froude_wave_resistance, multihull_wave_resistance_with, resistance, wave_resistance,
    wave_resistance_gradient, BSplineSurface, Conditions, Hull, Placement, WaveOptions,
    WaveResistance, STANDARD_GRAVITY,
};
use std::hint::black_box;
use std::time::{Duration, Instant};

fn sample_count() -> usize {
    std::env::var("MICHELL_BENCH_SAMPLES")
        .ok()
        .and_then(|value| value.parse().ok())
        .filter(|&value| value > 0)
        .unwrap_or(30)
}

struct SampleStats {
    median: Duration,
    q1: Duration,
    q3: Duration,
    best: Duration,
    checksum: f64,
}

fn sample_stats(mut times: Vec<Duration>, checksum: f64) -> SampleStats {
    times.sort_unstable();
    SampleStats {
        median: times[times.len() / 2],
        q1: times[times.len() / 4],
        q3: times[3 * times.len() / 4],
        best: times[0],
        checksum,
    }
}

fn measure(samples: usize, mut run: impl FnMut() -> f64) -> SampleStats {
    black_box(run());
    let mut times = Vec::with_capacity(samples);
    let mut checksum = 0.0;
    for _ in 0..samples {
        let start = Instant::now();
        checksum += black_box(run());
        times.push(start.elapsed());
    }
    sample_stats(times, checksum)
}

fn measure_pair(
    samples: usize,
    mut left: impl FnMut() -> f64,
    mut right: impl FnMut() -> f64,
) -> (SampleStats, SampleStats) {
    black_box(left());
    black_box(right());
    let mut left_times = Vec::with_capacity(samples);
    let mut right_times = Vec::with_capacity(samples);
    let mut left_checksum = 0.0;
    let mut right_checksum = 0.0;
    for sample in 0..samples {
        let mut run_left = || {
            let start = Instant::now();
            left_checksum += black_box(left());
            left_times.push(start.elapsed());
        };
        let mut run_right = || {
            let start = Instant::now();
            right_checksum += black_box(right());
            right_times.push(start.elapsed());
        };
        if sample.is_multiple_of(2) {
            run_left();
            run_right();
        } else {
            run_right();
            run_left();
        }
    }
    (
        sample_stats(left_times, left_checksum),
        sample_stats(right_times, right_checksum),
    )
}

fn millis(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1e3
}

fn print_stats(label: &str, stats: &SampleStats) {
    println!(
        "{label:<28} {:>9.3}  {:>9.3}..{:<9.3} {:>9.3}",
        millis(stats.median),
        millis(stats.q1),
        millis(stats.q3),
        millis(stats.best),
    );
}

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

fn marching_wave_resistance(hull: &Hull, conditions: &Conditions) -> WaveResistance {
    multihull_wave_resistance_with(
        &[(hull, Placement::default())],
        conditions,
        &WaveOptions::default(),
    )
    .unwrap()
}

fn main() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).expect("valid Wigley hull");
    let samples = sample_count();

    let sweep = measure(samples, || {
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
    let low = measure(samples, || {
        wave_resistance(&hull, &low_cond).unwrap().resistance
    });

    let very_low_fn = 0.02;
    let very_low_speed = very_low_fn * (STANDARD_GRAVITY * hull.length()).sqrt();
    let very_low_cond = Conditions::freshwater(very_low_speed);
    let very_low_direct = marching_wave_resistance(&hull, &very_low_cond);
    let very_low_reduced = low_froude_wave_resistance(&hull, &very_low_cond).unwrap();
    let (very_low_direct_samples, very_low_reduced_samples) = measure_pair(
        samples,
        || marching_wave_resistance(&hull, &very_low_cond).resistance,
        || {
            low_froude_wave_resistance(&hull, &very_low_cond)
                .unwrap()
                .resistance
        },
    );

    let gradient_hull = rebuild(
        hull.surface(),
        hull.surface()
            .control()
            .iter()
            .map(|value| value + 0.2)
            .collect(),
    );
    let gradient_speed = 0.35 * (STANDARD_GRAVITY * gradient_hull.length()).sqrt();
    let gradient_cond = Conditions::freshwater(gradient_speed);
    let (gradient_samples, finite_samples) = measure_pair(
        samples,
        || {
            wave_resistance_gradient(&gradient_hull, &gradient_cond)
                .unwrap()
                .control_gradient
                .iter()
                .map(|value| value.abs())
                .sum()
        },
        || {
            let step = 1e-5;
            let mut norm = 0.0;
            for index in 0..gradient_hull.surface().control().len() {
                let mut plus = gradient_hull.surface().control().to_vec();
                let mut minus = plus.clone();
                plus[index] += step;
                minus[index] -= step;
                let r_plus =
                    wave_resistance(&rebuild(gradient_hull.surface(), plus), &gradient_cond)
                        .unwrap()
                        .resistance;
                let r_minus =
                    wave_resistance(&rebuild(gradient_hull.surface(), minus), &gradient_cond)
                        .unwrap()
                        .resistance;
                norm += ((r_plus - r_minus) / (2.0 * step)).abs();
            }
            norm
        },
    );

    println!("michell Wigley benchmark ({samples} samples, default rel_tol=1e-5)");
    println!("case                         median_ms       IQR_ms          best_ms");
    print_stats("21 speeds Fn=0.10..0.50", &sweep);
    print_stats(&format!("low Froude Fn={low_fn:.2}"), &low);
    println!(
        "low-Fn diagnostics: evaluations={}, max_lambda={:.3}, est_rel_error={:.3e}",
        low_diagnostics.inner_evaluations,
        low_diagnostics.max_lambda,
        low_diagnostics.est_rel_error
    );
    println!(
        "checksums: sweep={:.12e}, low={:.12e}",
        sweep.checksum, low.checksum
    );
    println!("very-low-Froude method       median_ms       IQR_ms          best_ms");
    print_stats(
        &format!("marching Fn={very_low_fn:.2}"),
        &very_low_direct_samples,
    );
    print_stats(
        &format!("endpoint/NSD Fn={very_low_fn:.2}"),
        &very_low_reduced_samples,
    );
    println!(
        "very-low-Fn speedup: {:.2}x; evaluations {} -> {}; estimated endpoint error {:.3e}; checksums direct={:.12e}, endpoint={:.12e}",
        very_low_direct_samples.median.as_secs_f64()
            / very_low_reduced_samples.median.as_secs_f64(),
        very_low_direct.inner_evaluations,
        very_low_reduced.kernel_evaluations,
        very_low_reduced.est_rel_error,
        very_low_direct_samples.checksum,
        very_low_reduced_samples.checksum,
    );
    println!("gradient method              median_ms       IQR_ms          best_ms");
    print_stats("exact reverse (9 controls)", &gradient_samples);
    print_stats("centered finite differences", &finite_samples);
    println!(
        "gradient speedup: {:.2}x; checksums exact={:.12e}, finite={:.12e}",
        finite_samples.median.as_secs_f64() / gradient_samples.median.as_secs_f64(),
        gradient_samples.checksum,
        finite_samples.checksum
    );
}

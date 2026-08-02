//! Independent coverage check for the general marcher's error diagnostic.

use michell::{
    hulls, multihull_wave_resistance_with, Conditions, Placement, WaveMethod, WaveOptions,
    WaveOutcome, STANDARD_GRAVITY,
};

#[derive(Default)]
struct KahanSum {
    sum: f64,
    correction: f64,
}

impl KahanSum {
    fn add(&mut self, value: f64) {
        let corrected = value - self.correction;
        let next = self.sum + corrected;
        self.correction = (next - self.sum) - corrected;
        self.sum = next;
    }
}

fn wigley_amplitude(length: f64, beam: f64, draft: f64, nu: f64, lambda: f64) -> f64 {
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

/// Positive-real-axis reference using the analytic Wigley inner amplitude.
/// Panel boundaries advance the bow/stern phase by at most pi/2, 16-point
/// Gauss-Legendre integrates each panel, and Kahan summation preserves the
/// small high-lambda contributions through lambda=4000.
fn deep_wigley_reference(length: f64, beam: f64, draft: f64, conditions: &Conditions) -> f64 {
    const LAMBDA_MAX: f64 = 4_000.0;
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

    let nu = conditions.gravity / conditions.speed.powi(2);
    let half_length = length / 2.0;
    let du = std::f64::consts::FRAC_PI_2 / (nu * half_length);
    let panels = ((LAMBDA_MAX - 1.0) / du).ceil() as usize;
    let mut integral = KahanSum::default();

    for panel in 0..panels {
        let u0 = panel as f64 * du;
        let u1 = ((panel + 1) as f64 * du).min(LAMBDA_MAX - 1.0);
        let t0 = u0.sqrt();
        let t1 = u1.sqrt();
        let half = 0.5 * (t1 - t0);
        let mid = 0.5 * (t0 + t1);
        for (&node, &weight) in NODES.iter().zip(&WEIGHTS) {
            for sign in [-1.0, 1.0] {
                let t = mid + sign * half * node;
                let lambda = 1.0 + t * t;
                let amplitude = wigley_amplitude(length, beam, draft, nu, lambda);
                integral.add(
                    weight * half * amplitude.powi(2) * 2.0 * lambda.powi(2) / (2.0 + t * t).sqrt(),
                );
            }
        }
    }

    let coefficient = 4.0 * conditions.fluid.density * conditions.gravity.powi(2)
        / (std::f64::consts::PI * conditions.speed.powi(2));
    coefficient * integral.sum
}

#[test]
fn design_froude_tail_estimate_covers_actual_error() {
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(length, beam, draft).unwrap();
    for froude in [0.12, 0.20, 0.35] {
        let speed = froude * (STANDARD_GRAVITY * length).sqrt();
        let conditions = Conditions::freshwater(speed);
        let result = multihull_wave_resistance_with(
            &[(&hull, Placement::default())],
            &conditions,
            &WaveOptions {
                rel_tol: 1e-8,
                max_refinements: 8,
            },
        )
        .unwrap();
        let reference = deep_wigley_reference(length, beam, draft, &conditions);
        let actual_relative_error = (result.resistance - reference).abs() / reference;

        assert_eq!(result.method, WaveMethod::GeneralMarcher);
        assert_eq!(result.outcome, WaveOutcome::RefinementCap);
        eprintln!(
            "Fn={froude:.2}: reported={:.6e}, actual={actual_relative_error:.6e}, coverage={:.3}x",
            result.est_rel_error,
            result.est_rel_error / actual_relative_error,
        );
        assert!(
            actual_relative_error <= result.est_rel_error,
            "Fn={froude:.2}: reported {:.6e}, actual {actual_relative_error:.6e}; Rw={:.12e}, reference={reference:.12e}, lambda_max={:.3}, evaluations={}",
            result.est_rel_error,
            result.resistance,
            result.max_lambda,
            result.inner_evaluations,
        );
    }
}

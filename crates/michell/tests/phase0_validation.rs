//! Phase-0 validation harness for the numerical work in this repository.
//!
//! These tests deliberately sit above the unit tests: they pin one published
//! Wigley value, similarity laws, the asymmetric and multihull reductions, and
//! the relationship between requested tolerance, reported error, and a tighter
//! independently requested calculation.

use michell::{
    hulls, wave_resistance_with, BSplineSurface, Conditions, Hull, Placement, WaveOptions,
    STANDARD_GRAVITY,
};

fn rel_err(got: f64, want: f64) -> f64 {
    (got - want).abs() / want.abs().max(f64::MIN_POSITIVE)
}

fn scaled_surface(surface: &BSplineSurface, factor: f64) -> BSplineSurface {
    BSplineSurface::new(
        surface.degree_x(),
        surface.degree_z(),
        surface.knots_x().to_vec(),
        surface.knots_z().to_vec(),
        surface
            .control()
            .iter()
            .map(|value| factor * value)
            .collect(),
    )
    .unwrap()
}

/// Doctors & Beck, "Numerical Aspects of the Neumann-Kelvin Problem",
/// Journal of Ship Research 31(1), 1987, Table 1, reports the classical
/// thin-ship result 10^3 Cw = 1.2486 for the standard Wigley hull at Fn=0.35.
/// Hull proportions in the paper are B/L=0.1 and T/L=0.0625.
/// DOI: https://doi.org/10.5957/jsr.1987.31.1.1
#[test]
fn wigley_matches_published_thin_ship_value() {
    let length = 10.0;
    let hull = hulls::wigley(length, 1.0, 0.625).unwrap();
    let speed = 0.35 * (STANDARD_GRAVITY * length).sqrt();
    let cond = Conditions::freshwater(speed);
    let wave = wave_resistance_with(
        &hull,
        &cond,
        &WaveOptions {
            rel_tol: 1e-8,
            max_refinements: 7,
        },
    )
    .unwrap();
    let cw = wave.resistance / (0.5 * cond.fluid.density * speed * speed * hull.wetted_surface());
    let published_cw = 1.2486e-3;

    // The source prints five significant digits and uses its own physical
    // constants/geometric integration convention, so do not overfit it.
    assert!(
        rel_err(cw, published_cw) < 5e-3,
        "Cw={cw:.10e}, published={published_cw:.10e}"
    );
}

/// Geometrically similar hulls at equal length Froude number have invariant Cw
/// and wave resistance proportional to length cubed.
#[test]
fn froude_similarity_scales_wave_resistance_cubically() {
    let base = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let fn_ = 0.32;
    let base_speed = fn_ * (STANDARD_GRAVITY * base.length()).sqrt();
    let base_cond = Conditions::freshwater(base_speed);
    let base_wave = wave_resistance_with(&base, &base_cond, &WaveOptions::default())
        .unwrap()
        .resistance;

    for scale in [0.25_f64, 0.5, 2.0, 4.0] {
        let hull = hulls::wigley(10.0 * scale, scale, 0.625 * scale).unwrap();
        let speed = fn_ * (STANDARD_GRAVITY * hull.length()).sqrt();
        let cond = Conditions::freshwater(speed);
        let wave = wave_resistance_with(&hull, &cond, &WaveOptions::default())
            .unwrap()
            .resistance;
        assert!(
            rel_err(wave, base_wave * scale.powi(3)) < 2e-8,
            "scale={scale}: Rw={wave}, expected {}",
            base_wave * scale.powi(3)
        );
    }
}

/// Identical port and starboard control nets make the camber field exactly
/// zero, so the asymmetric entry path must reduce to classical Michell.
#[test]
fn symmetric_hull_has_zero_camber_contribution() {
    let symmetric = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let asymmetric_path = Hull::new_asymmetric(
        scaled_surface(symmetric.surface(), 1.0),
        scaled_surface(symmetric.surface(), 1.0),
    )
    .unwrap();

    for fn_ in [0.15, 0.25, 0.35, 0.50] {
        let speed = fn_ * (STANDARD_GRAVITY * symmetric.length()).sqrt();
        let cond = Conditions::freshwater(speed);
        let classical = wave_resistance_with(&symmetric, &cond, &WaveOptions::default())
            .unwrap()
            .resistance;
        let reduced = wave_resistance_with(&asymmetric_path, &cond, &WaveOptions::default())
            .unwrap()
            .resistance;
        assert!(
            rel_err(reduced, classical) < 1e-12,
            "Fn={fn_}: asymmetric path {reduced}, classical {classical}"
        );
    }
}

/// Analytic Wigley inner amplitude used only by the independent outer-integral
/// reference below.
fn wigley_j(length: f64, beam: f64, draft: f64, nu: f64, lambda: f64) -> f64 {
    let half_length = length / 2.0;
    let k = nu * lambda;
    let kappa = nu * lambda * lambda;
    let x = 2.0 * ((k * half_length).sin() / (k * k) - half_length * (k * half_length).cos() / k);
    let decay = (-kappa * draft).exp();
    let z = (1.0 - decay) / kappa
        - (2.0 - decay * (kappa * kappa * draft * draft + 2.0 * kappa * draft + 2.0))
            / (kappa.powi(3) * draft * draft);
    -(4.0 * beam / (length * length)) * x * z
}

fn catamaran_reference(
    length: f64,
    beam: f64,
    draft: f64,
    separation: f64,
    cond: &Conditions,
) -> f64 {
    // Independent, dense Simpson rule over theta, carrying the classical
    // identical-demihull multiplier 4 cos^2(ky*s/2).
    const PANELS: usize = 1_000_000;
    let nu = cond.gravity / (cond.speed * cond.speed);
    let theta_max = (1.0_f64 / 200.0).acos();
    let step = theta_max / PANELS as f64;
    let integrand = |theta: f64| {
        let sec = 1.0 / theta.cos();
        let j = wigley_j(length, beam, draft, nu, sec);
        let half_phase = 0.5 * nu * separation * sec * theta.tan();
        j * j * 4.0 * half_phase.cos().powi(2) * sec.powi(3)
    };
    let mut sum = integrand(0.0) + integrand(theta_max);
    for index in 1..PANELS {
        sum += (if index % 2 == 0 { 2.0 } else { 4.0 }) * integrand(index as f64 * step);
    }
    let integral = sum * step / 3.0;
    4.0 * cond.fluid.density * cond.gravity.powi(2) / (std::f64::consts::PI * cond.speed.powi(2))
        * integral
}

#[test]
fn identical_catamaran_reduces_to_classical_four_cosine_squared() {
    let (length, beam, draft) = (10.0, 1.0, 0.625);
    let separation = 2.0;
    let hull = hulls::wigley(length, beam, draft).unwrap();
    let speed = 0.35 * (STANDARD_GRAVITY * length).sqrt();
    let cond = Conditions::freshwater(speed);
    let members = [
        (
            &hull,
            Placement {
                x: 0.0,
                y: -separation / 2.0,
            },
        ),
        (
            &hull,
            Placement {
                x: 0.0,
                y: separation / 2.0,
            },
        ),
    ];
    let got = michell::multihull_wave_resistance_with(
        &members,
        &cond,
        &WaveOptions {
            rel_tol: 1e-7,
            max_refinements: 7,
        },
    )
    .unwrap()
    .resistance;
    let reference = catamaran_reference(length, beam, draft, separation, &cond);
    assert!(
        rel_err(got, reference) < 2e-4,
        "catamaran Rw={got}, 4cos^2 reference={reference}"
    );
}

#[test]
fn tolerance_tightening_self_converges() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let speed = 0.28 * (STANDARD_GRAVITY * hull.length()).sqrt();
    let cond = Conditions::freshwater(speed);
    let reference = wave_resistance_with(
        &hull,
        &cond,
        &WaveOptions {
            rel_tol: 1e-11,
            max_refinements: 9,
        },
    )
    .unwrap()
    .resistance;

    let mut previous_error = f64::INFINITY;
    for tolerance in [1e-3, 1e-5, 1e-7] {
        let result = wave_resistance_with(
            &hull,
            &cond,
            &WaveOptions {
                rel_tol: tolerance,
                max_refinements: 7,
            },
        )
        .unwrap();
        let actual = rel_err(result.resistance, reference);
        assert!(
            actual <= previous_error * 1.05,
            "tol={tolerance}: error grew"
        );
        assert!(
            actual < 5.0 * tolerance,
            "tol={tolerance}: actual relative error {actual:.3e}"
        );
        previous_error = actual;
    }
}

#[test]
fn reported_error_tracks_actual_error() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    for fn_ in [0.08, 0.20, 0.35, 0.50] {
        let speed = fn_ * (STANDARD_GRAVITY * hull.length()).sqrt();
        let cond = Conditions::freshwater(speed);
        let reference = wave_resistance_with(
            &hull,
            &cond,
            &WaveOptions {
                // Six forced refinements remain below the per-pass safety
                // cap even for Fn=0.08. Finer passes hit that cap and are not
                // valid references: their reported max_lambda retreats.
                rel_tol: 1e-14,
                max_refinements: 6,
            },
        )
        .unwrap();
        let result = wave_resistance_with(
            &hull,
            &cond,
            &WaveOptions {
                rel_tol: 1e-5,
                max_refinements: 7,
            },
        )
        .unwrap();
        assert!(
            reference.max_lambda >= 0.95 * result.max_lambda,
            "reference truncated earlier: reference λmax={:.3}, result λmax={:.3}",
            reference.max_lambda,
            result.max_lambda
        );
        let actual = rel_err(result.resistance, reference.resistance);
        let allowance = 10.0 * result.est_rel_error.max(1e-10);
        assert!(
            actual <= allowance,
            "Fn={fn_}: estimated {:.3e}, actual {actual:.3e}; result Rw={:.12e}, λmax={:.3}, evals={}; reference Rw={:.12e}, λmax={:.3}, evals={}",
            result.est_rel_error,
            result.resistance,
            result.max_lambda,
            result.inner_evaluations,
            reference.resistance,
            reference.max_lambda,
            reference.inner_evaluations,
        );
    }
}

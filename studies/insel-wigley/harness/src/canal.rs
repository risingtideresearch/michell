//! Independent finite-canal reference for the C2 Wigley hull.
//!
//! This module transcribes Insel's shallow-canal formulation rather than a
//! generic Sretensky kernel. The modal dispersion equations and finite-depth
//! vertical profile are equations (4.26), (4.29)--(4.31), printed page 48
//! (PDF 58). The monohull resistance weights are equations (4.40)--(4.41),
//! printed page 51 (PDF 61). A centered, symmetric catamaran multiplies each
//! demihull mode by `2 cos(n pi S/W)`, equations (4.42)--(4.48), printed pages
//! 52--53 (PDF 62--63), and its resistance is equation (4.50), printed page 54
//! (PDF 64).

use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
pub struct WigleyHull {
    pub length: f64,
    pub beam: f64,
    pub draft: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Flow {
    pub speed: f64,
    pub density: f64,
    pub gravity: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Canal {
    pub width: f64,
    pub depth: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ModalOptions {
    pub min_modes: usize,
    pub max_modes: usize,
    pub resistance_rel_tol: f64,
    pub interference_abs_tol: f64,
}

impl Default for ModalOptions {
    fn default() -> Self {
        Self {
            min_modes: 32,
            max_modes: 1_048_576,
            resistance_rel_tol: 5e-6,
            interference_abs_tol: 2e-6,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ModalResult {
    pub monohull_resistance: f64,
    pub catamaran_resistance: f64,
    pub interference: f64,
    pub modes: usize,
    pub resistance_rel_change: f64,
    pub interference_abs_change: f64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HistoricalVariant {
    Baseline,
    LiteralEquation429,
    DoubledCrossTerm,
    HalfNonzeroModeMultiplicity,
    DoubleNonzeroModeMultiplicity,
    CosineInsteadOfCosineSquared,
}

impl HistoricalVariant {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Baseline => "baseline_resolved_equations",
            Self::LiteralEquation429 => "literal_equation_4_29",
            Self::DoubledCrossTerm => "doubled_interference_cross_term",
            Self::HalfNonzeroModeMultiplicity => "half_nonzero_mode_multiplicity",
            Self::DoubleNonzeroModeMultiplicity => "double_nonzero_mode_multiplicity",
            Self::CosineInsteadOfCosineSquared => "cosine_instead_of_cosine_squared",
        }
    }
}

fn validate(hull: WigleyHull, flow: Flow, canal: Canal, separation: f64) -> Result<(), String> {
    let values = [
        hull.length,
        hull.beam,
        hull.draft,
        flow.speed,
        flow.density,
        flow.gravity,
        canal.width,
        canal.depth,
    ];
    if values
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err("hull, flow, and canal dimensions must be finite and positive".into());
    }
    if !separation.is_finite() || separation < 0.0 {
        return Err("separation must be finite and nonnegative".into());
    }
    if hull.draft >= canal.depth {
        return Err("the hull draft must be smaller than the canal depth".into());
    }
    let k0_h = flow.gravity * canal.depth / flow.speed.powi(2);
    if k0_h <= 1.0 {
        return Err("Insel equation (4.9) requires K0 H > 1".into());
    }
    Ok(())
}

fn sech_squared(value: f64) -> f64 {
    if value > 350.0 {
        0.0
    } else {
        1.0 / value.cosh().powi(2)
    }
}

fn two_x_over_sinh_two_x(value: f64) -> f64 {
    if value > 350.0 {
        0.0
    } else {
        2.0 * value / (2.0 * value).sinh()
    }
}

/// Solve Insel's equation (4.26) for the symmetric mode `m = 2n`.
fn mode_wavenumber(index: usize, k0: f64, width: f64, depth: f64) -> f64 {
    let transverse = 2.0 * PI * index as f64 / width;
    let residual = |wave_number: f64| {
        wave_number * wave_number
            - transverse * transverse
            - k0 * wave_number * (wave_number * depth).tanh()
    };
    let mut lower = if index == 0 {
        f64::EPSILON.sqrt() / depth
    } else {
        transverse
    };
    let mut upper = k0 + transverse;
    debug_assert!(residual(lower) <= 0.0);
    while residual(upper) <= 0.0 {
        upper *= 2.0;
    }
    for _ in 0..80 {
        let middle = 0.5 * (lower + upper);
        if residual(middle) <= 0.0 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    0.5 * (lower + upper)
}

/// `M_k(s) = integral_0^1 t^k exp(-s t) dt`, for k <= 2.
fn exponential_moments(s: f64) -> [f64; 3] {
    if s.abs() < 0.1 {
        let mut result = [0.0; 3];
        let mut factorial = 1.0;
        let mut power = 1.0;
        for order in 0..=28 {
            if order > 0 {
                factorial *= order as f64;
                power *= -s;
            }
            for (degree, value) in result.iter_mut().enumerate() {
                *value += power / (factorial * (degree + order + 1) as f64);
            }
        }
        return result;
    }
    let decay = (-s).exp();
    let m0 = -(-s).exp_m1() / s;
    let m1 = (m0 - decay) / s;
    let m2 = (2.0 * m1 - decay) / s;
    [m0, m1, m2]
}

/// Exact vertical integral of Insel's finite-depth factor over the Wigley
/// product parabola. Equation (4.29) contributes
/// `exp(-KH) cosh(K(H+z))`; this is not replaced by the deep-water kernel.
fn vertical_amplitude(wave_number: f64, depth: f64, draft: f64) -> f64 {
    let s = wave_number * draft;
    let moments = exponential_moments(s);
    let direct = moments[0] - moments[2];
    let image = (-wave_number * (2.0 * depth - draft)).exp() * (2.0 * moments[1] - moments[2]);
    0.5 * draft * (direct + image)
}

/// Exact longitudinal integral `integral f_x sin(w x) dx` for the Wigley
/// waterline parabola.
fn longitudinal_amplitude(wavenumber_x: f64, length: f64, beam: f64) -> f64 {
    let phase = 0.5 * wavenumber_x * length;
    if phase.abs() < 0.05 {
        let p2 = phase * phase;
        return 2.0
            * beam
            * phase
            * (-1.0 / 3.0 + p2 / 30.0 - p2 * p2 / 840.0 + p2.powi(3) / 45_360.0);
    }
    2.0 * beam * (phase * phase.cos() - phase.sin()) / phase.powi(2)
}

fn mode_contribution(
    hull: WigleyHull,
    flow: Flow,
    canal: Canal,
    separation: f64,
    index: usize,
    variant: HistoricalVariant,
) -> (f64, f64) {
    let k0 = flow.gravity / flow.speed.powi(2);
    let wave_number = mode_wavenumber(index, k0, canal.width, canal.depth);
    let transverse = 2.0 * PI * index as f64 / canal.width;
    let sin_theta = transverse / wave_number;
    let cos_squared = (1.0 - sin_theta * sin_theta).max(0.0);
    let wavenumber_x = wave_number * cos_squared.sqrt();
    let kh = wave_number * canal.depth;
    let denominator = 1.0 - k0 * canal.depth * sech_squared(kh) + sin_theta.powi(2);
    // Equation (4.25), printed page 47 (PDF 57), includes the dimensional
    // factor K0 + K cos^2(theta). It is absent from the printed definition of
    // tau_m in equation (4.29), but omitting it is dimensionally inconsistent
    // and fails the independently required wide/deep limit by 4 K0^2.
    let dimensional_factor = if variant == HistoricalVariant::LiteralEquation429 {
        1.0
    } else {
        k0 + wave_number * cos_squared
    };
    let source_factor = if index == 0 { -4.0 } else { -8.0 } * flow.speed.powi(2)
        / (canal.width * flow.gravity)
        * dimensional_factor;
    let eta = source_factor
        * longitudinal_amplitude(wavenumber_x, hull.length, hull.beam)
        * vertical_amplitude(wave_number, canal.depth, hull.draft)
        / denominator;
    let depth_ratio = two_x_over_sinh_two_x(kh);
    let resistance_weight = if index == 0 {
        1.0 - depth_ratio
    } else {
        1.0 - 0.5 * cos_squared * (1.0 + depth_ratio)
    };
    let prefactor = canal.width * flow.density * flow.gravity / 4.0;
    let multiplicity = match (variant, index) {
        (HistoricalVariant::HalfNonzeroModeMultiplicity, 1..) => 0.5,
        (HistoricalVariant::DoubleNonzeroModeMultiplicity, 1..) => 2.0,
        _ => 1.0,
    };
    let monohull = multiplicity * prefactor * eta.powi(2) * resistance_weight;
    let phase = PI * index as f64 * separation / canal.width;
    let catamaran_factor = match variant {
        HistoricalVariant::DoubledCrossTerm => 2.0 + 4.0 * (2.0 * phase).cos(),
        HistoricalVariant::CosineInsteadOfCosineSquared => 4.0 * phase.cos(),
        _ => 4.0 * phase.cos().powi(2),
    };
    (monohull, monohull * catamaran_factor)
}

pub fn fixed_mode_resistance(
    hull: WigleyHull,
    flow: Flow,
    canal: Canal,
    separation: f64,
    modes: usize,
) -> Result<ModalResult, String> {
    fixed_mode_resistance_variant(
        hull,
        flow,
        canal,
        separation,
        modes,
        HistoricalVariant::Baseline,
    )
}

pub fn fixed_mode_resistance_variant(
    hull: WigleyHull,
    flow: Flow,
    canal: Canal,
    separation: f64,
    modes: usize,
    variant: HistoricalVariant,
) -> Result<ModalResult, String> {
    validate(hull, flow, canal, separation)?;
    if modes == 0 {
        return Err("at least one mode is required".into());
    }
    let mut monohull = 0.0;
    let mut catamaran = 0.0;
    for index in 0..modes {
        let (mono_mode, cat_mode) =
            mode_contribution(hull, flow, canal, separation, index, variant);
        monohull += mono_mode;
        catamaran += cat_mode;
    }
    if !monohull.is_finite() || monohull <= 0.0 || !catamaran.is_finite() {
        return Err("finite-canal modal sum is not finite and positive".into());
    }
    Ok(ModalResult {
        monohull_resistance: monohull,
        catamaran_resistance: catamaran,
        interference: catamaran / (2.0 * monohull),
        modes,
        resistance_rel_change: f64::NAN,
        interference_abs_change: f64::NAN,
    })
}

/// Double the retained mode count until successive positive partial sums meet
/// both an absolute interference criterion and a relative monohull-resistance
/// criterion. Insel states only that high-angle contributions permit finite
/// truncation (printed page 51, PDF 61); these explicit criteria make that
/// qualitative rule reproducible.
pub fn converged_resistance(
    hull: WigleyHull,
    flow: Flow,
    canal: Canal,
    separation: f64,
    options: ModalOptions,
) -> Result<ModalResult, String> {
    converged_resistance_variant(
        hull,
        flow,
        canal,
        separation,
        options,
        HistoricalVariant::Baseline,
    )
}

pub fn converged_resistance_variant(
    hull: WigleyHull,
    flow: Flow,
    canal: Canal,
    separation: f64,
    options: ModalOptions,
    variant: HistoricalVariant,
) -> Result<ModalResult, String> {
    if options.min_modes == 0
        || options.max_modes < options.min_modes
        || !options.min_modes.is_power_of_two()
        || !options.max_modes.is_power_of_two()
        || options.resistance_rel_tol <= 0.0
        || options.interference_abs_tol <= 0.0
    {
        return Err("invalid modal convergence options".into());
    }
    let mut modes = options.min_modes;
    let mut previous =
        fixed_mode_resistance_variant(hull, flow, canal, separation, modes, variant)?;
    while modes < options.max_modes {
        modes *= 2;
        let mut current =
            fixed_mode_resistance_variant(hull, flow, canal, separation, modes, variant)?;
        current.resistance_rel_change =
            (current.monohull_resistance - previous.monohull_resistance).abs()
                / current.monohull_resistance;
        current.interference_abs_change = (current.interference - previous.interference).abs();
        if current.resistance_rel_change <= options.resistance_rel_tol
            && current.interference_abs_change <= options.interference_abs_tol
        {
            return Ok(current);
        }
        previous = current;
    }
    Err(format!(
        "modal sum did not converge by {} modes: resistance change {:.3e}, interference change {:.3e}",
        previous.modes, previous.resistance_rel_change, previous.interference_abs_change
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use michell::{
        hulls, multihull_wave_resistance_with, Conditions, Fluid, Placement, WaveOptions,
    };

    const LENGTH: f64 = 1.8;
    const BEAM: f64 = 0.18;
    const DRAFT: f64 = 0.1125;
    const DENSITY: f64 = 1000.0;
    const GRAVITY: f64 = 9.80665;

    const GL_NODES: [f64; 8] = [
        0.095_012_509_837_637_44,
        0.281_603_550_779_258_9,
        0.458_016_777_657_227_4,
        0.617_876_244_402_643_8,
        0.755_404_408_355_003,
        0.865_631_202_387_831_8,
        0.944_575_023_073_232_6,
        0.989_400_934_991_649_9,
    ];
    const GL_WEIGHTS: [f64; 8] = [
        0.189_450_610_455_068_5,
        0.182_603_415_044_923_6,
        0.169_156_519_395_002_54,
        0.149_595_988_816_576_73,
        0.124_628_971_255_533_87,
        0.095_158_511_682_492_78,
        0.062_253_523_938_647_89,
        0.027_152_459_411_754_095,
    ];

    fn flow(fn_: f64) -> Flow {
        Flow {
            speed: fn_ * (GRAVITY * LENGTH).sqrt(),
            density: DENSITY,
            gravity: GRAVITY,
        }
    }

    fn independent_deep_water(separation: f64, flow: Flow) -> (f64, f64) {
        let nu = flow.gravity / flow.speed.powi(2);
        let lambda_max = 500.0;
        let phase_step = PI / (2.0 * nu * LENGTH / 2.0);
        let panels = ((lambda_max - 1.0) / phase_step).ceil() as usize;
        let mut mono_integral = 0.0;
        let mut cat_integral = 0.0;
        for panel in 0..panels {
            let u0 = panel as f64 * phase_step;
            let u1 = ((panel + 1) as f64 * phase_step).min(lambda_max - 1.0);
            let t0 = u0.sqrt();
            let t1 = u1.sqrt();
            let half = 0.5 * (t1 - t0);
            let middle = 0.5 * (t1 + t0);
            for (&node, &weight) in GL_NODES.iter().zip(&GL_WEIGHTS) {
                for sign in [-1.0, 1.0] {
                    let t = middle + sign * half * node;
                    let lambda = 1.0 + t * t;
                    let kx = nu * lambda;
                    let vertical_k = nu * lambda * lambda;
                    let x = 2.0
                        * ((kx * LENGTH / 2.0).sin() / kx.powi(2)
                            - LENGTH / 2.0 * (kx * LENGTH / 2.0).cos() / kx);
                    let decay = (-vertical_k * DRAFT).exp();
                    let z = (1.0 - decay) / vertical_k
                        - (2.0
                            - decay
                                * (vertical_k.powi(2) * DRAFT.powi(2)
                                    + 2.0 * vertical_k * DRAFT
                                    + 2.0))
                            / (vertical_k.powi(3) * DRAFT.powi(2));
                    let amplitude = -(4.0 * BEAM / LENGTH.powi(2)) * x * z;
                    let jacobian = 2.0 * lambda.powi(2) / (2.0 + t * t).sqrt();
                    let contribution = weight * half * amplitude.powi(2) * jacobian;
                    mono_integral += contribution;
                    let half_phase = 0.5 * nu * separation * lambda * t * (2.0 + t * t).sqrt();
                    cat_integral += contribution * 4.0 * half_phase.cos().powi(2);
                }
            }
        }
        let prefactor = 4.0 * flow.density * flow.gravity.powi(2) / (PI * flow.speed.powi(2));
        let mono = prefactor * mono_integral;
        (mono, prefactor * cat_integral / (2.0 * mono))
    }

    fn library_reference(separation: f64, flow: Flow) -> (f64, f64) {
        let hull = hulls::wigley(LENGTH, BEAM, DRAFT).unwrap();
        let conditions = Conditions {
            speed: flow.speed,
            fluid: Fluid {
                density: flow.density,
                kinematic_viscosity: 1.141e-6,
            },
            gravity: flow.gravity,
        };
        let options = WaveOptions {
            rel_tol: 1e-8,
            max_refinements: 8,
        };
        let solo =
            multihull_wave_resistance_with(&[(&hull, Placement::default())], &conditions, &options)
                .unwrap()
                .resistance;
        let pair = multihull_wave_resistance_with(
            &[
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
            ],
            &conditions,
            &options,
        )
        .unwrap()
        .resistance;
        (solo, pair / (2.0 * solo))
    }

    fn relative_error(actual: f64, expected: f64) -> f64 {
        (actual - expected).abs() / expected.abs()
    }

    /// G1: as both canal dimensions grow, equations (4.26)--(4.41) must
    /// recover the independent analytic-Wigley Michell integral and the
    /// library implementation. Only convergence assertions are exposed; the
    /// W=3.7 m study result is intentionally not computed here.
    #[test]
    fn g1_wide_deep_canal_converges_to_unbounded_references() {
        let hull = WigleyHull {
            length: LENGTH,
            beam: BEAM,
            draft: DRAFT,
        };
        for fn_ in [0.25, 0.35, 0.50] {
            let flow = flow(fn_);
            for separation_ratio in [0.2, 0.3, 0.4, 0.5] {
                let separation = separation_ratio * LENGTH;
                let independent = independent_deep_water(separation, flow);
                let library = library_reference(separation, flow);
                assert!(relative_error(independent.0, library.0) < 5e-4);
                assert!((independent.1 - library.1).abs() < 2e-3);

                let mut errors = Vec::new();
                for scale in [5.0, 20.0, 80.0] {
                    let result = converged_resistance(
                        hull,
                        flow,
                        Canal {
                            width: 3.7 * scale,
                            depth: 1.85 * scale,
                        },
                        separation,
                        ModalOptions::default(),
                    )
                    .unwrap_or_else(|error| {
                        panic!("Fn={fn_}, S/L={separation_ratio}, scale={scale}: {error}")
                    });
                    errors.push((
                        relative_error(result.monohull_resistance, independent.0),
                        (result.interference - independent.1).abs(),
                    ));
                }
                assert!(
                    errors[2].0 < 5e-3,
                    "Fn={fn_}, S/L={separation_ratio}: {errors:?}"
                );
                assert!(
                    errors[2].1 < 0.02,
                    "Fn={fn_}, S/L={separation_ratio}: {errors:?}"
                );
                // Finite-width corrections oscillate with W, so convergence
                // is not pointwise monotone. The 80x endpoint is the gate.
            }
        }
    }

    /// G2: doubling a deliberately finite tank-dimension sum changes both
    /// observables by less than one tenth of the 0.02 comparison tolerance.
    #[test]
    fn g2_tank_dimension_mode_doubling_is_below_gate() {
        let hull = WigleyHull {
            length: LENGTH,
            beam: BEAM,
            draft: DRAFT,
        };
        let canal = Canal {
            width: 3.7,
            depth: 1.85,
        };
        for fn_ in [0.25, 0.35, 0.50] {
            for separation_ratio in [0.2, 0.3, 0.4, 0.5] {
                let coarse =
                    fixed_mode_resistance(hull, flow(fn_), canal, separation_ratio * LENGTH, 128)
                        .unwrap();
                let fine =
                    fixed_mode_resistance(hull, flow(fn_), canal, separation_ratio * LENGTH, 256)
                        .unwrap();
                assert!(
                    relative_error(coarse.monohull_resistance, fine.monohull_resistance) < 0.002
                );
                assert!((coarse.interference - fine.interference).abs() < 0.002);
            }
        }
    }
}

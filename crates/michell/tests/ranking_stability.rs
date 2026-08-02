//! Design-decision regressions: nearby hull variants must keep the same
//! ordering across solver tolerances, numerical routes, and spline refinement.

use michell::{
    multihull_wave_resistance_with, wave_resistance_with, BSplineSurface, Conditions, Hull,
    Placement, WaveMethod, WaveOptions, WaveOutcome, WaveResistance, STANDARD_GRAVITY,
};

const LENGTH: f64 = 10.0;
const BEAM: f64 = 1.0;
const DRAFT: f64 = 0.625;

struct Variant {
    name: &'static str,
    hull: Hull,
}

fn polynomial_product(left: &[f64], right: &[f64]) -> Vec<f64> {
    let mut product = vec![0.0; left.len() + right.len() - 1];
    for (i, &a) in left.iter().enumerate() {
        for (j, &b) in right.iter().enumerate() {
            product[i + j] += a * b;
        }
    }
    product
}

fn binomial(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    (0..k).fold(1.0, |value, index| {
        value * (n - index) as f64 / (index + 1) as f64
    })
}

/// Power-basis polynomial to degree-`degree` Bernstein controls.
fn bernstein_controls(power: &[f64], degree: usize) -> Vec<f64> {
    (0..=degree)
        .map(|i| {
            power
                .iter()
                .enumerate()
                .take(i + 1)
                .map(|(k, coefficient)| coefficient * binomial(i, k) / binomial(degree, k))
                .sum()
        })
        .collect()
}

fn wigley_variant(beam_scale: f64, fullness: f64, lcb_shift: f64) -> Hull {
    // t = (x + L/2)/L. The Wigley longitudinal factor is g = 4t(1-t).
    let wigley = vec![0.0, 4.0, -4.0];
    let squared = polynomial_product(&wigley, &wigley);
    let shifted = polynomial_product(&wigley, &[1.0 - lcb_shift, 2.0 * lcb_shift]);
    let mut longitudinal = vec![0.0; 5];
    for (index, value) in wigley.iter().enumerate() {
        longitudinal[index] += value;
    }
    for (index, value) in squared.iter().enumerate() {
        longitudinal[index] += fullness * value;
    }
    for (index, value) in shifted.iter().enumerate() {
        longitudinal[index] += value - wigley.get(index).copied().unwrap_or(0.0);
    }

    let x_control = bernstein_controls(&longitudinal, 4);
    let z_control = [1.0, 1.0, 0.0];
    let mut control = Vec::with_capacity(x_control.len() * z_control.len());
    for x_value in x_control {
        for z_value in z_control {
            control.push(0.5 * BEAM * beam_scale * x_value * z_value);
        }
    }
    let half_length = 0.5 * LENGTH;
    Hull::new(
        BSplineSurface::new(
            4,
            2,
            vec![
                -half_length,
                -half_length,
                -half_length,
                -half_length,
                -half_length,
                half_length,
                half_length,
                half_length,
                half_length,
                half_length,
            ],
            vec![0.0, 0.0, 0.0, DRAFT, DRAFT, DRAFT],
            control,
        )
        .unwrap(),
    )
    .unwrap()
}

fn variants() -> Vec<Variant> {
    vec![
        Variant {
            name: "base",
            hull: wigley_variant(1.0, 0.0, 0.0),
        },
        Variant {
            name: "beam_narrow",
            hull: wigley_variant(0.985, 0.0, 0.0),
        },
        Variant {
            name: "beam_wide",
            hull: wigley_variant(1.015, 0.0, 0.0),
        },
        Variant {
            name: "fullness_fine",
            hull: wigley_variant(1.0, -0.025, 0.0),
        },
        Variant {
            name: "fullness_full",
            hull: wigley_variant(1.0, 0.025, 0.0),
        },
        Variant {
            name: "lcb_forward",
            hull: wigley_variant(1.0, 0.0, 0.04),
        },
    ]
}

fn ordering(results: &[WaveResistance]) -> Vec<usize> {
    let mut indices: Vec<_> = (0..results.len()).collect();
    indices.sort_by(|&left, &right| {
        results[left]
            .resistance
            .total_cmp(&results[right].resistance)
    });
    indices
}

fn evaluate(
    variants: &[Variant],
    conditions: &Conditions,
    options: &WaveOptions,
    force_marcher: bool,
) -> Vec<WaveResistance> {
    variants
        .iter()
        .map(|variant| {
            if force_marcher {
                multihull_wave_resistance_with(
                    &[(&variant.hull, Placement::default())],
                    conditions,
                    options,
                )
                .unwrap()
            } else {
                wave_resistance_with(&variant.hull, conditions, options).unwrap()
            }
        })
        .collect()
}

fn insert_curve_knot(control: &[f64], knots: &[f64], degree: usize, knot: f64) -> Vec<f64> {
    let span = knots
        .windows(2)
        .position(|window| window[0] <= knot && knot < window[1])
        .unwrap();
    let multiplicity = knots.iter().filter(|&&value| value == knot).count();
    assert!(multiplicity < degree);
    let last_control = control.len() - 1;
    let mut refined = vec![0.0; control.len() + 1];
    refined[..=span - degree].copy_from_slice(&control[..=span - degree]);
    for index in span - multiplicity..=last_control {
        refined[index + 1] = control[index];
    }
    for index in span - degree + 1..=span - multiplicity {
        let alpha = (knot - knots[index]) / (knots[index + degree] - knots[index]);
        refined[index] = alpha * control[index] + (1.0 - alpha) * control[index - 1];
    }
    refined
}

fn insert_surface_knot_x(surface: &BSplineSurface, knot: f64) -> BSplineSurface {
    let old_nx = surface.n_ctrl_x();
    let nz = surface.n_ctrl_z();
    let mut control = vec![0.0; (old_nx + 1) * nz];
    for z in 0..nz {
        let column: Vec<_> = (0..old_nx)
            .map(|x| surface.control()[x * nz + z])
            .collect();
        let refined = insert_curve_knot(&column, surface.knots_x(), surface.degree_x(), knot);
        for (x, value) in refined.into_iter().enumerate() {
            control[x * nz + z] = value;
        }
    }
    let mut knots_x = surface.knots_x().to_vec();
    let insertion = knots_x.partition_point(|value| *value <= knot);
    knots_x.insert(insertion, knot);
    BSplineSurface::new(
        surface.degree_x(),
        surface.degree_z(),
        knots_x,
        surface.knots_z().to_vec(),
        control,
    )
    .unwrap()
}

fn insert_surface_knot_z(surface: &BSplineSurface, knot: f64) -> BSplineSurface {
    let nx = surface.n_ctrl_x();
    let old_nz = surface.n_ctrl_z();
    let mut control = Vec::with_capacity(nx * (old_nz + 1));
    for x in 0..nx {
        let row = &surface.control()[x * old_nz..(x + 1) * old_nz];
        control.extend(insert_curve_knot(
            row,
            surface.knots_z(),
            surface.degree_z(),
            knot,
        ));
    }
    let mut knots_z = surface.knots_z().to_vec();
    let insertion = knots_z.partition_point(|value| *value <= knot);
    knots_z.insert(insertion, knot);
    BSplineSurface::new(
        surface.degree_x(),
        surface.degree_z(),
        surface.knots_x().to_vec(),
        knots_z,
        control,
    )
    .unwrap()
}

fn refined_variant(variant: &Variant) -> Variant {
    let x_refined = insert_surface_knot_x(variant.hull.surface(), 0.0);
    let refined = insert_surface_knot_z(&x_refined, 0.5 * DRAFT);
    Variant {
        name: variant.name,
        hull: Hull::new(refined).unwrap(),
    }
}

#[test]
fn ordering_is_stable_across_requested_tolerances() {
    let variants = variants();
    let speed = 0.30 * (STANDARD_GRAVITY * LENGTH).sqrt();
    let conditions = Conditions::freshwater(speed);
    let mut reference_order = None;

    for rel_tol in [1e-4, 1e-6, 1e-8] {
        let options = WaveOptions {
            rel_tol,
            max_refinements: 8,
        };
        let results = evaluate(&variants, &conditions, &options, false);
        assert!(results.iter().all(|result| result.outcome.is_converged()));
        let order = ordering(&results);
        if let Some(reference) = &reference_order {
            assert_eq!(
                &order, reference,
                "ordering changed at rel_tol={rel_tol}: {}",
                order
                    .iter()
                    .map(|&index| variants[index].name)
                    .collect::<Vec<_>>()
                    .join(", "),
            );
        } else {
            reference_order = Some(order);
        }
    }
}

#[test]
fn endpoint_and_marcher_routes_agree_on_order_and_pairwise_margins() {
    let variants = variants();
    let speed = 0.05 * (STANDARD_GRAVITY * LENGTH).sqrt();
    let conditions = Conditions::freshwater(speed);
    let options = WaveOptions {
        rel_tol: 1e-6,
        max_refinements: 8,
    };
    let dispatched = evaluate(&variants, &conditions, &options, false);
    let marched = evaluate(&variants, &conditions, &options, true);
    assert!(dispatched.iter().all(|result| {
        result.method == WaveMethod::EndpointReduction && result.outcome == WaveOutcome::Converged
    }));
    assert!(marched.iter().all(|result| {
        result.method == WaveMethod::GeneralMarcher && result.outcome == WaveOutcome::Converged
    }));
    assert_eq!(ordering(&dispatched), ordering(&marched));

    for left in 0..variants.len() {
        for right in left + 1..variants.len() {
            let endpoint_margin = dispatched[left].resistance - dispatched[right].resistance;
            let marcher_margin = marched[left].resistance - marched[right].resistance;
            assert_eq!(
                endpoint_margin.is_sign_positive(),
                marcher_margin.is_sign_positive(),
                "{} vs {} changed order",
                variants[left].name,
                variants[right].name,
            );
            let combined_error = dispatched[left].est_rel_error
                * dispatched[left].resistance.abs()
                + dispatched[right].est_rel_error * dispatched[right].resistance.abs()
                + marched[left].est_rel_error * marched[left].resistance.abs()
                + marched[right].est_rel_error * marched[right].resistance.abs();
            let margin_difference = (endpoint_margin - marcher_margin).abs();
            assert!(
                margin_difference <= 1.05 * combined_error,
                "{} vs {}: endpoint margin={endpoint_margin:.12e}, marcher margin={marcher_margin:.12e}, difference={margin_difference:.3e}, combined estimate={combined_error:.3e}",
                variants[left].name,
                variants[right].name,
            );
        }
    }
}

#[test]
fn ordering_is_invariant_under_exact_knot_insertion() {
    let variants = variants();
    let refined: Vec<_> = variants.iter().map(refined_variant).collect();
    let speed = 0.30 * (STANDARD_GRAVITY * LENGTH).sqrt();
    let conditions = Conditions::freshwater(speed);
    let options = WaveOptions {
        rel_tol: 1e-8,
        max_refinements: 8,
    };

    for (original, refined) in variants.iter().zip(&refined) {
        for x_fraction in [0.0, 0.17, 0.5, 0.83, 1.0] {
            let x = -0.5 * LENGTH + x_fraction * LENGTH;
            for z_fraction in [0.0, 0.31, 0.72, 1.0] {
                let z = z_fraction * DRAFT;
                let scale = original.hull.surface().eval(x, z).abs().max(1.0);
                assert!(
                    (original.hull.surface().eval(x, z) - refined.hull.surface().eval(x, z))
                        .abs()
                        <= 2e-14 * scale,
                    "{} geometry changed at ({x}, {z})",
                    original.name,
                );
            }
        }
    }

    let original_results = evaluate(&variants, &conditions, &options, false);
    let refined_results = evaluate(&refined, &conditions, &options, false);
    assert_eq!(ordering(&original_results), ordering(&refined_results));
}

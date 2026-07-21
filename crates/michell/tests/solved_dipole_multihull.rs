//! Solved-dipole multihull: fleets whose asymmetric members carry a *solved*
//! centreplane-lifting dipole, superposed with placement phases
//! (`multihull_wave_resistance_lifting`).
//!
//! The care-point is phase-origin consistency — each member's dipole `G_j` is
//! re-centred at its hull midpoint so it shares the source's placement phase.
//! The tests pin that: a fleet of one reproduces the single-hull result, a
//! symmetric fleet reproduces the source-only multihull, a mirror-image fleet
//! is invariant, and well-separated members become additive.

use michell::{
    asymmetric_wave_resistance_lifting, multihull_wave_resistance_lifting,
    multihull_wave_resistance_with, BSplineSurface, Conditions, Hull, LiftingGrid, Placement,
};

fn scaled(s: &BSplineSurface, k: f64) -> BSplineSurface {
    BSplineSurface::new(
        s.degree_x(),
        s.degree_z(),
        s.knots_x().to_vec(),
        s.knots_z().to_vec(),
        s.control().iter().map(|c| k * c).collect(),
    )
    .unwrap()
}

/// `H` (port 0.5, stbd 1.0 ⇒ f_a > 0), its mirror `H_m` (port/stbd swapped),
/// and the symmetric mean hull `f_sym`.
fn variants() -> (Hull, Hull, Hull) {
    let base = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let stb = scaled(base.surface(), 1.0);
    let port = scaled(base.surface(), 0.5);
    let h = Hull::new_asymmetric(scaled(&port, 1.0), scaled(&stb, 1.0)).unwrap();
    let h_m = Hull::new_asymmetric(scaled(&stb, 1.0), scaled(&port, 1.0)).unwrap();
    let mean = Hull::new(scaled(h.surface(), 1.0)).unwrap();
    (h, h_m, mean)
}

const GRID: LiftingGrid = LiftingGrid { nx: 24, nz: 8 };

fn cond() -> Conditions {
    Conditions::seawater(3.0)
}

/// A fleet of one asymmetric hull is exactly the single-hull lifting result.
#[test]
fn fleet_of_one_matches_single() {
    let (h, _, _) = variants();
    let c = cond();
    let single = asymmetric_wave_resistance_lifting(&h, &c, &Default::default(), GRID)
        .unwrap()
        .resistance;
    let fleet = multihull_wave_resistance_lifting(
        &[(&h, Placement::default())],
        &c,
        &Default::default(),
        GRID,
    )
    .unwrap()
    .resistance;
    assert!((single - fleet).abs() < 1e-10 * single, "{single} vs {fleet}");
}

/// A fleet of purely symmetric hulls has no dipole, so the lifting path
/// reproduces the source-only strip multihull exactly.
#[test]
fn symmetric_fleet_matches_source_multihull() {
    let base = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let c = cond();
    let members = [
        (&base, Placement { x: 0.0, y: 1.5 }),
        (&base, Placement { x: 1.0, y: -1.5 }),
    ];
    let strip = multihull_wave_resistance_with(&members, &c, &Default::default())
        .unwrap()
        .resistance;
    let lift = multihull_wave_resistance_lifting(&members, &c, &Default::default(), GRID)
        .unwrap()
        .resistance;
    assert!((strip - lift).abs() < 1e-9 * strip, "{strip} vs {lift}");
}

/// A staggered asymmetric fleet and its mirror image about y = 0 are the same
/// physical flow: mirroring flips each hull's y **and** swaps port/starboard
/// (`H → H_m`). Equality confirms the dipole placement phases are consistent.
#[test]
fn mirror_image_fleet_is_invariant() {
    let (h, h_m, _) = variants();
    let c = cond();
    let s = 1.4;
    let original = [
        (&h, Placement { x: 0.0, y: s }),
        (&h, Placement { x: 1.7, y: -s }),
    ];
    let mirror = [
        (&h_m, Placement { x: 0.0, y: -s }),
        (&h_m, Placement { x: 1.7, y: s }),
    ];
    let ro = multihull_wave_resistance_lifting(&original, &c, &Default::default(), GRID)
        .unwrap()
        .resistance;
    let rm = multihull_wave_resistance_lifting(&mirror, &c, &Default::default(), GRID)
        .unwrap()
        .resistance;
    assert!((ro - rm).abs() < 1e-8 * ro, "{ro} vs mirror {rm}");
}

/// Well-separated members do not interfere: the combined resistance approaches
/// the sum of the standalone lifting resistances.
#[test]
fn far_apart_members_are_additive() {
    let (h, _, _) = variants();
    let c = cond();
    let single = asymmetric_wave_resistance_lifting(&h, &c, &Default::default(), GRID)
        .unwrap()
        .resistance;
    let far = multihull_wave_resistance_lifting(
        &[
            (&h, Placement { x: 0.0, y: 40.0 }),
            (&h, Placement { x: 0.0, y: -40.0 }),
        ],
        &c,
        &Default::default(),
        GRID,
    )
    .unwrap()
    .resistance;
    assert!(
        (far - 2.0 * single).abs() < 0.03 * 2.0 * single,
        "far {far} vs 2×single {}",
        2.0 * single
    );
}

/// An asymmetric-demihull catamaran differs from its symmetric-mean catamaran:
/// the solved dipole systems are active and (net) raise the resistance.
#[test]
fn asymmetric_catamaran_differs_from_symmetric() {
    let (h, h_m, mean) = variants();
    let c = cond();
    let s = 1.5;
    let asym = [
        (&h, Placement { x: 0.0, y: s }),
        (&h_m, Placement { x: 0.0, y: -s }),
    ];
    let sym = [
        (&mean, Placement { x: 0.0, y: s }),
        (&mean, Placement { x: 0.0, y: -s }),
    ];
    let ra = multihull_wave_resistance_lifting(&asym, &c, &Default::default(), GRID)
        .unwrap()
        .resistance;
    let rs = multihull_wave_resistance_lifting(&sym, &c, &Default::default(), GRID)
        .unwrap()
        .resistance;
    assert!(ra > rs * (1.0 + 1e-3), "asym {ra} should exceed sym {rs}");
}

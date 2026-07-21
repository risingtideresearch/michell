//! Asymmetric wave resistance with the dipole system driven by the *solved*
//! centreplane lifting distribution (`asymmetric_wave_resistance_lifting`)
//! rather than the prescribed strip closure `μ = 2U f_a`.
//!
//! The solved doublet density is a genuine lifting-surface potential jump, of
//! the same order as — but a different chordwise/vertical *shape* than — the
//! crude `2 f_a`, so the corrected dipole is positive, scales like `f_a²`, and
//! sits within an O(1) factor of the strip estimate (the ratio is
//! speed-dependent, since the two shapes weight the free-wave spectrum
//! differently). For a 2-D flat plate `μ_lift/μ_strip = π/2`, the archetype of
//! this correction.

use michell::{
    asymmetric_wave_resistance_lifting, wave_resistance, BSplineSurface, Conditions, Hull,
    LiftingGrid,
};

/// Rebuild a surface with its control net scaled by `k`.
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

fn clone_surface(s: &BSplineSurface) -> BSplineSurface {
    scaled(s, 1.0)
}

/// Build an asymmetric hull (full-beam starboard, `k`-scaled port) and return
/// it together with its symmetric mean hull `f_sym`.
fn asymmetric(k: f64) -> (Hull, Hull) {
    let base = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let starboard = clone_surface(base.surface());
    let port = scaled(base.surface(), k);
    let asym = Hull::new_asymmetric(clone_surface(&port), clone_surface(&starboard)).unwrap();
    // The asymmetric hull's own surface() is the symmetric mean f_sym.
    let mean = Hull::new(clone_surface(asym.surface())).unwrap();
    (asym, mean)
}

/// The solved lifting dipole is a positive addition to the source resistance,
/// of the same order of magnitude as the strip-closure dipole.
#[test]
fn lifting_dipole_is_positive_and_comparable_to_strip() {
    let (asym, mean) = asymmetric(0.5);
    let grid = LiftingGrid { nx: 36, nz: 10 };
    for u in [2.0, 3.0, 4.0] {
        let cond = Conditions::seawater(u);
        let r_source = wave_resistance(&mean, &cond).unwrap().resistance;
        let r_strip = wave_resistance(&asym, &cond).unwrap().resistance;
        let r_lift = asymmetric_wave_resistance_lifting(&asym, &cond, &Default::default(), grid)
            .unwrap()
            .resistance;
        let dip_lift = r_lift - r_source;
        let dip_strip = r_strip - r_source;
        assert!(dip_lift > 0.0, "u={u}: dipole {dip_lift} must be positive");
        // Same order as the strip estimate (ratio is speed-dependent).
        let ratio = dip_lift / dip_strip;
        assert!(
            (0.2..12.0).contains(&ratio),
            "u={u}: lift/strip dipole ratio {ratio} out of range"
        );
    }
}

/// The dipole resistance scales like `f_a²`: the centreplane solve is linear in
/// the camber forcing `∂f_a/∂x`, so μ ∝ f_a, G ∝ f_a and R_dipole ∝ f_a².
/// Halving the antisymmetric part quarters the dipole resistance.
#[test]
fn dipole_scales_quadratically_with_asymmetry() {
    // k = 0.5 ⇒ f_a = 0.25·base; k = 0.75 ⇒ f_a = 0.125·base (half).
    let cond = Conditions::seawater(3.0);
    let grid = LiftingGrid { nx: 40, nz: 12 };
    let dipole = |k: f64| {
        let (asym, mean) = asymmetric(k);
        let r_lift = asymmetric_wave_resistance_lifting(&asym, &cond, &Default::default(), grid)
            .unwrap()
            .resistance;
        let r_source = wave_resistance(&mean, &cond).unwrap().resistance;
        r_lift - r_source
    };
    let big = dipole(0.5);
    let small = dipole(0.75);
    assert!(big > 0.0 && small > 0.0);
    assert!(
        ((small / big) - 0.25).abs() < 0.02,
        "f_a² scaling: {small}/{big} = {} vs 0.25",
        small / big
    );
}

/// A symmetric hull has no camber system, so the lifting path is rejected.
#[test]
fn symmetric_hull_is_rejected() {
    let hull = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let cond = Conditions::seawater(3.0);
    assert!(
        asymmetric_wave_resistance_lifting(&hull, &cond, &Default::default(), LiftingGrid::default())
            .is_err()
    );
}

/// Mirroring the hull (`f_a → −f_a`) flips the sign of the solved μ but not
/// `|G|²`, so the resistance is unchanged — the solved dipole obeys the same
/// additive separation as the strip closure.
#[test]
fn mirror_image_has_equal_resistance() {
    let base = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let starboard = clone_surface(base.surface());
    let port = scaled(base.surface(), 0.4);
    let asym = Hull::new_asymmetric(clone_surface(&port), clone_surface(&starboard)).unwrap();
    let mirror = Hull::new_asymmetric(clone_surface(&starboard), clone_surface(&port)).unwrap();
    let cond = Conditions::seawater(3.0);
    let grid = LiftingGrid { nx: 28, nz: 8 };
    let ra = asymmetric_wave_resistance_lifting(&asym, &cond, &Default::default(), grid)
        .unwrap()
        .resistance;
    let rm = asymmetric_wave_resistance_lifting(&mirror, &cond, &Default::default(), grid)
        .unwrap()
        .resistance;
    assert!((ra - rm).abs() < 1e-7 * ra, "{ra} vs mirror {rm}");
}

/// The corrected resistance settles as the centreplane grid is refined.
#[test]
fn grid_convergence() {
    let (asym, _) = asymmetric(0.5);
    let cond = Conditions::seawater(3.0);
    let r = |nx, nz| {
        asymmetric_wave_resistance_lifting(&asym, &cond, &Default::default(), LiftingGrid { nx, nz })
            .unwrap()
            .resistance
    };
    let coarse = r(16, 5);
    let medium = r(28, 8);
    let fine = r(44, 12);
    assert!(
        (fine - medium).abs() < (medium - coarse).abs(),
        "not settling: {coarse}, {medium}, {fine}"
    );
    assert!((fine - medium).abs() < 0.05 * fine, "fine grid not converged");
}

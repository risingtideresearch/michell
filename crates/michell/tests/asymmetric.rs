//! Asymmetric (port ≠ starboard) thin-ship hulls: the antisymmetric half-beam
//! radiates an extra centreplane-dipole wave system on top of the classical
//! Michell source system. These tests pin the two guarantees the extension
//! must keep:
//!
//! 1. **Reduction.** An asymmetric hull whose two sides are identical is, to
//!    full floating-point precision, the symmetric [`Hull::new`] hull — the
//!    dipole term is exactly zero.
//! 2. **Additive separation.** Mirroring the hull (`f_a → −f_a`) leaves the
//!    resistance unchanged, proving the source–dipole cross term integrates to
//!    zero and `R_w = R_source(f_sym) + R_dipole(f_a)`.

use michell::{hulls, wave_resistance, BSplineSurface, Conditions, Hull};

/// Rebuild a surface with its control net scaled by `k` (a `k`-times narrower
/// or wider half-beam on the same parametrisation).
fn scaled(surface: &BSplineSurface, k: f64) -> BSplineSurface {
    BSplineSurface::new(
        surface.degree_x(),
        surface.degree_z(),
        surface.knots_x().to_vec(),
        surface.knots_z().to_vec(),
        surface.control().iter().map(|c| k * c).collect(),
    )
    .unwrap()
}

fn clone_surface(surface: &BSplineSurface) -> BSplineSurface {
    scaled(surface, 1.0)
}

/// Identical port and starboard sides ⇒ bit-for-bit the symmetric hull, in
/// both resistance and hydrostatics (the dipole net is all zeros).
#[test]
fn symmetric_sides_reduce_to_michell() {
    let base = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let asym =
        Hull::new_asymmetric(clone_surface(base.surface()), clone_surface(base.surface())).unwrap();

    // Geometry is unchanged.
    let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(1e-300);
    assert!(rel(asym.displaced_volume(), base.displaced_volume()) < 1e-12);
    assert!(rel(asym.wetted_surface(), base.wetted_surface()) < 1e-12);
    assert!(rel(asym.waterplane_area(), base.waterplane_area()) < 1e-12);
    assert!(
        rel(
            asym.waterplane_transverse_moment(),
            base.waterplane_transverse_moment()
        ) < 1e-12
    );

    // Wave resistance is unchanged across the whole Froude range.
    for u in [1.5, 2.5, 3.5, 5.0] {
        let cond = Conditions::seawater(u);
        let rb = wave_resistance(&base, &cond).unwrap().resistance;
        let ra = wave_resistance(&asym, &cond).unwrap().resistance;
        assert!(rel(ra, rb) < 1e-10, "u={u}: sym {rb} vs asym {ra}");
    }
}

/// A genuinely asymmetric hull: full-beam starboard, half-beam port. Its mean
/// (thickness) side is the 0.75-scaled symmetric hull; the antisymmetric part
/// is the 0.25-scaled half-beam. The dipole term must add strictly positive
/// resistance on top of the mean hull's source resistance.
#[test]
fn asymmetry_adds_wave_resistance() {
    let base = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let starboard = clone_surface(base.surface()); // f₊ = 1.00 × base
    let port = scaled(base.surface(), 0.5); //        f₋ = 0.50 × base
    let asym = Hull::new_asymmetric(clone_surface(&port), clone_surface(&starboard)).unwrap();

    // The source part alone is the symmetric hull with f_sym = 0.75 × base.
    let mean = Hull::new(scaled(base.surface(), 0.75)).unwrap();

    for u in [2.0, 3.0, 4.5] {
        let cond = Conditions::seawater(u);
        let r_mean = wave_resistance(&mean, &cond).unwrap().resistance;
        let r_asym = wave_resistance(&asym, &cond).unwrap().resistance;
        assert!(
            r_asym > r_mean * (1.0 + 1e-6),
            "u={u}: dipole should add resistance, mean {r_mean} vs asym {r_asym}"
        );
    }
}

/// Mirroring the hull sends `f_a → −f_a` while leaving `f_sym` fixed. If any
/// source–dipole cross term survived the θ-integration the resistance would
/// change sign-dependently; equality to full precision proves the additive
/// separation `R_w = R_source + R_dipole`.
#[test]
fn mirror_image_has_equal_resistance() {
    let base = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let starboard = clone_surface(base.surface());
    let port = scaled(base.surface(), 0.4);

    let asym = Hull::new_asymmetric(clone_surface(&port), clone_surface(&starboard)).unwrap();
    let mirror = Hull::new_asymmetric(clone_surface(&starboard), clone_surface(&port)).unwrap();

    let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(1e-300);
    for u in [2.0, 3.0, 4.5] {
        let cond = Conditions::seawater(u);
        let ra = wave_resistance(&asym, &cond).unwrap().resistance;
        let rm = wave_resistance(&mirror, &cond).unwrap().resistance;
        assert!(rel(ra, rm) < 1e-10, "u={u}: {ra} vs mirror {rm}");
    }
}

/// Contract violations are rejected rather than silently mis-decomposed.
#[test]
fn mismatched_parametrisations_are_rejected() {
    let a = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let b = hulls::wigley(12.0, 1.0, 0.625).unwrap(); // different x-knots
    assert!(Hull::new_asymmetric(clone_surface(a.surface()), clone_surface(b.surface())).is_err());
}

#[test]
fn negative_side_half_breadth_is_rejected() {
    let base = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let port = scaled(base.surface(), -0.1);
    let starboard = clone_surface(base.surface());

    assert!(
        Hull::new_asymmetric(port, starboard).is_err(),
        "each physical side must satisfy the non-negative half-breadth contract"
    );
}

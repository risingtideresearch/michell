//! Transom detection: does the half-breadth close at the aft end, and if not,
//! how much of a transom is it?

use michell::{hulls, BSplineSurface, Hull};

/// A wedge closing linearly toward the bow: `f = A (1 − x/L)(1 − z/T)`, so the
/// aft section is the full `f_T(z) = A(1 − z/T)` and every quantity the
/// detector reports has a closed form.
fn wedge(a: f64, length: f64, draft: f64) -> Hull {
    let knots_x = vec![0.0, 0.0, length, length];
    let knots_z = vec![0.0, 0.0, draft, draft];
    // Row-major, z fastest: (x0,z0) (x0,T) (x1,z0) (x1,T).
    let control = vec![a, 0.0, 0.0, 0.0];
    Hull::new(BSplineSurface::new(1, 1, knots_x, knots_z, control).unwrap()).unwrap()
}

#[test]
fn wigley_closes_at_both_ends() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    assert!(
        hull.transom().is_none(),
        "the Wigley half-breadth vanishes at x = ±L/2; it has no transom"
    );
}

#[test]
fn wigley_max_section_area_matches_the_closed_form() {
    let (l, b, t) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(l, b, t).unwrap();
    // A(x) = 2 ∫₀^T (B/2)(1 − (2x/L)²)(1 − (z/T)²) dz, maximal at x = 0:
    // A_X = B · (2T/3).
    let expect = 2.0 * b * t / 3.0;
    let got = hull.max_section_area();
    assert!(
        (got - expect).abs() < 1e-9 * expect,
        "A_X = {got}, expected {expect}"
    );
}

#[test]
fn wedge_transom_matches_the_closed_form() {
    let (a, l, t) = (0.4, 8.0, 0.25);
    let hull = wedge(a, l, t);
    let tr = hull.transom().expect("the wedge does not close aft");

    assert!((tr.x - 0.0).abs() < 1e-12, "transom station {}", tr.x);
    assert!(
        (tr.half_beam - a).abs() < 1e-12,
        "f_T(0) = {}, expected {a}",
        tr.half_beam
    );
    // A_T = 2 ∫₀^T A(1 − z/T) dz = A·T.
    let area = a * t;
    assert!(
        (tr.area - area).abs() < 1e-12 * area,
        "A_T = {}, expected {area}",
        tr.area
    );
    // A(x) = A·T·(1 − x/L) is maximal at the transom itself, so A_T/A_X = 1.
    assert!(
        (tr.area / hull.max_section_area() - 1.0).abs() < 1e-9,
        "A_T/A_X = {}",
        tr.area / hull.max_section_area()
    );
    // Equivalent rectangle of a linearly tapering transom: A_T/(2·A) = T/2.
    let depth = t / 2.0;
    assert!(
        (tr.depth - depth).abs() < 1e-12 * t,
        "depth = {}, expected {depth}",
        tr.depth
    );
}

#[test]
fn a_numerically_closed_stern_is_not_a_transom() {
    // A hull with a full midsection and a stern half-beam 1e-5 of it: the loft
    // wiggle case, which must not read as a transom.
    let (l, t) = (8.0, 0.25);
    let knots_x = vec![0.0, 0.0, 0.0, l, l, l];
    let knots_z = vec![0.0, 0.0, t, t];
    let control = vec![1e-5, 0.0, 1.0, 0.0, 0.0, 0.0];
    let hull = Hull::new(BSplineSurface::new(2, 1, knots_x, knots_z, control).unwrap()).unwrap();
    let tr = hull.transom();
    assert!(
        tr.is_none(),
        "stern area ratio {:?} should read as closed",
        tr.map(|t| t.area / hull.max_section_area())
    );
}

// ---------------------------------------------------------------------------
// The closure
// ---------------------------------------------------------------------------

use michell::{
    multihull_wave_resistance_with, wave_resistance_with, Conditions, Placement,
    TransomClosure, WaveOptions,
};

fn opts(transom: TransomClosure) -> WaveOptions {
    WaveOptions {
        transom,
        ..Default::default()
    }
}

#[test]
fn closure_is_inert_on_a_hull_that_closes_aft() {
    // Wigley has no transom, so every closure setting must give bit-for-bit
    // the same resistance — the guarantee that this feature cannot perturb
    // classical Michell.
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let cond = Conditions::seawater(3.0);
    let base = wave_resistance_with(&hull, &cond, &opts(TransomClosure::None))
        .unwrap()
        .resistance;
    for t in [
        TransomClosure::default(),
        TransomClosure::Ballistic { coeff: 4.0 },
        TransomClosure::Fixed { length: 2.0 },
    ] {
        let got = wave_resistance_with(&hull, &cond, &opts(t)).unwrap().resistance;
        assert_eq!(base.to_bits(), got.to_bits(), "{t:?} perturbed a closed hull");
    }
}

#[test]
fn zero_hollow_reproduces_the_analytic_step_term() {
    // As L_v -> 0 the virtual appendage collapses to a delta sheet at the
    // transom, whose amplitude is the closed form
    //     F_step = e^{i nu lambda x_T} * integral f_T(z) e^{-nu lambda^2 z} dz.
    // Here f_T(z) = A(1 - z/T), so the z-integral is elementary.
    let (a, l, t) = (0.4, 8.0, 0.25);
    let hull = wedge(a, l, t);
    let cond = Conditions::seawater(3.0);
    let nu = cond.gravity / (cond.speed * cond.speed);

    for lambda in [1.0, 1.7, 4.0] {
        let kappa = nu * lambda * lambda;
        // integral_0^T A(1 - z/T) e^{-kz} dz = A[(1-e^{-kT})/k - (1 - (1+kT)e^{-kT})/(k^2 T)]
        let e = (-kappa * t).exp();
        let zint =
            a * ((1.0 - e) / kappa - (1.0 - (1.0 + kappa * t) * e) / (kappa * kappa * t));
        // Phase is measured from the hull's x-centre, as the kernel does.
        let phase = nu * lambda * (0.0 - l / 2.0);
        let (want_re, want_im) = (zint * phase.cos(), zint * phase.sin());

        let closed = inner_integrals_with(&hull, &cond, lambda, TransomClosure::Fixed { length: 0.0 });
        let open = inner_integrals_with(&hull, &cond, lambda, TransomClosure::None);
        let (dre, dim) = (closed.0 - open.0, closed.1 - open.1);
        let scale = zint.abs().max(1e-12);
        assert!(
            (dre - want_re).abs() < 1e-9 * scale && (dim - want_im).abs() < 1e-9 * scale,
            "lambda {lambda}: closure added ({dre}, {dim}), analytic step ({want_re}, {want_im})"
        );
    }
}

/// `inner_integrals` with an explicit closure: one-hull fleet amplitude.
fn inner_integrals_with(
    hull: &michell::Hull,
    cond: &Conditions,
    lambda: f64,
    transom: TransomClosure,
) -> (f64, f64) {
    michell::inner_integrals_with(hull, cond, lambda, &opts(transom)).unwrap()
}

#[test]
fn a_longer_hollow_radiates_less() {
    // The physical content of the closure: stretching the hollow spreads the
    // same net change in half-beam over more length, so it makes less wave.
    // Monotone in the hollow length, and bracketed by the two limits.
    let hull = wedge(0.4, 8.0, 0.25);
    let cond = Conditions::seawater(3.0);
    let r = |t| {
        wave_resistance_with(&hull, &cond, &opts(t))
            .unwrap()
            .resistance
    };
    let step = r(TransomClosure::Fixed { length: 0.0 });
    let mut prev = step;
    for len in [0.5, 1.0, 2.0, 4.0] {
        let now = r(TransomClosure::Fixed { length: len });
        assert!(
            now < prev,
            "hollow {len} m gave {now} N, not less than {prev} N"
        );
        prev = now;
    }
    let open = r(TransomClosure::None);
    assert!(
        step > open,
        "an abrupt transom ({step} N) should make more wave than ignoring it ({open} N)"
    );
}

#[test]
fn ballistic_hollow_grows_with_speed() {
    // L_v = c*U*sqrt(d_T/g): the faster the boat, the longer the hollow, so
    // the transom's share of the wave-making falls away with speed.
    let hull = wedge(0.4, 8.0, 0.25);
    let tr = hull.transom().unwrap();
    let g = 9.80665;
    for u in [2.0, 4.0, 8.0] {
        let cond = Conditions::seawater(u);
        let want = std::f64::consts::SQRT_2 * u * (tr.depth / g).sqrt();
        // Match the ballistic default against an explicit Fixed hollow of the
        // length the formula predicts.
        let a = wave_resistance_with(&hull, &cond, &opts(TransomClosure::default()))
            .unwrap()
            .resistance;
        let b = wave_resistance_with(&hull, &cond, &opts(TransomClosure::Fixed { length: want }))
            .unwrap()
            .resistance;
        assert!(
            (a - b).abs() <= 1e-9 * a.abs().max(1.0),
            "U = {u}: ballistic {a} N vs explicit L_v = {want} m giving {b} N"
        );
    }
}

#[test]
fn closure_survives_placement_in_a_fleet() {
    // The transom term carries its own phase, so a fleet of two transom hulls
    // must still be translation-covariant: shifting the whole fleet in x
    // changes no resistance.
    let hull = wedge(0.4, 8.0, 0.25);
    let cond = Conditions::seawater(3.0);
    let o = opts(TransomClosure::default());
    let at = |dx: f64| {
        let m = [
            (&hull, Placement { x: dx, y: 1.6 }),
            (&hull, Placement { x: dx, y: -1.6 }),
        ];
        multihull_wave_resistance_with(&m, &cond, &o)
            .unwrap()
            .resistance
    };
    let (a, b) = (at(0.0), at(25.0));
    assert!(
        (a - b).abs() < 1e-9 * a,
        "fleet resistance moved with a rigid x-shift: {a} vs {b}"
    );
}

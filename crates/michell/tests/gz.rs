//! Righting-arm (GZ) validation: exact transverse hull properties on the
//! analytic Wigley hull, and rectangular barges pushed through the heel-aware
//! inclined-waterplane equilibrium.

use michell::body::{Body, BodyOptions};
use michell::float::{
    fleet_cg, solve_equilibrium_bodies, solve_equilibrium_heeled, HullLoad, LoadCase, PointLoad,
};
use michell::iges::HullPose;
use michell::inclined::InclinedGrid;
use michell::BSplineSurface;

#[test]
fn wigley_transverse_properties_are_exact() {
    let (l, b, t) = (10.0, 1.0, 0.625);
    let hull = michell::hulls::wigley(l, b, t).unwrap();
    // Waterline beam 2f = B·4ξ(1−ξ): I_T = ∫ (2/3) f³ dx = (4/105) B³ L.
    let want_it = 4.0 / 105.0 * b * b * b * l;
    let it = hull.waterplane_transverse_moment();
    assert!(
        (it - want_it).abs() < 1e-12 * want_it,
        "I_T {it} vs {want_it}"
    );
    // Parabolic sections 1 − (z/T)²: KB = 3T/8 below the waterline.
    let want_kb = 3.0 * t / 8.0;
    let kb = hull.vcb_z();
    assert!(
        (kb - want_kb).abs() < 1e-12 * want_kb,
        "KB {kb} vs {want_kb}"
    );
}

/// A rectangular barge as a full-band body: constant half-beam `b` over
/// length `l` and band depth `d`, design waterline `w` below the band top.
fn barge(l: f64, b: f64, d: f64, w: f64, centerplane: f64) -> Body {
    let surface =
        BSplineSurface::new(1, 1, vec![0.0, 0.0, l, l], vec![0.0, 0.0, d, d], vec![b; 4]).unwrap();
    Body::new(surface, w, centerplane).unwrap()
}

/// The fleet CG is the mass-weighted sum of the per-hull loads, carried through
/// each hull's pose: a symmetric pair cancels transversely; unequal masses
/// bias the CG; and dx/dy/dz translate a hull's CG one-for-one.
#[test]
fn fleet_cg_is_mass_weighted_and_tracks_pose() {
    let l = 10.0;
    let port = barge(l, 0.5, 1.2, 0.8, -1.5);
    let stbd = barge(l, 0.5, 1.2, 0.8, 1.5);
    let bodies = [&port, &stbd];
    let load = |mass| HullLoad {
        mass,
        lcg: 5.0,
        vcg: 0.3,
        points: vec![],
    };

    // Symmetric equal masses: transverse cancels, vertical/longitudinal shared.
    let cg = fleet_cg(
        &bodies,
        &[load(1000.0), load(1000.0)],
        &[HullPose::default(); 2],
    );
    assert!((cg.mass - 2000.0).abs() < 1e-9);
    assert!(cg.tcg.abs() < 1e-9, "tcg {}", cg.tcg);
    assert!((cg.lcg - 5.0).abs() < 1e-9 && (cg.vcg - 0.3).abs() < 1e-9);

    // Heavier port hull pulls the CG to port (−y): (3000·−1.5 + 1000·1.5)/4000.
    let cg = fleet_cg(
        &bodies,
        &[load(3000.0), load(1000.0)],
        &[HullPose::default(); 2],
    );
    assert!((cg.tcg - (-0.75)).abs() < 1e-9, "tcg {}", cg.tcg);

    // Pose translation: +dx raises lcg, +dy shifts tcg, +dz lowers vcg.
    let pose = HullPose {
        dx: 2.0,
        dy: 0.5,
        dz: 0.4,
        ..HullPose::default()
    };
    let center = barge(l, 0.5, 1.2, 0.8, 0.0);
    let cg = fleet_cg(&[&center], &[load(1000.0)], &[pose]);
    assert!((cg.lcg - 7.0).abs() < 1e-9, "lcg {}", cg.lcg);
    assert!((cg.tcg - 0.5).abs() < 1e-9, "tcg {}", cg.tcg);
    assert!((cg.vcg - (-0.1)).abs() < 1e-9, "vcg {}", cg.vcg);

    // Design trim rotates the CG about the hull midship (the pivot default):
    // a CG on the pivot station at height v maps to (5 − v·sinτ, v·cosτ).
    let tau = 0.2f64;
    let cg = fleet_cg(
        &[&center],
        &[load(1000.0)],
        &[HullPose {
            trim: tau,
            ..HullPose::default()
        }],
    );
    assert!(
        (cg.lcg - (5.0 - 0.3 * tau.sin())).abs() < 1e-9,
        "lcg {}",
        cg.lcg
    );
    assert!((cg.vcg - 0.3 * tau.cos()).abs() < 1e-9, "vcg {}", cg.vcg);

    // Massless hulls drop out entirely.
    let cg = fleet_cg(
        &bodies,
        &[load(0.0), load(1000.0)],
        &[HullPose::default(); 2],
    );
    assert!((cg.tcg - 1.5).abs() < 1e-9 && (cg.mass - 1000.0).abs() < 1e-9);

    // A point load adds a mass at an offset from the hull centerpoint (midship
    // x=5, centreplane y=0, floatplane): dz is +down, so it lowers the CG.
    let mut hl = load(1000.0); // structural: 1000 kg at (5, 0, +0.3)
    hl.points.push(PointLoad {
        mass: 1000.0,
        dx: 1.0, // → local x = 5 + 1 = 6
        dy: 0.4, // → transverse offset +0.4
        dz: 0.8, // +down → vcg = −0.8
    });
    let cg = fleet_cg(&[&center], &[hl.clone()], &[HullPose::default()]);
    assert!((cg.mass - 2000.0).abs() < 1e-9, "mass {}", cg.mass);
    assert!((cg.lcg - 5.5).abs() < 1e-9, "lcg {}", cg.lcg); // (5+6)/2
    assert!((cg.tcg - 0.2).abs() < 1e-9, "tcg {}", cg.tcg); // (0+0.4)/2
    assert!((cg.vcg - (-0.25)).abs() < 1e-9, "vcg {}", cg.vcg); // (0.3−0.8)/2

    // A point load rides with the hull's pose (dx shifts it too).
    let cg = fleet_cg(
        &[&center],
        &[hl],
        &[HullPose {
            dx: 2.0,
            ..HullPose::default()
        }],
    );
    assert!((cg.lcg - 7.5).abs() < 1e-9, "lcg {}", cg.lcg); // 5.5 + 2
}

/// At a small heel the inclined-cut righting arm reduces to the metacentric
/// `GZ = GM·sinφ` for a single centreline barge — the linear limit the exact
/// cut must recover (`GM = I_T/∇ − KB − vcg`, KB down / G up from the waterline;
/// design waterline == float waterline so the datums coincide).
#[test]
fn inclined_gz_matches_metacentric_gm_at_small_angle() {
    let (l, b, d) = (8.0f64, 0.5f64, 1.2f64);
    let t = 0.4f64;
    let rho = 1000.0;
    let mass = rho * 2.0 * b * l * t;
    let vcg = 0.1;
    let body = barge(l, b, d, d - t, 0.0);
    let bodies = [&body];
    let heel = 2.0f64.to_radians();
    let eq = solve_equilibrium_heeled(
        &bodies,
        0.0,
        &[HullPose::default()],
        &LoadCase { mass, lcg: None },
        rho,
        heel,
        vcg,
        0.0,
        &BodyOptions::default(),
        InclinedGrid::default(),
    )
    .unwrap();
    let gm = (2.0 / 3.0 * b * b * b * l) / (mass / rho) - t / 2.0 - vcg;
    let want = gm * heel.sin();
    assert!(
        (eq.gz - want).abs() < 3e-3 * want.abs(),
        "GZ {} vs GM·sinφ {want}",
        eq.gz
    );
}

// --- Heel-aware equilibrium with inclined-waterplane hydrostatics ---------

/// At zero heel the inclined-cut solver reproduces the horizontal-cut solver
/// (`solve_equilibrium_bodies`) to the loft/integration tolerance, and the
/// righting arm is zero.
#[test]
fn heeled_equilibrium_zero_heel_matches_upright() {
    let (l, b, d, w, s) = (8.0f64, 0.5f64, 1.2f64, 0.8f64, 1.5f64);
    let rho = 1000.0;
    let t_mean = 0.4f64;
    let mass = rho * 2.0 * (2.0 * b * l) * t_mean;
    let port = barge(l, b, d, w, -s);
    let stbd = barge(l, b, d, w, s);
    let bodies = [&port, &stbd];
    let poses = [HullPose::default(); 2];
    let load = LoadCase { mass, lcg: None };
    let opts = BodyOptions::default();
    let grid = InclinedGrid::default();

    let up = solve_equilibrium_bodies(&bodies, 0.0, &poses, &load, rho, &opts).unwrap();
    let he = solve_equilibrium_heeled(&bodies, 0.0, &poses, &load, rho, 0.0, 0.3, 0.0, &opts, grid)
        .unwrap();

    assert!(
        (he.volume - mass / rho).abs() < 1e-3 * mass / rho,
        "V {}",
        he.volume
    );
    assert!(
        (he.sinkage - up.sinkage).abs() < 2e-3,
        "sinkage {} vs {}",
        he.sinkage,
        up.sinkage
    );
    assert!(he.gz.abs() < 1e-6, "GZ(0) = {}", he.gz);
}

/// Heeled: the solver holds displacement, the righting arm is antisymmetric,
/// and (past the linear regime) it exceeds the metacentric straight line — the
/// inclined cut recovers the form stability the approximation drops.
#[test]
fn heeled_equilibrium_holds_displacement_and_beats_metacentric() {
    let (l, b, d, w, s) = (8.0f64, 0.5f64, 1.2f64, 0.8f64, 1.5f64);
    let rho = 1000.0;
    let t_mean = 0.4f64;
    let mass = rho * 2.0 * (2.0 * b * l) * t_mean;
    let vcg = 0.3;
    let i_t = 2.0 / 3.0 * b * b * b * l;
    let area = 2.0 * b * l;
    let port = barge(l, b, d, w, -s);
    let stbd = barge(l, b, d, w, s);
    let bodies = [&port, &stbd];
    let poses = [HullPose::default(); 2];
    let load = LoadCase { mass, lcg: None };
    let opts = BodyOptions::default();
    let grid = InclinedGrid::default();

    let gz_at = |deg: f64| {
        let he = solve_equilibrium_heeled(
            &bodies,
            0.0,
            &poses,
            &load,
            rho,
            (deg as f64).to_radians(),
            vcg,
            0.0,
            &opts,
            grid,
        )
        .unwrap();
        assert!(
            (he.volume - mass / rho).abs() < 1e-3 * mass / rho,
            "{deg}°: V {} vs {}",
            he.volume,
            mass / rho
        );
        he.gz
    };
    // Metacentric closed form (linear form stability).
    let metacentric = |deg: f64| {
        let (sin, cos) = (deg as f64).to_radians().sin_cos();
        let delta = s * sin;
        s * s * sin * cos / t_mean
            + sin * (i_t / (area * t_mean) - (t_mean * t_mean + delta * delta) / (2.0 * t_mean))
            - vcg * sin
    };

    assert!(gz_at(0.0).abs() < 1e-6);
    // Antisymmetric.
    assert!((gz_at(6.0) + gz_at(-6.0)).abs() < 1e-3 * gz_at(6.0).abs());
    // Beats the metacentric line by a margin that grows with heel.
    let e5 = gz_at(5.0) - metacentric(5.0);
    let e8 = gz_at(8.0) - metacentric(8.0);
    assert!(
        e5 > 0.0 && e8 > e5,
        "form-stability excess should grow: e5={e5} e8={e8}"
    );
}

/// The lcg path exercises the 2-D finite-difference Jacobian (sinkage + pitch).
/// A midship LCG must solve to ~zero trim; an LCG shifted forward must trim the
/// hull to bring the LCB there — both holding displacement.
#[test]
fn heeled_equilibrium_solves_trim() {
    // A single centred barge with topside band (so the design attitude has a
    // live waterplane) whose x-domain runs 0..l.
    let (l, b, d, w) = (10.0f64, 0.5f64, 1.2f64, 0.8f64);
    let rho = 1025.0;
    let t_mean = 0.4f64;
    let mass = rho * (2.0 * b * l) * t_mean;
    let body = barge(l, b, d, w, 0.0);
    let bodies = [&body];
    let poses = [HullPose::default()];
    let opts = BodyOptions::default();
    let grid = InclinedGrid::default();

    let solve = |lcg: f64| {
        solve_equilibrium_heeled(
            &bodies,
            0.0,
            &poses,
            &LoadCase {
                mass,
                lcg: Some(lcg),
            },
            rho,
            0.0,
            0.1,
            0.0,
            &opts,
            grid,
        )
        .unwrap()
    };

    // Midship LCG → ~zero trim.
    let mid = solve(l / 2.0);
    assert!(
        mid.volume_residual < 5e-3,
        "volume residual {}",
        mid.volume_residual
    );
    assert!(mid.lcb_residual < 1e-2, "lcb residual {}", mid.lcb_residual);
    assert!(
        mid.trim.abs() < 2e-3,
        "midship LCG should give ~zero trim, got {}",
        mid.trim
    );

    // LCG shifted forward → a definite trim that lands the LCB on the LCG.
    let fwd = solve(l / 2.0 + 0.5);
    assert!(
        fwd.volume_residual < 5e-3,
        "volume residual {}",
        fwd.volume_residual
    );
    assert!(fwd.lcb_residual < 1e-2, "lcb residual {}", fwd.lcb_residual);
    assert!(
        fwd.trim.abs() > 1e-3,
        "offset LCG should trim the hull, got {}",
        fwd.trim
    );
}

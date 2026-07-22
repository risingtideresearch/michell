//! Righting-arm (GZ) validation: exact transverse hull properties on the
//! analytic Wigley hull, and a rectangular-barge catamaran pushed through the
//! heel + equilibrium path and compared against the closed form.

use michell::body::{Body, BodyOptions};
use michell::float::{
    heel_poses, righting_arm, solve_equilibrium_bodies, solve_equilibrium_heeled, LoadCase,
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
    assert!((it - want_it).abs() < 1e-12 * want_it, "I_T {it} vs {want_it}");
    // Parabolic sections 1 − (z/T)²: KB = 3T/8 below the waterline.
    let want_kb = 3.0 * t / 8.0;
    let kb = hull.vcb_z();
    assert!((kb - want_kb).abs() < 1e-12 * want_kb, "KB {kb} vs {want_kb}");
}

/// A rectangular barge as a full-band body: constant half-beam `b` over
/// length `l` and band depth `d`, design waterline `w` below the band top.
fn barge(l: f64, b: f64, d: f64, w: f64, centerplane: f64) -> Body {
    let surface = BSplineSurface::new(
        1,
        1,
        vec![0.0, 0.0, l, l],
        vec![0.0, 0.0, d, d],
        vec![b; 4],
    )
    .unwrap();
    Body::new(surface, w, centerplane).unwrap()
}

#[test]
fn barge_catamaran_gz_matches_closed_form() {
    // Two barges (L = 8, beam 1) at y = ±1.5, floating at mean draft 0.4.
    let (l, b, d, w, s) = (8.0f64, 0.5f64, 1.2f64, 0.6f64, 1.5f64);
    let rho = 1000.0;
    let t_mean = 0.4f64;
    let area = 2.0 * b * l; // waterplane per barge
    let mass = rho * 2.0 * area * t_mean;
    let i_t = 2.0 / 3.0 * b * b * b * l; // per-barge transverse inertia
    let vcg = 0.3;

    let port = barge(l, b, d, w, -s);
    let stbd = barge(l, b, d, w, s);
    let bodies = [&port, &stbd];
    let opts = BodyOptions::default();

    let gz_at = |heel_deg: f64| -> f64 {
        let heel = heel_deg.to_radians();
        let poses = heel_poses(&bodies, &[HullPose::default(); 2], heel).unwrap();
        let eq = solve_equilibrium_bodies(
            &bodies,
            0.0,
            &poses,
            &LoadCase { mass, lcg: None },
            rho,
            &opts,
        )
        .unwrap();
        assert_eq!(eq.fleet.dry, 0, "a barge flew at {heel_deg} deg");
        assert!(
            (eq.volume - mass / rho).abs() < 1e-3 * (mass / rho),
            "volume {} vs {}",
            eq.volume,
            mass / rho
        );
        righting_arm(&eq.fleet, heel, vcg)
    };

    // Closed form for wall-sided barges with equal waterplanes (mean draft is
    // preserved under heel): drafts T ± δ with δ = s·sinφ,
    //   GZ = s² sinφ cosφ / T                      (transfer between hulls)
    //      + sinφ·(I_T/(A·T) − (T² + δ²)/(2T))     (per-hull BM − KB)
    //      − vcg·sinφ.
    for heel_deg in [2.0f64, 5.0, 8.0] {
        let (sin, cos) = heel_deg.to_radians().sin_cos();
        let delta = s * sin;
        let want = s * s * sin * cos / t_mean
            + sin * (i_t / (area * t_mean) - (t_mean * t_mean + delta * delta) / (2.0 * t_mean))
            - vcg * sin;
        let gz = gz_at(heel_deg);
        assert!(
            (gz - want).abs() < 2e-3 * want.abs(),
            "GZ({heel_deg} deg) = {gz} vs closed form {want}"
        );
    }

    // Antisymmetry and zero at upright.
    let up = gz_at(0.0);
    assert!(up.abs() < 1e-9, "GZ(0) = {up}");
    let (plus, minus) = (gz_at(5.0), gz_at(-5.0));
    assert!(
        (plus + minus).abs() < 1e-6 * plus.abs(),
        "GZ not antisymmetric: {plus} vs {minus}"
    );
}

#[test]
fn righting_arm_reduces_to_metacentric_gm_for_a_monohull() {
    // A single centerline barge heeled a small angle must give
    // GZ = GM·sinφ with GM = I_T/∇ − KB_depth − vcg (KB measured down, G
    // measured up from the waterplane).
    let (l, b, d, w) = (8.0f64, 0.5f64, 1.2f64, 0.6f64);
    let rho = 1000.0;
    let t = 0.5f64;
    let mass = rho * 2.0 * b * l * t;
    let vcg = 0.1;
    let body = barge(l, b, d, w, 0.0);
    let bodies = [&body];
    let heel = 2.0f64.to_radians();
    let poses = heel_poses(&bodies, &[HullPose::default()], heel).unwrap();
    let eq = solve_equilibrium_bodies(
        &bodies,
        0.0,
        &poses,
        &LoadCase { mass, lcg: None },
        rho,
        &BodyOptions::default(),
    )
    .unwrap();
    let gz = righting_arm(&eq.fleet, heel, vcg);
    let i_t = 2.0 / 3.0 * b * b * b * l;
    let gm = i_t / (mass / rho) - t / 2.0 - vcg;
    let want = gm * heel.sin();
    assert!(
        (gz - want).abs() < 2e-3 * want.abs(),
        "GZ {gz} vs GM·sinφ {want}"
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
    let he = solve_equilibrium_heeled(&bodies, 0.0, &poses, &load, rho, 0.0, 0.3, &opts, grid)
        .unwrap();

    assert!((he.volume - mass / rho).abs() < 1e-3 * mass / rho, "V {}", he.volume);
    assert!((he.sinkage - up.sinkage).abs() < 2e-3, "sinkage {} vs {}", he.sinkage, up.sinkage);
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
            &bodies, 0.0, &poses, &load, rho, (deg as f64).to_radians(), vcg, &opts, grid,
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
    assert!(e5 > 0.0 && e8 > e5, "form-stability excess should grow: e5={e5} e8={e8}");
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
            &LoadCase { mass, lcg: Some(lcg) },
            rho,
            0.0,
            0.1,
            &opts,
            grid,
        )
        .unwrap()
    };

    // Midship LCG → ~zero trim.
    let mid = solve(l / 2.0);
    assert!(mid.volume_residual < 5e-3, "volume residual {}", mid.volume_residual);
    assert!(mid.lcb_residual < 1e-2, "lcb residual {}", mid.lcb_residual);
    assert!(mid.trim.abs() < 2e-3, "midship LCG should give ~zero trim, got {}", mid.trim);

    // LCG shifted forward → a definite trim that lands the LCB on the LCG.
    let fwd = solve(l / 2.0 + 0.5);
    assert!(fwd.volume_residual < 5e-3, "volume residual {}", fwd.volume_residual);
    assert!(fwd.lcb_residual < 1e-2, "lcb residual {}", fwd.lcb_residual);
    assert!(fwd.trim.abs() > 1e-3, "offset LCG should trim the hull, got {}", fwd.trim);
}

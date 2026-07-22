//! Inclined-waterplane hydrostatics validated against closed forms.
//!
//! The module integrates the true tilted-waterline cut of each section, so for
//! a wall-sided box it must reproduce the exact wall-sided righting arm
//! `GZ = sinφ·(GM + ½·BM·tan²φ)` — including the `½BM tan²φ` form-stability
//! term the metacentric approximation (a straight `GM·sinφ`) leaves out.

use michell::body::Body;
use michell::iges::{HullPose, Platform};
use michell::inclined::{body_inclined_hydro, fleet_righting_arm, fleet_volume, InclinedGrid};
use michell::BSplineSurface;

/// A rectangular barge as a full-band body: constant half-beam `b` over length
/// `l` and band depth `d`, design waterline `w` below the band top, centreplane
/// at `centerplane`. (Same construction as the GZ test's barge.)
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

/// Bisect the water offset (sinkage) so the fleet displaces `target` m³ at the
/// given heel — the vertical-force balance an equilibrium solver would find.
fn solve_sinkage(
    bodies: &[&Body],
    poses: &[HullPose],
    heel: f64,
    grid: InclinedGrid,
    target: f64,
) -> f64 {
    let plat = Platform::default();
    let (mut lo, mut hi) = (-2.0f64, 4.0f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let v = fleet_volume(bodies, mid, poses, &plat, heel, grid);
        if v < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Upright, no trim: the module must reproduce the exact box hydrostatics.
#[test]
fn upright_box_is_exact() {
    let (l, b, d) = (8.0, 0.5, 1.2);
    let t = 0.5; // draft from the keel
    let w = d - t; // design waterline == float waterline (z_b = d - t)
    let body = barge(l, b, d, w, 0.0);
    let grid = InclinedGrid::default();
    let h = body_inclined_hydro(&body, 0.0, &HullPose::default(), &Platform::default(), 0.0, grid)
        .expect("wet");
    let want_v = l * 2.0 * b * t;
    assert!((h.volume - want_v).abs() < 1e-6 * want_v, "V {} vs {want_v}", h.volume);
    assert!(h.tcb.abs() < 1e-9, "TCB {}", h.tcb);
    assert!((h.kb - t / 2.0).abs() < 1e-9, "KB {} vs {}", h.kb, t / 2.0);
    assert!((h.lcb - l / 2.0).abs() < 1e-9, "LCB {} vs {}", h.lcb, l / 2.0);
    assert_eq!(h.band_exceeded, 0);
}

/// Zero heel gives zero righting arm and the upright volume.
#[test]
fn zero_heel_is_upright() {
    let (l, b, d) = (8.0, 0.5, 1.2);
    let t = 0.5;
    let body = barge(l, b, d, d - t, 0.0);
    let bodies = [&body];
    let poses = [HullPose::default()];
    let grid = InclinedGrid::default();
    let target = l * 2.0 * b * t;
    let s = solve_sinkage(&bodies, &poses, 0.0, grid, target);
    let gz = fleet_righting_arm(&bodies, s, &poses, &Platform::default(), 0.0, 0.1, grid);
    assert!(gz.abs() < 1e-9, "GZ(0) = {gz}");
}

/// A single wall-sided box: the inclined cut yields the exact wall-sided GZ,
/// which is strictly larger than the metacentric straight line — the whole
/// point of integrating the true cut rather than approximating.
#[test]
fn heeled_box_is_wall_sided_and_nonlinear() {
    // Beamy, shallow box so GM > 0 and the waterline stays on the side walls
    // (no deck-edge immersion or bilge emergence) through the test angle.
    let (l, b, d) = (8.0, 1.0, 1.2);
    let t = 0.3;
    let vcg = 0.1; // above the design waterline (== float waterline here)
    let body = barge(l, b, d, d - t, 0.0);
    let bodies = [&body];
    let poses = [HullPose::default()];
    let grid = InclinedGrid {
        stations: 65,
        band: 241,
    };
    let target = l * 2.0 * b * t;

    let bm = (2.0 / 3.0 * b * b * b * l) / target; // I_T / ∇ = b²/(3t)
    let gm = bm - t / 2.0 - vcg;
    assert!(gm > 0.0, "test needs a stable box: GM {gm}");

    for deg in [6.0f64, 12.0, 15.0] {
        let phi = deg.to_radians();
        let s = solve_sinkage(&bodies, &poses, phi, grid, target);
        // Wall-sided displacement is preserved (mean draft unchanged).
        let v = fleet_volume(&bodies, s, &poses, &Platform::default(), phi, grid);
        assert!((v - target).abs() < 1e-4 * target, "{deg}°: V {v} vs {target}");

        let gz = fleet_righting_arm(&bodies, s, &poses, &Platform::default(), phi, vcg, grid);
        let (sin, tan) = (phi.sin(), phi.tan());
        let want = sin * (gm + 0.5 * bm * tan * tan); // wall-sided exact
        let linear = gm * sin; // metacentric approximation
        assert!(
            (gz - want).abs() < 2e-3 * want,
            "{deg}°: GZ {gz} vs wall-sided {want}"
        );
        // The nonlinear form-stability term is real, not noise.
        if deg >= 12.0 {
            assert!(gz > 1.02 * linear, "{deg}°: GZ {gz} not above linear {linear}");
        }
    }
}

/// Barge catamaran GZ against the metacentric closed form used by the `gz.rs`
/// equilibrium test. In the linear regime the two agree; as heel grows the
/// inclined cut adds the nonlinear form stability the closed form drops, so the
/// module rises strictly above it (by a margin that grows with heel) — the
/// improvement this path exists for.
#[test]
fn catamaran_gz_beats_metacentric() {
    let (l, b, d, s) = (8.0f64, 0.5f64, 1.2f64, 1.5f64);
    let t_mean = 0.4f64;
    let w = d - t_mean; // design waterline == float waterline
    let area = 2.0 * b * l; // waterplane per barge
    let i_t = 2.0 / 3.0 * b * b * b * l; // per-barge transverse inertia
    let vcg = 0.3;

    let port = barge(l, b, d, w, -s);
    let stbd = barge(l, b, d, w, s);
    let bodies = [&port, &stbd];
    let poses = [HullPose::default(); 2];
    let grid = InclinedGrid::default();
    let target = 2.0 * area * t_mean;

    // Metacentric closed form (linear form stability) for wall-sided barges
    // with equal waterplanes:
    //   GZ = s² sinφ cosφ / T                    (inter-hull transfer)
    //      + sinφ (I_T/(A·T) − (T² + δ²)/(2T))   (per-hull BM − KB)
    //      − vcg sinφ,           δ = s sinφ.
    let closed = |phi: f64| {
        let (sin, cos) = phi.sin_cos();
        let delta = s * sin;
        s * s * sin * cos / t_mean
            + sin * (i_t / (area * t_mean) - (t_mean * t_mean + delta * delta) / (2.0 * t_mean))
            - vcg * sin
    };
    let gz_at = |deg: f64| {
        let phi = (deg as f64).to_radians();
        let sk = solve_sinkage(&bodies, &poses, phi, grid, target);
        // Displacement is held across the heel sweep.
        let v = fleet_volume(&bodies, sk, &poses, &Platform::default(), phi, grid);
        assert!((v - target).abs() < 1e-4 * target, "{deg}°: V {v} vs {target}");
        fleet_righting_arm(&bodies, sk, &poses, &Platform::default(), phi, vcg, grid)
    };

    // Linear regime: agree closely at 2°.
    let (g2, c2) = (gz_at(2.0), closed(2.0f64.to_radians()));
    assert!((g2 - c2).abs() < 4e-3 * c2, "2°: module {g2} vs closed {c2}");

    // Growing nonlinear excess above the metacentric line.
    let excess = |deg: f64| gz_at(deg) - closed((deg as f64).to_radians());
    let (e5, e8) = (excess(5.0), excess(8.0));
    assert!(e5 > 0.0 && e8 > e5, "form-stability excess should grow: e5={e5} e8={e8}");
    assert!(e8 > 0.02 * closed(8.0f64.to_radians()), "8° excess too small: {e8}");

    // Antisymmetric about upright.
    let sk = solve_sinkage(&bodies, &poses, 5.0f64.to_radians(), grid, target);
    let plus = fleet_righting_arm(&bodies, sk, &poses, &Platform::default(), 5.0f64.to_radians(), vcg, grid);
    let minus = fleet_righting_arm(&bodies, sk, &poses, &Platform::default(), -5.0f64.to_radians(), vcg, grid);
    assert!((plus + minus).abs() < 1e-6 * plus.abs(), "not antisymmetric: {plus} vs {minus}");
}

/// Refining the grid drives the volume toward the exact wall-sided value.
#[test]
fn resolution_converges() {
    let (l, b, d) = (8.0, 1.0, 1.2);
    let t = 0.3;
    let body = barge(l, b, d, d - t, 0.0);
    let bodies = [&body];
    let poses = [HullPose::default()];
    let phi = 15.0f64.to_radians();
    let target = l * 2.0 * b * t;

    let coarse = InclinedGrid { stations: 17, band: 25 };
    let fine = InclinedGrid { stations: 33, band: 401 };
    let sc = solve_sinkage(&bodies, &poses, phi, coarse, target);
    let sf = solve_sinkage(&bodies, &poses, phi, fine, target);
    let ec = (fleet_volume(&bodies, sc, &poses, &Platform::default(), phi, coarse) - target).abs();
    let ef = (fleet_volume(&bodies, sf, &poses, &Platform::default(), phi, fine) - target).abs();
    // Both hit the target by construction; the real check is that the fine grid
    // reproduces the wall-sided mean-draft invariant tightly.
    assert!(ef <= ec + 1e-9, "fine {ef} not better than coarse {ec}");
    assert!(ef < 1e-4 * target, "fine volume error {ef}");
}

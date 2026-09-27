//! Dynamic sinkage and trim: the hydrostatic equilibrium solver with a
//! speed-dependent vertical force and pitch moment in the balance, checked
//! against the hydrostatic solver and against closed-form shifts on a
//! fore-aft symmetric Wigley hull with freeboard.

use michell_geometry::float::{
    DynamicModel,
    solve_equilibrium_dynamic_with, solve_equilibrium_sectional,
    solve_equilibrium_sectional_dynamic, DynamicEquilibrium, DynamicLoad, FleetState, LoadCase,
};
use michell_geometry::iges::{self, HullPose, NurbsSurface3, Platform, SectionalOptions, SourceFleet};
use michell_geometry::SectionalHull;
use michell_geometry::source::SourceHull;
use michell::{Result, STANDARD_GRAVITY};

const RHO: f64 = 1000.0;
const G: f64 = STANDARD_GRAVITY;

/// The Wigley hull `(B/2)(1 − (2x/L)²)(1 − (z/T)²)` with freeboard: the
/// parabolic sections below the design waterline (CAD z = 0), carried
/// wall-sided up through a freeboard `fb` so the solver has topside to sink
/// into. Exactly a two-span biquadratic B-spline graph (control points at
/// the Greville abscissae); the C⁰ knot at the waterline is where the section
/// is flat anyway.
fn wigley_source(l: f64, b: f64, t: f64, fb: f64) -> SourceFleet {
    let a = l / 2.0;
    let knots_u = vec![-a, -a, -a, a, a, a];
    // v runs down from the freeboard top, as depth below it.
    let knots_v = vec![0.0, 0.0, 0.0, fb, fb, fb + t, fb + t, fb + t];
    let xs = [-a, 0.0, a];
    let depths = [0.0, fb / 2.0, fb, fb + t / 2.0, fb + t];
    // Bezier control values: 4t(1−t) in x; 1 over the freeboard, then 1−t².
    let gx = [0.0, 2.0, 0.0];
    let hz = [1.0, 1.0, 1.0, 1.0, 0.0];
    let side = |sgn: f64| {
        let mut ctrl = Vec::with_capacity(15);
        for i in 0..3 {
            for j in 0..5 {
                ctrl.push([xs[i], sgn * b / 2.0 * gx[i] * hz[j], fb - depths[j]]);
            }
        }
        NurbsSurface3 {
            degree_u: 2,
            degree_v: 2,
            knots_u: knots_u.clone(),
            knots_v: knots_v.clone(),
            n_ctrl_u: 3,
            n_ctrl_v: 5,
            ctrl,
            weights: vec![1.0; 15],
            trim_uv: None,
        }
    };
    iges::source_fleet_from_surfaces(vec![side(1.0), side(-1.0)], 1.0, 0.0).unwrap()
}

fn hulls(src: &SourceFleet) -> [SourceHull<'_>; 1] {
    [SourceHull {
        source: src,
        index: 0,
        waterline_z: 0.0,
        pose: HullPose::default(),
    }]
}

/// The standard-proportion Wigley with freeboard and the mass it floats at
/// its design waterline: `∇ = (4/9) B L T`.
fn setup() -> (SourceFleet, f64) {
    let (l, b, t, fb) = (10.0, 1.0, 0.625, 0.3);
    (wigley_source(l, b, t, fb), RHO * 4.0 / 9.0 * b * l * t)
}

fn solve_dynamic(
    src: &SourceFleet,
    load: &LoadCase,
    dynamic: impl DynamicModel<SectionalHull>,
    warm_start: Option<(f64, f64)>,
) -> DynamicEquilibrium {
    solve_equilibrium_sectional_dynamic(
        &hulls(src),
        load,
        RHO,
        G,
        &SectionalOptions::default(),
        dynamic,
        warm_start,
    )
    .unwrap()
}

/// Hydrostatic near-miss acceptance: 10× the fine-phase volume tolerance.
const V_TOL: f64 = 2e-3;

/// With a zero dynamic load the dynamic solver is the hydrostatic solver:
/// the same iterates and final attitude, bit for bit, both for a
/// sinkage-only balance and for one with trim. The final cut's hydrostatics
/// agree to rounding only: the dynamic solve's finite-difference probes
/// re-cut the hull at nearby attitudes too, and each cut is warm-started
/// from the last.
#[test]
fn zero_dynamic_load_reproduces_hydrostatic_solve_exactly() {
    let (src, mass) = setup();
    let opts = SectionalOptions::default();
    for lcg in [None, Some(0.2)] {
        let load = LoadCase { mass, lcg };
        let eq = solve_equilibrium_sectional(&hulls(&src), &load, RHO, &opts).unwrap();
        let dy = solve_dynamic(&src, &load, |_: &FleetState| Ok(DynamicLoad::default()), None);
        assert_eq!(
            dy.sinkage.to_bits(),
            eq.sinkage.to_bits(),
            "sinkage {lcg:?}"
        );
        assert_eq!(dy.trim.to_bits(), eq.trim.to_bits(), "trim {lcg:?}");
        let close = |a: f64, b: f64, scale: f64| (a - b).abs() <= 1e-9 * scale;
        assert!(close(dy.volume, eq.volume, eq.volume), "volume {lcg:?}");
        assert!(close(dy.lcb, eq.lcb, 10.0), "lcb {lcg:?}");
        assert_eq!(dy.iterations, eq.iterations, "iterations {lcg:?}");
        assert!(
            close(dy.volume_residual, eq.volume_residual, 1.0),
            "volume residual {lcg:?}"
        );
        assert!(
            close(dy.lcb_residual, eq.lcb_residual, 10.0),
            "lcb residual {lcg:?}"
        );
        assert_eq!(dy.fleet.members.len(), eq.fleet.members.len());
        assert_eq!(dy.dynamic, DynamicLoad::default());
        assert_eq!(dy.lift_fraction, 0.0);
    }
    // The offset-lcg case actually exercised the trim path.
    let eq = solve_equilibrium_sectional(
        &hulls(&src),
        &LoadCase {
            mass,
            lcg: Some(0.2),
        },
        RHO,
        &opts,
    )
    .unwrap();
    assert!(
        eq.trim != 0.0 && eq.iterations > 2,
        "trivial reference solve"
    );
}

/// A constant suction of 5 % of the weight must be carried by 5 % more
/// displacement: ρg(∇ − ∇₀) = 0.05 W ⇒ ∇ = 1.05 ∇₀, to the solver tolerance.
#[test]
fn constant_suction_adds_its_share_of_displacement() {
    let (src, mass) = setup();
    let v_target = mass / RHO;
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let weight = mass * G;
    let mut calls = 0usize;
    let dy = solve_dynamic(
        &src,
        &load,
        |_: &FleetState| {
            calls += 1;
            Ok(DynamicLoad {
                force_up: -0.05 * weight,
                moment_bow_up: 0.0,
            })
        },
        None,
    );
    let want = 1.05 * v_target;
    assert!(
        (dy.volume - want).abs() <= V_TOL * v_target,
        "V {} vs {want} (residual {})",
        dy.volume,
        dy.volume_residual
    );
    assert!(
        dy.volume_residual <= V_TOL,
        "residual {}",
        dy.volume_residual
    );
    assert!(
        dy.sinkage > 0.0,
        "suction must sink the hull: s = {}",
        dy.sinkage
    );
    assert!(
        (dy.lift_fraction + 0.05).abs() < 1e-12,
        "lift {}",
        dy.lift_fraction
    );
    assert!((dy.dynamic.force_up + 0.05 * weight).abs() < 1e-9);
    // Symmetric body, centred load, pure force: no trim to speak of.
    assert!(dy.trim.abs() < 1e-3, "trim {}", dy.trim);
    // The closure runs once per iteration for the base evaluation, plus up
    // to one finite-difference probe each for sinkage and trim sensitivity
    // (this case solves both, `lcg: Some(0.0)`) — the Jacobian the Newton
    // step now folds the dynamic load's own local sensitivity into, rather
    // than treating it as a constant added to a purely hydrostatic
    // residual — plus at most one more evaluation for the final report.
    assert!(
        calls <= 3 * dy.iterations + 1,
        "{calls} dynamic evaluations for {} iterations",
        dy.iterations
    );
}

/// A constant bow-up moment on a fore-aft symmetric body trims it bow-up,
/// shifting the LCB aft by M/(ρg∇) while leaving the displacement alone.
#[test]
fn bow_up_moment_trims_bow_up_without_changing_volume() {
    let (src, mass) = setup();
    let v_target = mass / RHO;
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let moment = 0.3 * mass * G; // 0.3 W·m ≈ 1.5° on this hull
    let dy = solve_dynamic(
        &src,
        &load,
        |_: &FleetState| {
            Ok(DynamicLoad {
                force_up: 0.0,
                moment_bow_up: moment,
            })
        },
        None,
    );
    assert!(dy.trim > 1e-3, "bow-up moment gave trim {}", dy.trim);
    assert!(
        (dy.volume - v_target).abs() <= V_TOL * v_target,
        "V {} vs {v_target}",
        dy.volume
    );
    assert!(dy.volume_residual <= V_TOL);
    // Moment balance: ρg(LCB·∇) = −M ⇒ LCB = −M/(ρg∇) = −0.3 m (to 10× the
    // 1e-4·L trim tolerance).
    let want_lcb = -moment / (RHO * G * dy.volume);
    assert!(
        (dy.lcb - want_lcb).abs() <= 1e-2,
        "LCB {} vs {want_lcb}",
        dy.lcb
    );
    assert!(dy.lcb_residual <= 1e-2, "lcb residual {}", dy.lcb_residual);
    assert_eq!(dy.lift_fraction, 0.0);
}

/// Warm-starting from a converged solution goes straight to the fine phase
/// and is accepted within a few iterations, reproducing the solution.
#[test]
fn warm_start_from_solution_converges_immediately() {
    let (src, mass) = setup();
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let weight = mass * G;
    let dynamic = |_: &FleetState| {
        Ok(DynamicLoad {
            force_up: -0.05 * weight,
            moment_bow_up: 0.2 * weight,
        })
    };
    let cold = solve_dynamic(&src, &load, dynamic, None);
    let warm = solve_dynamic(&src, &load, dynamic, Some((cold.sinkage, cold.trim)));
    assert!(
        warm.iterations <= 3,
        "warm start took {} iterations",
        warm.iterations
    );
    assert!(warm.iterations < cold.iterations);
    assert!((warm.sinkage - cold.sinkage).abs() < 1e-6);
    assert!((warm.trim - cold.trim).abs() < 1e-6);
    assert!(warm.volume_residual <= V_TOL && warm.lcb_residual <= 1e-2);
}

/// A dynamic force well past the thin-ship comfort zone is still solved and
/// reported, not refused: 40 % lift leaves 60 % of the displacement.
#[test]
fn large_lift_is_solved_and_reported() {
    let (src, mass) = setup();
    let v_target = mass / RHO;
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let weight = mass * G;
    let dy = solve_dynamic(
        &src,
        &load,
        |_: &FleetState| {
            Ok(DynamicLoad {
                force_up: 0.4 * weight,
                moment_bow_up: 0.0,
            })
        },
        None,
    );
    assert!(
        (dy.lift_fraction - 0.4).abs() < 1e-12,
        "lift {}",
        dy.lift_fraction
    );
    assert!(
        (dy.volume - 0.6 * v_target).abs() <= V_TOL * v_target,
        "V {} vs {}",
        dy.volume,
        0.6 * v_target
    );
    assert!(
        dy.sinkage < 0.0,
        "lift must raise the hull: s = {}",
        dy.sinkage
    );
}

// ---------------------------------------------------------------------------
// End-to-end: the real thin-ship squat closure (sectional), not a
// synthetic constant load — proving the full pipeline from hull geometry to
// solved dynamic attitude.
// ---------------------------------------------------------------------------

use michell::sectional::dynamic_load_closure;
use michell::squat::SquatOptions;
use michell::Conditions;

/// At speed, the real Wigley squat closure should sink the hull (negative
/// force pulls it down at these Froude numbers) relative to the purely
/// hydrostatic solve, and the solved volume should exceed the target
/// (buoyancy must overcome the extra downward pull).
#[test]
fn wigley_dynamic_equilibrium_sinks_relative_to_hydrostatic() {
    let (src, mass) = setup();
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let opts = SectionalOptions::default();

    let hydro = solve_equilibrium_sectional(&hulls(&src), &load, RHO, &opts).unwrap();

    // Fn 0.30 on L = 10 m: well inside where the closure's own tests pin
    // s/L against the Wigley literature values.
    let u = 0.30 * (G * 10.0).sqrt();
    let cond = Conditions::seawater(u);
    let squat_opts = SquatOptions::default();
    let closure = dynamic_load_closure(&cond, 0.0, &squat_opts);
    let dyn_eq = solve_dynamic(&src, &load, closure, None);

    assert!(
        dyn_eq.dynamic.force_up < 0.0,
        "Wigley at Fn 0.30 should feel a downward (suction) force, got {}",
        dyn_eq.dynamic.force_up
    );
    assert!(
        dyn_eq.sinkage > hydro.sinkage,
        "dynamic sinkage {} should exceed hydrostatic sinkage {} (deeper immersion)",
        dyn_eq.sinkage,
        hydro.sinkage
    );
    // Buoyant volume must make up for the downward dynamic force: the solved
    // hydrostatic displacement exceeds the target mass/density.
    assert!(
        dyn_eq.volume > mass / RHO,
        "solved volume {} should exceed target {} to balance the downward pull",
        dyn_eq.volume,
        mass / RHO
    );
    assert!(dyn_eq.lift_fraction < 0.0);
    assert!(dyn_eq.volume_residual < 1e-2);
}

/// The dynamic solve at a very low speed (negligible squat) should agree
/// with the hydrostatic solve closely.
#[test]
fn dynamic_equilibrium_reduces_to_hydrostatic_at_low_speed() {
    let (src, mass) = setup();
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let opts = SectionalOptions::default();

    let hydro = solve_equilibrium_sectional(&hulls(&src), &load, RHO, &opts).unwrap();

    let u = 0.02 * (G * 10.0).sqrt(); // Fn 0.02: squat negligible
    let cond = Conditions::seawater(u);
    let squat_opts = SquatOptions::default();
    let closure = dynamic_load_closure(&cond, 0.0, &squat_opts);
    let dyn_eq = solve_dynamic(&src, &load, closure, None);

    assert!(
        (dyn_eq.sinkage - hydro.sinkage).abs() < 1e-4,
        "at Fn 0.02, dynamic sinkage {} should barely differ from hydrostatic {}",
        dyn_eq.sinkage,
        hydro.sinkage
    );
    assert!(dyn_eq.lift_fraction.abs() < 1e-2);
}

// ---------------------------------------------------------------------------
// Robustness: a large coarse-to-fine geometry/force jump. A coarsened cut
// cannot resolve fine stern detail (a transom, a chine) the way the
// full-resolution one does, so `situate` can hand back a visibly different
// fleet at the very same (sinkage, trim) right at the coarse→fine boundary —
// and with a dynamic closure, the force jumps too. This is a synthetic stand
// -in (two very different Wigley beams for "coarse cut" vs "fine cut" of
// the same nominal hull), and it does exercise a real several-times jump in
// V, Aw, and the dynamic load at the handoff — but it converges cleanly
// within a handful of extra iterations with or without the coarse-phase
// exit damping below, so it does not by itself discriminate that fix.
//
// The fix was motivated and validated against the real case instead: an
// E12-catamaran-derived hull with a transom close to fully immersed, where
// the actual near-field force has a genuine nonlinear, possibly
// non-monotone dependence on attitude near the transom (its own hollow
// length depends on the current transom depth, which depends on attitude),
// not just a step change in magnitude — a limit cycle a synthetic geometry
// swap doesn't reproduce. There, the fix took a case that failed to
// converge (or stalled at several times the tolerance) to a clean solve.
// This test stays as a basic robustness check on the API — the solver
// should not error outright on a large handoff jump — while the sharper
// scenario remains open for a slower, real-geometry regression test.
// ---------------------------------------------------------------------------

#[test]
fn coarse_to_fine_handoff_does_not_stall_convergence() {
    let (l, t, fb) = (10.0, 0.625, 0.3);
    // Two meaningfully different cuts of the same nominal hull, standing in
    // for a coarse cut vs the full-resolution one.
    let coarse_src = wigley_source(l, 1.0, t, fb);
    let fine_src = wigley_source(l, 1.8, t, fb);
    let mass = RHO * 4.0 / 9.0 * 1.0 * l * t; // sized to the coarse beam

    let opts = SectionalOptions::default();
    let situate = |s: f64, tau: f64, coarse: bool| -> Result<FleetState> {
        let src = if coarse { &coarse_src } else { &fine_src };
        let platform = Platform {
            sinkage: s,
            trim: tau,
            pivot_x: 0.0,
        };
        match src.situate_sectional(0, 0.0, &HullPose::default(), &platform, &opts)? {
            Some(h) => Ok(FleetState {
                members: vec![(h.hull, h.placement)],
                dry: 0,
            }),
            None => Ok(FleetState {
                members: vec![],
                dry: 1,
            }),
        }
    };
    // A dynamic load that is a real consequence of the geometry it is
    // handed (the hull's own waterplane area), not a synthetic counter — so
    // it genuinely jumps when `situate`'s geometry jumps at the handoff.
    let dynamic = |fleet: &FleetState| -> Result<DynamicLoad> {
        let area: f64 = fleet.members.iter().map(|(h, _)| h.waterplane_area()).sum();
        let moment: f64 = fleet
            .members
            .iter()
            .map(|(h, p)| (h.lcb_x() + p.x) * h.displaced_volume())
            .sum();
        Ok(DynamicLoad {
            force_up: -400.0 * area,
            moment_bow_up: 300.0 * moment,
        })
    };

    let load = LoadCase {
        mass,
        lcg: Some(0.3),
    };
    let result = solve_equilibrium_dynamic_with(situate, dynamic, &load, RHO, G, None);
    assert!(
        result.is_ok(),
        "equilibrium should converge despite the coarse/fine geometry mismatch: {:?}",
        result.err()
    );
    let eq = result.unwrap();
    assert!(
        eq.volume_residual < 1e-2,
        "residual too loose after the handoff: {}",
        eq.volume_residual
    );
}

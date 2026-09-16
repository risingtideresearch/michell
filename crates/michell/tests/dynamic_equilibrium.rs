//! Dynamic sinkage and trim: the hydrostatic equilibrium solver with a
//! speed-dependent vertical force and pitch moment in the balance, checked
//! against the hydrostatic solver and against closed-form shifts on a
//! fore-aft symmetric Wigley body.

use michell::body::{Body, BodyOptions};
use michell::float::{
    solve_equilibrium_bodies, solve_equilibrium_bodies_dynamic, DynamicEquilibrium, DynamicLoad,
    FleetState, LoadCase,
};
use michell::iges::HullPose;
use michell::{BSplineSurface, Result, STANDARD_GRAVITY};

const RHO: f64 = 1000.0;
const G: f64 = STANDARD_GRAVITY;

/// The Wigley hull `(B/2)(1 − (2x/L)²)(1 − (z/T)²)` as a full-band body: the
/// parabolic sections below the design waterline, carried wall-sided up
/// through a freeboard `fb` so the solver has topside to sink into. Exactly
/// a two-span biquadratic B-spline; the C⁰ knot at the waterline is where
/// the section is flat anyway.
fn wigley_body(l: f64, b: f64, t: f64, fb: f64) -> Body {
    let a = l / 2.0;
    let knots_x = vec![-a, -a, -a, a, a, a];
    let knots_z = vec![0.0, 0.0, 0.0, fb, fb, fb + t, fb + t, fb + t];
    // Bezier control values: 4t(1−t) in x; 1 over the freeboard, then 1−t².
    let gx = [0.0, 2.0, 0.0];
    let hz = [1.0, 1.0, 1.0, 1.0, 0.0];
    let mut control = Vec::with_capacity(15);
    for gi in gx {
        for hj in hz {
            control.push(b / 2.0 * gi * hj);
        }
    }
    let surface = BSplineSurface::new(2, 2, knots_x, knots_z, control).unwrap();
    Body::new(surface, fb, 0.0).unwrap()
}

/// The standard-proportion Wigley body and the mass it floats at its design
/// waterline: `∇ = (4/9) B L T`.
fn setup() -> (Body, f64) {
    let (l, b, t, fb) = (10.0, 1.0, 0.625, 0.3);
    (wigley_body(l, b, t, fb), RHO * 4.0 / 9.0 * b * l * t)
}

fn solve_dynamic(
    body: &Body,
    load: &LoadCase,
    dynamic: impl FnMut(&FleetState) -> Result<DynamicLoad>,
    warm_start: Option<(f64, f64)>,
) -> DynamicEquilibrium {
    solve_equilibrium_bodies_dynamic(
        &[body],
        0.0,
        &[HullPose::default()],
        load,
        RHO,
        G,
        &BodyOptions::default(),
        dynamic,
        warm_start,
    )
    .unwrap()
}

/// Hydrostatic near-miss acceptance: 10× the fine-phase volume tolerance.
const V_TOL: f64 = 2e-3;

/// With a zero dynamic load the dynamic solver is the hydrostatic solver,
/// bit for bit — same iterates, same final state, same reported residuals —
/// both for a sinkage-only balance and for one with trim.
#[test]
fn zero_dynamic_load_reproduces_hydrostatic_solve_exactly() {
    let (body, mass) = setup();
    let bodies = [&body];
    let poses = [HullPose::default()];
    let opts = BodyOptions::default();
    for lcg in [None, Some(0.2)] {
        let load = LoadCase { mass, lcg };
        let eq = solve_equilibrium_bodies(&bodies, 0.0, &poses, &load, RHO, &opts).unwrap();
        let dy = solve_dynamic(&body, &load, |_| Ok(DynamicLoad::default()), None);
        assert_eq!(
            dy.sinkage.to_bits(),
            eq.sinkage.to_bits(),
            "sinkage {lcg:?}"
        );
        assert_eq!(dy.trim.to_bits(), eq.trim.to_bits(), "trim {lcg:?}");
        assert_eq!(dy.volume.to_bits(), eq.volume.to_bits(), "volume {lcg:?}");
        assert_eq!(dy.lcb.to_bits(), eq.lcb.to_bits(), "lcb {lcg:?}");
        assert_eq!(dy.iterations, eq.iterations, "iterations {lcg:?}");
        assert_eq!(
            dy.volume_residual.to_bits(),
            eq.volume_residual.to_bits(),
            "volume residual {lcg:?}"
        );
        assert_eq!(
            dy.lcb_residual.to_bits(),
            eq.lcb_residual.to_bits(),
            "lcb residual {lcg:?}"
        );
        assert_eq!(dy.fleet.members.len(), eq.fleet.members.len());
        assert_eq!(dy.dynamic, DynamicLoad::default());
        assert_eq!(dy.lift_fraction, 0.0);
    }
    // The offset-lcg case actually exercised the trim path.
    let eq = solve_equilibrium_bodies(
        &bodies,
        0.0,
        &poses,
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
    let (body, mass) = setup();
    let v_target = mass / RHO;
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let weight = mass * G;
    let mut calls = 0usize;
    let dy = solve_dynamic(
        &body,
        &load,
        |_| {
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
    // The expensive closure runs once per iteration, and at most once more
    // for the final report.
    assert!(
        calls <= dy.iterations + 1,
        "{calls} dynamic evaluations for {} iterations",
        dy.iterations
    );
}

/// A constant bow-up moment on a fore-aft symmetric body trims it bow-up,
/// shifting the LCB aft by M/(ρg∇) while leaving the displacement alone.
#[test]
fn bow_up_moment_trims_bow_up_without_changing_volume() {
    let (body, mass) = setup();
    let v_target = mass / RHO;
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let moment = 0.3 * mass * G; // 0.3 W·m ≈ 1.5° on this hull
    let dy = solve_dynamic(
        &body,
        &load,
        |_| {
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
    assert!(dy.fleet.band_exceeded == 0, "trim ran out of freeboard");
}

/// Warm-starting from a converged solution goes straight to the fine phase
/// and is accepted within a few iterations, reproducing the solution.
#[test]
fn warm_start_from_solution_converges_immediately() {
    let (body, mass) = setup();
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
    let cold = solve_dynamic(&body, &load, dynamic, None);
    let warm = solve_dynamic(&body, &load, dynamic, Some((cold.sinkage, cold.trim)));
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
    let (body, mass) = setup();
    let v_target = mass / RHO;
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let weight = mass * G;
    let dy = solve_dynamic(
        &body,
        &load,
        |_| {
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
// End-to-end: the real thin-ship squat closure (crate::squat), not a
// synthetic constant load — proving the full pipeline from hull geometry to
// solved dynamic attitude.
// ---------------------------------------------------------------------------

use michell::squat::{dynamic_load_closure, SquatOptions};
use michell::Conditions;

/// At speed, the real Wigley squat closure should sink the hull (negative
/// force pulls it down at these Froude numbers) relative to the purely
/// hydrostatic solve, and the solved volume should exceed the target
/// (buoyancy must overcome the extra downward pull).
#[test]
fn wigley_dynamic_equilibrium_sinks_relative_to_hydrostatic() {
    let (body, mass) = setup();
    let bodies = [&body];
    let poses = [HullPose::default()];
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let opts = BodyOptions::default();

    let hydro = solve_equilibrium_bodies(&bodies, 0.0, &poses, &load, RHO, &opts).unwrap();

    // Fn 0.30 on L = 10 m: well inside where the closure's own tests pin
    // s/L against the Wigley literature values.
    let u = 0.30 * (G * 10.0).sqrt();
    let cond = Conditions::seawater(u);
    let squat_opts = SquatOptions::default();
    let closure = dynamic_load_closure(&cond, 0.0, &squat_opts);
    let dyn_eq =
        solve_equilibrium_bodies_dynamic(&bodies, 0.0, &poses, &load, RHO, G, &opts, closure, None)
            .unwrap();

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
    let (body, mass) = setup();
    let bodies = [&body];
    let poses = [HullPose::default()];
    let load = LoadCase {
        mass,
        lcg: Some(0.0),
    };
    let opts = BodyOptions::default();

    let hydro = solve_equilibrium_bodies(&bodies, 0.0, &poses, &load, RHO, &opts).unwrap();

    let u = 0.02 * (G * 10.0).sqrt(); // Fn 0.02: squat negligible
    let cond = Conditions::seawater(u);
    let squat_opts = SquatOptions::default();
    let closure = dynamic_load_closure(&cond, 0.0, &squat_opts);
    let dyn_eq =
        solve_equilibrium_bodies_dynamic(&bodies, 0.0, &poses, &load, RHO, G, &opts, closure, None)
            .unwrap();

    assert!(
        (dyn_eq.sinkage - hydro.sinkage).abs() < 1e-4,
        "at Fn 0.02, dynamic sinkage {} should barely differ from hydrostatic {}",
        dyn_eq.sinkage,
        hydro.sinkage
    );
    assert!(dyn_eq.lift_fraction.abs() < 1e-2);
}

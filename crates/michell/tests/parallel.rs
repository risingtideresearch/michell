//! The thread fan-outs in the wave-resistance and near-field integrals must
//! not change a single bit of the answer: every parallel path evaluates the
//! same nodes and reduces them in the same order as the serial march. These
//! tests pin that, entry point by entry point, by running each once with a
//! budget of one worker and once with many and demanding exact equality.

use michell::float::{solve_equilibrium_bodies_dynamic, LoadCase};
use michell::iges::HullPose;
use michell::parallel::with_threads;
use michell::squat::{dynamic_load_closure, multihull_dynamic_force, SquatOptions};
use michell::{
    hulls, multihull_heel_wave_resistance, multihull_wave_resistance_lifting,
    multihull_wave_resistance_with, Conditions, Hull, LiftingGrid, Placement, WaveOptions,
};

fn catamaran() -> (Hull, Hull) {
    (
        hulls::wigley(10.0, 0.8, 0.6).unwrap(),
        hulls::wigley(9.0, 0.7, 0.55).unwrap(),
    )
}

fn placed<'a>(a: &'a Hull, b: &'a Hull) -> Vec<(&'a Hull, Placement)> {
    vec![
        (a, Placement { x: 0.0, y: -1.5 }),
        (b, Placement { x: 0.7, y: 1.5 }),
    ]
}

#[test]
fn wave_resistance_is_bitwise_independent_of_thread_count() {
    let (a, b) = catamaran();
    let members = placed(&a, &b);
    let opts = WaveOptions::default();
    for u in [2.0, 3.5, 5.0] {
        let cond = Conditions::seawater(u);
        let serial = with_threads(1, || multihull_wave_resistance_with(&members, &cond, &opts));
        let parallel = with_threads(8, || multihull_wave_resistance_with(&members, &cond, &opts));
        let (s, p) = (serial.unwrap(), parallel.unwrap());
        assert_eq!(s.resistance.to_bits(), p.resistance.to_bits(), "U = {u}");
        assert_eq!(s.est_rel_error.to_bits(), p.est_rel_error.to_bits());
        assert_eq!(s.max_lambda.to_bits(), p.max_lambda.to_bits());
        assert_eq!(s.inner_evaluations, p.inner_evaluations);
    }
}

#[test]
fn heeled_and_lifting_paths_are_bitwise_independent_of_thread_count() {
    let (a, b) = catamaran();
    let members = placed(&a, &b);
    let opts = WaveOptions::default();
    let cond = Conditions::seawater(3.0);

    let s = with_threads(1, || {
        multihull_heel_wave_resistance(&members, &cond, 0.2, &opts).unwrap()
    });
    let p = with_threads(8, || {
        multihull_heel_wave_resistance(&members, &cond, 0.2, &opts).unwrap()
    });
    assert_eq!(s.resistance.to_bits(), p.resistance.to_bits());

    let grid = LiftingGrid { nx: 8, nz: 4 };
    let s = with_threads(1, || {
        multihull_wave_resistance_lifting(&members, &cond, &opts, grid).unwrap()
    });
    let p = with_threads(8, || {
        multihull_wave_resistance_lifting(&members, &cond, &opts, grid).unwrap()
    });
    assert_eq!(s.resistance.to_bits(), p.resistance.to_bits());
}

#[test]
fn dynamic_force_is_bitwise_independent_of_thread_count() {
    let (a, b) = catamaran();
    let members = placed(&a, &b);
    let cond = Conditions::seawater(3.0);
    let opts = SquatOptions::default();
    let s = with_threads(1, || multihull_dynamic_force(&members, &cond, 0.3, &opts).unwrap());
    let p = with_threads(8, || multihull_dynamic_force(&members, &cond, 0.3, &opts).unwrap());
    assert_eq!(s.force_up.to_bits(), p.force_up.to_bits());
    assert_eq!(s.moment_bow_up.to_bits(), p.moment_bow_up.to_bits());
    assert_eq!(s.moment_wave.to_bits(), p.moment_wave.to_bits());
    assert_eq!(s.est_rel_error.to_bits(), p.est_rel_error.to_bits());
    assert_eq!(s.evaluations, p.evaluations);
}

#[test]
fn dynamic_equilibrium_is_bitwise_independent_of_thread_count() {
    // The whole Newton solve — situate, near-field force, resistance — through
    // the closure the manifest sweep uses, on a Wigley body with freeboard
    // (parabolic sections below the design waterline, wall-sided above).
    let (l, b, t, fb) = (10.0, 1.0, 0.625, 0.3);
    let a = l / 2.0;
    let knots_x = vec![-a, -a, -a, a, a, a];
    let knots_z = vec![0.0, 0.0, 0.0, fb, fb, fb + t, fb + t, fb + t];
    let mut control = Vec::with_capacity(15);
    for gi in [0.0, 2.0, 0.0] {
        for hj in [1.0, 1.0, 1.0, 1.0, 0.0] {
            control.push(b / 2.0 * gi * hj);
        }
    }
    let surface = michell::BSplineSurface::new(2, 2, knots_x, knots_z, control).unwrap();
    let body = michell::body::Body::new(surface, fb, 0.0).unwrap();
    let bodies = [&body];
    let poses = [HullPose::default()];
    let cond = Conditions::seawater(3.0);
    let load = LoadCase {
        mass: 2000.0,
        lcg: Some(0.0),
    };
    let squat = SquatOptions::default();
    let bopts = michell::body::BodyOptions::default();
    let solve = || {
        solve_equilibrium_bodies_dynamic(
            &bodies,
            0.0,
            &poses,
            &load,
            cond.fluid.density,
            cond.gravity,
            &bopts,
            dynamic_load_closure(&cond, 0.0, &squat),
            None,
        )
        .unwrap()
    };
    let s = with_threads(1, solve);
    let p = with_threads(8, solve);
    assert_eq!(s.sinkage.to_bits(), p.sinkage.to_bits());
    assert_eq!(s.trim.to_bits(), p.trim.to_bits());
    assert_eq!(s.dynamic.force_up.to_bits(), p.dynamic.force_up.to_bits());
    assert_eq!(s.iterations, p.iterations);
}

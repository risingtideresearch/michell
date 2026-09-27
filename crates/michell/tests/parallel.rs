//! The thread fan-outs in the wave-resistance and near-field integrals must
//! not change a single bit of the answer: every parallel path evaluates the
//! same nodes and reduces them in the same order as the serial march. These
//! tests pin that, entry point by entry point, by running each once with a
//! budget of one worker and once with many and demanding exact equality.

use michell_geometry::float::{solve_equilibrium_sectional_dynamic, LoadCase};
use michell_geometry::iges::{self, HullPose, Platform, SectionalOptions, SourceFleet};
use michell_geometry::parallel::with_threads;
use michell::sectional::{
    dynamic_load_closure, multihull_dynamic_force, multihull_wave_resistance,
};
use michell_geometry::source::SourceHull;
use michell::squat::SquatOptions;
use michell::{Conditions, Placement, SectionalHull, WaveOptions};

/// A Wigley hull's exact CAD surfaces, as a source to cut.
fn wigley_source(length: f64, beam: f64, draft: f64) -> SourceFleet {
    let surfaces = iges::wigley_surfaces(length, beam, draft).unwrap();
    iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap()
}

/// A coarse cut: bit-for-bit independence needs no resolution, only work.
fn opts() -> SectionalOptions {
    SectionalOptions {
        stations: 41,
        rays: 17,
        ..SectionalOptions::default()
    }
}

fn cut(src: &SourceFleet) -> SectionalHull {
    src.situate_sectional(0, 0.0, &HullPose::default(), &Platform::default(), &opts())
        .unwrap()
        .expect("wet")
        .hull
}

fn catamaran() -> (SectionalHull, SectionalHull) {
    (
        cut(&wigley_source(10.0, 0.8, 0.6)),
        cut(&wigley_source(9.0, 0.7, 0.55)),
    )
}

fn placed<'a>(a: &'a SectionalHull, b: &'a SectionalHull) -> Vec<(&'a SectionalHull, Placement)> {
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
    for u in [2.0, 5.0] {
        let cond = Conditions::seawater(u);
        let serial = with_threads(1, || multihull_wave_resistance(&members, &cond, &opts));
        let parallel = with_threads(8, || multihull_wave_resistance(&members, &cond, &opts));
        let (s, p) = (serial.unwrap(), parallel.unwrap());
        assert_eq!(s.resistance.to_bits(), p.resistance.to_bits(), "U = {u}");
        assert_eq!(s.est_rel_error.to_bits(), p.est_rel_error.to_bits());
        assert_eq!(s.max_lambda.to_bits(), p.max_lambda.to_bits());
        assert_eq!(s.inner_evaluations, p.inner_evaluations);
    }
}

#[test]
fn dynamic_force_is_bitwise_independent_of_thread_count() {
    let (a, b) = catamaran();
    let members = placed(&a, &b);
    let cond = Conditions::seawater(3.0);
    let opts = SquatOptions::default();
    let s = with_threads(1, || {
        multihull_dynamic_force(&members, &cond, 0.3, &opts).unwrap()
    });
    let p = with_threads(8, || {
        multihull_dynamic_force(&members, &cond, 0.3, &opts).unwrap()
    });
    assert_eq!(s.force_up.to_bits(), p.force_up.to_bits());
    assert_eq!(s.moment_bow_up.to_bits(), p.moment_bow_up.to_bits());
    assert_eq!(s.moment_wave.to_bits(), p.moment_wave.to_bits());
    assert_eq!(s.est_rel_error.to_bits(), p.est_rel_error.to_bits());
    assert_eq!(s.evaluations, p.evaluations);
}

#[test]
fn dynamic_equilibrium_is_bitwise_independent_of_thread_count() {
    // The whole Newton solve — re-cut, near-field force, resistance — through
    // the closure the manifest sweep uses, on a Wigley loaded light of its
    // design displacement (so it rises and never needs topside).
    let src = wigley_source(10.0, 1.0, 0.625);
    let hulls = [SourceHull {
        source: &src,
        index: 0,
        waterline_z: 0.0,
        pose: HullPose::default(),
    }];
    let cond = Conditions::seawater(3.0);
    let load = LoadCase {
        mass: 2000.0,
        lcg: Some(0.0),
    };
    let squat = SquatOptions::default();
    let opts = opts();
    let solve = || {
        solve_equilibrium_sectional_dynamic(
            &hulls,
            &load,
            cond.fluid.density,
            cond.gravity,
            &opts,
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

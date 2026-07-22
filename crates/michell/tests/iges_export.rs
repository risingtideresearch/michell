//! End-to-end tests for the IGES writer: geometry written by `iges::write`
//! must re-import (through the full sample-and-loft pipeline) as the same
//! hull, and posed exports must land where the pose says.

use michell::iges::{self, HullPose, ImportOptions, Platform};
use michell::{hulls, wave_resistance, Conditions};

#[test]
fn wigley_roundtrips_through_writer_and_importer() {
    let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    // The graph of the half-breadth spline, mirrored: exact geometry.
    let pair = iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
    let text = iges::write(&pair, "wigley").unwrap();
    let (imported, report) = iges::import_hull(&text, &ImportOptions::default()).unwrap();

    assert!(
        report.two_sided,
        "mirrored pair should fold as a full shell"
    );
    assert!(
        report.centerplane.abs() < 1e-6,
        "centerplane {}",
        report.centerplane
    );
    assert!(
        (report.draft - 0.625).abs() < 1e-6,
        "draft {}",
        report.draft
    );

    // The re-import goes through sampling + lofting, so allow loft-level
    // error in the resistance, not exactness.
    let cond = Conditions::seawater(3.0);
    let a = wave_resistance(&hull, &cond).unwrap().resistance;
    let b = wave_resistance(&imported, &cond).unwrap().resistance;
    assert!(
        (a - b).abs() < 0.02 * a,
        "Rw {a} N direct vs {b} N through IGES writer"
    );
}

#[test]
fn posed_pair_imports_as_fleet_at_the_posed_positions() {
    let hull = hulls::wigley(8.0, 0.8, 0.5).unwrap();
    let platform = Platform {
        sinkage: 0.1,
        ..Platform::default()
    };
    let mut surfaces = Vec::new();
    for dy in [2.0, -2.0] {
        let mut pair = iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
        let pose = HullPose {
            dy,
            ..HullPose::default()
        };
        iges::apply_pose(&mut pair, 0.0, &pose, &platform);
        surfaces.extend(pair);
    }
    let text = iges::write(&surfaces, "pair").unwrap();
    let fleet = iges::import_fleet(&text, &ImportOptions::default()).unwrap();
    assert_eq!(fleet.len(), 2, "expected two hulls");
    for (m, want_y) in fleet.iter().zip([-2.0, 2.0]) {
        assert!(
            (m.placement.y - want_y).abs() < 1e-6,
            "centerplane {} vs {want_y}",
            m.placement.y
        );
        // Sinkage was re-expressed as geometry moving down: the draft at the
        // fixed waterline grows by it.
        assert!(
            (m.report.draft - 0.6).abs() < 1e-6,
            "draft {} vs 0.6",
            m.report.draft
        );
    }
}

#[test]
fn source_fleet_posed_surfaces_roundtrip() {
    // Write a hull, read it back as a SourceFleet, re-pose and re-export:
    // the moved control net must be an exact translation of the original.
    let hull = hulls::wigley(6.0, 0.6, 0.4).unwrap();
    let pair = iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
    let text = iges::write(&pair, "one").unwrap();
    let fleet = iges::source_fleet(&text, 0.0).unwrap();
    assert_eq!(fleet.len(), 1);
    let pose = HullPose {
        dx: 1.5,
        dy: -0.75,
        dz: 0.05,
        ..HullPose::default()
    };
    let posed = fleet
        .posed_surfaces(0, 0.0, &pose, &Platform::default())
        .unwrap();
    assert_eq!(posed.len(), 2);
    for (orig, moved) in pair.iter().zip(&posed) {
        for (a, b) in orig.ctrl.iter().zip(&moved.ctrl) {
            assert!((b[0] - (a[0] + 1.5)).abs() < 1e-9);
            assert!((b[1] - (a[1] - 0.75)).abs() < 1e-9);
            assert!((b[2] - (a[2] - 0.05)).abs() < 1e-9);
        }
    }
    // And the re-export of the posed geometry parses.
    let text2 = iges::write(&posed, "one-posed").unwrap();
    assert_eq!(iges::source_fleet(&text2, 0.0).unwrap().len(), 1);
}

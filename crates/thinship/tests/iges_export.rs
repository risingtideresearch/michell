//! End-to-end tests for the IGES writer: geometry written by `iges::write`
//! must re-import (through the sectional cut) as the same hull, and posed
//! exports must land where the pose says.

use hullgeom::iges::{self, HullPose, Platform, SectionalImport, SectionalOptions, SourceFleet};
use thinship::{sectional, Conditions};

fn cut(src: &SourceFleet, idx: usize) -> SectionalImport {
    src.situate_sectional(
        idx,
        0.0,
        &HullPose::default(),
        &Platform::default(),
        &SectionalOptions::default(),
    )
    .unwrap()
    .expect("wet")
}

#[test]
fn wigley_roundtrips_through_writer_and_importer() {
    // The exact mirrored pair, cut directly and through a written file.
    let pair = iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
    let direct = cut(
        &iges::source_fleet_from_surfaces(pair.to_vec(), 1.0, 0.0).unwrap(),
        0,
    );
    let text = iges::write(&pair, "wigley").unwrap();
    let imported = cut(&iges::source_fleet(&text, 0.0).unwrap(), 0);

    let report = &imported.report;
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

    // The writer round-trips every real exactly, so the cut is the same.
    let cond = Conditions::seawater(3.0);
    let w = Default::default();
    let a = sectional::wave_resistance(&direct.hull, &cond, &w)
        .unwrap()
        .resistance;
    let b = sectional::wave_resistance(&imported.hull, &cond, &w)
        .unwrap()
        .resistance;
    assert!(
        (a - b).abs() < 1e-9 * a,
        "Rw {a} N direct vs {b} N through IGES writer"
    );
}

#[test]
fn posed_pair_imports_as_fleet_at_the_posed_positions() {
    let platform = Platform {
        sinkage: 0.1,
        ..Platform::default()
    };
    let mut surfaces = Vec::new();
    for dy in [2.0, -2.0] {
        let mut pair = iges::wigley_surfaces(8.0, 0.8, 0.5).unwrap();
        let pose = HullPose {
            dy,
            ..HullPose::default()
        };
        iges::apply_pose(&mut pair, 0.0, &pose, &platform);
        surfaces.extend(pair);
    }
    let text = iges::write(&surfaces, "pair").unwrap();
    let src = iges::source_fleet(&text, 0.0).unwrap();
    assert_eq!(src.len(), 2, "expected two hulls");
    for (idx, want_y) in [-2.0, 2.0].into_iter().enumerate() {
        let m = cut(&src, idx);
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
    let pair = iges::wigley_surfaces(6.0, 0.6, 0.4).unwrap();
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

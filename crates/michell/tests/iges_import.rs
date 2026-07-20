//! End-to-end IGES import: a synthetic file exercising units conversion (mm),
//! a 124 transformation matrix, port-half mirroring, waterline clipping, and
//! the sampled-inversion + loft pipeline — validated against the native
//! Wigley hull.

use michell::fit::FitOptions;
use michell::iges::{self, ImportOptions};
use michell::{hulls, Conditions};

fn line(content: &str, section: char, seq: usize) -> String {
    format!("{content:<72}{section}{seq:>7}\n")
}

/// Pack comma-separated parameters into 64-column P-section chunks, breaking
/// only at delimiters.
fn pack_params(params: &str, de_ptr: usize, seq_start: usize) -> (String, usize) {
    let mut out = String::new();
    let mut seq = seq_start;
    let mut rest = params;
    while !rest.is_empty() {
        let take = if rest.len() <= 64 {
            rest.len()
        } else {
            // Break after the last delimiter within 64 chars.
            rest[..64]
                .rfind([',', ';'])
                .map(|i| i + 1)
                .expect("a delimiter within 64 columns")
        };
        let (chunk, tail) = rest.split_at(take);
        out.push_str(&format!("{chunk:<64}{de_ptr:>8}P{seq:>7}\n"));
        seq += 1;
        rest = tail;
    }
    (out, seq - seq_start)
}

fn dir_entry(etype: i32, pd_ptr: usize, pd_count: usize, transform_de: usize, seq: usize) -> String {
    let l1 = format!(
        "{etype:>8}{pd_ptr:>8}{z:>8}{z:>8}{z:>8}{z:>8}{transform_de:>8}{z:>8}{z:>8}",
        z = 0
    );
    let l2 = format!(
        "{etype:>8}{z:>8}{z:>8}{pd_count:>8}{z:>8}{blank:>8}{blank:>8}{blank:>8}{z:>8}",
        z = 0,
        blank = ""
    );
    format!("{}{}", line(&l1, 'D', seq), line(&l2, 'D', seq + 1))
}

/// Wigley hull (L=10 m, B=1 m, T=0.625 m) as an IGES file in **millimetres**,
/// z up, keel at file z = 0, as the **port** half (y <= 0), moved up by a
/// 124 transform translating z by +2000 mm. The waterline in the final frame
/// sits at z = 625 + 2000 mm = 2.625 m.
fn wigley_iges(extra_dir: &str, weights_mid: f64) -> String {
    // Biquadratic Bezier: x = 10000 u, z = 625 v (v = 0 at keel),
    // y = -500 * 4u(1-u) * (2v - v^2)   [port half]
    let gx = [0.0, 2.0, 0.0];
    let hv = [0.0, 1.0, 1.0];
    let xs = [0.0, 5000.0, 10000.0];
    let zs = [0.0, 312.5, 625.0];

    let mut s = String::new();
    s.push_str(&line("michell test hull", 'S', 1));

    let global = ",,7Hmichell,10Hwigley.igs,7Hmichell,7Hmichell,32,38,6,308,15,\
                  7Hmichell,1.0,2,2HMM,1,0.01,15H20260719.000000,1E-08,100000.0,\
                  3Havi,7Hmichell,11,0,15H20260719.000000;";
    let mut gseq = 1;
    let mut rest: &str = global;
    while !rest.is_empty() {
        let take = rest.len().min(72);
        s.push_str(&line(&rest[..take], 'G', gseq));
        gseq += 1;
        rest = &rest[take..];
    }

    // P-section bodies.
    let p124 = "124,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,2000.0;";
    let mut p128 = String::from("128,2,2,2,2,0,0,1,0,0");
    for k in ["0.0", "0.0", "0.0", "1.0", "1.0", "1.0"] {
        p128.push_str(&format!(",{k}"));
    }
    for k in ["0.0", "0.0", "0.0", "1.0", "1.0", "1.0"] {
        p128.push_str(&format!(",{k}"));
    }
    // Weights, u index fastest.
    for j in 0..3 {
        for i in 0..3 {
            let w = if i == 1 && j == 1 { weights_mid } else { 1.0 };
            p128.push_str(&format!(",{w:.1}"));
        }
    }
    // Control points, u index fastest.
    for j in 0..3usize {
        for i in 0..3usize {
            let y = -500.0 * gx[i] * hv[j];
            p128.push_str(&format!(",{:.4},{y:.4},{:.4}", xs[i], zs[j]));
        }
    }
    p128.push_str(",0.0,1.0,0.0,1.0;");

    let (p124_text, p124_lines) = pack_params(p124, 1, 1);
    let (p128_text, _) = pack_params(&p128, 3, p124_lines + 1);

    // Directory: 124 at DE 1, 128 at DE 3 referencing the transform at DE 1.
    s.push_str(&dir_entry(124, 1, p124_lines, 0, 1));
    let p128_count = p128_text.lines().count();
    s.push_str(&dir_entry(128, p124_lines + 1, p128_count, 1, 3));
    s.push_str(extra_dir);

    s.push_str(&p124_text);
    s.push_str(&p128_text);
    s.push_str(&line("S      1G      2D      4P     40", 'T', 1));
    s
}

fn import_opts() -> ImportOptions {
    ImportOptions {
        waterline_z: 2.625,
        stations: 61,
        waterlines: 25,
        fit: FitOptions {
            degree_x: 2,
            degree_z: 2,
            n_ctrl_x: 9,
            n_ctrl_z: 7,
            ..FitOptions::default()
        },
        centerplane: None,
    }
}

#[test]
fn imports_wigley_and_reproduces_resistance() {
    let text = wigley_iges("", 1.0);
    let (hull, report) = iges::import_hull(&text, &import_opts()).unwrap();

    assert!((report.units_scale - 0.001).abs() < 1e-15);
    assert!(report.mirrored, "port half must be mirrored");
    assert_eq!(report.patches, 1);
    assert!(!report.two_sided);
    assert_eq!(report.centerplane, 0.0);
    assert!((report.draft - 0.625).abs() < 1e-9, "draft {}", report.draft);
    assert!((report.x_range.0 - 0.0).abs() < 1e-9);
    assert!((report.x_range.1 - 10.0).abs() < 1e-9);
    assert_eq!(report.failed_inversions, 0);
    assert!(
        report.fit.max_residual < 1e-8,
        "fit residual {}",
        report.fit.max_residual
    );

    let reference = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    assert!(
        (hull.displaced_volume() - reference.displaced_volume()).abs()
            < 1e-6 * reference.displaced_volume()
    );
    assert!(
        (hull.wetted_surface() - reference.wetted_surface()).abs()
            < 1e-6 * reference.wetted_surface()
    );
    for u in [2.0, 3.5] {
        let cond = Conditions::seawater(u);
        let rw = michell::wave_resistance(&hull, &cond).unwrap().resistance;
        let rw_ref = michell::wave_resistance(&reference, &cond)
            .unwrap()
            .resistance;
        assert!(
            (rw - rw_ref).abs() < 1e-5 * rw_ref,
            "U={u}: Rw {rw} vs {rw_ref}"
        );
    }
}

#[test]
fn rejects_rational_surface() {
    let text = wigley_iges("", 2.0); // one weight = 2.0 -> rational
    let err = iges::import_hull(&text, &import_opts()).unwrap_err();
    assert!(
        format!("{err}").contains("rational"),
        "unexpected error: {err}"
    );
}

#[test]
fn rejects_trimmed_surface_files() {
    // A 144 (trimmed surface) directory entry is enough to trigger rejection.
    let extra = dir_entry(144, 999, 0, 0, 5);
    let text = wigley_iges(&extra, 1.0);
    let err = iges::import_hull(&text, &import_opts()).unwrap_err();
    assert!(
        format!("{err}").contains("trimmed"),
        "unexpected error: {err}"
    );
}

#[test]
fn rejects_wrong_waterline() {
    // Waterline below the keel: nothing wetted.
    let text = wigley_iges("", 1.0);
    let mut opts = import_opts();
    opts.waterline_z = 1.5; // keel is at 2.0 m after the transform
    let err = iges::import_hull(&text, &opts).unwrap_err();
    assert!(
        format!("{err}").contains("above the specified waterline"),
        "unexpected error: {err}"
    );
}

/// The 4 patch bodies (fore/aft x starboard/port) of a Wigley full shell in
/// metres, z up with the DWL at z = 0.7, centred at y = y0.
fn wigley_shell_bodies(y0: f64) -> Vec<String> {
    let (b, t) = (1.0f64, 0.625f64);
    // Quadratic Bezier pieces of g(x) = 4(x/10)(1 - x/10) split at x = 5.
    let gx_fore = [0.0, 1.0, 1.0];
    let gx_aft = [1.0, 1.0, 0.0];
    let xs_fore = [0.0, 2.5, 5.0];
    let xs_aft = [5.0, 7.5, 10.0];
    // z' = t v downward: h = 1 - v^2 -> Bezier [1, 1, 0]; z_cad = 0.7 - t v.
    let hv = [1.0, 1.0, 0.0];
    let zs = [0.7, 0.7 - t / 2.0, 0.7 - t];

    let mut bodies: Vec<String> = Vec::new();
    for (xs, gx) in [(xs_fore, gx_fore), (xs_aft, gx_aft)] {
        for side in [1.0f64, -1.0] {
            let mut p = String::from("128,2,2,2,2,0,0,1,0,0");
            for _ in 0..2 {
                for k in ["0.0", "0.0", "0.0", "1.0", "1.0", "1.0"] {
                    p.push_str(&format!(",{k}"));
                }
            }
            for _ in 0..9 {
                p.push_str(",1.0");
            }
            for j in 0..3usize {
                for i in 0..3usize {
                    let y = side * b / 2.0 * gx[i] * hv[j] + y0;
                    p.push_str(&format!(",{:.6},{y:.6},{:.6}", xs[i], zs[j]));
                }
            }
            p.push_str(",0.0,1.0,0.0,1.0;");
            bodies.push(p);
        }
    }
    bodies
}

/// Assemble a complete metres-unit IGES file from 128-entity bodies.
fn iges_file_meters(bodies: &[String]) -> String {
    let mut s = String::new();
    s.push_str(&line("michell multipatch test", 'S', 1));
    let global = ",,7Hmichell,9Hmulti.igs,7Hmichell,7Hmichell,32,38,6,308,15,\
                  7Hmichell,1.0,6,1HM,1,0.01,15H20260719.000000,1E-08,100.0,\
                  3Havi,7Hmichell,11,0,15H20260719.000000;";
    let mut gseq = 1;
    let mut rest: &str = global;
    while !rest.is_empty() {
        let take = rest.len().min(72);
        s.push_str(&line(&rest[..take], 'G', gseq));
        gseq += 1;
        rest = &rest[take..];
    }
    // Pack all P bodies first to learn line counts, then directory entries.
    let mut packed = Vec::new();
    let mut p_at = 1usize;
    for (k, body) in bodies.iter().enumerate() {
        let de = 2 * k + 1;
        let (text, n) = pack_params(body, de, p_at);
        packed.push((text, p_at, n));
        p_at += n;
    }
    for (k, (_, ptr, n)) in packed.iter().enumerate() {
        s.push_str(&dir_entry(128, *ptr, *n, 0, 2 * k + 1));
    }
    for (text, _, _) in &packed {
        s.push_str(text);
    }
    s.push_str(&line("S      1G      2D      8P     99", 'T', 1));
    s
}

/// A single offset Wigley full shell (one detected hull at y = 7).
fn wigley_multipatch_iges() -> String {
    iges_file_meters(&wigley_shell_bodies(7.0))
}

/// Like [`wigley_shell_bodies`] but split into upper/lower z bands (8
/// patches), so a shallow waterline leaves the upper band entirely dry.
fn wigley_shell_bodies_zsplit(y0: f64) -> Vec<String> {
    let gx_fore = [0.0, 1.0, 1.0];
    let gx_aft = [1.0, 1.0, 0.0];
    let xs_fore = [0.0, 2.5, 5.0];
    let xs_aft = [5.0, 7.5, 10.0];
    // de Casteljau split at v = 1/2 of hv = [1,1,0] and zs (linear).
    let hv_upper = [1.0, 1.0, 0.75];
    let hv_lower = [0.75, 0.5, 0.0];
    let zs_upper = [0.7, 0.54375, 0.3875];
    let zs_lower = [0.3875, 0.23125, 0.075];
    let mut bodies = Vec::new();
    for (xs, gx) in [(xs_fore, gx_fore), (xs_aft, gx_aft)] {
        for side in [1.0f64, -1.0] {
            for (hv, zs) in [(hv_upper, zs_upper), (hv_lower, zs_lower)] {
                let mut p = String::from("128,2,2,2,2,0,0,1,0,0");
                for _ in 0..2 {
                    for k in ["0.0", "0.0", "0.0", "1.0", "1.0", "1.0"] {
                        p.push_str(&format!(",{k}"));
                    }
                }
                for _ in 0..9 {
                    p.push_str(",1.0");
                }
                for j in 0..3usize {
                    for i in 0..3usize {
                        let y = side * 0.5 * gx[i] * hv[j] + y0;
                        p.push_str(&format!(",{:.6},{y:.6},{:.6}", xs[i], zs[j]));
                    }
                }
                p.push_str(",0.0,1.0,0.0,1.0;");
                bodies.push(p);
            }
        }
    }
    bodies
}

#[test]
fn dry_patches_are_retained_for_deeper_poses() {
    use michell::iges::{HullPose, Platform};
    // Cluster at a shallow reference waterline (upper band of every quadrant
    // fully dry), then situate deeper: the volume must recover the full hull.
    let text = iges_file_meters(&wigley_shell_bodies_zsplit(0.0));
    let src = iges::source_fleet(&text, 0.3).unwrap();
    assert_eq!(src.len(), 1);
    let opts = import_opts();
    let fl = src
        .situate(0.7, &[HullPose::default()], &Platform::default(), &opts)
        .unwrap();
    let v = fl.members[0].hull.displaced_volume();
    let v_full = 4.0 * 1.0 * 10.0 * 0.625 / 9.0;
    assert!(
        (v - v_full).abs() < 1e-3 * v_full,
        "volume {v} vs full {v_full}"
    );
}

#[test]
fn multipatch_full_shell_detects_centerplane_and_matches_wigley() {
    let text = wigley_multipatch_iges();
    let mut opts = import_opts();
    opts.waterline_z = 0.7;
    let (hull, report) = iges::import_hull(&text, &opts).unwrap();

    assert_eq!(report.patches, 4);
    assert!(report.two_sided, "full shell must be detected as two-sided");
    assert!(
        (report.centerplane - 7.0).abs() < 1e-6,
        "centerplane {}",
        report.centerplane
    );
    assert!((report.draft - 0.625).abs() < 1e-9);
    assert!(report.max_asymmetry < 1e-8, "asymmetry {}", report.max_asymmetry);
    assert_eq!(report.failed_inversions, 0);
    assert!(
        report.fit.max_residual < 1e-8,
        "fit residual {}",
        report.fit.max_residual
    );

    let reference = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    assert!(
        (hull.displaced_volume() - reference.displaced_volume()).abs()
            < 1e-6 * reference.displaced_volume()
    );
    let cond = Conditions::seawater(3.0);
    let rw = michell::wave_resistance(&hull, &cond).unwrap().resistance;
    let rw_ref = michell::wave_resistance(&reference, &cond).unwrap().resistance;
    assert!(
        (rw - rw_ref).abs() < 1e-5 * rw_ref,
        "Rw {rw} vs {rw_ref}"
    );

    // Explicit centerplane gives the same answer.
    opts.centerplane = Some(7.0);
    let (_, report2) = iges::import_hull(&text, &opts).unwrap();
    assert_eq!(report2.centerplane, 7.0);
    assert!(report2.two_sided);
}

#[test]
fn trimaran_file_imports_as_fleet_with_detected_placements() {
    // Three full Wigley shells in one file: center hull at y = 0, amas at
    // y = ±7 — a whole trimaran modelled in position.
    let mut bodies = wigley_shell_bodies(0.0);
    bodies.extend(wigley_shell_bodies(7.0));
    bodies.extend(wigley_shell_bodies(-7.0));
    let text = iges_file_meters(&bodies);
    let mut opts = import_opts();
    opts.waterline_z = 0.7;

    let fleet = iges::import_fleet(&text, &opts).unwrap();
    assert_eq!(fleet.len(), 3, "expected 3 hulls");
    let ys: Vec<f64> = fleet.iter().map(|m| m.placement.y).collect();
    for (got, want) in ys.iter().zip([-7.0, 0.0, 7.0]) {
        assert!((got - want).abs() < 1e-6, "placements {ys:?}");
    }
    for m in &fleet {
        assert_eq!(m.report.patches, 4);
        assert!(m.report.two_sided);
        assert!((m.report.draft - 0.625).abs() < 1e-9);
        assert!(m.report.fit.max_residual < 1e-8);
    }

    // Resistance of the imported fleet matches a manually placed fleet of
    // reference Wigley hulls at the same transverse positions.
    let cond = Conditions::seawater(3.0);
    let members: Vec<(&michell::Hull, michell::Placement)> = fleet
        .iter()
        .map(|m| (&m.hull, m.placement))
        .collect();
    let got = michell::multihull_wave_resistance(&members, &cond)
        .unwrap()
        .resistance;
    let reference = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let refs = [
        (&reference, michell::Placement { x: 0.0, y: -7.0 }),
        (&reference, michell::Placement { x: 0.0, y: 0.0 }),
        (&reference, michell::Placement { x: 0.0, y: 7.0 }),
    ];
    let want = michell::multihull_wave_resistance(&refs, &cond)
        .unwrap()
        .resistance;
    assert!(
        (got - want).abs() < 1e-5 * want,
        "trimaran Rw {got} vs reference {want}"
    );

    // Single-hull import must refuse the multihull file with a clear count.
    let err = iges::import_hull(&text, &opts).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("3 separate hulls"), "unexpected error: {msg}");

    // A centerplane override is ambiguous for a multi-hull file.
    opts.centerplane = Some(0.0);
    let err = iges::import_fleet(&text, &opts).unwrap_err();
    assert!(
        format!("{err}").contains("ambiguous"),
        "unexpected error: {err}"
    );
}

#[test]
fn offset_one_sided_hull_is_rejected_without_centerplane() {
    // A one-sided hull far from y = 0 must error with advice (rather than
    // silently producing a ~7 m half-beam), and import cleanly once the
    // centerplane is supplied.
    let half = wigley_iges_offset_starboard();
    let mut opts = import_opts();
    opts.waterline_z = 0.7;
    let err = iges::import_hull(&half, &opts).unwrap_err();
    assert!(
        format!("{err}").contains("centerplane"),
        "unexpected error: {err}"
    );

    opts.centerplane = Some(7.0);
    let (hull, report) = iges::import_hull(&half, &opts).unwrap();
    assert!(!report.two_sided);
    assert_eq!(report.centerplane, 7.0);
    let reference = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let cond = Conditions::seawater(3.0);
    let rw = michell::wave_resistance(&hull, &cond).unwrap().resistance;
    let rw_ref = michell::wave_resistance(&reference, &cond).unwrap().resistance;
    assert!(
        (rw - rw_ref).abs() < 1e-5 * rw_ref,
        "Rw {rw} vs {rw_ref}"
    );
}

/// Single starboard-side Wigley patch offset to y ~ 7 (never reaches y = 0).
fn wigley_iges_offset_starboard() -> String {
    let (b, t, y0) = (1.0f64, 0.625f64, 7.0f64);
    let gx = [0.0, 2.0, 0.0];
    let hv = [1.0, 1.0, 0.0];
    let xs = [0.0, 5.0, 10.0];
    let zs = [0.7, 0.7 - t / 2.0, 0.7 - t];
    let mut p = String::from("128,2,2,2,2,0,0,1,0,0");
    for _ in 0..2 {
        for k in ["0.0", "0.0", "0.0", "1.0", "1.0", "1.0"] {
            p.push_str(&format!(",{k}"));
        }
    }
    for _ in 0..9 {
        p.push_str(",1.0");
    }
    for j in 0..3usize {
        for i in 0..3usize {
            let y = b / 2.0 * gx[i] * hv[j] + y0;
            p.push_str(&format!(",{:.6},{y:.6},{:.6}", xs[i], zs[j]));
        }
    }
    p.push_str(",0.0,1.0,0.0,1.0;");

    let mut s = String::new();
    s.push_str(&line("michell offset half test", 'S', 1));
    let global = ",,7Hmichell,8Hhalf.igs,7Hmichell,7Hmichell,32,38,6,308,15,\
                  7Hmichell,1.0,6,1HM,1,0.01,15H20260719.000000,1E-08,100.0,\
                  3Havi,7Hmichell,11,0,15H20260719.000000;";
    let mut gseq = 1;
    let mut rest: &str = global;
    while !rest.is_empty() {
        let take = rest.len().min(72);
        s.push_str(&line(&rest[..take], 'G', gseq));
        gseq += 1;
        rest = &rest[take..];
    }
    let (ptext, n) = pack_params(&p, 1, 1);
    s.push_str(&dir_entry(128, 1, n, 0, 1));
    s.push_str(&ptext);
    s.push_str(&line("S      1G      2D      2P     40", 'T', 1));
    s
}

#[test]
fn situate_dz_equals_waterline_shift() {
    use michell::iges::{HullPose, Platform};
    // Raising a hull by 0.1 m is the same wetted geometry as lowering the
    // waterline by 0.1 m.
    let text = iges_file_meters(&wigley_shell_bodies(7.0));
    let src = iges::source_fleet(&text, 0.7).unwrap();
    assert_eq!(src.len(), 1);
    let opts = import_opts();
    let raised = src
        .situate(
            0.7,
            &[HullPose {
                dz: -0.1,
                ..Default::default()
            }],
            &Platform::default(),
            &opts,
        )
        .unwrap();
    let lowered_wl = src
        .situate(0.6, &[HullPose::default()], &Platform::default(), &opts)
        .unwrap();
    let (a, b) = (&raised.members[0].hull, &lowered_wl.members[0].hull);
    assert!(
        (a.displaced_volume() - b.displaced_volume()).abs() < 1e-9 * b.displaced_volume(),
        "vol {} vs {}",
        a.displaced_volume(),
        b.displaced_volume()
    );
    assert!((a.draft() - b.draft()).abs() < 1e-9);
    let cond = Conditions::seawater(3.0);
    let rw_a = michell::wave_resistance(a, &cond).unwrap().resistance;
    let rw_b = michell::wave_resistance(b, &cond).unwrap().resistance;
    assert!((rw_a - rw_b).abs() < 1e-8 * rw_b, "Rw {rw_a} vs {rw_b}");
}

#[test]
fn situate_trim_is_symmetric_for_symmetric_hull() {
    use michell::iges::{HullPose, Platform};
    let text = iges_file_meters(&wigley_shell_bodies(0.0));
    let src = iges::source_fleet(&text, 0.7).unwrap();
    let opts = import_opts();
    let vol = |trim: f64| -> f64 {
        let fl = src
            .situate(
                0.7,
                &[HullPose {
                    trim,
                    ..Default::default()
                }],
                &Platform::default(),
                &opts,
            )
            .unwrap();
        fl.members[0].hull.displaced_volume()
    };
    let v0 = vol(0.0);
    let vp = vol(3.0f64.to_radians());
    let vm = vol((-3.0f64).to_radians());
    assert!(
        (vp - vm).abs() < 1e-5 * v0,
        "trim asymmetry: {vp} vs {vm} (v0 {v0})"
    );
    // Trimming a fore-aft symmetric hull about its midpoint is a second-order
    // volume effect — noticeable at 3 degrees on a 10 m hull, but bounded.
    assert!((vp - v0).abs() < 0.15 * v0, "vp {vp} vs v0 {v0}");
}

#[test]
fn equilibrium_matches_analytic_wigley() {
    use michell::float::{solve_equilibrium, LoadCase};
    use michell::iges::HullPose;
    // Wigley shell at design draft T0 = 0.625 under waterline 0.7. Target
    // immersion d = 0.5 -> analytic volume and sinkage = -0.125.
    let (l, b, t0, d) = (10.0f64, 1.0f64, 0.625f64, 0.5f64);
    let c = t0 - d;
    let v_analytic = b * (2.0 * l / 3.0) * (d - t0 / 3.0 + c.powi(3) / (3.0 * t0 * t0));
    let density = 1025.9;
    let mass = density * v_analytic;

    let text = iges_file_meters(&wigley_shell_bodies(0.0));
    let src = iges::source_fleet(&text, 0.7).unwrap();
    let opts = import_opts();

    // Weight-only balance (trim locked).
    let eq = solve_equilibrium(
        &src,
        0.7,
        &[HullPose::default()],
        &LoadCase { mass, lcg: None },
        density,
        &opts,
    )
    .unwrap();
    assert!(
        (eq.sinkage + 0.125).abs() < 2e-3,
        "sinkage {} (want -0.125)",
        eq.sinkage
    );
    assert!(eq.volume_residual < 5e-4, "vol residual {}", eq.volume_residual);
    assert_eq!(eq.trim, 0.0);
    assert_eq!(eq.fleet.dry, 0);

    // The shell spans x in [0, 10], so its symmetry plane is x = 5: with lcg
    // there the trim must stay ~0.
    let eq = solve_equilibrium(
        &src,
        0.7,
        &[HullPose::default()],
        &LoadCase {
            mass,
            lcg: Some(5.0),
        },
        density,
        &opts,
    )
    .unwrap();
    assert!(eq.trim.abs() < 1e-3, "trim {}", eq.trim);
    assert!(eq.lcb_residual < 1e-3, "lcb residual {}", eq.lcb_residual);

    // Shift the CG forward: solver must trim until LCB follows.
    let eq = solve_equilibrium(
        &src,
        0.7,
        &[HullPose::default()],
        &LoadCase {
            mass,
            lcg: Some(5.3),
        },
        density,
        &opts,
    )
    .unwrap();
    assert!(eq.volume_residual < 5e-4);
    assert!(eq.lcb_residual < 1.5e-3, "lcb residual {}", eq.lcb_residual);
    assert!((eq.lcb - 5.3).abs() < 1.5e-3, "lcb {}", eq.lcb);
    assert!(eq.trim.abs() > 1e-3, "expected nonzero trim, got {}", eq.trim);

    // An impossible CG (at the bow tip) must fail with a diagnosis, not hang.
    let err = solve_equilibrium(
        &src,
        0.7,
        &[HullPose::default()],
        &LoadCase {
            mass,
            lcg: Some(0.0),
        },
        density,
        &opts,
    )
    .unwrap_err();
    assert!(
        format!("{err}").contains("unreachable"),
        "unexpected error: {err}"
    );
}

#[test]
fn import_grid_carries_accurate_slopes() {
    // The synthetic file is the exact biquadratic Wigley (L=10, B=1,
    // T=0.625) with x ∈ [0, 10], so the sampled slope channels recovered
    // from the Newton Jacobian must match the analytic derivatives of
    // f(x, z) = (10x − x²)/50 · (1 − (z/T)²).
    let text = wigley_iges("", 1.0);
    let fleet = iges::import_fleet(&text, &import_opts()).unwrap();
    assert_eq!(fleet.len(), 1);
    let g = &fleet[0].grid;
    let (fx, fz) = (g.fx().expect("fx channel"), g.fz().expect("fz channel"));
    let t = 0.625f64;
    let mut checked = 0usize;
    for (i, &x) in g.stations().iter().enumerate() {
        for (j, &z) in g.waterlines().iter().enumerate() {
            let s = g.idx(i, j);
            // Interior samples only: away from the fold (f = 0) and with a
            // recovered slope.
            if g.half_beams()[s] < 1e-3 || !(fx[s].is_finite() && fz[s].is_finite()) {
                continue;
            }
            let want_fx = (10.0 - 2.0 * x) / 50.0 * (1.0 - (z / t).powi(2));
            let want_fz = -(10.0 * x - x * x) / 25.0 * z / (t * t);
            assert!(
                (fx[s] - want_fx).abs() < 1e-8,
                "x={x} z={z}: fx {} vs {want_fx}",
                fx[s]
            );
            assert!(
                (fz[s] - want_fz).abs() < 1e-8,
                "x={x} z={z}: fz {} vs {want_fz}",
                fz[s]
            );
            checked += 1;
        }
    }
    // Most of the wet grid must carry slopes (gaps allowed only at folds
    // and footprint edges).
    let total = g.stations().len() * g.waterlines().len();
    assert!(
        checked > total / 3,
        "only {checked} of {total} samples carried slopes"
    );
    assert!(fleet[0].report.fit.fx_residual.is_some());
}

#[test]
fn body_situate_matches_iges_situate() {
    use michell::body::{Body, BodyOptions};
    use michell::iges::{HullPose, Platform};
    // Loft the shell over its full band (this generator's hull tops out at
    // the DWL, so the band top is z_cad = 0.7 and the design waterline sits
    // at depth 0 below it), build a Body, and compare re-situating through
    // the body against re-situating through the IGES source.
    let text = iges_file_meters(&wigley_shell_bodies(0.0));
    let src = iges::source_fleet(&text, 0.7).unwrap();
    let opts = import_opts();
    let full = src
        .situate_one(0, 0.7, &HullPose::default(), &Platform::default(), &opts)
        .unwrap()
        .expect("wet");
    let body = Body::new(full.hull.surface().clone(), 0.0, full.report.centerplane).unwrap();
    let bopts = BodyOptions {
        stations: opts.stations,
        waterlines: opts.waterlines,
        fit: opts.fit,
    };

    // dz: raise by 0.1 — body vs IGES at the equivalent lowered waterline.
    let pose = HullPose {
        dz: -0.1,
        ..Default::default()
    };
    let via_body = body
        .situate(0.0, &pose, &Platform::default(), &bopts)
        .unwrap()
        .expect("wet");
    let via_iges = src
        .situate_one(0, 0.6, &HullPose::default(), &Platform::default(), &opts)
        .unwrap()
        .expect("wet");
    let (vb, vi) = (
        via_body.hull.displaced_volume(),
        via_iges.hull.displaced_volume(),
    );
    assert!((vb - vi).abs() < 1e-6 * vi, "dz: vol {vb} vs {vi}");
    let cond = Conditions::seawater(3.0);
    let rwb = michell::wave_resistance(&via_body.hull, &cond).unwrap().resistance;
    let rwi = michell::wave_resistance(&via_iges.hull, &cond).unwrap().resistance;
    assert!((rwb - rwi).abs() < 1e-4 * rwi, "dz: Rw {rwb} vs {rwi}");
    assert_eq!(via_body.band_exceeded, 0);

    // trim: 2 degrees through both paths.
    let pose = HullPose {
        trim: 2.0f64.to_radians(),
        ..Default::default()
    };
    let via_body = body
        .situate(0.0, &pose, &Platform::default(), &bopts)
        .unwrap()
        .expect("wet");
    let via_iges = src
        .situate_one(0, 0.7, &pose, &Platform::default(), &opts)
        .unwrap()
        .expect("wet");
    let (vb, vi) = (
        via_body.hull.displaced_volume(),
        via_iges.hull.displaced_volume(),
    );
    assert!((vb - vi).abs() < 1e-3 * vi, "trim: vol {vb} vs {vi}");
    let (db, di) = (via_body.hull.draft(), via_iges.hull.draft());
    assert!((db - di).abs() < 1e-3 * di, "trim: draft {db} vs {di}");
}

/// Tessellate Wigley full shells (L=10, B=1, T=0.625, DWL at z_cad = 0.7)
/// into ASCII STL, one shell per y offset.
fn wigley_stl_ascii(y0s: &[f64], nx: usize, nz: usize) -> String {
    let f = |x: f64, zp: f64| {
        0.5 * (4.0 * (x / 10.0) * (1.0 - x / 10.0)) * (1.0 - (zp / 0.625f64).powi(2))
    };
    let mut s = String::from("solid wigley\n");
    let mut tri = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
        s.push_str(" facet normal 0 0 0\n  outer loop\n");
        for v in [a, b, c] {
            s.push_str(&format!("   vertex {} {} {}\n", v[0], v[1], v[2]));
        }
        s.push_str("  endloop\n endfacet\n");
    };
    for &y0 in y0s {
        for side in [1.0f64, -1.0] {
            for i in 0..nx {
                for j in 0..nz {
                    let (x0, x1) = (
                        10.0 * i as f64 / nx as f64,
                        10.0 * (i + 1) as f64 / nx as f64,
                    );
                    let (z0, z1) = (
                        0.625 * j as f64 / nz as f64,
                        0.625 * (j + 1) as f64 / nz as f64,
                    );
                    let p = |x: f64, zp: f64| [x, side * f(x, zp) + y0, 0.7 - zp];
                    tri(p(x0, z0), p(x1, z0), p(x1, z1));
                    tri(p(x0, z0), p(x1, z1), p(x0, z1));
                }
            }
        }
    }
    s.push_str("endsolid wigley\n");
    s
}

#[test]
fn stl_import_matches_reference_wigley() {
    use michell::iges::{HullPose, Platform};
    let stl = wigley_stl_ascii(&[0.0], 160, 48);
    let mf = michell::stl::mesh_fleet(stl.as_bytes(), 1.0, 0.7).unwrap();
    assert_eq!(mf.len(), 1, "one hull expected");
    let opts = import_opts();
    let fl = mf
        .situate(0.7, &[HullPose::default()], &Platform::default(), &opts)
        .unwrap();
    let m = &fl.members[0];
    assert!(m.report.two_sided);
    assert!(m.report.centerplane.abs() < 1e-6);
    assert!((m.report.draft - 0.625).abs() < 1e-9);
    assert!(m.report.max_asymmetry < 1e-9);

    let reference = hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let (v, vr) = (m.hull.displaced_volume(), reference.displaced_volume());
    assert!((v - vr).abs() < 1e-3 * vr, "volume {v} vs {vr}");
    let cond = Conditions::seawater(3.0);
    let rw = michell::wave_resistance(&m.hull, &cond).unwrap().resistance;
    let rw_ref = michell::wave_resistance(&reference, &cond)
        .unwrap()
        .resistance;
    assert!(
        (rw - rw_ref).abs() < 1e-2 * rw_ref,
        "Rw {rw} vs {rw_ref}"
    );
}

#[test]
fn stl_catamaran_clusters_into_two_hulls() {
    use michell::iges::{HullPose, Platform};
    let stl = wigley_stl_ascii(&[3.0, -3.0], 60, 20);
    let mf = michell::stl::mesh_fleet(stl.as_bytes(), 1.0, 0.7).unwrap();
    assert_eq!(mf.len(), 2, "two hulls expected");
    let opts = import_opts();
    let fl = mf
        .situate(
            0.7,
            &[HullPose::default(), HullPose::default()],
            &Platform::default(),
            &opts,
        )
        .unwrap();
    assert_eq!(fl.members.len(), 2);
    assert!((fl.members[0].placement.y + 3.0).abs() < 1e-3);
    assert!((fl.members[1].placement.y - 3.0).abs() < 1e-3);
}

#[test]
fn parse_reports_inventory() {
    let text = wigley_iges("", 1.0);
    let file = iges::parse(&text).unwrap();
    assert_eq!(file.surfaces.len(), 1);
    assert!(file.entity_counts.contains(&(124, 1)));
    assert!(file.entity_counts.contains(&(128, 1)));
    let s = &file.surfaces[0];
    assert_eq!((s.degree_u, s.degree_v), (2, 2));
    // Units + transform applied: x in metres, z shifted by 2 m.
    let xs: Vec<f64> = s.ctrl.iter().map(|p| p[0]).collect();
    assert!(xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - 10.0 < 1e-12);
    let zmax = s.ctrl.iter().map(|p| p[2]).fold(f64::NEG_INFINITY, f64::max);
    assert!((zmax - 2.625).abs() < 1e-12, "zmax {zmax}");
}

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

/// Wigley hull as a 4-patch full shell (fore/aft x starboard/port), in
/// metres, z up with the DWL at z = 0.7, offset 7 m to starboard — the shape
/// of a real multi-patch export of an outrigger hull.
fn wigley_multipatch_iges() -> String {
    let (b, t, y0) = (1.0f64, 0.625f64, 7.0f64);
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

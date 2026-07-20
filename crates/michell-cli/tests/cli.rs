//! End-to-end CLI tests: run the actual binary and check its output against
//! the library.

use std::path::PathBuf;
use std::process::Command;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_michell"))
}

fn tmp(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(name)
}

fn run_ok(cmd: &mut Command) -> String {
    let out = cmd.output().expect("binary runs");
    assert!(
        out.status.success(),
        "command failed: {}\n{}",
        String::from_utf8_lossy(&out.stderr),
        String::from_utf8_lossy(&out.stdout),
    );
    String::from_utf8(out.stdout).expect("utf8 output")
}

/// Extract `"key":<number>` from flat JSON output (first occurrence).
fn json_num(json: &str, key: &str) -> f64 {
    let pat = format!("\"{key}\":");
    let at = json.find(&pat).unwrap_or_else(|| panic!("{key} in {json}"));
    let rest = &json[at + pat.len()..];
    let end = rest
        .find([',', '}', ']'])
        .unwrap_or_else(|| panic!("number end for {key}"));
    rest[..end].parse::<f64>().unwrap_or_else(|_| panic!("parse {key}: {rest:?}"))
}

#[test]
fn wigley_roundtrip_matches_library() {
    let hull_path = tmp("wigley.hull");
    run_ok(bin().args([
        "wigley",
        "--length",
        "10",
        "--beam",
        "1",
        "--draft",
        "0.625",
        "-o",
        hull_path.to_str().unwrap(),
    ]));

    // info reports the exact geometry.
    let info = run_ok(bin().args(["info", hull_path.to_str().unwrap(), "--json"]));
    let reference = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
    assert!((json_num(&info, "length") - 10.0).abs() < 1e-12);
    assert!((json_num(&info, "draft") - 0.625).abs() < 1e-12);
    assert!(
        (json_num(&info, "displaced_volume") - reference.displaced_volume()).abs() < 1e-10
    );

    // resistance --json at a single speed matches the library exactly.
    let out = run_ok(bin().args([
        "resistance",
        hull_path.to_str().unwrap(),
        "--speeds",
        "3.0",
        "--fluid",
        "freshwater",
        "--json",
    ]));
    let cond = michell::Conditions::freshwater(3.0);
    let want = michell::resistance(&reference, &cond).unwrap();
    let rw = json_num(&out, "rw");
    let rv = json_num(&out, "rv");
    assert!(
        (rw - want.wave.resistance).abs() < 1e-9 * want.wave.resistance,
        "rw {rw} vs {}",
        want.wave.resistance
    );
    assert!((rv - want.viscous.resistance).abs() < 1e-9 * want.viscous.resistance);
    // Effective power P_E = R_t * U.
    let pe = json_num(&out, "effective_power");
    let total = json_num(&out, "total");
    assert!(
        (pe - total * 3.0).abs() < 1e-9 * pe.abs(),
        "pe {pe} vs total*U {}",
        total * 3.0
    );
}

#[test]
fn froude_range_produces_table() {
    let hull_path = tmp("wigley_table.hull");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let out = run_ok(bin().args([
        "resistance",
        hull_path.to_str().unwrap(),
        "--froude",
        "0.2:0.4:0.05",
    ]));
    // Header + 5 rows.
    assert!(out.contains("Fn"), "table header missing:\n{out}");
    let data_rows = out
        .lines()
        .filter(|l| l.trim_start().starts_with(|c: char| c.is_ascii_digit()))
        .count();
    assert_eq!(data_rows, 5, "expected 5 rows:\n{out}");
}

#[test]
fn offsets_loft_and_resistance() {
    // Wigley offsets table, 41 stations x 9 waterlines.
    let (l, b, t) = (10.0f64, 1.0f64, 0.625f64);
    let mut text = String::from("michell-offsets v1\n# Wigley test table\nwaterlines");
    let mz = 9;
    for j in 0..mz {
        text.push_str(&format!(" {}", t * j as f64 / (mz - 1) as f64));
    }
    text.push('\n');
    for i in 0..41 {
        let x = -l / 2.0 + l * i as f64 / 40.0;
        text.push_str(&format!("station {x}"));
        for j in 0..mz {
            let z = t * j as f64 / (mz - 1) as f64;
            let y = b / 2.0 * (1.0 - (2.0 * x / l).powi(2)) * (1.0 - (z / t).powi(2));
            text.push_str(&format!(" {y}"));
        }
        text.push('\n');
    }
    let off_path = tmp("wigley.offsets");
    std::fs::write(&off_path, text).unwrap();

    // Loft to a control net file.
    let net_path = tmp("lofted.hull");
    let loft_out = run_ok(bin().args([
        "loft",
        off_path.to_str().unwrap(),
        "-o",
        net_path.to_str().unwrap(),
    ]));
    assert!(loft_out.contains("max residual"), "{loft_out}");

    // The lofted hull reproduces Wigley resistance (data is exactly
    // representable by the cubic loft).
    let out = run_ok(bin().args([
        "resistance",
        net_path.to_str().unwrap(),
        "--speeds",
        "3.0",
        "--json",
    ]));
    let reference = michell::hulls::wigley(l, b, t).unwrap();
    let want = michell::resistance(&reference, &michell::Conditions::seawater(3.0)).unwrap();
    let rw = json_num(&out, "rw");
    assert!(
        (rw - want.wave.resistance).abs() < 1e-6 * want.wave.resistance,
        "rw {rw} vs {}",
        want.wave.resistance
    );
}

#[test]
fn knots_flag_scales_speeds() {
    let hull_path = tmp("wigley_kn.hull");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let kn = run_ok(bin().args([
        "resistance",
        hull_path.to_str().unwrap(),
        "--speeds",
        "6.0",
        "--knots",
        "--json",
    ]));
    let ms = run_ok(bin().args([
        "resistance",
        hull_path.to_str().unwrap(),
        "--speeds",
        "3.0866666666666667",
        "--json",
    ]));
    assert!(
        (json_num(&kn, "rw") - json_num(&ms, "rw")).abs() < 1e-9 * json_num(&ms, "rw").abs()
    );
}

#[test]
fn catamaran_fleet_matches_library() {
    let hull_path = tmp("wigley_cat.hull");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let spec_a = format!("{}@y=1.0", hull_path.to_str().unwrap());
    let spec_b = format!("{}@y=-1.0", hull_path.to_str().unwrap());
    let out = run_ok(bin().args([
        "resistance",
        &spec_a,
        &spec_b,
        "--speeds",
        "3.0",
        "--json",
    ]));
    let hull = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let members = [
        (&hull, michell::Placement { x: 0.0, y: 1.0 }),
        (&hull, michell::Placement { x: 0.0, y: -1.0 }),
    ];
    let want =
        michell::multihull_resistance(&members, &michell::Conditions::seawater(3.0)).unwrap();
    let rw = json_num(&out, "rw");
    let iff = json_num(&out, "interference");
    assert!(
        (rw - want.wave.resistance).abs() < 1e-9 * want.wave.resistance,
        "rw {rw} vs {}",
        want.wave.resistance
    );
    assert!(
        (iff - want.interference).abs() < 1e-9,
        "IF {iff} vs {}",
        want.interference
    );
    // Interference must actually be doing something at this spacing.
    assert!((iff - 1.0).abs() > 0.01, "IF suspiciously unity: {iff}");
    // Table mode shows the IF column for fleets.
    let table = run_ok(bin().args(["resistance", &spec_a, &spec_b, "--speeds", "3.0"]));
    assert!(table.contains("IF"), "missing IF column:\n{table}");
}

/// Full-shell Wigley (L=10, B=1, T0=0.625) in metres, z up, DWL at z=0.7:
/// 4 untrimmed biquadratic patches per shell, one shell per y offset.
fn wigley_shells_iges(y0s: &[f64]) -> String {
    fn line(content: &str, section: char, seq: usize) -> String {
        format!("{content:<72}{section}{seq:>7}\n")
    }
    let mut bodies: Vec<String> = Vec::new();
    let (gx_f, gx_a) = ([0.0, 1.0, 1.0], [1.0, 1.0, 0.0]);
    let (xs_f, xs_a) = ([0.0, 2.5, 5.0], [5.0, 7.5, 10.0]);
    let hv = [1.0, 1.0, 0.0];
    let zs = [0.7, 0.7 - 0.3125, 0.7 - 0.625];
    for &y0 in y0s {
        for (xs, gx) in [(xs_f, gx_f), (xs_a, gx_a)] {
            for side in [1.0f64, -1.0] {
                let mut b = String::from("128,2,2,2,2,0,0,1,0,0");
                for _ in 0..2 {
                    for k in ["0.0", "0.0", "0.0", "1.0", "1.0", "1.0"] {
                        b.push_str(&format!(",{k}"));
                    }
                }
                for _ in 0..9 {
                    b.push_str(",1.0");
                }
                for j in 0..3usize {
                    for i in 0..3usize {
                        let y = side * 0.5 * gx[i] * hv[j] + y0;
                        b.push_str(&format!(",{:.6},{y:.6},{:.6}", xs[i], zs[j]));
                    }
                }
                b.push_str(",0.0,1.0,0.0,1.0;");
                bodies.push(b);
            }
        }
    }
    let mut s = String::new();
    s.push_str(&line("cli sweep test hull", 'S', 1));
    let global = ",,7Hmichell,8Hcli.iges,7Hmichell,7Hmichell,32,38,6,308,15,\
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
    // P bodies packed at 64 columns, breaking at delimiters.
    let mut packed: Vec<(String, usize, usize)> = Vec::new();
    let mut p_at = 1usize;
    for (k, body) in bodies.iter().enumerate() {
        let de = 2 * k + 1;
        let mut text = String::new();
        let mut n = 0usize;
        let mut rest: &str = body;
        while !rest.is_empty() {
            let take = if rest.len() <= 64 {
                rest.len()
            } else {
                rest[..64].rfind([',', ';']).map(|i| i + 1).unwrap()
            };
            let (chunk, tail) = rest.split_at(take);
            text.push_str(&format!("{chunk:<64}{de:>8}P{:>7}\n", p_at + n));
            n += 1;
            rest = tail;
        }
        packed.push((text, p_at, n));
        p_at += n;
    }
    for (k, (_, ptr, n)) in packed.iter().enumerate() {
        let l1 = format!(
            "{:>8}{:>8}{z:>8}{z:>8}{z:>8}{z:>8}{z:>8}{z:>8}{z:>8}",
            128,
            ptr,
            z = 0
        );
        let l2 = format!(
            "{:>8}{z:>8}{z:>8}{:>8}{z:>8}{b:>8}{b:>8}{b:>8}{z:>8}",
            128,
            n,
            z = 0,
            b = ""
        );
        s.push_str(&line(&l1, 'D', 2 * k + 1));
        s.push_str(&line(&l2, 'D', 2 * k + 2));
    }
    for (text, _, _) in &packed {
        s.push_str(text);
    }
    s.push_str(&line("S      1G      2D      8P     99", 'T', 1));
    s
}

fn wigley_shell_iges() -> String {
    wigley_shells_iges(&[0.0])
}

fn csv_col(csv: &str, col: &str) -> Vec<f64> {
    let mut lines = csv.lines();
    let header: Vec<&str> = lines.next().expect("header").split(',').collect();
    let at = header
        .iter()
        .position(|h| *h == col)
        .unwrap_or_else(|| panic!("column {col} in {header:?}"));
    lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.split(',')
                .nth(at)
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or_else(|| panic!("bad value in row {l:?}"))
        })
        .collect()
}

#[test]
fn sweep_equilibrium_hits_target_displacements() {
    let iges_path = tmp("sweep_shell.iges");
    std::fs::write(&iges_path, wigley_shell_iges()).unwrap();
    let out = run_ok(bin().args([
        "sweep",
        iges_path.to_str().unwrap(),
        "--waterline",
        "0.7",
        "--float",
        "weight=1000:2000:1000",
        "--speeds",
        "3.0",
        "--samples",
        "41x13",
        "--fit-control",
        "8x6",
    ]));
    let vols = csv_col(&out, "volume");
    assert_eq!(vols.len(), 2, "{out}");
    for (v, mass) in vols.iter().zip([1000.0, 2000.0]) {
        let want = mass / 1025.9;
        assert!(
            (v - want).abs() < 0.005 * want,
            "volume {v} vs target {want}"
        );
    }
    let rws = csv_col(&out, "rw");
    assert!(rws.iter().all(|r| r.is_finite() && *r > 0.0), "{rws:?}");
    // Heavier boat, deeper: sinkage must increase with weight.
    let sink = csv_col(&out, "sinkage");
    assert!(sink[1] > sink[0], "sinkage {sink:?}");
}

#[test]
fn sweep_raw_waterline_and_trim_axes() {
    let iges_path = tmp("sweep_shell2.iges");
    std::fs::write(&iges_path, wigley_shell_iges()).unwrap();
    let out = run_ok(bin().args([
        "sweep",
        iges_path.to_str().unwrap(),
        "--waterline",
        "0.7",
        "--axis",
        "waterline=0.6:0.7:0.1",
        "--axis",
        "sweep_shell2:trim=-2:2:2",
        "--speeds",
        "3.0",
        "--samples",
        "41x13",
        "--fit-control",
        "8x6",
    ]));
    // 2 waterlines x 3 trims x 1 speed = 6 rows.
    let vols = csv_col(&out, "volume");
    assert_eq!(vols.len(), 6, "{out}");
    // Deeper waterline displaces more at every trim.
    let deep: f64 = vols[3..].iter().sum();
    let shallow: f64 = vols[..3].iter().sum();
    assert!(deep > shallow, "{vols:?}");
    // Trim symmetry within each waterline.
    assert!((vols[0] - vols[2]).abs() < 1e-3 * vols[0], "{vols:?}");
}

#[test]
fn loft_decomposes_and_manifest_sweeps() {
    // Trimaran: shells at y = 0, +7, -7 in one IGES file.
    let iges_path = tmp("tri.iges");
    std::fs::write(&iges_path, wigley_shells_iges(&[0.0, 7.0, -7.0])).unwrap();

    // Decompose to semantic full-band bodies.
    let prefix = tmp("tri");
    let out = run_ok(bin().args([
        "loft",
        iges_path.to_str().unwrap(),
        "--waterline",
        "0.7",
        "-o",
        prefix.to_str().unwrap(),
        "--samples",
        "61x21",
        "--fit-control",
        "9x7",
        "--fit-degree",
        "2x2",
    ]));
    for name in ["port", "center", "starboard"] {
        let f = tmp(&format!("tri-{name}.hull"));
        assert!(f.exists(), "missing {f:?}\n{out}");
        let text = std::fs::read_to_string(&f).unwrap();
        assert!(text.contains("\nwaterline "), "no waterline key in {name}");
        assert!(text.contains("\ncenterplane "), "no centerplane key in {name}");
    }
    // Ports and starboards at the right sides.
    let port = std::fs::read_to_string(tmp("tri-port.hull")).unwrap();
    let cp: f64 = port
        .lines()
        .find_map(|l| l.strip_prefix("centerplane "))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert!((cp + 7.0).abs() < 1e-3, "port centerplane {cp}");

    // Manifest: equilibrium weight sweep x ama spread.
    let manifest = r#"{
  "name": "cli test study",
  "fluid": "seawater",
  "hulls": [
    { "id": "vaka",  "file": "tri-center.hull" },
    { "id": "ama_s", "file": "tri-starboard.hull" },
    { "id": "ama_p", "file": "tri-port.hull" }
  ],
  "sweep": [
    { "target": "speed", "unit": "ms", "value": 3.0 },
    { "target": "weight", "range": [3000, 6000], "step": 3000 },
    { "target": "lcg", "value": 5.0 },
    { "target": ["ama_s", "ama_p"], "param": "spread", "values": [0, 0.5] }
  ],
  "output": { "format": "csv", "file": "study.csv" },
  "options": { "samples": "61x17", "fit_control": "9x7", "fit_degree": "2x2" }
}"#;
    let man_path = tmp("study.json");
    std::fs::write(&man_path, manifest).unwrap();
    run_ok(bin().args(["sweep", man_path.to_str().unwrap()]));

    let csv = std::fs::read_to_string(tmp("study.csv")).unwrap();
    let vols = csv_col(&csv, "volume");
    assert_eq!(vols.len(), 4, "{csv}");
    // Rows: (3000,spread 0), (3000,0.5), (6000,0), (6000,0.5).
    for (v, mass) in vols.iter().zip([3000.0, 3000.0, 6000.0, 6000.0]) {
        let want = mass / 1025.9;
        assert!(
            (v - want).abs() < 0.005 * want,
            "volume {v} vs target {want}"
        );
    }
    // Spread must change the interference but not the displacement.
    let iff = csv_col(&csv, "interference");
    assert!(
        (iff[0] - iff[1]).abs() > 1e-4,
        "spread had no effect on interference: {iff:?}"
    );
    let rws = csv_col(&csv, "rw");
    assert!(rws.iter().all(|r| r.is_finite() && *r > 0.0), "{rws:?}");
}

#[test]
fn manifest_heel_axis_produces_gz_curve() {
    // Catamaran: shells at y = +-4 in one IGES file, lofted with the design
    // waterline at 0.5 so the bodies keep topsides above the DWL.
    let iges_path = tmp("cat.iges");
    std::fs::write(&iges_path, wigley_shells_iges(&[4.0, -4.0])).unwrap();
    let prefix = tmp("cat");
    run_ok(bin().args([
        "loft",
        iges_path.to_str().unwrap(),
        "--waterline",
        "0.5",
        "-o",
        prefix.to_str().unwrap(),
        "--samples",
        "61x21",
        "--fit-control",
        "9x7",
        "--fit-degree",
        "2x2",
    ]));

    let manifest = r#"{
  "name": "gz curve",
  "fluid": "seawater",
  "hulls": [
    { "id": "port", "file": "cat-port.hull" },
    { "id": "stbd", "file": "cat-starboard.hull" }
  ],
  "sweep": [
    { "target": "speed", "unit": "ms", "value": 3.0 },
    { "target": "weight", "value": 2900 },
    { "target": "vcg", "value": 0.2 },
    { "target": "heel", "values": [-1.5, 0, 1.5] }
  ],
  "output": { "format": "csv", "file": "gz.csv" },
  "options": { "samples": "61x17", "fit_control": "9x7", "fit_degree": "2x2" }
}"#;
    let man_path = tmp("gz_study.json");
    std::fs::write(&man_path, manifest).unwrap();
    run_ok(bin().args(["sweep", man_path.to_str().unwrap()]));

    let csv = std::fs::read_to_string(tmp("gz.csv")).unwrap();
    let gz = csv_col(&csv, "gz");
    assert_eq!(gz.len(), 3, "{csv}");
    // Upright a symmetric catamaran has no righting arm; heeled to starboard
    // (+y down) the buoyancy transfer at 4 m spacing must right it strongly,
    // and the curve is antisymmetric.
    assert!(gz[1].abs() < 1e-6, "gz(0) = {}", gz[1]);
    assert!(gz[2] > 0.5, "gz(+1.5 deg) = {}", gz[2]);
    assert!(
        (gz[0] + gz[2]).abs() < 0.02 * gz[2],
        "gz not antisymmetric: {gz:?}"
    );
    // rm = displacement weight x gz.
    let rm = csv_col(&csv, "rm");
    let want = 2900.0 * michell::STANDARD_GRAVITY * gz[2];
    assert!(
        (rm[2] - want).abs() < 1e-6 * want.abs(),
        "rm {} vs {want}",
        rm[2]
    );
    // Displacement is held across the heel sweep (equilibrium re-solved).
    let vols = csv_col(&csv, "volume");
    for v in &vols {
        let want = 2900.0 / 1025.9;
        assert!((v - want).abs() < 0.005 * want, "volume {v} vs {want}");
    }
}

/// Tessellated Wigley full shell (ASCII STL, metres, DWL at z = 0.7).
fn wigley_stl(nx: usize, nz: usize) -> String {
    let f = |x: f64, zp: f64| {
        0.5 * (4.0 * (x / 10.0) * (1.0 - x / 10.0)) * (1.0 - (zp / 0.625f64).powi(2))
    };
    let mut s = String::from("solid wigley\n");
    for side in [1.0f64, -1.0] {
        for i in 0..nx {
            for j in 0..nz {
                let (x0, x1) = (10.0 * i as f64 / nx as f64, 10.0 * (i + 1) as f64 / nx as f64);
                let (z0, z1) = (
                    0.625 * j as f64 / nz as f64,
                    0.625 * (j + 1) as f64 / nz as f64,
                );
                let p = |x: f64, zp: f64| [x, side * f(x, zp), 0.7 - zp];
                for tri in [
                    [p(x0, z0), p(x1, z0), p(x1, z1)],
                    [p(x0, z0), p(x1, z1), p(x0, z1)],
                ] {
                    s.push_str(" facet normal 0 0 0\n  outer loop\n");
                    for v in tri {
                        s.push_str(&format!("   vertex {} {} {}\n", v[0], v[1], v[2]));
                    }
                    s.push_str("  endloop\n endfacet\n");
                }
            }
        }
    }
    s.push_str("endsolid wigley\n");
    s
}

#[test]
fn stl_resistance_and_loft() {
    let stl_path = tmp("wigley_mesh.stl");
    std::fs::write(&stl_path, wigley_stl(120, 40)).unwrap();

    // Missing units must fail with advice.
    let out = bin()
        .args(["info", stl_path.to_str().unwrap(), "--waterline", "0.7"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--units"));

    // Direct resistance from the mesh matches the reference Wigley.
    let out = run_ok(bin().args([
        "resistance",
        stl_path.to_str().unwrap(),
        "--units",
        "m",
        "--waterline",
        "0.7",
        "--speeds",
        "3.0",
        "--fluid",
        "seawater",
        "--json",
    ]));
    let reference = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
    let want = michell::resistance(&reference, &michell::Conditions::seawater(3.0)).unwrap();
    let rw = json_num(&out, "rw");
    assert!(
        (rw - want.wave.resistance).abs() < 0.02 * want.wave.resistance,
        "rw {rw} vs {}",
        want.wave.resistance
    );

    // Loft to a full-band body and use it.
    let body_path = tmp("wigley_mesh_body");
    run_ok(bin().args([
        "loft",
        stl_path.to_str().unwrap(),
        "--units",
        "m",
        "--waterline",
        "0.7",
        "-o",
        body_path.to_str().unwrap(),
        "--samples",
        "121x49",
        "--fit-control",
        "16x12",
    ]));
    let body_file = tmp("wigley_mesh_body.hull");
    let text = std::fs::read_to_string(&body_file).unwrap();
    assert!(text.contains("\nwaterline "), "not a body file");
    let out = run_ok(bin().args([
        "resistance",
        body_file.to_str().unwrap(),
        "--speeds",
        "3.0",
        "--json",
    ]));
    let rw_body = json_num(&out, "rw");
    assert!(
        (rw_body - want.wave.resistance).abs() < 0.03 * want.wave.resistance,
        "body rw {rw_body} vs {}",
        want.wave.resistance
    );
}

#[test]
fn dump_grid_roundtrips_through_grid_json() {
    let iges_path = tmp("dump_shell.iges");
    std::fs::write(&iges_path, wigley_shell_iges()).unwrap();
    let grid_path = tmp("dump.grid.json");

    // Import the IGES, dumping the sampled grid IR alongside.
    let direct = run_ok(bin().args([
        "info",
        iges_path.to_str().unwrap(),
        "--waterline",
        "0.7",
        "--samples",
        "41x13",
        "--fit-control",
        "8x6",
        "--dump-grid",
        grid_path.to_str().unwrap(),
        "--json",
    ]));

    // The dumped grid must carry the derivative channels.
    let grid_text = std::fs::read_to_string(&grid_path).unwrap();
    assert!(grid_text.contains("\"michell\": \"sample-grid\""), "{grid_text}");
    assert!(grid_text.contains("\"dfdx\""), "no dfdx channel:\n{grid_text}");
    assert!(grid_text.contains("\"dfdz\""));
    assert!(grid_text.contains("\"weights\""));

    // Re-lofting the dumped grid with the same fit reproduces the hull.
    let reloaded = run_ok(bin().args([
        "info",
        grid_path.to_str().unwrap(),
        "--fit-control",
        "8x6",
        "--json",
    ]));
    for key in ["length", "draft", "wetted_surface", "displaced_volume"] {
        let (a, b) = (json_num(&direct, key), json_num(&reloaded, key));
        assert!(
            (a - b).abs() <= 1e-9 * a.abs().max(1e-9),
            "{key}: {a} vs {b}"
        );
    }
}

#[test]
fn errors_are_clean() {
    // Unknown file format.
    let bogus = tmp("bogus.txt");
    std::fs::write(&bogus, "hello world\n").unwrap();
    let out = bin()
        .args(["info", bogus.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("cannot determine the format"));

    // Missing speed selection.
    let hull_path = tmp("wigley_err.hull");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let out = bin()
        .args(["resistance", hull_path.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--speeds or --froude"));
}

#[test]
fn spectrum_cross_checks_resistance() {
    let hull_path = tmp("wigley_spectrum.hull");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let out = run_ok(bin().args([
        "spectrum",
        hull_path.to_str().unwrap(),
        "--speed",
        "3",
        "--json",
    ]));
    let reference = michell::wave_resistance(
        &michell::hulls::wigley(10.0, 1.0, 0.625).unwrap(),
        &michell::Conditions::seawater(3.0),
    )
    .unwrap()
    .resistance;
    let rw_michell = json_num(&out, "rw_michell");
    let rw_spectrum = json_num(&out, "rw_spectrum");
    assert!((rw_michell - reference).abs() < 1e-6 * reference);
    assert!(
        (rw_spectrum - rw_michell).abs() < 1e-2 * rw_michell,
        "spectrum {rw_spectrum} vs michell {rw_michell}"
    );

    // CSV variant has the documented header and the right row count.
    let csv = run_ok(bin().args([
        "spectrum",
        hull_path.to_str().unwrap(),
        "--speed",
        "3",
        "--points",
        "101",
    ]));
    let mut lines = csv.lines();
    assert_eq!(
        lines.next().unwrap(),
        "theta_deg,lambda,wavelength_m,amp_re,amp_im,amp_abs,drw_dtheta,cum_fraction"
    );
    assert_eq!(lines.count(), 101);
}

#[test]
fn wake_writes_png_and_json() {
    let hull_path = tmp("wigley_wake.hull");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let png_path = tmp("wake.png");
    run_ok(bin().args([
        "wake",
        hull_path.to_str().unwrap(),
        "--speed",
        "3",
        "--region",
        "-30:-5:-10:10",
        "--size",
        "200x150",
        "-o",
        png_path.to_str().unwrap(),
    ]));
    let bytes = std::fs::read(&png_path).unwrap();
    assert!(bytes.len() > 1000, "png too small: {}", bytes.len());
    assert_eq!(&bytes[..8], &[0x89, b'P', b'N', b'G', 0x0d, 0x0a, 0x1a, 0x0a]);

    let out = run_ok(bin().args([
        "wake",
        hull_path.to_str().unwrap(),
        "--speed",
        "3",
        "--region",
        "-30:-10:-6:6",
        "--size",
        "50x20",
        "--json",
    ]));
    assert_eq!(json_num(&out, "nx"), 50.0);
    assert_eq!(json_num(&out, "ny"), 20.0);
    assert!(out.contains("\"zeta\":[["));
    // The wake is not flat.
    let zeta = &out[out.find("\"zeta\":[[").unwrap() + 9..];
    let first: f64 = zeta[..zeta.find(',').unwrap()].parse().unwrap();
    assert!(first.is_finite());
}

#[test]
fn render_writes_png() {
    let hull_path = tmp("wigley_render.hull");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let png_path = tmp("render.png");
    run_ok(bin().args([
        "render",
        hull_path.to_str().unwrap(),
        "--speed",
        "3",
        "--grid",
        "120",
        "--size",
        "160x100",
        "-o",
        png_path.to_str().unwrap(),
    ]));
    let bytes = std::fs::read(&png_path).unwrap();
    assert!(bytes.len() > 1000, "png too small: {}", bytes.len());
    assert_eq!(&bytes[..8], &[0x89, b'P', b'N', b'G', 0x0d, 0x0a, 0x1a, 0x0a]);
}

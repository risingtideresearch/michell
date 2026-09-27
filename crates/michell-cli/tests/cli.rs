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
    rest[..end]
        .parse::<f64>()
        .unwrap_or_else(|_| panic!("parse {key}: {rest:?}"))
}

fn run_err(cmd: &mut Command) -> String {
    let out = cmd.output().expect("binary runs");
    assert!(!out.status.success(), "command unexpectedly succeeded");
    String::from_utf8_lossy(&out.stderr).into_owned()
}

/// A wigley control net to hang the viscous tests off.
/// The exact Wigley cut into sections through the public library API, as
/// the CLI cuts `michell wigley`'s IGES.
fn sectional_wigley() -> michell_geometry::SectionalHull {
    use michell_geometry::iges::{self, HullPose, Platform, SectionalOptions};
    let surfaces = iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
    let src = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
    src.situate_sectional(
        0,
        0.0,
        &HullPose::default(),
        &Platform::default(),
        &SectionalOptions::default(),
    )
    .unwrap()
    .unwrap()
    .hull
}

fn library_resistance(
    members: &[(&michell_geometry::SectionalHull, michell::Placement)],
    cond: &michell::Conditions,
) -> michell::MultihullResistance {
    michell::sectional::multihull_resistance(
        members,
        cond,
        &michell::WaveOptions::default(),
        &michell::ViscousOptions::default(),
    )
    .unwrap()
}

fn wigley_hull(name: &str) -> PathBuf {
    let path = tmp(name);
    run_ok(bin().args([
        "wigley",
        "--length",
        "10",
        "--beam",
        "1",
        "--draft",
        "0.625",
        "-o",
        path.to_str().unwrap(),
    ]));
    path
}

#[test]
fn roughness_is_separate_from_the_form_factor() {
    // C_V = (1+k)·C_F + ΔC_F. The two knobs must be independent, and the
    // roughness must sit *outside* the form factor: doubling k must not
    // change the roughness share of the total.
    let hull = wigley_hull("visc.igs");
    let h = hull.to_str().unwrap();
    let rv = |args: &[&str]| -> f64 {
        let mut c = bin();
        c.args(["resistance", h, "--speeds", "3", "--json"]);
        c.args(args);
        json_num(&run_ok(&mut c), "rv")
    };
    let plain = rv(&[]);
    let k_only = rv(&["--form-factor", "0.2"]);
    let cf_only = rv(&["--roughness", "cf=4e-4"]);
    let both = rv(&["--form-factor", "0.2", "--roughness", "cf=4e-4"]);

    assert!(k_only > plain && cf_only > plain);
    // Additivity: the ΔC_F contribution is the same whatever k is.
    let d0 = cf_only - plain;
    let dk = both - k_only;
    assert!(
        (d0 - dk).abs() < 1e-9 * d0,
        "roughness is not outside the form factor: {d0} vs {dk}"
    );
    // (1+k) multiplies only the friction part.
    assert!((k_only - plain - 0.2 * plain).abs() < 1e-9 * plain);

    // Reported, not silently folded in.
    let json = run_ok(bin().args([
        "resistance",
        h,
        "--speeds",
        "3",
        "--json",
        "--roughness",
        "cf=4e-4",
    ]));
    assert!((json_num(&json, "roughness_cf") - 4e-4).abs() < 1e-12);
    assert!(json.contains("\"delta_cf\":0.0004"));
    // Default is smooth and says so.
    let plain_json = run_ok(bin().args(["resistance", h, "--speeds", "3", "--json"]));
    assert_eq!(json_num(&plain_json, "roughness_cf"), 0.0);
    assert!(plain_json.contains("\"roughness\":null"));
}

#[test]
fn sand_grain_roughness_is_smooth_until_it_is_not() {
    // A finish inside the viscous sublayer costs exactly nothing; a coarse
    // one costs something; and the penalty grows with speed, because the
    // fully-rough branch is Re-independent while the smooth line falls.
    let hull = wigley_hull("visc2.igs");
    let h = hull.to_str().unwrap();
    let rv = |ks: &str, u: &str| -> f64 {
        let mut c = bin();
        c.args(["resistance", h, "--speeds", u, "--json"]);
        if !ks.is_empty() {
            c.args(["--roughness", ks]);
        }
        json_num(&run_ok(&mut c), "rv")
    };
    assert_eq!(rv("ks=1um", "3"), rv("", "3"));
    assert!(rv("ks=1mm", "3") > rv("", "3"));

    let excess = |u: &str| (rv("ks=1mm", u) - rv("", u)) / rv("", u);
    assert!(
        excess("6") > excess("2"),
        "roughness share did not grow with speed"
    );

    // Unit suffixes are equivalent ways of saying the same height.
    assert_eq!(rv("ks=1mm", "3"), rv("ks=1000um", "3"));
    assert_eq!(rv("ks=1mm", "3"), rv("ks=0.001", "3"));

    // The regime is reported, not left to the docs.
    let text = run_ok(bin().args(["resistance", h, "--speeds", "3", "--roughness", "ks=1mm"]));
    assert!(text.contains("k_s+"), "no regime note in:\n{text}");
}

#[test]
fn rejects_bad_roughness_specs() {
    let hull = wigley_hull("visc3.igs");
    let h = hull.to_str().unwrap();
    for spec in ["ks=-1um", "cf=-1e-4", "banana", "ks=", "cf=abc"] {
        let err = run_err(bin().args(["resistance", h, "--speeds", "3", "--roughness", spec]));
        assert!(
            err.contains("roughness") || err.contains("length"),
            "spec {spec:?} gave an unhelpful error: {err}"
        );
    }
    // `off` is accepted and means smooth.
    let a = run_ok(bin().args([
        "resistance",
        h,
        "--speeds",
        "3",
        "--json",
        "--roughness",
        "off",
    ]));
    let b = run_ok(bin().args(["resistance", h, "--speeds", "3", "--json"]));
    assert_eq!(json_num(&a, "rv"), json_num(&b, "rv"));
}

#[test]
fn wigley_roundtrip_matches_library() {
    let hull_path = tmp("wigley.igs");
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
    let reference = sectional_wigley();
    assert!((json_num(&info, "length") - 10.0).abs() < 1e-9);
    assert!((json_num(&info, "draft") - 0.625).abs() < 1e-9);
    let vol = reference.displaced_volume();
    assert!((json_num(&info, "displaced_volume") - vol).abs() < 1e-8 * vol);

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
    let want = library_resistance(&[(&reference, michell::Placement::default())], &cond);
    let rw = json_num(&out, "rw");
    let rv = json_num(&out, "rv");
    assert!(
        (rw - want.wave.resistance).abs() < 1e-6 * want.wave.resistance,
        "rw {rw} vs {}",
        want.wave.resistance
    );
    // Viscous: the wetted surface is the one the CLI reports, within its
    // strip construction of the exact Wigley shell 2∬√(1 + f_x² + f_z²).
    let (n, (l, b, t)) = (800, (10.0f64, 1.0f64, 0.625f64));
    let mut exact = 0.0;
    for i in 0..n {
        for j in 0..n {
            let (u, v) = ((i as f64 + 0.5) / n as f64, (j as f64 + 0.5) / n as f64);
            let (x, z) = (-l / 2.0 + l * u, t * v);
            let fx = b / 2.0 * (-8.0 * x / (l * l)) * (1.0 - (z / t).powi(2));
            let fz = b / 2.0 * (1.0 - (2.0 * x / l).powi(2)) * (-2.0 * z / (t * t));
            exact += 2.0 * (1.0 + fx * fx + fz * fz).sqrt() * (l / n as f64) * (t / n as f64);
        }
    }
    let s_cli = json_num(&info, "wetted_surface");
    assert!((s_cli - exact).abs() < 2e-3 * exact, "S {s_cli} vs {exact}");
    let rv_want = michell::viscous_resistance_for(10.0, s_cli, &cond, &Default::default())
        .unwrap()
        .resistance;
    assert!((rv - rv_want).abs() < 1e-9 * rv_want, "rv {rv} vs {rv_want}");
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
    let hull_path = tmp("wigley_table.igs");
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
fn offsets_and_loft_are_rejected() {
    let off_path = tmp("wigley.offsets");
    std::fs::write(&off_path, "michell-offsets v1\nwaterlines 0 0.5\nstation 0 1 0\n").unwrap();
    let err = run_err(bin().args(["info", off_path.to_str().unwrap()]));
    assert!(err.contains("no longer read"), "{err}");
    let err = run_err(bin().args(["loft", off_path.to_str().unwrap(), "-o", "x"]));
    assert!(err.contains("was removed"), "{err}");
}

#[test]
fn knots_flag_scales_speeds() {
    let hull_path = tmp("wigley_kn.igs");
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
    assert!((json_num(&kn, "rw") - json_num(&ms, "rw")).abs() < 1e-9 * json_num(&ms, "rw").abs());
}

#[test]
fn catamaran_fleet_matches_library() {
    let hull_path = tmp("wigley_cat.igs");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let spec_a = format!("{}@y=1.0", hull_path.to_str().unwrap());
    let spec_b = format!("{}@y=-1.0", hull_path.to_str().unwrap());
    let out = run_ok(bin().args(["resistance", &spec_a, &spec_b, "--speeds", "3.0", "--json"]));
    let hull = sectional_wigley();
    let members = [
        (&hull, michell::Placement { x: 0.0, y: 1.0 }),
        (&hull, michell::Placement { x: 0.0, y: -1.0 }),
    ];
    let want = library_resistance(&members, &michell::Conditions::seawater(3.0));
    let rw = json_num(&out, "rw");
    let iff = json_num(&out, "interference");
    assert!(
        (rw - want.wave.resistance).abs() < 1e-6 * want.wave.resistance,
        "rw {rw} vs {}",
        want.wave.resistance
    );
    assert!(
        (iff - want.interference).abs() < 1e-6,
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
        "--stations",
        "41",
        "--rays",
        "17",
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
        "--stations",
        "41",
        "--rays",
        "17",
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
fn manifest_sweeps_iges_hulls() {
    // Trimaran: shells at y = 0, +7, -7 in one IGES file.
    let iges_path = tmp("tri.iges");
    std::fs::write(&iges_path, wigley_shells_iges(&[0.0, 7.0, -7.0])).unwrap();

    // The file's hulls, ordered by transverse position.
    let info = run_ok(bin().args([
        "info",
        iges_path.to_str().unwrap(),
        "--waterline",
        "0.7",
        "--json",
    ]));
    let ys: Vec<f64> = info
        .split("\"y\":")
        .skip(1)
        .map(|t| t.split(['}', ',']).next().unwrap().parse().unwrap())
        .collect();
    assert_eq!(ys.len(), 3, "{info}");
    for (y, want) in ys.iter().zip([-7.0, 0.0, 7.0]) {
        assert!((y - want).abs() < 1e-3, "centerplanes {ys:?}");
    }

    // Manifest: equilibrium weight sweep x ama spread, each hull picked out
    // of the one file.
    let manifest = r#"{
  "name": "cli test study",
  "fluid": "seawater",
  "hulls": [
    { "id": "vaka",  "file": "tri.iges", "hull": 1, "load": { "mass": 3000, "lcg": 5.0 } },
    { "id": "ama_s", "file": "tri.iges", "hull": 2 },
    { "id": "ama_p", "file": "tri.iges", "hull": 0 }
  ],
  "sweep": [
    { "target": "speed", "unit": "ms", "value": 3.0 },
    { "target": "vaka", "param": "mass", "values": [0, 3000] },
    { "target": ["ama_s", "ama_p"], "param": "spread", "values": [0, 0.5] }
  ],
  "output": { "format": "csv", "file": "study.csv" },
  "options": { "waterline": 0.7, "stations": 41, "rays": 17 }
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

    // A multihull file needs its hull picked.
    let bad = manifest.replace(r#", "hull": 1"#, "");
    std::fs::write(&man_path, bad).unwrap();
    let err = run_err(bin().args(["sweep", man_path.to_str().unwrap()]));
    assert!(err.contains("pick one with"), "{err}");
}

/// Heel was removed. Old inputs that still ask for it — the `--heel` flag, a
/// manifest `options.heel` block, a point load's transverse `dy` (base value or
/// swept axis) — fail loudly rather than being silently ignored.
#[test]
fn heel_inputs_are_rejected() {
    let hull = wigley_hull("noheel.igs");
    let err = run_err(bin().args([
        "resistance",
        hull.to_str().unwrap(),
        "--speed",
        "3",
        "--heel",
        "10",
    ]));
    assert!(err.contains("--heel was removed"), "{err}");

    let iges_path = tmp("noheel.iges");
    std::fs::write(&iges_path, wigley_shells_iges(&[0.0])).unwrap();
    let cases = [
        (
            r#"{ "id": "vaka", "file": "noheel.iges", "load": { "mass": 1500 } }"#,
            r#"{ "target": "speed", "unit": "ms", "value": 3.0 }"#,
            r#", "heel": { "gz_step": 5 }"#,
            "heel was removed",
        ),
        (
            r#"{ "id": "vaka", "file": "noheel.iges", "load": { "mass": 1500 },
                 "points": [ { "id": "crew", "mass": 80, "dy": 0.5 } ] }"#,
            r#"{ "target": "speed", "unit": "ms", "value": 3.0 }"#,
            "",
            "\"dy\" was removed",
        ),
        (
            r#"{ "id": "vaka", "file": "noheel.iges", "load": { "mass": 1500 },
                 "points": [ { "id": "crew", "mass": 80 } ] }"#,
            r#"{ "target": "speed", "unit": "ms", "value": 3.0 },
               { "target": "crew", "param": "dy", "values": [0.0, 0.5] }"#,
            "",
            "\"dy\" was removed",
        ),
    ];
    for (k, (hull, sweep, opts, want)) in cases.iter().enumerate() {
        let manifest = format!(
            r#"{{ "name": "no heel", "fluid": "seawater", "hulls": [ {hull} ],
  "sweep": [ {sweep} ], "output": {{ "format": "csv", "file": "noheel.csv" }},
  "options": {{ "waterline": 0.5{opts} }} }}"#
        );
        let man_path = tmp(&format!("noheel{k}.json"));
        std::fs::write(&man_path, manifest).unwrap();
        let err = run_err(bin().args(["sweep", man_path.to_str().unwrap()]));
        assert!(err.contains(want), "case {k}: {err}");
    }
}

// --- Minimal `.msw` archive reader for the binary-output test. Mirrors the
// format documented in `archive.rs` (the module itself is private to the
// binary crate, so integration tests re-implement the walk). ---

struct Blob {
    kind: u32,
    name: String,
    data: Vec<u8>,
}

fn parse_msw(bytes: &[u8]) -> Vec<Blob> {
    assert_eq!(&bytes[0..4], b"MSWP", "bad magic");
    assert_eq!(
        u32::from_le_bytes(bytes[4..8].try_into().unwrap()),
        2,
        "version"
    );
    let mut blobs = Vec::new();
    let mut p = 8usize;
    while p < bytes.len() {
        let kind = u32::from_le_bytes(bytes[p..p + 4].try_into().unwrap());
        p += 4;
        let nl = u32::from_le_bytes(bytes[p..p + 4].try_into().unwrap()) as usize;
        p += 4;
        let name = String::from_utf8(bytes[p..p + nl].to_vec()).unwrap();
        p += nl;
        let dl = u64::from_le_bytes(bytes[p..p + 8].try_into().unwrap()) as usize;
        p += 8;
        let data = bytes[p..p + dl].to_vec();
        p += dl;
        blobs.push(Blob { kind, name, data });
    }
    blobs
}

/// The binary output format bundles the manifest, the referenced hull files,
/// and per-row parameters + metrics + spectrum into one `.msw`
/// file that round-trips without re-running the study.
#[test]
fn manifest_binary_archive_bundles_everything() {
    let iges_path = tmp("bin_arc.iges");
    std::fs::write(&iges_path, wigley_shells_iges(&[0.0])).unwrap();
    let hull_bytes = std::fs::read(&iges_path).unwrap();

    let manifest = r#"{
  "name": "binary archive test",
  "fluid": "seawater",
  "hulls": [ { "id": "vaka", "file": "bin_arc.iges", "load": { "mass": 1500, "vcg": 0.0 } } ],
  "sweep": [ { "target": "speed", "unit": "ms", "values": [2.5, 3.5] } ],
  "output": { "format": "binary", "file": "study.msw", "spectrum": { "points": 129 } },
  "options": { "waterline": 0.5, "stations": 41, "rays": 17 }
}"#;
    let man_path = tmp("bin_arc.json");
    std::fs::write(&man_path, manifest).unwrap();
    run_ok(bin().args(["sweep", man_path.to_str().unwrap()]));

    let bytes = std::fs::read(tmp("study.msw")).unwrap();
    let blobs = parse_msw(&bytes);

    // Manifest blob preserved verbatim.
    let man = blobs.iter().find(|b| b.kind == 1).expect("manifest blob");
    assert_eq!(man.name, "bin_arc.json");
    assert_eq!(man.data, manifest.as_bytes());

    // Hull file bundled byte-for-byte.
    let hull = blobs.iter().find(|b| b.kind == 2).expect("hull blob");
    assert_eq!(hull.name, "bin_arc.iges");
    assert_eq!(hull.data, hull_bytes);

    // Meta names the columns.
    let meta = blobs.iter().find(|b| b.kind == 3).expect("meta blob");
    let meta_txt = String::from_utf8(meta.data.clone()).unwrap();
    assert!(meta_txt.contains("\"metric_labels\""), "{meta_txt}");
    assert!(meta_txt.contains("\"vcg\""), "{meta_txt}");
    assert!(
        !meta_txt.contains("gz"),
        "no GZ in the archive any more: {meta_txt}"
    );
    assert!(meta_txt.contains("\"speeds_ms\":[2.5,3.5]"), "{meta_txt}");

    // Rows: two speeds → two rows, each with a spectrum.
    let rows = blobs.iter().find(|b| b.kind == 4).expect("rows blob");
    let d = &rows.data;
    let mut p = 0usize;
    let u32_at = |p: &mut usize| {
        let v = u32::from_le_bytes(d[*p..*p + 4].try_into().unwrap());
        *p += 4;
        v
    };
    let n_rows = u32_at(&mut p);
    let n_axes = u32_at(&mut p);
    let n_metrics = u32_at(&mut p);
    assert_eq!(n_rows, 2, "one row per speed");
    assert_eq!(n_axes, 0, "no pose/load axes in this study");
    assert!(n_metrics > 10, "metrics present: {n_metrics}");

    // Walk row 0 and confirm the spectrum is non-empty.
    p += (n_axes as usize + n_metrics as usize) * 8;
    p += 16; // wavenumber + transverse wavelength
    let spec_n = u32_at(&mut p);
    assert_eq!(spec_n, 129, "spectrum sampled at the requested resolution");
}

/// A point load mounted on a hull adds to the derived fleet CG and rides with
/// its own swept offset: lowering it (`dz` +down) pulls `vcg` down, and its
/// mass shows up in the derived `mass` column — all pure CG arithmetic, so the
/// derived columns are exact.
#[test]
fn manifest_point_load_moves_derived_cg() {
    let iges_path = tmp("ptload.iges");
    std::fs::write(&iges_path, wigley_shells_iges(&[0.0])).unwrap();

    let manifest = r#"{
  "name": "point load cg",
  "fluid": "seawater",
  "hulls": [
    { "id": "vaka", "file": "ptload.iges",
      "load": { "mass": 1000, "vcg": 0.5 },
      "points": [ { "id": "keel", "mass": 500, "dz": 1.0 } ] }
  ],
  "sweep": [
    { "target": "speed", "unit": "ms", "value": 3.0 },
    { "target": "keel", "param": "dz", "values": [0.0, 2.0] }
  ],
  "output": { "format": "csv", "file": "ptload.csv" },
  "options": { "waterline": 0.5, "stations": 41, "rays": 17 }
}"#;
    let man_path = tmp("ptload.json");
    std::fs::write(&man_path, manifest).unwrap();
    run_ok(bin().args(["sweep", man_path.to_str().unwrap()]));

    let csv = std::fs::read_to_string(tmp("ptload.csv")).unwrap();
    // Derived fleet mass = structural 1000 + point 500, at both points.
    let mass = csv_col(&csv, "mass");
    assert_eq!(mass.len(), 2, "{csv}");
    assert!(mass.iter().all(|m| (m - 1500.0).abs() < 1e-6), "{mass:?}");
    // vcg = (1000·0.5 + 500·(−dz)) / 1500. dz base 1.0: row0 dz=1.0 → 0;
    // row1 dz=3.0 → (500 − 1500)/1500 = −2/3.
    let vcg = csv_col(&csv, "vcg");
    assert!(vcg[0].abs() < 1e-6, "vcg(dz=1) = {}", vcg[0]);
    assert!((vcg[1] + 2.0 / 3.0).abs() < 1e-6, "vcg(dz=3) = {}", vcg[1]);
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
                let (x0, x1) = (
                    10.0 * i as f64 / nx as f64,
                    10.0 * (i + 1) as f64 / nx as f64,
                );
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
fn stl_resistance() {
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
    let reference = sectional_wigley();
    let want = library_resistance(
        &[(&reference, michell::Placement::default())],
        &michell::Conditions::seawater(3.0),
    );
    let rw = json_num(&out, "rw");
    assert!(
        (rw - want.wave.resistance).abs() < 0.02 * want.wave.resistance,
        "rw {rw} vs {}",
        want.wave.resistance
    );
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
    let hull_path = tmp("wigley_err.igs");
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
    let hull_path = tmp("wigley_spectrum.igs");
    run_ok(bin().args(["wigley", "-o", hull_path.to_str().unwrap()]));
    let out = run_ok(bin().args([
        "spectrum",
        hull_path.to_str().unwrap(),
        "--speed",
        "3",
        "--json",
    ]));
    let reference = library_resistance(
        &[(&sectional_wigley(), michell::Placement::default())],
        &michell::Conditions::seawater(3.0),
    )
    .wave
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
    let hull_path = tmp("wigley_wake.igs");
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
    assert_eq!(
        &bytes[..8],
        &[0x89, b'P', b'N', b'G', 0x0d, 0x0a, 0x1a, 0x0a]
    );

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
    let hull_path = tmp("wigley_render.igs");
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
    assert_eq!(
        &bytes[..8],
        &[0x89, b'P', b'N', b'G', 0x0d, 0x0a, 0x1a, 0x0a]
    );
}

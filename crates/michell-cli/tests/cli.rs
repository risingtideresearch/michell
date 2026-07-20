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

/// Full-shell Wigley (L=10, B=1, T0=0.625) in metres, z up, DWL at z=0.7,
/// centred at y=0: 4 untrimmed biquadratic patches.
fn wigley_shell_iges() -> String {
    fn line(content: &str, section: char, seq: usize) -> String {
        format!("{content:<72}{section}{seq:>7}\n")
    }
    let mut bodies: Vec<String> = Vec::new();
    let (gx_f, gx_a) = ([0.0, 1.0, 1.0], [1.0, 1.0, 0.0]);
    let (xs_f, xs_a) = ([0.0, 2.5, 5.0], [5.0, 7.5, 10.0]);
    let hv = [1.0, 1.0, 0.0];
    let zs = [0.7, 0.7 - 0.3125, 0.7 - 0.625];
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
                    let y = side * 0.5 * gx[i] * hv[j];
                    b.push_str(&format!(",{:.6},{y:.6},{:.6}", xs[i], zs[j]));
                }
            }
            b.push_str(",0.0,1.0,0.0,1.0;");
            bodies.push(b);
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

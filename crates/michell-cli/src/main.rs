//! `michell` — thin-ship wave resistance from hull files.

mod formats;

use formats::{load_hull, parse_pair, parse_range, write_hull_file, LoadSettings, Source};
use michell::{Conditions, Fluid, Hull, Placement, WaveOptions, STANDARD_GRAVITY};
use std::collections::HashMap;

const KNOT: f64 = 1852.0 / 3600.0; // m/s

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("resistance") => cmd_resistance(&args[1..]),
        Some("info") => cmd_info(&args[1..]),
        Some("loft") => cmd_loft(&args[1..]),
        Some("wigley") => cmd_wigley(&args[1..]),
        Some("help") | Some("-h") | Some("--help") | None => {
            print!("{HELP}");
            Ok(())
        }
        Some(other) => Err(format!("unknown command {other:?}; run `michell help`")),
    }
}

const HELP: &str = "\
michell — thin-ship wave resistance (Michell's integral) + ITTC-57 friction

USAGE
  michell resistance <hull>... --speeds A[:B:STEP] [options]  resistance curve
  michell info <hull>... [options]                            geometry & diagnostics
  michell loft <offsets|iges> -o OUT.hull [options]           convert to a control net
  michell wigley [-o OUT.hull] [--length L --beam B --draft T]

HULL INPUTS (sniffed by header / extension)
  *.hull            canonical B-spline control net (exact)
  offsets table     `michell-offsets v1` station x waterline half-beams (lofted)
  *.igs, *.iges     untrimmed NURBS surface(s) (sampled and lofted)

MULTIHULLS
  Pass several hulls; each may carry a placement suffix `@x=DX,y=Y`:
      michell resistance vaka.hull ama.igs@y=1.9 ama.igs@y=-1.9 --speeds 3:8:0.5
  y positions the hull's centerplane; x is added to the hull file's own x
  coordinates (0 keeps them, so hulls modelled in position line up).
  Wave interference is computed exactly (thin-ship superposition); the IF
  column reports combined Rw / sum of standalone Rw. Froude numbers use the
  longest hull's length. Demihulls are assumed symmetric about their own
  centerplanes.

SPEED SELECTION (resistance)
  --speeds A[:B:STEP]   speeds in m/s (inclusive range)
  --froude A[:B:STEP]   length Froude numbers instead of speeds
  --knots               interpret and display speeds in knots

IMPORT / LOFT OPTIONS (offsets and IGES inputs)
  --waterline Z         IGES: design waterline height in the file frame,
                        metres after unit conversion (z up; default 0)
  --centerplane Y       IGES: transverse position of the hull centerplane
                        (default: auto-detect; full shells fold about their
                        midplane, half hulls measure from y = 0)
  --samples NxM         IGES: sample grid, stations x waterlines (default 121x33)
  --fit-degree PxQ      spline degrees for the loft (default 3x3)
  --fit-control NxM     loft control net (default: 12x8 offsets, 20x12 IGES)

PHYSICS OPTIONS
  --fluid NAME          seawater | freshwater, at 15 C (default seawater)
  --rho R, --nu V       override density [kg/m3] / kinematic viscosity [m2/s]
  --gravity G           override g [m/s2]
  --form-factor K       viscous form factor (1+k), default 0
  --rel-tol T           wave-integral relative tolerance (default 1e-5)

OUTPUT
  --json                machine-readable output (resistance, info)
";

// ---------------------------------------------------------------------------
// Argument handling
// ---------------------------------------------------------------------------

struct Parsed {
    positional: Vec<String>,
    flags: HashMap<String, String>,
    switches: Vec<String>,
}

const SWITCHES: &[&str] = &["--json", "--knots"];

fn parse_args(args: &[String]) -> Result<Parsed, String> {
    let mut p = Parsed {
        positional: Vec::new(),
        flags: HashMap::new(),
        switches: Vec::new(),
    };
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if a == "-o" {
            let val = args.get(i + 1).ok_or("-o requires a value")?;
            p.flags.insert("output".to_string(), val.clone());
            i += 1;
        } else if let Some(name) = a.strip_prefix("--") {
            if SWITCHES.contains(&a.as_str()) {
                p.switches.push(a.clone());
            } else {
                let val = args
                    .get(i + 1)
                    .ok_or_else(|| format!("--{name} requires a value"))?;
                p.flags.insert(name.to_string(), val.clone());
                i += 1;
            }
        } else {
            p.positional.push(a.clone());
        }
        i += 1;
    }
    Ok(p)
}

impl Parsed {
    fn switch(&self, name: &str) -> bool {
        self.switches.iter().any(|s| s == name)
    }

    fn f64_flag(&self, name: &str) -> Result<Option<f64>, String> {
        match self.flags.get(name) {
            None => Ok(None),
            Some(v) => v
                .parse::<f64>()
                .map(Some)
                .map_err(|_| format!("--{name}: cannot parse number {v:?}")),
        }
    }

    fn load_settings(&self) -> Result<LoadSettings, String> {
        let mut s = LoadSettings::default();
        if let Some(w) = self.f64_flag("waterline")? {
            s.waterline_z = w;
        }
        s.centerplane = self.f64_flag("centerplane")?;
        if let Some(v) = self.flags.get("samples") {
            s.samples = parse_pair(v)?;
        }
        if let Some(v) = self.flags.get("fit-degree") {
            let (px, pz) = parse_pair(v)?;
            s.fit.degree_x = px;
            s.fit.degree_z = pz;
            s.fit_explicit = true;
        }
        if let Some(v) = self.flags.get("fit-control") {
            let (nx, nz) = parse_pair(v)?;
            s.fit.n_ctrl_x = nx;
            s.fit.n_ctrl_z = nz;
            s.fit_explicit = true;
        }
        Ok(s)
    }

    fn conditions(&self, speed: f64) -> Result<Conditions, String> {
        let mut cond = match self.flags.get("fluid").map(String::as_str) {
            None | Some("seawater") => Conditions::seawater(speed),
            Some("freshwater") => Conditions::freshwater(speed),
            Some(other) => {
                return Err(format!(
                    "--fluid {other:?}: expected seawater or freshwater"
                ))
            }
        };
        if let Some(rho) = self.f64_flag("rho")? {
            cond.fluid.density = rho;
        }
        if let Some(nu) = self.f64_flag("nu")? {
            cond.fluid.kinematic_viscosity = nu;
        }
        if let Some(g) = self.f64_flag("gravity")? {
            cond.gravity = g;
        }
        Ok(cond)
    }
}

/// Parse `path` or `path@x=DX,y=Y` into a file and a placement.
fn parse_hull_spec(spec: &str) -> Result<(String, Placement), String> {
    let Some((path, rest)) = spec.split_once('@') else {
        return Ok((spec.to_string(), Placement::default()));
    };
    let mut place = Placement::default();
    for part in rest.split(',') {
        let (k, v) = part
            .split_once('=')
            .ok_or_else(|| format!("bad placement {rest:?}: expected x=DX,y=Y"))?;
        let val: f64 = v
            .trim()
            .parse()
            .map_err(|_| format!("bad placement value {v:?} in {spec:?}"))?;
        match k.trim() {
            "x" => place.x = val,
            "y" => place.y = val,
            other => return Err(format!("unknown placement key {other:?} (use x, y)")),
        }
    }
    Ok((path.to_string(), place))
}

/// Fleet placements plus the loaded hulls, keyed by file path.
type Fleet = (Vec<(String, Placement)>, HashMap<String, (Hull, Source)>);

/// Load a fleet of hull specs, reading each unique file once.
fn load_fleet(specs: &[String], settings: &LoadSettings) -> Result<Fleet, String> {
    let mut cache: HashMap<String, (Hull, Source)> = HashMap::new();
    let mut fleet = Vec::new();
    for spec in specs {
        let (path, place) = parse_hull_spec(spec)?;
        if !cache.contains_key(&path) {
            let loaded = load_hull(&path, settings)?;
            cache.insert(path.clone(), loaded);
        }
        fleet.push((path, place));
    }
    Ok((fleet, cache))
}

fn describe_source(source: &Source) -> Vec<String> {
    match source {
        Source::Native => vec!["source: native control net (exact)".into()],
        Source::Offsets(r) => vec![format!(
            "source: offsets table, lofted (max residual {:.3e} m at x={:.3} z={:.3}, rms {:.3e} m)",
            r.max_residual, r.max_residual_at.0, r.max_residual_at.1, r.rms_residual
        )],
        Source::Iges(r) => {
            let sides = if r.two_sided {
                format!("full shell folded about y = {:.4} m", r.centerplane)
            } else if r.mirrored {
                format!("port half mirrored about y = {:.4} m", r.centerplane)
            } else {
                format!("half hull about y = {:.4} m", r.centerplane)
            };
            let mut v = vec![format!(
                "source: IGES ({} patches, units scale {}, {sides})",
                r.patches, r.units_scale
            )];
            v.push(format!(
                "geometry: draft {:.4} m, x {:.4}..{:.4} m, fold asymmetry {:.3e} m",
                r.draft, r.x_range.0, r.x_range.1, r.max_asymmetry
            ));
            v.push(format!(
                "loft: max residual {:.3e} m at x={:.3} z={:.3}, rms {:.3e} m, \
                 {} failed inversions, {} ambiguous samples",
                r.fit.max_residual,
                r.fit.max_residual_at.0,
                r.fit.max_residual_at.1,
                r.fit.rms_residual,
                r.failed_inversions,
                r.ambiguous_samples
            ));
            v
        }
    }
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

fn cmd_info(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err("usage: michell info <hull>... [options]".into());
    }
    let (fleet, cache) = load_fleet(&p.positional, &p.load_settings()?)?;
    if p.switch("--json") {
        let mut out = String::from("[");
        for (i, (path, _)) in fleet.iter().enumerate() {
            let (hull, _) = &cache[path];
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"path\":{path:?},\"length\":{},\"draft\":{},\"wetted_surface\":{},\
                 \"displaced_volume\":{}}}",
                hull.length(),
                hull.draft(),
                hull.wetted_surface(),
                hull.displaced_volume()
            ));
        }
        out.push(']');
        println!("{out}");
        return Ok(());
    }
    for (i, (path, place)) in fleet.iter().enumerate() {
        let (hull, source) = &cache[path];
        if i > 0 {
            println!();
        }
        if *place == Placement::default() {
            println!("hull: {path}");
        } else {
            println!("hull: {path} (placed at x{:+}, y = {})", place.x, place.y);
        }
        for line in describe_source(source) {
            println!("{line}");
        }
        println!("length          {:>10.4} m", hull.length());
        println!("draft           {:>10.4} m", hull.draft());
        println!("wetted surface  {:>10.4} m^2", hull.wetted_surface());
        println!("displaced vol   {:>10.4} m^3", hull.displaced_volume());
        let s = hull.surface();
        println!(
            "spline          degree {}x{}, control net {}x{}",
            s.degree_x(),
            s.degree_z(),
            s.n_ctrl_x(),
            s.n_ctrl_z()
        );
    }
    if fleet.len() > 1 {
        let total_s: f64 = fleet.iter().map(|(p, _)| cache[p].0.wetted_surface()).sum();
        let total_v: f64 = fleet
            .iter()
            .map(|(p, _)| cache[p].0.displaced_volume())
            .sum();
        println!();
        println!(
            "fleet: {} hulls, wetted surface {total_s:.4} m^2, displaced vol {total_v:.4} m^3",
            fleet.len()
        );
    }
    Ok(())
}

fn cmd_resistance(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err(
            "usage: michell resistance <hull>[@x=DX,y=Y]... --speeds A[:B:STEP] [options]".into(),
        );
    }
    let (fleet, cache) = load_fleet(&p.positional, &p.load_settings()?)?;
    let members: Vec<(&Hull, Placement)> = fleet
        .iter()
        .map(|(path, place)| (&cache[path].0, *place))
        .collect();
    let multi = members.len() > 1;
    // Reference length for Froude number: the longest hull.
    let l_ref = members
        .iter()
        .map(|(h, _)| h.length())
        .fold(0.0f64, f64::max);
    let total_s: f64 = members.iter().map(|(h, _)| h.wetted_surface()).sum();
    let total_v: f64 = members.iter().map(|(h, _)| h.displaced_volume()).sum();

    let knots = p.switch("--knots");
    let g = p.f64_flag("gravity")?.unwrap_or(STANDARD_GRAVITY);
    let speeds: Vec<f64> = match (p.flags.get("speeds"), p.flags.get("froude")) {
        (Some(_), Some(_)) => return Err("give either --speeds or --froude, not both".into()),
        (Some(s), None) => {
            let v = parse_range(s)?;
            if knots {
                v.into_iter().map(|u| u * KNOT).collect()
            } else {
                v
            }
        }
        (None, Some(f)) => parse_range(f)?
            .into_iter()
            .map(|fr| fr * (g * l_ref).sqrt())
            .collect(),
        (None, None) => return Err("select speeds with --speeds or --froude".into()),
    };

    let form_factor = p.f64_flag("form-factor")?.unwrap_or(0.0);
    let mut wave_opts = WaveOptions::default();
    if let Some(t) = p.f64_flag("rel-tol")? {
        wave_opts.rel_tol = t;
    }

    let mut rows = Vec::new();
    for &u in &speeds {
        let cond = p.conditions(u)?;
        let r = michell::multihull_resistance_with(&members, &cond, &wave_opts, form_factor)
            .map_err(|e| format!("at U = {u} m/s: {e}"))?;
        rows.push((u, cond, r));
    }

    if p.switch("--json") {
        let fluid = rows
            .first()
            .map(|(_, c, _)| c.fluid)
            .unwrap_or(Fluid::SEAWATER_15C);
        let mut out = String::from("{\"hulls\":[");
        for (i, (path, place)) in fleet.iter().enumerate() {
            let hull = &cache[path].0;
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"path\":{path:?},\"placement\":{{\"x\":{},\"y\":{}}},\"length\":{},\
                 \"draft\":{},\"wetted_surface\":{},\"displaced_volume\":{}}}",
                place.x,
                place.y,
                hull.length(),
                hull.draft(),
                hull.wetted_surface(),
                hull.displaced_volume()
            ));
        }
        out.push_str("],");
        out.push_str(&format!(
            "\"fluid\":{{\"density\":{},\"kinematic_viscosity\":{}}},\"form_factor\":{},",
            fluid.density, fluid.kinematic_viscosity, form_factor
        ));
        out.push_str("\"points\":[");
        for (i, (u, cond, r)) in rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"speed\":{u},\"froude\":{},\"rw\":{},\"rv\":{},\"total\":{},\
                 \"effective_power\":{},\"interference\":{},\"cw\":{},\"cv\":{},\"ct\":{},\
                 \"wave_est_rel_error\":{}}}",
                cond.froude_number(l_ref),
                r.wave.resistance,
                r.viscous_total,
                r.total,
                r.effective_power,
                r.interference,
                r.cw,
                r.cv,
                r.ct,
                r.wave.est_rel_error
            ));
        }
        out.push_str("]}");
        println!("{out}");
        return Ok(());
    }

    let mut seen: Vec<&str> = Vec::new();
    for (path, _) in &fleet {
        if seen.contains(&path.as_str()) {
            continue;
        }
        seen.push(path);
        println!("hull: {path}");
        for line in describe_source(&cache[path].1) {
            println!("{line}");
        }
    }
    if multi {
        let placements: Vec<String> = fleet
            .iter()
            .map(|(path, pl)| format!("{path}@x={},y={}", pl.x, pl.y))
            .collect();
        println!("fleet: {}", placements.join("  "));
    }
    println!(
        "L(ref) = {l_ref:.3} m, S = {total_s:.3} m^2, vol = {total_v:.3} m^3, form factor {form_factor}",
    );
    let u_label = if knots { "U[kn]" } else { "U[m/s]" };
    if multi {
        println!(
            "{:>8} {:>6} {:>11} {:>11} {:>11} {:>9} {:>7} {:>9} {:>9}",
            u_label, "Fn", "Rw[N]", "Rv[N]", "Rt[N]", "Pe[kW]", "IF", "Cw*1e3", "Ct*1e3"
        );
    } else {
        println!(
            "{:>8} {:>6} {:>11} {:>11} {:>11} {:>9} {:>9} {:>9}",
            u_label, "Fn", "Rw[N]", "Rv[N]", "Rt[N]", "Pe[kW]", "Cw*1e3", "Ct*1e3"
        );
    }
    for (u, cond, r) in &rows {
        let shown = if knots { u / KNOT } else { *u };
        if multi {
            println!(
                "{shown:>8.3} {:>6.3} {:>11.2} {:>11.2} {:>11.2} {:>9.3} {:>7.3} {:>9.4} {:>9.4}",
                cond.froude_number(l_ref),
                r.wave.resistance,
                r.viscous_total,
                r.total,
                r.effective_power / 1000.0,
                r.interference,
                r.cw * 1e3,
                r.ct * 1e3
            );
        } else {
            println!(
                "{shown:>8.3} {:>6.3} {:>11.2} {:>11.2} {:>11.2} {:>9.3} {:>9.4} {:>9.4}",
                cond.froude_number(l_ref),
                r.wave.resistance,
                r.viscous_total,
                r.total,
                r.effective_power / 1000.0,
                r.cw * 1e3,
                r.ct * 1e3
            );
        }
    }
    Ok(())
}

fn cmd_loft(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    let [path] = p.positional.as_slice() else {
        return Err("usage: michell loft <offsets|iges> -o OUT.hull [options]".into());
    };
    let out_path = p
        .flags
        .get("output")
        .ok_or("loft requires an output path: -o OUT.hull")?;
    let (hull, source) = load_hull(path, &p.load_settings()?)?;
    if matches!(source, Source::Native) {
        return Err(format!("{path} is already a control-net file"));
    }
    std::fs::write(out_path, write_hull_file(&hull))
        .map_err(|e| format!("cannot write {out_path}: {e}"))?;
    for line in describe_source(&source) {
        println!("{line}");
    }
    println!("wrote {out_path}");
    Ok(())
}

fn cmd_wigley(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    let l = p.f64_flag("length")?.unwrap_or(10.0);
    let b = p.f64_flag("beam")?.unwrap_or(l / 10.0);
    let t = p.f64_flag("draft")?.unwrap_or(b * 0.625);
    let hull = michell::hulls::wigley(l, b, t).map_err(|e| format!("{e}"))?;
    let text = write_hull_file(&hull);
    match p.flags.get("output") {
        Some(path) => {
            std::fs::write(path, text).map_err(|e| format!("cannot write {path}: {e}"))?;
            println!("wrote {path} (Wigley L={l} B={b} T={t})");
        }
        None => print!("{text}"),
    }
    Ok(())
}

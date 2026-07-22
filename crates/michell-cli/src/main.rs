//! `michell` — thin-ship wave resistance from hull files.

mod formats;
mod gridio;
mod json;
mod manifest;
mod png;
mod render;
mod view;

use formats::{load_hulls, parse_pair, parse_range, write_hull_file, LoadSettings, Source};
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
        Some("sweep") => {
            // A JSON manifest is the preferred sweep interface.
            if let Some(path) = args.get(1).filter(|a| a.ends_with(".json")) {
                if args.len() > 2 {
                    return Err("manifest sweeps take no further arguments; put \
                                everything in the JSON file"
                        .into());
                }
                manifest::run(path)
            } else {
                cmd_sweep(&args[1..])
            }
        }
        Some("info") => cmd_info(&args[1..]),
        Some("spectrum") => cmd_spectrum(&args[1..]),
        Some("wake") => cmd_wake(&args[1..]),
        Some("render") => cmd_render(&args[1..]),
        Some("view") => view::cmd_view(&args[1..]),
        Some("loft") => cmd_loft(&args[1..]),
        Some("place") => cmd_place(&args[1..]),
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
  michell spectrum <hull>... --speed U [options]              free-wave spectrum
  michell wake <hull>... --speed U [-o wake.png] [options]    Kelvin wake heatmap
  michell view <hull>... --speed U [--port N] [options]       interactive fleet viewer
  michell loft <offsets|iges> -o OUT.hull [options]           convert to a control net
  michell place <hull>[@dx=..,dy=..,dz=..]... -o OUT.igs      write posed CAD geometry
  michell wigley [-o OUT.hull] [--length L --beam B --draft T]

HULL INPUTS (sniffed by header / extension)
  *.hull            canonical B-spline control net (exact); with a
                    `waterline` key, a re-situatable full-band body
  offsets table     `michell-offsets v1` station x waterline half-beams (lofted)
  *.grid.json       derivative-augmented sample grid (the IR written by
                    --dump-grid), lofted on load
  *.igs, *.iges     untrimmed NURBS surface(s) (sampled and lofted, with
                    surface slopes recovered from the CAD geometry)
  *.stl             triangle mesh, binary or ASCII (ray-sampled and lofted);
                    requires --units; quality tracks the export's chord
                    tolerance — use fine tessellations

MULTIHULLS
  Pass several hulls; each may carry a placement suffix:
      michell resistance vaka.hull ama.igs@y=1.9 ama.igs@y=-1.9 --speeds 3:8:0.5
  Suffix keys: y=Y  places the hull's centerplane absolutely (single-hull
  files only); dy=S shifts transversely; x=DX / dx=DX shifts longitudinally
  (added to the file's own x coordinates).
  An IGES file containing a whole multihull imports as a fleet: hulls are
  detected by clustering wetted patches, each at its detected centerplane
  (dry structure like beams and decks is dropped); dy/dx then shift all of
  them together. Wave interference is computed exactly (thin-ship
  superposition); the IF column reports combined Rw / sum of standalone Rw.
  Froude numbers use the longest hull's length. Demihulls are assumed
  symmetric about their own centerplanes.

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
  --units U             STL: scale to metres (mm|cm|m|in|ft or a number);
                        required for STL, which has no units field
  --samples NxM         IGES/STL: sample grid, stations x waterlines
                        (default 121x33)
  --fit-degree PxQ      spline degrees for the loft (default 3x3)
  --fit-control NxM     loft control net (default: 12x8 offsets, 20x12 IGES)
  --fit-deriv-weight W  relative weight of sampled surface slopes in the
                        loft (default 1; 0 fits values only)
  --dump-grid PATH      write the sampled grid(s) as *.grid.json before
                        lofting (multihull files get -0, -1, ... suffixes)

PHYSICS OPTIONS
  --fluid NAME          seawater | freshwater, at 15 C (default seawater)
  --rho R, --nu V       override density [kg/m3] / kinematic viscosity [m2/s]
  --gravity G           override g [m/s2]
  --form-factor K       viscous form factor (1+k), default 0
  --rel-tol T           wave-integral relative tolerance (default 1e-5)
  --heel DEG            (resistance) heel the whole fleet DEG degrees about the
                        platform's longitudinal axis; adds the tilted-thickness
                        wave-making (|DEG| < 90). Viscous resistance is
                        unchanged: no waterline re-clip, no yaw side-force

WAVE FIELD (spectrum, wake)
  Both take one speed: --speed U (m/s; knots with --knots) or --froude F.
  The ship advances toward +x (the bow is the high-x end of the hull file).

  spectrum: the far-field free-wave spectrum by propagation angle theta —
  where the wave energy goes. CSV to stdout (or --json): amplitude density
  |A| [m/rad], phase, dRw/dtheta [N/rad], cumulative resistance fraction;
  integrating dRw/dtheta recovers Rw (cross-check on stderr).
    --points N          theta samples (default 721)

  wake: the Kelvin wave pattern zeta(x, y) reconstructed from the spectrum.
    -o OUT.png|.csv|.json  output (default wake.png); the PNG is a heatmap
                        (blue trough, red crest, gray hull waterplanes;
                        faded where not astern of every hull)
    --region X0:X1:Y0:Y1   window in fleet coordinates [m]
                        (default: ~3 hull lengths of wake, auto width)
    --size WxH          grid points (default 900 wide, aspect-matched)
    --zmax M            color saturation elevation [m] (default: 99.5th
                        percentile of |zeta|)
  The pattern is the far-field free-wave part of the linear solution: it is
  physical astern of each hull, not on or ahead of it.

  render: a 3D shot of the fleet sitting in its wake (software-rendered).
    -o OUT.png          output (default render.png)
    --region X0:X1:Y0:Y1   water extent [m] (default: as wake)
    --size WxH          image pixels (default 1200x800)
    --grid N            water mesh columns (default 700)
    --camera AZ:EL[:D]  degrees off dead astern (positive to starboard),
                        elevation degrees, distance m (default 35:18, auto)
    --z-scale S         vertical exaggeration of the water (default 1)
    --freeboard F       topsides height above waterline [m] (default T/2)
    --zmax M            tint saturation elevation [m] (default: 99.5th pct)
  Same caveat as wake (paled where not astern of every hull); hulls sit at
  the static waterline, without dynamic sinkage or trim.

SWEEPS
  michell sweep study.json      preferred: a JSON manifest referencing
                                full-band .hull bodies (from `michell loft`),
                                with speed/waterline axes, per-hull load
                                (mass/lcg/vcg), point loads (mass at an
                                offset from a hull), pose axes, and a derived
                                fleet CG — see the README for the schema; when
                                the fleet carries mass each row also rolls up
                                the heel behaviour (GZ-curve summaries + a
                                per-angle resistance rise)
  michell loft boat.igs --waterline Z -o boat
                                decompose an IGES multihull into full-band
                                body files (boat-port.hull, ...); --wetted
                                keeps the old single-hull wetted output

PLACE (reconstruct CAD geometry from a studied configuration)
  michell place ama.igs@dy=1.7,dz=0.05 ama.igs@dy=-1.7,dz=0.05 -o boat.igs
  Writes the input geometry, posed, as a new IGES file (untrimmed 128
  surfaces, metres) for import back into CAD — e.g. amas at the dx/dy/dz a
  sweep found good. Inputs: IGES files (surfaces pass through exactly, with
  each spec's pose applied to every hull in that file) and .hull control
  nets or bodies (the half-breadth spline converts exactly to a mirrored
  pair of surfaces). Suffix keys, all optional:
      dx / x    longitudinal shift [m]        dy  transverse shift [m]
      dz        immersion [m], + is deeper    y   absolute centerplane
                                                  (.hull inputs only)
      trim      design trim [deg], + raises the +x end, about pivot=X
                (default: the hull's x mid) at the waterline
  --waterline Z         CAD height of the design waterline (default 0);
                        poses and the sinkage re-expression are relative to it
  --sinkage S           platform sinkage [m] from a solved equilibrium row;
                        + moves every hull deeper (the DWL stays at Z)
  --platform-trim DEG   platform pitch about (--pivot-x, the waterline)
  --pivot-x X           platform trim pivot station (default 0)
  Bounded (143/141) source patches are exported as full base surfaces with
  their parameter range restricted to the bounded box.

SWEEPS (flag form, IGES inputs; hulls modelled in position)
  michell sweep boat.igs --speeds 3:8:1 [axes...]         long-form CSV/JSON
  --axis waterline=A:B:S        raw waterline sweep
  --axis SEL:PARAM=A[:B:S]      design-pose sweep; SEL = file stem, or
                                stem#0, stem#0+2 for specific hulls (indexed
                                by transverse position); PARAM one of
                                dz (immersion, +down), dx, dy,
                                spread (outboard shift, sign follows side),
                                trim (degrees, + raises the +x end),
                                scale (uniform size factor, >0; 1 = unchanged)
  --float weight=A[:B:S]        solve sinkage (and pitch, with lcg) so the
  --float lcg=A[:B:S]           fleet floats each load; excludes a waterline
                                axis; single input file only
  Rows are the Cartesian product of all axes x speeds; each row carries the
  solved sinkage/trim, displacement, LCB, and the resistance breakdown.
  Output is CSV on stdout (use --json for JSON); progress goes to stderr.
  All poses are static (hydrostatic) attitudes: no dynamic sinkage/trim.

OUTPUT
  --json                machine-readable output (resistance, info, sweep)
  --csv                 sweep: comma-separated output (the default)
";

// ---------------------------------------------------------------------------
// Argument handling
// ---------------------------------------------------------------------------

struct Parsed {
    positional: Vec<String>,
    /// Every `--flag value` occurrence, in order (flags may repeat).
    pairs: Vec<(String, String)>,
    switches: Vec<String>,
}

const SWITCHES: &[&str] = &["--json", "--knots", "--csv", "--wetted"];

fn parse_args(args: &[String]) -> Result<Parsed, String> {
    let mut p = Parsed {
        positional: Vec::new(),
        pairs: Vec::new(),
        switches: Vec::new(),
    };
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if a == "-o" {
            let val = args.get(i + 1).ok_or("-o requires a value")?;
            p.pairs.push(("output".to_string(), val.clone()));
            i += 1;
        } else if let Some(name) = a.strip_prefix("--") {
            if SWITCHES.contains(&a.as_str()) {
                p.switches.push(a.clone());
            } else {
                let val = args
                    .get(i + 1)
                    .ok_or_else(|| format!("--{name} requires a value"))?;
                p.pairs.push((name.to_string(), val.clone()));
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

    /// Last occurrence of a flag.
    fn flag(&self, name: &str) -> Option<&String> {
        self.pairs
            .iter()
            .rev()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v)
    }

    /// All occurrences of a flag, in order.
    fn all(&self, name: &str) -> Vec<&String> {
        self.pairs
            .iter()
            .filter(|(k, _)| k == name)
            .map(|(_, v)| v)
            .collect()
    }

    fn f64_flag(&self, name: &str) -> Result<Option<f64>, String> {
        match self.flag(name) {
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
        if let Some(u) = self.flag("units") {
            s.units = Some(formats::parse_units(u)?);
        }
        if let Some(v) = self.flag("samples") {
            s.samples = parse_pair(v)?;
        }
        if let Some(v) = self.flag("fit-degree") {
            let (px, pz) = parse_pair(v)?;
            s.fit.degree_x = px;
            s.fit.degree_z = pz;
            s.fit_explicit = true;
        }
        if let Some(v) = self.flag("fit-control") {
            let (nx, nz) = parse_pair(v)?;
            s.fit.n_ctrl_x = nx;
            s.fit.n_ctrl_z = nz;
            s.fit_explicit = true;
        }
        if let Some(w) = self.f64_flag("fit-deriv-weight")? {
            s.fit.derivative_weight = w;
        }
        if let Some(p) = self.flag("dump-grid") {
            s.dump_grid = Some(p.clone());
        }
        Ok(s)
    }

    fn conditions(&self, speed: f64) -> Result<Conditions, String> {
        let mut cond = match self.flag("fluid").map(String::as_str) {
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

/// Placement request from a hull spec suffix.
#[derive(Default, Clone, Copy)]
struct SpecPlacement {
    /// Absolute centerplane position (single-hull files only).
    y_abs: Option<f64>,
    /// Transverse shift applied to the file's detected placements.
    dy: f64,
    /// Longitudinal shift added to the file's x coordinates.
    dx: f64,
}

/// Parse `path` or `path@key=V,...` (keys: `y` absolute centerplane,
/// `dy` transverse shift, `x`/`dx` longitudinal shift).
fn parse_hull_spec(spec: &str) -> Result<(String, SpecPlacement), String> {
    let Some((path, rest)) = spec.split_once('@') else {
        return Ok((spec.to_string(), SpecPlacement::default()));
    };
    let mut place = SpecPlacement::default();
    for part in rest.split(',') {
        let (k, v) = part
            .split_once('=')
            .ok_or_else(|| format!("bad placement {rest:?}: expected key=value pairs"))?;
        let val: f64 = v
            .trim()
            .parse()
            .map_err(|_| format!("bad placement value {v:?} in {spec:?}"))?;
        match k.trim() {
            "y" => place.y_abs = Some(val),
            "dy" => place.dy = val,
            "x" | "dx" => place.dx = val,
            other => return Err(format!("unknown placement key {other:?} (use y, dy, x/dx)")),
        }
    }
    if place.y_abs.is_some() && place.dy != 0.0 {
        return Err(format!(
            "{spec:?}: give either y (absolute) or dy (shift), not both"
        ));
    }
    Ok((path.to_string(), place))
}

/// One hull of the working fleet: a file may contribute several (a multihull
/// IGES export), and a spec's `@x,y` offset shifts everything from that file.
struct Member {
    path: String,
    hull: Hull,
    placement: Placement,
    source: Source,
}

/// Load a fleet of hull specs, reading each unique file once. A file
/// containing several hulls contributes all of them, each at its detected
/// placement plus the spec's offset.
fn load_fleet(specs: &[String], settings: &LoadSettings) -> Result<Vec<Member>, String> {
    let mut cache: HashMap<String, Vec<(Hull, Placement, Source)>> = HashMap::new();
    let mut members = Vec::new();
    for spec in specs {
        let (path, sp) = parse_hull_spec(spec)?;
        if !cache.contains_key(&path) {
            cache.insert(path.clone(), load_hulls(&path, settings)?);
        }
        let hulls = &cache[&path];
        if sp.y_abs.is_some() && hulls.len() > 1 {
            return Err(format!(
                "{spec:?}: absolute y placement is ambiguous for a file with {} hulls; \
                 use dy=SHIFT instead",
                hulls.len()
            ));
        }
        for (hull, detected, source) in hulls {
            members.push(Member {
                path: path.clone(),
                hull: hull.clone(),
                placement: Placement {
                    x: detected.x + sp.dx,
                    y: sp.y_abs.unwrap_or(detected.y + sp.dy),
                },
                source: source.clone(),
            });
        }
    }
    Ok(members)
}

/// Slope-channel residual summary, when the loft fit any.
fn slope_line(r: &michell::fit::FitReport) -> Option<String> {
    match (r.fx_residual, r.fz_residual) {
        (None, None) => None,
        (fx, fz) => {
            let one = |name: &str, c: Option<michell::fit::ChannelResiduals>| {
                c.map(|c| format!("{name} max {:.3e} rms {:.3e}", c.max, c.rms))
            };
            let parts: Vec<String> = [one("dfdx", fx), one("dfdz", fz)]
                .into_iter()
                .flatten()
                .collect();
            Some(format!("slopes: {}", parts.join(", ")))
        }
    }
}

fn describe_source(source: &Source) -> Vec<String> {
    match source {
        Source::Native => vec!["source: native control net (exact)".into()],
        Source::Body(r) => {
            let mut v = vec![format!(
                "source: full-band body, situated at its design waterline \
                 (loft max residual {:.3e} m, rms {:.3e} m)",
                r.max_residual, r.rms_residual
            )];
            v.extend(slope_line(r));
            v
        }
        Source::Offsets(r) => vec![format!(
            "source: offsets table, lofted (max residual {:.3e} m at x={:.3} z={:.3}, rms {:.3e} m)",
            r.max_residual, r.max_residual_at.0, r.max_residual_at.1, r.rms_residual
        )],
        Source::Grid(r) => {
            let mut v = vec![format!(
                "source: sample grid, lofted (max residual {:.3e} m at x={:.3} z={:.3}, rms {:.3e} m)",
                r.max_residual, r.max_residual_at.0, r.max_residual_at.1, r.rms_residual
            )];
            v.extend(slope_line(r));
            v
        }
        Source::Stl(r) => {
            let sides = if r.two_sided {
                format!("full shell folded about y = {:.4} m", r.centerplane)
            } else {
                format!("one-sided about y = {:.4} m", r.centerplane)
            };
            vec![
                format!(
                    "source: STL mesh ({} triangles, units scale {}, {sides})",
                    r.patches, r.units_scale
                ),
                format!(
                    "geometry: draft {:.4} m, x {:.4}..{:.4} m, fold asymmetry {:.3e} m",
                    r.draft, r.x_range.0, r.x_range.1, r.max_asymmetry
                ),
                format!(
                    "loft: max residual {:.3e} m at x={:.3} z={:.3}, rms {:.3e} m, \
                     {} ambiguous samples",
                    r.fit.max_residual,
                    r.fit.max_residual_at.0,
                    r.fit.max_residual_at.1,
                    r.fit.rms_residual,
                    r.ambiguous_samples
                ),
            ]
        }
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
                 {} failed inversions, {} ambiguous samples, {} slope gaps",
                r.fit.max_residual,
                r.fit.max_residual_at.0,
                r.fit.max_residual_at.1,
                r.fit.rms_residual,
                r.failed_inversions,
                r.ambiguous_samples,
                r.derivative_gaps
            ));
            v.extend(slope_line(&r.fit));
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
    let members = load_fleet(&p.positional, &p.load_settings()?)?;
    if p.switch("--json") {
        let mut out = String::from("[");
        for (i, m) in members.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"path\":{:?},\"placement\":{{\"x\":{},\"y\":{}}},\"length\":{},\
                 \"draft\":{},\"wetted_surface\":{},\"displaced_volume\":{}}}",
                m.path,
                m.placement.x,
                m.placement.y,
                m.hull.length(),
                m.hull.draft(),
                m.hull.wetted_surface(),
                m.hull.displaced_volume()
            ));
        }
        out.push(']');
        println!("{out}");
        return Ok(());
    }
    // Several hulls (a multihull file, or several specs): lead with a fleet
    // overview so the shape of the fleet is visible before the per-hull
    // detail blocks.
    if members.len() > 1 {
        println!("fleet: {} hulls", members.len());
        for (i, m) in members.iter().enumerate() {
            println!(
                "  [{}] {}  y {:+.4} m  length {:.4} m  draft {:.4} m  \
                 displaced vol {:.4} m^3",
                i + 1,
                m.path,
                m.placement.y,
                m.hull.length(),
                m.hull.draft(),
                m.hull.displaced_volume()
            );
        }
        let total_s: f64 = members.iter().map(|m| m.hull.wetted_surface()).sum();
        let total_v: f64 = members.iter().map(|m| m.hull.displaced_volume()).sum();
        println!("  total: wetted surface {total_s:.4} m^2, displaced vol {total_v:.4} m^3");
    }
    for (i, m) in members.iter().enumerate() {
        if i > 0 || members.len() > 1 {
            println!();
        }
        let name = if members.len() > 1 {
            format!("hull [{}]: {}", i + 1, m.path)
        } else {
            format!("hull: {}", m.path)
        };
        if m.placement == Placement::default() {
            println!("{name}");
        } else {
            println!(
                "{name} (placed at dx {:+.3} m, y {:.4} m)",
                m.placement.x, m.placement.y
            );
        }
        for line in describe_source(&m.source) {
            println!("{line}");
        }
        println!("length          {:>10.4} m", m.hull.length());
        println!("draft           {:>10.4} m", m.hull.draft());
        println!("wetted surface  {:>10.4} m^2", m.hull.wetted_surface());
        println!("displaced vol   {:>10.4} m^3", m.hull.displaced_volume());
        let s = m.hull.surface();
        println!(
            "spline          degree {}x{}, control net {}x{}",
            s.degree_x(),
            s.degree_z(),
            s.n_ctrl_x(),
            s.n_ctrl_z()
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
    let loaded = load_fleet(&p.positional, &p.load_settings()?)?;
    let members: Vec<(&Hull, Placement)> = loaded.iter().map(|m| (&m.hull, m.placement)).collect();
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
    let speeds: Vec<f64> = match (p.flag("speeds"), p.flag("froude")) {
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
    let heel_deg = p.f64_flag("heel")?.unwrap_or(0.0);
    let heel = heel_deg.to_radians();
    let mut wave_opts = WaveOptions::default();
    if let Some(t) = p.f64_flag("rel-tol")? {
        wave_opts.rel_tol = t;
    }

    let mut rows = Vec::new();
    for &u in &speeds {
        let cond = p.conditions(u)?;
        let r = if heel != 0.0 {
            michell::multihull_resistance_heeled(&members, &cond, &wave_opts, form_factor, heel)
        } else {
            michell::multihull_resistance_with(&members, &cond, &wave_opts, form_factor)
        }
        .map_err(|e| format!("at U = {u} m/s: {e}"))?;
        rows.push((u, cond, r));
    }

    if p.switch("--json") {
        let fluid = rows
            .first()
            .map(|(_, c, _)| c.fluid)
            .unwrap_or(Fluid::SEAWATER_15C);
        let mut out = String::from("{\"hulls\":[");
        for (i, m) in loaded.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"path\":{:?},\"placement\":{{\"x\":{},\"y\":{}}},\"length\":{},\
                 \"draft\":{},\"wetted_surface\":{},\"displaced_volume\":{}}}",
                m.path,
                m.placement.x,
                m.placement.y,
                m.hull.length(),
                m.hull.draft(),
                m.hull.wetted_surface(),
                m.hull.displaced_volume()
            ));
        }
        out.push_str("],");
        out.push_str(&format!(
            "\"fluid\":{{\"density\":{},\"kinematic_viscosity\":{}}},\"form_factor\":{},\
             \"heel_deg\":{},",
            fluid.density, fluid.kinematic_viscosity, form_factor, heel_deg
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
    for m in &loaded {
        if seen.contains(&m.path.as_str()) {
            continue;
        }
        seen.push(&m.path);
        println!("hull: {}", m.path);
        for line in describe_source(&m.source) {
            println!("{line}");
        }
    }
    if multi {
        let placements: Vec<String> = loaded
            .iter()
            .map(|m| format!("{}@x={},y={}", m.path, m.placement.x, m.placement.y))
            .collect();
        println!("fleet: {}", placements.join("  "));
    }
    println!(
        "L(ref) = {l_ref:.3} m, S = {total_s:.3} m^2, vol = {total_v:.3} m^3, form factor {form_factor}",
    );
    if heel != 0.0 {
        println!(
            "heeled {heel_deg:.1} deg about the platform axis (tilted-thickness \
             wave-making only; no waterline re-clip, no yaw side-force)"
        );
    }
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

// ---------------------------------------------------------------------------
// Sweep
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
enum PoseParam {
    Dz,
    Dx,
    Dy,
    Spread,
    TrimDeg,
    Scale,
}

enum Target {
    Waterline,
    Weight,
    Lcg,
    Pose {
        file: usize,
        hulls: Vec<usize>,
        param: PoseParam,
    },
}

struct Axis {
    label: String,
    values: Vec<f64>,
    target: Target,
}

fn cmd_sweep(args: &[String]) -> Result<(), String> {
    use michell::float::{solve_equilibrium, LoadCase};
    use michell::iges::{self, HullPose, ImportOptions, Platform, SourceFleet};

    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err(
            "usage: michell sweep <boat.igs>... --speeds A[:B:S] [--float weight=... \
             [--float lcg=...]] [--axis KEY=A:B:S]..."
                .into(),
        );
    }
    if p.positional.iter().any(|s| s.contains('@')) {
        return Err(
            "sweep does not accept @ placement suffixes; use --axis with dz/dx/dy/spread/trim/scale"
                .into(),
        );
    }
    let settings = p.load_settings()?;
    if settings.centerplane.is_some() {
        return Err("--centerplane is not supported by sweep".into());
    }
    if settings.dump_grid.is_some() {
        return Err("--dump-grid is not supported by sweep (it re-lofts one \
                    grid per pose); use `michell info` or `michell loft`"
            .into());
    }
    let base_wl = settings.waterline_z;
    let opts = ImportOptions {
        waterline_z: base_wl,
        stations: settings.samples.0,
        waterlines: settings.samples.1,
        fit: if settings.fit_explicit {
            settings.fit
        } else {
            ImportOptions::default().fit
        },
        centerplane: None,
    };

    // Load the source fleets and learn each hull's base transverse position.
    struct File {
        stem: String,
        src: SourceFleet,
        base_y: Vec<f64>,
        n: usize,
    }
    let mut files: Vec<File> = Vec::new();
    let mut l_ref = 0.0f64;
    for path in &p.positional {
        let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
        let src = iges::source_fleet(&text, base_wl).map_err(|e| format!("{path}: {e}"))?;
        let poses = vec![HullPose::default(); src.len()];
        let fl = src
            .situate(base_wl, &poses, &Platform::default(), &opts)
            .map_err(|e| format!("{path}: {e}"))?;
        let mut base_y = vec![0.0; src.len()];
        let mut mi = 0;
        for (hi, y) in base_y.iter_mut().enumerate() {
            if fl.dry.contains(&hi) {
                continue;
            }
            *y = fl.members[mi].placement.y;
            l_ref = l_ref.max(fl.members[mi].hull.length());
            mi += 1;
        }
        let stem = std::path::Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(path)
            .to_string();
        eprintln!(
            "loaded {path}: {} hull(s) at y = {:?}",
            src.len(),
            base_y
                .iter()
                .map(|y| (y * 1e4).round() / 1e4)
                .collect::<Vec<_>>()
        );
        files.push(File {
            stem,
            n: src.len(),
            src,
            base_y,
        });
    }

    // Axes: --float entries first, then --axis entries, in CLI order.
    let mut axes: Vec<Axis> = Vec::new();
    for v in p.all("float") {
        let (key, range) = v
            .split_once('=')
            .ok_or_else(|| format!("--float {v:?}: expected weight=... or lcg=..."))?;
        let target = match key.trim() {
            "weight" => Target::Weight,
            "lcg" => Target::Lcg,
            other => return Err(format!("--float key {other:?}: expected weight or lcg")),
        };
        axes.push(Axis {
            label: key.trim().to_string(),
            values: parse_range(range)?,
            target,
        });
    }
    for v in p.all("axis") {
        let (key, range) = v
            .split_once('=')
            .ok_or_else(|| format!("--axis {v:?}: expected KEY=A[:B:S]"))?;
        let values = parse_range(range)?;
        let target = if key.trim() == "waterline" {
            Target::Waterline
        } else {
            let (sel, pname) = key
                .rsplit_once(':')
                .ok_or_else(|| format!("--axis key {key:?}: expected SEL:PARAM or waterline"))?;
            let param = match pname.trim() {
                "dz" => PoseParam::Dz,
                "dx" => PoseParam::Dx,
                "dy" => PoseParam::Dy,
                "spread" => PoseParam::Spread,
                "trim" => PoseParam::TrimDeg,
                "scale" => PoseParam::Scale,
                other => {
                    return Err(format!(
                        "unknown pose parameter {other:?} (dz, dx, dy, spread, trim, scale)"
                    ))
                }
            };
            if pname.trim() == "scale" && values.iter().any(|&v| !(v > 0.0 && v.is_finite())) {
                return Err(format!(
                    "--axis {key:?}: scale factors must be positive and finite"
                ));
            }
            let (stem, idx_spec) = match sel.split_once('#') {
                None => (sel.trim(), None),
                Some((st, is)) => (st.trim(), Some(is)),
            };
            let file = files
                .iter()
                .position(|f| f.stem == stem)
                .ok_or_else(|| format!("no input file with stem {stem:?}"))?;
            let hulls: Vec<usize> = match idx_spec {
                None => (0..files[file].n).collect(),
                Some(is) => is
                    .split('+')
                    .map(|t| {
                        t.trim()
                            .parse::<usize>()
                            .map_err(|_| format!("bad hull index {t:?} in {key:?}"))
                    })
                    .collect::<Result<_, _>>()?,
            };
            if hulls.iter().any(|&h| h >= files[file].n) {
                return Err(format!(
                    "hull index out of range in {key:?}: file has {} hulls",
                    files[file].n
                ));
            }
            Target::Pose { file, hulls, param }
        };
        axes.push(Axis {
            label: key.trim().to_string(),
            values,
            target,
        });
    }

    let float_mode = axes.iter().any(|a| matches!(a.target, Target::Weight));
    if axes.iter().any(|a| matches!(a.target, Target::Lcg)) && !float_mode {
        return Err("--float lcg=... requires --float weight=...".into());
    }
    if float_mode {
        if axes.iter().any(|a| matches!(a.target, Target::Waterline)) {
            return Err("a waterline axis cannot be combined with --float (the \
                        waterline is solved)"
                .into());
        }
        if files.len() != 1 {
            return Err(
                "--float requires a single input file (the whole platform in one \
                 IGES) so the rigid-body equilibrium is well defined"
                    .into(),
            );
        }
    }

    // Speeds.
    let knots = p.switch("--knots");
    let g = p.f64_flag("gravity")?.unwrap_or(STANDARD_GRAVITY);
    let speeds: Vec<f64> = match (p.flag("speeds"), p.flag("froude")) {
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
    let density = p.conditions(1.0)?.fluid.density;

    let points: usize = axes
        .iter()
        .map(|a| a.values.len())
        .product::<usize>()
        .max(1);
    if points * speeds.len() > 100_000 {
        return Err(format!(
            "sweep would produce {} rows; narrow the axes",
            points * speeds.len()
        ));
    }
    eprintln!(
        "sweep: {points} point(s) x {} speed(s){}",
        speeds.len(),
        if float_mode { ", equilibrium mode" } else { "" }
    );

    // Output assembly.
    let header: Vec<String> = axes
        .iter()
        .map(|a| a.label.clone())
        .chain(
            [
                "sinkage",
                "trim_deg",
                "volume",
                "lcb",
                "dry",
                "speed",
                "froude",
                "rw",
                "rv",
                "rt",
                "pe",
                "interference",
                "cw",
                "ct",
            ]
            .iter()
            .map(|s| s.to_string()),
        )
        .collect();
    let json = p.switch("--json");
    let mut out = String::new();
    if json {
        out.push('[');
    } else {
        out.push_str(&header.join(","));
        out.push('\n');
    }
    let mut first_row = true;

    // Mixed-radix strides so a flat point index decomposes into its per-axis
    // values independently — the last axis varies fastest, matching the
    // odometer the serial version walked. This lets any point be evaluated in
    // isolation, which is what makes the fan-out below safe.
    let mut strides = vec![1usize; axes.len()];
    for i in (0..axes.len()).rev() {
        if i + 1 < axes.len() {
            strides[i] = strides[i + 1] * axes[i + 1].values.len();
        }
    }

    // Evaluate one grid point into its speed rows. Everything it reads
    // (`axes`, `files`, `opts`, `wave_opts`, `speeds`, `p`, scalars) is
    // immutable, so points are independent and evaluate in any order/thread.
    let eval_point = |point: usize| -> Result<Vec<Vec<f64>>, String> {
        let vals: Vec<f64> = axes
            .iter()
            .enumerate()
            .map(|(i, a)| a.values[(point / strides[i]) % a.values.len()])
            .collect();

        // Assemble poses and load for this point.
        let mut poses: Vec<Vec<HullPose>> = files
            .iter()
            .map(|f| vec![HullPose::default(); f.n])
            .collect();
        let mut waterline = base_wl;
        let mut weight = None;
        let mut lcg = None;
        for (a, &v) in axes.iter().zip(&vals) {
            match &a.target {
                Target::Waterline => waterline = v,
                Target::Weight => weight = Some(v),
                Target::Lcg => lcg = Some(v),
                Target::Pose { file, hulls, param } => {
                    for &h in hulls {
                        let pose = &mut poses[*file][h];
                        match param {
                            PoseParam::Dz => pose.dz = v,
                            PoseParam::Dx => pose.dx = v,
                            PoseParam::Dy => pose.dy = v,
                            PoseParam::Spread => {
                                pose.dy = if files[*file].base_y[h] < 0.0 { -v } else { v }
                            }
                            PoseParam::TrimDeg => pose.trim = v.to_radians(),
                            PoseParam::Scale => pose.scale = v,
                        }
                    }
                }
            }
        }

        // Situate (raw) or solve (float).
        let mut fleets: Vec<michell::float::FleetState> = Vec::new();
        let (sinkage, trim_deg, volume, lcb, dry) = if let Some(mass) = weight {
            let eq = solve_equilibrium(
                &files[0].src,
                base_wl,
                &poses[0],
                &LoadCase { mass, lcg },
                density,
                &opts,
            )
            .map_err(|e| format!("point {}: {e}", point + 1))?;
            let out = (
                eq.sinkage,
                eq.trim.to_degrees(),
                eq.volume,
                eq.lcb,
                eq.fleet.dry,
            );
            fleets.push(eq.fleet);
            out
        } else {
            let mut volume = 0.0;
            let mut moment = 0.0;
            let mut dry = 0usize;
            for (f, fp) in files.iter().zip(&poses) {
                let fl = f
                    .src
                    .situate(waterline, fp, &Platform::default(), &opts)
                    .map_err(|e| format!("point {}: {e}", point + 1))?;
                dry += fl.dry.len();
                let mut members = Vec::new();
                for m in fl.members {
                    volume += m.hull.displaced_volume();
                    moment += m.hull.lcb_x() * m.hull.displaced_volume();
                    members.push((m.hull, m.placement));
                }
                fleets.push(michell::float::FleetState {
                    members,
                    dry: 0,
                    band_exceeded: 0,
                });
            }
            let lcb = if volume > 0.0 { moment / volume } else { 0.0 };
            (0.0, 0.0, volume, lcb, dry)
        };
        let members: Vec<(&Hull, Placement)> = fleets
            .iter()
            .flat_map(|fl| fl.members.iter().map(|(h, p)| (h, *p)))
            .collect();

        let mut rows = Vec::with_capacity(speeds.len());
        for &u in &speeds {
            let cond = p.conditions(u)?;
            let (rw, rv, rt, pe, iff, cw, ct) = if members.is_empty() {
                (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
            } else {
                let r =
                    michell::multihull_resistance_with(&members, &cond, &wave_opts, form_factor)
                        .map_err(|e| format!("point {} U={u}: {e}", point + 1))?;
                (
                    r.wave.resistance,
                    r.viscous_total,
                    r.total,
                    r.effective_power,
                    r.interference,
                    r.cw,
                    r.ct,
                )
            };
            let froude = u / (g * l_ref).sqrt();
            rows.push(
                vals.iter()
                    .cloned()
                    .chain([
                        sinkage, trim_deg, volume, lcb, dry as f64, u, froude, rw, rv, rt, pe, iff,
                        cw, ct,
                    ])
                    .collect(),
            );
        }
        Ok(rows)
    };

    // Fan out over points, one worker per core. Contiguous chunks mean the
    // results concatenate back into grid order without a sort, and a shared
    // counter reports completion (order is nondeterministic, count is not).
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    let nworkers = cores.min(points).max(1);
    let chunks: Vec<(usize, usize)> = (0..nworkers)
        .map(|w| (w * points / nworkers, (w + 1) * points / nworkers))
        .filter(|(a, b)| b > a)
        .collect();
    let done = std::sync::atomic::AtomicUsize::new(0);
    let done_ref = &done;
    let eval_ref = &eval_point;
    let chunk_results: Vec<Result<Vec<Vec<f64>>, String>> = std::thread::scope(|scope| {
        let handles: Vec<_> = chunks
            .iter()
            .map(|&(a, b)| {
                scope.spawn(move || {
                    let mut local: Vec<Vec<f64>> = Vec::new();
                    for point in a..b {
                        local.extend(eval_ref(point)?);
                        let n = done_ref.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                        eprintln!("point {n}/{points} done");
                    }
                    Ok(local)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    // Serialize in grid order: the chunks are contiguous and already ordered.
    for chunk in chunk_results {
        for nums in chunk? {
            if json {
                if !first_row {
                    out.push(',');
                }
                first_row = false;
                out.push('{');
                for (i, (k, v)) in header.iter().zip(&nums).enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    out.push_str(&format!("{k:?}:{v}"));
                }
                out.push('}');
            } else {
                let row: Vec<String> = nums.iter().map(|v| format!("{v}")).collect();
                out.push_str(&row.join(","));
                out.push('\n');
            }
        }
    }
    if json {
        out.push(']');
        println!("{out}");
    } else {
        print!("{out}");
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Wave field
// ---------------------------------------------------------------------------

/// One `--speed U` or `--froude F` (with `--knots` applying to `--speed`).
fn single_speed(p: &Parsed, l_ref: f64) -> Result<f64, String> {
    let g = p.f64_flag("gravity")?.unwrap_or(STANDARD_GRAVITY);
    match (p.f64_flag("speed")?, p.f64_flag("froude")?) {
        (Some(_), Some(_)) => Err("give either --speed or --froude, not both".into()),
        (Some(u), None) => Ok(if p.switch("--knots") { u * KNOT } else { u }),
        (None, Some(f)) => Ok(f * (g * l_ref).sqrt()),
        (None, None) => Err("select a speed with --speed U or --froude F".into()),
    }
}

fn parse_region(s: &str) -> Result<[f64; 4], String> {
    let parts: Vec<&str> = s.split(':').collect();
    let [x0, x1, y0, y1] = parts.as_slice() else {
        return Err(format!("--region {s:?}: expected X0:X1:Y0:Y1"));
    };
    let mut out = [0.0; 4];
    for (slot, v) in out.iter_mut().zip([x0, x1, y0, y1]) {
        *slot = v
            .trim()
            .parse::<f64>()
            .map_err(|_| format!("--region: cannot parse number {v:?}"))?;
    }
    Ok(out)
}

fn cmd_spectrum(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err("usage: michell spectrum <hull>... --speed U [options]".into());
    }
    let loaded = load_fleet(&p.positional, &p.load_settings()?)?;
    let members: Vec<(&Hull, Placement)> = loaded.iter().map(|m| (&m.hull, m.placement)).collect();
    let l_ref = members
        .iter()
        .map(|(h, _)| h.length())
        .fold(0.0f64, f64::max);
    let u = single_speed(&p, l_ref)?;
    let cond = p.conditions(u)?;
    let n = match p.flag("points") {
        None => 721,
        Some(v) => v
            .parse::<usize>()
            .map_err(|_| format!("--points: cannot parse {v:?}"))?
            .max(9),
    };

    let mut spec = michell::FreeWaveSpectrum::new(&members, &cond).map_err(|e| format!("{e}"))?;

    // Significant angular range: where dRw/dθ still matters.
    let lim = 89.5f64.to_radians();
    let scan = 4096;
    let mut peak = 0.0f64;
    let mut theta_max = 0.0f64;
    for i in 0..=scan {
        let theta = -lim + 2.0 * lim * i as f64 / scan as f64;
        let d = spec.resistance_density(theta);
        if d > peak {
            peak = d;
        }
    }
    for i in 0..=scan {
        let theta = -lim + 2.0 * lim * i as f64 / scan as f64;
        if spec.resistance_density(theta) > 1e-5 * peak {
            theta_max = theta_max.max(theta.abs());
        }
    }
    let theta_max = (theta_max * 1.05).min(lim).max(1e-3);

    // Sample the output rows and the cumulative resistance integral.
    let mut rows = Vec::with_capacity(n);
    let h = 2.0 * theta_max / (n - 1) as f64;
    let mut cum = Vec::with_capacity(n);
    let mut total = 0.0f64;
    let mut prev_d = 0.0f64;
    for i in 0..n {
        let theta = -theta_max + h * i as f64;
        let a = spec.amplitude(theta);
        let d = spec.resistance_density(theta);
        if i > 0 {
            total += 0.5 * (prev_d + d) * h;
        }
        prev_d = d;
        cum.push(total);
        rows.push((theta, a, d));
    }
    let rw_spectrum = total;
    let rw = michell::multihull_wave_resistance(&members, &cond)
        .map_err(|e| format!("{e}"))?
        .resistance;
    eprintln!(
        "U = {u:.3} m/s (Fn {:.3}): Rw = {rw:.4} N (spectrum integral {rw_spectrum:.4} N), \
         transverse wavelength {:.3} m, theta range +-{:.2} deg",
        cond.froude_number(l_ref),
        spec.transverse_wavelength(),
        theta_max.to_degrees()
    );

    if p.switch("--json") {
        let mut out = String::from("{");
        out.push_str(&format!(
            "\"speed\":{u},\"froude\":{},\"k0\":{},\"transverse_wavelength\":{},\
             \"rw_michell\":{rw},\"rw_spectrum\":{rw_spectrum},\"points\":[",
            cond.froude_number(l_ref),
            spec.wavenumber(),
            spec.transverse_wavelength()
        ));
        for (i, (theta, a, d)) in rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"theta_deg\":{},\"lambda\":{},\"wavelength\":{},\"amp_re\":{},\
                 \"amp_im\":{},\"amp_abs\":{},\"drw_dtheta\":{},\"cum_fraction\":{}}}",
                theta.to_degrees(),
                1.0 / theta.cos(),
                spec.transverse_wavelength() * theta.cos() * theta.cos(),
                a.re,
                a.im,
                a.abs(),
                d,
                if rw_spectrum > 0.0 {
                    cum[i] / rw_spectrum
                } else {
                    0.0
                }
            ));
        }
        out.push_str("]}");
        println!("{out}");
        return Ok(());
    }

    println!("theta_deg,lambda,wavelength_m,amp_re,amp_im,amp_abs,drw_dtheta,cum_fraction");
    for (i, (theta, a, d)) in rows.iter().enumerate() {
        println!(
            "{},{},{},{},{},{},{},{}",
            theta.to_degrees(),
            1.0 / theta.cos(),
            spec.transverse_wavelength() * theta.cos() * theta.cos(),
            a.re,
            a.im,
            a.abs(),
            d,
            if rw_spectrum > 0.0 {
                cum[i] / rw_spectrum
            } else {
                0.0
            }
        );
    }
    Ok(())
}

fn cmd_wake(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err("usage: michell wake <hull>... --speed U [-o wake.png] [options]".into());
    }
    let loaded = load_fleet(&p.positional, &p.load_settings()?)?;
    let members: Vec<(&Hull, Placement)> = loaded.iter().map(|m| (&m.hull, m.placement)).collect();
    let l_ref = members
        .iter()
        .map(|(h, _)| h.length())
        .fold(0.0f64, f64::max);
    let u = single_speed(&p, l_ref)?;
    let cond = p.conditions(u)?;

    // Fleet extents in fleet coordinates.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_abs = 0.0f64;
    for m in &loaded {
        let (h0, h1) = m.hull.surface().x_domain();
        x_lo = x_lo.min(h0 + m.placement.x);
        x_hi = x_hi.max(h1 + m.placement.x);
        y_abs = y_abs.max(m.placement.y.abs());
    }

    let [x0, x1, y0, y1] = match p.flag("region") {
        Some(s) => parse_region(s)?,
        None => {
            let x1 = x_hi + 0.35 * l_ref;
            let x0 = x_lo - 3.0 * l_ref;
            let yh = (0.42 * (x1 - x0)).max(y_abs + 0.8 * l_ref);
            [x0, x1, -yh, yh]
        }
    };
    let (nx, ny) = match p.flag("size") {
        Some(s) => parse_pair(s)?,
        None => {
            let nx = 900usize;
            let ny = ((nx as f64) * (y1 - y0) / (x1 - x0)).round() as usize;
            (nx, ny.clamp(64, 1200))
        }
    };
    if nx < 2 || ny < 2 {
        return Err("--size: need at least 2x2 grid points".into());
    }

    let mut spec = michell::FreeWaveSpectrum::new(&members, &cond).map_err(|e| format!("{e}"))?;
    let grid = spec
        .elevation_grid(x0, x1, y0, y1, nx, ny)
        .map_err(|e| format!("{e}"))?;

    let out_path = p.flag("output").cloned().unwrap_or_else(|| {
        if p.switch("--json") || p.switch("--csv") {
            String::new()
        } else {
            "wake.png".to_string()
        }
    });

    let mut lo = 0.0f64;
    let mut hi = 0.0f64;
    for &v in &grid.zeta {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    let summary = format!(
        "wake: {nx}x{ny} over x {x0:.2}..{x1:.2} m, y {y0:.2}..{y1:.2} m at U = {u:.3} m/s \
         (Fn {:.3})\nzeta {:.4}..{:.4} m, transverse wavelength {:.3} m, max lambda {:.1}{}\n\
         ship advances toward +x; the pattern is physical astern of each hull",
        cond.froude_number(l_ref),
        lo,
        hi,
        spec.transverse_wavelength(),
        grid.max_lambda,
        if grid.resolution_limited {
            " (grid-resolution limited; finer --size reveals shorter diverging waves)"
        } else {
            ""
        }
    );

    // Data outputs.
    let json_out = out_path.ends_with(".json") || (out_path.is_empty() && p.switch("--json"));
    let csv_out = out_path.ends_with(".csv") || (out_path.is_empty() && p.switch("--csv"));
    if json_out || csv_out {
        let mut out = String::new();
        if json_out {
            out.push_str(&format!(
                "{{\"speed\":{u},\"x0\":{x0},\"x1\":{x1},\"y0\":{y0},\"y1\":{y1},\
                 \"nx\":{nx},\"ny\":{ny},\"transverse_wavelength\":{},\"max_lambda\":{},\
                 \"resolution_limited\":{},\"zeta\":[",
                spec.transverse_wavelength(),
                grid.max_lambda,
                grid.resolution_limited
            ));
            for iy in 0..ny {
                if iy > 0 {
                    out.push(',');
                }
                out.push('[');
                for ix in 0..nx {
                    if ix > 0 {
                        out.push(',');
                    }
                    out.push_str(&format!("{}", grid.get(ix, iy)));
                }
                out.push(']');
            }
            out.push_str("]}");
        } else {
            out.push_str("x,y,zeta\n");
            for iy in 0..ny {
                for ix in 0..nx {
                    out.push_str(&format!(
                        "{},{},{}\n",
                        grid.x(ix),
                        grid.y(iy),
                        grid.get(ix, iy)
                    ));
                }
            }
        }
        eprintln!("{summary}");
        if out_path.is_empty() {
            print!("{out}");
            if json_out {
                println!();
            }
        } else {
            std::fs::write(&out_path, out).map_err(|e| format!("cannot write {out_path}: {e}"))?;
            println!("wrote {out_path}");
        }
        return Ok(());
    }
    if !out_path.ends_with(".png") {
        return Err(format!(
            "wake output {out_path:?}: expected a .png, .csv, or .json path"
        ));
    }

    // Color scale: saturate at the 99.5th percentile of |zeta| so a single
    // extreme pixel does not wash out the pattern.
    let vmax = match p.f64_flag("zmax")? {
        Some(v) if v > 0.0 => v,
        Some(v) => return Err(format!("--zmax must be positive, got {v}")),
        None => {
            let mut abs: Vec<f64> = grid.zeta.iter().map(|v| v.abs()).collect();
            abs.sort_by(|a, b| a.total_cmp(b));
            abs[((abs.len() - 1) as f64 * 0.995) as usize].max(1e-12)
        }
    };

    // The free-wave field is physical only astern of every hull: fade the
    // columns from the aft-most stern forward so the trustworthy wake reads
    // at full strength and the rest is visibly indicative.
    const FADE_TOWARD: [u8; 3] = [0xf0, 0xef, 0xec];
    const FADE_FRACTION: f64 = 0.55;
    let mut rgb = vec![0u8; 3 * nx * ny];
    for iy in 0..ny {
        // PNG row 0 is the top of the image = the +y edge.
        let row = ny - 1 - iy;
        for ix in 0..nx {
            let t = grid.get(ix, iy) / vmax;
            let mut c = png::diverging(t);
            if grid.x(ix) > x_lo {
                c = png::fade(c, FADE_TOWARD, FADE_FRACTION);
            }
            rgb[3 * (row * nx + ix)..3 * (row * nx + ix) + 3].copy_from_slice(&c);
        }
    }
    // Hull waterplane footprints in neutral dark gray.
    const HULL_GRAY: [u8; 3] = [0x52, 0x51, 0x4e];
    for m in &loaded {
        let (h0, h1) = m.hull.surface().x_domain();
        for ix in 0..nx {
            let x = grid.x(ix) - m.placement.x;
            if x < h0 || x > h1 {
                continue;
            }
            let half_beam = m.hull.surface().eval(x, 0.0);
            for iy in 0..ny {
                if (grid.y(iy) - m.placement.y).abs() <= half_beam {
                    let row = ny - 1 - iy;
                    rgb[3 * (row * nx + ix)..3 * (row * nx + ix) + 3].copy_from_slice(&HULL_GRAY);
                }
            }
        }
    }
    std::fs::write(&out_path, png::encode_rgb(nx, ny, &rgb))
        .map_err(|e| format!("cannot write {out_path}: {e}"))?;
    println!("{summary}");
    if x1 > x_lo {
        println!("faded ahead of x = {x_lo:.2} m (aft-most stern): not physical there");
    }
    println!("color scale: +-{vmax:.4} m; wrote {out_path}");
    Ok(())
}

fn cmd_render(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err(
            "usage: michell render <hull>... --speed U [-o render.png] [options]".into(),
        );
    }
    let loaded = load_fleet(&p.positional, &p.load_settings()?)?;
    let members: Vec<(&Hull, Placement)> =
        loaded.iter().map(|m| (&m.hull, m.placement)).collect();
    let l_ref = members
        .iter()
        .map(|(h, _)| h.length())
        .fold(0.0f64, f64::max);
    let u = single_speed(&p, l_ref)?;
    let cond = p.conditions(u)?;

    // Fleet extents in fleet coordinates.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_abs = 0.0f64;
    for m in &loaded {
        let (h0, h1) = m.hull.surface().x_domain();
        x_lo = x_lo.min(h0 + m.placement.x);
        x_hi = x_hi.max(h1 + m.placement.x);
        y_abs = y_abs.max(m.placement.y.abs());
    }

    let [x0, x1, y0, y1] = match p.flag("region") {
        Some(s) => parse_region(s)?,
        None => {
            let x1 = x_hi + 0.35 * l_ref;
            let x0 = x_lo - 3.0 * l_ref;
            let yh = (0.42 * (x1 - x0)).max(y_abs + 0.8 * l_ref);
            [x0, x1, -yh, yh]
        }
    };
    let (gx, gy) = {
        let gx = match p.flag("grid") {
            Some(s) => s
                .parse::<usize>()
                .map_err(|_| format!("--grid: cannot parse count {s:?}"))?,
            None => 700,
        };
        let gy = ((gx as f64) * (y1 - y0) / (x1 - x0)).round() as usize;
        (gx.max(2), gy.clamp(64, 1200))
    };
    let (iw, ih) = match p.flag("size") {
        Some(s) => parse_pair(s)?,
        None => (1200, 800),
    };
    if iw < 16 || ih < 16 {
        return Err("--size: image must be at least 16x16 pixels".into());
    }

    let mut spec =
        michell::FreeWaveSpectrum::new(&members, &cond).map_err(|e| format!("{e}"))?;
    let grid = spec
        .elevation_grid(x0, x1, y0, y1, gx, gy)
        .map_err(|e| format!("{e}"))?;

    let vmax = match p.f64_flag("zmax")? {
        Some(v) if v > 0.0 => v,
        Some(v) => return Err(format!("--zmax must be positive, got {v}")),
        None => {
            let mut abs: Vec<f64> = grid.zeta.iter().map(|v| v.abs()).collect();
            abs.sort_by(|a, b| a.total_cmp(b));
            abs[((abs.len() - 1) as f64 * 0.995) as usize].max(1e-12)
        }
    };
    let z_scale = match p.f64_flag("z-scale")? {
        Some(s) if s > 0.0 => s,
        Some(s) => return Err(format!("--z-scale must be positive, got {s}")),
        None => 1.0,
    };

    let mut scene = render::Scene::default();
    let fade = render::PhysFade {
        x_phys: x_lo,
        feather: 0.5 * l_ref,
    };
    render::add_water(&mut scene, &grid, z_scale, vmax, fade);
    render::add_skirt(&mut scene, &grid, z_scale, vmax, fade);
    for m in &loaded {
        let fb = match p.f64_flag("freeboard")? {
            Some(f) if f >= 0.0 => f,
            Some(f) => return Err(format!("--freeboard must be non-negative, got {f}")),
            // Slender hulls (amas) have tiny drafts; keep a visible sheer.
            None => (0.5 * m.hull.draft()).max(0.02 * m.hull.length()),
        };
        render::add_hull(
            &mut scene,
            m.hull.surface(),
            m.placement.x,
            m.placement.y,
            fb,
        );
    }

    // Camera: degrees off dead astern (positive toward +y), elevation, and
    // an optional distance in metres (default frames the fleet).
    let (az, el, dist) = match p.flag("camera") {
        Some(s) => {
            let parts: Vec<&str> = s.split(':').collect();
            if parts.len() < 2 || parts.len() > 3 {
                return Err(format!("--camera {s:?}: expected AZ:EL[:DIST]"));
            }
            let mut nums = [0.0f64; 3];
            for (slot, v) in nums.iter_mut().zip(&parts) {
                *slot = v
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("--camera: cannot parse number {v:?}"))?;
            }
            let dist = if parts.len() == 3 { Some(nums[2]) } else { None };
            (nums[0], nums[1], dist)
        }
        None => (35.0, 18.0, None),
    };
    let l_fleet = x_hi - x_lo;
    let dist = dist.unwrap_or(2.1 * (l_fleet + y_abs));
    let target = [0.5 * (x_lo + x_hi), 0.0, 0.0];
    let (azr, elr) = (az.to_radians(), el.to_radians());
    let eye = render::add(
        target,
        [
            -dist * azr.cos() * elr.cos(),
            dist * azr.sin() * elr.cos(),
            dist * elr.sin(),
        ],
    );
    let cam = render::Camera {
        eye,
        target,
        fov_deg: 32.0,
    };
    let light = render::Light {
        dir: render::normalize([0.4, -0.4, 0.83]),
        ambient: 0.45,
        diffuse: 0.6,
    };

    let rgb = render::render(&scene, &cam, &light, iw, ih);
    let out_path = p
        .flag("output")
        .cloned()
        .unwrap_or_else(|| "render.png".to_string());
    if !out_path.ends_with(".png") {
        return Err(format!("render output {out_path:?}: expected a .png path"));
    }
    std::fs::write(&out_path, png::encode_rgb(iw, ih, &rgb))
        .map_err(|e| format!("cannot write {out_path}: {e}"))?;
    println!(
        "render: {iw}x{ih} px, water {gx}x{gy} over x {x0:.2}..{x1:.2} m, y {y0:.2}..{y1:.2} m \
         at U = {u:.3} m/s (Fn {:.3})\ncamera {az:.0} deg off astern, {el:.0} deg up, {dist:.1} m; \
         tint +-{vmax:.4} m{}",
        cond.froude_number(l_ref),
        if z_scale != 1.0 {
            format!(", water z x{z_scale}")
        } else {
            String::new()
        }
    );
    println!(
        "free-wave surface: physical astern of each hull; paled ahead of x = {x_lo:.2} m; \
         hulls at static waterline (no sinkage/trim); wrote {out_path}"
    );
    Ok(())
}

fn cmd_loft(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    let [path] = p.positional.as_slice() else {
        return Err("usage: michell loft <offsets|iges> -o OUT.hull [options]".into());
    };
    let out_path = p
        .flag("output")
        .ok_or("loft requires an output path: -o PREFIX (or OUT.hull)")?;
    let settings = p.load_settings()?;

    // IGES and STL files loft to full-band bodies by default, decomposing
    // multihull files into one body per hull. --wetted keeps the old
    // single-hull wetted-only output; offsets tables always loft wetted.
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let lower = path.to_ascii_lowercase();
    let is_text = std::str::from_utf8(&bytes).is_ok();
    let is_stl = lower.ends_with(".stl") || formats::looks_binary_stl(&bytes) || !is_text;
    let first = if is_stl {
        ""
    } else {
        std::str::from_utf8(&bytes)
            .expect("checked utf8")
            .lines()
            .find(|l| !l.trim().is_empty())
            .unwrap_or("")
            .trim_end()
    };
    let is_iges = !is_stl
        && (lower.ends_with(".igs")
            || lower.ends_with(".iges")
            || first.len() >= 73 && matches!(first.as_bytes()[72], b'S' | b'G'));

    if !(is_iges || is_stl) || p.switch("--wetted") {
        let mut hulls = load_hulls(path, &settings)?;
        if hulls.len() > 1 {
            return Err(format!(
                "{path} contains {} hulls; --wetted lofts exactly one — drop \
                 --wetted to decompose into full-band bodies",
                hulls.len()
            ));
        }
        let (hull, _, source) = hulls.pop().expect("one hull");
        if matches!(source, Source::Native | Source::Body(_)) {
            return Err(format!("{path} is already a control-net file"));
        }
        std::fs::write(out_path, write_hull_file(&hull))
            .map_err(|e| format!("cannot write {out_path}: {e}"))?;
        for line in describe_source(&source) {
            println!("{line}");
        }
        println!("wrote {out_path}");
        return Ok(());
    }

    // Full-band decomposition. A body is lofted once and re-situated many
    // times, so default to a much denser sampling and control net than the
    // one-shot import path — wave resistance is sensitive to loft resolution
    // near support boundaries (keel rocker, stem), and the band is taller
    // than the wetted zone.
    use michell::iges::{self, HullPose, ImportOptions, Platform};
    let design_wl = settings.waterline_z;
    let opts = ImportOptions {
        waterline_z: design_wl,
        stations: if p.flag("samples").is_some() {
            settings.samples.0
        } else {
            241
        },
        waterlines: if p.flag("samples").is_some() {
            settings.samples.1
        } else {
            97
        },
        fit: if settings.fit_explicit {
            settings.fit
        } else {
            michell::fit::FitOptions {
                degree_x: 3,
                degree_z: 3,
                n_ctrl_x: 28,
                n_ctrl_z: 32,
                // Honour --fit-deriv-weight even with the default net.
                derivative_weight: settings.fit.derivative_weight,
            }
        },
        centerplane: settings.centerplane,
    };
    enum LoftSrc {
        Iges(michell::iges::SourceFleet),
        Mesh(michell::stl::MeshFleet),
    }
    impl LoftSrc {
        fn len(&self) -> usize {
            match self {
                LoftSrc::Iges(s) => s.len(),
                LoftSrc::Mesh(s) => s.len(),
            }
        }
        fn hull_z_top(&self, i: usize) -> f64 {
            match self {
                LoftSrc::Iges(s) => s.hull_z_top(i),
                LoftSrc::Mesh(s) => s.hull_z_top(i),
            }
        }
        fn hull_z_bottom(&self, i: usize) -> f64 {
            match self {
                LoftSrc::Iges(s) => s.hull_z_bottom(i),
                LoftSrc::Mesh(s) => s.hull_z_bottom(i),
            }
        }
        fn situate_one(
            &self,
            i: usize,
            wl: f64,
            pose: &HullPose,
            platform: &Platform,
            opts: &ImportOptions,
        ) -> Result<Option<michell::iges::ImportedHull>, michell::Error> {
            match self {
                LoftSrc::Iges(s) => s.situate_one(i, wl, pose, platform, opts),
                LoftSrc::Mesh(s) => s.situate_one(i, wl, pose, platform, opts),
            }
        }
    }
    let src = if is_stl {
        let scale = settings.units.ok_or_else(|| {
            format!(
                "{path}: STL files carry no units; pass --units mm|cm|m|in|ft \
                 (or a scale to metres)"
            )
        })?;
        LoftSrc::Mesh(
            michell::stl::mesh_fleet(&bytes, scale, design_wl)
                .map_err(|e| format!("{path}: {e}"))?,
        )
    } else {
        let text = std::str::from_utf8(&bytes).expect("checked utf8");
        LoftSrc::Iges(iges::source_fleet(text, design_wl).map_err(|e| format!("{path}: {e}"))?)
    };
    let n = src.len();
    // The band reaches from the keel to `--band` metres above the design
    // waterline (default: half the design draft). Including the deck in the
    // fit would distort the wetted geometry — a deck is a cliff for a
    // height-field loft — so the band should stay below it; a sweep that
    // rises past the band is reported per-pose as band_exceeded.
    let band_flag = p.f64_flag("band")?;
    let mut lofted = Vec::new();
    for idx in 0..n {
        let top = src.hull_z_top(idx);
        let bottom = src.hull_z_bottom(idx);
        let draft_est = design_wl - bottom;
        if draft_est <= 0.0 {
            return Err(format!(
                "{path} hull {idx}: design waterline (z = {design_wl}) is below \
                 the hull (keel bound z = {bottom:.3})"
            ));
        }
        let margin = band_flag.unwrap_or(0.5 * draft_est).max(0.0);
        let band_top = (design_wl + margin).min(top);
        // Detect the centerplane at the *design* waterline, where the hull is
        // symmetric and the fold physically matters — band-top probes sample
        // topsides, where fittings can skew the detection — then hold it
        // fixed for the band loft.
        let mut detect_opts = opts;
        detect_opts.stations = 61;
        detect_opts.waterlines = 17;
        detect_opts.fit.n_ctrl_x = detect_opts.fit.n_ctrl_x.min(10);
        detect_opts.fit.n_ctrl_z = detect_opts.fit.n_ctrl_z.min(7);
        let mut hull_opts = opts;
        if hull_opts.centerplane.is_none() {
            hull_opts.centerplane = src
                .situate_one(
                    idx,
                    design_wl,
                    &HullPose::default(),
                    &Platform::default(),
                    &detect_opts,
                )
                .map_err(|e| format!("{path} hull {idx}: {e}"))?
                .map(|m| m.report.centerplane);
        }
        let m = src
            .situate_one(
                idx,
                band_top,
                &HullPose::default(),
                &Platform::default(),
                &hull_opts,
            )
            .map_err(|e| format!("{path} hull {idx}: {e}"))?
            .ok_or_else(|| format!("{path} hull {idx}: nothing below the band top?"))?;
        let wl_depth = band_top - design_wl;
        lofted.push((m, wl_depth));
    }
    let grids: Vec<(&michell::SampleGrid, Option<f64>)> = lofted
        .iter()
        .map(|(m, _)| (&m.grid, Some(m.report.centerplane)))
        .collect();
    formats::dump_grids(&settings, &grids)?;
    let ys: Vec<f64> = lofted.iter().map(|(m, _)| m.report.centerplane).collect();
    let span = ys.last().unwrap_or(&0.0) - ys.first().unwrap_or(&0.0);
    let names: Vec<String> = if n == 1 {
        vec![String::new()]
    } else if n == 2 && (ys[0] + ys[1]).abs() < 0.1 * span {
        vec!["port".into(), "starboard".into()]
    } else if n == 3
        && (ys[0] + ys[2]).abs() < 0.1 * span
        && (ys[1] - (ys[0] + ys[2]) / 2.0).abs() < 0.25 * span
    {
        vec!["port".into(), "center".into(), "starboard".into()]
    } else {
        (0..n).map(|i| i.to_string()).collect()
    };
    let prefix = out_path.strip_suffix(".hull").unwrap_or(out_path);
    println!(
        "{:<28} {:>12} {:>9} {:>9} {:>11} {:>18}",
        "file", "centerplane", "band[m]", "WL depth", "residual", "at (x, z)"
    );
    for ((m, wl_depth), name) in lofted.iter().zip(&names) {
        let file = if name.is_empty() {
            format!("{prefix}.hull")
        } else {
            format!("{prefix}-{name}.hull")
        };
        let body_text = formats::write_body_file(m.hull.surface(), *wl_depth, m.report.centerplane);
        std::fs::write(&file, body_text).map_err(|e| format!("cannot write {file}: {e}"))?;
        println!(
            "{file:<28} {:>12.4} {:>9.4} {:>9.4} {:>11.3e} {:>9.3},{:>7.3}",
            m.report.centerplane,
            m.hull.draft(),
            wl_depth,
            m.report.fit.max_residual,
            m.report.fit.max_residual_at.0,
            m.report.fit.max_residual_at.1,
        );
    }
    Ok(())
}

/// One `place` input: a path plus the pose to apply to every hull in it.
struct PlaceSpec {
    path: String,
    /// Absolute centerplane (`y=`); `.hull` inputs only.
    y_abs: Option<f64>,
    pose: michell::iges::HullPose,
}

/// Parse `path` or `path@key=V,...` (keys: dx/x, dy, y, dz, trim [deg],
/// pivot [m]).
fn parse_place_spec(spec: &str) -> Result<PlaceSpec, String> {
    let (path, rest) = match spec.split_once('@') {
        Some((p, r)) => (p, r),
        None => (spec, ""),
    };
    let mut out = PlaceSpec {
        path: path.to_string(),
        y_abs: None,
        pose: michell::iges::HullPose::default(),
    };
    for part in rest.split(',').filter(|s| !s.trim().is_empty()) {
        let (k, v) = part
            .split_once('=')
            .ok_or_else(|| format!("bad pose {rest:?}: expected key=value pairs"))?;
        let val: f64 = v
            .trim()
            .parse()
            .map_err(|_| format!("bad pose value {v:?} in {spec:?}"))?;
        match k.trim() {
            "dx" | "x" => out.pose.dx = val,
            "dy" => out.pose.dy = val,
            "y" => out.y_abs = Some(val),
            "dz" => out.pose.dz = val,
            "trim" => out.pose.trim = val.to_radians(),
            "pivot" => out.pose.pivot_x = Some(val),
            other => {
                return Err(format!(
                    "unknown pose key {other:?} (use dx/x, dy, y, dz, trim, pivot)"
                ))
            }
        }
    }
    if out.y_abs.is_some() && out.pose.dy != 0.0 {
        return Err(format!(
            "{spec:?}: give either y (absolute) or dy (shift), not both"
        ));
    }
    Ok(out)
}

fn cmd_place(args: &[String]) -> Result<(), String> {
    use michell::iges::{self, Platform};
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err(
            "usage: michell place <hull>[@dx=..,dy=..,dz=..,trim=..]... -o OUT.igs \
             [--waterline Z] [--sinkage S] [--platform-trim DEG] [--pivot-x X]"
                .into(),
        );
    }
    let out_path = p
        .flag("output")
        .ok_or("place requires an output path: -o OUT.igs")?;
    let waterline_z = p.f64_flag("waterline")?.unwrap_or(0.0);
    let platform = Platform {
        sinkage: p.f64_flag("sinkage")?.unwrap_or(0.0),
        trim: p.f64_flag("platform-trim")?.unwrap_or(0.0).to_radians(),
        pivot_x: p.f64_flag("pivot-x")?.unwrap_or(0.0),
    };

    // Each input file is parsed once; a spec's pose applies rigidly to every
    // hull the file contains.
    enum Loaded {
        Fleet(iges::SourceFleet),
        Spline(formats::HullFileData),
    }
    let mut cache: HashMap<String, Loaded> = HashMap::new();
    let mut surfaces: Vec<michell::iges::NurbsSurface3> = Vec::new();
    for raw in &p.positional {
        let spec = parse_place_spec(raw)?;
        if !cache.contains_key(&spec.path) {
            let path = &spec.path;
            let text =
                std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
            let first = text
                .lines()
                .find(|l| !l.trim().is_empty())
                .unwrap_or("")
                .trim_end();
            let lower = path.to_ascii_lowercase();
            let loaded = if first.starts_with("michell-hull") {
                Loaded::Spline(formats::parse_hull_data(&text).map_err(|e| format!("{path}: {e}"))?)
            } else if lower.ends_with(".igs")
                || lower.ends_with(".iges")
                || first.len() >= 73 && matches!(first.as_bytes()[72], b'S' | b'G')
            {
                Loaded::Fleet(
                    iges::source_fleet(&text, waterline_z).map_err(|e| format!("{path}: {e}"))?,
                )
            } else {
                return Err(format!(
                    "{path}: place takes IGES files and .hull control nets or \
                     bodies; loft other formats first (michell loft)"
                ));
            };
            cache.insert(spec.path.clone(), loaded);
        }
        let before = surfaces.len();
        let hulls = match &cache[&spec.path] {
            Loaded::Fleet(fleet) => {
                if spec.y_abs.is_some() {
                    return Err(format!(
                        "{raw:?}: absolute y placement needs a recorded centerplane, \
                         which IGES inputs don't carry; use dy=SHIFT"
                    ));
                }
                for i in 0..fleet.len() {
                    surfaces.extend(
                        fleet
                            .posed_surfaces(i, waterline_z, &spec.pose, &platform)
                            .map_err(|e| format!("{}: {e}", spec.path))?,
                    );
                }
                fleet.len()
            }
            Loaded::Spline(data) => {
                let centerplane = data.centerplane.unwrap_or(0.0);
                let z_top_cad = waterline_z + data.waterline.unwrap_or(0.0);
                let mut pose = spec.pose;
                if let Some(y) = spec.y_abs {
                    pose.dy = y - centerplane;
                }
                let mut pair = iges::halfbreadth_surfaces(&data.surface, centerplane, z_top_cad);
                iges::apply_pose(&mut pair, waterline_z, &pose, &platform);
                surfaces.extend(pair);
                1
            }
        };
        eprintln!(
            "{}: {} hull{}, {} patch{}",
            raw,
            hulls,
            if hulls == 1 { "" } else { "s" },
            surfaces.len() - before,
            if surfaces.len() - before == 1 {
                ""
            } else {
                "es"
            },
        );
    }

    let stem = std::path::Path::new(out_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("michell");
    let text = iges::write(&surfaces, stem).map_err(|e| format!("{e}"))?;
    std::fs::write(out_path, text).map_err(|e| format!("cannot write {out_path}: {e}"))?;
    println!(
        "wrote {out_path}: {} surface patches, metres, design waterline at z = {waterline_z}",
        surfaces.len()
    );
    Ok(())
}

fn cmd_wigley(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    let l = p.f64_flag("length")?.unwrap_or(10.0);
    let b = p.f64_flag("beam")?.unwrap_or(l / 10.0);
    let t = p.f64_flag("draft")?.unwrap_or(b * 0.625);
    let hull = michell::hulls::wigley(l, b, t).map_err(|e| format!("{e}"))?;
    let text = write_hull_file(&hull);
    match p.flag("output") {
        Some(path) => {
            std::fs::write(path, text).map_err(|e| format!("cannot write {path}: {e}"))?;
            println!("wrote {path} (Wigley L={l} B={b} T={t})");
        }
        None => print!("{text}"),
    }
    Ok(())
}

//! `michell` — thin-ship wave resistance from hull files.
//!
//! The CLI lives here as a library so front-ends (the web server) can load
//! hulls and run sweeps in-process rather than shelling out to the `michell`
//! binary. [`info`] and [`run_manifest`] report progress through a
//! [`Reporter`] and return the text they would otherwise have printed to
//! stdout; every other command still prints directly.

mod archive;
pub mod fleet;
mod formats;
mod json;
pub mod manifest;
mod pdf;
mod png;
mod render;
mod report;
mod scene;
mod seakeeping;

use fleet::{describe, max_beam, LoadSettings};
pub use formats::parse_units;
use formats::{parse_pair, parse_range};
use michell::{Conditions, Fluid, Placement, WaveOptions, STANDARD_GRAVITY};
use std::collections::HashMap;
use std::fmt::Write as _;

const KNOT: f64 = 1852.0 / 3600.0; // m/s

/// A progress sink. Each call carries one human-readable line (what the CLI
/// writes to stderr) and, when the line marks progress, a machine-readable
/// `(done, total)` fraction a front-end can draw as a bar.
pub type Reporter<'a> = dyn FnMut(&str, Option<(usize, usize)>) + 'a;

/// Preferred sweep entry point: run a JSON manifest, returning the CSV/JSON it
/// would have printed to stdout (empty when the manifest names an output file,
/// which is written directly). Progress flows through `report`.
pub fn run_manifest(manifest_path: &str, report: &mut Reporter) -> Result<String, String> {
    manifest::run(manifest_path, report)
}

/// Dispatch a full CLI argument vector (without the program name), reporting
/// progress through `report`. Commands that produce stdout do so directly here;
/// the in-process entry points build their output as a string which this prints.
pub fn run(args: &[String], report: &mut Reporter) -> Result<(), String> {
    match args.first().map(String::as_str) {
        Some("resistance") => cmd_resistance(&args[1..]),
        Some("report") => cmd_report(&args[1..]),
        Some("squat") => cmd_squat(&args[1..]),
        Some("seakeeping") => seakeeping::cmd_seakeeping(&args[1..]),
        Some("sweep") => {
            // A JSON manifest is the preferred sweep interface.
            if let Some(path) = args.get(1).filter(|a| a.ends_with(".json")) {
                if args.len() > 2 {
                    return Err("manifest sweeps take no further arguments; put \
                                everything in the JSON file"
                        .into());
                }
                let out = manifest::run(path, report)?;
                print!("{out}");
                Ok(())
            } else {
                cmd_sweep(&args[1..])
            }
        }
        Some("info") => {
            let out = info(&args[1..])?;
            print!("{out}");
            Ok(())
        }
        Some("spectrum") => cmd_spectrum(&args[1..]),
        Some("wake") => cmd_wake(&args[1..]),
        Some("field") => cmd_field(&args[1..]),
        Some("render") => cmd_render(&args[1..]),
        Some("view") => Err("`michell view` was removed; use the web viewer (michell-web)".into()),
        Some("loft") => Err("`michell loft` was removed: every command cuts hulls into \
                             sections straight from the IGES or STL file"
            .into()),
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
  michell squat <hull>... --speeds A[:B:STEP] [options]       dynamic sinkage/trim force
  michell seakeeping <hull>... --froude F [options]           heave/pitch RAOs, added resistance
  michell info <hull>... [options]                            geometry & diagnostics
  michell spectrum <hull>... --speed U [options]              free-wave spectrum
  michell wake <hull>... --speed U [-o wake.png] [options]    Kelvin wake heatmap
  michell field <hull>... --speed U -o scene.json [options]   3-D scene: sections, closure, wave field
  michell place <hull>[@dx=..,dy=..,dz=..]... -o OUT.igs      write posed CAD geometry
  michell wigley [-o OUT.igs] [--length L --beam B --draft T]   exact Wigley hull as IGES

HULL INPUTS (sniffed by header / extension)
  Every hull is cut into sections (stations along x, each integrated along
  rays from its top centreplane point) straight from its source geometry,
  and re-cut at every pose a sweep or equilibrium solve visits.
  *.igs, *.iges     NURBS surfaces (types 128/143/141), clustered into hulls
  *.stl             triangle mesh, binary or ASCII; requires --units

SEAKEEPING
  michell seakeeping <hull>... (--speed U | --froude F) [--heading DEG]
      [--lambda A:B:STEP] [--kyy FRAC] [--mass KG] [--lcg X] [--panels N]
      [--sea hs=H,tp=T[,gamma=G]]
  Heave and pitch RAOs and added resistance by strip theory (Salvesen–
  Tuck–Faltinsen; Gerritsma–Beukelman) over wave lengths λ/L (default
  0.5:3:0.125), head seas (180) by default; mass defaults to the
  displacement at the loaded waterline with the LCG over the LCB, and
  k_yy to 0.25 L. --sea adds significant motions, bow and LCG vertical
  accelerations and mean added resistance in a Bretschneider (or, with
  gamma, JONSWAP) sea. Multihulls move as one rigid platform, without
  hull-to-hull wave interaction. --dynamic first floats the platform at
  its thin-ship dynamic sinkage and trim at that speed, and takes the
  motions about that attitude.

MULTIHULLS
  Pass several hulls; each may carry a placement suffix:
      michell resistance vaka.igs ama.igs@y=1.9 ama.igs@y=-1.9 --speeds 3:8:0.5
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

SQUAT (thin-ship dynamic sinkage and trim)
  michell squat <hull>... --speeds A[:B:STEP] [--knots] [--pivot X]
  The near-field pressure's vertical force Fz (+ up) and pitch moment M
  (+ bow-up, about --pivot at the waterline; default the fleet's LCF) on
  the hulls held at their current attitude, from the same exact hull
  transforms as the wave integral (sinkage is the local field; trim is
  mostly the wave part). Reports the lift fraction Fz/(rho g V) and the
  first-order equivalent sinkage -Fz/(rho g Aw) and trim M/(rho g I_L),
  I_L about the LCF. Thin-ship overstates both ~20-40% vs surface-panel
  linear theory at Fn 0.3-0.4; the linearisation expires once |lift| is a
  real share of the weight.

SPEED SELECTION (resistance, squat)
  --speeds A[:B:STEP]   speeds in m/s (inclusive range)
  --froude A[:B:STEP]   length Froude numbers instead of speeds
  --knots               interpret and display speeds in knots

IMPORT OPTIONS
  --waterline Z         IGES/STL: design waterline height in the file frame,
                        metres after unit conversion (z up; default 0)
  --centerplane Y       transverse position of the hull centerplane
                        (default: auto-detect; full shells fold about their
                        midplane, half hulls measure from y = 0)
  --units U             STL: scale to metres (mm|cm|m|in|ft or a number);
                        required for STL, which has no units field
  --stations N          stations along each hull (default 121, cosine-spaced)
  --rays M              rays across each section (default 33)

PHYSICS OPTIONS
  --fluid NAME          seawater | freshwater, at 15 C (default seawater)
  --rho R, --nu V       override density [kg/m3] / kinematic viscosity [m2/s]
  --gravity G           override g [m/s2]
  --form-factor K       form factor k, applied as (1+k)*C_F; default 0 (a
                        bare flat plate). A property of the SHAPE: streamline
                        curvature and the stern's viscous pressure defect
  --roughness SPEC      surface-roughness allowance, added OUTSIDE (1+k) as
                        C_V = (1+k)*C_F + dC_F, per ITTC-78. A property of the
                        SKIN, so it is kept separate from the form factor:
                          off                 hydraulically smooth (default)
                          cf=DELTA_CF         a prescribed dC_F you trust
                          ks=HEIGHT           equivalent sand-grain height
                                              (100um, 0.1mm, 1e-4); dC_F is
                                              then speed-dependent
                        Guide: sprayed topcoat ~30um, rolled antifouling
                        ~100-150um, light slime a few hundred um. On a small
                        hull this term can exceed the form factor
  --rel-tol T           wave-integral relative tolerance (default 1e-5)
  --transom SPEC        transom-stern closure: off | ballistic[=COEFF] |
                        hollow=METRES. A transom leaves the half-breadth open
                        at the stern; the closure appends a virtual hollow of
                        length L_v = COEFF*U*sqrt(d_T/g) (default COEFF=sqrt2,
                        the ballistic free-fall value) so the body closes.
                        Inert on a hull that closes aft

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
    --zmax M            tint saturation elevation [m] (default: 99.5th pct)
  Same caveat as wake (paled where not astern of every hull); hulls are drawn
  from their source geometry at the static waterline, without dynamic
  sinkage or trim.

SWEEPS
  michell sweep study.json      preferred: a JSON manifest referencing hull
                                files (IGES/STL; `hull: N` picks one
                                of a multihull file's hulls), with
                                speed/waterline axes, per-hull load
                                (mass/lcg/vcg), point loads (mass at an
                                offset from a hull), pose axes, and a derived
                                fleet CG — see the README for the schema

PLACE (reconstruct CAD geometry from a studied configuration)
  michell place ama.igs@dy=1.7,dz=0.05 ama.igs@dy=-1.7,dz=0.05 -o boat.igs
  Writes the input geometry, posed, as a new IGES file (untrimmed 128
  surfaces, metres) for import back into CAD — e.g. amas at the dx/dy/dz a
  sweep found good. Inputs: IGES files (surfaces pass through exactly, with
  each spec's pose applied to every hull in that file). Suffix keys, all
  optional:
      dx / x    longitudinal shift [m]        dy  transverse shift [m]
      dz        immersion [m], + is deeper
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

SWEEPS (flag form; hulls modelled in position)
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

// `--sections` is accepted and ignored: every hull is sectional now.
const SWITCHES: &[&str] = &["--json", "--knots", "--csv", "--sections", "--dynamic"];

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

    /// The viscous knobs: form factor and roughness allowance. They are read
    /// together because they compose as `(1+k)·C_F + ΔC_F` and are easy to
    /// conflate — a roughness folded into `k` is invisible in the report.
    fn viscous_options(&self) -> Result<michell::ViscousOptions, String> {
        Ok(michell::ViscousOptions {
            form_factor: self.f64_flag("form-factor")?.unwrap_or(0.0),
            roughness: match self.flag("roughness") {
                Some(spec) => formats::parse_roughness(spec)?,
                None => michell::Roughness::None,
            },
        })
    }

    /// The import settings: `--waterline`, `--centerplane`, and the
    /// sectioning resolution `--stations N` / `--rays M`.
    fn load_settings(&self) -> Result<LoadSettings, String> {
        for gone in [
            "samples",
            "fit-degree",
            "fit-control",
            "fit-deriv-weight",
            "fit-fairing",
            "dump-grid",
        ] {
            if self.flag(gone).is_some() {
                return Err(format!(
                    "--{gone} was removed with lofting; hulls are cut into sections \
                     (resolution: --stations N --rays M)"
                ));
            }
        }
        let d = LoadSettings::default();
        let count = |name: &str, default: usize, min: usize| -> Result<usize, String> {
            match self.flag(name) {
                None => Ok(default),
                Some(v) => match v.trim().parse::<usize>() {
                    Ok(n) if n >= min => Ok(n),
                    _ => Err(format!(
                        "--{name} {v:?}: expected a count of at least {min}"
                    )),
                },
            }
        };
        Ok(LoadSettings {
            waterline_z: self.f64_flag("waterline")?.unwrap_or(0.0),
            centerplane: self.f64_flag("centerplane")?,
            stations: count("stations", d.stations, 8)?,
            rays: count("rays", d.rays, 5)?,
            units: self
                .flag("units")
                .map(|u| formats::parse_units(u))
                .transpose()?,
        })
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

/// One-line description of a roughness allowance, empty when smooth, for the
/// header where the form factor is reported. Kept visible: the allowance can
/// exceed the form factor on a small hull, and it is speed-dependent when it
/// comes from a sand-grain height, so it must not look like a constant.
fn describe_roughness(r: &michell::Roughness) -> String {
    match r {
        michell::Roughness::None => String::new(),
        michell::Roughness::DeltaCf(c) => format!(", roughness dCf {c:.3e}"),
        michell::Roughness::SandGrain(k) => {
            format!(", roughness k_s {:.0} um (dCf varies with speed)", k * 1e6)
        }
    }
}

/// Where a sand-grain roughness spec lands on the smooth / transitional /
/// fully-rough scale, over the speeds actually run.
///
/// `k_s⁺ = k_s·u_τ/ν` is the only thing that says whether a given finish
/// matters at a given speed, and it is what bounds how far to trust the
/// allowance: the estimator bridges the transitional band by the smooth /
/// fully-rough crossover rather than fitting it, so it reads high in there.
/// Printing the range keeps that visible instead of leaving it in the docs.
fn roughness_regime_note<'a>(
    mut viscous: impl Iterator<Item = &'a michell::ViscousResistance>,
) -> Option<String> {
    let first = viscous.find_map(|v| v.roughness_reynolds)?;
    let (mut lo, mut hi) = (first, first);
    for v in viscous.filter_map(|v| v.roughness_reynolds) {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    let regime = if hi < 5.0 {
        "hydraulically smooth over this speed range — the allowance is ~0 and \
         the finish is not costing you anything"
    } else if lo > 70.0 {
        "fully rough: the allowance is on its firm asymptote"
    } else {
        "transitional (5 < k_s+ < 70), where the estimator bridges the \
         crossover rather than fitting it and so reads high; treat dC_F as an \
         upper bound, or pass --roughness cf=... if you have a better number"
    };
    Some(format!(
        "                roughness k_s+ {lo:.1}..{hi:.1}, {regime}"
    ))
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

/// `michell info` in-process: geometry and diagnostics for a fleet of hull
/// specs (arguments without the leading `info`). Returns the text the CLI would
/// have printed to stdout (a JSON array with `--json`, otherwise a report).
pub fn info(args: &[String]) -> Result<String, String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err("usage: michell info <hull>... [options]".into());
    }
    let settings = p.load_settings()?;
    let fleet = fleet::load(&p.positional, &settings)?;
    let members = &fleet.members;
    let mut out = String::new();
    if p.switch("--json") {
        out.push('[');
        for (i, m) in members.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"path\":{:?},\"placement\":{{\"x\":{},\"y\":{}}},\"length\":{},\
                 \"beam\":{},\"draft\":{},\"wetted_surface\":{},\"displaced_volume\":{}}}",
                m.path,
                m.placement.x,
                m.placement.y,
                m.hull.length(),
                max_beam(&m.hull),
                m.hull.draft(),
                m.hull.wetted_surface(),
                m.hull.displaced_volume()
            ));
        }
        out.push_str("]\n");
        return Ok(out);
    }
    // Several hulls (a multihull file, or several specs): lead with a fleet
    // overview so the shape of the fleet is visible before the per-hull
    // detail blocks.
    if members.len() > 1 {
        let _ = writeln!(out, "fleet: {} hulls", members.len());
        for (i, m) in members.iter().enumerate() {
            let _ = writeln!(
                out,
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
        let _ = writeln!(
            out,
            "  total: wetted surface {total_s:.4} m^2, displaced vol {total_v:.4} m^3"
        );
    }
    for (i, m) in members.iter().enumerate() {
        if i > 0 || members.len() > 1 {
            out.push('\n');
        }
        let name = if members.len() > 1 {
            format!("hull [{}]: {}", i + 1, m.path)
        } else {
            format!("hull: {}", m.path)
        };
        if m.placement.x.abs() < 1e-9 && m.placement.y.abs() < 1e-9 {
            let _ = writeln!(out, "{name}");
        } else {
            let _ = writeln!(
                out,
                "{name} (placed at dx {:+.3} m, y {:.4} m)",
                m.placement.x, m.placement.y
            );
        }
        for line in describe(&fleet, m) {
            let _ = writeln!(out, "{line}");
        }
        let _ = writeln!(out, "length          {:>10.4} m", m.hull.length());
        let _ = writeln!(out, "beam            {:>10.4} m", max_beam(&m.hull));
        let _ = writeln!(out, "draft           {:>10.4} m", m.hull.draft());
        let _ = writeln!(out, "wetted surface  {:>10.4} m^2", m.hull.wetted_surface());
        let _ = writeln!(
            out,
            "displaced vol   {:>10.4} m^3",
            m.hull.displaced_volume()
        );
        let _ = writeln!(out, "LCB             {:>10.4} m", m.hull.lcb_x());
        let _ = writeln!(
            out,
            "waterplane      {:>10.4} m^2",
            m.hull.waterplane_area()
        );
        let _ = writeln!(
            out,
            "sections        {} stations, {} rays each",
            m.report.stations, settings.rays
        );
    }
    Ok(out)
}

/// Parse `--transom`: `off`, `ballistic[=COEFF]`, or `hollow=METRES`.
///
/// Default (flag absent) is the ballistic hollow — inert on a hull that closes
/// aft, so only transom-sterned geometry is affected.
pub(crate) fn parse_transom(spec: Option<&str>) -> Result<michell::TransomClosure, String> {
    use michell::TransomClosure;
    let Some(spec) = spec else {
        return Ok(TransomClosure::default());
    };
    let (key, val) = match spec.split_once('=') {
        Some((k, v)) => (k, Some(v)),
        None => (spec, None),
    };
    let num = |v: Option<&str>, what: &str| -> Result<f64, String> {
        v.ok_or_else(|| format!("--transom {key}: expected {key}=<{what}>"))?
            .parse::<f64>()
            .map_err(|_| format!("--transom: cannot parse number {:?}", v.unwrap_or("")))
            .and_then(|x| {
                if x.is_finite() && x >= 0.0 {
                    Ok(x)
                } else {
                    Err(format!("--transom: {what} must be finite and >= 0"))
                }
            })
    };
    match key {
        "off" | "none" => Ok(TransomClosure::None),
        "ballistic" => Ok(match val {
            None => TransomClosure::default(),
            Some(_) => TransomClosure::Ballistic {
                coeff: num(val, "coefficient")?,
            },
        }),
        "hollow" | "length" => Ok(TransomClosure::Fixed {
            length: num(val, "metres")?,
        }),
        other => Err(format!(
            "--transom {other:?}: expected off | ballistic[=COEFF] | hollow=METRES"
        )),
    }
}

fn cmd_squat(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err(
            "usage: michell squat <hull>[@x=DX,y=Y]... --speeds A[:B:STEP] [options]".into(),
        );
    }
    let fleet = fleet::load(&p.positional, &p.load_settings()?)?;
    let cut_members = fleet.hulls();
    // Per hull: (length, waterplane area, its first and second moments about
    // x = 0, volume, placement).
    let hydro: Vec<(f64, f64, f64, f64, f64, Placement)> = cut_members
        .iter()
        .map(|(h, pl)| {
            (
                h.length(),
                h.waterplane_area(),
                h.waterplane_moment(),
                h.waterplane_second_moment(),
                h.displaced_volume(),
                *pl,
            )
        })
        .collect();
    let l_ref = hydro.iter().map(|h| h.0).fold(0.0f64, f64::max);
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
    let mut opts = michell::squat::SquatOptions::default();
    if let Some(t) = p.f64_flag("rel-tol")? {
        opts.rel_tol = t;
    }
    opts.wave.transom = parse_transom(p.flag("transom").map(|s| s.as_str()))?;

    // Fleet waterplane in fleet coordinates: area, first and second moments,
    // and the LCF, which is the default pivot and the axis I_L is taken about.
    let (mut aw, mut mw, mut iw, mut vol) = (0.0, 0.0, 0.0, 0.0);
    for &(_, a, m1, m2, v, pl) in &hydro {
        aw += a;
        mw += m1 + pl.x * a;
        iw += m2 + 2.0 * pl.x * m1 + pl.x * pl.x * a;
        vol += v;
    }
    if aw <= 0.0 {
        return Err("fleet has no waterplane".into());
    }
    let lcf = mw / aw;
    let pivot = p.f64_flag("pivot")?.unwrap_or(lcf);
    // Trimming inertia about the LCF (the parallel-axis correction of I_w).
    let i_l = iw - mw * mw / aw;

    for (h, _) in &cut_members {
        if let Some(t) = h.transom() {
            eprintln!(
                "note: transom immersed ({:.1}% of max section); squat is evaluated on the \
                 closed composite body ({:?})",
                100.0 * t.area / h.max_section_area(),
                opts.wave.transom
            );
        }
    }
    println!(
        "fleet: {} hull(s), L(ref) {:.3} m, Aw {:.3} m^2, LCF {:.3} m, I_L {:.3} m^4, vol {:.3} m^3; pivot {:.3} m",
        hydro.len(),
        l_ref,
        aw,
        lcf,
        i_l,
        vol,
        pivot
    );
    println!(
        "{:>8} {:>7} {:>11} {:>8} {:>12} {:>12} {:>12} {:>9} {:>8}",
        if knots { "U[kn]" } else { "U[m/s]" },
        "Fn",
        "Fz[N]",
        "lift%",
        "M[N m]",
        "s_eq[mm]",
        "trim_eq[deg]",
        "M_wave%",
        "err"
    );
    for &u in &speeds {
        let cond = p.conditions(u)?;
        let d = michell::sectional::multihull_dynamic_force(&cut_members, &cond, pivot, &opts)
            .map_err(|e| format!("{e}"))?;
        let rho = cond.fluid.density;
        let s_eq = -d.force_up / (rho * g * aw);
        let trim_eq = (d.moment_bow_up / (rho * g * i_l)).to_degrees();
        let wave_share = if d.moment_bow_up.abs() > 0.0 {
            100.0 * d.moment_wave / d.moment_bow_up
        } else {
            0.0
        };
        println!(
            "{:8.3} {:7.3} {:11.2} {:8.2} {:12.2} {:12.2} {:12.4} {:9.1} {:8.1e}",
            if knots { u / KNOT } else { u },
            u / (g * l_ref).sqrt(),
            d.force_up,
            100.0 * d.lift_fraction,
            d.moment_bow_up,
            1000.0 * s_eq,
            trim_eq,
            wave_share,
            d.est_rel_error
        );
    }
    Ok(())
}

fn cmd_resistance(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.flag("heel").is_some() {
        return Err("--heel was removed: michell no longer models heel".into());
    }
    if p.positional.is_empty() {
        return Err(
            "usage: michell resistance <hull>[@x=DX,y=Y]... --speeds A[:B:STEP] [options]".into(),
        );
    }
    let fleet = fleet::load(&p.positional, &p.load_settings()?)?;
    let summary = &fleet.members;
    let cut_members = fleet.hulls();
    let multi = summary.len() > 1;
    // Reference length for Froude number: the longest hull.
    let l_ref = fleet.l_ref();
    let total_s: f64 = summary.iter().map(|m| m.hull.wetted_surface()).sum();
    let total_v: f64 = summary.iter().map(|m| m.hull.displaced_volume()).sum();

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

    let viscous_opts = p.viscous_options()?;
    let form_factor = viscous_opts.form_factor;
    let mut wave_opts = WaveOptions::default();
    if let Some(t) = p.f64_flag("rel-tol")? {
        wave_opts.rel_tol = t;
    }
    wave_opts.transom = parse_transom(p.flag("transom").map(|s| s.as_str()))?;

    let mut rows = Vec::new();
    for &u in &speeds {
        let cond = p.conditions(u)?;
        let r = michell::sectional::multihull_resistance(
            &cut_members,
            &cond,
            &wave_opts,
            &viscous_opts,
        )
        .map_err(|e| format!("at U = {u} m/s: {e}"))?;
        rows.push((u, cond, r));
    }

    if p.switch("--json") {
        let fluid = rows
            .first()
            .map(|(_, c, _)| c.fluid)
            .unwrap_or(Fluid::SEAWATER_15C);
        let mut out = String::from("{\"hulls\":[");
        for (i, m) in summary.iter().enumerate() {
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
             \"roughness\":{},",
            fluid.density,
            fluid.kinematic_viscosity,
            form_factor,
            match viscous_opts.roughness {
                michell::Roughness::None => "null".to_string(),
                michell::Roughness::DeltaCf(c) => format!("{{\"delta_cf\":{c}}}"),
                michell::Roughness::SandGrain(k) => format!("{{\"sand_grain_m\":{k}}}"),
            },
        ));
        out.push_str("\"points\":[");
        for (i, (u, cond, r)) in rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"speed\":{u},\"froude\":{},\"rw\":{},\"rv\":{},\"total\":{},\
                 \"effective_power\":{},\"interference\":{},\"cw\":{},\"cv\":{},\"ct\":{},\
                 \"roughness_cf\":{},\"roughness_reynolds\":{},\
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
                // Per-member roughness is the same allowance applied to each,
                // so the fleet's value is any member's; the k_s+ regime is
                // per-member (it scales with each hull's own C_F).
                r.viscous.first().map_or(0.0, |v| v.roughness_cf),
                r.viscous
                    .first()
                    .and_then(|v| v.roughness_reynolds)
                    .map_or("null".to_string(), |k| format!("{k}")),
                r.wave.est_rel_error
            ));
        }
        out.push_str("]}");
        println!("{out}");
        return Ok(());
    }

    let mut seen: Vec<&str> = Vec::new();
    for m in summary {
        if seen.contains(&m.path.as_str()) {
            continue;
        }
        seen.push(&m.path);
        println!("hull: {}", m.path);
        for line in describe(&fleet, m) {
            println!("{line}");
        }
    }
    if multi {
        let placements: Vec<String> = summary
            .iter()
            .map(|m| format!("{}@x={},y={}", m.path, m.placement.x, m.placement.y))
            .collect();
        println!("fleet: {}", placements.join("  "));
    }
    println!(
        "L(ref) = {l_ref:.3} m, S = {total_s:.3} m^2, vol = {total_v:.3} m^3, \
         form factor {form_factor}{}",
        describe_roughness(&viscous_opts.roughness),
    );
    if let Some(note) = roughness_regime_note(rows.iter().flat_map(|(_, _, r)| r.viscous.iter())) {
        println!("{note}");
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

fn cmd_report(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    let Some(manifest_path) = p.positional.first() else {
        return Err(
            "usage: michell report <study.json> [-o report.pdf] [--cache solved.json]".into(),
        );
    };
    if !manifest_path.ends_with(".json") {
        return Err("michell report takes a JSON manifest, same schema as `michell sweep`".into());
    }
    if p.positional.len() > 1 {
        return Err("report takes a single manifest; put every axis in the JSON file".into());
    }
    let out_path = p
        .flag("output")
        .cloned()
        .unwrap_or_else(|| "report.pdf".to_string());
    report::run(
        manifest_path,
        &out_path,
        p.flag("cache").map(String::as_str),
    )
}

fn cmd_sweep(args: &[String]) -> Result<(), String> {
    use michell_geometry::float::{solve_equilibrium_sectional, LoadCase};
    use michell_geometry::iges::{HullPose, Platform};
    use michell_geometry::source::SourceHull;

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
    let base_wl = settings.waterline_z;
    let opts = settings.sectional(base_wl);

    // Load the source geometry and learn each hull's base transverse position.
    struct File {
        stem: String,
        src: fleet::SourceFile,
        base_y: Vec<f64>,
        n: usize,
    }
    let mut files: Vec<File> = Vec::new();
    let mut l_ref = 0.0f64;
    for path in &p.positional {
        let src = fleet::open_source(path, &settings)?;
        let n = src.source.len();
        let mut base_y = vec![0.0; n];
        for (hi, y) in base_y.iter_mut().enumerate() {
            let cut = src
                .source
                .situate_sectional(
                    hi,
                    src.waterline_z,
                    &HullPose::default(),
                    &Platform::default(),
                    &settings.sectional(src.waterline_z),
                )
                .map_err(|e| format!("{path}: hull {}: {e}", hi + 1))?;
            if let Some(h) = cut {
                *y = h.placement.y;
                l_ref = l_ref.max(h.hull.length());
            }
        }
        let stem = std::path::Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(path)
            .to_string();
        eprintln!(
            "loaded {path}: {n} hull(s) at y = {:?}",
            base_y
                .iter()
                .map(|y| (y * 1e4).round() / 1e4)
                .collect::<Vec<_>>()
        );
        files.push(File {
            stem,
            n,
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
    let viscous_opts = p.viscous_options()?;
    let mut wave_opts = WaveOptions::default();
    if let Some(t) = p.f64_flag("rel-tol")? {
        wave_opts.rel_tol = t;
    }
    wave_opts.transom = parse_transom(p.flag("transom").map(|s| s.as_str()))?;
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
        let mut members: Vec<(michell_geometry::SectionalHull, Placement)> = Vec::new();
        let (sinkage, trim_deg, volume, lcb, dry) = if let Some(mass) = weight {
            let f = &files[0];
            let hulls: Vec<SourceHull> = poses[0]
                .iter()
                .enumerate()
                .map(|(index, pose)| SourceHull {
                    source: f.src.source.as_ref(),
                    index,
                    waterline_z: f.src.waterline_z,
                    pose: *pose,
                })
                .collect();
            let eq = solve_equilibrium_sectional(&hulls, &LoadCase { mass, lcg }, density, &opts)
                .map_err(|e| format!("point {}: {e}", point + 1))?;
            let out = (
                eq.sinkage,
                eq.trim.to_degrees(),
                eq.volume,
                eq.lcb,
                eq.fleet.dry,
            );
            members = eq.fleet.members;
            out
        } else {
            let mut volume = 0.0;
            let mut moment = 0.0;
            let mut dry = 0usize;
            let o = settings.sectional(waterline);
            for (f, fp) in files.iter().zip(&poses) {
                for (hi, pose) in fp.iter().enumerate() {
                    match f
                        .src
                        .source
                        .situate_sectional(hi, waterline, pose, &Platform::default(), &o)
                        .map_err(|e| format!("point {}: {e}", point + 1))?
                    {
                        Some(m) => {
                            volume += m.hull.displaced_volume();
                            moment += m.hull.lcb_x() * m.hull.displaced_volume();
                            members.push((m.hull, m.placement));
                        }
                        None => dry += 1,
                    }
                }
            }
            let lcb = if volume > 0.0 { moment / volume } else { 0.0 };
            (0.0, 0.0, volume, lcb, dry)
        };
        let members: Vec<(&michell_geometry::SectionalHull, Placement)> =
            members.iter().map(|(h, p)| (h, *p)).collect();

        let mut rows = Vec::with_capacity(speeds.len());
        for &u in &speeds {
            let cond = p.conditions(u)?;
            let (rw, rv, rt, pe, iff, cw, ct) = if members.is_empty() {
                (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
            } else {
                let r = michell::sectional::multihull_resistance(
                    &members,
                    &cond,
                    &wave_opts,
                    &viscous_opts,
                )
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
    let fleet = fleet::load(&p.positional, &p.load_settings()?)?;
    let members = fleet.hulls();
    let l_ref = fleet.l_ref();
    let transom = parse_transom(p.flag("transom").map(|s| s.as_str()))?;
    let u = single_speed(&p, l_ref)?;
    let cond = p.conditions(u)?;
    let n = match p.flag("points") {
        None => 721,
        Some(v) => v
            .parse::<usize>()
            .map_err(|_| format!("--points: cannot parse {v:?}"))?
            .max(9),
    };

    let mut spec = michell::FreeWaveSpectrum::new_sectional(&members, &cond, transom)
        .map_err(|e| format!("{e}"))?;

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
    let wave_opts = WaveOptions {
        transom,
        ..WaveOptions::default()
    };
    let rw = michell::sectional::multihull_wave_resistance(&members, &cond, &wave_opts)
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

/// `michell field`: a viewer-neutral 3-D scene of IGES hulls cut into
/// sections at one speed (see [`scene`]).
fn cmd_field(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err(
            "usage: michell field <hull>[@x=DX,y=Y]... --speed U | --froude F -o scene.json \
             [--waterline Z --transom SPEC --region X0,X1,Y0,Y1 --size NX,NY --stations N --rays M]"
                .into(),
        );
    }
    let out_path = p
        .flag("o")
        .or_else(|| p.flag("output"))
        .ok_or("michell field needs -o scene.json")?
        .clone();
    let fleet = fleet::load(&p.positional, &p.load_settings()?)?;
    let l_ref = fleet.l_ref();
    let u = single_speed(&p, l_ref)?;
    let cond = p.conditions(u)?;
    let closure = parse_transom(p.flag("transom").map(|s| s.as_str()))?;

    let (mut x_lo, mut x_hi, mut y_abs) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
    for m in &fleet.members {
        let (a, b) = m.report.x_range;
        x_lo = x_lo.min(a + m.placement.x);
        x_hi = x_hi.max(b + m.placement.x);
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
            let nx = 400usize;
            let ny = ((nx as f64) * (y1 - y0) / (x1 - x0)).round() as usize;
            (nx, ny.clamp(32, 600))
        }
    };
    if nx < 2 || ny < 2 {
        return Err("--size: need at least 2x2 grid points".into());
    }
    let hulls: Vec<scene::SceneHull> = fleet
        .members
        .iter()
        .map(|m| {
            let file = &fleet.files[m.file];
            scene::SceneHull {
                name: if file.source.len() > 1 {
                    format!("{}#{}", m.path, m.index + 1)
                } else {
                    m.path.clone()
                },
                hull: &m.hull,
                placement: m.placement,
                shift: m.shift,
                source: file.source.as_ref(),
                index: m.index,
                waterline_z: file.waterline_z,
            }
        })
        .collect();
    let surface = scene::Surface {
        x0,
        x1,
        y0,
        y1,
        nx,
        ny,
    };
    let t = std::time::Instant::now();
    let text = scene::build(&hulls, &cond, closure, &surface, l_ref)?;
    std::fs::write(&out_path, text).map_err(|e| format!("cannot write {out_path}: {e}"))?;
    eprintln!(
        "wrote {out_path}: {} hull(s) at U = {:.3} m/s (Fn {:.3}), free surface {nx}x{ny} over \
         x {x0:.2}..{x1:.2}, y {y0:.2}..{y1:.2} m ({:.1} s)",
        hulls.len(),
        u,
        u / (cond.gravity * l_ref).sqrt(),
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

fn cmd_wake(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err("usage: michell wake <hull>... --speed U [-o wake.png] [options]".into());
    }
    let fleet = fleet::load(&p.positional, &p.load_settings()?)?;
    let loaded = &fleet.members;
    let members = fleet.hulls();
    let l_ref = fleet.l_ref();
    let transom = parse_transom(p.flag("transom").map(|s| s.as_str()))?;
    let u = single_speed(&p, l_ref)?;
    let cond = p.conditions(u)?;

    // Fleet extents in fleet coordinates.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_abs = 0.0f64;
    for m in loaded {
        let (h0, h1) = m.hull.x_range();
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

    let mut spec = michell::FreeWaveSpectrum::new_sectional(&members, &cond, transom)
        .map_err(|e| format!("{e}"))?;
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
    for m in loaded {
        let (h0, h1) = m.hull.x_range();
        for ix in 0..nx {
            let x = grid.x(ix) - m.placement.x;
            if x < h0 || x > h1 {
                continue;
            }
            let half_beam = m.hull.waterline_half_beam(x);
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
        return Err("usage: michell render <hull>... --speed U [-o render.png] [options]".into());
    }
    let fleet = fleet::load(&p.positional, &p.load_settings()?)?;
    let loaded = &fleet.members;
    let members = fleet.hulls();
    let l_ref = fleet.l_ref();
    let transom = parse_transom(p.flag("transom").map(|s| s.as_str()))?;
    let u = single_speed(&p, l_ref)?;
    let cond = p.conditions(u)?;

    // Fleet extents in fleet coordinates.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_abs = 0.0f64;
    for m in loaded {
        let (h0, h1) = m.hull.x_range();
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

    let mut spec = michell::FreeWaveSpectrum::new_sectional(&members, &cond, transom)
        .map_err(|e| format!("{e}"))?;
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
    for m in loaded {
        let file = &fleet.files[m.file];
        let (verts, tris) = file
            .source
            .posed_tessellation(
                m.index,
                file.waterline_z,
                &michell_geometry::iges::HullPose::default(),
                &michell_geometry::iges::Platform::default(),
            )
            .map_err(|e| format!("{}: {e}", m.path))?;
        render::add_mesh(&mut scene, &verts, &tris, m.shift.x, m.shift.y);
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
            let dist = if parts.len() == 3 {
                Some(nums[2])
            } else {
                None
            };
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

/// One `place` input: a path plus the pose to apply to every hull in it.
struct PlaceSpec {
    path: String,
    /// Absolute centerplane (`y=`): rejected, IGES carries none.
    y_abs: Option<f64>,
    pose: michell_geometry::iges::HullPose,
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
        pose: michell_geometry::iges::HullPose::default(),
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
    use michell_geometry::iges::{self, Platform};
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
    let mut cache: HashMap<String, iges::SourceFleet> = HashMap::new();
    let mut surfaces: Vec<michell_geometry::iges::NurbsSurface3> = Vec::new();
    for raw in &p.positional {
        let spec = parse_place_spec(raw)?;
        if spec.y_abs.is_some() {
            return Err(format!(
                "{raw:?}: absolute y placement needs a recorded centerplane, \
                 which IGES inputs don't carry; use dy=SHIFT"
            ));
        }
        if !cache.contains_key(&spec.path) {
            let path = &spec.path;
            let text =
                std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
            let fleet =
                iges::source_fleet(&text, waterline_z).map_err(|e| format!("{path}: {e}"))?;
            cache.insert(spec.path.clone(), fleet);
        }
        let before = surfaces.len();
        let fleet = &cache[&spec.path];
        for i in 0..fleet.len() {
            surfaces.extend(
                fleet
                    .posed_surfaces(i, waterline_z, &spec.pose, &platform)
                    .map_err(|e| format!("{}: {e}", spec.path))?,
            );
        }
        let hulls = fleet.len();
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
    let surfaces = michell_geometry::iges::wigley_surfaces(l, b, t).map_err(|e| format!("{e}"))?;
    let text = michell_geometry::iges::write(&surfaces, "wigley").map_err(|e| format!("{e}"))?;
    match p.flag("output") {
        Some(path) => {
            std::fs::write(path, text).map_err(|e| format!("cannot write {path}: {e}"))?;
            println!("wrote {path} (Wigley L={l} B={b} T={t}, IGES, design waterline at z = 0)");
        }
        None => print!("{text}"),
    }
    Ok(())
}

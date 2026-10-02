//! `boatmath` — hulls, cases, studies and their results as JSON records,
//! one per line, made and read by small commands in a pipe:
//!
//! ```sh
//! boatmath hull guillemot.igs --waterline 0.12 \
//!   | boatmath case --mass 150,200 \
//!   | boatmath study --froude 0.2:0.6:0.05 \
//!   | boatmath run -j 8 \
//!   | boatmath table study.case.params.mass study.params.froude forces.rt
//! ```
//!
//! Every record is saved in the store (`$BOATMATH_HOME`, default
//! `~/.boatmath`) and refers to its parents by id. See docs/boatmath-cli.md.

mod list;
mod path;
mod plot;
mod props;
mod records;
mod run;
mod store;
mod units;
mod views;

use boatmath::params::{CaseParams, Closure, Sea, StudyParams, Waves};
use boatmath::{LoftRequest, MassBy};
use clap::{Parser, Subcommand};
use serde_json::Value;
use std::io::{Read, Write};
use std::path::PathBuf;
use store::Store;

#[derive(Parser)]
#[command(
    name = "boatmath",
    version,
    about = "Hulls, cases and studies as JSON records in a pipe"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Hulls from files (IGES, STL): a hull record per file and setting, its
    /// geometry inline as JSON (B-spline patches, or triangles from an STL),
    /// in metres with the design waterline at z = 0.
    Hull {
        files: Vec<PathBuf>,
        #[arg(long)]
        name: Option<String>,
        /// Design waterline height in the file's frame [m] (IGES); a LIST.
        #[arg(long, allow_hyphen_values = true)]
        waterline: Option<String>,
        /// Centreplane override [m] (IGES).
        #[arg(long, allow_hyphen_values = true)]
        centerplane: Option<f64>,
        /// Stations along the hull (IGES).
        #[arg(long)]
        stations: Option<usize>,
        /// Rays across each section (IGES).
        #[arg(long)]
        rays: Option<usize>,
        /// Units of an STL: mm, cm, m, in, ft, or metres per unit.
        #[arg(long)]
        units: Option<String>,
    },
    /// Hulls made by scaling hulls on stdin, about the design waterline.
    Scale {
        /// Scale the whole hull by this factor; a LIST.
        #[arg(long)]
        by: Option<String>,
        /// Scale beam and draft only by this factor; a LIST.
        #[arg(long)]
        beam: Option<String>,
        #[arg(long)]
        name: Option<String>,
    },
    /// Load hulls on stdin into cases, with their statics: a case per hull
    /// and combination of the LISTs.
    Case {
        /// Make a catamaran with this centre span [m]; a LIST.
        #[arg(long)]
        span: Option<String>,
        /// Load [kg]; default the design displacement; a LIST.
        #[arg(long)]
        mass: Option<String>,
        /// How a given mass is carried: sinking, scale, or scale_yz.
        #[arg(long, default_value = "sinking")]
        mass_by: String,
        /// Longitudinal centre of gravity [m]; default the LCB; a LIST.
        #[arg(long, allow_hyphen_values = true)]
        lcg: Option<String>,
        /// Centre of gravity above the design waterline [m]; a LIST.
        #[arg(long, allow_hyphen_values = true)]
        vcg: Option<String>,
        /// Roll radius of gyration [m]; a LIST.
        #[arg(long)]
        kxx: Option<String>,
        /// Pitch radius of gyration, a fraction of the length; a LIST.
        #[arg(long)]
        kyy: Option<String>,
        /// Yaw radius of gyration [m]; a LIST.
        #[arg(long)]
        kzz: Option<String>,
        /// Roll damping, a fraction of critical; a LIST.
        #[arg(long)]
        roll_damping: Option<String>,
        #[arg(long)]
        name: Option<String>,
        /// Compute the statics again even if the store has them.
        #[arg(long)]
        force: bool,
    },
    /// Studies on cases on stdin: one per case and combination of the LISTs.
    Study {
        /// Length Froude numbers; a LIST.
        #[arg(long)]
        froude: String,
        /// Hold the attitude at rest instead of floating at speed.
        #[arg(long)]
        hold: bool,
        /// Transom closure: ballistic[:COEFF], fixed:LENGTH, or off.
        #[arg(long)]
        closure: Option<String>,
        /// Free-surface grid columns (calm water).
        #[arg(long)]
        grid: Option<usize>,
        /// In waves from these headings [deg, 180 head seas]; a LIST.
        #[arg(long, allow_hyphen_values = true)]
        waves: Option<String>,
        /// Wavelengths over the length (with --waves); a LIST.
        #[arg(long)]
        lambdas: Option<String>,
        /// An irregular sea (with --waves): bretschneider:hs=H,tp=T or
        /// jonswap:hs=H,tp=T[,gamma=G].
        #[arg(long)]
        sea: Option<String>,
    },
    /// Compute studies on stdin; each result is written as it is ready.
    Run {
        /// Studies at once; default the number of cores.
        #[arg(short, long)]
        jobs: Option<usize>,
        /// Compute again even if the store has a result.
        #[arg(long)]
        force: bool,
        /// No progress on stderr.
        #[arg(short, long)]
        quiet: bool,
    },
    /// Print stored records by id (or a prefix of one).
    Get { ids: Vec<String> },
    /// The wake of calm-water results on stdin, from above, as SVG.
    Wake {
        /// Write here; with several results, a pattern naming each by
        /// `{id8}` or `{froude}` (wake-{froude}.svg).
        #[arg(short)]
        o: Option<String>,
        /// Colour scale ±RANGE [m]; default the 99th percentile of |ζ|.
        #[arg(long)]
        range: Option<f64>,
        #[arg(long)]
        title: Option<String>,
    },
    /// The pressure on the hulls of calm-water results on stdin, from
    /// below, as SVG.
    Pressure {
        /// As for `wake`.
        #[arg(short)]
        o: Option<String>,
        /// Colour scale ±RANGE in C_p; default the 99th percentile of |C_p|.
        #[arg(long)]
        range: Option<f64>,
        #[arg(long)]
        title: Option<String>,
    },
    /// The hull of calm-water results on stdin in profile at its attitude,
    /// with the wave along its side, as SVG.
    Profile {
        /// As for `wake`.
        #[arg(short)]
        o: Option<String>,
        /// Draw the wave this many times its height.
        #[arg(long, default_value_t = 1.0)]
        wave_scale: f64,
        #[arg(long)]
        title: Option<String>,
    },
    /// The best B-series propeller for each calm-water result on stdin (its
    /// speed, and the thrust R_t / (1 − t)), or for --speed and --thrust.
    Prop {
        /// Total thrust (N, kN, kgf, lbf), instead of reading results.
        #[arg(long)]
        thrust: Option<String>,
        /// Ship speed (m/s, kn, mph, km/h), with --thrust.
        #[arg(long)]
        speed: Option<String>,
        /// Shafts the thrust is split across; default one per hull.
        #[arg(long)]
        shafts: Option<u32>,
        /// Wake fraction w: the propeller sees V (1 − w).
        #[arg(long, default_value_t = 0.0)]
        wake: f64,
        /// Thrust deduction t.
        #[arg(long, default_value_t = 0.0)]
        thrust_deduction: f64,
        /// Largest diameter (m, mm, in, …).
        #[arg(long)]
        d_max: String,
        #[arg(long, default_value = "40mm")]
        d_min: String,
        /// Blade counts to consider; a LIST.
        #[arg(long, default_value = "2,3,4,5,6,7")]
        blades: String,
        /// Shaft immersion, for cavitation.
        #[arg(long, default_value = "0.3")]
        depth: String,
        /// Keller's margin k: 0.2 single screw, 0.1 twin, 0 fast.
        #[arg(long, default_value_t = 0.2)]
        keller_k: f64,
        /// Allow cavitating designs.
        #[arg(long)]
        no_cavitation: bool,
        /// Leave the series at its Re = 2e6 fit.
        #[arg(long)]
        no_re_correct: bool,
        /// Keep each blade count to the series' own blade-area range.
        #[arg(long)]
        strict_ear: bool,
        /// A second point the same propeller must reach: the result of the
        /// same study at this Froude number (already run).
        #[arg(long)]
        top_froude: Option<f64>,
        /// Or the second point's speed and total thrust.
        #[arg(long)]
        top_speed: Option<String>,
        #[arg(long)]
        top_thrust: Option<String>,
        /// Keep the whole feasible curve, not just to 1.6× the least power.
        #[arg(long)]
        full_range: bool,
    },
    /// The motors that can drive each prop on stdin, each at its best
    /// reduction and point on the prop's curve, ranked: a drive record each.
    Match {
        /// Rank by electrical power, mass or price.
        #[arg(long, default_value = "power")]
        rank_by: String,
        /// Every winding, not just the best of each family.
        #[arg(long)]
        all_windings: bool,
        /// Direct drive only (a ratio of 1).
        #[arg(long)]
        direct_only: bool,
        /// A fixed reduction ratio, instead of each motor's best.
        #[arg(long)]
        ratio: Option<f64>,
        /// The largest reduction considered.
        #[arg(long, default_value_t = 25.0)]
        ratio_max: f64,
        /// Efficiency of one reduction stage.
        #[arg(long, default_value_t = propeller::motor::GEAR_ETA_DEFAULT)]
        gear_eta: f64,
        /// Leave out the controller's loss (motor efficiency only).
        #[arg(long)]
        no_controller_loss: bool,
        /// Hold the design point on peak ratings, not continuous.
        #[arg(long)]
        peak: bool,
        /// Don't require the prop's second point.
        #[arg(long)]
        no_top: bool,
        #[command(flatten)]
        filter: MotorFilter,
        /// Motors from this file instead of the vendored database: a
        /// motors.js, its JSON, or motor records (as `boatmath motors`
        /// writes, filtered with jq, say).
        #[arg(long)]
        motors: Option<PathBuf>,
    },
    /// The motor database, a record per motor.
    Motors {
        #[command(flatten)]
        filter: MotorFilter,
    },
    /// Print a stored blob (a result's `field`, say) by its id.
    Blob { id: String },
    /// Records on stdin as a table, a column per FIELD path (`forces.rt`,
    /// `study.case.params.span`, `|heave|`, ...); tab-separated.
    Table {
        #[arg(required = true)]
        fields: Vec<String>,
        /// A row per element of the array at this path, its fields read
        /// from the element first and then the record.
        #[arg(long)]
        explode: Option<String>,
        /// Comma-separated, quoted where needed.
        #[arg(long)]
        csv: bool,
        /// No header row.
        #[arg(long)]
        no_header: bool,
    },
    /// Records on stdin as an SVG line plot.
    Plot {
        /// The x field.
        #[arg(short)]
        x: String,
        /// A y field; repeat for several.
        #[arg(short, required = true)]
        y: Vec<String>,
        /// A series per value of this field (or fields, repeated).
        #[arg(long)]
        by: Vec<String>,
        /// A point per element of the array at this path (see `table`).
        #[arg(long)]
        explode: Option<String>,
        #[arg(long)]
        title: Option<String>,
        #[arg(long)]
        xlabel: Option<String>,
        #[arg(long)]
        ylabel: Option<String>,
        /// Write here instead of stdout.
        #[arg(short)]
        o: Option<PathBuf>,
    },
}

#[derive(clap::Args)]
struct MotorFilter {
    /// Only this maker's motors; repeat for several.
    #[arg(long)]
    vendor: Vec<String>,
    /// Only motors with a published efficiency map (tiers A, A-).
    #[arg(long)]
    mapped_only: bool,
    /// No stacked (multi-machine) motors.
    #[arg(long)]
    single_only: bool,
    /// Heaviest motor [kg].
    #[arg(long)]
    max_mass: Option<f64>,
    /// Largest outside diameter (m, mm, in, …).
    #[arg(long)]
    max_od: Option<String>,
    /// Dearest motor [USD].
    #[arg(long)]
    max_price: Option<f64>,
}

impl MotorFilter {
    fn filter(&self) -> Result<props::Filter, String> {
        Ok(props::Filter {
            vendors: self.vendor.clone(),
            mapped_only: self.mapped_only,
            single_only: self.single_only,
            max_mass: self.max_mass,
            max_od: self
                .max_od
                .as_deref()
                .map(units::length)
                .transpose()?
                .map(|m| m * 1e3),
            max_price: self.max_price,
        })
    }
}

/// A motor database from a file: propopt's motors.js or its JSON, or motor
/// records (their maps and controller efficiency the vendored database's).
fn motor_database(path: &PathBuf) -> Result<propeller::motor::Database, String> {
    use propeller::motor::Database;
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let start = text.find('{').unwrap_or(0);
    if let Ok(v) = serde_json::from_str::<Value>(text[start..].trim().trim_end_matches(';')) {
        if v.get("motors").is_some() {
            return Database::from_value(v);
        }
    }
    if text.trim_start().starts_with("const") || text.contains("MOTOR_DB") {
        return Database::from_js(&text);
    }
    let motors: Vec<Value> = serde_json::Deserializer::from_str(&text)
        .into_iter::<Value>()
        .collect::<Result<_, _>>()
        .map_err(|e| format!("{}: {e}", path.display()))?;
    let base = &Database::vendored().raw;
    Database::from_value(serde_json::json!({
        "motors": motors,
        "maps": base["maps"],
        "controller_eta": base["controller_eta"],
    }))
}

/// The views a picture command draws.
enum View {
    Wake(Option<f64>),
    Pressure(Option<f64>),
    Profile(f64),
}

/// Draw `view` of each calm-water result on stdin, to `o` (or stdout).
fn draw(
    store: &Store,
    out: &mut impl Write,
    view: View,
    o: Option<String>,
    title: Option<String>,
) -> Result<usize, String> {
    let records = read_records()?;
    if records.len() > 1 && !o.as_deref().is_some_and(|p| p.contains('{')) {
        return Err(format!(
            "{} results: give -o a pattern with {{id8}} or {{froude}} to name each picture",
            records.len()
        ));
    }
    let res = path::Resolver::new(store);
    let mut failed = 0;
    for r in &records {
        let id = records::expect(r, "result")?;
        let picture = (|| -> Result<String, String> {
            if r["kind"] != "calm" {
                return Err("not a calm-water result".into());
            }
            let blob = r["field"]["blob"].as_str().ok_or("no field")?;
            let v: Value =
                serde_json::from_slice(&store.blob(blob)?).map_err(|e| format!("field: {e}"))?;
            let field = views::Field::parse(&v)?;
            let name = path::cell(&res.get(r, "study.case.hull.name"));
            let froude = r["froude"].as_f64().unwrap_or(0.0);
            let f = &r["forces"];
            let num = |k: &str| f[k].as_f64().map_or("?".into(), |v| format!("{v:.1}"));
            let title = title
                .clone()
                .unwrap_or_else(|| format!("{name} · Fn {froude}"));
            let sub = format!(
                "R_t {} N · R_w {} N · {}",
                num("rt"),
                num("rw"),
                store::short(id)
            );
            Ok(match view {
                View::Wake(range) => views::wake(&field, &title, &sub, range),
                View::Pressure(range) => views::pressure(&field, &title, &sub, range),
                View::Profile(scale) => {
                    let sec = res
                        .get(r, "sections.attitude")
                        .ok_or("no sections at its attitude")?;
                    let attitude = (
                        sec["sinkage"].as_f64().unwrap_or(0.0),
                        sec["trim_rad"].as_f64().unwrap_or(0.0),
                    );
                    let case = store.need(
                        "case",
                        res.get(r, "study.case.id")
                            .and_then(|v| v.as_str().map(String::from))
                            .as_deref()
                            .unwrap_or(""),
                    )?;
                    let hull = store.need("hull", case["hull"].as_str().unwrap_or(""))?;
                    let src = records::hull_source(&hull)?;
                    let params = serde_json::from_value(case["params"].clone())
                        .map_err(|e| format!("case params: {e}"))?;
                    let meshes = boatmath::sections::meshes(
                        &src.file_name,
                        src.bytes,
                        &src.import,
                        &params,
                        attitude,
                    )?;
                    views::profile(&field, &meshes[0].0, &title, &sub, attitude, scale)?
                }
            })
        })();
        match picture {
            Ok(svg) => match &o {
                Some(p) => {
                    let froude = r["froude"].as_f64().unwrap_or(0.0);
                    let path = p
                        .replace("{id8}", &id[..8.min(id.len())])
                        .replace("{froude}", &format!("{froude}"));
                    std::fs::write(&path, svg).map_err(|e| format!("{path}: {e}"))?;
                }
                None => out.write_all(svg.as_bytes()).map_err(stdout_err)?,
            },
            Err(e) => {
                eprintln!("boatmath: result {}: {e}", store::short(id));
                failed += 1;
            }
        }
    }
    Ok(failed)
}

/// A field path's last segment, to name a series by: `span` for
/// `study.case.params.span`, `|heave|` for `|seakeeping.heave|`.
fn short_field(f: &str) -> String {
    match f.strip_prefix('|').and_then(|p| p.strip_suffix('|')) {
        Some(inner) => format!("|{}|", short_field(inner)),
        None => f.rsplit('.').next().unwrap_or(f).to_string(),
    }
}

/// A failed write to stdout. A reader that stops early (`| head`) closes
/// the pipe; like any Unix filter, stop quietly then.
fn stdout_err(e: std::io::Error) -> String {
    if e.kind() == std::io::ErrorKind::BrokenPipe {
        std::process::exit(0);
    }
    format!("stdout: {e}")
}

fn csv_cell(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Records on stdin: JSON values, one per line or simply one after another.
fn read_records() -> Result<Vec<Value>, String> {
    let mut s = String::new();
    std::io::stdin()
        .read_to_string(&mut s)
        .map_err(|e| format!("stdin: {e}"))?;
    serde_json::Deserializer::from_str(&s)
        .into_iter::<Value>()
        .map(|v| v.map_err(|e| format!("stdin: {e}")))
        .collect()
}

fn list(s: &Option<String>, what: &str) -> Result<Option<Vec<f64>>, String> {
    s.as_deref()
        .map(|s| list::parse(s).map_err(|e| format!("--{what}: {e}")))
        .transpose()
}

fn mass_by(s: &str) -> Result<MassBy, String> {
    serde_json::from_value(Value::String(s.into()))
        .map_err(|_| format!("--mass-by {s:?}: expected sinking, scale or scale_yz"))
}

fn closure(s: &str) -> Result<Closure, String> {
    let (kind, arg) = s.split_once(':').unwrap_or((s, ""));
    let num = |a: &str| {
        a.parse::<f64>()
            .map_err(|_| format!("--closure {s:?}: expected a number after the colon"))
    };
    Ok(match kind {
        "ballistic" if arg.is_empty() => Closure::Ballistic { coeff: None },
        "ballistic" => Closure::Ballistic {
            coeff: Some(num(arg)?),
        },
        "fixed" => Closure::Fixed { length: num(arg)? },
        "off" => Closure::Off,
        _ => {
            return Err(format!(
                "--closure {s:?}: expected ballistic[:COEFF], fixed:LENGTH or off"
            ))
        }
    })
}

fn sea(s: &str) -> Result<Sea, String> {
    let (kind, args) = s.split_once(':').unwrap_or((s, ""));
    let (mut hs, mut tp, mut gamma) = (None, None, None);
    for kv in args.split(',').filter(|kv| !kv.is_empty()) {
        let (k, v) = kv
            .split_once('=')
            .ok_or_else(|| format!("--sea: {kv:?}: expected key=value"))?;
        let v: f64 = v
            .parse()
            .map_err(|_| format!("--sea: {kv:?}: expected a number"))?;
        match k {
            "hs" => hs = Some(v),
            "tp" => tp = Some(v),
            "gamma" => gamma = Some(v),
            _ => return Err(format!("--sea: {k:?}: expected hs, tp or gamma")),
        }
    }
    let (hs, tp) = match (hs, tp) {
        (Some(h), Some(t)) => (h, t),
        _ => return Err(format!("--sea {s:?}: needs hs= and tp=")),
    };
    match kind {
        "bretschneider" | "ittc" if gamma.is_none() => Ok(Sea::Bretschneider { hs, tp }),
        "jonswap" => Ok(Sea::Jonswap { hs, tp, gamma }),
        _ => Err(format!(
            "--sea {s:?}: expected bretschneider:… or jonswap:…"
        )),
    }
}

fn main() {
    let cli = Cli::parse();
    match go(cli) {
        Ok(0) => {}
        Ok(_) => std::process::exit(1),
        Err(e) => {
            eprintln!("boatmath: {e}");
            std::process::exit(1);
        }
    }
}

/// Run a command; returns how many of its records failed.
fn go(cli: Cli) -> Result<usize, String> {
    let store = Store::open()?;
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    let mut emit = |r: &Value| -> Result<(), String> { writeln!(out, "{r}").map_err(stdout_err) };
    let mut failed = 0;
    match cli.command {
        Command::Hull {
            files,
            name,
            waterline,
            centerplane,
            stations,
            rays,
            units,
        } => {
            if files.is_empty() {
                return Err("hull: name a hull file".into());
            }
            let units = units.as_deref().map(boatmath::parse_units).transpose()?;
            for f in &files {
                for w in list::product(&[list(&waterline, "waterline")?]) {
                    let import = LoftRequest {
                        waterline: w[0],
                        centerplane,
                        stations,
                        rays,
                        units,
                        ..LoftRequest::default()
                    };
                    match records::hull_from_file(&store, f, name.as_deref(), import) {
                        Ok(r) => emit(&r)?,
                        Err(e) => {
                            eprintln!("boatmath: {e}");
                            failed += 1;
                        }
                    }
                }
            }
        }
        Command::Scale { by, beam, name } => {
            for h in read_records()? {
                for k in list::product(&[list(&by, "by")?, list(&beam, "beam")?]) {
                    emit(&records::scaled_hull(
                        &store,
                        &h,
                        k[0],
                        k[1],
                        name.as_deref(),
                    )?)?;
                }
            }
        }
        Command::Case {
            span,
            mass,
            mass_by: by,
            lcg,
            vcg,
            kxx,
            kyy,
            kzz,
            roll_damping,
            name,
            force,
        } => {
            let by = mass_by(&by)?;
            let axes = [
                list(&span, "span")?,
                list(&mass, "mass")?,
                list(&lcg, "lcg")?,
                list(&vcg, "vcg")?,
                list(&kxx, "kxx")?,
                list(&kyy, "kyy")?,
                list(&kzz, "kzz")?,
                list(&roll_damping, "roll-damping")?,
            ];
            for h in read_records()? {
                for v in list::product(&axes) {
                    let params = CaseParams {
                        span: v[0],
                        mass: v[1],
                        mass_by: by,
                        lcg: v[2],
                        vcg: v[3],
                        kxx: v[4],
                        kyy: v[5],
                        kzz: v[6],
                        roll_damping: v[7].unwrap_or(0.0),
                    };
                    let r = records::case(&store, &h, params, name.as_deref(), force)?;
                    if let Some(e) = r["error"].as_str() {
                        eprintln!(
                            "boatmath: case {}: {e}",
                            store::short(r["id"].as_str().unwrap_or(""))
                        );
                        failed += 1;
                    }
                    emit(&r)?;
                }
            }
        }
        Command::Study {
            froude,
            hold,
            closure: cl,
            grid,
            waves,
            lambdas,
            sea: sea_arg,
        } => {
            let froudes = list::parse(&froude).map_err(|e| format!("--froude: {e}"))?;
            let closure = cl.as_deref().map(closure).transpose()?.unwrap_or_default();
            let lambdas = list(&lambdas, "lambdas")?;
            let sea = sea_arg.as_deref().map(sea).transpose()?;
            if waves.is_none() && (lambdas.is_some() || sea.is_some()) {
                return Err("--lambdas and --sea need --waves".into());
            }
            let headings = list(&waves, "waves")?;
            for c in read_records()? {
                for v in list::product(&[headings.clone(), Some(froudes.clone())]) {
                    let params = StudyParams {
                        froude: v[1].expect("a froude number"),
                        dynamic: !hold,
                        closure,
                        grid: grid.unwrap_or_else(boatmath::params::default_grid),
                        waves: v[0].map(|heading| Waves {
                            heading,
                            lambdas: lambdas.clone(),
                            sea,
                        }),
                    };
                    emit(&records::study(&store, &c, params)?)?;
                }
            }
        }
        Command::Run { jobs, force, quiet } => {
            let jobs =
                jobs.unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
            // The workers write to stdout themselves: let go of its lock.
            drop(out);
            failed += run::run(
                &store,
                &read_records()?,
                &run::Options { jobs, force, quiet },
            )?;
        }
        Command::Wake { o, range, title } => {
            failed += draw(&store, &mut out, View::Wake(range), o, title)?;
        }
        Command::Pressure { o, range, title } => {
            failed += draw(&store, &mut out, View::Pressure(range), o, title)?;
        }
        Command::Profile {
            o,
            wave_scale,
            title,
        } => {
            failed += draw(&store, &mut out, View::Profile(wave_scale), o, title)?;
        }
        Command::Prop {
            thrust,
            speed,
            shafts,
            wake,
            thrust_deduction,
            d_max,
            d_min,
            blades,
            depth,
            keller_k,
            no_cavitation,
            no_re_correct,
            strict_ear,
            top_froude,
            top_speed,
            top_thrust,
            full_range,
        } => {
            let blades = list::parse(&blades)
                .map_err(|e| format!("--blades: {e}"))?
                .into_iter()
                .map(|z| {
                    let z = z.round();
                    if (2.0..=7.0).contains(&z) {
                        Ok(z as u32)
                    } else {
                        Err(format!("--blades: {z}: the series has 2 to 7 blades"))
                    }
                })
                .collect::<Result<Vec<u32>, String>>()?;
            if !(0.0..1.0).contains(&wake) || !(0.0..1.0).contains(&thrust_deduction) {
                return Err("--wake and --thrust-deduction are fractions, 0 to 1".into());
            }
            let explicit_top = match (&top_speed, &top_thrust) {
                (Some(s), Some(t)) => Some(propeller::bseries::TopInputs {
                    speed: units::speed(s)?,
                    thrust: units::force(t)?,
                }),
                (None, None) => None,
                _ => return Err("--top-speed and --top-thrust go together".into()),
            };
            let base = |speed: f64, thrust: f64, shafts: u32, top| propeller::bseries::Inputs {
                speed,
                thrust,
                shafts,
                wake,
                thrust_deduction,
                d_min: 0.0,
                d_max: 0.0,
                blades: blades.clone(),
                depth: 0.0,
                keller_k,
                cavitation: !no_cavitation,
                re_correct: !no_re_correct,
                strict_ear,
                top,
            };
            let (d_min, d_max, depth) = (
                units::length(&d_min)?,
                units::length(&d_max)?,
                units::length(&depth)?,
            );
            let cap = if full_range { f64::INFINITY } else { 1.6 };
            let mut cases: Vec<(propeller::bseries::Inputs, Option<String>)> = Vec::new();
            match (&thrust, &speed) {
                (Some(t), Some(s)) => {
                    if top_froude.is_some() {
                        return Err("--top-froude needs results on stdin; give --top-speed and --top-thrust".into());
                    }
                    cases.push((
                        base(
                            units::speed(s)?,
                            units::force(t)?,
                            shafts.unwrap_or(1),
                            explicit_top,
                        ),
                        None,
                    ));
                }
                (None, None) => {
                    for r in read_records()? {
                        let one = (|| -> Result<_, String> {
                            let id = records::expect(&r, "result")?.to_string();
                            let (v, t) = props::operating_point(&r, thrust_deduction)?;
                            let top = match top_froude {
                                Some(f) => {
                                    let o = props::result_at_froude(&store, &r, f)?;
                                    let (vt, tt) = props::operating_point(&o, thrust_deduction)?;
                                    Some(propeller::bseries::TopInputs {
                                        speed: vt,
                                        thrust: tt,
                                    })
                                }
                                None => explicit_top,
                            };
                            let n = shafts.unwrap_or_else(|| props::hull_count(&store, &r) as u32);
                            Ok((base(v, t, n, top), Some(id)))
                        })();
                        match one {
                            Ok(c) => cases.push(c),
                            Err(e) => {
                                eprintln!("boatmath: {e}");
                                failed += 1;
                            }
                        }
                    }
                }
                _ => return Err("--thrust and --speed go together".into()),
            }
            for (mut inputs, result) in cases {
                inputs.d_min = d_min;
                inputs.d_max = d_max;
                inputs.depth = depth;
                let r = props::prop(&store, &inputs, result.as_deref(), cap)?;
                let b = &r["best"];
                if b.is_null() {
                    eprintln!(
                        "boatmath: prop {}: no B-series propeller can do it (try a larger --d-max, or --no-cavitation)",
                        store::short(r["id"].as_str().unwrap_or(""))
                    );
                    failed += 1;
                } else {
                    eprintln!(
                        "boatmath: prop {}: {:.0} W at {:.0} rpm, {}-blade, D {:.0} mm, P/D {:.2}, EAR {:.2}, η₀ {:.3}",
                        store::short(r["id"].as_str().unwrap_or("")),
                        b["P_shaft"].as_f64().unwrap_or(0.0),
                        b["rpm"].as_f64().unwrap_or(0.0),
                        b["Z"].as_f64().unwrap_or(0.0),
                        1e3 * b["D"].as_f64().unwrap_or(0.0),
                        b["PD"].as_f64().unwrap_or(0.0),
                        b["EAR"].as_f64().unwrap_or(0.0),
                        b["eta0"].as_f64().unwrap_or(0.0),
                    );
                }
                emit(&r)?;
            }
        }
        Command::Match {
            rank_by,
            all_windings,
            direct_only,
            ratio,
            ratio_max,
            gear_eta,
            no_controller_loss,
            peak,
            no_top,
            filter,
            motors,
        } => {
            use propeller::motor::RankBy;
            let rank_by = match rank_by.as_str() {
                "power" => RankBy::Power,
                "mass" => RankBy::Mass,
                "price" => RankBy::Price,
                o => return Err(format!("--rank-by {o:?}: expected power, mass or price")),
            };
            let owned;
            let db = match &motors {
                Some(p) => {
                    owned = motor_database(p)?;
                    &owned
                }
                None => propeller::motor::Database::vendored(),
            };
            let o = props::MatchOptions {
                filter: filter.filter()?,
                rank_by,
                all_windings,
                ratio: if direct_only { Some(1.0) } else { ratio },
                ratio_max,
                gear_eta,
                controller_loss: !no_controller_loss,
                peak,
                check_top: !no_top,
            };
            for p in read_records()? {
                let (drives, unreachable) = props::drives(&store, db, &p, &o)?;
                eprintln!(
                    "boatmath: prop {}: {} motors can drive it, {} can't",
                    store::short(p["id"].as_str().unwrap_or("")),
                    drives.len(),
                    unreachable.len()
                );
                for d in &drives {
                    emit(d)?;
                }
            }
        }
        Command::Motors { filter } => {
            let f = filter.filter()?;
            for m in &propeller::motor::Database::vendored().motors {
                if f.admits(m) {
                    emit(&props::motor_record(m))?;
                }
            }
        }
        Command::Blob { id } => {
            let bytes = store.blob(&id)?;
            out.write_all(&bytes).map_err(stdout_err)?;
        }
        Command::Get { ids } => {
            for id in &ids {
                let found = store.find(id)?;
                if found.is_empty() {
                    eprintln!("boatmath: no record {id}");
                    failed += 1;
                }
                for r in found {
                    emit(&r)?;
                }
            }
        }
        Command::Table {
            fields,
            explode,
            csv,
            no_header,
        } => {
            let res = path::Resolver::new(&store);
            let line = |cells: Vec<String>| {
                if csv {
                    cells
                        .iter()
                        .map(|c| csv_cell(c))
                        .collect::<Vec<_>>()
                        .join(",")
                } else {
                    cells
                        .iter()
                        .map(|c| c.replace(['\t', '\n'], " "))
                        .collect::<Vec<_>>()
                        .join("\t")
                }
            };
            if !no_header {
                writeln!(out, "{}", line(fields.clone())).map_err(stdout_err)?;
            }
            for r in read_records()? {
                for row in res.rows(&r, explode.as_deref()) {
                    let cells = fields
                        .iter()
                        .map(|f| path::cell(&res.get_in(&row, &r, f)))
                        .collect();
                    writeln!(out, "{}", line(cells)).map_err(stdout_err)?;
                }
            }
        }
        Command::Plot {
            x,
            y,
            by,
            explode,
            title,
            xlabel,
            ylabel,
            o,
        } => {
            let res = path::Resolver::new(&store);
            let mut series: Vec<plot::Series> = Vec::new();
            for r in read_records()? {
                for row in res.rows(&r, explode.as_deref()) {
                    let Some(xv) = res.get_in(&row, &r, &x).and_then(|v| v.as_f64()) else {
                        continue;
                    };
                    let group: Vec<String> = by
                        .iter()
                        .map(|f| {
                            format!(
                                "{} {}",
                                short_field(f),
                                path::cell(&res.get_in(&row, &r, f))
                            )
                        })
                        .collect();
                    for yf in &y {
                        let Some(yv) = res.get_in(&row, &r, yf).and_then(|v| v.as_f64()) else {
                            continue;
                        };
                        let mut name = group.clone();
                        if y.len() > 1 {
                            name.insert(0, short_field(yf));
                        }
                        let name = name.join(", ");
                        match series.iter_mut().find(|s| s.name == name) {
                            Some(s) => s.points.push((xv, yv)),
                            None => series.push(plot::Series {
                                name,
                                points: vec![(xv, yv)],
                            }),
                        }
                    }
                }
            }
            let svg = plot::svg(&plot::Plot {
                title,
                x_label: xlabel.unwrap_or_else(|| x.clone()),
                y_label: ylabel.unwrap_or_else(|| y.join(", ")),
                series,
            })?;
            match o {
                Some(p) => std::fs::write(&p, svg).map_err(|e| format!("{}: {e}", p.display()))?,
                None => out.write_all(svg.as_bytes()).map_err(stdout_err)?,
            }
        }
    }
    Ok(failed)
}

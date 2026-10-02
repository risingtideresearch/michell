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
mod records;
mod run;
mod store;

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
    /// Cut hull files (IGES, STL) into sections: a hull record per file and
    /// setting.
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

/// A field path's last segment, to name a series by: `span` for
/// `study.case.params.span`, `|heave|` for `|seakeeping.heave|`.
fn short_field(f: &str) -> String {
    match f.strip_prefix('|').and_then(|p| p.strip_suffix('|')) {
        Some(inner) => format!("|{}|", short_field(inner)),
        None => f.rsplit('.').next().unwrap_or(f).to_string(),
    }
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
    let mut emit = |r: &Value| -> Result<(), String> {
        writeln!(out, "{r}").map_err(|e| format!("stdout: {e}"))
    };
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
                writeln!(out, "{}", line(fields.clone())).map_err(|e| e.to_string())?;
            }
            for r in read_records()? {
                for row in res.rows(&r, explode.as_deref()) {
                    let cells = fields
                        .iter()
                        .map(|f| path::cell(&res.get_in(&row, &r, f)))
                        .collect();
                    writeln!(out, "{}", line(cells)).map_err(|e| e.to_string())?;
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
                None => out.write_all(svg.as_bytes()).map_err(|e| e.to_string())?,
            }
        }
    }
    Ok(failed)
}

//! Hull file formats and input sniffing.
//!
//! Two line-oriented text formats (SI units, `#` comments):
//!
//! **`.hull` — canonical control net** (what the core consumes):
//! ```text
//! michell-hull v1
//! degree-x 2
//! degree-z 2
//! knots-x -5 -5 -5 5 5 5
//! knots-z 0 0 0 0.625 0.625 0.625
//! row 0 0 0        # one line per x control index, n_ctrl_z values each
//! row 1 1 0
//! row 0 0 0
//! ```
//!
//! **offsets — station × waterline half-beam table** (lofted on load):
//! ```text
//! michell-offsets v1
//! waterlines 0 0.125 0.25 0.375 0.5 0.625
//! station -5.0   0 0 0 0 0 0
//! station -4.0   0.36 0.34 0.30 0.24 0.15 0
//! ...
//! ```
//!
//! IGES files (`.igs`/`.iges`, or sniffed by the section letter in column
//! 73) are imported via `michell::iges`.

use michell::body::{Body, BodyOptions};
use michell::fit::{fit_offsets, FitOptions, FitReport};
use michell::iges::{self, HullPose, ImportOptions, ImportReport, Platform};
use michell::{BSplineSurface, Hull, Placement};

/// Where a hull came from, with any conversion diagnostics.
#[derive(Clone)]
pub enum Source {
    /// Native wetted control-net file: exact, no fit involved.
    Native,
    /// A full-band body file, situated at its design waterline.
    Body(FitReport),
    /// Lofted from an offsets table.
    Offsets(FitReport),
    /// Imported from IGES.
    Iges(ImportReport),
}

/// Import/loft settings shared by every command that reads a hull.
pub struct LoadSettings {
    pub waterline_z: f64,
    pub centerplane: Option<f64>,
    pub samples: (usize, usize),
    pub fit: FitOptions,
    /// True when the user set fit options explicitly (otherwise offsets
    /// lofting adapts the control count to the grid).
    pub fit_explicit: bool,
}

impl Default for LoadSettings {
    fn default() -> Self {
        LoadSettings {
            waterline_z: 0.0,
            centerplane: None,
            samples: (121, 33),
            fit: FitOptions::default(),
            fit_explicit: false,
        }
    }
}

/// Load every hull contained in a file. Native control nets and offsets
/// tables hold one hull at the default placement; an IGES file may contain a
/// whole multihull, each member carrying its detected placement.
pub fn load_hulls(
    path: &str,
    settings: &LoadSettings,
) -> Result<Vec<(Hull, Placement, Source)>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let first = text
        .lines()
        .find(|l| !l.trim().is_empty())
        .unwrap_or("")
        .trim_end();
    if first.starts_with("michell-hull") {
        let data = parse_hull_data(&text)?;
        let y = data.centerplane.unwrap_or(0.0);
        return match data.waterline {
            // Full-band body: situate at its design waterline.
            Some(wl) => {
                let body = Body::new(data.surface, wl, y).map_err(|e| format!("{e}"))?;
                let bopts = body_options(settings);
                let situated = body
                    .situate(0.0, &HullPose::default(), &Platform::default(), &bopts)
                    .map_err(|e| format!("{e}"))?
                    .ok_or_else(|| format!("{path}: body is dry at its design waterline"))?;
                Ok(vec![(situated.hull, situated.placement, Source::Body(situated.fit))])
            }
            None => {
                let hull = Hull::new(data.surface).map_err(|e| format!("{e}"))?;
                Ok(vec![(
                    hull,
                    Placement { x: 0.0, y },
                    Source::Native,
                )])
            }
        };
    }
    if first.starts_with("michell-offsets") {
        let (st, wl, y) = parse_offsets_file(&text)?;
        let mut fit = settings.fit;
        if !settings.fit_explicit {
            // Adapt the default control count to the grid so small tables
            // still loft.
            fit.n_ctrl_x = fit.n_ctrl_x.min(st.len().saturating_sub(2)).max(fit.degree_x + 1);
            fit.n_ctrl_z = fit.n_ctrl_z.min(wl.len().saturating_sub(2)).max(fit.degree_z + 1);
        }
        let (hull, report) =
            fit_offsets(&st, &wl, &y, &fit).map_err(|e| format!("loft failed: {e}"))?;
        return Ok(vec![(hull, Placement::default(), Source::Offsets(report))]);
    }
    let lower = path.to_ascii_lowercase();
    let looks_iges = lower.ends_with(".igs")
        || lower.ends_with(".iges")
        || first.len() >= 73 && matches!(first.as_bytes()[72], b'S' | b'G');
    if looks_iges {
        let opts = ImportOptions {
            waterline_z: settings.waterline_z,
            stations: settings.samples.0,
            waterlines: settings.samples.1,
            // Unless set explicitly, use the denser IGES default net rather
            // than the offsets-table default.
            fit: if settings.fit_explicit {
                settings.fit
            } else {
                ImportOptions::default().fit
            },
            centerplane: settings.centerplane,
        };
        let fleet =
            iges::import_fleet(&text, &opts).map_err(|e| format!("IGES import failed: {e}"))?;
        return Ok(fleet
            .into_iter()
            .map(|m| (m.hull, m.placement, Source::Iges(m.report)))
            .collect());
    }
    Err(format!(
        "cannot determine the format of {path}: expected a `michell-hull v1` or \
         `michell-offsets v1` header, or an IGES file"
    ))
}

fn strip_comment(line: &str) -> &str {
    match line.find('#') {
        Some(i) => &line[..i],
        None => line,
    }
}

fn parse_floats(s: &str, what: &str) -> Result<Vec<f64>, String> {
    s.split_whitespace()
        .map(|t| {
            t.parse::<f64>()
                .map_err(|_| format!("{what}: cannot parse number {t:?}"))
        })
        .collect()
}

/// Contents of a `.hull` file: the spline plus optional body metadata.
pub struct HullFileData {
    pub surface: BSplineSurface,
    /// Present = full-band body: depth of the design WL below the band top.
    pub waterline: Option<f64>,
    pub centerplane: Option<f64>,
}

pub fn parse_hull_data(text: &str) -> Result<HullFileData, String> {
    let mut degree_x: Option<usize> = None;
    let mut degree_z: Option<usize> = None;
    let mut knots_x: Option<Vec<f64>> = None;
    let mut knots_z: Option<Vec<f64>> = None;
    let mut waterline: Option<f64> = None;
    let mut centerplane: Option<f64> = None;
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let mut saw_header = false;
    for (ln, raw) in text.lines().enumerate() {
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        let at = |msg: String| format!("line {}: {msg}", ln + 1);
        let (key, rest) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
        match key {
            "michell-hull" => {
                if rest.trim() != "v1" {
                    return Err(at(format!("unsupported version {:?}", rest.trim())));
                }
                saw_header = true;
            }
            "degree-x" => degree_x = Some(parse_usize(rest).map_err(&at)?),
            "degree-z" => degree_z = Some(parse_usize(rest).map_err(&at)?),
            "knots-x" => knots_x = Some(parse_floats(rest, "knots-x").map_err(&at)?),
            "knots-z" => knots_z = Some(parse_floats(rest, "knots-z").map_err(&at)?),
            "waterline" => waterline = Some(parse_f64(rest).map_err(&at)?),
            "centerplane" => centerplane = Some(parse_f64(rest).map_err(&at)?),
            "row" => rows.push(parse_floats(rest, "row").map_err(&at)?),
            other => return Err(at(format!("unknown key {other:?}"))),
        }
    }
    if !saw_header {
        return Err("missing `michell-hull v1` header".into());
    }
    let (Some(px), Some(pz), Some(kx), Some(kz)) = (degree_x, degree_z, knots_x, knots_z) else {
        return Err("missing one of: degree-x, degree-z, knots-x, knots-z".into());
    };
    let nx = kx.len().saturating_sub(px + 1);
    let nz = kz.len().saturating_sub(pz + 1);
    if rows.len() != nx {
        return Err(format!(
            "expected {nx} `row` lines (knots-x implies {nx} control rows), got {}",
            rows.len()
        ));
    }
    let mut control = Vec::with_capacity(nx * nz);
    for (i, r) in rows.iter().enumerate() {
        if r.len() != nz {
            return Err(format!(
                "row {}: expected {nz} values (knots-z implies {nz} control columns), got {}",
                i + 1,
                r.len()
            ));
        }
        control.extend_from_slice(r);
    }
    let surface =
        BSplineSurface::new(px, pz, kx, kz, control).map_err(|e| format!("{e}"))?;
    Ok(HullFileData {
        surface,
        waterline,
        centerplane,
    })
}

/// Load a full-band body from a `.hull` file (requires the `waterline` key).
pub fn load_body(path: &str) -> Result<Body, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let data = parse_hull_data(&text).map_err(|e| format!("{path}: {e}"))?;
    let wl = data.waterline.ok_or_else(|| {
        format!(
            "{path} is a wetted-only hull (no `waterline` key); re-run \
             `michell loft` on the source geometry to produce a full-band body"
        )
    })?;
    Body::new(data.surface, wl, data.centerplane.unwrap_or(0.0)).map_err(|e| format!("{path}: {e}"))
}

/// Body sampling options derived from the CLI load settings.
pub fn body_options(settings: &LoadSettings) -> BodyOptions {
    BodyOptions {
        stations: settings.samples.0,
        waterlines: settings.samples.1,
        fit: if settings.fit_explicit {
            settings.fit
        } else {
            BodyOptions::default().fit
        },
    }
}

fn parse_f64(s: &str) -> Result<f64, String> {
    s.trim()
        .parse::<f64>()
        .map_err(|_| format!("cannot parse number {:?}", s.trim()))
}

fn parse_usize(s: &str) -> Result<usize, String> {
    s.trim()
        .parse::<usize>()
        .map_err(|_| format!("cannot parse integer {:?}", s.trim()))
}

pub fn write_hull_file(hull: &Hull) -> String {
    write_spline_file(hull.surface(), None, None)
}

/// Write a full-band body file (`waterline` = design WL depth below the band
/// top, `centerplane` = transverse position).
pub fn write_body_file(surface: &BSplineSurface, waterline: f64, centerplane: f64) -> String {
    write_spline_file(surface, Some(waterline), Some(centerplane))
}

fn write_spline_file(
    s: &BSplineSurface,
    waterline: Option<f64>,
    centerplane: Option<f64>,
) -> String {
    let join = |v: &[f64]| {
        v.iter()
            .map(|x| format!("{x}"))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let mut out = String::from("michell-hull v1\n");
    if let Some(w) = waterline {
        out.push_str(&format!("waterline {w}\n"));
    }
    if let Some(c) = centerplane {
        out.push_str(&format!("centerplane {c}\n"));
    }
    out.push_str(&format!("degree-x {}\n", s.degree_x()));
    out.push_str(&format!("degree-z {}\n", s.degree_z()));
    out.push_str(&format!("knots-x {}\n", join(s.knots_x())));
    out.push_str(&format!("knots-z {}\n", join(s.knots_z())));
    let nz = s.n_ctrl_z();
    for i in 0..s.n_ctrl_x() {
        out.push_str(&format!("row {}\n", join(&s.control()[i * nz..(i + 1) * nz])));
    }
    out
}

/// Stations, waterlines, and half-beams (row-major, stations-major).
pub type OffsetsTable = (Vec<f64>, Vec<f64>, Vec<f64>);

pub fn parse_offsets_file(text: &str) -> Result<OffsetsTable, String> {
    let mut waterlines: Option<Vec<f64>> = None;
    let mut stations: Vec<f64> = Vec::new();
    let mut grid: Vec<f64> = Vec::new();
    let mut saw_header = false;
    for (ln, raw) in text.lines().enumerate() {
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        let at = |msg: String| format!("line {}: {msg}", ln + 1);
        let (key, rest) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
        match key {
            "michell-offsets" => {
                if rest.trim() != "v1" {
                    return Err(at(format!("unsupported version {:?}", rest.trim())));
                }
                saw_header = true;
            }
            "waterlines" => waterlines = Some(parse_floats(rest, "waterlines").map_err(&at)?),
            "station" => {
                let vals = parse_floats(rest, "station").map_err(&at)?;
                let Some(wl) = &waterlines else {
                    return Err(at("`waterlines` must come before `station` lines".into()));
                };
                if vals.len() != wl.len() + 1 {
                    return Err(at(format!(
                        "expected x plus {} half-beams, got {} values",
                        wl.len(),
                        vals.len()
                    )));
                }
                stations.push(vals[0]);
                grid.extend_from_slice(&vals[1..]);
            }
            other => return Err(at(format!("unknown key {other:?}"))),
        }
    }
    if !saw_header {
        return Err("missing `michell-offsets v1` header".into());
    }
    let wl = waterlines.ok_or("missing `waterlines` line")?;
    if stations.len() < 2 {
        return Err("need at least two `station` lines".into());
    }
    Ok((stations, wl, grid))
}

/// Parse `a`, or an inclusive range `a:b:step`.
pub fn parse_range(s: &str) -> Result<Vec<f64>, String> {
    let parts: Vec<&str> = s.split(':').collect();
    let num = |t: &str| {
        t.parse::<f64>()
            .map_err(|_| format!("cannot parse number {t:?} in {s:?}"))
    };
    match parts.as_slice() {
        [one] => Ok(vec![num(one)?]),
        [a, b, step] => {
            let (a, b, step) = (num(a)?, num(b)?, num(step)?);
            if !(step > 0.0 && b >= a) {
                return Err(format!("bad range {s:?}: need start <= end and step > 0"));
            }
            let n = ((b - a) / step + 1e-9).floor() as usize;
            Ok((0..=n).map(|i| a + i as f64 * step).collect())
        }
        _ => Err(format!("bad range {s:?}: expected `value` or `start:end:step`")),
    }
}

/// Parse `NxM` (e.g. `12x8`).
pub fn parse_pair(s: &str) -> Result<(usize, usize), String> {
    let (a, b) = s
        .split_once(['x', 'X'])
        .ok_or_else(|| format!("expected NxM, got {s:?}"))?;
    Ok((
        a.parse().map_err(|_| format!("bad integer {a:?}"))?,
        b.parse().map_err(|_| format!("bad integer {b:?}"))?,
    ))
}

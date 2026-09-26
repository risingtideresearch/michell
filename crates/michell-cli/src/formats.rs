//! Hull file formats and argument parsing helpers.
//!
//! **`.hull` — B-spline half-breadth control net** (SI units, `#` comments):
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
//! An optional `waterline D` key marks a full-band body: its top is a band
//! top, with the design waterline `D` below it; `centerplane Y` places it.
//! Loaded hulls are converted to their exact surfaces and cut into sections
//! ([`crate::fleet`]).

use michell::Roughness;
use michell::{BSplineSurface, Hull};

pub fn parse_length(s: &str) -> Result<f64, String> {
    let t = s.trim();
    let (num, scale) = [
        ("um", 1e-6),
        ("µm", 1e-6),
        ("mm", 1e-3),
        ("cm", 1e-2),
        ("ft", 0.3048),
        ("in", 0.0254),
        ("m", 1.0),
    ]
    .iter()
    .find_map(|(suf, k)| t.strip_suffix(suf).map(|n| (n, *k)))
    .unwrap_or((t, 1.0));
    num.trim()
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v >= 0.0)
        .map(|v| v * scale)
        .ok_or_else(|| format!("cannot parse length {t:?} (try 100um, 0.1mm, or 1e-4)"))
}

/// Parse a `--roughness` spec: `off`, `cf=VALUE` (a ΔC_F added outside the
/// form factor), or `ks=LENGTH` (equivalent sand-grain height).
pub fn parse_roughness(spec: &str) -> Result<Roughness, String> {
    let t = spec.trim();
    if matches!(t, "off" | "none" | "smooth") {
        return Ok(Roughness::None);
    }
    if let Some(v) = t.strip_prefix("cf=") {
        let c = v
            .trim()
            .parse::<f64>()
            .ok()
            .filter(|x| x.is_finite() && *x >= 0.0)
            .ok_or_else(|| format!("--roughness cf={v:?}: expected a non-negative number"))?;
        return Ok(Roughness::DeltaCf(c));
    }
    if let Some(v) = t.strip_prefix("ks=") {
        return Ok(Roughness::SandGrain(parse_length(v)?));
    }
    Err(format!(
        "--roughness {t:?}: expected off | cf=DELTA_CF | ks=HEIGHT (e.g. ks=100um)"
    ))
}

/// Parse a units name (or raw scale) to metres-per-unit.
pub fn parse_units(s: &str) -> Result<f64, String> {
    match s.trim() {
        "mm" => Ok(0.001),
        "cm" => Ok(0.01),
        "m" => Ok(1.0),
        "in" => Ok(0.0254),
        "ft" => Ok(0.3048),
        other => other
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite() && *v > 0.0)
            .ok_or_else(|| format!("--units {other:?}: expected mm|cm|m|in|ft or a scale")),
    }
}

/// Binary-STL detection: exact size match on the triangle count.
pub fn looks_binary_stl(bytes: &[u8]) -> bool {
    bytes.len() >= 84 && {
        let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        bytes.len() == 84 + 50 * n
    }
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
    let surface = BSplineSurface::new(px, pz, kx, kz, control).map_err(|e| format!("{e}"))?;
    Ok(HullFileData {
        surface,
        waterline,
        centerplane,
    })
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
        out.push_str(&format!(
            "row {}\n",
            join(&s.control()[i * nz..(i + 1) * nz])
        ));
    }
    out
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
        _ => Err(format!(
            "bad range {s:?}: expected `value` or `start:end:step`"
        )),
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

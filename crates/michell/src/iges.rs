//! Minimal, dependency-free IGES reader and hull importer.
//!
//! Scope: **untrimmed rational B-spline surfaces** (entity 128) — one or many
//! patches — with optional transformation matrices (entity 124) and unit
//! conversion from the global section. Bounded-surface wrappers (143/141), as
//! produced by SubD/T-spline → NURBS exports, are tolerated on the assumption
//! that the boundaries are the natural patch rectangles; genuinely *trimmed*
//! surfaces (entities 142/144) are rejected. This is a hull-surface importer,
//! not a CAD kernel.
//!
//! ## Import pipeline
//!
//! A CAD surface is parametric, `S(u,v) → (x, y, z)`; the crate's core wants a
//! height field `y = f(x, z)`. The importer therefore:
//!
//! 1. parses every 128 entity (all weights must be uniform — the polynomial
//!    contract; see crate docs),
//! 2. maps patches to the hull frame (`x` longitudinal as in the file, `z'`
//!    downward from the user-specified waterline),
//! 3. finds the hull's **centerplane**: given via [`ImportOptions::centerplane`],
//!    or auto-detected — interior probes count shell intersections along y;
//!    a two-sided (full) shell folds about the midplane of its intersections,
//!    a one-sided (half-hull) file measures from y = 0,
//! 4. clips to the wetted region `z' ∈ [0, T]` and samples half-beams on a
//!    station × waterline grid by per-patch 2-D Newton inversion of
//!    `(x(u,v), z'(u,v))`, taking the outermost fold `max |y - y_c|`; the
//!    surface slopes `∂y/∂x`, `∂y/∂z` come for free from the converged
//!    Newton Jacobian (implicit function theorem),
//! 5. lofts the resulting [`crate::grid::SampleGrid`] — values, slopes, and
//!    per-sample weights (failed inversions are excluded, not zeroed) — with
//!    [`crate::fit::fit_grid`].
//!
//! Expected CAD frame: `x` longitudinal, `z` **up**, `y` transverse.
//! `waterline_z` gives the design waterline height in the file's frame, **in
//! metres** (after unit conversion). Diagnostics (fold asymmetry, ambiguous
//! samples, failed inversions, loft residuals) are reported so a bad import
//! is visible.

use crate::bspline::{ders_basis, find_span};
use crate::error::{Error, Result};
use crate::fit::{fit_grid, FitOptions, FitReport};
use crate::grid::SampleGrid;
use crate::hull::Hull;
use crate::michell::Placement;

// ---------------------------------------------------------------------------
// Parsed geometry
// ---------------------------------------------------------------------------

/// An untrimmed NURBS surface from an IGES entity 128, with any 124 transform
/// applied and coordinates scaled to metres.
#[derive(Debug, Clone)]
pub struct NurbsSurface3 {
    pub degree_u: usize,
    pub degree_v: usize,
    pub knots_u: Vec<f64>,
    pub knots_v: Vec<f64>,
    pub n_ctrl_u: usize,
    pub n_ctrl_v: usize,
    /// Row-major, v fastest: `ctrl[iu * n_ctrl_v + iv]`.
    pub ctrl: Vec<[f64; 3]>,
    /// Same layout as `ctrl`.
    pub weights: Vec<f64>,
}

impl NurbsSurface3 {
    /// True when all weights are equal (the surface is polynomial).
    pub fn is_polynomial(&self) -> bool {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &w in &self.weights {
            lo = lo.min(w);
            hi = hi.max(w);
        }
        hi - lo <= 1e-9 * hi.abs().max(1.0)
    }

    pub fn u_domain(&self) -> (f64, f64) {
        (self.knots_u[self.degree_u], self.knots_u[self.n_ctrl_u])
    }

    pub fn v_domain(&self) -> (f64, f64) {
        (self.knots_v[self.degree_v], self.knots_v[self.n_ctrl_v])
    }

    /// Point and first partials. Assumes uniform weights (polynomial).
    #[allow(clippy::needless_range_loop)]
    fn eval1(&self, u: f64, v: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let (pu, pv) = (self.degree_u, self.degree_v);
        let su = find_span(&self.knots_u, pu, self.n_ctrl_u, u);
        let sv = find_span(&self.knots_v, pv, self.n_ctrl_v, v);
        let du = ders_basis(&self.knots_u, pu, su, u, 1);
        let dv = ders_basis(&self.knots_v, pv, sv, v, 1);
        let mut out = [[0.0f64; 3]; 3]; // S, Su, Sv
        for i in 0..=pu {
            let ci = su - pu + i;
            for j in 0..=pv {
                let p = self.ctrl[ci * self.n_ctrl_v + (sv - pv + j)];
                let b00 = du[0][i] * dv[0][j];
                let b10 = du[1][i] * dv[0][j];
                let b01 = du[0][i] * dv[1][j];
                for c in 0..3 {
                    out[0][c] += b00 * p[c];
                    out[1][c] += b10 * p[c];
                    out[2][c] += b01 * p[c];
                }
            }
        }
        (out[0], out[1], out[2])
    }
}

/// Result of parsing an IGES file.
#[derive(Debug, Clone)]
pub struct IgesFile {
    pub surfaces: Vec<NurbsSurface3>,
    /// Multiplier applied to convert file units to metres.
    pub units_scale: f64,
    /// Entity type -> count, for everything in the directory section.
    pub entity_counts: Vec<(i64, usize)>,
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

/// Parse an IGES file, extracting all untrimmed 128 surfaces.
pub fn parse(text: &str) -> Result<IgesFile> {
    let mut g_lines: Vec<String> = Vec::new();
    let mut d_lines: Vec<String> = Vec::new();
    let mut p_lines: Vec<String> = Vec::new();
    for raw in text.lines() {
        if raw.trim().is_empty() {
            continue;
        }
        let line: String = if raw.len() < 80 {
            format!("{raw:<80}")
        } else {
            raw.to_string()
        };
        match line.as_bytes().get(72) {
            Some(b'S') | Some(b'T') => {}
            Some(b'G') => g_lines.push(line),
            Some(b'D') => d_lines.push(line),
            Some(b'P') => p_lines.push(line),
            _ => {
                return Err(Error::Parse(format!(
                    "line without a valid section letter in column 73: {:?}",
                    raw.get(..73.min(raw.len())).unwrap_or(raw)
                )))
            }
        }
    }
    if d_lines.is_empty() || p_lines.is_empty() {
        return Err(Error::Parse(
            "missing directory or parameter section; not an IGES file?".into(),
        ));
    }
    if !d_lines.len().is_multiple_of(2) {
        return Err(Error::Parse(
            "directory section must contain an even number of lines".into(),
        ));
    }

    let units_scale = parse_units(&g_lines)?;

    // Directory entries: two 80-column lines per entity, 8-column fields.
    struct Dir {
        etype: i64,
        pd_ptr: usize,
        pd_count: usize,
        transform_de: usize, // DE sequence number (0 = none)
    }
    fn int_field(line: &str, k: usize) -> i64 {
        line[8 * k..8 * (k + 1)].trim().parse::<i64>().unwrap_or(0)
    }
    let mut dirs = Vec::new();
    for pair in d_lines.chunks(2) {
        dirs.push(Dir {
            etype: int_field(&pair[0], 0),
            pd_ptr: int_field(&pair[0], 1).max(0) as usize,
            transform_de: int_field(&pair[0], 6).max(0) as usize,
            pd_count: int_field(&pair[1], 3).max(0) as usize,
        });
    }

    let mut counts: Vec<(i64, usize)> = Vec::new();
    for d in &dirs {
        match counts.iter_mut().find(|(t, _)| *t == d.etype) {
            Some((_, c)) => *c += 1,
            None => counts.push((d.etype, 1)),
        }
    }

    // Parameter data for one entity: concatenate data columns (1..=64) of
    // pd_count lines starting at pd_ptr (1-based line number in P section).
    let entity_params = |pd_ptr: usize, pd_count: usize| -> Result<Vec<f64>> {
        if pd_ptr == 0 || pd_ptr + pd_count - 1 > p_lines.len() {
            return Err(Error::Parse(format!(
                "parameter data pointer {pd_ptr}+{pd_count} out of range"
            )));
        }
        let mut text = String::new();
        for line in &p_lines[pd_ptr - 1..pd_ptr - 1 + pd_count] {
            text.push_str(&line[..64]);
        }
        let body = text.split(';').next().unwrap_or("");
        body.split(',').map(parse_number).collect()
    };

    // Transformation matrices (entity 124), keyed by DE sequence number.
    let mut transforms: Vec<(usize, [f64; 12])> = Vec::new();
    for (idx, d) in dirs.iter().enumerate() {
        if d.etype == 124 {
            let params = entity_params(d.pd_ptr, d.pd_count)?;
            if params.len() < 13 {
                return Err(Error::Parse(
                    "transformation matrix (124) with fewer than 12 parameters".into(),
                ));
            }
            let mut m = [0.0f64; 12];
            m.copy_from_slice(&params[1..13]);
            transforms.push((2 * idx + 1, m));
        }
    }

    let mut surfaces = Vec::new();
    for d in dirs.iter().filter(|d| d.etype == 128) {
        let p = entity_params(d.pd_ptr, d.pd_count)?;
        let mut surf = parse_surface_128(&p)?;
        if d.transform_de != 0 {
            let m = transforms
                .iter()
                .find(|(de, _)| *de == d.transform_de)
                .map(|(_, m)| *m)
                .ok_or_else(|| {
                    Error::Parse(format!(
                        "entity references transformation matrix DE {} which is not a 124 entity",
                        d.transform_de
                    ))
                })?;
            for p in surf.ctrl.iter_mut() {
                let q = *p;
                for r in 0..3 {
                    p[r] = m[4 * r] * q[0] + m[4 * r + 1] * q[1] + m[4 * r + 2] * q[2]
                        + m[4 * r + 3];
                }
            }
        }
        for p in surf.ctrl.iter_mut() {
            for c in p.iter_mut() {
                *c *= units_scale;
            }
        }
        surfaces.push(surf);
    }

    Ok(IgesFile {
        surfaces,
        units_scale,
        entity_counts: counts,
    })
}

/// Entity 128 parameter list -> surface (file units, no transform).
fn parse_surface_128(p: &[f64]) -> Result<NurbsSurface3> {
    let need = |cond: bool, what: &str| -> Result<()> {
        if cond {
            Ok(())
        } else {
            Err(Error::Parse(format!("malformed 128 entity: {what}")))
        }
    };
    need(p.len() >= 10, "fewer than 10 parameters")?;
    need(p[0] as i64 == 128, "entity type mismatch")?;
    let k1 = p[1] as i64;
    let k2 = p[2] as i64;
    let m1 = p[3] as i64;
    let m2 = p[4] as i64;
    need(k1 >= 1 && k2 >= 1 && m1 >= 1 && m2 >= 1, "bad indices")?;
    need(k1 >= m1 && k2 >= m2, "fewer control points than degree + 1")?;
    let (nu, nv) = (k1 as usize + 1, k2 as usize + 1);
    let (pu, pv) = (m1 as usize, m2 as usize);
    let nku = nu + pu + 1;
    let nkv = nv + pv + 1;
    // 1 entity type + 9 header parameters (K1, K2, M1, M2, PROP1..PROP5).
    let total = 10 + nku + nkv + nu * nv + 3 * nu * nv + 4;
    need(
        p.len() >= total,
        &format!("expected at least {total} parameters, got {}", p.len()),
    )?;
    let mut at = 10usize;
    let knots_u = p[at..at + nku].to_vec();
    at += nku;
    let knots_v = p[at..at + nkv].to_vec();
    at += nkv;
    // Weights and points are stored with the u index varying fastest.
    let mut weights = vec![0.0f64; nu * nv];
    for j in 0..nv {
        for i in 0..nu {
            weights[i * nv + j] = p[at];
            at += 1;
        }
    }
    let mut ctrl = vec![[0.0f64; 3]; nu * nv];
    for j in 0..nv {
        for i in 0..nu {
            ctrl[i * nv + j] = [p[at], p[at + 1], p[at + 2]];
            at += 3;
        }
    }
    if weights.iter().any(|&w| !(w.is_finite() && w > 0.0)) {
        return Err(Error::Parse("128 entity has non-positive weights".into()));
    }
    if ctrl.iter().flatten().any(|v| !v.is_finite()) {
        return Err(Error::Parse(
            "128 entity has non-finite control points".into(),
        ));
    }
    for (knots, deg, n, dir) in [(&knots_u, pu, nu, "u"), (&knots_v, pv, nv, "v")] {
        if knots.iter().any(|k| !k.is_finite())
            || knots.windows(2).any(|w| w[1] < w[0])
            || knots[deg] >= knots[n]
        {
            return Err(Error::Parse(format!(
                "128 entity has an invalid knot vector in {dir}"
            )));
        }
    }
    Ok(NurbsSurface3 {
        degree_u: pu,
        degree_v: pv,
        knots_u,
        knots_v,
        n_ctrl_u: nu,
        n_ctrl_v: nv,
        ctrl,
        weights,
    })
}

/// IGES numbers may use FORTRAN D-exponents; blank fields read as 0.
fn parse_number(tok: &str) -> Result<f64> {
    let t = tok.trim();
    if t.is_empty() {
        return Ok(0.0);
    }
    let t = t.replace(['D', 'd'], "E");
    t.parse::<f64>()
        .map_err(|_| Error::Parse(format!("cannot parse number {tok:?}")))
}

/// Units scale (to metres) from the global section, field 14 (units flag)
/// with field 15 (units name) as fallback for flag 3.
fn parse_units(g_lines: &[String]) -> Result<f64> {
    let mut text = String::new();
    for line in g_lines {
        text.push_str(&line[..72]);
    }
    let fields = split_global(&text);
    let flag = fields
        .get(13)
        .and_then(|s| s.trim().parse::<i64>().ok())
        .ok_or_else(|| Error::Parse("global section lacks a units flag (field 14)".into()))?;
    let name = fields.get(14).cloned().unwrap_or_default();
    let scale = match flag {
        1 => 0.0254,       // inches
        2 => 0.001,        // millimetres
        3 => match name.trim().to_ascii_uppercase().as_str() {
            "M" | "METER" | "METERS" | "METRE" | "METRES" => 1.0,
            "MM" | "MILLIMETER" | "MILLIMETERS" => 0.001,
            "CM" | "CENTIMETER" | "CENTIMETERS" => 0.01,
            "IN" | "INCH" | "INCHES" => 0.0254,
            "FT" | "FOOT" | "FEET" => 0.3048,
            other => {
                return Err(Error::Unsupported(format!(
                    "units specified by name {other:?} are not recognised"
                )))
            }
        },
        4 => 0.3048,       // feet
        5 => 1609.344,     // miles
        6 => 1.0,          // metres
        7 => 1000.0,       // kilometres
        8 => 2.54e-5,      // mils
        9 => 1e-6,         // microns
        10 => 0.01,        // centimetres
        11 => 2.54e-8,     // microinches
        other => {
            return Err(Error::Unsupported(format!(
                "unknown IGES units flag {other}"
            )))
        }
    };
    Ok(scale)
}

/// Split the global section on the parameter delimiter, honouring Hollerith
/// strings (`nHxxxx`). The delimiters themselves are read from the leading
/// fields (defaults `,` and `;`).
fn split_global(text: &str) -> Vec<String> {
    let bytes = text.as_bytes();
    // Leading delimiter definitions: either ",," (defaults) or "1Hx,1Hy,".
    let mut param_delim = b',';
    let mut i = 0usize;
    if text.starts_with("1H") && bytes.len() > 3 {
        param_delim = bytes[2];
        i = 4; // skip "1Hx" + delimiter
    } else if bytes.first() == Some(&b',') {
        i = 1;
    }
    // Record delimiter field: "1Hx" or empty.
    if bytes.get(i) == Some(&b'1') && bytes.get(i + 1) == Some(&b'H') {
        i += 3;
        if bytes.get(i) == Some(&param_delim) {
            i += 1;
        }
    } else if bytes.get(i) == Some(&param_delim) {
        i += 1;
    }
    let mut fields = vec![String::new(), String::new()]; // the two delim fields
    let mut cur = String::new();
    while i < bytes.len() {
        // Hollerith: digits followed by 'H'.
        let mut j = i;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
        if j > i && bytes.get(j) == Some(&b'H') {
            if let Ok(len) = text[i..j].parse::<usize>() {
                let start = j + 1;
                let end = (start + len).min(bytes.len());
                cur.push_str(&text[start..end]);
                i = end;
                continue;
            }
        }
        let c = bytes[i];
        if c == param_delim || c == b';' {
            fields.push(std::mem::take(&mut cur));
            if c == b';' {
                return fields;
            }
        } else {
            cur.push(c as char);
        }
        i += 1;
    }
    fields.push(cur);
    fields
}

// ---------------------------------------------------------------------------
// Hull import
// ---------------------------------------------------------------------------

/// Options controlling the IGES → hull conversion.
#[derive(Debug, Clone, Copy)]
pub struct ImportOptions {
    /// Height of the design waterline in the file's frame (z up), in
    /// **metres** (i.e. after unit conversion).
    pub waterline_z: f64,
    /// Number of sample stations along x.
    pub stations: usize,
    /// Number of sample waterlines over the draft.
    pub waterlines: usize,
    pub fit: FitOptions,
    /// Transverse position of the hull's centerplane in the file frame [m].
    /// `None` auto-detects: a two-sided (full) shell folds about the midplane
    /// of its shell intersections; a one-sided (half-hull) file measures from
    /// y = 0.
    pub centerplane: Option<f64>,
}

impl Default for ImportOptions {
    fn default() -> Self {
        ImportOptions {
            waterline_z: 0.0,
            stations: 121,
            waterlines: 33,
            // Real CAD hulls carry more shape than hand-typed offset tables;
            // a denser net keeps low-Froude wave resistance converged.
            fit: FitOptions {
                degree_x: 3,
                degree_z: 3,
                n_ctrl_x: 20,
                n_ctrl_z: 12,
                ..FitOptions::default()
            },
            centerplane: None,
        }
    }
}

/// Import diagnostics.
#[derive(Debug, Clone, Copy)]
pub struct ImportReport {
    pub units_scale: f64,
    /// Number of surface patches used.
    pub patches: usize,
    /// True if the file contained a full (both-sided) shell that was folded
    /// about the centerplane.
    pub two_sided: bool,
    /// Centerplane y used for folding (auto-detected or user-supplied) [m].
    pub centerplane: f64,
    /// True if a one-sided file contained the port half (y < y_c).
    pub mirrored: bool,
    /// Draft found below the waterline [m].
    pub draft: f64,
    /// Longitudinal extent of the wetted surface [m].
    pub x_range: (f64, f64),
    /// Largest port/starboard half-beam disagreement among folded sample
    /// pairs [m]. Large values mean the shell is not symmetric about the
    /// centerplane — or that patch rectangles extend past real trims.
    pub max_asymmetry: f64,
    /// Samples with more than two shell intersections (overlapping or
    /// interior geometry; a few near seams are harmless).
    pub ambiguous_samples: usize,
    /// Grid samples inside a waterline's footprint that could not be
    /// inverted onto any patch (excluded from the loft by zero weight).
    pub failed_inversions: usize,
    /// Converged samples whose surface slope was unavailable (tangent-vertical
    /// shell, boundary-clamped near-misses, the centerplane fold): their
    /// half-beam still constrains the loft, their slope simply does not.
    pub derivative_gaps: usize,
    pub fit: FitReport,
}

/// One patch, presampled in the hull frame (x, y, z-downward).
struct Patch {
    surf: NurbsSurface3,
    /// (u, v, x, y, z-downward), grid_n × grid_n row-major (u-major).
    pts: Vec<(f64, f64, f64, f64, f64)>,
    grid_n: usize,
    /// (x_lo, x_hi, z_lo, z_hi) over the whole patch presample.
    bbox: (f64, f64, f64, f64),
    /// (x_lo, x_hi, y_lo, y_hi, z_lo, z_hi) over the wetted presample only;
    /// `None` when the patch is entirely above the waterline.
    wet_box: Option<[f64; 6]>,
}

/// One hull detected in an IGES file, with its recovered placement.
#[derive(Debug, Clone)]
pub struct ImportedHull {
    pub hull: Hull,
    /// Placement recovered from the file: `x = 0` (the file's longitudinal
    /// coordinates are kept) and `y` = the hull's detected centerplane.
    pub placement: Placement,
    pub report: ImportReport,
    /// The derivative-augmented sample grid the hull was lofted from.
    pub grid: SampleGrid,
}

/// Per-hull **design** pose: how a hull is mounted relative to the platform.
/// Applied to the source geometry before the waterline clip, so all fields
/// change the wetted shape exactly (affine maps of the control nets).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct HullPose {
    /// Longitudinal shift [m].
    pub dx: f64,
    /// Transverse shift [m].
    pub dy: f64,
    /// Immersion shift [m]; positive lowers the hull (deeper).
    pub dz: f64,
    /// Pitch rotation [rad]; positive raises the hull's +x end.
    pub trim: f64,
    /// Pivot station for `trim` (default: the hull's x mid); the pivot height
    /// is the base waterline.
    pub pivot_x: Option<f64>,
}

/// Whole-platform **state**: rigid-body sinkage and pitch, normally solved
/// from a load case rather than chosen.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Platform {
    /// Additional immersion of the whole platform [m]; positive = deeper.
    pub sinkage: f64,
    /// Pitch rotation [rad]; positive raises the +x end.
    pub trim: f64,
    /// Pivot station for `trim`, at the effective waterline height.
    pub pivot_x: f64,
}

/// A situated fleet: the wetted hulls plus which source hulls were dry.
#[derive(Debug)]
pub struct SituatedFleet {
    /// Wetted members, in source-hull order.
    pub members: Vec<ImportedHull>,
    /// Indices (into the [`SourceFleet`]) of hulls entirely above the water.
    pub dry: Vec<usize>,
}

/// Parsed and clustered source geometry, kept in CAD coordinates so hulls can
/// be re-situated (waterline, immersion, trim, position) repeatedly.
pub struct SourceFleet {
    units_scale: f64,
    /// Patches of each detected hull, CAD frame, metres, sorted by y.
    hulls: Vec<Vec<NurbsSurface3>>,
}

/// Parse an IGES file and cluster its patches into hulls at a reference
/// waterline (use the deepest waterline you intend to sweep, so cluster
/// membership stays fixed).
pub fn source_fleet(text: &str, reference_waterline: f64) -> Result<SourceFleet> {
    let file = validated_file(text)?;
    let patches = presample_surfaces(&file.surfaces, reference_waterline);
    let mut global_wet: Option<[f64; 6]> = None;
    for p in &patches {
        if let Some(b) = p.wet_box {
            let g = global_wet.get_or_insert(b);
            for k in 0..3 {
                g[2 * k] = g[2 * k].min(b[2 * k]);
                g[2 * k + 1] = g[2 * k + 1].max(b[2 * k + 1]);
            }
        }
    }
    let Some(g) = global_wet else {
        return Err(Error::InvalidGeometry(
            "the surface lies entirely above the specified waterline".into(),
        ));
    };
    let scale = (g[1] - g[0]).max(g[3] - g[2]).max(g[5] - g[4]);
    if scale <= 0.0 || !scale.is_finite() {
        return Err(Error::InvalidGeometry(
            "the wetted part of the surface is degenerate".into(),
        ));
    }
    let mut clusters = cluster_patches(&patches, 0.01 * scale);
    // Deterministic order: by wetted-y midpoint.
    let key = |idxs: &Vec<usize>| -> f64 {
        let mids: Vec<f64> = idxs
            .iter()
            .filter_map(|&i| patches[i].wet_box.map(|b| (b[2] + b[3]) / 2.0))
            .collect();
        mids.iter().sum::<f64>() / mids.len().max(1) as f64
    };
    clusters.sort_by(|a, b| key(a).total_cmp(&key(b)));

    // Dry patches (decks, topsides above the reference waterline) belong to
    // *some* hull and must be retained — a deeper pose may wet them. Attach
    // each to the nearest cluster by box distance. Clusters themselves are
    // still formed from wetted geometry only, so dry structure cannot merge
    // two hulls.
    let cluster_boxes: Vec<[f64; 6]> = clusters
        .iter()
        .map(|idxs| {
            let mut b = [
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ];
            for &i in idxs {
                if let Some(w) = patches[i].wet_box {
                    for k in 0..3 {
                        b[2 * k] = b[2 * k].min(w[2 * k]);
                        b[2 * k + 1] = b[2 * k + 1].max(w[2 * k + 1]);
                    }
                }
            }
            b
        })
        .collect();
    let assigned: Vec<bool> = {
        let mut a = vec![false; patches.len()];
        for idxs in &clusters {
            for &i in idxs {
                a[i] = true;
            }
        }
        a
    };
    for (pi, patch) in patches.iter().enumerate() {
        if assigned[pi] {
            continue;
        }
        // Full 3-D bbox of the dry patch from its presample.
        let mut pb = [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ];
        for q in &patch.pts {
            pb[0] = pb[0].min(q.2);
            pb[1] = pb[1].max(q.2);
            pb[2] = pb[2].min(q.3);
            pb[3] = pb[3].max(q.3);
            pb[4] = pb[4].min(q.4);
            pb[5] = pb[5].max(q.4);
        }
        let dist = |a: &[f64; 6], b: &[f64; 6]| -> f64 {
            (0..3)
                .map(|k| {
                    let gap = (a[2 * k] - b[2 * k + 1]).max(b[2 * k] - a[2 * k + 1]).max(0.0);
                    gap * gap
                })
                .sum::<f64>()
        };
        let nearest = cluster_boxes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| dist(&pb, a).total_cmp(&dist(&pb, b)))
            .map(|(i, _)| i)
            .expect("at least one cluster");
        clusters[nearest].push(pi);
    }

    Ok(SourceFleet {
        units_scale: file.units_scale,
        hulls: clusters
            .iter()
            .map(|idxs| idxs.iter().map(|&i| file.surfaces[i].clone()).collect())
            .collect(),
    })
}

impl SourceFleet {
    /// Number of hulls detected.
    pub fn len(&self) -> usize {
        self.hulls.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hulls.is_empty()
    }

    pub fn units_scale(&self) -> f64 {
        self.units_scale
    }

    /// Situate the fleet: apply each hull's design pose and the platform
    /// state, clip at the effective waterline `waterline_z + sinkage`, and
    /// loft. Hulls that end up entirely dry are reported, not errors.
    pub fn situate(
        &self,
        waterline_z: f64,
        poses: &[HullPose],
        platform: &Platform,
        opts: &ImportOptions,
    ) -> Result<SituatedFleet> {
        if poses.len() != self.hulls.len() {
            return Err(Error::InvalidInput(format!(
                "{} poses supplied for {} hulls",
                poses.len(),
                self.hulls.len()
            )));
        }
        if self.hulls.len() > 1 && opts.centerplane.is_some() {
            return Err(Error::InvalidConditions(format!(
                "a centerplane override is ambiguous: {} separate hulls detected",
                self.hulls.len()
            )));
        }
        if opts.stations < 8 || opts.waterlines < 6 {
            return Err(Error::InvalidInput(
                "need at least 8 stations and 6 waterlines to sample".into(),
            ));
        }
        let mut members = Vec::new();
        let mut dry = Vec::new();
        for (hi, pose) in poses.iter().enumerate() {
            match self.situate_hull(hi, waterline_z, pose, platform, opts)? {
                Some(m) => members.push(m),
                None => dry.push(hi),
            }
        }
        Ok(SituatedFleet { members, dry })
    }

    /// Situate a single hull of the fleet; `Ok(None)` when it is dry.
    pub fn situate_one(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &ImportOptions,
    ) -> Result<Option<ImportedHull>> {
        if idx >= self.hulls.len() {
            return Err(Error::InvalidInput(format!(
                "hull index {idx} out of range ({} hulls)",
                self.hulls.len()
            )));
        }
        self.situate_hull(idx, waterline_z, pose, platform, opts)
    }

    /// Highest z (CAD frame, up) of a hull's control net — an upper bound on
    /// its geometry, used to loft full bands.
    pub fn hull_z_top(&self, idx: usize) -> f64 {
        self.hulls[idx]
            .iter()
            .flat_map(|s| s.ctrl.iter())
            .fold(f64::NEG_INFINITY, |m, p| m.max(p[2]))
    }

    /// Lowest z (CAD frame, up) of a hull's control net — a lower bound on
    /// its keel.
    pub fn hull_z_bottom(&self, idx: usize) -> f64 {
        self.hulls[idx]
            .iter()
            .flat_map(|s| s.ctrl.iter())
            .fold(f64::INFINITY, |m, p| m.min(p[2]))
    }

    fn situate_hull(
        &self,
        hi: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &ImportOptions,
    ) -> Result<Option<ImportedHull>> {
        let surfs = &self.hulls[hi];
        let wl = waterline_z + platform.sinkage;
        let mut moved = surfs.clone();
        // Design trim about (pivot_x, base waterline), then shifts.
        if pose.trim != 0.0 {
            let px = pose.pivot_x.unwrap_or_else(|| ctrl_x_mid(surfs));
            let (sin, cos) = pose.trim.sin_cos();
            for s in &mut moved {
                for p in s.ctrl.iter_mut() {
                    rotate_xz(p, px, waterline_z, cos, sin);
                }
            }
        }
        if pose.dx != 0.0 || pose.dy != 0.0 || pose.dz != 0.0 {
            for s in &mut moved {
                for p in s.ctrl.iter_mut() {
                    p[0] += pose.dx;
                    p[1] += pose.dy;
                    p[2] -= pose.dz;
                }
            }
        }
        // Platform pitch about (pivot_x, effective waterline).
        if platform.trim != 0.0 {
            let (sin, cos) = platform.trim.sin_cos();
            for s in &mut moved {
                for p in s.ctrl.iter_mut() {
                    rotate_xz(p, platform.pivot_x, wl, cos, sin);
                }
            }
        }
        let patches = presample_surfaces(&moved, wl);
        if patches.iter().all(|p| p.wet_box.is_none()) {
            return Ok(None);
        }
        let (hull, report, grid) = import_cluster(patches, opts, self.units_scale)?;
        Ok(Some(ImportedHull {
            placement: Placement {
                x: 0.0,
                y: report.centerplane,
            },
            hull,
            report,
            grid,
        }))
    }
}

/// Pitch rotation of a point in the x–z plane about `(px, pz)` (z up);
/// positive angle raises the +x side.
#[inline]
pub(crate) fn rotate_xz(p: &mut [f64; 3], px: f64, pz: f64, cos: f64, sin: f64) {
    let (dx, dz) = (p[0] - px, p[2] - pz);
    p[0] = px + dx * cos - dz * sin;
    p[2] = pz + dz * cos + dx * sin;
}

fn ctrl_x_mid(surfs: &[NurbsSurface3]) -> f64 {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for s in surfs {
        for p in &s.ctrl {
            lo = lo.min(p[0]);
            hi = hi.max(p[0]);
        }
    }
    0.5 * (lo + hi)
}

/// Parse + reject unsupported content (shared by every entry point).
fn validated_file(text: &str) -> Result<IgesFile> {
    let file = parse(text)?;
    if file.entity_counts.iter().any(|&(t, _)| t == 144 || t == 142) {
        return Err(Error::Unsupported(
            "the file contains trimmed surfaces (entities 142/144); export the \
             hull as untrimmed (or naturally bounded) surfaces"
                .into(),
        ));
    }
    if file.surfaces.is_empty() {
        let inventory: Vec<String> = file
            .entity_counts
            .iter()
            .map(|(t, c)| format!("{c} x type {t}"))
            .collect();
        return Err(Error::Unsupported(format!(
            "no B-spline surface (entity 128) found; file contains: {}",
            inventory.join(", ")
        )));
    }
    let rational = file.surfaces.iter().filter(|s| !s.is_polynomial()).count();
    if rational > 0 {
        return Err(Error::Unsupported(format!(
            "{rational} of {} surfaces are rational (non-uniform NURBS weights); \
             this crate's geometry contract is polynomial B-splines — re-export \
             with unit weights or refit",
            file.surfaces.len()
        )));
    }
    Ok(file)
}

/// Presample CAD-frame surfaces into hull-frame patches
/// (z' = waterline_z - z, downward).
fn presample_surfaces(surfaces: &[NurbsSurface3], waterline_z: f64) -> Vec<Patch> {
    const GRID_N: usize = 21;
    let mut patches = Vec::with_capacity(surfaces.len());
    for s in surfaces {
        let mut hs = s.clone();
        for p in hs.ctrl.iter_mut() {
            p[2] = waterline_z - p[2];
        }
        let (u0, u1) = hs.u_domain();
        let (v0, v1) = hs.v_domain();
        let mut pts = Vec::with_capacity(GRID_N * GRID_N);
        let mut bbox = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        let mut wet_box: Option<[f64; 6]> = None;
        for i in 0..GRID_N {
            let u = u0 + (u1 - u0) * i as f64 / (GRID_N - 1) as f64;
            for j in 0..GRID_N {
                let v = v0 + (v1 - v0) * j as f64 / (GRID_N - 1) as f64;
                let (p3, _, _) = hs.eval1(u, v);
                let (x, y, zd) = (p3[0], p3[1], p3[2]);
                bbox.0 = bbox.0.min(x);
                bbox.1 = bbox.1.max(x);
                bbox.2 = bbox.2.min(zd);
                bbox.3 = bbox.3.max(zd);
                if zd >= -1e-12 {
                    let b = wet_box.get_or_insert([x, x, y, y, zd, zd]);
                    b[0] = b[0].min(x);
                    b[1] = b[1].max(x);
                    b[2] = b[2].min(y);
                    b[3] = b[3].max(y);
                    b[4] = b[4].min(zd);
                    b[5] = b[5].max(zd);
                }
                pts.push((u, v, x, y, zd));
            }
        }
        patches.push(Patch {
            surf: hs,
            pts,
            grid_n: GRID_N,
            bbox,
            wet_box,
        });
    }
    patches
}

/// Import every hull found in an IGES file at the fixed waterline in `opts`.
///
/// Patches are clustered by wetted-geometry proximity: seams between patches
/// of one hull touch to CAD tolerance, while distinct hulls of a multihull
/// are far apart, so a whole boat modelled in position imports as a fleet
/// with its true placements. Only wetted geometry clusters, so dry structure
/// (cross-beams, decks) cannot bridge two hulls and is dropped. Results are
/// sorted by transverse position.
///
/// [`ImportOptions::centerplane`] may only be set when the file contains a
/// single hull. For repeated re-situating (sweeps over waterline, immersion,
/// trim), use [`source_fleet`] + [`SourceFleet::situate`].
pub fn import_fleet(text: &str, opts: &ImportOptions) -> Result<Vec<ImportedHull>> {
    let src = source_fleet(text, opts.waterline_z)?;
    let poses = vec![HullPose::default(); src.len()];
    let fl = src.situate(opts.waterline_z, &poses, &Platform::default(), opts)?;
    Ok(fl.members)
}

/// Import a hull from an IGES file that contains exactly one; errors (listing
/// the detected centerplanes) when the file holds several hulls — use
/// [`import_fleet`] for whole-multihull files.
pub fn import_hull(text: &str, opts: &ImportOptions) -> Result<(Hull, ImportReport)> {
    let mut fleet = import_fleet(text, opts)?;
    if fleet.len() == 1 {
        let m = fleet.pop().expect("one member");
        return Ok((m.hull, m.report));
    }
    let ys: Vec<String> = fleet
        .iter()
        .map(|m| format!("{:.3}", m.placement.y))
        .collect();
    Err(Error::Unsupported(format!(
        "the file contains {} separate hulls (centerplanes at y = {}); import \
         them as a fleet (import_fleet) or export hulls separately",
        fleet.len(),
        ys.join(", ")
    )))
}

/// Union-find clustering of patches whose wetted bounding boxes come within
/// `eps` of touching; dry patches are excluded.
#[allow(clippy::needless_range_loop)]
fn cluster_patches(patches: &[Patch], eps: f64) -> Vec<Vec<usize>> {
    let n = patches.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for i in 0..n {
        let Some(a) = patches[i].wet_box else { continue };
        for j in (i + 1)..n {
            let Some(b) = patches[j].wet_box else { continue };
            let touch = (0..3).all(|k| {
                a[2 * k] - eps <= b[2 * k + 1] && b[2 * k] - eps <= a[2 * k + 1]
            });
            if touch {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                parent[ri] = rj;
            }
        }
    }
    let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
    for i in 0..n {
        if patches[i].wet_box.is_none() {
            continue;
        }
        let r = find(&mut parent, i);
        match groups.iter_mut().find(|(root, _)| *root == r) {
            Some((_, v)) => v.push(i),
            None => groups.push((r, vec![i])),
        }
    }
    groups.into_iter().map(|(_, v)| v).collect()
}

/// Run the single-hull pipeline on one cluster of patches.
fn import_cluster(
    patches: Vec<Patch>,
    opts: &ImportOptions,
    units_scale: f64,
) -> Result<(Hull, ImportReport, SampleGrid)> {
    // Wetted statistics of this cluster.
    let mut draft = 0.0f64;
    let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut y_lo, mut y_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut y_sum = 0.0f64;
    let mut wet_count = 0usize;
    for p in &patches {
        for q in &p.pts {
            if q.4 >= -1e-12 {
                draft = draft.max(q.4);
                x_min = x_min.min(q.2);
                x_max = x_max.max(q.2);
                y_lo = y_lo.min(q.3);
                y_hi = y_hi.max(q.3);
                y_sum += q.3;
                wet_count += 1;
            }
        }
    }
    if wet_count == 0 || !(draft > 0.0 && x_max > x_min) {
        return Err(Error::InvalidGeometry(
            "a detected hull's wetted geometry is degenerate (zero draft or length)".into(),
        ));
    }
    let length = x_max - x_min;
    let scale = length.max(draft).max(y_hi - y_lo);

    // Centerplane: probe interior depths and count shell intersections.
    let mut probe_counts: Vec<usize> = Vec::new();
    let mut probe_mids: Vec<f64> = Vec::new();
    {
        let targets: Vec<(f64, f64)> = patches
            .iter()
            .flat_map(|p| p.pts.iter())
            .filter(|q| q.4 >= 0.3 * draft && q.4 <= 0.7 * draft)
            .map(|q| (q.2, q.4))
            .collect();
        let step = (targets.len() / 48).max(1);
        for t in targets.iter().step_by(step) {
            let ys = shell_intersections(&patches, t.0, t.1, scale);
            if !ys.is_empty() {
                probe_counts.push(ys.len());
                if ys.len() >= 2 {
                    probe_mids.push((ys[0].y + ys[ys.len() - 1].y) / 2.0);
                }
            }
        }
    }
    probe_counts.sort_unstable();
    let two_sided = !probe_counts.is_empty()
        && probe_counts[probe_counts.len() / 2] >= 2
        && !probe_mids.is_empty();
    let y_c = opts.centerplane.unwrap_or(if two_sided {
        probe_mids.iter().sum::<f64>() / probe_mids.len() as f64
    } else {
        0.0
    });
    let mirrored = !two_sided && y_sum / wet_count as f64 <= y_c;
    if !two_sided {
        // A half hull must actually reach its centerplane (keel/stem lines).
        let nearest = patches
            .iter()
            .flat_map(|p| p.pts.iter())
            .filter(|q| q.4 >= -1e-12)
            .fold(f64::INFINITY, |m, q| m.min((q.3 - y_c).abs()));
        if nearest > 0.2 * (y_hi - y_lo).max(1e-12) {
            return Err(Error::InvalidGeometry(format!(
                "the surface is one-sided but never approaches the centerplane \
                 y = {y_c}; if this is an offset hull, supply the centerplane \
                 position explicitly"
            )));
        }
    }

    // Sample grid: cosine-spaced stations, uniform waterlines to the draft.
    let ns = opts.stations;
    let nw = opts.waterlines;
    let stations: Vec<f64> = (0..ns)
        .map(|i| {
            let c = (std::f64::consts::PI * i as f64 / (ns - 1) as f64).cos();
            x_min + length * (1.0 - c) / 2.0
        })
        .collect();
    let waterlines: Vec<f64> = (0..nw).map(|j| draft * j as f64 / (nw - 1) as f64).collect();

    let mut grid = vec![0.0f64; ns * nw];
    let mut fx = vec![f64::NAN; ns * nw];
    let mut fz = vec![f64::NAN; ns * nw];
    let mut weight = vec![1.0f64; ns * nw];
    let mut failed = 0usize;
    let mut ambiguous = 0usize;
    let mut deriv_gaps = 0usize;
    let mut max_asym = 0.0f64;
    for (j, &zj) in waterlines.iter().enumerate() {
        let Some((flo, fhi)) = footprint(&patches, zj, draft) else {
            continue; // no hull at this depth
        };
        for (i, &xi) in stations.iter().enumerate() {
            let s = i * nw + j;
            if xi < flo || xi > fhi {
                continue;
            }
            let ys = shell_intersections(&patches, xi, zj, scale);
            if ys.is_empty() {
                // Unknown geometry, not zero beam: exclude from the loft.
                failed += 1;
                weight[s] = 0.0;
                continue;
            }
            if ys.len() > 2 {
                ambiguous += 1;
            }
            let outer = ys
                .iter()
                .max_by(|a, b| (a.y - y_c).abs().total_cmp(&(b.y - y_c).abs()))
                .expect("non-empty");
            let folded = (outer.y - y_c).abs();
            if two_sided && ys.len() >= 2 {
                let stb = (ys[ys.len() - 1].y - y_c).abs();
                let prt = (ys[0].y - y_c).abs();
                max_asym = max_asym.max((stb - prt).abs());
            }
            grid[s] = folded;
            // Half-beam is |y - y_c|: fold the slope's sign with y. At the
            // fold itself (keel/stem lines) the derivative is one-sided;
            // leave it unconstrained there.
            if outer.deriv_ok && folded > 1e-9 * scale {
                let sign = if outer.y >= y_c { 1.0 } else { -1.0 };
                fx[s] = sign * outer.y_x;
                fz[s] = sign * outer.y_z;
            } else {
                deriv_gaps += 1;
            }
        }
    }

    let sample_grid = SampleGrid::new(stations, waterlines, grid)?
        .with_fx(fx)?
        .with_fz(fz)?
        .with_weights(weight)?;
    let (hull, fit_report) = fit_grid(&sample_grid, &opts.fit)?;
    Ok((
        hull,
        ImportReport {
            units_scale,
            patches: patches.len(),
            two_sided,
            centerplane: y_c,
            mirrored,
            draft,
            x_range: (x_min, x_max),
            max_asymmetry: max_asym,
            ambiguous_samples: ambiguous,
            failed_inversions: failed,
            derivative_gaps: deriv_gaps,
            fit: fit_report,
        },
        sample_grid,
    ))
}

/// Longitudinal footprint [x_lo, x_hi] of the hull at depth `z`: union over
/// patches of level-set crossings of the presample grids, plus exact-plateau
/// points (so flat bottoms at maximum draft and sheer lines at z = 0 keep
/// their extent).
fn footprint(patches: &[Patch], z: f64, draft: f64) -> Option<(f64, f64)> {
    let tol = 1e-9 * draft.max(1e-300);
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut found = false;
    for p in patches {
        for q in &p.pts {
            if (q.4 - z).abs() <= tol {
                lo = lo.min(q.2);
                hi = hi.max(q.2);
                found = true;
            }
        }
        let n = p.grid_n;
        let mut cross = |a: &(f64, f64, f64, f64, f64), b: &(f64, f64, f64, f64, f64)| {
            if (a.4 - z) * (b.4 - z) < 0.0 {
                let t = (z - a.4) / (b.4 - a.4);
                let x = a.2 + t * (b.2 - a.2);
                lo = lo.min(x);
                hi = hi.max(x);
                found = true;
            }
        };
        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                if j + 1 < n {
                    cross(&p.pts[idx], &p.pts[idx + 1]);
                }
                if i + 1 < n {
                    cross(&p.pts[idx], &p.pts[idx + n]);
                }
            }
        }
    }
    found.then_some((lo, hi))
}

/// One shell intersection along y at a target `(x, z)`.
#[derive(Debug, Clone, Copy)]
struct ShellHit {
    y: f64,
    /// Surface slopes `∂y/∂x`, `∂y/∂z` (hull frame, z down) at the
    /// intersection; meaningful only when `deriv_ok`.
    y_x: f64,
    y_z: f64,
    deriv_ok: bool,
}

/// Slopes beyond this magnitude carry no loftable information at realistic
/// sample densities (the shell is effectively tangent-vertical there); the
/// sample keeps its value but drops its derivative.
const MAX_USEFUL_SLOPE: f64 = 1e2;

/// Build a hit at a converged parameter point, recovering `∂y/∂x`, `∂y/∂z`
/// from the surface partials by the implicit function theorem.
fn hit_at(surf: &NurbsSurface3, u: f64, v: f64) -> ShellHit {
    let (s, du, dv) = surf.eval1(u, v);
    // (x, z)(u, v) Jacobian determinant — the same one Newton inverts.
    let det = du[0] * dv[2] - dv[0] * du[2];
    let det_scale = du[0].abs() * dv[2].abs() + dv[0].abs() * du[2].abs();
    if det.abs() <= 1e-9 * det_scale.max(1e-300) {
        return ShellHit {
            y: s[1],
            y_x: 0.0,
            y_z: 0.0,
            deriv_ok: false,
        };
    }
    let y_x = (du[1] * dv[2] - dv[1] * du[2]) / det;
    let y_z = (dv[1] * du[0] - du[1] * dv[0]) / det;
    let ok = y_x.is_finite()
        && y_z.is_finite()
        && y_x.abs() < MAX_USEFUL_SLOPE
        && y_z.abs() < MAX_USEFUL_SLOPE;
    ShellHit {
        y: s[1],
        y_x,
        y_z,
        deriv_ok: ok,
    }
}

/// All shell intersections along y at `(x, z)`: per patch, Newton from the
/// nearest presample seeds; converged hits deduped across patches and
/// returned sorted by y.
fn shell_intersections(patches: &[Patch], x_t: f64, z_t: f64, scale: f64) -> Vec<ShellHit> {
    let tol = 1e-11 * scale;
    let loose = 1e-7 * scale;
    let margin = 0.05 * scale;
    let mut ys: Vec<ShellHit> = Vec::new();
    let mut best_loose: Option<(f64, f64)> = None; // (residual, y)
    for p in patches {
        if x_t < p.bbox.0 - margin
            || x_t > p.bbox.1 + margin
            || z_t < p.bbox.2 - margin
            || z_t > p.bbox.3 + margin
        {
            continue;
        }
        // Three nearest seeds within this patch.
        let mut seeds = [(f64::INFINITY, 0.0f64, 0.0f64); 3];
        for q in &p.pts {
            let d = (q.2 - x_t) * (q.2 - x_t) + (q.4 - z_t) * (q.4 - z_t);
            if d < seeds[2].0 {
                seeds[2] = (d, q.0, q.1);
                seeds.sort_by(|a, b| a.0.total_cmp(&b.0));
            }
        }
        for &(_, su, sv) in &seeds {
            let (res, y, u, v) = newton_on(&p.surf, su, sv, x_t, z_t, tol);
            if res < tol {
                ys.push(hit_at(&p.surf, u, v));
                break;
            }
            if best_loose.is_none_or(|(r, _)| res < r) {
                best_loose = Some((res, y));
            }
        }
    }
    if ys.is_empty() {
        // Accept a boundary-clamped near-miss (footprint edge) as a single
        // intersection; its parameters don't hit the target, so no slope.
        if let Some((r, y)) = best_loose {
            if r < loose {
                ys.push(ShellHit {
                    y,
                    y_x: 0.0,
                    y_z: 0.0,
                    deriv_ok: false,
                });
            }
        }
    }
    ys.sort_by(|a, b| a.y.total_cmp(&b.y));
    ys.dedup_by(|a, b| (a.y - b.y).abs() <= 1e-6 * scale);
    ys
}

/// 2-D Newton on `(x(u,v), z(u,v)) = (x_t, z_t)` from one seed; returns the
/// best residual reached and the y value and parameters there. Stops early
/// below `tol`.
fn newton_on(
    surf: &NurbsSurface3,
    mut u: f64,
    mut v: f64,
    x_t: f64,
    z_t: f64,
    tol: f64,
) -> (f64, f64, f64, f64) {
    let (u0, u1) = surf.u_domain();
    let (v0, v1) = surf.v_domain();
    let mut best = (f64::INFINITY, 0.0f64, u, v);
    for _ in 0..60 {
        let (s, du, dv) = surf.eval1(u, v);
        let (fx, fz) = (s[0] - x_t, s[2] - z_t);
        let res = fx.abs() + fz.abs();
        if res < best.0 {
            best = (res, s[1], u, v);
        }
        if res < tol {
            break;
        }
        let det = du[0] * dv[2] - dv[0] * du[2];
        if det.abs() < 1e-30 {
            break;
        }
        let step_u = (fx * dv[2] - fz * dv[0]) / det;
        let step_v = (-fx * du[2] + fz * du[0]) / det;
        u = (u - step_u).clamp(u0, u1);
        v = (v - step_v).clamp(v0, v1);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Nonlinearly parametrised test surface:
    /// x = L(0.7u + 0.3u³), z = T·v, y = (1 + u)(2 - v)/4.
    fn nonlinear_surface() -> NurbsSurface3 {
        let l = 12.0;
        let t = 2.0;
        // Cubic Bezier of x(u): power coeffs a0=0, a1=0.7L, a2=0, a3=0.3L.
        let bx = [0.0, 0.7 * l / 3.0, 2.0 * 0.7 * l / 3.0, l];
        // y is bilinear; represent at degree (3,1). A cubic Bezier of a linear
        // function has equally spaced control values.
        let yu = [1.0, 1.0 + 1.0 / 3.0, 1.0 + 2.0 / 3.0, 2.0];
        let yv = [2.0 / 4.0, 1.0 / 4.0];
        let zv = [0.0, t];
        let mut ctrl = Vec::new();
        let mut weights = Vec::new();
        for i in 0..4 {
            for j in 0..2 {
                ctrl.push([bx[i], yu[i] * yv[j], zv[j]]);
                weights.push(1.0);
            }
        }
        NurbsSurface3 {
            degree_u: 3,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            n_ctrl_u: 4,
            n_ctrl_v: 2,
            ctrl,
            weights,
        }
    }

    fn as_patch(surf: NurbsSurface3, n: usize) -> Patch {
        let (u0, u1) = surf.u_domain();
        let (v0, v1) = surf.v_domain();
        let mut pts = Vec::new();
        let mut bbox = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for i in 0..n {
            let u = u0 + (u1 - u0) * i as f64 / (n - 1) as f64;
            for j in 0..n {
                let v = v0 + (v1 - v0) * j as f64 / (n - 1) as f64;
                let (s, _, _) = surf.eval1(u, v);
                bbox.0 = bbox.0.min(s[0]);
                bbox.1 = bbox.1.max(s[0]);
                bbox.2 = bbox.2.min(s[2]);
                bbox.3 = bbox.3.max(s[2]);
                pts.push((u, v, s[0], s[1], s[2]));
            }
        }
        Patch {
            surf,
            pts,
            grid_n: n,
            bbox,
            wet_box: None,
        }
    }

    #[test]
    fn newton_inversion_recovers_y_under_nonlinear_parametrisation() {
        let surf = nonlinear_surface();
        let patches = vec![as_patch(surf.clone(), 33)];
        for &fu in &[0.03, 0.31, 0.5, 0.77, 0.99] {
            for &fv in &[0.02, 0.4, 0.96] {
                let (s, _, _) = surf.eval1(fu, fv);
                let ys = shell_intersections(&patches, s[0], s[2], 12.0);
                assert_eq!(ys.len(), 1, "u={fu} v={fv}: {ys:?}");
                assert!(
                    (ys[0].y - s[1]).abs() < 1e-9,
                    "u={fu} v={fv}: y {} vs {}",
                    ys[0].y,
                    s[1]
                );
            }
        }
    }

    #[test]
    fn inversion_recovers_surface_slopes() {
        // x = L(0.7u + 0.3u³), z = T·v, y = (1 + u)(2 − v)/4, so
        // ∂y/∂x = ((2 − v)/4) / (L(0.7 + 0.9u²)) and ∂y/∂z = −(1 + u)/(4T).
        let (l, t) = (12.0, 2.0);
        let surf = nonlinear_surface();
        let patches = vec![as_patch(surf.clone(), 33)];
        for &fu in &[0.03, 0.31, 0.5, 0.77, 0.95] {
            for &fv in &[0.05, 0.4, 0.9] {
                let (s, _, _) = surf.eval1(fu, fv);
                let ys = shell_intersections(&patches, s[0], s[2], 12.0);
                assert_eq!(ys.len(), 1, "u={fu} v={fv}: {ys:?}");
                let hit = ys[0];
                assert!(hit.deriv_ok, "u={fu} v={fv}");
                let want_yx = (2.0 - fv) / 4.0 / (l * (0.7 + 0.9 * fu * fu));
                let want_yz = -(1.0 + fu) / (4.0 * t);
                assert!(
                    (hit.y_x - want_yx).abs() < 1e-8,
                    "u={fu} v={fv}: y_x {} vs {want_yx}",
                    hit.y_x
                );
                assert!(
                    (hit.y_z - want_yz).abs() < 1e-8,
                    "u={fu} v={fv}: y_z {} vs {want_yz}",
                    hit.y_z
                );
            }
        }
    }

    #[test]
    fn two_patches_yield_both_intersections() {
        // Two parallel walls: y = +1 and y = -1, both spanning the same
        // (x, z) rectangle — a degenerate "full shell".
        let wall = |y: f64| -> NurbsSurface3 {
            NurbsSurface3 {
                degree_u: 1,
                degree_v: 1,
                knots_u: vec![0.0, 0.0, 10.0, 10.0],
                knots_v: vec![0.0, 0.0, 2.0, 2.0],
                n_ctrl_u: 2,
                n_ctrl_v: 2,
                ctrl: vec![[0.0, y, 0.0], [0.0, y, 2.0], [10.0, y, 0.0], [10.0, y, 2.0]],
                weights: vec![1.0; 4],
            }
        };
        let patches = vec![as_patch(wall(1.0), 9), as_patch(wall(-1.0), 9)];
        let ys = shell_intersections(&patches, 5.0, 1.0, 10.0);
        assert_eq!(ys.len(), 2, "{ys:?}");
        assert!((ys[0].y + 1.0).abs() < 1e-9 && (ys[1].y - 1.0).abs() < 1e-9);
        // A vertical wall has zero slope in both directions.
        assert!(ys.iter().all(|h| h.deriv_ok));
        assert!(ys.iter().all(|h| h.y_x.abs() < 1e-9 && h.y_z.abs() < 1e-9));
    }

    #[test]
    fn units_parsing() {
        // Global text with default delimiters, units flag field 14 = 2 (mm).
        let g = ",,7Hmichell,4Htest,7Hmichell,7Hmichell,32,38,6,308,15,4Htest,1.0,2,2HMM,1,0.01,15H20260719.000000,1E-08,10000.0,3Havi,7Hmichell,11,0;";
        let fields = split_global(g);
        assert_eq!(fields[13].trim(), "2");
        assert_eq!(fields[14].trim(), "MM");
    }

    #[test]
    fn number_forms() {
        assert_eq!(parse_number(" 1.5 ").unwrap(), 1.5);
        assert_eq!(parse_number("1.0D2").unwrap(), 100.0);
        assert_eq!(parse_number("").unwrap(), 0.0);
        assert!(parse_number("abc").is_err());
    }
}

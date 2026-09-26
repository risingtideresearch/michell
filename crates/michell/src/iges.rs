//! Minimal, dependency-free IGES reader and hull importer.
//!
//! Scope: **untrimmed rational B-spline surfaces** (entity 128) — one or many
//! patches — with optional transformation matrices (entity 124) and unit
//! conversion from the global section. Bounded-surface wrappers (143/141), as
//! produced by SubD/T-spline → NURBS exports, are supported approximately:
//! each base surface is restricted to the parameter-space **bounding box** of
//! its boundary curves, which is exact for rectangular boundaries (the common
//! natural-patch case) and conservative otherwise. Genuinely *trimmed*
//! surfaces (entities 142/144) are rejected. This is a hull-surface importer,
//! not a CAD kernel.
//!
//! ## Import pipeline
//!
//! A CAD surface is parametric, `S(u,v) → (x, y, z)`. The importer
//!
//! 1. parses every 128 entity (all weights must be uniform — the polynomial
//!    contract; see crate docs),
//! 2. clusters patches into hulls by wetted-geometry proximity
//!    ([`source_fleet`]), keeping each hull in the CAD frame so it can be
//!    re-posed ([`HullPose`], [`Platform`]) and re-cut repeatedly,
//! 3. finds each hull's **centerplane**: given via
//!    [`SectionalOptions::centerplane`], or auto-detected — interior probes
//!    count shell intersections along y; a two-sided (full) shell folds about
//!    the midplane of its intersections, a one-sided (half-hull) file measures
//!    from y = 0,
//! 4. cuts the posed hull at stations below the waterline and integrates
//!    each section along a fan of rays, taking the outermost fold, into a
//!    [`SectionalHull`] ([`SourceFleet::situate_sectional`]).
//!
//! Expected CAD frame: `x` longitudinal, `z` **up**, `y` transverse.
//! `waterline_z` gives the design waterline height in the file's frame, **in
//! metres** (after unit conversion). Diagnostics (ambiguous rays, the
//! detected transom) are reported so a bad import is visible.

use crate::bspline::{ders_basis, find_span, BSplineSurface};
use crate::error::{Error, Result};
use crate::michell::Placement;
use crate::sectional::{DepthQuadrature, SectionNodes, SectionalHull};

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
    /// `[u0, u1, v0, v1]` bounding box of this surface's boundary in
    /// parameter space, when the file marks it as the base of a bounded
    /// surface (entity 143). Sampling is restricted to this box: the region
    /// outside the boundary is construction geometry, not shell.
    pub trim_uv: Option<[f64; 4]>,
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
        let full = (self.knots_u[self.degree_u], self.knots_u[self.n_ctrl_u]);
        clip_domain(full, self.trim_uv.map(|t| (t[0], t[1])))
    }

    pub fn v_domain(&self) -> (f64, f64) {
        let full = (self.knots_v[self.degree_v], self.knots_v[self.n_ctrl_v]);
        clip_domain(full, self.trim_uv.map(|t| (t[2], t[3])))
    }

    /// Surface point at parameter `(u, v)` (see [`NurbsSurface3::u_domain`] /
    /// [`NurbsSurface3::v_domain`] for the valid range). Assumes uniform
    /// weights, like the rest of the sampler.
    pub fn point(&self, u: f64, v: f64) -> [f64; 3] {
        self.eval1(u, v).0
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

/// Intersect a knot domain with an optional trim interval; a trim that
/// leaves no proper interval is ignored rather than producing an empty or
/// inverted domain.
fn clip_domain(full: (f64, f64), trim: Option<(f64, f64)>) -> (f64, f64) {
    let Some((t0, t1)) = trim else { return full };
    let (a, b) = (full.0.max(t0), full.1.min(t1));
    if a < b && a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        full
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

    // Bounded surfaces (143): the untrimmed base surface can extend far past
    // the shell it carries (e.g. a full-width plane trimmed down to a keel
    // plank). Restrict each base surface to the parameter-space bounding box
    // of its boundary (141) so the phantom untrimmed region is never
    // sampled — left in, it can bridge the hulls of a multihull into one
    // cluster. Boundaries that cannot be resolved leave the surface
    // unrestricted (the pre-existing behaviour).
    let de_dir = |de: i64| -> Option<&Dir> {
        let de = de.unsigned_abs() as usize;
        if de == 0 || de.is_multiple_of(2) {
            return None;
        }
        dirs.get((de - 1) / 2)
    };
    // Bounding box over the control points of a parameter-space curve
    // (x, y) = (u, v); entities 126, 110, and 102 composites thereof.
    let curve_uv_bbox = |de0: i64| -> Option<[f64; 4]> {
        let mut stack = vec![de0];
        let mut b = [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ];
        let mut grow = |u: f64, v: f64| {
            b[0] = b[0].min(u);
            b[1] = b[1].max(u);
            b[2] = b[2].min(v);
            b[3] = b[3].max(v);
        };
        while let Some(de) = stack.pop() {
            let d = de_dir(de)?;
            let p = entity_params(d.pd_ptr, d.pd_count).ok()?;
            match d.etype {
                126 => {
                    let k = *p.get(1)? as usize;
                    let m = *p.get(2)? as usize;
                    let start = 7 + (k + m + 2) + (k + 1);
                    for i in 0..=k {
                        grow(*p.get(start + 3 * i)?, *p.get(start + 3 * i + 1)?);
                    }
                }
                110 => {
                    grow(*p.get(1)?, *p.get(2)?);
                    grow(*p.get(4)?, *p.get(5)?);
                }
                102 => {
                    let n = *p.get(1)? as usize;
                    for i in 0..n {
                        stack.push(*p.get(2 + i)? as i64);
                    }
                }
                _ => return None,
            }
        }
        // Iso-parameter edges are degenerate in one direction; only require
        // that some points were seen. Degenerate *unions* are discarded when
        // the trim is applied (`clip_domain` ignores improper intervals).
        (b[0] <= b[1] && b[2] <= b[3]).then_some(b)
    };
    // Boundary (141) -> parameter-space bbox, when it carries pcurves.
    let boundary_uv_bbox = |de: i64| -> Option<[f64; 4]> {
        let d = de_dir(de)?;
        if d.etype != 141 {
            return None;
        }
        let p = entity_params(d.pd_ptr, d.pd_count).ok()?;
        if *p.get(1)? as i64 != 1 {
            return None; // model-space representation only
        }
        let n = *p.get(4)? as usize;
        let mut b: Option<[f64; 4]> = None;
        let mut at = 5;
        for _ in 0..n {
            let k = *p.get(at + 2)? as usize;
            for j in 0..k {
                let c = curve_uv_bbox(*p.get(at + 3 + j)? as i64)?;
                let u = b.get_or_insert(c);
                u[0] = u[0].min(c[0]);
                u[1] = u[1].max(c[1]);
                u[2] = u[2].min(c[2]);
                u[3] = u[3].max(c[3]);
            }
            at += 3 + k;
        }
        b
    };
    // Base-surface DE -> union of its boundaries' bboxes; a 143 with any
    // unresolvable boundary imposes no restriction.
    let mut trims: Vec<(usize, [f64; 4])> = Vec::new();
    for d in dirs.iter().filter(|d| d.etype == 143) {
        let Ok(p) = entity_params(d.pd_ptr, d.pd_count) else {
            continue;
        };
        let (Some(&sptr), Some(&n)) = (p.get(2), p.get(3)) else {
            continue;
        };
        let boxes: Option<Vec<[f64; 4]>> = (0..n as usize)
            .map(|i| p.get(4 + i).and_then(|&b| boundary_uv_bbox(b as i64)))
            .collect();
        let Some(boxes) = boxes else { continue };
        let Some(joined) = boxes.into_iter().reduce(|mut a, c| {
            a[0] = a[0].min(c[0]);
            a[1] = a[1].max(c[1]);
            a[2] = a[2].min(c[2]);
            a[3] = a[3].max(c[3]);
            a
        }) else {
            continue;
        };
        trims.push((sptr.abs() as usize, joined));
    }

    let mut surfaces = Vec::new();
    for (idx, d) in dirs.iter().enumerate().filter(|(_, d)| d.etype == 128) {
        let p = entity_params(d.pd_ptr, d.pd_count)?;
        let mut surf = parse_surface_128(&p)?;
        surf.trim_uv = trims
            .iter()
            .filter(|(de, _)| *de == 2 * idx + 1)
            .map(|(_, b)| *b)
            .reduce(|mut a, c| {
                a[0] = a[0].min(c[0]);
                a[1] = a[1].max(c[1]);
                a[2] = a[2].min(c[2]);
                a[3] = a[3].max(c[3]);
                a
            });
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
                    p[r] =
                        m[4 * r] * q[0] + m[4 * r + 1] * q[1] + m[4 * r + 2] * q[2] + m[4 * r + 3];
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
        if deg > crate::bspline::MAX_DEGREE {
            return Err(Error::Parse(format!(
                "128 entity has degree {deg} in {dir}; the supported maximum is {}",
                crate::bspline::MAX_DEGREE
            )));
        }
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
        trim_uv: None,
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
        1 => 0.0254, // inches
        2 => 0.001,  // millimetres
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
        4 => 0.3048,   // feet
        5 => 1609.344, // miles
        6 => 1.0,      // metres
        7 => 1000.0,   // kilometres
        8 => 2.54e-5,  // mils
        9 => 1e-6,     // microns
        10 => 0.01,    // centimetres
        11 => 2.54e-8, // microinches
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

/// One patch, presampled in the hull frame (x, y, z-downward).
struct Patch {
    surf: NurbsSurface3,
    /// (u, v, x, y, z-downward), grid_n × grid_n row-major (u-major).
    pts: Vec<(f64, f64, f64, f64, f64)>,
    grid_n: usize,
    /// (x_lo, x_hi, z_lo, z_hi) over the whole patch presample.
    bbox: (f64, f64, f64, f64),
    /// (y_lo, y_hi) over the whole patch presample.
    ybox: (f64, f64),
    /// (x_lo, x_hi, y_lo, y_hi, z_lo, z_hi) over the wetted presample only;
    /// `None` when the patch is entirely above the waterline.
    wet_box: Option<[f64; 6]>,
}

/// Per-hull **design** pose: how a hull is mounted relative to the platform.
/// Applied to the source geometry before the waterline clip, so all fields
/// change the wetted shape exactly (affine maps of the control nets).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HullPose {
    /// Longitudinal shift [m].
    pub dx: f64,
    /// Transverse shift [m].
    pub dy: f64,
    /// Immersion shift [m]; positive lowers the hull (deeper).
    pub dz: f64,
    /// Pitch rotation [rad]; positive raises the hull's +x end.
    pub trim: f64,
    /// Uniform geometric scale factor (default `1.0`). Applied before every
    /// other field, about the hull's design waterline and transverse centre
    /// and the `pivot_x` station, so it grows or shrinks the whole hull in
    /// place (length, beam, and draft all scale together). Must be positive.
    pub scale: f64,
    /// Pivot station for `trim` (default: the hull's x mid); the pivot height
    /// is the base waterline.
    pub pivot_x: Option<f64>,
}

impl Default for HullPose {
    fn default() -> Self {
        HullPose {
            dx: 0.0,
            dy: 0.0,
            dz: 0.0,
            trim: 0.0,
            scale: 1.0,
            pivot_x: None,
        }
    }
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

/// Parsed and clustered source geometry, kept in CAD coordinates so hulls can
/// be re-situated (waterline, immersion, trim, position) repeatedly.
pub struct SourceFleet {
    units_scale: f64,
    /// Patches of each detected hull, CAD frame, metres, sorted by y.
    hulls: Vec<Vec<NurbsSurface3>>,
    /// Each hull's patches tessellated once, for fast sectioning at any pose.
    meshes: Vec<Tessellation>,
}

/// Parse an IGES file and cluster its patches into hulls at a reference
/// waterline (use the deepest waterline you intend to sweep, so cluster
/// membership stays fixed).
pub fn source_fleet(text: &str, reference_waterline: f64) -> Result<SourceFleet> {
    let file = validated_file(text)?;
    source_fleet_from_surfaces(file.surfaces, file.units_scale, reference_waterline)
}

/// [`source_fleet`] from surfaces already in hand (CAD frame, z up, metres),
/// e.g. [`wigley_surfaces`].
pub fn source_fleet_from_surfaces(
    surfaces: Vec<NurbsSurface3>,
    units_scale: f64,
    reference_waterline: f64,
) -> Result<SourceFleet> {
    let patches = presample_surfaces(&surfaces, reference_waterline);
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
                    let gap = (a[2 * k] - b[2 * k + 1])
                        .max(b[2 * k] - a[2 * k + 1])
                        .max(0.0);
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

    let hulls: Vec<Vec<NurbsSurface3>> = clusters
        .iter()
        .map(|idxs| idxs.iter().map(|&i| surfaces[i].clone()).collect())
        .collect();
    Ok(SourceFleet {
        units_scale,
        meshes: hulls.iter().map(|h| Tessellation::new(h)).collect(),
        hulls,
    })
}

impl SourceFleet {
    /// A B-spline half-breadth surface `y = f(x, z')` as source geometry:
    /// its exact mirrored pair of surfaces about `centerplane`, with the
    /// spline's top (z' = 0) at CAD height `top_z`: 0 for a wetted surface,
    /// or the freeboard for one carried up above the design waterline —
    /// either way the design waterline lands at CAD height 0, and the hull
    /// re-poses like any CAD hull. Test support.
    #[cfg(test)]
    pub(crate) fn from_halfbreadth(
        surface: &BSplineSurface,
        centerplane: f64,
        top_z: f64,
    ) -> Result<SourceFleet> {
        let surfaces = halfbreadth_surfaces(surface, centerplane, top_z).to_vec();
        source_fleet_from_surfaces(surfaces, 1.0, 0.0)
    }

    /// The x mid of a hull's control net: the default pivot of its design
    /// trim ([`HullPose::pivot_x`]).
    pub fn x_mid(&self, idx: usize) -> f64 {
        ctrl_x_mid(&self.hulls[idx])
    }

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

    /// Highest z (CAD frame, up) of a hull's control net — an upper bound on
    /// its geometry.
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

    /// A hull's patches with a pose applied — the geometry
    /// [`SourceFleet::situate_sectional`] would cut, expressed in the CAD frame
    /// (z up, metres) with the water surface back at `waterline_z`: platform
    /// sinkage moves the *hulls* down rather than the waterline up, so the
    /// result drops into a CAD model whose waterplane is fixed. Intended for
    /// re-exporting a studied configuration (see [`write`]).
    pub fn posed_surfaces(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> Result<Vec<NurbsSurface3>> {
        if idx >= self.hulls.len() {
            return Err(Error::InvalidInput(format!(
                "hull index {idx} out of range ({} hulls)",
                self.hulls.len()
            )));
        }
        let mut surfs = self.hulls[idx].clone();
        apply_pose(&mut surfs, waterline_z, pose, platform);
        Ok(surfs)
    }
}

/// Apply a design pose and platform state to CAD-frame surfaces (z up,
/// metres), by moving their control nets, and re-express platform sinkage
/// as a geometry shift so the water surface stays at `waterline_z`. This is
/// the configuration [`SourceFleet::situate_sectional`] evaluates (it instead
/// raises the waterline to `waterline_z + sinkage`, which is equivalent),
/// in a frame that drops into a CAD model whose waterplane is fixed.
pub fn apply_pose(
    surfs: &mut [NurbsSurface3],
    waterline_z: f64,
    pose: &HullPose,
    platform: &Platform,
) {
    pose_ctrl(surfs, waterline_z, pose, platform);
    if platform.sinkage != 0.0 {
        for s in surfs.iter_mut() {
            for p in s.ctrl.iter_mut() {
                p[2] -= platform.sinkage;
            }
        }
    }
}

/// The transform [`SourceFleet::situate_sectional`] applies before cutting at the
/// effective waterline `waterline_z + sinkage`: uniform `scale` about
/// `(pose.pivot_x, y-mid, waterline_z)`, then design trim about
/// `(pose.pivot_x, waterline_z)`, then the `dx`/`dy`/`dz` shifts
/// (`dz` positive lowers the hull), then platform pitch about
/// `(platform.pivot_x, waterline_z + sinkage)`.
fn pose_ctrl(surfs: &mut [NurbsSurface3], waterline_z: f64, pose: &HullPose, platform: &Platform) {
    let map = PoseMap::new(surfs, waterline_z, pose, platform);
    for p in surfs.iter_mut().flat_map(|s| s.ctrl.iter_mut()) {
        map.apply(p);
    }
}

/// [`pose_ctrl`]'s transform as a map on points, with its pivots taken from
/// the unposed control nets — so control points and points on the surface
/// (a tessellation's vertices) move by exactly the same affine map, under
/// which a B-spline surface and its points move together.
pub(crate) struct PoseMap {
    waterline_z: f64,
    wl: f64,
    scale: Option<(f64, f64, f64)>,
    trim: Option<(f64, f64, f64)>,
    shift: [f64; 3],
    platform_trim: Option<(f64, f64, f64)>,
}

impl PoseMap {
    fn new(
        surfs: &[NurbsSurface3],
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> Self {
        Self::with_mids(
            waterline_z,
            pose,
            platform,
            || ctrl_x_mid(surfs),
            || ctrl_y_mid(surfs),
        )
    }

    /// The same map for any geometry, given its longitudinal and transverse
    /// mids (the default trim and scale pivots), evaluated only if needed.
    pub(crate) fn with_mids(
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        x_mid: impl Fn() -> f64,
        y_mid: impl Fn() -> f64,
    ) -> Self {
        let px = || pose.pivot_x.unwrap_or_else(&x_mid);
        PoseMap {
            waterline_z,
            wl: waterline_z + platform.sinkage,
            scale: (pose.scale != 1.0).then(|| (pose.scale, px(), y_mid())),
            trim: (pose.trim != 0.0).then(|| {
                let (sin, cos) = pose.trim.sin_cos();
                (px(), cos, sin)
            }),
            shift: [pose.dx, pose.dy, -pose.dz],
            platform_trim: (platform.trim != 0.0).then(|| {
                let (sin, cos) = platform.trim.sin_cos();
                (platform.pivot_x, cos, sin)
            }),
        }
    }

    pub(crate) fn apply(&self, p: &mut [f64; 3]) {
        if let Some((s, px, py)) = self.scale {
            p[0] = px + s * (p[0] - px);
            p[1] = py + s * (p[1] - py);
            p[2] = self.waterline_z + s * (p[2] - self.waterline_z);
        }
        if let Some((px, cos, sin)) = self.trim {
            rotate_xz(p, px, self.waterline_z, cos, sin);
        }
        for k in 0..3 {
            p[k] += self.shift[k];
        }
        if let Some((px, cos, sin)) = self.platform_trim {
            rotate_xz(p, px, self.wl, cos, sin);
        }
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

/// Transverse mid of a hull's control net — the pivot `scale` shrinks toward,
/// so a symmetric hull scales about its centreplane and stays in place.
fn ctrl_y_mid(surfs: &[NurbsSurface3]) -> f64 {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for s in surfs {
        for p in &s.ctrl {
            lo = lo.min(p[1]);
            hi = hi.max(p[1]);
        }
    }
    0.5 * (lo + hi)
}

/// Parse + reject unsupported content (shared by every entry point).
fn validated_file(text: &str) -> Result<IgesFile> {
    let file = parse(text)?;
    if file
        .entity_counts
        .iter()
        .any(|&(t, _)| t == 144 || t == 142)
    {
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
        let mut bbox = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
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
                pts.push((u, v, x, y, zd));
            }
        }
        // Wetted box: every presample cell with a wet corner, all four of its
        // corners included. Wet points alone understate it by up to a cell —
        // on a long patch that is decimetres of bow (enough, on one CAD file,
        // to split the stem off as a separate "hull" and truncate the rest).
        for i in 0..GRID_N - 1 {
            for j in 0..GRID_N - 1 {
                let corners = [
                    pts[i * GRID_N + j],
                    pts[i * GRID_N + j + 1],
                    pts[(i + 1) * GRID_N + j],
                    pts[(i + 1) * GRID_N + j + 1],
                ];
                if corners.iter().all(|q| q.4 < -1e-12) {
                    continue;
                }
                for q in corners {
                    let b = wet_box.get_or_insert([q.2, q.2, q.3, q.3, q.4, q.4]);
                    b[0] = b[0].min(q.2);
                    b[1] = b[1].max(q.2);
                    b[2] = b[2].min(q.3);
                    b[3] = b[3].max(q.3);
                    b[4] = b[4].min(q.4.max(0.0));
                    b[5] = b[5].max(q.4);
                }
            }
        }
        let ybox = pts
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), q| {
                (lo.min(q.3), hi.max(q.3))
            });
        patches.push(Patch {
            surf: hs,
            pts,
            grid_n: GRID_N,
            bbox,
            ybox,
            wet_box,
        });
    }
    patches
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
        let Some(a) = patches[i].wet_box else {
            continue;
        };
        for j in (i + 1)..n {
            let Some(b) = patches[j].wet_box else {
                continue;
            };
            let touch =
                (0..3).all(|k| a[2 * k] - eps <= b[2 * k + 1] && b[2 * k] - eps <= a[2 * k + 1]);
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

/// A hull cluster's wetted extent and centreplane.
struct ClusterFrame {
    draft: f64,
    x_min: f64,
    x_max: f64,
    /// Largest wetted extent, for tolerances.
    scale: f64,
    y_c: f64,
    two_sided: bool,
    /// One-sided shell on the port side of its centreplane.
    mirrored: bool,
}

/// What a previous pose's frame lets the next one skip: its centreplane
/// (unchanged by sinkage and trim) and its ends (to re-bracket, not search).
#[derive(Debug, Clone, Copy)]
struct FrameMemo {
    y_c: f64,
    two_sided: bool,
    mirrored: bool,
    x_lo: f64,
    x_hi: f64,
}

/// Wetted statistics of a cluster and its centreplane: given, or detected
/// by probing interior depths for shell intersections — a two-sided shell
/// folds about the midplane of its intersections, a one-sided one measures
/// from y = 0 (and must reach it). Warm-started from a nearby pose's frame
/// when given.
fn cluster_frame_with(
    patches: &[Patch],
    centerplane: Option<f64>,
    prev: Option<&FrameMemo>,
    mesh: Option<&PosedMesh>,
) -> Result<ClusterFrame> {
    // Wetted statistics of this cluster.
    let mut draft = 0.0f64;
    let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut y_lo, mut y_hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut y_sum = 0.0f64;
    let mut wet_count = 0usize;
    for q in frame_points(patches, mesh) {
        if q[2] >= -1e-12 {
            draft = draft.max(q[2]);
            x_min = x_min.min(q[0]);
            x_max = x_max.max(q[0]);
            y_lo = y_lo.min(q[1]);
            y_hi = y_hi.max(q[1]);
            y_sum += q[1];
            wet_count += 1;
        }
    }
    if wet_count == 0 || !(draft > 0.0 && x_max > x_min) {
        return Err(Error::InvalidGeometry(
            "a detected hull's wetted geometry is degenerate (zero draft or length)".into(),
        ));
    }
    let length = x_max - x_min;
    let scale = length.max(draft).max(y_hi - y_lo);

    if let Some(m) = prev {
        // Same centreplane; ends re-bracketed from where they were.
        let sides = frame_sides(m.two_sided, m.mirrored);
        let step = 1e-3 * length;
        let x_lo = hull_end(
            patches,
            m.x_lo + step,
            -1.0,
            m.y_c,
            sides,
            length,
            scale,
            step,
            mesh,
        );
        let x_hi = hull_end(
            patches,
            m.x_hi - step,
            1.0,
            m.y_c,
            sides,
            length,
            scale,
            step,
            mesh,
        );
        let (x_min, x_max) = if x_hi > x_lo {
            (x_lo, x_hi)
        } else {
            (x_min, x_max)
        };
        return Ok(ClusterFrame {
            draft,
            x_min,
            x_max,
            scale,
            y_c: m.y_c,
            two_sided: m.two_sided,
            mirrored: m.mirrored,
        });
    }
    // Centerplane: probe interior depths and count shell intersections.
    let mut probe_counts: Vec<usize> = Vec::new();
    let mut probe_mids: Vec<f64> = Vec::new();
    {
        let band = |z: f64| z >= 0.3 * draft && z <= 0.7 * draft;
        // Probe points on the shell: the presample's, or — for a bare mesh —
        // triangle centroids (a vertex would put the probe exactly on the
        // facets' edges).
        let targets: Vec<(f64, f64)> = match (patches.is_empty(), mesh) {
            (true, Some(m)) => m
                .tess
                .tris
                .iter()
                .map(|t| {
                    let c = t.map(|k| m.verts[k as usize]);
                    (
                        (c[0][0] + c[1][0] + c[2][0]) / 3.0,
                        (c[0][2] + c[1][2] + c[2][2]) / 3.0,
                    )
                })
                .filter(|&(_, z)| band(z))
                .collect(),
            _ => patches
                .iter()
                .flat_map(|p| p.pts.iter())
                .filter(|q| band(q.4))
                .map(|q| (q.2, q.4))
                .collect(),
        };
        let step = (targets.len() / 48).max(1);
        for t in targets.iter().step_by(step) {
            let ys: Vec<f64> = match (patches.is_empty(), mesh) {
                (true, Some(m)) => m.transverse(t.0, t.1, scale),
                _ => shell_intersections(patches, t.0, t.1, scale),
            };
            if !ys.is_empty() {
                probe_counts.push(ys.len());
                if ys.len() >= 2 {
                    probe_mids.push((ys[0] + ys[ys.len() - 1]) / 2.0);
                }
            }
        }
    }
    probe_counts.sort_unstable();
    let two_sided = !probe_counts.is_empty()
        && probe_counts[probe_counts.len() / 2] >= 2
        && !probe_mids.is_empty();
    let y_c = centerplane.unwrap_or(if two_sided {
        probe_mids.iter().sum::<f64>() / probe_mids.len() as f64
    } else {
        0.0
    });
    let mirrored = !two_sided && y_sum / wet_count as f64 <= y_c;
    if !two_sided {
        // A half hull must actually reach its centerplane (keel/stem lines).
        let nearest = frame_points(patches, mesh)
            .filter(|q| q[2] >= -1e-12)
            .fold(f64::INFINITY, |m, q| m.min((q[1] - y_c).abs()));
        if nearest > 0.2 * (y_hi - y_lo).max(1e-12) {
            return Err(Error::InvalidGeometry(format!(
                "the surface is one-sided but never approaches the centerplane \
                 y = {y_c}; if this is an offset hull, supply the centerplane \
                 position explicitly"
            )));
        }
    }

    // The presample only brackets the ends (by up to a cell, which on a long
    // patch is decimetres); find where closed sections actually stop.
    let sides = frame_sides(two_sided, mirrored);
    let step = 0.02 * length;
    let x_lo = hull_end(
        patches,
        x_min + step,
        -1.0,
        y_c,
        sides,
        length,
        scale,
        step,
        mesh,
    );
    let x_hi = hull_end(
        patches,
        x_max - step,
        1.0,
        y_c,
        sides,
        length,
        scale,
        step,
        mesh,
    );
    let (x_min, x_max) = if x_hi > x_lo {
        (x_lo, x_hi)
    } else {
        (x_min, x_max)
    };
    Ok(ClusterFrame {
        draft,
        x_min,
        x_max,
        scale,
        y_c,
        two_sided,
        mirrored,
    })
}

/// Points on a cluster's shell (hull frame) for its wetted statistics: the
/// patches' presamples, or a bare mesh's vertices when there are no patches
/// (an STL import, where the mesh is the geometry).
fn frame_points<'a>(
    patches: &'a [Patch],
    mesh: Option<&'a PosedMesh>,
) -> impl Iterator<Item = [f64; 3]> + 'a {
    let bare = mesh.filter(|_| patches.is_empty());
    patches
        .iter()
        .flat_map(|p| p.pts.iter().map(|q| [q.2, q.3, q.4]))
        .chain(bare.into_iter().flat_map(|m| m.verts.iter().copied()))
}

/// Which side(s) of the centreplane carry shell: both for a full shell,
/// else the one a half hull lies on.
fn frame_sides(two_sided: bool, mirrored: bool) -> &'static [f64] {
    if two_sided {
        &[1.0, -1.0]
    } else if mirrored {
        &[-1.0]
    } else {
        &[1.0]
    }
}

/// All shell intersections along y at `(x, z)`: per patch, Newton from the
/// nearest presample seeds; converged hits deduped across patches and
/// returned sorted by y.
fn shell_intersections(patches: &[Patch], x_t: f64, z_t: f64, scale: f64) -> Vec<f64> {
    let tol = 1e-11 * scale;
    let loose = 1e-7 * scale;
    let margin = 0.05 * scale;
    let mut ys: Vec<f64> = Vec::new();
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
            let (res, y) = newton_on(&p.surf, su, sv, x_t, z_t, tol);
            if res < tol {
                ys.push(y);
                break;
            }
            if best_loose.is_none_or(|(r, _)| res < r) {
                best_loose = Some((res, y));
            }
        }
    }
    if ys.is_empty() {
        // Accept a boundary-clamped near-miss (footprint edge) as a single
        // intersection.
        if let Some((r, y)) = best_loose {
            if r < loose {
                ys.push(y);
            }
        }
    }
    ys.sort_by(|a, b| a.total_cmp(b));
    ys.dedup_by(|a, b| (*a - *b).abs() <= 1e-6 * scale);
    ys
}

/// 2-D Newton on `(x(u,v), z(u,v)) = (x_t, z_t)` from one seed; returns the
/// best residual reached and the y value there. Stops early below `tol`.
fn newton_on(
    surf: &NurbsSurface3,
    mut u: f64,
    mut v: f64,
    x_t: f64,
    z_t: f64,
    tol: f64,
) -> (f64, f64) {
    let (u0, u1) = surf.u_domain();
    let (v0, v1) = surf.v_domain();
    let mut best = (f64::INFINITY, 0.0f64);
    for _ in 0..60 {
        let (s, du, dv) = surf.eval1(u, v);
        let (fx, fz) = (s[0] - x_t, s[2] - z_t);
        let res = fx.abs() + fz.abs();
        if res < best.0 {
            best = (res, s[1]);
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

// ---------------------------------------------------------------------------
// Sectional import
// ---------------------------------------------------------------------------

/// Options for a sectional import: the hull is cut at stations and each
/// section integrated along rays.
#[derive(Debug, Clone, Copy)]
pub struct SectionalOptions {
    /// Height of the design waterline in the file's frame (z up), in
    /// **metres** (i.e. after unit conversion).
    pub waterline_z: f64,
    /// Transverse position of the hull's centreplane in the file frame [m].
    /// `None` auto-detects: a two-sided (full) shell folds about the midplane
    /// of its shell intersections; a one-sided (half-hull) file measures from
    /// y = 0.
    pub centerplane: Option<f64>,
    /// Stations along the hull, cosine-spaced with both ends included.
    pub stations: usize,
    /// Rays sampled across each section (Chebyshev–Lobatto in angle); the
    /// depth quadrature evaluates their interpolant.
    pub rays: usize,
    pub quadrature: DepthQuadrature,
}

impl Default for SectionalOptions {
    fn default() -> Self {
        SectionalOptions {
            waterline_z: 0.0,
            centerplane: None,
            stations: 121,
            rays: 33,
            quadrature: DepthQuadrature::default(),
        }
    }
}

/// Diagnostics of a sectional import.
#[derive(Debug, Clone)]
pub struct SectionalReport {
    pub units_scale: f64,
    /// Surface patches the hull was cut from (triangles, for a mesh import).
    pub patches: usize,
    pub two_sided: bool,
    pub centerplane: f64,
    pub mirrored: bool,
    pub draft: f64,
    pub x_range: (f64, f64),
    /// Stations the hull was built from.
    pub stations: usize,
    /// Interior stations dropped because a ray found no shell.
    pub dropped_stations: usize,
    /// Rays that met the shell more than once (the section is not
    /// star-shaped about its top centreplane point there, or other shell
    /// lies beyond it; the nearest hit is used).
    pub ambiguous_rays: usize,
    /// Largest port/starboard disagreement of a ray's reach [m], for a
    /// two-sided shell (the sides are averaged: the symmetric thickness).
    pub max_asymmetry: f64,
    /// The transom the aft end presents, if immersed (see
    /// [`SectionalHull::transom`]).
    pub transom: Option<crate::sectional::Transom>,
}

/// One hull imported by sections.
#[derive(Debug, Clone)]
pub struct SectionalImport {
    pub hull: SectionalHull,
    pub placement: Placement,
    pub report: SectionalReport,
    /// Each station's x and its sampled outline, `(half-beam, depth)` from
    /// the waterline (or the section's top) round to the keel — for display.
    pub sections: Vec<(f64, Vec<(f64, f64)>)>,
}

/// Every hull of an IGES file imported by sections: those that could be, and
/// why the others could not (a sliver of appendage geometry, say), so one
/// bad fragment does not sink the import of the rest.
#[derive(Debug, Clone)]
pub struct SectionalFleet {
    pub hulls: Vec<SectionalImport>,
    /// `(source hull index, reason)` for each hull that failed.
    pub failed: Vec<(usize, String)>,
}

/// Import every hull in an IGES file by sections, at the fixed waterline in
/// `opts`.
pub fn import_sectional(text: &str, opts: &SectionalOptions) -> Result<SectionalFleet> {
    let fleet = source_fleet(text, opts.waterline_z)?;
    collect_sectional(fleet.len(), |i| {
        fleet.situate_sectional(
            i,
            opts.waterline_z,
            &HullPose::default(),
            &Platform::default(),
            opts,
        )
    })
}

/// Every hull of a fleet by sections, unposed: the ones that could be, the
/// reasons the others could not.
pub(crate) fn collect_sectional(
    n: usize,
    situate: impl Fn(usize) -> Result<Option<SectionalImport>>,
) -> Result<SectionalFleet> {
    let mut out = SectionalFleet {
        hulls: Vec::new(),
        failed: Vec::new(),
    };
    for i in 0..n {
        match situate(i) {
            Ok(Some(h)) => out.hulls.push(h),
            Ok(None) => {}
            Err(e) => out.failed.push((i, e.to_string())),
        }
    }
    if out.hulls.is_empty() {
        return Err(Error::InvalidGeometry(match out.failed.first() {
            Some((_, e)) => format!("no hull could be sectioned: {e}"),
            None => "no hull is wetted at this waterline".into(),
        }));
    }
    Ok(out)
}

/// What one sectional situate leaves for the next, so a sweep that re-poses
/// the same hull many times (sinkage and trim in an equilibrium solve)
/// skips the centreplane detection and re-brackets the hull's ends rather
/// than searching for them. Carry one per hull between calls to
/// [`SourceFleet::situate_sectional_warm`]; a fresh one is a cold start.
/// Only valid between poses of the same hull with the same options,
/// differing by sinkage and trim.
#[derive(Debug, Clone, Default)]
pub struct SectionalState {
    frame: Option<FrameMemo>,
}

impl SourceFleet {
    /// Situate one hull — apply its design pose and the platform state, clip
    /// at the effective waterline `waterline_z + sinkage` — and build it by
    /// sections; `Ok(None)` when it is dry.
    pub fn situate_sectional(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &SectionalOptions,
    ) -> Result<Option<SectionalImport>> {
        self.situate_sectional_warm(
            idx,
            waterline_z,
            pose,
            platform,
            opts,
            &mut SectionalState::default(),
        )
    }

    /// A hull's tessellation at a pose, in the water frame: `x` forward, `y`
    /// transverse, `z` **up** from the effective waterline — vertices and
    /// triangles, the whole hull (above water too), for display.
    pub fn posed_tessellation(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>)> {
        if idx >= self.hulls.len() {
            return Err(Error::InvalidInput(format!(
                "hull index {idx} out of range ({} hulls)",
                self.hulls.len()
            )));
        }
        let map = PoseMap::new(&self.hulls[idx], waterline_z, pose, platform);
        let wl = waterline_z + platform.sinkage;
        let tess = &self.meshes[idx];
        let verts = tess
            .verts
            .iter()
            .map(|&p| {
                let mut q = p;
                map.apply(&mut q);
                [q[0], q[1], q[2] - wl]
            })
            .collect();
        Ok((verts, tess.tris.clone()))
    }

    /// [`SourceFleet::situate_sectional`], warm-started from (and updating)
    /// `state` — the cheap way to re-pose.
    pub fn situate_sectional_warm(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &SectionalOptions,
        state: &mut SectionalState,
    ) -> Result<Option<SectionalImport>> {
        if idx >= self.hulls.len() {
            return Err(Error::InvalidInput(format!(
                "hull index {idx} out of range ({} hulls)",
                self.hulls.len()
            )));
        }
        if opts.stations < 8 || opts.rays < 5 {
            return Err(Error::InvalidInput(
                "need at least 8 stations and 5 rays per section".into(),
            ));
        }
        let wl = waterline_z + platform.sinkage;
        let mut moved = self.hulls[idx].clone();
        let map = PoseMap::new(&self.hulls[idx], waterline_z, pose, platform);
        for p in moved.iter_mut().flat_map(|s| s.ctrl.iter_mut()) {
            map.apply(p);
        }
        let patches = presample_surfaces(&moved, wl);
        if patches.iter().all(|p| p.wet_box.is_none()) {
            return Ok(None);
        }
        let mesh = PosedMesh::new(&self.meshes[idx], |mut p| {
            map.apply(&mut p);
            [p[0], p[1], wl - p[2]]
        });
        sectional_cluster(&patches, Some(&mesh), opts, self.units_scale, state).map(Some)
    }
}

/// Build one posed hull by sections from a bare mesh that is the geometry
/// itself (an STL): the same stations, rays, reach rule, end caps, ends and
/// frame as for patches, with the hits taken straight from the facets.
/// `Ok(None)` when the mesh is dry.
pub(crate) fn sectional_mesh(
    mesh: &PosedMesh,
    opts: &SectionalOptions,
    units_scale: f64,
    state: &mut SectionalState,
) -> Result<Option<SectionalImport>> {
    if mesh.verts.iter().all(|p| p[2] < -1e-12) {
        return Ok(None);
    }
    let mut imp = sectional_cluster(&[], Some(mesh), opts, units_scale, state)?;
    imp.report.patches = mesh.tess.tris.len();
    Ok(Some(imp))
}

/// Build one posed hull by sections. `patches` are its exact surfaces (the
/// mesh, when given, guides the search and is polished onto them); with no
/// patches, `mesh` is required and is the geometry itself.
fn sectional_cluster(
    patches: &[Patch],
    mesh: Option<&PosedMesh>,
    opts: &SectionalOptions,
    units_scale: f64,
    state: &mut SectionalState,
) -> Result<SectionalImport> {
    let frame = cluster_frame_with(patches, opts.centerplane, state.frame.as_ref(), mesh)?;
    let (scale, y_c) = (frame.scale, frame.y_c);
    // Which side(s) of the centreplane carry shell.
    let sides = frame_sides(frame.two_sided, frame.mirrored);
    let half_pi = std::f64::consts::FRAC_PI_2;
    let nr = opts.rays;
    let thetas: Vec<f64> = (0..nr)
        .map(|k| 0.5 * half_pi * (1.0 - (std::f64::consts::PI * k as f64 / (nr - 1) as f64).cos()))
        .collect();

    let ns = opts.stations;
    // End stations exactly where closed sections stop (see
    // [`cluster_frame`]): a station a hair inside a pointed end would carry a
    // small section the interpolant then drops to nothing — a spurious
    // microscopic transom, visible at low speed, where the amplitude is a
    // small residue of bow–stern cancellation.
    let (lo, hi) = (frame.x_min, frame.x_max);
    let xs: Vec<f64> = (0..ns)
        .map(|i| {
            let c = (std::f64::consts::PI * i as f64 / (ns - 1) as f64).cos();
            lo + (hi - lo) * (1.0 - c) / 2.0
        })
        .collect();

    let (mut ambiguous, mut dropped, mut max_asym) = (0usize, 0usize, 0.0f64);
    let mut kept_x = Vec::with_capacity(ns);
    let mut sections = Vec::with_capacity(ns);
    let mut outlines = Vec::with_capacity(ns);
    // Stations are independent: sample them across cores, then assemble in
    // order (the result does not depend on the thread count).
    let sampled_all = crate::parallel::map_indexed(
        ns,
        || (),
        |_, i| {
            let mut amb = 0;
            let s = sample_section(patches, mesh, xs[i], y_c, sides, &thetas, scale, &mut amb);
            (s, amb)
        },
    );
    for (i, (&x, (sampled, amb))) in xs.iter().zip(sampled_all).enumerate() {
        ambiguous += amb;
        let end = i == 0 || i + 1 == ns;
        match sampled {
            Some((z0, beam, depth, radii, asym)) => {
                max_asym = max_asym.max(asym);
                let interp = Lobatto::new(&thetas, &radii);
                sections.push(SectionNodes::from_polar(
                    z0,
                    beam,
                    depth,
                    |t| interp.eval(t),
                    &opts.quadrature,
                ));
                outlines.push((
                    x,
                    thetas
                        .iter()
                        .zip(&radii)
                        .map(|(&t, &r)| (beam * r * t.cos(), z0 + depth * r * t.sin()))
                        .collect(),
                ));
                kept_x.push(x);
            }
            None if end => {
                // Past the hull's tip: no section.
                sections.push(SectionNodes::empty());
                outlines.push((x, Vec::new()));
                kept_x.push(x);
            }
            None => dropped += 1,
        }
    }
    let n = kept_x.len();
    if n < 8 {
        return Err(Error::InvalidGeometry(format!(
            "only {n} of {ns} stations could be sectioned"
        )));
    }
    // Cubic interpolation along x, not-a-knot: interior knots at every
    // station but the second and second-to-last.
    let p = 3;
    let mut knots = vec![kept_x[0]; p + 1];
    knots.extend_from_slice(&kept_x[2..n - 2]);
    knots.extend(std::iter::repeat_n(kept_x[n - 1], p + 1));
    let hull = SectionalHull::new(p, knots, &kept_x, sections)?;
    let transom = hull.transom().cloned();
    state.frame = Some(FrameMemo {
        y_c,
        two_sided: frame.two_sided,
        mirrored: frame.mirrored,
        x_lo: lo,
        x_hi: hi,
    });
    Ok(SectionalImport {
        hull,
        placement: Placement { x: 0.0, y: y_c },
        report: SectionalReport {
            units_scale,
            patches: patches.len(),
            two_sided: frame.two_sided,
            centerplane: y_c,
            mirrored: frame.mirrored,
            draft: frame.draft,
            x_range: (lo, hi),
            stations: n,
            dropped_stations: dropped,
            ambiguous_rays: ambiguous,
            max_asymmetry: max_asym,
            transom,
        },
        sections: outlines,
    })
}

/// One station plane `x = const` through a cluster's patches. Where a posed
/// tessellation of the patches is available ([`PosedMesh`]), the plane's cut
/// through it gives every place a ray can meet the shell — with the patch and
/// parameters already known, so each hit is one or two Newton steps on the
/// exact surface away. Without one (or when a polish strays), hits are searched
/// for from presample seeds. Patch-edge crossings are found up front either
/// way, because a section's defining points often lie *on* an edge — a keel
/// or stem line where the half-breadth is exactly zero, a seam, an open rim —
/// where a Newton solve in the parameter domain has nowhere to go but the
/// boundary.
struct Station<'a> {
    patches: &'a [Patch],
    x: f64,
    y_c: f64,
    scale: f64,
    /// `(y, z)` (hull frame, z down) of every patch-edge crossing.
    edges: Vec<(f64, f64)>,
    /// The posed tessellation cut by the plane: segments between points on
    /// the shell, with where they lie on their patch. `None` without a mesh.
    cut: Option<Vec<[Cut; 2]>>,
    /// The mesh is the geometry (no patches to polish on): hits are the
    /// cut's own.
    exact: bool,
}

/// A point where the station plane cuts a tessellation edge: its position
/// (hull frame) and its patch and parameters (interpolated along the edge).
#[derive(Debug, Clone, Copy)]
struct Cut {
    y: f64,
    z: f64,
    patch: usize,
    u: f64,
    v: f64,
}

/// One place a ray meets the shell: its reach and, for a patch hit, the
/// patch and parameters.
type Hit = (f64, Option<(usize, f64, f64)>);

impl<'a> Station<'a> {
    fn new(patches: &'a [Patch], mesh: Option<&PosedMesh>, x: f64, y_c: f64, scale: f64) -> Self {
        let mut edges = Vec::new();
        for p in patches {
            if x < p.bbox.0 || x > p.bbox.1 {
                continue;
            }
            let g = p.grid_n;
            let runs: [Vec<usize>; 4] = [
                (0..g).collect(),
                (0..g).map(|j| (g - 1) * g + j).collect(),
                (0..g).map(|i| i * g).collect(),
                (0..g).map(|i| i * g + g - 1).collect(),
            ];
            for run in &runs {
                for w in run.windows(2) {
                    let (a, b) = (p.pts[w[0]], p.pts[w[1]]);
                    if (a.2 - x) * (b.2 - x) > 0.0 {
                        continue;
                    }
                    let at = |t: f64| p.surf.point(a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1));
                    let (mut ta, mut tb) = (0.0f64, 1.0f64);
                    let fa = at(ta)[0] - x;
                    for _ in 0..60 {
                        let tm = 0.5 * (ta + tb);
                        if (at(tm)[0] - x) * fa > 0.0 {
                            ta = tm;
                        } else {
                            tb = tm;
                        }
                    }
                    let t = 0.5 * (ta + tb);
                    let (u, v) = (a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1));
                    if is_end_cap(&p.surf, u, v) {
                        continue;
                    }
                    let q = p.surf.point(u, v);
                    edges.push((q[1], q[2]));
                }
            }
        }
        Station {
            patches,
            x,
            y_c,
            scale,
            edges,
            cut: mesh.map(|m| m.cut(x)),
            exact: mesh.is_some_and(|m| m.tess.params.is_empty()),
        }
    }

    /// Whether the plane cuts a **closed** wetted section: one its rays can
    /// sweep — reaching the shell along the waterline (or hanging from the
    /// section's top below it) and down the centreplane to a keel. Bare
    /// skins with nothing between them (side shells running on past a
    /// recessed transom, say) enclose no hull and do not count.
    fn has_section(&self, sides: &[f64]) -> bool {
        let mut amb = 0;
        let half_pi = std::f64::consts::FRAC_PI_2;
        let sweeps = |z0: f64, amb: &mut usize| {
            (z0 > 0.0
                || sides
                    .iter()
                    .all(|&sd| self.reach(sd, z0, 0.0, amb).is_some()))
                && self.reach(sides[0], z0, half_pi, amb).is_some()
        };
        if sweeps(0.0, &mut amb) {
            return true;
        }
        self.top()
            .is_some_and(|z0| z0 > 0.0 && sweeps(z0, &mut amb))
    }

    /// Distance from `(y_c, z0)` along the ray at angle `θ` below the
    /// horizontal, toward side `sd = ±1`, to the boundary of the hull's
    /// section — the region defined by the **outermost fold**: at each depth the section reaches as far from
    /// the centreplane as the farthest shell there. The reach is where the ray
    /// first *leaves* that region: at each hit in turn, the ray's point
    /// midway to the next hit is tested against the fold at its depth. Past
    /// internal geometry (a floor, a beam, a bulkhead) the ray is still inside
    /// the fold and goes on; past a boundary with water beyond it — a raked
    /// stem face with nothing but patch slivers below it, say — it is not,
    /// and stops, even where the boundary is no wider than the fold at its
    /// own depth. (Taking instead the farthest hit that lies on the fold, as
    /// if the region were star-shaped, carried rays through such a stem face
    /// to the slivers and hung a spike off every section under it.) More than
    /// one distinct hit counts as ambiguous.
    fn reach(&self, sd: f64, z0: f64, theta: f64, ambiguous: &mut usize) -> Option<f64> {
        let hits = self.hits(sd, z0, theta);
        if hits.len() > 1 {
            *ambiguous += 1;
        }
        if hits.len() <= 1 {
            return hits.last().map(|h| h.0);
        }
        let (st, ct) = theta.sin_cos();
        let tol = 1e-9 * self.scale;
        for w in hits.windows(2) {
            let r = 0.5 * (w[0].0 + w[1].0);
            let (yb, z) = (r * ct, z0 + r * st);
            // The fold's outermost reach at this depth: the ray is inside
            // the section there only if some shell lies farther out.
            let outer = self.hits(sd, z, 0.0).last().map_or(0.0, |h| h.0);
            if yb + tol >= outer {
                return Some(w[0].0);
            }
        }
        hits.last().map(|h| h.0)
    }

    /// Newton for the ray's hit on one patch from `(u, v)`: the reach and
    /// where it landed, if it converged ahead of the origin — `Some(None)`
    /// when it landed on an end cap (see [`is_end_cap`]), which is no hit.
    #[allow(clippy::option_option)]
    fn solve_on(
        &self,
        pi: usize,
        u: f64,
        v: f64,
        sd: f64,
        z0: f64,
        theta: f64,
    ) -> Option<Option<(f64, (usize, f64, f64))>> {
        let (st, ct) = theta.sin_cos();
        let (x_t, y_c, scale) = (self.x, self.y_c, self.scale);
        let off = |y: f64, z: f64| sd * (y - y_c) * st - (z - z0) * ct;
        let p = self.patches.get(pi)?;
        let (s, uv) = newton_in_domain_uv(&p.surf, u, v, 1e-11 * scale, |s, du, dv| {
            (
                [s[0] - x_t, off(s[1], s[2])],
                [
                    [du[0], dv[0]],
                    [sd * du[1] * st - du[2] * ct, sd * dv[1] * st - dv[2] * ct],
                ],
            )
        })?;
        if is_end_cap(&p.surf, uv.0, uv.1) {
            return Some(None);
        }
        let r = sd * (s[1] - y_c) * ct + (s[2] - z0) * st;
        (r > 1e-12 * scale).then_some(Some((r, (pi, uv.0, uv.1))))
    }

    /// Every distinct distance along the ray from `(y_c, z0)` at angle `θ`
    /// toward side `sd` at which it meets the shell, nearest first,
    /// including the edge crossings that lie on the ray. From the mesh cut
    /// when there is one: each cut segment the ray crosses, polished on the
    /// exact patch from its interpolated parameters (a polish that fails, or
    /// wanders off to another hit, falls back to the seed search for the
    /// whole ray).
    fn hits(&self, sd: f64, z0: f64, theta: f64) -> Vec<Hit> {
        let Some(cut) = &self.cut else {
            return self.search_hits(sd, z0, theta);
        };
        if self.exact {
            return self.mesh_hits(cut, sd, z0, theta);
        }
        let (st, ct) = theta.sin_cos();
        let (dy, dz) = (sd * ct, st);
        let mut hits = self.edge_hits(sd, z0, theta);
        for [a, b] in cut {
            // Ray (y_c, z0) + r·(dy, dz) against segment a + s·(b − a).
            let (ey, ez) = (b.y - a.y, b.z - a.z);
            let den = dy * ez - dz * ey;
            if den.abs() < 1e-300 {
                continue;
            }
            let (qy, qz) = (a.y - self.y_c, a.z - z0);
            let r = (qy * ez - qz * ey) / den;
            let s = (qy * dz - qz * dy) / den;
            if !(-1e-9..=1.0 + 1e-9).contains(&s) || r <= 0.0 {
                continue;
            }
            let s = s.clamp(0.0, 1.0);
            let (u, v) = (a.u + s * (b.u - a.u), a.v + s * (b.v - a.v));
            match self.solve_on(a.patch, u, v, sd, z0, theta) {
                // A polish lands within a facet's chord of where the mesh
                // said; one that went much farther found some other hit.
                Some(Some((rr, at))) if (rr - r).abs() <= 0.02 * self.scale => {
                    hits.push((rr, Some(at)))
                }
                // Onto an end cap: not a hit.
                Some(None) => {}
                _ => return self.search_hits(sd, z0, theta),
            }
        }
        self.tidy(hits)
    }

    /// [`Station::hits`] on a mesh that is itself the geometry (an STL): the
    /// ray's crossings of the cut segments, as they are.
    fn mesh_hits(&self, cut: &[[Cut; 2]], sd: f64, z0: f64, theta: f64) -> Vec<Hit> {
        let (st, ct) = theta.sin_cos();
        let (dy, dz) = (sd * ct, st);
        let min_reach = 1e-12 * self.scale;
        let mut hits = Vec::new();
        for [a, b] in cut {
            let (ey, ez) = (b.y - a.y, b.z - a.z);
            let den = dy * ez - dz * ey;
            if den.abs() < 1e-300 {
                continue;
            }
            let (qy, qz) = (a.y - self.y_c, a.z - z0);
            let r = (qy * ez - qz * ey) / den;
            let s = (qy * dz - qz * dy) / den;
            if (-1e-9..=1.0 + 1e-9).contains(&s) && r > min_reach {
                hits.push((r, None));
            }
        }
        self.tidy(hits)
    }

    /// [`Station::hits`] by search: patch hits by Newton on `(u, v)` for
    /// `x(u,v) = x` and the point lying on the ray, from the nearest
    /// presample seeds ahead of the origin.
    fn search_hits(&self, sd: f64, z0: f64, theta: f64) -> Vec<Hit> {
        let (st, ct) = theta.sin_cos();
        let (x_t, y_c, scale) = (self.x, self.y_c, self.scale);
        let tol = 1e-11 * scale;
        let margin = 0.02 * scale;
        let min_reach = 1e-12 * scale;
        let off = |y: f64, z: f64| sd * (y - y_c) * st - (z - z0) * ct;
        let along = |y: f64, z: f64| sd * (y - y_c) * ct + (z - z0) * st;
        let mut hits = self.edge_hits(sd, z0, theta);
        for (pi, p) in self.patches.iter().enumerate() {
            if x_t < p.bbox.0 - margin || x_t > p.bbox.1 + margin {
                continue;
            }
            // Skip patches the ray passes well clear of: every corner of the
            // patch's (y, z) box on one side of the ray, or behind its origin.
            let corners = [
                (p.ybox.0, p.bbox.2),
                (p.ybox.0, p.bbox.3),
                (p.ybox.1, p.bbox.2),
                (p.ybox.1, p.bbox.3),
            ];
            let offs = corners.map(|(y, z)| off(y, z));
            if offs.iter().all(|&o| o > margin)
                || offs.iter().all(|&o| o < -margin)
                || corners.iter().all(|&(y, z)| along(y, z) < -margin)
            {
                continue;
            }
            // Seeds ahead of the origin only: the shell behind it (a stem
            // head above the water on a trimmed hull, say) converges to hits
            // the ray never reaches.
            let mut seeds = [(f64::INFINITY, 0.0f64, 0.0f64); 3];
            for q in &p.pts {
                if along(q.3, q.4) <= 0.0 {
                    continue;
                }
                let d = (q.2 - x_t).powi(2) + off(q.3, q.4).powi(2);
                if d < seeds[2].0 {
                    seeds[2] = (d, q.0, q.1);
                    seeds.sort_by(|a, b| a.0.total_cmp(&b.0));
                }
            }
            for &(d, u, v) in &seeds {
                if !d.is_finite() {
                    break;
                }
                let solved = newton_in_domain_uv(&p.surf, u, v, tol, |s, du, dv| {
                    (
                        [s[0] - x_t, off(s[1], s[2])],
                        [
                            [du[0], dv[0]],
                            [sd * du[1] * st - du[2] * ct, sd * dv[1] * st - dv[2] * ct],
                        ],
                    )
                });
                if let Some((s, (su, sv))) = solved {
                    let r = along(s[1], s[2]);
                    if r > min_reach && !is_end_cap(&p.surf, su, sv) {
                        hits.push((r, Some((pi, su, sv))));
                    }
                }
            }
        }
        self.tidy(hits)
    }

    /// The edge crossings that lie on the ray, as hits.
    fn edge_hits(&self, sd: f64, z0: f64, theta: f64) -> Vec<Hit> {
        let (st, ct) = theta.sin_cos();
        let (y_c, scale) = (self.y_c, self.scale);
        let off = |y: f64, z: f64| sd * (y - y_c) * st - (z - z0) * ct;
        let along = |y: f64, z: f64| sd * (y - y_c) * ct + (z - z0) * st;
        self.edges
            .iter()
            .filter(|&&(y, z)| off(y, z).abs() <= 1e-9 * scale && along(y, z) > 1e-12 * scale)
            .map(|&(y, z)| (along(y, z), None))
            .collect()
    }

    /// Sort hits by reach and merge coincident ones, keeping a patch hit
    /// over an edge crossing at the same place.
    fn tidy(&self, mut hits: Vec<Hit>) -> Vec<Hit> {
        let scale = self.scale;
        hits.sort_by(|a, b| a.0.total_cmp(&b.0));
        hits.dedup_by(|a, b| {
            let same = (a.0 - b.0).abs() <= 1e-6 * scale;
            if same && b.1.is_none() {
                b.1 = a.1;
            }
            same
        });
        hits
    }

    /// Shallowest wetted point of the section: its top on the centreplane (a
    /// stem face, a bulb top) or an edge crossing (an open shell's rim that
    /// has gone under, on a trimmed or deeply immersed hull with no deck).
    /// Rays hang from `(y_c, top)`, which for an open rim lies on the lid
    /// closing it — as outermost-fold half-breadths imply.
    fn top(&self) -> Option<f64> {
        let mut best = self.centerplane_top();
        for &(_, z) in &self.edges {
            if z >= 0.0 {
                best = Some(best.map_or(z, |m: f64| m.min(z)));
            }
        }
        best
    }

    /// Shallowest wetted point on the centreplane: by Newton on
    /// `x(u,v) = x`, `y(u,v) = y_c` from where the mesh cut crosses the
    /// centreplane (or, without a mesh, from presample seeds), and among
    /// edge crossings there (keel and stem lines are usually patch edges).
    fn centerplane_top(&self) -> Option<f64> {
        let (x_t, y_c, scale) = (self.x, self.y_c, self.scale);
        let tol = 1e-11 * scale;
        let margin = 0.05 * scale;
        let mut best: Option<f64> = None;
        let mut take = |z: f64| {
            if z >= 0.0 {
                best = Some(best.map_or(z, |b: f64| b.min(z)));
            }
        };
        for &(y, z) in &self.edges {
            if (y - y_c).abs() <= 1e-9 * scale {
                take(z);
            }
        }
        let system = |s: [f64; 3], du: [f64; 3], dv: [f64; 3]| {
            ([s[0] - x_t, s[1] - y_c], [[du[0], dv[0]], [du[1], dv[1]]])
        };
        if let Some(cut) = &self.cut {
            for [a, b] in cut {
                if (a.y - y_c) * (b.y - y_c) > 0.0 {
                    continue;
                }
                let s = if (b.y - a.y).abs() > 1e-300 {
                    ((y_c - a.y) / (b.y - a.y)).clamp(0.0, 1.0)
                } else {
                    0.5
                };
                if self.exact {
                    take(a.z + s * (b.z - a.z));
                    continue;
                }
                let (u, v) = (a.u + s * (b.u - a.u), a.v + s * (b.v - a.v));
                if let Some(p) = self.patches.get(a.patch) {
                    if let Some((q, (qu, qv))) = newton_in_domain_uv(&p.surf, u, v, tol, system) {
                        if !is_end_cap(&p.surf, qu, qv) {
                            take(q[2]);
                        }
                    }
                }
            }
            return best;
        }
        for p in self.patches {
            if x_t < p.bbox.0 - margin || x_t > p.bbox.1 + margin {
                continue;
            }
            let mut seeds: Vec<(f64, f64, f64)> = p
                .pts
                .iter()
                .filter(|q| q.4 >= -margin)
                .map(|q| ((q.2 - x_t).powi(2) + (q.3 - y_c).powi(2), q.0, q.1))
                .collect();
            seeds.sort_by(|a, b| a.0.total_cmp(&b.0));
            for &(_, u, v) in seeds.iter().take(4) {
                if let Some((s, (su, sv))) = newton_in_domain_uv(&p.surf, u, v, tol, system) {
                    if !is_end_cap(&p.surf, su, sv) {
                        take(s[2]);
                    }
                }
            }
        }
        best
    }
}

/// A hull's patches tessellated once, in the CAD frame, each vertex
/// remembering its patch and parameters — so a pose moves vertices (by the
/// same affine map as the control nets, under which a B-spline surface and
/// its points move together) rather than re-evaluating surfaces, and a
/// station's cut is found among the few triangles spanning it.
#[derive(Debug, Clone)]
pub(crate) struct Tessellation {
    pub(crate) verts: Vec<[f64; 3]>,
    /// `(patch, u, v)` of each vertex; empty for a mesh that is itself the
    /// geometry (an STL), whose cuts are then taken as exact.
    params: Vec<(usize, f64, f64)>,
    pub(crate) tris: Vec<[u32; 3]>,
}

impl Tessellation {
    /// A bare triangle mesh (CAD frame) that is the geometry itself: no
    /// patches behind it, so sections come straight from its facets.
    pub(crate) fn from_mesh(verts: Vec<[f64; 3]>, tris: Vec<[u32; 3]>) -> Self {
        Tessellation {
            verts,
            params: Vec::new(),
            tris,
        }
    }

    /// Grid-tessellate every patch finely enough that a facet is ~1/400 of
    /// the hull's largest extent: fine enough to catch every crossing
    /// (thin plates are their own patches, so still get cells), coarse
    /// enough to slice in microseconds. Accuracy comes from the Newton
    /// polish on the exact surface, not from the facets.
    pub(crate) fn new(surfs: &[NurbsSurface3]) -> Self {
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in surfs.iter().flat_map(|s| s.ctrl.iter()) {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let extent = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
        let h = (extent / 400.0).max(1e-12);
        let (mut verts, mut params, mut tris) = (Vec::new(), Vec::new(), Vec::new());
        for (pi, s) in surfs.iter().enumerate() {
            // Control-polygon lengths along u and v bound the patch's.
            let (nu, nv) = (s.n_ctrl_u, s.n_ctrl_v);
            let dist = |a: [f64; 3], b: [f64; 3]| {
                ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
            };
            let len_u = (0..nv)
                .map(|j| {
                    (1..nu)
                        .map(|i| dist(s.ctrl[i * nv + j], s.ctrl[(i - 1) * nv + j]))
                        .sum::<f64>()
                })
                .fold(0.0, f64::max);
            let len_v = (0..nu)
                .map(|i| {
                    (1..nv)
                        .map(|j| dist(s.ctrl[i * nv + j], s.ctrl[i * nv + j - 1]))
                        .sum::<f64>()
                })
                .fold(0.0, f64::max);
            let cells = |len: f64| ((len / h).ceil() as usize).clamp(8, 256);
            let (cu, cv) = (cells(len_u), cells(len_v));
            let (u0, u1) = s.u_domain();
            let (v0, v1) = s.v_domain();
            let base = verts.len() as u32;
            for i in 0..=cu {
                let u = u0 + (u1 - u0) * i as f64 / cu as f64;
                for j in 0..=cv {
                    let v = v0 + (v1 - v0) * j as f64 / cv as f64;
                    verts.push(s.point(u, v));
                    params.push((pi, u, v));
                }
            }
            let idx = |i: usize, j: usize| base + (i * (cv + 1) + j) as u32;
            for i in 0..cu {
                for j in 0..cv {
                    tris.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
                    tris.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
                }
            }
        }
        Tessellation {
            verts,
            params,
            tris,
        }
    }
}

/// A [`Tessellation`] at one pose, in the hull frame (z down from the
/// waterline), with its triangles bucketed by x so a station's cut touches
/// only the triangles that span it.
pub(crate) struct PosedMesh<'t> {
    tess: &'t Tessellation,
    verts: Vec<[f64; 3]>,
    x0: f64,
    bucket_w: f64,
    buckets: Vec<Vec<u32>>,
}

impl<'t> PosedMesh<'t> {
    /// Pose the vertices with `to_hull` (CAD frame → hull frame) and bucket
    /// the triangles.
    pub(crate) fn new(tess: &'t Tessellation, to_hull: impl Fn([f64; 3]) -> [f64; 3]) -> Self {
        let verts: Vec<[f64; 3]> = tess.verts.iter().map(|&p| to_hull(p)).collect();
        let (x0, x1) = verts
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), p| {
                (a.min(p[0]), b.max(p[0]))
            });
        let nb = (tess.tris.len() / 64).clamp(1, 4096);
        let bucket_w = ((x1 - x0) / nb as f64).max(1e-300);
        let mut buckets = vec![Vec::new(); nb];
        let slot = |x: f64| (((x - x0) / bucket_w) as usize).min(nb - 1);
        for (ti, t) in tess.tris.iter().enumerate() {
            let xs = t.map(|k| verts[k as usize][0]);
            let (lo, hi) = (xs[0].min(xs[1]).min(xs[2]), xs[0].max(xs[1]).max(xs[2]));
            for b in slot(lo)..=slot(hi) {
                buckets[b].push(ti as u32);
            }
        }
        PosedMesh {
            tess,
            verts,
            x0,
            bucket_w,
            buckets,
        }
    }

    /// Every y at which the transverse line through `(x, z)` (hull frame)
    /// crosses the mesh, sorted, coincident crossings merged.
    fn transverse(&self, x: f64, z: f64, scale: f64) -> Vec<f64> {
        let mut ys: Vec<f64> = Vec::new();
        for [a, b] in self.cut(x) {
            if (a.z - z) * (b.z - z) > 0.0 {
                continue;
            }
            let dz = b.z - a.z;
            let y = if dz.abs() > 1e-300 {
                a.y + (z - a.z) / dz * (b.y - a.y)
            } else {
                0.5 * (a.y + b.y)
            };
            ys.push(y);
        }
        ys.sort_by(f64::total_cmp);
        ys.dedup_by(|a, b| (*a - *b).abs() <= 1e-6 * scale);
        ys
    }

    /// The plane `x = const` through the mesh: one segment per triangle it
    /// crosses, its ends where it crosses the triangle's edges.
    fn cut(&self, x: f64) -> Vec<[Cut; 2]> {
        let b = (x - self.x0) / self.bucket_w;
        if b < 0.0 || b >= self.buckets.len() as f64 + 1e-9 {
            return Vec::new();
        }
        let b = (b as usize).min(self.buckets.len() - 1);
        let mut out = Vec::new();
        for &ti in &self.buckets[b] {
            let t = self.tess.tris[ti as usize];
            let [p0, p1, p2] = t.map(|k| self.verts[k as usize]);
            if normal_is_cap(cross3(sub3(p1, p0), sub3(p2, p0))) {
                continue;
            }
            let mut ends = [None, None];
            let mut n = 0;
            for (a, c) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let (pa, pc) = (self.verts[a as usize], self.verts[c as usize]);
                let (da, dc) = (pa[0] - x, pc[0] - x);
                // Half-open, so a plane through a vertex is counted once.
                if (da >= 0.0) == (dc >= 0.0) {
                    continue;
                }
                let s = da / (da - dc);
                let param = |k: u32| {
                    self.tess
                        .params
                        .get(k as usize)
                        .copied()
                        .unwrap_or((usize::MAX, 0.0, 0.0))
                };
                let (qa, qc) = (param(a), param(c));
                if n < 2 {
                    ends[n] = Some(Cut {
                        y: pa[1] + s * (pc[1] - pa[1]),
                        z: pa[2] + s * (pc[2] - pa[2]),
                        patch: qa.0,
                        u: qa.1 + s * (qc.1 - qa.1),
                        v: qa.2 + s * (qc.2 - qa.2),
                    });
                    n += 1;
                }
            }
            if let [Some(a), Some(c)] = ends {
                out.push([a, c]);
            }
        }
        out
    }
}

/// A surface within ~18° of parallel to the station planes — facing fore or
/// aft, like a transom face or a flat stem face — is an **end cap**, not
/// part of the shell a section is cut from. A station plane meets it, if at
/// all, along an ill-conditioned line (the whole face lies in the plane of an
/// upright transom; a trimmed one crosses it at a grazing angle), and taking
/// that line as section boundary turns the end section into a degenerate
/// sliver — so the hull would close over one span, as if a transom wall were
/// there, the moment it trims. Sections are cut from the side and bottom
/// shell alone; where that shell stops, the hull ends with its full end
/// section, open, at any trim.
const END_CAP_NX: f64 = 0.95;

fn normal_is_cap(n: [f64; 3]) -> bool {
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    len > 0.0 && n[0].abs() > END_CAP_NX * len
}

fn is_end_cap(surf: &NurbsSurface3, u: f64, v: f64) -> bool {
    let (_, du, dv) = surf.eval1(u, v);
    normal_is_cap(cross3(du, dv))
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// 2-D Newton on a patch for `F(S(u,v)) = 0`, with `F` and its Jacobian with
/// respect to `(u, v)` supplied from the point and partials. Steps are cut
/// back along their own direction to stay inside the parameter domain —
/// clamping each coordinate separately bends the step off course, and
/// cannot settle on a solution that lies on the domain's edge. Returns the
/// converged point and its parameters.
fn newton_in_domain_uv(
    surf: &NurbsSurface3,
    mut u: f64,
    mut v: f64,
    tol: f64,
    f: impl Fn([f64; 3], [f64; 3], [f64; 3]) -> ([f64; 2], [[f64; 2]; 2]),
) -> Option<([f64; 3], (f64, f64))> {
    let (u0, u1) = surf.u_domain();
    let (v0, v1) = surf.v_domain();
    for _ in 0..60 {
        let (s, du, dv) = surf.eval1(u, v);
        let (r, j) = f(s, du, dv);
        if r[0].abs() + r[1].abs() < tol {
            return Some((s, (u, v)));
        }
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det.abs() < 1e-30 {
            return None;
        }
        let su = (r[0] * j[1][1] - r[1] * j[0][1]) / det;
        let sv = (-r[0] * j[1][0] + r[1] * j[0][0]) / det;
        // Largest fraction of the step that stays in the domain.
        let mut t = 1.0f64;
        for (x, dx, lo, hi) in [(u, -su, u0, u1), (v, -sv, v0, v1)] {
            if x + dx > hi {
                t = t.min((hi - x) / dx);
            } else if x + dx < lo {
                t = t.min((lo - x) / dx);
            }
        }
        if t <= 0.0 {
            // Pinned against the boundary: slide along it.
            u = (u - su).clamp(u0, u1);
            v = (v - sv).clamp(v0, v1);
        } else {
            u = (u - t * su).clamp(u0, u1);
            v = (v - t * sv).clamp(v0, v1);
        }
    }
    None
}

/// The hull's end in direction `dir = ±1` from a station `x_in` known to
/// have a section: step outward (from `first_step`, growing) until sections
/// stop, then bisect the transition to ~1e-12 of the length. Returns the **inside** end of the
/// final bracket, so the end station always has its section — which at a
/// transom is the whole transom, and must not be read as empty (that would
/// close the hull over one span: a spurious cliff in `Z(x)`).
#[allow(clippy::too_many_arguments)]
fn hull_end(
    patches: &[Patch],
    x_in: f64,
    dir: f64,
    y_c: f64,
    sides: &[f64],
    length: f64,
    scale: f64,
    first_step: f64,
    mesh: Option<&PosedMesh>,
) -> f64 {
    // Start from a station that has a section (the presample's bracket can
    // sit in a recess or an open skin), moving inward if need be.
    let probe = |x: f64| Station::new(patches, mesh, x, y_c, scale).has_section(sides);
    let mut inside = x_in;
    for k in 1..=20 {
        if probe(inside) {
            break;
        }
        inside = x_in - dir * 0.01 * length * k as f64;
    }
    let mut step = first_step;
    let mut outside = None;
    for _ in 0..40 {
        let x = inside + dir * step;
        if probe(x) {
            inside = x;
        } else {
            outside = Some(x);
            break;
        }
        step *= 1.5;
    }
    let Some(mut outside) = outside else {
        return inside;
    };
    for _ in 0..60 {
        if (outside - inside).abs() <= 1e-12 * length {
            break;
        }
        let mid = 0.5 * (inside + outside);
        if probe(mid) {
            inside = mid;
        } else {
            outside = mid;
        }
    }
    inside
}

/// One station's section, sampled for [`SectionNodes::from_polar`]: the ray
/// origin depth `z₀` (0 when the section reaches the waterline, else its top
/// centreplane point), the section's beam and depth scales, each scaled
/// ray's scaled reach averaged over the shell's sides, and the largest
/// side-to-side difference [m]. `None` when the plane misses the hull or
/// some ray finds no shell.
#[allow(clippy::too_many_arguments)]
fn sample_section(
    patches: &[Patch],
    mesh: Option<&PosedMesh>,
    x: f64,
    y_c: f64,
    sides: &[f64],
    thetas: &[f64],
    scale: f64,
    ambiguous: &mut usize,
) -> Option<(f64, f64, f64, Vec<f64>, f64)> {
    let station = Station::new(patches, mesh, x, y_c, scale);
    sample_station(&station, sides, thetas, scale, ambiguous)
}

fn sample_station(
    station: &Station,
    sides: &[f64],
    thetas: &[f64],
    scale: f64,
    ambiguous: &mut usize,
) -> Option<(f64, f64, f64, Vec<f64>, f64)> {
    let half_pi = std::f64::consts::FRAC_PI_2;
    // Physical-angle reach, averaged over the sides.
    let reach = |z0: f64, t: f64, amb: &mut usize, asym: &mut f64| -> Option<f64> {
        let mut rs = [0.0f64; 2];
        for (k, &sd) in sides.iter().enumerate() {
            rs[k] = station.reach(sd, z0, t, amb)?;
        }
        if sides.len() == 2 {
            *asym = asym.max((rs[0] - rs[1]).abs());
        }
        Some(rs[..sides.len()].iter().sum::<f64>() / sides.len() as f64)
    };
    let mut asym = 0.0f64;
    // A section that reaches the waterline is swept from the centreplane at
    // the surface; one that doesn't (forefoot ahead of the waterline entry,
    // a bulb, an open shell whose rims are under) from its shallowest point.
    let wl = reach(0.0, 0.0, ambiguous, &mut asym);
    let z0 = if wl.is_some() { 0.0 } else { station.top()? };
    let depth = reach(z0, half_pi, ambiguous, &mut asym)?;
    // Beam scale: the widest of a coarse fan (the waterline, usually).
    let mut beam = wl.unwrap_or(0.0);
    for k in 1..8 {
        let t = half_pi * k as f64 / 8.0;
        if let Some(r) = reach(z0, t, ambiguous, &mut asym) {
            beam = beam.max(r * t.cos());
        }
    }
    // A section with no depth or no beam to speak of (a hull's very tip) is
    // no section: its fan, scaled to nothing in one direction, would sweep
    // along the other and pick up whatever lies there. Its area is nil either
    // way; an end station without one is empty, which is what it is.
    if !(beam > 1e-6 * scale && depth > 1e-6 * scale) {
        return None;
    }
    let mut radii = Vec::with_capacity(thetas.len());
    for &t in thetas {
        let (st, ct) = t.sin_cos();
        let (dy, dz) = (beam * ct, depth * st);
        let phys = dz.atan2(dy);
        let r = reach(z0, phys, ambiguous, &mut asym)?;
        radii.push(r / dy.hypot(dz));
    }
    Some((z0, beam, depth, radii, asym))
}

/// Barycentric interpolation on Chebyshev–Lobatto points (as `thetas` are
/// laid out in [`sectional_cluster`]): spectrally accurate for the smooth
/// ray reach of a fair section.
struct Lobatto<'a> {
    t: &'a [f64],
    f: &'a [f64],
    w: Vec<f64>,
}

impl<'a> Lobatto<'a> {
    fn new(t: &'a [f64], f: &'a [f64]) -> Self {
        let n = t.len();
        let w = (0..n)
            .map(|k| {
                let s = if k % 2 == 0 { 1.0 } else { -1.0 };
                if k == 0 || k + 1 == n {
                    0.5 * s
                } else {
                    s
                }
            })
            .collect();
        Lobatto { t, f, w }
    }

    fn eval(&self, x: f64) -> f64 {
        let (mut num, mut den) = (0.0, 0.0);
        for ((&tk, &fk), &wk) in self.t.iter().zip(self.f).zip(&self.w) {
            let d = x - tk;
            if d == 0.0 {
                return fk;
            }
            num += wk * fk / d;
            den += wk / d;
        }
        num / den
    }
}

// ---------------------------------------------------------------------------
// IGES export
// ---------------------------------------------------------------------------

/// Serialize surfaces to a minimal IGES 5.3 file: one B-spline surface
/// (entity 128) per patch, coordinates in **metres** (units flag 6). The
/// output round-trips through [`parse`] and imports into CAD systems.
///
/// Surfaces carrying a bounded-surface restriction (`trim_uv`) are written as
/// plain 128 entities whose parameter range (`U0..U1`, `V0..V1`) is the
/// restricted box; readers that honour the range see the bounded patch,
/// others see the full base surface.
///
/// `product` names the model in the global section (product ID / file name).
pub fn write(surfaces: &[NurbsSurface3], product: &str) -> Result<String> {
    if surfaces.is_empty() {
        return Err(Error::InvalidInput("no surfaces to write".into()));
    }
    for (i, s) in surfaces.iter().enumerate() {
        let n = s.n_ctrl_u * s.n_ctrl_v;
        if s.ctrl.len() != n
            || s.weights.len() != n
            || s.knots_u.len() != s.n_ctrl_u + s.degree_u + 1
            || s.knots_v.len() != s.n_ctrl_v + s.degree_v + 1
            || s.ctrl.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput(format!(
                "surface {i} is inconsistent (control/weight/knot counts or \
                 non-finite coordinates)"
            )));
        }
    }
    let product: String = product
        .chars()
        .map(|c| {
            if c.is_ascii_graphic() || c == ' ' {
                c
            } else {
                '_'
            }
        })
        .take(60)
        .collect();
    let product = if product.trim().is_empty() {
        "michell".to_string()
    } else {
        product
    };
    let dt = iges_datetime();
    let max_coord = surfaces
        .iter()
        .flat_map(|s| s.ctrl.iter())
        .flat_map(|p| p.iter())
        .fold(1.0f64, |m, &c| m.max(c.abs()));

    // Global section (IGES 5.3 field order).
    let holl = |s: &str| format!("{}H{}", s.len(), s);
    let g_tokens: Vec<String> = vec![
        "1H,".into(),                                            //  1 parameter delimiter
        "1H;".into(),                                            //  2 record delimiter
        holl(&product),                                          //  3 product ID (sender)
        holl(&product),                                          //  4 file name
        holl("michell"),                                         //  5 native system ID
        holl(&format!("michell {}", env!("CARGO_PKG_VERSION"))), // 6 preprocessor
        "32".into(),                                             //  7 integer bits
        "38".into(),                                             //  8 single max exponent
        "6".into(),                                              //  9 single sig digits
        "308".into(),                                            // 10 double max exponent
        "15".into(),                                             // 11 double sig digits
        holl(&product),                                          // 12 product ID (receiver)
        "1.".into(),                                             // 13 model space scale
        "6".into(),                                              // 14 units flag: metres
        holl("M"),                                               // 15 units name
        "1".into(),                                              // 16 line weight gradations
        "0.01".into(),                                           // 17 max line weight
        holl(&dt),                                               // 18 file generation time
        fmt_real(1e-9),                                          // 19 minimum resolution
        fmt_real(max_coord),                                     // 20 approx max coordinate
        holl("michell"),                                         // 21 author
        holl("michell"),                                         // 22 organization
        "11".into(),                                             // 23 version flag (5.3)
        "0".into(),                                              // 24 drafting standard
        holl(&dt),                                               // 25 last modified
    ];
    let g_lines = pack_tokens(&g_tokens, 72);

    // Parameter + directory sections.
    let mut p_lines: Vec<(String, usize)> = Vec::new(); // (data, owner DE)
    let mut d_lines: Vec<String> = Vec::new();
    for (k, s) in surfaces.iter().enumerate() {
        let (nu, nv) = (s.n_ctrl_u, s.n_ctrl_v);
        let mut t: Vec<String> = vec![
            "128".into(),
            (nu - 1).to_string(),
            (nv - 1).to_string(),
            s.degree_u.to_string(),
            s.degree_v.to_string(),
            "0".into(), // PROP1: not closed in u
            "0".into(), // PROP2: not closed in v
            if s.is_polynomial() { "1" } else { "0" }.into(),
            "0".into(), // PROP4: not periodic in u
            "0".into(), // PROP5: not periodic in v
        ];
        t.extend(s.knots_u.iter().map(|&v| fmt_real(v)));
        t.extend(s.knots_v.iter().map(|&v| fmt_real(v)));
        // Weights and points with the u index varying fastest (see parse).
        for j in 0..nv {
            for i in 0..nu {
                t.push(fmt_real(s.weights[i * nv + j]));
            }
        }
        for j in 0..nv {
            for i in 0..nu {
                let p = s.ctrl[i * nv + j];
                t.extend(p.iter().map(|&v| fmt_real(v)));
            }
        }
        let (u0, u1) = s.u_domain();
        let (v0, v1) = s.v_domain();
        t.extend([u0, u1, v0, v1].iter().map(|&v| fmt_real(v)));

        let de_seq = 2 * k + 1;
        let pd_ptr = p_lines.len() + 1;
        let chunk = pack_tokens(&t, 64);
        let pd_count = chunk.len();
        p_lines.extend(chunk.into_iter().map(|l| (l, de_seq)));

        let f = |v: usize| format!("{v:>8}");
        d_lines.push(format!(
            "{}{}{}{}{}{}{}{}{:>8}",
            f(128),
            f(pd_ptr),
            f(0),
            f(0),
            f(0),
            f(0),
            f(0),
            f(0),
            "00000000"
        ));
        d_lines.push(format!(
            "{}{}{}{}{}{:8}{:8}{:8}{}",
            f(128),
            f(0),
            f(0),
            f(pd_count),
            f(0),
            "",
            "",
            "",
            f(0)
        ));
    }

    let mut out = String::new();
    let mut line = |data: &str, letter: char, seq: usize| {
        out.push_str(&format!("{data:<72}{letter}{seq:>7}\n"));
    };
    line(&format!("{product} - michell IGES export"), 'S', 1);
    for (i, g) in g_lines.iter().enumerate() {
        line(g, 'G', i + 1);
    }
    for (i, d) in d_lines.iter().enumerate() {
        line(d, 'D', i + 1);
    }
    for (i, (p, de)) in p_lines.iter().enumerate() {
        line(&format!("{p:<64}{de:>8}"), 'P', i + 1);
    }
    let totals = format!(
        "S{:>7}G{:>7}D{:>7}P{:>7}",
        1,
        g_lines.len(),
        d_lines.len(),
        p_lines.len()
    );
    line(&totals, 'T', 1);
    Ok(out)
}

/// The Wigley parabolic hull `y = ±(B/2)(1 − (2x/L)²)(1 − (z/T)²)` as its
/// **exact** mirrored pair of CAD patches (x along the hull, `x ∈ [−L/2,
/// L/2]`; z up, the design waterline at `z = 0` and the keel at
/// `z = −draft`; centreplane `y = 0`), returned as `[starboard (+y), port
/// (−y)]`. Write it out with [`write`], or cut it directly through
/// [`source_fleet_from_surfaces`]. The classic proportions are `L/B = 10`,
/// `B/T = 1.6`, e.g. `wigley_surfaces(10.0, 1.0, 0.625)`.
pub fn wigley_surfaces(length: f64, beam: f64, draft: f64) -> Result<[NurbsSurface3; 2]> {
    let s = crate::hulls::wigley_surface(length, beam, draft)?;
    Ok(halfbreadth_surfaces(&s, 0.0, 0.0))
}

/// Convert a half-breadth spline `y = f(x, z')` (the crate's hull frame,
/// z' downward) into the mirrored pair of CAD-frame (z up) patches its graph
/// describes — **exact**, because linear precision puts the Greville
/// abscissae of each knot vector on the graph's coordinate lines.
/// `z_top_cad` is the CAD height of z' = 0 (the design waterline for a
/// wetted hull; the band top for a full-band body); `centerplane` is the
/// transverse position the two sides mirror about. Returned as
/// `[starboard (+y), port (−y)]`.
pub(crate) fn halfbreadth_surfaces(
    s: &BSplineSurface,
    centerplane: f64,
    z_top_cad: f64,
) -> [NurbsSurface3; 2] {
    let greville = |knots: &[f64], p: usize, n: usize| -> Vec<f64> {
        (0..n)
            .map(|i| knots[i + 1..i + 1 + p].iter().sum::<f64>() / p as f64)
            .collect()
    };
    let gx = greville(s.knots_x(), s.degree_x(), s.n_ctrl_x());
    let gz = greville(s.knots_z(), s.degree_z(), s.n_ctrl_z());
    let (nx, nz) = (s.n_ctrl_x(), s.n_ctrl_z());
    #[allow(clippy::needless_range_loop)]
    let build = |side: f64| -> NurbsSurface3 {
        let mut ctrl = Vec::with_capacity(nx * nz);
        for i in 0..nx {
            for j in 0..nz {
                ctrl.push([
                    gx[i],
                    centerplane + side * s.control()[i * nz + j],
                    z_top_cad - gz[j],
                ]);
            }
        }
        NurbsSurface3 {
            degree_u: s.degree_x(),
            degree_v: s.degree_z(),
            knots_u: s.knots_x().to_vec(),
            knots_v: s.knots_z().to_vec(),
            n_ctrl_u: nx,
            n_ctrl_v: nz,
            ctrl,
            weights: vec![1.0; nx * nz],
            trim_uv: None,
        }
    };
    [build(1.0), build(-1.0)]
}

/// Format a real for an IGES parameter field: Rust's shortest round-trip
/// representation, with the decimal point IGES requires; extreme magnitudes
/// switch to `E` exponents so no token can outgrow a parameter line.
fn fmt_real(v: f64) -> String {
    let a = v.abs();
    let s = if a != 0.0 && !(1e-4..1e7).contains(&a) {
        format!("{v:E}")
    } else {
        format!("{v}")
    };
    if let Some(e) = s.find('E') {
        let (m, ex) = s.split_at(e);
        let m = if m.contains('.') {
            m.to_string()
        } else {
            format!("{m}.0")
        };
        format!("{m}E{}", &ex[1..])
    } else if s.contains('.') {
        s
    } else {
        format!("{s}.")
    }
}

/// Join parameter tokens with `,` (record-terminated by `;`) and pack them
/// into lines of at most `width` columns, breaking only between tokens.
fn pack_tokens(tokens: &[String], width: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut cur = String::new();
    for (i, t) in tokens.iter().enumerate() {
        let delim = if i + 1 == tokens.len() { ';' } else { ',' };
        let piece = format!("{t}{delim}");
        if !cur.is_empty() && cur.len() + piece.len() > width {
            lines.push(std::mem::take(&mut cur));
        }
        cur.push_str(&piece);
    }
    if !cur.is_empty() {
        lines.push(cur);
    }
    lines
}

/// Current UTC time as the IGES `YYYYMMDD.HHMMSS` timestamp.
fn iges_datetime() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let (days, rem) = (secs / 86400, secs % 86400);
    // Civil date from day count (Howard Hinnant's algorithm).
    let z = days as i64 + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe as i64 + era * 400 + i64::from(m <= 2);
    format!(
        "{y:04}{m:02}{d:02}.{:02}{:02}{:02}",
        rem / 3600,
        (rem % 3600) / 60,
        rem % 60
    )
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
            trim_uv: None,
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
        let ybox = pts
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), q| {
                (lo.min(q.3), hi.max(q.3))
            });
        Patch {
            surf,
            pts,
            grid_n: n,
            bbox,
            ybox,
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
                    (ys[0] - s[1]).abs() < 1e-9,
                    "u={fu} v={fv}: y {} vs {}",
                    ys[0],
                    s[1]
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
                trim_uv: None,
            }
        };
        let patches = vec![as_patch(wall(1.0), 9), as_patch(wall(-1.0), 9)];
        let ys = shell_intersections(&patches, 5.0, 1.0, 10.0);
        assert_eq!(ys.len(), 2, "{ys:?}");
        assert!((ys[0] + 1.0).abs() < 1e-9 && (ys[1] - 1.0).abs() < 1e-9);
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

    #[test]
    fn real_formatting() {
        assert_eq!(fmt_real(1.5), "1.5");
        assert_eq!(fmt_real(1.0), "1.");
        assert_eq!(fmt_real(-3.0), "-3.");
        assert_eq!(fmt_real(1e-9), "1.0E-9");
        assert_eq!(fmt_real(-2.5e-7), "-2.5E-7");
        assert_eq!(fmt_real(0.0), "0.");
        for &v in &[0.1, -123.456, 1e-9, 6.02e23, 12.0] {
            assert_eq!(parse_number(&fmt_real(v)).unwrap(), v);
        }
    }

    /// The Wigley helper written out and read back: one two-sided hull on
    /// the centreplane, the draft below CAD z = 0, and the analytic volume
    /// `∇ = (4/9) L B T`.
    #[test]
    fn wigley_surfaces_roundtrip_through_write_and_source_fleet() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let text = write(&wigley_surfaces(l, b, t).unwrap(), "wigley").unwrap();
        let src = source_fleet(&text, 0.0).unwrap();
        assert_eq!(src.len(), 1);
        assert!((src.hull_z_bottom(0) + t).abs() < 1e-12);
        let cut = src
            .situate_sectional(
                0,
                0.0,
                &HullPose::default(),
                &Platform::default(),
                &SectionalOptions::default(),
            )
            .unwrap()
            .expect("the Wigley is wet at its design waterline");
        assert!(cut.report.two_sided);
        assert!(
            cut.report.centerplane.abs() < 1e-9,
            "{}",
            cut.report.centerplane
        );
        assert!(
            (cut.report.draft - t).abs() < 1e-9,
            "draft {}",
            cut.report.draft
        );
        let want = 4.0 / 9.0 * l * b * t;
        let got = cut.hull.displaced_volume();
        assert!((got - want).abs() < 1e-6 * want, "volume {got} vs {want}");
        assert!(wigley_surfaces(0.0, 1.0, 1.0).is_err());
    }

    #[test]
    fn write_roundtrips_through_parse() {
        let a = nonlinear_surface();
        let mut b = nonlinear_surface();
        for p in b.ctrl.iter_mut() {
            p[0] += 3.25;
            p[1] = -p[1];
            p[2] -= 0.75;
        }
        let text = write(&[a.clone(), b.clone()], "roundtrip").unwrap();
        // Fixed-width layout: every line is 80 columns with the section
        // letter in column 73.
        for line in text.lines() {
            assert_eq!(line.len(), 80, "{line:?}");
            assert!(matches!(
                line.as_bytes()[72],
                b'S' | b'G' | b'D' | b'P' | b'T'
            ));
        }
        let file = parse(&text).unwrap();
        assert_eq!(file.units_scale, 1.0);
        assert_eq!(file.surfaces.len(), 2);
        for (got, want) in file.surfaces.iter().zip([&a, &b]) {
            assert_eq!(got.degree_u, want.degree_u);
            assert_eq!(got.degree_v, want.degree_v);
            assert_eq!(got.knots_u, want.knots_u);
            assert_eq!(got.knots_v, want.knots_v);
            assert_eq!(got.weights, want.weights);
            for (g, w) in got.ctrl.iter().zip(&want.ctrl) {
                for c in 0..3 {
                    assert_eq!(g[c], w[c], "control point mismatch");
                }
            }
        }
    }

    #[test]
    fn halfbreadth_graph_is_exact() {
        // An arbitrary half-breadth spline; its graph surfaces must evaluate
        // to (x, ±f(x, z) + yc, z_top − z) exactly.
        let s = BSplineSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 4.0, 10.0, 10.0, 10.0],
            vec![0.0, 0.0, 0.0, 1.5, 1.5, 1.5],
            vec![
                0.0, 0.0, 0.0, //
                0.8, 0.6, 0.1, //
                1.0, 0.9, 0.3, //
                0.2, 0.1, 0.0,
            ],
        )
        .unwrap();
        let (yc, ztop) = (1.9, 0.42);
        let [stbd, port] = halfbreadth_surfaces(&s, yc, ztop);
        for &(u, v) in &[(0.0, 0.0), (2.7, 0.3), (5.0, 1.1), (9.3, 1.5), (10.0, 0.7)] {
            let f = s.eval(u, v);
            let (p, _, _) = stbd.eval1(u, v);
            assert!((p[0] - u).abs() < 1e-12, "x at ({u},{v}): {}", p[0]);
            assert!((p[1] - (yc + f)).abs() < 1e-12, "y at ({u},{v})");
            assert!((p[2] - (ztop - v)).abs() < 1e-12, "z at ({u},{v})");
            let (q, _, _) = port.eval1(u, v);
            assert!((q[1] - (yc - f)).abs() < 1e-12, "port y at ({u},{v})");
        }
    }

    #[test]
    fn posed_export_matches_situate_transform() {
        // apply_pose + sinkage re-expression: a point on the DWL moves down
        // by dz + sinkage relative to the fixed CAD waterline.
        let mut surfs = vec![nonlinear_surface()];
        let pose = HullPose {
            dx: 1.0,
            dy: -0.5,
            dz: 0.2,
            ..HullPose::default()
        };
        let platform = Platform {
            sinkage: 0.1,
            ..Platform::default()
        };
        let before = surfs[0].ctrl[3];
        apply_pose(&mut surfs, 0.0, &pose, &platform);
        let after = surfs[0].ctrl[3];
        assert!((after[0] - (before[0] + 1.0)).abs() < 1e-12);
        assert!((after[1] - (before[1] - 0.5)).abs() < 1e-12);
        assert!((after[2] - (before[2] - 0.3)).abs() < 1e-12);
    }
}

#[cfg(test)]
mod sectional_ends {
    use super::*;

    fn e12() -> Option<SourceFleet> {
        let text = crate::cad_fixture("e12.igs")?;
        Some(source_fleet(&text, -0.95).unwrap())
    }

    /// e12's transom is a flat face its side skins run on past. Trimmed bow
    /// up, the face tilts across the aft station plane; read as section
    /// boundary it collapsed the end section to a sliver, closing the hull
    /// over one span (the near-field force doubled from trim 0 to +0.02°).
    /// Faces are end caps, not shell: the aft station keeps the full
    /// transom section at every trim.
    #[test]
    fn a_trimmed_transom_keeps_its_section() {
        let Some(fleet) = e12() else {
            return;
        };
        let idx = (0..fleet.len())
            .max_by_key(|&i| fleet.hulls[i].len())
            .unwrap();
        let so = SectionalOptions {
            waterline_z: -0.95,
            stations: 61,
            rays: 17,
            ..Default::default()
        };
        let mut last_vol = 0.0;
        for trim_deg in [-0.1f64, 0.0, 0.02, 0.05, 0.1] {
            let plat = Platform {
                sinkage: 0.012,
                trim: trim_deg.to_radians(),
                pivot_x: 4.5,
            };
            let h = fleet
                .situate_sectional(idx, -0.95, &HullPose::default(), &plat, &so)
                .unwrap()
                .unwrap();
            let (st, _) = h.hull.depth_integral_curve(0.0, 1);
            assert!(
                st[0].1 > 0.9 * st[1].1,
                "trim {trim_deg}°: end section {:.3e} against its neighbour's {:.3e}",
                st[0].1,
                st[1].1
            );
            let vol = h.hull.displaced_volume();
            assert!(vol > last_vol, "bow-up trim sinks this stern deeper: {vol}");
            last_vol = vol;
        }
    }
}

#[cfg(test)]
mod pose_timing {
    use super::*;

    /// How long a re-pose takes: the per-pose work (control-net pose,
    /// presample, frame and ends, sections) on real CAD at a few nearby
    /// poses. A report, not a check.
    #[test]
    #[ignore = "timing report"]
    fn sectional_repose_cost() {
        for (file, wl) in [("ama.igs", 0.0), ("e12.igs", -0.95)] {
            let Some(text) = crate::cad_fixture(file) else {
                continue;
            };
            let t = std::time::Instant::now();
            let fleet = source_fleet(&text, wl).unwrap();
            let tris: usize = fleet.meshes.iter().map(|m| m.tris.len()).sum();
            eprintln!(
                "{file}: load (parse, cluster, tessellate) {:.1} ms, {tris} triangles",
                1e3 * t.elapsed().as_secs_f64()
            );
            let idx = (0..fleet.len())
                .max_by_key(|&i| fleet.hulls[i].len())
                .unwrap();
            let opts = SectionalOptions {
                waterline_z: wl,
                ..Default::default()
            };
            // A cold start, then a run of nearby poses as an equilibrium
            // solve would visit them, each warm-started from the last and
            // checked against a cold import of the same pose.
            let mut state = SectionalState::default();
            for (sink, trim_deg) in [
                (0.0, 0.0),
                (0.005, 0.0),
                (0.01, 0.0),
                (0.01, 0.1),
                (0.01, 0.2),
                (0.012, 0.3),
            ] {
                let plat = Platform {
                    sinkage: sink,
                    trim: f64::to_radians(trim_deg),
                    pivot_x: 0.0,
                };
                let t = std::time::Instant::now();
                let warm = fleet
                    .situate_sectional_warm(idx, wl, &HullPose::default(), &plat, &opts, &mut state)
                    .unwrap()
                    .unwrap();
                let t_warm = t.elapsed().as_secs_f64();
                let t = std::time::Instant::now();
                let cold = fleet
                    .situate_sectional(idx, wl, &HullPose::default(), &plat, &opts)
                    .unwrap()
                    .unwrap();
                let t_cold = t.elapsed().as_secs_f64();
                let dv = (warm.hull.displaced_volume() / cold.hull.displaced_volume() - 1.0).abs();
                eprintln!(
                    "{file} sinkage {sink} trim {trim_deg}°: warm {:.1} ms, cold {:.1} ms, volume {:.5} (warm vs cold {dv:.1e})",
                    1e3 * t_warm,
                    1e3 * t_cold,
                    warm.hull.displaced_volume()
                );
            }
        }
    }
}

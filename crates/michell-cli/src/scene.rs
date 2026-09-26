//! `michell field` — a viewer-neutral **scene** of a sectional hull (or
//! fleet) at one speed: named geometries with named quantities on them, for
//! any 3-D viewer to draw (the polyscope script in `python/tools/`, the web
//! app).
//!
//! Frame: `x` forward (the ship advances toward +x), `y` transverse, `z`
//! **up** from the waterline, metres. The file is JSON:
//!
//! ```json
//! { "michell": "scene", "version": 1,
//!   "meta": { "speed": U, "froude": Fn, "nu": ν, "closure": "...", ... },
//!   "objects": [
//!     { "name": "...", "kind": "mesh",   "vertices": [[x,y,z],...], "faces": [[i,j,k],...],
//!       "quantities": [{ "name": "...", "on": "vertices", "values": [...] }], "note": "..." },
//!     { "name": "...", "kind": "curves", "vertices": [...], "edges": [[i,j],...], "quantities": [...] }
//!   ] }
//! ```
//!
//! Per hull: the CAD hull as its posed tessellation; the station curves the
//! physics integrates, carrying each station's section half-area, depth
//! integral and source strength `∂Z/∂x` at the transverse wave's decay rate;
//! and, for a closed transom, the virtual appendage. For the fleet: the
//! far-field free-wave elevation ζ(x, y) as a height-field mesh.

use michell::iges::{HullPose, Platform, SectionalImport, SourceFleet};
use michell::sectional::SectionalHull;
use michell::{Conditions, FreeWaveSpectrum, Placement, TransomClosure};
use std::fmt::Write as _;

/// A geometry and the quantities on its vertices.
struct Object {
    name: String,
    kind: &'static str,
    vertices: Vec<[f64; 3]>,
    /// Triangles for a mesh, segments for curves.
    cells: Vec<Vec<usize>>,
    quantities: Vec<(String, Vec<f64>)>,
    note: Option<String>,
}

/// Everything a scene needs about one hull.
pub(crate) struct SceneHull<'a> {
    pub name: String,
    pub import: &'a SectionalImport,
    pub placement: Placement,
    /// Where the hull came from, for its CAD tessellation.
    pub source: &'a SourceFleet,
    pub index: usize,
    pub waterline_z: f64,
}

/// Region and resolution of the free-surface grid.
pub(crate) struct Surface {
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub nx: usize,
    pub ny: usize,
}

pub(crate) fn build(
    hulls: &[SceneHull],
    cond: &Conditions,
    closure: TransomClosure,
    surface: &Surface,
    l_ref: f64,
) -> Result<String, String> {
    let nu = cond.gravity / (cond.speed * cond.speed);
    let mut objects = Vec::new();
    for h in hulls {
        let pl = h.placement;
        let hull = &h.import.hull;
        // The CAD hull itself, posed, whole.
        let (verts, tris) = h
            .source
            .posed_tessellation(h.index, h.waterline_z, &HullPose::default(), &Platform::default())
            .map_err(|e| e.to_string())?;
        let dy = pl.y - h.import.placement.y;
        let verts: Vec<[f64; 3]> = verts
            .iter()
            .map(|p| [p[0] + pl.x, p[1] + dy, p[2]])
            .collect();
        let height = verts.iter().map(|p| p[2]).collect();
        objects.push(Object {
            name: format!("{} hull (CAD)", h.name),
            kind: "mesh",
            cells: tris.iter().map(|t| t.iter().map(|&i| i as usize).collect()).collect(),
            vertices: verts,
            quantities: vec![("height above waterline [m]".into(), height)],
            note: Some("the posed CAD tessellation (display only; the physics uses the stations)".into()),
        });
        objects.push(stations(h, hull, pl, nu));
        if let Some(app) = appendage(h, hull, pl, nu, closure) {
            objects.push(app);
        }
    }
    objects.push(free_surface(hulls, cond, closure, surface)?);

    let mut out = String::from("{\"michell\":\"scene\",\"version\":1,\"meta\":{");
    let _ = write!(
        out,
        "\"speed\":{},\"froude\":{},\"nu\":{},\"transverse_wavelength\":{},\"closure\":{:?},\"frame\":\"x forward, y transverse, z up from the waterline [m]\"",
        num(cond.speed),
        num(cond.speed / (cond.gravity * l_ref).sqrt()),
        num(nu),
        num(2.0 * std::f64::consts::PI / nu),
        format!("{closure:?}")
    );
    out.push_str("},\"objects\":[");
    for (i, o) in objects.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        write_object(&mut out, o);
    }
    out.push_str("]}\n");
    Ok(out)
}

/// The station curves as the physics integrates them, both sides, each
/// vertex carrying its station's section half-area, depth integral and
/// source strength at the transverse wave's decay rate κ = ν.
fn stations(h: &SceneHull, hull: &SectionalHull, pl: Placement, nu: f64) -> Object {
    let (area, _) = hull.depth_integral_curve(0.0, 1);
    let (z_nu, _) = hull.depth_integral_curve(nu, 1);
    let slope = hull.depth_integral_slope_curve(nu, 8);
    let at = |curve: &[(f64, f64)], x: f64| -> f64 {
        // Linear interpolation on a dense (x, value) curve.
        let j = curve.partition_point(|p| p.0 < x).clamp(1, curve.len().max(2) - 1);
        if curve.len() < 2 {
            return curve.first().map_or(0.0, |p| p.1);
        }
        let (a, b) = (curve[j - 1], curve[j]);
        let t = if b.0 > a.0 { (x - a.0) / (b.0 - a.0) } else { 0.0 };
        a.1 + t * (b.1 - a.1)
    };
    let (mut vertices, mut cells) = (Vec::new(), Vec::new());
    let (mut qa, mut qz, mut qs) = (Vec::new(), Vec::new(), Vec::new());
    for (i, (x, c)) in hull.curves().enumerate() {
        for side in [1.0, -1.0] {
            let base = vertices.len();
            for &(yb, zd) in c {
                vertices.push([x + pl.x, pl.y + side * yb, -zd]);
                qa.push(area[i].1);
                qz.push(z_nu[i].1);
                qs.push(at(&slope, x));
            }
            for k in 1..c.len() {
                cells.push(vec![base + k - 1, base + k]);
            }
        }
    }
    Object {
        name: format!("{} stations", h.name),
        kind: "curves",
        vertices,
        cells,
        quantities: vec![
            ("section half-area Z(x;0) [m^2]".into(), qa),
            ("depth integral Z(x;nu) [m^2]".into(), qz),
            ("source strength dZ/dx(x;nu) [m]".into(), qs),
        ],
        note: Some("each station's section curve, as the depth integral integrates it".into()),
    }
}

/// The transom closure's virtual appendage: the transom section with its
/// half-beam faded by φ(s) = 1 − 3s² + 2s³ over the hollow length.
fn appendage(
    h: &SceneHull,
    hull: &SectionalHull,
    pl: Placement,
    nu: f64,
    closure: TransomClosure,
) -> Option<Object> {
    let t = hull.transom()?;
    let lv = closure.hollow_length(t.depth, nu)?;
    if lv <= 0.0 {
        return None;
    }
    let (_, sec) = hull
        .curves()
        .min_by(|a, b| (a.0 - t.x).abs().total_cmp(&(b.0 - t.x).abs()))?;
    let phi = |s: f64| 1.0 - 3.0 * s * s + 2.0 * s * s * s;
    let ns = 25;
    let m = sec.len();
    let (mut vertices, mut cells, mut fade) = (Vec::new(), Vec::new(), Vec::new());
    for side in [1.0, -1.0] {
        let base = vertices.len();
        for i in 0..ns {
            let s = i as f64 / (ns - 1) as f64;
            for &(yb, zd) in sec {
                vertices.push([t.x - s * lv + pl.x, pl.y + side * yb * phi(s), -zd]);
                fade.push(phi(s));
            }
        }
        for i in 0..ns - 1 {
            for j in 0..m - 1 {
                let a = base + i * m + j;
                cells.push(vec![a, a + m, a + m + 1]);
                cells.push(vec![a, a + m + 1, a + 1]);
            }
        }
    }
    Some(Object {
        name: format!("{} transom closure", h.name),
        kind: "mesh",
        vertices,
        cells,
        quantities: vec![("phi(s)".into(), fade)],
        note: Some(format!(
            "virtual appendage, hollow length {lv:.4} m ({closure:?})"
        )),
    })
}

/// The fleet's far-field free-wave elevation as a height-field mesh.
fn free_surface(
    hulls: &[SceneHull],
    cond: &Conditions,
    closure: TransomClosure,
    s: &Surface,
) -> Result<Object, String> {
    let members: Vec<(&SectionalHull, Placement)> =
        hulls.iter().map(|h| (&h.import.hull, h.placement)).collect();
    let mut spec =
        FreeWaveSpectrum::new_sectional(&members, cond, closure).map_err(|e| e.to_string())?;
    let g = spec
        .elevation_grid(s.x0, s.x1, s.y0, s.y1, s.nx, s.ny)
        .map_err(|e| e.to_string())?;
    let mut vertices = Vec::with_capacity(g.nx * g.ny);
    for iy in 0..g.ny {
        for ix in 0..g.nx {
            vertices.push([g.x(ix), g.y(iy), g.get(ix, iy)]);
        }
    }
    let mut cells = Vec::new();
    for iy in 0..g.ny.saturating_sub(1) {
        for ix in 0..g.nx.saturating_sub(1) {
            let a = iy * g.nx + ix;
            cells.push(vec![a, a + 1, a + g.nx + 1]);
            cells.push(vec![a, a + g.nx + 1, a + g.nx]);
        }
    }
    Ok(Object {
        name: "free surface".into(),
        kind: "mesh",
        vertices,
        cells,
        quantities: vec![("wave elevation zeta [m]".into(), g.zeta.clone())],
        note: Some(format!(
            "far-field free waves only: the full linear answer aft of the sterns, not \
             abreast of or ahead of the hulls{}",
            if g.resolution_limited {
                "; waves shorter than the grid resolves were tapered away"
            } else {
                ""
            }
        )),
    })
}

/// Seven significant figures: far below what any viewer resolves, and a
/// third of the size of full round-trip precision.
fn num(v: f64) -> String {
    if v == 0.0 {
        "0".into()
    } else if v.is_finite() {
        let s = format!("{v:.6e}");
        let (m, e) = s.split_once('e').unwrap_or((&s, "0"));
        let m = if m.contains('.') { m.trim_end_matches('0').trim_end_matches('.') } else { m };
        if e == "0" { m.to_string() } else { format!("{m}e{e}") }
    } else {
        "null".into()
    }
}

fn write_object(out: &mut String, o: &Object) {
    let _ = write!(out, "{{\"name\":{:?},\"kind\":{:?},\"vertices\":[", o.name, o.kind);
    for (i, v) in o.vertices.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(out, "[{},{},{}]", num(v[0]), num(v[1]), num(v[2]));
    }
    out.push_str(if o.kind == "curves" { "],\"edges\":[" } else { "],\"faces\":[" });
    for (i, c) in o.cells.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('[');
        for (k, v) in c.iter().enumerate() {
            if k > 0 {
                out.push(',');
            }
            let _ = write!(out, "{v}");
        }
        out.push(']');
    }
    out.push_str("],\"quantities\":[");
    for (i, (name, vals)) in o.quantities.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(out, "{{\"name\":{name:?},\"on\":\"vertices\",\"values\":[");
        for (k, v) in vals.iter().enumerate() {
            if k > 0 {
                out.push(',');
            }
            out.push_str(&num(*v));
        }
        out.push_str("]}");
    }
    out.push(']');
    if let Some(n) = &o.note {
        let _ = write!(out, ",\"note\":{n:?}");
    }
    out.push('}');
}

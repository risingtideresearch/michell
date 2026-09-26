//! `michell-web` — a browser front end for the `michell` tools.
//!
//! Upload a hull and see it the way the physics sees it: cut into sections.
//! IGES and STL hulls are cut straight from their patches or triangles, and
//! exact B-spline `.hull` files from their exact surfaces — the same loader
//! the CLI uses ([`michell_cli::fleet`]). The page is sent each
//! station's section curve (what the depth integral integrates), the CAD
//! ray hits it was interpolated from, the
//! depth-integral curve the kernel interpolates along x, the hydrostatics,
//! and the transom with what the closure needs to draw its virtual appendage
//! at any speed.

use michell::iges::{HullPose, Platform, SectionalImport};
use michell::sectional::SectionalHull;
use michell::Placement;
use michell_cli::fleet::{open_source_bytes, Kind, LoadSettings};
use serde_json::{json, Value};

/// Largest upload accepted [bytes]. Big enough for a finely tessellated STL.
pub const MAX_UPLOAD: usize = 128 << 20;

/// Import options the page can set.
#[derive(Default)]
pub struct LoftRequest {
    /// Design waterline height in the file's frame [m] (IGES).
    pub waterline: Option<f64>,
    /// Centreplane override [m] (IGES).
    pub centerplane: Option<f64>,
    /// Stations along the hull (IGES).
    pub stations: Option<usize>,
    /// Rays across each section (IGES).
    pub rays: Option<usize>,
    /// Scale to metres (STL, which carries no units): mm, m, in, ... or a
    /// number.
    pub units: Option<f64>,
}

impl LoftRequest {
    /// Parse from `key=value` query pairs; unknown keys are ignored.
    pub fn from_query(pairs: &[(String, String)]) -> Result<LoftRequest, String> {
        let mut r = LoftRequest::default();
        let num = |k: &str, v: &str| {
            v.trim()
                .parse::<f64>()
                .map_err(|_| format!("{k}: expected a number, got {v:?}"))
        };
        let count = |k: &str, v: &str, min: usize| match v.trim().parse::<usize>() {
            Ok(n) if n >= min => Ok(n),
            _ => Err(format!(
                "{k}: expected a count of at least {min}, got {v:?}"
            )),
        };
        for (k, v) in pairs {
            if v.trim().is_empty() {
                continue;
            }
            match k.as_str() {
                "waterline" => r.waterline = Some(num(k, v)?),
                "centerplane" => r.centerplane = Some(num(k, v)?),
                "stations" => r.stations = Some(count(k, v, 8)?),
                "rays" => r.rays = Some(count(k, v, 5)?),
                "units" => r.units = Some(michell_cli::parse_units(v.trim())?),
                _ => {}
            }
        }
        Ok(r)
    }
}

/// Cut an uploaded hull into sections and describe it as the JSON the page
/// draws.
pub fn loft(name: &str, bytes: Vec<u8>, req: &LoftRequest) -> Result<Value, String> {
    let t0 = std::time::Instant::now();
    let d = LoadSettings::default();
    let settings = LoadSettings {
        waterline_z: req.waterline.unwrap_or(0.0),
        centerplane: req.centerplane,
        stations: req.stations.unwrap_or(d.stations),
        rays: req.rays.unwrap_or(d.rays),
        units: req.units,
    };
    let file = open_source_bytes(name, bytes, &settings)?;
    let opts = settings.sectional(file.waterline_z);
    let mut notes = Vec::new();
    let mut hulls = Vec::new();
    for i in 0..file.source.len() {
        match file.source.situate_sectional(
            i,
            file.waterline_z,
            &HullPose::default(),
            &Platform::default(),
            &opts,
        ) {
            Ok(Some(h)) => hulls.push(sectioned_json(&h, file.kind)),
            Ok(None) => notes.push(format!("hull {} is dry at this waterline", i + 1)),
            Err(e) => notes.push(format!("hull {} not sectioned: {e}", i + 1)),
        }
    }
    if hulls.is_empty() {
        return Err(notes
            .first()
            .cloned()
            .unwrap_or_else(|| "no hull found".into()));
    }
    Ok(json!({
        "name": name,
        "seconds": t0.elapsed().as_secs_f64(),
        "notes": notes,
        "hulls": hulls,
    }))
}

fn sectioned_json(imp: &SectionalImport, kind: Kind) -> Value {
    let r = &imp.report;
    let sides = if r.two_sided {
        format!("both sides averaged about y = {:.4} m", r.centerplane)
    } else {
        format!("one side about y = {:.4} m", r.centerplane)
    };
    let what = match kind {
        Kind::Iges => format!("IGES, {} patches", r.patches),
        Kind::Stl => format!("STL, {} triangles", r.patches),
        Kind::Spline => "exact B-spline control net".to_string(),
    };
    let mut lines = vec![format!(
        "source: {what}, cut at {} stations over x {:.4}..{:.4} m ({sides})",
        r.stations, r.x_range.0, r.x_range.1
    )];
    if r.dropped_stations > 0 {
        lines.push(format!(
            "{} interior stations could not be sectioned (bridged by the interpolant)",
            r.dropped_stations
        ));
    }
    if r.max_asymmetry > 1e-3 * r.draft.max(1e-9) {
        lines.push(format!(
            "port and starboard differ by up to {:.3e} m (averaged)",
            r.max_asymmetry
        ));
    }
    hull_json(&imp.hull, imp.placement, Some(&imp.sections), lines)
}

/// One hull, drawn as the physics uses it (see the module docs).
fn hull_json(
    hull: &SectionalHull,
    placement: Placement,
    rays: Option<&Vec<(f64, Vec<(f64, f64)>)>>,
    lines: Vec<String>,
) -> Value {
    // Each station as the curve the quadrature integrates (not a polyline
    // through its quadrature nodes, which are graded toward the waterline).
    let stations: Vec<(f64, Vec<(f64, f64)>)> =
        hull.curves().map(|(x, o)| (x, o.to_vec())).collect();
    // A see-through surface between stations, for orientation only (the
    // kernel interpolates each station's depth integral along x, not a
    // surface): rows joined at equal fractions of girth, so neighbouring
    // sections with different node spacing still meet cleanly.
    let rows = 48usize;
    let (mut mx, mut mz, mut my) = (Vec::new(), Vec::new(), Vec::new());
    for (x, o) in &stations {
        for (y, z) in by_girth(o, rows) {
            mx.push(*x);
            mz.push(z);
            my.push(y);
        }
    }
    let keel: Vec<(f64, f64)> = stations
        .iter()
        .map(|(x, o)| (*x, o.last().map_or(0.0, |p| p.1)))
        .collect();
    let beam = 2.0
        * stations
            .iter()
            .flat_map(|(_, o)| o.first().map(|p| p.0))
            .fold(0.0, f64::max);
    // The depth-integral curves at κ = 0 (sectional area) and at a short
    // wave's decay rate: λ = 2 at Fn 0.15, κ = νλ², ν = g/U² = 1/(0.0225 L).
    let kappas = [0.0, 4.0 / (0.0225 * hull.length())];
    let curves: Vec<Value> = kappas
        .iter()
        .map(|&kappa| {
            let (st, c) = hull.depth_integral_curve(kappa, 8);
            json!({ "kappa": kappa, "stations": st, "curve": c })
        })
        .collect();
    // The transom and its section — the aft end station's outline, whose
    // depth integral is the closure's depth factor — so the page can draw
    // the virtual appendage at whatever speed and closure it is asked about.
    let transom = hull.transom().map(|t| {
        let aft = stations
            .iter()
            .min_by(|a, b| (a.0 - t.x).abs().total_cmp(&(b.0 - t.x).abs()));
        json!({
            "x": t.x,
            "depth": t.depth,
            "half_beam": t.half_beam,
            "area": t.area,
            "area_ratio": t.area / hull.max_section_area().max(1e-300),
            "outline": aft.map(|s| by_girth(&s.1, rows)).unwrap_or_default(),
            // Z_T at each curve's κ: the kernel extends Z over the hollow as
            // Z_T · φ(s).
            "z_t": kappas
                .iter()
                .map(|&k| hull.depth_integral_curve(k, 1).0.first().map_or(0.0, |p| p.1))
                .collect::<Vec<_>>(),
        })
    });
    json!({
        "placement": { "x": placement.x, "y": placement.y },
        "diagnostics": lines,
        "length": hull.length(),
        "beam": beam,
        "draft": hull.draft(),
        "displaced_volume": hull.displaced_volume(),
        "wetted_surface": hull.wetted_surface(),
        "lcb_x": hull.lcb_x(),
        "waterplane_area": hull.waterplane_area(),
        "stations": stations,
        "rays": rays,
        "mesh": { "nx": stations.len(), "nz": rows + 1, "x": mx, "z": mz, "y": my },
        "keel": keel,
        "area": curves,
        "transom": transom,
    })
}

/// A section outline resampled at `rows + 1` points evenly spaced in girth
/// (a point for an empty section).
fn by_girth(o: &[(f64, f64)], rows: usize) -> Vec<(f64, f64)> {
    if o.len() < 2 {
        return vec![(0.0, 0.0); rows + 1];
    }
    let mut cum = vec![0.0];
    for w in o.windows(2) {
        let d = ((w[1].0 - w[0].0).powi(2) + (w[1].1 - w[0].1).powi(2)).sqrt();
        cum.push(cum.last().unwrap() + d);
    }
    let total = *cum.last().unwrap();
    (0..=rows)
        .map(|k| {
            let g = total * k as f64 / rows as f64;
            let j = cum.partition_point(|&c| c < g).clamp(1, o.len() - 1);
            let t = if cum[j] > cum[j - 1] {
                (g - cum[j - 1]) / (cum[j] - cum[j - 1])
            } else {
                0.0
            };
            (
                o[j - 1].0 + t * (o[j].0 - o[j - 1].0),
                o[j - 1].1 + t * (o[j].1 - o[j - 1].1),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The Wigley as a native `.hull` file, shown by sections.
    fn wigley_file(hull: &michell::Hull) -> String {
        let s = hull.surface();
        let join = |v: &[f64]| {
            v.iter()
                .map(|x| format!("{x}"))
                .collect::<Vec<_>>()
                .join(" ")
        };
        let mut out = format!(
            "michell-hull v1\ndegree-x {}\ndegree-z {}\nknots-x {}\nknots-z {}\n",
            s.degree_x(),
            s.degree_z(),
            join(s.knots_x()),
            join(s.knots_z())
        );
        let nz = s.n_ctrl_z();
        for i in 0..s.n_ctrl_x() {
            out.push_str(&format!(
                "row {}\n",
                join(&s.control()[i * nz..(i + 1) * nz])
            ));
        }
        out
    }

    #[test]
    fn a_native_wigley_is_shown_in_sections() {
        let hull = michell::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let v = loft(
            "w.hull",
            wigley_file(&hull).into_bytes(),
            &LoftRequest::default(),
        )
        .unwrap();
        let h = &v["hulls"][0];
        let vol = h["displaced_volume"].as_f64().unwrap();
        // The Wigley's exact volume, 4/9 L B T.
        let exact = 4.0 / 9.0 * 10.0 * 0.625;
        assert!((vol - exact).abs() < 1e-8 * exact, "{vol}");
        assert!(h["transom"].is_null());
        assert!(!h["stations"].as_array().unwrap().is_empty());
    }

    #[test]
    fn an_stl_needs_its_units() {
        let e = loft("x.stl", vec![0; 200], &LoftRequest::default()).unwrap_err();
        assert!(e.contains("units"), "{e}");
    }
}

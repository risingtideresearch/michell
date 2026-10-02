//! A calm-water result (or a prop on one) as an IGES file for CAD: the
//! hulls at their solved attitude, the free surface, and with a prop its
//! discs. Each group is on its own level, coloured and labelled:
//!
//! | level | colour | label | what |
//! |---|---|---|---|
//! | 1 | white | `HULL1`, `HULL2` | each hull's patches |
//! | 2 | cyan | `WATER` | the free surface |
//! | 3 | red | `PROP1`, … | each propeller disc |
//!
//! Frame: x forward, y to port, z up, the still water at z = 0, metres.

use crate::records::{expect, hull_source};
use crate::stream::{short, Stream};
use boatmath::params::CaseParams;
use boatmath::propulsion::Position;
use michell::propulsion::Disc;
use michell_geometry::iges::{write_labelled, Color, Label, NurbsSurface3};
use serde_json::Value;

pub struct Options {
    pub water: bool,
    pub wave_scale: f64,
    /// Where a prop without a placement of its own puts its discs: forward
    /// of each hull's aft end, out from its centreplane, and down (default
    /// the aft end, on the centreplane, at the prop's shaft depth).
    pub prop_x: Option<f64>,
    pub prop_y: Option<f64>,
    pub depth: Option<f64>,
}

/// The IGES file for one record: a calm-water result, or a prop.
pub fn model(st: &Stream, rec: &Value, o: &Options) -> Result<String, String> {
    let (result, prop) = match rec["type"].as_str() {
        Some("prop") => {
            let id = rec["result"]
                .as_str()
                .ok_or("a prop given outright (no result) has no hull to draw")?;
            (st.need("result", id)?.clone(), Some(rec))
        }
        _ => (rec.clone(), None),
    };
    let rid = expect(&result, "result")?;
    if result["kind"] != "calm" {
        return Err(format!("result {}: not a calm-water result", short(rid)));
    }
    let study = st.follow(&result, "study")?;
    let case = st.follow(study, "case")?;
    let hull = st.follow(case, "hull")?;
    let src = hull_source(hull)?;
    let cp: CaseParams = crate::records::case_params(case)?;
    let f = &result["forces"];
    let attitude = (
        f["sinkage"].as_f64().unwrap_or(0.0),
        f["trim_rad"].as_f64().unwrap_or(0.0),
    );

    let mut surfaces: Vec<NurbsSurface3> = Vec::new();
    let mut labels: Vec<Label> = Vec::new();
    let hulls = boatmath::cad::hull_surfaces(
        &src.file_name,
        src.bytes.clone(),
        &src.import,
        &cp,
        attitude,
    )?;
    for (k, patches) in hulls.iter().enumerate() {
        for p in patches {
            surfaces.push(p.clone());
            labels.push(Label {
                level: 1,
                color: Color::White,
                name: format!("HULL{}", k + 1),
            });
        }
    }

    if o.water {
        let v = st.follow(&result, "field")?;
        let s = &v["surface"];
        let n = |k: &str| s[k].as_f64().ok_or_else(|| format!("field: surface {k}"));
        let zeta: Vec<f64> = crate::views::f32s(&s["zeta"], "zeta")?
            .into_iter()
            .map(f64::from)
            .collect();
        surfaces.push(boatmath::cad::water_surface(
            n("x0")?,
            n("x1")?,
            n("nx")? as usize,
            n("y0")?,
            n("y1")?,
            n("ny")? as usize,
            &zeta,
            o.wave_scale,
        )?);
        labels.push(Label {
            level: 2,
            color: Color::Cyan,
            name: "WATER".into(),
        });
    }

    if let Some(p) = prop {
        for (k, d) in discs(&src, &cp, attitude, p, o)?.iter().enumerate() {
            surfaces.push(boatmath::cad::disc_surface(d));
            labels.push(Label {
                level: 3,
                color: Color::Red,
                name: format!("PROP{}", k + 1),
            });
        }
    }

    let product = format!(
        "{} Fn {}",
        hull["name"].as_str().unwrap_or("hull"),
        result["froude"].as_f64().unwrap_or(0.0)
    );
    write_labelled(&surfaces, &labels, &product).map_err(|e| e.to_string())
}

/// A prop's discs: where its own wake estimate placed them, or at `at`.
fn discs(
    src: &crate::records::HullSource,
    cp: &CaseParams,
    attitude: (f64, f64),
    p: &Value,
    o: &Options,
) -> Result<Vec<Disc>, String> {
    let radius = 0.5 * p["best"]["D"].as_f64().ok_or("the prop has no propeller")?;
    if let Some(ds) = p["interaction"]["discs"].as_array() {
        return Ok(ds
            .iter()
            .map(|d| Disc {
                x: d["x"].as_f64().unwrap_or(0.0),
                y: d["y"].as_f64().unwrap_or(0.0),
                depth: d["depth"].as_f64().unwrap_or(0.0),
                radius: d["radius"].as_f64().unwrap_or(radius),
                hub: 0.2,
            })
            .collect());
    }
    let at = Position {
        forward_of_aft: o.prop_x.unwrap_or(0.0),
        outboard: o.prop_y.unwrap_or(0.0),
        depth: o.depth.or(p["inputs"]["depth"].as_f64()).unwrap_or(0.3),
    };
    let m =
        boatmath::propulsion::mounts(&src.file_name, src.bytes.clone(), &src.import, cp, attitude)?;
    let shafts = p["inputs"]["shafts"].as_u64().unwrap_or(1) as usize;
    let per_hull = (shafts / m.hulls.len().max(1)).max(1);
    Ok(m.discs(&at, per_hull, radius))
}

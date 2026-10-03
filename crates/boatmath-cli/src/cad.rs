//! A calm-water result (or a prop on one) as an IGES file for CAD: the
//! hulls at their solved attitude with their drives, the free surface, and
//! with a prop its discs. Each group is on its own level, coloured and
//! labelled:
//!
//! | level | colour | label | what |
//! |---|---|---|---|
//! | 1 | white | `HULL1`, `HULL2` | each hull's patches |
//! | 2 | cyan | `WATER` | the free surface |
//! | 3 | red | `PROP1`, … | each propeller: its B-series blades and hub (or a disc) |
//! | 4 | yellow | `MOUNT1`, … | each drive's leg and pod |
//!
//! Frame: x forward, y to port, z up, the still water at z = 0, metres.

use crate::records::{case_source, expect, CaseSource};
use crate::stream::{short, Stream};
use boatmath::params::CaseParams;
use hullgeom::iges::{write_labelled, Color, Label, NurbsSurface3};
use serde_json::Value;
use thinship::propulsion::Disc;

pub struct Options {
    pub water: bool,
    pub wave_scale: f64,
    /// A prop's propellers as their blades (else plain discs).
    pub blades: bool,
    pub right_handed: bool,
}

/// A member's level, colour and label: hulls on 1, a drive's parts on 4.
fn member_label(role: &str, hull_no: usize, part_no: usize) -> Label {
    if role == "hull" {
        Label {
            level: 1,
            color: Color::White,
            name: format!("HULL{hull_no}"),
        }
    } else {
        Label {
            level: 4,
            color: Color::Yellow,
            name: format!("MOUNT{part_no}"),
        }
    }
}

/// Number the members: hulls in order, and parts in order.
fn labelled<T>(members: &[(String, T)]) -> Vec<Label> {
    let (mut h, mut p) = (0, 0);
    members
        .iter()
        .map(|(role, _)| {
            if role == "hull" {
                h += 1;
            } else {
                p += 1;
            }
            member_label(role, h, p)
        })
        .collect()
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
    let CaseSource {
        hull,
        params: cp,
        src,
        thrusts,
        ..
    } = case_source(st, case)?;
    let f = &result["forces"];
    let attitude = (
        f["sinkage"].as_f64().unwrap_or(0.0),
        f["trim_rad"].as_f64().unwrap_or(0.0),
    );

    let mut surfaces: Vec<NurbsSurface3> = Vec::new();
    let mut labels: Vec<Label> = Vec::new();
    let members = boatmath::cad::hull_surfaces(
        &src.file_name,
        src.bytes.clone(),
        &src.import,
        &cp,
        attitude,
    )?;
    for ((_, patches), label) in members.iter().zip(labelled(&members)) {
        for p in patches {
            surfaces.push(p.clone());
            labels.push(label.clone());
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
        let ds = discs(&src, &cp, attitude, &thrusts, p)?;
        if ds.is_empty() {
            eprintln!(
                "boatmath: prop {}: its case has no mount, so no disc is drawn (`boatmath mount`)",
                short(p["id"].as_str().unwrap_or(""))
            );
        }
        // The shaft's angle to the still water: its own and the trim (the
        // drive's thrust line keeps both), else the trim.
        let angle = p["thrust_line"]["angle_deg"]
            .as_f64()
            .map_or(attitude.1, f64::to_radians);
        let b = &p["best"];
        let chosen = match (
            b["Z"].as_f64(),
            b["D"].as_f64(),
            b["PD"].as_f64(),
            b["EAR"].as_f64(),
        ) {
            (Some(z), Some(dia), Some(pd), Some(ear)) if o.blades => {
                Some(boatmath::cad::Propeller {
                    blades: z as usize,
                    diameter: dia,
                    pitch_ratio: pd,
                    area_ratio: ear,
                    right_handed: o.right_handed,
                })
            }
            _ => None,
        };
        for (k, d) in ds.iter().enumerate() {
            let label = Label {
                level: 3,
                color: Color::Red,
                name: format!("PROP{}", k + 1),
            };
            let parts = match &chosen {
                Some(prop) => boatmath::cad::propeller_surfaces(d, angle, prop)?,
                None => vec![boatmath::cad::disc_surface(d)],
            };
            for s in parts {
                surfaces.push(s);
                labels.push(label.clone());
            }
        }
    }

    let product = format!(
        "{} Fn {}",
        hull["name"].as_str().unwrap_or("hull"),
        result["froude"].as_f64().unwrap_or(0.0)
    );
    write_labelled(&surfaces, &labels, &product).map_err(|e| e.to_string())
}

/// A prop's discs: where its wake estimate placed them, else where the
/// case's mount puts its propellers; none without a mount.
fn discs(
    src: &crate::records::HullSource,
    cp: &CaseParams,
    attitude: (f64, f64),
    thrusts: &[boatmath::mount::Thrust],
    p: &Value,
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
                hub: d["hub"].as_f64().unwrap_or(0.2),
            })
            .collect());
    }
    if thrusts.is_empty() {
        return Ok(Vec::new());
    }
    let m =
        boatmath::propulsion::mounts(&src.file_name, src.bytes.clone(), &src.import, cp, attitude)?;
    Ok(m.mount_discs(thrusts, radius))
}

// ------------------------------------------------------------- statics

/// A pose of the hull: heel φ about x (lifting +y), trim θ bow up about
/// x = 0, then down by the sinkage — from the hull's own axes (z up from
/// the design waterline) to the water's.
#[derive(Clone, Copy)]
struct Pose {
    heel: f64,
    trim: f64,
    sinkage: f64,
}

impl Pose {
    fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        let (sf, cf) = self.heel.sin_cos();
        let (st, ct) = self.trim.sin_cos();
        let (y1, z1) = (p[1] * cf - p[2] * sf, p[1] * sf + p[2] * cf);
        [p[0] * ct - z1 * st, y1, p[0] * st + z1 * ct - self.sinkage]
    }
}

/// A flat patch of still water over `[x0, x1] × [y0, y1]`.
fn water_plane(x0: f64, x1: f64, y0: f64, y1: f64) -> NurbsSurface3 {
    NurbsSurface3 {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        n_ctrl_u: 2,
        n_ctrl_v: 2,
        ctrl: vec![[x0, y0, 0.0], [x0, y1, 0.0], [x1, y0, 0.0], [x1, y1, 0.0]],
        weights: vec![1.0; 4],
        trim_uv: None,
    }
}

fn point3(v: &Value) -> Option<[f64; 3]> {
    let a = v.as_array()?;
    Some([
        a.first()?.as_f64()?,
        a.get(1)?.as_f64()?,
        a.get(2)?.as_f64()?,
    ])
}

/// A statics record as IGES: a case floated at rest, at its maximum GZ, at
/// its angle of vanishing stability and at each of `heels` [deg], each pose
/// on its own level with its G, B and righting arm; a hull's statics, the
/// hull upright at its design waterline. The still water on level 2.
pub fn statics_model(st: &Stream, rec: &Value, heels: &[f64]) -> Result<String, String> {
    use hullgeom::iges::{write_entities, Entity};
    let id = expect(rec, "statics")?;
    if rec.get("error").is_some_and(|e| !e.is_null()) {
        return Err(format!("statics {} failed: {}", short(id), rec["error"]));
    }
    let (hull, cp, is_case, src) = match rec.get("case").and_then(|c| c.as_str()) {
        Some(_) => {
            let case = st.follow(rec, "case")?;
            let cs = case_source(st, case)?;
            (cs.hull, cs.params, true, cs.src)
        }
        None => {
            let h = st.follow(rec, "hull")?;
            (
                h,
                CaseParams::default(),
                false,
                crate::records::hull_source(h)?,
            )
        }
    };
    // The hulls in their own axes: the design pose, a catamaran's demihulls
    // at their span.
    let design = boatmath::cad::hull_surfaces(
        &src.file_name,
        src.bytes.clone(),
        &src.import,
        &cp,
        (0.0, 0.0),
    )?;
    let label = |level: u32, color: Color, name: &str| Label {
        level,
        color,
        name: name.to_string(),
    };
    let mut items: Vec<(Entity, Label)> = Vec::new();
    if is_case {
        // The poses: (level, name, heel [deg]).
        let gz = &rec["gz"];
        let mut poses: Vec<(u32, String, f64)> = vec![(1, "REST".into(), 0.0)];
        if let Some(h) = gz["heel_at_max_deg"].as_f64() {
            poses.push((11, "GZMAX".into(), h));
        }
        // A vanishing angle of 0 is a platform with no positive stability;
        // of 180, one that never loses it: neither is a pose of its own.
        if let Some(h) = gz["vanishing_deg"]
            .as_f64()
            .filter(|h| *h > 0.5 && *h < 179.9)
        {
            poses.push((12, "VANISH".into(), h));
        }
        for (k, &h) in heels.iter().enumerate() {
            poses.push((13 + k as u32, format!("HEEL{}", h.round()), h));
        }
        let want: Vec<f64> = poses.iter().map(|p| p.2).collect();
        let floated = boatmath::heeled(&src.file_name, src.bytes.clone(), &src.import, &cp, &want)?;
        let g = point3(&floated["g"]).ok_or("no centre of gravity")?;
        let fl = floated["poses"].as_array().ok_or("no poses")?;
        if fl.len() != poses.len() {
            return Err(format!(
                "{} poses asked for, {} floated",
                poses.len(),
                fl.len()
            ));
        }
        for ((level, name, _), f) in poses.iter().zip(fl) {
            if let Some(e) = f["error"].as_str() {
                eprintln!(
                    "boatmath: statics {}: {name} at {}°: {e}; left out",
                    short(id),
                    f["asked_deg"]
                );
                continue;
            }
            let (asked, got) = (
                f["asked_deg"].as_f64().unwrap_or(0.0),
                f["heel_deg"].as_f64().unwrap_or(0.0),
            );
            if (asked - got).abs() > 1e-6 {
                eprintln!(
                    "boatmath: statics {}: {name} drawn at {got:.2}°, just short of {asked:.2}° where the float fails",
                    short(id)
                );
            }
            let pose = Pose {
                heel: f["heel_deg"].as_f64().unwrap_or(0.0).to_radians(),
                trim: f["trim_rad"].as_f64().unwrap_or(0.0),
                sinkage: f["sinkage"].as_f64().unwrap_or(0.0),
            };
            for (role, patches) in &design {
                let color = if role == "hull" {
                    Color::White
                } else {
                    Color::Yellow
                };
                for p in patches {
                    let mut q = p.clone();
                    for c in q.ctrl.iter_mut() {
                        *c = pose.apply(*c);
                    }
                    items.push((Entity::Surface(q), label(*level, color, name)));
                }
            }
            let cb = point3(&f["cb"]).ok_or("no centre of buoyancy")?;
            let (gw, bw) = (pose.apply(g), pose.apply(cb));
            items.push((
                Entity::Point(gw),
                label(*level, Color::Red, &format!("{name}-G")),
            ));
            items.push((
                Entity::Point(bw),
                label(*level, Color::Blue, &format!("{name}-B")),
            ));
            // The righting arm: from G across to the vertical through B.
            items.push((
                Entity::Line(gw, [gw[0], bw[1], gw[2]]),
                label(*level, Color::Yellow, &format!("{name}-A")),
            ));
        }
    } else {
        for ((_, patches), l) in design.iter().zip(labelled(&design)) {
            for p in patches {
                items.push((Entity::Surface(p.clone()), l.clone()));
            }
        }
    }
    let (mut xa, mut xb, mut yb) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
    for c in design
        .iter()
        .flat_map(|(_, ps)| ps.iter())
        .flat_map(|p| p.ctrl.iter())
    {
        xa = xa.min(c[0]);
        xb = xb.max(c[0]);
        yb = yb.max(c[1].abs());
    }
    let (l, pad) = (xb - xa, 0.15 * (xb - xa));
    items.push((
        Entity::Surface(water_plane(
            xa - pad,
            xb + pad,
            -(yb + 0.3 * l),
            yb + 0.3 * l,
        )),
        label(2, Color::Cyan, "WATER"),
    ));
    let product = format!("{} statics", hull["name"].as_str().unwrap_or("hull"));
    write_entities(&items, &product).map_err(|e| e.to_string())
}

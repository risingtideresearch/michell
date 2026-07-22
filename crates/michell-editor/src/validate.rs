//! A lightweight validator that mirrors the constraints `manifest.rs` enforces
//! at run time, so the editor can flag a broken study before it is saved. It is
//! deliberately geometry-free: anything that needs the loaded `.hull` bodies
//! (e.g. a `spread` axis on a centreline hull, or the wetted reference length)
//! is left to the CLI.

use crate::model::*;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Level {
    Error,
    Warning,
}

pub struct Issue {
    pub level: Level,
    pub msg: String,
}

fn err(msg: impl Into<String>) -> Issue {
    Issue {
        level: Level::Error,
        msg: msg.into(),
    }
}

fn warn(msg: impl Into<String>) -> Issue {
    Issue {
        level: Level::Warning,
        msg: msg.into(),
    }
}

/// The number of values an axis expands to, or `None` if its spec cannot be
/// parsed yet (the value fields are surfaced separately, so we stay quiet).
fn axis_count(v: &ValueSpec) -> Option<usize> {
    match v.mode {
        ValueMode::Scalar => v.scalar.trim().parse::<f64>().ok().map(|_| 1),
        ValueMode::List => {
            let items: Vec<&str> = v
                .list
                .split([',', ' ', '\t', '\n'])
                .filter(|t| !t.trim().is_empty())
                .collect();
            if items.iter().all(|t| t.trim().parse::<f64>().is_ok()) && !items.is_empty() {
                Some(items.len())
            } else {
                None
            }
        }
        ValueMode::Range => {
            let a = v.range_start.trim().parse::<f64>().ok()?;
            let b = v.range_stop.trim().parse::<f64>().ok()?;
            if b < a {
                return None;
            }
            if b == a {
                return Some(1);
            }
            let step = if v.step.trim().is_empty() {
                (b - a) / 5.0
            } else {
                let s = v.step.trim().parse::<f64>().ok()?;
                if s <= 0.0 {
                    return None;
                }
                s
            };
            Some(((b - a) / step + 1e-9).floor() as usize + 1)
        }
    }
}

/// Parsed axis values, for range/scale checks (empty if unparseable).
fn axis_values(v: &ValueSpec) -> Vec<f64> {
    match v.mode {
        ValueMode::Scalar => v.scalar.trim().parse().into_iter().collect(),
        ValueMode::List => v
            .list
            .split([',', ' ', '\t', '\n'])
            .filter_map(|t| t.trim().parse().ok())
            .collect(),
        ValueMode::Range => [&v.range_start, &v.range_stop]
            .iter()
            .filter_map(|s| s.trim().parse().ok())
            .collect(),
    }
}

pub fn validate(m: &Manifest) -> Vec<Issue> {
    let mut out = Vec::new();

    // Hulls.
    if m.hulls.is_empty() {
        out.push(err("manifest has no hulls"));
    }
    for (i, h) in m.hulls.iter().enumerate() {
        let name = if h.id.trim().is_empty() {
            format!("#{}", i + 1)
        } else {
            h.id.clone()
        };
        if h.id.trim().is_empty() {
            out.push(err(format!("hull #{} needs an id", i + 1)));
        }
        if h.file.trim().is_empty() {
            out.push(err(format!("hull {name:?} needs a file")));
        }
        if h.pose.enabled && !(h.pose.scale > 0.0 && h.pose.scale.is_finite()) {
            out.push(err(format!("hull {name:?}: pose scale must be positive")));
        }
        if h.load.enabled && h.load.mass < 0.0 {
            out.push(err(format!("hull {name:?}: load mass must be >= 0")));
        }
        for p in &h.points {
            if p.id.trim().is_empty() {
                out.push(err(format!("hull {name:?}: a point load needs an id")));
            }
            if p.mass < 0.0 {
                out.push(err(format!(
                    "hull {name:?} point {:?}: mass must be >= 0",
                    p.id
                )));
            }
        }
    }

    // Global id uniqueness (hull ids and point ids share one namespace).
    let mut seen: Vec<&str> = Vec::new();
    for h in &m.hulls {
        for id in std::iter::once(&h.id).chain(h.points.iter().map(|p| &p.id)) {
            let id = id.as_str();
            if id.trim().is_empty() {
                continue;
            }
            if seen.contains(&id) {
                out.push(err(format!(
                    "duplicate id {id:?} (hull and point ids must be unique)"
                )));
            } else {
                seen.push(id);
            }
        }
    }

    // Speed axis.
    let speed_count = m.axes.iter().filter(|a| a.kind == AxisKind::Speed).count();
    if speed_count == 0 {
        out.push(err("the sweep needs a speed axis"));
    } else if speed_count > 1 {
        out.push(err("only one speed axis is allowed"));
    }

    // Targets resolve; scale-axis values positive.
    let hull_ids: Vec<&str> = m.hulls.iter().map(|h| h.id.as_str()).collect();
    let point_ids: Vec<String> = m.point_ids();
    for a in &m.axes {
        match a.kind {
            AxisKind::Hull => {
                if a.targets.is_empty() {
                    out.push(err("a hull-param axis needs at least one target hull"));
                }
                for t in &a.targets {
                    if !hull_ids.contains(&t.as_str()) {
                        out.push(err(format!("sweep target {t:?} is not a hull id")));
                    }
                }
                if a.hull_param == HullParam::Scale
                    && axis_values(&a.values)
                        .iter()
                        .any(|&v| !(v > 0.0 && v.is_finite()))
                {
                    out.push(err("scale factors must be positive and finite"));
                }
            }
            AxisKind::Point => {
                if a.targets.is_empty() {
                    out.push(err("a point-load axis needs at least one target point"));
                }
                for t in &a.targets {
                    if !point_ids.iter().any(|p| p == t) {
                        out.push(err(format!("sweep target {t:?} is not a point-load id")));
                    }
                }
            }
            _ => {}
        }
    }

    // Equilibrium (float) mode: the fleet carries mass, either a base load or a
    // swept mass axis. The CG is derived, so lcg/vcg only make sense there.
    let base_mass: f64 = m
        .hulls
        .iter()
        .map(|h| {
            (if h.load.enabled { h.load.mass } else { 0.0 })
                + h.points.iter().map(|p| p.mass).sum::<f64>()
        })
        .sum();
    let mass_axis = m.axes.iter().any(|a| {
        (a.kind == AxisKind::Hull && a.hull_param == HullParam::Mass)
            || (a.kind == AxisKind::Point && a.point_param == PointParam::Mass)
    });
    let float_mode = base_mass > 0.0 || mass_axis;

    let has_waterline = m.axes.iter().any(|a| a.kind == AxisKind::Waterline);
    if float_mode && has_waterline {
        out.push(err(
            "a waterline axis cannot be combined with hull mass (the waterline is solved)",
        ));
    }
    let has_cg_axis = m.axes.iter().any(|a| {
        a.kind == AxisKind::Hull && matches!(a.hull_param, HullParam::Lcg | HullParam::Vcg)
    });
    if has_cg_axis && !float_mode {
        out.push(err(
            "an lcg/vcg axis needs the fleet to carry mass — give a hull a load mass or sweep a mass axis",
        ));
    }

    // Heel roll-up options.
    if m.options.heel.enabled {
        for a in m
            .options
            .heel
            .resistance_angles
            .split([',', ' '])
            .filter(|t| !t.trim().is_empty())
        {
            match a.trim().parse::<f64>() {
                Ok(v) if v.is_finite() && v.abs() < 90.0 => {}
                _ => out.push(err(format!(
                    "options.heel.resistance_angles: {a:?} must be within ±90°"
                ))),
            }
        }
        for (name, f) in [
            ("gz_step", &m.options.heel.gz_step),
            ("gz_max", &m.options.heel.gz_max),
        ] {
            if !f.trim().is_empty() && !f.trim().parse::<f64>().map(|v| v > 0.0).unwrap_or(false) {
                out.push(err(format!(
                    "options.heel.{name} must be a positive number"
                )));
            }
        }
        if !float_mode {
            out.push(warn(
                "heel metrics only appear when the fleet carries mass (equilibrium mode)",
            ));
        }
    }

    // Row-count ceiling (the CLI rejects sweeps over 100_000 rows).
    let points: Option<usize> = m
        .axes
        .iter()
        .filter(|a| a.kind != AxisKind::Speed)
        .try_fold(1usize, |acc, a| axis_count(&a.values).map(|c| acc * c));
    let speeds: Option<usize> = m
        .axes
        .iter()
        .find(|a| a.kind == AxisKind::Speed)
        .and_then(|a| axis_count(&a.values));
    if let (Some(p), Some(s)) = (points, speeds) {
        let rows = p.max(1) * s;
        if rows > 100_000 {
            out.push(err(format!(
                "sweep would produce {rows} rows; narrow the axes (cap 100000)"
            )));
        } else if rows > 5_000 {
            out.push(warn(format!(
                "sweep produces {rows} rows — this may take a while"
            )));
        }
    }

    out
}

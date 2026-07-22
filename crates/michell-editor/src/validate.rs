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

pub fn validate(m: &Manifest) -> Vec<Issue> {
    let mut out = Vec::new();

    // Hulls.
    if m.hulls.is_empty() {
        out.push(err("manifest has no hulls"));
    }
    for (i, h) in m.hulls.iter().enumerate() {
        if h.id.trim().is_empty() {
            out.push(err(format!("hull #{} needs an id", i + 1)));
        }
        if h.file.trim().is_empty() {
            out.push(err(format!(
                "hull {:?} needs a file",
                if h.id.is_empty() {
                    format!("#{}", i + 1)
                } else {
                    h.id.clone()
                }
            )));
        }
    }
    for i in 0..m.hulls.len() {
        for j in (i + 1)..m.hulls.len() {
            if !m.hulls[i].id.trim().is_empty() && m.hulls[i].id == m.hulls[j].id {
                out.push(err(format!("duplicate hull id {:?}", m.hulls[i].id)));
            }
        }
    }

    // Axis kinds present.
    let count = |k: AxisKind| m.axes.iter().filter(|a| a.kind == k).count();
    let has = |k: AxisKind| count(k) > 0;

    if count(AxisKind::Speed) == 0 {
        out.push(err("the sweep needs a speed axis"));
    } else if count(AxisKind::Speed) > 1 {
        out.push(err("only one speed axis is allowed"));
    }

    let float_mode = has(AxisKind::Weight);
    if has(AxisKind::Lcg) && !float_mode {
        out.push(err("an lcg axis requires a weight axis"));
    }
    if float_mode && has(AxisKind::Waterline) {
        out.push(err(
            "a waterline axis cannot be combined with a weight axis (the waterline is solved)",
        ));
    }
    let vcg_mode = has(AxisKind::Vcg);
    if vcg_mode && !float_mode {
        out.push(err(
            "a vcg axis requires a weight axis (gz is computed at a solved equilibrium)",
        ));
    }
    if has(AxisKind::Heel) && !vcg_mode {
        out.push(err("a heel axis requires a vcg axis so gz is well defined"));
    }

    // Pose axes reference real hulls; a coupled list must be non-empty.
    let ids: Vec<&str> = m.hulls.iter().map(|h| h.id.as_str()).collect();
    for a in &m.axes {
        if a.kind == AxisKind::Pose {
            if a.targets.is_empty() {
                out.push(err("a hull-pose axis needs at least one target hull"));
            }
            for t in &a.targets {
                if !ids.contains(&t.as_str()) {
                    out.push(err(format!("sweep target {t:?} is not a hull id")));
                }
            }
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

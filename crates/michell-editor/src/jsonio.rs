//! Load a manifest into the editing model and write it back out. `serde_json`
//! is used only as the JSON substrate — the mapping to and from the schema is
//! by hand, because axes are a tagged union and their hull-vs-point kind is
//! resolved against the declared ids.

use crate::model::*;
use serde_json::{json, Map, Value};

// -------------------------------------------------------------------------
// Load
// -------------------------------------------------------------------------

pub fn from_str(text: &str) -> Result<Manifest, String> {
    let doc: Value = serde_json::from_str(text).map_err(|e| format!("JSON: {e}"))?;
    from_value(&doc)
}

fn from_value(doc: &Value) -> Result<Manifest, String> {
    let obj = doc.as_object().ok_or("top level must be an object")?;

    let name = obj
        .get("name")
        .and_then(Value::as_str)
        .unwrap_or("")
        .to_string();
    let fluid = match obj
        .get("fluid")
        .and_then(Value::as_str)
        .unwrap_or("seawater")
    {
        "seawater" => Fluid::Seawater,
        "freshwater" => Fluid::Freshwater,
        other => return Err(format!("fluid {other:?}: expected seawater or freshwater")),
    };

    let mut hulls = Vec::new();
    if let Some(arr) = obj.get("hulls").and_then(Value::as_array) {
        for h in arr {
            hulls.push(hull_from(h)?);
        }
    }

    // Axis kind (hull vs point) is resolved against the declared ids.
    let hull_ids: Vec<String> = hulls.iter().map(|h| h.id.clone()).collect();
    let point_ids: Vec<String> = hulls
        .iter()
        .flat_map(|h| h.points.iter().map(|p| p.id.clone()))
        .collect();

    let mut axes = Vec::new();
    if let Some(arr) = obj.get("sweep").and_then(Value::as_array) {
        for a in arr {
            axes.push(axis_from(a, &hull_ids, &point_ids)?);
        }
    }

    let output = match obj.get("output") {
        Some(o) => Output {
            format: match o.get("format").and_then(Value::as_str).unwrap_or("csv") {
                "csv" => OutputFormat::Csv,
                "json" => OutputFormat::Json,
                "binary" => OutputFormat::Binary,
                other => {
                    return Err(format!(
                        "output format {other:?}: expected csv, json, or binary"
                    ))
                }
            },
            file: o
                .get("file")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
        },
        None => Output::default(),
    };

    let options = options_from(obj.get("options"));

    Ok(Manifest {
        name,
        fluid,
        hulls,
        axes,
        output,
        options,
    })
}

fn hull_from(h: &Value) -> Result<HullSpec, String> {
    let id = h
        .get("id")
        .and_then(Value::as_str)
        .ok_or("every hull needs an \"id\"")?;
    let file = h
        .get("file")
        .and_then(Value::as_str)
        .ok_or_else(|| format!("hull {id:?} needs a \"file\""))?;

    let mut pose = Pose::default();
    if let Some(p) = h.get("pose") {
        pose.enabled = true;
        pose.dx = p.get("dx").and_then(Value::as_f64).unwrap_or(0.0);
        pose.dy = p.get("dy").and_then(Value::as_f64).unwrap_or(0.0);
        pose.dz = p.get("dz").and_then(Value::as_f64).unwrap_or(0.0);
        pose.trim_deg = p.get("trim").and_then(Value::as_f64).unwrap_or(0.0);
        pose.scale = p.get("scale").and_then(Value::as_f64).unwrap_or(1.0);
    }

    let mut load = Load::default();
    if let Some(l) = h.get("load") {
        load.enabled = true;
        load.mass = l.get("mass").and_then(Value::as_f64).unwrap_or(0.0);
        load.vcg = l.get("vcg").and_then(Value::as_f64).unwrap_or(0.0);
        if let Some(lcg) = l.get("lcg").and_then(Value::as_f64) {
            load.lcg_set = true;
            load.lcg = lcg;
        }
    }

    let mut points = Vec::new();
    if let Some(arr) = h.get("points").and_then(Value::as_array) {
        for p in arr {
            let pid = p
                .get("id")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("hull {id:?}: every point load needs an \"id\""))?;
            points.push(PointLoad {
                id: pid.to_string(),
                mass: p.get("mass").and_then(Value::as_f64).unwrap_or(0.0),
                dx: p.get("dx").and_then(Value::as_f64).unwrap_or(0.0),
                dy: p.get("dy").and_then(Value::as_f64).unwrap_or(0.0),
                dz: p.get("dz").and_then(Value::as_f64).unwrap_or(0.0),
            });
        }
    }

    Ok(HullSpec {
        id: id.to_string(),
        file: file.to_string(),
        pose,
        load,
        points,
    })
}

fn axis_from(a: &Value, hull_ids: &[String], point_ids: &[String]) -> Result<AxisSpec, String> {
    let targets: Vec<String> = match a.get("target") {
        Some(Value::String(s)) => vec![s.clone()],
        Some(Value::Array(arr)) => arr
            .iter()
            .map(|t| {
                t.as_str()
                    .map(str::to_string)
                    .ok_or_else(|| "target list entries must be strings".to_string())
            })
            .collect::<Result<_, _>>()?,
        _ => return Err("every sweep entry needs a \"target\"".into()),
    };
    let param = a.get("param").and_then(Value::as_str);

    // Reserved single-target global axes.
    if targets.len() == 1 && param.is_none() {
        match targets[0].as_str() {
            "speed" => {
                let mut axis = AxisSpec::new(AxisKind::Speed);
                axis.unit = match a.get("unit").and_then(Value::as_str).unwrap_or("ms") {
                    "ms" | "m/s" => SpeedUnit::Ms,
                    "knots" | "kn" => SpeedUnit::Knots,
                    "froude" | "fn" => SpeedUnit::Froude,
                    other => return Err(format!("unknown speed unit {other:?}")),
                };
                axis.values = values_from(a)?;
                return Ok(axis);
            }
            "waterline" => {
                let mut axis = AxisSpec::new(AxisKind::Waterline);
                axis.values = values_from(a)?;
                return Ok(axis);
            }
            _ => {}
        }
    }

    // Otherwise a hull- or point-targeted axis: classify by the first target's
    // id (unknown ids default to a hull axis so the user can fix them).
    let is_point = point_ids.iter().any(|p| p == &targets[0]);
    let mut axis = AxisSpec::new(if is_point {
        AxisKind::Point
    } else {
        AxisKind::Hull
    });
    axis.targets = targets;
    if let Some(p) = param {
        if is_point {
            axis.point_param = match p {
                "mass" => PointParam::Mass,
                "dx" => PointParam::Dx,
                "dy" => PointParam::Dy,
                "dz" => PointParam::Dz,
                other => return Err(format!("unknown point-load param {other:?}")),
            };
        } else {
            axis.hull_param = match p {
                "dx" => HullParam::Dx,
                "dy" => HullParam::Dy,
                "dz" => HullParam::Dz,
                "trim" => HullParam::Trim,
                "spread" => HullParam::Spread,
                "scale" => HullParam::Scale,
                "mass" => HullParam::Mass,
                "lcg" => HullParam::Lcg,
                "vcg" => HullParam::Vcg,
                other => return Err(format!("unknown hull param {other:?}")),
            };
        }
    }
    let _ = hull_ids;
    axis.values = values_from(a)?;
    Ok(axis)
}

fn values_from(a: &Value) -> Result<ValueSpec, String> {
    let mut v = ValueSpec::default();
    if let Some(x) = a.get("value").and_then(Value::as_f64) {
        v.mode = ValueMode::Scalar;
        v.scalar = fmt_num(x);
    } else if let Some(arr) = a.get("values").and_then(Value::as_array) {
        v.mode = ValueMode::List;
        let nums: Option<Vec<f64>> = arr.iter().map(Value::as_f64).collect();
        let nums = nums.ok_or("\"values\" must be an array of numbers")?;
        v.list = nums
            .iter()
            .map(|n| fmt_num(*n))
            .collect::<Vec<_>>()
            .join(", ");
    } else if let Some(arr) = a.get("range").and_then(Value::as_array) {
        v.mode = ValueMode::Range;
        if arr.len() != 2 {
            return Err("\"range\" must be [start, stop]".into());
        }
        v.range_start = fmt_num(arr[0].as_f64().ok_or("range start must be a number")?);
        v.range_stop = fmt_num(arr[1].as_f64().ok_or("range stop must be a number")?);
        if let Some(s) = a.get("step").and_then(Value::as_f64) {
            v.step = fmt_num(s);
        }
    } else {
        return Err("sweep entry needs one of: value, values, range".into());
    }
    Ok(v)
}

fn options_from(o: Option<&Value>) -> Options {
    let mut opts = Options::default();
    let Some(o) = o else { return opts };
    let set_str = |f: &mut OptField, v: Option<&str>| {
        if let Some(s) = v {
            f.enabled = true;
            f.text = s.to_string();
        }
    };
    let set_num = |f: &mut OptField, v: Option<f64>| {
        if let Some(n) = v {
            f.enabled = true;
            f.text = fmt_num(n);
        }
    };
    set_str(&mut opts.samples, o.get("samples").and_then(Value::as_str));
    set_str(
        &mut opts.fit_degree,
        o.get("fit_degree").and_then(Value::as_str),
    );
    set_str(
        &mut opts.fit_control,
        o.get("fit_control").and_then(Value::as_str),
    );
    set_num(&mut opts.rel_tol, o.get("rel_tol").and_then(Value::as_f64));
    set_num(
        &mut opts.form_factor,
        o.get("form_factor").and_then(Value::as_f64),
    );
    set_num(&mut opts.gravity, o.get("gravity").and_then(Value::as_f64));
    set_num(&mut opts.rho, o.get("rho").and_then(Value::as_f64));
    set_num(&mut opts.nu, o.get("nu").and_then(Value::as_f64));

    if let Some(h) = o.get("heel") {
        opts.heel.enabled = true;
        if let Some(arr) = h.get("resistance_angles").and_then(Value::as_array) {
            let nums: Vec<String> = arr.iter().filter_map(Value::as_f64).map(fmt_num).collect();
            opts.heel.resistance_angles = nums.join(", ");
        }
        if let Some(s) = h.get("gz_step").and_then(Value::as_f64) {
            opts.heel.gz_step = fmt_num(s);
        }
        if let Some(s) = h.get("gz_max").and_then(Value::as_f64) {
            opts.heel.gz_max = fmt_num(s);
        }
    }
    opts
}

// -------------------------------------------------------------------------
// Save
// -------------------------------------------------------------------------

/// Serialize the model to pretty JSON. Returns an error if any numeric field
/// fails to parse, so a save never silently drops a typo.
pub fn to_string(m: &Manifest) -> Result<String, String> {
    let v = to_value(m)?;
    serde_json::to_string_pretty(&v).map_err(|e| e.to_string())
}

fn to_value(m: &Manifest) -> Result<Value, String> {
    let mut root = Map::new();
    if !m.name.is_empty() {
        root.insert("name".into(), json!(m.name));
    }
    root.insert("fluid".into(), json!(m.fluid.as_str()));

    let hulls: Vec<Value> = m.hulls.iter().map(hull_to).collect();
    root.insert("hulls".into(), Value::Array(hulls));

    let mut axes = Vec::new();
    for a in &m.axes {
        axes.push(axis_to(a)?);
    }
    root.insert("sweep".into(), Value::Array(axes));

    let mut out = Map::new();
    out.insert("format".into(), json!(m.output.format.as_str()));
    if !m.output.file.is_empty() {
        out.insert("file".into(), json!(m.output.file));
    }
    root.insert("output".into(), Value::Object(out));

    if !m.options.is_empty() {
        root.insert("options".into(), options_to(&m.options)?);
    }

    Ok(Value::Object(root))
}

fn hull_to(h: &HullSpec) -> Value {
    let mut o = Map::new();
    o.insert("id".into(), json!(h.id));
    o.insert("file".into(), json!(h.file));

    if h.pose.enabled {
        let mut p = Map::new();
        if h.pose.dx != 0.0 {
            p.insert("dx".into(), json!(h.pose.dx));
        }
        if h.pose.dy != 0.0 {
            p.insert("dy".into(), json!(h.pose.dy));
        }
        if h.pose.dz != 0.0 {
            p.insert("dz".into(), json!(h.pose.dz));
        }
        if h.pose.trim_deg != 0.0 {
            p.insert("trim".into(), json!(h.pose.trim_deg));
        }
        if h.pose.scale != 1.0 {
            p.insert("scale".into(), json!(h.pose.scale));
        }
        o.insert("pose".into(), Value::Object(p));
    }

    if h.load.enabled {
        let mut l = Map::new();
        l.insert("mass".into(), num_value(h.load.mass));
        if h.load.lcg_set {
            l.insert("lcg".into(), num_value(h.load.lcg));
        }
        if h.load.vcg != 0.0 {
            l.insert("vcg".into(), num_value(h.load.vcg));
        }
        o.insert("load".into(), Value::Object(l));
    }

    if !h.points.is_empty() {
        let pts: Vec<Value> = h
            .points
            .iter()
            .map(|p| {
                let mut po = Map::new();
                po.insert("id".into(), json!(p.id));
                po.insert("mass".into(), num_value(p.mass));
                if p.dx != 0.0 {
                    po.insert("dx".into(), json!(p.dx));
                }
                if p.dy != 0.0 {
                    po.insert("dy".into(), json!(p.dy));
                }
                if p.dz != 0.0 {
                    po.insert("dz".into(), json!(p.dz));
                }
                Value::Object(po)
            })
            .collect();
        o.insert("points".into(), Value::Array(pts));
    }

    Value::Object(o)
}

fn axis_to(a: &AxisSpec) -> Result<Value, String> {
    let mut o = Map::new();
    let label: String = match a.kind {
        AxisKind::Speed => {
            o.insert("target".into(), json!("speed"));
            o.insert("unit".into(), json!(a.unit.as_str()));
            "speed".into()
        }
        AxisKind::Waterline => {
            o.insert("target".into(), json!("waterline"));
            "waterline".into()
        }
        AxisKind::Hull => {
            if a.targets.is_empty() {
                return Err("a hull-param axis needs at least one target hull".into());
            }
            insert_targets(&mut o, &a.targets);
            o.insert("param".into(), json!(a.hull_param.as_str()));
            a.hull_param.as_str().into()
        }
        AxisKind::Point => {
            if a.targets.is_empty() {
                return Err("a point-load axis needs at least one target point".into());
            }
            insert_targets(&mut o, &a.targets);
            o.insert("param".into(), json!(a.point_param.as_str()));
            a.point_param.as_str().into()
        }
    };
    write_values(&mut o, &a.values, &label)?;
    Ok(Value::Object(o))
}

fn insert_targets(o: &mut Map<String, Value>, targets: &[String]) {
    if targets.len() == 1 {
        o.insert("target".into(), json!(targets[0]));
    } else {
        o.insert("target".into(), json!(targets));
    }
}

fn write_values(o: &mut Map<String, Value>, v: &ValueSpec, label: &str) -> Result<(), String> {
    let ctx = |what: &str| format!("{label} axis: {what}");
    match v.mode {
        ValueMode::Scalar => {
            let n = parse_num(&v.scalar).map_err(|e| ctx(&format!("value {e}")))?;
            o.insert("value".into(), num_value(n));
        }
        ValueMode::List => {
            let nums = parse_list(&v.list).map_err(|e| ctx(&format!("values: {e}")))?;
            if nums.is_empty() {
                return Err(ctx("values list is empty"));
            }
            o.insert(
                "values".into(),
                Value::Array(nums.into_iter().map(num_value).collect()),
            );
        }
        ValueMode::Range => {
            let a = parse_num(&v.range_start).map_err(|e| ctx(&format!("range start {e}")))?;
            let b = parse_num(&v.range_stop).map_err(|e| ctx(&format!("range stop {e}")))?;
            o.insert(
                "range".into(),
                Value::Array(vec![num_value(a), num_value(b)]),
            );
            if !v.step.trim().is_empty() {
                let s = parse_num(&v.step).map_err(|e| ctx(&format!("step {e}")))?;
                o.insert("step".into(), num_value(s));
            }
        }
    }
    Ok(())
}

fn options_to(opts: &Options) -> Result<Value, String> {
    let mut o = Map::new();
    let put_str = |o: &mut Map<String, Value>, key: &str, f: &OptField| {
        if f.enabled {
            o.insert(key.into(), json!(f.text.trim()));
        }
    };
    let put_num = |o: &mut Map<String, Value>, key: &str, f: &OptField| -> Result<(), String> {
        if f.enabled {
            let n = parse_num(&f.text).map_err(|e| format!("options.{key} {e}"))?;
            o.insert(key.into(), num_value(n));
        }
        Ok(())
    };
    put_str(&mut o, "samples", &opts.samples);
    put_str(&mut o, "fit_degree", &opts.fit_degree);
    put_str(&mut o, "fit_control", &opts.fit_control);
    put_num(&mut o, "rel_tol", &opts.rel_tol)?;
    put_num(&mut o, "form_factor", &opts.form_factor)?;
    put_num(&mut o, "gravity", &opts.gravity)?;
    put_num(&mut o, "rho", &opts.rho)?;
    put_num(&mut o, "nu", &opts.nu)?;

    if opts.heel.enabled {
        let mut h = Map::new();
        if !opts.heel.resistance_angles.trim().is_empty() {
            let angles = parse_list(&opts.heel.resistance_angles)
                .map_err(|e| format!("options.heel.resistance_angles: {e}"))?;
            h.insert(
                "resistance_angles".into(),
                Value::Array(angles.into_iter().map(num_value).collect()),
            );
        }
        if !opts.heel.gz_step.trim().is_empty() {
            let s =
                parse_num(&opts.heel.gz_step).map_err(|e| format!("options.heel.gz_step {e}"))?;
            h.insert("gz_step".into(), num_value(s));
        }
        if !opts.heel.gz_max.trim().is_empty() {
            let s = parse_num(&opts.heel.gz_max).map_err(|e| format!("options.heel.gz_max {e}"))?;
            h.insert("gz_max".into(), num_value(s));
        }
        o.insert("heel".into(), Value::Object(h));
    }
    Ok(Value::Object(o))
}

// -------------------------------------------------------------------------
// Number helpers
// -------------------------------------------------------------------------

fn parse_num(s: &str) -> Result<f64, String> {
    let t = s.trim();
    if t.is_empty() {
        return Err("is empty".into());
    }
    t.parse::<f64>()
        .map_err(|_| format!("{t:?} is not a number"))
}

fn parse_list(s: &str) -> Result<Vec<f64>, String> {
    s.split([',', ' ', '\t', '\n'])
        .filter(|t| !t.trim().is_empty())
        .map(|t| {
            t.trim()
                .parse::<f64>()
                .map_err(|_| format!("{:?} is not a number", t.trim()))
        })
        .collect()
}

/// Prefer an integer JSON token when the value is whole, so `700.0` round-trips
/// as `700` — matching how the hand-written manifests read.
fn num_value(n: f64) -> Value {
    if n.is_finite() && n.fract() == 0.0 && n.abs() < 1e15 {
        json!(n as i64)
    } else {
        json!(n)
    }
}

/// Render a number for a text field: whole numbers without a trailing `.0`.
fn fmt_num(n: f64) -> String {
    if n.is_finite() && n.fract() == 0.0 && n.abs() < 1e15 {
        format!("{}", n as i64)
    } else {
        format!("{n}")
    }
}

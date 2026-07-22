//! Load a manifest into the editing model and write it back out. `serde_json`
//! is used only as the JSON substrate — the mapping to and from the schema is
//! by hand, because the axis shape is a tagged union that no derive expresses
//! cleanly.

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

    let mut axes = Vec::new();
    if let Some(arr) = obj.get("sweep").and_then(Value::as_array) {
        for a in arr {
            axes.push(axis_from(a)?);
        }
    }

    let output = match obj.get("output") {
        Some(o) => Output {
            format: match o.get("format").and_then(Value::as_str).unwrap_or("csv") {
                "csv" => OutputFormat::Csv,
                "json" => OutputFormat::Json,
                other => return Err(format!("output format {other:?}: expected csv or json")),
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
    }
    Ok(HullSpec {
        id: id.to_string(),
        file: file.to_string(),
        pose,
    })
}

fn axis_from(a: &Value) -> Result<AxisSpec, String> {
    // Target is a string or a list of strings.
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

    let mut axis = if targets.len() == 1 && param.is_none() {
        match targets[0].as_str() {
            "speed" => AxisSpec::new(AxisKind::Speed),
            "weight" => AxisSpec::new(AxisKind::Weight),
            "lcg" => AxisSpec::new(AxisKind::Lcg),
            "vcg" => AxisSpec::new(AxisKind::Vcg),
            "heel" => AxisSpec::new(AxisKind::Heel),
            "waterline" => AxisSpec::new(AxisKind::Waterline),
            // A single hull id with no param is still a pose axis (but then it
            // needs a param — leave that to the validator).
            _ => pose_axis(targets, param)?,
        }
    } else {
        pose_axis(targets, param)?
    };

    if axis.kind == AxisKind::Speed {
        axis.unit = match a.get("unit").and_then(Value::as_str).unwrap_or("ms") {
            "ms" | "m/s" => SpeedUnit::Ms,
            "knots" | "kn" => SpeedUnit::Knots,
            "froude" | "fn" => SpeedUnit::Froude,
            other => return Err(format!("unknown speed unit {other:?}")),
        };
    }

    axis.values = values_from(a)?;
    Ok(axis)
}

fn pose_axis(targets: Vec<String>, param: Option<&str>) -> Result<AxisSpec, String> {
    let mut axis = AxisSpec::new(AxisKind::Pose);
    axis.targets = targets;
    if let Some(p) = param {
        axis.param = match p {
            "dx" => PoseParam::Dx,
            "dy" => PoseParam::Dy,
            "dz" => PoseParam::Dz,
            "spread" => PoseParam::Spread,
            "trim" => PoseParam::Trim,
            other => return Err(format!("unknown pose param {other:?}")),
        };
    }
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
        // Emit only the non-zero components to keep the file tidy.
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
        o.insert("pose".into(), Value::Object(p));
    }
    Value::Object(o)
}

fn axis_to(a: &AxisSpec) -> Result<Value, String> {
    let mut o = Map::new();
    match a.kind {
        AxisKind::Pose => {
            if a.targets.is_empty() {
                return Err("a hull-pose axis needs at least one target hull".into());
            }
            if a.targets.len() == 1 {
                o.insert("target".into(), json!(a.targets[0]));
            } else {
                o.insert("target".into(), json!(a.targets));
            }
            o.insert("param".into(), json!(a.param.as_str()));
        }
        AxisKind::Speed => {
            o.insert("target".into(), json!("speed"));
            o.insert("unit".into(), json!(a.unit.as_str()));
        }
        kind => {
            o.insert("target".into(), json!(kind.label()));
        }
    }
    write_values(&mut o, &a.values, a.kind)?;
    Ok(Value::Object(o))
}

fn write_values(o: &mut Map<String, Value>, v: &ValueSpec, kind: AxisKind) -> Result<(), String> {
    let ctx = |what: &str| format!("{} axis: {}", kind.label(), what);
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

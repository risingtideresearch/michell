//! Making records: hulls from files (their geometry inline), cases on hulls, studies on cases. Each
//! is saved in the store as it is made, and its id is that of its inputs,
//! so making the same thing twice finds the first.

use crate::store::{id_of, sha256, short, Store};
use boatmath::params::{CaseParams, StudyParams};
use boatmath::{hull_summary, loft, statics, LoftRequest, SOLVER_VERSION};
use serde_json::{json, Value};
use std::path::Path;

/// The type a record must be, or an error naming what it is instead.
pub fn expect<'a>(r: &'a Value, kind: &str) -> Result<&'a str, String> {
    let id = r["id"].as_str().unwrap_or("?");
    match r["type"].as_str() {
        Some(t) if t == kind => Ok(id),
        Some(t) => Err(format!("expected a {kind}, got a {t} ({})", short(id))),
        None => Err(format!("expected a {kind} record, got {}", truncate(r))),
    }
}

fn truncate(v: &Value) -> String {
    let s = v.to_string();
    if s.len() > 60 {
        format!("{}…", &s[..60])
    } else {
        s
    }
}

/// How a hull is cut into sections: settings of the cut, not of the
/// geometry, so they stay with the hull record.
#[derive(Clone, Copy, Debug, Default, serde::Serialize, serde::Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CutSettings {
    /// Stations along the hull.
    pub stations: Option<usize>,
    /// Rays across each section.
    pub rays: Option<usize>,
    /// Centreplane override [m, y].
    pub centerplane: Option<f64>,
}

/// What the solver needs of a hull: its geometry, as the bytes `boatmath`
/// reads, and how it is cut.
pub struct HullSource {
    pub file_name: String,
    pub bytes: Vec<u8>,
    pub import: LoftRequest,
}

pub fn hull_source(hull: &Value) -> Result<HullSource, String> {
    let cut: CutSettings =
        serde_json::from_value(hull["cut"].clone()).map_err(|e| format!("hull cut: {e}"))?;
    if !hull["geometry"].is_object() {
        return Err("a hull without geometry".into());
    }
    Ok(HullSource {
        file_name: "geometry.json".into(),
        bytes: serde_json::to_vec(&hull["geometry"]).map_err(|e| e.to_string())?,
        import: LoftRequest {
            stations: cut.stations,
            rays: cut.rays,
            centerplane: cut.centerplane,
            ..LoftRequest::default()
        },
    })
}

fn hull_record(
    store: &Store,
    name: &str,
    source: Value,
    geometry: Value,
    cut: CutSettings,
    parent: Option<&str>,
) -> Result<Value, String> {
    let mut r = json!({
        "type": "hull",
        "id": id_of(&json!({ "geometry": geometry, "cut": cut })),
        "name": name,
        "source": source,
        "cut": cut,
        "parent": parent,
        "geometry": geometry,
        "solver_version": SOLVER_VERSION,
    });
    let src = hull_source(&r)?;
    let sections = loft(&src.file_name, src.bytes, &src.import)?;
    r["summary"] = hull_summary(&sections);
    store.put(&r)?;
    Ok(r)
}

/// A hull from a file: its geometry, with the file's import settings
/// (waterline, units) applied once. `source` records where it came from.
pub fn hull_from_file(
    store: &Store,
    path: &Path,
    name: Option<&str>,
    import: LoftRequest,
) -> Result<Value, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let file_name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "hull".into());
    let stem = path
        .file_stem()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| file_name.clone());
    let abs = std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf());
    let source = json!({
        "path": abs.to_string_lossy(),
        "sha256": sha256(&bytes),
        "waterline": import.waterline,
        "units": import.units,
    });
    let geometry = boatmath::native::from_file(&file_name, bytes, &import)
        .map_err(|e| format!("{file_name}: {e}"))?;
    let cut = CutSettings {
        stations: import.stations,
        rays: import.rays,
        centerplane: import.centerplane,
    };
    hull_record(store, name.unwrap_or(&stem), source, geometry, cut, None)
}

/// A hull made by scaling another: wholly by `by`, and/or its beam and
/// draft by `yz`, about its design waterline. The new geometry is the old
/// with its control points (or vertices) scaled.
pub fn scaled_hull(
    store: &Store,
    hull: &Value,
    by: Option<f64>,
    yz: Option<f64>,
    name: Option<&str>,
) -> Result<Value, String> {
    let id = expect(hull, "hull")?;
    for (k, v) in [("--by", by), ("--beam", yz)] {
        if let Some(x) = v {
            if !(x > 0.0 && x.is_finite()) {
                return Err(format!("{k} {x}: expected a positive factor"));
            }
        }
    }
    let pose = LoftRequest {
        scale: by,
        scale_yz: yz,
        ..LoftRequest::default()
    }
    .pose();
    let geometry = boatmath::native::posed(&hull["geometry"], &pose)?;
    let cut: CutSettings =
        serde_json::from_value(hull["cut"].clone()).map_err(|e| format!("hull cut: {e}"))?;
    let base = hull["name"].as_str().unwrap_or("hull");
    let default_name = match (by, yz) {
        (Some(k), None) => format!("{base} ×{k}"),
        (None, Some(k)) => format!("{base} ×{k} beam"),
        (Some(k), Some(j)) => format!("{base} ×{k}, ×{j} beam"),
        (None, None) => base.to_string(),
    };
    let source = json!({ "scaled": { "by": by, "beam": yz } });
    hull_record(
        store,
        name.unwrap_or(&default_name),
        source,
        geometry,
        cut,
        Some(id),
    )
}

/// A case on a hull, with its statics. One already in the store, from this
/// solver, is reused unless `force`.
pub fn case(
    store: &Store,
    hull: &Value,
    params: CaseParams,
    name: Option<&str>,
    force: bool,
) -> Result<Value, String> {
    let hull_id = expect(hull, "hull")?;
    let params = params.canonical()?;
    let id = id_of(&json!({ "hull": hull_id, "params": params }));
    if let Some(mut old) = store.get("case", &id)? {
        if !force && old["solver_version"] == SOLVER_VERSION {
            if let Some(n) = name {
                if old["name"] != n {
                    old["name"] = json!(n);
                    store.put(&old)?;
                }
            }
            return Ok(old);
        }
    }
    let src = hull_source(hull)?;
    let mut r = json!({
        "type": "case",
        "id": id,
        "name": name.unwrap_or(""),
        "hull": hull_id,
        "params": params,
        "solver_version": SOLVER_VERSION,
    });
    match statics(&src.file_name, src.bytes, &src.import, &params) {
        Ok(mut s) => {
            let o = s.as_object_mut().expect("statics are an object");
            // The meshes are for drawing; `show` makes them again.
            o.remove("meshes");
            o.remove("body_meshes");
            r["seconds"] = o.remove("seconds").unwrap_or(Value::Null);
            let rest = (
                s["at_rest"]["sinkage"].as_f64(),
                s["at_rest"]["trim_rad"].as_f64(),
            );
            r["statics"] = s;
            if let (Some(z), Some(t)) = rest {
                match sections(store, &r, (z, t)) {
                    Ok(id) => r["sections"] = json!(id),
                    Err(e) => eprintln!("boatmath: case {}: sections: {e}", short(&id)),
                }
            }
        }
        Err(e) => r["error"] = json!(e),
    }
    store.put(&r)?;
    Ok(r)
}

/// The sections record of a case cut at `(sinkage [m], trim [rad])`,
/// made if the store hasn't it; returns its id. The id is that of the case
/// and the attitude, so a case at rest and the studies held there share one.
pub fn sections(store: &Store, case: &Value, attitude: (f64, f64)) -> Result<String, String> {
    let case_id = expect(case, "case")?;
    let id = id_of(&json!({ "case": case_id, "attitude": [attitude.0, attitude.1] }));
    if store
        .get("sections", &id)?
        .is_some_and(|s| s["solver_version"] == SOLVER_VERSION)
    {
        return Ok(id);
    }
    let hull = store.need("hull", case["hull"].as_str().unwrap_or(""))?;
    let src = hull_source(&hull)?;
    let params: CaseParams =
        serde_json::from_value(case["params"].clone()).map_err(|e| format!("case params: {e}"))?;
    let mut r =
        boatmath::sections::sections(&src.file_name, src.bytes, &src.import, &params, attitude)?;
    let o = r.as_object_mut().expect("an object");
    o.insert("type".into(), json!("sections"));
    o.insert("id".into(), json!(id));
    o.insert("case".into(), json!(case_id));
    o.insert("solver_version".into(), json!(SOLVER_VERSION));
    store.put(&r)?;
    Ok(id)
}

/// A study on a case: only the request, saved; `run` computes it.
pub fn study(store: &Store, case: &Value, params: StudyParams) -> Result<Value, String> {
    let case_id = expect(case, "case")?;
    if let Some(e) = case["error"].as_str() {
        return Err(format!("case {} failed: {e}", short(case_id)));
    }
    let params = params.canonical()?;
    let r = study_record(case_id, params);
    store.put(&r)?;
    Ok(r)
}

pub fn study_record(case_id: &str, params: StudyParams) -> Value {
    json!({
        "type": "study",
        "id": id_of(&json!({ "case": case_id, "params": params })),
        "case": case_id,
        "kind": params.kind(),
        "params": params,
    })
}

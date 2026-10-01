//! Making records: hulls from files, cases on hulls, studies on cases. Each
//! is saved in the store as it is made, and its id is that of its inputs,
//! so making the same thing twice finds the first.

use crate::store::{id_of, short, Store};
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

/// A hull's file and how it is cut.
pub struct HullSource {
    pub file_name: String,
    pub bytes: Vec<u8>,
    pub import: LoftRequest,
}

pub fn hull_source(store: &Store, hull: &Value) -> Result<HullSource, String> {
    let import: LoftRequest =
        serde_json::from_value(hull["import"].clone()).map_err(|e| format!("hull import: {e}"))?;
    let blob = hull["source"]["blob"]
        .as_str()
        .ok_or("hull without a source blob")?;
    Ok(HullSource {
        file_name: hull["source"]["file_name"]
            .as_str()
            .unwrap_or("hull")
            .to_string(),
        bytes: store.blob(blob)?,
        import,
    })
}

fn hull_record(
    store: &Store,
    name: &str,
    source: Value,
    file_name: &str,
    bytes: Vec<u8>,
    import: LoftRequest,
    parent: Option<&str>,
) -> Result<Value, String> {
    let blob = store.put_blob(&bytes)?;
    let sections = loft(file_name, bytes, &import).map_err(|e| format!("{file_name}: {e}"))?;
    let mut source = source;
    source["blob"] = json!(blob);
    source["file_name"] = json!(file_name);
    let id = id_of(&json!({ "file": blob, "import": import }));
    let r = json!({
        "type": "hull",
        "id": id,
        "name": name,
        "source": source,
        "import": import,
        "parent": parent,
        "summary": hull_summary(&sections),
        "solver_version": SOLVER_VERSION,
    });
    store.put(&r)?;
    Ok(r)
}

/// A hull from a file: its bytes go in the store, so the record stays good
/// when the file moves or changes.
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
    hull_record(
        store,
        name.unwrap_or(&stem),
        json!({ "path": abs.to_string_lossy() }),
        &file_name,
        bytes,
        import,
        None,
    )
}

/// A hull made by scaling another: wholly by `by`, and/or its beam and
/// draft by `yz`, about its design waterline.
pub fn scaled_hull(
    store: &Store,
    hull: &Value,
    by: Option<f64>,
    yz: Option<f64>,
    name: Option<&str>,
) -> Result<Value, String> {
    let id = expect(hull, "hull")?;
    let src = hull_source(store, hull)?;
    let mut import = src.import;
    for (k, v) in [("--by", by), ("--beam", yz)] {
        if let Some(x) = v {
            if !(x > 0.0 && x.is_finite()) {
                return Err(format!("{k} {x}: expected a positive factor"));
            }
        }
    }
    if let Some(k) = by {
        import.scale = Some(import.scale.unwrap_or(1.0) * k);
    }
    if let Some(k) = yz {
        import.scale_yz = Some(import.scale_yz.unwrap_or(1.0) * k);
    }
    let base = hull["name"].as_str().unwrap_or("hull");
    let default_name = match (by, yz) {
        (Some(k), None) => format!("{base} ×{k}"),
        (None, Some(k)) => format!("{base} ×{k} beam"),
        (Some(k), Some(j)) => format!("{base} ×{k}, ×{j} beam"),
        (None, None) => base.to_string(),
    };
    hull_record(
        store,
        name.unwrap_or(&default_name),
        hull["source"].clone(),
        &src.file_name,
        src.bytes,
        import,
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
    let src = hull_source(store, hull)?;
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
            r["statics"] = s;
        }
        Err(e) => r["error"] = json!(e),
    }
    store.put(&r)?;
    Ok(r)
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

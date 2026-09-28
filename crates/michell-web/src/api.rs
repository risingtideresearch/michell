//! The queue's JSON API: hulls, studies, cases, their results, and the
//! queue itself. Bodies and answers are JSON unless noted.
//!
//! ```text
//! GET  /api/version                     {"solver_version": ..}
//! GET  /api/hulls                       every hull, newest first
//! POST /api/hulls?file=F&name=&notes=&by=&parent=&waterline=&centerplane=&stations=&rays=&units=
//!                                       body = the file's bytes; cut, saved →
//!                                       {"id", "created", "hull"}
//! GET  /api/hulls/:id                   one hull
//! POST /api/hulls/:id                   {"name"?, "notes"?}
//! GET  /api/hulls/:id/sections          the cut (as /api/loft answers)
//! GET  /api/hulls/:id/geometry          the display geometry (as /api/geometry)
//! GET  /api/hulls/:id/file              the uploaded file
//! GET  /api/studies                     every study, with case counts
//! POST /api/cases                       {"hull", "cases": [CaseParams..],
//!                                        "study"?: {"name", "notes"?} | "study_id"?,
//!                                        "priority"?, "by"?}
//!                                       → {"study_id", "cases": [{"id", "created"}..]}
//! GET  /api/cases?hull=&study=&status=&limit=
//!                                       cases with their results' scalars
//! GET  /api/cases/:id                   one case
//! GET  /api/cases/:id/results/:idx      that result's whole answer (fields and all)
//! POST /api/cases/:id/cancel            stop it (queued, or running)
//! POST /api/cases/:id/requeue           run it again
//! POST /api/cases/:id/priority?p=N      higher runs sooner
//! POST /api/requeue_stale?hull=         requeue every case with a stale result
//! GET  /api/queue                       {"running", "cases", "counts"}
//! ```

use michell_web::case::CaseParams;
use michell_web::store::{CaseFilter, NewHull, Store, SOLVER_VERSION};
use michell_web::worker::Worker;
use michell_web::{geometry, loft, LoftRequest};
use serde_json::{json, Value};
use std::sync::Arc;
use tiny_http::Method;

pub struct App {
    pub store: Arc<Store>,
    pub worker: Arc<Worker>,
}

pub enum Reply {
    Json(Value),
    /// A stored blob, gzipped, and its content type.
    Gz(Vec<u8>, &'static str),
}

/// An answer that is not a success: a status and a message.
pub type Fail = (u16, String);

fn bad(e: impl std::fmt::Display) -> Fail {
    (400, e.to_string())
}

fn missing(what: &str) -> Fail {
    (404, format!("no such {what}"))
}

fn internal(e: String) -> Fail {
    (500, e)
}

/// Route an `/api/…` request; `None` if no route matches.
pub fn route(
    app: &App,
    method: &Method,
    path: &str,
    q: &[(String, String)],
    body: &mut dyn FnMut() -> Result<Vec<u8>, String>,
) -> Option<Result<Reply, Fail>> {
    let seg: Vec<&str> = path.trim_matches('/').split('/').collect();
    let get = |k: &str| {
        q.iter()
            .find(|(n, v)| n == k && !v.trim().is_empty())
            .map(|(_, v)| v.trim())
    };
    let id = |s: &str| s.parse::<i64>().map_err(|_| bad(format!("bad id {s:?}")));
    let opt_id = |k: &str| get(k).map(id).transpose();
    let json_body = |body: &mut dyn FnMut() -> Result<Vec<u8>, String>| -> Result<Value, Fail> {
        serde_json::from_slice(&body().map_err(bad)?).map_err(|e| bad(format!("body: {e}")))
    };
    let s = &app.store;
    use Method::{Get, Post};
    let r = match (method, &seg[1..]) {
        (Get, ["version"]) => Ok(Reply::Json(json!({ "solver_version": SOLVER_VERSION }))),
        (Get, ["hulls"]) => s.hulls().map(|v| Reply::Json(json!(v))).map_err(internal),
        (Post, ["hulls"]) => upload(app, q, body),
        (Get, ["hulls", h]) => id(h).and_then(|h| {
            s.hull(h)
                .map_err(internal)?
                .map(Reply::Json)
                .ok_or(missing("hull"))
        }),
        (Post, ["hulls", h]) => id(h).and_then(|h| {
            let v = json_body(body)?;
            let text = |k| v[k].as_str().map(str::trim);
            if s.update_hull(h, text("name"), text("notes"))
                .map_err(internal)?
            {
                Ok(Reply::Json(
                    s.hull(h).map_err(internal)?.unwrap_or_default(),
                ))
            } else {
                Err(missing("hull"))
            }
        }),
        (Get, ["hulls", h, which @ ("file" | "sections" | "geometry")]) => id(h).and_then(|h| {
            let blob = s
                .hull_blob(h, which)
                .map_err(internal)?
                .ok_or(missing("hull"))?;
            let kind = if *which == "file" {
                "application/octet-stream"
            } else {
                "application/json"
            };
            Ok(Reply::Gz(s.blob_gz(&blob).map_err(internal)?, kind))
        }),
        (Get, ["studies"]) => s.studies().map(|v| Reply::Json(json!(v))).map_err(internal),
        (Post, ["cases"]) => json_body(body).and_then(|v| add_cases(app, &v)),
        (Get, ["cases"]) => (|| {
            let f = CaseFilter {
                hull: opt_id("hull")?,
                study: opt_id("study")?,
                status: get("status").map(String::from),
                limit: opt_id("limit")?,
            };
            s.cases(&f).map(|v| Reply::Json(json!(v))).map_err(internal)
        })(),
        (Get, ["cases", c]) => id(c).and_then(|c| {
            s.case(c)
                .map_err(internal)?
                .map(Reply::Json)
                .ok_or(missing("case"))
        }),
        (Get, ["cases", c, "results", i]) => id(c).and_then(|c| {
            let blob = s
                .result_blob(c, id(i)?)
                .map_err(internal)?
                .ok_or(missing("result"))?;
            Ok(Reply::Gz(
                s.blob_gz(&blob).map_err(internal)?,
                "application/json",
            ))
        }),
        (Post, ["cases", c, "cancel"]) => id(c).and_then(|c| {
            let before = s
                .cancel_queued(c)
                .map_err(internal)?
                .ok_or(missing("case"))?;
            let stopping = before == "running" && app.worker.cancel(c);
            Ok(Reply::Json(json!({ "was": before, "stopping": stopping })))
        }),
        (Post, ["cases", c, "requeue"]) => id(c).and_then(|c| {
            let n = s.requeue(Some(&[c]), None).map_err(internal)?;
            app.worker.poke();
            Ok(Reply::Json(json!({ "requeued": n })))
        }),
        (Post, ["cases", c, "priority"]) => id(c).and_then(|c| {
            let p = get("p").ok_or(bad("p is required"))?;
            let p = p
                .parse::<i64>()
                .map_err(|_| bad(format!("p: bad number {p:?}")))?;
            if s.set_priority(c, p).map_err(internal)? {
                Ok(Reply::Json(json!({ "id": c, "priority": p })))
            } else {
                Err(missing("case"))
            }
        }),
        (Post, ["requeue_stale"]) => opt_id("hull").and_then(|h| {
            let n = s.requeue(None, h).map_err(internal)?;
            app.worker.poke();
            Ok(Reply::Json(json!({ "requeued": n })))
        }),
        (Get, ["queue"]) => (|| {
            Ok(Reply::Json(json!({
                "running": app.worker.status(),
                "cases": s.queue().map_err(internal)?,
                "counts": s.counts().map_err(internal)?,
                "solver_version": SOLVER_VERSION,
            })))
        })(),
        _ => return None,
    };
    Some(r)
}

/// `POST /api/hulls`: cut the upload (refusing what cannot be cut) and save
/// it with its sections and display geometry.
fn upload(
    app: &App,
    q: &[(String, String)],
    body: &mut dyn FnMut() -> Result<Vec<u8>, String>,
) -> Result<Reply, Fail> {
    let get = |k: &str| q.iter().find(|(n, _)| n == k).map_or("", |(_, v)| v.trim());
    let file = get("file");
    if file.is_empty() {
        return Err(bad("file (the file's name, for its format) is required"));
    }
    let import = LoftRequest::from_query(q).map_err(bad)?;
    let parent_id = match get("parent") {
        "" => None,
        p => Some(
            p.parse::<i64>()
                .map_err(|_| bad(format!("parent: bad id {p:?}")))?,
        ),
    };
    let bytes = body().map_err(bad)?;
    let t0 = std::time::Instant::now();
    let sections = loft(file, bytes.clone(), &import).map_err(bad)?;
    let geom = geometry(file, bytes.clone(), import.units).map_err(bad)?;
    let stem = file.rsplit_once('.').map_or(file, |(s, _)| s);
    let name = match get("name") {
        "" => stem,
        n => n,
    };
    let (id, created) = app
        .store
        .add_hull(&NewHull {
            name,
            notes: get("notes"),
            uploaded_by: get("by"),
            file_name: file,
            bytes: &bytes,
            import,
            parent_id,
            sections: &sections,
            geometry: &geom,
        })
        .map_err(internal)?;
    eprintln!(
        "hull {id} ({file}): {} ({:.2} s)",
        if created { "saved" } else { "already saved" },
        t0.elapsed().as_secs_f64()
    );
    Ok(Reply::Json(json!({
        "id": id,
        "created": created,
        "hull": app.store.hull(id).map_err(internal)?,
    })))
}

/// `POST /api/cases`: every case checked before any is queued.
fn add_cases(app: &App, v: &Value) -> Result<Reply, Fail> {
    let s = &app.store;
    let hull = v["hull"].as_i64().ok_or(bad("hull (an id) is required"))?;
    if s.hull(hull).map_err(internal)?.is_none() {
        return Err(missing("hull"));
    }
    let cases: Vec<CaseParams> = v["cases"]
        .as_array()
        .filter(|a| !a.is_empty())
        .ok_or(bad("cases: expected a list of cases"))?
        .iter()
        .enumerate()
        .map(|(i, c)| {
            serde_json::from_value::<CaseParams>(c.clone())
                .map_err(|e| e.to_string())
                .and_then(CaseParams::canonical)
                .map_err(|e| bad(format!("case {}: {e}", i + 1)))
        })
        .collect::<Result<_, _>>()?;
    let by = v["by"].as_str().unwrap_or("").trim();
    let priority = v["priority"].as_i64().unwrap_or(0);
    let study_id = match (&v["study"], v["study_id"].as_i64()) {
        (Value::Object(st), None) => {
            let name = st
                .get("name")
                .and_then(Value::as_str)
                .map(str::trim)
                .filter(|n| !n.is_empty())
                .ok_or(bad("study: a name is required"))?;
            let notes = st.get("notes").and_then(Value::as_str).unwrap_or("");
            Some(s.add_study(name, notes, by).map_err(internal)?)
        }
        (Value::Null, id) => id,
        _ => return Err(bad("give study or study_id, not both")),
    };
    let mut out = Vec::new();
    for c in &cases {
        let (id, created) = s
            .add_case(hull, c, study_id, priority, by)
            .map_err(internal)?;
        out.push(json!({ "id": id, "created": created }));
    }
    app.worker.poke();
    Ok(Reply::Json(json!({ "study_id": study_id, "cases": out })))
}

//! The JSON API: hulls, their cases, studies on those,
//! results, and the queue. Bodies and answers are JSON unless noted.
//!
//! ```text
//! GET  /api/version                     {"solver_version": ..}
//! GET  /api/whoami                      {"login", "name"} from the tailnet, or null
//!
//! Behind `tailscale serve` (the server bound to loopback), uploads and
//! studies are labelled with the tailnet's name for the asker, and a `by`
//! given in the request is ignored.
//!
//! GET  /api/hulls                       every hull, newest first
//! POST /api/hulls?file=F&name=&notes=&by=&parent=&waterline=&centerplane=&stations=&rays=&units=
//!                                       body = the file's bytes; cut, saved →
//!                                       {"id", "created", "hull"}
//! GET  /api/hulls/:id                   one hull
//! POST /api/hulls/:id                   {"name"?, "notes"?}
//! GET  /api/hulls/:id/sections          the cut (as /api/loft answers)
//! GET  /api/hulls/:id/geometry          the display geometry (as /api/geometry)
//! GET  /api/hulls/:id/file              the uploaded file
//!
//! GET  /api/cases?hull=               cases, with their statics' summary
//! POST /api/cases                     {"hull", "params": CaseParams, "name"?, "by"?,
//!                                        "dry_run"?}
//!                                       → {"id", "created", "case"}; the statics
//!                                       are computed before it is saved. A dry study
//!                                       checks it and answers {"params", "existing"}
//! GET  /api/cases/:id                 one (its statics computed again if stale)
//! POST /api/cases/:id                 {"name"?, "notes"?}
//! GET  /api/cases/:id/statics         the whole statics (GZ curve, meshes at rest)
//!
//! POST /api/studies                        {"case_ids": [..], "studies": [StudyParams..],
//!                                        "priority"?, "by"?, "dry_run"?}
//!                                       every study on every case (a study in
//!                                       waves brings its calm-water study) →
//!                                       {"studies": [..]}; a dry run
//!                                       queues nothing and says which exist
//! GET  /api/studies?hull=&case=&ids=a,b&status=a,b&kind=&limit=&sort=finished
//!                                       studies with their case and result
//! GET  /api/studies/:id                    one study
//! GET  /api/studies/:id/result             its whole answer (fields and all)
//! POST /api/studies/:id/cancel             stop it (queued, or running)
//! POST /api/studies/:id/requeue            run it again
//! POST /api/studies/:id/priority?p=N       higher runs sooner
//! POST /api/requeue_stale?hull=         requeue every study with a stale result
//! GET  /api/queue                       {"running", "studies", "counts"}
//! ```

use boatmath_web::params::{CaseParams, StudyParams};
use boatmath_web::store::{NewHull, Store, StudyFilter, SOLVER_VERSION};
use boatmath_web::worker::Worker;
use boatmath_web::{geometry, loft, statics, LoftRequest};
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

/// Who is asking, as `tailscale serve` says (its `Tailscale-User-*`
/// headers): what uploads and studies are labelled with, in place of the name
/// the page sends.
pub struct Who {
    pub login: String,
    pub name: String,
}

/// The label for work asked for: the tailnet's name for the asker when
/// there is one, else what the request gave.
fn by<'a>(who: Option<&'a Who>, given: &'a str) -> &'a str {
    match who {
        Some(w) if !w.name.is_empty() => &w.name,
        Some(w) => &w.login,
        None => given.trim(),
    }
}

/// Route an `/api/…` request; `None` if no route matches.
pub fn route(
    app: &App,
    method: &Method,
    path: &str,
    q: &[(String, String)],
    who: Option<&Who>,
    body: &mut dyn FnMut() -> Result<Vec<u8>, String>,
) -> Option<Result<Reply, Fail>> {
    let seg: Vec<&str> = path.trim_matches('/').split('/').collect();
    if seg.first() != Some(&"api") {
        return None;
    }
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
    let list = |v: Result<Vec<Value>, String>| v.map(|v| Reply::Json(json!(v))).map_err(internal);
    let one = |v: Result<Option<Value>, String>, what: &str| {
        v.map_err(internal)?.map(Reply::Json).ok_or(missing(what))
    };
    use Method::{Get, Post};
    let r = match (method, &seg[1..]) {
        (Get, ["version"]) => Ok(Reply::Json(json!({ "solver_version": SOLVER_VERSION }))),
        (Get, ["whoami"]) => Ok(Reply::Json(match who {
            Some(w) => json!({ "login": w.login, "name": w.name }),
            None => Value::Null,
        })),

        (Get, ["hulls"]) => list(s.hulls()),
        (Post, ["hulls"]) => upload(app, q, who, body),
        (Get, ["hulls", h]) => id(h).and_then(|h| one(s.hull(h), "hull")),
        (Post, ["hulls", h]) => id(h).and_then(|h| {
            let v = json_body(body)?;
            let text = |k| v[k].as_str().map(str::trim);
            if s.update_hull(h, text("name"), text("notes"))
                .map_err(internal)?
            {
                one(s.hull(h), "hull")
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

        (Get, ["cases"]) => opt_id("hull").and_then(|h| list(s.cases(h))),
        (Post, ["cases"]) => json_body(body).and_then(|v| {
            let hull = v["hull"].as_i64().ok_or(bad("hull (an id) is required"))?;
            let p = case_params(&v["params"])?;
            if v["dry_run"].as_bool() == Some(true) {
                // Checked, and found if it exists; nothing made.
                check_case(app, hull, &p)?;
                return Ok(Reply::Json(json!({
                    "params": p,
                    "existing": s.find_case(hull, &p).map_err(internal)?,
                })));
            }
            let by = by(who, v["by"].as_str().unwrap_or(""));
            let (id, created, _) = make_case(app, hull, &p, v["name"].as_str(), by)?;
            Ok(Reply::Json(json!({
                "id": id,
                "created": created,
                "case": s.case(id).map_err(internal)?,
            })))
        }),
        (Get, ["cases", c]) => id(c).and_then(|c| {
            let cfg = s.case(c).map_err(internal)?.ok_or(missing("case"))?;
            if cfg["stale"] == true {
                // Statics are cheap: computed again when looked at.
                recompute(app, c)?;
                return one(s.case(c), "case");
            }
            Ok(Reply::Json(cfg))
        }),
        (Post, ["cases", c]) => id(c).and_then(|c| {
            let v = json_body(body)?;
            let text = |k| v[k].as_str().map(str::trim);
            if s.update_case(c, text("name"), text("notes"))
                .map_err(internal)?
            {
                one(s.case(c), "case")
            } else {
                Err(missing("case"))
            }
        }),
        (Get, ["cases", c, "statics"]) => id(c).and_then(|c| {
            let (_, _, blob) = s.case_source(c).map_err(internal)?.ok_or(missing("case"))?;
            let blob = blob.ok_or((409, "its statics could not be computed".to_string()))?;
            Ok(Reply::Gz(
                s.blob_gz(&blob).map_err(internal)?,
                "application/json",
            ))
        }),

        (Post, ["studies"]) => json_body(body).and_then(|v| add_studies(app, &v, who)),
        (Get, ["studies"]) => (|| {
            let f = StudyFilter {
                hull: opt_id("hull")?,
                case: opt_id("case")?,
                ids: get("ids").map(String::from),
                status: get("status").map(String::from),
                kind: get("kind").map(String::from),
                limit: opt_id("limit")?,
                by_finish: get("sort") == Some("finished"),
            };
            list(s.studies(&f))
        })(),
        (Get, ["studies", r]) => id(r).and_then(|r| one(s.study(r), "study")),
        (Get, ["studies", r, "result"]) => id(r).and_then(|r| {
            let blob = s
                .result_blob(r)
                .map_err(internal)?
                .ok_or(missing("result"))?;
            Ok(Reply::Gz(
                s.blob_gz(&blob).map_err(internal)?,
                "application/json",
            ))
        }),
        (Post, ["studies", r, "cancel"]) => id(r).and_then(|r| {
            let before = s
                .cancel_queued(r)
                .map_err(internal)?
                .ok_or(missing("study"))?;
            let stopping = before == "running" && app.worker.cancel(r);
            Ok(Reply::Json(json!({ "was": before, "stopping": stopping })))
        }),
        (Post, ["studies", r, "requeue"]) => id(r).and_then(|r| {
            let n = s.requeue(Some(&[r]), None).map_err(internal)?;
            app.worker.poke();
            Ok(Reply::Json(json!({ "requeued": n })))
        }),
        (Post, ["studies", r, "priority"]) => id(r).and_then(|r| {
            let p = get("p").ok_or(bad("p is required"))?;
            let p = p
                .parse::<i64>()
                .map_err(|_| bad(format!("p: bad number {p:?}")))?;
            if s.set_priority(r, p).map_err(internal)? {
                Ok(Reply::Json(json!({ "id": r, "priority": p })))
            } else {
                Err(missing("study"))
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
                "studies": s.queue().map_err(internal)?,
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
    who: Option<&Who>,
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
            uploaded_by: by(who, get("by")),
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

fn case_params(v: &Value) -> Result<CaseParams, Fail> {
    let v = if v.is_null() { json!({}) } else { v.clone() };
    serde_json::from_value::<CaseParams>(v)
        .map_err(|e| e.to_string())
        .and_then(CaseParams::canonical)
        .map_err(bad)
}

/// The checks a case must pass before its statics are tried: a
/// catamaran doubles one hull, its demihulls clear of each other.
fn check_case(app: &App, hull: i64, p: &CaseParams) -> Result<(), Fail> {
    let h = app
        .store
        .hull(hull)
        .map_err(internal)?
        .ok_or(missing("hull"))?;
    let hulls = h["summary"]["hulls"].as_array().map_or(0, Vec::len);
    if let Some(span) = p.span {
        if hulls != 1 {
            return Err(bad(format!(
                "a catamaran doubles a single hull; this file holds {hulls}"
            )));
        }
        let beam = h["summary"]["hulls"][0]["beam"].as_f64().unwrap_or(0.0);
        if span <= beam {
            return Err(bad(format!(
                "span {span} m: the demihulls overlap (beam {beam:.3} m)"
            )));
        }
    }
    Ok(())
}

/// A case's statics, from its hull's file.
fn compute_statics(app: &App, hull: i64, p: &CaseParams) -> Result<Result<Value, String>, Fail> {
    let src = app
        .store
        .hull_source(hull)
        .map_err(internal)?
        .ok_or(missing("hull"))?;
    let bytes = app.store.blob(&src.file_blob).map_err(internal)?;
    let t0 = std::time::Instant::now();
    let r = statics(&src.file_name, bytes, &src.import, p);
    eprintln!(
        "statics on hull {hull}: {} ({:.2} s)",
        r.as_ref().map_or_else(|e| e.as_str(), |_| "ok"),
        t0.elapsed().as_secs_f64()
    );
    Ok(r)
}

/// The case on `hull` with these parameters: found, or made (its
/// statics computed first). `(id, created, statics error)`.
fn make_case(
    app: &App,
    hull: i64,
    p: &CaseParams,
    name: Option<&str>,
    by: &str,
) -> Result<(i64, bool, Option<String>), Fail> {
    if let Some(id) = app.store.find_case(hull, p).map_err(internal)? {
        let error = app
            .store
            .case(id)
            .map_err(internal)?
            .and_then(|c| c["error"].as_str().map(String::from));
        return Ok((id, false, error));
    }
    check_case(app, hull, p)?;
    let st = compute_statics(app, hull, p)?;
    let name = name.map(str::trim).unwrap_or("");
    let (id, created) = app
        .store
        .add_case(hull, p, name, by, st.as_ref().map_err(String::as_str))
        .map_err(internal)?;
    Ok((id, created, st.err()))
}

fn recompute(app: &App, id: i64) -> Result<(), Fail> {
    let (hull, p, _) = app
        .store
        .case_source(id)
        .map_err(internal)?
        .ok_or(missing("case"))?;
    let st = compute_statics(app, hull, &p)?;
    app.store
        .set_statics(id, st.as_ref().map_err(String::as_str))
        .map_err(internal)
}

/// `POST /api/studies`: every study on every case named (each made
/// already, its statics computed), all checked before anything is queued; a
/// study in waves brings the calm-water study it is taken about (with the grid
/// of a calm study asked for alongside, if any).
fn add_studies(app: &App, v: &Value, who: Option<&Who>) -> Result<Reply, Fail> {
    let s = &app.store;
    let cases: Vec<i64> = v["case_ids"]
        .as_array()
        .filter(|a| !a.is_empty())
        .ok_or(bad("case_ids: expected a list of case ids"))?
        .iter()
        .map(|c| {
            let id = c.as_i64().ok_or(bad(format!("case_ids: bad id {c}")))?;
            let cfg = s.case(id).map_err(internal)?.ok_or(missing("case"))?;
            if let Some(e) = cfg["error"].as_str() {
                return Err(bad(format!(
                    "case #{id} has no statics, so nothing can float on it: {e}"
                )));
            }
            Ok(id)
        })
        .collect::<Result<_, Fail>>()?;
    let asked: Vec<StudyParams> = v["studies"]
        .as_array()
        .filter(|a| !a.is_empty())
        .ok_or(bad("studies: expected a list of studies"))?
        .iter()
        .enumerate()
        .map(|(i, r)| {
            serde_json::from_value::<StudyParams>(r.clone())
                .map_err(|e| e.to_string())
                .and_then(StudyParams::canonical)
                .map_err(|e| bad(format!("study {}: {e}", i + 1)))
        })
        .collect::<Result<_, _>>()?;
    // The studies to make on each case: those asked for, each study in
    // waves after the calm-water study it needs (`implied` when not asked).
    let mut studies: Vec<(StudyParams, bool)> = Vec::new();
    for r in &asked {
        if r.waves.is_some() {
            let grid = asked
                .iter()
                .find(|o| o.waves.is_none() && o.attitude_key() == r.attitude_key())
                .map_or(640, |o| o.grid);
            let calm = r.calm(grid);
            if !studies
                .iter()
                .any(|(o, _)| o.waves.is_none() && o.attitude_key() == calm.attitude_key())
            {
                studies.push((calm, !asked.contains(&r.calm(grid))));
            }
        }
        if !studies.iter().any(|(o, _)| o == r) {
            studies.push((r.clone(), false));
        }
    }
    let total = cases.len() * studies.len();
    if total > 400 {
        return Err(bad(format!(
            "that is {total} studies; ask for at most 400 at a time"
        )));
    }

    if v["dry_run"].as_bool() == Some(true) {
        let mut rs = Vec::new();
        for &c in &cases {
            for (r, implied) in &studies {
                let found = s.find_study(c, r).map_err(internal)?;
                rs.push(json!({
                    "case_id": c,
                    "params": r,
                    "kind": r.kind(),
                    "implied": implied,
                    "existing": found.map(|(id, status, stale)| json!({
                        "id": id, "status": status, "stale": stale,
                    })),
                }));
            }
        }
        return Ok(Reply::Json(json!({ "studies": rs })));
    }

    let by = by(who, v["by"].as_str().unwrap_or(""));
    let priority = v["priority"].as_i64().unwrap_or(0);
    let mut rs = Vec::new();
    for &cid in &cases {
        for (r, implied) in &studies {
            let (rid, created) = s.add_study(cid, r, priority, by).map_err(internal)?;
            rs.push(json!({ "id": rid, "created": created, "case_id": cid, "implied": implied }));
        }
    }
    app.worker.poke();
    Ok(Reply::Json(json!({ "studies": rs })))
}

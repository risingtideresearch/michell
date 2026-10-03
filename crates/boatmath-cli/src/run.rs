//! Running studies. A study whose result is already in the stream (or the
//! cache) isn't computed again. A study in waves is held at the attitude of
//! the calm-water study at its speed (as the web app's queue does), so
//! those calm-water studies run first, and join the stream if they weren't
//! in it.
//!
//! Each calm-water result is written with its sections and its field just
//! before it, so the stream stays whole.

use crate::cache::Cache;
use crate::records::{case_source, expect, sections, study_record, CaseSource};
use crate::stream::{short, Out, Stream};
use boatmath::params::{default_grid, StudyParams};
use boatmath::SOLVER_VERSION;
use boatmath::{
    flow_with_progress, wave_scalars, waves_with_progress, FlowRequest, Progress, ThrustLine,
};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};
use std::io::Stdout;
use std::sync::Mutex;
use std::time::Instant;

pub struct Options {
    pub jobs: usize,
    pub quiet: bool,
}

/// Run `f` on each item (with its index) on up to `jobs` threads.
fn pool<T: Sync>(items: &[T], jobs: usize, f: impl Fn(usize, &T) + Sync) {
    let next = Mutex::new(0usize);
    std::thread::scope(|s| {
        for _ in 0..jobs.clamp(1, items.len().max(1)) {
            s.spawn(|| loop {
                let i = {
                    let mut n = next.lock().unwrap();
                    *n += 1;
                    *n - 1
                };
                match items.get(i) {
                    Some(item) => f(i, item),
                    None => break,
                }
            });
        }
    });
}

fn params(study: &Value) -> Result<StudyParams, String> {
    serde_json::from_value(study["params"].clone()).map_err(|e| format!("study params: {e}"))
}

/// The calm-water study a study in waves is held at the attitude of.
fn calm_of(study: &Value) -> Result<Option<Value>, String> {
    let p = params(study)?;
    Ok(p.waves
        .is_some()
        .then(|| study_record(study["case"].as_str().unwrap_or(""), p.calm(default_grid()))))
}

/// What's known as the run goes: results computed (or found) so far, and
/// each case's attitude at rest.
#[derive(Default)]
struct Known {
    results: HashMap<String, Value>,
    at_rest: HashMap<String, (f64, f64)>,
}

/// A result already in hand: in the stream, or in the cache with its field
/// and sections. Returns the records to write, the result last.
fn found(s: &Stream, cache: &Cache, id: &str) -> Option<Vec<Value>> {
    if let Some(r) = s
        .get("result", id)
        .filter(|r| r["solver_version"] == SOLVER_VERSION)
    {
        return Some(vec![r.clone()]);
    }
    let r = cache.get("result", id)?;
    let mut out = Vec::new();
    for key in ["sections", "field"] {
        if let Some(pid) = r[key].as_str() {
            out.push(cache.get(key, pid)?);
        }
    }
    out.push(r);
    Some(out)
}

/// Run the studies in the stream, writing it through and then each result
/// (after the records it refers to) as it's ready. Returns how many failed.
pub fn run(
    s: &Stream,
    out: &mut Out<Stdout>,
    cache: &Cache,
    opts: &Options,
) -> Result<usize, String> {
    out.pass(s)?;
    // Results are written in study order, each as soon as every one before
    // it is done, so a stream comes out the same however the jobs fall.
    struct Queue<'o> {
        out: &'o mut Out<Stdout>,
        next: usize,
        ready: HashMap<usize, Vec<Value>>,
    }
    let queue = Mutex::new(Queue {
        out,
        next: 0,
        ready: HashMap::new(),
    });
    let emit = |rs: &[Value]| -> Result<(), String> {
        let mut q = queue.lock().unwrap();
        for r in rs {
            q.out.emit(r)?;
        }
        Ok(())
    };
    let finish = |i: usize, rs: Vec<Value>| -> Result<(), String> {
        let mut q = queue.lock().unwrap();
        q.ready.insert(i, rs);
        loop {
            let n = q.next;
            let Some(rs) = q.ready.remove(&n) else { break };
            for r in &rs {
                q.out.emit(r)?;
            }
            q.next += 1;
        }
        Ok(())
    };
    let known = Mutex::new(Known::default());
    let failed = Mutex::new(0usize);

    // First the calm-water studies, those asked for and those the studies
    // in waves are held about; then the studies in waves.
    let mut seen: HashSet<String> = HashSet::new();
    let mut first: Vec<Value> = Vec::new();
    let mut later: Vec<Value> = Vec::new();
    for study in s.of_type("study") {
        let id = study["id"].as_str().unwrap_or("").to_string();
        if !seen.insert(id.clone()) {
            continue;
        }
        match calm_of(study)? {
            None => first.push(study.clone()),
            Some(calm) => {
                let cid = calm["id"].as_str().unwrap_or("").to_string();
                if seen.insert(cid) {
                    // Into the stream, before its result.
                    emit(std::slice::from_ref(&calm))?;
                    first.push(calm);
                }
                later.push(study.clone());
            }
        }
    }
    let attempt = |seq: usize, study: &Value| {
        let id = study["id"].as_str().unwrap_or("");
        let records = match found(s, cache, id) {
            Some(rs) => Ok(rs),
            None => compute(s, cache, &known, study, opts.quiet),
        };
        let rs = match records {
            Ok(rs) => {
                if let Some(r) = rs.last() {
                    known
                        .lock()
                        .unwrap()
                        .results
                        .insert(id.to_string(), r.clone());
                }
                rs
            }
            Err(e) => {
                eprintln!("boatmath: study {}: {e}", short(id));
                *failed.lock().unwrap() += 1;
                Vec::new()
            }
        };
        if let Err(e) = finish(seq, rs) {
            eprintln!("boatmath: {e}");
        }
    };
    pool(&first, opts.jobs, |i, st| attempt(i, st));
    let n = first.len();
    pool(&later, opts.jobs, |i, st| attempt(n + i, st));
    Ok(failed.into_inner().unwrap())
}

/// Each hull's entry less its display mesh.
fn without_meshes(hulls: &Value) -> Value {
    let mut h = hulls.clone();
    for e in h.as_array_mut().into_iter().flatten() {
        if let Some(o) = e.as_object_mut() {
            o.remove("mesh");
        }
    }
    h
}

fn pair(v: &Value, a: &str, b: &str) -> Option<(f64, f64)> {
    Some((v[a].as_f64()?, v[b].as_f64()?))
}

/// A case's attitude at rest: from its statics in the stream, or solved.
fn at_rest(s: &Stream, known: &Mutex<Known>, case: &Value) -> Result<(f64, f64), String> {
    let id = case["id"].as_str().unwrap_or("");
    if let Some(&a) = known.lock().unwrap().at_rest.get(id) {
        return Ok(a);
    }
    let a = match s
        .get("statics", id)
        .and_then(|st| pair(&st["at_rest"], "sinkage", "trim_rad"))
    {
        Some(a) => a,
        None => {
            let CaseSource { params, src, .. } = case_source(s, case)?;
            boatmath::at_rest(&src.file_name, src.bytes, &src.import, &params)?
        }
    };
    known.lock().unwrap().at_rest.insert(id.to_string(), a);
    Ok(a)
}

/// Compute one study: the records it makes, its result last.
fn compute(
    s: &Stream,
    cache: &Cache,
    known: &Mutex<Known>,
    study: &Value,
    quiet: bool,
) -> Result<Vec<Value>, String> {
    let id = expect(study, "study")?;
    let p = params(study)?;
    let case = s.follow(study, "case")?;
    let CaseSource {
        hull,
        params: cp,
        src,
        thrusts,
        ..
    } = case_source(s, case)?;

    let t0 = Instant::now();
    let label = format!(
        "{} {} Fn {}{}",
        short(id),
        hull["name"].as_str().unwrap_or(""),
        p.froude,
        p.waves
            .as_ref()
            .map_or(String::new(), |w| format!(" waves {}°", w.heading))
    );
    let mut last = "";
    let mut report = |pr: &Progress| {
        if !quiet && pr.stage != last {
            eprintln!(
                "boatmath: {label}: {} ({:.0}%)",
                pr.stage.to_lowercase(),
                100.0 * pr.fraction
            );
            last = pr.stage;
        }
        true
    };

    let mut records = Vec::new();
    let mut r = match &p.waves {
        None => {
            let hold = if p.dynamic {
                None
            } else {
                Some(at_rest(s, known, case)?)
            };
            // With a drive and floating at speed, the attitude is the
            // self-propelled one: towed first, for the resistance the
            // thrust must meet, then with that thrust along the drives'
            // line. (Held at rest, the attitude is given.)
            let angle = thrusts.first().map_or(0.0, |t| t.angle);
            let mut propulsion = None;
            let thrust = if p.dynamic && !thrusts.is_empty() {
                let towed = FlowRequest {
                    cut: src.import,
                    case: cp.clone(),
                    study: p.clone(),
                    warm: None,
                    hold: None,
                    thrust: None,
                    field: false,
                };
                let v0 =
                    flow_with_progress(&src.file_name, src.bytes.clone(), &towed, &mut report)?;
                let rt0 = v0["forces"]["rt"]
                    .as_f64()
                    .ok_or("the towed pass has no R_t")?;
                let t = rt0 / angle.cos();
                propulsion = Some(json!({
                    "thrust": t,
                    "shaft_angle_deg": angle.to_degrees(),
                    "drives": thrusts,
                    "towed": {
                        "rt": rt0,
                        "sinkage": v0["forces"]["sinkage"],
                        "trim_deg": v0["forces"]["trim_deg"],
                    },
                }));
                Some(ThrustLine {
                    thrust: t,
                    at: thrusts.iter().map(|d| (d.x, d.z)).collect(),
                    angle,
                })
            } else {
                None
            };
            let req = FlowRequest {
                cut: src.import,
                case: cp,
                study: p.clone(),
                warm: None,
                hold,
                thrust,
                field: true,
            };
            let v = flow_with_progress(&src.file_name, src.bytes, &req, &mut report)?;
            // The meshes are for drawing; the pictures re-pose the hull's
            // geometry at the attitude instead of keeping them.
            let field = json!({
                "type": "field",
                "id": id,
                "study": id,
                "hulls": without_meshes(&v["hulls"]),
                "surface": v["surface"],
                "solver_version": SOLVER_VERSION,
            });
            let sec = match pair(&v["forces"], "sinkage", "trim_rad") {
                Some(a) => Some(sections(s, case, a, cache)?),
                None => None,
            };
            let r = json!({
                "kind": "calm",
                "froude": v["froude"],
                "speed": v["speed"],
                "transverse_wavelength": v["transverse_wavelength"],
                "seconds": v["seconds"],
                "forces": v["forces"],
                "self_propelled": propulsion,
                "field": id,
                "sections": sec.as_ref().map(|x| x["id"].clone()),
            });
            if let Some(sec) = sec {
                records.push(sec);
            }
            cache.put(&field);
            records.push(field);
            r
        }
        Some(_) => {
            let calm = calm_of(study)?.expect("a study in waves");
            let calm_id = calm["id"].as_str().unwrap_or("");
            let cr = known
                .lock()
                .unwrap()
                .results
                .get(calm_id)
                .cloned()
                .ok_or_else(|| format!("its calm-water study {} has no result", short(calm_id)))?;
            let hold = pair(&cr["forces"], "sinkage", "trim_rad")
                .ok_or("its calm-water study has no attitude")?;
            let v = waves_with_progress(
                &src.file_name,
                src.bytes,
                &src.import,
                &cp,
                &p,
                hold,
                &mut report,
            )?;
            json!({
                "kind": "waves",
                "froude": v["froude"],
                "speed": v["speed"],
                "seconds": v["seconds"],
                "attitude": v["attitude"],
                "calm": calm_id,
                "sections": cr["sections"],
                "peaks": wave_scalars(&v),
                "seakeeping": v["seakeeping"],
            })
        }
    };
    let o = r.as_object_mut().expect("an object");
    o.insert("type".into(), json!("result"));
    o.insert("id".into(), json!(id));
    o.insert("study".into(), json!(id));
    o.insert("solver_version".into(), json!(SOLVER_VERSION));
    cache.put(&r);
    records.push(r);
    if !quiet {
        eprintln!(
            "boatmath: {label}: done in {:.1} s",
            t0.elapsed().as_secs_f64()
        );
    }
    Ok(records)
}

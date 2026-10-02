//! Running studies. A result is saved under its study's id; one already
//! there from this solver is reused. A study in waves is held at the
//! attitude of the calm-water study at its speed (as the web app's queue
//! does), so those calm-water studies run first.

use crate::records::{expect, hull_source, study_record};
use crate::store::{short, Store};
use boatmath::params::{default_grid, CaseParams, StudyParams};
use boatmath::SOLVER_VERSION;
use boatmath::{flow_with_progress, wave_scalars, waves_with_progress, FlowRequest, Progress};
use serde_json::{json, Value};
use std::collections::HashSet;
use std::io::Write;
use std::sync::Mutex;
use std::time::Instant;

pub struct Options {
    pub jobs: usize,
    pub force: bool,
    pub quiet: bool,
}

fn fresh(store: &Store, id: &str) -> Result<Option<Value>, String> {
    Ok(store
        .get("result", id)?
        .filter(|r| r["solver_version"] == SOLVER_VERSION))
}

/// Run `f` on each item on up to `jobs` threads.
fn pool<T: Sync>(items: &[T], jobs: usize, f: impl Fn(&T) + Sync) {
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
                    Some(item) => f(item),
                    None => break,
                }
            });
        }
    });
}

/// The calm-water study a study in waves is held at the attitude of.
fn calm_of(study: &Value) -> Result<Option<Value>, String> {
    let p: StudyParams = serde_json::from_value(study["params"].clone())
        .map_err(|e| format!("study params: {e}"))?;
    Ok(p.waves
        .is_some()
        .then(|| study_record(study["case"].as_str().unwrap_or(""), p.calm(default_grid()))))
}

/// Run the studies, writing each one's result to stdout as it is ready.
/// Returns how many failed.
pub fn run(store: &Store, studies: &[Value], opts: &Options) -> Result<usize, String> {
    for s in studies {
        expect(s, "study")?;
    }
    let out = Mutex::new(std::io::stdout());
    let emit = |r: &Value| {
        let mut o = out.lock().unwrap();
        let _ = writeln!(o, "{r}");
        let _ = o.flush();
    };
    let failed = Mutex::new(0usize);
    let attempt = |s: &Value, show: bool| {
        let id = s["id"].as_str().unwrap_or("");
        match compute(store, s, opts.quiet) {
            Ok(r) => {
                if show {
                    emit(&r)
                }
            }
            Err(e) => {
                eprintln!("boatmath: study {}: {e}", short(id));
                *failed.lock().unwrap() += 1;
            }
        }
    };

    // First the calm-water studies: those asked for, and those the studies
    // in waves are held about.
    let mut wanted: HashSet<String> = HashSet::new();
    let mut first: Vec<(Value, bool)> = Vec::new();
    let mut later: Vec<Value> = Vec::new();
    for s in studies {
        let id = s["id"].as_str().unwrap_or("").to_string();
        if !wanted.insert(id.clone()) {
            continue;
        }
        if !opts.force {
            if let Some(r) = fresh(store, &id)? {
                emit(&r);
                continue;
            }
        }
        match calm_of(s)? {
            None => first.push((s.clone(), true)),
            Some(calm) => {
                let cid = calm["id"].as_str().unwrap_or("").to_string();
                if !studies.iter().any(|t| t["id"] == calm["id"])
                    && wanted.insert(cid.clone())
                    && fresh(store, &cid)?.is_none()
                {
                    store.put(&calm)?;
                    first.push((calm, false));
                }
                later.push(s.clone());
            }
        }
    }
    pool(&first, opts.jobs, |(s, show)| attempt(s, *show));
    pool(&later, opts.jobs, |s| attempt(s, true));
    Ok(failed.into_inner().unwrap())
}

/// Compute one study and save its result.
fn compute(store: &Store, study: &Value, quiet: bool) -> Result<Value, String> {
    let id = study["id"].as_str().unwrap_or("");
    let p: StudyParams = serde_json::from_value(study["params"].clone())
        .map_err(|e| format!("study params: {e}"))?;
    let case = store.need("case", study["case"].as_str().unwrap_or(""))?;
    if let Some(e) = case["error"].as_str() {
        return Err(format!("its case failed: {e}"));
    }
    let hull = store.need("hull", case["hull"].as_str().unwrap_or(""))?;
    let cp: CaseParams =
        serde_json::from_value(case["params"].clone()).map_err(|e| format!("case params: {e}"))?;
    let src = hull_source(&hull)?;

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

    let pair = |v: &Value, a: &str, b: &str| -> Option<(f64, f64)> {
        Some((v[a].as_f64()?, v[b].as_f64()?))
    };
    let mut r = match &p.waves {
        None => {
            let hold = if p.dynamic {
                None
            } else {
                Some(
                    pair(&case["statics"]["at_rest"], "sinkage", "trim_rad")
                        .ok_or("the case has no attitude at rest")?,
                )
            };
            let req = FlowRequest {
                cut: src.import,
                case: cp,
                study: p.clone(),
                warm: None,
                hold,
            };
            let v = flow_with_progress(&src.file_name, src.bytes, &req, &mut report)?;
            json!({
                "kind": "calm",
                "froude": v["froude"],
                "speed": v["speed"],
                "transverse_wavelength": v["transverse_wavelength"],
                "seconds": v["seconds"],
                "forces": v["forces"],
                "field": store.put_json_blob(&json!({ "hulls": v["hulls"], "surface": v["surface"] }))?,
            })
        }
        Some(_) => {
            let calm = calm_of(study)?.expect("a study in waves");
            let calm_id = calm["id"].as_str().unwrap_or("");
            let cr = fresh(store, calm_id)?
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
                "peaks": wave_scalars(&v),
                "seakeeping": v["seakeeping"],
                "field": store.put_json_blob(&json!({ "hulls": v["hulls"] }))?,
            })
        }
    };
    let o = r.as_object_mut().expect("an object");
    o.insert("type".into(), json!("result"));
    o.insert("id".into(), json!(id));
    o.insert("study".into(), json!(id));
    o.insert("solver_version".into(), json!(SOLVER_VERSION));
    store.put(&r)?;
    if !quiet {
        eprintln!(
            "boatmath: {label}: done in {:.1} s",
            t0.elapsed().as_secs_f64()
        );
    }
    Ok(r)
}

//! The worker: takes queued studies one at a time (the solver already uses
//! every core), runs them, and saves their results. The running study's
//! progress is kept here, in memory, for the queue page; a cancel request
//! for it is answered at its next progress report.

use crate::params::StudyParams;
use crate::store::{status, Claimed, NewResult, Store};
use crate::{flow_with_progress, waves_with_progress, FlowRequest, Progress, CANCELLED};
use serde_json::{json, Value};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant};

pub struct Worker {
    store: Arc<Store>,
    running: Mutex<Option<Running>>,
    wake: Condvar,
}

struct Running {
    study_id: i64,
    started: Instant,
    progress: Value,
    cancel: bool,
}

impl Worker {
    pub fn new(store: Arc<Store>) -> Arc<Worker> {
        Arc::new(Worker {
            store,
            running: Mutex::new(None),
            wake: Condvar::new(),
        })
    }

    /// Work through the queue on a thread of its own, forever.
    pub fn spawn(self: &Arc<Self>) -> std::thread::JoinHandle<()> {
        let w = Arc::clone(self);
        std::thread::spawn(move || loop {
            match w.run_next() {
                Ok(true) => {}
                Ok(false) => {
                    // Nothing to take: sleep until a study is added (or a
                    // while, in case one was added some other way).
                    let g = w.lock();
                    let _ = w.wake.wait_timeout(g, Duration::from_secs(10));
                }
                Err(e) => {
                    eprintln!("worker: {e}");
                    std::thread::sleep(Duration::from_secs(5));
                }
            }
        })
    }

    /// Say that studies have been queued.
    pub fn poke(&self) {
        self.wake.notify_all();
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, Option<Running>> {
        self.running.lock().unwrap_or_else(|p| p.into_inner())
    }

    /// The running study and how far it has got, if any.
    pub fn status(&self) -> Value {
        match &*self.lock() {
            None => Value::Null,
            Some(r) => json!({
                "study_id": r.study_id,
                "seconds": r.started.elapsed().as_secs_f64(),
                "progress": r.progress,
                "cancelling": r.cancel,
            }),
        }
    }

    /// Ask the running study to stop; `false` if it is not the one running.
    pub fn cancel(&self, study_id: i64) -> bool {
        match &mut *self.lock() {
            Some(r) if r.study_id == study_id => {
                r.cancel = true;
                true
            }
            _ => false,
        }
    }

    /// Run the next study that can go, if there is one, to its end.
    pub fn run_next(&self) -> Result<bool, String> {
        let Some(c) = self.store.claim()? else {
            return Ok(false);
        };
        let id = c.id;
        *self.lock() = Some(Running {
            study_id: id,
            started: Instant::now(),
            progress: Value::Null,
            cancel: false,
        });
        let t0 = Instant::now();
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| self.compute(&c)))
            .unwrap_or_else(|p| {
                let what = p
                    .downcast_ref::<&str>()
                    .map(|s| s.to_string())
                    .or_else(|| p.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "unknown".into());
                Err(format!("the solver panicked: {what}"))
            });
        // A cancel that came after the last progress report (the resistance
        // integral and the meshes report none) still stands.
        let cancelled = self.lock().take().is_some_and(|r| r.cancel);
        let outcome = match outcome {
            Ok(_) if cancelled => Err(CANCELLED.to_string()),
            o => o,
        };
        let secs = t0.elapsed().as_secs_f64();
        match outcome {
            Ok(result) => {
                eprintln!("study {id}: done ({secs:.1} s)");
                self.store.finish(id, &result)?;
            }
            Err(e) if e == CANCELLED => {
                eprintln!("study {id}: cancelled after {secs:.1} s");
                self.store.stop(id, status::CANCELLED, None)?;
            }
            Err(e) => {
                eprintln!("study {id}: {e}");
                self.store.stop(id, status::FAILED, Some(&e))?;
            }
        }
        Ok(true)
    }

    /// The progress callback: record where the study has got, and whether to
    /// go on.
    fn report(&self, t0: Instant, p: &Progress) -> bool {
        match &mut *self.lock() {
            Some(r) => {
                r.progress = json!({
                    "stage": p.stage,
                    "step": p.step,
                    "steps": p.steps,
                    "detail": p.detail,
                    "fraction": p.fraction,
                    "seconds": t0.elapsed().as_secs_f64(),
                });
                !r.cancel
            }
            None => true,
        }
    }

    fn compute(&self, c: &Claimed) -> Result<NewResult, String> {
        let bytes = self.store.blob(&c.hull.file_blob)?;
        let t0 = Instant::now();
        let mut report = |p: &Progress| self.report(t0, p);
        let p = &c.params;
        match &c.calm {
            None => {
                let (warm, warm_from) = match self.warm_start(c.case_id, p)? {
                    Some((w, id)) => (Some(w), Some(id)),
                    None => (None, None),
                };
                let hold = if p.dynamic {
                    None
                } else {
                    Some(
                        c.at_rest
                            .ok_or("the case has no attitude at rest (its statics failed)")?,
                    )
                };
                let req = FlowRequest {
                    cut: c.hull.import,
                    case: c.case.clone(),
                    study: p.clone(),
                    warm,
                    hold,
                };
                let v = flow_with_progress(&c.hull.file_name, bytes, &req, &mut report)?;
                let mut scalars = v["forces"].clone();
                for k in ["froude", "speed", "transverse_wavelength"] {
                    scalars[k] = v[k].clone();
                }
                Ok(NewResult {
                    scalars,
                    seconds: v["seconds"].as_f64().unwrap_or(0.0),
                    field: v,
                    warm_from,
                    attitude_from: None,
                })
            }
            Some((calm_id, hold)) => {
                let v = waves_with_progress(
                    &c.hull.file_name,
                    bytes,
                    &c.hull.import,
                    &c.case,
                    p,
                    *hold,
                    &mut report,
                )?;
                Ok(NewResult {
                    scalars: wave_scalars(&v),
                    seconds: v["seconds"].as_f64().unwrap_or(0.0),
                    field: v,
                    warm_from: None,
                    attitude_from: Some(*calm_id),
                })
            }
        }
    }

    /// The equilibrium of the nearest speed already solved on this
    /// case otherwise alike, to start the solve from — and which
    /// study it was.
    fn warm_start(&self, case_id: i64, p: &StudyParams) -> Result<Option<WarmStart>, String> {
        if !p.dynamic {
            return Ok(None);
        }
        Ok(self
            .store
            .solved_neighbours(case_id, p)?
            .into_iter()
            .min_by(|a, b| (a.1 - p.froude).abs().total_cmp(&(b.1 - p.froude).abs()))
            .map(|(id, _, z, t)| ((z, t), id)))
    }
}

/// Where to start an equilibrium solve, `(sinkage, trim)`, and the study it
/// was solved for.
type WarmStart = ((f64, f64), i64);

/// A response measure: its name, and how to read it off a wavelength's point.
type Measure<'a> = (&'a str, &'a dyn Fn(&Value) -> Option<f64>);

/// What a study in waves is plotted by: each response's peak over the
/// wavelengths and where it is, the added resistance's, the roll stability
/// at the attitude, and the irregular sea's statistics.
fn wave_scalars(v: &Value) -> Value {
    let sk = &v["seakeeping"];
    let h = &sk["headings"][0];
    let pts = h["points"].as_array().cloned().unwrap_or_default();
    let abs = |z: &Value| {
        z.as_array().map(|a| {
            a[0].as_f64()
                .unwrap_or(0.0)
                .hypot(a[1].as_f64().unwrap_or(0.0))
        })
    };
    let peak = |f: &dyn Fn(&Value) -> Option<f64>| {
        pts.iter()
            .filter_map(|p| Some((f(p)?, p["lambda"].as_f64()?)))
            .fold(None, |best: Option<(f64, f64)>, (y, l)| match best {
                Some((b, _)) if b >= y => best,
                _ => Some((y, l)),
            })
    };
    let mut out = json!({
        "froude": v["froude"],
        "speed": v["speed"],
        "heading": h["heading"],
        "sinkage": v["attitude"]["sinkage"],
        "trim_rad": v["attitude"]["trim_rad"],
        "trim_deg": v["attitude"]["trim_deg"],
        "gm_t": sk["gm_t"],
        "roll_period": sk["roll_period"],
        "failed_points": pts.iter().filter(|p| !p["error"].is_null()).count(),
    });
    let heave = |p: &Value| abs(&p["heave"]);
    let pitch = |p: &Value| abs(&p["pitch"]);
    let roll = |p: &Value| abs(&p["roll"]);
    let added = |p: &Value| p["raw_gb"].as_f64();
    let measures: [Measure; 4] = [
        ("heave", &heave),
        ("pitch", &pitch),
        ("roll", &roll),
        ("added_resistance", &added),
    ];
    for (k, f) in measures {
        if let Some((y, l)) = peak(f) {
            out[format!("{k}_peak")] = json!(y);
            out[format!("{k}_peak_lambda")] = json!(l);
        }
    }
    let sea = &h["sea"];
    if sea.is_object() && sea["error"].is_null() {
        for k in [
            "heave",
            "pitch_deg",
            "accel_bow",
            "accel_lcg",
            "raw_gb",
            "raw_far",
        ] {
            out[format!("sea_{k}")] = sea[k].clone();
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::CaseParams;
    use crate::store::{NewHull, StudyFilter};
    use crate::LoftRequest;

    fn temp_store(tag: &str) -> Arc<Store> {
        let dir = std::env::temp_dir().join(format!(
            "boatmath-queue-{tag}-{}-{}",
            std::process::id(),
            crate::store::now()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        Arc::new(Store::open(&dir).unwrap())
    }

    /// A Wigley uploaded, a case made, studies asked for (one twice),
    /// run through the queue and their results kept: the second request
    /// finds the first study, a cancelled study does not run, and a study in waves
    /// waits for its calm-water study and is then taken about its attitude.
    #[test]
    fn studies_go_through_the_queue_and_are_kept() {
        let store = temp_store("e2e");
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let bytes = michell_geometry::iges::write(&surfaces, "wigley")
            .unwrap()
            .into_bytes();
        let import = LoftRequest::default();
        let sections = crate::loft("w.igs", bytes.clone(), &import).unwrap();
        let geometry = crate::geometry("w.igs", bytes.clone(), &import).unwrap();
        let new = NewHull {
            name: "wigley",
            notes: "",
            uploaded_by: "test",
            file_name: "w.igs",
            bytes: &bytes,
            import,
            parent_id: None,
            sections: &sections,
            geometry: &geometry,
        };
        let (hull, created) = store.add_hull(&new).unwrap();
        assert!(created);
        assert_eq!(store.add_hull(&new).unwrap(), (hull, false));

        // This Wigley ends at its waterline, so held at rest's attitude —
        // the design one, at the design load.
        let case = CaseParams::default().canonical().unwrap();
        let statics = crate::statics("w.igs", bytes.clone(), &import, &case).unwrap();
        assert!(
            statics["at_rest"]["sinkage"].as_f64().unwrap().abs() < 1e-3,
            "{}",
            statics["at_rest"]
        );
        let (cfg, _) = store
            .add_case(hull, &case, "design", "test", Ok(&statics))
            .unwrap();
        let study = |s: &str| {
            serde_json::from_str::<StudyParams>(s)
                .unwrap()
                .canonical()
                .unwrap()
        };
        let a = study(r#"{"froude": 0.35, "grid": 40, "dynamic": false}"#);
        let b = study(r#"{"froude": 0.4, "grid": 40, "dynamic": false}"#);
        let w = study(
            r#"{"froude": 0.35, "dynamic": false, "waves": {"heading": 180, "lambdas": [1, 2]}}"#,
        );
        // The study in waves first, and ahead: it waits for its calm study anyway.
        let (iw, _) = store.add_study(cfg, &w, 5, "test").unwrap();
        let (ia, new_a) = store.add_study(cfg, &a, 0, "test").unwrap();
        let (ib, _) = store.add_study(cfg, &b, 0, "test").unwrap();
        assert!(new_a);
        assert_eq!(store.add_study(cfg, &a, 0, "test").unwrap(), (ia, false));
        assert_eq!(store.cancel_queued(ib).unwrap().as_deref(), Some("queued"));

        let worker = Worker::new(Arc::clone(&store));
        assert!(worker.run_next().unwrap());
        assert_eq!(
            store.study(ia).unwrap().unwrap()["status"],
            "done",
            "the calm study goes first"
        );
        assert!(worker.run_next().unwrap());
        assert!(!worker.run_next().unwrap(), "the cancelled study ran");

        let done = store
            .studies(&StudyFilter {
                case: Some(cfg),
                status: Some("done".into()),
                ..StudyFilter::default()
            })
            .unwrap();
        assert_eq!(done.len(), 2);
        let calm = store.study(ia).unwrap().unwrap();
        let r = &calm["result"];
        assert!(r["scalars"]["rw"].as_f64().unwrap() > 0.0, "{r}");
        assert_eq!(r["scalars"]["froude"], 0.35);
        assert_eq!(r["stale"], false);
        let blob = store.result_blob(ia).unwrap().unwrap();
        let field: Value = serde_json::from_slice(&store.blob(&blob).unwrap()).unwrap();
        assert_eq!(field["surface"]["nx"], 40);

        let waves = store.study(iw).unwrap().unwrap();
        let s = &waves["result"]["scalars"];
        assert_eq!(waves["result"]["attitude_from"], ia);
        assert!(s["heave_peak"].as_f64().unwrap() > 0.0, "{s}");
        assert_eq!(s["heading"], 180.0);

        // Asking again for the cancelled study puts it back on the queue.
        assert_eq!(store.add_study(cfg, &b, 0, "test").unwrap(), (ib, false));
        assert_eq!(store.study(ib).unwrap().unwrap()["status"], "queued");
    }
}

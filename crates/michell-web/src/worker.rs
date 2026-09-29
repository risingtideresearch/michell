//! The worker: takes queued cases one at a time (the solver already uses
//! every core), runs them, and saves their results. The running case's
//! progress is kept here, in memory, for the queue page; a cancel request
//! for it is answered at its next progress report.

use crate::case::CaseParams;
use crate::store::{status, Claimed, NewResult, Store};
use crate::{flow_with_progress, span_sweep_with_progress, Progress, CANCELLED};
use serde_json::{json, Value};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant};

pub struct Worker {
    store: Arc<Store>,
    running: Mutex<Option<Running>>,
    wake: Condvar,
}

struct Running {
    case_id: i64,
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
                    // Nothing queued: sleep until a case is added (or a
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

    /// Say that cases have been queued.
    pub fn poke(&self) {
        self.wake.notify_all();
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, Option<Running>> {
        self.running.lock().unwrap_or_else(|p| p.into_inner())
    }

    /// The running case and how far it has got, if any.
    pub fn status(&self) -> Value {
        match &*self.lock() {
            None => Value::Null,
            Some(r) => json!({
                "case_id": r.case_id,
                "seconds": r.started.elapsed().as_secs_f64(),
                "progress": r.progress,
                "cancelling": r.cancel,
            }),
        }
    }

    /// Ask the running case to stop; `false` if it is not the one running.
    pub fn cancel(&self, case_id: i64) -> bool {
        match &mut *self.lock() {
            Some(r) if r.case_id == case_id => {
                r.cancel = true;
                true
            }
            _ => false,
        }
    }

    /// Run the next queued case, if there is one, to its end.
    pub fn run_next(&self) -> Result<bool, String> {
        let Some(c) = self.store.claim()? else {
            return Ok(false);
        };
        let id = c.id;
        *self.lock() = Some(Running {
            case_id: id,
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
        *self.lock() = None;
        let secs = t0.elapsed().as_secs_f64();
        match outcome {
            Ok((results, warm_from)) => {
                eprintln!("case {id}: done ({secs:.1} s)");
                self.store.finish(id, &results, warm_from)?;
            }
            Err(e) if e == CANCELLED => {
                eprintln!("case {id}: cancelled after {secs:.1} s");
                self.store.stop(id, status::CANCELLED, None)?;
            }
            Err(e) => {
                eprintln!("case {id}: {e}");
                self.store.stop(id, status::FAILED, Some(&e))?;
            }
        }
        Ok(true)
    }

    /// The progress callback: record where the case has got, and whether
    /// to go on.
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

    fn compute(&self, c: &Claimed) -> Result<(Vec<NewResult>, Option<i64>), String> {
        let bytes = self.store.blob(&c.file_blob)?;
        let t0 = Instant::now();
        let mut report = |p: &Progress| self.report(t0, p);
        let p = &c.params;
        match &p.spans {
            None => {
                let start = self.warm_start(c.hull_id, p)?;
                let (warm, warm_from) = (start.map(|s| s.0), start.map(|s| s.1));
                let req = p.flow_request(c.import, warm);
                let v = flow_with_progress(&c.file_name, bytes, &req, &mut report)?;
                Ok((vec![result(json!({}), v)], warm_from))
            }
            Some(spans) => {
                let req = p.flow_request(c.import, None);
                let mut out: Vec<(usize, Value)> = Vec::new();
                span_sweep_with_progress(
                    &c.file_name,
                    bytes,
                    &req,
                    spans,
                    &mut |i, v| {
                        out.push((i, v));
                        true
                    },
                    &mut report,
                )?;
                out.sort_by_key(|(i, _)| *i);
                Ok((
                    out.into_iter()
                        .map(|(i, v)| result(json!({ "span": spans[i] }), v))
                        .collect(),
                    None,
                ))
            }
        }
    }

    /// The equilibrium of the nearest speed already solved for this case
    /// otherwise, to start the solve from — and which case it was.
    fn warm_start(&self, hull_id: i64, p: &CaseParams) -> Result<Option<WarmStart>, String> {
        if !p.dynamic {
            return Ok(None);
        }
        Ok(self
            .store
            .solved_neighbours(hull_id, p)?
            .into_iter()
            .min_by(|a, b| (a.1 - p.froude).abs().total_cmp(&(b.1 - p.froude).abs()))
            .map(|(id, _, z, t)| ((z, t), id)))
    }
}

/// Where to start an equilibrium solve, `(sinkage, trim)`, and the case it
/// was solved for.
type WarmStart = ((f64, f64), i64);

/// A computed point: the forces and speed as scalars, the whole answer as
/// its field.
fn result(point: Value, v: Value) -> NewResult {
    let mut scalars = v["forces"].clone();
    for k in ["froude", "speed", "transverse_wavelength"] {
        scalars[k] = v[k].clone();
    }
    NewResult {
        point,
        scalars,
        seconds: v["seconds"].as_f64().unwrap_or(0.0),
        field: v,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::store::{CaseFilter, NewHull};
    use crate::LoftRequest;

    fn temp_store(tag: &str) -> Arc<Store> {
        let dir = std::env::temp_dir().join(format!(
            "michell-queue-{tag}-{}-{}",
            std::process::id(),
            crate::store::now()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        Arc::new(Store::open(&dir).unwrap())
    }

    /// A Wigley uploaded, two cases asked for (one twice), run through the
    /// queue, and their results kept: the second request finds the first
    /// case, a cancelled case does not run, and the saved answer is whole.
    #[test]
    fn cases_run_through_the_queue_and_are_kept() {
        let store = temp_store("e2e");
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let bytes = michell_geometry::iges::write(&surfaces, "wigley")
            .unwrap()
            .into_bytes();
        let import = LoftRequest::default();
        let sections = crate::loft("w.igs", bytes.clone(), &import).unwrap();
        let geometry = crate::geometry("w.igs", bytes.clone(), None).unwrap();
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

        // This Wigley ends at its waterline, so at the design attitude.
        let case = |s: &str| {
            serde_json::from_str::<CaseParams>(s)
                .unwrap()
                .canonical()
                .unwrap()
        };
        let a = case(r#"{"froude": 0.35, "grid": 40, "dynamic": false}"#);
        let b = case(r#"{"froude": 0.4, "grid": 40, "dynamic": false}"#);
        let study = store.add_study("speeds", "", "test").unwrap();
        let (ia, new_a) = store.add_case(hull, &a, Some(study), 0, "test").unwrap();
        let (ib, _) = store.add_case(hull, &b, Some(study), 0, "test").unwrap();
        assert!(new_a);
        assert_eq!(
            store.add_case(hull, &a, None, 0, "test").unwrap(),
            (ia, false)
        );
        assert_eq!(store.cancel_queued(ib).unwrap().as_deref(), Some("queued"));

        let worker = Worker::new(Arc::clone(&store));
        assert!(worker.run_next().unwrap());
        assert!(!worker.run_next().unwrap(), "the cancelled case ran");

        let done = store
            .cases(&CaseFilter {
                study: Some(study),
                status: Some("done".into()),
                ..CaseFilter::default()
            })
            .unwrap();
        assert_eq!(done.len(), 1);
        let r = &done[0]["results"][0];
        assert!(r["scalars"]["rw"].as_f64().unwrap() > 0.0, "{r}");
        assert_eq!(r["scalars"]["froude"], 0.35);
        assert_eq!(r["stale"], false);
        let blob = store.result_blob(ia, 0).unwrap().unwrap();
        let field: Value = serde_json::from_slice(&store.blob(&blob).unwrap()).unwrap();
        assert_eq!(field["surface"]["nx"], 40);

        // Asking again for the cancelled case puts it back on the queue.
        assert_eq!(
            store.add_case(hull, &b, None, 0, "test").unwrap(),
            (ib, false)
        );
        assert_eq!(store.case(ib).unwrap().unwrap()["status"], "queued");
    }
}

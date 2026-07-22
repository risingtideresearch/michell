//! Runs a loft or sweep on a background thread by calling the `michell_cli`
//! library in-process — the editor is a standalone binary and does not shell
//! out to the `michell` CLI. The library reports progress through a callback
//! (a human-readable line plus an optional `(done, total)` fraction), which we
//! forward to the UI over a channel; its return value is the text the CLI would
//! have written to stdout (the loft summary table, or the sweep CSV/JSON).

use std::path::PathBuf;
use std::sync::mpsc::{channel, Receiver};
use std::thread;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobKind {
    Loft,
    Sweep,
}

/// A successful run's captured stdout (the loft summary table / sweep CSV) plus
/// the working directory it ran in, so the caller can resolve any files it
/// produced.
pub struct Success {
    pub stdout: String,
    pub cwd: PathBuf,
}

enum Update {
    Log(String),
    Progress { done: usize, total: usize },
    Finished(Result<Success, String>),
}

/// A live job. Poll it each frame; when `outcome` becomes `Some`, it is done.
pub struct Job {
    pub kind: JobKind,
    pub title: String,
    pub log: Vec<String>,
    pub progress: Option<(usize, usize)>,
    pub outcome: Option<Result<Success, String>>,
    rx: Receiver<Update>,
}

impl Job {
    /// Drain any pending updates. Returns true if anything changed (so the app
    /// can refresh), and leaves `outcome` set once the run finishes.
    pub fn poll(&mut self) -> bool {
        let mut changed = false;
        while let Ok(u) = self.rx.try_recv() {
            changed = true;
            match u {
                Update::Log(l) => self.log.push(l),
                Update::Progress { done, total } => self.progress = Some((done, total)),
                Update::Finished(r) => self.outcome = Some(r),
            }
        }
        changed
    }

    pub fn is_running(&self) -> bool {
        self.outcome.is_none()
    }

    /// A 0..1 bar fraction, if a total is known.
    pub fn fraction(&self) -> Option<f32> {
        self.progress
            .filter(|(_, total)| *total > 0)
            .map(|(done, total)| (done as f32 / total as f32).clamp(0.0, 1.0))
    }
}

/// Start a loft or sweep on a background thread, streaming progress. `args` is a
/// CLI-style argument vector: for a loft, `["loft", <source>, ...]` (the `-o`
/// output path should be absolute, since the work runs in-process with no
/// per-job working directory); for a sweep, `["sweep", <manifest.json>]`. `cwd`
/// is reported back on success so the caller can resolve produced files.
pub fn spawn(args: Vec<String>, cwd: PathBuf, kind: JobKind, title: String) -> Job {
    let (tx, rx) = channel();
    let cwd_report = cwd;

    thread::spawn(move || {
        // The library's progress callback: forward the fraction (if any) and
        // the human-readable line to the UI channel.
        let mut report = |line: &str, frac: Option<(usize, usize)>| {
            if let Some((done, total)) = frac {
                let _ = tx.send(Update::Progress { done, total });
            }
            let _ = tx.send(Update::Log(line.to_string()));
        };

        let result = match kind {
            JobKind::Loft => michell_cli::loft(&args[1..], &mut report),
            JobKind::Sweep => match args.get(1) {
                Some(path) => michell_cli::run_manifest(path, &mut report),
                None => Err("sweep: no manifest path".to_string()),
            },
        };

        let msg = result.map(|stdout| Success {
            stdout,
            cwd: cwd_report,
        });
        let _ = tx.send(Update::Finished(msg));
    });

    Job {
        kind,
        title,
        log: Vec::new(),
        progress: None,
        outcome: None,
        rx,
    }
}

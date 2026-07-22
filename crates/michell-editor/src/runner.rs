//! Runs the `michell` CLI as a subprocess (loft, sweep) and streams its
//! progress back to the UI over a channel. The CLI already prints
//! machine-parseable progress to stderr — `sweep: N point(s) …` and
//! `point k/N done` for sweeps, `lofting hull i/n` for lofts — so we parse
//! those lines into a `(done, total)` fraction rather than reimplementing any
//! of the numerics.

use std::ffi::OsString;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::mpsc::{channel, Receiver};
use std::thread;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum JobKind {
    Loft,
    Sweep,
}

/// A successful run's captured stdout (the CLI writes its result tables / CSV
/// there) plus the working directory it ran in, so the caller can resolve any
/// files it produced.
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
    /// can refresh), and leaves `outcome` set once the process exits.
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

/// Interpret a CLI stderr line as progress, if it is one.
fn parse_progress(line: &str) -> Option<(usize, usize)> {
    // "point 3/9 done" — 3 of 9 points complete.
    if let Some(rest) = line.strip_prefix("point ") {
        let frac = rest.strip_suffix(" done")?;
        let (a, b) = frac.split_once('/')?;
        return Some((a.trim().parse().ok()?, b.trim().parse().ok()?));
    }
    // "sweep: 9 point(s) x 5 speed(s)" — total known up front, none done yet.
    if let Some(rest) = line.strip_prefix("sweep: ") {
        let n = rest.split_whitespace().next()?.parse().ok()?;
        return Some((0, n));
    }
    // "lofting hull 1/3" — hull i is starting, so i-1 are complete.
    if let Some(rest) = line.strip_prefix("lofting hull ") {
        let (a, b) = rest.trim().split_once('/')?;
        let i: usize = a.trim().parse().ok()?;
        return Some((i.saturating_sub(1), b.trim().parse().ok()?));
    }
    None
}

/// Spawn `bin` with `args` in `cwd`, streaming progress. Returns immediately;
/// the process runs on a background thread.
pub fn spawn(bin: &Path, args: Vec<OsString>, cwd: &Path, kind: JobKind, title: String) -> Job {
    let (tx, rx) = channel();
    let bin = bin.to_owned();
    let cwd = cwd.to_owned();
    let cwd_report = cwd.clone();

    thread::spawn(move || {
        let mut child = match Command::new(&bin)
            .args(&args)
            .current_dir(&cwd)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
        {
            Ok(c) => c,
            Err(e) => {
                let _ = tx.send(Update::Finished(Err(format!(
                    "cannot run {}: {e}",
                    bin.display()
                ))));
                return;
            }
        };

        // Drain stdout on its own thread so a large CSV (sweep with no output
        // file) can't deadlock the pipe while we read stderr.
        let mut stdout = child.stdout.take().expect("piped stdout");
        let out_handle = thread::spawn(move || {
            let mut s = String::new();
            let _ = stdout.read_to_string(&mut s);
            s
        });

        let stderr = child.stderr.take().expect("piped stderr");
        let mut last_err = String::new();
        for line in BufReader::new(stderr).lines().map_while(Result::ok) {
            if let Some((done, total)) = parse_progress(&line) {
                let _ = tx.send(Update::Progress { done, total });
            }
            if line.starts_with("error:") {
                last_err = line.clone();
            }
            let _ = tx.send(Update::Log(line));
        }

        let status = child.wait();
        let stdout_str = out_handle.join().unwrap_or_default();
        let msg = match status {
            Ok(s) if s.success() => Ok(Success {
                stdout: stdout_str,
                cwd: cwd_report,
            }),
            Ok(s) => Err(if last_err.is_empty() {
                format!("michell exited with {s}")
            } else {
                last_err.trim_start_matches("error:").trim().to_string()
            }),
            Err(e) => Err(e.to_string()),
        };
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

/// Locate the `michell` binary: prefer the one shipped next to this editor
/// (same `target/…` dir), then fall back to the bare name on `PATH`.
pub fn default_michell_bin() -> PathBuf {
    let exe_name = if cfg!(windows) {
        "michell.exe"
    } else {
        "michell"
    };
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let sibling = dir.join(exe_name);
            if sibling.exists() {
                return sibling;
            }
        }
    }
    PathBuf::from("michell")
}

#[cfg(test)]
mod tests {
    use super::parse_progress;

    #[test]
    fn parses_sweep_point_lines() {
        assert_eq!(parse_progress("point 3/9 done"), Some((3, 9)));
        assert_eq!(
            parse_progress("sweep: 9 point(s) x 5 speed(s)"),
            Some((0, 9))
        );
        assert_eq!(
            parse_progress("sweep: 12 point(s) x 3 speed(s), equilibrium mode"),
            Some((0, 12))
        );
    }

    #[test]
    fn parses_loft_lines() {
        assert_eq!(parse_progress("lofting hull 1/3"), Some((0, 3)));
        assert_eq!(parse_progress("lofting hull 3/3"), Some((2, 3)));
    }

    #[test]
    fn ignores_other_lines() {
        assert_eq!(parse_progress("study: ama placement"), None);
        assert_eq!(parse_progress("hull h: ama.hull (centerplane 0.0)"), None);
        assert_eq!(parse_progress(""), None);
    }
}

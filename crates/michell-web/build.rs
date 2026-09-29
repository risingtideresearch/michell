//! Stamps `MICHELL_SOLVER_VERSION`: the last commit to touch the solver's
//! code, `-dirty` if that code has uncommitted changes. Only the solver's
//! code counts, so a change to the pages does not make saved results stale.
//! A build without git (a container) passes the version in the environment.

use std::process::Command;

const SOLVER: &[&str] = &[
    "../michell-geometry/src",
    "../michell/src",
    "../michell-seakeeping/src",
    "../michell-cli/src",
    "src/lib.rs",
];

fn git(args: &[&str]) -> Option<(bool, String)> {
    let out = Command::new("git").args(args).output().ok()?;
    Some((
        out.status.success(),
        String::from_utf8_lossy(&out.stdout).trim().to_string(),
    ))
}

fn main() {
    println!("cargo:rerun-if-env-changed=MICHELL_SOLVER_VERSION");
    if let Ok(v) = std::env::var("MICHELL_SOLVER_VERSION") {
        println!("cargo:rustc-env=MICHELL_SOLVER_VERSION={v}");
        return;
    }
    for p in SOLVER {
        println!("cargo:rerun-if-changed={p}");
    }
    // A commit or checkout moves HEAD and rewrites the index.
    if let Some((true, dir)) = git(&["rev-parse", "--git-dir"]) {
        for f in ["HEAD", "index"] {
            println!("cargo:rerun-if-changed={dir}/{f}");
        }
    }
    let mut args = vec!["log", "-1", "--format=%h", "--abbrev=12", "--"];
    args.extend(SOLVER);
    let version = match git(&args) {
        Some((true, h)) if !h.is_empty() => {
            let mut diff = vec!["diff", "--quiet", "HEAD", "--"];
            diff.extend(SOLVER);
            let clean = Command::new("git")
                .args(&diff)
                .status()
                .is_ok_and(|s| s.success());
            if clean {
                h
            } else {
                format!("{h}-dirty")
            }
        }
        _ => "unknown".to_string(),
    };
    println!("cargo:rustc-env=MICHELL_SOLVER_VERSION={version}");
}

//! The cache: an opt-in memo of expensive records (statics, sections,
//! results, fields, props), by type and id, from this solver. A record
//! found here is the record the computation would make, so the output is
//! the same with or without it; deleting it costs only recomputation.
//!
//! ```text
//! $BOATMATH_CACHE/<type>/ab/cdef….json
//! ```

use crate::stream::short;
use boatmath::SOLVER_VERSION;
use serde_json::Value;
use std::path::{Path, PathBuf};

#[derive(Default, Clone)]
pub struct Cache {
    dir: Option<PathBuf>,
}

fn path(dir: &Path, kind: &str, id: &str) -> PathBuf {
    let (a, b) = id.split_at(2.min(id.len()));
    dir.join(kind).join(a).join(format!("{b}.json"))
}

impl Cache {
    /// `--cache DIR`, else `$BOATMATH_CACHE`, else none.
    pub fn open(dir: Option<PathBuf>) -> Cache {
        Cache {
            dir: dir.or_else(|| std::env::var_os("BOATMATH_CACHE").map(PathBuf::from)),
        }
    }

    /// The cached record of this type and id, if it's from this solver.
    pub fn get(&self, kind: &str, id: &str) -> Option<Value> {
        let p = path(self.dir.as_ref()?, kind, id);
        let v: Value = serde_json::from_str(&std::fs::read_to_string(p).ok()?).ok()?;
        (v["solver_version"] == SOLVER_VERSION).then_some(v)
    }

    /// Keep a record, written to a temporary file and renamed into place so
    /// concurrent runs never see half of one. Failing to cache is only a
    /// warning: the computation stands.
    pub fn put(&self, r: &Value) {
        let Some(dir) = &self.dir else { return };
        let (Some(kind), Some(id)) = (r["type"].as_str(), r["id"].as_str()) else {
            return;
        };
        let p = path(dir, kind, id);
        let write = || -> std::io::Result<()> {
            std::fs::create_dir_all(p.parent().expect("a parent"))?;
            let tmp = p.with_extension(format!("tmp{}", std::process::id()));
            std::fs::write(&tmp, r.to_string())?;
            std::fs::rename(&tmp, &p)
        };
        if let Err(e) = write() {
            eprintln!("boatmath: cache: {kind} {}: {e}", short(id));
        }
    }
}

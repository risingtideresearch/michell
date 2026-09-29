//! The permanent record, in SQLite, with the bulky data — uploaded files, a
//! hull's sections and display geometry, a configuration's statics, a
//! run's fields — in a content-addressed directory of gzipped blobs beside
//! it.
//!
//! ```text
//! $DATA_DIR/michell.db
//! $DATA_DIR/blobs/ab/cdef….gz    SHA-256 of the uncompressed bytes
//! ```
//!
//! Three levels:
//!
//! - a **hull** is a file together with its import settings (the cut
//!   depends on both);
//! - a **configuration** is a hull and [`ConfigParams`] (unique on the
//!   pair), with its statics, computed when it is made;
//! - a **run** is a configuration and [`RunParams`] (unique on the pair),
//!   queued, and its **result** once done. A run in waves waits for the
//!   calm-water run it is taken about (the same `attitude_key`).
//!
//! Studies are named groups of runs asked for together. Statics and
//! results record the solver version they were computed with; one from
//! another version is stale, kept until it is computed again.

use crate::case::{hex, ConfigParams, RunParams};
use crate::LoftRequest;
use rusqlite::{params, Connection, OptionalExtension, Row};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// The solver version results are stamped with: the last commit to touch
/// the solver's code (the geometry, thin-ship, seakeeping and CLI crates,
/// and the web crate's computations), `-dirty` if they have uncommitted
/// changes. Set by `build.rs`; `MICHELL_SOLVER_VERSION` at build time
/// overrides it.
pub const SOLVER_VERSION: &str = env!("MICHELL_SOLVER_VERSION");

const SCHEMA: &str = "
CREATE TABLE IF NOT EXISTS hulls (
    id            INTEGER PRIMARY KEY,
    name          TEXT NOT NULL,
    notes         TEXT NOT NULL DEFAULT '',
    uploaded_by   TEXT NOT NULL DEFAULT '',
    file_name     TEXT NOT NULL,
    file_blob     TEXT NOT NULL,
    import        TEXT NOT NULL,
    parent_id     INTEGER REFERENCES hulls(id),
    summary       TEXT NOT NULL,
    sections_blob TEXT NOT NULL,
    geometry_blob TEXT NOT NULL,
    created_at    INTEGER NOT NULL,
    UNIQUE (file_blob, import)
);
CREATE TABLE IF NOT EXISTS configs (
    id             INTEGER PRIMARY KEY,
    hull_id        INTEGER NOT NULL REFERENCES hulls(id),
    name           TEXT NOT NULL DEFAULT '',
    notes          TEXT NOT NULL DEFAULT '',
    params         TEXT NOT NULL,
    params_hash    TEXT NOT NULL,
    summary        TEXT,
    statics_blob   TEXT,
    error          TEXT,
    solver_version TEXT NOT NULL,
    created_by     TEXT NOT NULL DEFAULT '',
    created_at     INTEGER NOT NULL,
    UNIQUE (hull_id, params_hash)
);
CREATE TABLE IF NOT EXISTS studies (
    id         INTEGER PRIMARY KEY,
    name       TEXT NOT NULL,
    notes      TEXT NOT NULL DEFAULT '',
    created_by TEXT NOT NULL DEFAULT '',
    created_at INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS runs (
    id           INTEGER PRIMARY KEY,
    config_id    INTEGER NOT NULL REFERENCES configs(id),
    kind         TEXT NOT NULL,
    params       TEXT NOT NULL,
    params_hash  TEXT NOT NULL,
    attitude_key TEXT NOT NULL,
    status       TEXT NOT NULL,
    priority     INTEGER NOT NULL DEFAULT 0,
    requested_by TEXT NOT NULL DEFAULT '',
    error        TEXT,
    created_at   INTEGER NOT NULL,
    queued_at    INTEGER,
    started_at   INTEGER,
    finished_at  INTEGER,
    UNIQUE (config_id, params_hash)
);
CREATE INDEX IF NOT EXISTS runs_queue ON runs (status, priority, id);
CREATE INDEX IF NOT EXISTS runs_attitude ON runs (config_id, attitude_key);
CREATE TABLE IF NOT EXISTS study_runs (
    study_id INTEGER NOT NULL REFERENCES studies(id),
    run_id   INTEGER NOT NULL REFERENCES runs(id),
    PRIMARY KEY (study_id, run_id)
);
CREATE TABLE IF NOT EXISTS results (
    run_id         INTEGER PRIMARY KEY REFERENCES runs(id),
    scalars        TEXT NOT NULL,
    field_blob     TEXT NOT NULL,
    solver_version TEXT NOT NULL,
    seconds        REAL NOT NULL,
    warm_from      INTEGER REFERENCES runs(id),
    attitude_from  INTEGER REFERENCES runs(id),
    created_at     INTEGER NOT NULL
);
";

/// A run's state in the queue.
pub mod status {
    pub const QUEUED: &str = "queued";
    pub const RUNNING: &str = "running";
    pub const DONE: &str = "done";
    pub const FAILED: &str = "failed";
    pub const CANCELLED: &str = "cancelled";
}

pub struct Store {
    db: Mutex<Connection>,
    blobs: PathBuf,
}

/// A new hull, as uploaded and cut.
pub struct NewHull<'a> {
    pub name: &'a str,
    pub notes: &'a str,
    pub uploaded_by: &'a str,
    pub file_name: &'a str,
    pub bytes: &'a [u8],
    pub import: LoftRequest,
    pub parent_id: Option<i64>,
    /// The cut, as [`crate::loft`] describes it.
    pub sections: &'a Value,
    /// The display geometry, as [`crate::geometry`] describes it.
    pub geometry: &'a Value,
}

/// What a computation needs of a hull: its file and how it is cut.
pub struct HullSource {
    pub file_name: String,
    pub file_blob: String,
    pub import: LoftRequest,
}

/// A run the worker has claimed, with what it needs to run it.
pub struct Claimed {
    pub id: i64,
    pub config_id: i64,
    pub params: RunParams,
    pub config: ConfigParams,
    pub hull: HullSource,
    /// The configuration's attitude at rest `(sinkage, trim)`, for a run
    /// held there.
    pub at_rest: Option<(f64, f64)>,
    /// A run in waves: its calm-water run and that run's attitude.
    pub calm: Option<(i64, (f64, f64))>,
}

/// A computed run, ready to save.
pub struct NewResult {
    pub scalars: Value,
    /// The full answer (fields and all), stored as a blob.
    pub field: Value,
    pub seconds: f64,
    pub warm_from: Option<i64>,
    pub attitude_from: Option<i64>,
}

pub fn now() -> i64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs() as i64)
}

fn err(e: impl std::fmt::Display) -> String {
    e.to_string()
}

fn parse(s: Option<String>) -> Value {
    s.and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or(Value::Null)
}

impl Store {
    /// Open (or create) the store in `dir`. Runs left running by a server
    /// that stopped go back on the queue.
    pub fn open(dir: &Path) -> Result<Store, String> {
        let blobs = dir.join("blobs");
        std::fs::create_dir_all(&blobs).map_err(|e| format!("{}: {e}", blobs.display()))?;
        let db = Connection::open(dir.join("michell.db")).map_err(err)?;
        db.execute_batch("PRAGMA journal_mode = WAL; PRAGMA foreign_keys = ON;")
            .map_err(err)?;
        let old: bool = db
            .query_row(
                "SELECT EXISTS (SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'cases')",
                [],
                |r| r.get(0),
            )
            .map_err(err)?;
        if old {
            return Err(format!(
                "{} holds cases from before configurations and runs; start a new store",
                dir.display()
            ));
        }
        db.execute_batch(SCHEMA).map_err(err)?;
        db.execute(
            "UPDATE runs SET status = ?1, started_at = NULL WHERE status = ?2",
            params![status::QUEUED, status::RUNNING],
        )
        .map_err(err)?;
        Ok(Store {
            db: Mutex::new(db),
            blobs,
        })
    }

    fn db(&self) -> std::sync::MutexGuard<'_, Connection> {
        self.db.lock().unwrap_or_else(|p| p.into_inner())
    }

    // --- blobs ---------------------------------------------------------

    fn blob_path(&self, hash: &str) -> Result<PathBuf, String> {
        if hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(format!("bad blob id {hash:?}"));
        }
        Ok(self
            .blobs
            .join(&hash[..2])
            .join(format!("{}.gz", &hash[2..])))
    }

    /// Store `bytes`, returning their SHA-256; storing the same bytes again
    /// is free.
    pub fn put_blob(&self, bytes: &[u8]) -> Result<String, String> {
        let hash = hex(&Sha256::digest(bytes));
        let path = self.blob_path(&hash)?;
        if path.exists() {
            return Ok(hash);
        }
        let dir = path.parent().expect("blob paths have a parent");
        std::fs::create_dir_all(dir).map_err(err)?;
        // Written aside and renamed, so a blob that exists is whole.
        static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let tmp = dir.join(format!(
            "{}.{}.{}.tmp",
            &hash[2..],
            std::process::id(),
            SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        let mut gz = flate2::write::GzEncoder::new(
            std::fs::File::create(&tmp).map_err(err)?,
            flate2::Compression::default(),
        );
        gz.write_all(bytes).map_err(err)?;
        gz.finish().map_err(err)?.sync_all().map_err(err)?;
        std::fs::rename(&tmp, &path).map_err(err)?;
        Ok(hash)
    }

    /// A blob as stored, gzipped: what to send with `Content-Encoding: gzip`.
    pub fn blob_gz(&self, hash: &str) -> Result<Vec<u8>, String> {
        std::fs::read(self.blob_path(hash)?).map_err(|e| format!("blob {hash}: {e}"))
    }

    pub fn blob(&self, hash: &str) -> Result<Vec<u8>, String> {
        let mut out = Vec::new();
        flate2::read::GzDecoder::new(&self.blob_gz(hash)?[..])
            .read_to_end(&mut out)
            .map_err(|e| format!("blob {hash}: {e}"))?;
        Ok(out)
    }

    fn put_json(&self, v: &Value) -> Result<String, String> {
        self.put_blob(v.to_string().as_bytes())
    }

    // --- hulls ---------------------------------------------------------

    /// Save a hull; the same file with the same import settings is the hull
    /// already saved, returned as `(id, false)`.
    pub fn add_hull(&self, h: &NewHull) -> Result<(i64, bool), String> {
        let file_blob = self.put_blob(h.bytes)?;
        let import = serde_json::to_string(&h.import).map_err(err)?;
        if let Some(id) = self
            .db()
            .query_row(
                "SELECT id FROM hulls WHERE file_blob = ?1 AND import = ?2",
                params![file_blob, import],
                |r| r.get(0),
            )
            .optional()
            .map_err(err)?
        {
            return Ok((id, false));
        }
        let summary = hull_summary(h.sections);
        let sections_blob = self.put_json(h.sections)?;
        let geometry_blob = self.put_json(h.geometry)?;
        let db = self.db();
        db.execute(
            "INSERT INTO hulls (name, notes, uploaded_by, file_name, file_blob, import,
                                parent_id, summary, sections_blob, geometry_blob, created_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            params![
                h.name,
                h.notes,
                h.uploaded_by,
                h.file_name,
                file_blob,
                import,
                h.parent_id,
                summary.to_string(),
                sections_blob,
                geometry_blob,
                now()
            ],
        )
        .map_err(err)?;
        Ok((db.last_insert_rowid(), true))
    }

    /// Every hull, newest first, with its configuration count and its runs'
    /// counts by status.
    pub fn hulls(&self) -> Result<Vec<Value>, String> {
        let db = self.db();
        let mut q = db
            .prepare(&format!("SELECT {HULL_COLS} FROM hulls h ORDER BY id DESC"))
            .map_err(err)?;
        let rows = q.query_map([], hull_row).map_err(err)?;
        rows.collect::<Result<_, _>>().map_err(err)
    }

    pub fn hull(&self, id: i64) -> Result<Option<Value>, String> {
        self.db()
            .query_row(
                &format!("SELECT {HULL_COLS} FROM hulls h WHERE id = ?1"),
                [id],
                hull_row,
            )
            .optional()
            .map_err(err)
    }

    /// What a computation on hull `id` needs: its file and import settings.
    pub fn hull_source(&self, id: i64) -> Result<Option<HullSource>, String> {
        let row: Option<(String, String, String)> = self
            .db()
            .query_row(
                "SELECT file_name, file_blob, import FROM hulls WHERE id = ?1",
                [id],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
            )
            .optional()
            .map_err(err)?;
        row.map(|(file_name, file_blob, import)| {
            Ok(HullSource {
                file_name,
                file_blob,
                import: serde_json::from_str(&import).map_err(err)?,
            })
        })
        .transpose()
    }

    /// A hull's blob (`file`, `sections` or `geometry`), by the hull's id.
    pub fn hull_blob(&self, id: i64, which: &str) -> Result<Option<String>, String> {
        let col = match which {
            "file" => "file_blob",
            "sections" => "sections_blob",
            "geometry" => "geometry_blob",
            _ => return Ok(None),
        };
        self.db()
            .query_row(
                &format!("SELECT {col} FROM hulls WHERE id = ?1"),
                [id],
                |r| r.get(0),
            )
            .optional()
            .map_err(err)
    }

    pub fn update_hull(
        &self,
        id: i64,
        name: Option<&str>,
        notes: Option<&str>,
    ) -> Result<bool, String> {
        let n = self
            .db()
            .execute(
                "UPDATE hulls SET name = COALESCE(?2, name), notes = COALESCE(?3, notes)
                 WHERE id = ?1",
                params![id, name, notes],
            )
            .map_err(err)?;
        Ok(n > 0)
    }

    // --- configurations ------------------------------------------------

    /// The configuration already made on `hull_id` with these parameters.
    pub fn find_config(&self, hull_id: i64, p: &ConfigParams) -> Result<Option<i64>, String> {
        self.db()
            .query_row(
                "SELECT id FROM configs WHERE hull_id = ?1 AND params_hash = ?2",
                params![hull_id, p.hash()],
                |r| r.get(0),
            )
            .optional()
            .map_err(err)
    }

    /// Save a configuration with its statics (or the error computing them);
    /// the same parameters on the same hull are the configuration already
    /// saved, returned as `(id, false)`.
    pub fn add_config(
        &self,
        hull_id: i64,
        p: &ConfigParams,
        name: &str,
        by: &str,
        statics: Result<&Value, &str>,
    ) -> Result<(i64, bool), String> {
        if let Some(id) = self.find_config(hull_id, p)? {
            return Ok((id, false));
        }
        let (summary, blob, error) = match statics {
            Ok(v) => (
                Some(statics_summary(v).to_string()),
                Some(self.put_json(v)?),
                None,
            ),
            Err(e) => (None, None, Some(e)),
        };
        let db = self.db();
        db.execute(
            "INSERT INTO configs (hull_id, name, params, params_hash, summary, statics_blob,
                                  error, solver_version, created_by, created_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)",
            params![
                hull_id,
                name,
                p.to_json(),
                p.hash(),
                summary,
                blob,
                error,
                SOLVER_VERSION,
                by,
                now()
            ],
        )
        .map_err(|e| match e {
            rusqlite::Error::SqliteFailure(f, _)
                if f.code == rusqlite::ErrorCode::ConstraintViolation =>
            {
                format!("no hull {hull_id}")
            }
            e => err(e),
        })?;
        Ok((db.last_insert_rowid(), true))
    }

    /// Replace a configuration's statics (computed again, by this solver).
    pub fn set_statics(&self, id: i64, statics: Result<&Value, &str>) -> Result<(), String> {
        let (summary, blob, error) = match statics {
            Ok(v) => (
                Some(statics_summary(v).to_string()),
                Some(self.put_json(v)?),
                None,
            ),
            Err(e) => (None, None, Some(e)),
        };
        self.db()
            .execute(
                "UPDATE configs SET summary = ?2, statics_blob = ?3, error = ?4,
                                    solver_version = ?5 WHERE id = ?1",
                params![id, summary, blob, error, SOLVER_VERSION],
            )
            .map(|_| ())
            .map_err(err)
    }

    /// Configurations (of one hull, or all), oldest first, each with its
    /// statics' summary and its runs' counts by status.
    pub fn configs(&self, hull: Option<i64>) -> Result<Vec<Value>, String> {
        let db = self.db();
        let mut q = db
            .prepare(&format!(
                "SELECT {CONFIG_COLS} FROM configs g WHERE (?1 IS NULL OR g.hull_id = ?1)
                 ORDER BY g.id"
            ))
            .map_err(err)?;
        let rows = q.query_map([hull], config_row).map_err(err)?;
        rows.collect::<Result<_, _>>().map_err(err)
    }

    pub fn config(&self, id: i64) -> Result<Option<Value>, String> {
        self.db()
            .query_row(
                &format!("SELECT {CONFIG_COLS} FROM configs g WHERE g.id = ?1"),
                [id],
                config_row,
            )
            .optional()
            .map_err(err)
    }

    /// A configuration's hull, parameters and statics blob.
    pub fn config_source(
        &self,
        id: i64,
    ) -> Result<Option<(i64, ConfigParams, Option<String>)>, String> {
        let row: Option<(i64, String, Option<String>)> = self
            .db()
            .query_row(
                "SELECT hull_id, params, statics_blob FROM configs WHERE id = ?1",
                [id],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
            )
            .optional()
            .map_err(err)?;
        row.map(|(h, p, b)| Ok((h, serde_json::from_str(&p).map_err(err)?, b)))
            .transpose()
    }

    pub fn update_config(
        &self,
        id: i64,
        name: Option<&str>,
        notes: Option<&str>,
    ) -> Result<bool, String> {
        let n = self
            .db()
            .execute(
                "UPDATE configs SET name = COALESCE(?2, name), notes = COALESCE(?3, notes)
                 WHERE id = ?1",
                params![id, name, notes],
            )
            .map_err(err)?;
        Ok(n > 0)
    }

    // --- studies and runs ----------------------------------------------

    pub fn add_study(&self, name: &str, notes: &str, by: &str) -> Result<i64, String> {
        let db = self.db();
        db.execute(
            "INSERT INTO studies (name, notes, created_by, created_at) VALUES (?1, ?2, ?3, ?4)",
            params![name, notes, by, now()],
        )
        .map_err(err)?;
        Ok(db.last_insert_rowid())
    }

    pub fn studies(&self) -> Result<Vec<Value>, String> {
        let db = self.db();
        let mut q = db
            .prepare(
                "SELECT s.id, s.name, s.notes, s.created_by, s.created_at,
                        (SELECT COUNT(*) FROM study_runs sr WHERE sr.study_id = s.id),
                        (SELECT COUNT(*) FROM study_runs sr JOIN runs r ON r.id = sr.run_id
                          WHERE sr.study_id = s.id AND r.status = 'done')
                 FROM studies s ORDER BY s.id DESC",
            )
            .map_err(err)?;
        let rows = q
            .query_map([], |r| {
                Ok(json!({
                    "id": r.get::<_, i64>(0)?,
                    "name": r.get::<_, String>(1)?,
                    "notes": r.get::<_, String>(2)?,
                    "created_by": r.get::<_, String>(3)?,
                    "created_at": r.get::<_, i64>(4)?,
                    "runs": r.get::<_, i64>(5)?,
                    "done": r.get::<_, i64>(6)?,
                }))
            })
            .map_err(err)?;
        rows.collect::<Result<_, _>>().map_err(err)
    }

    /// The run already asked for on `config_id` with these parameters, if
    /// any: its id, status, and whether its result is stale.
    pub fn find_run(
        &self,
        config_id: i64,
        p: &RunParams,
    ) -> Result<Option<(i64, String, bool)>, String> {
        self.db()
            .query_row(
                "SELECT r.id, r.status, EXISTS (SELECT 1 FROM results x
                         WHERE x.run_id = r.id AND x.solver_version != ?3)
                 FROM runs r WHERE r.config_id = ?1 AND r.params_hash = ?2",
                params![config_id, p.hash(), SOLVER_VERSION],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
            )
            .optional()
            .map_err(err)
    }

    /// Ask for a run on a configuration, returning `(id, created)`. A run
    /// already asked for is not duplicated: it joins the study, and if it
    /// failed or was cancelled it goes back on the queue.
    pub fn add_run(
        &self,
        config_id: i64,
        p: &RunParams,
        study_id: Option<i64>,
        priority: i64,
        by: &str,
    ) -> Result<(i64, bool), String> {
        let mut db = self.db();
        let tx = db.transaction().map_err(err)?;
        let hash = p.hash();
        let existing: Option<(i64, String)> = tx
            .query_row(
                "SELECT id, status FROM runs WHERE config_id = ?1 AND params_hash = ?2",
                params![config_id, hash],
                |r| Ok((r.get(0)?, r.get(1)?)),
            )
            .optional()
            .map_err(err)?;
        let (id, created) = match existing {
            Some((id, st)) => {
                if st == status::FAILED || st == status::CANCELLED {
                    tx.execute(
                        "UPDATE runs SET status = 'queued', error = NULL, queued_at = ?2,
                                         priority = MAX(priority, ?3)
                         WHERE id = ?1",
                        params![id, now(), priority],
                    )
                    .map_err(err)?;
                } else if st == status::QUEUED {
                    tx.execute(
                        "UPDATE runs SET priority = MAX(priority, ?2) WHERE id = ?1",
                        params![id, priority],
                    )
                    .map_err(err)?;
                }
                (id, false)
            }
            None => {
                tx.execute(
                    "INSERT INTO runs (config_id, kind, params, params_hash, attitude_key, status,
                                       priority, requested_by, created_at, queued_at)
                     VALUES (?1, ?2, ?3, ?4, ?5, 'queued', ?6, ?7, ?8, ?8)",
                    params![
                        config_id,
                        p.kind(),
                        p.to_json(),
                        hash,
                        p.attitude_key(),
                        priority,
                        by,
                        now()
                    ],
                )
                .map_err(|e| match e {
                    rusqlite::Error::SqliteFailure(f, _)
                        if f.code == rusqlite::ErrorCode::ConstraintViolation =>
                    {
                        format!("no configuration {config_id}")
                    }
                    e => err(e),
                })?;
                (tx.last_insert_rowid(), true)
            }
        };
        if let Some(s) = study_id {
            tx.execute(
                "INSERT OR IGNORE INTO study_runs (study_id, run_id) VALUES (?1, ?2)",
                params![s, id],
            )
            .map_err(err)?;
        }
        tx.commit().map_err(err)?;
        Ok((id, created))
    }

    /// Runs matching the filters (each optional), newest first — or most
    /// recently finished first — each with its configuration, hull and
    /// result's scalars.
    pub fn runs(&self, f: &RunFilter) -> Result<Vec<Value>, String> {
        let db = self.db();
        let mut q = db
            .prepare(&format!(
                "SELECT {RUN_COLS} FROM runs r JOIN configs g ON g.id = r.config_id
                 JOIN hulls h ON h.id = g.hull_id LEFT JOIN results x ON x.run_id = r.id
                 WHERE (?1 IS NULL OR g.hull_id = ?1)
                   AND (?2 IS NULL OR r.config_id = ?2)
                   AND (?3 IS NULL OR r.id IN (SELECT run_id FROM study_runs WHERE study_id = ?3))
                   AND (?4 IS NULL OR instr(',' || ?4 || ',', ',' || r.status || ',') > 0)
                   AND (?5 IS NULL OR r.kind = ?5)
                 ORDER BY CASE WHEN ?7 THEN r.finished_at END DESC, r.id DESC LIMIT ?6"
            ))
            .map_err(err)?;
        let rows = q
            .query_map(
                params![
                    f.hull,
                    f.config,
                    f.study,
                    f.status,
                    f.kind,
                    f.limit.unwrap_or(10_000),
                    f.by_finish
                ],
                run_row,
            )
            .map_err(err)?;
        rows.collect::<Result<_, _>>().map_err(err)
    }

    pub fn run(&self, id: i64) -> Result<Option<Value>, String> {
        let db = self.db();
        let r = db
            .query_row(
                &format!(
                    "SELECT {RUN_COLS} FROM runs r JOIN configs g ON g.id = r.config_id
                     JOIN hulls h ON h.id = g.hull_id LEFT JOIN results x ON x.run_id = r.id
                     WHERE r.id = ?1"
                ),
                [id],
                run_row,
            )
            .optional()
            .map_err(err)?;
        let Some(mut r) = r else { return Ok(None) };
        let mut q = db
            .prepare("SELECT study_id FROM study_runs WHERE run_id = ?1")
            .map_err(err)?;
        r["studies"] = json!(q
            .query_map([id], |r| r.get::<_, i64>(0))
            .map_err(err)?
            .collect::<Result<Vec<_>, _>>()
            .map_err(err)?);
        Ok(Some(r))
    }

    /// The field blob of a run's result.
    pub fn result_blob(&self, run_id: i64) -> Result<Option<String>, String> {
        self.db()
            .query_row(
                "SELECT field_blob FROM results WHERE run_id = ?1",
                [run_id],
                |r| r.get(0),
            )
            .optional()
            .map_err(err)
    }

    /// The queue: running then queued runs, in the order they will be
    /// taken (a run in waves waits, besides, for its calm-water run).
    pub fn queue(&self) -> Result<Vec<Value>, String> {
        let db = self.db();
        let mut q = db
            .prepare(&format!(
                "SELECT {RUN_COLS} FROM runs r JOIN configs g ON g.id = r.config_id
                 JOIN hulls h ON h.id = g.hull_id LEFT JOIN results x ON x.run_id = r.id
                 WHERE r.status IN ('running', 'queued')
                 ORDER BY r.status = 'running' DESC, r.priority DESC, r.id"
            ))
            .map_err(err)?;
        let rows = q.query_map([], run_row).map_err(err)?;
        rows.collect::<Result<_, _>>().map_err(err)
    }

    /// Counts of runs by status, and of stale results and statics.
    pub fn counts(&self) -> Result<Value, String> {
        let db = self.db();
        let mut q = db
            .prepare("SELECT status, COUNT(*) FROM runs GROUP BY status")
            .map_err(err)?;
        let mut out = json!({});
        for r in q
            .query_map([], |r| Ok((r.get::<_, String>(0)?, r.get::<_, i64>(1)?)))
            .map_err(err)?
        {
            let (s, n) = r.map_err(err)?;
            out[s] = json!(n);
        }
        let stale: i64 = db
            .query_row(
                "SELECT COUNT(*) FROM results x JOIN runs r ON r.id = x.run_id
                 WHERE x.solver_version != ?1 AND r.status = 'done'",
                [SOLVER_VERSION],
                |r| r.get(0),
            )
            .map_err(err)?;
        out["stale"] = json!(stale);
        Ok(out)
    }

    /// The calm-water run a run in waves is taken about: the done one if
    /// there is one, else the one asked for (queued, running, failed …).
    fn calm_for(
        db: &Connection,
        config_id: i64,
        attitude_key: &str,
    ) -> rusqlite::Result<Option<(i64, String, Option<String>)>> {
        db.query_row(
            "SELECT r.id, r.status, x.scalars FROM runs r LEFT JOIN results x ON x.run_id = r.id
             WHERE r.config_id = ?1 AND r.attitude_key = ?2 AND r.kind = 'calm'
             ORDER BY r.status = 'done' DESC, r.status IN ('queued', 'running') DESC, r.id
             LIMIT 1",
            params![config_id, attitude_key],
            |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
        )
        .optional()
    }

    /// Take the next run that can go (highest priority, then oldest; a run
    /// in waves once its calm-water run is done) and mark it running. A run
    /// in waves whose calm-water run failed, or was cancelled, fails too.
    pub fn claim(&self) -> Result<Option<Claimed>, String> {
        let mut db = self.db();
        let tx = db.transaction().map_err(err)?;
        let candidates: Vec<Candidate> = {
            let mut q = tx
                .prepare(
                    "SELECT r.id, r.config_id, r.kind, r.params, r.attitude_key, g.hull_id,
                            g.params, g.statics_blob
                     FROM runs r JOIN configs g ON g.id = r.config_id
                     WHERE r.status = 'queued' ORDER BY r.priority DESC, r.id LIMIT 200",
                )
                .map_err(err)?;
            let rows = q
                .query_map([], |r| {
                    Ok((
                        r.get(0)?,
                        r.get(1)?,
                        r.get(2)?,
                        r.get(3)?,
                        r.get(4)?,
                        r.get(5)?,
                        r.get(6)?,
                        r.get(7)?,
                    ))
                })
                .map_err(err)?;
            rows.collect::<Result<_, _>>().map_err(err)?
        };
        let mut chosen = None;
        for (id, config_id, kind, params_json, key, hull_id, config_json, statics) in candidates {
            let calm = if kind == "waves" {
                match Self::calm_for(&tx, config_id, &key).map_err(err)? {
                    Some((cid, st, Some(scalars))) if st == status::DONE => {
                        let s: Value = serde_json::from_str(&scalars).unwrap_or_default();
                        match (s["sinkage"].as_f64(), s["trim_rad"].as_f64()) {
                            (Some(z), Some(t)) => Some((cid, (z, t))),
                            _ => {
                                fail(&tx, id, &format!("calm-water run #{cid} has no attitude"))?;
                                continue;
                            }
                        }
                    }
                    Some((_, st, _)) if st == status::QUEUED || st == status::RUNNING => continue,
                    Some((cid, st, _)) => {
                        fail(&tx, id, &format!("its calm-water run #{cid} is {st}"))?;
                        continue;
                    }
                    None => {
                        fail(&tx, id, "it has no calm-water run to be taken about")?;
                        continue;
                    }
                }
            } else {
                None
            };
            chosen = Some((
                id,
                config_id,
                params_json,
                hull_id,
                config_json,
                statics,
                calm,
            ));
            break;
        }
        let Some((id, config_id, params_json, hull_id, config_json, statics, calm)) = chosen else {
            tx.commit().map_err(err)?;
            return Ok(None);
        };
        tx.execute(
            "UPDATE runs SET status = 'running', started_at = ?2, error = NULL WHERE id = ?1",
            params![id, now()],
        )
        .map_err(err)?;
        let hull: (String, String, String) = tx
            .query_row(
                "SELECT file_name, file_blob, import FROM hulls WHERE id = ?1",
                [hull_id],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
            )
            .map_err(err)?;
        tx.commit().map_err(err)?;
        drop(db);
        // The attitude at rest, from the statics (outside the lock: a blob).
        let at_rest = match statics {
            Some(b) => {
                let v: Value = serde_json::from_slice(&self.blob(&b)?).map_err(err)?;
                match (
                    v["at_rest"]["sinkage"].as_f64(),
                    v["at_rest"]["trim_rad"].as_f64(),
                ) {
                    (Some(z), Some(t)) => Some((z, t)),
                    _ => None,
                }
            }
            None => None,
        };
        Ok(Some(Claimed {
            id,
            config_id,
            params: serde_json::from_str(&params_json).map_err(err)?,
            config: serde_json::from_str(&config_json).map_err(err)?,
            hull: HullSource {
                file_name: hull.0,
                file_blob: hull.1,
                import: serde_json::from_str(&hull.2).map_err(err)?,
            },
            at_rest,
            calm,
        }))
    }

    /// A running run finished: its result replaces any it had.
    pub fn finish(&self, id: i64, r: &NewResult) -> Result<(), String> {
        let blob = self.put_json(&r.field)?;
        let mut db = self.db();
        let tx = db.transaction().map_err(err)?;
        let t = now();
        tx.execute(
            "INSERT OR REPLACE INTO results (run_id, scalars, field_blob, solver_version, seconds,
                                             warm_from, attitude_from, created_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            params![
                id,
                r.scalars.to_string(),
                blob,
                SOLVER_VERSION,
                r.seconds,
                r.warm_from,
                r.attitude_from,
                t
            ],
        )
        .map_err(err)?;
        tx.execute(
            "UPDATE runs SET status = 'done', finished_at = ?2, error = NULL WHERE id = ?1",
            params![id, t],
        )
        .map_err(err)?;
        tx.commit().map_err(err)
    }

    /// A running run stopped without a result (`failed` or `cancelled`).
    /// A result from an earlier run is kept.
    pub fn stop(&self, id: i64, st: &str, error: Option<&str>) -> Result<(), String> {
        self.db()
            .execute(
                "UPDATE runs SET status = ?2, finished_at = ?3, error = ?4 WHERE id = ?1",
                params![id, st, now(), error],
            )
            .map(|_| ())
            .map_err(err)
    }

    /// Cancel a queued run, returning its status before (the worker stops a
    /// running one itself).
    pub fn cancel_queued(&self, id: i64) -> Result<Option<String>, String> {
        let db = self.db();
        let st: Option<String> = db
            .query_row("SELECT status FROM runs WHERE id = ?1", [id], |r| r.get(0))
            .optional()
            .map_err(err)?;
        if st.as_deref() == Some(status::QUEUED) {
            db.execute(
                "UPDATE runs SET status = 'cancelled', finished_at = ?2 WHERE id = ?1",
                params![id, now()],
            )
            .map_err(err)?;
        }
        Ok(st)
    }

    /// Put finished runs back on the queue — the given ones (with, for a run
    /// in waves, its calm-water run if that failed or was cancelled), or
    /// every one with a stale result (optionally on one hull). Returns how
    /// many.
    pub fn requeue(
        &self,
        ids: Option<&[i64]>,
        stale_on_hull: Option<i64>,
    ) -> Result<usize, String> {
        let db = self.db();
        let t = now();
        match ids {
            Some(ids) => {
                let mut n = 0;
                for id in ids {
                    n += db
                        .execute(
                            "UPDATE runs SET status = 'queued', error = NULL, queued_at = ?2
                             WHERE id = ?1 AND status IN ('done', 'failed', 'cancelled')",
                            params![id, t],
                        )
                        .map_err(err)?;
                    n += db
                        .execute(
                            "UPDATE runs SET status = 'queued', error = NULL, queued_at = ?2
                             WHERE status IN ('failed', 'cancelled') AND kind = 'calm'
                               AND (config_id, attitude_key) IN
                                   (SELECT config_id, attitude_key FROM runs
                                     WHERE id = ?1 AND kind = 'waves')",
                            params![id, t],
                        )
                        .map_err(err)?;
                }
                Ok(n)
            }
            None => db
                .execute(
                    "UPDATE runs SET status = 'queued', error = NULL, queued_at = ?1
                     WHERE status = 'done'
                       AND (?2 IS NULL OR config_id IN (SELECT id FROM configs WHERE hull_id = ?2))
                       AND id IN (SELECT run_id FROM results WHERE solver_version != ?3)",
                    params![t, stale_on_hull, SOLVER_VERSION],
                )
                .map_err(err),
        }
    }

    pub fn set_priority(&self, id: i64, priority: i64) -> Result<bool, String> {
        let n = self
            .db()
            .execute(
                "UPDATE runs SET priority = ?2 WHERE id = ?1",
                params![id, priority],
            )
            .map_err(err)?;
        Ok(n > 0)
    }

    /// For a warm start: the done calm-water runs on `config_id` that differ
    /// from `p` only in speed (and grid), as `(run id, froude, sinkage,
    /// trim)` from their solved equilibria.
    pub fn solved_neighbours(
        &self,
        config_id: i64,
        p: &RunParams,
    ) -> Result<Vec<(i64, f64, f64, f64)>, String> {
        let db = self.db();
        let mut q = db
            .prepare(
                "SELECT r.id, r.params, x.scalars FROM runs r JOIN results x ON x.run_id = r.id
                 WHERE r.config_id = ?1 AND r.kind = 'calm' AND r.status = 'done'",
            )
            .map_err(err)?;
        let rows = q
            .query_map([config_id], |r| {
                Ok((
                    r.get::<_, i64>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, String>(2)?,
                ))
            })
            .map_err(err)?;
        let mut out = Vec::new();
        for row in rows {
            let (id, q, s) = row.map_err(err)?;
            let (Ok(q), Ok(s)) = (
                serde_json::from_str::<RunParams>(&q),
                serde_json::from_str::<Value>(&s),
            ) else {
                continue;
            };
            if !p.same_but_speed(&q) || s["solved"].as_bool() != Some(true) {
                continue;
            }
            if let (Some(z), Some(t)) = (s["sinkage"].as_f64(), s["trim_rad"].as_f64()) {
                out.push((id, q.froude, z, t));
            }
        }
        Ok(out)
    }
}

/// A queued run the worker might take: its id, configuration, kind, params
/// and attitude key, and its configuration's hull, params and statics blob.
type Candidate = (
    i64,
    i64,
    String,
    String,
    String,
    i64,
    String,
    Option<String>,
);

fn fail(tx: &rusqlite::Transaction, id: i64, why: &str) -> Result<(), String> {
    tx.execute(
        "UPDATE runs SET status = 'failed', error = ?2, finished_at = ?3 WHERE id = ?1",
        params![id, why, now()],
    )
    .map(|_| ())
    .map_err(err)
}

/// Filters for [`Store::runs`].
#[derive(Default)]
pub struct RunFilter {
    pub hull: Option<i64>,
    pub config: Option<i64>,
    pub study: Option<i64>,
    /// One status, or several separated by commas.
    pub status: Option<String>,
    /// `calm` or `waves`.
    pub kind: Option<String>,
    pub limit: Option<i64>,
    /// Most recently finished first, rather than newest.
    pub by_finish: bool,
}

/// What the hull list shows of a cut: per hull, the principal dimensions and
/// hydrostatics.
fn hull_summary(sections: &Value) -> Value {
    let hulls: Vec<Value> = sections["hulls"]
        .as_array()
        .map(|a| {
            a.iter()
                .map(|h| {
                    let mut s = json!({});
                    for k in [
                        "length",
                        "beam",
                        "draft",
                        "displaced_volume",
                        "wetted_surface",
                        "lcb_x",
                        "waterplane_area",
                    ] {
                        s[k] = h[k].clone();
                    }
                    s["transom"] = json!(!h["transom"].is_null());
                    s
                })
                .collect()
        })
        .unwrap_or_default();
    json!({ "hulls": hulls, "notes": sections["notes"] })
}

/// What lists show of a configuration's statics: everything but the meshes
/// and the GZ curve's points.
fn statics_summary(v: &Value) -> Value {
    let mut s = v.clone();
    if let Some(o) = s.as_object_mut() {
        o.remove("meshes");
    }
    if let Some(g) = s["gz"].as_object_mut() {
        for k in ["heel_deg", "gz", "sinkage", "trim_deg"] {
            g.remove(k);
        }
    }
    s
}

const HULL_COLS: &str = "h.id, h.name, h.notes, h.uploaded_by, h.file_name, h.import, h.parent_id,
    h.summary, h.created_at,
    (SELECT COUNT(*) FROM configs WHERE hull_id = h.id),
    (SELECT json_group_object(status, n) FROM
        (SELECT r.status, COUNT(*) AS n FROM runs r JOIN configs g ON g.id = r.config_id
          WHERE g.hull_id = h.id GROUP BY r.status))";

fn hull_row(r: &Row) -> rusqlite::Result<Value> {
    Ok(json!({
        "id": r.get::<_, i64>(0)?,
        "name": r.get::<_, String>(1)?,
        "notes": r.get::<_, String>(2)?,
        "uploaded_by": r.get::<_, String>(3)?,
        "file_name": r.get::<_, String>(4)?,
        "import": parse(r.get(5)?),
        "parent_id": r.get::<_, Option<i64>>(6)?,
        "summary": parse(r.get(7)?),
        "created_at": r.get::<_, i64>(8)?,
        "configs": r.get::<_, i64>(9)?,
        "runs": parse(r.get(10)?),
    }))
}

const CONFIG_COLS: &str = "g.id, g.hull_id, g.name, g.notes, g.params, g.summary, g.error,
    g.solver_version, g.created_by, g.created_at,
    (SELECT name FROM hulls WHERE id = g.hull_id),
    (SELECT json_group_object(status, n) FROM
        (SELECT status, COUNT(*) AS n FROM runs WHERE config_id = g.id GROUP BY status))";

fn config_row(r: &Row) -> rusqlite::Result<Value> {
    let version: String = r.get(7)?;
    Ok(json!({
        "id": r.get::<_, i64>(0)?,
        "hull_id": r.get::<_, i64>(1)?,
        "name": r.get::<_, String>(2)?,
        "notes": r.get::<_, String>(3)?,
        "params": parse(r.get(4)?),
        "statics": parse(r.get(5)?),
        "error": r.get::<_, Option<String>>(6)?,
        "stale": version != SOLVER_VERSION,
        "solver_version": version,
        "created_by": r.get::<_, String>(8)?,
        "created_at": r.get::<_, i64>(9)?,
        "hull_name": r.get::<_, Option<String>>(10)?,
        "runs": parse(r.get(11)?),
    }))
}

const RUN_COLS: &str = "r.id, r.config_id, r.kind, r.params, r.status, r.priority, r.requested_by,
    r.error, r.created_at, r.queued_at, r.started_at, r.finished_at,
    g.hull_id, h.name, g.name, g.params,
    x.scalars, x.solver_version, x.seconds, x.warm_from, x.attitude_from, x.created_at";

fn run_row(r: &Row) -> rusqlite::Result<Value> {
    let version: Option<String> = r.get(17)?;
    let result = match version {
        None => Value::Null,
        Some(v) => json!({
            "scalars": parse(r.get(16)?),
            "stale": v != SOLVER_VERSION,
            "solver_version": v,
            "seconds": r.get::<_, f64>(18)?,
            "warm_from": r.get::<_, Option<i64>>(19)?,
            "attitude_from": r.get::<_, Option<i64>>(20)?,
            "created_at": r.get::<_, i64>(21)?,
        }),
    };
    Ok(json!({
        "id": r.get::<_, i64>(0)?,
        "config_id": r.get::<_, i64>(1)?,
        "kind": r.get::<_, String>(2)?,
        "params": parse(r.get(3)?),
        "status": r.get::<_, String>(4)?,
        "priority": r.get::<_, i64>(5)?,
        "requested_by": r.get::<_, String>(6)?,
        "error": r.get::<_, Option<String>>(7)?,
        "created_at": r.get::<_, i64>(8)?,
        "queued_at": r.get::<_, Option<i64>>(9)?,
        "started_at": r.get::<_, Option<i64>>(10)?,
        "finished_at": r.get::<_, Option<i64>>(11)?,
        "hull_id": r.get::<_, i64>(12)?,
        "hull_name": r.get::<_, String>(13)?,
        "config_name": r.get::<_, String>(14)?,
        "config": parse(r.get(15)?),
        "stale": result["stale"] == true,
        "result": result,
    }))
}

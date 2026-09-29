//! The permanent record: hulls, cases, their results, and studies (named
//! groups of cases asked for together), in SQLite, with the bulky data —
//! uploaded files, a hull's sections and display geometry, a result's
//! fields — in a content-addressed directory of gzipped blobs beside it.
//!
//! ```text
//! $DATA_DIR/michell.db
//! $DATA_DIR/blobs/ab/cdef….gz    SHA-256 of the uncompressed bytes
//! ```
//!
//! A hull is a file together with its import settings (the cut depends on
//! both); a case is a hull and [`CaseParams`], unique on the pair; a result
//! is one computed point of a case (a flow has one, a span sweep one per
//! span) and records the solver version it was computed with — a result
//! from another version is stale, kept until its case is run again.

use crate::case::{hex, CaseParams};
use crate::LoftRequest;
use rusqlite::{params, Connection, OptionalExtension, Row};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// The solver version results are stamped with: the last commit to touch
/// the solver's code (the geometry, thin-ship, seakeeping and CLI crates,
/// and the web crate's flow wrapper), `-dirty` if they have uncommitted changes. Set by
/// `build.rs`; `MICHELL_SOLVER_VERSION` at build time overrides it.
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
CREATE TABLE IF NOT EXISTS studies (
    id         INTEGER PRIMARY KEY,
    name       TEXT NOT NULL,
    notes      TEXT NOT NULL DEFAULT '',
    created_by TEXT NOT NULL DEFAULT '',
    created_at INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS cases (
    id           INTEGER PRIMARY KEY,
    hull_id      INTEGER NOT NULL REFERENCES hulls(id),
    kind         TEXT NOT NULL,
    params       TEXT NOT NULL,
    params_hash  TEXT NOT NULL,
    status       TEXT NOT NULL,
    priority     INTEGER NOT NULL DEFAULT 0,
    requested_by TEXT NOT NULL DEFAULT '',
    error        TEXT,
    created_at   INTEGER NOT NULL,
    queued_at    INTEGER,
    started_at   INTEGER,
    finished_at  INTEGER,
    UNIQUE (hull_id, params_hash)
);
CREATE INDEX IF NOT EXISTS cases_queue ON cases (status, priority, id);
CREATE TABLE IF NOT EXISTS study_cases (
    study_id INTEGER NOT NULL REFERENCES studies(id),
    case_id  INTEGER NOT NULL REFERENCES cases(id),
    PRIMARY KEY (study_id, case_id)
);
CREATE TABLE IF NOT EXISTS results (
    case_id        INTEGER NOT NULL REFERENCES cases(id),
    idx            INTEGER NOT NULL,
    point          TEXT NOT NULL,
    scalars        TEXT NOT NULL,
    field_blob     TEXT NOT NULL,
    solver_version TEXT NOT NULL,
    seconds        REAL NOT NULL,
    warm_from      INTEGER REFERENCES cases(id),
    created_at     INTEGER NOT NULL,
    PRIMARY KEY (case_id, idx)
);
";

/// A case's state in the queue.
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

/// A case the worker has claimed, with what it needs to run it.
pub struct Claimed {
    pub id: i64,
    pub hull_id: i64,
    pub params: CaseParams,
    pub file_name: String,
    pub file_blob: String,
    pub import: LoftRequest,
}

/// One computed point, ready to save.
pub struct NewResult {
    pub point: Value,
    pub scalars: Value,
    /// The full answer (fields and all), stored as a blob.
    pub field: Value,
    pub seconds: f64,
}

pub fn now() -> i64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs() as i64)
}

fn err(e: impl std::fmt::Display) -> String {
    e.to_string()
}

impl Store {
    /// Open (or create) the store in `dir`. Cases left running by a server
    /// that stopped go back on the queue.
    pub fn open(dir: &Path) -> Result<Store, String> {
        let blobs = dir.join("blobs");
        std::fs::create_dir_all(&blobs).map_err(|e| format!("{}: {e}", blobs.display()))?;
        let db = Connection::open(dir.join("michell.db")).map_err(err)?;
        db.execute_batch("PRAGMA journal_mode = WAL; PRAGMA foreign_keys = ON;")
            .map_err(err)?;
        db.execute_batch(SCHEMA).map_err(err)?;
        db.execute(
            "UPDATE cases SET status = ?1, started_at = NULL WHERE status = ?2",
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

    /// Every hull, newest first, with its case counts by status.
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

    // --- studies and cases ---------------------------------------------

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
                        (SELECT COUNT(*) FROM study_cases sc WHERE sc.study_id = s.id),
                        (SELECT COUNT(*) FROM study_cases sc JOIN cases c ON c.id = sc.case_id
                          WHERE sc.study_id = s.id AND c.status = 'done')
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
                    "cases": r.get::<_, i64>(5)?,
                    "done": r.get::<_, i64>(6)?,
                }))
            })
            .map_err(err)?;
        rows.collect::<Result<_, _>>().map_err(err)
    }

    /// Ask for a case on a hull, returning `(id, created)`. A case already
    /// asked for is not duplicated: it joins the study, and if it failed or
    /// was cancelled it goes back on the queue.
    pub fn add_case(
        &self,
        hull_id: i64,
        params: &CaseParams,
        study_id: Option<i64>,
        priority: i64,
        by: &str,
    ) -> Result<(i64, bool), String> {
        let mut db = self.db();
        let tx = db.transaction().map_err(err)?;
        let hash = params.hash();
        let existing: Option<(i64, String)> = tx
            .query_row(
                "SELECT id, status FROM cases WHERE hull_id = ?1 AND params_hash = ?2",
                params![hull_id, hash],
                |r| Ok((r.get(0)?, r.get(1)?)),
            )
            .optional()
            .map_err(err)?;
        let (id, created) = match existing {
            Some((id, st)) => {
                if st == status::FAILED || st == status::CANCELLED {
                    tx.execute(
                        "UPDATE cases SET status = ?2, error = NULL, queued_at = ?3,
                                          priority = MAX(priority, ?4)
                         WHERE id = ?1",
                        params![id, status::QUEUED, now(), priority],
                    )
                    .map_err(err)?;
                }
                (id, false)
            }
            None => {
                tx.execute(
                    "INSERT INTO cases (hull_id, kind, params, params_hash, status, priority,
                                        requested_by, created_at, queued_at)
                     VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?8)",
                    params![
                        hull_id,
                        params.kind(),
                        params.to_json(),
                        hash,
                        status::QUEUED,
                        priority,
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
                (tx.last_insert_rowid(), true)
            }
        };
        if let Some(s) = study_id {
            tx.execute(
                "INSERT OR IGNORE INTO study_cases (study_id, case_id) VALUES (?1, ?2)",
                params![s, id],
            )
            .map_err(err)?;
        }
        tx.commit().map_err(err)?;
        Ok((id, created))
    }

    /// The case already asked for on `hull_id` with these parameters, if
    /// any: its id, status, and whether its results are stale.
    pub fn find_case(
        &self,
        hull_id: i64,
        params: &CaseParams,
    ) -> Result<Option<(i64, String, bool)>, String> {
        self.db()
            .query_row(
                "SELECT c.id, c.status, EXISTS (SELECT 1 FROM results r
                         WHERE r.case_id = c.id AND r.solver_version != ?3)
                 FROM cases c WHERE c.hull_id = ?1 AND c.params_hash = ?2",
                params![hull_id, params.hash(), SOLVER_VERSION],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)),
            )
            .optional()
            .map_err(err)
    }

    /// Cases matching the filters (each optional), newest first — or most
    /// recently finished first — each with its results' scalars.
    pub fn cases(&self, f: &CaseFilter) -> Result<Vec<Value>, String> {
        let db = self.db();
        let mut q = db
            .prepare(&format!(
                "SELECT {CASE_COLS} FROM cases c
                 WHERE (?1 IS NULL OR c.hull_id = ?1)
                   AND (?2 IS NULL OR c.id IN (SELECT case_id FROM study_cases WHERE study_id = ?2))
                   AND (?3 IS NULL OR instr(',' || ?3 || ',', ',' || c.status || ',') > 0)
                 ORDER BY CASE WHEN ?5 THEN c.finished_at END DESC, c.id DESC LIMIT ?4"
            ))
            .map_err(err)?;
        let mut cases: Vec<Value> = q
            .query_map(
                params![
                    f.hull,
                    f.study,
                    f.status,
                    f.limit.unwrap_or(10_000),
                    f.by_finish
                ],
                case_row,
            )
            .map_err(err)?
            .collect::<Result<_, _>>()
            .map_err(err)?;
        attach_results(&db, &mut cases)?;
        Ok(cases)
    }

    pub fn case(&self, id: i64) -> Result<Option<Value>, String> {
        let db = self.db();
        let c = db
            .query_row(
                &format!("SELECT {CASE_COLS} FROM cases c WHERE c.id = ?1"),
                [id],
                case_row,
            )
            .optional()
            .map_err(err)?;
        let Some(c) = c else { return Ok(None) };
        let mut v = vec![c];
        attach_results(&db, &mut v)?;
        let mut c = v.pop().expect("one case");
        let mut q = db
            .prepare("SELECT study_id FROM study_cases WHERE case_id = ?1")
            .map_err(err)?;
        c["studies"] = json!(q
            .query_map([id], |r| r.get::<_, i64>(0))
            .map_err(err)?
            .collect::<Result<Vec<_>, _>>()
            .map_err(err)?);
        Ok(Some(c))
    }

    /// The field blob of one of a case's results.
    pub fn result_blob(&self, case_id: i64, idx: i64) -> Result<Option<String>, String> {
        self.db()
            .query_row(
                "SELECT field_blob FROM results WHERE case_id = ?1 AND idx = ?2",
                params![case_id, idx],
                |r| r.get(0),
            )
            .optional()
            .map_err(err)
    }

    /// The queue: running then queued cases, in the order they will run.
    pub fn queue(&self) -> Result<Vec<Value>, String> {
        let db = self.db();
        let mut q = db
            .prepare(&format!(
                "SELECT {CASE_COLS} FROM cases c WHERE c.status IN ('running', 'queued')
                 ORDER BY c.status = 'running' DESC, c.priority DESC, c.id"
            ))
            .map_err(err)?;
        let mut cases: Vec<Value> = q
            .query_map([], case_row)
            .map_err(err)?
            .collect::<Result<_, _>>()
            .map_err(err)?;
        // A requeued case keeps its earlier results until it runs again.
        attach_results(&db, &mut cases)?;
        Ok(cases)
    }

    /// Counts of cases by status.
    pub fn counts(&self) -> Result<Value, String> {
        let db = self.db();
        let mut q = db
            .prepare("SELECT status, COUNT(*) FROM cases GROUP BY status")
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
                "SELECT COUNT(DISTINCT case_id) FROM results WHERE solver_version != ?1",
                [SOLVER_VERSION],
                |r| r.get(0),
            )
            .map_err(err)?;
        out["stale"] = json!(stale);
        Ok(out)
    }

    /// Take the next queued case (highest priority, then oldest) and mark it
    /// running.
    pub fn claim(&self) -> Result<Option<Claimed>, String> {
        let mut db = self.db();
        let tx = db.transaction().map_err(err)?;
        let c = tx
            .query_row(
                "SELECT c.id, c.hull_id, c.params, h.file_name, h.file_blob, h.import
                 FROM cases c JOIN hulls h ON h.id = c.hull_id
                 WHERE c.status = 'queued' ORDER BY c.priority DESC, c.id LIMIT 1",
                [],
                |r| {
                    Ok((
                        r.get::<_, i64>(0)?,
                        r.get::<_, i64>(1)?,
                        r.get::<_, String>(2)?,
                        r.get::<_, String>(3)?,
                        r.get::<_, String>(4)?,
                        r.get::<_, String>(5)?,
                    ))
                },
            )
            .optional()
            .map_err(err)?;
        let Some((id, hull_id, params, file_name, file_blob, import)) = c else {
            return Ok(None);
        };
        tx.execute(
            "UPDATE cases SET status = ?2, started_at = ?3, error = NULL WHERE id = ?1",
            params![id, status::RUNNING, now()],
        )
        .map_err(err)?;
        tx.commit().map_err(err)?;
        Ok(Some(Claimed {
            id,
            hull_id,
            params: serde_json::from_str(&params).map_err(err)?,
            file_name,
            file_blob,
            import: serde_json::from_str(&import).map_err(err)?,
        }))
    }

    /// A running case finished: its results replace any it had.
    pub fn finish(
        &self,
        id: i64,
        results: &[NewResult],
        warm_from: Option<i64>,
    ) -> Result<(), String> {
        let blobs: Vec<String> = results
            .iter()
            .map(|r| self.put_json(&r.field))
            .collect::<Result<_, _>>()?;
        let mut db = self.db();
        let tx = db.transaction().map_err(err)?;
        tx.execute("DELETE FROM results WHERE case_id = ?1", [id])
            .map_err(err)?;
        let t = now();
        for (i, (r, blob)) in results.iter().zip(&blobs).enumerate() {
            tx.execute(
                "INSERT INTO results (case_id, idx, point, scalars, field_blob, solver_version,
                                      seconds, warm_from, created_at)
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
                params![
                    id,
                    i as i64,
                    r.point.to_string(),
                    r.scalars.to_string(),
                    blob,
                    SOLVER_VERSION,
                    r.seconds,
                    warm_from,
                    t
                ],
            )
            .map_err(err)?;
        }
        tx.execute(
            "UPDATE cases SET status = ?2, finished_at = ?3, error = NULL WHERE id = ?1",
            params![id, status::DONE, t],
        )
        .map_err(err)?;
        tx.commit().map_err(err)
    }

    /// A running case stopped without results (`failed` or `cancelled`).
    /// Results from an earlier run are kept.
    pub fn stop(&self, id: i64, st: &str, error: Option<&str>) -> Result<(), String> {
        self.db()
            .execute(
                "UPDATE cases SET status = ?2, finished_at = ?3, error = ?4 WHERE id = ?1",
                params![id, st, now(), error],
            )
            .map(|_| ())
            .map_err(err)
    }

    /// Cancel a queued case, returning its status before (the worker stops
    /// a running one itself).
    pub fn cancel_queued(&self, id: i64) -> Result<Option<String>, String> {
        let db = self.db();
        let st: Option<String> = db
            .query_row("SELECT status FROM cases WHERE id = ?1", [id], |r| r.get(0))
            .optional()
            .map_err(err)?;
        if st.as_deref() == Some(status::QUEUED) {
            db.execute(
                "UPDATE cases SET status = ?2, finished_at = ?3 WHERE id = ?1",
                params![id, status::CANCELLED, now()],
            )
            .map_err(err)?;
        }
        Ok(st)
    }

    /// Put finished cases back on the queue — the given ones, or every one
    /// with a stale result (optionally on one hull). Returns how many.
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
                            "UPDATE cases SET status = 'queued', error = NULL, queued_at = ?2
                             WHERE id = ?1 AND status IN ('done', 'failed', 'cancelled')",
                            params![id, t],
                        )
                        .map_err(err)?;
                }
                Ok(n)
            }
            None => db
                .execute(
                    "UPDATE cases SET status = 'queued', error = NULL, queued_at = ?1
                     WHERE status = 'done' AND (?2 IS NULL OR hull_id = ?2)
                       AND id IN (SELECT case_id FROM results WHERE solver_version != ?3)",
                    params![t, stale_on_hull, SOLVER_VERSION],
                )
                .map_err(err),
        }
    }

    pub fn set_priority(&self, id: i64, priority: i64) -> Result<bool, String> {
        let n = self
            .db()
            .execute(
                "UPDATE cases SET priority = ?2 WHERE id = ?1",
                params![id, priority],
            )
            .map_err(err)?;
        Ok(n > 0)
    }

    /// For a warm start: the done flow cases on `hull_id` that differ from
    /// `params` only in speed (and grid), as `(case id, froude, sinkage,
    /// trim)` from their solved equilibria.
    pub fn solved_neighbours(
        &self,
        hull_id: i64,
        params: &CaseParams,
    ) -> Result<Vec<(i64, f64, f64, f64)>, String> {
        let db = self.db();
        let mut q = db
            .prepare(
                "SELECT c.id, c.params, r.scalars FROM cases c
                 JOIN results r ON r.case_id = c.id AND r.idx = 0
                 WHERE c.hull_id = ?1 AND c.kind = 'flow' AND c.status = 'done'",
            )
            .map_err(err)?;
        let rows = q
            .query_map([hull_id], |r| {
                Ok((
                    r.get::<_, i64>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, String>(2)?,
                ))
            })
            .map_err(err)?;
        let mut out = Vec::new();
        for row in rows {
            let (id, p, s) = row.map_err(err)?;
            let (Ok(p), Ok(s)) = (
                serde_json::from_str::<CaseParams>(&p),
                serde_json::from_str::<Value>(&s),
            ) else {
                continue;
            };
            if !params.same_but_speed(&p) || s["solved"].as_bool() != Some(true) {
                continue;
            }
            if let (Some(z), Some(t)) = (s["sinkage"].as_f64(), s["trim_rad"].as_f64()) {
                out.push((id, p.froude, z, t));
            }
        }
        Ok(out)
    }
}

/// Filters for [`Store::cases`].
#[derive(Default)]
pub struct CaseFilter {
    pub hull: Option<i64>,
    pub study: Option<i64>,
    /// One status, or several separated by commas.
    pub status: Option<String>,
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

const HULL_COLS: &str = "h.id, h.name, h.notes, h.uploaded_by, h.file_name, h.import, h.parent_id,
    h.summary, h.created_at,
    (SELECT json_group_object(status, n) FROM
        (SELECT status, COUNT(*) AS n FROM cases WHERE hull_id = h.id GROUP BY status))";

fn hull_row(r: &Row) -> rusqlite::Result<Value> {
    let parse = |s: String| serde_json::from_str::<Value>(&s).unwrap_or(Value::Null);
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
        "cases": parse(r.get(9)?),
    }))
}

const CASE_COLS: &str = "c.id, c.hull_id, c.kind, c.params, c.status, c.priority, c.requested_by,
    c.error, c.created_at, c.queued_at, c.started_at, c.finished_at,
    (SELECT name FROM hulls WHERE id = c.hull_id)";

fn case_row(r: &Row) -> rusqlite::Result<Value> {
    Ok(json!({
        "id": r.get::<_, i64>(0)?,
        "hull_id": r.get::<_, i64>(1)?,
        "kind": r.get::<_, String>(2)?,
        "params": serde_json::from_str::<Value>(&r.get::<_, String>(3)?).unwrap_or(Value::Null),
        "status": r.get::<_, String>(4)?,
        "priority": r.get::<_, i64>(5)?,
        "requested_by": r.get::<_, String>(6)?,
        "error": r.get::<_, Option<String>>(7)?,
        "created_at": r.get::<_, i64>(8)?,
        "queued_at": r.get::<_, Option<i64>>(9)?,
        "started_at": r.get::<_, Option<i64>>(10)?,
        "finished_at": r.get::<_, Option<i64>>(11)?,
        "hull_name": r.get::<_, Option<String>>(12)?,
    }))
}

/// Each case's results (scalars, not fields), as `results`, and whether any
/// is stale.
fn attach_results(db: &Connection, cases: &mut [Value]) -> Result<(), String> {
    let mut q = db
        .prepare(
            "SELECT idx, point, scalars, solver_version, seconds, warm_from, created_at
             FROM results WHERE case_id = ?1 ORDER BY idx",
        )
        .map_err(err)?;
    for c in cases.iter_mut() {
        let id = c["id"].as_i64().expect("a case id");
        let results: Vec<Value> = q
            .query_map([id], |r| {
                let version: String = r.get(3)?;
                Ok(json!({
                    "idx": r.get::<_, i64>(0)?,
                    "point": serde_json::from_str::<Value>(&r.get::<_, String>(1)?)
                        .unwrap_or(Value::Null),
                    "scalars": serde_json::from_str::<Value>(&r.get::<_, String>(2)?)
                        .unwrap_or(Value::Null),
                    "stale": version != SOLVER_VERSION,
                    "solver_version": version,
                    "seconds": r.get::<_, f64>(4)?,
                    "warm_from": r.get::<_, Option<i64>>(5)?,
                    "created_at": r.get::<_, i64>(6)?,
                }))
            })
            .map_err(err)?
            .collect::<Result<_, _>>()
            .map_err(err)?;
        c["stale"] = json!(results.iter().any(|r| r["stale"] == true));
        c["results"] = json!(results);
    }
    Ok(())
}

//! The store: every record a command writes, by type and id, and the bulky
//! data records refer to, by content.
//!
//! ```text
//! $BOATMATH_HOME/records/hull/ab/cdef….json
//! $BOATMATH_HOME/records/case/…  study/…  result/…  sections/…
//! $BOATMATH_HOME/blobs/ab/cdef….gz        results' fields, by SHA-256 of the uncompressed bytes
//! ```
//!
//! `$BOATMATH_HOME` defaults to `~/.boatmath`. A result's id is its study's,
//! so the two live under their own types. Writes go to a temporary file and
//! are renamed into place, so concurrent runs never see half a record.

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

/// The record types, in ancestry order.
pub const TYPES: [&str; 5] = ["hull", "case", "study", "result", "sections"];

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

pub fn sha256(bytes: &[u8]) -> String {
    hex(&Sha256::digest(bytes))
}

/// The id of a record's inputs: the SHA-256 of their JSON.
pub fn id_of(inputs: &Value) -> String {
    sha256(inputs.to_string().as_bytes())
}

pub struct Store {
    root: PathBuf,
}

fn split(dir: &Path, id: &str, ext: &str) -> PathBuf {
    let (a, b) = id.split_at(2.min(id.len()));
    dir.join(a).join(format!("{b}.{ext}"))
}

fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let dir = path.parent().expect("store paths have a parent");
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let tmp = path.with_extension(format!("tmp{}", std::process::id()));
    std::fs::write(&tmp, bytes).map_err(|e| format!("{}: {e}", tmp.display()))?;
    std::fs::rename(&tmp, path).map_err(|e| format!("{}: {e}", path.display()))
}

impl Store {
    pub fn open() -> Result<Store, String> {
        let root = match std::env::var_os("BOATMATH_HOME") {
            Some(d) => PathBuf::from(d),
            None => {
                let home =
                    std::env::var_os("HOME").ok_or("neither BOATMATH_HOME nor HOME is set")?;
                PathBuf::from(home).join(".boatmath")
            }
        };
        Ok(Store { root })
    }

    fn record_path(&self, kind: &str, id: &str) -> PathBuf {
        split(&self.root.join("records").join(kind), id, "json")
    }

    /// Save a record under its `type` and `id`.
    pub fn put(&self, record: &Value) -> Result<(), String> {
        let kind = record["type"].as_str().ok_or("a record without a type")?;
        let id = record["id"].as_str().ok_or("a record without an id")?;
        let text = serde_json::to_string(record).map_err(|e| e.to_string())?;
        write_atomic(&self.record_path(kind, id), text.as_bytes())
    }

    /// The record of this type and id, if the store has it.
    pub fn get(&self, kind: &str, id: &str) -> Result<Option<Value>, String> {
        let path = self.record_path(kind, id);
        match std::fs::read_to_string(&path) {
            Ok(s) => serde_json::from_str(&s)
                .map(Some)
                .map_err(|e| format!("{}: {e}", path.display())),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(format!("{}: {e}", path.display())),
        }
    }

    /// The record of this type and id, or an error saying it is missing.
    pub fn need(&self, kind: &str, id: &str) -> Result<Value, String> {
        self.get(kind, id)?
            .ok_or_else(|| format!("no {kind} {} in the store", short(id)))
    }

    /// Every stored record whose id starts with `prefix`, of any type.
    pub fn find(&self, prefix: &str) -> Result<Vec<Value>, String> {
        if prefix.len() < 2 {
            return Err(format!("{prefix:?}: give at least two characters of an id"));
        }
        let (a, b) = prefix.split_at(2);
        let mut out = Vec::new();
        for kind in TYPES {
            let dir = self.root.join("records").join(kind).join(a);
            let Ok(entries) = std::fs::read_dir(&dir) else {
                continue;
            };
            let mut names: Vec<String> = entries
                .filter_map(|e| e.ok()?.file_name().into_string().ok())
                .filter(|n| n.starts_with(b) && n.ends_with(".json"))
                .collect();
            names.sort();
            for n in names {
                let id = format!("{a}{}", n.trim_end_matches(".json"));
                out.extend(self.get(kind, &id)?);
            }
        }
        Ok(out)
    }

    /// Store bytes, gzipped, by the SHA-256 of the bytes; returns it.
    pub fn put_blob(&self, bytes: &[u8]) -> Result<String, String> {
        let id = sha256(bytes);
        let path = split(&self.root.join("blobs"), &id, "gz");
        if !path.exists() {
            let mut gz = GzEncoder::new(Vec::new(), flate2::Compression::default());
            gz.write_all(bytes).map_err(|e| e.to_string())?;
            write_atomic(&path, &gz.finish().map_err(|e| e.to_string())?)?;
        }
        Ok(id)
    }

    pub fn blob(&self, id: &str) -> Result<Vec<u8>, String> {
        let path = split(&self.root.join("blobs"), id, "gz");
        let f = std::fs::File::open(&path).map_err(|e| format!("blob {}: {e}", short(id)))?;
        let mut out = Vec::new();
        GzDecoder::new(f)
            .read_to_end(&mut out)
            .map_err(|e| format!("blob {}: {e}", short(id)))?;
        Ok(out)
    }

    /// Store a JSON value as a blob; returns `{"blob": id}`.
    pub fn put_json_blob(&self, v: &Value) -> Result<Value, String> {
        let id = self.put_blob(v.to_string().as_bytes())?;
        Ok(serde_json::json!({ "blob": id }))
    }
}

/// An id's first twelve characters, for messages.
pub fn short(id: &str) -> &str {
    &id[..12.min(id.len())]
}

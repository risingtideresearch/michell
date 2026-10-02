//! Streams: JSON records, one after another, each stream carrying every
//! record its records refer to, before them. Records are keyed by type and
//! id (a result's id is its study's); a record seen twice is the same
//! record, and the first is kept.

use serde_json::Value;
use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};
use std::io::{Read, Write};
use std::path::PathBuf;

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

/// An id's first twelve characters, for messages.
pub fn short(id: &str) -> &str {
    &id[..12.min(id.len())]
}

/// The keys that hold a parent's id, and the parent's type.
pub const LINKS: [(&str, &str); 9] = [
    ("hull", "hull"),
    ("parent", "*"),
    ("case", "case"),
    ("study", "study"),
    ("calm", "study"),
    ("result", "result"),
    ("sections", "sections"),
    ("field", "field"),
    ("prop", "prop"),
];

/// The parent type a link key names in a record of type `kind` (`parent`
/// is a hull's parent hull or a case's parent case).
pub fn link_type<'a>(key: &str, kind: &'a str) -> Option<&'a str> {
    LINKS
        .iter()
        .find(|(k, _)| *k == key)
        .map(|(_, t)| match *t {
            "*" => kind,
            t => t,
        })
}

/// Every `(type, id)` a record refers to.
pub fn parents(r: &Value) -> Vec<(String, String)> {
    let kind = r["type"].as_str().unwrap_or("");
    LINKS
        .iter()
        .filter_map(|(key, _)| {
            let id = r.get(*key)?.as_str()?;
            Some((link_type(key, kind)?.to_string(), id.to_string()))
        })
        .collect()
}

fn key(r: &Value) -> Option<(String, String)> {
    Some((
        r["type"].as_str()?.to_string(),
        r["id"].as_str()?.to_string(),
    ))
}

#[derive(Default)]
pub struct Stream {
    pub records: Vec<Value>,
    at: HashMap<(String, String), usize>,
}

impl Stream {
    /// The records of these files, or of stdin with none.
    pub fn read(files: &[PathBuf]) -> Result<Stream, String> {
        let mut s = Stream::default();
        let mut add = |text: &str, what: &str| -> Result<(), String> {
            for v in serde_json::Deserializer::from_str(text).into_iter::<Value>() {
                s.push(v.map_err(|e| format!("{what}: {e}"))?);
            }
            Ok(())
        };
        if files.is_empty() {
            let mut text = String::new();
            std::io::stdin()
                .read_to_string(&mut text)
                .map_err(|e| format!("stdin: {e}"))?;
            add(&text, "stdin")?;
        } else {
            for f in files {
                let text =
                    std::fs::read_to_string(f).map_err(|e| format!("{}: {e}", f.display()))?;
                add(&text, &f.display().to_string())?;
            }
        }
        Ok(s)
    }

    /// Add a record unless it's already here; true if it was new.
    pub fn push(&mut self, r: Value) -> bool {
        let Some(k) = key(&r) else {
            // A record without a type or id is kept, but can't be found.
            self.records.push(r);
            return true;
        };
        if self.at.contains_key(&k) {
            return false;
        }
        self.at.insert(k, self.records.len());
        self.records.push(r);
        true
    }

    pub fn get(&self, kind: &str, id: &str) -> Option<&Value> {
        self.at
            .get(&(kind.to_string(), id.to_string()))
            .map(|&i| &self.records[i])
    }

    /// The record, or an error saying the stream lacks it.
    pub fn need(&self, kind: &str, id: &str) -> Result<&Value, String> {
        self.get(kind, id).ok_or_else(|| {
            format!(
                "{kind} {} is not in the stream (a filter may have dropped it: \
                 `boatmath pick` keeps ancestors)",
                short(id)
            )
        })
    }

    pub fn of_type(&self, kind: &str) -> impl Iterator<Item = &Value> + '_ {
        let kind = kind.to_string();
        self.records
            .iter()
            .filter(move |r| r["type"] == kind.as_str())
    }

    /// The record a link key in `r` names.
    pub fn follow(&self, r: &Value, key: &str) -> Result<&Value, String> {
        let kind = r["type"].as_str().unwrap_or("");
        let t = link_type(key, kind).ok_or_else(|| format!("{key}: not a link"))?;
        let id = r[key]
            .as_str()
            .ok_or_else(|| format!("a {kind} without a {key}"))?;
        self.need(t, id)
    }

    /// `ids` (as `(type, id)`) and all their ancestors, in stream order.
    pub fn with_ancestors(&self, picked: &[usize]) -> Vec<usize> {
        let mut keep: HashSet<usize> = HashSet::new();
        let mut todo: Vec<usize> = picked.to_vec();
        while let Some(i) = todo.pop() {
            if !keep.insert(i) {
                continue;
            }
            for (t, id) in parents(&self.records[i]) {
                if let Some(&j) = self.at.get(&(t, id)) {
                    todo.push(j);
                }
            }
        }
        let mut out: Vec<usize> = keep.into_iter().collect();
        out.sort_unstable();
        out
    }
}

/// Where a command's stream goes: every record once, the input first.
pub struct Out<W: Write> {
    w: W,
    seen: HashSet<(String, String)>,
}

impl<W: Write> Out<W> {
    pub fn new(w: W) -> Self {
        Out {
            w,
            seen: HashSet::new(),
        }
    }

    /// Write a record unless it's been written already.
    pub fn emit(&mut self, r: &Value) -> Result<(), String> {
        if let Some(k) = key(r) {
            if !self.seen.insert(k) {
                return Ok(());
            }
        }
        writeln!(self.w, "{r}").map_err(crate::stdout_err)?;
        self.w.flush().map_err(crate::stdout_err)
    }

    /// Write a whole stream through.
    pub fn pass(&mut self, s: &Stream) -> Result<(), String> {
        for r in &s.records {
            self.emit(r)?;
        }
        Ok(())
    }

    pub fn raw(&mut self) -> &mut W {
        &mut self.w
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn stream() -> Stream {
        let mut s = Stream::default();
        for r in [
            json!({ "type": "hull", "id": "h1" }),
            json!({ "type": "case", "id": "c1", "hull": "h1" }),
            json!({ "type": "case", "id": "c2", "hull": "h1" }),
            json!({ "type": "study", "id": "s1", "case": "c2" }),
            json!({ "type": "result", "id": "s1", "study": "s1", "field": "s1" }),
            json!({ "type": "field", "id": "s1", "study": "s1" }),
        ] {
            assert!(s.push(r));
        }
        s
    }

    /// The same record twice is one; a result shares its study's id but not
    /// its type.
    #[test]
    fn records_are_keyed_by_type_and_id() {
        let mut s = stream();
        assert!(!s.push(json!({ "type": "hull", "id": "h1", "name": "again" })));
        assert_eq!(s.records.len(), 6);
        assert_eq!(s.get("result", "s1").unwrap()["field"], "s1");
        assert_eq!(s.get("study", "s1").unwrap()["case"], "c2");
        assert!(s
            .need("case", "c9")
            .unwrap_err()
            .contains("not in the stream"));
    }

    /// A result's ancestors are its study, field, case and hull, not the
    /// other case.
    #[test]
    fn ancestors_follow_the_links() {
        let s = stream();
        let kept = s.with_ancestors(&[4]);
        let ids: Vec<String> = kept
            .iter()
            .map(|&i| {
                let r = &s.records[i];
                format!(
                    "{}:{}",
                    r["type"].as_str().unwrap(),
                    r["id"].as_str().unwrap()
                )
            })
            .collect();
        assert_eq!(
            ids,
            ["hull:h1", "case:c2", "study:s1", "result:s1", "field:s1"]
        );
    }
}

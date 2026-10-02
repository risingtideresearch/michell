//! Field paths into records, as `table` and `plot` name their columns:
//!
//! - `forces.rt`, `seakeeping.headings.0.points` — keys and array indices;
//! - a parent's id is followed into its record, so on a result
//!   `study.case.params.span` reads the span of the result's case
//!   (`hull`, `case`, `study`, `parent`, `calm` and `sections` are
//!   followed);
//! - `|heave|` — the modulus of a complex `[re, im]` pair;
//! - `id8` — a record's id, shortened to eight characters, as a label.

use crate::store::Store;
use serde_json::{json, Value};
use std::cell::RefCell;
use std::collections::HashMap;

/// The record type a key holding an id refers to.
fn parent_type(key: &str) -> Option<&'static str> {
    match key {
        "hull" | "parent" => Some("hull"),
        "case" => Some("case"),
        "study" | "calm" => Some("study"),
        "sections" => Some("sections"),
        _ => None,
    }
}

pub struct Resolver<'a> {
    store: &'a Store,
    cache: RefCell<HashMap<(String, String), Option<Value>>>,
}

impl<'a> Resolver<'a> {
    pub fn new(store: &'a Store) -> Self {
        Resolver {
            store,
            cache: RefCell::new(HashMap::new()),
        }
    }

    fn load(&self, kind: &str, id: &str) -> Option<Value> {
        let key = (kind.to_string(), id.to_string());
        if let Some(v) = self.cache.borrow().get(&key) {
            return v.clone();
        }
        let v = self.store.get(kind, id).ok().flatten();
        self.cache.borrow_mut().insert(key, v.clone());
        v
    }

    /// The value at `path` in `root`, or `None` if there is none.
    pub fn get(&self, root: &Value, path: &str) -> Option<Value> {
        if let Some(inner) = path.strip_prefix('|').and_then(|p| p.strip_suffix('|')) {
            let v = self.get(root, inner)?;
            let a = v.as_array()?;
            return match a[..] {
                [ref re, ref im] => Some(json!(re.as_f64()?.hypot(im.as_f64()?))),
                _ => None,
            };
        }
        if path == "id8" {
            return root["id"].as_str().map(|s| json!(&s[..8.min(s.len())]));
        }
        let mut cur = root.clone();
        let segments: Vec<&str> = path.split('.').collect();
        for (i, seg) in segments.iter().enumerate() {
            let next = match &cur {
                Value::Array(a) => seg.parse::<usize>().ok().and_then(|i| a.get(i).cloned()),
                Value::Object(o) => o.get(*seg).cloned(),
                _ => None,
            }?;
            // An id with more path after it is followed into its record.
            cur = match (&next, parent_type(seg)) {
                (Value::String(id), Some(kind)) if i + 1 < segments.len() => self.load(kind, id)?,
                _ => next,
            };
        }
        Some(cur)
    }

    /// The value at `path` in `row`, else in `record` (the record a row was
    /// exploded from).
    pub fn get_in(&self, row: &Value, record: &Value, path: &str) -> Option<Value> {
        self.get(row, path)
            .filter(|v| !v.is_null())
            .or_else(|| self.get(record, path))
    }

    /// The rows of `record`: itself, or each element of the array at
    /// `explode`.
    pub fn rows(&self, record: &Value, explode: Option<&str>) -> Vec<Value> {
        match explode {
            None => vec![record.clone()],
            Some(p) => match self.get(record, p) {
                Some(Value::Array(a)) => a,
                _ => Vec::new(),
            },
        }
    }
}

/// A value as a table cell: numbers and strings plain, nothing empty, the
/// rest as JSON.
pub fn cell(v: &Option<Value>) -> String {
    match v {
        None | Some(Value::Null) => String::new(),
        Some(Value::String(s)) => s.clone(),
        Some(Value::Number(n)) => n.to_string(),
        Some(v) => v.to_string(),
    }
}

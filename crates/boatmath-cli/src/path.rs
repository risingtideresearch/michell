//! Field paths into records, as `table` and `plot` name their columns:
//!
//! - `forces.rt`, `seakeeping.headings.0.points` — keys and array indices;
//! - a parent's id is followed into its record in the stream, so on a
//!   result `study.case.params.span` reads the span of the result's case
//!   (the links of [`crate::stream::LINKS`] are followed, and `motor` into
//!   the motor database);
//! - `|heave|` — the modulus of a complex `[re, im]` pair;
//! - `id8` — a record's id, shortened to eight characters, as a label.

use crate::stream::{link_type, Stream};
use serde_json::{json, Value};

pub struct Resolver<'a> {
    stream: &'a Stream,
}

impl<'a> Resolver<'a> {
    pub fn new(stream: &'a Stream) -> Self {
        Resolver { stream }
    }

    fn load(&self, kind: &str, id: &str) -> Option<Value> {
        // Motors are the database's, not the stream's.
        if kind == "motor" {
            return propeller::motor::Database::vendored()
                .motor(id)
                .map(crate::props::motor_record);
        }
        self.stream.get(kind, id).cloned()
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
            let kind = cur["type"].as_str().unwrap_or("").to_string();
            let next = match &cur {
                Value::Array(a) => seg.parse::<usize>().ok().and_then(|i| a.get(i).cloned()),
                Value::Object(o) => o.get(*seg).cloned(),
                _ => None,
            }?;
            // An id with more path after it is followed into its record.
            let parent = if *seg == "motor" {
                Some("motor")
            } else {
                link_type(seg, &kind)
            };
            cur = match (&next, parent) {
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
        // A whole number written as a float (a blade count, say) as an integer.
        Some(Value::Number(n)) => match n.as_f64() {
            Some(x) if n.is_f64() && x.fract() == 0.0 && x.abs() < 1e15 => format!("{}", x as i64),
            _ => n.to_string(),
        },
        Some(v) => v.to_string(),
    }
}

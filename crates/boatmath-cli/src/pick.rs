//! `pick`: filter a stream by record, keeping every record the kept ones
//! refer to, so what comes out is still a whole stream.

use crate::path::{cell, Resolver};
use crate::stream::Stream;

pub struct Pick {
    /// Keep records of this type (with their ancestors); none: every type.
    pub kind: Option<String>,
    /// Keep only the first N of them, in stream order.
    pub first: Option<usize>,
    /// `PATH=VALUE`: keep a record only where each path reads that value.
    pub wheres: Vec<(String, String)>,
    /// Leave records of these types out altogether (ancestors included).
    pub drop: Vec<String>,
}

pub fn parse_where(s: &str) -> Result<(String, String), String> {
    s.split_once('=')
        .map(|(p, v)| (p.trim().to_string(), v.trim().to_string()))
        .ok_or_else(|| format!("--where {s:?}: expected PATH=VALUE"))
}

/// A value as the text a `--where` compares against: numbers equal as
/// numbers (`0.4` matches `0.40`), anything else as its table cell.
fn matches(v: &Option<serde_json::Value>, want: &str) -> bool {
    match (v.as_ref().and_then(|v| v.as_f64()), want.parse::<f64>()) {
        (Some(a), Ok(b)) => (a - b).abs() <= 1e-9 * b.abs().max(1.0),
        _ => cell(v) == want,
    }
}

/// The indices of the records `p` keeps, in stream order.
pub fn pick(s: &Stream, p: &Pick) -> Vec<usize> {
    let res = Resolver::new(s);
    let mut chosen: Vec<usize> = s
        .records
        .iter()
        .enumerate()
        .filter(|(_, r)| {
            let t = r["type"].as_str().unwrap_or("");
            p.kind.as_deref().is_none_or(|k| k == t)
                && !p.drop.iter().any(|d| d == t)
                && p.wheres
                    .iter()
                    .all(|(path, want)| matches(&res.get(r, path), want))
        })
        .map(|(i, _)| i)
        .collect();
    if let Some(n) = p.first {
        chosen.truncate(n);
    }
    let kept = if p.kind.is_some() || !p.wheres.is_empty() || p.first.is_some() {
        s.with_ancestors(&chosen)
    } else {
        chosen
    };
    kept.into_iter()
        .filter(|&i| !p.drop.iter().any(|d| s.records[i]["type"] == d.as_str()))
        .collect()
}

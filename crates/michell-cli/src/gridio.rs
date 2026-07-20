//! Sample-grid JSON interchange (`*.grid.json`).
//!
//! Serializes the library's [`SampleGrid`] IR — the station × waterline
//! half-beam samples every importer produces, with optional `∂f/∂x`, `∂f/∂z`
//! and weight channels — so sampled geometry can be inspected, diffed, and
//! re-lofted without the source CAD file. Unknown derivative entries (NaN)
//! are written as JSON `null`.
//!
//! ```json
//! {
//!   "michell": "sample-grid",
//!   "version": 1,
//!   "centerplane": 0.0,
//!   "stations": [ ... ],
//!   "waterlines": [ ... ],
//!   "half_beams": [ ... ],          // row-major, waterline index fastest
//!   "dfdx": [ ..., null, ... ],     // optional
//!   "dfdz": [ ... ],                // optional
//!   "weights": [ ... ]              // optional; 0 excludes a sample
//! }
//! ```

use crate::json::{self, Json};
use michell::SampleGrid;

fn arr(v: &[f64]) -> String {
    let mut s = String::from("[");
    for (i, x) in v.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        if x.is_finite() {
            s.push_str(&format!("{x}"));
        } else {
            s.push_str("null");
        }
    }
    s.push(']');
    s
}

/// Serialize a grid (and the hull's transverse centerplane position, so a
/// multihull member re-lofts at its detected placement).
pub fn write_grid_json(g: &SampleGrid, centerplane: Option<f64>) -> String {
    let mut out = String::from("{\n  \"michell\": \"sample-grid\",\n  \"version\": 1");
    if let Some(y) = centerplane {
        out.push_str(&format!(",\n  \"centerplane\": {y}"));
    }
    out.push_str(&format!(",\n  \"stations\": {}", arr(g.stations())));
    out.push_str(&format!(",\n  \"waterlines\": {}", arr(g.waterlines())));
    out.push_str(&format!(",\n  \"half_beams\": {}", arr(g.half_beams())));
    if let Some(fx) = g.fx() {
        out.push_str(&format!(",\n  \"dfdx\": {}", arr(fx)));
    }
    if let Some(fz) = g.fz() {
        out.push_str(&format!(",\n  \"dfdz\": {}", arr(fz)));
    }
    if let Some(w) = g.weights() {
        out.push_str(&format!(",\n  \"weights\": {}", arr(w)));
    }
    out.push_str("\n}\n");
    out
}

/// Parse a sample-grid JSON file: the grid plus the stored centerplane
/// position, if any.
pub fn parse_grid_json(text: &str) -> Result<(SampleGrid, Option<f64>), String> {
    let j = json::parse(text)?;
    if j.get("michell").and_then(Json::as_str) != Some("sample-grid") {
        return Err(
            "not a michell sample-grid file (expected \"michell\": \"sample-grid\")".into(),
        );
    }
    match j.get("version").and_then(Json::as_f64) {
        Some(1.0) => {}
        other => return Err(format!("unsupported sample-grid version {other:?}")),
    }
    let numbers = |key: &str| -> Result<Vec<f64>, String> {
        let arr = j
            .get(key)
            .and_then(Json::as_arr)
            .ok_or_else(|| format!("missing array {key:?}"))?;
        arr.iter()
            .map(|v| match v {
                Json::Num(n) => Ok(*n),
                Json::Null => Ok(f64::NAN),
                _ => Err(format!("{key:?} must contain only numbers or null")),
            })
            .collect()
    };
    let mut grid = SampleGrid::new(
        numbers("stations")?,
        numbers("waterlines")?,
        numbers("half_beams")?,
    )
    .map_err(|e| format!("{e}"))?;
    if j.get("dfdx").is_some() {
        grid = grid.with_fx(numbers("dfdx")?).map_err(|e| format!("{e}"))?;
    }
    if j.get("dfdz").is_some() {
        grid = grid.with_fz(numbers("dfdz")?).map_err(|e| format!("{e}"))?;
    }
    if j.get("weights").is_some() {
        grid = grid
            .with_weights(numbers("weights")?)
            .map_err(|e| format!("{e}"))?;
    }
    let centerplane = match j.get("centerplane") {
        None | Some(Json::Null) => None,
        Some(v) => Some(
            v.as_f64()
                .ok_or_else(|| "\"centerplane\" must be a number".to_string())?,
        ),
    };
    Ok((grid, centerplane))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_all_channels() {
        let grid = SampleGrid::new(
            vec![0.0, 1.0, 2.0],
            vec![0.0, 0.5],
            vec![0.1, 0.0, 0.4, 0.2, 0.3, 0.0],
        )
        .unwrap()
        .with_fx(vec![0.5, f64::NAN, 0.25, -0.5, f64::NAN, 0.0])
        .unwrap()
        .with_weights(vec![1.0, 1.0, 0.0, 1.0, 1.0, 1.0])
        .unwrap();
        let text = write_grid_json(&grid, Some(-1.25));
        let (back, y) = parse_grid_json(&text).unwrap();
        assert_eq!(y, Some(-1.25));
        assert_eq!(back.stations(), grid.stations());
        assert_eq!(back.waterlines(), grid.waterlines());
        assert_eq!(back.half_beams(), grid.half_beams());
        assert_eq!(back.weights(), grid.weights());
        assert!(back.fz().is_none());
        let (a, b) = (back.fx().unwrap(), grid.fx().unwrap());
        for (x, y) in a.iter().zip(b) {
            assert!((x.is_nan() && y.is_nan()) || x == y);
        }
    }

    #[test]
    fn rejects_wrong_kind_and_shape() {
        assert!(parse_grid_json("{\"michell\":\"manifest\"}").is_err());
        assert!(parse_grid_json(
            "{\"michell\":\"sample-grid\",\"version\":2,\"stations\":[],\
             \"waterlines\":[],\"half_beams\":[]}"
        )
        .is_err());
        // Null in a required value channel is rejected by grid validation.
        assert!(parse_grid_json(
            "{\"michell\":\"sample-grid\",\"version\":1,\"stations\":[0,1],\
             \"waterlines\":[0,1],\"half_beams\":[0,null,0,0]}"
        )
        .is_err());
    }
}

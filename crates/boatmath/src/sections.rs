//! A case's hulls cut into sections at one attitude, as JSON: each station
//! exactly as the solver integrates it (see
//! [`michell_geometry::iges::PolarSection`]).
//!
//! ```json
//! { "attitude": { "sinkage": 0.008, "trim_rad": 0.0015, "trim_deg": 0.087 },
//!   "hulls": [ { "placement": { "x": 0.0, "y": 0.0 }, "transom": true,
//!                "stations": [ { "x": 0.31, "z0": 0.0, "beam": 0.21, "depth": 0.18,
//!                                "radii": [1.0, 0.998, …] }, … ] } ] }
//! ```
//!
//! A station's section is `y = beam·R(θ)·cos θ`, `z = z0 + depth·R(θ)·sin θ`
//! (`y` the half-breadth from the hull's centreplane, `z` the depth below
//! the water), `R` the polynomial through `radii` at the Chebyshev–Lobatto
//! angles `θ_k = ¼π(1 − cos(πk/(n−1)))`. A station with no radii lies past
//! the hull's tip. With `transom`, the aft station is a transom.

use crate::params::CaseParams;
use crate::{setup, LoftRequest};
use michell_geometry::iges::PolarSection;
use serde_json::{json, Value};

fn station(p: &PolarSection) -> Value {
    json!({ "x": p.x, "z0": p.z0, "beam": p.beam, "depth": p.depth, "radii": p.radii })
}

/// The case's hulls cut at `(sinkage [m], trim [rad])`, as their sections.
pub fn sections(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
    attitude: (f64, f64),
) -> Result<Value, String> {
    let s = setup(name, bytes, cut, c)?;
    let (sinkage, trim) = attitude;
    let platform = s.platform(sinkage, trim);
    let mut hulls = Vec::new();
    for (i, pose) in &s.layout {
        let h = s
            .cut
            .file
            .source
            .situate_sectional(*i, s.cut.file.waterline_z, pose, &platform, &s.cut.opts)
            .map_err(|e| e.to_string())?
            .ok_or("a hull is dry at this attitude")?;
        hulls.push(json!({
            "placement": { "x": h.placement.x, "y": h.placement.y },
            "transom": h.hull.transom().is_some(),
            "stations": h.polar.iter().map(station).collect::<Vec<_>>(),
        }));
    }
    Ok(json!({
        "attitude": { "sinkage": sinkage, "trim_rad": trim, "trim_deg": trim.to_degrees() },
        "hulls": hulls,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wigley() -> Vec<u8> {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        michell_geometry::iges::write(&surfaces, "wigley")
            .unwrap()
            .into_bytes()
    }

    /// Each station's half-area from its polar form, `½·beam·depth·∫R² dθ`,
    /// is the one the solver integrates (its depth integral at κ = 0).
    #[test]
    fn the_polar_form_is_the_section_integrated() {
        let cut = crate::cut("w.igs", wigley(), &LoftRequest::default()).unwrap();
        let imp = &cut.hulls[0];
        let (at_stations, _) = imp.hull.depth_integral_curve(0.0, 1);
        assert_eq!(at_stations.len(), imp.polar.len());
        let (gx, gw) = (64, std::f64::consts::FRAC_PI_2);
        for (p, &(x, z)) in imp.polar.iter().zip(&at_stations) {
            assert!((p.x - x).abs() < 1e-12);
            if p.radii.is_empty() {
                assert!(z.abs() < 1e-9);
                continue;
            }
            let t = p.angles();
            // Midpoint rule on R(θ)² through the polar points.
            let r2: f64 = (0..gx)
                .map(|k| {
                    let th = (k as f64 + 0.5) / gx as f64 * gw;
                    let (y, zz) = p.point(th);
                    let (dy, dz) = (y / p.beam, (zz - p.z0) / p.depth);
                    dy * dy + dz * dz
                })
                .sum::<f64>()
                * gw
                / gx as f64;
            let area = 0.5 * p.beam * p.depth * r2;
            assert!(t.len() == p.radii.len());
            assert!(
                (area - z).abs() < 1e-3 * z.abs().max(1e-6),
                "x {x}: {area} vs {z}"
            );
        }
        let v = sections(
            "w.igs",
            wigley(),
            &LoftRequest::default(),
            &CaseParams::default(),
            (0.0, 0.0),
        )
        .unwrap();
        assert_eq!(v["hulls"].as_array().unwrap().len(), 1);
        assert_eq!(
            v["hulls"][0]["stations"].as_array().unwrap().len(),
            imp.polar.len()
        );
    }
}

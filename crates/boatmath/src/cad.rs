//! A calm-water result as CAD surfaces: the hulls at their attitude (their
//! own patches, re-posed exactly), the free surface through every point of
//! the wave grid, and propeller discs. Frame: x forward, y to port, z up,
//! the still water at z = 0, metres.

use crate::params::CaseParams;
use crate::{setup, LoftRequest};
use michell::propulsion::Disc;
use michell_geometry::iges::NurbsSurface3;

/// Each of the case's hulls as its surface patches at `(sinkage [m],
/// trim [rad])`.
pub fn hull_surfaces(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
    attitude: (f64, f64),
) -> Result<Vec<Vec<NurbsSurface3>>, String> {
    let s = setup(name, bytes, cut, c)?;
    let platform = s.platform(attitude.0, attitude.1);
    s.layout
        .iter()
        .map(|(i, pose)| {
            s.cut
                .file
                .source
                .posed_surfaces(*i, s.cut.file.waterline_z, pose, &platform)
                .ok_or("a mesh hull (from an STL) has no surfaces to write")?
                .map_err(|e| e.to_string())
        })
        .collect()
}

/// The free surface `z = scale · ζ(x, y)` through every point of its grid.
#[allow(clippy::too_many_arguments)]
pub fn water_surface(
    x0: f64,
    x1: f64,
    nx: usize,
    y0: f64,
    y1: f64,
    ny: usize,
    zeta: &[f64],
    scale: f64,
) -> Result<NurbsSurface3, String> {
    let z: Vec<f64> = zeta.iter().map(|v| scale * v).collect();
    michell_geometry::iges::graph_surface(x0, x1, nx, y0, y1, ny, &z).map_err(|e| e.to_string())
}

/// A propeller disc as a flat annulus, hub to tip, square to the shaft:
/// a rational quadratic circle swept radially, exact.
pub fn disc_surface(d: &Disc) -> NurbsSurface3 {
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let square = [
        (1.0, 0.0, 1.0),
        (1.0, 1.0, h),
        (0.0, 1.0, 1.0),
        (-1.0, 1.0, h),
        (-1.0, 0.0, 1.0),
        (-1.0, -1.0, h),
        (0.0, -1.0, 1.0),
        (1.0, -1.0, h),
        (1.0, 0.0, 1.0),
    ];
    let (yc, zc) = (d.y, -d.depth);
    let radii = [d.hub * d.radius, d.radius];
    let mut ctrl = Vec::new();
    let mut weights = Vec::new();
    for &(cy, cz, w) in &square {
        for &r in &radii {
            ctrl.push([d.x, yc + r * cy, zc + r * cz]);
            weights.push(w);
        }
    }
    NurbsSurface3 {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        n_ctrl_u: 9,
        n_ctrl_v: 2,
        ctrl,
        weights,
        trim_uv: None,
    }
}

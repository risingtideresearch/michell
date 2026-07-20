//! Reference hull constructors.

use crate::bspline::BSplineSurface;
use crate::error::Result;
use crate::hull::Hull;

/// The standard Wigley parabolic test hull,
/// `f(x, z) = (B/2) (1 - (2x/L)²) (1 - (z/T)²)`,
/// represented **exactly** as a single-span biquadratic B-spline.
///
/// `x ∈ [-L/2, L/2]`, `z ∈ [0, T]`. The classic proportions are
/// `L/B = 10`, `B/T = 1.6` (e.g. `wigley(10.0, 1.0, 0.625)`).
pub fn wigley(length: f64, beam: f64, draft: f64) -> Result<Hull> {
    let a = length / 2.0;
    let knots_x = vec![-a, -a, -a, a, a, a];
    let knots_z = vec![0.0, 0.0, 0.0, draft, draft, draft];
    // Quadratic Bezier control values reproducing 4t(1-t) in x and 1-t² in z.
    let gx = [0.0, 2.0, 0.0];
    let hz = [1.0, 1.0, 0.0];
    let mut control = Vec::with_capacity(9);
    for gi in gx {
        for hj in hz {
            control.push(beam / 2.0 * gi * hj);
        }
    }
    Hull::new(BSplineSurface::new(2, 2, knots_x, knots_z, control)?)
}

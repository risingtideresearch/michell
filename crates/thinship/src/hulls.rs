//! Reference hulls for the exact B-spline oracle.

use crate::hull::Hull;
use hullgeom::hulls::wigley_surface;
use hullgeom::Result;

/// The Wigley hull as the exact B-spline oracle's [`Hull`].
pub(crate) fn wigley(length: f64, beam: f64, draft: f64) -> Result<Hull> {
    Hull::new(wigley_surface(length, beam, draft)?)
}

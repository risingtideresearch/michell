//! Hydrostatic restoring coefficients for heave and pitch.
//!
//! A heave `η₃` (up) and pitch `η₅` (bow up, about `x_ref`) lift the hull
//! section at `x` by `η₃ + (x − x_ref) η₅`, which the waterplane resists:
//!
//! ```text
//! F₃ = −C₃₃ η₃ − C₃₅ η₅,     M₅ = −C₅₃ η₃ − C₅₅ η₅,
//! C₃₃ = ρ g A_w,  C₃₅ = C₅₃ = ρ g ∫ (x − x_ref) b dx,  C₅₅ = ρ g ∫ (x − x_ref)² b dx
//! ```
//!
//! with `b = 2 f(x, 0)` the waterline beam. `C₅₅` here is the waterplane's
//! part, `ρ g I_L`; the full pitch stiffness adds `ρ g ∇ (z_B − z_G)`
//! (`z` up), which is ∇·BG/I_L — a percent or so — of it on a slender hull,
//! and is left to the caller until the geometry reports a vertical centre of
//! buoyancy.

use michell_geometry::SectionalHull;

/// Heave–pitch hydrostatic restoring coefficients about a pivot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Restoring {
    /// `C₃₃` [N/m].
    pub c33: f64,
    /// `C₃₅ = C₅₃` [N/rad = N·m/m].
    pub c35: f64,
    /// `C₅₅` [N·m/rad] (waterplane part; see the module docs).
    pub c55: f64,
}

/// The restoring coefficients of `hull` about `x_ref` (in the hull's own x),
/// in a fluid of `density` under `gravity`.
pub fn restoring(hull: &SectionalHull, density: f64, gravity: f64, x_ref: f64) -> Restoring {
    let rg = density * gravity;
    let (a, m, i) = (
        hull.waterplane_area(),
        hull.waterplane_moment(),
        hull.waterplane_second_moment(),
    );
    Restoring {
        c33: rg * a,
        c35: rg * (m - x_ref * a),
        c55: rg * (i - 2.0 * x_ref * m + x_ref * x_ref * a),
    }
}

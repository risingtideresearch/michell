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

use michell_geometry::{Placement, SectionalHull};

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

/// The platform's transverse metacentric height `GM_T = BM_T − BG` [m] for
/// roll about its centreplane (`y = 0`): `BM_T = I_T/∇`, the waterplanes'
/// transverse second moment taken about the platform centreplane — each
/// hull's own `∫⅔f³dx` plus `A_w y²` for its offset — so a catamaran's comes
/// almost entirely from its spacing. `bg` is the height of the centre of
/// gravity above the centre of buoyancy.
pub fn transverse_metacentric_height(members: &[(&SectionalHull, Placement)], bg: f64) -> f64 {
    let (mut it, mut vol) = (0.0, 0.0);
    for (hull, pl) in members {
        it += hull.transverse_waterplane_inertia() + hull.waterplane_area() * pl.y * pl.y;
        vol += hull.displaced_volume();
    }
    if vol > 0.0 {
        it / vol - bg
    } else {
        0.0
    }
}

/// Depth of the platform's centre of buoyancy below the waterline [m]
/// (volume-weighted over its hulls).
pub fn buoyancy_depth(members: &[(&SectionalHull, Placement)]) -> f64 {
    let (mut m, mut vol) = (0.0, 0.0);
    for (hull, _) in members {
        let v = hull.displaced_volume();
        m += v * hull.buoyancy_depth();
        vol += v;
    }
    if vol > 0.0 {
        m / vol
    } else {
        0.0
    }
}

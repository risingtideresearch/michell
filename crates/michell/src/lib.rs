//! # michell
//!
//! Thin-ship **wave resistance** (Michell's integral) and **viscous
//! resistance** (ITTC-57) for hulls described as clamped, polynomial
//! (non-rational) tensor-product B-spline half-breadth surfaces.
//!
//! ## Geometry contract
//!
//! The hull is port/starboard symmetric and given by its half-beam
//! `y = f(x, z) >= 0` on the centerplane:
//!
//! - `x` runs along the hull (arbitrary origin), in metres;
//! - `z` runs vertically **downward** from the undisturbed waterline, so the
//!   surface domain is `z ∈ [0, T]` with `T` the draft;
//! - all quantities are SI.
//!
//! ## Theory
//!
//! With `ν = g/U²` (Tuck 1989; Dambrine, Pierre & Rousseaux 2016):
//!
//! ```text
//! R_w = (4 ρ g²)/(π U²) ∫_1^∞ (I² + J²) λ²/√(λ² − 1) dλ
//! I(λ) + i J(λ) = ∬ (∂f/∂x) exp(−ν λ² z) exp(i ν λ x) dx dz
//! ```
//!
//! Because `f` is piecewise polynomial, the inner integrals are evaluated in
//! closed form on every knot span (polynomial × oscillatory / exponential
//! moments); only the smooth outer integral is quadratured, with panels sized
//! to the local oscillation rate and refined to a requested tolerance.
//!
//! Viscous resistance uses the ITTC-57 correlation line with an optional form
//! factor, referencing the thin-ship wetted surface
//! `S = 2 ∬ √(1 + fx² + fz²) dx dz`.
//!
//! ## Assumptions and limits
//!
//! - Thin-ship (Michell) linearisation: slender hull, `|∂f/∂x| ≪ 1`; no
//!   sinkage, trim, or wave-breaking; deep water; infinite fluid extent.
//! - The half-breadth should close at both ends (`f = 0` at bow and stern).
//!   Transom sterns are not yet modelled.
//! - The control net must be non-negative (sufficient condition for
//!   `f >= 0`).
//!
//! ## Example
//!
//! ```
//! use michell::{hulls, Conditions};
//!
//! let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
//! let cond = Conditions::freshwater(3.0); // 3 m/s
//! let r = michell::resistance(&hull, &cond).unwrap();
//! assert!(r.wave.resistance > 0.0 && r.viscous.resistance > 0.0);
//! assert!(r.total > r.wave.resistance);
//! ```

pub mod body;
mod bspline;
mod conditions;
mod error;
pub mod fit;
pub mod float;
mod friction;
mod hull;
pub mod hulls;
pub mod iges;
mod michell;
mod moments;
mod quadrature;
pub mod spectrum;
pub mod stl;

pub use bspline::BSplineSurface;
pub use conditions::{Conditions, Fluid, STANDARD_GRAVITY};
pub use error::{Error, Result};
pub use friction::{ittc57_cf, viscous_resistance, viscous_resistance_with, ViscousResistance};
pub use hull::Hull;
pub use michell::{
    inner_integrals, multihull_wave_resistance, multihull_wave_resistance_with, wave_resistance,
    wave_resistance_with, Placement, WaveOptions, WaveResistance,
};
pub use moments::C64;
pub use spectrum::{FreeWaveSpectrum, WaveGrid};

/// Combined resistance breakdown.
#[derive(Debug, Clone, Copy)]
pub struct Resistance {
    pub wave: WaveResistance,
    pub viscous: ViscousResistance,
    /// R_w + R_v [N].
    pub total: f64,
    /// Effective (towing) power P_E = (R_w + R_v) · U [W].
    pub effective_power: f64,
    /// Wave resistance coefficient C_w = R_w / (½ ρ U² S).
    pub cw: f64,
    /// Viscous resistance coefficient C_v = R_v / (½ ρ U² S).
    pub cv: f64,
    /// Total resistance coefficient C_t = C_w + C_v.
    pub ct: f64,
}

/// Wave + viscous resistance with default options and zero form factor.
pub fn resistance(hull: &Hull, cond: &Conditions) -> Result<Resistance> {
    resistance_with(hull, cond, &WaveOptions::default(), 0.0)
}

/// Combined resistance breakdown for a multihull.
#[derive(Debug, Clone)]
pub struct MultihullResistance {
    /// Combined wave resistance (with interference).
    pub wave: WaveResistance,
    /// Per-member viscous resistance, in input order.
    pub viscous: Vec<ViscousResistance>,
    /// Σ of the members' viscous resistances [N].
    pub viscous_total: f64,
    /// Combined wave + viscous resistance [N].
    pub total: f64,
    /// Effective (towing) power P_E = R_t · U [W].
    pub effective_power: f64,
    /// Σ of the members' wetted surfaces [m²] (reference area for Cw/Cv/Ct).
    pub wetted_surface: f64,
    pub cw: f64,
    pub cv: f64,
    pub ct: f64,
    /// Σ of each member's standalone wave resistance [N].
    pub solo_wave_total: f64,
    /// Wave interference factor: combined R_w / Σ standalone R_w
    /// (< 1 favourable, > 1 unfavourable, 1 when far apart).
    pub interference: f64,
}

/// Multihull resistance with default options and zero form factor.
pub fn multihull_resistance(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
) -> Result<MultihullResistance> {
    multihull_resistance_with(members, cond, &WaveOptions::default(), 0.0)
}

/// Multihull resistance with explicit quadrature options and form factor
/// (applied to every member).
pub fn multihull_resistance_with(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    wave_opts: &WaveOptions,
    form_factor: f64,
) -> Result<MultihullResistance> {
    let wave = multihull_wave_resistance_with(members, cond, wave_opts)?;
    let mut solo_wave_total = 0.0;
    for m in members {
        solo_wave_total +=
            multihull_wave_resistance_with(&[*m], cond, wave_opts)?.resistance;
    }
    let viscous: Vec<ViscousResistance> = members
        .iter()
        .map(|(h, _)| viscous_resistance_with(h, cond, form_factor))
        .collect::<Result<_>>()?;
    let viscous_total: f64 = viscous.iter().map(|v| v.resistance).sum();
    let wetted_surface: f64 = members.iter().map(|(h, _)| h.wetted_surface()).sum();
    let total = wave.resistance + viscous_total;
    let q = 0.5 * cond.fluid.density * cond.speed * cond.speed * wetted_surface;
    Ok(MultihullResistance {
        wave,
        viscous,
        viscous_total,
        total,
        effective_power: total * cond.speed,
        wetted_surface,
        cw: wave.resistance / q,
        cv: viscous_total / q,
        ct: total / q,
        solo_wave_total,
        interference: if solo_wave_total > f64::MIN_POSITIVE {
            wave.resistance / solo_wave_total
        } else {
            1.0
        },
    })
}

/// Wave + viscous resistance with explicit quadrature options and form factor.
pub fn resistance_with(
    hull: &Hull,
    cond: &Conditions,
    wave_opts: &WaveOptions,
    form_factor: f64,
) -> Result<Resistance> {
    let wave = wave_resistance_with(hull, cond, wave_opts)?;
    let viscous = viscous_resistance_with(hull, cond, form_factor)?;
    let q = 0.5 * cond.fluid.density * cond.speed * cond.speed * hull.wetted_surface();
    let total = wave.resistance + viscous.resistance;
    Ok(Resistance {
        wave,
        viscous,
        total,
        effective_power: total * cond.speed,
        cw: wave.resistance / q,
        cv: viscous.resistance / q,
        ct: total / q,
    })
}

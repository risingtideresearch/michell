//! Viscous (skin-friction) resistance via the ITTC-57 correlation line.

use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::hull::Hull;

/// Viscous resistance result.
#[derive(Debug, Clone, Copy)]
pub struct ViscousResistance {
    /// Viscous resistance R_v = (1 + k) · ½ ρ U² S C_F [N].
    pub resistance: f64,
    /// ITTC-57 flat-plate friction coefficient C_F.
    pub cf: f64,
    /// Reynolds number U·L/ν based on the hull length.
    pub reynolds: f64,
    /// Wetted surface S used as reference area [m²].
    pub wetted_surface: f64,
    /// Form factor k applied as (1 + k).
    pub form_factor: f64,
}

/// ITTC-57 correlation line: `C_F = 0.075 / (log10(Re) - 2)²`.
pub fn ittc57_cf(reynolds: f64) -> f64 {
    let d = reynolds.log10() - 2.0;
    0.075 / (d * d)
}

/// Viscous resistance with zero form factor (bare ITTC-57 flat plate).
pub fn viscous_resistance(hull: &Hull, cond: &Conditions) -> Result<ViscousResistance> {
    viscous_resistance_with(hull, cond, 0.0)
}

/// Viscous resistance with a form factor k, i.e. `R_v = (1 + k) R_F`.
pub fn viscous_resistance_with(
    hull: &Hull,
    cond: &Conditions,
    form_factor: f64,
) -> Result<ViscousResistance> {
    cond.validate()?;
    if !form_factor.is_finite() || form_factor < 0.0 {
        return Err(Error::InvalidConditions(format!(
            "form factor must be finite and non-negative, got {form_factor}"
        )));
    }
    let re = cond.speed * hull.length() / cond.fluid.kinematic_viscosity;
    if re <= 1e3 {
        return Err(Error::InvalidConditions(format!(
            "Reynolds number {re:.3e} is outside the ITTC-57 line's sensible range"
        )));
    }
    let cf = ittc57_cf(re);
    let s = hull.wetted_surface();
    let resistance =
        (1.0 + form_factor) * 0.5 * cond.fluid.density * cond.speed * cond.speed * s * cf;
    Ok(ViscousResistance {
        resistance,
        cf,
        reynolds: re,
        wetted_surface: s,
        form_factor,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ittc_line_known_value() {
        // At Re = 1e8: C_F = 0.075 / (8 - 2)^2 = 0.075 / 36.
        assert!((ittc57_cf(1e8) - 0.075 / 36.0).abs() < 1e-15);
    }
}

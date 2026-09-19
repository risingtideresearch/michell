//! Viscous resistance: the ITTC-57 correlation line, a form factor, and a
//! surface-roughness allowance.
//!
//! Following the ITTC-78 shape, the three pieces compose as
//!
//! ```text
//! C_V = (1 + k)·C_F(Re)  +  ΔC_F
//! ```
//!
//! — the form factor multiplies flat-plate friction (it is a *property of the
//! shape*: streamline curvature speeds the flow over most of the hull, and the
//! stern boundary layer costs a viscous pressure defect), while the roughness
//! allowance is added **outside** it, because surface finish is a property of
//! the skin, not of the form. Keeping them separate matters: they are
//! independently sourced (one from geometry, one from the paint) and on a small
//! craft the roughness term can be the larger of the two, so folding it into
//! `k` hides it.

use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::hull::Hull;

/// How the hull's surface finish is charged on top of flat-plate friction.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Roughness {
    /// Hydraulically smooth: no allowance.
    #[default]
    None,
    /// A prescribed `ΔC_F`, added to `C_F` outside `(1 + k)`.
    ///
    /// This is the unambiguous knob. ITTC-78's own correlation allowance
    /// `C_A` enters here; so does anything else you want charged per unit
    /// wetted area at `½ρU²`.
    DeltaCf(f64),
    /// Equivalent sand-grain roughness height [m]; `ΔC_F` is then estimated
    /// by [`roughness_delta_cf`], which is speed-dependent.
    ///
    /// Rough guide for a hull: a well-sprayed topcoat is ~30 µm, a rolled
    /// antifouling ~100–150 µm, light slime a few hundred µm, and anything
    /// calcareous is millimetres. These are equivalent sand-grain heights,
    /// not paint film thicknesses or profilometer readings.
    SandGrain(f64),
}

/// Viscous knobs, grouped so adding one does not churn every signature.
#[derive(Debug, Clone, Copy, Default)]
pub struct ViscousOptions {
    /// Form factor `k`, applied as `(1 + k)·C_F`. `0` is the bare flat plate.
    pub form_factor: f64,
    /// Surface-roughness allowance, added outside `(1 + k)`.
    pub roughness: Roughness,
}

impl ViscousOptions {
    /// Just a form factor, no roughness allowance.
    pub fn form_factor(k: f64) -> Self {
        ViscousOptions {
            form_factor: k,
            roughness: Roughness::None,
        }
    }

    fn validate(&self) -> Result<()> {
        if !self.form_factor.is_finite() || self.form_factor < 0.0 {
            return Err(Error::InvalidConditions(format!(
                "form factor must be finite and non-negative, got {}",
                self.form_factor
            )));
        }
        match self.roughness {
            Roughness::None => Ok(()),
            Roughness::DeltaCf(c) if c.is_finite() && c >= 0.0 => Ok(()),
            Roughness::DeltaCf(c) => Err(Error::InvalidConditions(format!(
                "roughness ΔC_F must be finite and non-negative, got {c}"
            ))),
            Roughness::SandGrain(k) if k.is_finite() && k >= 0.0 => Ok(()),
            Roughness::SandGrain(k) => Err(Error::InvalidConditions(format!(
                "sand-grain roughness height must be finite and non-negative, got {k} m"
            ))),
        }
    }
}

/// Viscous resistance result.
#[derive(Debug, Clone, Copy)]
pub struct ViscousResistance {
    /// Viscous resistance `R_v = [(1 + k)·C_F + ΔC_F] · ½ ρ U² S` [N].
    pub resistance: f64,
    /// ITTC-57 flat-plate friction coefficient C_F.
    pub cf: f64,
    /// Reynolds number U·L/ν based on the hull length.
    pub reynolds: f64,
    /// Wetted surface S used as reference area [m²].
    pub wetted_surface: f64,
    /// Form factor k applied as (1 + k).
    pub form_factor: f64,
    /// Roughness allowance `ΔC_F` actually applied (0 when smooth).
    pub roughness_cf: f64,
    /// Roughness Reynolds number `k_s⁺ = k_s·u_τ/ν`, when the allowance came
    /// from a sand-grain height. Below ~5 the surface is hydraulically
    /// smooth; above ~70 it is fully rough; between, transitional — and that
    /// is the band where [`roughness_delta_cf`] is least trustworthy.
    pub roughness_reynolds: Option<f64>,
}

impl ViscousResistance {
    /// The roughness share of the viscous coefficient,
    /// `ΔC_F / [(1+k)·C_F + ΔC_F]`.
    pub fn roughness_fraction(&self) -> f64 {
        let cv = (1.0 + self.form_factor) * self.cf + self.roughness_cf;
        if cv > 0.0 {
            self.roughness_cf / cv
        } else {
            0.0
        }
    }
}

/// ITTC-57 correlation line: `C_F = 0.075 / (log10(Re) - 2)²`.
pub fn ittc57_cf(reynolds: f64) -> f64 {
    let d = reynolds.log10() - 2.0;
    0.075 / (d * d)
}

/// Prandtl–Schlichting **fully-rough** flat-plate friction,
/// `C_F = (1.89 + 1.62·log10(L/k_s))^(−2.5)`.
///
/// Reynolds-independent, as fully-rough flow is. Schlichting fits it for
/// `10² < L/k_s < 10⁶`; outside that the argument is clamped, since the
/// formula is an interpolation of plate data, not a law.
pub fn schlichting_rough_cf(length_over_ks: f64) -> f64 {
    let r = length_over_ks.clamp(1e2, 1e6);
    (1.89 + 1.62 * r.log10()).powf(-2.5)
}

/// Roughness Reynolds number `k_s⁺ = k_s·u_τ/ν` with `u_τ = U·√(C_F/2)`.
///
/// The regime indicator: `≲ 5` hydraulically smooth (the roughness sits
/// inside the viscous sublayer and costs nothing), `≳ 70` fully rough,
/// transitional between.
pub fn roughness_reynolds(k_s: f64, speed: f64, cf: f64, kinematic_viscosity: f64) -> f64 {
    k_s * speed * (cf / 2.0).sqrt() / kinematic_viscosity
}

/// `ΔC_F` for an equivalent sand-grain height `k_s` on a plate of length `L`
/// at Reynolds number `Re`: the excess of the fully-rough friction over the
/// smooth line, floored at zero.
///
/// This is the textbook Moody construction — take whichever of the smooth and
/// fully-rough curves is higher — and it has the two properties that matter:
/// it returns **exactly zero** while the surface is hydraulically smooth at
/// this Reynolds number, and it grows with speed as the smooth line falls
/// away beneath the (Re-independent) rough one.
///
/// **What it is not.** The transitional band is *bridged by the crossover*,
/// not fitted to Nikuradse or Colebrook data, so between roughly
/// `k_s⁺ = 5` and `70` it reads high — check
/// [`ViscousResistance::roughness_reynolds`] to see whether you are in it. Nor
/// is it a ship *correlation* allowance: the published ones (Bowden–Davison,
/// Townsin) are regressed at ship scale and fall apart on small craft, where
/// the `(k_s/L)^{1/3}` term runs away as `L` shrinks. On an 8 m hull at
/// `Re ≈ 9×10⁶` Bowden–Davison returns `1.8×10⁻³` (about 60% of `C_F`, which
/// is not credible for paint) while Townsin returns a negative number; that
/// spread, not either value, is the honest state of the art at this size. If
/// you have a number you trust, pass [`Roughness::DeltaCf`] instead.
pub fn roughness_delta_cf(k_s: f64, length: f64, reynolds: f64) -> f64 {
    if !(k_s > 0.0 && length > 0.0 && reynolds > 1e3) {
        return 0.0;
    }
    (schlichting_rough_cf(length / k_s) - ittc57_cf(reynolds)).max(0.0)
}

/// Viscous resistance with zero form factor and a smooth hull (bare ITTC-57).
pub fn viscous_resistance(hull: &Hull, cond: &Conditions) -> Result<ViscousResistance> {
    viscous_resistance_with_options(hull, cond, &ViscousOptions::default())
}

/// Viscous resistance with a form factor k, i.e. `R_v = (1 + k) R_F`, and no
/// roughness allowance.
pub fn viscous_resistance_with(
    hull: &Hull,
    cond: &Conditions,
    form_factor: f64,
) -> Result<ViscousResistance> {
    viscous_resistance_with_options(hull, cond, &ViscousOptions::form_factor(form_factor))
}

/// Viscous resistance with an explicit form factor and roughness allowance.
pub fn viscous_resistance_with_options(
    hull: &Hull,
    cond: &Conditions,
    opts: &ViscousOptions,
) -> Result<ViscousResistance> {
    cond.validate()?;
    opts.validate()?;
    let re = cond.speed * hull.length() / cond.fluid.kinematic_viscosity;
    if re <= 1e3 {
        return Err(Error::InvalidConditions(format!(
            "Reynolds number {re:.3e} is outside the ITTC-57 line's sensible range"
        )));
    }
    let cf = ittc57_cf(re);
    let (roughness_cf, roughness_reynolds) = match opts.roughness {
        Roughness::None => (0.0, None),
        Roughness::DeltaCf(c) => (c, None),
        Roughness::SandGrain(k_s) => (
            roughness_delta_cf(k_s, hull.length(), re),
            Some(self::roughness_reynolds(
                k_s,
                cond.speed,
                cf,
                cond.fluid.kinematic_viscosity,
            )),
        ),
    };
    let s = hull.wetted_surface();
    let cv = (1.0 + opts.form_factor) * cf + roughness_cf;
    let resistance = cv * 0.5 * cond.fluid.density * cond.speed * cond.speed * s;
    Ok(ViscousResistance {
        resistance,
        cf,
        reynolds: re,
        wetted_surface: s,
        form_factor: opts.form_factor,
        roughness_cf,
        roughness_reynolds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hulls;

    #[test]
    fn ittc_line_known_value() {
        // At Re = 1e8: C_F = 0.075 / (8 - 2)^2 = 0.075 / 36.
        assert!((ittc57_cf(1e8) - 0.075 / 36.0).abs() < 1e-15);
    }

    #[test]
    fn smooth_surface_costs_nothing() {
        // A finish fine enough to sit inside the viscous sublayer must be
        // free, not merely cheap — the crossover has to land on zero exactly.
        let l = 8.234;
        let re = 9.3e6;
        assert_eq!(roughness_delta_cf(1e-6, l, re), 0.0);
        assert_eq!(roughness_delta_cf(0.0, l, re), 0.0);
        // ... and a rough one must not be.
        assert!(roughness_delta_cf(1e-3, l, re) > 0.0);
    }

    #[test]
    fn roughness_grows_with_height_and_with_speed() {
        let l = 8.234;
        let re = 9.3e6;
        let mut prev = 0.0;
        for k_s in [1e-4, 3e-4, 1e-3, 3e-3] {
            let d = roughness_delta_cf(k_s, l, re);
            assert!(d > prev, "k_s {k_s}: {d} not above {prev}");
            prev = d;
        }
        // The rough branch is Re-independent while the smooth line falls, so
        // the penalty grows with speed.
        let slow = roughness_delta_cf(1e-3, l, 9.3e6);
        let fast = roughness_delta_cf(1e-3, l, 2.8e7);
        assert!(fast > slow, "{fast} not above {slow}");
    }

    #[test]
    fn roughness_is_added_outside_the_form_factor() {
        // C_V = (1+k)·C_F + ΔC_F, not (1+k)(C_F + ΔC_F): the ITTC-78 shape.
        let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let cond = Conditions::seawater(3.0);
        let dcf = 4e-4;
        let opts = ViscousOptions {
            form_factor: 0.2,
            roughness: Roughness::DeltaCf(dcf),
        };
        let v = viscous_resistance_with_options(&hull, &cond, &opts).unwrap();
        let q = 0.5 * cond.fluid.density * cond.speed * cond.speed * hull.wetted_surface();
        let want = (1.2 * v.cf + dcf) * q;
        assert!((v.resistance - want).abs() < 1e-9 * want);
        assert_eq!(v.roughness_cf, dcf);
        assert!((v.roughness_fraction() - dcf / (1.2 * v.cf + dcf)).abs() < 1e-12);
        // And the two knobs are genuinely separate.
        let k_only = viscous_resistance_with(&hull, &cond, 0.2).unwrap();
        assert_eq!(k_only.roughness_cf, 0.0);
        assert!((v.resistance - k_only.resistance - dcf * q).abs() < 1e-9 * want);
    }

    #[test]
    fn defaults_reproduce_the_bare_flat_plate() {
        let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let cond = Conditions::seawater(3.0);
        let a = viscous_resistance(&hull, &cond).unwrap();
        let b = viscous_resistance_with_options(&hull, &cond, &ViscousOptions::default()).unwrap();
        assert_eq!(a.resistance, b.resistance);
        assert_eq!(a.roughness_cf, 0.0);
        assert!(a.roughness_reynolds.is_none());
        let q = 0.5 * cond.fluid.density * cond.speed * cond.speed * hull.wetted_surface();
        assert!((a.resistance - a.cf * q).abs() < 1e-12 * a.resistance);
    }

    #[test]
    fn sand_grain_reports_its_regime() {
        let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let cond = Conditions::seawater(3.0);
        let opts = ViscousOptions {
            form_factor: 0.0,
            roughness: Roughness::SandGrain(2e-4),
        };
        let v = viscous_resistance_with_options(&hull, &cond, &opts).unwrap();
        let ksp = v.roughness_reynolds.expect("k_s+ reported");
        // k_s+ = k_s·U·sqrt(C_F/2)/ν, a few units for 0.2 mm at 3 m/s.
        let want = 2e-4 * 3.0 * (v.cf / 2.0).sqrt() / cond.fluid.kinematic_viscosity;
        assert!((ksp - want).abs() < 1e-9 * want);
        assert!(ksp > 1.0 && ksp < 100.0, "k_s+ {ksp} implausible");
    }

    #[test]
    fn rejects_negative_allowances() {
        let hull = hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let cond = Conditions::seawater(3.0);
        for r in [Roughness::DeltaCf(-1e-4), Roughness::SandGrain(-1e-4)] {
            let opts = ViscousOptions {
                form_factor: 0.0,
                roughness: r,
            };
            assert!(viscous_resistance_with_options(&hull, &cond, &opts).is_err());
        }
        assert!(viscous_resistance_with(&hull, &cond, -0.1).is_err());
    }
}

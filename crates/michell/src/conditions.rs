//! Fluid properties and operating conditions.

use crate::error::{Error, Result};

/// Fluid properties. ITTC recommended values at 15 °C are provided as
/// constants; supply your own for other temperatures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fluid {
    /// Density ρ [kg/m³].
    pub density: f64,
    /// Kinematic viscosity ν [m²/s].
    pub kinematic_viscosity: f64,
}

impl Fluid {
    /// Salt water at 15 °C (ITTC): ρ = 1025.9 kg/m³, ν = 1.1892e-6 m²/s.
    pub const SEAWATER_15C: Fluid = Fluid {
        density: 1025.9,
        kinematic_viscosity: 1.1892e-6,
    };

    /// Fresh water at 15 °C (ITTC): ρ = 999.1 kg/m³, ν = 1.1386e-6 m²/s.
    pub const FRESHWATER_15C: Fluid = Fluid {
        density: 999.1,
        kinematic_viscosity: 1.1386e-6,
    };
}

/// Steady-ahead operating conditions in calm, deep water.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conditions {
    /// Ship speed U [m/s].
    pub speed: f64,
    pub fluid: Fluid,
    /// Gravitational acceleration g [m/s²].
    pub gravity: f64,
}

pub const STANDARD_GRAVITY: f64 = 9.80665;

impl Conditions {
    pub fn seawater(speed: f64) -> Conditions {
        Conditions {
            speed,
            fluid: Fluid::SEAWATER_15C,
            gravity: STANDARD_GRAVITY,
        }
    }

    pub fn freshwater(speed: f64) -> Conditions {
        Conditions {
            speed,
            fluid: Fluid::FRESHWATER_15C,
            gravity: STANDARD_GRAVITY,
        }
    }

    /// Length Froude number U/√(gL).
    pub fn froude_number(&self, length: f64) -> f64 {
        self.speed / (self.gravity * length).sqrt()
    }

    pub(crate) fn validate(&self) -> Result<()> {
        let ok = |v: f64| v.is_finite() && v > 0.0;
        if !ok(self.speed) {
            return Err(Error::InvalidConditions(format!(
                "speed must be finite and positive, got {}",
                self.speed
            )));
        }
        if !ok(self.fluid.density) || !ok(self.fluid.kinematic_viscosity) {
            return Err(Error::InvalidConditions(
                "fluid density and kinematic viscosity must be finite and positive".into(),
            ));
        }
        if !ok(self.gravity) {
            return Err(Error::InvalidConditions(
                "gravity must be finite and positive".into(),
            ));
        }
        Ok(())
    }
}

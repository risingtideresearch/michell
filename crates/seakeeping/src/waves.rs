//! Deep-water regular-wave kinematics.

/// Deep-water wavenumber `k = ω²/g` [rad/m] of a wave of absolute
/// frequency `omega` [rad/s].
pub fn wavenumber(omega: f64, gravity: f64) -> f64 {
    omega * omega / gravity
}

/// Encounter frequency `ω_e = ω − k U cos β` [rad/s] of a deep-water wave of
/// absolute frequency `omega` met at `speed` and `heading` (`β = π` head
/// seas, `0` following). Negative in following seas when the ship overtakes
/// the wave (`U cos β` beyond the phase speed): the wave is then met from
/// astern, and its sign says so.
pub fn encounter_frequency(omega: f64, speed: f64, heading: f64, gravity: f64) -> f64 {
    omega - wavenumber(omega, gravity) * speed * heading.cos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn head_seas_raise_and_following_seas_lower_the_encounter_frequency() {
        let (g, w, u) = (9.81, 1.2, 5.0);
        let k = wavenumber(w, g);
        assert!((encounter_frequency(w, u, PI, g) - (w + k * u)).abs() < 1e-12);
        assert!((encounter_frequency(w, u, 0.0, g) - (w - k * u)).abs() < 1e-12);
        assert!((encounter_frequency(w, 0.0, 0.7, g) - w).abs() < 1e-12);
        // Phase speed ω/k = g/ω.
        assert!((w / k - g / w).abs() < 1e-12);
    }
}

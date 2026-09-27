//! **Irregular seas**: wave spectra and the statistics of the response.
//!
//! A long-crested sea of spectrum `S(ω)` met at one heading drives each
//! response `r` with spectral moments `m_n = ∫ ω_e^n |R(ω)|² S(ω) dω`
//! (`R` the response per unit wave amplitude). The integral is taken over
//! the **absolute** wave frequency, so no encounter-spectrum Jacobian is
//! needed; for a narrow-banded response the significant (mean of the
//! highest third) single amplitude is `2√m₀`. The mean added resistance is
//! `2∫ (R_aw/ζ_a²) S(ω) dω`.
//!
//! Following and stern-quartering seas whose components are overtaken
//! (encounter frequency ≤ 0) are outside the strip solver; those components
//! are skipped and counted in [`SeaResponse::skipped_energy`].

use crate::strip::{added_resistance_fleet, MassProperties, StripOptions, Wave};
use michell_geometry::{Placement, Result, SectionalHull};
use std::f64::consts::PI;

/// A one-dimensional wave spectrum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Spectrum {
    /// Bretschneider / ITTC two-parameter spectrum (fully developed sea
    /// shape) of significant height `hs` [m] and peak period `tp` [s].
    Bretschneider { hs: f64, tp: f64 },
    /// JONSWAP (fetch-limited) of significant height `hs`, peak period
    /// `tp` and peak enhancement `gamma` (3.3 is the mean North Sea value),
    /// normalised to `hs` by the usual `1 − 0.287 ln γ`.
    Jonswap { hs: f64, tp: f64, gamma: f64 },
}

impl Spectrum {
    /// Peak frequency `ω_p` [rad/s].
    pub fn peak_frequency(&self) -> f64 {
        let tp = match *self {
            Spectrum::Bretschneider { tp, .. } | Spectrum::Jonswap { tp, .. } => tp,
        };
        2.0 * PI / tp
    }

    /// Spectral density `S(ω)` [m²·s/rad].
    pub fn density(&self, omega: f64) -> f64 {
        if !(omega > 0.0) {
            return 0.0;
        }
        let wp = self.peak_frequency();
        let pm = |hs: f64| {
            let r = wp / omega;
            5.0 / 16.0 * hs * hs * wp.powi(4) / omega.powi(5) * (-1.25 * r.powi(4)).exp()
        };
        match *self {
            Spectrum::Bretschneider { hs, .. } => pm(hs),
            Spectrum::Jonswap { hs, gamma, .. } => {
                let sigma = if omega <= wp { 0.07 } else { 0.09 };
                let peak = (-(omega - wp).powi(2) / (2.0 * sigma * sigma * wp * wp)).exp();
                (1.0 - 0.287 * gamma.ln()) * pm(hs) * gamma.powf(peak)
            }
        }
    }
}

/// Statistics of the heave–pitch response and added resistance in a sea.
#[derive(Debug, Clone, PartialEq)]
pub struct SeaResponse {
    /// Significant heave single amplitude `2√m₀` [m].
    pub heave: f64,
    /// Significant pitch single amplitude [rad].
    pub pitch: f64,
    /// Significant vertical acceleration single amplitude at each requested
    /// station [m/s²].
    pub accelerations: Vec<f64>,
    /// Mean added resistance [N].
    pub added_resistance: f64,
    /// Share of the sea's energy `m₀` the solver could not take (overtaken
    /// components in following seas).
    pub skipped_energy: f64,
}

/// The response of `hull` at `speed` and `heading` to the sea `spectrum`,
/// sampled at `samples` frequencies from 0.3 to 4 times the peak frequency.
/// `stations` are hull-x positions for vertical accelerations (e.g. bow and
/// helm).
#[allow(clippy::too_many_arguments)]
pub fn sea_response(
    hull: &SectionalHull,
    mass: &MassProperties,
    spectrum: &Spectrum,
    heading: f64,
    speed: f64,
    stations: &[f64],
    samples: usize,
    opts: &StripOptions,
) -> Result<SeaResponse> {
    sea_response_fleet(
        &[(hull, Placement::default())],
        mass,
        spectrum,
        heading,
        speed,
        stations,
        samples,
        opts,
    )
}

/// [`sea_response`] for a platform of several placed hulls; `stations` are
/// platform-x positions.
#[allow(clippy::too_many_arguments)]
pub fn sea_response_fleet(
    members: &[(&SectionalHull, Placement)],
    mass: &MassProperties,
    spectrum: &Spectrum,
    heading: f64,
    speed: f64,
    stations: &[f64],
    samples: usize,
    opts: &StripOptions,
) -> Result<SeaResponse> {
    let wp = spectrum.peak_frequency();
    let (lo, hi) = (0.3 * wp, 4.0 * wp);
    let n = samples.max(3) | 1; // odd, for Simpson's rule
    let h = (hi - lo) / (n - 1) as f64;
    let (mut m_heave, mut m_pitch, mut raw, mut total, mut skipped) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut m_acc = vec![0.0; stations.len()];
    for i in 0..n {
        let omega = lo + h * i as f64;
        let w = h / 3.0
            * if i == 0 || i == n - 1 {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
        let s = spectrum.density(omega);
        total += w * s;
        let wave = Wave {
            omega,
            heading,
            speed,
        };
        let r = match added_resistance_fleet(members, mass, &wave, opts) {
            Ok(r) => r,
            Err(_) => {
                skipped += w * s;
                continue;
            }
        };
        let resp = r.response;
        m_heave += w * s * resp.heave.abs_sq();
        m_pitch += w * s * resp.pitch.abs_sq();
        for (m, &x) in m_acc.iter_mut().zip(stations) {
            *m += w * s * resp.omega_e.powi(4) * resp.vertical_motion(x, mass.lcg).abs_sq();
        }
        raw += w * s * 2.0 * r.per_amplitude_sq;
    }
    Ok(SeaResponse {
        heave: 2.0 * m_heave.sqrt(),
        pitch: 2.0 * m_pitch.sqrt(),
        accelerations: m_acc.iter().map(|m| 2.0 * m.sqrt()).collect(),
        added_resistance: raw,
        skipped_energy: if total > 0.0 { skipped / total } else { 0.0 },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m0(s: &Spectrum) -> f64 {
        let n = 20000;
        let (lo, hi) = (0.01, 20.0 * s.peak_frequency());
        let h = (hi - lo) / n as f64;
        (0..n)
            .map(|i| s.density(lo + h * (i as f64 + 0.5)) * h)
            .sum()
    }

    #[test]
    fn spectra_carry_their_significant_height() {
        for s in [
            Spectrum::Bretschneider { hs: 2.0, tp: 8.0 },
            Spectrum::Jonswap {
                hs: 2.0,
                tp: 8.0,
                gamma: 3.3,
            },
            Spectrum::Jonswap {
                hs: 1.0,
                tp: 5.0,
                gamma: 1.0,
            },
        ] {
            let hs = 4.0 * m0(&s).sqrt();
            let want = match s {
                Spectrum::Bretschneider { hs, .. } | Spectrum::Jonswap { hs, .. } => hs,
            };
            // The JONSWAP normalisation is an approximation (~1–2%).
            assert!((hs - want).abs() < 0.02 * want, "{s:?}: Hs {hs}");
        }
    }

    fn wigley(l: f64) -> SectionalHull {
        use michell_geometry::iges::{self, HullPose, Platform, SectionalOptions};
        let surfaces = iges::wigley_surfaces(l, 0.1 * l, 0.0625 * l).unwrap();
        let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
        let so = SectionalOptions {
            stations: 31,
            ..Default::default()
        };
        source
            .situate_sectional(0, 0.0, &HullPose::default(), &Platform::default(), &so)
            .unwrap()
            .unwrap()
            .hull
    }

    /// A swell far longer than the hull lifts it with the water: its
    /// significant heave is the sea's significant amplitude `Hs/2`.
    #[test]
    fn a_long_swell_is_ridden() {
        let l = 3.0;
        let hull = wigley(l);
        let opts = StripOptions {
            panels: 12,
            density: 1000.0,
            gravity: 9.81,
        };
        let mass = MassProperties::floating(&hull, 1000.0, 0.25 * l);
        let sea = Spectrum::Bretschneider { hs: 0.2, tp: 12.0 };
        let r = sea_response(&hull, &mass, &sea, PI, 0.0, &[0.0], 21, &opts).unwrap();
        // Simpson over 0.3–4 ω_p holds ~97% of m₀.
        assert!((r.heave - 0.1).abs() < 0.05 * 0.1, "heave {}", r.heave);
        assert!(r.added_resistance >= 0.0 && r.skipped_energy == 0.0);
    }

    #[test]
    fn jonswap_with_unit_gamma_is_bretschneider() {
        let b = Spectrum::Bretschneider { hs: 1.5, tp: 7.0 };
        let j = Spectrum::Jonswap {
            hs: 1.5,
            tp: 7.0,
            gamma: 1.0,
        };
        for w in [0.5, 0.9, 1.3, 2.0] {
            assert!((b.density(w) - j.density(w)).abs() < 1e-12);
        }
    }
}

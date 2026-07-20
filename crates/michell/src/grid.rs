//! The sampling intermediate representation: a station × waterline grid of
//! half-beam samples, optionally augmented with first derivatives and
//! per-sample weights.
//!
//! Every input front-end reduces to a `SampleGrid` — a hand-typed offset
//! table carries values only, the IGES importer adds `∂f/∂x` and `∂f/∂z`
//! recovered from the CAD surface, and a re-situated [`crate::body::Body`]
//! adds exact derivatives from its source spline. [`crate::fit::fit_grid`]
//! lofts a grid to the canonical B-spline hull, using whichever channels are
//! present.
//!
//! Channel conventions:
//! - the value channel is required, finite, and non-negative (small negative
//!   samples are treated as measurement noise and clamped to zero);
//! - derivative channels are optional; a `NaN` entry means "unknown at this
//!   sample" and is simply not used as an observation;
//! - the weight channel is optional; weight `0` excludes a sample entirely
//!   (e.g. a failed CAD inversion — *unknown* geometry, as opposed to a dry
//!   sample where a half-beam of `0` is correct data).

use crate::error::{Error, Result};

/// A gridded set of half-beam samples with optional derivative and weight
/// channels. See the module docs for channel conventions.
#[derive(Debug, Clone)]
pub struct SampleGrid {
    stations: Vec<f64>,
    waterlines: Vec<f64>,
    /// Row-major, waterline index fastest: `f[i * waterlines.len() + j]`.
    f: Vec<f64>,
    fx: Option<Vec<f64>>,
    fz: Option<Vec<f64>>,
    weight: Option<Vec<f64>>,
}

impl SampleGrid {
    /// Build a value-only grid. `stations` (strictly increasing x, metres)
    /// and `waterlines` (strictly increasing z downward from the waterline,
    /// starting at 0) define the grid; `half_beams[i * waterlines.len() + j]`
    /// is the half-beam at `(stations[i], waterlines[j])`. Small negative
    /// samples (measurement noise) are clamped to zero; clearly negative
    /// values are rejected.
    pub fn new(stations: Vec<f64>, waterlines: Vec<f64>, half_beams: Vec<f64>) -> Result<Self> {
        let mx = stations.len();
        let mz = waterlines.len();
        if mx < 2 || mz < 2 {
            return Err(Error::InvalidInput(format!(
                "a sample grid needs at least 2 stations and 2 waterlines; \
                 got {mx} x {mz}"
            )));
        }
        if half_beams.len() != mx * mz {
            return Err(Error::InvalidInput(format!(
                "half_beams has {} entries, expected {} stations x {} waterlines = {}",
                half_beams.len(),
                mx,
                mz,
                mx * mz
            )));
        }
        check_strictly_increasing(&stations, "stations")?;
        check_strictly_increasing(&waterlines, "waterlines")?;
        let z_span = waterlines[mz - 1] - waterlines[0];
        if waterlines[0] < 0.0 || waterlines[0] > 1e-9 * z_span {
            return Err(Error::InvalidInput(format!(
                "waterlines must start at the design waterline z = 0 (z downward); got {}",
                waterlines[0]
            )));
        }
        if half_beams.iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput(
                "half_beams contains a non-finite value".into(),
            ));
        }
        let y_scale = half_beams.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        if half_beams.iter().any(|&v| v < -1e-6 * y_scale.max(1.0)) {
            return Err(Error::InvalidInput(
                "half_beams contains clearly negative values; half-beams must be >= 0".into(),
            ));
        }
        Ok(SampleGrid {
            stations,
            waterlines,
            f: half_beams.into_iter().map(|v| v.max(0.0)).collect(),
            fx: None,
            fz: None,
            weight: None,
        })
    }

    /// Attach a `∂f/∂x` channel (same layout as the values; `NaN` = unknown
    /// at that sample).
    pub fn with_fx(mut self, fx: Vec<f64>) -> Result<Self> {
        self.check_deriv_channel(&fx, "fx")?;
        self.fx = Some(fx);
        Ok(self)
    }

    /// Attach a `∂f/∂z` channel (same layout as the values; `NaN` = unknown
    /// at that sample).
    pub fn with_fz(mut self, fz: Vec<f64>) -> Result<Self> {
        self.check_deriv_channel(&fz, "fz")?;
        self.fz = Some(fz);
        Ok(self)
    }

    /// Attach per-sample weights (same layout; finite, `>= 0`; weight `0`
    /// excludes the sample from every channel).
    pub fn with_weights(mut self, weight: Vec<f64>) -> Result<Self> {
        if weight.len() != self.f.len() {
            return Err(Error::InvalidInput(format!(
                "weight channel has {} entries, expected {}",
                weight.len(),
                self.f.len()
            )));
        }
        if weight.iter().any(|w| !(w.is_finite() && *w >= 0.0)) {
            return Err(Error::InvalidInput(
                "weights must be finite and non-negative".into(),
            ));
        }
        self.weight = Some(weight);
        Ok(self)
    }

    fn check_deriv_channel(&self, chan: &[f64], name: &str) -> Result<()> {
        if chan.len() != self.f.len() {
            return Err(Error::InvalidInput(format!(
                "{name} channel has {} entries, expected {}",
                chan.len(),
                self.f.len()
            )));
        }
        // NaN is the "unknown here" marker; infinities are garbage.
        if chan.iter().any(|v| v.is_infinite()) {
            return Err(Error::InvalidInput(format!(
                "{name} channel contains an infinite value"
            )));
        }
        Ok(())
    }

    pub fn stations(&self) -> &[f64] {
        &self.stations
    }

    pub fn waterlines(&self) -> &[f64] {
        &self.waterlines
    }

    pub fn half_beams(&self) -> &[f64] {
        &self.f
    }

    pub fn fx(&self) -> Option<&[f64]> {
        self.fx.as_deref()
    }

    pub fn fz(&self) -> Option<&[f64]> {
        self.fz.as_deref()
    }

    pub fn weights(&self) -> Option<&[f64]> {
        self.weight.as_deref()
    }

    /// Flat index of sample `(station i, waterline j)`.
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        i * self.waterlines.len() + j
    }
}

fn check_strictly_increasing(t: &[f64], name: &str) -> Result<()> {
    if t.iter().any(|v| !v.is_finite()) {
        return Err(Error::InvalidInput(format!(
            "{name} contains a non-finite value"
        )));
    }
    if t.windows(2).any(|w| w[1] <= w[0]) {
        return Err(Error::InvalidInput(format!(
            "{name} must be strictly increasing"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid_2x2() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (vec![0.0, 1.0], vec![0.0, 0.5], vec![0.1, 0.2, 0.3, 0.4])
    }

    #[test]
    fn accepts_valid_channels() {
        let (st, wl, f) = grid_2x2();
        let g = SampleGrid::new(st, wl, f)
            .unwrap()
            .with_fx(vec![0.0, f64::NAN, 1.0, -1.0])
            .unwrap()
            .with_fz(vec![f64::NAN; 4])
            .unwrap()
            .with_weights(vec![1.0, 0.0, 2.0, 1.0])
            .unwrap();
        assert!(g.fx().unwrap()[1].is_nan());
        assert_eq!(g.weights().unwrap()[1], 0.0);
        assert_eq!(g.idx(1, 1), 3);
    }

    #[test]
    fn clamps_noise_rejects_negative() {
        let (st, wl, mut f) = grid_2x2();
        f[0] = -1e-12;
        let g = SampleGrid::new(st.clone(), wl.clone(), f.clone()).unwrap();
        assert_eq!(g.half_beams()[0], 0.0);
        f[0] = -0.5;
        assert!(SampleGrid::new(st, wl, f).is_err());
    }

    #[test]
    fn rejects_bad_input() {
        let (st, wl, f) = grid_2x2();
        // Wrong lengths.
        assert!(SampleGrid::new(st.clone(), wl.clone(), vec![0.0; 3]).is_err());
        assert!(SampleGrid::new(st.clone(), wl.clone(), f.clone())
            .unwrap()
            .with_fx(vec![0.0; 3])
            .is_err());
        assert!(SampleGrid::new(st.clone(), wl.clone(), f.clone())
            .unwrap()
            .with_weights(vec![1.0; 3])
            .is_err());
        // Non-increasing stations.
        assert!(SampleGrid::new(vec![0.0, 0.0], wl.clone(), f.clone()).is_err());
        // Waterlines not starting at 0.
        assert!(SampleGrid::new(st.clone(), vec![0.1, 0.5], f.clone()).is_err());
        // NaN value; infinite derivative; negative weight.
        assert!(SampleGrid::new(st.clone(), wl.clone(), vec![0.0, f64::NAN, 0.0, 0.0]).is_err());
        assert!(SampleGrid::new(st.clone(), wl.clone(), f.clone())
            .unwrap()
            .with_fx(vec![0.0, f64::INFINITY, 0.0, 0.0])
            .is_err());
        assert!(SampleGrid::new(st, wl, f)
            .unwrap()
            .with_weights(vec![1.0, -1.0, 1.0, 1.0])
            .is_err());
    }
}

//! Least-squares lofting of gridded half-breadth samples to the canonical
//! B-spline hull representation.
//!
//! This is the shared back-end for every input front-end: a hand-typed
//! station × waterline offset table and a CAD import both reduce to a
//! [`SampleGrid`], which is fit here by weighted tensor-product least
//! squares. When the grid carries derivative channels (`∂f/∂x`, `∂f/∂z`),
//! they join the fit as additional observations — pinning slopes as well as
//! values, which matters most where the sampling is coarse relative to the
//! shape (stem, sterns, re-lofted poses). The fit is deterministic: interior
//! knots are placed at sample quantiles (which keeps the normal equations
//! well-conditioned whenever the sample grid is denser than the control
//! net), and the achieved residuals are reported per channel so callers can
//! judge fit quality.

use crate::bspline::{self, BSplineSurface};
use crate::error::{Error, Result};
use crate::grid::SampleGrid;
use crate::hull::Hull;

/// Fit configuration.
#[derive(Debug, Clone, Copy)]
pub struct FitOptions {
    pub degree_x: usize,
    pub degree_z: usize,
    /// Control points along x (stations direction).
    pub n_ctrl_x: usize,
    /// Control points along z (waterlines direction).
    pub n_ctrl_z: usize,
    /// Relative weight of the derivative channels. Derivative residuals are
    /// scaled by the local sample spacing (making a slope error over one
    /// sample cell commensurate with a value error of the same size), then
    /// by this factor. `0` ignores derivative channels entirely.
    pub derivative_weight: f64,
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            degree_x: 3,
            degree_z: 3,
            n_ctrl_x: 12,
            n_ctrl_z: 8,
            derivative_weight: 1.0,
        }
    }
}

/// Residual summary of one derivative channel, in the channel's units.
#[derive(Debug, Clone, Copy)]
pub struct ChannelResiduals {
    /// Largest |fitted − sample| over the observed samples.
    pub max: f64,
    /// Root-mean-square residual over the observed samples.
    pub rms: f64,
}

/// Fit quality diagnostics, in the units of the input half-beams.
#[derive(Debug, Clone, Copy)]
pub struct FitReport {
    /// Largest |fitted − sample| over the (weight > 0) grid samples.
    pub max_residual: f64,
    /// Sample location (x, z) of the largest residual.
    pub max_residual_at: (f64, f64),
    /// Root-mean-square residual over the (weight > 0) grid samples.
    pub rms_residual: f64,
    /// Most negative control value floored to zero to enforce `f >= 0`
    /// (0.0 if the unconstrained fit was already non-negative).
    pub floored: f64,
    /// Residuals of the `∂f/∂x` channel, when the grid carried one (and at
    /// least one sample of it was observed).
    pub fx_residual: Option<ChannelResiduals>,
    /// Residuals of the `∂f/∂z` channel, likewise.
    pub fz_residual: Option<ChannelResiduals>,
}

/// Fit a hull to gridded half-beam samples (value channel only). Convenience
/// wrapper over [`fit_grid`]; see [`SampleGrid::new`] for the grid contract.
///
/// Requires a few more samples than control points in each direction
/// (`len >= n_ctrl + 2`).
pub fn fit_offsets(
    stations: &[f64],
    waterlines: &[f64],
    half_beams: &[f64],
    opts: &FitOptions,
) -> Result<(Hull, FitReport)> {
    let grid = SampleGrid::new(stations.to_vec(), waterlines.to_vec(), half_beams.to_vec())?;
    fit_grid(&grid, opts)
}

/// Fit a hull to a [`SampleGrid`], using every channel the grid carries.
///
/// All observations — values, and derivatives where present and known — are
/// assembled into one weighted least-squares system over the tensor-product
/// control net and solved by dense Cholesky (the net is small; a few hundred
/// unknowns). With a value-only, unit-weight grid this reproduces the
/// classical separable fit exactly (the normal equations are identical).
pub fn fit_grid(grid: &SampleGrid, opts: &FitOptions) -> Result<(Hull, FitReport)> {
    let st = grid.stations();
    let wl = grid.waterlines();
    let (mx, mz) = (st.len(), wl.len());
    let (px, pz) = (opts.degree_x, opts.degree_z);
    let (nx, nz) = (opts.n_ctrl_x, opts.n_ctrl_z);
    if !(opts.derivative_weight.is_finite() && opts.derivative_weight >= 0.0) {
        return Err(Error::InvalidInput(
            "derivative_weight must be finite and non-negative".into(),
        ));
    }
    let kx = approx_knots(st, px, nx)?;
    let kz = approx_knots(wl, pz, nz)?;

    // Per-sample basis rows: (first control index, [values, derivatives]).
    let rows_x: Vec<(usize, Vec<Vec<f64>>)> = st
        .iter()
        .map(|&u| bspline::basis_rows1(&kx, px, nx, u))
        .collect();
    let rows_z: Vec<(usize, Vec<Vec<f64>>)> = wl
        .iter()
        .map(|&u| bspline::basis_rows1(&kz, pz, nz, u))
        .collect();

    // Derivative observations are commensurated with value observations by
    // the local sample spacing (a slope error over one sample cell ≈ a value
    // error of the same size — so where samples cluster, e.g. cosine-spaced
    // stations at the ends, slope observations claim proportionally less of
    // the fit), times the user multiplier; weights enter the normal
    // equations squared.
    let local = |t: &[f64], i: usize| -> f64 {
        (t[(i + 1).min(t.len() - 1)] - t[i.saturating_sub(1)])
            / (((i + 1).min(t.len() - 1) - i.saturating_sub(1)) as f64)
    };
    let dw2 = opts.derivative_weight * opts.derivative_weight;
    let hx: Vec<f64> = (0..mx).map(|i| local(st, i)).collect();
    let hz: Vec<f64> = (0..mz).map(|j| local(wl, j)).collect();

    let f = grid.half_beams();
    let fx = grid.fx().filter(|_| dw2 > 0.0);
    let fz = grid.fz().filter(|_| dw2 > 0.0);
    let weights = grid.weights();

    // A slope observation predicting more change across one sample cell than
    // the half-beam anywhere on its stencil has left the geometry's
    // resolvable scale (near-vertical shell at a bilge turn, the keel fold);
    // no loft at this sampling can honour it, and fitting it only distorts
    // the values. Skip it — the value observation stays. Comparing against
    // the neighbourhood, not the sample alone, keeps the informative slopes
    // where f itself runs out smoothly (stem, stern, waterline endings):
    // there the adjacent half-beam is ≈ |slope|·h by construction.
    let usable = |slope: f64, h: f64, fmax: f64| slope.is_finite() && slope.abs() * h <= fmax;
    let f_near_x = |s: usize| -> f64 {
        let mut m = f[s];
        if s >= mz {
            m = m.max(f[s - mz]);
        }
        if s + mz < f.len() {
            m = m.max(f[s + mz]);
        }
        m
    };
    let f_near_z = |s: usize, j: usize| -> f64 {
        let mut m = f[s];
        if j > 0 {
            m = m.max(f[s - 1]);
        }
        if j + 1 < mz {
            m = m.max(f[s + 1]);
        }
        m
    };

    let n = nx * nz;
    let mut a = vec![0.0f64; n * n];
    let mut rhs = vec![0.0f64; n];
    let mut scratch: Vec<(usize, f64)> = Vec::with_capacity((px + 1) * (pz + 1));
    for (i, (x0, brx)) in rows_x.iter().enumerate() {
        for (j, (z0, brz)) in rows_z.iter().enumerate() {
            let s = i * mz + j;
            let w = weights.map_or(1.0, |ws| ws[s]);
            if w <= 0.0 {
                continue;
            }
            add_observation(
                &mut a,
                &mut rhs,
                nz,
                (*x0, &brx[0]),
                (*z0, &brz[0]),
                w,
                f[s],
                &mut scratch,
            );
            if let Some(fx) = fx {
                if usable(fx[s], hx[i], f_near_x(s)) {
                    add_observation(
                        &mut a,
                        &mut rhs,
                        nz,
                        (*x0, &brx[1]),
                        (*z0, &brz[0]),
                        w * dw2 * hx[i] * hx[i],
                        fx[s],
                        &mut scratch,
                    );
                }
            }
            if let Some(fz) = fz {
                if usable(fz[s], hz[j], f_near_z(s, j)) {
                    add_observation(
                        &mut a,
                        &mut rhs,
                        nz,
                        (*x0, &brx[0]),
                        (*z0, &brz[1]),
                        w * dw2 * hz[j] * hz[j],
                        fz[s],
                        &mut scratch,
                    );
                }
            }
        }
    }

    let l = cholesky(a, n).ok_or_else(ill_conditioned)?;
    let mut control = rhs;
    chol_solve(&l, n, &mut control);

    // Enforce the hull contract f >= 0: snap numerical dust to zero, floor
    // genuine ringing (reported so the caller can add control points if it
    // is significant).
    let c_scale = control.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    let mut floored = 0.0f64;
    for v in control.iter_mut() {
        if v.abs() < 1e-12 * c_scale {
            *v = 0.0;
        } else if *v < 0.0 {
            floored = floored.min(*v);
            *v = 0.0;
        }
    }

    // Residuals of the final (floored) control net, per channel, over the
    // observed (weight > 0) samples.
    let mut max_res = 0.0f64;
    let mut max_at = (st[0], wl[0]);
    let mut sum_sq = 0.0f64;
    let mut n_val = 0usize;
    let mut acc_fx = (0.0f64, 0.0f64, 0usize); // (max, sum_sq, count)
    let mut acc_fz = (0.0f64, 0.0f64, 0usize);
    for (i, (x0, brx)) in rows_x.iter().enumerate() {
        for (j, (z0, brz)) in rows_z.iter().enumerate() {
            let s = i * mz + j;
            let w = weights.map_or(1.0, |ws| ws[s]);
            if w <= 0.0 {
                continue;
            }
            let r = (tensor_dot((*x0, &brx[0]), (*z0, &brz[0]), &control, nz) - f[s]).abs();
            if r > max_res {
                max_res = r;
                max_at = (st[i], wl[j]);
            }
            sum_sq += r * r;
            n_val += 1;
            if let Some(fx) = fx {
                if usable(fx[s], hx[i], f_near_x(s)) {
                    let r =
                        (tensor_dot((*x0, &brx[1]), (*z0, &brz[0]), &control, nz) - fx[s]).abs();
                    acc_fx = (acc_fx.0.max(r), acc_fx.1 + r * r, acc_fx.2 + 1);
                }
            }
            if let Some(fz) = fz {
                if usable(fz[s], hz[j], f_near_z(s, j)) {
                    let r =
                        (tensor_dot((*x0, &brx[0]), (*z0, &brz[1]), &control, nz) - fz[s]).abs();
                    acc_fz = (acc_fz.0.max(r), acc_fz.1 + r * r, acc_fz.2 + 1);
                }
            }
        }
    }
    let channel = |(max, sum_sq, count): (f64, f64, usize)| {
        (count > 0).then(|| ChannelResiduals {
            max,
            rms: (sum_sq / count as f64).sqrt(),
        })
    };
    let report = FitReport {
        max_residual: max_res,
        max_residual_at: max_at,
        rms_residual: (sum_sq / n_val.max(1) as f64).sqrt(),
        floored,
        fx_residual: channel(acc_fx),
        fz_residual: channel(acc_fz),
    };

    let surface = BSplineSurface::new(px, pz, kx, kz, control)?;
    Ok((Hull::new(surface)?, report))
}

/// Accumulate one observation row (the tensor product of two basis rows)
/// into the normal equations `A += w r rᵀ`, `rhs += w r t`.
#[allow(clippy::too_many_arguments)]
fn add_observation(
    a: &mut [f64],
    rhs: &mut [f64],
    nz: usize,
    bx: (usize, &[f64]),
    bz: (usize, &[f64]),
    w: f64,
    target: f64,
    scratch: &mut Vec<(usize, f64)>,
) {
    let n = rhs.len();
    scratch.clear();
    for (i, &vx) in bx.1.iter().enumerate() {
        for (j, &vz) in bz.1.iter().enumerate() {
            scratch.push(((bx.0 + i) * nz + (bz.0 + j), vx * vz));
        }
    }
    for &(ci, c1) in scratch.iter() {
        rhs[ci] += w * c1 * target;
        let row = &mut a[ci * n..(ci + 1) * n];
        for &(cj, c2) in scratch.iter() {
            row[cj] += w * c1 * c2;
        }
    }
}

/// Evaluate a tensor-product row against the control net.
fn tensor_dot(bx: (usize, &[f64]), bz: (usize, &[f64]), control: &[f64], nz: usize) -> f64 {
    let mut s = 0.0;
    for (i, &vx) in bx.1.iter().enumerate() {
        let row = &control[(bx.0 + i) * nz..];
        let mut inner = 0.0;
        for (j, &vz) in bz.1.iter().enumerate() {
            inner += vz * row[bz.0 + j];
        }
        s += vx * inner;
    }
    s
}

fn ill_conditioned() -> Error {
    Error::InvalidInput(
        "least-squares normal equations are singular or ill-conditioned; \
         increase sample density or reduce the number of control points"
            .into(),
    )
}

/// Clamped knot vector for least-squares approximation: interior knots at
/// strictly increasing sample quantiles, so every knot span contains samples
/// (Schoenberg–Whitney condition in practice).
fn approx_knots(t: &[f64], degree: usize, n_ctrl: usize) -> Result<Vec<f64>> {
    let m = t.len();
    if n_ctrl < degree + 1 {
        return Err(Error::InvalidInput(format!(
            "need at least degree + 1 = {} control points, got {n_ctrl}",
            degree + 1
        )));
    }
    if m < n_ctrl + 2 {
        return Err(Error::InvalidInput(format!(
            "need at least n_ctrl + 2 = {} samples for {n_ctrl} control points, got {m}",
            n_ctrl + 2
        )));
    }
    let n_interior = n_ctrl - degree - 1;
    let mut knots = vec![t[0]; degree + 1];
    let mut prev_idx = 0usize;
    for k in 1..=n_interior {
        let mut idx = (k as f64 * (m - 1) as f64 / (n_ctrl - degree) as f64).round() as usize;
        if idx <= prev_idx {
            idx = prev_idx + 1;
        }
        if idx > m - 2 {
            return Err(ill_conditioned());
        }
        knots.push(t[idx]);
        prev_idx = idx;
    }
    knots.extend(std::iter::repeat_n(t[m - 1], degree + 1));
    Ok(knots)
}

/// In-place Cholesky A = L L'; returns the lower factor, or None if the
/// matrix is not (numerically) positive definite.
fn cholesky(mut a: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    for k in 0..n {
        let mut d = a[k * n + k];
        for j in 0..k {
            d -= a[k * n + j] * a[k * n + j];
        }
        if d <= 0.0 || !d.is_finite() {
            return None;
        }
        let d = d.sqrt();
        a[k * n + k] = d;
        for i in (k + 1)..n {
            let mut s = a[i * n + k];
            for j in 0..k {
                s -= a[i * n + j] * a[k * n + j];
            }
            a[i * n + k] = s / d;
        }
    }
    Some(a)
}

/// Solve L L' x = rhs in place.
fn chol_solve(l: &[f64], n: usize, rhs: &mut [f64]) {
    for i in 0..n {
        let mut s = rhs[i];
        for j in 0..i {
            s -= l[i * n + j] * rhs[j];
        }
        rhs[i] = s / l[i * n + i];
    }
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..n {
            s -= l[j * n + i] * rhs[j];
        }
        rhs[i] = s / l[i * n + i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wigley_value(l: f64, b: f64, t: f64, x: f64, z: f64) -> f64 {
        b / 2.0 * (1.0 - (2.0 * x / l).powi(2)) * (1.0 - (z / t).powi(2))
    }

    fn wigley_fx(l: f64, b: f64, t: f64, x: f64, z: f64) -> f64 {
        b / 2.0 * (-8.0 * x / (l * l)) * (1.0 - (z / t).powi(2))
    }

    fn wigley_fz(l: f64, b: f64, t: f64, x: f64, z: f64) -> f64 {
        b / 2.0 * (1.0 - (2.0 * x / l).powi(2)) * (-2.0 * z / (t * t))
    }

    fn wigley_grid(l: f64, b: f64, t: f64, mx: usize, mz: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let stations: Vec<f64> = (0..mx)
            .map(|i| -l / 2.0 + l * i as f64 / (mx - 1) as f64)
            .collect();
        let waterlines: Vec<f64> = (0..mz).map(|j| t * j as f64 / (mz - 1) as f64).collect();
        let mut y = vec![0.0; mx * mz];
        for (i, &x) in stations.iter().enumerate() {
            for (j, &z) in waterlines.iter().enumerate() {
                y[i * mz + j] = wigley_value(l, b, t, x, z);
            }
        }
        (stations, waterlines, y)
    }

    #[test]
    fn recovers_exactly_representable_surface() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let (st, wl, y) = wigley_grid(l, b, t, 41, 17);
        // Cubic fit of a biquadratic: representable exactly.
        let opts = FitOptions {
            degree_x: 3,
            degree_z: 3,
            n_ctrl_x: 10,
            n_ctrl_z: 7,
            ..FitOptions::default()
        };
        let (hull, report) = fit_offsets(&st, &wl, &y, &opts).unwrap();
        assert!(
            report.max_residual < 1e-10 * b,
            "max residual {}",
            report.max_residual
        );
        assert_eq!(report.floored, 0.0);
        assert!(report.fx_residual.is_none() && report.fz_residual.is_none());
        // The fitted hull must reproduce Wigley resistance.
        let reference = crate::hulls::wigley(l, b, t).unwrap();
        let cond = crate::Conditions::freshwater(3.0);
        let rw_fit = crate::wave_resistance(&hull, &cond).unwrap().resistance;
        let rw_ref = crate::wave_resistance(&reference, &cond)
            .unwrap()
            .resistance;
        assert!(
            (rw_fit - rw_ref).abs() < 1e-6 * rw_ref,
            "Rw {rw_fit} vs {rw_ref}"
        );
        // Geometry too.
        assert!((hull.displaced_volume() - reference.displaced_volume()).abs() < 1e-8);
    }

    #[test]
    fn hermite_recovers_biquadratic_from_coarse_grid() {
        // A grid barely denser than the control net, but with exact
        // derivative channels: the fit must still be exact in every channel.
        let (l, b, t) = (10.0, 1.0, 0.625);
        let (st, wl, y) = wigley_grid(l, b, t, 13, 7);
        let (mut fx, mut fz) = (vec![0.0; 13 * 7], vec![0.0; 13 * 7]);
        for (i, &x) in st.iter().enumerate() {
            for (j, &z) in wl.iter().enumerate() {
                fx[i * 7 + j] = wigley_fx(l, b, t, x, z);
                fz[i * 7 + j] = wigley_fz(l, b, t, x, z);
            }
        }
        let grid = SampleGrid::new(st, wl, y)
            .unwrap()
            .with_fx(fx)
            .unwrap()
            .with_fz(fz)
            .unwrap();
        let opts = FitOptions {
            degree_x: 3,
            degree_z: 3,
            n_ctrl_x: 10,
            n_ctrl_z: 5,
            ..FitOptions::default()
        };
        let (hull, report) = fit_grid(&grid, &opts).unwrap();
        assert!(
            report.max_residual < 1e-10 * b,
            "max residual {}",
            report.max_residual
        );
        let fxr = report.fx_residual.expect("fx channel observed");
        let fzr = report.fz_residual.expect("fz channel observed");
        assert!(fxr.max < 1e-10, "fx residual {}", fxr.max);
        assert!(fzr.max < 1e-10, "fz residual {}", fzr.max);
        let reference = crate::hulls::wigley(l, b, t).unwrap();
        assert!((hull.displaced_volume() - reference.displaced_volume()).abs() < 1e-8);
    }

    #[test]
    fn derivative_channels_improve_a_coarse_fit() {
        // Non-representable hull sampled coarsely: derivative observations
        // must reduce the true (off-grid) error, not just the residuals.
        let (l, t) = (10.0, 0.8);
        let truth = |x: f64, z: f64| -> f64 {
            let g = (std::f64::consts::PI * x / l).sin();
            0.5 * g * g * (1.0 - (z / t).powi(3))
        };
        let truth_fx = |x: f64, z: f64| -> f64 {
            let p = std::f64::consts::PI / l;
            (p * x).sin() * (p * x).cos() * p * (1.0 - (z / t).powi(3))
        };
        let truth_fz = |x: f64, z: f64| -> f64 {
            let g = (std::f64::consts::PI * x / l).sin();
            0.5 * g * g * (-3.0 * z * z / (t * t * t))
        };
        let (mx, mz) = (16, 9);
        let st: Vec<f64> = (0..mx).map(|i| l * i as f64 / (mx - 1) as f64).collect();
        let wl: Vec<f64> = (0..mz).map(|j| t * j as f64 / (mz - 1) as f64).collect();
        let mut y = vec![0.0; mx * mz];
        let mut fx = vec![0.0; mx * mz];
        let mut fz = vec![0.0; mx * mz];
        for (i, &x) in st.iter().enumerate() {
            for (j, &z) in wl.iter().enumerate() {
                y[i * mz + j] = truth(x, z);
                fx[i * mz + j] = truth_fx(x, z);
                fz[i * mz + j] = truth_fz(x, z);
            }
        }
        let grid = SampleGrid::new(st, wl, y)
            .unwrap()
            .with_fx(fx)
            .unwrap()
            .with_fz(fz)
            .unwrap();
        let opts = FitOptions {
            degree_x: 3,
            degree_z: 3,
            n_ctrl_x: 12,
            n_ctrl_z: 6,
            ..FitOptions::default()
        };
        let value_only = FitOptions {
            derivative_weight: 0.0,
            ..opts
        };
        let (h1, _) = fit_grid(&grid, &value_only).unwrap();
        let (h2, r2) = fit_grid(&grid, &opts).unwrap();
        assert!(r2.fx_residual.is_some());
        let err = |h: &Hull| -> f64 {
            let mut e = 0.0f64;
            for i in 0..97 {
                for j in 0..33 {
                    let x = l * i as f64 / 96.0;
                    let z = t * j as f64 / 32.0;
                    e = e.max((h.surface().eval(x, z) - truth(x, z)).abs());
                }
            }
            e
        };
        let (e1, e2) = (err(&h1), err(&h2));
        assert!(
            e2 < e1,
            "hermite fit error {e2} not below value-only error {e1}"
        );
        assert!(e2 < 2e-3, "hermite fit error {e2}");
    }

    #[test]
    fn zero_weight_excludes_a_poisoned_sample() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let (st, wl, y) = wigley_grid(l, b, t, 41, 17);
        let opts = FitOptions::default();
        let (clean, _) = fit_offsets(&st, &wl, &y, &opts).unwrap();
        let mut bad = y.clone();
        let mut w = vec![1.0; y.len()];
        bad[20 * 17 + 8] = 1e3; // garbage value at an interior sample
        w[20 * 17 + 8] = 0.0; // ... excluded by weight
        let grid = SampleGrid::new(st, wl, bad)
            .unwrap()
            .with_weights(w)
            .unwrap();
        let (fitted, _) = fit_grid(&grid, &opts).unwrap();
        let (ca, cb) = (clean.surface().control(), fitted.surface().control());
        for (a, b) in ca.iter().zip(cb) {
            assert!((a - b).abs() < 1e-9, "{a} vs {b}");
        }
    }

    #[test]
    fn fits_non_representable_data_with_small_residual() {
        // A hull with a sine-modulated waterline: not polynomial, so the fit
        // approximates; residuals must be small and reported.
        let (l, t) = (10.0, 0.8);
        let mx = 61;
        let mz = 21;
        let stations: Vec<f64> = (0..mx).map(|i| l * i as f64 / (mx - 1) as f64).collect();
        let waterlines: Vec<f64> = (0..mz).map(|j| t * j as f64 / (mz - 1) as f64).collect();
        let mut y = vec![0.0; mx * mz];
        for (i, &x) in stations.iter().enumerate() {
            for (j, &z) in waterlines.iter().enumerate() {
                let g = (std::f64::consts::PI * x / l).sin();
                y[i * mz + j] = 0.5 * g * g * (1.0 - (z / t).powi(3)).max(0.0);
            }
        }
        let opts = FitOptions {
            degree_x: 3,
            degree_z: 3,
            n_ctrl_x: 14,
            n_ctrl_z: 9,
            ..FitOptions::default()
        };
        let (hull, report) = fit_offsets(&stations, &waterlines, &y, &opts).unwrap();
        assert!(
            report.max_residual < 2e-3,
            "max residual {}",
            report.max_residual
        );
        assert!(report.rms_residual < 5e-4);
        assert!(hull.displaced_volume() > 0.0);
    }

    #[test]
    fn rejects_bad_grids() {
        let ok_st: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let ok_wl: Vec<f64> = (0..12).map(|j| j as f64 * 0.1).collect();
        let y = vec![1.0; 20 * 12];
        let opts = FitOptions::default();
        // Non-increasing stations.
        let mut bad = ok_st.clone();
        bad[3] = bad[2];
        assert!(fit_offsets(&bad, &ok_wl, &y, &opts).is_err());
        // Waterlines not starting at 0.
        let bad_wl: Vec<f64> = ok_wl.iter().map(|v| v + 0.5).collect();
        assert!(fit_offsets(&ok_st, &bad_wl, &y, &opts).is_err());
        // Wrong grid size.
        assert!(fit_offsets(&ok_st, &ok_wl, &y[..100], &opts).is_err());
        // Clearly negative half-beam.
        let mut bad_y = y.clone();
        bad_y[5] = -0.5;
        assert!(fit_offsets(&ok_st, &ok_wl, &bad_y, &opts).is_err());
        // Too few samples for the control net.
        let few: Vec<f64> = (0..5).map(|i| i as f64).collect();
        assert!(fit_offsets(&few, &ok_wl, &vec![1.0; 5 * 12], &opts).is_err());
        // Bad derivative weight.
        let bad_opts = FitOptions {
            derivative_weight: f64::NAN,
            ..FitOptions::default()
        };
        assert!(fit_offsets(&ok_st, &ok_wl, &y, &bad_opts).is_err());
    }
}

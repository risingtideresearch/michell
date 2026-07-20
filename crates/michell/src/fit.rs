//! Least-squares lofting of gridded half-breadth samples to the canonical
//! B-spline hull representation.
//!
//! This is the shared back-end for every input front-end: a hand-typed
//! station × waterline offset table and a CAD import both reduce to "a grid of
//! half-beam samples", which is fit here by separable tensor-product least
//! squares. The fit is deterministic: interior knots are placed at sample
//! quantiles (which keeps the normal equations well-conditioned whenever the
//! sample grid is denser than the control net), and the achieved residuals
//! are reported so callers can judge fit quality.

use crate::bspline::{self, BSplineSurface};
use crate::error::{Error, Result};
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
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            degree_x: 3,
            degree_z: 3,
            n_ctrl_x: 12,
            n_ctrl_z: 8,
        }
    }
}

/// Fit quality diagnostics, in the units of the input half-beams.
#[derive(Debug, Clone, Copy)]
pub struct FitReport {
    /// Largest |fitted − sample| over the grid.
    pub max_residual: f64,
    /// Sample location (x, z) of the largest residual.
    pub max_residual_at: (f64, f64),
    /// Root-mean-square residual over the grid.
    pub rms_residual: f64,
    /// Most negative control value floored to zero to enforce `f >= 0`
    /// (0.0 if the unconstrained fit was already non-negative).
    pub floored: f64,
}

/// Fit a hull to gridded half-beam samples.
///
/// `stations` (strictly increasing x, metres) and `waterlines` (strictly
/// increasing z downward from the waterline, starting at 0) define a grid;
/// `half_beams[i * waterlines.len() + j]` is the half-beam at
/// `(stations[i], waterlines[j])`. Small negative samples (measurement noise)
/// are clamped to zero; clearly negative values are rejected.
///
/// Requires a few more samples than control points in each direction
/// (`len >= n_ctrl + 2`).
pub fn fit_offsets(
    stations: &[f64],
    waterlines: &[f64],
    half_beams: &[f64],
    opts: &FitOptions,
) -> Result<(Hull, FitReport)> {
    let mx = stations.len();
    let mz = waterlines.len();
    if half_beams.len() != mx * mz {
        return Err(Error::InvalidInput(format!(
            "half_beams has {} entries, expected {} stations x {} waterlines = {}",
            half_beams.len(),
            mx,
            mz,
            mx * mz
        )));
    }
    check_strictly_increasing(stations, "stations")?;
    check_strictly_increasing(waterlines, "waterlines")?;
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
    let y: Vec<f64> = half_beams.iter().map(|&v| v.max(0.0)).collect();

    let (px, pz) = (opts.degree_x, opts.degree_z);
    let (nx, nz) = (opts.n_ctrl_x, opts.n_ctrl_z);
    let kx = approx_knots(stations, px, nx)?;
    let kz = approx_knots(waterlines, pz, nz)?;

    // Dense basis matrices (small: samples x control points).
    let bx = basis_matrix(&kx, px, nx, stations);
    let bz = basis_matrix(&kz, pz, nz, waterlines);

    // Separable normal equations:
    //   C = (Bx' Bx)^{-1} Bx'  Y  Bz (Bz' Bz)^{-1}
    let ax = normal_matrix(&bx, mx, nx);
    let az = normal_matrix(&bz, mz, nz);
    let lx = cholesky(ax, nx).ok_or_else(ill_conditioned)?;
    let lz = cholesky(az, nz).ok_or_else(ill_conditioned)?;

    // R1 = Bx' Y  (nx x mz), then solve Ax C1 = R1 column-wise.
    let mut c1 = vec![0.0f64; nx * mz];
    for k in 0..nx {
        for j in 0..mz {
            let mut s = 0.0;
            for i in 0..mx {
                s += bx[i * nx + k] * y[i * mz + j];
            }
            c1[k * mz + j] = s;
        }
    }
    let mut col = vec![0.0f64; nx];
    for j in 0..mz {
        for k in 0..nx {
            col[k] = c1[k * mz + j];
        }
        chol_solve(&lx, nx, &mut col);
        for k in 0..nx {
            c1[k * mz + j] = col[k];
        }
    }

    // R2 = C1 Bz  (nx x nz), then solve Az c_row = r2_row for each row.
    let mut control = vec![0.0f64; nx * nz];
    let mut row = vec![0.0f64; nz];
    for k in 0..nx {
        for l in 0..nz {
            let mut s = 0.0;
            for j in 0..mz {
                s += c1[k * mz + j] * bz[j * nz + l];
            }
            row[l] = s;
        }
        chol_solve(&lz, nz, &mut row);
        control[k * nz..(k + 1) * nz].copy_from_slice(&row);
    }

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

    // Residuals of the final (floored) control net over the sample grid.
    let mut max_res = 0.0f64;
    let mut max_at = (stations[0], waterlines[0]);
    let mut sum_sq = 0.0f64;
    let mut tmp = vec![0.0f64; nz];
    for i in 0..mx {
        for l in 0..nz {
            let mut s = 0.0;
            for k in 0..nx {
                s += bx[i * nx + k] * control[k * nz + l];
            }
            tmp[l] = s;
        }
        for j in 0..mz {
            let mut fitted = 0.0;
            for l in 0..nz {
                fitted += tmp[l] * bz[j * nz + l];
            }
            let r = (fitted - y[i * mz + j]).abs();
            if r > max_res {
                max_res = r;
                max_at = (stations[i], waterlines[j]);
            }
            sum_sq += r * r;
        }
    }
    let report = FitReport {
        max_residual: max_res,
        max_residual_at: max_at,
        rms_residual: (sum_sq / (mx * mz) as f64).sqrt(),
        floored,
    };

    let surface = BSplineSurface::new(px, pz, kx, kz, control)?;
    Ok((Hull::new(surface)?, report))
}

fn ill_conditioned() -> Error {
    Error::InvalidInput(
        "least-squares normal equations are singular or ill-conditioned; \
         increase sample density or reduce the number of control points"
            .into(),
    )
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
        let mut idx =
            (k as f64 * (m - 1) as f64 / (n_ctrl - degree) as f64).round() as usize;
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

/// Dense basis matrix B (samples x control points, row-major).
fn basis_matrix(knots: &[f64], degree: usize, n_ctrl: usize, samples: &[f64]) -> Vec<f64> {
    let m = samples.len();
    let mut b = vec![0.0f64; m * n_ctrl];
    for (i, &u) in samples.iter().enumerate() {
        let (first, vals) = bspline::basis_row(knots, degree, n_ctrl, u);
        for (j, &v) in vals.iter().enumerate() {
            b[i * n_ctrl + first + j] = v;
        }
    }
    b
}

/// A = B'B (n x n, row-major).
fn normal_matrix(b: &[f64], m: usize, n: usize) -> Vec<f64> {
    let mut a = vec![0.0f64; n * n];
    for i in 0..m {
        let row = &b[i * n..(i + 1) * n];
        for k in 0..n {
            if row[k] == 0.0 {
                continue;
            }
            for l in k..n {
                a[k * n + l] += row[k] * row[l];
            }
        }
    }
    // Mirror to the lower triangle.
    for k in 0..n {
        for l in 0..k {
            a[k * n + l] = a[l * n + k];
        }
    }
    a
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

    fn wigley_grid(l: f64, b: f64, t: f64, mx: usize, mz: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let stations: Vec<f64> = (0..mx)
            .map(|i| -l / 2.0 + l * i as f64 / (mx - 1) as f64)
            .collect();
        let waterlines: Vec<f64> = (0..mz).map(|j| t * j as f64 / (mz - 1) as f64).collect();
        let mut y = vec![0.0; mx * mz];
        for (i, &x) in stations.iter().enumerate() {
            for (j, &z) in waterlines.iter().enumerate() {
                y[i * mz + j] =
                    b / 2.0 * (1.0 - (2.0 * x / l).powi(2)) * (1.0 - (z / t).powi(2));
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
        };
        let (hull, report) = fit_offsets(&st, &wl, &y, &opts).unwrap();
        assert!(
            report.max_residual < 1e-10 * b,
            "max residual {}",
            report.max_residual
        );
        assert_eq!(report.floored, 0.0);
        // The fitted hull must reproduce Wigley resistance.
        let reference = crate::hulls::wigley(l, b, t).unwrap();
        let cond = crate::Conditions::freshwater(3.0);
        let rw_fit = crate::wave_resistance(&hull, &cond).unwrap().resistance;
        let rw_ref = crate::wave_resistance(&reference, &cond).unwrap().resistance;
        assert!(
            (rw_fit - rw_ref).abs() < 1e-6 * rw_ref,
            "Rw {rw_fit} vs {rw_ref}"
        );
        // Geometry too.
        assert!((hull.displaced_volume() - reference.displaced_volume()).abs() < 1e-8);
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
    }
}

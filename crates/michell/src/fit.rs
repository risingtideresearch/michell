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
    /// Fairing weight `λ` (dimensionless; `0` = plain least squares). Adds
    /// `λ · ∬ (f_uu² + 2 f_uv² + f_vv²) du dv` to the mean squared residual,
    /// with `(u, v)` the loft domain mapped to the unit square so that length
    /// and draft are penalised against their own extents. The energy is what
    /// keeps a tensor-product loft from ringing where the data has a crease
    /// no knot line can follow — above all the keel/profile line, where the
    /// half-beam falls to zero at an angle to both parameter directions.
    pub fairing: f64,
}

impl Default for FitOptions {
    fn default() -> Self {
        FitOptions {
            degree_x: 3,
            degree_z: 3,
            n_ctrl_x: 12,
            n_ctrl_z: 8,
            derivative_weight: 1.0,
            fairing: 0.0,
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
    /// How many control values were floored.
    pub held_at_zero: usize,
    /// Residuals of the `∂f/∂x` channel, when the grid carried one (and at
    /// least one sample of it was observed).
    pub fx_residual: Option<ChannelResiduals>,
    /// Residuals of the `∂f/∂z` channel, likewise.
    pub fz_residual: Option<ChannelResiduals>,
    /// Largest half-beam among the fitted samples [m] — the scale the
    /// residuals are meaningful against. `0.0` if nothing was observed.
    pub sample_scale: f64,
}

impl FitReport {
    /// RMS residual as a fraction of the hull's own half-beam scale.
    pub fn relative_rms(&self) -> f64 {
        if self.sample_scale > 0.0 {
            self.rms_residual / self.sample_scale
        } else {
            0.0
        }
    }

    /// Largest residual as a fraction of the hull's own half-beam scale.
    pub fn relative_max(&self) -> f64 {
        if self.sample_scale > 0.0 {
            self.max_residual / self.sample_scale
        } else {
            0.0
        }
    }

    /// Whether the control net is too coarse to hold the sampled geometry.
    ///
    /// A loft that cannot reach its samples is not a cosmetic problem: the
    /// least-squares fit removes exactly the short-scale content of `∂f/∂x`
    /// that feeds the **diverging** (large-`λ`) end of the free-wave
    /// spectrum, so the first thing it biases is `R_w` at low Froude number,
    /// where that end of the spectrum carries the resistance. The cure is
    /// more control points (and enough samples to support them), which the
    /// banded normal equations make cheap.
    ///
    /// The test is on the **RMS** residual, not the largest one. A single
    /// bad sample at a stem tip or a transom corner drives `max_residual` to
    /// a large fraction of a locally tiny half-beam without saying anything
    /// about the hull as a whole, and calibrating against real imports (a
    /// converged ama sits at ~1.2% max-normalised RMS; a visibly
    /// unconverged hull at ~2.3%, and the same hull on a net three times too
    /// coarse at ~8%) puts the useful line just above 2%.
    ///
    /// This is a proxy, and a conservative one: it measures the error in `f`,
    /// while what the wave integral actually sees is the error in `∂f/∂x`. A
    /// grid that carries observed slopes reports those separately in
    /// [`Self::fx_residual`], which is the sharper signal where it exists.
    pub fn under_resolved(&self) -> bool {
        self.relative_rms() > 0.02
    }
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
    if !(opts.fairing.is_finite() && opts.fairing >= 0.0) {
        return Err(Error::InvalidInput(
            "fairing must be finite and non-negative".into(),
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
    // The normal equations are **banded**, not dense: a sample's basis row
    // touches control indices `(x0 + a) * nz + (z0 + c)` for `a <= px`,
    // `c <= pz`, so `A[i][j]` can only be non-zero for `|i - j| <= b`. Storing
    // and factoring the band turns an O(n³) / O(n²)-memory solve into
    // O(n·b²) / O(n·b), which is what makes a finely-resolved control net
    // affordable (a 160x30 net drops from tens of seconds to milliseconds).
    let b = (px * nz + pz).min(n.saturating_sub(1));
    let mut a = vec![0.0f64; n * (b + 1)];
    let mut rhs = vec![0.0f64; n];
    let mut scratch: Vec<(usize, f64)> = Vec::with_capacity((px + 1) * (pz + 1));
    let mut w_total = 0.0f64;
    for (i, (x0, brx)) in rows_x.iter().enumerate() {
        for (j, (z0, brz)) in rows_z.iter().enumerate() {
            let s = i * mz + j;
            let w = weights.map_or(1.0, |ws| ws[s]);
            if w <= 0.0 {
                continue;
            }
            w_total += w;
            add_observation(
                &mut a,
                &mut rhs,
                b,
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
                        b,
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
                        b,
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

    if opts.fairing > 0.0 {
        // The data term is a weighted *sum* of squares; scaling the energy
        // by the total weight makes λ trade it against the *mean*.
        add_fairing(
            &mut a,
            b,
            (&kx, px, nx),
            (&kz, pz, nz),
            opts.fairing * w_total,
        );
    }

    let l = band_cholesky(a, n, b).ok_or_else(ill_conditioned)?;
    let mut control = rhs;
    band_solve(&l, n, b, &mut control);

    // Enforce the hull contract f >= 0: snap numerical dust to zero, floor
    // genuine ringing (reported so the caller can add control points if it
    // is significant). Solving with `c >= 0` as a constraint instead
    // ([`solve_nonneg`]) gives the true constrained optimum, but it is the
    // wrong fix: "every control >= 0" is far stronger than `f >= 0`, the
    // closure onto the keel wants controls below zero, and holding them at
    // zero *during* the solve makes the loft ring more than clipping after
    // it does (on the ama import, curvature sign changes 1,156 -> 1,426).
    let c_scale = control.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    let mut floored = 0.0f64;
    let mut held_at_zero = 0usize;
    for v in control.iter_mut() {
        if v.abs() < 1e-12 * c_scale {
            *v = 0.0;
        } else if *v < 0.0 {
            floored = floored.min(*v);
            held_at_zero += 1;
            *v = 0.0;
        }
    }

    // Residuals of the final (floored) control net, per channel, over the
    // observed (weight > 0) samples.
    let mut max_res = 0.0f64;
    let mut max_at = (st[0], wl[0]);
    let mut sum_sq = 0.0f64;
    let mut n_val = 0usize;
    let mut sample_scale = 0.0f64;
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
            sample_scale = sample_scale.max(f[s].abs());
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
        held_at_zero,
        fx_residual: channel(acc_fx),
        fz_residual: channel(acc_fz),
        sample_scale,
    };

    let surface = BSplineSurface::new(px, pz, kx, kz, control)?;
    Ok((Hull::new(surface)?, report))
}

/// A least-squares surface through **scattered** samples `(u, v, value)` on
/// caller-chosen clamped knots — the building block for lofts in a
/// coordinate other than the hull's own (x, z), e.g. a keel-following depth.
///
/// `fairing` is as in [`FitOptions::fairing`]; some fairing is needed
/// whenever the samples leave parts of the net unconstrained. With
/// `zero_at_v_end` the last row of control values in `v` is pinned to zero,
/// which — the knots being clamped — makes the surface exactly zero along
/// the `v = v_max` edge. The net is constrained non-negative, as in
/// [`fit_grid`], except for the controls `exempt` marks (index `i·n_v + j`):
/// those a trimmed loft never integrates over, where `f >= 0` does not apply. Returns the surface and the RMS residual over the samples.
pub fn fit_scattered(
    samples: &[[f64; 3]],
    (px, pz): (usize, usize),
    knots_u: Vec<f64>,
    knots_v: Vec<f64>,
    fairing: f64,
    zero_at_v_end: bool,
    exempt: Option<&[bool]>,
) -> Result<(BSplineSurface, f64)> {
    if !(fairing.is_finite() && fairing >= 0.0) {
        return Err(Error::InvalidInput(
            "fairing must be finite and non-negative".into(),
        ));
    }
    let (nx, nz) = (knots_u.len() - px - 1, knots_v.len() - pz - 1);
    let n = nx * nz;
    let b = (px * nz + pz).min(n.saturating_sub(1));
    let mut a = vec![0.0f64; n * (b + 1)];
    let mut rhs = vec![0.0f64; n];
    let mut scratch = Vec::with_capacity((px + 1) * (pz + 1));
    let rows: Vec<_> = samples
        .iter()
        .map(|&[u, v, _]| {
            (
                bspline::basis_rows1(&knots_u, px, nx, u),
                bspline::basis_rows1(&knots_v, pz, nz, v),
            )
        })
        .collect();
    for (((x0, bx), (z0, bz)), &[_, _, f]) in rows.iter().zip(samples) {
        add_observation(
            &mut a,
            &mut rhs,
            b,
            nz,
            (*x0, &bx[0]),
            (*z0, &bz[0]),
            1.0,
            f,
            &mut scratch,
        );
    }
    if fairing > 0.0 {
        add_fairing(
            &mut a,
            b,
            (&knots_u, px, nx),
            (&knots_v, pz, nz),
            fairing * samples.len() as f64,
        );
    }
    let mut pinned = vec![false; n];
    if zero_at_v_end {
        for i in 0..nx {
            pinned[i * nz + nz - 1] = true;
        }
    }
    let free = exempt.map_or_else(|| vec![false; n], |e| e.to_vec());
    let control = solve_nonneg(a, rhs, n, b, &pinned, &free)?;
    let mut sum_sq = 0.0;
    for (((x0, bx), (z0, bz)), &[_, _, f]) in rows.iter().zip(samples) {
        let r = tensor_dot((*x0, &bx[0]), (*z0, &bz[0]), &control, nz) - f;
        sum_sq += r * r;
    }
    let rms = (sum_sq / samples.len().max(1) as f64).sqrt();
    Ok((BSplineSurface::new(px, pz, knots_u, knots_v, control)?, rms))
}

/// Minimise `½ cᵀ H c − rᵀ c` subject to `c >= 0` (`c = 0` where
/// `pinned`, and no bound where `unconstrained`), `H` given as its lower
/// band. This is the least-squares fit
/// with the hull's `f >= 0` contract as a constraint, rather than an
/// unconstrained fit with its negative controls clipped afterwards — which
/// is not the solution of any fitting problem, and leaves a bump wherever
/// it clips.
///
/// Primal–dual active set (Hintermüller, Ito & Kunisch 2002): guess which
/// controls sit on the constraint, solve exactly for the rest (one banded
/// Cholesky), then move a control into the set if it came out negative and
/// out of it if holding it at zero costs objective (its multiplier
/// `μ = (Hc − r)ᵢ` went negative). Stops when the set repeats, at which
/// point the KKT conditions hold exactly. It typically settles in a handful
/// of iterations; the cap is a safeguard, after which the last iterate is
/// clipped as before.
fn solve_nonneg(
    h: Vec<f64>,
    r: Vec<f64>,
    n: usize,
    b: usize,
    pinned: &[bool],
    unconstrained: &[bool],
) -> Result<Vec<f64>> {
    const MAX_ITER: usize = 100;
    let diag: Vec<f64> = (0..n).map(|i| h[band_idx(i, i, b)]).collect();
    let mut active: Vec<bool> = pinned.to_vec();
    let mut c = Vec::new();
    for _ in 0..MAX_ITER {
        // Solve with the active controls held at zero: their rows and
        // columns become the identity.
        let mut a = h.clone();
        let mut rhs = r.clone();
        for k in (0..n).filter(|&k| active[k]) {
            for j in k.saturating_sub(b)..=k {
                a[band_idx(k, j, b)] = 0.0;
            }
            for i in k + 1..=(k + b).min(n - 1) {
                a[band_idx(i, k, b)] = 0.0;
            }
            a[band_idx(k, k, b)] = 1.0;
            rhs[k] = 0.0;
        }
        let l = band_cholesky(a, n, b).ok_or_else(ill_conditioned)?;
        band_solve(&l, n, b, &mut rhs);
        c = rhs;
        // Multipliers μ = H c − r (zero on the free set, up to rounding).
        let mut mu: Vec<f64> = r.iter().map(|v| -v).collect();
        for i in 0..n {
            for j in i.saturating_sub(b)..=i {
                let v = h[band_idx(i, j, b)];
                mu[i] += v * c[j];
                if j != i {
                    mu[j] += v * c[i];
                }
            }
        }
        let c_scale = c.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        let next: Vec<bool> = (0..n)
            .map(|k| {
                pinned[k] || !unconstrained[k] && mu[k] - diag[k] * c[k] > 1e-12 * diag[k] * c_scale
            })
            .collect();
        if next == active {
            break;
        }
        active = next;
    }
    // Rounding dust on the free set; a no-op at convergence otherwise.
    for (v, &u) in c.iter_mut().zip(unconstrained) {
        if !u {
            *v = v.max(0.0);
        }
    }
    Ok(c)
}

/// Accumulate one observation row (the tensor product of two basis rows)
/// into the normal equations `A += w r rᵀ`, `rhs += w r t`.
#[allow(clippy::too_many_arguments)]
fn add_observation(
    a: &mut [f64],
    rhs: &mut [f64],
    b: usize,
    nz: usize,
    bx: (usize, &[f64]),
    bz: (usize, &[f64]),
    w: f64,
    target: f64,
    scratch: &mut Vec<(usize, f64)>,
) {
    scratch.clear();
    for (i, &vx) in bx.1.iter().enumerate() {
        for (j, &vz) in bz.1.iter().enumerate() {
            scratch.push(((bx.0 + i) * nz + (bz.0 + j), vx * vz));
        }
    }
    for &(ci, c1) in scratch.iter() {
        rhs[ci] += w * c1 * target;
        for &(cj, c2) in scratch.iter() {
            // Lower band only; the factorisation reads the symmetric half.
            if cj <= ci {
                a[band_idx(ci, cj, b)] += w * c1 * c2;
            }
        }
    }
}

/// Add `scale · ∬ (f_uu² + 2 f_uv² + f_vv²) du dv` to the normal equations,
/// `(u, v)` being the domain mapped to the unit square. The energy of a
/// tensor-product surface is a sum of Kronecker products of 1-D Gram matrices
/// of basis derivatives, each banded to the degree, so it lands entirely
/// inside the band the data term already occupies.
fn add_fairing(
    a: &mut [f64],
    b: usize,
    (kx, px, nx): (&[f64], usize, usize),
    (kz, pz, nz): (&[f64], usize, usize),
    scale: f64,
) {
    let gx: Vec<Vec<f64>> = (0..=2).map(|k| derivative_gram(kx, px, nx, k)).collect();
    let gz: Vec<Vec<f64>> = (0..=2).map(|k| derivative_gram(kz, pz, nz, k)).collect();
    let gx_at = |k: usize, i: usize, j: usize| gx[k][i * (2 * px + 1) + (j + px - i)];
    let gz_at = |k: usize, i: usize, j: usize| gz[k][i * (2 * pz + 1) + (j + pz - i)];
    for i in 0..nx {
        for j in 0..nz {
            let ci = i * nz + j;
            for k in i.saturating_sub(px)..=(i + px).min(nx - 1) {
                for l in j.saturating_sub(pz)..=(j + pz).min(nz - 1) {
                    let cj = k * nz + l;
                    if cj > ci {
                        continue;
                    }
                    let e = gx_at(2, i, k) * gz_at(0, j, l)
                        + 2.0 * gx_at(1, i, k) * gz_at(1, j, l)
                        + gx_at(0, i, k) * gz_at(2, j, l);
                    a[band_idx(ci, cj, b)] += scale * e;
                }
            }
        }
    }
}

/// Banded Gram matrix `G[i][j] = ∫ N_i^(k)(u) N_j^(k)(u) du` over the knot
/// vector's domain mapped to `[0, 1]`, stored `n × (2p + 1)` with `G[i][j]`
/// at `i * (2p + 1) + (j + p - i)`. Exact: `p + 1` Gauss points per span
/// integrate the degree-`2(p - k)` products without error.
fn derivative_gram(knots: &[f64], p: usize, n: usize, k: usize) -> Vec<f64> {
    let w = 2 * p + 1;
    let mut g = vec![0.0f64; n * w];
    let (t0, t1) = (knots[p], knots[n]);
    let len = t1 - t0;
    // d/du = len · d/dx and du = dx / len.
    let factor = len.powi(2 * k as i32 - 1);
    let (nodes, weights) = crate::quadrature::gauss_legendre(p + 1);
    for span in p..n {
        let (a, bnd) = (knots[span], knots[span + 1]);
        if bnd <= a {
            continue;
        }
        let half = 0.5 * (bnd - a);
        for (&t, &wt) in nodes.iter().zip(&weights) {
            let u = a + half * (t + 1.0);
            let d = bspline::ders_basis(knots, p, span, u, k);
            let row = &d[k];
            for (r, &vr) in row.iter().enumerate() {
                for (c, &vc) in row.iter().enumerate() {
                    let (i, j) = (span - p + r, span - p + c);
                    g[i * w + (j + p - i)] += factor * half * wt * vr * vc;
                }
            }
        }
    }
    g
}

/// Index of `A[i][j]` (`i - b <= j <= i`) in lower-band storage.
#[inline]
fn band_idx(i: usize, j: usize, b: usize) -> usize {
    i * (b + 1) + (j + b - i)
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
/// (Schoenberg–Whitney condition in practice). This is the placement
/// [`fit_grid`] uses for each axis of its sample grid.
pub fn quantile_knots(t: &[f64], degree: usize, n_ctrl: usize) -> Result<Vec<f64>> {
    approx_knots(t, degree, n_ctrl)
}

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
fn band_cholesky(mut a: Vec<f64>, n: usize, b: usize) -> Option<Vec<f64>> {
    for k in 0..n {
        let lo_k = k.saturating_sub(b);
        let mut d = a[band_idx(k, k, b)];
        for j in lo_k..k {
            let v = a[band_idx(k, j, b)];
            d -= v * v;
        }
        if d <= 0.0 || !d.is_finite() {
            return None;
        }
        let d = d.sqrt();
        a[band_idx(k, k, b)] = d;
        for i in (k + 1)..(k + b + 1).min(n) {
            // L[i][j] and L[k][j] are both in band only for j >= i - b.
            let lo = i.saturating_sub(b);
            let mut sacc = a[band_idx(i, k, b)];
            for j in lo..k {
                sacc -= a[band_idx(i, j, b)] * a[band_idx(k, j, b)];
            }
            a[band_idx(i, k, b)] = sacc / d;
        }
    }
    Some(a)
}

/// Solve `L Lᵀ x = rhs` in place, `L` in lower-band storage.
fn band_solve(l: &[f64], n: usize, b: usize, rhs: &mut [f64]) {
    for i in 0..n {
        let mut s = rhs[i];
        for j in i.saturating_sub(b)..i {
            s -= l[band_idx(i, j, b)] * rhs[j];
        }
        rhs[i] = s / l[band_idx(i, i, b)];
    }
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..(i + b + 1).min(n) {
            s -= l[band_idx(j, i, b)] * rhs[j];
        }
        rhs[i] = s / l[band_idx(i, i, b)];
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
    fn band_solve_is_exact_on_a_dense_net() {
        // The normal equations are factored in band storage; a wrong band
        // index shows up as a wrong fit, not a crash, so pin the one case
        // with a known answer: a biquadratic is exactly representable by a
        // bicubic net of any size.
        let (l, b, t) = (10.0, 1.0, 0.625);
        let (st, wl, y) = wigley_grid(l, b, t, 401, 81);
        for (nx, nz) in [(20usize, 12usize), (80, 24), (160, 40)] {
            let opts = FitOptions {
                degree_x: 3,
                degree_z: 3,
                n_ctrl_x: nx,
                n_ctrl_z: nz,
                ..FitOptions::default()
            };
            let (hull, rep) = fit_offsets(&st, &wl, &y, &opts).unwrap();
            assert!(
                rep.max_residual < 1e-9,
                "{nx}x{nz}: max residual {} on an exactly representable surface",
                rep.max_residual
            );
            assert!(!rep.under_resolved(), "{nx}x{nz} wrongly flagged");
            // And the geometry that comes out of it is still the Wigley one.
            let exact = 4.0 * b * l * t / 9.0;
            let vol = hull.displaced_volume();
            assert!(
                (vol - exact).abs() < 1e-6 * exact,
                "{nx}x{nz}: volume {vol} vs {exact}"
            );
        }
    }

    #[test]
    fn residuals_fall_with_the_net_and_a_good_fit_is_not_flagged() {
        // A hard slope break in x at 0.7L: the short-scale content of df/dx
        // that a coarse net cannot hold. Refining must reduce the residual
        // monotonically, and a well-resolved fit must not raise the flag.
        let (l, t) = (10.0, 0.6);
        let (mx, mz) = (401, 41);
        let stations: Vec<f64> = (0..mx).map(|i| l * i as f64 / (mx - 1) as f64).collect();
        let waterlines: Vec<f64> = (0..mz).map(|j| t * j as f64 / (mz - 1) as f64).collect();
        let mut y = vec![0.0; mx * mz];
        for (i, &x) in stations.iter().enumerate() {
            let xi = 2.0 * x / l - 1.0;
            let kink = if x < 0.7 * l {
                1.0
            } else {
                1.0 + 2.5 * (x - 0.7 * l) / l
            };
            for (j, &z) in waterlines.iter().enumerate() {
                y[i * mz + j] = 0.5 * (1.0 - xi * xi) * kink * (1.0 - (z / t) * (z / t));
            }
        }
        let fit = |nx, nz| {
            let opts = FitOptions {
                degree_x: 3,
                degree_z: 3,
                n_ctrl_x: nx,
                n_ctrl_z: nz,
                ..FitOptions::default()
            };
            fit_offsets(&stations, &waterlines, &y, &opts).unwrap().1
        };
        let r: Vec<_> = [(6usize, 5usize), (20, 12), (80, 24), (160, 32)]
            .iter()
            .map(|&(nx, nz)| fit(nx, nz))
            .collect();
        for w in r.windows(2) {
            assert!(
                w[1].rms_residual < w[0].rms_residual,
                "refining did not help: {} -> {}",
                w[0].rms_residual,
                w[1].rms_residual
            );
        }
        let fine = r.last().unwrap();
        assert!(fine.sample_scale > 0.0);
        assert!(
            !fine.under_resolved(),
            "well-resolved fit flagged (rms {:.3}%)",
            100.0 * fine.relative_rms()
        );
    }

    #[test]
    fn under_resolved_reads_residuals_against_the_hull_scale() {
        // The flag is pure arithmetic on the report, so pin it directly
        // rather than through a hull that happens to sit near the line.
        // Calibration: a converged real import lands near 1% RMS-on-scale, a
        // visibly unconverged one above 2%.
        let report = |rms: f64, scale: f64| FitReport {
            max_residual: 10.0 * rms,
            max_residual_at: (0.0, 0.0),
            rms_residual: rms,
            floored: 0.0,
            held_at_zero: 0,
            fx_residual: None,
            fz_residual: None,
            sample_scale: scale,
        };
        assert!(!report(0.010, 1.0).under_resolved());
        assert!(!report(0.019, 1.0).under_resolved());
        assert!(report(0.021, 1.0).under_resolved());
        assert!(report(0.21, 10.0).under_resolved());
        // Scale-free: the same relative error at a different hull size.
        assert!(report(0.0021, 0.1).under_resolved());
        // A degenerate fit with nothing observed must not warn.
        assert!(!report(0.5, 0.0).under_resolved());
        assert_eq!(report(0.5, 0.0).relative_rms(), 0.0);
        assert_eq!(report(0.5, 0.0).relative_max(), 0.0);
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

    /// Bending energy `∬ (f_uu² + 2 f_uv² + f_vv²) du dv` of a surface over
    /// its domain mapped to the unit square, by brute-force quadrature on
    /// `eval_deriv` — independent of the Gram-matrix assembly it checks.
    fn bending_energy(s: &BSplineSurface) -> f64 {
        let ((x0, x1), (z0, z1)) = (s.x_domain(), s.z_domain());
        let (lx, lz) = (x1 - x0, z1 - z0);
        let (gx, gw) = crate::quadrature::gauss_legendre(8);
        let spans = |k: &[f64]| -> Vec<(f64, f64)> {
            k.windows(2)
                .filter(|w| w[1] > w[0])
                .map(|w| (w[0], w[1]))
                .collect()
        };
        let mut e = 0.0;
        for (xa, xb) in spans(s.knots_x()) {
            for (za, zb) in spans(s.knots_z()) {
                for (&tx, &wx) in gx.iter().zip(&gw) {
                    let x = xa + 0.5 * (xb - xa) * (tx + 1.0);
                    for (&tz, &wz) in gx.iter().zip(&gw) {
                        let z = za + 0.5 * (zb - za) * (tz + 1.0);
                        let fuu = lx * lx * s.eval_deriv(x, z, 2, 0);
                        let fuv = lx * lz * s.eval_deriv(x, z, 1, 1);
                        let fvv = lz * lz * s.eval_deriv(x, z, 0, 2);
                        let jac = 0.25 * (xb - xa) * (zb - za) / (lx * lz);
                        e += wx * wz * jac * (fuu * fuu + 2.0 * fuv * fuv + fvv * fvv);
                    }
                }
            }
        }
        e
    }

    #[test]
    fn fairing_energy_matches_direct_quadrature() {
        // Non-uniform knots, unequal degrees, an arbitrary control net.
        let (px, pz, nx, nz) = (3, 2, 7, 5);
        let kx = vec![-2.0, -2.0, -2.0, -2.0, -1.2, 0.1, 0.4, 3.0, 3.0, 3.0, 3.0];
        let kz = vec![0.0, 0.0, 0.0, 0.2, 0.5, 0.7, 0.7, 0.7];
        let c: Vec<f64> = (0..nx * nz)
            .map(|k| 1.0 + ((k * 7919) % 13) as f64 * 0.1)
            .collect();
        let n = nx * nz;
        let b = px * nz + pz;
        let mut a = vec![0.0; n * (b + 1)];
        add_fairing(&mut a, b, (&kx, px, nx), (&kz, pz, nz), 1.0);
        // cᵀ A c from the symmetric lower band.
        let mut quad = 0.0;
        for i in 0..n {
            for j in i.saturating_sub(b)..=i {
                let v = a[band_idx(i, j, b)] * c[i] * c[j];
                quad += if i == j { v } else { 2.0 * v };
            }
        }
        let s = BSplineSurface::new(px, pz, kx, kz, c).unwrap();
        let direct = bending_energy(&s);
        assert!(
            (quad - direct).abs() < 1e-10 * direct,
            "banded {quad} vs direct {direct}"
        );
    }

    /// A Wigley-like hull with a rockered keel: `f = B/2 (1 − ξ²)(1 − (z/d)²)`
    /// above the keel depth `d(x) = T(1 − 0.8ξ²)` and zero below it. The
    /// half-beam is creased along the keel line, which runs diagonally across
    /// the parameter grid where no knot line can follow it.
    fn rockered_keel() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let (mx, mz) = (121, 31);
        let st: Vec<f64> = (0..mx)
            .map(|i| -5.0 + 10.0 * i as f64 / (mx - 1) as f64)
            .collect();
        let wl: Vec<f64> = (0..mz).map(|j| 0.5 * j as f64 / (mz - 1) as f64).collect();
        let mut y = Vec::with_capacity(mx * mz);
        for &x in &st {
            let xi = x / 5.0;
            let d = 0.5 * (1.0 - 0.8 * xi * xi);
            for &z in &wl {
                y.push(if z < d {
                    0.5 * (1.0 - xi * xi) * (1.0 - (z / d).powi(2))
                } else {
                    0.0
                });
            }
        }
        (st, wl, y)
    }

    #[test]
    fn fairing_trades_residual_for_bending_energy_monotonically() {
        // A smooth, strictly positive surface the net cannot hold exactly,
        // so the f >= 0 floor never acts and the loft is exactly the
        // penalised optimum: mean-square residual + λ·energy. Its residual
        // must then rise and its energy fall as λ grows.
        let st: Vec<f64> = (0..61).map(|i| 4.0 * i as f64 / 60.0).collect();
        let wl: Vec<f64> = (0..21).map(|j| 0.5 * j as f64 / 20.0).collect();
        let mut y = Vec::new();
        for &x in &st {
            for &z in &wl {
                y.push(1.0 + 0.3 * (3.0 * x).sin() * (4.0 * z).cos());
            }
        }
        let mut last: Option<(f64, f64)> = None;
        for lam in [0.0, 1e-7, 1e-5, 1e-3] {
            let opts = FitOptions {
                n_ctrl_x: 10,
                n_ctrl_z: 6,
                fairing: lam,
                ..FitOptions::default()
            };
            let (hull, rep) = fit_offsets(&st, &wl, &y, &opts).unwrap();
            assert_eq!(rep.floored, 0.0);
            let e = bending_energy(hull.surface());
            if let Some((e0, r0)) = last {
                assert!(e < e0, "λ {lam}: energy {e} not below {e0}");
                assert!(rep.rms_residual > r0, "λ {lam}: rms fell");
            }
            last = Some((e, rep.rms_residual));
        }
    }

    #[test]
    fn fairing_stops_a_loft_ringing_along_a_diagonal_keel() {
        let (st, wl, y) = rockered_keel();
        let mz = wl.len();
        // Curvature sign changes down each section, over the samples the
        // hull actually occupies. The true sections are convex throughout.
        let inflections = |hull: &Hull| -> usize {
            let s = hull.surface();
            let mut count = 0;
            for (i, &x) in st.iter().enumerate() {
                let mut last = 0.0f64;
                for (j, &z) in wl.iter().enumerate() {
                    if y[i * mz + j] <= 1e-3 {
                        break;
                    }
                    let c = s.eval_deriv(x, z, 0, 2);
                    if c.abs() > 1e-2 {
                        if last != 0.0 && c.signum() != last.signum() {
                            count += 1;
                        }
                        last = c;
                    }
                }
            }
            count
        };
        let opts = FitOptions {
            n_ctrl_x: 40,
            n_ctrl_z: 12,
            ..FitOptions::default()
        };
        let (plain, plain_rep) = fit_offsets(&st, &wl, &y, &opts).unwrap();
        let (faired, faired_rep) = fit_offsets(
            &st,
            &wl,
            &y,
            &FitOptions {
                fairing: 1e-7,
                ..opts
            },
        )
        .unwrap();
        let (n_plain, n_faired) = (inflections(&plain), inflections(&faired));
        eprintln!(
            "section inflections {n_plain} -> {n_faired}, rms {:.3e} -> {:.3e}",
            plain_rep.rms_residual, faired_rep.rms_residual
        );
        // Rounding the crease costs each section about one inflection where
        // it turns into the keel; the plain loft rings well beyond that.
        let n_st = st.len();
        assert!(2 * n_plain > 3 * n_st, "plain loft should ring ({n_plain})");
        assert!(5 * n_faired < 6 * n_st, "{n_plain} -> {n_faired}");
        // At this weight the fit is no worse: fairing mostly removes ringing
        // the f >= 0 floor would otherwise have clipped.
        assert!(faired_rep.rms_residual <= plain_rep.rms_residual);
    }

    #[test]
    fn fairing_leaves_a_representable_surface_alone() {
        // A bicubic the net can hold exactly, whose bending energy is tiny
        // next to the data term at a sane λ: the faired loft is still exact
        // to well within the fit tolerance.
        let st: Vec<f64> = (0..21).map(|i| i as f64 / 20.0).collect();
        let wl: Vec<f64> = (0..11).map(|j| 0.5 * j as f64 / 10.0).collect();
        let mut y = Vec::new();
        for &x in &st {
            for &z in &wl {
                y.push(1.0 + 0.2 * x * (1.0 - x) - 0.3 * z * z);
            }
        }
        let opts = FitOptions {
            n_ctrl_x: 6,
            n_ctrl_z: 5,
            fairing: 1e-9,
            ..FitOptions::default()
        };
        let (_, rep) = fit_offsets(&st, &wl, &y, &opts).unwrap();
        assert!(rep.max_residual < 1e-6, "{}", rep.max_residual);
        let bad = FitOptions {
            fairing: -1.0,
            ..opts
        };
        assert!(fit_offsets(&st, &wl, &y, &bad).is_err());
    }

    #[test]
    fn scattered_fit_recovers_a_representable_surface_and_pins_its_edge() {
        // A bicubic that vanishes at v = 1, sampled at scattered points: the
        // fit with a pinned edge reproduces it and is exactly zero there.
        let f = |u: f64, v: f64| (1.0 + 0.3 * u - 0.1 * u * u) * (1.0 - v) * (1.0 + v);
        let mut pts = Vec::new();
        for k in 0..600 {
            // Deterministic scatter over [0, 2] x [0, 1].
            let u = 2.0 * ((k * 7919) % 997) as f64 / 996.0;
            let v = ((k * 104_729) % 991) as f64 / 990.0;
            pts.push([u, v, f(u, v)]);
        }
        let ku = vec![0.0, 0.0, 0.0, 0.0, 0.7, 1.3, 2.0, 2.0, 2.0, 2.0];
        let kv = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        let (surf, rms) = fit_scattered(&pts, (3, 3), ku, kv, 0.0, true, None).unwrap();
        assert!(rms < 1e-10, "rms {rms}");
        for k in 0..=20 {
            let u = 2.0 * k as f64 / 20.0;
            assert_eq!(surf.eval(u, 1.0), 0.0);
            assert!((surf.eval(u, 0.3) - f(u, 0.3)).abs() < 1e-10);
        }
    }
}

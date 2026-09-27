//! The steady **near-field** of a fleet of sectional hulls: the linearised
//! dynamic pressure on the hulls and the free-surface elevation around them,
//! local (non-wave) part included — what [`crate::FreeWaveSpectrum`], which
//! holds only the far-field waves, leaves out abreast of and ahead of the
//! hulls.
//!
//! ## Method
//!
//! Thin-ship theory puts a source sheet of density `σ = −2U ∂f/∂x` on each
//! hull's centreplane (stream toward −x, `z` down). The Kelvin source splits
//! (Havelock) into a Rankine source minus its negative image and a regular
//! wave part:
//!
//! ```text
//! G = e^{−k|z−ζ|} − e^{−k(z+ζ)} + (A+1)·e^{−k(z+ζ)},
//! A + 1 = 2νk/(νk − k_x²) = −2ν sec²θ/(k − k₀),   k₀ = ν sec²θ
//! ```
//!
//! * The **Rankine pair** (`1/r − 1/r'`) is evaluated in physical space, in
//!   closed form, over the sheet sampled as piecewise-constant `σ` cells
//!   between stations and depth rows: `∂/∂x` of a cell reduces to `asinh`s
//!   at its corners. It vanishes on `z = 0`, so it carries none of the
//!   free-surface elevation, and is only needed on the hulls.
//! * The **wave part** is a 2-D wavenumber integral of the sectional hull's
//!   own transforms `q(k_x, κ) = ∬ ∂f/∂x e^{−κζ} e^{−ik_x(ξ−x_c)}`, exactly
//!   as [`crate::squat`] integrates the force: `k` on a shared log grid (one
//!   z-contraction per node), `θ` by Gauss panels, the pole `k₀` by
//!   subtraction (principal value), and the radiation condition as the
//!   half-residue on the dispersion curve, which puts the waves astern.
//!
//! The field at a point is then a sum over wavenumber nodes of
//! `Re(c·e^{ik_x X})·cos(k_y Y)·e^{−kz}`, so many field points share one
//! node set; on the free surface each node's term is separable in `x` and
//! `y` and the grid is a blocked matrix product.
//!
//! The pressure is `p = ρUφ_x` (reported as `C_p = 2φ_x/U`) on the
//! centreplane, and the elevation `ζ = (U/g)φ_x` at `z = 0`, positive up.

use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::michell::{Placement, TransomClosure};
use crate::moments::C64;
use crate::quadrature::gauss_legendre;
use crate::sectional::{SectionalContracted, SectionalHull};
use crate::spectrum::WaveGrid;
use std::f64::consts::{FRAC_PI_2, PI};

/// Sign of the half-residue that makes the waves trail the hull (checked by
/// the upstream-decay test below).
const RADIATION_SIGN: f64 = -1.0;

/// How the near field is resolved.
#[derive(Debug, Clone, Copy)]
pub struct NearFieldOptions {
    /// Transom closure, applied to the wave part and to the sheet alike.
    pub closure: TransomClosure,
    /// Depth rows of each hull's source sheet and pressure mesh.
    pub depth_rows: usize,
    /// Quadrature refinement (panel counts scale with it; the free surface
    /// runs at twice this).
    pub level: usize,
}

impl Default for NearFieldOptions {
    fn default() -> Self {
        NearFieldOptions {
            closure: TransomClosure::default(),
            depth_rows: 40,
            level: 2,
        }
    }
}

/// The pressure on one hull, on a mesh of columns midway between its
/// stations and rows at its sheet's depths.
#[derive(Debug, Clone)]
pub struct HullPressure {
    /// Fleet-frame x of each column [m].
    pub x: Vec<f64>,
    /// Depth of each row below the waterline [m].
    pub depth: Vec<f64>,
    /// Half-beam at each mesh node, `[ix * depth.len() + j]` [m].
    pub half_beam: Vec<f64>,
    /// Pressure coefficient `p/(½ρU²)` at each mesh node, same layout.
    pub cp: Vec<f64>,
    /// The hull's centreplane in the fleet frame [m].
    pub y: f64,
    /// Vertical force of this pressure on the hull, `−2∬ p ∂f/∂z dx dz`
    /// [N], positive up.
    pub force_up: f64,
}

/// One hull's centreplane sheet: `f` sampled at stations × depth nodes, and
/// the corner coefficients of its piecewise-constant `σ` cells.
struct Sheet {
    /// Station x in the fleet frame, increasing.
    xs: Vec<f64>,
    /// Depth nodes, uniform from 0 to the draft.
    depth: Vec<f64>,
    /// Half-beam `[i * nz1 + j]`.
    f: Vec<f64>,
    /// `∂_x` Rankine coefficients at each (station, depth node).
    d: Vec<f64>,
    /// Stations `0..real_from` are the closure's virtual appendage.
    real_from: usize,
    y: f64,
}

impl Sheet {
    fn nz1(&self) -> usize {
        self.depth.len()
    }
}

/// Half-beam of a section curve (`(half-beam, depth)` from top to keel) at
/// depth `z`: the widest crossing, zero below the keel.
fn half_beam_at(curve: &[(f64, f64)], z: f64) -> f64 {
    let Some(&(y0, z0)) = curve.first() else {
        return 0.0;
    };
    if z <= z0 {
        return y0;
    }
    let mut best: f64 = 0.0;
    for w in curve.windows(2) {
        let ((ya, za), (yb, zb)) = (w[0], w[1]);
        let (lo, hi) = (za.min(zb), za.max(zb));
        if z >= lo && z <= hi && hi > lo {
            let t = (z - za) / (zb - za);
            best = best.max(ya + t * (yb - ya));
        }
    }
    best
}

fn build_sheet(
    hull: &SectionalHull,
    pl: Placement,
    u: f64,
    nu: f64,
    opts: &NearFieldOptions,
) -> Sheet {
    let nz = opts.depth_rows.max(4);
    let t = hull.draft().max(1e-9);
    let depth: Vec<f64> = (0..=nz).map(|j| t * j as f64 / nz as f64).collect();
    let mut xs = Vec::new();
    let mut f = Vec::new();
    for (x, c) in hull.curves() {
        xs.push(x + pl.x);
        f.extend(depth.iter().map(|&z| half_beam_at(c, z)));
    }
    // Stations must increase; the aft (low-x) end carries any transom.
    let mut real_from = 0;
    if let Some(tr) = hull.transom() {
        if let Some(lv) = opts
            .closure
            .hollow_length(tr.depth, nu)
            .filter(|&l| l > 0.0)
        {
            let nz1 = depth.len();
            let f_t: Vec<f64> = f[..nz1].to_vec();
            let na = 12;
            let mut app_x = Vec::new();
            let mut app_f = Vec::new();
            for i in (1..=na).rev() {
                let s = i as f64 / na as f64;
                let phi = 1.0 - 3.0 * s * s + 2.0 * s * s * s;
                app_x.push(xs[0] - s * lv);
                app_f.extend(f_t.iter().map(|v| v * phi));
            }
            real_from = na;
            app_x.extend(xs);
            app_f.extend(f);
            xs = app_x;
            f = app_f;
        }
    }
    let nx = xs.len();
    let nz1 = depth.len();
    // σ per cell (i, j) between stations i, i+1 and depth nodes j, j+1.
    let sigma = |i: usize, j: usize| -> f64 {
        if i + 1 >= nx || j + 1 >= nz1 {
            return 0.0;
        }
        let dx = xs[i + 1] - xs[i];
        if dx <= 0.0 {
            return 0.0;
        }
        let fx =
            (f[(i + 1) * nz1 + j] + f[(i + 1) * nz1 + j + 1] - f[i * nz1 + j] - f[i * nz1 + j + 1])
                / (2.0 * dx);
        -2.0 * u * fx
    };
    // Δσ across station edge i in row j (σ is 0 outside the sheet).
    let dsig = |i: usize, j: usize| -> f64 {
        let left = if i == 0 { 0.0 } else { sigma(i - 1, j) };
        left - sigma(i, j)
    };
    let mut d = vec![0.0; nx * nz1];
    for i in 0..nx {
        for j in 0..nz1 {
            let below = if j + 1 < nz1 { dsig(i, j) } else { 0.0 };
            let above = if j > 0 { dsig(i, j - 1) } else { 0.0 };
            d[i * nz1 + j] = above - below;
        }
    }
    Sheet {
        xs,
        depth,
        f,
        d,
        real_from,
        y: pl.y,
    }
}

/// `φ_x` of the Rankine pairs (source minus negative image) of every sheet
/// at `(x, y, z)`, `z` down.
fn rankine_phi_x(sheets: &[Sheet], x: f64, y: f64, z: f64) -> f64 {
    let mut total = 0.0;
    for s in sheets {
        let yy = (y - s.y) * (y - s.y);
        let nz1 = s.nz1();
        for (i, &xe) in s.xs.iter().enumerate() {
            let dx = x - xe;
            let rho = (dx * dx + yy).sqrt();
            if rho < 1e-12 {
                continue;
            }
            let row = &s.d[i * nz1..(i + 1) * nz1];
            for (j, &dij) in row.iter().enumerate() {
                if dij == 0.0 {
                    continue;
                }
                let zeta = s.depth[j];
                total += dij * (((zeta - z) / rho).asinh() - ((zeta + z) / rho).asinh());
            }
        }
    }
    total / (4.0 * PI)
}

/// Wavenumber nodes of the wave part: the field is
/// `Σ_n Σ_m Re(c[n·nm+m]·e^{ik_x(x − cx_m)})·cos(k_y(y − y_m))·e^{−kz}`.
struct Nodes {
    kx: Vec<f64>,
    ky: Vec<f64>,
    k: Vec<f64>,
    c: Vec<C64>,
    nm: usize,
}

#[derive(Clone, Copy)]
enum Mode {
    /// On the hulls: the full range the force integral uses.
    Hull { level: usize },
    /// On a free-surface grid: wavenumbers tapered away above `k_cap`.
    Surface { level: usize, k_cap: f64 },
}

fn build_nodes(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    closure: TransomClosure,
    mode: Mode,
) -> Nodes {
    let u = cond.speed;
    let nu = cond.gravity / (u * u);
    let nm = members.len();
    let l = members.iter().map(|(h, _)| h.length()).fold(0.0, f64::max);
    let t = members
        .iter()
        .map(|(h, _)| h.draft())
        .fold(0.0, f64::max)
        .max(1e-6);
    let (level, k_max, cap) = match mode {
        Mode::Hull { level } => (level, (1.0e4 / l).max(60.0 / t), None),
        Mode::Surface { level, k_cap } => (level, k_cap, Some(k_cap)),
    };
    let taper = |k: f64| -> f64 {
        match cap {
            None => 1.0,
            Some(c) => {
                let a = 0.6 * c;
                if k <= a {
                    1.0
                } else if k >= c {
                    0.0
                } else {
                    let s = (k - a) / (c - a);
                    let v = (0.5 * PI * s).cos();
                    v * v
                }
            }
        }
    };
    let (gx, gw) = gauss_legendre(8);
    let panel = |ka: f64, kb: f64, out: &mut Vec<(f64, f64)>| {
        for (g, &xk) in gx.iter().enumerate() {
            out.push((
                0.5 * (kb - ka) * xk + 0.5 * (ka + kb),
                0.5 * (kb - ka) * gw[g],
            ));
        }
    };
    let k_lo = 0.05 / l;
    let n_log = 12 * level;
    let mut knodes = Vec::new();
    panel(0.0, k_lo, &mut knodes);
    let ratio = (k_max / k_lo).powf(1.0 / n_log as f64);
    let mut ka = k_lo;
    for _ in 0..n_log {
        panel(ka, ka * ratio, &mut knodes);
        ka *= ratio;
    }
    let contract_all = |k: f64| -> Vec<SectionalContracted> {
        members
            .iter()
            .map(|(h, _)| {
                let mut zc = SectionalContracted::default();
                h.contract(k, &mut zc);
                zc
            })
            .collect()
    };
    let zcs: Vec<Vec<SectionalContracted>> =
        crate::parallel::map_indexed(knodes.len(), || (), |_, i| contract_all(knodes[i].0));
    let q_of = |zc: &[SectionalContracted], kx: f64| -> Vec<C64> {
        let mut xm = Vec::with_capacity(8);
        members
            .iter()
            .zip(zc)
            .map(|((h, _), z)| h.q_closed(z, kx, nu, closure, &mut xm))
            .collect()
    };

    // θ nodes for the principal value.
    let theta_top = match mode {
        Mode::Hull { .. } => ((400.0 / l) / k_max).min(1.0).acos().min(FRAC_PI_2 - 1e-7),
        Mode::Surface { .. } => FRAC_PI_2,
    };
    let n_theta = 24 * level;
    let mut thetas = Vec::new();
    for it in 0..n_theta {
        let (a, b) = (
            theta_top * it as f64 / n_theta as f64,
            theta_top * (it + 1) as f64 / n_theta as f64,
        );
        panel(a, b, &mut thetas);
    }
    let pref = 2.0 * u * nu / (PI * PI);
    let per_theta: Vec<Vec<(f64, f64, f64, Vec<C64>)>> = crate::parallel::map_indexed(
        thetas.len(),
        || (),
        |_, it| {
            let (theta, wth) = thetas[it];
            let (s, c) = theta.sin_cos();
            let sec = 1.0 / c;
            let k0 = nu * sec * sec;
            let mut out = Vec::new();
            // The grid, extended past a pole that lies beyond it (on the
            // hulls only: on the surface the taper has killed it).
            let mut grid: Vec<(f64, f64, Option<usize>)> = knodes
                .iter()
                .enumerate()
                .map(|(i, &(k, w))| (k, w, Some(i)))
                .collect();
            let mut k_hi = k_max;
            if matches!(mode, Mode::Hull { .. }) && k0 >= k_max {
                k_hi = 4.0 * k0;
                let n_ext = 4 * level;
                let r = (k_hi / k_max).powf(1.0 / n_ext as f64);
                let mut ka = k_max;
                let mut ext = Vec::new();
                for _ in 0..n_ext {
                    panel(ka, ka * r, &mut ext);
                    ka *= r;
                }
                grid.extend(ext.into_iter().map(|(k, w)| (k, w, None)));
            }
            let h0_live = taper(k0) > 0.0 && k0 < k_hi;
            let mut sum_w = 0.0;
            for (k, wk, zc) in grid {
                let d = k - k0;
                if d.abs() < 1e-9 * k0 {
                    continue;
                }
                sum_w += wk / d;
                let tk = taper(k);
                if tk == 0.0 {
                    continue;
                }
                let q = match zc {
                    Some(i) => q_of(&zcs[i], k * c),
                    None => q_of(&contract_all(k), k * c),
                };
                let wn = wth * pref * sec * k * tk * wk / d;
                // Im(E q) = Re(E · (−i q)).
                out.push((
                    k * c,
                    k * s,
                    k,
                    q.iter().map(|q| C64::new(q.im, -q.re).scale(wn)).collect(),
                ));
            }
            if h0_live {
                let w0 = wth * pref * sec * k0 * taper(k0) * (((k_hi - k0) / k0).ln() - sum_w);
                let q = q_of(&contract_all(k0), k0 * c);
                out.push((
                    k0 * c,
                    k0 * s,
                    k0,
                    q.iter().map(|q| C64::new(q.im, -q.re).scale(w0)).collect(),
                ));
            }
            out
        },
    );

    // Half-residue on the dispersion curve k = ν sec²θ, k_x = ν sec θ.
    let theta_r = match mode {
        Mode::Hull { .. } => {
            let kx_cap = (400.0 / l).max(2.0 * nu);
            (nu / kx_cap).min(1.0).acos()
        }
        Mode::Surface { k_cap, .. } => (nu / k_cap).min(1.0).sqrt().acos(),
    }
    .min(FRAC_PI_2 - 1e-7);
    let n_res = 48 * level;
    let mut rthetas = Vec::new();
    for it in 0..n_res {
        let (a, b) = (
            theta_r * it as f64 / n_res as f64,
            theta_r * (it + 1) as f64 / n_res as f64,
        );
        panel(a, b, &mut rthetas);
    }
    let res_pref = RADIATION_SIGN * 2.0 * u * nu / PI;
    let res: Vec<(f64, f64, f64, Vec<C64>)> = crate::parallel::map_indexed(
        rthetas.len(),
        || (),
        |_, it| {
            let (theta, wth) = rthetas[it];
            let (s, c) = theta.sin_cos();
            let sec = 1.0 / c;
            let k0 = nu * sec * sec;
            let w = wth * res_pref * sec * k0 * taper(k0);
            let q = q_of(&contract_all(k0), k0 * c);
            (k0 * c, k0 * s, k0, q.iter().map(|q| q.scale(w)).collect())
        },
    );

    let mut nodes = Nodes {
        kx: Vec::new(),
        ky: Vec::new(),
        k: Vec::new(),
        c: Vec::new(),
        nm,
    };
    for (kx, ky, k, c) in per_theta.into_iter().flatten().chain(res) {
        if c.iter().all(|v| v.re == 0.0 && v.im == 0.0) {
            continue;
        }
        nodes.kx.push(kx);
        nodes.ky.push(ky);
        nodes.k.push(k);
        nodes.c.extend(c);
    }
    nodes
}

fn check(members: &[(&SectionalHull, Placement)], cond: &Conditions) -> Result<()> {
    cond.validate()?;
    if members.is_empty() {
        return Err(Error::InvalidGeometry("empty fleet".into()));
    }
    Ok(())
}

/// The near-field pressure on every hull of a fleet (see the module docs).
///
/// Only the real hull is covered: with a transom closure the virtual
/// appendage shapes the field but carries no pressure of its own here, so
/// [`HullPressure::force_up`] is the force on the hull. (The thin-ship force
/// of [`crate::sectional::dynamic_force`] instead includes the appendage's
/// share, and on an open transom it omits the end term the integration by
/// parts along the hull drops; see the CAD report test.)
pub fn hull_pressure(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    opts: &NearFieldOptions,
) -> Result<Vec<HullPressure>> {
    hull_pressure_with(members, cond, opts, false)
}

/// [`hull_pressure`], optionally over the closure's virtual appendage too.
fn hull_pressure_with(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    opts: &NearFieldOptions,
    with_appendage: bool,
) -> Result<Vec<HullPressure>> {
    check(members, cond)?;
    let u = cond.speed;
    let nu = cond.gravity / (u * u);
    let sheets: Vec<Sheet> = members
        .iter()
        .map(|(h, p)| build_sheet(h, *p, u, nu, opts))
        .collect();
    let nodes = build_nodes(
        members,
        cond,
        opts.closure,
        Mode::Hull {
            level: opts.level.max(1),
        },
    );
    let cxs: Vec<f64> = members.iter().map(|(h, p)| h.x_center() + p.x).collect();
    let q_dyn = 0.5 * cond.fluid.density * u * u;
    let mut out = Vec::new();
    for sh in &sheets {
        let nz1 = sh.nz1();
        let dz = sh.depth[1] - sh.depth[0];
        // Columns midway between the real hull's stations.
        let first = if with_appendage { 0 } else { sh.real_from };
        let cols: Vec<usize> = (first..sh.xs.len() - 1).collect();
        let per_col: Vec<Vec<f64>> = crate::parallel::map_indexed(
            cols.len(),
            || vec![0.0f64; nz1],
            |acc, ci| {
                let i = cols[ci];
                let x = 0.5 * (sh.xs[i] + sh.xs[i + 1]);
                acc.iter_mut().for_each(|v| *v = 0.0);
                for n in 0..nodes.k.len() {
                    let mut e = 0.0;
                    for m in 0..nodes.nm {
                        let c = nodes.c[n * nodes.nm + m];
                        let ph = C64::cis(nodes.kx[n] * (x - cxs[m]));
                        let re = c.re * ph.re - c.im * ph.im;
                        e += re * (nodes.ky[n] * (sh.y - members[m].1.y)).cos();
                    }
                    let r = (-nodes.k[n] * dz).exp();
                    let mut p = 1.0;
                    for v in acc.iter_mut() {
                        *v += e * p;
                        p *= r;
                    }
                }
                (0..nz1)
                    .map(|j| {
                        let phi_x = acc[j] + rankine_phi_x(&sheets, x, sh.y, sh.depth[j]);
                        2.0 * phi_x / u
                    })
                    .collect()
            },
        );
        let mut hp = HullPressure {
            x: Vec::new(),
            depth: sh.depth.clone(),
            half_beam: Vec::new(),
            cp: Vec::new(),
            y: sh.y,
            force_up: 0.0,
        };
        for (ci, &i) in cols.iter().enumerate() {
            let dx = sh.xs[i + 1] - sh.xs[i];
            hp.x.push(0.5 * (sh.xs[i] + sh.xs[i + 1]));
            let fm: Vec<f64> = (0..nz1)
                .map(|j| 0.5 * (sh.f[i * nz1 + j] + sh.f[(i + 1) * nz1 + j]))
                .collect();
            for j in 0..nz1 - 1 {
                let pm = 0.5 * (per_col[ci][j] + per_col[ci][j + 1]) * q_dyn;
                hp.force_up -= 2.0 * pm * (fm[j + 1] - fm[j]) * dx;
            }
            hp.half_beam.extend(fm);
            hp.cp.extend(per_col[ci].iter().copied());
        }
        out.push(hp);
    }
    Ok(out)
}

/// The free-surface elevation `ζ = (U/g)φ_x` on a grid, local field and waves
/// together (see the module docs). Inside a hull's waterplane the value is
/// the (unphysical) field continued there.
#[allow(clippy::too_many_arguments)]
pub fn free_surface(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    opts: &NearFieldOptions,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    nx: usize,
    ny: usize,
) -> Result<WaveGrid> {
    check(members, cond)?;
    if nx < 2 || ny < 2 {
        return Err(Error::InvalidInput("need at least a 2 x 2 grid".into()));
    }
    let u = cond.speed;
    let dxg = (x1 - x0) / (nx - 1) as f64;
    let dyg = (y1 - y0) / (ny - 1) as f64;
    // Two points per shortest wave the grid can show.
    let k_cap = PI / dxg.max(dyg).max(1e-9);
    let t_nodes = std::time::Instant::now();
    let nodes = build_nodes(
        members,
        cond,
        opts.closure,
        // The oscillatory θ integrals at grid range need twice the panels.
        Mode::Surface {
            level: 2 * opts.level.max(1),
            k_cap,
        },
    );
    let cxs: Vec<f64> = members.iter().map(|(h, p)| h.x_center() + p.x).collect();
    let ys: Vec<f64> = members.iter().map(|(_, p)| p.y).collect();
    let nn = nodes.k.len();
    if std::env::var("NF_TIME").is_ok() {
        eprintln!(
            "nodes {nn} built in {:.2}s",
            t_nodes.elapsed().as_secs_f64()
        );
    }
    // Per node, a_c(x)·cos(k_y y) + a_s(x)·sin(k_y y), where
    // a_c = Σ_m Re(c_m E_m(x)) cos(k_y y_m), a_s the same with sin.
    const BLOCK: usize = 2048;
    let mut zeta = vec![0.0f64; nx * ny];
    let mut start = 0;
    while start < nn {
        let end = (start + BLOCK).min(nn);
        let nb = end - start;
        let mut ac = vec![0.0f64; nb * nx];
        let mut asn = vec![0.0f64; nb * nx];
        for b in 0..nb {
            let n = start + b;
            for m in 0..nodes.nm {
                let c = nodes.c[n * nodes.nm + m];
                let (sy, cy) = (nodes.ky[n] * ys[m]).sin_cos();
                let step = C64::cis(nodes.kx[n] * dxg);
                let mut e = C64::cis(nodes.kx[n] * (x0 - cxs[m]));
                for ix in 0..nx {
                    let re = c.re * e.re - c.im * e.im;
                    ac[b * nx + ix] += re * cy;
                    asn[b * nx + ix] += re * sy;
                    e = C64::new(
                        e.re * step.re - e.im * step.im,
                        e.re * step.im + e.im * step.re,
                    );
                }
            }
        }
        let rows: Vec<Vec<f64>> = crate::parallel::map_indexed(
            ny,
            || (),
            |_, iy| {
                let y = y0 + dyg * iy as f64;
                let mut row = vec![0.0f64; nx];
                for b in 0..nb {
                    let (s, c) = (nodes.ky[start + b] * y).sin_cos();
                    let (a1, a2) = (&ac[b * nx..(b + 1) * nx], &asn[b * nx..(b + 1) * nx]);
                    for ix in 0..nx {
                        row[ix] += a1[ix] * c + a2[ix] * s;
                    }
                }
                row
            },
        );
        for (iy, row) in rows.into_iter().enumerate() {
            for (ix, v) in row.into_iter().enumerate() {
                zeta[iy * nx + ix] += v;
            }
        }
        start = end;
    }
    let scale = u / cond.gravity;
    zeta.iter_mut().for_each(|v| *v *= scale);
    Ok(WaveGrid {
        x0,
        x1,
        y0,
        y1,
        nx,
        ny,
        zeta,
        max_lambda: (k_cap / (cond.gravity / (u * u))).sqrt(),
        theta_samples: nn,
        resolution_limited: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sectional::dynamic_force;
    use crate::squat::SquatOptions;
    use crate::FreeWaveSpectrum;

    fn wigley() -> SectionalHull {
        let w = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let text = crate::iges::write(
            &crate::iges::halfbreadth_surfaces(w.surface(), 0.0, 0.0),
            "wigley",
        )
        .unwrap();
        let f = crate::iges::import_sectional(&text, &Default::default()).unwrap();
        f.hulls[0].hull.clone()
    }

    /// The q-only transform is the full one's `q`.
    #[test]
    fn q_closed_is_the_transform_q() {
        let h = wigley();
        let mut zc = SectionalContracted::default();
        let mut xm = Vec::new();
        for (kx, k) in [(0.3, 0.5), (2.0, 7.0), (11.0, 40.0)] {
            h.contract(k, &mut zc);
            for cl in [TransomClosure::None, TransomClosure::Fixed { length: 0.4 }] {
                let a = h.transforms_at_closed(&zc, kx, 0.8, cl).q;
                let b = h.q_closed(&zc, kx, 0.8, cl, &mut xm);
                assert!((a.re - b.re).abs() + (a.im - b.im).abs() < 1e-12 * (1.0 + a.abs()));
            }
        }
    }

    /// The hull pressure integrates to the force `squat` computes, and the
    /// free surface is quiet upstream and matches the free waves astern.
    #[test]
    fn near_field_matches_force_and_far_field() {
        let h = wigley();
        let cond = Conditions::seawater(0.35 * (9.81f64 * 10.0).sqrt());
        let opts = NearFieldOptions {
            closure: TransomClosure::None,
            ..Default::default()
        };
        let t = std::time::Instant::now();
        let hp = hull_pressure(&[(&h, Placement::default())], &cond, &opts).unwrap();
        let t_hull = t.elapsed().as_secs_f64();
        let so = SquatOptions {
            wave: crate::WaveOptions {
                transom: TransomClosure::None,
                ..Default::default()
            },
            ..Default::default()
        };
        let d = dynamic_force(&h, &cond, 0.0, &so).unwrap();
        eprintln!(
            "force: pressure {:.3} N vs squat {:.3} N ({:.2}%), {t_hull:.2}s",
            hp[0].force_up,
            d.force_up,
            100.0 * (hp[0].force_up / d.force_up - 1.0)
        );
        assert!((hp[0].force_up - d.force_up).abs() < 0.05 * d.force_up.abs());

        let t = std::time::Instant::now();
        let fs = free_surface(
            &[(&h, Placement::default())],
            &cond,
            &opts,
            -35.0,
            25.0,
            -12.0,
            12.0,
            241,
            97,
        )
        .unwrap();
        let t_fs = t.elapsed().as_secs_f64();
        let mut spec = FreeWaveSpectrum::new_sectional(
            &[(&h, Placement::default())],
            &cond,
            TransomClosure::None,
        )
        .unwrap();
        let far = spec
            .elevation_grid(-35.0, 25.0, -12.0, 12.0, 241, 97)
            .unwrap();
        let (mut num, mut den, mut up, mut peak) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for iy in 0..fs.ny {
            for ix in 0..fs.nx {
                let (x, v) = (fs.x(ix), fs.get(ix, iy));
                peak = peak.max(v.abs());
                if x < -25.0 {
                    num += (v - far.get(ix, iy)).powi(2);
                    den += far.get(ix, iy).powi(2);
                }
                if x > 20.0 {
                    up = up.max(v.abs());
                }
            }
        }
        let rel = (num / den).sqrt();
        eprintln!(
            "surface {}x{} with {} nodes in {t_fs:.2}s: astern rms rel {rel:.3e}, \
             upstream max {up:.3e} of peak {peak:.3e}",
            fs.nx, fs.ny, fs.theta_samples
        );
        assert!(rel < 0.05, "astern mismatch {rel}");
        assert!(up < 0.05 * peak, "waves upstream {up} vs {peak}");
    }

    /// On a CAD hull (e12, a shallow transom-sterned canoe body): the hull
    /// pressure's force, its convergence, and how it relates to the `squat`
    /// force with the transom open and closed. A report.
    #[test]
    #[ignore = "CAD convergence report"]
    fn e12_pressure_force_report() {
        let Some(text) = crate::cad_fixture("e12.igs") else {
            return;
        };
        let so = crate::iges::SectionalOptions {
            waterline_z: -0.95,
            ..Default::default()
        };
        let f = crate::iges::import_sectional(&text, &so).unwrap();
        let h = &f.hulls[0].hull;
        let cond = Conditions::seawater(0.3 * (9.81f64 * h.length()).sqrt());
        let sq = SquatOptions {
            wave: crate::WaveOptions {
                transom: TransomClosure::None,
                ..Default::default()
            },
            ..Default::default()
        };
        let d = dynamic_force(h, &cond, 0.0, &sq).unwrap();
        eprintln!("squat force {:.2} N", d.force_up);
        for len in [0.5, 1.0] {
            let cl = TransomClosure::Fixed { length: len };
            let sq = SquatOptions {
                wave: crate::WaveOptions {
                    transom: cl,
                    ..Default::default()
                },
                ..Default::default()
            };
            let d = dynamic_force(h, &cond, 0.0, &sq).unwrap();
            let o = NearFieldOptions {
                closure: cl,
                ..Default::default()
            };
            let hull = hull_pressure(&[(h, Placement::default())], &cond, &o).unwrap();
            let all = hull_pressure_with(&[(h, Placement::default())], &cond, &o, true).unwrap();
            eprintln!(
                "closed {len} m: squat {:.2} N; pressure on hull {:.2} N, with appendage {:.2} N",
                d.force_up, hull[0].force_up, all[0].force_up
            );
        }
        for (rows, level) in [(40, 2), (80, 2), (40, 4)] {
            let o = NearFieldOptions {
                closure: TransomClosure::None,
                depth_rows: rows,
                level,
            };
            let hp = hull_pressure(&[(h, Placement::default())], &cond, &o).unwrap();
            eprintln!("rows {rows} level {level}: {:.2} N", hp[0].force_up);
        }
    }
}

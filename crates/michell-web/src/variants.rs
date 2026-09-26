//! Experimental loft representations, lofted side by side with the current
//! one so they can be compared by eye and by number.
//!
//! The current loft is a tensor-product spline `f(x, z)` over the whole
//! rectangle, zero below the keel: the keel line is a crease no knot line can
//! follow where it runs diagonally. Both alternatives take the keel as a
//! **piecewise-linear** line `d(x)` with breakpoints at the loft's x knots —
//! linear per span being what would keep the Michell inner integrals closed
//! form — fitted to the keel depth each station's samples imply.
//!
//! - **Trimmed**: the same `f(x, z)` spline, fitted only to samples above
//!   `d(x)` and integrated only there; what lies below is "don't care"
//!   (filled smoothly by the fairing term), with zeros observed along the
//!   keel so the trimmed sections close. No crease, no phantom fin, but
//!   the square-root closure of a round bilge onto the keel still has to be
//!   resolved on a diagonal.
//! - **Keel-following**: a spline `f(x, s)` in the normalised depth
//!   `s = z / d(x)`, pinned to zero along the keel `s = 1`, with s knots
//!   graded toward it. The keel is a domain edge and the bilge closure a
//!   parameter line, so a small net can hold both.

use michell::fit::{fit_scattered, quantile_knots};
use michell::sectional::{DepthQuadrature, SectionNodes, SectionalHull};
use michell::{BSplineSurface, Hull, SampleGrid};
use serde_json::{json, Value};

/// Knobs for the experimental variants.
pub struct VariantOptions {
    /// Fairing weight for all three (the don't-care region of the trimmed
    /// loft needs some; a floor of 1e-9 is applied there).
    pub fairing: f64,
    /// Control net of the keel-following loft (x, s).
    pub keel_net: (usize, usize),
    /// Keel segments per x knot span (1 = breakpoints at the knots). More
    /// segments track the profile more closely; the kernel would then split
    /// each span at them, at proportional cost.
    pub keel_subdiv: usize,
}

/// What every variant reports: a structured display mesh in the hull frame,
/// the control net in the hull frame, the loft's residual at every grid
/// sample, and summary numbers.
struct Variant {
    key: &'static str,
    label: &'static str,
    net: (usize, usize),
    /// Display mesh, `mx × mz` row-major (z fastest): (x, z, half-beam).
    mesh: (usize, usize, Vec<[f64; 3]>),
    /// Drop mesh quads whose corners all have zero beam (the rectangle
    /// below the keel, for the untrimmed loft).
    cull_dry: bool,
    control: (usize, usize, Vec<[f64; 3]>),
    /// Loft − sample at each grid sample (None where the sample is unused).
    residual: Vec<Option<f64>>,
    /// Piecewise-linear keel (x, d) breakpoints, if the variant has one.
    keel: Option<Vec<(f64, f64)>>,
    volume: f64,
    seconds: f64,
}

/// Loft the grid three ways. `current` is the hull as loaded (the first
/// variant, unchanged).
pub fn variants(grid: &SampleGrid, current: &Hull, opts: &VariantOptions) -> Value {
    let samples = Samples::new(grid);
    let s = current.surface();
    let net = (s.n_ctrl_x(), s.n_ctrl_z());
    let tag = |v: Result<Variant, String>, key: &str| match v {
        Ok(v) => variant_json(&v, &samples),
        Err(e) => json!({ "key": key, "error": e }),
    };
    let vals = vec![
        variant_json(&current_variant(&samples, current), &samples),
        tag(trimmed_variant(&samples, net, opts), "trimmed"),
        tag(keel_variant(&samples, opts), "keel"),
    ];
    json!({
        "variants": vals,
        "samples": {
            "inflections": samples.inflections(|i, j| samples.y(i, j)),
            "volume": samples.volume(),
        },
    })
}

/// The sample grid with its wetted mask, per-station keel estimates, and
/// the fairness measure every variant is scored by.
struct Samples<'a> {
    g: &'a SampleGrid,
    scale: f64,
    /// Estimated keel depth per station (0 where the station is dry).
    keel: Vec<f64>,
}

impl<'a> Samples<'a> {
    fn new(g: &'a SampleGrid) -> Self {
        let (st, wl) = (g.stations(), g.waterlines());
        let mz = wl.len();
        let mut scale = 0.0f64;
        for i in 0..st.len() {
            for j in 0..mz {
                if Self::used(g, i, j) {
                    scale = scale.max(g.half_beams()[g.idx(i, j)]);
                }
            }
        }
        let wet = |i: usize, j: usize| Self::used(g, i, j) && g.half_beams()[g.idx(i, j)] > 0.0;
        // The keel lies between the deepest wetted sample and the next one
        // down. Extrapolate f² linearly to zero from the two deepest wetted
        // samples — exact for the square-root closure of a round bilge, and
        // clipped into the bracketing interval for anything else.
        let keel = (0..st.len())
            .map(|i| {
                let Some(j) = (0..mz).rev().find(|&j| wet(i, j)) else {
                    return 0.0;
                };
                let hi = if j + 1 < mz { wl[j + 1] } else { wl[j] };
                if j == 0 || !wet(i, j - 1) {
                    return 0.5 * (wl[j] + hi);
                }
                let f1 = g.half_beams()[g.idx(i, j)].powi(2);
                let f0 = g.half_beams()[g.idx(i, j - 1)].powi(2);
                let d = if f0 > f1 {
                    wl[j] + f1 * (wl[j] - wl[j - 1]) / (f0 - f1)
                } else {
                    hi
                };
                d.clamp(wl[j], hi)
            })
            .collect();
        Samples { g, scale, keel }
    }

    fn used(g: &SampleGrid, i: usize, j: usize) -> bool {
        let k = g.idx(i, j);
        g.weights().is_none_or(|w| w[k] > 0.0) && g.half_beams()[k].is_finite()
    }

    fn y(&self, i: usize, j: usize) -> f64 {
        self.g.half_beams()[self.g.idx(i, j)]
    }

    fn wetted(&self, i: usize, j: usize) -> bool {
        Self::used(self.g, i, j) && self.y(i, j) > 1e-3 * self.scale
    }

    /// Curvature sign changes of `val` along waterlines and down sections,
    /// over the wetted samples (second differences, ignoring curvature
    /// below a small fraction of the hull's own).
    fn inflections(&self, val: impl Fn(usize, usize) -> f64) -> [usize; 2] {
        let (st, wl) = (self.g.stations(), self.g.waterlines());
        let (mx, mz) = (st.len(), wl.len());
        let len = st[mx - 1] - st[0];
        let depth = wl[mz - 1] - wl[0];
        let d2 = |t: &[f64], k: usize, a: f64, b: f64, c: f64| {
            let (h0, h1) = (t[k] - t[k - 1], t[k + 1] - t[k]);
            2.0 * ((c - b) / h1 - (b - a) / h0) / (h0 + h1)
        };
        let count = |curv: &mut dyn Iterator<Item = Option<f64>>, tol: f64| {
            let (mut n, mut last) = (0, 0.0f64);
            for c in curv {
                match c {
                    None => last = 0.0,
                    Some(c) if c.abs() >= tol => {
                        if last != 0.0 && c.signum() != last.signum() {
                            n += 1;
                        }
                        last = c;
                    }
                    Some(_) => {}
                }
            }
            n
        };
        let tol_x = 1e-3 * self.scale / (len / 10.0).powi(2);
        let tol_z = 1e-3 * self.scale / (depth / 5.0).powi(2);
        let mut along = 0;
        for j in 0..mz {
            let mut it = (1..mx - 1).map(|i| {
                (self.wetted(i - 1, j) && self.wetted(i, j) && self.wetted(i + 1, j))
                    .then(|| d2(st, i, val(i - 1, j), val(i, j), val(i + 1, j)))
            });
            along += count(&mut it, tol_x);
        }
        let mut down = 0;
        for i in 0..mx {
            let mut it = (1..mz - 1).map(|j| {
                (self.wetted(i, j - 1) && self.wetted(i, j) && self.wetted(i, j + 1))
                    .then(|| d2(wl, j, val(i, j - 1), val(i, j), val(i, j + 1)))
            });
            down += count(&mut it, tol_z);
        }
        [along, down]
    }

    /// Displaced volume of the samples themselves (trapezoidal).
    fn volume(&self) -> f64 {
        2.0 * self.integrate(|i, j| {
            if Self::used(self.g, i, j) {
                self.y(i, j)
            } else {
                0.0
            }
        })
    }

    /// Trapezoidal ∬ over the sample grid of a per-sample value.
    fn integrate(&self, val: impl Fn(usize, usize) -> f64) -> f64 {
        let (st, wl) = (self.g.stations(), self.g.waterlines());
        let w = |t: &[f64], k: usize| {
            let lo = if k > 0 { t[k] - t[k - 1] } else { 0.0 };
            let hi = if k + 1 < t.len() {
                t[k + 1] - t[k]
            } else {
                0.0
            };
            0.5 * (lo + hi)
        };
        let mut v = 0.0;
        for i in 0..st.len() {
            for j in 0..wl.len() {
                v += val(i, j) * w(st, i) * w(wl, j);
            }
        }
        v
    }

    /// Least-squares piecewise-linear keel with breakpoints at `xs` (the
    /// distinct knots of a loft) through the per-station estimates, kept
    /// within [0, draft].
    fn pl_keel(&self, xs: &[f64]) -> Vec<(f64, f64)> {
        let st = self.g.stations();
        let n = xs.len();
        // Normal equations of the hat-function basis: tridiagonal.
        let (mut diag, mut off, mut rhs) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        for (i, &x) in st.iter().enumerate() {
            let k = xs.partition_point(|&b| b <= x).clamp(1, n - 1) - 1;
            let t = ((x - xs[k]) / (xs[k + 1] - xs[k])).clamp(0.0, 1.0);
            let (a, b) = (1.0 - t, t);
            diag[k] += a * a;
            diag[k + 1] += b * b;
            off[k] += a * b;
            rhs[k] += a * self.keel[i];
            rhs[k + 1] += b * self.keel[i];
        }
        // A whisker of smoothing keeps sample-free spans determined.
        for d in diag.iter_mut() {
            *d += 1e-9;
        }
        let d = solve_tridiagonal(&diag, &off, &rhs);
        let draft = *self.g.waterlines().last().unwrap();
        xs.iter()
            .zip(d)
            .map(|(&x, d)| (x, d.clamp(0.0, draft)))
            .collect()
    }
}

fn solve_tridiagonal(diag: &[f64], off: &[f64], rhs: &[f64]) -> Vec<f64> {
    // Symmetric: sub- and super-diagonal are both `off[..n-1]`.
    let n = diag.len();
    let (mut c, mut d) = (vec![0.0; n], vec![0.0; n]);
    c[0] = off[0] / diag[0];
    d[0] = rhs[0] / diag[0];
    for k in 1..n {
        let m = diag[k] - off[k - 1] * c[k - 1];
        c[k] = if k + 1 < n { off[k] / m } else { 0.0 };
        d[k] = (rhs[k] - off[k - 1] * d[k - 1]) / m;
    }
    for k in (0..n - 1).rev() {
        d[k] -= c[k] * d[k + 1];
    }
    d
}

fn keel_at(keel: &[(f64, f64)], x: f64) -> f64 {
    let k = keel
        .partition_point(|&(b, _)| b <= x)
        .clamp(1, keel.len() - 1)
        - 1;
    let ((x0, d0), (x1, d1)) = (keel[k], keel[k + 1]);
    let t = ((x - x0) / (x1 - x0)).clamp(0.0, 1.0);
    d0 + (d1 - d0) * t
}

fn subdivide(breaks: &[f64], m: usize) -> Vec<f64> {
    let m = m.max(1);
    let mut v = Vec::with_capacity(breaks.len() * m);
    for w in breaks.windows(2) {
        for k in 0..m {
            v.push(w[0] + (w[1] - w[0]) * k as f64 / m as f64);
        }
    }
    v.extend(breaks.last());
    v
}

fn distinct(knots: &[f64]) -> Vec<f64> {
    let mut v = knots.to_vec();
    v.dedup();
    v
}

fn display_axis(lo: f64, hi: f64, extra: &[f64], n: usize) -> Vec<f64> {
    let mut v: Vec<f64> = (0..=n)
        .map(|i| lo + (hi - lo) * i as f64 / n as f64)
        .chain(extra.iter().copied().filter(|t| (lo..=hi).contains(t)))
        .collect();
    v.sort_by(f64::total_cmp);
    v.dedup_by(|a, b| (*a - *b).abs() < 1e-9 * (hi - lo));
    v
}

fn greville(knots: &[f64], p: usize, n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| knots[i + 1..=i + p].iter().sum::<f64>() / p.max(1) as f64)
        .collect()
}

/// Mesh and volume over a keel-conforming parameterisation: z = t·d(x),
/// with `eval(x, t, z)` the half-beam there.
fn conforming_mesh(
    keel: &[(f64, f64)],
    xs: &[f64],
    ts: &[f64],
    eval: &dyn Fn(f64, f64, f64) -> f64,
) -> ((usize, usize, Vec<[f64; 3]>), f64) {
    let mut pts = Vec::with_capacity(xs.len() * ts.len());
    for &x in xs {
        let d = keel_at(keel, x);
        for &t in ts {
            pts.push([x, t * d, eval(x, t, t * d)]);
        }
    }
    // Volume 2∫ d(x) ∫ f(x, t d) dt dx by Gauss on each keel segment × t panels.
    const G: [(f64, f64); 4] = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    let mut vol = 0.0;
    for w in keel.windows(2) {
        let (a, b) = (w[0].0, w[1].0);
        let sub = 4;
        for q in 0..sub {
            let (qa, qb) = (
                a + (b - a) * q as f64 / sub as f64,
                a + (b - a) * (q + 1) as f64 / sub as f64,
            );
            for &(gx, wx) in &G {
                let x = 0.5 * (qa + qb) + 0.5 * (qb - qa) * gx;
                let d = keel_at(keel, x);
                let mut inner = 0.0;
                let panels = 16;
                for pnl in 0..panels {
                    let (ta, tb) = (pnl as f64 / panels as f64, (pnl + 1) as f64 / panels as f64);
                    for &(gt, wt) in &G {
                        let t = 0.5 * (ta + tb) + 0.5 * (tb - ta) * gt;
                        inner += 0.5 * (tb - ta) * wt * eval(x, t, t * d);
                    }
                }
                vol += 0.5 * (qb - qa) * wx * d * inner;
            }
        }
    }
    ((xs.len(), ts.len(), pts), 2.0 * vol)
}

fn current_variant(s: &Samples, hull: &Hull) -> Variant {
    let surf = hull.surface();
    let (x0, x1) = surf.x_domain();
    let (_, z1) = surf.z_domain();
    let xs = display_axis(x0, x1, surf.knots_x(), 240);
    let zs = display_axis(0.0, z1, surf.knots_z(), 48);
    let mesh: Vec<[f64; 3]> = xs
        .iter()
        .flat_map(|&x| zs.iter().map(move |&z| [x, z, surf.eval(x, z)]))
        .collect();
    let (gx, gz) = (
        greville(surf.knots_x(), surf.degree_x(), surf.n_ctrl_x()),
        greville(surf.knots_z(), surf.degree_z(), surf.n_ctrl_z()),
    );
    let ctrl = control_points(surf, &gx, &gz, &|_, z| z);
    Variant {
        key: "current",
        label: "Current (x, z)",
        net: (surf.n_ctrl_x(), surf.n_ctrl_z()),
        mesh: (xs.len(), zs.len(), mesh),
        cull_dry: true,
        control: ctrl,
        residual: residuals(s, &|x, z| surf.eval(x, z.min(z1))),
        keel: None,
        volume: hull.displaced_volume(),
        seconds: 0.0,
    }
}

fn control_points(
    surf: &BSplineSurface,
    gx: &[f64],
    gz: &[f64],
    z_of: &dyn Fn(f64, f64) -> f64,
) -> (usize, usize, Vec<[f64; 3]>) {
    let nz = gz.len();
    let c = surf.control();
    let pts = gx
        .iter()
        .enumerate()
        .flat_map(|(i, &x)| {
            gz.iter()
                .enumerate()
                .map(move |(j, &v)| [x, z_of(x, v), c[i * nz + j]])
        })
        .collect();
    (gx.len(), nz, pts)
}

/// Loft − sample at every used grid sample, given the loft's half-beam at
/// (x, z) in the hull frame (zero outside its own domain).
fn residuals(s: &Samples, loft: &dyn Fn(f64, f64) -> f64) -> Vec<Option<f64>> {
    let (st, wl) = (s.g.stations(), s.g.waterlines());
    let mut out = Vec::with_capacity(st.len() * wl.len());
    for (i, &x) in st.iter().enumerate() {
        for (j, &z) in wl.iter().enumerate() {
            out.push(Samples::used(s.g, i, j).then(|| loft(x, z) - s.y(i, j)));
        }
    }
    out
}

fn trimmed_variant(
    s: &Samples,
    net: (usize, usize),
    opts: &VariantOptions,
) -> Result<Variant, String> {
    let t0 = std::time::Instant::now();
    let g = s.g;
    let (st, wl) = (g.stations(), g.waterlines());
    let (px, pz) = (3usize, 3usize);
    let kx = quantile_knots(st, px, net.0).map_err(|e| e.to_string())?;
    let kz = quantile_knots(wl, pz, net.1).map_err(|e| e.to_string())?;
    let keel = s.pl_keel(&subdivide(&distinct(&kx), opts.keel_subdiv));
    // Fit what lies inside the trimmed domain z <= d(x), values only (as the
    // keel-following loft is), plus zeros along the keel line itself: a
    // diagonal trim can't be pinned through the control net, so it is
    // pinned through observations instead, twice per station.
    let mut pts = Vec::new();
    for (i, &x) in st.iter().enumerate() {
        let d = keel_at(&keel, x);
        for (j, &z) in wl.iter().enumerate() {
            if Samples::used(g, i, j) && z <= d {
                pts.push([x, z, s.y(i, j)]);
            }
        }
    }
    let z_max = wl[wl.len() - 1];
    let n_keel = 2 * st.len();
    for k in 0..=n_keel {
        let x = st[0] + (st[st.len() - 1] - st[0]) * k as f64 / n_keel as f64;
        pts.push([x, keel_at(&keel, x).min(z_max), 0.0]);
    }
    // No bound on the control values at all. "Every control >= 0" is only a
    // sufficient condition for f >= 0, and far stronger than the hull
    // needs: the closure onto the keel wants controls below zero, and
    // denying them is what makes a constrained loft ring. Unconstrained, the
    // loft stays non-negative inside the trimmed domain to within a few
    // millimetres (reported as `min_beam`); a production version would
    // bound f itself at points inside the domain instead.
    let (nx, nz) = (kx.len() - px - 1, kz.len() - pz - 1);
    let exempt = vec![true; nx * nz];
    let (surf, _) = fit_scattered(
        &pts,
        (px, pz),
        kx,
        kz,
        opts.fairing.max(1e-9),
        false,
        Some(&exempt),
    )
    .map_err(|e| e.to_string())?;

    let (x0, x1) = surf.x_domain();
    let z1 = surf.z_domain().1;
    let xs = display_axis(x0, x1, &distinct(surf.knots_x()), 240);
    let ts = display_axis(0.0, 1.0, &[], 48);
    let (mesh, volume) = conforming_mesh(&keel, &xs, &ts, &|x, _, z| surf.eval(x, z.min(z1)));
    let (gx, gz) = (
        greville(surf.knots_x(), surf.degree_x(), surf.n_ctrl_x()),
        greville(surf.knots_z(), surf.degree_z(), surf.n_ctrl_z()),
    );
    let control = control_points(&surf, &gx, &gz, &|_, z| z);
    let residual = residuals(s, &|x, z| {
        if z <= keel_at(&keel, x) {
            surf.eval(x, z.min(z1))
        } else {
            0.0
        }
    });
    Ok(Variant {
        key: "trimmed",
        label: "Trimmed at PL keel",
        net,
        mesh,
        cull_dry: false,
        control,
        residual,
        keel: Some(keel),
        volume,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

/// The keel-following surface `f(x, s)`, `s = z/d(x)`, and the
/// piecewise-linear keel `d(x)` it is lofted against.
fn keel_fit(
    s: &Samples,
    opts: &VariantOptions,
) -> Result<(BSplineSurface, Vec<(f64, f64)>), String> {
    let g = s.g;
    let (st, wl) = (g.stations(), g.waterlines());
    let (nx, ns) = opts.keel_net;
    let (px, ps) = (3usize, 3usize);
    if ns < ps + 2 {
        return Err(format!("keel-following net needs at least {} in s", ps + 2));
    }
    let kx = quantile_knots(st, px, nx).map_err(|e| e.to_string())?;
    let keel = s.pl_keel(&subdivide(&distinct(&kx), opts.keel_subdiv));
    // s knots graded toward the keel, where a round bilge closes like a
    // square root: uniform in 1 − (1 − s)^½.
    let n_int = ns - ps - 1;
    let mut ks = vec![0.0; ps + 1];
    for k in 1..=n_int {
        let t = k as f64 / (n_int + 1) as f64;
        ks.push(1.0 - (1.0 - t).powi(2));
    }
    ks.extend(std::iter::repeat_n(1.0, ps + 1));
    // Samples in (x, s): everything inside the keel line, zeros included.
    let mut pts = Vec::new();
    for (i, &x) in st.iter().enumerate() {
        let d = keel_at(&keel, x);
        if d <= 1e-9 * wl[wl.len() - 1] {
            continue;
        }
        for (j, &z) in wl.iter().enumerate() {
            if Samples::used(g, i, j) && z <= d {
                pts.push([x, z / d, s.y(i, j)]);
            }
        }
    }
    let (surf, _) = fit_scattered(&pts, (px, ps), kx, ks, opts.fairing.max(1e-9), true, None)
        .map_err(|e| e.to_string())?;
    Ok((surf, keel))
}

/// The keel-following loft as a [`SectionalHull`], for the physics: each of
/// `stations` cosine-spaced stations integrates its section in `s` (where
/// `f(x, s)` is smooth, the keel being `s = 1`), and the sectional kernel
/// interpolates along x — so it can be compared with the sectional import
/// cut straight from CAD on equal terms.
pub fn keel_following_hull(
    grid: &SampleGrid,
    opts: &VariantOptions,
    stations: usize,
) -> Result<SectionalHull, String> {
    let s = Samples::new(grid);
    let (surf, keel) = keel_fit(&s, opts)?;
    let (x0, x1) = surf.x_domain();
    let xs: Vec<f64> = (0..stations)
        .map(|i| {
            let c = (std::f64::consts::PI * i as f64 / (stations - 1) as f64).cos();
            x0 + (x1 - x0) * (1.0 - c) / 2.0
        })
        .collect();
    let quad = DepthQuadrature::default();
    let s_breaks: Vec<f64> = distinct(surf.knots_z());
    let sections = xs
        .iter()
        .map(|&x| {
            let d = keel_at(&keel, x);
            if d <= 0.0 {
                return SectionNodes::empty();
            }
            let breaks: Vec<f64> = s_breaks.iter().map(|b| b * d).collect();
            SectionNodes::from_depth_function(d, &breaks, |z| surf.eval(x, (z / d).min(1.0)), &quad)
        })
        .collect();
    let n = xs.len();
    let mut knots = vec![xs[0]; 4];
    knots.extend_from_slice(&xs[2..n - 2]);
    knots.extend(std::iter::repeat_n(xs[n - 1], 4));
    SectionalHull::new(3, knots, &xs, sections).map_err(|e| e.to_string())
}

fn keel_variant(s: &Samples, opts: &VariantOptions) -> Result<Variant, String> {
    let t0 = std::time::Instant::now();
    let (nx, ns) = opts.keel_net;
    let (px, ps) = (3usize, 3usize);
    let (surf, keel) = keel_fit(s, opts)?;
    let (x0, x1) = surf.x_domain();
    let xs = display_axis(x0, x1, &distinct(surf.knots_x()), 240);
    let ts = display_axis(0.0, 1.0, &distinct(surf.knots_z()), 48);
    let (mesh, volume) = conforming_mesh(&keel, &xs, &ts, &|x, t, _| surf.eval(x, t));
    let (gx, gs) = (
        greville(surf.knots_x(), px, surf.n_ctrl_x()),
        greville(surf.knots_z(), ps, surf.n_ctrl_z()),
    );
    let control = control_points(&surf, &gx, &gs, &|x, v| v * keel_at(&keel, x));
    let residual = residuals(s, &|x, z| {
        let d = keel_at(&keel, x);
        if d > 0.0 && z <= d {
            surf.eval(x, z / d)
        } else {
            0.0
        }
    });
    Ok(Variant {
        key: "keel",
        label: "Keel-following (x, s)",
        net: (nx, ns),
        mesh,
        cull_dry: false,
        control,
        residual,
        keel: Some(keel),
        volume,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

fn variant_json(v: &Variant, s: &Samples) -> Value {
    let flat = |pts: &[[f64; 3]], k: usize| pts.iter().map(|p| p[k]).collect::<Vec<_>>();
    // Score fairness at the samples: the loft value there is sample + residual.
    let mz = s.g.waterlines().len();
    let inflections = s.inflections(|i, j| s.y(i, j) + v.residual[i * mz + j].unwrap_or(0.0));
    let (mut sum_sq, mut n) = (0.0, 0usize);
    for (i, _) in s.g.stations().iter().enumerate() {
        for j in 0..mz {
            if s.wetted(i, j) {
                if let Some(r) = v.residual[i * mz + j] {
                    sum_sq += r * r;
                    n += 1;
                }
            }
        }
    }
    let wetted_rms = (sum_sq / n.max(1) as f64).sqrt() / s.scale;
    // Volume the loft puts where the samples have no beam, as a fraction of
    // the samples' own volume (one-cell disagreements at a step, e.g. a
    // shallow transom, are real but small; a fin below the keel is not).
    let phantom =
        2.0 * s.integrate(|i, j| {
            if Samples::used(s.g, i, j) && s.y(i, j) == 0.0 {
                v.residual[i * mz + j].unwrap_or(0.0).max(0.0)
            } else {
                0.0
            }
        }) / s.volume();
    json!({
        "key": v.key,
        "label": v.label,
        "net": [v.net.0, v.net.1],
        "mesh": { "nx": v.mesh.0, "nz": v.mesh.1, "x": flat(&v.mesh.2, 0), "z": flat(&v.mesh.2, 1), "y": flat(&v.mesh.2, 2) },
        "cull_dry": v.cull_dry,
        "control": { "nx": v.control.0, "nz": v.control.1, "x": flat(&v.control.2, 0), "z": flat(&v.control.2, 1), "y": flat(&v.control.2, 2) },
        "residual": v.residual,
        "keel": v.keel,
        "stats": {
            "wetted_rms": wetted_rms,
            "inflections": inflections,
            "volume": v.volume,
            "phantom_volume": phantom,
            "min_beam": v.mesh.2.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min),
            "seconds": v.seconds,
        },
    })
}

/// The sectional import as a comparison variant, drawn as the physics uses
/// it: each station at its depth-quadrature nodes (the curve the integral
/// sees), the CAD ray hits those were interpolated from, a see-through
/// surface between stations for orientation only (the kernel interpolates
/// each station's depth integral along x, not a surface), and the
/// interpolated depth-integral curve itself at κ = 0 (sectional area) and at
/// a short wave's decay rate. Not a loft of the samples, so no residuals.
pub fn sectional_variant(imp: &michell::iges::SectionalImport, seconds: f64) -> Value {
    let hull = &imp.hull;
    let stations: Vec<(f64, Vec<(f64, f64)>)> =
        hull.sections().map(|(x, o)| (x, o.to_vec())).collect();
    let nodes = stations.iter().map(|(_, o)| o.len()).max().unwrap_or(0);
    let (mut x, mut z, mut y) = (Vec::new(), Vec::new(), Vec::new());
    for (xs, outline) in &stations {
        for k in 0..nodes {
            // An empty end station collapses to a point on the waterline.
            let (hb, depth) = outline.get(k).copied().unwrap_or((0.0, 0.0));
            x.push(*xs);
            z.push(depth);
            y.push(hb);
        }
    }
    let keel: Vec<(f64, f64)> = stations
        .iter()
        .map(|(x, o)| (*x, o.last().map_or(0.0, |p| p.1)))
        .collect();
    let curve = |kappa: f64| {
        let (st, c) = hull.depth_integral_curve(kappa, 8);
        json!({ "kappa": kappa, "stations": st, "curve": c })
    };
    // A short wave at low speed: λ = 2 at Fn 0.15, κ = νλ² with
    // ν = g/U² = 1/(0.0225 L) — where the weight has moved toward the surface.
    let kappa_short = 4.0 / (0.0225 * hull.length());
    let r = &imp.report;
    json!({
        "key": "sectional",
        "label": "Sectional (CAD)",
        "net": [stations.len(), imp.sections.first().map_or(0, |s| s.1.len())],
        "mesh": { "nx": stations.len(), "nz": nodes, "x": x, "z": z, "y": y },
        "ghost": true,
        "cull_dry": false,
        "control": { "nx": 0, "nz": 0, "x": [], "z": [], "y": [] },
        "residual": null,
        "keel": keel,
        "stations": stations,
        "rays": imp.sections,
        "area": [curve(0.0), curve(kappa_short)],
        "stats": {
            "wetted_rms": null,
            "inflections": null,
            "volume": hull.displaced_volume(),
            "phantom_volume": 0.0,
            "min_beam": 0.0,
            "seconds": seconds,
        },
        "report": {
            "stations": r.stations,
            "dropped_stations": r.dropped_stations,
            "ambiguous_rays": r.ambiguous_rays,
            "max_asymmetry": r.max_asymmetry,
            "x_range": [r.x_range.0, r.x_range.1],
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use michell::{Conditions, TransomClosure, WaveOptions};

    /// Keel-following loft against the sectional import on the same CAD
    /// files: resistance per Froude number. A report, not a check.
    #[test]
    #[ignore = "comparison report"]
    fn keel_following_against_cad_sections() {
        let wave = WaveOptions {
            transom: TransomClosure::None,
            ..WaveOptions::default()
        };
        for (file, wl) in [("ama.igs", 0.0), ("e12.igs", -0.95)] {
            let text =
                std::fs::read_to_string(format!("{}/../../{file}", env!("CARGO_MANIFEST_DIR")))
                    .unwrap();
            let so = michell::iges::SectionalOptions {
                waterline_z: wl,
                ..Default::default()
            };
            let cad = michell::iges::import_sectional(&text, &so)
                .unwrap()
                .hulls
                .remove(0)
                .hull;
            let io = michell::iges::ImportOptions {
                waterline_z: wl,
                ..Default::default()
            };
            let fleet = michell::iges::import_fleet(&text, &io).unwrap();
            let lofted = fleet
                .iter()
                .max_by(|a, b| a.hull.length().total_cmp(&b.hull.length()))
                .unwrap();
            let mut kfs = Vec::new();
            for (net, fairing) in [((20, 12), 1e-7), ((40, 16), 1e-7), ((20, 12), 1e-9)] {
                let o = VariantOptions {
                    fairing,
                    keel_net: net,
                    keel_subdiv: 1,
                };
                kfs.push((
                    format!("KF {}x{} λ{fairing:e}", net.0, net.1),
                    keel_following_hull(&lofted.grid, &o, 121).unwrap(),
                ));
            }
            eprintln!(
                "{file}: volume CAD-sectional {:.5}  {}",
                cad.displaced_volume(),
                kfs.iter()
                    .map(|(n, h)| format!("{n} {:.5}", h.displaced_volume()))
                    .collect::<Vec<_>>()
                    .join("  ")
            );
            let l = cad.length();
            for fnum in [0.15, 0.2, 0.3, 0.4, 0.5] {
                let cond = Conditions::seawater(fnum * (9.81 * l).sqrt());
                let r0 = michell::sectional::wave_resistance(&cad, &cond, &wave)
                    .unwrap()
                    .resistance;
                let rs: Vec<String> = kfs
                    .iter()
                    .map(|(n, h)| {
                        let r = michell::sectional::wave_resistance(h, &cond, &wave)
                            .unwrap()
                            .resistance;
                        format!("{n} {r:.4} ({:+.1}%)", 100.0 * (r / r0 - 1.0))
                    })
                    .collect();
                eprintln!("{file} Fn {fnum}: CAD-sectional {r0:.4}  {}", rs.join("  "));
            }
        }
    }
}

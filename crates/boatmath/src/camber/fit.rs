//! A camber document fitted to a hull: the inverse of [`super::geometry`].
//!
//! camber's hulls are a family: a section swept along a sheer plan,
//! blended between a few stations and trimmed. A hull from CAD is fitted
//! into that family by least squares, through camber's own sweep (the port
//! in [`super`]):
//!
//! 1. **A first guess**, read off the hull. The deck datum is the sheer's
//!    highest point. The sheer plan (a clamped B-spline) is fitted to the
//!    sheer's half-breadth, the trim (camber's PCHIP) to its heights, and the
//!    transom to the hull's aft edge. Each station is the hull cut in
//!    camber's own plane there (vertical, square to the plan), its points
//!    spaced down the girth, the last one carried just past the centreline so
//!    the keel emerges where it should.
//! 2. **Levenberg–Marquardt** over every number at once. At planes spread
//!    along the plan, camber's swept section and the hull's slice in the same
//!    plane are compared both ways: each point of camber's section to the
//!    hull's slice (camber where the hull isn't), and points of the hull's
//!    slice to camber's section (hull camber misses).
//!
//! The numbers are kept to what camber's editor writes: control points in
//! order along x, each section's points stepping down (as positive
//! increments), stations apart in u. Knuckles are smooth, but for a hard
//! chine: if most of the hull's sections turn sharply, the middle point
//! rides the chine as a knuckle (k = 1) all along.

use super::{Document, Hull, Station, V2, V3};
use serde_json::{json, Value};

/// How many of each a fitted document has: few, so it stays editable.
#[derive(Clone, Copy, Debug)]
pub struct Shape {
    pub stations: usize,
    pub points: usize,
    pub plan: usize,
    pub trim: usize,
}

impl Default for Shape {
    fn default() -> Self {
        Shape {
            stations: 3,
            points: 5,
            plan: 4,
            trim: 4,
        }
    }
}

/// How well a fit fits, in metres.
#[derive(Debug, Clone)]
pub struct Report {
    /// camber's sections to the hull's slices: RMS and largest.
    pub rms_to_hull: f64,
    pub max_to_hull: f64,
    /// The hull's slices to camber's sections.
    pub rms_to_camber: f64,
    pub max_to_camber: f64,
    pub iterations: usize,
    /// How far the stations' last points were reached inboard to close
    /// every section.
    pub reach: f64,
    /// Where the fitted document's x = 0 is in the hull's frame.
    pub x_origin: f64,
}

// ---------------------------------------------------------------------------
// The hull to fit: triangles, sliced by camber's station planes
// ---------------------------------------------------------------------------

/// The hull's starboard surface as triangles, in its own frame (z up, the
/// design waterline at 0).
pub struct Target {
    tris: Vec<[V3; 3]>,
    /// Triangle indices by their least x, and each one's x range.
    order: Vec<usize>,
    x_lo: Vec<f64>,
    x_hi: Vec<f64>,
    pub x_min: f64,
    pub x_max: f64,
    pub z_max: f64,
    pub y_max: f64,
}

impl Target {
    /// The first hull of a boatmath geometry (patches or triangles).
    pub fn from_geometry(g: &Value) -> Result<Target, String> {
        let hull = &g["hulls"][0];
        let mut tris: Vec<[V3; 3]> = Vec::new();
        match g["kind"].as_str() {
            Some("nurbs") => {
                for p in hull["patches"].as_array().ok_or("geometry: no patches")? {
                    let s = crate::native::patch_from_value(p)?;
                    let ((u0, u1), (v0, v1)) = match s.trim_uv {
                        Some([a, b, c, d]) => ((a, b), (c, d)),
                        None => (s.u_domain(), s.v_domain()),
                    };
                    let nu = (3 * s.n_ctrl_u).clamp(16, 240);
                    let nv = (3 * s.n_ctrl_v).clamp(16, 120);
                    let grid: Vec<Vec<V3>> = (0..=nu)
                        .map(|i| {
                            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
                            (0..=nv)
                                .map(|j| s.point(u, v0 + (v1 - v0) * j as f64 / nv as f64))
                                .collect()
                        })
                        .collect();
                    for i in 0..nu {
                        for j in 0..nv {
                            let (a, b, c, d) = (
                                grid[i][j],
                                grid[i + 1][j],
                                grid[i + 1][j + 1],
                                grid[i][j + 1],
                            );
                            tris.push([a, b, c]);
                            tris.push([a, c, d]);
                        }
                    }
                }
            }
            Some("mesh") => {
                let v: Vec<V3> = serde_json::from_value(hull["vertices"].clone())
                    .map_err(|e| format!("geometry: {e}"))?;
                let t: Vec<[usize; 3]> = serde_json::from_value(hull["triangles"].clone())
                    .map_err(|e| format!("geometry: {e}"))?;
                tris.extend(t.iter().map(|t| [v[t[0]], v[t[1]], v[t[2]]]));
            }
            k => return Err(format!("geometry kind {k:?}")),
        }
        // Starboard: the hull's side at y ≥ 0 (a centred hull).
        tris.retain(|t| t.iter().any(|p| p[1] > 1e-9));
        if tris.is_empty() {
            return Err("the hull has no surface at y > 0".into());
        }
        let x_lo: Vec<f64> = tris
            .iter()
            .map(|t| t.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min))
            .collect();
        let x_hi: Vec<f64> = tris
            .iter()
            .map(|t| t.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max))
            .collect();
        let mut order: Vec<usize> = (0..tris.len()).collect();
        order.sort_by(|&a, &b| x_lo[a].total_cmp(&x_lo[b]));
        let all = || tris.iter().flatten();
        let (x_min, x_max) = all().fold((f64::INFINITY, f64::NEG_INFINITY), |m, p| {
            (m.0.min(p[0]), m.1.max(p[0]))
        });
        let z_max = all().map(|p| p[2]).fold(f64::NEG_INFINITY, f64::max);
        let y_max = all().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
        Ok(Target {
            tris,
            order,
            x_lo,
            x_hi,
            x_min,
            x_max,
            z_max,
            y_max,
        })
    }

    /// The hull cut by the vertical plane through `p` (in plan) square to the
    /// unit heading `t`: segments, as (n, z) with n along `(t.y, −t.x)` (in
    /// toward the centreline for a heading forward) from `p`, z as the hull's.
    pub fn slice(&self, p: V2, t: V2) -> Vec<[V2; 2]> {
        let n_hat = [t[1], -t[0]];
        // The plane's x over the hull's breadth.
        let (xa, xb) = if t[0].abs() > 1e-3 {
            let at = |y: f64| p[0] - (y - p[1]) * t[1] / t[0];
            let (a, b) = (at(0.0), at(self.y_max));
            (a.min(b), a.max(b))
        } else {
            (f64::NEG_INFINITY, f64::INFINITY)
        };
        let dist = |q: V3| (q[0] - p[0]) * t[0] + (q[1] - p[1]) * t[1];
        let local = |q: V3| -> V2 { [(q[0] - p[0]) * n_hat[0] + (q[1] - p[1]) * n_hat[1], q[2]] };
        let mut out = Vec::new();
        let end = self.order.partition_point(|&i| self.x_lo[i] <= xb);
        for &i in &self.order[..end] {
            if self.x_hi[i] < xa {
                continue;
            }
            let tri = self.tris[i];
            let d = tri.map(dist);
            let mut pts: Vec<V2> = Vec::new();
            for k in 0..3 {
                let (a, b) = (k, (k + 1) % 3);
                if (d[a] < 0.0) != (d[b] < 0.0) {
                    let s = d[a] / (d[a] - d[b]);
                    let q = [
                        tri[a][0] + s * (tri[b][0] - tri[a][0]),
                        tri[a][1] + s * (tri[b][1] - tri[a][1]),
                        tri[a][2] + s * (tri[b][2] - tri[a][2]),
                    ];
                    if q[1] >= -1e-9 {
                        pts.push(local(q));
                    }
                }
            }
            if pts.len() == 2 {
                out.push([pts[0], pts[1]]);
            }
        }
        out
    }
}

fn seg_dist(q: V2, a: V2, b: V2) -> f64 {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let l2 = dx * dx + dy * dy;
    let s = if l2 > 0.0 {
        (((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (q[0] - a[0] - s * dx).hypot(q[1] - a[1] - s * dy)
}

/// A positive length from a free parameter, its logarithm (kept finite).
fn pos(a: f64) -> f64 {
    a.clamp(-30.0, 5.0).exp()
}

// ---------------------------------------------------------------------------
// The parameters: a camber document as a free vector
// ---------------------------------------------------------------------------

/// What the fit holds fixed: the deck datum, the plan's start, the
/// transom, the stations' u and the knuckles.
struct Frame {
    shape: Shape,
    z0: f64,
    x0: f64,
    /// Where the trim's points may go, in x.
    trim_range: (f64, f64),
    /// Each section point's knuckle, the same on every station.
    knuckles: Vec<f64>,
    transom: [V2; 2],
    us: Vec<f64>,
}

impl Frame {
    /// The document a parameter vector makes, in the deck frame (x from the
    /// plan's start, z from the deck datum).
    fn document(&self, q: &[f64]) -> Document {
        let s = self.shape;
        let mut it = q.iter().copied();
        let mut plan = Vec::with_capacity(s.plan);
        let mut x = self.x0;
        for i in 0..s.plan {
            if i > 0 {
                x += pos(it.next().unwrap());
            }
            plan.push([x, it.next().unwrap()]);
        }
        // The trim's x: positive gaps (one more than its points), shared
        // out over its range, so they stay in order and in range.
        let gaps: Vec<f64> = (0..=s.trim)
            .map(|_| pos(it.next().unwrap()) + 0.1)
            .collect();
        let total: f64 = gaps.iter().sum();
        let (lo, hi) = self.trim_range;
        let mut acc = 0.0;
        let trim: Vec<V2> = (0..s.trim)
            .map(|i| {
                acc += gaps[i];
                [lo + (hi - lo) * acc / total, it.next().unwrap()]
            })
            .collect();
        let stations = self
            .us
            .iter()
            .map(|&u| {
                let mut pts = vec![[0.0, 0.0, 1.0]];
                let mut z = 0.0;
                for i in 1..s.points {
                    let n = it.next().unwrap();
                    z -= pos(it.next().unwrap());
                    let k = if i == s.points - 1 {
                        1.0
                    } else {
                        self.knuckles[i]
                    };
                    pts.push([n, z, k]);
                }
                Station { u, pts }
            })
            .collect();
        let shift = |p: V2| [p[0] - self.x0, p[1]];
        Document {
            name: String::new(),
            unit_m: 1.0,
            waterline: self.z0,
            deck_trim: 0.0,
            plan: plan.into_iter().map(shift).collect::<Vec<_>>(),
            trim: trim.into_iter().map(shift).collect(),
            transom: self.transom.map(shift),
            stations,
        }
    }
}

/// A document's parameter vector (the inverse of [`Frame::document`]).
fn params(shape: Shape, plan: &[V2], trim: &[V2], range: V2, stations: &[Vec<V2>]) -> Vec<f64> {
    let mut q = Vec::new();
    for (i, p) in plan.iter().enumerate() {
        if i > 0 {
            q.push((p[0] - plan[i - 1][0]).max(1e-6).ln());
        }
        q.push(p[1]);
    }
    let mut x = range[0];
    for p in trim {
        q.push((p[0] - x).max(1e-9).ln());
        x = p[0];
    }
    q.push((range[1] - x).max(1e-9).ln());
    q.extend(trim.iter().map(|p| p[1]));
    for st in stations {
        let mut z = 0.0;
        for p in st.iter().skip(1).take(shape.points - 1) {
            q.push(p[0]);
            q.push((z - p[1]).max(1e-6).ln());
            z = p[1];
        }
    }
    q
}

// ---------------------------------------------------------------------------
// The misfit
// ---------------------------------------------------------------------------

/// The u of comparison plane m: cosine-spaced, so closer toward the ends.
fn plane_u(m: usize) -> f64 {
    0.5 * (1.0 - (std::f64::consts::PI * (m as f64 + 0.5) / PLANES as f64).cos())
}

/// Planes along the plan, and samples of the hull's slice in each.
const PLANES: usize = 40;
const SAMPLES: usize = 12;
const ROWS: usize = 4;
/// How much an open section counts against a fit, per metre short.
const OPEN_WEIGHT: f64 = 10.0;
/// How sharply a section must turn to be a hard chine [deg].
const CHINE_TURN: f64 = 20.0;
const OPEN_CHECKS: usize = 40;

struct Misfit {
    residuals: Vec<f64>,
    to_hull: Vec<f64>,
    to_camber: Vec<f64>,
}

fn misfit(target: &Target, frame: &Frame, q: &[f64]) -> Misfit {
    misfit_of(target, &frame.document(q), frame.z0, frame.x0)
}

/// The misfit of a document whose deck datum is at height `z0` and whose
/// x = 0 is at `x0` in the hull's frame.
fn misfit_of(target: &Target, doc: &Document, z0: f64, x0: f64) -> Misfit {
    let hull = Hull::new(doc);
    let n_cam = ROWS * (doc.stations[0].pts.len() - 1) + 1;
    let mut residuals = Vec::with_capacity(PLANES * (n_cam + SAMPLES + 1));
    let (mut to_hull, mut to_camber) = (Vec::new(), Vec::new());
    for m in 0..PLANES {
        // Closer toward the ends, where a stem or a transom turns.
        let u = plane_u(m);
        let sec = hull.section(u);
        // The plane in the hull's frame.
        let p = [sec.p[0] + x0, sec.p[1]];
        let t = [-sec.n_hat[1], sec.n_hat[0]];
        let slice = target.slice(p, t);
        // camber's section in the plane's (n, z), z in the hull's frame.
        let col = hull.column(u, ROWS);
        let cam: Vec<V2> = col
            .map(|c| {
                c.pts
                    .iter()
                    .map(|w| {
                        let d = [w[0] - sec.p[0], w[1] - sec.p[1]];
                        [d[0] * sec.n_hat[0] + d[1] * sec.n_hat[1], w[2] + z0]
                    })
                    .collect()
            })
            .unwrap_or_default();
        // camber to the hull: each of its points to the slice, or, with no
        // slice here, its distance from its own top (a section that
        // shouldn't be there, shrinking to nothing).
        for k in 0..n_cam {
            let r = match cam.get(k) {
                None => 0.0,
                Some(&c) if slice.is_empty() => (c[0] - cam[0][0]).hypot(c[1] - cam[0][1]),
                Some(&c) => {
                    let d = slice
                        .iter()
                        .map(|s| seg_dist(c, s[0], s[1]))
                        .fold(f64::INFINITY, f64::min);
                    to_hull.push(d);
                    d
                }
            };
            residuals.push(r);
        }
        // The hull to camber: its slice's points at heights spread from top
        // to bottom, each to camber's section, or with none, to the plan.
        // camber's section is closed as the hull's is: level across to the
        // centreline from where it ends (on the transom, the transom's face;
        // on the keel, nothing). In this vertical plane the centreline is at
        // one n.
        let mut closed = cam.clone();
        if let (Some(&e), true) = (cam.last(), sec.n_hat[1].abs() > 1e-9) {
            closed.push([-sec.p[1] / sec.n_hat[1], e[1]]);
        }
        // The hull's end face, which camber's sections stop at rather than
        // include, is left out: points on or near the transom's plane (an
        // end face can stand a little forward of it, where the first
        // section closes it). camber's sections still answer to the hull
        // there.
        let on_transom = |h: &V2| {
            let x = sec.p[0] + h[0] * sec.n_hat[0];
            x - hull.x_transom(h[1] - z0) < 0.02 * doc.loa()
        };
        let mut pts: Vec<V2> = slice
            .iter()
            .flatten()
            .copied()
            .filter(|h| !on_transom(h))
            .collect();
        pts.sort_by(|a, b| b[1].total_cmp(&a[1]));
        for k in 0..SAMPLES {
            let r = if pts.is_empty() {
                0.0
            } else {
                let h = pts[((k as f64 + 0.5) / SAMPLES as f64 * pts.len() as f64) as usize];
                let d = if closed.len() < 2 {
                    h[0].hypot(h[1] - z0)
                } else {
                    closed
                        .windows(2)
                        .map(|w| seg_dist(h, w[0], w[1]))
                        .fold(f64::INFINITY, f64::min)
                };
                to_camber.push(d);
                d
            };
            residuals.push(r);
        }
    }
    // Every section must close, on the centreline or the transom: one that
    // stops short is an open bottom, as far short as its end is off the
    // centreline. Checked finer than the planes, where it's cheap.
    for m in 0..=OPEN_CHECKS {
        let u = m as f64 / OPEN_CHECKS as f64;
        residuals.push(match hull.column(u, 1) {
            Some(c) if c.bottom == super::End::Sheet => OPEN_WEIGHT * c.pts[c.pts.len() - 1][1],
            _ => 0.0,
        });
    }
    Misfit {
        residuals,
        to_hull,
        to_camber,
    }
}

// ---------------------------------------------------------------------------
// The first guess
// ---------------------------------------------------------------------------

/// The hull cut square to x at `x`: its top (sheer) and bottom points, and
/// the top's half-breadth.
fn x_slice(target: &Target, x: f64) -> Option<(V2, f64)> {
    let s = target.slice([x, 0.0], [1.0, 0.0]);
    // n here is −y (the heading is +x, so n̂ = (0, −1)).
    let top = s.iter().flatten().max_by(|a, b| a[1].total_cmp(&b[1]))?;
    let bot = s.iter().flatten().min_by(|a, b| a[1].total_cmp(&b[1]))?;
    Some(([-top[0], top[1]], bot[1]))
}

/// A starting point: how a station's points spread down its girth (point
/// i at the fraction (i/(S−1))^spread of it from the top, so a spread
/// under 1 puts them low, where a section turns, below a straight
/// topside).
#[derive(Clone, Copy, Debug)]
struct Start {
    spread: f64,
}

/// Planes along the plan the first guess's stations are fitted through.
const GUESS_PLANES: usize = 60;

const STARTS: [Start; 3] = [
    Start { spread: 1.0 },
    Start { spread: 0.6 },
    Start { spread: 0.35 },
];

fn first_guess(target: &Target, shape: Shape, start: Start) -> Result<(Frame, Vec<f64>), String> {
    let l = target.x_max - target.x_min;
    let z0 = target.z_max;
    // The sheer, along the hull.
    let sheer: Vec<(f64, f64, f64)> = (0..=60)
        .filter_map(|i| {
            let x = target.x_min + l * (0.005 + 0.99 * i as f64 / 60.0);
            x_slice(target, x).map(|(top, _)| (x, top[0], top[1]))
        })
        .collect();
    if sheer.len() < 8 {
        return Err("the hull is too small to fit".into());
    }
    let depth_at = |x: f64| x_slice(target, x).map_or(0.0, |(top, bot)| top[1] - bot);
    let max_depth = sheer.iter().map(|s| depth_at(s.0)).fold(0.0, f64::max);

    // The transom: the hull's aft edge where the hull has breadth there (a
    // transom has width; a sternpost doesn't), as a line x(z) by least
    // squares, raked or upright. With none, the transom goes out of the
    // way, aft.
    let z_bot = target
        .tris
        .iter()
        .flatten()
        .map(|p| p[2])
        .fold(f64::INFINITY, f64::min);
    let band = 0.02 * (z0 - z_bot);
    let edge: Vec<V2> = (1..24)
        .filter_map(|k| {
            let z = z_bot + (z0 - z_bot) * k as f64 / 24.0;
            let x = target
                .tris
                .iter()
                .flatten()
                .filter(|p| (p[2] - z).abs() < band && p[0] < target.x_min + 0.3 * l)
                .map(|p| p[0])
                .min_by(f64::total_cmp)?;
            // The half-breadth just forward of it, at this height.
            let hb = target
                .slice([x + 0.01 * l, 0.0], [1.0, 0.0])
                .iter()
                .flatten()
                .filter(|h| (h[1] - z).abs() < band)
                .map(|h| -h[0])
                .fold(0.0, f64::max);
            (hb > 0.15 * target.y_max).then_some([x, z])
        })
        .collect();
    let has_transom = edge.len() >= 3;
    let transom = if has_transom {
        let n = edge.len() as f64;
        let (mx, mz) = (
            edge.iter().map(|p| p[0]).sum::<f64>() / n,
            edge.iter().map(|p| p[1]).sum::<f64>() / n,
        );
        let szz: f64 = edge.iter().map(|p| (p[1] - mz).powi(2)).sum();
        let szx: f64 = edge.iter().map(|p| (p[1] - mz) * (p[0] - mx)).sum();
        let slope = if szz > 0.0 { szx / szz } else { 0.0 };
        let x_at = |z: f64| mx + slope * (z - mz);
        // Past the hull's top and bottom, so it cuts all the way.
        let (zt, zb) = (z0 + 0.1 * (z0 - z_bot), z_bot - 0.5 * (z0 - z_bot));
        [[x_at(zt), zt - z0], [x_at(zb), zb - z0]]
    } else {
        let x = target.x_min - 0.05 * l;
        [[x, 0.5 * max_depth], [x, -2.0 * max_depth]]
    };
    let x0 = if has_transom {
        target.x_min.min(transom[0][0]).min(transom[1][0]) - 0.01 * l
    } else {
        target.x_min
    };
    let x_bow = target.x_max;

    // The plan: control points evenly in x, their half-breadths by least
    // squares to the sheer's (the B-spline is linear in them once its x are
    // fixed), the bow's at the centreline.
    let n = shape.plan;
    let xs: Vec<f64> = (0..n)
        .map(|i| x0 + (x_bow - x0) * i as f64 / (n - 1) as f64)
        .collect();
    let basis = |x: f64| -> Vec<f64> {
        let pts: Vec<V2> = xs.iter().map(|&x| [x, 0.0]).collect();
        let plan = super::Plan::new(pts.clone());
        let (mut lo, mut hi) = (0.0, 1.0);
        for _ in 0..50 {
            let m = 0.5 * (lo + hi);
            if plan.at(m)[0] < x {
                lo = m;
            } else {
                hi = m;
            }
        }
        let u = 0.5 * (lo + hi);
        (0..n)
            .map(|i| {
                let mut p = pts.clone();
                p[i][1] = 1.0;
                super::Plan::new(p).at(u)[1]
            })
            .collect()
    };
    let mut ata = vec![vec![0.0; n]; n];
    let mut atb = vec![0.0; n];
    let mut add = |b: &[f64], y: f64, w: f64| {
        for i in 0..n {
            atb[i] += w * b[i] * y;
            for j in 0..n {
                ata[i][j] += w * b[i] * b[j];
            }
        }
    };
    for &(x, y, _) in &sheer {
        add(&basis(x), y, 1.0);
    }
    add(&basis(x_bow), 0.0, 10.0);
    let ys = super::Lu::new(ata).solve(&atb);
    let plan: Vec<V2> = xs.iter().zip(&ys).map(|(&x, &y)| [x, y]).collect();

    // The trim: the sheer's height below the deck datum, at evenly spaced x.
    let trim_x: Vec<f64> = (0..shape.trim)
        .map(|i| target.x_min + (x_bow - target.x_min) * i as f64 / (shape.trim - 1) as f64)
        .collect();
    let sheer_z = |x: f64| {
        let i = sheer
            .iter()
            .position(|s| s.0 >= x)
            .unwrap_or(sheer.len() - 1)
            .max(1);
        let (a, b) = (sheer[i - 1], sheer[i]);
        let f = ((x - a.0) / (b.0 - a.0)).clamp(0.0, 1.0);
        a.2 + f * (b.2 - a.2)
    };
    let trim: Vec<V2> = trim_x.iter().map(|&x| [x, sheer_z(x) - z0]).collect();

    // The stations. camber's loft is linear in its stations' points (a
    // Catmull–Rom in u), so they're solved for by least squares from the
    // hull's sections all along it: each cut in camber's plane, its points
    // spread down its girth, the last carried past the centreline. That
    // reaches a station the hull itself doesn't (one at the plan's end,
    // behind a transom or past the stem) from the sections it shapes.
    let us: Vec<f64> = (0..shape.stations)
        .map(|j| j as f64 / (shape.stations - 1).max(1) as f64)
        .collect();
    let plan_curve = super::Plan::new(plan.clone());
    let samples: Vec<(f64, Vec<V2>)> = (0..=GUESS_PLANES)
        .filter_map(|m| {
            let u = m as f64 / GUESS_PLANES as f64;
            let (p, d) = (plan_curve.at(u), plan_curve.d(u));
            let l = d[0].hypot(d[1]);
            let slice = target.slice(p, [d[0] / l, d[1] / l]);
            let pts: Vec<V2> = slice
                .iter()
                .flatten()
                .copied()
                .filter(|h| p[0] + h[0] * d[1] / l - x_transom_at(&transom, h[1] - z0) > 1e-3 * l)
                .collect();
            (pts.len() >= 16).then_some((u, pts))
        })
        .collect();
    // A hard chine, if most sections have one: the middle point rides it,
    // a knuckle all along.
    let chines: Vec<Option<f64>> = samples.iter().map(|(_, p)| chine_at(p, z0)).collect();
    let found: Vec<f64> = chines.iter().flatten().copied().collect();
    let chined = shape.points >= 5 && 2 * found.len() > samples.len();
    let typical = if found.is_empty() {
        0.5
    } else {
        let mut f = found.clone();
        f.sort_by(f64::total_cmp);
        f[f.len() / 2]
    };
    let samples: Vec<(f64, Vec<V2>)> = samples
        .iter()
        .zip(&chines)
        .map(|((u, pts), c)| {
            let chine = chined.then(|| c.unwrap_or(typical));
            (
                *u,
                section_points(pts, shape.points, z0, start.spread, chine),
            )
        })
        .collect();
    let mut knuckles = vec![0.0; shape.points];
    if chined {
        knuckles[chine_index(shape.points)] = 1.0;
    }
    if samples.len() < shape.stations + 1 {
        return Err("too few of the hull's sections to fit stations to".into());
    }
    // Each station's weight in the loft at each sample's u.
    let weights: Vec<Vec<f64>> = samples
        .iter()
        .map(|(u, _)| {
            (0..shape.stations)
                .map(|j| {
                    let sts: Vec<Station> = us
                        .iter()
                        .enumerate()
                        .map(|(k, &uk)| Station {
                            u: uk,
                            pts: vec![[if k == j { 1.0 } else { 0.0 }, 0.0, 0.0]; 2],
                        })
                        .collect();
                    super::Loft::new(&sts).at(*u).0[0][0]
                })
                .collect()
        })
        .collect();
    let k = shape.stations;
    let mut ata = vec![vec![0.0; k]; k];
    for w in &weights {
        for a in 0..k {
            for b in 0..k {
                ata[a][b] += w[a] * w[b];
            }
        }
    }
    // A little pull toward the neighbours' values keeps a station the
    // samples hardly see from running off.
    for a in 0..k {
        ata[a][a] += 1e-6 * samples.len() as f64;
    }
    let lu = super::Lu::new(ata);
    let mut stations: Vec<Vec<V2>> = vec![vec![[0.0, 0.0]; shape.points]; k];
    for i in 1..shape.points {
        for c in 0..2 {
            let mut atb = vec![0.0; k];
            for (w, (_, pts)) in weights.iter().zip(&samples) {
                for a in 0..k {
                    atb[a] += w[a] * pts[i][c];
                }
            }
            for (j, v) in lu.solve(&atb).into_iter().enumerate() {
                stations[j][i][c] = v;
            }
        }
    }
    // Each point a little below the last, as camber's editor keeps them.
    for st in stations.iter_mut() {
        for i in 1..st.len() {
            st[i][1] = st[i][1].min(st[i - 1][1] - 1e-4);
        }
    }
    let trim_range = (x0 - 0.02 * l, x_bow + 0.02 * l);
    let q = params(shape, &plan, &trim, [trim_range.0, trim_range.1], &stations);
    Ok((
        Frame {
            shape,
            z0,
            x0,
            trim_range,
            knuckles,
            transom,
            us,
        },
        q,
    ))
}

/// The transom's x at height z (the deck frame), as camber reads it.
fn x_transom_at(t: &[V2; 2], z: f64) -> f64 {
    let dz = t[1][1] - t[0][1];
    t[0][0] + (t[1][0] - t[0][0]) * ((z - t[0][1]) / if dz != 0.0 { dz } else { 1.0 })
}

/// A section's points from the hull's slice in its plane (points (n, z),
/// z in the hull's frame): the deck point, points spread down the girth,
/// and the last just past the centreline.
fn section_points(slice: &[V2], s: usize, z0: f64, spread: f64, chine: Option<f64>) -> Vec<V2> {
    let (pts, g) = girth(slice, z0);
    let total = g[g.len() - 1].max(1e-9);
    let at = |f: f64| {
        pts[g
            .iter()
            .position(|&x| x >= f * total)
            .unwrap_or(pts.len() - 1)]
    };
    // Where each interior point goes down the girth: spread from the top,
    // or with a chine, the middle one on it and the rest either side.
    let fracs: Vec<f64> = match chine {
        None => (1..s - 1)
            .map(|i| (i as f64 / (s - 1) as f64).powf(spread))
            .collect(),
        Some(fc) => {
            let c = chine_index(s);
            (1..s - 1)
                .map(|i| {
                    if i < c {
                        fc * (i as f64 / c as f64).powf(spread)
                    } else if i == c {
                        fc
                    } else {
                        fc + (1.0 - fc) * (i - c) as f64 / (s - 1 - c) as f64
                    }
                })
                .collect()
        }
    };
    let mut out = vec![[0.0, 0.0]];
    out.extend(fracs.iter().map(|&f| at(f)));
    let keel = pts[pts.len() - 1];
    let prev = at(0.9);
    let (dn, dz) = (keel[0] - prev[0], keel[1] - prev[1]);
    let l = dn.hypot(dz).max(1e-9);
    let reach = 0.1 * total;
    out.push([
        keel[0] + reach * dn / l,
        keel[1] + reach * dz.min(-0.2 * l) / l,
    ]);
    // Each a little below the last.
    for i in 1..out.len() {
        out[i][1] = out[i][1].min(out[i - 1][1] - 1e-4);
    }
    out
}

/// A slice's points from the top down (z in the deck frame), and the girth
/// to each.
fn girth(slice: &[V2], z0: f64) -> (Vec<V2>, Vec<f64>) {
    let mut pts: Vec<V2> = slice.iter().map(|p| [p[0], p[1] - z0]).collect();
    pts.sort_by(|a, b| b[1].total_cmp(&a[1]));
    pts.dedup_by(|a, b| (a[0] - b[0]).hypot(a[1] - b[1]) < 1e-9);
    let mut g = vec![0.0];
    for w in pts.windows(2) {
        g.push(g[g.len() - 1] + (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1]));
    }
    (pts, g)
}

/// A hard chine on a slice: where down its girth it turns most sharply,
/// as a fraction of the girth, if it turns more than `CHINE_TURN`.
fn chine_at(slice: &[V2], z0: f64) -> Option<f64> {
    let (pts, g) = girth(slice, z0);
    let total = g[g.len() - 1];
    if total <= 0.0 || pts.len() < 8 {
        return None;
    }
    let w = 0.05 * total;
    let at = |x: f64| {
        let i = g.partition_point(|&y| y < x).min(pts.len() - 1);
        pts[i]
    };
    let mut best: Option<(f64, f64)> = None;
    for (k, &gk) in g.iter().enumerate() {
        if gk < w || gk > total - w {
            continue;
        }
        let (a, b, c) = (at(gk - w), pts[k], at(gk + w));
        let (d1, d2) = ([b[0] - a[0], b[1] - a[1]], [c[0] - b[0], c[1] - b[1]]);
        let turn = (d1[0] * d2[1] - d1[1] * d2[0])
            .atan2(d1[0] * d2[0] + d1[1] * d2[1])
            .abs();
        if best.is_none_or(|b| turn > b.0) {
            best = Some((turn, gk / total));
        }
    }
    best.filter(|b| b.0 > CHINE_TURN.to_radians()).map(|b| b.1)
}

/// The point index a chine gets: the middle one.
fn chine_index(s: usize) -> usize {
    (s - 1) / 2
}

// ---------------------------------------------------------------------------
// The fit
// ---------------------------------------------------------------------------

/// How well a document (in metres) fits a hull, its deck datum at `z0` and
/// its x = 0 at `x0` in the hull's frame: RMS both ways.
pub fn score(geometry: &Value, doc: &Document, z0: f64, x0: f64) -> Result<(f64, f64), String> {
    let target = Target::from_geometry(geometry)?;
    let m = misfit_of(&target, &doc.in_metres(), z0, x0);
    let rms = |v: &[f64]| (v.iter().map(|x| x * x).sum::<f64>() / v.len().max(1) as f64).sqrt();
    Ok((rms(&m.to_hull), rms(&m.to_camber)))
}

/// How many sections of a document stop short of the centreline.
fn open_sections(doc: &Document) -> usize {
    let hull = Hull::new(doc);
    (0..=1000)
        .filter(|&m| {
            hull.column(m as f64 / 1000.0, 1)
                .is_some_and(|c| c.bottom == super::End::Sheet)
        })
        .count()
}

/// The parameters holding each station's last point's n.
fn last_points(frame: &Frame) -> Vec<usize> {
    let s = frame.shape;
    let start = 2 * s.plan - 1 + (s.trim + 1) + s.trim;
    let per = 2 * (s.points - 1);
    (0..s.stations).map(|j| start + j * per + per - 2).collect()
}

/// Levenberg–Marquardt from `q`, the Jacobian by forward differences, until
/// a step gains less than a part in 10⁵: the parameters, and the steps.
fn levenberg_marquardt(target: &Target, frame: &Frame, mut q: Vec<f64>) -> (Vec<f64>, usize) {
    let cost = |r: &[f64]| r.iter().map(|x| x * x).sum::<f64>();
    let mut m = misfit(target, frame, &q);
    let mut c = cost(&m.residuals);
    let mut lambda = 1e-3;
    let mut iterations = 0;
    for it in 0..300 {
        iterations = it + 1;
        let n = q.len();
        let cols: Vec<Vec<f64>> = (0..n)
            .map(|j| {
                let h = 1e-4 * q[j].abs().max(1e-2);
                let mut qp = q.clone();
                qp[j] += h;
                let r = misfit(target, frame, &qp).residuals;
                r.iter()
                    .zip(&m.residuals)
                    .map(|(a, b)| (a - b) / h)
                    .collect()
            })
            .collect();
        let mut jtj = vec![vec![0.0; n]; n];
        let mut jtr = vec![0.0; n];
        for i in 0..n {
            jtr[i] = cols[i].iter().zip(&m.residuals).map(|(a, b)| a * b).sum();
            for j in i..n {
                let v: f64 = cols[i].iter().zip(&cols[j]).map(|(a, b)| a * b).sum();
                jtj[i][j] = v;
                jtj[j][i] = v;
            }
        }
        let mut improved = false;
        for _ in 0..10 {
            let mut a = jtj.clone();
            for i in 0..n {
                a[i][i] += lambda * jtj[i][i].max(1e-12);
            }
            let step = super::Lu::new(a).solve(&jtr.iter().map(|x| -x).collect::<Vec<_>>());
            let qn: Vec<f64> = q.iter().zip(&step).map(|(a, b)| a + b).collect();
            let mn = misfit(target, frame, &qn);
            let cn = cost(&mn.residuals);
            if cn < c {
                let gain = (c - cn) / c;
                q = qn;
                m = mn;
                c = cn;
                lambda = (lambda * 0.3).max(1e-9);
                improved = gain > 1e-8;
                break;
            }
            lambda *= 10.0;
        }
        if !improved {
            break;
        }
    }
    (q, iterations)
}

/// A camber document fitted to a hull's geometry, and how well it fits.
pub fn fit(geometry: &Value, shape: Shape) -> Result<(Value, Report), String> {
    if shape.stations < 2 || shape.points < 3 || shape.plan < 3 || shape.trim < 3 {
        return Err(
            "a fit needs at least 2 stations, 3 points each, and 3 plan and 3 trim points".into(),
        );
    }
    let target = Target::from_geometry(geometry)?;
    // From each start, the best fit.
    let cost = |r: &[f64]| r.iter().map(|x| x * x).sum::<f64>();
    let mut best: Option<(f64, Frame, Vec<f64>, usize)> = None;
    for start in STARTS {
        let (frame, q) = first_guess(&target, shape, start)?;
        let (q, steps) = levenberg_marquardt(&target, &frame, q);
        let c = cost(&misfit(&target, &frame, &q).residuals);
        if std::env::var("FIT_DEBUG").is_ok() {
            eprintln!("debug: start {start:?}: cost {c:.6e} after {steps} steps");
        }
        if best.as_ref().is_none_or(|b| c < b.0) {
            best = Some((c, frame, q, steps));
        }
    }
    let (_, frame, mut q, iterations) = best.unwrap();

    // A section the fit leaves just short of the centreline (near a fine
    // bow, say) is closed by reaching each station's last point further
    // inboard: that point lies past the centreline, in what camber trims
    // away, so the hull hardly moves.
    let step = 1e-4 * (target.x_max - target.x_min);
    let mut reach = 0.0;
    while open_sections(&frame.document(&q)) > 0 {
        if reach > 0.1 * (target.x_max - target.x_min) {
            return Err("the fitted sections don't close on the centreline".into());
        }
        for j in last_points(&frame) {
            q[j] += step;
        }
        reach += step;
    }
    let m = misfit(&target, &frame, &q);
    if std::env::var("FIT_DEBUG").is_ok() {
        let per = m.residuals.len() / PLANES;
        let _ = per;
        let n_cam = ROWS * (shape.points - 1) + 1;
        let stride = n_cam + SAMPLES;
        for k in 0..PLANES {
            let r = &m.residuals[k * stride..(k + 1) * stride];
            let (a, b) = r.split_at(n_cam);
            let mx = |v: &[f64]| v.iter().copied().fold(0.0f64, f64::max);
            eprintln!("debug: plane {k:2} u {:.3}  camber->hull max {:5.1} mm  hull->camber max {:5.1} mm", plane_u(k), 1e3 * mx(a), 1e3 * mx(b));
        }
    }
    let rms = |v: &[f64]| (v.iter().map(|x| x * x).sum::<f64>() / v.len().max(1) as f64).sqrt();
    let max = |v: &[f64]| v.iter().copied().fold(0.0, f64::max);
    let report = Report {
        rms_to_hull: rms(&m.to_hull),
        max_to_hull: max(&m.to_hull),
        rms_to_camber: rms(&m.to_camber),
        max_to_camber: max(&m.to_camber),
        iterations,
        reach,
        x_origin: frame.x0,
    };
    Ok((frame.document(&q).to_json(), report))
}

impl Document {
    /// The same document with its lengths in metres.
    pub fn in_metres(&self) -> Document {
        let k = self.unit_m;
        let m = |p: &V2| [p[0] * k, p[1] * k];
        Document {
            name: self.name.clone(),
            unit_m: 1.0,
            waterline: self.waterline * k,
            deck_trim: self.deck_trim,
            plan: self.plan.iter().map(m).collect(),
            trim: self.trim.iter().map(m).collect(),
            transom: [m(&self.transom[0]), m(&self.transom[1])],
            stations: self
                .stations
                .iter()
                .map(|s| Station {
                    u: s.u,
                    pts: s.pts.iter().map(|p| [p[0] * k, p[1] * k, p[2]]).collect(),
                })
                .collect(),
        }
    }

    /// The document as camber saves it (format version 2).
    pub fn to_json(&self) -> Value {
        let round = |x: f64| (x * 1e9).round() / 1e9;
        json!({
            "version": 2,
            "name": self.name,
            "unit": "m",
            "waterline": round(self.waterline * self.unit_m),
            "deckTrimDeg": self.deck_trim.to_degrees(),
            "sheerPlan": self.plan.iter().map(|p| json!({ "x": round(p[0] * self.unit_m), "y": round(p[1] * self.unit_m) })).collect::<Vec<_>>(),
            "sheerTrim": self.trim.iter().map(|p| json!({ "x": round(p[0] * self.unit_m), "z": round(p[1] * self.unit_m), "k": 0 })).collect::<Vec<_>>(),
            "transom": self.transom.iter().map(|p| json!({ "x": round(p[0] * self.unit_m), "z": round(p[1] * self.unit_m) })).collect::<Vec<_>>(),
            "stations": self.stations.iter().map(|s| json!({
                "u": s.u,
                "keelK": 0,
                "points": s.pts.iter().map(|p| json!({ "n": round(p[0] * self.unit_m), "z": round(p[1] * self.unit_m), "k": p[2] })).collect::<Vec<_>>(),
            })).collect::<Vec<_>>(),
        })
    }
}

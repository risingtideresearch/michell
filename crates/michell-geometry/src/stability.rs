//! Large-angle transverse stability: the righting-arm (GZ) curve of a
//! platform of hulls, from their whole sections.
//!
//! The hydrodynamic sections ([`crate::SectionalHull`]) stop at the
//! waterline; a heeled hull immerses what was above it and lifts out what
//! was below. So this works on each station's **whole** outline, keel to
//! sheer, closed across the top (a flush deck at the sheer — deck openings
//! and downflooding are not modelled), and on both sides of its centreplane.
//!
//! ## Method
//!
//! Heel `φ` is a rotation about the longitudinal axis, trim `θ` one about
//! the transverse axis, and sinkage `s` a drop; together a rigid map, under
//! which a body point's height above the still water is
//!
//! ```text
//! Z = x sinθ + cosθ (y sinφ + z cosφ) − s        (x fore, y port, z up)
//! ```
//!
//! — at a station (`x` fixed) a straight line in its `(y, z)` plane. Each
//! section polygon is clipped by the half-plane `Z ≤ 0` (Sutherland–Hodgman,
//! exact at the waterline, so a wall-sided box comes out exact), and the
//! clipped polygon's area and centroid (shoelace) are integrated along `x`
//! (trapezoidal between stations) for the immersed volume and the centre of
//! buoyancy `B` in body axes. With the centre of gravity `G`, the righting
//! arm is the horizontal transverse distance between them in the water's
//! frame,
//!
//! ```text
//! GZ = Y(G) − Y(B),        Y(p) = y cosφ − z sinφ,
//! ```
//!
//! positive when the couple rights the platform. At each heel the platform
//! is balanced first: its sinkage carries the mass, and with free trim its
//! trim puts `B` under `G` fore and aft as well (a Newton solve on `(s, θ)`,
//! from the previous heel's answer).
//!
//! Positive `φ` lifts the `+y` (port) side: the platform heels to
//! starboard. The curve is odd in `φ` for a platform symmetric about `y = 0`.

use crate::error::{Error, Result};

/// One station's whole section: `(half-beam ≥ 0, height above the design
/// waterline)` from the top of its side round to the keel, on one side of
/// the hull's centreplane.
#[derive(Debug, Clone)]
pub struct FullSection {
    pub x: f64,
    pub half: Vec<(f64, f64)>,
}

/// A hull of the platform: its whole sections, and where its centreplane
/// is (`y`) — the sections' `x` are already in the platform's frame.
#[derive(Debug, Clone)]
pub struct StabilityHull {
    pub sections: Vec<FullSection>,
    pub y: f64,
}

/// What the platform carries: its mass and where its centre of gravity is.
#[derive(Debug, Clone, Copy)]
pub struct StabilityLoad {
    /// [kg]
    pub mass: f64,
    /// Longitudinal centre of gravity [m, platform x].
    pub lcg: f64,
    /// Height of the centre of gravity above the design waterline [m].
    pub vcg: f64,
}

/// Options for [`gz_curve`].
#[derive(Debug, Clone)]
pub struct GzOptions {
    /// The heel angles [rad], each solved in turn (ascending from 0 is
    /// fastest: each starts from the last).
    pub heels: Vec<f64>,
    /// Let the platform trim at each heel (else it keeps its upright trim).
    pub free_trim: bool,
    pub density: f64,
    /// A wave frozen around the platform, in place of calm water (see
    /// [`Wave`]); `None` in calm water.
    pub wave: Option<Wave>,
}

/// A regular wave held still around the platform: the quasi-static GZ in
/// waves, where the water's surface is the wave's rather than a plane. The
/// elevation (up, from the mean water) is
///
/// ```text
/// ζ(X, Y) = ½H cos(k (X cos μ + Y sin μ) − ψ),    k = 2π/λ,
/// ```
///
/// in the water's axes (X fore, Y port), for a heading `μ` — 90° beam
/// seas, crests along the platform; 45° stern quartering; 135° bow
/// quartering; 0° or 180° crests across it — and a crest on the line
/// `k (X cos μ + Y sin μ) = ψ`. The wave's orbital pressures (the Smith
/// effect) and its dynamics are not modelled: the hull floats on the
/// frozen surface as it would on still water shaped so.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wave {
    /// Wavelength λ [m].
    pub length: f64,
    /// Crest-to-trough height H [m].
    pub height: f64,
    /// μ [rad].
    pub heading: f64,
    /// ψ [rad].
    pub phase: f64,
}

impl Wave {
    pub fn elevation(&self, x: f64, y: f64) -> f64 {
        let k = 2.0 * std::f64::consts::PI / self.length;
        0.5 * self.height
            * (k * (x * self.heading.cos() + y * self.heading.sin()) - self.phase).cos()
    }

    /// The phase that puts a crest through `(x, y)`.
    pub fn crest_at(length: f64, heading: f64, x: f64, y: f64) -> f64 {
        2.0 * std::f64::consts::PI / length * (x * heading.cos() + y * heading.sin())
    }
}

impl Default for GzOptions {
    /// 0° to 180° in 2° steps, free trim, sea water.
    fn default() -> Self {
        GzOptions {
            heels: (0..=90).map(|i| (2.0 * i as f64).to_radians()).collect(),
            free_trim: true,
            density: crate::Fluid::SEAWATER_15C.density,
            wave: None,
        }
    }
}

/// The balanced platform at one heel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GzPoint {
    /// [rad]
    pub heel: f64,
    /// Righting arm [m].
    pub gz: f64,
    /// Drop of the design waterline's origin below the water [m].
    pub sinkage: f64,
    /// Bow up [rad].
    pub trim: f64,
    /// Immersed volume [m³] (the mass over the density, when balanced).
    pub volume: f64,
    /// Centre of buoyancy, body axes [m].
    pub cb: [f64; 3],
}

/// A GZ curve with the figures read from it.
#[derive(Debug, Clone, PartialEq)]
pub struct GzCurve {
    pub points: Vec<GzPoint>,
    /// The curve's slope at upright, `dGZ/dφ` [m/rad] — the metacentric
    /// height the large-angle integration implies.
    pub gm: f64,
    /// The largest righting arm [m] and its heel [rad].
    pub max_gz: f64,
    pub heel_at_max: f64,
    /// The first heel past upright where the arm turns negative (the angle
    /// of vanishing stability) [rad], interpolated; `None` if it stays
    /// positive over the heels asked for. Zero when the platform does not
    /// right itself at all from upright (GM ≤ 0: it lolls, or capsizes) —
    /// any later crossing is past an angle of loll, not a range of
    /// stability.
    pub vanishing: Option<f64>,
    /// `∫ GZ dφ` from upright to 30°, to 40° and to the vanishing angle
    /// (or the last heel) [m·rad], over the heels asked for.
    pub area_30: f64,
    pub area_40: f64,
    pub area_total: f64,
}

/// A section polygon from its half outline: down the `−y` side from the
/// top to the keel, back up the `+y` side, shut across the top. (Either
/// orientation integrates alike: [`clipped`] divides it out.)
fn polygon(half: &[(f64, f64)], y0: f64) -> Vec<(f64, f64)> {
    let mut p: Vec<(f64, f64)> = Vec::with_capacity(2 * half.len());
    for &(h, z) in half {
        p.push((y0 - h, z));
    }
    for &(h, z) in half.iter().rev() {
        p.push((y0 + h, z));
    }
    p.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12);
    p
}

/// Area and centroid of the part of `poly` where `a y + b z + c ≤ 0`.
fn clipped(poly: &[(f64, f64)], a: f64, b: f64, c: f64) -> (f64, f64, f64) {
    let n = poly.len();
    if n < 3 {
        return (0.0, 0.0, 0.0);
    }
    let f = |p: (f64, f64)| a * p.0 + b * p.1 + c;
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(n + 2);
    for i in 0..n {
        let (p, q) = (poly[i], poly[(i + 1) % n]);
        let (fp, fq) = (f(p), f(q));
        if fp <= 0.0 {
            out.push(p);
        }
        if (fp <= 0.0) != (fq <= 0.0) {
            let t = fp / (fp - fq);
            out.push((p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1)));
        }
    }
    let m = out.len();
    if m < 3 {
        return (0.0, 0.0, 0.0);
    }
    let (mut area2, mut cy, mut cz) = (0.0, 0.0, 0.0);
    for i in 0..m {
        let (p, q) = (out[i], out[(i + 1) % m]);
        let cross = p.0 * q.1 - q.0 * p.1;
        area2 += cross;
        cy += (p.0 + q.0) * cross;
        cz += (p.1 + q.1) * cross;
    }
    if area2.abs() < 1e-300 {
        return (0.0, 0.0, 0.0);
    }
    // The polygon's orientation sets the sign of the shoelace sum.
    let area = 0.5 * area2.abs();
    (area, cy / (3.0 * area2), cz / (3.0 * area2))
}

/// The part of `poly` (at station `x`) under a wave's surface, the platform
/// at `(heel, trim, sinkage)`: the polygon clipped by `Z − ζ(X, Y) ≤ 0`,
/// its crossings found linearly along each edge. The polygon is densified
/// beforehand (see [`Prepared`]) so that the surface is near straight along
/// an edge; between crossings the waterline is the chord, within
/// `H k² w²/16` of the wave over a width `w` of section.
fn clipped_wave(
    poly: &[(f64, f64)],
    x: f64,
    heel: f64,
    trim: f64,
    sinkage: f64,
    wave: &Wave,
) -> (f64, f64, f64) {
    let (cp, sp, ct, st) = (heel.cos(), heel.sin(), trim.cos(), trim.sin());
    let f = |(y, z): (f64, f64)| {
        let zr = y * sp + z * cp;
        let (xe, ye, ze) = (
            x * ct - st * zr,
            y * cp - z * sp,
            x * st + ct * zr - sinkage,
        );
        ze - wave.elevation(xe, ye)
    };
    clip_by(poly, f)
}

/// Area and centroid of the part of `poly` where `f ≤ 0`, `f` linear
/// enough along each edge to find its crossing there.
fn clip_by(poly: &[(f64, f64)], f: impl Fn((f64, f64)) -> f64) -> (f64, f64, f64) {
    let n = poly.len();
    if n < 3 {
        return (0.0, 0.0, 0.0);
    }
    let fs: Vec<f64> = poly.iter().map(|&q| f(q)).collect();
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(n + 4);
    for i in 0..n {
        let j = (i + 1) % n;
        let (p, q, fp, fq) = (poly[i], poly[j], fs[i], fs[j]);
        if fp <= 0.0 {
            out.push(p);
        }
        if (fp <= 0.0) != (fq <= 0.0) {
            let t = fp / (fp - fq);
            out.push((p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1)));
        }
    }
    shoelace(&out)
}

fn shoelace(out: &[(f64, f64)]) -> (f64, f64, f64) {
    let m = out.len();
    if m < 3 {
        return (0.0, 0.0, 0.0);
    }
    let (mut area2, mut cy, mut cz) = (0.0, 0.0, 0.0);
    for i in 0..m {
        let (p, q) = (out[i], out[(i + 1) % m]);
        let cross = p.0 * q.1 - q.0 * p.1;
        area2 += cross;
        cy += (p.0 + q.0) * cross;
        cz += (p.1 + q.1) * cross;
    }
    if area2.abs() < 1e-300 {
        return (0.0, 0.0, 0.0);
    }
    let area = 0.5 * area2.abs();
    (area, cy / (3.0 * area2), cz / (3.0 * area2))
}

/// `poly` with no edge longer than `d`.
fn densified(poly: &[(f64, f64)], d: f64) -> Vec<(f64, f64)> {
    let n = poly.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let (p, q) = (poly[i], poly[(i + 1) % n]);
        let len = (q.0 - p.0).hypot(q.1 - p.1);
        let k = (len / d).ceil().max(1.0) as usize;
        for j in 0..k {
            let t = j as f64 / k as f64;
            out.push((p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1)));
        }
    }
    out
}

/// A station: its x and its section polygon.
type Station = (f64, Vec<(f64, f64)>);

/// Section polygons of the platform, per hull per station, with their x —
/// densified for a wave's curved surface when there is one.
struct Prepared {
    hulls: Vec<Vec<Station>>,
    wave: Option<Wave>,
}

impl Prepared {
    fn new(hulls: &[StabilityHull], wave: Option<Wave>) -> Prepared {
        Prepared {
            hulls: hulls
                .iter()
                .map(|h| {
                    let mut s: Vec<(f64, Vec<(f64, f64)>)> = h
                        .sections
                        .iter()
                        .filter(|s| s.half.len() >= 2)
                        .map(|s| {
                            let p = polygon(&s.half, h.y);
                            let p = match wave {
                                Some(w) => densified(&p, w.length / 32.0),
                                None => p,
                            };
                            (s.x, p)
                        })
                        .collect();
                    s.sort_by(|a, b| a.0.total_cmp(&b.0));
                    s
                })
                .collect(),
            wave,
        }
    }

    /// Immersed volume and centre of buoyancy (body axes) at an attitude.
    fn immersed(&self, heel: f64, trim: f64, sinkage: f64) -> (f64, [f64; 3]) {
        let (a, b) = (trim.cos() * heel.sin(), trim.cos() * heel.cos());
        let (mut v, mut mx, mut my, mut mz) = (0.0, 0.0, 0.0, 0.0);
        for stations in &self.hulls {
            let cut: Vec<(f64, f64, f64, f64)> = stations
                .iter()
                .map(|(x, poly)| {
                    let (area, cy, cz) = match &self.wave {
                        Some(w) => clipped_wave(poly, *x, heel, trim, sinkage, w),
                        None => clipped(poly, a, b, x * trim.sin() - sinkage),
                    };
                    (*x, area, cy, cz)
                })
                .collect();
            for w in cut.windows(2) {
                let (x0, a0, y0, z0) = w[0];
                let (x1, a1, y1, z1) = w[1];
                let dx = x1 - x0;
                // Trapezoidal in x for the area and each first moment.
                v += 0.5 * dx * (a0 + a1);
                mx += 0.5 * dx * (a0 * x0 + a1 * x1);
                my += 0.5 * dx * (a0 * y0 + a1 * y1);
                mz += 0.5 * dx * (a0 * z0 + a1 * z1);
            }
        }
        if v <= 0.0 {
            return (0.0, [0.0; 3]);
        }
        (v, [mx / v, my / v, mz / v])
    }
}

/// The residuals of the balance at `(s, θ)`: volume short of the load's,
/// and (for free trim) how far B lies ahead of G along the water.
fn residuals(
    p: &Prepared,
    load: &StabilityLoad,
    target: f64,
    heel: f64,
    s: f64,
    t: f64,
) -> (f64, f64, f64, [f64; 3]) {
    let (v, cb) = p.immersed(heel, t, s);
    let along = |q: [f64; 3]| q[0] * t.cos() - t.sin() * (q[1] * heel.sin() + q[2] * heel.cos());
    let dx = if v > 0.0 {
        along(cb) - along([load.lcg, 0.0, load.vcg])
    } else {
        0.0
    };
    (v - target, dx, v, cb)
}

/// Balance the platform at one heel from `(s, t)`: the sinkage (and, for
/// free trim, the trim) that floats the load.
#[allow(clippy::too_many_arguments)]
fn balance(
    p: &Prepared,
    load: &StabilityLoad,
    target: f64,
    heel: f64,
    mut s: f64,
    mut t: f64,
    free_trim: bool,
    scale: (f64, f64),
) -> Result<(f64, f64, f64, [f64; 3])> {
    let (ls, lx) = scale; // a length for sinkage steps, and for trim moments
    let tol_v = 1e-9 * target;
    let tol_x = 1e-7 * lx;
    for _ in 0..60 {
        let (rv, rx, v, cb) = residuals(p, load, target, heel, s, t);
        if rv.abs() <= tol_v && (!free_trim || rx.abs() <= tol_x) && v > 0.0 {
            return Ok((s, t, v, cb));
        }
        // A dry platform: sink it until it floats.
        if v <= 0.0 {
            s += 0.25 * ls;
            continue;
        }
        // Central differences: at a deck or keel exactly at the water a
        // one-sided slope can be zero (nothing more to immerse on that side).
        let hs = 1e-6 * ls;
        let (rv_p, rx_p, _, _) = residuals(p, load, target, heel, s + hs, t);
        let (rv_m, rx_m, _, _) = residuals(p, load, target, heel, s - hs, t);
        let (dvs, dxs) = ((rv_p - rv_m) / (2.0 * hs), (rx_p - rx_m) / (2.0 * hs));
        let (ds, dt) = if free_trim {
            let ht = 1e-6;
            let (rv_p, rx_p, _, _) = residuals(p, load, target, heel, s, t + ht);
            let (rv_m, rx_m, _, _) = residuals(p, load, target, heel, s, t - ht);
            let (dvt, dxt) = ((rv_p - rv_m) / (2.0 * ht), (rx_p - rx_m) / (2.0 * ht));
            let det = dvs * dxt - dvt * dxs;
            if det.abs() < 1e-300 {
                return Err(Error::InvalidInput(format!(
                    "GZ at {:.1}°: no attitude floats the load — the platform is fully \
                     immersed (no freeboard to sink into) or clear of the water",
                    heel.to_degrees()
                )));
            }
            ((-rv * dxt + rx * dvt) / det, (-dvs * rx + dxs * rv) / det)
        } else {
            if dvs.abs() < 1e-300 {
                return Err(Error::InvalidInput(format!(
                    "GZ at {:.1}°: the waterplane has vanished",
                    heel.to_degrees()
                )));
            }
            (-rv / dvs, 0.0)
        };
        // Damped: no step past half the depth scale, or 10° of trim.
        let k = (0.5 * ls / ds.abs().max(1e-300))
            .min(0.17 / dt.abs().max(1e-300))
            .min(1.0);
        s += k * ds;
        t += k * dt;
    }
    Err(Error::InvalidInput(format!(
        "GZ at {:.1}°: the balance did not converge",
        heel.to_degrees()
    )))
}

/// The platform's GZ curve under `load`, over `opts.heels`.
pub fn gz_curve(
    hulls: &[StabilityHull],
    load: &StabilityLoad,
    opts: &GzOptions,
) -> Result<GzCurve> {
    if load.mass.is_nan() || load.mass <= 0.0 {
        return Err(Error::InvalidInput("GZ: the mass must be positive".into()));
    }
    let p = Prepared::new(hulls, opts.wave);
    if p.hulls.iter().all(|h| h.len() < 2) {
        return Err(Error::InvalidInput("GZ: no sections".into()));
    }
    // Length scales for the solve: the platform's height and length.
    let (mut zlo, mut zhi, mut xlo, mut xhi) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for h in hulls {
        for s in &h.sections {
            xlo = xlo.min(s.x);
            xhi = xhi.max(s.x);
            for &(_, z) in &s.half {
                zlo = zlo.min(z);
                zhi = zhi.max(z);
            }
        }
    }
    let scale = ((zhi - zlo).max(1e-6), (xhi - xlo).max(1e-6));
    let target = load.mass / opts.density;
    let (mut s, mut t) = (0.0, 0.0);
    // Upright first: its trim is the one a fixed-trim curve keeps.
    let (s0, t0, _, _) = balance(&p, load, target, 0.0, s, t, true, scale)?;
    let upright_trim = t0;
    s = s0;
    t = t0;
    let mut points = Vec::with_capacity(opts.heels.len());
    for &heel in &opts.heels {
        let t_start = if opts.free_trim { t } else { upright_trim };
        let (s1, t1, v, cb) = balance(&p, load, target, heel, s, t_start, opts.free_trim, scale)?;
        s = s1;
        t = t1;
        let across = |q: [f64; 3]| q[1] * heel.cos() - q[2] * heel.sin();
        let gz = across([load.lcg, 0.0, load.vcg]) - across(cb);
        points.push(GzPoint {
            heel,
            gz,
            sinkage: s1,
            trim: t1,
            volume: v,
            cb,
        });
    }
    Ok(summarise(
        points,
        hulls,
        load,
        opts,
        &p,
        target,
        scale,
        upright_trim,
    ))
}

#[allow(clippy::too_many_arguments)]
fn summarise(
    points: Vec<GzPoint>,
    _hulls: &[StabilityHull],
    load: &StabilityLoad,
    opts: &GzOptions,
    p: &Prepared,
    target: f64,
    scale: (f64, f64),
    upright_trim: f64,
) -> GzCurve {
    // The slope at upright from a small heel either side of it (central,
    // so a platform asymmetric in y still reads right).
    let h = 0.5f64.to_radians();
    let arm = |heel: f64| {
        balance(
            p,
            load,
            target,
            heel,
            0.0,
            upright_trim,
            opts.free_trim,
            scale,
        )
        .map(|(_, _, _, cb)| {
            let across = |q: [f64; 3]| q[1] * heel.cos() - q[2] * heel.sin();
            across([load.lcg, 0.0, load.vcg]) - across(cb)
        })
        .unwrap_or(f64::NAN)
    };
    let gm = (arm(h) - arm(-h)) / (2.0 * h);
    let (mut max_gz, mut heel_at_max) = (f64::NEG_INFINITY, 0.0);
    for q in &points {
        if q.gz > max_gz {
            max_gz = q.gz;
            heel_at_max = q.heel;
        }
    }
    let mut vanishing = None;
    // Unstable upright: no range of positive stability to vanish at.
    let upright_unstable = gm.is_nan() || gm <= 0.0
        || points
            .iter()
            .find(|q| q.heel > 0.0)
            .is_some_and(|q| q.gz <= 0.0);
    for w in points.windows(2).filter(|_| !upright_unstable) {
        if w[0].heel > 0.0 && w[0].gz > 0.0 && w[1].gz <= 0.0 {
            let f = w[0].gz / (w[0].gz - w[1].gz);
            vanishing = Some(w[0].heel + f * (w[1].heel - w[0].heel));
            break;
        }
    }
    // ∫ GZ dφ from 0 to `upto`, trapezoidal over the points (linear within
    // the last interval).
    let area = |upto: f64| {
        let mut a = 0.0;
        for w in points.windows(2) {
            let (h0, h1) = (w[0].heel.max(0.0), w[1].heel.min(upto));
            if h1 <= h0 || w[1].heel <= w[0].heel {
                continue;
            }
            let at =
                |x: f64| w[0].gz + (w[1].gz - w[0].gz) * (x - w[0].heel) / (w[1].heel - w[0].heel);
            a += 0.5 * (h1 - h0) * (at(h0) + at(h1));
        }
        a
    };
    let last = points.last().map_or(0.0, |q| q.heel);
    GzCurve {
        area_30: area(30f64.to_radians()),
        area_40: area(40f64.to_radians()),
        area_total: area(vanishing.unwrap_or(last)),
        points,
        gm,
        max_gz: if max_gz.is_finite() { max_gz } else { 0.0 },
        heel_at_max,
        vanishing: if upright_unstable {
            Some(0.0)
        } else {
            vanishing
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A box hull, `l × b × depth`, its keel `draft` below the design
    /// waterline, cut into `n` stations.
    fn boxy(l: f64, b: f64, depth: f64, draft: f64, y: f64, n: usize) -> StabilityHull {
        let half = vec![(0.5 * b, depth - draft), (0.5 * b, -draft), (0.0, -draft)];
        StabilityHull {
            sections: (0..n)
                .map(|i| FullSection {
                    x: l * i as f64 / (n - 1) as f64,
                    half: half.clone(),
                })
                .collect(),
            y,
        }
    }

    fn opts(deg: &[f64]) -> GzOptions {
        GzOptions {
            heels: deg.iter().map(|d| d.to_radians()).collect(),
            ..GzOptions::default()
        }
    }

    /// A wall-sided box, heeled short of immersing its deck or lifting its
    /// bilge: GZ = sin φ (GM + ½ BM tan² φ) exactly.
    #[test]
    fn a_box_follows_the_wall_sided_formula() {
        let (l, b, d, t) = (10.0, 2.0, 1.5, 0.5);
        let vcg = 0.2; // KG = t + vcg
        let rho = GzOptions::default().density;
        let load = StabilityLoad {
            mass: rho * l * b * t,
            lcg: 0.5 * l,
            vcg,
        };
        let hull = boxy(l, b, d, t, 0.0, 11);
        let heels = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0];
        let c = gz_curve(&[hull], &load, &opts(&heels)).unwrap();
        let bm = b * b / (12.0 * t);
        let gm = 0.5 * t + bm - (t + vcg);
        for q in &c.points {
            let f = q.heel;
            let exact = f.sin() * (gm + 0.5 * bm * f.tan().powi(2));
            assert!(
                (q.gz - exact).abs() < 1e-7,
                "{:.0}°: {} vs {exact}",
                f.to_degrees(),
                q.gz
            );
            assert!(
                q.trim.abs() < 1e-9,
                "a box symmetric fore and aft stays level"
            );
        }
        assert!((c.gm - gm).abs() < 1e-4, "{} vs {gm}", c.gm);
        // Deck edge in at atan(2(d−t)/b) = 45°; the arm keeps rising a while.
        assert!(c.max_gz > 0.0);
    }

    /// A box with G above its metacentre does not right itself: its range
    /// of stability is nothing, however its curve runs later.
    #[test]
    fn an_unstable_box_has_no_range_of_stability() {
        let (l, b, d, t) = (10.0, 2.0, 1.5, 0.5);
        let rho = GzOptions::default().density;
        // KM = T/2 + B²/12T = 0.917 m; G at 0.5 + 0.6 = 1.1 m above the keel.
        let load = StabilityLoad {
            mass: rho * l * b * t,
            lcg: 0.5 * l,
            vcg: 0.6,
        };
        let c = gz_curve(
            &[boxy(l, b, d, t, 0.0, 11)],
            &load,
            &opts(&[0.0, 2.0, 10.0, 30.0, 60.0]),
        )
        .unwrap();
        assert!(c.gm < 0.0, "{}", c.gm);
        assert_eq!(c.vanishing, Some(0.0));
    }

    /// A wave too low to matter is calm water; so is a wave so long that
    /// the platform sits on its crest, which only lifts the water.
    #[test]
    fn a_negligible_wave_is_calm_water() {
        let hull = boxy(10.0, 2.0, 1.5, 0.5, 0.0, 11);
        let rho = GzOptions::default().density;
        let load = StabilityLoad {
            mass: rho * 10.0,
            lcg: 5.0,
            vcg: 0.1,
        };
        let heels = [0.0, 10.0, 30.0, 60.0];
        let calm = gz_curve(std::slice::from_ref(&hull), &load, &opts(&heels)).unwrap();
        for wave in [
            Wave {
                length: 20.0,
                height: 1e-9,
                heading: 0.5 * std::f64::consts::PI,
                phase: 0.3,
            },
            // λ 10 km, the crest over the platform: the water raised 0.2 m.
            Wave {
                length: 1e4,
                height: 0.4,
                heading: 0.5 * std::f64::consts::PI,
                phase: 0.0,
            },
        ] {
            let c = gz_curve(
                std::slice::from_ref(&hull),
                &load,
                &GzOptions {
                    wave: Some(wave),
                    ..opts(&heels)
                },
            )
            .unwrap();
            for (a, b) in calm.points.iter().zip(&c.points) {
                assert!(
                    (a.gz - b.gz).abs() < 2e-5,
                    "{wave:?} at {:.0}°: {} vs {}",
                    a.heel.to_degrees(),
                    b.gz,
                    a.gz
                );
            }
        }
    }

    /// A catamaran in beam seas, λ twice its span: with a crest under the
    /// windward (port, lifting) hull and a trough under the lee, the water
    /// slopes down to leeward — the way the platform heels — so the
    /// platform heels less against the water than against the vertical,
    /// both hulls stay in, and the arm is far smaller than in calm water;
    /// with the crest under the lee hull the slope runs the other way and
    /// the windward hull lifts sooner.
    #[test]
    fn a_beam_sea_crest_moves_the_arm() {
        let (l, b, d, t, span) = (10.0, 0.8, 1.2, 0.4, 4.0);
        let rho = GzOptions::default().density;
        let load = StabilityLoad {
            mass: rho * 2.0 * l * b * t,
            lcg: 0.5 * l,
            vcg: 1.0,
        };
        let hulls = [
            boxy(l, b, d, t, 0.5 * span, 11),
            boxy(l, b, d, t, -0.5 * span, 11),
        ];
        let heels = [0.0, 5.0, 15.0, 25.0];
        let calm = gz_curve(&hulls, &load, &opts(&heels)).unwrap();
        let beam = 0.5 * std::f64::consts::PI;
        let (lambda, h) = (2.0 * span, 0.8);
        let at = |y: f64| Wave {
            length: lambda,
            height: h,
            heading: beam,
            phase: Wave::crest_at(lambda, beam, 0.0, y),
        };
        let windward = gz_curve(
            &hulls,
            &load,
            &GzOptions {
                wave: Some(at(0.5 * span)),
                ..opts(&heels)
            },
        )
        .unwrap();
        let lee = gz_curve(
            &hulls,
            &load,
            &GzOptions {
                wave: Some(at(-0.5 * span)),
                ..opts(&heels)
            },
        )
        .unwrap();
        let k = 2; // 15°: the calm-water windward hull is clear.
        assert!(
            windward.points[k].gz < 0.5 * calm.points[k].gz,
            "{} vs {}",
            windward.points[k].gz,
            calm.points[k].gz
        );
        // With the crest to leeward the windward hull lifts sooner: at 5°,
        // with both hulls still in in calm water, a larger arm; by 15°, the
        // lee hull carrying all either way, about the same.
        assert!(
            lee.points[1].gz > calm.points[1].gz,
            "{} vs {}",
            lee.points[1].gz,
            calm.points[1].gz
        );
        assert!((lee.points[k].gz - calm.points[k].gz).abs() < 0.02 * calm.points[k].gz);
        // Upright, the slope to leeward heels it: a negative arm.
        assert!(windward.points[0].gz < 0.0 && lee.points[0].gz > 0.0);
        for c in [&windward, &lee] {
            for q in &c.points {
                assert!((q.volume * rho - load.mass).abs() < 1e-6 * load.mass);
            }
        }
    }

    /// Heeling the other way mirrors the curve.
    #[test]
    fn the_curve_is_odd_in_heel() {
        let hull = boxy(8.0, 1.2, 1.0, 0.4, 0.0, 9);
        let load = StabilityLoad {
            mass: 3500.0,
            lcg: 3.5,
            vcg: 0.1,
        };
        let c = gz_curve(
            &[hull],
            &load,
            &opts(&[0.0, 20.0, 40.0, 60.0, -20.0, -40.0, -60.0]),
        )
        .unwrap();
        for k in 1..4 {
            let (a, b) = (c.points[k].gz, c.points[k + 3].gz);
            assert!((a + b).abs() < 1e-6 * a.abs().max(1e-3), "{a} {b}");
        }
    }

    /// Two boxes as a catamaran: at small heel GM = I_T/∇ − BG with
    /// I_T = 2 (l b³/12 + l b (s/2)²); later a hull flies and the arm peaks.
    #[test]
    fn a_catamaran_of_boxes() {
        let (l, b, d, t, span) = (10.0, 0.8, 1.2, 0.4, 4.0);
        let vcg = 1.0;
        let rho = GzOptions::default().density;
        let vol = 2.0 * l * b * t;
        let load = StabilityLoad {
            mass: rho * vol,
            lcg: 0.5 * l,
            vcg,
        };
        let hulls = [
            boxy(l, b, d, t, 0.5 * span, 11),
            boxy(l, b, d, t, -0.5 * span, 11),
        ];
        let heels: Vec<f64> = (0..=45).map(|i| 2.0 * i as f64).collect();
        let c = gz_curve(&hulls, &load, &opts(&heels)).unwrap();
        let it = 2.0 * (l * b.powi(3) / 12.0 + l * b * (0.5 * span).powi(2));
        let gm = it / vol - (t + vcg - 0.5 * t);
        assert!((c.gm - gm).abs() < 1e-3 * gm, "{} vs {gm}", c.gm);
        // The windward hull lifts clear near atan(2t/span) ≈ 11°: the arm
        // peaks near there, then falls as G swings over the lee hull.
        let peak = c.heel_at_max.to_degrees();
        assert!(peak > 6.0 && peak < 20.0, "{peak}");
        assert!(c.vanishing.is_some());
        assert!(c.area_30 > 0.0 && c.area_40 >= c.area_30);
    }

    /// A load aft of the hull's middle trims it by the stern, upright and
    /// heeled; with fixed trim the heeled curve keeps the upright trim.
    #[test]
    fn free_trim_follows_the_load() {
        let hull = boxy(10.0, 2.0, 1.5, 0.5, 0.0, 11);
        let rho = GzOptions::default().density;
        let load = StabilityLoad {
            mass: rho * 10.0,
            lcg: 4.5,
            vcg: 0.0,
        };
        let free = gz_curve(std::slice::from_ref(&hull), &load, &opts(&[0.0, 30.0])).unwrap();
        assert!(free.points[0].trim > 0.0, "bow up: {}", free.points[0].trim);
        let fixed = gz_curve(
            &[hull],
            &load,
            &GzOptions {
                free_trim: false,
                ..opts(&[0.0, 30.0])
            },
        )
        .unwrap();
        assert!((fixed.points[1].trim - fixed.points[0].trim).abs() < 1e-12);
        for q in &free.points {
            assert!((q.volume * rho - load.mass).abs() < 1e-6 * load.mass);
        }
    }
}

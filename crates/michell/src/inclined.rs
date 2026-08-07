//! Inclined-waterplane hydrostatics.
//!
//! The [`Hull`](crate::Hull) integrals (`displaced_volume`, `vcb_z`,
//! `waterplane_transverse_moment`, …) are thin-ship forms that assume a
//! **horizontal** cut at `z = 0` and a section symmetric about the centreplane
//! (`|y| ≤ f`). They cannot describe a hull heeled relative to the free
//! surface: the keel swings off the earth-vertical centreplane, so the immersed
//! part of each section is the region under a **tilted** waterline, an
//! asymmetric shape the closed forms do not cover.
//!
//! This module integrates that true cut directly on a full-band [`Body`], for
//! an arbitrary rigid attitude — sinkage, trim, **and heel**. It is the
//! geometry behind an inclined-waterplane righting arm (a genuinely nonlinear
//! GZ curve, form stability included), replacing the metacentric approximation
//! `sinφ·(I_T/∇ − KB)` that [`crate::float::righting_arm`] adds for the
//! hull-local heel a half-breadth surface cannot rotate.
//!
//! ## Method
//!
//! Heel is a rotation about the longitudinal axis, so it leaves `x` untouched
//! and only tilts each section in its own `(y, z)` plane. The whole
//! body→earth map (design pose, heel, platform trim, sinkage) is a rigid
//! transform, hence **affine** in the body coordinates `(x_b, z_b, y_b)`, so
//! the earth-frame depth below the still-water surface is
//!
//! ```text
//! D(x_b, z_b, y_b) = a·x_b + b·z_b + c·y_b + d.
//! ```
//!
//! At each station `x_b` the hull section is the polygon `|y_b| ≤ f(x_b, z_b)`,
//! `z_b ∈ [0, depth]`; the immersed part is that polygon **clipped** by the
//! half-plane `D ≥ 0` — a straight line in `(y_b, z_b)`. Clipping
//! (Sutherland–Hodgman) intersects the waterline with the section edges
//! *exactly*, so there is no step/kink error at the free surface (the volume of
//! a wall-sided box comes out exact). The clipped polygon's area and centroid
//! (shoelace) give the station's immersed area `A(x_b)` and body-frame
//! centroid; the volume is `∫ A dx_b` (the rigid map preserves volume) and the
//! earth-frame centre of buoyancy weights each station's centroid, mapped
//! through the affine transform, by `A`. Exact for the spline geometry up to
//! the polygonal sampling of the section outline (exact for straight frames)
//! and the station quadrature.
//!
//! A symmetric body (`|y| ≤ f`) is assumed; the immersed **asymmetry** comes
//! from the tilt, not the hull. (Port/starboard bands for cambered hulls are
//! future work, alongside the asymmetric-hull front-end.)

use crate::body::Body;
use crate::iges::{HullPose, Platform};

/// Integration resolution: body stations × band samples (a tensor grid over
/// `x_b ∈ [x0, x1]` and `z_b ∈ [0, depth]`, trapezoidal in both directions).
#[derive(Debug, Clone, Copy)]
pub struct InclinedGrid {
    pub stations: usize,
    pub band: usize,
}

impl Default for InclinedGrid {
    fn default() -> Self {
        InclinedGrid {
            stations: 161,
            band: 81,
        }
    }
}

/// Immersed hydrostatics of one body at one attitude, in the earth frame.
/// `z` is depth below the still-water surface (positive **down**).
#[derive(Debug, Clone, Copy)]
pub struct InclinedHydro {
    /// Displaced volume `∇` [m³].
    pub volume: f64,
    /// Longitudinal centre of buoyancy (earth `x`) [m].
    pub lcb: f64,
    /// Transverse centre of buoyancy (earth `y`, measured from the heel axis on
    /// the platform centerline) [m]. This is the buoyancy lever the righting
    /// arm is built from.
    pub tcb: f64,
    /// Vertical centre of buoyancy: depth of `B` below the still-water surface
    /// (`KB`, positive down) [m].
    pub kb: f64,
    /// Band samples at the deck (band top, `z_b = 0`) found immersed — the
    /// modelled band was exceeded, i.e. the deck edge is under water
    /// (down-flooding / geometry beyond the file). A non-zero count means the
    /// volume is missing whatever sits above the modelled sheer.
    pub band_exceeded: usize,
}

/// Affine body→earth map `(x_b, z_b, y_b) → (X, Y, Z)` with `Z` the depth below
/// the still-water surface. Recovered from four evaluations of the rigid
/// transform (constant Jacobian).
struct Affine {
    base: [f64; 3],
    /// Partial derivatives w.r.t. `x_b`, `z_b`, `y_b`.
    dx: [f64; 3],
    dz: [f64; 3],
    dy: [f64; 3],
}

impl Affine {
    /// Earth coordinates of a body point.
    #[inline]
    fn at(&self, xb: f64, zb: f64, yb: f64) -> [f64; 3] {
        [
            self.base[0] + xb * self.dx[0] + zb * self.dz[0] + yb * self.dy[0],
            self.base[1] + xb * self.dx[1] + zb * self.dz[1] + yb * self.dy[1],
            self.base[2] + xb * self.dx[2] + zb * self.dz[2] + yb * self.dy[2],
        ]
    }
}

/// Build the body→earth affine map for the given attitude. Mirrors the frame
/// sequence of [`Body::situate`], with heel inserted as a rotation about the
/// platform centerline at the design floatplane (`Y = 0`, `z_a = 0`, positive
/// heel puts the `+y` side down) — matching [`crate::float::heel_poses`].
fn attitude_affine(
    body: &Body,
    water_offset: f64,
    pose: &HullPose,
    platform: &Platform,
    heel: f64,
) -> Affine {
    let (x0, x1) = body.surface().x_domain();
    let waterline = body.waterline();
    let cp = body.centerplane();
    let zw = -(water_offset + platform.sinkage);
    let pose_pivot = pose.pivot_x.unwrap_or(0.5 * (x0 + x1));
    let (pts, ptc) = pose.trim.sin_cos();
    let (hs, hc) = heel.sin_cos();
    let (fts, ftc) = platform.trim.sin_cos();

    let to_earth = |xb: f64, zb: f64, yb: f64| -> [f64; 3] {
        // 1. align design waterlines and lay in the transverse position.
        let mut x = xb;
        let mut z = zb - waterline;
        let mut y = cp + pose.dy + yb;
        // 2. design pitch about (pose_pivot, z_a = 0), then +dx, +dz.
        if pose.trim != 0.0 {
            let (dx, dz) = (x - pose_pivot, z);
            x = pose_pivot + dx * ptc + dz * pts;
            z = dz * ptc - dx * pts;
        }
        x += pose.dx;
        z += pose.dz;
        // 3. platform heel about (Y = 0, z_a = 0): +y side goes down.
        if heel != 0.0 {
            let (y0, z0) = (y, z);
            y = y0 * hc - z0 * hs;
            z = z0 * hc + y0 * hs;
        }
        // 4. platform pitch about (pivot_x, z_a = zw), a point on the surface.
        if platform.trim != 0.0 {
            let (dx, dz) = (x - platform.pivot_x, z - zw);
            x = platform.pivot_x + dx * ftc + dz * fts;
            z = zw + dz * ftc - dx * fts;
        }
        // 5. depth below the water surface.
        [x, y, z - zw]
    };

    let base = to_earth(0.0, 0.0, 0.0);
    let ex = to_earth(1.0, 0.0, 0.0);
    let ez = to_earth(0.0, 1.0, 0.0);
    let ey = to_earth(0.0, 0.0, 1.0);
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    Affine {
        base,
        dx: sub(ex, base),
        dz: sub(ez, base),
        dy: sub(ey, base),
    }
}

/// Signed area and (area) centroid of a polygon in `(y, z)`, by the shoelace
/// formula. Returns `None` for a degenerate (zero-area) polygon.
fn polygon_area_centroid(poly: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
    if poly.len() < 3 {
        return None;
    }
    let mut a2 = 0.0; // 2·area
    let mut cy = 0.0;
    let mut cz = 0.0;
    for i in 0..poly.len() {
        let (y0, z0) = poly[i];
        let (y1, z1) = poly[(i + 1) % poly.len()];
        let cross = y0 * z1 - y1 * z0;
        a2 += cross;
        cy += (y0 + y1) * cross;
        cz += (z0 + z1) * cross;
    }
    if a2.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    // area = a2/2; centroid = (1/(6·area))·Σ… = (1/(3·a2))·Σ.
    Some((0.5 * a2.abs(), cy / (3.0 * a2), cz / (3.0 * a2)))
}

/// Clip a polygon to the half-plane `cy·y + cz·z + c0 ≥ 0` (Sutherland–Hodgman
/// against one edge). Vertices on the boundary are kept.
fn clip_halfplane(poly: &[(f64, f64)], cy: f64, cz: f64, c0: f64) -> Vec<(f64, f64)> {
    let inside = |p: &(f64, f64)| cy * p.0 + cz * p.1 + c0 >= 0.0;
    let intersect = |a: &(f64, f64), b: &(f64, f64)| {
        let da = cy * a.0 + cz * a.1 + c0;
        let db = cy * b.0 + cz * b.1 + c0;
        let t = da / (da - db);
        (a.0 + t * (b.0 - a.0), a.1 + t * (b.1 - a.1))
    };
    let mut out = Vec::with_capacity(poly.len() + 4);
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        let (ain, bin) = (inside(&a), inside(&b));
        if ain {
            out.push(a);
            if !bin {
                out.push(intersect(&a, &b));
            }
        } else if bin {
            out.push(intersect(&a, &b));
        }
    }
    out
}

/// Immersed hydrostatics of a body at the given attitude, or `None` if it is
/// entirely dry. `water_offset + platform.sinkage` is the water surface depth
/// below the design floatplane (positive sinks the hull deeper); `heel` is the
/// platform heel angle [rad] (positive `+y` side down).
pub fn body_inclined_hydro(
    body: &Body,
    water_offset: f64,
    pose: &HullPose,
    platform: &Platform,
    heel: f64,
    grid: InclinedGrid,
) -> Option<InclinedHydro> {
    let ns = grid.stations.max(2);
    let nb = grid.band.max(2);
    let (x0, x1) = body.surface().x_domain();
    let (_, depth) = body.surface().z_domain();
    let aff = attitude_affine(body, water_offset, pose, platform, heel);

    let hx = (x1 - x0) / (ns - 1) as f64;
    let trap = |i: usize, n: usize| if i == 0 || i == n - 1 { 0.5 } else { 1.0 };
    let ztol = 1e-9 * depth.max(1.0);

    let mut volume = 0.0;
    let mut mx = 0.0;
    let mut my = 0.0;
    let mut mz = 0.0;
    let mut band_exceeded = 0usize;

    // Reusable section-outline buffer (starboard deck→keel, then port keel→deck).
    let mut section = Vec::with_capacity(2 * nb);
    for i in 0..ns {
        let xb = if i == ns - 1 { x1 } else { x0 + hx * i as f64 };
        let wx = trap(i, ns) * hx;

        // Section polygon at this station, in (y_b, z_b).
        section.clear();
        for j in 0..nb {
            let zb = depth * j as f64 / (nb - 1) as f64;
            section.push((body.surface().eval(xb, zb).max(0.0), zb));
        }
        for j in (0..nb).rev() {
            let zb = depth * j as f64 / (nb - 1) as f64;
            section.push((-body.surface().eval(xb, zb).max(0.0), zb));
        }

        // Immersed part: clip by D(x_b, z_b, y_b) ≥ 0, i.e.
        // dy·y_b + dz·z_b + (dx·x_b + base) ≥ 0.
        let clipped = clip_halfplane(&section, aff.dy[2], aff.dz[2], aff.dx[2] * xb + aff.base[2]);
        let Some((area, ybar, zbar)) = polygon_area_centroid(&clipped) else {
            continue;
        };
        let p = aff.at(xb, zbar, ybar);
        volume += wx * area;
        mx += wx * area * p[0];
        my += wx * area * p[1];
        mz += wx * area * p[2];
        // The band top (deck) is under water here: geometry beyond the file.
        if clipped.iter().any(|&(_, zb)| zb <= ztol) {
            band_exceeded += 1;
        }
    }

    if volume <= 1e-12 {
        return None;
    }
    Some(InclinedHydro {
        volume,
        lcb: mx / volume,
        tcb: my / volume,
        kb: mz / volume,
        band_exceeded,
    })
}

/// Aggregate inclined hydrostatics of a whole fleet at one attitude — the
/// quantities an equilibrium solve and a righting-arm read off.
#[derive(Debug, Clone, Copy, Default)]
pub struct FleetInclined {
    /// Total displaced volume [m³].
    pub volume: f64,
    /// Longitudinal buoyancy moment `Σ v·LCB` (earth `x`) [m⁴].
    pub moment_x: f64,
    /// Transverse buoyancy moment `Σ v·TCB` (earth `y`, from the heel axis)
    /// [m⁴].
    pub moment_y: f64,
    /// Source bodies found entirely dry at this attitude.
    pub dry: usize,
    /// Total band-top-immersed samples across the fleet (deck under water).
    pub band_exceeded: usize,
}

/// Inclined hydrostatics summed over a fleet of bodies at one attitude.
pub fn fleet_inclined(
    bodies: &[&Body],
    water_offset: f64,
    poses: &[HullPose],
    platform: &Platform,
    heel: f64,
    grid: InclinedGrid,
) -> FleetInclined {
    let mut f = FleetInclined::default();
    for (body, pose) in bodies.iter().zip(poses) {
        match body_inclined_hydro(body, water_offset, pose, platform, heel, grid) {
            Some(h) => {
                f.volume += h.volume;
                f.moment_x += h.volume * h.lcb;
                f.moment_y += h.volume * h.tcb;
                f.band_exceeded += h.band_exceeded;
            }
            None => f.dry += 1,
        }
    }
    f
}

/// Righting arm `GZ` [m] of a heeled fleet by inclined-waterplane
/// hydrostatics: the earth-frame transverse separation of the buoyancy and
/// gravity lines of action. `vcg` is the centre of gravity metres **above** the
/// design floatplane; `tcg` is its transverse offset in the fleet frame (0 for
/// a laterally symmetric load). The arm is `GZ = TCB − vcg·sin φ − tcg·cos φ`;
/// positive `GZ` rights the platform.
///
/// Unlike [`crate::float::righting_arm`] the per-hull form stability is
/// integrated from the true tilted cut (nonlinear in `φ`), not added
/// metacentrically. Only meaningful when the fleet is in vertical balance
/// (buoyancy = weight) at this `water_offset`.
#[allow(clippy::too_many_arguments)]
pub fn fleet_righting_arm(
    bodies: &[&Body],
    water_offset: f64,
    poses: &[HullPose],
    platform: &Platform,
    heel: f64,
    vcg: f64,
    tcg: f64,
    grid: InclinedGrid,
) -> f64 {
    let f = fleet_inclined(bodies, water_offset, poses, platform, heel, grid);
    if f.volume <= 0.0 {
        return 0.0;
    }
    f.moment_y / f.volume - vcg * heel.sin() - tcg * heel.cos()
}

/// Total displaced volume of a fleet at the given attitude — the vertical-force
/// residual an equilibrium solver drives to the target displacement.
pub fn fleet_volume(
    bodies: &[&Body],
    water_offset: f64,
    poses: &[HullPose],
    platform: &Platform,
    heel: f64,
    grid: InclinedGrid,
) -> f64 {
    fleet_inclined(bodies, water_offset, poses, platform, heel, grid).volume
}

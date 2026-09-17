//! Full-band hull bodies: a half-breadth spline over the hull's entire
//! modelled height, re-situatable (waterline, immersion, trim, position)
//! without the source CAD file.
//!
//! ## Frames
//!
//! A body stores `Y(x, z_b)` with `z_b` measured **downward from the top of
//! the modelled band** (`z_b = 0` at the sheer/deck, `z_b = depth` at the
//! keel), plus the depth of the **design waterline** below the band top.
//! Fleets assemble on the *design floatplane*: every body is aligned so its
//! design waterline sits at assembly height 0, and poses/sinkage move it from
//! there (positive z is down / deeper throughout).
//!
//! ## Why bodies are fast
//!
//! Situating never needs Newton inversion: pitch rotation of a height field
//! `y = f(x, z)` leaves y untouched, so the rotated surface is exactly
//! `f(R⁻¹(x, z))` — a domain lookup. Sampling the wetted band and re-lofting
//! is the whole job.

use crate::bspline::BSplineSurface;
use crate::error::{Error, Result};
use crate::fit::{fit_grid, FitOptions, FitReport};
use crate::grid::SampleGrid;
use crate::hull::Hull;
use crate::iges::{HullPose, Platform};
use crate::michell::Placement;

/// Sampling / loft options for situating a body.
#[derive(Debug, Clone, Copy)]
pub struct BodyOptions {
    pub stations: usize,
    pub waterlines: usize,
    pub fit: FitOptions,
}

impl Default for BodyOptions {
    fn default() -> Self {
        BodyOptions {
            stations: 121,
            waterlines: 33,
            fit: FitOptions {
                degree_x: 3,
                degree_z: 3,
                n_ctrl_x: 20,
                n_ctrl_z: 12,
                ..FitOptions::default()
            },
        }
    }
}

/// A body situated into the water.
#[derive(Debug)]
pub struct SituatedBody {
    pub hull: Hull,
    pub placement: Placement,
    pub fit: FitReport,
    /// The derivative-augmented sample grid the wetted hull was lofted from.
    pub grid: SampleGrid,
    /// Wetted samples that mapped **above** the stored band (the water rose
    /// past the modelled sheer): their geometry is unknown and was taken as
    /// zero. A non-zero count means the pose exceeds what the file covers.
    pub band_exceeded: usize,
}

/// A full-band hull body.
#[derive(Debug, Clone)]
pub struct Body {
    surface: BSplineSurface,
    waterline: f64,
    centerplane: f64,
}

impl Body {
    /// `surface` is the half-breadth over the full band with the z domain
    /// starting at 0 (band top); `waterline` is the design waterline's depth
    /// below the band top; `centerplane` the hull's transverse position.
    pub fn new(surface: BSplineSurface, waterline: f64, centerplane: f64) -> Result<Body> {
        let (z0, z1) = surface.z_domain();
        if z0 != 0.0 {
            return Err(Error::InvalidGeometry(format!(
                "a body's z domain must start at 0 (the band top); got {z0}"
            )));
        }
        if !(waterline.is_finite() && (0.0..=z1).contains(&waterline)) {
            return Err(Error::InvalidGeometry(format!(
                "design waterline {waterline} outside the body's band [0, {z1}]"
            )));
        }
        if !centerplane.is_finite() {
            return Err(Error::InvalidGeometry("centerplane must be finite".into()));
        }
        let scale = surface
            .control()
            .iter()
            .fold(0.0f64, |m, &v| m.max(v.abs()));
        if surface
            .control()
            .iter()
            .any(|&v| v < -1e-12 * scale.max(1.0))
        {
            return Err(Error::InvalidGeometry(
                "body control net contains negative half-beam values".into(),
            ));
        }
        Ok(Body {
            surface,
            waterline,
            centerplane,
        })
    }

    pub fn surface(&self) -> &BSplineSurface {
        &self.surface
    }

    /// Depth of the design waterline below the band top [m].
    pub fn waterline(&self) -> f64 {
        self.waterline
    }

    /// Transverse position of the hull's centerplane [m].
    pub fn centerplane(&self) -> f64 {
        self.centerplane
    }

    /// Situate the body: apply the design pose and platform state, clip at
    /// the water surface `water_offset + platform.sinkage` below the design
    /// floatplane, and loft the wetted part. `Ok(None)` means the body is
    /// entirely dry at this pose.
    pub fn situate(
        &self,
        water_offset: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &BodyOptions,
    ) -> Result<Option<SituatedBody>> {
        if opts.stations < 8 || opts.waterlines < 6 {
            return Err(Error::InvalidInput(
                "need at least 8 stations and 6 waterlines to sample".into(),
            ));
        }
        if !(pose.scale > 0.0 && pose.scale.is_finite()) {
            return Err(Error::InvalidInput(format!(
                "hull scale must be a positive, finite factor; got {}",
                pose.scale
            )));
        }
        // Water surface position in assembly coordinates (down-positive):
        // sinking the platform (positive) puts the water ABOVE the design
        // floatplane, i.e. at negative z_a.
        let zw = -(water_offset + platform.sinkage);
        let (x0, x1) = self.surface.x_domain();
        let (_, depth) = self.surface.z_domain();
        let px = pose.pivot_x.unwrap_or(0.5 * (x0 + x1));
        let map = FrameMap {
            waterline: self.waterline,
            pose: *pose,
            pose_pivot_x: px,
            platform: *platform,
            zw,
        };

        // Forward-scan the body per body-fixed x column: the column's deepest
        // immersion m(x) is a max of finitely many functions continuous in
        // the pose/state, so extents derived from its zero crossings vary
        // CONTINUOUSLY with sinkage and trim. This matters: snapping extents
        // to scan cells makes V(sinkage) a staircase whose steps can exceed
        // the equilibrium tolerance, trapping the solver between treads.
        const SCAN: usize = 97;
        let mut col_x = [0.0f64; SCAN]; // water-frame x at the column's deepest point
        let mut col_m = [f64::NEG_INFINITY; SCAN]; // deepest immersion of the column
        for i in 0..SCAN {
            let xb = x0 + (x1 - x0) * i as f64 / (SCAN - 1) as f64;
            for j in 0..SCAN {
                let zb = depth * j as f64 / (SCAN - 1) as f64;
                let (xw, zdepth) = map.body_to_water(xb, zb);
                if zdepth > col_m[i] {
                    col_m[i] = zdepth;
                    col_x[i] = xw;
                }
            }
        }
        let draft = col_m.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        if !(draft > 1e-12) {
            return Ok(None);
        }
        let first = col_m.iter().position(|&m| m >= 0.0).expect("draft > 0");
        let last = col_m.iter().rposition(|&m| m >= 0.0).expect("draft > 0");
        // Interpolate the wet/dry crossing into the neighbouring dry column.
        let mut wx_lo = col_x[first];
        let mut wx_hi = col_x[last];
        if first > 0 {
            let (ma, mb) = (col_m[first - 1], col_m[first]);
            let t = ma / (ma - mb);
            wx_lo = col_x[first - 1] + (col_x[first] - col_x[first - 1]) * t;
        }
        if last + 1 < SCAN {
            let (ma, mb) = (col_m[last + 1], col_m[last]);
            let t = ma / (ma - mb);
            wx_hi = col_x[last + 1] + (col_x[last] - col_x[last + 1]) * t;
        }
        if wx_lo > wx_hi {
            std::mem::swap(&mut wx_lo, &mut wx_hi);
        }
        if !(wx_hi > wx_lo) {
            return Ok(None);
        }

        // Water-frame sample grid: cosine stations, uniform waterlines.
        let ns = opts.stations;
        let nw = opts.waterlines;
        let stations: Vec<f64> = (0..ns)
            .map(|i| {
                let c = (std::f64::consts::PI * i as f64 / (ns - 1) as f64).cos();
                wx_lo + (wx_hi - wx_lo) * (1.0 - c) / 2.0
            })
            .collect();
        let waterlines: Vec<f64> = (0..nw)
            .map(|j| draft * j as f64 / (nw - 1) as f64)
            .collect();

        // The water → body map is affine, so its Jacobian is constant and
        // exactly recovered from three evaluations; the sampled slopes are
        // the body slopes pushed through it by the chain rule.
        let o = map.water_to_body(0.0, 0.0);
        let jx = map.water_to_body(1.0, 0.0); // column d(xb, zb)/dxw + o
        let jz = map.water_to_body(0.0, 1.0); // column d(xb, zb)/dzw + o
        let (jxx, jzx) = (jx.0 - o.0, jx.1 - o.1);
        let (jxz, jzz) = (jz.0 - o.0, jz.1 - o.1);

        let mut grid = vec![0.0f64; ns * nw];
        let mut fx = vec![f64::NAN; ns * nw];
        let mut fz = vec![f64::NAN; ns * nw];
        let mut band_exceeded = 0usize;
        // The forward scan and this inverse map round-trip through the same
        // rotations, so grid rows meant to land exactly on a domain edge (the
        // keel row, the end stations) can overshoot by roundoff; clamp a
        // whisker rather than silently zeroing them.
        let xtol = 1e-9 * (x1 - x0);
        let ztol = 1e-9 * depth;
        for (i, &xw) in stations.iter().enumerate() {
            for (j, &zw) in waterlines.iter().enumerate() {
                let (xb, zb) = map.water_to_body(xw, zw);
                let s = i * nw + j;
                if xb < x0 - xtol || xb > x1 + xtol {
                    continue; // beyond the ends: no hull, half-beam 0 is data
                }
                if zb > depth + ztol {
                    continue; // below the keel
                }
                if zb < -ztol {
                    // Above the modelled band while under water: unknown
                    // geometry, taken as zero so the loft still covers the
                    // wetted rectangle (excluding the whole wedge would
                    // leave control points unconstrained); the slope stays
                    // unconstrained and the count reports that the pose
                    // exceeds what the file models.
                    band_exceeded += 1;
                    continue;
                }
                let (xb, zb) = (xb.clamp(x0, x1), zb.clamp(0.0, depth));
                // The half-beam scales with the hull; the water→body map (hence
                // the Jacobian below) already carries the reciprocal, so the
                // chain rule leaves one factor of `scale` on both slopes.
                let sc = pose.scale;
                grid[s] = self.surface.eval(xb, zb).max(0.0) * sc;
                let fxb = self.surface.eval_deriv(xb, zb, 1, 0);
                let fzb = self.surface.eval_deriv(xb, zb, 0, 1);
                fx[s] = (fxb * jxx + fzb * jzx) * sc;
                fz[s] = (fxb * jxz + fzb * jzz) * sc;
            }
        }

        let grid = SampleGrid::new(stations, waterlines, grid)?
            .with_fx(fx)?
            .with_fz(fz)?;
        let (hull, fit) = fit_grid(&grid, &opts.fit)?;
        Ok(Some(SituatedBody {
            hull,
            placement: Placement {
                x: 0.0,
                y: self.centerplane + pose.dy,
            },
            fit,
            grid,
            band_exceeded,
        }))
    }
}

/// The body → water frame map (z positive down everywhere):
/// 1. align design waterlines: `z_a = z_b - waterline`;
/// 2. design trim about `(pose_pivot_x, z_a = 0)`, then `+dx`, `+dz`;
/// 3. platform trim about `(platform.pivot_x, z_a = zw)` — a point on the
///    water surface;
/// 4. depth below water: `z' = z_a - zw`.
struct FrameMap {
    waterline: f64,
    pose: HullPose,
    pose_pivot_x: f64,
    platform: Platform,
    /// Water surface in assembly coordinates: `-(water_offset + sinkage)`.
    zw: f64,
}

/// Rotation in down-positive coordinates; positive angle raises the +x side.
#[inline]
fn rot_down(x: f64, zd: f64, px: f64, pzd: f64, sin: f64, cos: f64) -> (f64, f64) {
    let (dx, dz) = (x - px, zd - pzd);
    (px + dx * cos + dz * sin, pzd + dz * cos - dx * sin)
}

impl FrameMap {
    fn body_to_water(&self, xb: f64, zb: f64) -> (f64, f64) {
        // Uniform scale first, about the pivot station and the design
        // waterline (z = 0 in this design frame): the hull grows or shrinks
        // in place, then the pose/state transforms position it.
        let sc = self.pose.scale;
        let mut x = self.pose_pivot_x + sc * (xb - self.pose_pivot_x);
        let mut z = sc * (zb - self.waterline);
        if self.pose.trim != 0.0 {
            let (s, c) = self.pose.trim.sin_cos();
            (x, z) = rot_down(x, z, self.pose_pivot_x, 0.0, s, c);
        }
        x += self.pose.dx;
        z += self.pose.dz;
        if self.platform.trim != 0.0 {
            let (s, c) = self.platform.trim.sin_cos();
            (x, z) = rot_down(x, z, self.platform.pivot_x, self.zw, s, c);
        }
        (x, z - self.zw)
    }

    fn water_to_body(&self, xw: f64, zdepth: f64) -> (f64, f64) {
        let mut x = xw;
        let mut z = zdepth + self.zw;
        if self.platform.trim != 0.0 {
            let (s, c) = (-self.platform.trim).sin_cos();
            (x, z) = rot_down(x, z, self.platform.pivot_x, self.zw, s, c);
        }
        x -= self.pose.dx;
        z -= self.pose.dz;
        if self.pose.trim != 0.0 {
            let (s, c) = (-self.pose.trim).sin_cos();
            (x, z) = rot_down(x, z, self.pose_pivot_x, 0.0, s, c);
        }
        // Undo the uniform scale (inverse of `body_to_water`'s first step).
        let sc = self.pose.scale;
        let x = self.pose_pivot_x + (x - self.pose_pivot_x) / sc;
        let z = z / sc;
        (x, z + self.waterline)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frame_map_roundtrips() {
        let map = FrameMap {
            waterline: 0.4,
            pose: HullPose {
                dx: 1.2,
                dy: 0.0,
                dz: -0.07,
                trim: 0.03,
                pivot_x: Some(3.0),
                ..Default::default()
            },
            pose_pivot_x: 3.0,
            platform: Platform {
                sinkage: 0.05,
                trim: -0.02,
                pivot_x: 5.5,
            },
            zw: -0.11,
        };
        for &(xb, zb) in &[(0.0, 0.0), (2.7, 0.31), (10.0, 0.65), (-4.0, 1.0)] {
            let (xw, zw) = map.body_to_water(xb, zb);
            let (xb2, zb2) = map.water_to_body(xw, zw);
            assert!((xb - xb2).abs() < 1e-12 && (zb - zb2).abs() < 1e-12);
        }
    }

    #[test]
    fn positive_trim_raises_the_positive_x_end() {
        let map = FrameMap {
            waterline: 0.0,
            pose: HullPose {
                trim: 0.1,
                ..Default::default()
            },
            pose_pivot_x: 0.0,
            platform: Platform::default(),
            zw: 0.0,
        };
        // A point forward of the pivot at the waterline must move up
        // (smaller depth) under positive trim.
        let (_, z) = map.body_to_water(5.0, 0.0);
        assert!(z < 0.0, "z {z}");
    }

    #[test]
    fn situate_grid_slopes_match_value_samples() {
        // A posed Wigley body: the sampled slope channels must agree with
        // finite differences of the sampled values. The surface is
        // biquadratic and the pose map affine, so a 3-point non-uniform
        // central difference is exact up to rounding wherever the stencil
        // stays inside the hull.
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let body = Body::new(hull.surface().clone(), 0.2, 0.0).unwrap();
        let pose = HullPose {
            dz: 0.05,
            trim: 0.04,
            ..Default::default()
        };
        let sb = body
            .situate(0.0, &pose, &Platform::default(), &BodyOptions::default())
            .unwrap()
            .expect("wet");
        assert!(sb.fit.fx_residual.is_some() && sb.fit.fz_residual.is_some());
        let g = &sb.grid;
        let (st, wl) = (g.stations(), g.waterlines());
        let (f, fx, fz) = (g.half_beams(), g.fx().unwrap(), g.fz().unwrap());
        let (mx, mz) = (st.len(), wl.len());
        let mut checked = 0usize;
        for i in 1..mx - 1 {
            for j in 1..mz - 1 {
                let s = i * mz + j;
                let stencil = [s, s - mz, s + mz, s - 1, s + 1];
                if stencil.iter().any(|&k| f[k] <= 1e-3)
                    || !(fx[s].is_finite() && fz[s].is_finite())
                {
                    continue;
                }
                let (h0, h1) = (st[i] - st[i - 1], st[i + 1] - st[i]);
                let dfdx = (h0 * h0 * f[s + mz] + (h1 * h1 - h0 * h0) * f[s] - h1 * h1 * f[s - mz])
                    / (h0 * h1 * (h0 + h1));
                let (k0, k1) = (wl[j] - wl[j - 1], wl[j + 1] - wl[j]);
                let dfdz = (k0 * k0 * f[s + 1] + (k1 * k1 - k0 * k0) * f[s] - k1 * k1 * f[s - 1])
                    / (k0 * k1 * (k0 + k1));
                assert!(
                    (fx[s] - dfdx).abs() < 1e-4,
                    "i={i} j={j}: fx {} vs FD {dfdx}",
                    fx[s]
                );
                assert!(
                    (fz[s] - dfdz).abs() < 1e-4,
                    "i={i} j={j}: fz {} vs FD {dfdz}",
                    fz[s]
                );
                checked += 1;
            }
        }
        assert!(checked > 100, "only {checked} interior samples checked");
    }

    #[test]
    fn platform_sinkage_equals_pose_dz() {
        // Sinking the platform by s must immerse a hull exactly like
        // lowering it by dz = s.
        let sunk = FrameMap {
            waterline: 0.3,
            pose: HullPose::default(),
            pose_pivot_x: 0.0,
            platform: Platform {
                sinkage: 0.17,
                trim: 0.0,
                pivot_x: 0.0,
            },
            zw: -0.17,
        };
        let lowered = FrameMap {
            waterline: 0.3,
            pose: HullPose {
                dz: 0.17,
                ..Default::default()
            },
            pose_pivot_x: 0.0,
            platform: Platform::default(),
            zw: 0.0,
        };
        for &(xb, zb) in &[(0.0, 0.0), (4.0, 0.5)] {
            let a = sunk.body_to_water(xb, zb);
            let b = lowered.body_to_water(xb, zb);
            assert!((a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12);
        }
    }

    #[test]
    fn scale_maps_roundtrip() {
        // A scaled frame map must still invert exactly.
        let map = FrameMap {
            waterline: 0.4,
            pose: HullPose {
                dx: 1.2,
                dz: -0.07,
                trim: 0.03,
                scale: 1.7,
                pivot_x: Some(3.0),
                ..Default::default()
            },
            pose_pivot_x: 3.0,
            platform: Platform {
                sinkage: 0.05,
                trim: -0.02,
                pivot_x: 5.5,
            },
            zw: -0.11,
        };
        for &(xb, zb) in &[(0.0, 0.0), (2.7, 0.31), (10.0, 0.65), (-4.0, 1.0)] {
            let (xw, zw) = map.body_to_water(xb, zb);
            let (xb2, zb2) = map.water_to_body(xw, zw);
            assert!((xb - xb2).abs() < 1e-12 && (zb - zb2).abs() < 1e-12);
        }
    }

    #[test]
    fn scale_is_geometrically_similar() {
        // Uniform scale must reproduce a geometrically similar hull: length
        // and draft go as s, wetted area as s^2, displaced volume as s^3.
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let body = Body::new(hull.surface().clone(), 0.4, 0.0).unwrap();
        let opts = BodyOptions::default();
        let base = body
            .situate(0.0, &HullPose::default(), &Platform::default(), &opts)
            .unwrap()
            .expect("wet");
        for &s in &[0.5, 1.5, 2.0] {
            let scaled = body
                .situate(
                    0.0,
                    &HullPose {
                        scale: s,
                        ..Default::default()
                    },
                    &Platform::default(),
                    &opts,
                )
                .unwrap()
                .expect("wet");
            assert_eq!(scaled.band_exceeded, 0);
            let (b, c) = (&base.hull, &scaled.hull);
            let rel = |got: f64, want: f64| (got - want).abs() <= 1e-4 * want.abs().max(1e-9);
            assert!(
                rel(c.length(), s * b.length()),
                "s={s}: length {} vs {}",
                c.length(),
                s * b.length()
            );
            assert!(
                rel(c.draft(), s * b.draft()),
                "s={s}: draft {} vs {}",
                c.draft(),
                s * b.draft()
            );
            assert!(
                rel(c.wetted_surface(), s * s * b.wetted_surface()),
                "s={s}: wetted {} vs {}",
                c.wetted_surface(),
                s * s * b.wetted_surface()
            );
            assert!(
                rel(c.displaced_volume(), s * s * s * b.displaced_volume()),
                "s={s}: volume {} vs {}",
                c.displaced_volume(),
                s * s * s * b.displaced_volume()
            );
        }
    }

    #[test]
    fn nonpositive_scale_is_rejected() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let body = Body::new(hull.surface().clone(), 0.4, 0.0).unwrap();
        for bad in [0.0, -1.0, f64::NAN] {
            let r = body.situate(
                0.0,
                &HullPose {
                    scale: bad,
                    ..Default::default()
                },
                &Platform::default(),
                &BodyOptions::default(),
            );
            assert!(r.is_err(), "scale {bad} should be rejected");
        }
    }
}

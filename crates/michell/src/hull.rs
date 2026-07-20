//! Validated hull wrapper: a B-spline half-breadth surface plus everything
//! the resistance computations need, precomputed once.

use crate::bspline::BSplineSurface;
use crate::error::{Error, Result};
use crate::quadrature::gauss_legendre;

/// One non-empty knot span in one direction.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Span {
    /// Physical coordinate of the span start.
    pub start: f64,
    /// Span length (> 0).
    pub len: f64,
}

/// A validated hull.
///
/// Geometry contract (see crate docs): `y = f(x, z) >= 0` is the local
/// half-beam of a port/starboard-symmetric hull; `x` runs along the hull
/// (arbitrary origin), `z` runs vertically **downward** from the undisturbed
/// waterline at `z = 0`; the surface domain is `x ∈ [x0, x1]`, `z ∈ [0, T]`.
#[derive(Debug, Clone)]
pub struct Hull {
    surface: BSplineSurface,
    spans_x: Vec<Span>,
    spans_z: Vec<Span>,
    /// Local polynomial coefficients of ∂f/∂x per span pair:
    /// `fx(x, z) = Σ_{a,b} c[a][b] (x - x0)^a (z - z0)^b`,
    /// flattened as `[(sx * spans_z.len() + sz) * p + a][(q+1) elements over b]`,
    /// i.e. index `((sx * nsz + sz) * p + a) * (q + 1) + b`,
    /// with `a = 0..p` (p = degree_x) and `b = 0..=q` (q = degree_z).
    fx_coeff: Vec<f64>,
    length: f64,
    draft: f64,
    x_center: f64,
    wetted_surface: f64,
    displaced_volume: f64,
    lcb_x: f64,
    vcb_z: f64,
    waterplane_area: f64,
    waterplane_moment: f64,
    waterplane_second_moment: f64,
    waterplane_transverse_moment: f64,
}

impl Hull {
    /// Validate the surface as a hull and precompute span polynomials and
    /// geometric integrals.
    ///
    /// Beyond spline validity this requires:
    /// - the z-domain starts at the waterline: `z0 = 0` (small negative
    ///   tolerance rejected) and has positive draft;
    /// - a non-negative control net (a conservative sufficient condition for
    ///   `f >= 0`, by the B-spline convex-hull property).
    pub fn new(surface: BSplineSurface) -> Result<Hull> {
        let (z0, z1) = surface.z_domain();
        let (x0, x1) = surface.x_domain();
        let draft = z1;
        if z0 < 0.0 || z0 > 1e-9 * (z1 - z0) {
            return Err(Error::InvalidGeometry(format!(
                "z-domain must start at the waterline z = 0 (z measured downward); got z0 = {z0}"
            )));
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
                "control net contains negative half-beam values; the half-breadth \
                 surface must satisfy f(x, z) >= 0"
                    .into(),
            ));
        }

        let p = surface.degree_x();
        let q = surface.degree_z();
        let xs: Vec<usize> = surface.x_span_indices();
        let zs: Vec<usize> = surface.z_span_indices();
        let spans_x: Vec<Span> = xs
            .iter()
            .map(|&s| Span {
                start: surface.knots_x()[s],
                len: surface.knots_x()[s + 1] - surface.knots_x()[s],
            })
            .collect();
        let spans_z: Vec<Span> = zs
            .iter()
            .map(|&s| Span {
                start: surface.knots_z()[s],
                len: surface.knots_z()[s + 1] - surface.knots_z()[s],
            })
            .collect();

        // Factorials up to max degree (degrees are small).
        let mut fact = vec![1.0f64; p.max(q) + 2];
        for i in 1..fact.len() {
            fact[i] = fact[i - 1] * i as f64;
        }

        // Local polynomial coefficients of fx = ∂f/∂x on every span pair.
        let nsz = zs.len();
        let mut fx_coeff = vec![0.0f64; xs.len() * nsz * p * (q + 1)];
        for (isx, &sx) in xs.iter().enumerate() {
            for (isz, &sz) in zs.iter().enumerate() {
                let d = surface.corner_partials(sx, sz);
                for a in 0..p {
                    for b in 0..=q {
                        // f = Σ D[a][b]/(a! b!) X^a Z^b  =>
                        // fx coefficient of X^a Z^b is D[a+1][b]/(a! b!).
                        fx_coeff[((isx * nsz + isz) * p + a) * (q + 1) + b] =
                            d[a + 1][b] / (fact[a] * fact[b]);
                    }
                }
            }
        }

        // Geometric integrals by per-span Gauss-Legendre.
        // Volume: integrand is polynomial of degree (p, q) => exact.
        // Wetted surface: smooth integrand, use a generous rule.
        // x·f and z·f raise one degree by one; +2 keeps the rule exact.
        let n_vol = (p.max(q) + 2) / 2 + 2;
        // Wetted-surface integrand is smooth but non-polynomial; a hull may be
        // a single span, so use a high-order rule per span.
        let n_wet = 24;
        let (xv, wv) = gauss_legendre(n_vol);
        let (xw, ww) = gauss_legendre(n_wet);
        let mut volume = 0.0;
        let mut volume_mx = 0.0;
        let mut volume_mz = 0.0;
        let mut wetted = 0.0;
        for sx in &spans_x {
            for sz in &spans_z {
                let jac = sx.len / 2.0 * (sz.len / 2.0);
                for (i, &xi) in xv.iter().enumerate() {
                    let x = sx.start + sx.len * (xi + 1.0) / 2.0;
                    for (j, &zj) in xv.iter().enumerate() {
                        let z = zj_map(sz, zj);
                        let f = surface.eval(x, z);
                        volume += wv[i] * wv[j] * jac * f;
                        volume_mx += wv[i] * wv[j] * jac * x * f;
                        volume_mz += wv[i] * wv[j] * jac * z * f;
                    }
                }
                for (i, &xi) in xw.iter().enumerate() {
                    let x = sx.start + sx.len * (xi + 1.0) / 2.0;
                    for (j, &zj) in xw.iter().enumerate() {
                        let z = zj_map(sz, zj);
                        let fx = surface.eval_deriv(x, z, 1, 0);
                        let fz = surface.eval_deriv(x, z, 0, 1);
                        wetted += ww[i] * ww[j] * jac * (1.0 + fx * fx + fz * fz).sqrt();
                    }
                }
            }
        }
        // Both sides of the hull.
        volume *= 2.0;
        volume_mx *= 2.0;
        volume_mz *= 2.0;
        wetted *= 2.0;

        // Waterplane properties: 1-D integrals of the beam b(x) = 2 f(x, 0)
        // (the transverse inertia integrand f³ has degree 3p; the rule is
        // sized for it, which also covers x²·f at degree p + 2).
        let (xq, wq) = gauss_legendre(3 * p / 2 + 2);
        let mut wp_area = 0.0;
        let mut wp_mx = 0.0;
        let mut wp_ixx = 0.0;
        let mut wp_iyy = 0.0;
        for sx in &spans_x {
            let jac = sx.len / 2.0;
            for (i, &xi) in xq.iter().enumerate() {
                let x = sx.start + sx.len * (xi + 1.0) / 2.0;
                let f = surface.eval(x, z0.max(0.0));
                let b = 2.0 * f;
                wp_area += wq[i] * jac * b;
                wp_mx += wq[i] * jac * x * b;
                wp_ixx += wq[i] * jac * x * x * b;
                wp_iyy += wq[i] * jac * (2.0 / 3.0) * f * f * f;
            }
        }

        Ok(Hull {
            surface,
            spans_x,
            spans_z,
            fx_coeff,
            length: x1 - x0,
            draft,
            x_center: 0.5 * (x0 + x1),
            wetted_surface: wetted,
            displaced_volume: volume,
            lcb_x: if volume > 0.0 { volume_mx / volume } else { 0.0 },
            vcb_z: if volume > 0.0 { volume_mz / volume } else { 0.0 },
            waterplane_area: wp_area,
            waterplane_moment: wp_mx,
            waterplane_second_moment: wp_ixx,
            waterplane_transverse_moment: wp_iyy,
        })
    }

    pub fn surface(&self) -> &BSplineSurface {
        &self.surface
    }

    /// Length of the x-domain [m]. Used as the reference length.
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Depth of the z-domain below the waterline [m].
    pub fn draft(&self) -> f64 {
        self.draft
    }

    /// Wetted surface area of both sides [m²], computed on the centerplane
    /// projection: `S = 2 ∬ √(1 + fx² + fz²) dx dz` (consistent with
    /// thin-ship theory; no girth correction).
    pub fn wetted_surface(&self) -> f64 {
        self.wetted_surface
    }

    /// Displaced volume `∇ = 2 ∬ f dx dz` [m³] (thin-ship approximation).
    pub fn displaced_volume(&self) -> f64 {
        self.displaced_volume
    }

    /// Longitudinal centre of buoyancy `x_B = ∬ x f / ∬ f` [m], in the hull's
    /// x coordinates.
    pub fn lcb_x(&self) -> f64 {
        self.lcb_x
    }

    /// Vertical centre of buoyancy `z_B = ∬ z f / ∬ f` [m], measured
    /// **downward** from the waterline (KB below the water surface).
    pub fn vcb_z(&self) -> f64 {
        self.vcb_z
    }

    /// Waterplane area `A_w = ∫ 2 f(x, 0) dx` [m²].
    pub fn waterplane_area(&self) -> f64 {
        self.waterplane_area
    }

    /// First moment of the waterplane about x = 0: `∫ 2 x f(x, 0) dx` [m³].
    pub fn waterplane_moment(&self) -> f64 {
        self.waterplane_moment
    }

    /// Second moment of the waterplane about x = 0: `∫ 2 x² f(x, 0) dx` [m⁴].
    pub fn waterplane_second_moment(&self) -> f64 {
        self.waterplane_second_moment
    }

    /// Second moment of the waterplane about the hull's own centerplane:
    /// `I_T = ∫ (2/3) f(x, 0)³ dx` [m⁴] — the transverse inertia that sets
    /// the hull's individual metacentric radius `BM_T = I_T / ∇`.
    pub fn waterplane_transverse_moment(&self) -> f64 {
        self.waterplane_transverse_moment
    }

    /// Longitudinal centre of flotation (waterplane centroid) [m].
    pub fn lcf_x(&self) -> f64 {
        if self.waterplane_area > 0.0 {
            self.waterplane_moment / self.waterplane_area
        } else {
            0.0
        }
    }

    pub(crate) fn spans_x(&self) -> &[Span] {
        &self.spans_x
    }

    pub(crate) fn spans_z(&self) -> &[Span] {
        &self.spans_z
    }

    pub(crate) fn fx_coeff(&self) -> &[f64] {
        &self.fx_coeff
    }

    /// Half-extent of the hull about its x-midpoint (bandwidth of the
    /// oscillatory inner integral).
    pub(crate) fn x_half_extent(&self) -> f64 {
        self.length / 2.0
    }

    pub(crate) fn x_center(&self) -> f64 {
        self.x_center
    }
}

#[inline]
fn zj_map(sz: &Span, node: f64) -> f64 {
    sz.start + sz.len * (node + 1.0) / 2.0
}

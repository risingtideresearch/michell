//! Validated hull wrapper: a B-spline half-breadth surface plus everything
//! the resistance computations need, precomputed once.

use crate::bspline::{ders_basis, BSplineSurface};
use crate::error::{Error, Result};
use crate::quadrature::gauss_legendre;

/// Largest spline degree accepted by [`Hull`] in either parametric direction.
///
/// The resistance kernels use local power-basis coefficients. Independent
/// degree-elevation and high-precision tests validate that representation
/// through degree 16. Higher degrees are rejected because the former
/// corner-derivative/Taylor conversion becomes ill-conditioned and can return
/// finite but catastrophically wrong resistance values.
pub const MAX_SUPPORTED_SPLINE_DEGREE: usize = 16;

// The independent high-precision degree-(16, 12) stress case differed by
// 1.12e-9 in resistance. This 2e-8 validated floor is intentionally about 18
// times larger; the permanent exact degree-elevation sweep enforces it across
// the supported envelope, and every accepted result is charged for it.
const VALIDATED_COEFFICIENT_RESISTANCE_REL_ERROR: f64 = 2.0e-8;
const MAX_SPAN_RECONSTRUCTION_REL_ERROR: f64 = 1.0e-8;

/// One non-empty knot span in one direction.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Span {
    /// Physical coordinate of the span start.
    pub start: f64,
    /// Span length (> 0).
    pub len: f64,
}

/// Constraint derivatives with respect to one B-spline control net.
#[derive(Debug, Clone)]
pub struct ConstraintGradient {
    /// `∂∇/∂Pᵢ` \[m²] for displaced volume `∇`.
    pub displaced_volume: Vec<f64>,
    /// `∂x_B/∂Pᵢ` for the longitudinal centre of buoyancy.
    pub lcb_x: Vec<f64>,
    /// `∂S/∂Pᵢ` \[m] for wetted surface area `S`.
    pub wetted_surface: Vec<f64>,
}

/// Hull-constraint derivatives for the physical control net or nets.
#[derive(Debug, Clone)]
pub enum HullConstraintGradients {
    /// A symmetric half-breadth control net, representing both physical sides.
    Symmetric(ConstraintGradient),
    /// Independent derivatives for the two physical half-breadth control nets.
    Asymmetric {
        port: ConstraintGradient,
        starboard: ConstraintGradient,
    },
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
    /// Validated resistance-level error floor for the shared local coefficient
    /// representation. Both numerical routes charge this same contribution.
    fx_coeff_rel_error_bound: f64,
    /// `∂f_a/∂x` coefficients for the **antisymmetric** half-beam
    /// `f_a = (f₊ − f₋)/2` of an asymmetric hull, same layout as `fx_coeff`.
    /// `None` for a port/starboard-symmetric hull (the default contract), in
    /// which case every wave computation reduces exactly to classical Michell.
    fx_a_coeff: Option<Vec<f64>>,
    /// The antisymmetric half-beam surface `f_a = (f₊ − f₋)/2` itself, kept for
    /// evaluating `∂f_a/∂x` at arbitrary `(x, z)` (the camber forcing of the
    /// centreplane lifting solve). `Some` exactly when `fx_a_coeff` is.
    a_surface: Option<BSplineSurface>,
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
        validate_resistance_degree(&surface)?;
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

        // Local polynomial coefficients of fx = ∂f/∂x on every span pair.
        let (fx_coeff, fx_coeff_rel_error_bound) = compute_fx_coeff(&surface, &xs, &zs)?;

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
            fx_coeff_rel_error_bound,
            fx_a_coeff: None,
            a_surface: None,
            length: x1 - x0,
            draft,
            x_center: 0.5 * (x0 + x1),
            wetted_surface: wetted,
            displaced_volume: volume,
            lcb_x: if volume > 0.0 {
                volume_mx / volume
            } else {
                0.0
            },
            vcb_z: if volume > 0.0 {
                volume_mz / volume
            } else {
                0.0
            },
            waterplane_area: wp_area,
            waterplane_moment: wp_mx,
            waterplane_second_moment: wp_ixx,
            waterplane_transverse_moment: wp_iyy,
        })
    }

    /// Validate an **asymmetric** hull from its two half-breadth surfaces.
    ///
    /// `starboard` is the half-beam `y = +f₊(x, z) >= 0` and `port` is
    /// `y = −f₋(x, z) <= 0` (both surfaces store the *magnitude* `f ≥ 0`). The
    /// hull is split into a symmetric thickness part `f_sym = (f₊ + f₋)/2` and
    /// an antisymmetric camber part `f_a = (f₊ − f₋)/2`; the wave field then
    /// superposes a source system from `f_sym` (classical Michell) and a
    /// centreplane y-dipole system from `f_a` (see the `michell` module docs).
    ///
    /// Contract: both surfaces must share **identical degrees and knot
    /// vectors** (they may differ only in their control nets), so the two
    /// decomposed surfaces live on the same spans. Each side must be a valid
    /// non-negative half-breadth. When `port == starboard` the result is bit
    /// -for-bit the symmetric [`Hull::new`] with `f_a ≡ 0`.
    pub fn new_asymmetric(port: BSplineSurface, starboard: BSplineSurface) -> Result<Hull> {
        if starboard.degree_x() != port.degree_x() || starboard.degree_z() != port.degree_z() {
            return Err(Error::InvalidGeometry(
                "asymmetric hull: port and starboard surfaces must share degrees".into(),
            ));
        }
        let knots_match = |a: &[f64], b: &[f64]| {
            a.len() == b.len()
                && a.iter()
                    .zip(b)
                    .all(|(x, y)| (x - y).abs() <= 1e-12 * (1.0 + x.abs()))
        };
        if !knots_match(starboard.knots_x(), port.knots_x())
            || !knots_match(starboard.knots_z(), port.knots_z())
        {
            return Err(Error::InvalidGeometry(
                "asymmetric hull: port and starboard surfaces must share knot vectors \
                 (only the control nets may differ)"
                    .into(),
            ));
        }
        if starboard.control().len() != port.control().len() {
            return Err(Error::InvalidGeometry(
                "asymmetric hull: control nets differ in length".into(),
            ));
        }
        for (side, surface) in [("port", &port), ("starboard", &starboard)] {
            let scale = surface
                .control()
                .iter()
                .fold(0.0f64, |maximum, &value| maximum.max(value.abs()));
            if surface
                .control()
                .iter()
                .any(|&value| value < -1e-12 * scale.max(1.0))
            {
                return Err(Error::InvalidGeometry(format!(
                    "asymmetric hull: {side} control net contains negative half-beam values"
                )));
            }
        }

        // Symmetric and antisymmetric control nets on the shared parametrisation.
        let sym_ctrl: Vec<f64> = starboard
            .control()
            .iter()
            .zip(port.control())
            .map(|(s, p)| 0.5 * (s + p))
            .collect();
        let a_ctrl: Vec<f64> = starboard
            .control()
            .iter()
            .zip(port.control())
            .map(|(s, p)| 0.5 * (s - p))
            .collect();
        let sym_surface = BSplineSurface::new(
            starboard.degree_x(),
            starboard.degree_z(),
            starboard.knots_x().to_vec(),
            starboard.knots_z().to_vec(),
            sym_ctrl,
        )?;
        let a_surface = BSplineSurface::new(
            starboard.degree_x(),
            starboard.degree_z(),
            starboard.knots_x().to_vec(),
            starboard.knots_z().to_vec(),
            a_ctrl,
        )?;

        // The symmetric mean is itself a valid non-negative half-breadth
        // (mean of two non-negative nets), so build the base hull from it —
        // this reuses every validated symmetric computation unchanged.
        let mut hull = Hull::new(sym_surface)?;

        // Antisymmetric ∂f_a/∂x on the same spans.
        let xs = a_surface.x_span_indices();
        let zs = a_surface.z_span_indices();
        let (fx_a_coeff, fx_a_rel_error_bound) = compute_fx_coeff(&a_surface, &xs, &zs)?;
        hull.fx_a_coeff = Some(fx_a_coeff);
        hull.fx_coeff_rel_error_bound = hull.fx_coeff_rel_error_bound.max(fx_a_rel_error_bound);
        hull.a_surface = Some(a_surface);

        // Two-sided geometry corrections: wetted surface and the centreplane
        // transverse inertia see each side separately, not twice the mean.
        // (Volume, LCB/VCB and waterplane area/moment depend only on the sum
        // f₊ + f₋ = 2 f_sym and are already correct from the base hull.)
        hull.wetted_surface = wetted_surface_of(&starboard) + wetted_surface_of(&port);
        hull.waterplane_transverse_moment =
            waterplane_transverse_moment_of(&starboard) + waterplane_transverse_moment_of(&port);
        Ok(hull)
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

    /// Derivatives of displacement volume, LCB, and wetted surface with
    /// respect to every physical B-spline control value.
    ///
    /// Volume and first-moment derivatives are polynomial basis integrals and
    /// are evaluated exactly to floating-point roundoff. Wetted-surface
    /// derivatives differentiate the same 24-point-per-span Gauss–Legendre
    /// rule used by [`Self::wetted_surface`]. An asymmetric hull returns
    /// separate port and starboard derivatives. LCB is undefined for a hull
    /// with zero displaced volume, which is reported as an error.
    pub fn constraint_gradients(&self) -> Result<HullConstraintGradients> {
        if self.displaced_volume <= 0.0 {
            return Err(Error::InvalidGeometry(
                "constraint gradients require positive displaced volume".into(),
            ));
        }
        match self.a_surface.as_ref() {
            None => Ok(HullConstraintGradients::Symmetric(
                self.constraint_gradient_for_side(0.0, 2.0),
            )),
            Some(_) => Ok(HullConstraintGradients::Asymmetric {
                port: self.constraint_gradient_for_side(-1.0, 1.0),
                starboard: self.constraint_gradient_for_side(1.0, 1.0),
            }),
        }
    }

    /// Constraint derivative for `f_side = f_sym + camber_sign · f_a`.
    /// `multiplicity` is two for a symmetric net and one for a physical side.
    fn constraint_gradient_for_side(
        &self,
        camber_sign: f64,
        multiplicity: f64,
    ) -> ConstraintGradient {
        let surface = &self.surface;
        let p = surface.degree_x();
        let q = surface.degree_z();
        let nx = surface.n_ctrl_x();
        let nz = surface.n_ctrl_z();
        let xs = surface.x_span_indices();
        let zs = surface.z_span_indices();
        let mut displaced_volume = vec![0.0; nx * nz];
        let mut volume_moment = vec![0.0; nx * nz];
        let mut wetted_surface = vec![0.0; nx * nz];

        let n_volume = (p.max(q) + 2) / 2 + 2;
        let (volume_nodes, volume_weights) = gauss_legendre(n_volume);
        for &sx in &xs {
            let x_start = surface.knots_x()[sx];
            let x_len = surface.knots_x()[sx + 1] - x_start;
            for &sz in &zs {
                let z_start = surface.knots_z()[sz];
                let z_len = surface.knots_z()[sz + 1] - z_start;
                let jacobian = x_len * z_len / 4.0;
                for (ix, &node_x) in volume_nodes.iter().enumerate() {
                    let x = x_start + x_len * (node_x + 1.0) / 2.0;
                    let basis_x = ders_basis(surface.knots_x(), p, sx, x, 0);
                    for (iz, &node_z) in volume_nodes.iter().enumerate() {
                        let z = z_start + z_len * (node_z + 1.0) / 2.0;
                        let basis_z = ders_basis(surface.knots_z(), q, sz, z, 0);
                        let weighted_jacobian =
                            multiplicity * volume_weights[ix] * volume_weights[iz] * jacobian;
                        for (local_x, &value_x) in basis_x[0].iter().enumerate() {
                            let control_x = sx - p + local_x;
                            for (local_z, &value_z) in basis_z[0].iter().enumerate() {
                                let control_z = sz - q + local_z;
                                let index = control_x * nz + control_z;
                                let derivative = weighted_jacobian * value_x * value_z;
                                displaced_volume[index] += derivative;
                                volume_moment[index] += x * derivative;
                            }
                        }
                    }
                }
            }
        }

        let (wetted_nodes, wetted_weights) = gauss_legendre(24);
        for &sx in &xs {
            let x_start = surface.knots_x()[sx];
            let x_len = surface.knots_x()[sx + 1] - x_start;
            for &sz in &zs {
                let z_start = surface.knots_z()[sz];
                let z_len = surface.knots_z()[sz + 1] - z_start;
                let jacobian = x_len * z_len / 4.0;
                for (ix, &node_x) in wetted_nodes.iter().enumerate() {
                    let x = x_start + x_len * (node_x + 1.0) / 2.0;
                    let basis_x = ders_basis(surface.knots_x(), p, sx, x, 1);
                    for (iz, &node_z) in wetted_nodes.iter().enumerate() {
                        let z = z_start + z_len * (node_z + 1.0) / 2.0;
                        let basis_z = ders_basis(surface.knots_z(), q, sz, z, 1);
                        let mean_fx = surface.eval_deriv(x, z, 1, 0);
                        let mean_fz = surface.eval_deriv(x, z, 0, 1);
                        let camber_fx = self
                            .a_surface
                            .as_ref()
                            .map_or(0.0, |camber| camber.eval_deriv(x, z, 1, 0));
                        let camber_fz = self
                            .a_surface
                            .as_ref()
                            .map_or(0.0, |camber| camber.eval_deriv(x, z, 0, 1));
                        let fx = mean_fx + camber_sign * camber_fx;
                        let fz = mean_fz + camber_sign * camber_fz;
                        let area_scale = (1.0 + fx * fx + fz * fz).sqrt();
                        let weighted_jacobian =
                            multiplicity * wetted_weights[ix] * wetted_weights[iz] * jacobian;
                        for (local_x, &value_x) in basis_x[0].iter().enumerate() {
                            let control_x = sx - p + local_x;
                            for (local_z, &value_z) in basis_z[0].iter().enumerate() {
                                let control_z = sz - q + local_z;
                                let index = control_x * nz + control_z;
                                let derivative_fx = basis_x[1][local_x] * value_z;
                                let derivative_fz = value_x * basis_z[1][local_z];
                                wetted_surface[index] += weighted_jacobian
                                    * (fx * derivative_fx + fz * derivative_fz)
                                    / area_scale;
                            }
                        }
                    }
                }
            }
        }

        let lcb_x = volume_moment
            .iter()
            .zip(&displaced_volume)
            .map(|(moment, volume)| (moment - self.lcb_x * volume) / self.displaced_volume)
            .collect();
        ConstraintGradient {
            displaced_volume,
            lcb_x,
            wetted_surface,
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

    pub(crate) fn fx_coeff_rel_error_bound(&self) -> f64 {
        self.fx_coeff_rel_error_bound
    }

    /// Reverse the linear map from surface controls to the local polynomial
    /// coefficients of `∂f/∂x` used by the exact Michell inner integral.
    pub(crate) fn fx_control_adjoint(&self, coeff_adjoint: &[f64]) -> Vec<f64> {
        let surface = &self.surface;
        let (p, q) = (surface.degree_x(), surface.degree_z());
        let xs = surface.x_span_indices();
        let zs = surface.z_span_indices();
        assert_eq!(coeff_adjoint.len(), xs.len() * zs.len() * p * (q + 1));

        let mut factorial = vec![1.0f64; p.max(q) + 2];
        for index in 1..factorial.len() {
            factorial[index] = factorial[index - 1] * index as f64;
        }
        let mut control_adjoint = vec![0.0; surface.control().len()];
        for (isx, &sx) in xs.iter().enumerate() {
            let ndx = ders_basis(surface.knots_x(), p, sx, surface.knots_x()[sx], p);
            for (isz, &sz) in zs.iter().enumerate() {
                let ndz = ders_basis(surface.knots_z(), q, sz, surface.knots_z()[sz], q);
                for a in 0..p {
                    for b in 0..=q {
                        let coefficient = coeff_adjoint
                            [((isx * zs.len() + isz) * p + a) * (q + 1) + b]
                            / (factorial[a] * factorial[b]);
                        for (i, &bx) in ndx[a + 1].iter().enumerate() {
                            let ci = sx - p + i;
                            for (j, &bz) in ndz[b].iter().enumerate() {
                                let cj = sz - q + j;
                                control_adjoint[ci * surface.n_ctrl_z() + cj] +=
                                    coefficient * bx * bz;
                            }
                        }
                    }
                }
            }
        }
        control_adjoint
    }

    /// `∂f_a/∂x` coefficients for the antisymmetric half-beam, or `None` if the
    /// hull is port/starboard symmetric.
    pub(crate) fn fx_a_coeff(&self) -> Option<&[f64]> {
        self.fx_a_coeff.as_deref()
    }

    /// Whether the hull is asymmetric (has an antisymmetric camber part).
    pub fn is_asymmetric(&self) -> bool {
        self.a_surface.is_some()
    }

    /// `∂f_a/∂x` at `(x, z)` in the hull's own coordinates, evaluated from the
    /// antisymmetric surface. Returns 0 for a symmetric hull. This is the
    /// camber-slope forcing of the centreplane lifting solve.
    pub(crate) fn eval_fx_a(&self, x: f64, z: f64) -> f64 {
        self.a_surface
            .as_ref()
            .map_or(0.0, |s| s.eval_deriv(x, z, 1, 0))
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

/// Local polynomial coefficients of `∂f/∂x` per span pair, in the flattened
/// layout documented on [`Hull::fx_coeff`]. Shared by the symmetric and
/// asymmetric constructors so both paths use identical arithmetic.
fn compute_fx_coeff(
    surface: &BSplineSurface,
    xs: &[usize],
    zs: &[usize],
) -> Result<(Vec<f64>, f64)> {
    let p = surface.degree_x();
    let q = surface.degree_z();
    // Factorials up to max degree (degrees are small).
    let mut fact = vec![1.0f64; p.max(q) + 2];
    for i in 1..fact.len() {
        fact[i] = fact[i - 1] * i as f64;
    }
    let nsz = zs.len();
    let mut fx_coeff = vec![0.0f64; xs.len() * nsz * p * (q + 1)];
    for (isx, &sx) in xs.iter().enumerate() {
        for (isz, &sz) in zs.iter().enumerate() {
            let d = surface.corner_partials(sx, sz);
            if d.iter().flatten().any(|value| !value.is_finite()) {
                return Err(Error::Unsupported(
                    "non-finite span derivative encountered while constructing resistance coefficients"
                        .into(),
                ));
            }
            for a in 0..p {
                for b in 0..=q {
                    // f = Σ D[a][b]/(a! b!) X^a Z^b  =>
                    // fx coefficient of X^a Z^b is D[a+1][b]/(a! b!).
                    let coefficient = d[a + 1][b] / (fact[a] * fact[b]);
                    if !coefficient.is_finite() {
                        return Err(Error::Unsupported(
                            "non-finite local resistance coefficient encountered".into(),
                        ));
                    }
                    fx_coeff[((isx * nsz + isz) * p + a) * (q + 1) + b] = coefficient;
                }
            }
        }
    }
    let reconstruction_rel_error = span_reconstruction_rel_error(surface, xs, zs, &fx_coeff);
    if !reconstruction_rel_error.is_finite()
        || reconstruction_rel_error > MAX_SPAN_RECONSTRUCTION_REL_ERROR
    {
        return Err(Error::Unsupported(format!(
            "local resistance coefficients failed the independent span-reconstruction gate: \
             relative residual {reconstruction_rel_error:.3e} exceeds \
             {MAX_SPAN_RECONSTRUCTION_REL_ERROR:.3e}"
        )));
    }
    Ok((
        fx_coeff,
        VALIDATED_COEFFICIENT_RESISTANCE_REL_ERROR.max(4.0 * reconstruction_rel_error),
    ))
}

fn span_reconstruction_rel_error(
    surface: &BSplineSurface,
    xs: &[usize],
    zs: &[usize],
    coefficients: &[f64],
) -> f64 {
    let p = surface.degree_x();
    let q = surface.degree_z();
    let (nodes_x, _) = gauss_legendre(p + 1);
    let (nodes_z, _) = gauss_legendre(q + 1);
    let mut max_reference = 0.0f64;
    let mut max_residual = 0.0f64;
    let control_scale = surface
        .control()
        .iter()
        .fold(0.0f64, |scale, value| scale.max(value.abs()));
    let min_x_span = xs
        .iter()
        .map(|&span| surface.knots_x()[span + 1] - surface.knots_x()[span])
        .fold(f64::INFINITY, f64::min);
    let characteristic_slope = control_scale / min_x_span;
    for (span_x_index, &span_x) in xs.iter().enumerate() {
        let x0 = surface.knots_x()[span_x];
        let dx = surface.knots_x()[span_x + 1] - x0;
        for (span_z_index, &span_z) in zs.iter().enumerate() {
            let z0 = surface.knots_z()[span_z];
            let dz = surface.knots_z()[span_z + 1] - z0;
            let offset = (span_x_index * zs.len() + span_z_index) * p * (q + 1);
            let span_coefficients = &coefficients[offset..offset + p * (q + 1)];
            for &node_x in &nodes_x {
                let local_x = dx * (node_x + 1.0) / 2.0;
                let x = x0 + local_x;
                for &node_z in &nodes_z {
                    let local_z = dz * (node_z + 1.0) / 2.0;
                    let z = z0 + local_z;
                    let reference = surface.eval_deriv(x, z, 1, 0);
                    let reconstructed =
                        evaluate_local_polynomial(span_coefficients, p, q, local_x, local_z);
                    max_reference = max_reference.max(reference.abs());
                    max_residual = max_residual.max((reconstructed - reference).abs());
                }
            }
        }
    }
    let scale = max_reference.max(characteristic_slope);
    if scale == 0.0 {
        if max_residual == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        max_residual / scale
    }
}

fn evaluate_local_polynomial(coefficients: &[f64], p: usize, q: usize, x: f64, z: f64) -> f64 {
    (0..p).rev().fold(0.0, |sum_x, a| {
        let row = &coefficients[a * (q + 1)..(a + 1) * (q + 1)];
        let sum_z = row
            .iter()
            .rev()
            .fold(0.0, |sum_z, coefficient| sum_z * z + coefficient);
        sum_x * x + sum_z
    })
}

fn validate_resistance_degree(surface: &BSplineSurface) -> Result<()> {
    let degree_x = surface.degree_x();
    let degree_z = surface.degree_z();
    if degree_x > MAX_SUPPORTED_SPLINE_DEGREE || degree_z > MAX_SUPPORTED_SPLINE_DEGREE {
        return Err(Error::Unsupported(format!(
            "resistance evaluation supports spline degrees up to {MAX_SUPPORTED_SPLINE_DEGREE} \
             in each direction; got degree ({degree_x}, {degree_z})"
        )));
    }
    Ok(())
}

/// One-sided wetted surface `∬ √(1 + fx² + fz²) dx dz` of a half-breadth
/// surface (thin-ship projection, no girth correction).
fn wetted_surface_of(surface: &BSplineSurface) -> f64 {
    let xs = surface.x_span_indices();
    let zs = surface.z_span_indices();
    let (xw, ww) = gauss_legendre(24);
    let mut wetted = 0.0;
    for &sxi in &xs {
        let x_start = surface.knots_x()[sxi];
        let x_len = surface.knots_x()[sxi + 1] - x_start;
        for &szi in &zs {
            let z_start = surface.knots_z()[szi];
            let z_len = surface.knots_z()[szi + 1] - z_start;
            let jac = x_len / 2.0 * (z_len / 2.0);
            for (i, &xi) in xw.iter().enumerate() {
                let x = x_start + x_len * (xi + 1.0) / 2.0;
                for (j, &zj) in xw.iter().enumerate() {
                    let z = z_start + z_len * (zj + 1.0) / 2.0;
                    let fx = surface.eval_deriv(x, z, 1, 0);
                    let fz = surface.eval_deriv(x, z, 0, 1);
                    wetted += ww[i] * ww[j] * jac * (1.0 + fx * fx + fz * fz).sqrt();
                }
            }
        }
    }
    wetted
}

/// One-sided waterplane transverse inertia about the centreplane,
/// `∫ (1/3) f(x, 0)³ dx`.
fn waterplane_transverse_moment_of(surface: &BSplineSurface) -> f64 {
    let p = surface.degree_x();
    let (z0, _) = surface.z_domain();
    let xs = surface.x_span_indices();
    let (xq, wq) = gauss_legendre(3 * p / 2 + 2);
    let mut iyy = 0.0;
    for &sxi in &xs {
        let x_start = surface.knots_x()[sxi];
        let x_len = surface.knots_x()[sxi + 1] - x_start;
        let jac = x_len / 2.0;
        for (i, &xi) in xq.iter().enumerate() {
            let x = x_start + x_len * (xi + 1.0) / 2.0;
            let f = surface.eval(x, z0.max(0.0));
            iyy += wq[i] * jac * (1.0 / 3.0) * f * f * f;
        }
    }
    iyy
}

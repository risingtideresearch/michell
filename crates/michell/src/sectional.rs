//! **Sectional** hull transforms: the Michell and near-field wavenumber
//! transforms evaluated from a hull's stations rather than from a lofted
//! half-breadth spline `f(x, z)`.
//!
//! Every quantity the wave and sinkage/trim integrals need is a transform of
//! the form `∬ f(x, z) e^{−κz} e^{ik_x x} dx dz` (or its `∂f/∂x` and
//! `(x − x_c) ∂f/∂x` variants). Write it as
//!
//! ```text
//! ∫ Z(x; κ) e^{ik_x x} dx,      Z(x; κ) = ∫ f(x, z) e^{−κz} dz,
//! ```
//!
//! and do the depth integral **per station**, along the station's own section
//! curve. The keel is then just the end of that curve: the crease a lofted
//! `f(x, z)` has along a rockered keel line, and the square-root closure of a
//! round bilge onto it, are both inside a 1-D integral that absorbs them, and
//! `Z(x; κ)` is smooth in `x` wherever the hull's sections vary smoothly.
//!
//! What keeps this affordable is that the split preserves the separability
//! the near-field quadrature depends on. Everything that depends on `κ` is
//! one pass over the stations' depth nodes (`Σ_j w_j e^{−κ z_j}` per station)
//! followed by one banded solve for the B-spline interpolant of `Z(·; κ)`
//! along `x`; each `k_x` at that `κ` then costs one closed-form
//! oscillatory-moment call per x-span, exactly as the exact B-spline kernel's does.
//!
//! Hulls are cut from CAD through [`crate::source::HullSource`]. In tests,
//! `SectionalHull::from_hull` builds one from a B-spline half-breadth
//! surface — stations at the Greville points of its own x knots, where the
//! interpolant reproduces `Z(x; κ)` exactly — which is the reference
//! harness: any difference from the exact B-spline kernel
//! (`crate::michell::InnerIntegral`, a test oracle) is the depth
//! quadrature's.

use crate::bspline::{ders_basis, find_span};
use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::float::{DynamicLoad, FleetState};
use crate::friction::{viscous_resistance_for, ViscousOptions, ViscousResistance};
#[cfg(test)]
use crate::hull::Hull;
use crate::michell::Placement;
use crate::michell::{
    add_transom_step, appendage_source, hollow_shape_moment, run_outer, MemberWave, NearFieldKernel, OuterParams, SquatTransforms,
    TransomClosure, WaveOptions, WaveResistance,
};
use crate::moments::{osc_moments, C64};
use crate::quadrature::gauss_legendre;
use crate::squat::{integrate_force, DynamicForce, Fleet, Member, SquatOptions};
use crate::MultihullResistance;
use std::f64::consts::PI;

/// Below this, `e^{−κz}` cannot matter against the shallowest node's weight
/// (the same floor the exact B-spline kernel drops z-spans at).
const DECAY_EXPONENT_FLOOR: f64 = 46.0; // e^{-46} ≈ 1e-20

/// One station's depth integral `Z(κ) = ∫ f(x_i, z) e^{−κz} dz`, as a
/// quadrature that costs one exponential per node at each `κ`.
#[derive(Debug, Clone)]
pub struct SectionNodes {
    nodes: Nodes,
    /// Half-beam at the waterline, `f(x_i, 0)`.
    waterline: f64,
    /// Depth of the section's lowest point.
    depth: f64,
    /// `(half-beam, depth)` at each quadrature node, in node order: where
    /// the depth integral evaluates the section.
    outline: Vec<(f64, f64)>,
    /// The section curve the quadrature integrates, sampled densely and
    /// evenly from the waterline (or the section's top) round to the keel —
    /// for drawing it and for measuring its girth, where joining the
    /// quadrature nodes (graded toward the waterline, sparse at the keel)
    /// would cut corners.
    curve: Vec<(f64, f64)>,
}

/// Points per section in [`SectionNodes::curve`].
const CURVE_POINTS: usize = 129;

#[derive(Debug, Clone)]
enum Nodes {
    /// `Σ_j w_j e^{−κ z_j}`, nodes sorted by depth: a half-beam sampled as
    /// a function of depth.
    Depth { z: Vec<f64>, w: Vec<f64> },
    /// `e^{−κ z₀} Σ_k w_k F(κ d sin θ_k, R_k)` with `F(a, R) = ∫_0^R e^{−ar} r dr`
    /// (`sin` stores `d sin θ_k`, `w` carries the scales `b·d`):
    /// the section as a region swept by rays from its top centreplane point
    /// `(0, z₀)`, the ray at angle `θ` below the horizontal reaching the shell
    /// at distance `R(θ)`. The depth integral is then an area integral in
    /// polar coordinates whose radial part is closed form, and `R(θ)` is
    /// smooth right into the keel — where a half-beam as a function of depth
    /// closes like a square root on any round bilge.
    Polar {
        z0: f64,
        sin: Vec<f64>,
        r: Vec<f64>,
        w: Vec<f64>,
    },
}

/// `∫_0^R e^{−ar} r dr = R²·(1 − e^{−x}(1 + x))/x²`, `x = aR`, by series where
/// the closed form cancels.
fn polar_radial(a: f64, r: f64) -> f64 {
    let x = a * r;
    if x < 0.1 {
        // Σ_n (−1)ⁿ (n+1) xⁿ/(n+2)!
        let (mut term, mut sum, mut fact) = (1.0f64, 0.5f64, 2.0f64);
        for n in 1..12 {
            term *= -x;
            fact *= (n + 2) as f64;
            sum += (n + 1) as f64 * term / fact;
        }
        r * r * sum
    } else {
        r * r * (1.0 - (-x).exp() * (1.0 + x)) / (x * x)
    }
}

impl SectionNodes {
    /// Nodes for a section given as a half-beam function of depth on
    /// `[0, depth]`: Gauss–Legendre panels between the `breaks` (kinks,
    /// chines, knot lines — anywhere the half-beam is not smooth), further
    /// split geometrically toward the waterline so that `e^{−κz}` stays
    /// resolved when a large `κ` confines it to a sliver under the surface.
    pub fn from_depth_function(
        depth: f64,
        breaks: &[f64],
        half_beam: impl Fn(f64) -> f64,
        opts: &DepthQuadrature,
    ) -> SectionNodes {
        let cuts = graded_cuts(depth, breaks, opts);
        let (gx, gw) = gauss_legendre(opts.points);
        let (mut z, mut w) = (Vec::new(), Vec::new());
        let mut outline = vec![(half_beam(0.0), 0.0)];
        for c in cuts.windows(2) {
            let (a, b) = (c[0], c[1]);
            let half = 0.5 * (b - a);
            for (&t, &wt) in gx.iter().zip(&gw) {
                let zj = a + half * (t + 1.0);
                let h = half_beam(zj);
                z.push(zj);
                w.push(half * wt * h);
                outline.push((h, zj));
            }
        }
        outline.push((half_beam(depth), depth));
        // Evenly in depth, with every break (a chine, a knot line) on it.
        let mut zs: Vec<f64> = (0..CURVE_POINTS)
            .map(|k| depth * k as f64 / (CURVE_POINTS - 1) as f64)
            .chain(breaks.iter().copied().filter(|&b| b > 0.0 && b < depth))
            .collect();
        zs.sort_by(f64::total_cmp);
        let curve = zs.iter().map(|&z| (half_beam(z), z)).collect();
        SectionNodes {
            nodes: Nodes::Depth { z, w },
            waterline: half_beam(0.0),
            depth,
            outline,
            curve,
        }
    }

    /// Nodes for a section swept by rays from `(0, z0)` on the centreplane,
    /// in the section's own proportions: with half-beam scaled by `beam` and
    /// depth by `depth`, `radius(θ)` is the scaled distance to the shell
    /// along the scaled ray `θ` below the horizontal, `θ ∈ [0, π/2]` (`π/2`
    /// runs down the centreplane to the keel). Sweeping the scaled plane
    /// keeps rays spread over a thin section (a fine entry is millimetres
    /// wide and the full draft deep), where rays even in physical angle
    /// would crowd its whole outline into a sliver of angle beside the
    /// keel. Requires the section to be star-shaped about `(0, z0)`. Graded
    /// toward `θ = 0`, where a large `κ` concentrates the integrand against
    /// the surface.
    pub fn from_polar(
        z0: f64,
        beam: f64,
        depth: f64,
        radius: impl Fn(f64) -> f64,
        opts: &DepthQuadrature,
    ) -> SectionNodes {
        let half_pi = std::f64::consts::FRAC_PI_2;
        let cuts = graded_cuts(half_pi, &[], opts);
        let (gx, gw) = gauss_legendre(opts.points);
        let (mut sin, mut r, mut w) = (Vec::new(), Vec::new(), Vec::new());
        let at = |th: f64, rr: f64| (beam * rr * th.cos(), z0 + depth * rr * th.sin());
        let mut outline = vec![at(0.0, radius(0.0))];
        for c in cuts.windows(2) {
            let (a, b) = (c[0], c[1]);
            let half = 0.5 * (b - a);
            for (&t, &wt) in gx.iter().zip(&gw) {
                let th = a + half * (t + 1.0);
                let rr = radius(th);
                // Physical depth along the scaled ray is depth·r·sin θ.
                sin.push(depth * th.sin());
                r.push(rr);
                w.push(half * wt * beam * depth);
                outline.push(at(th, rr));
            }
        }
        outline.push(at(half_pi, radius(half_pi)));
        // Evenly in the scaled ray angle: evenly round a section of any
        // proportions.
        let curve = (0..CURVE_POINTS)
            .map(|k| {
                let th = half_pi * k as f64 / (CURVE_POINTS - 1) as f64;
                at(th, radius(th))
            })
            .collect();
        SectionNodes {
            nodes: Nodes::Polar { z0, sin, r, w },
            waterline: if z0 == 0.0 { beam * radius(0.0) } else { 0.0 },
            depth: z0 + depth * radius(half_pi),
            outline,
            curve,
        }
    }

    /// A station with no wetted section (beyond the hull's ends).
    pub fn empty() -> SectionNodes {
        SectionNodes {
            nodes: Nodes::Depth {
                z: Vec::new(),
                w: Vec::new(),
            },
            waterline: 0.0,
            depth: 0.0,
            outline: Vec::new(),
            curve: Vec::new(),
        }
    }

    /// The section at its quadrature nodes, `(half-beam, depth)` from the
    /// waterline (or the section's top) round to the keel: exactly what the
    /// depth integral integrates, for display.
    pub fn outline(&self) -> &[(f64, f64)] {
        &self.outline
    }

    /// The section curve the quadrature integrates, sampled densely and
    /// evenly, `(half-beam, depth)` from the waterline (or top) to the keel.
    pub fn curve(&self) -> &[(f64, f64)] {
        &self.curve
    }

    /// `Z(κ)`.
    fn integrate(&self, kappa: f64) -> f64 {
        match &self.nodes {
            Nodes::Depth { z, w } => {
                let mut s = 0.0;
                for (&z, &w) in z.iter().zip(w) {
                    let e = kappa * z;
                    if e > DECAY_EXPONENT_FLOOR {
                        break;
                    }
                    s += w * (-e).exp();
                }
                s
            }
            Nodes::Polar { z0, sin, r, w } => {
                let e0 = kappa * z0;
                if e0 > DECAY_EXPONENT_FLOOR {
                    return 0.0;
                }
                let s: f64 = sin
                    .iter()
                    .zip(r)
                    .zip(w)
                    .map(|((&sn, &r), &w)| w * polar_radial(kappa * sn, r))
                    .sum();
                s * (-e0).exp()
            }
        }
    }

    /// Number of quadrature nodes.
    pub fn len(&self) -> usize {
        match &self.nodes {
            Nodes::Depth { z, .. } => z.len(),
            Nodes::Polar { r, .. } => r.len(),
        }
    }

    /// True when the section has no nodes.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Panel edges on `[0, extent]`: the `breaks`, plus a geometric sequence
/// toward 0 with ratio `opts.grading` down to `opts.finest · extent`.
fn graded_cuts(extent: f64, breaks: &[f64], opts: &DepthQuadrature) -> Vec<f64> {
    let mut cuts = vec![0.0, extent];
    let mut t = extent;
    while t > opts.finest * extent {
        t *= opts.grading;
        cuts.push(t);
    }
    cuts.extend(breaks.iter().copied().filter(|&b| b > 0.0 && b < extent));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|a, b| (*a - *b).abs() <= 1e-14 * extent);
    cuts
}

/// How finely each station's depth integral is resolved.
#[derive(Debug, Clone, Copy)]
pub struct DepthQuadrature {
    /// Gauss–Legendre points per panel.
    pub points: usize,
    /// Ratio between successive panel edges toward the waterline.
    pub grading: f64,
    /// Shallowest panel edge, as a fraction of the section depth.
    pub finest: f64,
}

impl Default for DepthQuadrature {
    fn default() -> Self {
        // 8 points on panels graded by 0.35 toward the waterline hold the
        // transforms to ~1e-11 of the exact B-spline kernel's exact values on the
        // Wigley (one z-span: the deepest panel is ~0.4 T long) and ~1e-12
        // on a CAD import, for κ up to ~1e4·ν. Six points leave ~1e-8.
        DepthQuadrature {
            points: 8,
            grading: 0.35,
            finest: 1e-7,
        }
    }
}

/// One non-empty knot span in one direction.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Span {
    /// Physical coordinate of the span start.
    pub start: f64,
    /// Span length (> 0).
    pub len: f64,
}

/// A transom whose area is under this fraction of the hull's maximum section
/// area is not reported: it is hydrodynamically negligible, and at a closing
/// stern indistinguishable from the geometry's own noise.
pub(crate) const TRANSOM_AREA_REL: f64 = 1e-3;

/// The aft-end section of a hull whose half-breadth does not close there — a
/// **transom**.
///
/// The geometry contract puts the bow at the high-`x` end, so the transom, if
/// there is one, is the section at the hull's lowest `x`. `f_T(z) = f(x_T, z)`
/// is the aft end station's section.
///
/// Presence alone is not a warning: a transom clear of the water is simply a
/// closed hull as far as the wave integral is concerned. What matters is
/// [`Transom::depth`] and the area ratio against
/// [`SectionalHull::max_section_area`].
#[derive(Debug, Clone)]
pub struct Transom {
    /// Station of the transom — the aft end of the hull's x-domain [m].
    pub x: f64,
    /// Immersion depth [m], as the **equivalent rectangle**:
    /// `A_T / (2 · max_z f_T)` — the depth of a rectangle of the transom's
    /// widest beam carrying the same immersed area. Exact for a rectangular
    /// transom; `T/2` for one tapering linearly to the keel.
    ///
    /// Deliberately not a level crossing of `f_T`: an area measure is
    /// insensitive to a thin tail of beam running down the draft, and it is
    /// also the scale a closure wants, since what sets the hollow is how much
    /// water has to fill in behind the transom, not where its edge sits.
    pub depth: f64,
    /// Immersed transom area, `2∫₀^T f_T(z) dz` [m²] — both sides.
    pub area: f64,
    /// Half-beam at the waterline, `f_T(0)` [m].
    pub half_beam: f64,
    /// The exact B-spline oracle's transom section: per-z-span polynomial
    /// coefficients `f_T(z) = Σ_b c[b] (z − z0_sz)^b`, flattened as
    /// `[sz * (q + 1) + b]`. Empty on a sectional hull cut from geometry.
    #[cfg(test)]
    pub(crate) coeff: Vec<f64>,
}

/// A hull as stations along `x`, each with its section's depth quadrature,
/// and the B-spline space along `x` its depth integrals are interpolated in.
#[derive(Debug, Clone)]
pub struct SectionalHull {
    p: usize,
    n: usize,
    sections: Vec<SectionNodes>,
    /// Banded LU of the collocation matrix (interpolation at the stations).
    lu: BandLu,
    spans: Vec<Span>,
    /// First basis index active on each span.
    span_first: Vec<usize>,
    /// Per span, `(p+1) × (p+1)`: `N_{first+r}^{(a)}(x_s) / a!` at `[a][r]` —
    /// maps B-spline coefficients to the span's local power basis.
    taylor: Vec<f64>,
    /// Waterline `f(x, 0)` and `∂f/∂x(x, 0)` in each span's power basis.
    wl_f: Vec<f64>,
    wl_fx: Vec<f64>,
    x_center: f64,
    length: f64,
    draft: f64,
    volume: f64,
    xs: Vec<f64>,
    lcb_x: f64,
    waterplane: [f64; 3],
    wetted_surface: f64,
    /// The aft end station's section, when it is a real (immersed) one: a
    /// transom, closed by [`TransomClosure`] as on any hull.
    transom: Option<Transom>,
}

/// A sectional hull's depth integrals at one `κ`, in each x-span's local
/// power basis (the layout the exact B-spline kernel's `ZContracted` uses).
#[derive(Debug, Clone, Default)]
pub struct SectionalContracted {
    g_f: Vec<f64>,
    g_fx: Vec<f64>,
    any: bool,
    /// The aft end station's depth integral: the transom section's
    /// z-factor for the closure.
    z_t: f64,
}

impl SectionalHull {
    /// A sectional hull from stations `xs` (strictly increasing) and their
    /// sections, interpolated along `x` by a clamped B-spline of degree `p`
    /// on `knots` (which must satisfy Schoenberg–Whitney against `xs`: one
    /// station per basis function, each inside its support).
    pub fn new(p: usize, knots: Vec<f64>, xs: &[f64], sections: Vec<SectionNodes>) -> Result<Self> {
        let n = knots.len() - p - 1;
        if xs.len() != n || sections.len() != n {
            return Err(Error::InvalidInput(format!(
                "sectional hull: {} stations and {} sections for {n} basis functions",
                xs.len(),
                sections.len()
            )));
        }
        // Collocation: row i holds the basis functions alive at station i.
        let mut lu = BandLu::new(n, p);
        for (i, &x) in xs.iter().enumerate() {
            let span = find_span(&knots, p, n, x);
            let d = ders_basis(&knots, p, span, x, 0);
            for r in 0..=p {
                lu.set(i, span - p + r, d[0][r]);
            }
        }
        if !lu.factor() {
            return Err(Error::InvalidInput(
                "sectional hull: stations do not interpolate the x knots \
                 (Schoenberg–Whitney violated)"
                    .into(),
            ));
        }
        let mut spans = Vec::new();
        let mut span_first = Vec::new();
        let mut taylor = Vec::new();
        for s in p..n {
            let (a, b) = (knots[s], knots[s + 1]);
            if b <= a {
                continue;
            }
            spans.push(Span {
                start: a,
                len: b - a,
            });
            span_first.push(s - p);
            let d = ders_basis(&knots, p, s, a, p);
            let mut fact = 1.0;
            for k in 0..=p {
                if k > 0 {
                    fact *= k as f64;
                }
                for r in 0..=p {
                    taylor.push(d[k][r] / fact);
                }
            }
        }
        let mut hull = SectionalHull {
            p,
            n,
            sections,
            lu,
            spans,
            span_first,
            taylor,
            wl_f: Vec::new(),
            wl_fx: Vec::new(),
            x_center: 0.5 * (knots[p] + knots[n]),
            length: knots[n] - knots[p],
            draft: 0.0,
            volume: 0.0,
            xs: xs.to_vec(),
            lcb_x: 0.0,
            waterplane: [0.0; 3],
            wetted_surface: 0.0,
            transom: None,
        };
        hull.draft = hull.sections.iter().map(|s| s.depth).fold(0.0, f64::max);
        // Displaced volume 2∬ f = 2∫ Z(x; 0) dx and its x-moment, exactly
        // from the interpolant of the section areas.
        let mut z0 = SectionalContracted::default();
        hull.contract(0.0, &mut z0);
        hull.volume = 2.0 * hull.x_moment(&z0.g_f, 0);
        hull.lcb_x = if hull.volume > 0.0 {
            2.0 * hull.x_moment(&z0.g_f, 1) / hull.volume
        } else {
            0.0
        };
        let wl: Vec<f64> = hull.sections.iter().map(|s| s.waterline).collect();
        let (f, fx) = hull.to_power_basis(wl);
        hull.wl_f = f;
        hull.wl_fx = fx;
        // Waterplane: A_w = ∫ 2 f(x, 0) dx and its first and second moments
        // about x = 0, from the waterline interpolant.
        hull.waterplane = [0, 1, 2].map(|k| 2.0 * hull.x_moment(&hull.wl_f, k));
        hull.wetted_surface = hull.strip_area();
        hull.transom = hull.detect_transom();
        Ok(hull)
    }

    /// The reference construction: sections of a B-spline hull at the
    /// Greville points of its own x knots, depth-integrated panel by panel
    /// between its z knots. `Z(x; κ)` is then exactly in the interpolating
    /// space, so the only approximation is the depth quadrature.
    #[cfg(test)]
    pub(crate) fn from_hull(hull: &Hull, opts: &DepthQuadrature) -> Result<Self> {
        let s = hull.surface();
        let (p, n) = (s.degree_x(), s.n_ctrl_x());
        let kx = s.knots_x();
        let xs: Vec<f64> = (0..n)
            .map(|i| kx[i + 1..=i + p].iter().sum::<f64>() / p.max(1) as f64)
            .collect();
        let depth = s.z_domain().1;
        let sections = xs
            .iter()
            .map(|&x| SectionNodes::from_depth_function(depth, s.knots_z(), |z| s.eval(x, z), opts))
            .collect();
        // The spline's own transom description (its depth measure samples
        // the spline, not the quadrature nodes), so the harness checks the
        // closure machinery exactly.
        let mut sec = SectionalHull::new(p, kx.to_vec(), &xs, sections)?;
        sec.transom = hull.transom().cloned();
        Ok(sec)
    }

    /// The transom, if the aft end station carries a real section: immersed area
    /// above `TRANSOM_AREA_REL` of the largest section's — with its depth as
    /// the same equivalent rectangle, `A_T / (2 · max f_T)`. The geometry
    /// contract puts the bow at high x, so this is the lowest-x station.
    ///
    /// Only the aft end is closed. A truncated **bow**
    /// (a forward end station with a real section) is left as it is: its
    /// section steps to nothing past the end of the interpolant, which the
    /// wave kernel reads as an open, bluff forward face — no virtual
    /// appendage is added ahead of it.
    fn detect_transom(&self) -> Option<Transom> {
        let aft = self.sections.first()?;
        let area = 2.0 * aft.integrate(0.0);
        let max_section = self
            .sections
            .iter()
            .map(|s| 2.0 * s.integrate(0.0))
            .fold(0.0, f64::max);
        if !(area > TRANSOM_AREA_REL * max_section) {
            return None;
        }
        let peak = aft
            .outline()
            .iter()
            .map(|p| p.0)
            .fold(aft.waterline, f64::max);
        let depth = if peak > 0.0 {
            (area / (2.0 * peak)).min(self.draft)
        } else {
            0.0
        };
        Some(Transom {
            x: self.xs[0],
            depth,
            area,
            half_beam: aft.waterline,
            #[cfg(test)]
            coeff: Vec::new(),
        })
    }

    /// The transom the aft end presents, if any.
    pub fn transom(&self) -> Option<&Transom> {
        self.transom.as_ref()
    }

    /// Largest immersed section area `2∫ f dz` over the stations [m²] — the
    /// reference a transom's area is judged against.
    pub fn max_section_area(&self) -> f64 {
        self.sections
            .iter()
            .map(|s| 2.0 * s.integrate(0.0))
            .fold(0.0, f64::max)
    }

    /// `∫ x^k g(x) dx` over the hull for `g` given in each span's local power
    /// basis (`[s·(p+1) + a]`, local coordinate `t = x − start`): exact, by
    /// expanding `(start + t)^k` binomially.
    fn x_moment(&self, g: &[f64], k: usize) -> f64 {
        let p = self.p;
        let binom = |n: usize, r: usize| -> f64 {
            (0..r).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
        };
        self.spans
            .iter()
            .enumerate()
            .map(|(s, sx)| {
                let mut total = 0.0;
                for a in 0..=p {
                    let c = g[s * (p + 1) + a];
                    for j in 0..=k {
                        // ∫_0^len start^(k−j) t^(j+a) dt
                        let e = j + a + 1;
                        total +=
                            c * binom(k, j) * sx.start.powi((k - j) as i32) * sx.len.powi(e as i32)
                                / e as f64;
                    }
                }
                total
            })
            .sum()
    }

    /// Wetted surface of both sides [m²]: the shell between neighbouring
    /// stations as a strip of triangles joining their outlines node for
    /// node, which carries the slope along x a projected area would miss.
    /// An end station without a section closes its strip to a point.
    fn strip_area(&self) -> f64 {
        let tri = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
            let (u, v) = (
                [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
            );
            let w = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt()
        };
        let mut area = 0.0;
        for i in 0..self.xs.len().saturating_sub(1) {
            let (xa, xb) = (self.xs[i], self.xs[i + 1]);
            let (oa, ob) = (self.sections[i].curve(), self.sections[i + 1].curve());
            let n = oa.len().max(ob.len());
            if n < 2 {
                continue;
            }
            // Node k of an outline, stretched over n nodes (a point if empty).
            let at = |o: &[(f64, f64)], x: f64, k: usize| -> [f64; 3] {
                if o.is_empty() {
                    return [x, 0.0, 0.0];
                }
                let j = (k * (o.len() - 1)) / (n - 1);
                [x, o[j].0, o[j].1]
            };
            for k in 0..n - 1 {
                let (a0, a1) = (at(oa, xa, k), at(oa, xa, k + 1));
                let (b0, b1) = (at(ob, xb, k), at(ob, xb, k + 1));
                area += tri(a0, b0, b1) + tri(a0, b1, a1);
            }
        }
        2.0 * area
    }

    /// B-spline interpolant of per-station values, as power coefficients of
    /// the value (`[s·(p+1) + a]`) and its x-derivative (`[s·p + a]`) on
    /// each span.
    fn to_power_basis(&self, mut v: Vec<f64>) -> (Vec<f64>, Vec<f64>) {
        let mut f = vec![0.0; self.spans.len() * (self.p + 1)];
        let mut fx = vec![0.0; self.spans.len() * self.p];
        self.lu.solve(&mut v);
        self.power_into(&v, &mut f, &mut fx);
        (f, fx)
    }

    fn power_into(&self, c: &[f64], f: &mut [f64], fx: &mut [f64]) {
        let p = self.p;
        for (s, &first) in self.span_first.iter().enumerate() {
            let t = &self.taylor[s * (p + 1) * (p + 1)..(s + 1) * (p + 1) * (p + 1)];
            for a in 0..=p {
                let row = &t[a * (p + 1)..(a + 1) * (p + 1)];
                f[s * (p + 1) + a] = row
                    .iter()
                    .zip(&c[first..=first + p])
                    .map(|(m, c)| m * c)
                    .sum();
            }
            for a in 0..p {
                fx[s * p + a] = (a + 1) as f64 * f[s * (p + 1) + a + 1];
            }
        }
    }

    /// Everything that depends on `κ`: each station's depth integral, and
    /// their interpolant along `x` in each span's power basis.
    pub fn contract(&self, kappa: f64, out: &mut SectionalContracted) {
        let mut v: Vec<f64> = self.sections.iter().map(|s| s.integrate(kappa)).collect();
        out.any = v.iter().any(|&x| x != 0.0);
        out.z_t = v.first().copied().unwrap_or(0.0);
        out.g_f.resize(self.spans.len() * (self.p + 1), 0.0);
        out.g_fx.resize(self.spans.len() * self.p, 0.0);
        if !out.any {
            out.g_f.fill(0.0);
            out.g_fx.fill(0.0);
            return;
        }
        self.lu.solve(&mut v);
        self.power_into(&v, &mut out.g_f, &mut out.g_fx);
    }

    /// The source free-wave amplitude `∬ ∂f/∂x e^{−κz} e^{ik_x(x−x_c)}` at
    /// `λ` (the exact B-spline kernel's `InnerIntegral::eval` convention, no transom
    /// closure).
    pub fn amplitude(&self, nu: f64, lambda: f64, scratch: &mut SectionalContracted) -> C64 {
        self.amplitude_closed(nu, lambda, TransomClosure::None, scratch)
    }

    /// [`SectionalHull::amplitude`] with the transom closed by `closure`: the
    /// virtual appendage `f_v = f_T(z)·φ(s)` over the hollow behind the
    /// transom, exactly as the exact B-spline kernel adds it — its z-factor is the
    /// aft end station's own depth integral.
    pub fn amplitude_closed(
        &self,
        nu: f64,
        lambda: f64,
        closure: TransomClosure,
        scratch: &mut SectionalContracted,
    ) -> C64 {
        let kx = nu * lambda;
        self.contract(nu * lambda * lambda, scratch);
        if !scratch.any {
            return C64::ZERO;
        }
        let closing = self.transom_term(kx, nu, closure, scratch.z_t);
        let p = self.p;
        let mut xm = Vec::with_capacity(p + 1);
        let mut f = C64::ZERO;
        for (s, sx) in self.spans.iter().enumerate() {
            osc_moments(kx, sx.len, p - 1, &mut xm);
            let phase = C64::cis(kx * (sx.start - self.x_center));
            let mut sum = C64::ZERO;
            for (a, &m) in xm.iter().enumerate() {
                sum = sum + m.scale(scratch.g_fx[s * p + a]);
            }
            f = f + phase * sum;
        }
        f + closing
    }

    /// The transom appendage's free-wave amplitude (kernel convention): the
    /// exact `InnerIntegral::transom_term` with the transom section's
    /// z-factor `z_t` taken from the aft end station.
    fn transom_term(&self, kx: f64, nu: f64, closure: TransomClosure, z_t: f64) -> C64 {
        let Some(tr) = &self.transom else {
            return C64::ZERO;
        };
        let Some(lv) = closure.hollow_length(tr.depth, nu) else {
            return C64::ZERO;
        };
        let mut sm = Vec::with_capacity(3);
        let shape = hollow_shape_moment(-kx * lv, &mut sm);
        let phase = C64::cis(kx * (tr.x - self.x_center));
        C64::ZERO - (phase * shape).scale(z_t)
    }

    /// The transom's share of the near-field transforms (kernel convention),
    /// as the exact `InnerIntegral::add_transom_transforms`: the step in the
    /// weights, the closing appendage in the sources (see
    /// [`SquatTransforms`]).
    fn add_transom_transforms(
        &self,
        kx: f64,
        nu: f64,
        closure: TransomClosure,
        z_t: f64,
        t: &mut SquatTransforms,
    ) {
        let Some(tr) = &self.transom else {
            return;
        };
        let dx_t = tr.x - self.x_center;
        let phase = C64::cis(kx * dx_t);
        add_transom_step(phase, dx_t, z_t, tr.half_beam, t);
        let Some(lv) = closure.hollow_length(tr.depth, nu) else {
            return;
        };
        let mut m = Vec::with_capacity(4);
        osc_moments(-kx * lv, 1.0, 3, &mut m);
        t.q_src = t.q_src + appendage_source(phase, &m).scale(z_t);
    }

    /// The six near-field transforms at `k_x` from a contraction at some `κ`
    /// — the same quantities, convention and x origin as the exact B-spline kernel's
    /// `InnerIntegral::transforms_at` (no transom closure).
    #[allow(
        dead_code,
        reason = "the open-transom form, exercised by the parity tests"
    )]
    pub(crate) fn transforms_at(&self, zc: &SectionalContracted, kx: f64) -> SquatTransforms {
        self.transforms_at_closed(zc, kx, 0.0, TransomClosure::None)
    }

    /// [`SectionalHull::transforms_at`] with the transom closed by `closure`
    /// at `ν` (the hollow length depends on it).
    pub(crate) fn transforms_at_closed(
        &self,
        zc: &SectionalContracted,
        kx: f64,
        nu: f64,
        closure: TransomClosure,
    ) -> SquatTransforms {
        let p = self.p;
        let mut t = SquatTransforms::default();
        if !zc.any {
            return t;
        }
        let mut xm = Vec::with_capacity(p + 1);
        for (s, sx) in self.spans.iter().enumerate() {
            osc_moments(kx, sx.len, p, &mut xm);
            let d = sx.start - self.x_center;
            let phase = C64::cis(kx * d);
            let (mut q_s, mut q1_s, mut w_s, mut q1w_s) =
                (C64::ZERO, C64::ZERO, C64::ZERO, C64::ZERO);
            for a in 0..p {
                let (g, gw) = (zc.g_fx[s * p + a], self.wl_fx[s * p + a]);
                q_s = q_s + xm[a].scale(g);
                w_s = w_s + xm[a].scale(gw);
                let shifted = xm[a + 1] + xm[a].scale(d);
                q1_s = q1_s + shifted.scale(g);
                q1w_s = q1w_s + shifted.scale(gw);
            }
            let (mut p_s, mut pw_s) = (C64::ZERO, C64::ZERO);
            for a in 0..=p {
                p_s = p_s + xm[a].scale(zc.g_f[s * (p + 1) + a]);
                pw_s = pw_s + xm[a].scale(self.wl_f[s * (p + 1) + a]);
            }
            t.q = t.q + phase * q_s;
            t.q1 = t.q1 + phase * q1_s;
            t.w = t.w + phase * w_s;
            t.q1_wl = t.q1_wl + phase * q1w_s;
            t.p = t.p + phase * p_s;
            t.p_wl = t.p_wl + phase * pw_s;
        }
        t.q_src = t.q;
        self.add_transom_transforms(kx, nu, closure, zc.z_t, &mut t);
        let conj = |v: C64| C64::new(v.re, -v.im);
        t.q_src = conj(t.q_src);
        t.q = conj(t.q);
        t.p = conj(t.p);
        t.q1 = conj(t.q1);
        t.w = conj(t.w);
        t.p_wl = conj(t.p_wl);
        t.q1_wl = conj(t.q1_wl);
        t
    }

    /// Length between the end stations' knots [m].
    pub fn length(&self) -> f64 {
        self.length
    }

    /// Deepest depth node's section depth [m] (the deepest section).
    pub fn draft(&self) -> f64 {
        self.draft
    }

    /// Displaced volume `2∬ f dx dz` [m³].
    pub fn displaced_volume(&self) -> f64 {
        self.volume
    }

    /// Longitudinal centre of buoyancy [m], in the hull's x coordinates.
    pub fn lcb_x(&self) -> f64 {
        self.lcb_x
    }

    /// Waterplane area `∫ 2 f(x, 0) dx` [m²].
    pub fn waterplane_area(&self) -> f64 {
        self.waterplane[0]
    }

    /// First moment of the waterplane about x = 0 [m³].
    pub fn waterplane_moment(&self) -> f64 {
        self.waterplane[1]
    }

    /// Second moment of the waterplane about x = 0 [m⁴].
    pub fn waterplane_second_moment(&self) -> f64 {
        self.waterplane[2]
    }

    /// Wetted surface of both sides [m²] (see [`SectionalHull`]'s strip
    /// construction): the shell itself, not its centreplane projection.
    pub fn wetted_surface(&self) -> f64 {
        self.wetted_surface
    }

    /// Centre of the x domain, the phase reference of the transforms [m].
    pub fn x_center(&self) -> f64 {
        self.x_center
    }

    /// The x domain: the end stations [m].
    pub fn x_range(&self) -> (f64, f64) {
        (
            self.spans.first().map_or(0.0, |s| s.start),
            self.spans.last().map_or(0.0, |s| s.start + s.len),
        )
    }

    /// Waterline half-beam `f(x, 0)` [m] as the kernel interpolates it; 0
    /// outside the hull.
    pub fn waterline_half_beam(&self, x: f64) -> f64 {
        let p = self.p;
        let Some(s) = self
            .spans
            .iter()
            .position(|sp| x >= sp.start && x <= sp.start + sp.len)
        else {
            return 0.0;
        };
        let t = x - self.spans[s].start;
        (0..=p)
            .rev()
            .fold(0.0, |acc, a| acc * t + self.wl_f[s * (p + 1) + a])
            .max(0.0)
    }

    /// The stations: each one's x and its section at the quadrature nodes
    /// (see [`SectionNodes::outline`]).
    pub fn sections(&self) -> impl Iterator<Item = (f64, &[(f64, f64)])> {
        self.xs
            .iter()
            .copied()
            .zip(self.sections.iter().map(|s| s.outline()))
    }

    /// The stations: each one's x and its section curve, sampled densely
    /// (see [`SectionNodes::curve`]).
    pub fn curves(&self) -> impl Iterator<Item = (f64, &[(f64, f64)])> {
        self.xs
            .iter()
            .copied()
            .zip(self.sections.iter().map(|s| s.curve()))
    }

    /// Each station's depth integral `Z(xᵢ; κ)`, and the interpolant the
    /// kernel actually integrates along x, sampled `per_span` times per
    /// x-span: `(stations, curve)` as `(x, Z)` pairs. At `κ = 0` this is the
    /// sectional-area curve (half-areas).
    pub fn depth_integral_curve(
        &self,
        kappa: f64,
        per_span: usize,
    ) -> (Vec<(f64, f64)>, Vec<(f64, f64)>) {
        let at_st = self
            .xs
            .iter()
            .zip(&self.sections)
            .map(|(&x, s)| (x, s.integrate(kappa)))
            .collect();
        let mut zc = SectionalContracted::default();
        self.contract(kappa, &mut zc);
        let p = self.p;
        let mut curve = Vec::new();
        for (s, sx) in self.spans.iter().enumerate() {
            for k in 0..per_span {
                let t = sx.len * k as f64 / per_span as f64;
                let v = if zc.any {
                    (0..=p)
                        .rev()
                        .fold(0.0, |acc, a| acc * t + zc.g_f[s * (p + 1) + a])
                } else {
                    0.0
                };
                curve.push((sx.start + t, v));
            }
        }
        // Close the curve at the last station.
        if let (Some(&xe), Some(se)) = (self.xs.last(), self.sections.last()) {
            curve.push((xe, se.integrate(kappa)));
        }
        (at_st, curve)
    }

    /// The x-derivative of the interpolated depth integral, `∂Z/∂x(x; κ)`,
    /// sampled `per_span` times per x-span as `(x, ∂Z/∂x)` — the source
    /// strength the kernel integrates (the Michell amplitude is its Fourier
    /// transform along x).
    pub fn depth_integral_slope_curve(&self, kappa: f64, per_span: usize) -> Vec<(f64, f64)> {
        let mut zc = SectionalContracted::default();
        self.contract(kappa, &mut zc);
        let p = self.p;
        let mut curve = Vec::new();
        for (s, sx) in self.spans.iter().enumerate() {
            for k in 0..=per_span {
                if k == per_span && s + 1 < self.spans.len() {
                    continue;
                }
                let t = sx.len * k as f64 / per_span as f64;
                let v = if zc.any {
                    (0..p)
                        .rev()
                        .fold(0.0, |acc, a| acc * t + zc.g_fx[s * p + a])
                } else {
                    0.0
                };
                curve.push((sx.start + t, v));
            }
        }
        curve
    }

    /// Number of stations.
    pub fn stations(&self) -> usize {
        self.n
    }

    /// Total depth nodes over all stations.
    pub fn depth_nodes(&self) -> usize {
        self.sections.iter().map(SectionNodes::len).sum()
    }

    /// Number of x-spans (what each `k_x` evaluation walks).
    pub fn x_spans(&self) -> usize {
        self.spans.len()
    }
}

/// The sectional hull as a near-field kernel, for [`crate::squat`].
#[derive(Clone)]
pub(crate) struct SectionalKernel<'h> {
    hull: &'h SectionalHull,
    nu: f64,
    closure: TransomClosure,
}

impl NearFieldKernel for SectionalKernel<'_> {
    type Contracted = SectionalContracted;
    fn contract_z(&mut self, kappa: f64, out: &mut SectionalContracted) {
        self.hull.contract(kappa, out)
    }
    fn transforms_at(&mut self, zc: &SectionalContracted, kx: f64) -> SquatTransforms {
        self.hull
            .transforms_at_closed(zc, kx, self.nu, self.closure)
    }
}

/// The sectional hull as a wave-resistance member, for the shared outer
/// quadrature, with its placement phase offsets.
#[derive(Clone)]
struct SectionalMember<'h> {
    hull: &'h SectionalHull,
    scratch: SectionalContracted,
    dx: f64,
    dy: f64,
    closure: TransomClosure,
}

impl MemberWave for SectionalMember<'_> {
    fn amps(&mut self, nu: f64, lambda: f64) -> (C64, C64) {
        let f = self
            .hull
            .amplitude_closed(nu, lambda, self.closure, &mut self.scratch);
        (f, f)
    }
    fn offsets(&self) -> (f64, f64) {
        (self.dx, self.dy)
    }
}

/// Michell wave resistance of a single sectional hull: an adaptive outer
/// quadrature over `λ`, the transom (if any) closed by `opts.transom`.
pub fn wave_resistance(
    hull: &SectionalHull,
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    multihull_wave_resistance(&[(hull, Placement::default())], cond, opts)
}

/// Wave resistance of a fleet of sectional hulls, placements and interference
/// included: each member's amplitude carries its placement phase
/// `exp(i ν (λ Δx_j ± λ √(λ²−1) y_j))`, and the integrand is the mean of the
/// two (±θ) wave systems. For two identical hulls separated by `s` this
/// reduces to the classical catamaran factor `4 cos²(½ ν s λ √(λ²−1))`.
pub fn multihull_wave_resistance(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    cond.validate()?;
    if members.is_empty() {
        return Err(Error::InvalidConditions("empty fleet".into()));
    }
    let (u, g) = (cond.speed, cond.gravity);
    let nu = g / (u * u);
    let n = members.len() as f64;
    let cx_ref = members.iter().map(|(h, p)| h.x_center + p.x).sum::<f64>() / n;
    let y_ref = members.iter().map(|(_, p)| p.y).sum::<f64>() / n;
    let params = OuterParams {
        nu,
        x_half: members
            .iter()
            .map(|(h, p)| (h.x_center + p.x - cx_ref).abs() + 0.5 * h.length)
            .fold(0.0, f64::max),
        y_half: members
            .iter()
            .map(|(_, p)| (p.y - y_ref).abs())
            .fold(0.0, f64::max),
        t_max: members.iter().map(|(h, _)| h.draft).fold(0.0, f64::max),
    };
    let coeff = 4.0 * cond.fluid.density * g * g / (PI * u * u);
    let mem = members
        .iter()
        .map(|(h, p)| SectionalMember {
            hull: h,
            scratch: SectionalContracted::default(),
            dx: h.x_center + p.x - cx_ref,
            dy: p.y - y_ref,
            closure: opts.transom,
        })
        .collect();
    Ok(run_outer(&params, opts, coeff, mem))
}

/// Wave + viscous resistance of a fleet of sectional hulls, broken down into
/// [`MultihullResistance`]. The viscous part uses each hull's own shell area
/// ([`SectionalHull::wetted_surface`]).
pub fn multihull_resistance(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    wave_opts: &WaveOptions,
    viscous_opts: &ViscousOptions,
) -> Result<MultihullResistance> {
    let wave = multihull_wave_resistance(members, cond, wave_opts)?;
    let mut solo_wave_total = 0.0;
    for m in members {
        solo_wave_total += multihull_wave_resistance(&[*m], cond, wave_opts)?.resistance;
    }
    let viscous: Vec<ViscousResistance> = members
        .iter()
        .map(|(h, _)| viscous_resistance_for(h.length, h.wetted_surface, cond, viscous_opts))
        .collect::<Result<_>>()?;
    let viscous_total: f64 = viscous.iter().map(|v| v.resistance).sum();
    let wetted_surface: f64 = members.iter().map(|(h, _)| h.wetted_surface).sum();
    let total = wave.resistance + viscous_total;
    let q = 0.5 * cond.fluid.density * cond.speed * cond.speed * wetted_surface;
    Ok(MultihullResistance {
        wave,
        viscous,
        viscous_total,
        total,
        effective_power: total * cond.speed,
        wetted_surface,
        cw: wave.resistance / q,
        cv: viscous_total / q,
        ct: total / q,
        solo_wave_total,
        interference: if solo_wave_total > f64::MIN_POSITIVE {
            wave.resistance / solo_wave_total
        } else {
            1.0
        },
    })
}

/// Near-field vertical force and pitch moment on a single sectional hull
/// about `x_ref` (see [`crate::squat`] for the theory), the transom closed by
/// `opts.wave.transom`.
pub fn dynamic_force(
    hull: &SectionalHull,
    cond: &Conditions,
    x_ref: f64,
    opts: &SquatOptions,
) -> Result<DynamicForce> {
    multihull_dynamic_force(&[(hull, Placement::default())], cond, x_ref, opts)
}

/// Near-field force and moment on a fleet of sectional hulls, demihull
/// interaction included: every member's pressure includes what every other
/// member's source sheet induces on it, through the placement phases.
/// `x_ref` is the pitch pivot in fleet coordinates.
pub fn multihull_dynamic_force(
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    x_ref: f64,
    opts: &SquatOptions,
) -> Result<DynamicForce> {
    cond.validate()?;
    if members.is_empty() {
        return Err(Error::InvalidGeometry("empty fleet".into()));
    }
    let nu = cond.gravity / (cond.speed * cond.speed);
    let fleet = Fleet {
        members: members
            .iter()
            .map(|(h, p)| Member {
                inner: SectionalKernel {
                    hull: h,
                    nu,
                    closure: opts.wave.transom,
                },
                cx: h.x_center + p.x,
                y: p.y,
            })
            .collect(),
        nu,
        x_ref,
        scratch: vec![SquatTransforms::default(); members.len()],
        zc_tmp: vec![SectionalContracted::default(); members.len()],
    };
    let l_max = members.iter().map(|(h, _)| h.length).fold(0.0, f64::max);
    let t_max = members.iter().map(|(h, _)| h.draft).fold(0.0, f64::max);
    let volume = members.iter().map(|(h, _)| h.volume).sum();
    Ok(integrate_force(fleet, cond, l_max, t_max, volume, opts))
}

/// The dynamic load a sectional fleet carries at speed, as the equilibrium
/// solver's closure for [`crate::float::solve_equilibrium_sectional_dynamic`].
/// A dry fleet reports zero load: the hydrostatic side of the solver already
/// sinks it until something gets wet.
pub fn dynamic_load_closure<'a>(
    cond: &'a Conditions,
    x_ref: f64,
    opts: &'a SquatOptions,
) -> impl FnMut(&FleetState<SectionalHull>) -> Result<DynamicLoad> + 'a {
    move |fleet: &FleetState<SectionalHull>| {
        if fleet.members.is_empty() {
            return Ok(DynamicLoad::default());
        }
        let members: Vec<(&SectionalHull, Placement)> =
            fleet.members.iter().map(|(h, p)| (h, *p)).collect();
        let d = multihull_dynamic_force(&members, cond, x_ref, opts)?;
        Ok(DynamicLoad {
            force_up: d.force_up,
            moment_bow_up: d.moment_bow_up,
        })
    }
}

/// LU of a banded matrix without pivoting — B-spline collocation matrices
/// are totally positive, for which elimination without pivoting is stable
/// (de Boor). Dense storage, band-limited loops: stations number in the
/// hundreds.
#[derive(Debug, Clone)]
struct BandLu {
    n: usize,
    /// Half-bandwidth: entries with `|i − j| > w` are zero.
    w: usize,
    a: Vec<f64>,
}

impl BandLu {
    fn new(n: usize, p: usize) -> Self {
        BandLu {
            n,
            w: p + 1,
            a: vec![0.0; n * n],
        }
    }

    fn set(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.n + j] = v;
    }

    fn factor(&mut self) -> bool {
        let (n, w) = (self.n, self.w);
        for k in 0..n {
            let piv = self.a[k * n + k];
            if piv.abs() < 1e-14 {
                return false;
            }
            for i in k + 1..(k + w + 1).min(n) {
                let l = self.a[i * n + k] / piv;
                if l == 0.0 {
                    continue;
                }
                self.a[i * n + k] = l;
                for j in k + 1..(k + w + 1).min(n) {
                    self.a[i * n + j] -= l * self.a[k * n + j];
                }
            }
        }
        true
    }

    fn solve(&self, b: &mut [f64]) {
        let (n, w) = (self.n, self.w);
        for i in 0..n {
            let mut s = b[i];
            for k in i.saturating_sub(w)..i {
                s -= self.a[i * n + k] * b[k];
            }
            b[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = b[i];
            for j in i + 1..(i + w + 1).min(n) {
                s -= self.a[i * n + j] * b[j];
            }
            b[i] = s / self.a[i * n + i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::michell::{InnerIntegral, TransomClosure};
    use crate::squat::SquatOptions;

    /// Largest relative disagreement of the six transforms over a (k_x, κ)
    /// grid spanning the near-field quadrature's range, each transform
    /// against its own largest magnitude on the grid.
    fn transform_error(hull: &Hull, sec: &SectionalHull, nu: f64) -> f64 {
        let mut inner = InnerIntegral::new(hull, nu, TransomClosure::None);
        let mut zc = SectionalContracted::default();
        let parts = |t: &SquatTransforms| [t.q, t.p, t.q1, t.w, t.p_wl, t.q1_wl];
        let (mut worst_abs, mut scale) = ([0.0f64; 6], [0.0f64; 6]);
        for ik in 0..40 {
            let kappa = nu * 1e-3 * 10f64.powf(ik as f64 * 0.15); // up to ~1e3·ν
            sec.contract(kappa, &mut zc);
            for ix in 0..40 {
                let kx = nu * (ix as f64 * 0.5 - 10.0);
                let a = parts(&inner.eval_transforms(kx, kappa));
                let b = parts(&sec.transforms_at(&zc, kx));
                for k in 0..6 {
                    worst_abs[k] = worst_abs[k].max((a[k] - b[k]).abs());
                    scale[k] = scale[k].max(a[k].abs());
                }
            }
        }
        (0..6).map(|k| worst_abs[k] / scale[k]).fold(0.0, f64::max)
    }

    fn amplitude_error(hull: &Hull, sec: &SectionalHull, nu: f64) -> f64 {
        let mut inner = InnerIntegral::new(hull, nu, TransomClosure::None);
        let mut zc = SectionalContracted::default();
        let (mut worst, mut scale) = (0.0f64, 0.0f64);
        for i in 0..400 {
            let lambda = 1.0 + i as f64 * 0.25; // κ up to ~1e4·ν
            let a = inner.eval(lambda);
            let b = sec.amplitude(nu, lambda, &mut zc);
            worst = worst.max((a - b).abs());
            scale = scale.max(a.abs());
        }
        worst / scale
    }

    #[test]
    fn reproduces_the_lofted_kernel_on_wigley() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        for nu in [0.3, 1.1, 4.0] {
            let (t, a) = (
                transform_error(&hull, &sec, nu),
                amplitude_error(&hull, &sec, nu),
            );
            assert!(
                t < 1e-9 && a < 1e-9,
                "ν {nu}: transforms {t:.1e}, amplitude {a:.1e}"
            );
        }
    }

    #[test]
    #[ignore = "timing report, not a check"]
    fn cost_against_the_lofted_kernel() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        let nu = 9.81 / 25.0;
        let mut inner = InnerIntegral::new(&hull, nu, TransomClosure::None);
        let n = 2000;
        let mut zl = crate::michell::ZContracted::default();
        let t = std::time::Instant::now();
        for i in 0..n {
            inner.contract_z(nu * (1.0 + i as f64 * 0.01), &mut zl);
        }
        let lofted_k = t.elapsed().as_secs_f64() / n as f64;
        let mut zs = SectionalContracted::default();
        let t = std::time::Instant::now();
        for i in 0..n {
            sec.contract(nu * (1.0 + i as f64 * 0.01), &mut zs);
        }
        let sec_k = t.elapsed().as_secs_f64() / n as f64;
        inner.contract_z(nu * 2.0, &mut zl);
        sec.contract(nu * 2.0, &mut zs);
        let m = 20_000;
        let mut acc = C64::ZERO;
        let t = std::time::Instant::now();
        for i in 0..m {
            acc = acc + inner.transforms_at(&zl, nu * (1.0 + i as f64 * 1e-3)).q;
        }
        let lofted_x = t.elapsed().as_secs_f64() / m as f64;
        let t = std::time::Instant::now();
        for i in 0..m {
            acc = acc + sec.transforms_at(&zs, nu * (1.0 + i as f64 * 1e-3)).q;
        }
        let sec_x = t.elapsed().as_secs_f64() / m as f64;
        eprintln!(
            "lofted {}x{} spans: {:.1} us/κ, {:.1} us/kx\n\
             sectional {} stations, {} depth nodes, {} x-spans: {:.1} us/κ, {:.1} us/kx ({acc:?})",
            hull.spans_x().len(),
            hull.spans_z().len(),
            lofted_k * 1e6,
            lofted_x * 1e6,
            sec.stations(),
            sec.depth_nodes(),
            sec.x_spans(),
            sec_k * 1e6,
            sec_x * 1e6
        );
    }

    /// Speeds spanning Fn ≈ 0.15–0.5 for a hull of length `l`.
    fn speeds(l: f64) -> Vec<f64> {
        [0.15, 0.25, 0.35, 0.5]
            .iter()
            .map(|fn_| fn_ * (9.81 * l).sqrt())
            .collect()
    }

    fn untransomed() -> WaveOptions {
        WaveOptions {
            transom: TransomClosure::None,
            ..WaveOptions::default()
        }
    }

    fn compare_end_to_end(hull: &Hull, name: &str, tol: f64) {
        let sec = SectionalHull::from_hull(hull, &DepthQuadrature::default()).unwrap();
        let wave = untransomed();
        let squat = SquatOptions {
            wave,
            ..SquatOptions::default()
        };
        let x_ref = hull.x_center();
        for u in speeds(hull.length()) {
            let cond = Conditions::seawater(u);
            let t = std::time::Instant::now();
            let a = crate::michell::wave_resistance_with(hull, &cond, &wave).unwrap();
            let fa = crate::squat::dynamic_force(hull, &cond, x_ref, &squat).unwrap();
            let t_lofted = t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            let b = wave_resistance(&sec, &cond, &wave).unwrap();
            let fb = dynamic_force(&sec, &cond, x_ref, &squat).unwrap();
            let t_sec = t.elapsed().as_secs_f64();
            let rel = |x: f64, y: f64| (x - y).abs() / x.abs().max(1e-300);
            let (er, ef, em) = (
                rel(a.resistance, b.resistance),
                rel(fa.force_up, fb.force_up),
                rel(fa.moment_bow_up, fb.moment_bow_up),
            );
            eprintln!(
                "{name} U {u:.2}: Rw {:.4e}/{:.4e} ({er:.1e}), Fz {:.4e} ({ef:.1e}), M {:.4e} ({em:.1e}); \
                 {t_lofted:.3} s lofted, {t_sec:.3} s sectional",
                a.resistance, b.resistance, fa.force_up, fa.moment_bow_up
            );
            assert!(er < tol && ef < tol && em < tol, "{name} U {u}");
        }
        let vol = rel_vol(hull, &sec);
        assert!(vol < 1e-10, "{name} volume {vol:.1e}");
    }

    fn rel_vol(hull: &Hull, sec: &SectionalHull) -> f64 {
        (hull.displaced_volume() - sec.displaced_volume()).abs() / hull.displaced_volume()
    }

    #[test]
    fn wave_resistance_and_squat_match_the_lofted_hull_on_wigley() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        compare_end_to_end(&hull, "wigley", 1e-7);
    }

    /// The sectional importer end to end on geometry with an exact answer:
    /// the Wigley hull written out as CAD patches and cut back into
    /// sections. Its sections are smooth all the way to the keel and its
    /// depth integral is quadratic in x, so the rays and the x interpolant
    /// should both be essentially exact.
    #[test]
    fn iges_sections_of_a_wigley_reproduce_the_exact_hull() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let surfs = crate::iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
        let text = crate::iges::write(&surfs, "wigley").unwrap();
        let imp = crate::iges::import_sectional(&text, &crate::iges::SectionalOptions::default())
            .unwrap()
            .hulls
            .remove(0);
        eprintln!("{:?}", imp.report);
        let sec = &imp.hull;
        let vol = rel_vol(&hull, sec);
        let wave = untransomed();
        let squat = SquatOptions {
            wave,
            ..SquatOptions::default()
        };
        let rel = |x: f64, y: f64| (x - y).abs() / x.abs();
        for u in speeds(hull.length()) {
            let cond = Conditions::seawater(u);
            let a = crate::michell::wave_resistance_with(&hull, &cond, &wave).unwrap();
            let b = wave_resistance(sec, &cond, &wave).unwrap();
            let fa = crate::squat::dynamic_force(&hull, &cond, 0.0, &squat).unwrap();
            let fb = dynamic_force(sec, &cond, 0.0, &squat).unwrap();
            let (er, ef, em) = (
                rel(a.resistance, b.resistance),
                rel(fa.force_up, fb.force_up),
                rel(fa.moment_bow_up, fb.moment_bow_up),
            );
            eprintln!("U {u:.2}: Rw {er:.1e}, Fz {ef:.1e}, M {em:.1e}");
            assert!(
                er < 1e-8 && ef < 1e-9 && em < 1e-7,
                "U {u}: {er:.1e} {ef:.1e} {em:.1e}"
            );
        }
        eprintln!("volume {vol:.1e}");
        assert!(vol < 1e-8);
    }

    /// A real CAD hull with every awkward end at once — a wet transom whose
    /// side skins run on past it, a separately-patched stem, a fine entry —
    /// imported by sections at two resolutions: the transom is carried as a
    /// section (not closed to nothing), and resistance has converged.
    #[test]
    fn cad_sections_converge_and_keep_the_transom() {
        let Some(text) = crate::cad_fixture("e12.igs") else {
            return;
        };
        let wave = untransomed();
        let import = |stations: usize, rays: usize| {
            let so = crate::iges::SectionalOptions {
                waterline_z: -0.95,
                stations,
                rays,
                ..Default::default()
            };
            let mut fleet = crate::iges::import_sectional(&text, &so).unwrap();
            assert!(fleet.failed.is_empty(), "{:?}", fleet.failed);
            assert_eq!(fleet.hulls.len(), 1, "the stem is part of the hull");
            fleet.hulls.remove(0)
        };
        let (coarse, fine) = (import(61, 17), import(121, 33));
        // The aft end station is the transom, 8–9 mm deep and 0.37 m wide.
        let (x_aft, aft) = &fine.sections[0];
        assert!((x_aft - 0.1799).abs() < 1e-3, "aft end at {x_aft}");
        let wl_half = aft.first().unwrap().0;
        let depth = aft.last().unwrap().1;
        assert!(
            wl_half > 0.18 && depth > 0.008,
            "transom {wl_half} x {depth}"
        );
        let rel = (coarse.hull.displaced_volume() - fine.hull.displaced_volume()).abs()
            / fine.hull.displaced_volume();
        assert!(rel < 1e-4, "volume {rel:.1e}");
        let l = fine.hull.length();
        for fnum in [0.2, 0.3] {
            let cond = Conditions::seawater(fnum * (9.81 * l).sqrt());
            let a = wave_resistance(&coarse.hull, &cond, &wave)
                .unwrap()
                .resistance;
            let b = wave_resistance(&fine.hull, &cond, &wave)
                .unwrap()
                .resistance;
            assert!((a - b).abs() / b < 2e-3, "Fn {fnum}: {a} vs {b}");
        }
    }

    /// Hydrostatics from sections against the exact Wigley's: the
    /// waterplane and buoyancy come from the interpolants exactly; the
    /// wetted surface from strips between station curves converges to the
    /// graph's true area.
    #[test]
    fn sectional_hydrostatics_match_the_exact_hull() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let surfs = crate::iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
        let text = crate::iges::write(&surfs, "wigley").unwrap();
        let sec = crate::iges::import_sectional(&text, &crate::iges::SectionalOptions::default())
            .unwrap()
            .hulls
            .remove(0)
            .hull;
        let rel = |a: f64, b: f64| (a - b).abs() / b.abs().max(1e-300);
        let area = rel(sec.waterplane_area(), hull.waterplane_area());
        let m2 = rel(
            sec.waterplane_second_moment(),
            hull.waterplane_second_moment(),
        );
        let m1 = (sec.waterplane_moment() - hull.waterplane_moment()).abs();
        let lcb = (sec.lcb_x() - hull.lcb_x()).abs();
        let s = rel(sec.wetted_surface(), hull.wetted_surface());
        eprintln!("A_w {area:.1e}, I {m2:.1e}, M {m1:.1e}, lcb {lcb:.1e}, S {s:.1e}");
        assert!(
            area < 1e-9 && m2 < 1e-9 && m1 < 1e-9 && lcb < 1e-9,
            "{area} {m2} {m1} {lcb}"
        );
        assert!(s < 1e-3, "wetted surface {s:.1e}");
    }
}

#[cfg(test)]
mod transom_tests {
    use super::*;
    use crate::michell::{InnerIntegral, TransomClosure};

    const E12_WL: f64 = -0.95;

    fn e12_text() -> Option<String> {
        crate::cad_fixture("e12.igs")
    }

    /// A smooth hull with a wet transom: cubic along x, full aft (half-beam
    /// 0.3 m at x = 0) and closed at the bow (x = 8), quadratic in depth and
    /// closed at the keel (T = 0.5).
    fn transom_hull() -> Hull {
        let knots_x = vec![0.0, 0.0, 0.0, 0.0, 4.0, 8.0, 8.0, 8.0, 8.0];
        let knots_z = vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5];
        let control: Vec<f64> = [0.3, 0.45, 0.5, 0.35, 0.0]
            .iter()
            .flat_map(|&c| [c, 0.9 * c, 0.0])
            .collect();
        let surface = crate::bspline::BSplineSurface::new(3, 2, knots_x, knots_z, control).unwrap();
        Hull::new(surface).unwrap()
    }

    fn closures() -> [TransomClosure; 3] {
        [
            TransomClosure::default(),
            TransomClosure::Fixed { length: 0.3 },
            TransomClosure::None,
        ]
    }

    /// The closure machinery against the exact B-spline kernel's: on the
    /// harness (the spline's own stations and transom), the amplitude and all
    /// six near-field transforms with each closure, at two speeds.
    #[test]
    fn transom_closure_matches_the_lofted_kernel() {
        let hull = transom_hull();
        assert!(
            hull.transom().is_some(),
            "the hull should carry a wet transom"
        );
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        let l = hull.length();
        for fnum in [0.2, 0.4] {
            let nu = 1.0 / (fnum * fnum * l);
            for closure in closures() {
                let mut inner = InnerIntegral::new(&hull, nu, closure);
                let mut zc = SectionalContracted::default();
                let (mut worst, mut scale) = (0.0f64, 0.0f64);
                for i in 0..400 {
                    let lambda = 1.0 + i as f64 * 0.25;
                    let a = inner.eval(lambda);
                    let b = sec.amplitude_closed(nu, lambda, closure, &mut zc);
                    worst = worst.max((a - b).abs());
                    scale = scale.max(a.abs());
                }
                let amp = worst / scale;
                let parts = |t: &SquatTransforms| [t.q, t.p, t.q1, t.w, t.p_wl, t.q1_wl];
                let (mut wa, mut sa) = ([0.0f64; 6], [0.0f64; 6]);
                for ik in 0..30 {
                    let kappa = nu * 1e-3 * 10f64.powf(ik as f64 * 0.2);
                    sec.contract(kappa, &mut zc);
                    for ix in 0..30 {
                        let kx = nu * (ix as f64 * 0.7 - 10.0);
                        let a = parts(&inner.eval_transforms(kx, kappa));
                        let b = parts(&sec.transforms_at_closed(&zc, kx, nu, closure));
                        for k in 0..6 {
                            wa[k] = wa[k].max((a[k] - b[k]).abs());
                            sa[k] = sa[k].max(a[k].abs());
                        }
                    }
                }
                let tr = (0..6).map(|k| wa[k] / sa[k]).fold(0.0, f64::max);
                eprintln!("Fn {fnum} {closure:?}: amplitude {amp:.1e}, transforms {tr:.1e}");
                assert!(
                    amp < 1e-9 && tr < 1e-9,
                    "Fn {fnum} {closure:?}: {amp:.1e} {tr:.1e}"
                );
            }
        }
    }

    /// End to end through the outer quadrature and the near-field
    /// integrals: resistance, force and moment with the closure on the
    /// harness equal the lofted hull's.
    #[test]
    fn closed_transom_resistance_and_squat_match_the_lofted_hull() {
        let hull = transom_hull();
        let sec = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        let wave = WaveOptions::default();
        let squat = SquatOptions {
            wave,
            ..SquatOptions::default()
        };
        let x_ref = hull.x_center();
        let rel = |a: f64, b: f64| (a - b).abs() / a.abs();
        for fnum in [0.2, 0.4] {
            let cond = Conditions::seawater(fnum * (9.81 * hull.length()).sqrt());
            let a = crate::michell::wave_resistance_with(&hull, &cond, &wave).unwrap();
            let b = wave_resistance(&sec, &cond, &wave).unwrap();
            let fa = crate::squat::dynamic_force(&hull, &cond, x_ref, &squat).unwrap();
            let fb = dynamic_force(&sec, &cond, x_ref, &squat).unwrap();
            let (er, ef, em) = (
                rel(a.resistance, b.resistance),
                rel(fa.force_up, fb.force_up),
                rel(fa.moment_bow_up, fb.moment_bow_up),
            );
            eprintln!("Fn {fnum}: Rw {er:.1e}, Fz {ef:.1e}, M {em:.1e}");
            assert!(
                er < 1e-8 && ef < 1e-8 && em < 1e-8,
                "Fn {fnum}: {er:.1e} {ef:.1e} {em:.1e}"
            );
        }
    }

    /// Sections cut from CAD find e12's transom: the aft end station, about
    /// 0.37 m wide at the waterline and a centimetre deep.
    #[test]
    fn a_cad_transom_is_detected_from_its_end_section() {
        let so = crate::iges::SectionalOptions {
            waterline_z: E12_WL,
            ..Default::default()
        };
        let Some(text) = e12_text() else {
            return;
        };
        let imp = crate::iges::import_sectional(&text, &so)
            .unwrap()
            .hulls
            .remove(0);
        let tr = imp.hull.transom().expect("e12 has a wet transom");
        eprintln!("transom {tr:?}");
        assert!((tr.x - 0.18).abs() < 5e-3, "x {}", tr.x);
        assert!(
            tr.half_beam > 0.17 && tr.half_beam < 0.2,
            "half-beam {}",
            tr.half_beam
        );
        assert!(tr.depth > 0.004 && tr.depth < 0.012, "depth {}", tr.depth);
        assert_eq!(imp.report.transom.as_ref().map(|t| t.x), Some(tr.x));
    }

    /// The transom's closure is a hollow in the water, not hull: it adds
    /// sources but no pressure-bearing surface. So on e12, whose transom is
    /// barely immersed (6.5 mm equivalent depth, but 0.37 m wide at the
    /// waterline), the near-field force and moment hardly depend on the
    /// hollow's length — they once moved by 35% and changed sign when the
    /// appendage's waterplane was counted as hull.
    #[test]
    fn a_shallow_transom_closure_barely_moves_the_squat_force() {
        let Some(text) = e12_text() else {
            return;
        };
        let opts = crate::iges::SectionalOptions {
            waterline_z: -0.95,
            ..Default::default()
        };
        let fleet = crate::iges::import_sectional(&text, &opts).unwrap();
        let h = &fleet.hulls[0].hull;
        assert!(h.transom().is_some());
        let cond = Conditions::seawater(0.3 * (9.81f64 * h.length()).sqrt());
        let force = |closure| {
            let so = crate::squat::SquatOptions {
                wave: WaveOptions {
                    transom: closure,
                    ..Default::default()
                },
                ..Default::default()
            };
            dynamic_force(h, &cond, h.lcb_x(), &so).unwrap()
        };
        let open = force(TransomClosure::None);
        for closure in [TransomClosure::default(), TransomClosure::Fixed { length: 1.0 }] {
            let d = force(closure);
            let (df, dm) = (
                (d.force_up - open.force_up).abs() / open.force_up.abs(),
                (d.moment_bow_up - open.moment_bow_up).abs() / open.moment_bow_up.abs(),
            );
            assert!(df < 0.05 && dm < 0.05, "{closure:?}: force {df:.3}, moment {dm:.3} off");
        }
    }
}

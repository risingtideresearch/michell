//! # 3D lifting-surface solve (horseshoe vortex lattice)
//!
//! The engine for the rigorous "level 2" asymmetric-hull camber solve. A
//! lifting surface is discretised into panels; each carries a **horseshoe
//! vortex** (a bound segment at the panel quarter-chord plus two trailing legs
//! running downstream to infinity), and flow tangency is collocated at each
//! panel's three-quarter-chord point. The trailing legs make the Kutta
//! condition automatic, exactly as the 1/4–3/4 rule does in the 2D
//! [`crate::lifting`] solver — this module is that solver's 3D generalisation.
//!
//! Solving `A Γ = b` (A = normal-wash influence coefficients, b = the freestream
//! normal component to cancel) gives the bound circulation `Γ` on every panel,
//! hence the loading and the doublet/potential-jump distribution.
//!
//! ## Why this is the hull engine
//!
//! For an asymmetric hull the antisymmetric half-beam `f_a(x, z)` is a camber
//! surface on the centreplane, and its wave resistance needs the centreplane
//! doublet density — the solution of a hypersingular lifting-surface equation
//! (non-local; no pointwise closure). This vortex lattice *is* the numerical
//! solution of that equation: the bound circulation `Γ` integrates to the
//! doublet density `μ`, which feeds the free-wave amplitude
//! `G(λ) = ∬ μ e^{−κz} e^{iνλx}`. This first increment builds and validates the
//! engine on a **flat wing**, where the aspect-ratio dependence of the
//! lift-curve slope brackets it between two exact asymptotic limits:
//!
//! - `AR → ∞` (2-D limit): `dC_L/dα → 2π`;
//! - `AR → 0` (slender-wing / R.T. Jones limit): `dC_L/dα → (π/2)·AR`.
//!
//! The hull-specific wiring — a vertical centreplane, the free-surface rigid-
//! wall image (double-body model), the `∂f_a/∂x` forcing, and passing the
//! solved `μ` to [`crate::spectrum`] — layers on top of this engine and is the
//! next increment.

use crate::lifting::solve_dense;
use std::f64::consts::PI;

type V3 = [f64; 3];

#[inline]
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Squared cut-off below which a field point is treated as lying on a filament
/// (the Biot–Savart kernel is singular there and the contribution is dropped).
const CORE_SQ: f64 = 1e-20;

/// Induced velocity at `p` of a **finite** straight vortex segment from `p1` to
/// `p2` with unit circulation (Katz & Plotkin VORTXL).
fn seg_velocity(p: V3, p1: V3, p2: V3) -> V3 {
    let r1 = sub(p, p1);
    let r2 = sub(p, p2);
    let r1x2 = cross(r1, r2);
    let denom = dot(r1x2, r1x2);
    let n1 = norm(r1);
    let n2 = norm(r2);
    if denom < CORE_SQ || n1 < 1e-10 || n2 < 1e-10 {
        return [0.0, 0.0, 0.0];
    }
    let r0 = sub(p2, p1);
    let k = (dot(r0, r1) / n1 - dot(r0, r2) / n2) / (4.0 * PI * denom);
    scale(r1x2, k)
}

/// Induced velocity at `p` of a **semi-infinite** straight vortex filament that
/// starts at `b` and runs to `+∞` along `+x̂`, carrying unit circulation in the
/// `b → +∞` sense. Closed form (no huge far-point coordinates):
/// `V = (1/4π) (d×r₀)/(|r₀|² − (r₀·d)²) (1 + (r₀·d)/|r₀|)`, `d = x̂`.
fn semi_inf_x_velocity(p: V3, b: V3) -> V3 {
    let r0 = sub(p, b);
    let d: V3 = [1.0, 0.0, 0.0];
    let dxr = cross(d, r0);
    let h2 = dot(dxr, dxr); // |d×r0|² = |r0|² − (r0·d)²  (d unit)
    if h2 < CORE_SQ {
        return [0.0, 0.0, 0.0];
    }
    let n0 = norm(r0);
    let a = dot(r0, d);
    scale(dxr, (1.0 + a / n0) / (4.0 * PI * h2))
}

/// One lifting-surface panel with its horseshoe geometry.
#[derive(Debug, Clone, Copy)]
pub struct Panel {
    /// Bound-vortex endpoints (at the panel quarter-chord).
    pub bound: (V3, V3),
    /// Collocation point (panel three-quarter-chord, mid-span).
    pub collocation: V3,
    /// Unit surface normal at the collocation point.
    pub normal: V3,
    /// Spanwise width used in the Kutta–Joukowski lift sum.
    pub span_width: f64,
}

impl Panel {
    /// Normal-wash influence: the component along this panel's normal of the
    /// velocity a unit-circulation horseshoe on `src` induces here.
    fn influence_from(&self, src: &Panel) -> f64 {
        dot(self.horseshoe_velocity(src), self.normal)
    }

    /// Velocity at this panel's collocation point from a unit horseshoe on
    /// `src`: trailing leg into `b1` (from downstream), bound `b1 → b2`, and
    /// trailing leg out of `b2` (to downstream).
    fn horseshoe_velocity(&self, src: &Panel) -> V3 {
        let (b1, b2) = src.bound;
        let leg_in = scale(semi_inf_x_velocity(self.collocation, b1), -1.0);
        let bound = seg_velocity(self.collocation, b1, b2);
        let leg_out = semi_inf_x_velocity(self.collocation, b2);
        [
            leg_in[0] + bound[0] + leg_out[0],
            leg_in[1] + bound[1] + leg_out[1],
            leg_in[2] + bound[2] + leg_out[2],
        ]
    }
}

/// Result of a 3D lifting-surface solve (freestream speed normalised to 1).
#[derive(Debug, Clone)]
pub struct WingSolution {
    /// Lift coefficient `C_L`.
    pub cl: f64,
    /// Bound circulation on every panel, in input order.
    pub gamma: Vec<f64>,
    /// Reference area used to normalise `C_L`.
    pub area: f64,
}

/// Solve the lifting-surface problem for a set of panels at incidence `alpha`
/// (radians). The freestream is `U∞ = (cos α, 0, sin α)`; the panel normals
/// encode any camber. Returns the bound circulations and `C_L`.
pub fn solve(panels: &[Panel], alpha: f64) -> WingSolution {
    let n = panels.len();
    assert!(n >= 1, "need at least one panel");
    let u_inf: V3 = [alpha.cos(), 0.0, alpha.sin()];

    let mut a = vec![vec![0.0f64; n]; n];
    let mut b = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = panels[i].influence_from(&panels[j]);
        }
        // (U∞ + Σ Γ v)·n = 0  ⇒  Σ Γ (v·n) = −U∞·n.
        b[i] = -dot(u_inf, panels[i].normal);
    }
    let gamma = solve_dense(a, b);

    let area: f64 = panels.iter().map(|p| panel_area(p)).sum();
    // Kutta–Joukowski, freestream ≈ x̂: L = ρU Σ Γ_j Δy_j ⇒ C_L = 2 Σ Γ Δy / S.
    let lift: f64 = gamma
        .iter()
        .zip(panels)
        .map(|(g, p)| g * p.span_width)
        .sum();
    WingSolution {
        cl: 2.0 * lift / area,
        gamma,
        area,
    }
}

/// Panel area from its bound width and the chordwise extent implied by the
/// 1/4–3/4 geometry (bound at 1/4, collocation at 3/4 ⇒ chord = 2×(x_col −
/// x_bound)). Robust for the planar builders here.
fn panel_area(p: &Panel) -> f64 {
    let x_bound = 0.5 * (p.bound.0[0] + p.bound.1[0]);
    let chord = 2.0 * (p.collocation[0] - x_bound);
    chord.abs() * p.span_width
}

/// Build a flat rectangular wing of aspect ratio `ar` (span `b = ar`, chord
/// `c = 1`, area `S = ar`), discretised into `nc` chordwise × `ns` spanwise
/// panels. The wing lies in the `z = 0` plane; incidence is applied through the
/// freestream in [`solve`].
pub fn rectangular_wing(ar: f64, nc: usize, ns: usize) -> Vec<Panel> {
    let chord = 1.0;
    let span = ar * chord;
    let dy = span / ns as f64;
    let dc = chord / nc as f64;
    let mut panels = Vec::with_capacity(nc * ns);
    for j in 0..ns {
        let y_left = -0.5 * span + j as f64 * dy;
        let y_right = y_left + dy;
        let y_mid = 0.5 * (y_left + y_right);
        for i in 0..nc {
            let x_le = i as f64 * dc;
            let x_bound = x_le + 0.25 * dc;
            let x_col = x_le + 0.75 * dc;
            panels.push(Panel {
                bound: ([x_bound, y_left, 0.0], [x_bound, y_right, 0.0]),
                collocation: [x_col, y_mid, 0.0],
                normal: [0.0, 0.0, 1.0],
                span_width: dy,
            });
        }
    }
    panels
}

/// Spanwise loading: total bound circulation summed over the chordwise strip at
/// each spanwise station, returned as `(y_mid, Γ_strip)` for an `nc × ns` wing
/// built by [`rectangular_wing`] (panels in that row-major order).
pub fn spanwise_loading(panels: &[Panel], sol: &WingSolution, nc: usize, ns: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(ns);
    for j in 0..ns {
        let mut g = 0.0;
        let mut y = 0.0;
        for i in 0..nc {
            let idx = j * nc + i;
            g += sol.gamma[idx];
            y = panels[idx].collocation[1];
        }
        out.push((y, g));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const DEG: f64 = PI / 180.0;

    fn lift_slope(ar: f64, nc: usize, ns: usize) -> f64 {
        // C_L is linear in α; take the slope from a small angle.
        let a = 3.0 * DEG;
        solve(&rectangular_wing(ar, nc, ns), a).cl / a
    }

    /// Zero incidence ⇒ zero lift, and C_L is exactly proportional to sin α
    /// (the freestream is `(cos α, 0, sin α)`, so the whole linear system scales
    /// with sin α).
    #[test]
    fn zero_and_proportional_to_sin_alpha() {
        let w = rectangular_wing(6.0, 4, 20);
        assert!(solve(&w, 0.0).cl.abs() < 1e-12);
        let c1 = solve(&w, 2.0 * DEG).cl;
        let c2 = solve(&w, 4.0 * DEG).cl;
        assert!(c1 > 0.0, "positive incidence must give positive lift");
        let k1 = c1 / (2.0 * DEG).sin();
        let k2 = c2 / (4.0 * DEG).sin();
        assert!((k1 - k2).abs() < 1e-9 * k1, "C_L not ∝ sin α: {k1} vs {k2}");
    }

    /// The lift-curve slope stays below the 2-D value and climbs monotonically
    /// toward it as the aspect ratio grows (finite-span downwash).
    #[test]
    fn slope_approaches_two_pi_with_aspect_ratio() {
        let mut prev = 0.0;
        for &ar in &[4.0, 8.0, 16.0, 32.0] {
            let a = lift_slope(ar, 4, 40);
            assert!(a < 2.0 * PI, "AR={ar}: slope {a} must stay below 2π");
            assert!(a > prev, "AR={ar}: slope {a} did not increase from {prev}");
            prev = a;
        }
        // By AR = 32 a rectangular wing is ~92% of 2π (lifting-line gives
        // 2π/(1+2/32) = 5.91; a rectangular planform sits a little under it).
        assert!(prev > 0.90 * 2.0 * PI, "AR=32 slope {prev} too far below 2π");
    }

    /// Prandtl lifting-line sanity: for a high-AR rectangular wing the slope is
    /// close to a = 2π/(1 + 2/AR) (elliptic-loading estimate; a rectangular
    /// wing sits a little under it).
    #[test]
    fn slope_matches_lifting_line_estimate() {
        let ar = 12.0;
        let a = lift_slope(ar, 5, 48);
        let llt = 2.0 * PI / (1.0 + 2.0 / ar);
        assert!(
            a < llt && a > 0.9 * llt,
            "AR={ar}: VLM slope {a}, LLT estimate {llt}"
        );
    }

    /// Slender-wing (R.T. Jones) limit: for small AR, dC_L/dα → (π/2)·AR, so
    /// the slope-to-AR ratio approaches π/2 from... below, staying well under
    /// the 2-D value.
    #[test]
    fn slender_wing_limit() {
        let ar = 0.5;
        let a = lift_slope(ar, 6, 24);
        let jones = 0.5 * PI * ar;
        // Within ~25% of the leading-order slender estimate, and far below 2π.
        assert!(
            (a - jones).abs() < 0.25 * jones,
            "AR={ar}: slope {a}, Jones {jones}"
        );
        assert!(a < 0.5 * 2.0 * PI);
    }

    /// The spanwise loading is symmetric about the root and falls to (near)
    /// zero at the tips — the qualitative shape any correct 3-D solve produces.
    #[test]
    fn loading_is_symmetric_and_tip_relieved() {
        let (nc, ns) = (4, 40);
        let w = rectangular_wing(8.0, nc, ns);
        let sol = solve(&w, 5.0 * DEG);
        let load = spanwise_loading(&w, &sol, nc, ns);
        // Symmetry: station j mirrors station ns-1-j.
        for j in 0..ns / 2 {
            let (yl, gl) = load[j];
            let (yr, gr) = load[ns - 1 - j];
            assert!((yl + yr).abs() < 1e-9, "stations not mirrored in y");
            assert!((gl - gr).abs() < 1e-9 * gl.abs().max(1e-12), "loading not symmetric");
        }
        // Tip relief: the outermost strip carries much less than the root.
        let root = load[ns / 2].1.abs();
        let tip = load[0].1.abs();
        assert!(tip < 0.5 * root, "tip loading {tip} not relieved vs root {root}");
        assert!(tip > 0.0);
    }

    /// Grid convergence: the slope settles as the lattice is refined.
    #[test]
    fn grid_convergence() {
        let ar = 8.0;
        let coarse = lift_slope(ar, 3, 24);
        let medium = lift_slope(ar, 5, 40);
        let fine = lift_slope(ar, 8, 64);
        assert!((fine - medium).abs() < (medium - coarse).abs(), "not settling");
        assert!((fine - medium).abs() < 0.03 * fine, "fine grid not converged");
    }
}

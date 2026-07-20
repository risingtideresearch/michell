//! # 2D thin-airfoil lifting solve — a proving ground for level-2 asymmetry
//!
//! Rigorous wave resistance of a **port/starboard-asymmetric** hull needs the
//! antisymmetric (camber) part treated as a *lifting* problem: the centreplane
//! doublet density is fixed by the *mean* of the two-sided normal velocities,
//! which is a **non-local** functional of the density — unlike the source
//! (thickness) strength, which the boundary condition fixes pointwise. The
//! crate's asymmetric wave-resistance path sidesteps this with a crude
//! prescribed strip closure `μ = 2U f_a`; the honest fix ("level 2") is to
//! *solve* a hypersingular lifting-surface integral equation for the doublet
//! density.
//!
//! Before committing to the delicate 3D free-surface solve, this module builds
//! and validates its **2D ancestor** — classical thin-airfoil theory — which
//! exercises the same machinery that makes the 3D version risky:
//!
//! - discretising a lifting sheet into elements and assembling an influence
//!   matrix,
//! - imposing the **Kutta condition** (smooth flow-off at the trailing edge)
//!   so the circulation is uniquely determined,
//! - recovering the loading (bound vorticity / doublet density) and integrating
//!   it to force and moment coefficients.
//!
//! Two independent solvers are provided and cross-checked against each other
//! and against closed-form thin-airfoil results:
//!
//! 1. [`solve_vortex_lattice`] — the lumped-vortex (vortex-lattice) method with
//!    the Pistolesi 1/4–3/4 rule, the direct 2D form of the 3D vortex/doublet
//!    lattice. The 3/4-chord collocation enforces the Kutta condition
//!    automatically, and the point-vortex kernel is the desingularised form of
//!    the hypersingular kernel the 3D solve must integrate.
//! 2. [`glauert`] — the closed-form Fourier-series (Glauert) solution, with the
//!    series coefficients obtained by numerical projection of the camber slope.
//!    Independent of the lattice; used as the analytic reference.
//!
//! ## Theory
//!
//! A thin airfoil of chord `c` (here normalised to `c = 1`) at incidence `α`
//! with camber-line slope `dη/dx(x)` is represented by a bound-vortex sheet
//! `γ(x)` on the chord line. Flow tangency (the linearised body condition) is
//!
//! ```text
//! (1/2π) ⨍₀¹ γ(ξ)/(x − ξ) dξ = U (α − dη/dx(x)),
//! ```
//!
//! a Cauchy singular integral equation of the first kind. Its Glauert solution
//! with `x = ½(1 − cos θ)` is
//!
//! ```text
//! γ(θ)/U = 2[ A₀ (1+cosθ)/sinθ + Σ_{k≥1} Aₖ sin(kθ) ],
//! A₀ = α − (1/π) ∫₀^π (dη/dx) dθ,   Aₖ = (2/π) ∫₀^π (dη/dx) cos(kθ) dθ,
//! C_L = π(2A₀ + A₁),   C_{m,c/4} = −(π/4)(A₁ − A₂).
//! ```
//!
//! The `(1+cosθ)/sinθ` term is the physical `√`-type leading-edge singularity;
//! the `sin(kθ)` terms vanish at both edges, so the Kutta condition (zero
//! loading at the trailing edge, `θ = π`) is built in. A flat plate
//! (`dη/dx ≡ 0`) gives `A₀ = α`, hence `C_L = 2πα` — the canonical check.
//!
//! ## Mapping back to the hull problem
//!
//! In hull terms the camber line is the antisymmetric half-beam `f_a(x)` at a
//! waterline strip and the forcing `U(α − dη/dx)` is `−U ∂f_a/∂x` (no separate
//! incidence). The solved loading `γ` integrates to the doublet density
//! `μ(x) = ∫ γ` (the potential jump), which is exactly the quantity the 3D
//! free-wave amplitude `G(λ) = ∬ μ e^{−κz} e^{iνλx}` carries. This module
//! returns `μ` alongside the coefficients so that bridge is explicit.

use crate::quadrature::gauss_legendre;
use std::f64::consts::PI;

/// Result of a 2D thin-airfoil lifting solve (chord normalised to 1, loading
/// scaled by the free-stream speed `U` so the results are dimensionless).
#[derive(Debug, Clone)]
pub struct AirfoilSolution {
    /// Lift coefficient `C_L`.
    pub cl: f64,
    /// Pitching-moment coefficient about the leading edge (nose-up positive).
    pub cm_le: f64,
    /// Pitching-moment coefficient about the quarter-chord.
    pub cm_quarter: f64,
    /// Centre of pressure as a chord fraction `x_cp/c = −C_{m,LE}/C_L`
    /// (`NaN` when `C_L = 0`).
    pub center_of_pressure: f64,
    /// Bound-vortex loading `γ(x)/U` sampled at `(x/c, γ/U)`.
    pub loading: Vec<(f64, f64)>,
    /// Doublet density / potential jump `μ(x)/U = ∫₀ˣ γ/U dξ` sampled at
    /// `(x/c, μ/U)` — the 3D free-wave-amplitude integrand's 2D analogue.
    pub doublet: Vec<(f64, f64)>,
}

/// Solve the thin-airfoil problem by the lumped-vortex (vortex-lattice) method.
///
/// The chord `[0, 1]` is split into `n` equal panels; a point vortex sits at
/// each panel's quarter-point and flow tangency is collocated at each panel's
/// three-quarter-point (the Pistolesi 1/4–3/4 rule, which enforces the Kutta
/// condition automatically). `dydx(x)` is the camber-line slope at chord
/// fraction `x ∈ [0, 1]`; `alpha` is the angle of attack in radians.
///
/// A flat plate (`dydx` returning 0) yields `C_L = 2πα` exactly for any `n ≥ 1`.
pub fn solve_vortex_lattice(alpha: f64, dydx: impl Fn(f64) -> f64, n: usize) -> AirfoilSolution {
    assert!(n >= 1, "need at least one panel");
    let h = 1.0 / n as f64;
    // Vortex points (quarter) and collocation points (three-quarter).
    let vortex: Vec<f64> = (0..n).map(|j| (j as f64 + 0.25) * h).collect();
    let collo: Vec<f64> = (0..n).map(|i| (i as f64 + 0.75) * h).collect();

    // Influence matrix A_ij = 1/(2π (x_i − ξ_j)); RHS w_i = α − dη/dx(x_i).
    let mut a = vec![vec![0.0f64; n]; n];
    let mut b = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = 1.0 / (2.0 * PI * (collo[i] - vortex[j]));
        }
        b[i] = alpha - dydx(collo[i]);
    }
    let g = solve_dense(a, b); // g_j = Γ_j / U

    // Coefficients (c = 1): C_L = 2 Σ g_j; moments from the vortex positions.
    let cl: f64 = 2.0 * g.iter().sum::<f64>();
    let cm_le: f64 = -2.0 * g.iter().zip(&vortex).map(|(gj, xj)| gj * xj).sum::<f64>();
    let cm_quarter: f64 =
        -2.0 * g.iter().zip(&vortex).map(|(gj, xj)| gj * (xj - 0.25)).sum::<f64>();

    // Loading γ/U is a set of point circulations Γ_j/U; report the equivalent
    // sheet strength γ ≈ Γ_j / h at each vortex station.
    let loading: Vec<(f64, f64)> = vortex.iter().zip(&g).map(|(x, gj)| (*x, gj / h)).collect();
    // Doublet μ(x)/U = cumulative bound circulation ahead of x.
    let mut cum = 0.0;
    let doublet: Vec<(f64, f64)> = vortex
        .iter()
        .zip(&g)
        .map(|(x, gj)| {
            cum += gj;
            (*x, cum)
        })
        .collect();

    AirfoilSolution {
        cl,
        cm_le,
        cm_quarter,
        center_of_pressure: if cl != 0.0 { -cm_le / cl } else { f64::NAN },
        loading,
        doublet,
    }
}

/// Closed-form thin-airfoil coefficients via the Glauert Fourier series, with
/// the coefficients `A₀ … A_{n_terms}` obtained by `n_quad`-point Gauss
/// projection of the camber slope. Independent of [`solve_vortex_lattice`];
/// used as the analytic reference.
pub fn glauert(
    alpha: f64,
    dydx: impl Fn(f64) -> f64,
    n_terms: usize,
    n_quad: usize,
) -> AirfoilSolution {
    assert!(n_terms >= 2, "need A0..A2 for the moment");
    let (nodes, weights) = gauss_legendre(n_quad);
    // θ ∈ [0, π]; x = ½(1 − cos θ).
    let slope_at = |theta: f64| dydx(0.5 * (1.0 - theta.cos()));

    // A0 = α − (1/π)∫₀^π s dθ; Ak = (2/π)∫₀^π s cos(kθ) dθ.
    let mut acoef = vec![0.0f64; n_terms + 1];
    for (k, coef) in acoef.iter_mut().enumerate() {
        let mut integral = 0.0;
        for (t, w) in nodes.iter().zip(&weights) {
            let theta = 0.5 * PI * (t + 1.0); // map [-1,1] → [0,π]
            let jac = 0.5 * PI;
            let s = slope_at(theta);
            integral += w * jac * s * (k as f64 * theta).cos();
        }
        *coef = if k == 0 {
            alpha - integral / PI
        } else {
            2.0 * integral / PI
        };
    }

    let cl = PI * (2.0 * acoef[0] + acoef[1]);
    let cm_quarter = -(PI / 4.0) * (acoef[1] - acoef[2]);
    let cm_le = cm_quarter - 0.25 * cl;

    // Sample γ(θ)/U and μ(x)/U at safe interior stations (avoid the LE
    // singularity at θ = 0).
    let m = 200usize;
    let mut loading = Vec::with_capacity(m);
    for i in 0..m {
        let theta = PI * (i as f64 + 0.5) / m as f64;
        let x = 0.5 * (1.0 - theta.cos());
        let mut g = acoef[0] * (1.0 + theta.cos()) / theta.sin();
        for (k, coef) in acoef.iter().enumerate().skip(1) {
            g += coef * (k as f64 * theta).sin();
        }
        loading.push((x, 2.0 * g));
    }
    // μ(x)/U by trapezoidal integration of the sampled loading.
    let mut doublet = Vec::with_capacity(m);
    let mut cum = 0.0;
    for i in 0..m {
        if i > 0 {
            let (x0, g0) = loading[i - 1];
            let (x1, g1) = loading[i];
            cum += 0.5 * (g0 + g1) * (x1 - x0);
        }
        doublet.push((loading[i].0, cum));
    }

    AirfoilSolution {
        cl,
        cm_le,
        cm_quarter,
        center_of_pressure: if cl != 0.0 { -cm_le / cl } else { f64::NAN },
        loading,
        doublet,
    }
}

/// Camber-line slope of a parabolic (circular-arc) mean line with maximum
/// camber ratio `m` at mid-chord: `η/c = 4m (x/c)(1 − x/c)`, so
/// `dη/dx = 4m(1 − 2x/c)`. Analytic result: `C_L = 2π(α + 2m)`,
/// `C_{m,c/4} = −π m`. A convenience for tests and demos.
pub fn parabolic_camber_slope(m: f64) -> impl Fn(f64) -> f64 {
    move |x: f64| 4.0 * m * (1.0 - 2.0 * x)
}

/// Dense linear solve `A x = b` by Gaussian elimination with partial pivoting.
/// `a` is consumed row-major; systems here are small (a few hundred rows).
#[allow(clippy::needless_range_loop)] // index loops read clearest for elimination
fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot.
        let mut piv = col;
        let mut best = a[col][col].abs();
        for r in (col + 1)..n {
            let v = a[r][col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let d = a[col][col];
        debug_assert!(d != 0.0, "singular lifting influence matrix");
        for r in (col + 1)..n {
            let factor = a[r][col] / d;
            if factor != 0.0 {
                for c in col..n {
                    a[r][c] -= factor * a[col][c];
                }
                b[r] -= factor * b[col];
            }
        }
    }
    // Back-substitution.
    let mut x = vec![0.0f64; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for c in (row + 1)..n {
            s -= a[row][c] * x[c];
        }
        x[row] = s / a[row][row];
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    const DEG: f64 = PI / 180.0;

    /// Flat plate: C_L = 2πα exactly, for any panel count (the 1/4–3/4 rule is
    /// exact here), with the centre of pressure at the quarter-chord.
    #[test]
    fn flat_plate_is_exact() {
        let alpha = 5.0 * DEG;
        for n in [1usize, 2, 8, 64] {
            let s = solve_vortex_lattice(alpha, |_| 0.0, n);
            assert!(
                (s.cl - 2.0 * PI * alpha).abs() < 1e-12,
                "n={n}: C_L={} vs {}",
                s.cl,
                2.0 * PI * alpha
            );
            assert!(s.cm_quarter.abs() < 1e-12, "n={n}: C_m,c/4={}", s.cm_quarter);
            assert!(
                (s.center_of_pressure - 0.25).abs() < 1e-12,
                "n={n}: x_cp={}",
                s.center_of_pressure
            );
        }
    }

    /// Glauert reference reproduces the analytic flat-plate result.
    #[test]
    fn glauert_flat_plate() {
        let alpha = 4.0 * DEG;
        let s = glauert(alpha, |_| 0.0, 8, 64);
        assert!((s.cl - 2.0 * PI * alpha).abs() < 1e-10);
        assert!(s.cm_quarter.abs() < 1e-10);
    }

    /// Parabolic camber: both solvers hit the closed-form `C_L = 2π(α + 2m)`
    /// and `C_{m,c/4} = −πm`, independent of α for the moment.
    #[test]
    fn parabolic_camber_matches_closed_form() {
        let m = 0.05;
        let slope = parabolic_camber_slope(m);
        for &alpha in &[0.0, 3.0 * DEG, 6.0 * DEG] {
            let cl_exact = 2.0 * PI * (alpha + 2.0 * m);
            let cm_exact = -PI * m;

            let gl = glauert(alpha, &slope, 8, 128);
            assert!(
                (gl.cl - cl_exact).abs() < 1e-9,
                "glauert C_L={} vs {cl_exact}",
                gl.cl
            );
            assert!(
                (gl.cm_quarter - cm_exact).abs() < 1e-9,
                "glauert C_m={} vs {cm_exact}",
                gl.cm_quarter
            );

            // The 1/4–3/4 rule captures the total lift of a linear-slope
            // (parabolic) mean line to machine precision; the moment of the
            // distributed load converges as O(1/n²), so allow it more room.
            let vl = solve_vortex_lattice(alpha, &slope, 64);
            assert!(
                (vl.cl - cl_exact).abs() < 1e-9,
                "VLM C_L={} vs {cl_exact}",
                vl.cl
            );
            assert!(
                (vl.cm_quarter - cm_exact).abs() < 5e-4,
                "VLM C_m={} vs {cm_exact}",
                vl.cm_quarter
            );
        }
    }

    /// The zero-lift angle of the parabolic mean line is α₀ = −2m.
    #[test]
    fn zero_lift_angle() {
        let m = 0.04;
        let s = glauert(-2.0 * m, parabolic_camber_slope(m), 8, 128);
        assert!(s.cl.abs() < 1e-9, "C_L at α₀ = {}", s.cl);
    }

    /// For a genuinely curved (quadratic-slope) mean line the lattice is *not*
    /// exact, so its error against the analytic reference must shrink as the
    /// lattice is refined. A quadratic slope is spanned by `A₀, A₁, A₂`, so the
    /// Glauert result is itself exact and a clean reference.
    #[test]
    fn vortex_lattice_converges_to_glauert() {
        let slope = |x: f64| 0.3 * x * x;
        let alpha = 4.0 * DEG;
        let gl = glauert(alpha, &slope, 4, 128);
        let mut prev = f64::INFINITY;
        for &n in &[10usize, 40, 160] {
            let vl = solve_vortex_lattice(alpha, &slope, n);
            let err = (vl.cl - gl.cl).abs();
            assert!(
                err < prev,
                "n={n}: error {err:e} did not decrease from {prev:e}"
            );
            prev = err;
        }
        assert!(prev < 1e-4, "finest lattice error {prev:e}");
    }

    /// Loading vanishes at the trailing edge (Kutta) and the doublet density
    /// grows monotonically from the leading edge — the qualitative shape the
    /// 3D solve must also produce.
    #[test]
    fn loading_satisfies_kutta_and_doublet_is_monotone() {
        let s = solve_vortex_lattice(5.0 * DEG, parabolic_camber_slope(0.05), 80);
        // Trailing-edge loading (last panel) is small relative to the peak.
        let peak = s.loading.iter().map(|(_, g)| g.abs()).fold(0.0, f64::max);
        let te = s.loading.last().unwrap().1.abs();
        assert!(te < 0.15 * peak, "TE loading {te} not small vs peak {peak}");
        // Positive-lift doublet accumulates monotonically.
        for w in s.doublet.windows(2) {
            assert!(w[1].1 >= w[0].1 - 1e-12, "doublet not monotone");
        }
        assert!(s.doublet.last().unwrap().1 > 0.0);
    }
}

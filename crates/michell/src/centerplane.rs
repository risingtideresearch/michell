//! # Hull centreplane lifting solve (double-body model)
//!
//! Applies the [`crate::lifting3d`] horseshoe vortex lattice to the **vertical
//! centreplane** of an asymmetric hull — the step that replaces the crude
//! prescribed strip closure `μ = 2U f_a` with an actual solution of the
//! lifting-surface equation.
//!
//! ## Geometry
//!
//! The centreplane occupies `x ∈ [0, L]` (streamwise, along the hull) and
//! `z ∈ [0, T]` (vertical, downward from the waterline to the draft `T`). The
//! antisymmetric half-beam `f_a(x, z)` is a **camber surface** on this plane;
//! its streamwise slope `∂f_a/∂x` is the downwash the transverse (y-normal)
//! doublet sheet must produce — the direct analogue of a wing's camber slope.
//!
//! ## Free surface as a rigid wall (double body)
//!
//! In the low-Froude near field the free surface `z = 0` behaves as a rigid
//! wall (`∂φ/∂z = 0`), consistent with how Michell's source strength is set on
//! the double-body flow. We enforce it by **reflecting the centreplane about
//! `z = 0`** with the same circulation sign, so the bound vortices run
//! continuously through the waterline: the waterline acts as a wing **root**
//! (maximum loading) and the keel `z = T` as a free **tip** (loading relieved
//! to zero). A surface-piercing lifting surface of draft `T` therefore behaves
//! as a wing of span `2T` — the classic effective-aspect-ratio doubling, which
//! this module's tests check against [`crate::lifting3d::rectangular_wing`].
//!
//! ## Output and the next step
//!
//! Solving for the bound circulation gives the **doublet density**
//! `μ(x, z) = ∫ γ` (the potential jump) on the physical half `z ≥ 0`. That is
//! exactly the field the free-wave amplitude `G(λ) = ∬ μ e^{−κz} e^{iνλx}`
//! integrates ([`doublet_free_wave_amplitude`] computes it here). Feeding this
//! `G` into [`crate::spectrum`] in place of the strip closure — once the
//! asymmetric-`Hull` decomposition and this solver share a branch — is the
//! final wiring step.

use crate::lifting3d::{solve, Panel};

/// A sample of the doublet density `μ(x, z)/U` (potential jump) at a physical
/// centreplane collocation point.
#[derive(Debug, Clone, Copy)]
pub struct Doublet {
    /// Streamwise position [m].
    pub x: f64,
    /// Depth below the waterline [m] (z downward, `0 ≤ z ≤ T`).
    pub z: f64,
    /// Doublet density / potential jump `μ/U` [m].
    pub mu: f64,
}

/// Result of a centreplane lifting solve.
#[derive(Debug, Clone)]
pub struct CenterplaneSolution {
    /// Side-force coefficient `C_Y` (reference area `L·T`, the physical
    /// centreplane).
    pub side_force: f64,
    /// Doublet density `μ(x, z)/U` on the physical half `z ≥ 0`, one sample per
    /// physical panel, ordered spanwise-major (increasing `z`) then streamwise.
    pub doublet: Vec<Doublet>,
    /// Streamwise / vertical panel counts on the physical half.
    pub nx: usize,
    pub nz: usize,
    length: f64,
    draft: f64,
}

/// Solve the centreplane lifting problem for an asymmetric hull whose
/// antisymmetric half-beam has streamwise slope `dfadx(x, z)` on
/// `x ∈ [0, length]`, `z ∈ [0, draft]` (z downward). The free surface is a
/// rigid wall (double body): the centreplane is reflected about `z = 0` and the
/// physical solution is the `z ≥ 0` half.
pub fn solve_centerplane(
    length: f64,
    draft: f64,
    dfadx: impl Fn(f64, f64) -> f64,
    nx: usize,
    nz: usize,
) -> CenterplaneSolution {
    assert!(nx >= 1 && nz >= 1, "need at least one panel each way");
    let dx = length / nx as f64;
    let dz = draft / nz as f64;

    // Panels over z ∈ [−T, T]: 2·nz spanwise strips, image half (z < 0) first.
    // Bound vortices run along z at each panel's quarter-chord; trailing legs
    // run downstream (+x). The camber enters through the panel normal, so the
    // freestream stays (1, 0, 0) — call `solve` with zero incidence.
    let mut panels = Vec::with_capacity(nx * 2 * nz);
    for kz in 0..(2 * nz) {
        let z_a = -draft + kz as f64 * dz;
        let z_b = z_a + dz;
        let z_mid = 0.5 * (z_a + z_b);
        let z_phys = z_mid.abs(); // forcing mirrored about the waterline
        for ix in 0..nx {
            let x_le = ix as f64 * dx;
            let x_bound = x_le + 0.25 * dx;
            let x_col = x_le + 0.75 * dx;
            let slope = dfadx(x_col, z_phys);
            let nlen = (1.0 + slope * slope).sqrt();
            panels.push(Panel {
                bound: ([x_bound, 0.0, z_a], [x_bound, 0.0, z_b]),
                collocation: [x_col, 0.0, z_mid],
                normal: [-slope / nlen, 1.0 / nlen, 0.0],
                span_width: dz,
            });
        }
    }
    let sol = solve(&panels, 0.0);

    // Doublet density μ(x,z) = cumulative bound circulation from the leading
    // edge, per physical spanwise strip (z > 0).
    let mut doublet = Vec::with_capacity(nx * nz);
    let mut sum_phys = 0.0;
    for kz in nz..(2 * nz) {
        let mut cum = 0.0;
        for ix in 0..nx {
            let idx = kz * nx + ix;
            cum += sol.gamma[idx];
            let c = panels[idx].collocation;
            doublet.push(Doublet {
                x: c[0],
                z: c[2],
                mu: cum,
            });
            sum_phys += sol.gamma[idx] * dz;
        }
    }

    // Side-force coefficient on the physical reference area L·T. Only the
    // physical (z > 0) half is a real surface; C_Y = 2 Σ_phys Γ Δz / (L·T).
    let side_force = 2.0 * sum_phys / (length * draft);

    CenterplaneSolution {
        side_force,
        doublet,
        nx,
        nz,
        length,
        draft,
    }
}

impl CenterplaneSolution {
    /// The free-wave amplitude `G(λ) = ∬ μ(x,z) e^{−νλ²z} e^{iνλx} dx dz` over
    /// the physical centreplane, returned as `(Re, Im)`. This is the dipole
    /// integrand of the asymmetric wave resistance; the dipole free-wave
    /// amplitude is `A_d = i·νλ√(λ²−1)·G` and the resistance follows as in
    /// [`crate::spectrum`]. Evaluated by midpoint quadrature over the panel
    /// grid (each `μ` sample carries cell area `Δx·Δz`).
    pub fn doublet_free_wave_amplitude(&self, nu: f64, lambda: f64) -> (f64, f64) {
        let dx = self.length / self.nx as f64;
        let dz = self.draft / self.nz as f64;
        let cell = dx * dz;
        let kx = nu * lambda;
        let kappa = nu * lambda * lambda;
        let (mut re, mut im) = (0.0, 0.0);
        for d in &self.doublet {
            let w = d.mu * (-kappa * d.z).exp() * cell;
            re += w * (kx * d.x).cos();
            im += w * (kx * d.x).sin();
        }
        (re, im)
    }
}

/// Camber-slope field of a hull that is straight but at a small **drift angle**
/// `beta` (yaw): `f_a = β·x`, so `∂f_a/∂x = β` everywhere — the surface-piercing
/// flat-plate case, used to validate the solver against a wing of the
/// equivalent doubled aspect ratio. Returned as a closure of `(x, z)`.
pub fn drift_slope(beta: f64) -> impl Fn(f64, f64) -> f64 {
    move |_x: f64, _z: f64| beta
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lifting3d::{rectangular_wing, solve};

    const DEG: f64 = std::f64::consts::PI / 180.0;

    /// A surface-piercing flat plate of draft `T` and chord `L` must reproduce
    /// the lift of a full wing of span `2T` and chord `L`: the free-surface
    /// rigid-wall image **doubles the effective aspect ratio**. A geometric
    /// aspect ratio `T/L = 3` plate behaves as an `AR = 2T/L = 6` wing.
    ///
    /// The two code paths agree to *linear order*: `rectangular_wing` tilts the
    /// freestream (exact normal), while the centreplane tilts the panel normal
    /// to encode the camber slope (the correct formulation for real camber), so
    /// they differ only at O(α²). At α = 3° that is a ~0.1% residual.
    #[test]
    fn surface_piercing_plate_doubles_aspect_ratio() {
        let alpha = 3.0 * DEG;
        let ar = 6.0; // wing span 2T = ar ⇒ T = ar/2, chord L = 1
        let (nc, ns) = (5, 30);
        let wing = solve(&rectangular_wing(ar, nc, ns), alpha);
        // Draft ar/2, ns = 2·nz spanwise strips over the doubled span.
        let cp = solve_centerplane(1.0, ar / 2.0, drift_slope(alpha.tan()), nc, ns / 2);
        let rel = (cp.side_force.abs() - wing.cl.abs()).abs() / wing.cl.abs();
        assert!(
            rel < 3e-3,
            "effective-AR doubling: C_Y {} vs wing C_L {} (rel {rel:e})",
            cp.side_force,
            wing.cl
        );
    }

    /// Zero camber ⇒ no doublet and no side force.
    #[test]
    fn zero_camber_is_inert() {
        let cp = solve_centerplane(1.0, 2.0, |_, _| 0.0, 4, 8);
        assert!(cp.side_force.abs() < 1e-12);
        assert!(cp.doublet.iter().all(|d| d.mu.abs() < 1e-12));
    }

    /// Side force is proportional to the drift angle (linear theory).
    #[test]
    fn side_force_is_linear_in_drift() {
        let s1 = solve_centerplane(1.0, 2.0, drift_slope(0.02), 5, 12).side_force;
        let s2 = solve_centerplane(1.0, 2.0, drift_slope(0.04), 5, 12).side_force;
        assert!(s1.abs() > 0.0);
        assert!((s2 - 2.0 * s1).abs() < 1e-9 * s1.abs(), "not linear: {s1}, {s2}");
    }

    /// The doublet density accumulates monotonically from the leading edge
    /// (positive-loading flat plate), is largest near the waterline root, and
    /// relieves toward the keel tip.
    #[test]
    fn doublet_is_root_loaded_and_tip_relieved() {
        let (nx, nz) = (6, 16);
        let cp = solve_centerplane(1.0, 3.0, drift_slope(0.05), nx, nz);
        // Trailing-edge μ per strip = last sample of each spanwise block.
        let te_mu = |strip: usize| cp.doublet[strip * nx + (nx - 1)].mu.abs();
        // Within a strip, μ grows monotonically along x.
        for strip in 0..nz {
            for ix in 1..nx {
                let a = cp.doublet[strip * nx + ix - 1].mu.abs();
                let b = cp.doublet[strip * nx + ix].mu.abs();
                assert!(b >= a - 1e-12, "μ not monotone in x at strip {strip}");
            }
        }
        // Root strip (nearest waterline, strip 0) carries more than the keel
        // tip strip (strip nz-1).
        assert!(
            te_mu(0) > te_mu(nz - 1),
            "root {} not > tip {}",
            te_mu(0),
            te_mu(nz - 1)
        );
        // Keel tip is substantially relieved.
        assert!(te_mu(nz - 1) < 0.6 * te_mu(0), "keel not relieved");
    }

    /// The doublet field feeds a finite free-wave amplitude whose dipole weight
    /// `√(λ²−1)` vanishes at λ = 1 (a wave running dead ahead is blind to
    /// asymmetry), matching the analysis behind the strip closure.
    #[test]
    fn free_wave_amplitude_is_finite_and_grows_from_head_sea() {
        let cp = solve_centerplane(1.0, 2.0, drift_slope(0.05), 8, 16);
        let nu = 1.0;
        let g1 = cp.doublet_free_wave_amplitude(nu, 1.0);
        let g2 = cp.doublet_free_wave_amplitude(nu, 2.5);
        assert!(g1.0.is_finite() && g1.1.is_finite());
        // Dipole *contribution* ∝ √(λ²−1)·|G|; zero at λ=1, positive beyond.
        let contrib = |lambda: f64, g: (f64, f64)| {
            (lambda * lambda - 1.0).max(0.0).sqrt() * (g.0 * g.0 + g.1 * g.1).sqrt()
        };
        assert!(contrib(1.0, g1).abs() < 1e-12);
        assert!(contrib(2.5, g2) > 0.0);
    }
}

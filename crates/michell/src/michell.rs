//! Michell's thin-ship wave resistance integral.
//!
//! With ν = g/U², half-beam f(x, z), z downward from the waterline:
//!
//! ```text
//! R_w = (4 ρ g²)/(π U²) ∫_1^∞ (I² + J²) λ²/√(λ² − 1) dλ
//! I(λ) + i J(λ) = ∬ (∂f/∂x) exp(−ν λ² z) exp(i ν λ x) dx dz
//! ```
//!
//! (Tuck 1989; Dambrine, Pierre & Rousseaux 2016.) Substituting λ = sec θ
//! turns the outer integral into `∫_0^{π/2} (I² + J²) sec³θ dθ`, removing the
//! integrable singularity at λ = 1.
//!
//! Because the hull is piecewise polynomial, I and J are evaluated **exactly**
//! per knot span via the closed-form moments in [`crate::moments`]; the only
//! numerical error lives in the smooth outer θ-integral, which is integrated
//! with Gauss–Legendre panels sized to the local oscillation rate and then
//! refined until the requested tolerance is met.
//!
//! ## Asymmetric hulls
//!
//! A hull built with [`Hull::new_asymmetric`] is split into a symmetric
//! thickness part `f_sym = (f₊ + f₋)/2` and an antisymmetric camber part
//! `f_a = (f₊ − f₋)/2`. The thickness part is the classical **source** system
//! above; the camber part adds a centreplane **y-dipole** system whose
//! free-wave amplitude is the same inner integral over `∂f_a/∂x` weighted by
//! the transverse wavenumber (see [`dipole_weight`]). Because the source
//! amplitude is even in θ and the dipole amplitude is odd, their cross term
//! integrates to zero over the Kelvin fan, so
//!
//! ```text
//! R_w = R_source(f_sym) + R_dipole(f_a)
//! ```
//!
//! with the symmetric hull (`f_a ≡ 0`) recovering classical Michell exactly.
//!
//! The [`multihull_wave_resistance_with`] path fixes the dipole *magnitude* with
//! the crude prescribed strip closure `μ = 2U f_a` (read qualitatively — see
//! [`DIPOLE_WEIGHT_C`]). [`asymmetric_wave_resistance_lifting`] instead *solves*
//! the centreplane lifting-surface problem ([`crate::centerplane`]) for the
//! doublet density and forms the same dipole from the solved `μ` — the
//! physically grounded magnitude, sharing this term's normalisation exactly (the
//! strip closure is what it reduces to when `μ = 2U f_a`).

use crate::centerplane::CenterplaneSolution;
use crate::conditions::Conditions;
use crate::error::{Error, Result};
use crate::hull::Hull;
use crate::moments::{exp_moments, exp_moments_complex, osc_moments, C64};
use crate::quadrature::gauss_legendre;
use std::f64::consts::{FRAC_PI_2, PI};

/// Options for the outer-integral quadrature.
#[derive(Debug, Clone, Copy)]
pub struct WaveOptions {
    /// Target relative tolerance on the resistance.
    pub rel_tol: f64,
    /// Maximum number of panel-halving refinement passes.
    pub max_refinements: usize,
}

impl Default for WaveOptions {
    fn default() -> Self {
        WaveOptions {
            rel_tol: 1e-5,
            max_refinements: 4,
        }
    }
}

/// Numerical method used to produce a [`WaveResistance`] result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WaveMethod {
    /// General real-axis marching quadrature, valid for every supported hull.
    GeneralMarcher,
    /// Low-Froude endpoint reduction on a steepest-descent contour.
    EndpointReduction,
}

impl WaveMethod {
    /// Stable lower-case name for diagnostics and file formats.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::GeneralMarcher => "general_marcher",
            Self::EndpointReduction => "endpoint_reduction",
        }
    }

    /// Stable numeric code for numeric-only sweep tables.
    pub const fn code(self) -> u32 {
        match self {
            Self::GeneralMarcher => 0,
            Self::EndpointReduction => 1,
        }
    }
}

/// Convergence outcome for a [`WaveResistance`] result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WaveOutcome {
    /// The method's stopping and requested-refinement criteria were met.
    Converged,
    /// The marcher reached its hard λ tail cap before a quiet tail was found.
    TailCap,
    /// The marcher exhausted its per-pass evaluation budget.
    EvalCap,
    /// Panel refinement was exhausted before `rel_tol` was met.
    RefinementCap,
}

impl WaveOutcome {
    /// Stable lower-case name for diagnostics and file formats.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Converged => "converged",
            Self::TailCap => "tail_cap",
            Self::EvalCap => "eval_cap",
            Self::RefinementCap => "refinement_cap",
        }
    }

    /// Stable numeric code for numeric-only sweep tables.
    pub const fn code(self) -> u32 {
        match self {
            Self::Converged => 0,
            Self::TailCap => 1,
            Self::EvalCap => 2,
            Self::RefinementCap => 3,
        }
    }

    /// Whether the method met all of its convergence criteria.
    pub const fn is_converged(self) -> bool {
        matches!(self, Self::Converged)
    }
}

/// Wave resistance result with quadrature diagnostics.
#[derive(Debug, Clone, Copy)]
pub struct WaveResistance {
    /// Wave resistance R_w [N].
    pub resistance: f64,
    /// Method-specific relative error diagnostic.
    ///
    /// For [`WaveMethod::GeneralMarcher`] this is the larger of the last-pass
    /// refinement difference and a power-law extrapolation of the terminating
    /// quiet window. It remains a heuristic rather than a rigorous bound, but
    /// accounts for the long aggregate tail that a small local window alone
    /// can hide at low Froude number.
    /// For [`WaveMethod::EndpointReduction`] it combines a rigorous bound on
    /// omitted submerged endpoints with an empirical contour-quadrature
    /// estimate.
    pub est_rel_error: f64,
    /// Total number of inner-integral evaluations performed, or transformed
    /// kernel nodes for an accepted low-Froude endpoint reduction.
    pub inner_evaluations: usize,
    /// Largest λ = sec θ reached before truncation. This is infinity when the
    /// low-Froude steepest-descent contour evaluates the infinite interval
    /// without real-axis truncation.
    pub max_lambda: f64,
    /// Numerical route used to produce this result.
    pub method: WaveMethod,
    /// Whether that route converged or stopped at a safety/refinement cap.
    pub outcome: WaveOutcome,
}

/// Michell resistance and its exact first derivative with respect to the
/// symmetric hull's B-spline control net (row-major, matching
/// [`crate::BSplineSurface::control`]).
#[derive(Debug, Clone)]
pub struct WaveResistanceGradient {
    pub wave: WaveResistance,
    /// `∂R_w/∂Pᵢⱼ` [N/m] for each half-breadth control value.
    pub control_gradient: Vec<f64>,
    /// Inner-integral evaluations in the one reverse pass, separate from the
    /// primal evaluations reported by [`WaveResistance::inner_evaluations`].
    pub gradient_evaluations: usize,
}

/// Control-net derivative for one fleet member.
#[derive(Debug, Clone)]
pub enum ControlNetGradient {
    /// One control net represents a port/starboard-symmetric half-breadth.
    Symmetric(Vec<f64>),
    /// Independent physical half-breadth control nets for an asymmetric hull.
    Asymmetric {
        /// Derivatives with respect to the port half-breadth controls.
        port: Vec<f64>,
        /// Derivatives with respect to the starboard half-breadth controls.
        starboard: Vec<f64>,
    },
}

/// Derivatives with respect to one member's rigid placement [N/m].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PlacementGradient {
    /// ∂R_w/∂x for the longitudinal offset.
    pub longitudinal: f64,
    /// ∂R_w/∂y for the transverse offset.
    pub transverse: f64,
}

/// Exact derivatives for one member of a multihull resistance result.
#[derive(Debug, Clone)]
pub struct MemberWaveResistanceGradient {
    pub placement: PlacementGradient,
    pub control: ControlNetGradient,
}

/// Michell resistance and exact derivatives for every multihull member.
#[derive(Debug, Clone)]
pub struct MultihullWaveResistanceGradient {
    pub wave: WaveResistance,
    pub members: Vec<MemberWaveResistanceGradient>,
    /// Outer nodes in the final reverse pass, separate from the primal count.
    pub gradient_evaluations: usize,
}

/// Position of one hull of a multihull, in the global (fleet) frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Placement {
    /// Longitudinal shift **added** to the hull's own x coordinates [m]
    /// (0 keeps the coordinates from the hull's file).
    pub x: f64,
    /// Transverse position of the hull's centerplane [m].
    pub y: f64,
}

/// Compute Michell wave resistance with default options.
pub fn wave_resistance(hull: &Hull, cond: &Conditions) -> Result<WaveResistance> {
    wave_resistance_with(hull, cond, &WaveOptions::default())
}

/// Compute Michell wave resistance with explicit quadrature options.
///
/// For a symmetric hull at sufficiently low Froude number this first tries the
/// waterline-endpoint/steepest-descent reduction. It is accepted only when its
/// omitted-endpoint bound and contour estimate meet `opts.rel_tol`; all other
/// cases retain the general marching quadrature.
pub fn wave_resistance_with(
    hull: &Hull,
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    if opts.rel_tol.is_finite() && opts.rel_tol > 0.0 {
        if let Ok(reduced) = crate::low_froude::low_froude_wave_resistance(hull, cond) {
            if reduced.est_rel_error <= opts.rel_tol {
                return Ok(WaveResistance {
                    resistance: reduced.resistance,
                    est_rel_error: reduced.est_rel_error,
                    inner_evaluations: reduced.kernel_evaluations,
                    max_lambda: f64::INFINITY,
                    method: WaveMethod::EndpointReduction,
                    outcome: WaveOutcome::Converged,
                });
            }
        }
    }
    multihull_wave_resistance_with(&[(hull, Placement::default())], cond, opts)
}

/// Wave resistance and exact B-spline control-net gradient with default
/// quadrature options. See [`wave_resistance_gradient_with`].
pub fn wave_resistance_gradient(
    hull: &Hull,
    cond: &Conditions,
) -> Result<WaveResistanceGradient> {
    wave_resistance_gradient_with(hull, cond, &WaveOptions::default())
}

/// Wave resistance and exact B-spline control-net gradient.
///
/// The inner amplitude is linear in every control value and the resistance is
/// quadratic in that amplitude. This routine reverse-accumulates
/// `2 Re(conj(F) · ∂F/∂Pᵢⱼ)` on the final adaptive outer-quadrature pass, then
/// applies the transpose of the exact control-to-span-polynomial map. Its cost
/// is one primal convergence plus one reverse pass, independent of the number
/// of controls; no finite differencing is used.
///
/// The derivative holds knots, degrees, speed, and fluid conditions fixed. It
/// is currently defined for a symmetric [`Hull`]; an asymmetric hull has two
/// physical control nets and should use
/// [`multihull_wave_resistance_gradient_with`] to receive both.
///
/// This gradient deliberately uses [`WaveMethod::GeneralMarcher`] even where
/// [`wave_resistance_with`] would dispatch the primal value to the endpoint
/// reduction. At low Froude number the returned primal can therefore differ
/// from that separate API; differentiating the endpoint reduction is not yet
/// implemented.
pub fn wave_resistance_gradient_with(
    hull: &Hull,
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistanceGradient> {
    if hull.is_asymmetric() {
        return Err(Error::InvalidGeometry(
            "wave_resistance_gradient requires a symmetric hull".into(),
        ));
    }
    let fleet = multihull_wave_resistance_gradient_with(
        &[(hull, Placement::default())],
        cond,
        opts,
    )?;
    let member = fleet.members.into_iter().next().unwrap();
    let ControlNetGradient::Symmetric(control_gradient) = member.control else {
        unreachable!("symmetric hull returned asymmetric controls")
    };
    Ok(WaveResistanceGradient {
        wave: fleet.wave,
        control_gradient,
        gradient_evaluations: fleet.gradient_evaluations,
    })
}

/// Exact multihull control-net and placement gradients with default options.
pub fn multihull_wave_resistance_gradient(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
) -> Result<MultihullWaveResistanceGradient> {
    multihull_wave_resistance_gradient_with(members, cond, &WaveOptions::default())
}

/// Exact multihull control-net and placement gradients.
///
/// Member phases are differentiated analytically, so the placement entries are
/// derivatives with respect to each raw [`Placement::x`] and [`Placement::y`].
/// Asymmetric members cover the same source plus prescribed strip-closure
/// camber/dipole terms as [`multihull_wave_resistance_with`]; the chain rule
/// returns separate port and starboard control nets. This API does not
/// differentiate the optional solved-lifting closure.
///
/// Like [`wave_resistance_gradient_with`], the primal and reverse pass both
/// use [`WaveMethod::GeneralMarcher`] at every Froude number.
pub fn multihull_wave_resistance_gradient_with(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<MultihullWaveResistanceGradient> {
    validate_fleet(members, cond, opts)?;
    let u = cond.speed;
    let g = cond.gravity;
    let nu = g / (u * u);
    let coeff = 4.0 * cond.fluid.density * g * g / (PI * u * u);
    let (cx_ref, y_ref) = fleet_phase_refs(members);
    let params = fleet_outer_params(members, nu, cx_ref, y_ref, 0.0);

    let make_members = || {
        members
            .iter()
            .map(|(hull, placement)| SourceMember {
                inner: InnerIntegral::new(hull, nu),
                dx: hull.x_center() + placement.x - cx_ref,
                dy: placement.y - y_ref,
            })
            .collect::<Vec<_>>()
    };
    let mut primal_members = make_members();
    let (wave, frac) = run_outer_with_frac(&params, opts, coeff, |lambda| {
        superpose(&mut primal_members, nu, lambda)
    });

    let mut reverse_members = make_members();
    let mut source_coeff_adjoint: Vec<Vec<f64>> = members
        .iter()
        .map(|(hull, _)| vec![0.0; hull.fx_coeff().len()])
        .collect();
    let mut camber_coeff_adjoint: Vec<Option<Vec<f64>>> = members
        .iter()
        .map(|(hull, _)| hull.fx_a_coeff().map(|coeff| vec![0.0; coeff.len()]))
        .collect();
    let mut placement_gradient = vec![PlacementGradient::default(); members.len()];
    let mut contributions = vec![SourceContribution::default(); members.len()];

    const GL_N: usize = 16;
    let (gx, gw) = gauss_legendre(GL_N);
    let mut gradient_evaluations = 0usize;
    let reverse = integrate_outer(
        &params,
        frac,
        &gx,
        &gw,
        &mut |lambda, weight| {
            source_gradient_integrand(
                &mut reverse_members,
                nu,
                lambda,
                weight,
                &mut contributions,
                &mut source_coeff_adjoint,
                &mut camber_coeff_adjoint,
                &mut placement_gradient,
            )
        },
        &mut gradient_evaluations,
    );
    debug_assert!(
        (coeff * reverse.integral - wave.resistance).abs()
            <= 1e-12 * wave.resistance.abs().max(1.0)
    );
    debug_assert_eq!(reverse.max_lambda, wave.max_lambda);

    let member_gradients = members
        .iter()
        .enumerate()
        .map(|(index, (hull, _))| {
            let source = hull.fx_control_adjoint(&source_coeff_adjoint[index]);
            let control = match &camber_coeff_adjoint[index] {
                None => ControlNetGradient::Symmetric(
                    source.into_iter().map(|value| coeff * value).collect(),
                ),
                Some(camber_coeff) => {
                    let camber = hull.fx_control_adjoint(camber_coeff);
                    let port = source
                        .iter()
                        .zip(&camber)
                        .map(|(source, camber)| 0.5 * coeff * (source - camber))
                        .collect();
                    let starboard = source
                        .iter()
                        .zip(camber)
                        .map(|(source, camber)| 0.5 * coeff * (source + camber))
                        .collect();
                    ControlNetGradient::Asymmetric { port, starboard }
                }
            };
            MemberWaveResistanceGradient {
                placement: PlacementGradient {
                    longitudinal: coeff * placement_gradient[index].longitudinal,
                    transverse: coeff * placement_gradient[index].transverse,
                },
                control,
            }
        })
        .collect();

    Ok(MultihullWaveResistanceGradient {
        wave,
        members: member_gradients,
        gradient_evaluations,
    })
}

/// Wave resistance of a hull **heeled** by `heel` radians about its
/// longitudinal axis.
///
/// A heeled hull is asymmetric relative to the horizontal free surface, but its
/// keel swings off the earth-vertical centreplane, so it cannot be written as
/// port/starboard half-beams there. Thin-ship theory instead keeps the sources
/// on the ship's own (now tilted) centreplane: a strip at ship-depth `z` sits
/// at earth depth `z·cosφ` and transverse offset `−z·sinφ`, which turns the
/// vertical decay **complex**,
///
/// ```text
/// κ = νλ²·cosφ + i·νλ√(λ²−1)·sinφ,
/// ```
/// the imaginary part being the transverse-wavenumber phase of the tilt (the
/// same dipole coupling as an asymmetric hull, arising here from geometry). The
/// upright kernel is untouched; heel just swaps in a complex-`κ` variant of the
/// inner integral, so `heel = 0` reproduces [`wave_resistance`] exactly. The
/// result is even in `heel` (port and starboard heel are mirror images).
///
/// This captures the asymmetric **wave-making** of the tilted thickness
/// distribution — the leading heel effect. It does not include the lifting
/// side-force a heeled-and-yawed (drifting) hull develops; that is a separate
/// forcing into the centreplane lifting solve.
pub fn heel_wave_resistance(
    hull: &Hull,
    cond: &Conditions,
    heel: f64,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    // The single-hull case of the fleet path: one member at the origin, so the
    // placement phase is unity and the result is bit-for-bit the standalone
    // heeled amplitude averaged over the ±θ systems.
    multihull_heel_wave_resistance(&[(hull, Placement::default())], cond, heel, opts)
}

/// Combined wave resistance of a **heeled multihull** — a rigid platform (e.g.
/// a catamaran or trimaran) heeled by `heel` radians about its own longitudinal
/// axis.
///
/// A rigid heel rotates *every* demihull's centreplane by the same angle `φ`,
/// so the single scalar `heel` is applied to each member's tilted-centreplane
/// (complex-`κ`) kernel — the same one [`heel_wave_resistance`] uses. Rotating
/// about the distant platform axis decomposes into a rotation about each
/// demihull's own axis (captured by the kernel) plus a rigid translation of the
/// centreplane; the transverse part of that translation is exactly the
/// member's transverse offset `y_j`, carried by the usual placement phase. So a
/// caller who has already repositioned the demihulls to the heeled attitude
/// (e.g. via [`crate::float::heel_poses`], which also supplies the immersion)
/// and passes the same `φ` here gets both the per-hull tilt and the demihull
/// interference of the heeled fleet.
///
/// The demihull free-wave systems (each tilted, each at its heeled transverse
/// offset) superpose exactly as in [`multihull_wave_resistance_with`]:
///
/// ```text
/// A_± = Σ_j F_j^φ(λ, ±) · exp(iν(λ Δx_j ± λ√(λ²−1) y_j)),
/// R_w = (4ρg²/πU²) ∫₁^∞ ½(|A₊|² + |A₋|²) λ²/√(λ²−1) dλ.
/// ```
///
/// Unlike the single-hull case the result is **not** even in `φ`: for a fleet
/// that is not mirror-symmetric about its mean centreplane, or once the demihull
/// interference is included, port-down and starboard-down heel differ, so both
/// Kelvin half-systems are carried explicitly (they are not folded together).
/// `heel = 0` reproduces [`multihull_wave_resistance_with`] for a symmetric
/// fleet exactly.
///
/// Scope (inherited from [`heel_wave_resistance`]): this captures the asymmetric
/// **wave-making** of the tilted thickness distribution. It does not re-clip
/// each hull to the true tilted waterline (the emerging/submerging wedges — the
/// same approximation the hydrostatic heel sweep makes), nor include the lifting
/// side-force of a heeled-and-yawed (drifting) fleet.
pub fn multihull_heel_wave_resistance(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    heel: f64,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    validate_fleet(members, cond, opts)?;
    if !heel.is_finite() || heel.abs() >= FRAC_PI_2 {
        return Err(Error::InvalidGeometry(
            "heel angle must be finite with |heel| < 90°".into(),
        ));
    }
    let u = cond.speed;
    let g = cond.gravity;
    let nu = g / (u * u);
    let rho = cond.fluid.density;

    let (cx_ref, y_ref) = fleet_phase_refs(members);
    // The tilt spreads each centreplane transversely by up to T·|sinφ| beyond
    // its placement offset; widen the panel-sizing extent so the outer integral
    // resolves that heel-induced oscillation.
    let t_max = members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max);
    let params = fleet_outer_params(members, nu, cx_ref, y_ref, t_max * heel.sin().abs());

    let mut mem: Vec<HeelMember> = members
        .iter()
        .map(|(h, p)| HeelMember {
            inner: HeelInner::new(h, nu, heel),
            dx: h.x_center() + p.x - cx_ref,
            dy: p.y - y_ref,
        })
        .collect();

    let coeff = 4.0 * rho * g * g / (PI * u * u);
    Ok(run_outer(&params, opts, coeff, |lambda| {
        superpose(&mut mem, nu, lambda)
    }))
}

/// Combined wave resistance of several thin hulls (multihull), with default
/// options. See [`multihull_wave_resistance_with`].
pub fn multihull_wave_resistance(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
) -> Result<WaveResistance> {
    multihull_wave_resistance_with(members, cond, &WaveOptions::default())
}

/// Combined wave resistance of several thin hulls.
///
/// In thin-ship theory the far-field free-wave amplitudes superpose: hull `j`
/// at longitudinal offset `Δx_j` and transverse position `y_j` contributes
/// `F_j(λ) · exp(i ν (λ Δx_j ± λ √(λ²−1) y_j))` to the wave system
/// propagating at ±θ off the track, and the resistance integrand is the mean
/// of the two systems, `½ (|Σ_j …₊|² + |Σ_j …₋|²)` (the θ < 0 half of the
/// free-wave spectrum; the two differ only for fleets that are not
/// mirror-symmetric about their mean centerplane, e.g. staggered pairs). For
/// two identical hulls separated by `s` this reduces to the classical
/// catamaran interference factor `4 cos²(½ ν s λ √(λ²−1))`.
///
/// A member built with [`Hull::new_asymmetric`] additionally contributes a
/// centreplane y-dipole wave system from its antisymmetric half-beam `f_a`
/// (weighted by `dipole_weight`); this superposes on the source system exactly
/// like any other free-wave amplitude, so demihull asymmetry and demihull
/// interference are handled together. The dipole *magnitude* rests on an
/// approximate strip closure — see the `DIPOLE_WEIGHT_C` constant.
pub fn multihull_wave_resistance_with(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<WaveResistance> {
    validate_fleet(members, cond, opts)?;
    let u = cond.speed;
    let g = cond.gravity;
    let nu = g / (u * u);
    let rho = cond.fluid.density;

    let (cx_ref, y_ref) = fleet_phase_refs(members);
    let params = fleet_outer_params(members, nu, cx_ref, y_ref, 0.0);

    let mut mem: Vec<SourceMember> = members
        .iter()
        .map(|(h, p)| SourceMember {
            inner: InnerIntegral::new(h, nu),
            dx: h.x_center() + p.x - cx_ref,
            dy: p.y - y_ref,
        })
        .collect();

    let coeff = 4.0 * rho * g * g / (PI * u * u);
    Ok(run_outer(&params, opts, coeff, |lambda| {
        superpose(&mut mem, nu, lambda)
    }))
}

/// Grid resolution for the centreplane lifting solve behind
/// [`asymmetric_wave_resistance_lifting`]: `nx` streamwise panels along the
/// hull, `nz` vertical panels over the draft (physical half).
#[derive(Debug, Clone, Copy)]
pub struct LiftingGrid {
    pub nx: usize,
    pub nz: usize,
}

impl Default for LiftingGrid {
    fn default() -> Self {
        LiftingGrid { nx: 48, nz: 16 }
    }
}

/// Wave resistance of a single **asymmetric** hull with the camber (dipole)
/// system taken from a *solved* centreplane lifting distribution rather than
/// the prescribed strip closure `μ = 2U f_a` that
/// [`multihull_wave_resistance_with`] uses.
///
/// The symmetric thickness part is the usual Michell source integral. The
/// antisymmetric part runs the [`crate::centerplane`] vortex-lattice solve for
/// the doublet density `μ(x, z)` and forms the dipole free-wave amplitude
///
/// ```text
/// A_d(λ) = i · νλ√(λ²−1) · ½ · G(λ),   G(λ) = ∬ (μ/U) e^{−νλ²z} e^{iνλx} dx dz.
/// ```
///
/// Substituting the strip closure `μ/U = 2 f_a` into this (via integration by
/// parts in x) reproduces the `dipole_weight` term exactly, so the two paths
/// share one normalisation; the solved `μ` has a different chordwise/vertical
/// shape than `2 f_a`, giving a physically grounded dipole of the same order
/// (their ratio is speed-dependent; for a 2-D flat plate `μ_lift/μ_strip =
/// π/2`). The source amplitude is even in θ and the dipole odd, so they add with
/// no cross term: `R_w = R_source(f_sym) + R_dipole(μ)`.
///
/// Errors if `hull` is symmetric (use [`wave_resistance`] instead). This is the
/// single-hull case of [`multihull_wave_resistance_lifting`].
pub fn asymmetric_wave_resistance_lifting(
    hull: &Hull,
    cond: &Conditions,
    opts: &WaveOptions,
    grid: LiftingGrid,
) -> Result<WaveResistance> {
    if !hull.is_asymmetric() {
        return Err(Error::InvalidGeometry(
            "asymmetric_wave_resistance_lifting requires a hull built with \
             Hull::new_asymmetric"
                .into(),
        ));
    }
    multihull_wave_resistance_lifting(&[(hull, Placement::default())], cond, opts, grid)
}

/// One fleet member's precomputed inner integral and (for an asymmetric hull)
/// its solved centreplane dipole distribution.
struct LiftMember<'h> {
    inner: InnerIntegral<'h>,
    /// The solved doublet distribution and the hull length (to re-centre the
    /// dipole phase at the hull midpoint). `None` for a symmetric member.
    dipole: Option<(CenterplaneSolution, f64)>,
    dx: f64,
    dy: f64,
}

/// Combined wave resistance of a fleet whose asymmetric members carry a
/// **solved** centreplane-lifting dipole (rather than the strip closure of
/// [`multihull_wave_resistance_with`]).
///
/// Each asymmetric member's centreplane lifting problem is solved once for its
/// doublet density `μ_j`; the dipole amplitude `A_{d,j} = i·νλ√(λ²−1)·½·G_j`
/// then superposes on the source system with the member's placement phase,
/// exactly like any other free-wave amplitude:
///
/// ```text
/// A₊ = Σ_j (F_j + A_{d,j}) e^{iν(λΔx_j + λ√(λ²−1) y_j)},   A₋ likewise with −y_j,
/// R_w = (4ρg²/πU²) ∫₁^∞ ½(|A₊|² + |A₋|²) λ²/√(λ²−1) dλ.
/// ```
///
/// This is the far-field coupling: each member's lifting problem is solved
/// independently (its own free-surface image), and their dipole **wave**
/// systems interfere through the superposition — so an asymmetric-demihull
/// catamaran's demihull interference is captured. Near-field lifting
/// cross-induction between close demihulls is not (that needs one coupled
/// lifting solve over all centreplanes).
///
/// Symmetric members contribute source-only (no dipole); a fleet of one
/// asymmetric hull is [`asymmetric_wave_resistance_lifting`].
pub fn multihull_wave_resistance_lifting(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
    grid: LiftingGrid,
) -> Result<WaveResistance> {
    validate_fleet(members, cond, opts)?;
    if grid.nx == 0 || grid.nz == 0 {
        return Err(Error::InvalidConditions(
            "lifting grid must have at least one panel each way".into(),
        ));
    }
    let u = cond.speed;
    let g = cond.gravity;
    let nu = g / (u * u);
    let rho = cond.fluid.density;

    let (cx_ref, y_ref) = fleet_phase_refs(members);
    let params = fleet_outer_params(members, nu, cx_ref, y_ref, 0.0);

    // Solve each asymmetric member's centreplane once (up front, not per λ).
    let mut mem: Vec<LiftMember> = members
        .iter()
        .map(|(h, p)| {
            let dipole = if h.is_asymmetric() {
                let x0 = h.surface().x_domain().0;
                let sol = crate::centerplane::solve_centerplane(
                    h.length(),
                    h.draft(),
                    |xl, z| h.eval_fx_a(x0 + xl, z),
                    grid.nx,
                    grid.nz,
                );
                Some((sol, h.length()))
            } else {
                None
            };
            LiftMember {
                inner: InnerIntegral::new(h, nu),
                dipole,
                dx: h.x_center() + p.x - cx_ref,
                dy: p.y - y_ref,
            }
        })
        .collect();

    let coeff = 4.0 * rho * g * g / (PI * u * u);
    Ok(run_outer(&params, opts, coeff, |lambda| {
        superpose(&mut mem, nu, lambda)
    }))
}

/// The Michell inner integrals `(I(λ), J(λ))` — the free-wave amplitude
/// functions — evaluated exactly (per-span closed forms) for `λ >= 1`.
///
/// Phases are taken relative to the hull's x-midpoint, so I and J individually
/// depend on that (physically irrelevant) choice of origin while `I² + J²`
/// does not.
pub fn inner_integrals(hull: &Hull, cond: &Conditions, lambda: f64) -> Result<(f64, f64)> {
    cond.validate()?;
    if !(lambda.is_finite() && lambda >= 1.0) {
        return Err(Error::InvalidConditions(format!(
            "lambda must be >= 1, got {lambda}"
        )));
    }
    let nu = cond.gravity / (cond.speed * cond.speed);
    let mut inner = InnerIntegral::new(hull, nu);
    let f = inner.eval(lambda);
    Ok((f.re, f.im))
}

/// Closure constant for the antisymmetric (dipole) free-wave amplitude — see
/// [`dipole_weight`]. `1.0` is the value self-consistent with this crate's
/// source normalisation *under the prescribed strip closure* `μ = 2U f_a` (the
/// naive antisymmetric analogue of Michell's `σ = 2U ∂f_sym/∂x`).
///
/// **This magnitude is approximate.** Unlike the source strength — which the
/// boundary condition fixes *pointwise*, because a source sheet's normal-
/// velocity *jump* equals its local strength — the dipole density is fixed by
/// the *mean* of the two-sided normal velocities, and a doublet sheet's mean
/// normal velocity is a hypersingular (finite-part) integral of `μ`. So `μ` is
/// genuinely **non-local** in `∂f_a/∂x`; no exact pointwise closure exists, and
/// the rigorous density solves a hypersingular Fredholm equation of the first
/// kind (Kaklis & Papanikolaou; 21st Symp. Naval Hydro., Appendix A). `μ = 2U
/// f_a` is the crudest strip estimate of that solve — the *structure* below
/// (weight, θ-parity, additive separation) is exact, but the overall constant
/// should be validated against a reference before the dipole magnitude is
/// trusted quantitatively.
const DIPOLE_WEIGHT_C: f64 = 1.0;

/// Real spectral weight relating the antisymmetric dipole amplitude to the
/// companion source integral `G(λ) = ∬ (∂f_a/∂x) e^{−κz} e^{iνλx} dx dz`.
///
/// A y-normal centreplane doublet of density `μ(x,z)` radiates a free wave
/// carrying one extra factor of the **transverse wavenumber**
/// `k_y = νλ√(λ²−1)` relative to a source:
/// `A_d(θ) = i k_y ∬ μ e^{−κz} e^{iνλx} dx dz`. With the prescribed closure
/// `μ = 2U f_a`, `f_a` closing at bow and stern, integration by parts in `x`
/// rewrites `∬ f_a e^{iνλx} = (i/νλ) G`, and the two factors of `i` collapse to
/// a **real** weight on `G` in this crate's normalisation:
///
/// ```text
/// A_d(λ) = i · k_y · (2U) · (i/νλ) G(λ) / (2U)   (÷2U → crate amplitude units)
///        = −√(λ²−1) · DIPOLE_WEIGHT_C · G(λ).
/// ```
///
/// The sign is carried in [`multihull_wave_resistance_with`] (the +θ system
/// gets `−w·G`); this function returns the magnitude `w = √(λ²−1)·C`. It
/// vanishes as `λ → 1` (θ → 0): a wave running dead ahead carries no transverse
/// wavenumber and is blind to side-to-side asymmetry, so `R_dipole` is fed
/// entirely by the diverging (large-λ) part of the spectrum.
#[inline]
pub(crate) fn dipole_weight(lambda: f64) -> f64 {
    DIPOLE_WEIGHT_C * (lambda * lambda - 1.0).max(0.0).sqrt()
}

/// Geometry-derived phase-rate parameters for the outer quadrature.
struct OuterParams {
    nu: f64,
    /// Half-extent of the whole fleet about its longitudinal phase centre.
    x_half: f64,
    /// Largest transverse offset from the fleet's phase centre.
    y_half: f64,
    /// Deepest draft in the fleet.
    t_max: f64,
}

#[derive(Debug, Clone, Copy)]
struct OuterLimits {
    lambda_hard_cap: f64,
    max_evals_per_pass: usize,
}

const DEFAULT_OUTER_LIMITS: OuterLimits = OuterLimits {
    lambda_hard_cap: 1e4,
    max_evals_per_pass: 4_000_000,
};

struct OuterPass {
    integral: f64,
    /// Estimated integral beyond `max_lambda` from the terminating quiet
    /// window. Infinite when the pass stopped at a safety cap.
    tail_abs_estimate: f64,
    max_lambda: f64,
    outcome: WaveOutcome,
}

/// Marching-panel Gauss–Legendre integration of
/// `∫_0^{π/2} |A(sec θ)|² sec³θ dθ`, where `amp_sq` supplies the combined
/// `|A|²` of the fleet's two (±θ) wave systems.
fn integrate_outer(
    params: &OuterParams,
    frac: f64,
    gx: &[f64],
    gw: &[f64],
    amp_sq: &mut impl FnMut(f64, f64) -> f64,
    evals: &mut usize,
) -> OuterPass {
    integrate_outer_with_limits(
        params,
        frac,
        gx,
        gw,
        amp_sq,
        evals,
        DEFAULT_OUTER_LIMITS,
    )
}

fn integrate_outer_with_limits(
    params: &OuterParams,
    frac: f64,
    gx: &[f64],
    gw: &[f64],
    amp_sq: &mut impl FnMut(f64, f64) -> f64,
    evals: &mut usize,
    limits: OuterLimits,
) -> OuterPass {
    const GL_N: usize = 16;
    /// Truncate once a full quiet window contributes below this fraction.
    const STOP_REL: f64 = 1e-9;
    /// Width of the quiet window in accumulated phase — several full periods
    /// of the oscillating integrand, so a trough of cos²-type oscillation
    /// (width < π) can never trigger truncation on its own. Phase-based, so
    /// the criterion is independent of the panel-refinement level.
    const STOP_WINDOW_PHASE: f64 = 8.0 * PI;

    let nu = params.nu;
    let (x_half, y_half, t_max) = (params.x_half, params.y_half, params.t_max);

    // Local phase rate of |A|² in θ: the x-oscillation contributes
    // 2 ν x_half d(sec θ)/dθ, the z-decay envelope 2 ν T d(sec²θ)/dθ, and the
    // transverse separation phase ν y λ√(λ²−1) = ν y sec θ tan θ contributes
    // 2 ν y_half d(sec θ tan θ)/dθ = 2 ν y_half sec θ (sec²θ + tan²θ).
    let rate = |sec: f64, tan: f64| -> f64 {
        2.0 * nu * sec * tan * (x_half + t_max * sec)
            + 2.0 * nu * y_half * sec * (sec * sec + tan * tan)
            + 4.0
    };
    // Near θ = 0 the longitudinal phase grows like ν x_half θ², so cap the
    // first panels at one period of that quadratic phase (the transverse
    // phase is linear near 0 and already covered by rate(0)).
    let cap = (2.0 * PI / (2.0 * nu * x_half).sqrt().max(1.0)).min(0.12);

    let mut theta = 0.0f64;
    let mut total = 0.0f64;
    let mut window_sum = 0.0f64;
    let mut window_phase = 0.0f64;
    let mut window_lambda_start = 1.0f64;
    let mut pass_evals = 0usize;
    let mut lambda = 1.0f64;
    let mut outcome = WaveOutcome::TailCap;
    let mut tail_abs_estimate = f64::INFINITY;
    while theta < FRAC_PI_2 - 1e-12 {
        let (sin_theta, cos_theta) = theta.sin_cos();
        let sec_theta = 1.0 / cos_theta;
        let local_rate = rate(sec_theta, sin_theta * sec_theta);
        let dt = (frac * 2.0 * PI / local_rate)
            .min(frac * cap)
            .min(FRAC_PI_2 - theta)
            .max(1e-15);
        let half = dt / 2.0;
        let mid = theta + half;
        let mut panel = 0.0;
        for (i, &xi) in gx.iter().enumerate() {
            let th = mid + half * xi;
            let sec = 1.0 / th.cos();
            let sec_cubed = sec * sec * sec;
            panel += gw[i] * amp_sq(sec, half * gw[i] * sec_cubed) * sec_cubed;
        }
        panel *= half;
        total += panel;
        pass_evals += GL_N;
        theta += dt;
        lambda = 1.0 / theta.cos();

        // Truncation: only past λ = 2, and only when an entire window of
        // accumulated oscillation phase contributed negligibly.
        if lambda > 2.0 {
            if window_phase == 0.0 {
                window_lambda_start = sec_theta;
            }
            window_sum += panel;
            window_phase += local_rate * dt;
            if window_phase >= STOP_WINDOW_PHASE {
                // Every piecewise-polynomial hull amplitude is O(λ⁻³) or
                // faster at the waterline, so the transformed resistance
                // density is O(λ⁻⁵) and its remaining integral is
                // asymptotically one quarter of the local density times λ.
                // Recover that density from the full phase window rather
                // than treating one tiny low-Froude window as the tail.
                let lambda_width = (lambda - window_lambda_start).max(f64::MIN_POSITIVE);
                let extrapolation = (lambda / (4.0 * lambda_width)).max(1.0);
                const TAIL_SAFETY: f64 = 1.25;
                let candidate_tail = TAIL_SAFETY * window_sum.abs() * extrapolation;
                if window_sum.abs() <= STOP_REL * total.abs() + f64::MIN_POSITIVE {
                    tail_abs_estimate = candidate_tail;
                    outcome = WaveOutcome::Converged;
                    break;
                }
                window_sum = 0.0;
                window_phase = 0.0;
            }
        }
        if lambda > limits.lambda_hard_cap {
            outcome = WaveOutcome::TailCap;
            break;
        }
        if pass_evals >= limits.max_evals_per_pass {
            outcome = WaveOutcome::EvalCap;
            break;
        }
    }
    *evals += pass_evals;
    OuterPass {
        integral: total,
        tail_abs_estimate,
        max_lambda: lambda,
        outcome,
    }
}

// ---------------------------------------------------------------------------
// Shared multihull driver
//
// Every wave-resistance path integrates the same outer integral of the fleet's
// combined free-wave amplitude; they differ only in how each member's ±θ
// amplitude is formed. `MemberWave` captures that per-member amplitude;
// `superpose` shares the placement-phase superposition and `run_outer` the
// panel refinement and result assembly.
// ---------------------------------------------------------------------------

/// One fleet member's contribution to the far-field free-wave amplitude.
trait MemberWave {
    /// `(A₊, A₋)` — the amplitudes this member carries into the +θ and −θ wave
    /// systems at `λ` (before the placement phase).
    fn amps(&mut self, nu: f64, lambda: f64) -> (C64, C64);
    /// `(dx, dy)`: longitudinal and transverse offsets from the fleet phase
    /// centre. The placement phase is `exp(iν(λ dx ± λ√(λ²−1) dy))`.
    fn offsets(&self) -> (f64, f64);
}

/// Combined `½(|A₊|² + |A₋|²)` of the fleet's two Kelvin half-systems at `λ`,
/// summing each member's amplitude with its placement phase. Shared by every
/// multihull path — the members differ only in [`MemberWave::amps`].
fn superpose<M: MemberWave>(members: &mut [M], nu: f64, lambda: f64) -> f64 {
    let kx = nu * lambda;
    let ky = nu * lambda * (lambda * lambda - 1.0).max(0.0).sqrt();
    let mut plus = C64::ZERO;
    let mut minus = C64::ZERO;
    for m in members.iter_mut() {
        let (a_plus, a_minus) = m.amps(nu, lambda);
        if a_plus == C64::ZERO && a_minus == C64::ZERO {
            continue;
        }
        let (dx, dy) = m.offsets();
        plus = plus + a_plus * C64::cis(kx * dx + ky * dy);
        minus = minus + a_minus * C64::cis(kx * dx - ky * dy);
    }
    0.5 * (plus.abs_sq() + minus.abs_sq())
}

/// March the outer integral to the requested tolerance and assemble the
/// [`WaveResistance`]; `coeff = 4ρg²/(πU²)` is the Michell prefactor. Shared by
/// every wave-resistance entry point.
fn run_outer_with_frac(
    params: &OuterParams,
    opts: &WaveOptions,
    coeff: f64,
    amp_sq: impl FnMut(f64) -> f64,
) -> (WaveResistance, f64) {
    run_outer_with_limits(params, opts, coeff, amp_sq, DEFAULT_OUTER_LIMITS)
}

fn run_outer_with_limits(
    params: &OuterParams,
    opts: &WaveOptions,
    coeff: f64,
    mut amp_sq: impl FnMut(f64) -> f64,
    limits: OuterLimits,
) -> (WaveResistance, f64) {
    const GL_N: usize = 16;
    let (gx, gw) = gauss_legendre(GL_N);
    let mut sample = |lambda: f64, _weight: f64| amp_sq(lambda);
    let mut evals_total = 0usize;
    let mut frac = 1.0;
    let mut evals = 0usize;
    let mut pass = integrate_outer_with_limits(
        params,
        frac,
        &gx,
        &gw,
        &mut sample,
        &mut evals,
        limits,
    );
    evals_total += evals;
    let mut est_rel = f64::INFINITY;
    for _ in 0..opts.max_refinements {
        frac *= 0.5;
        let mut evals = 0usize;
        let refined = integrate_outer_with_limits(
            params,
            frac,
            &gx,
            &gw,
            &mut sample,
            &mut evals,
            limits,
        );
        evals_total += evals;
        let scale = refined.integral.abs().max(f64::MIN_POSITIVE);
        let refinement_rel = (refined.integral - pass.integral).abs() / scale;
        let tail_rel = refined.tail_abs_estimate / scale;
        est_rel = refinement_rel.max(tail_rel);
        // Panel halving cannot reduce the fixed quiet-window truncation. Once
        // its discretisation change is already below that tail floor, further
        // refinement is pure cost; stop and report RefinementCap below.
        let tail_limited = tail_rel > opts.rel_tol && refinement_rel <= tail_rel;
        pass = refined;
        if pass.outcome.is_converged() && (est_rel <= opts.rel_tol || tail_limited) {
            break;
        }
    }
    let outcome = if pass.outcome.is_converged() && est_rel > opts.rel_tol {
        WaveOutcome::RefinementCap
    } else {
        pass.outcome
    };
    (
        WaveResistance {
            resistance: coeff * pass.integral,
            est_rel_error: est_rel,
            inner_evaluations: evals_total,
            max_lambda: pass.max_lambda,
            method: WaveMethod::GeneralMarcher,
            outcome,
        },
        frac,
    )
}

fn run_outer(
    params: &OuterParams,
    opts: &WaveOptions,
    coeff: f64,
    amp_sq: impl FnMut(f64) -> f64,
) -> WaveResistance {
    run_outer_with_frac(params, opts, coeff, amp_sq).0
}

/// Shared validation for the multihull entry points.
fn validate_fleet(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    opts: &WaveOptions,
) -> Result<()> {
    cond.validate()?;
    if !(opts.rel_tol.is_finite() && opts.rel_tol > 0.0) {
        return Err(Error::InvalidConditions(
            "rel_tol must be finite and positive".into(),
        ));
    }
    if members.is_empty() {
        return Err(Error::InvalidConditions(
            "at least one hull is required".into(),
        ));
    }
    if members.iter().any(|(_, p)| !(p.x.is_finite() && p.y.is_finite())) {
        return Err(Error::InvalidConditions(
            "hull placements must be finite".into(),
        ));
    }
    Ok(())
}

/// Fleet phase references: the mean hull x-centre and mean transverse position.
/// Constant overall phase is irrelevant; centring the oscillatory arguments
/// keeps them small.
fn fleet_phase_refs(members: &[(&Hull, Placement)]) -> (f64, f64) {
    let n = members.len() as f64;
    let cx_ref = members.iter().map(|(h, p)| h.x_center() + p.x).sum::<f64>() / n;
    let y_ref = members.iter().map(|(_, p)| p.y).sum::<f64>() / n;
    (cx_ref, y_ref)
}

/// Outer-integral panel-sizing parameters for a fleet. `y_tilt` widens the
/// transverse extent for a heeled fleet (the tilt spreads each centreplane by
/// up to `T·|sinφ|` beyond its placement offset); pass 0 for an upright fleet.
fn fleet_outer_params(
    members: &[(&Hull, Placement)],
    nu: f64,
    cx_ref: f64,
    y_ref: f64,
    y_tilt: f64,
) -> OuterParams {
    OuterParams {
        nu,
        x_half: members
            .iter()
            .map(|(h, p)| (h.x_center() + p.x - cx_ref).abs() + h.x_half_extent())
            .fold(0.0, f64::max),
        y_half: members
            .iter()
            .map(|(_, p)| (p.y - y_ref).abs())
            .fold(0.0, f64::max)
            + y_tilt,
        t_max: members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max),
    }
}

/// Upright source (thickness) member carrying the strip-closure camber dipole.
struct SourceMember<'h> {
    inner: InnerIntegral<'h>,
    dx: f64,
    dy: f64,
}

#[derive(Clone, Copy)]
struct SourceContribution {
    phase_plus: C64,
    phase_minus: C64,
    carried_plus: C64,
    carried_minus: C64,
    camber_weight: f64,
}

impl Default for SourceContribution {
    fn default() -> Self {
        Self {
            phase_plus: C64::ZERO,
            phase_minus: C64::ZERO,
            carried_plus: C64::ZERO,
            carried_minus: C64::ZERO,
            camber_weight: 0.0,
        }
    }
}

#[inline]
fn conjugate(value: C64) -> C64 {
    C64::new(value.re, -value.im)
}

#[inline]
fn complex_dot(left: C64, right: C64) -> f64 {
    left.re * right.re + left.im * right.im
}

#[inline]
fn multiply_i(value: C64) -> C64 {
    C64::new(-value.im, value.re)
}

#[allow(clippy::too_many_arguments)]
fn source_gradient_integrand(
    members: &mut [SourceMember<'_>],
    nu: f64,
    lambda: f64,
    weight: f64,
    contributions: &mut [SourceContribution],
    source_coeff_adjoint: &mut [Vec<f64>],
    camber_coeff_adjoint: &mut [Option<Vec<f64>>],
    placement_gradient: &mut [PlacementGradient],
) -> f64 {
    let kx = nu * lambda;
    let ky = nu * lambda * (lambda * lambda - 1.0).max(0.0).sqrt();
    let camber_weight = dipole_weight(lambda);
    let mut plus = C64::ZERO;
    let mut minus = C64::ZERO;

    for (member, contribution) in members.iter_mut().zip(contributions.iter_mut()) {
        let (source, camber) = member.inner.eval_pair(lambda);
        let weighted_camber =
            camber.map_or(C64::ZERO, |value| value.scale(camber_weight));
        let phase_plus = C64::cis(kx * member.dx + ky * member.dy);
        let phase_minus = C64::cis(kx * member.dx - ky * member.dy);
        let carried_plus = (source - weighted_camber) * phase_plus;
        let carried_minus = (source + weighted_camber) * phase_minus;
        *contribution = SourceContribution {
            phase_plus,
            phase_minus,
            carried_plus,
            carried_minus,
            camber_weight,
        };
        plus = plus + carried_plus;
        minus = minus + carried_minus;
    }

    for index in 0..members.len() {
        let contribution = contributions[index];
        let dx_plus = multiply_i(contribution.carried_plus).scale(kx);
        let dx_minus = multiply_i(contribution.carried_minus).scale(kx);
        let dy_plus = multiply_i(contribution.carried_plus).scale(ky);
        let dy_minus = multiply_i(contribution.carried_minus).scale(-ky);
        placement_gradient[index].longitudinal +=
            weight * (complex_dot(plus, dx_plus) + complex_dot(minus, dx_minus));
        placement_gradient[index].transverse +=
            weight * (complex_dot(plus, dy_plus) + complex_dot(minus, dy_minus));

        // accumulate_coeff_adjoint differentiates |A|² and therefore carries
        // a factor of two. The half-system average supplies the compensating
        // 1/2 in these effective local amplitudes.
        let local_plus = plus * conjugate(contribution.phase_plus);
        let local_minus = minus * conjugate(contribution.phase_minus);
        let source_effective = (local_plus + local_minus).scale(0.5);
        let camber_effective = (local_minus - local_plus)
            .scale(0.5 * contribution.camber_weight);
        members[index].inner.accumulate_pair_coeff_adjoint(
            lambda,
            source_effective,
            camber_effective,
            weight,
            &mut source_coeff_adjoint[index],
            camber_coeff_adjoint[index].as_deref_mut(),
        );
    }

    0.5 * (plus.abs_sq() + minus.abs_sq())
}

impl MemberWave for SourceMember<'_> {
    fn amps(&mut self, _nu: f64, lambda: f64) -> (C64, C64) {
        // Source amplitude plus the strip-closure y-dipole. The dipole spectral
        // weight scales the source integral by the transverse wavenumber (see
        // [`dipole_weight`]); it is *odd* in θ (the +θ system carries −w·G, the
        // −θ system +w·G), so the ½(|A₊|²+|A₋|²) average in `superpose` cancels
        // the source–dipole cross term, leaving R_w = R_source + R_dipole.
        let (f, g) = self.inner.eval_pair(lambda);
        let wg = g.map_or(C64::ZERO, |g| g.scale(dipole_weight(lambda)));
        (f - wg, f + wg)
    }
    fn offsets(&self) -> (f64, f64) {
        (self.dx, self.dy)
    }
}

impl MemberWave for LiftMember<'_> {
    fn amps(&mut self, nu: f64, lambda: f64) -> (C64, C64) {
        let f = self.inner.eval(lambda); // source (thickness) amplitude
        let a_d = match &self.dipole {
            Some((sol, len)) => {
                let kx = nu * lambda;
                let wd = 0.5 * nu * lambda * (lambda * lambda - 1.0).max(0.0).sqrt();
                let (gre, gim) = sol.doublet_free_wave_amplitude(nu, lambda);
                // Re-centre the dipole phase at the hull midpoint so it shares
                // the source's origin (the solver uses x ∈ [0, L]).
                let gc = C64::new(gre, gim) * C64::cis(-kx * len * 0.5);
                // A_d = i · wd · G (odd in θ: +θ gets +A_d, −θ gets −A_d).
                C64::new(-gc.im, gc.re).scale(wd)
            }
            None => C64::ZERO,
        };
        (f + a_d, f - a_d)
    }
    fn offsets(&self) -> (f64, f64) {
        (self.dx, self.dy)
    }
}

/// Heeled member: the tilted-centreplane (complex-`κ`) source amplitude for the
/// ±θ systems, at a shared platform heel angle.
struct HeelMember<'h> {
    inner: HeelInner<'h>,
    dx: f64,
    dy: f64,
}

impl MemberWave for HeelMember<'_> {
    fn amps(&mut self, _nu: f64, lambda: f64) -> (C64, C64) {
        (self.inner.eval(lambda, 1.0), self.inner.eval(lambda, -1.0))
    }
    fn offsets(&self) -> (f64, f64) {
        (self.dx, self.dy)
    }
}

/// Exact (per-span closed-form) evaluation of I + iJ at λ.
pub(crate) struct InnerIntegral<'h> {
    hull: &'h Hull,
    nu: f64,
    /// Scratch: z-moments including the e^{−κ z0} shift, [n_spans_z][q+1].
    zm: Vec<f64>,
    /// Scratch: raw x-moments for one span.
    xm: Vec<C64>,
    /// Scratch: raw z-moments for one span.
    zm_raw: Vec<f64>,
    p: usize,
    q: usize,
}

impl<'h> InnerIntegral<'h> {
    pub(crate) fn new(hull: &'h Hull, nu: f64) -> Self {
        let p = hull.surface().degree_x();
        let q = hull.surface().degree_z();
        InnerIntegral {
            hull,
            nu,
            zm: vec![0.0; hull.spans_z().len() * (q + 1)],
            xm: Vec::with_capacity(p),
            zm_raw: Vec::with_capacity(q + 1),
            p,
            q,
        }
    }

    /// Source free-wave amplitude `I + iJ` at λ = sec θ. Phases use x relative
    /// to the hull midpoint (a pure phase factor on I + iJ that leaves
    /// |I + iJ|² unchanged, but keeps the oscillatory arguments as small as
    /// possible).
    pub(crate) fn eval(&mut self, lambda: f64) -> C64 {
        self.eval_pair(lambda).0
    }

    /// The source amplitude `F(λ)` and, for an asymmetric hull, the companion
    /// **dipole** amplitude `G(λ) = ∬ (∂f_a/∂x) e^{−κz} e^{iνλx} dx dz` from the
    /// antisymmetric half-beam. `G` is `None` for a symmetric hull. Both share
    /// the single z-moment pass; only the x-span accumulation is repeated with
    /// the second coefficient array, so an asymmetric evaluation costs little
    /// more than a symmetric one.
    pub(crate) fn eval_pair(&mut self, lambda: f64) -> (C64, Option<C64>) {
        let hull = self.hull;
        let kx = self.nu * lambda;
        let kappa = self.nu * lambda * lambda;
        if !self.fill_zm(kappa) {
            return (C64::ZERO, hull.fx_a_coeff().map(|_| C64::ZERO));
        }
        let f = self.accumulate(kx, hull.fx_coeff());
        let g = hull.fx_a_coeff().map(|c| self.accumulate(kx, c));
        (f, g)
    }

    /// Reverse-accumulate source and camber coefficient derivatives after an
    /// [`Self::eval_pair`] call at the same `lambda`. The effective amplitudes
    /// already contain the fleet superposition and the ±θ system average.
    fn accumulate_pair_coeff_adjoint(
        &mut self,
        lambda: f64,
        source_effective: C64,
        camber_effective: C64,
        weight: f64,
        source_coeff_adjoint: &mut [f64],
        camber_coeff_adjoint: Option<&mut [f64]>,
    ) {
        let kx = self.nu * lambda;
        self.accumulate_coeff_adjoint(
            kx,
            source_effective,
            weight,
            source_coeff_adjoint,
        );
        if let Some(adjoint) = camber_coeff_adjoint {
            debug_assert!(self.hull.fx_a_coeff().is_some());
            self.accumulate_coeff_adjoint(kx, camber_effective, weight, adjoint);
        }
    }

    /// Fill the per-span z-moments (shifted by the decay to each span start) at
    /// decay rate `kappa = ν λ²`. Returns `false` if every span has underflowed
    /// (the whole hull is below the exponential's support), in which case the
    /// amplitudes are zero.
    fn fill_zm(&mut self, kappa: f64) -> bool {
        let hull = self.hull;
        let q = self.q;
        let mut any = false;
        for (t, sz) in hull.spans_z().iter().enumerate() {
            let decay = (-kappa * sz.start).exp();
            if decay == 0.0 {
                for b in 0..=q {
                    self.zm[t * (q + 1) + b] = 0.0;
                }
                continue;
            }
            any = true;
            exp_moments(kappa, sz.len, q, &mut self.zm_raw);
            for b in 0..=q {
                self.zm[t * (q + 1) + b] = decay * self.zm_raw[b];
            }
        }
        any
    }

    /// Accumulate `∬ (Σ coeff·XᵃZᵇ) e^{−κz} e^{iνλx} dx dz` over all span pairs
    /// using the already-filled z-moments and the per-span x-oscillation
    /// moments. `coeff` is either the source `∂f/∂x` net or the antisymmetric
    /// `∂f_a/∂x` net (identical layout), so source and dipole amplitudes reuse
    /// this one kernel.
    fn accumulate(&mut self, kx: f64, coeff: &[f64]) -> C64 {
        let hull = self.hull;
        let (p, q) = (self.p, self.q);
        let nsz = hull.spans_z().len();
        let x_center = hull.x_center();
        let mut f = C64::ZERO;
        for (s, sx) in hull.spans_x().iter().enumerate() {
            osc_moments(kx, sx.len, p - 1, &mut self.xm);
            let phase = C64::cis(kx * (sx.start - x_center));
            let mut span_sum = C64::ZERO;
            for (a, &xma) in self.xm.iter().enumerate() {
                // g_a = Σ_t Σ_b c[s,t,a,b] · Zm[t][b]
                let mut g_a = 0.0;
                for t in 0..nsz {
                    let base = ((s * nsz + t) * p + a) * (q + 1);
                    let zrow = &self.zm[t * (q + 1)..(t + 1) * (q + 1)];
                    let crow = &coeff[base..base + q + 1];
                    for b in 0..=q {
                        g_a += crow[b] * zrow[b];
                    }
                }
                span_sum = span_sum + xma.scale(g_a);
            }
            f = f + phase * span_sum;
        }
        f
    }

    fn accumulate_coeff_adjoint(
        &mut self,
        kx: f64,
        amplitude: C64,
        weight: f64,
        coeff_adjoint: &mut [f64],
    ) {
        let hull = self.hull;
        let (p, q) = (self.p, self.q);
        let nsz = hull.spans_z().len();
        let x_center = hull.x_center();
        assert_eq!(coeff_adjoint.len(), hull.fx_coeff().len());
        for (s, sx) in hull.spans_x().iter().enumerate() {
            osc_moments(kx, sx.len, p - 1, &mut self.xm);
            let phase = C64::cis(kx * (sx.start - x_center));
            for (a, &xma) in self.xm.iter().enumerate() {
                let x_basis = phase * xma;
                for t in 0..nsz {
                    for b in 0..=q {
                        let index = ((s * nsz + t) * p + a) * (q + 1) + b;
                        let basis = x_basis.scale(self.zm[t * (q + 1) + b]);
                        coeff_adjoint[index] += weight
                            * 2.0
                            * (amplitude.re * basis.re + amplitude.im * basis.im);
                    }
                }
            }
        }
    }
}

/// Heeled-hull free-wave amplitude kernel — the upright [`InnerIntegral`] with a
/// **complex** vertical decay `κ = νλ²cosφ + i·νλ√(λ²−1)·sinφ` (see
/// [`heel_wave_resistance`]). The x-oscillation moments and the per-span
/// accumulate are identical to the upright kernel; only the z-moment is complex
/// (via [`exp_moments_complex`]), so the real path is entirely untouched.
pub(crate) struct HeelInner<'h> {
    hull: &'h Hull,
    nu: f64,
    cos_phi: f64,
    sin_phi: f64,
    /// Scratch: complex z-moments including the e^{−κ z0} shift, [n_spans_z][q+1].
    zm: Vec<C64>,
    /// Scratch: x-moments for one span.
    xm: Vec<C64>,
    /// Scratch: raw complex z-moments for one span.
    zm_raw: Vec<C64>,
    p: usize,
    q: usize,
}

impl<'h> HeelInner<'h> {
    fn new(hull: &'h Hull, nu: f64, heel: f64) -> Self {
        let p = hull.surface().degree_x();
        let q = hull.surface().degree_z();
        HeelInner {
            hull,
            nu,
            cos_phi: heel.cos(),
            sin_phi: heel.sin(),
            zm: vec![C64::ZERO; hull.spans_z().len() * (q + 1)],
            xm: Vec::with_capacity(p),
            zm_raw: Vec::with_capacity(q + 1),
            p,
            q,
        }
    }

    /// Free-wave amplitude of the tilted source distribution at λ, for the wave
    /// system with transverse-wavenumber sign `ky_sign` (±1 ⇒ the ±θ system).
    fn eval(&mut self, lambda: f64, ky_sign: f64) -> C64 {
        let kx = self.nu * lambda;
        let ky = self.nu * lambda * (lambda * lambda - 1.0).max(0.0).sqrt();
        let kappa = C64::new(
            self.nu * lambda * lambda * self.cos_phi,
            ky_sign * ky * self.sin_phi,
        );
        if !self.fill_zm(kappa) {
            return C64::ZERO;
        }
        self.accumulate(kx)
    }

    fn fill_zm(&mut self, kappa: C64) -> bool {
        let hull = self.hull;
        let q = self.q;
        let mut any = false;
        for (t, sz) in hull.spans_z().iter().enumerate() {
            let decay = kappa.scale(-sz.start).exp();
            if decay == C64::ZERO {
                for b in 0..=q {
                    self.zm[t * (q + 1) + b] = C64::ZERO;
                }
                continue;
            }
            any = true;
            exp_moments_complex(kappa, sz.len, q, &mut self.zm_raw);
            for b in 0..=q {
                self.zm[t * (q + 1) + b] = decay * self.zm_raw[b];
            }
        }
        any
    }

    fn accumulate(&mut self, kx: f64) -> C64 {
        let hull = self.hull;
        let (p, q) = (self.p, self.q);
        let nsz = hull.spans_z().len();
        let x_center = hull.x_center();
        let coeff = hull.fx_coeff();
        let mut f = C64::ZERO;
        for (s, sx) in hull.spans_x().iter().enumerate() {
            osc_moments(kx, sx.len, p - 1, &mut self.xm);
            let phase = C64::cis(kx * (sx.start - x_center));
            let mut span_sum = C64::ZERO;
            for (a, &xma) in self.xm.iter().enumerate() {
                let mut g_a = C64::ZERO;
                for t in 0..nsz {
                    let base = ((s * nsz + t) * p + a) * (q + 1);
                    let zrow = &self.zm[t * (q + 1)..(t + 1) * (q + 1)];
                    let crow = &coeff[base..base + q + 1];
                    for b in 0..=q {
                        g_a = g_a + zrow[b].scale(crow[b]);
                    }
                }
                span_sum = span_sum + xma * g_a;
            }
            f = f + phase * span_sum;
        }
        f
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hulls::wigley;

    fn synthetic_outer_with_limits(limits: OuterLimits) -> WaveResistance {
        let params = OuterParams {
            nu: 1.0,
            x_half: 0.0,
            y_half: 0.0,
            t_max: 0.0,
        };
        run_outer_with_limits(
            &params,
            &WaveOptions {
                rel_tol: 1e-5,
                max_refinements: 1,
            },
            1.0,
            |_| 1.0,
            limits,
        )
        .0
    }

    #[test]
    fn lambda_cap_is_reported_as_not_converged() {
        let wave = synthetic_outer_with_limits(OuterLimits {
            lambda_hard_cap: 1.5,
            max_evals_per_pass: usize::MAX,
        });
        assert_eq!(wave.method, WaveMethod::GeneralMarcher);
        assert_eq!(wave.outcome, WaveOutcome::TailCap);
    }

    #[test]
    fn evaluation_cap_is_reported_as_not_converged() {
        let wave = synthetic_outer_with_limits(OuterLimits {
            lambda_hard_cap: f64::INFINITY,
            max_evals_per_pass: 16,
        });
        assert_eq!(wave.method, WaveMethod::GeneralMarcher);
        assert_eq!(wave.outcome, WaveOutcome::EvalCap);
    }

    #[test]
    fn staggered_fleet_resistance_is_mirror_symmetric() {
        // A staggered pair and its mirror image about y = 0 are the same
        // physical system; the resistance must not change. (Regression: the
        // integrand once used only the +θ wave system, which broke this.)
        let hull = wigley(8.0, 0.8, 0.5).unwrap();
        let cond = crate::Conditions::seawater(2.5);
        let a = [
            (&hull, Placement { x: 0.0, y: 1.4 }),
            (&hull, Placement { x: 1.7, y: -1.4 }),
        ];
        let b = [
            (&hull, Placement { x: 0.0, y: -1.4 }),
            (&hull, Placement { x: 1.7, y: 1.4 }),
        ];
        let ra = multihull_wave_resistance(&a, &cond).unwrap().resistance;
        let rb = multihull_wave_resistance(&b, &cond).unwrap().resistance;
        assert!(
            (ra - rb).abs() <= 1e-9 * ra,
            "staggered {ra} vs mirrored {rb}"
        );
    }

    #[test]
    fn inner_integral_is_translation_invariant_in_modulus() {
        let l = 10.0;
        let hull = wigley(l, 1.0, 0.625).unwrap();
        // Same hull, shifted +37 m in x.
        let s = hull.surface();
        let shifted = crate::bspline::BSplineSurface::new(
            s.degree_x(),
            s.degree_z(),
            s.knots_x().iter().map(|k| k + 37.0).collect(),
            s.knots_z().to_vec(),
            s.control().to_vec(),
        )
        .unwrap();
        let hull2 = Hull::new(shifted).unwrap();
        let nu = 9.80665 / (3.0f64 * 3.0);
        let mut i1 = InnerIntegral::new(&hull, nu);
        let mut i2 = InnerIntegral::new(&hull2, nu);
        for lambda in [1.0, 1.4, 2.5, 6.0] {
            let a = i1.eval(lambda).abs_sq();
            let b = i2.eval(lambda).abs_sq();
            assert!(
                (a - b).abs() <= 1e-9 * a.abs().max(1e-300),
                "lambda={lambda}: {a} vs {b}"
            );
        }
    }
}

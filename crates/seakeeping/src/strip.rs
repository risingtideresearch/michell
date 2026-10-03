//! **Strip theory**: heave and pitch of a sectional hull in regular waves at
//! forward speed, after Salvesen, Tuck & Faltinsen (1970).
//!
//! Each station's section is solved in two dimensions at the **encounter**
//! frequency ([`crate::section2d`]) for its added mass `a(x)`, damping
//! `b(x)` and heave radiation potential; the hull's coefficients are their
//! integrals along the length with the forward-speed corrections, and the
//! end terms of a **transom** (the aft station's section, when it has one):
//!
//! ```text
//! A33 = ∫a − (U/ω²) b_A                      B33 = ∫b + U a_A
//! A35 = −∫xa − (U/ω²) B⁰ + (U/ω²) x_A b_A    B35 = −∫xb + U A⁰ − U x_A a_A
//! A53 = −∫xa + (U/ω²) B⁰ + (U/ω²) x_A b_A    B53 = −∫xb − U A⁰ − U x_A a_A
//! A55 = ∫x²a + (U²/ω²) A⁰ − (U/ω²) x_A² b_A + (U²/ω²) x_A a_A
//! B55 = ∫x²b + (U²/ω²) B⁰ + U x_A² a_A + (U²/ω²) x_A b_A
//! ```
//!
//! (`A⁰ = ∫a`, `B⁰ = ∫b`; `x` from the pivot, the centre of gravity), and the
//! exciting force is the Froude–Krylov part ([`crate::froude_krylov`]) plus
//! the diffraction part from the Haskind form of each section's radiation
//! potential, with its own speed terms:
//!
//! ```text
//! F3_D = ∫h + (U/iω) h_A,     F5_D = −∫xh − (U/iω) ∫h − (U/iω) x_A h_A
//! ```
//!
//! These are STF's own expressions, in their convention — time `e^{iωt}`
//! with `ω` the encounter frequency, pitch positive **bow down** — which is
//! how they are assembled here; results are converted to this crate's
//! convention (`e^{−iωt}`, pitch **bow up**, see the crate docs) on the way
//! out. Following STF, the diffraction potential is taken at the encounter
//! frequency: `h(x) = ρ ω₀ ω_e ∫ ψ e^{kz} n_z ds` (per unit wave amplitude;
//! with the transverse factors of oblique seas), the classical
//! relative-motion form `−e^{−kT}(ω₀ω_e a + iω₀ b)` in the limit of a
//! shallow section.
//!
//! Strip theory is a slender-body, moderate-speed approximation: trust it
//! to Fn ≈ 0.4 and for waves not much shorter than the hull; the section
//! solver's irregular frequencies (see [`crate::section2d`]) set a further
//! upper limit on the encounter frequency.

use crate::froude_krylov::froude_krylov;
use crate::restoring::restoring;
use crate::section2d::{HeaveSolution, LateralSolution, Section};
use hullgeom::parallel::map_indexed;
use hullgeom::quadrature::gauss_legendre;
use hullgeom::{Error, Placement, Result, SectionalHull, C64};
use std::f64::consts::PI;

/// The ship's mass properties, about its centre of gravity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassProperties {
    /// Mass [kg].
    pub mass: f64,
    /// Longitudinal centre of gravity [m], in the hull's x.
    pub lcg: f64,
    /// Pitch radius of gyration about the centre of gravity [m] (commonly
    /// ≈ 0.25 L).
    pub radius_of_gyration: f64,
    /// Height of the centre of gravity above the centre of buoyancy [m]:
    /// subtracts `ρg∇·BG` from the pitch and roll stiffnesses, and places
    /// the roll axis (through G) above the waterline.
    pub bg: f64,
    /// Roll radius of gyration about the centre of gravity [m] (commonly
    /// 0.35–0.40 of the beam for a monohull; a multihull's is set by its
    /// hull spacing).
    pub roll_radius_of_gyration: f64,
    /// Yaw radius of gyration about the centre of gravity [m] (commonly
    /// close to the pitch one).
    pub yaw_radius_of_gyration: f64,
}

impl MassProperties {
    /// A hull floating freely at its cut attitude: mass `ρ∇`, centre of
    /// gravity on the centre of buoyancy, pitch radius of gyration `k_yy`
    /// (and yaw the same), roll radius of gyration 0.35 of the waterline
    /// beam.
    pub fn floating(hull: &SectionalHull, density: f64, radius_of_gyration: f64) -> Self {
        let (a, b) = hull.x_range();
        let beam = (0..=200)
            .map(|i| 2.0 * hull.waterline_half_beam(a + (b - a) * i as f64 / 200.0))
            .fold(0.0f64, f64::max);
        MassProperties {
            mass: density * hull.displaced_volume(),
            lcg: hull.lcb_x(),
            radius_of_gyration,
            bg: 0.0,
            roll_radius_of_gyration: 0.35 * beam,
            yaw_radius_of_gyration: radius_of_gyration,
        }
    }
}

/// Sea state and speed of one regular-wave case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wave {
    /// Absolute (earth-frame) wave frequency `ω₀` [rad/s].
    pub omega: f64,
    /// Heading `β` [rad]: `π` head seas, `0` following.
    pub heading: f64,
    /// Ship speed `U` [m/s].
    pub speed: f64,
}

/// Solver settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StripOptions {
    /// Panels on each section's half-girth.
    pub panels: usize,
    pub density: f64,
    pub gravity: f64,
    /// Roll damping added to the potential flow's, as a fraction of
    /// critical (`B₄₄ += 2ζ√(C₄₄(I₄₄ + A₄₄))`). Potential flow alone
    /// leaves a monohull's roll almost undamped — the real damping is
    /// viscous (friction, eddies off the bilges, keel lift); a few percent
    /// stands in for it. 0 by default.
    pub roll_damping: f64,
}

impl Default for StripOptions {
    fn default() -> Self {
        StripOptions {
            panels: 20,
            density: 1025.0,
            gravity: hullgeom::STANDARD_GRAVITY,
            roll_damping: 0.0,
        }
    }
}

/// The heave–pitch system at one encounter frequency, in this crate's
/// convention: pitch bow up about the centre of gravity, forces and motions
/// as `e^{−iω_e t}` amplitudes per unit wave amplitude, phases relative to a
/// crest at the centre of gravity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coefficients {
    /// Added mass `[[A33, A35], [A53, A55]]`.
    pub added_mass: [[f64; 2]; 2],
    /// Damping `[[B33, B35], [B53, B55]]`.
    pub damping: [[f64; 2]; 2],
    /// Restoring `[[C33, C35], [C53, C55]]`.
    pub restoring: [[f64; 2]; 2],
    /// Mass `[M, M k_yy²]`.
    pub mass: [f64; 2],
    /// Exciting force and moment: Froude–Krylov.
    pub froude_krylov: [C64; 2],
    /// Exciting force and moment: diffraction.
    pub diffraction: [C64; 2],
    /// The platform's full five-mode system — sway, heave, roll, pitch,
    /// yaw about G, in that order, roll right-handed about +x (port side
    /// up), pitch here **bow down** and yaw bow to port (Salvesen–Tuck–
    /// Faltinsen's axes) — as added mass, damping and restoring.
    pub full: FullCoefficients,
}

/// The five-mode added mass, damping and restoring matrices (see
/// [`Coefficients::full`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FullCoefficients {
    pub added_mass: [[f64; 5]; 5],
    pub damping: [[f64; 5]; 5],
    pub restoring: [[f64; 5]; 5],
    pub mass: [f64; 5],
}

/// The response to one regular wave.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Response {
    pub wave: Wave,
    /// Encounter frequency `ω_e` [rad/s].
    pub omega_e: f64,
    /// Wavenumber `k` [rad/m].
    pub k: f64,
    /// Heave per unit wave amplitude (complex; up).
    pub heave: C64,
    /// Pitch per unit wave amplitude [rad/m] (complex; bow up).
    pub pitch: C64,
    /// Sway of G per unit wave amplitude (complex; to port, +y).
    pub sway: C64,
    /// Roll per unit wave amplitude [rad/m] (complex; right-handed about
    /// +x, the port side rising).
    pub roll: C64,
    /// Yaw per unit wave amplitude [rad/m] (complex; bow to port).
    pub yaw: C64,
    pub coefficients: Coefficients,
}

impl Response {
    /// Heave RAO `|η₃|/ζ_a`.
    pub fn heave_rao(&self) -> f64 {
        self.heave.abs()
    }

    /// Pitch RAO in the usual non-dimensional form `|η₅|/(k ζ_a)`.
    pub fn pitch_rao(&self) -> f64 {
        self.pitch.abs() / self.k
    }

    /// Vertical motion (up) at longitudinal position `x` (hull x) per unit
    /// wave amplitude: `η₃ + (x − x_G) η₅`.
    pub fn vertical_motion(&self, x: f64, lcg: f64) -> C64 {
        self.heave + self.pitch.scale(x - lcg)
    }

    /// Sway RAO `|η₂|/ζ_a`.
    pub fn sway_rao(&self) -> f64 {
        self.sway.abs()
    }

    /// Roll RAO in the usual non-dimensional form `|η₄|/(k ζ_a)`.
    pub fn roll_rao(&self) -> f64 {
        self.roll.abs() / self.k
    }

    /// Yaw RAO `|η₆|/(k ζ_a)`.
    pub fn yaw_rao(&self) -> f64 {
        self.yaw.abs() / self.k
    }
}

/// One station's 2-D solution (or none, for a dry or degenerate station).
struct Strip {
    x: f64,
    sol: Option<HeaveSolution>,
    /// Sway and roll, when the platform's lateral motions are wanted.
    lat: Option<LateralSolution>,
}

impl Strip {
    fn a(&self) -> f64 {
        self.sol.as_ref().map_or(0.0, |s| s.added_mass)
    }
    fn b(&self) -> f64 {
        self.sol.as_ref().map_or(0.0, |s| s.damping)
    }
}

/// Trapezoidal `∫ f dx` over the stations.
fn trapz(strips: &[Strip], f: impl Fn(&Strip) -> C64) -> C64 {
    strips.windows(2).fold(C64::ZERO, |acc, w| {
        acc + (f(&w[0]) + f(&w[1])).scale(0.5 * (w[1].x - w[0].x))
    })
}

/// The heave–pitch response of `hull` with mass properties `mass` to the
/// regular wave `wave`.
pub fn response(
    hull: &SectionalHull,
    mass: &MassProperties,
    wave: &Wave,
    opts: &StripOptions,
) -> Result<Response> {
    response_fleet(&[(hull, Placement::default())], mass, wave, opts)
}

/// The heave–pitch response of a rigid platform of several hulls (a
/// catamaran, trimaran or proa) — each placed in the platform frame, `mass`
/// the whole platform's with its centre of gravity in platform x. The hulls
/// act independently, coupled only through the platform: strip theory has
/// no hull-to-hull wave interaction, which is fair for well-spaced hulls at
/// speed and least trustworthy near the gap's own resonances. In oblique
/// seas each hull feels the wave at its own transverse position.
pub fn response_fleet(
    members: &[(&SectionalHull, Placement)],
    mass: &MassProperties,
    wave: &Wave,
    opts: &StripOptions,
) -> Result<Response> {
    solve(members, mass, wave, opts).map(|(r, _)| r)
}

/// Mean added resistance in a regular wave, with the response it comes
/// from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AddedResistance {
    pub response: Response,
    /// Mean added resistance per unit wave amplitude squared,
    /// `R_aw/ζ_a²` [N/m²].
    pub per_amplitude_sq: f64,
}

impl AddedResistance {
    /// The usual non-dimensional form `R_aw / (ρ g ζ_a² B²/L)`.
    pub fn coefficient(&self, density: f64, gravity: f64, beam: f64, length: f64) -> f64 {
        self.per_amplitude_sq / (density * gravity * beam * beam / length)
    }
}

/// Mean added resistance by the **radiated-energy** method of Gerritsma &
/// Beukelman (1972): the work the ship does radiating waves, through each
/// section's damping, as it moves relative to the local water,
///
/// ```text
/// R_aw = −(k cos β / 2ω_e) ∫ b*(x) |V(x)|² dx,     b* = b − U da/dx,
/// V = −iω_e [η₃ + (x − x_G) η₅ − ζ*(x)]
/// ```
///
/// with `ζ*` the incident wave reduced by the section-mean Smith factor
/// `Z(x; k)/Z(x; 0)` (the mean of `e^{−kz}` over the section's area).
/// Radiated energy vanishes with the motions, so in waves much shorter than
/// the hull — where the real added resistance is the bow's reflection — the
/// method falls to zero; it has no short-wave correction.
pub fn added_resistance(
    hull: &SectionalHull,
    mass: &MassProperties,
    wave: &Wave,
    opts: &StripOptions,
) -> Result<AddedResistance> {
    added_resistance_fleet(&[(hull, Placement::default())], mass, wave, opts)
}

/// [`added_resistance`] for a platform of several hulls (see
/// [`response_fleet`]): each hull's radiated energy, summed.
pub fn added_resistance_fleet(
    members: &[(&SectionalHull, Placement)],
    mass: &MassProperties,
    wave: &Wave,
    opts: &StripOptions,
) -> Result<AddedResistance> {
    let (response, fleet) = solve(members, mass, wave, opts)?;
    let per_amplitude_sq = radiated_energy(members, &fleet, &response, wave);
    Ok(AddedResistance {
        response,
        per_amplitude_sq,
    })
}

/// Both added-resistance estimates from one strip solution:
/// `(response, Gerritsma–Beukelman, Maruo)`, the last two per unit wave
/// amplitude squared [N/m²]. They bracket the tank on Journée's Wigley
/// hulls from different sides (see the crate docs), so report both.
pub fn added_resistance_both(
    members: &[(&SectionalHull, Placement)],
    mass: &MassProperties,
    wave: &Wave,
    opts: &StripOptions,
) -> Result<(Response, f64, f64)> {
    let (response, fleet) = solve(members, mass, wave, opts)?;
    let gb = radiated_energy(members, &fleet, &response, wave);
    let maruo = far_field(members, &fleet, &response, wave, opts);
    Ok((response, gb, maruo))
}

/// Gerritsma–Beukelman's integral (see [`added_resistance`]).
fn radiated_energy(
    members: &[(&SectionalHull, Placement)],
    fleet: &[Vec<Strip>],
    response: &Response,
    wave: &Wave,
) -> f64 {
    let (k, we, u) = (response.k, response.omega_e, wave.speed);
    let (cb, sb) = (wave.heading.cos(), wave.heading.sin());
    let mut integral = 0.0;
    for ((hull, pl), strips) in members.iter().zip(fleet) {
        let smith: Vec<f64> = hull
            .section_nodes()
            .iter()
            .map(|n| {
                let z0 = n.depth_integral(0.0);
                if z0 > 0.0 {
                    n.depth_integral(k) / z0
                } else {
                    0.0
                }
            })
            .collect();
        let n = strips.len();
        let dadx = |i: usize| -> f64 {
            let (lo, hi) = (i.saturating_sub(1), (i + 1).min(n - 1));
            if hi == lo {
                return 0.0;
            }
            (strips[hi].a() - strips[lo].a()) / (strips[hi].x - strips[lo].x)
        };
        let side = C64::cis(k * sb * pl.y);
        let integrand: Vec<f64> = (0..n)
            .map(|i| {
                let s = &strips[i];
                let zeta =
                    (C64::cis(k * cb * s.x) * side).scale(smith.get(i).copied().unwrap_or(0.0));
                let rel = response.heave + response.pitch.scale(s.x) - zeta;
                let b_star = s.b() - u * dadx(i);
                b_star * we * we * rel.abs_sq()
            })
            .collect();
        integral += (0..n.saturating_sub(1))
            .map(|i| 0.5 * (integrand[i] + integrand[i + 1]) * (strips[i + 1].x - strips[i].x))
            .sum::<f64>();
    }
    -k * cb / (2.0 * we) * integral
}

/// Mean added resistance by the **far-field** method of Maruo (1960): the
/// longitudinal momentum the ship's own waves — radiated by its motions and
/// scattered off it — carry away, which in the frame moving with the ship
/// is (as stated by Liu, Liang & Chen, IWWWFB 2025, after Maruo)
///
/// ```text
/// R_aw = ρ/8π { ∫_{−π/2}^{−α₀} + ∫_{α₀}^{π/2} + ∫_{π/2}^{3π/2} } |H(k₁,θ)|² k₁(k₁cosθ − k cosβ)/√(1−4τcosθ) dθ
///      + ρ/8π ∫_{α₀}^{2π−α₀} |H(k₂,θ)|² k₂(k₂cosθ − k cosβ)/√(1−4τcosθ) dθ
/// k₁,₂(θ) = (K₀/2)(1 − 2τcosθ ± √(1−4τcosθ))/cos²θ,   τ = ω_e U/g,   K₀ = g/U²
/// ```
///
/// (`α₀ = acos(1/4τ)` beyond `τ = 1/4`, else 0), which at zero speed is
/// `ρk²/8π ∫ |H(k,θ)|²(cosθ − cosβ) dθ`. The Kochin function is the hull's
/// outward source flux summed against each far-field wave,
/// `H(k,θ) = ∬ σ e^{kz} e^{−ik(x cosθ + y sinθ)} dS` ([`kochin`]), from the
/// strip solution: each station's heave sources at its relative vertical
/// velocity `V(x) = −iω_e(η₃ + xη₅) − Uη₅` plus its diffraction sources.
///
/// Waves shorter than twice the station spacing are left out (see the
/// source).
///
/// Unlike [`added_resistance`] (Gerritsma–Beukelman), which sums each
/// strip's radiated energy, the waves of all the sections interfere here
/// before their momentum is counted, and the scattered incident wave is in
/// the balance too. Head and following seas only: the sections carry only
/// the part of the diffraction problem symmetric about their centreplanes.
pub fn added_resistance_maruo(
    hull: &SectionalHull,
    mass: &MassProperties,
    wave: &Wave,
    opts: &StripOptions,
) -> Result<AddedResistance> {
    added_resistance_maruo_fleet(&[(hull, Placement::default())], mass, wave, opts)
}

/// [`added_resistance_maruo`] for a platform of several placed hulls, whose
/// waves interfere in the one Kochin function.
pub fn added_resistance_maruo_fleet(
    members: &[(&SectionalHull, Placement)],
    mass: &MassProperties,
    wave: &Wave,
    opts: &StripOptions,
) -> Result<AddedResistance> {
    let (response, fleet) = solve(members, mass, wave, opts)?;
    let per_amplitude_sq = far_field(members, &fleet, &response, wave, opts);
    Ok(AddedResistance {
        response,
        per_amplitude_sq,
    })
}

/// Maruo's integral (see [`added_resistance_maruo`]).
fn far_field(
    members: &[(&SectionalHull, Placement)],
    fleet: &[Vec<Strip>],
    response: &Response,
    wave: &Wave,
    opts: &StripOptions,
) -> f64 {
    let (k, we, u, g) = (response.k, response.omega_e, wave.speed, opts.gravity);
    let cb = wave.heading.cos();
    let h = |kk: f64, theta: f64| kochin(members, fleet, response, wave, kk, theta).abs_sq();
    // The Kochin function is assembled from stations; waves shorter than
    // twice their spacing only alias it (the hull's own, smooth sources
    // radiate next to nothing there), so the integrand stops at that
    // Nyquist wavenumber — which the k₁ system's short divergent waves
    // would otherwise reach near θ = ±π/2.
    let spacing = fleet
        .iter()
        .flat_map(|st| st.windows(2).map(|w| w[1].x - w[0].x))
        .fold(0.0f64, f64::max);
    let k_max = PI / spacing.max(1e-9);
    let integral = if u == 0.0 {
        // One wave system, k₂ = k all round.
        integrate(0.0, 2.0 * PI, |t| h(k, t) * k * (k * t.cos() - k * cb))
    } else {
        let tau = we * u / g;
        let k0 = g / (u * u);
        let a0 = if tau > 0.25 {
            (1.0 / (4.0 * tau)).acos()
        } else {
            0.0
        };
        let root = |c: f64| (1.0 - 4.0 * tau * c).max(0.0).sqrt();
        let k1 = |c: f64| 0.5 * k0 * (1.0 - 2.0 * tau * c + root(c)) / (c * c);
        // The conjugate form, free of cancellation as cos θ → 0.
        let k2 = |c: f64| 2.0 * k0 * tau * tau / (1.0 - 2.0 * tau * c + root(c));
        let f = |kj: f64, t: f64| {
            let c = t.cos();
            let r = root(c);
            if !(kj.is_finite()) || kj > k_max || r == 0.0 {
                return 0.0;
            }
            h(kj, t) * kj * (kj * c - k * cb) / r
        };
        let branch1 = |t: f64| f(k1(t.cos()), t);
        let branch2 = |t: f64| f(k2(t.cos()), t);
        integrate(-0.5 * PI, -a0, branch1)
            + integrate(a0, 0.5 * PI, branch1)
            + integrate(0.5 * PI, 1.5 * PI, branch1)
            + integrate(a0, 2.0 * PI - a0, branch2)
    };
    opts.density / (8.0 * PI) * integral
}

/// The Kochin function `H(k,θ) = ∬ σ e^{kz} e^{−ik(x cosθ + y sinθ)} dS` of
/// the platform's disturbance — each station's heave sources at its
/// relative vertical velocity, plus its diffraction sources — per unit
/// incident wave amplitude, `x` from the centre of gravity.
fn kochin(
    members: &[(&SectionalHull, Placement)],
    fleet: &[Vec<Strip>],
    response: &Response,
    wave: &Wave,
    kk: f64,
    theta: f64,
) -> C64 {
    let (st, ct) = theta.sin_cos();
    let (k, we, u) = (response.k, response.omega_e, wave.speed);
    let (cb, sb) = (wave.heading.cos(), wave.heading.sin());
    let mut total = C64::ZERO;
    for ((_, pl), strips) in members.iter().zip(fleet) {
        let side = C64::cis(-kk * pl.y * st);
        let incident_side = C64::cis(k * sb * pl.y);
        let term = |s: &Strip| -> C64 {
            s.sol.as_ref().map_or(C64::ZERO, |sol| {
                let v = C64::new(0.0, -we) * (response.heave + response.pitch.scale(s.x))
                    - response.pitch.scale(u);
                let d = C64::cis(k * cb * s.x) * incident_side;
                sol.source_spectrum(kk, st, v, d) * C64::cis(-kk * s.x * ct)
            })
        };
        total = total + trapz(strips, term) * side;
    }
    total
}

/// `∫_a^b f(θ) dθ` for integrands with inverse-square-root singularities
/// at either end: `θ = a + (b − a)(1 − cos πt)/2` makes them bounded, then
/// Gauss–Legendre panels in `t`.
fn integrate(a: f64, b: f64, f: impl Fn(f64) -> f64) -> f64 {
    if !(b > a) {
        return 0.0;
    }
    const PANELS: usize = 96;
    let (gx, gw) = gauss_legendre(8);
    let mut total = 0.0;
    for p in 0..PANELS {
        let (t0, t1) = (p as f64 / PANELS as f64, (p + 1) as f64 / PANELS as f64);
        for (&x, &w) in gx.iter().zip(&gw) {
            let t = 0.5 * (t0 + t1) + 0.5 * (t1 - t0) * x;
            let theta = a + (b - a) * 0.5 * (1.0 - (PI * t).cos());
            let jac = (b - a) * 0.5 * PI * (PI * t).sin();
            total += 0.5 * (t1 - t0) * w * jac * f(theta);
        }
    }
    total
}

/// A 5×5 complex matrix and helpers for the five-mode system.
type M5 = [[C64; 5]; 5];

fn zero5() -> M5 {
    [[C64::ZERO; 5]; 5]
}

/// The section kinematics: platform motions `(sway, heave, roll, pitch
/// bow-down, yaw)` about G to the section's local sway, heave and roll
/// (about its waterline centre), at `X` from G, hull offset `y`, G `zg`
/// above the waterline; and its x-derivative.
fn kinematics(x: f64, y: f64, zg: f64) -> ([[f64; 5]; 3], [[f64; 5]; 3]) {
    let t = [
        [1.0, 0.0, zg, 0.0, x],
        [0.0, 1.0, y, -x, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
    ];
    let tp = [
        [0.0, 0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, -1.0, 0.0],
        [0.0; 5],
    ];
    (t, tp)
}

/// `Lᵀ c R` for 3×5 real `L`, `R` given as complex combinations
/// `(l0 T + l1 T')` and `(r0 T + r1 T')`.
fn sandwich(
    t: &[[f64; 5]; 3],
    tp: &[[f64; 5]; 3],
    l: (C64, C64),
    c: &[[C64; 3]; 3],
    r: (C64, C64),
) -> M5 {
    let mut left = [[C64::ZERO; 5]; 3];
    let mut right = [[C64::ZERO; 5]; 3];
    for i in 0..3 {
        for j in 0..5 {
            left[i][j] = l.0.scale(t[i][j]) + l.1.scale(tp[i][j]);
            right[i][j] = r.0.scale(t[i][j]) + r.1.scale(tp[i][j]);
        }
    }
    let mut out = zero5();
    for a in 0..5 {
        for b in 0..5 {
            let mut acc = C64::ZERO;
            for i in 0..3 {
                for j in 0..3 {
                    acc = acc + left[i][a] * c[i][j] * right[j][b];
                }
            }
            out[a][b] = acc;
        }
    }
    out
}

/// A station's complex added mass `c = a − ib/ω` (the `e^{iωt}` form),
/// local `(sway, heave, roll)`.
fn station_c(s: &Strip, we: f64) -> [[C64; 3]; 3] {
    let cx = |a: f64, b: f64| C64::new(a, -b / we);
    let mut c = [[C64::ZERO; 3]; 3];
    if let Some(h) = &s.sol {
        c[1][1] = cx(h.added_mass, h.damping);
    }
    if let Some(l) = &s.lat {
        let (a, b) = (l.added_mass, l.damping);
        c[0][0] = cx(a[0][0], b[0][0]);
        c[0][2] = cx(a[0][1], b[0][1]);
        c[2][0] = cx(a[1][0], b[1][0]);
        c[2][2] = cx(a[1][1], b[1][1]);
    }
    c
}

/// Solve the platform's five-mode system — sway, heave, roll, pitch, yaw
/// about G — by strip theory; also returns each member's stations (x from
/// the centre of gravity) and their 2-D solutions.
///
/// Each station's local motion (sway, heave and roll about its waterline
/// centre) follows from the platform's through a kinematic map `T(x)`
/// (the hull's offset and G's height enter here, so multihulls and their
/// roll need nothing special), and its hydrodynamic force is the
/// forward-speed operator `D = iω − U∂ₓ` wrapped round its complex added
/// mass, `f = −D[c D(Tη)]`. Integrated along the hull (by parts, the aft
/// station's section left as the transom's end term) this gives
///
/// ```text
/// K = ∫ (iωT + UT′)ᵀ c (iωT − UT′) dx + U T_Aᵀ c_A (iωT_A − UT′_A),   F_rad = −Kη
/// ```
///
/// — every coefficient of Salvesen, Tuck & Faltinsen (1970), sway, roll and
/// yaw with heave and pitch, their transom terms included. The exciting
/// force is Froude–Krylov (closed form in the vertical, per station in the
/// lateral) plus Haskind diffraction with STF's speed terms,
/// `∫Tᵀh dx + (U/iω)(∫T′ᵀh dx + T_Aᵀh_A)`. Lateral sections are solved only
/// when the waves are oblique or the platform asymmetric; otherwise the
/// lateral motions are unforced and uncoupled, and zero.
fn solve(
    members: &[(&SectionalHull, Placement)],
    mass: &MassProperties,
    wave: &Wave,
    opts: &StripOptions,
) -> Result<(Response, Vec<Vec<Strip>>)> {
    if members.is_empty() {
        return Err(Error::InvalidGeometry("empty fleet".into()));
    }
    let (rho, g) = (opts.density, opts.gravity);
    let w0 = wave.omega;
    let k = w0 * w0 / g;
    let (cb, sb) = (wave.heading.cos(), wave.heading.sin());
    let we = w0 - k * wave.speed * cb;
    if !(we > 1e-6 * w0) {
        return Err(Error::InvalidInput(format!(
            "encounter frequency {we:.4} rad/s is not positive: strip theory here \
             needs the waves met from ahead of the stern (ω₀ = {w0}, U cos β = {})",
            wave.speed * cb
        )));
    }
    let u = wave.speed;
    let xg = mass.lcg;
    // The platform's waterplane: its transverse first moment says whether it
    // is symmetric (roll uncoupled from heave and pitch).
    let area: f64 = members.iter().map(|(h, _)| h.waterplane_area()).sum();
    let first_y: f64 = members
        .iter()
        .map(|(h, pl)| h.waterplane_area() * pl.y)
        .sum();
    let length = members
        .iter()
        .map(|(h, _)| h.length())
        .fold(0.0f64, f64::max);
    let lateral = sb.abs() > 1e-9 || first_y.abs() > 1e-9 * area * length.max(1.0);
    let volume: f64 = members.iter().map(|(h, _)| h.displaced_volume()).sum();
    let zg = mass.bg - crate::restoring::buoyancy_depth(members);
    // Stations and their 2-D solutions at the encounter frequency, per hull.
    let fleet: Vec<Vec<Strip>> = members
        .iter()
        .map(|(hull, pl)| {
            let curves: Vec<(f64, Vec<(f64, f64)>)> =
                hull.curves().map(|(x, c)| (x, c.to_vec())).collect();
            map_indexed(
                curves.len(),
                || (),
                |_, i| {
                    let (x, c) = &curves[i];
                    let sec = Section::from_curve(c, opts.panels)
                        .filter(|s| s.half_beam() > 0.0 && s.draft() > 0.0);
                    let sol = sec
                        .as_ref()
                        .and_then(|s| s.heave_with_diffraction(we, k, wave.heading, g, rho));
                    let lat = if lateral {
                        sec.as_ref()
                            .and_then(|s| s.lateral(we, g, rho, Some((k, wave.heading))))
                    } else {
                        None
                    };
                    Strip {
                        x: x + pl.x - xg,
                        sol,
                        lat,
                    }
                },
            )
        })
        .collect();
    let iw = C64::new(0.0, we);
    let uc = C64::new(u, 0.0);
    // Radiation.
    let mut kmat = zero5();
    let add = |acc: &mut M5, m: &M5, w: f64| {
        for a in 0..5 {
            for b in 0..5 {
                acc[a][b] = acc[a][b] + m[a][b].scale(w);
            }
        }
    };
    for ((_, pl), strips) in members.iter().zip(&fleet) {
        let station = |s: &Strip| {
            let (t, tp) = kinematics(s.x, pl.y, zg);
            sandwich(&t, &tp, (iw, uc), &station_c(s, we), (iw, -uc))
        };
        for w in strips.windows(2) {
            let dx = w[1].x - w[0].x;
            add(&mut kmat, &station(&w[0]), 0.5 * dx);
            add(&mut kmat, &station(&w[1]), 0.5 * dx);
        }
        if let Some(a) = strips.first().filter(|s| s.sol.is_some()) {
            let (t, tp) = kinematics(a.x, pl.y, zg);
            add(
                &mut kmat,
                &sandwich(&t, &tp, (uc, C64::ZERO), &station_c(a, we), (iw, -uc)),
                1.0,
            );
        }
    }
    // Restoring about G (STF axes: pitch bow down).
    let mut cmat = [[0.0; 5]; 5];
    let (mut it, mut c35s, mut c55s, mut yc35) = (0.0, 0.0, 0.0, 0.0);
    for (hull, pl) in members {
        let c = restoring(hull, rho, g, xg - pl.x);
        c35s += c.c35;
        c55s += c.c55;
        yc35 += pl.y * c.c35;
        it += hull.transverse_waterplane_inertia() + hull.waterplane_area() * pl.y * pl.y;
    }
    let rg = rho * g;
    cmat[1][1] = rg * area;
    cmat[1][2] = rg * first_y;
    cmat[2][1] = cmat[1][2];
    cmat[1][3] = -c35s;
    cmat[3][1] = -c35s;
    cmat[2][2] = rg * (it - volume * mass.bg);
    cmat[2][3] = -yc35;
    cmat[3][2] = -yc35;
    cmat[3][3] = c55s - rg * volume * mass.bg;
    let m = mass.mass;
    let inertia = [
        m,
        m,
        m * mass.roll_radius_of_gyration.powi(2),
        m * mass.radius_of_gyration.powi(2),
        m * mass.yaw_radius_of_gyration.powi(2),
    ];
    // Extra roll damping as a fraction of critical.
    if opts.roll_damping > 0.0 && cmat[2][2] > 0.0 {
        let a44 = -kmat[2][2].re / (we * we);
        let b = 2.0 * opts.roll_damping * (cmat[2][2] * (inertia[2] + a44).max(0.0)).sqrt();
        kmat[2][2] = kmat[2][2] + iw.scale(b);
    }
    // Excitation (STF convention: conjugates of this crate's amplitudes).
    let mut f_fk = [C64::ZERO; 5];
    let mut f_d = [C64::ZERO; 5];
    for ((hull, pl), strips) in members.iter().zip(&fleet) {
        let side = C64::cis(k * sb * pl.y);
        // Vertical Froude–Krylov: the closed form, about G.
        let fk = froude_krylov(hull, rho, g, k, wave.heading, xg - pl.x);
        let (f3, f5) = ((fk.heave * side).conj(), -(fk.pitch * side).conj());
        f_fk[1] = f_fk[1] + f3;
        f_fk[2] = f_fk[2] + f3.scale(pl.y);
        f_fk[3] = f_fk[3] + f5;
        // Per station: lateral Froude–Krylov, and diffraction (with STF's
        // encounter-frequency Haskind scaling), as local (sway, heave, roll).
        let local = |s: &Strip| -> ([C64; 3], [C64; 3]) {
            let ph = C64::cis(k * cb * s.x) * side;
            let mut fk_l = [C64::ZERO; 3];
            let mut h_l = [C64::ZERO; 3];
            if let Some(sol) = &s.sol {
                h_l[1] = (sol.diffraction(k, wave.heading, g, rho).scale(we / w0) * ph).conj();
            }
            if let Some(lat) = &s.lat {
                let fkl = lat.froude_krylov(k, wave.heading, g, rho);
                let hl = lat.diffraction(k, wave.heading, g, rho);
                fk_l[0] = (fkl[0] * ph).conj();
                fk_l[2] = (fkl[1] * ph).conj();
                h_l[0] = (hl[0].scale(we / w0) * ph).conj();
                h_l[2] = (hl[1].scale(we / w0) * ph).conj();
            }
            (fk_l, h_l)
        };
        let project = |t: &[[f64; 5]; 3], v: &[C64; 3]| -> [C64; 5] {
            let mut out = [C64::ZERO; 5];
            for (a, o) in out.iter_mut().enumerate() {
                for i in 0..3 {
                    *o = *o + v[i].scale(t[i][a]);
                }
            }
            out
        };
        let u_iw = C64::new(0.0, -u / we); // U/(iω)
        for w in strips.windows(2) {
            let dx = w[1].x - w[0].x;
            for s in [&w[0], &w[1]] {
                let (t, tp) = kinematics(s.x, pl.y, zg);
                let (fk_l, h_l) = local(s);
                let (a, b, c) = (project(&t, &fk_l), project(&t, &h_l), project(&tp, &h_l));
                for i in 0..5 {
                    f_fk[i] = f_fk[i] + a[i].scale(0.5 * dx);
                    f_d[i] = f_d[i] + (b[i] + u_iw * c[i]).scale(0.5 * dx);
                }
            }
        }
        if let Some(a) = strips.first().filter(|s| s.sol.is_some()) {
            let (t, _) = kinematics(a.x, pl.y, zg);
            let (_, h_l) = local(a);
            let e = project(&t, &h_l);
            for i in 0..5 {
                f_d[i] = f_d[i] + u_iw * e[i];
            }
        }
    }
    // Solve (−ω²M + K + C) η = F.
    let mut z = vec![C64::ZERO; 25];
    for a in 0..5 {
        for b in 0..5 {
            let mut v = kmat[a][b] + C64::new(cmat[a][b], 0.0);
            if a == b {
                v = v - C64::new(we * we * inertia[a], 0.0);
            }
            z[a * 5 + b] = v;
        }
    }
    let f: Vec<C64> = (0..5).map(|i| f_fk[i] + f_d[i]).collect();
    let eta = crate::linalg::solve(z, f)
        .ok_or_else(|| Error::InvalidInput("singular motion system".into()))?;
    // Coefficients: A = −Re K/ω², B = Im K/ω.
    let mut full = FullCoefficients {
        added_mass: [[0.0; 5]; 5],
        damping: [[0.0; 5]; 5],
        restoring: cmat,
        mass: inertia,
    };
    for a in 0..5 {
        for b in 0..5 {
            full.added_mass[a][b] = -kmat[a][b].re / (we * we);
            full.damping[a][b] = kmat[a][b].im / we;
        }
    }
    // To this crate's convention: conjugate (time), negate pitch (bow up).
    let (am, dm) = (full.added_mass, full.damping);
    let coefficients = Coefficients {
        added_mass: [[am[1][1], -am[1][3]], [-am[3][1], am[3][3]]],
        damping: [[dm[1][1], -dm[1][3]], [-dm[3][1], dm[3][3]]],
        restoring: [[cmat[1][1], -cmat[1][3]], [-cmat[3][1], cmat[3][3]]],
        mass: [inertia[1], inertia[3]],
        froude_krylov: [f_fk[1].conj(), -f_fk[3].conj()],
        diffraction: [f_d[1].conj(), -f_d[3].conj()],
        full,
    };
    let response = Response {
        wave: *wave,
        omega_e: we,
        k,
        heave: eta[1].conj(),
        pitch: -eta[3].conj(),
        sway: eta[0].conj(),
        roll: eta[2].conj(),
        yaw: eta[4].conj(),
        coefficients,
    };
    Ok((response, fleet))
}

#[cfg(test)]
mod tests {
    use super::*;
    use hullgeom::iges::{self, HullPose, Platform, SectionalOptions};
    use std::f64::consts::PI;

    const RHO: f64 = 1000.0;
    const G: f64 = 9.81;

    /// Journée's Wigley models are 3 m long; this is the classic
    /// parabolic Wigley at those proportions (L/B = 10, B/T = 1.6).
    fn wigley(l: f64) -> SectionalHull {
        let surfaces = iges::wigley_surfaces(l, 0.1 * l, 0.0625 * l).unwrap();
        let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
        let opts = SectionalOptions {
            stations: 41,
            ..SectionalOptions::default()
        };
        source
            .situate_sectional(0, 0.0, &HullPose::default(), &Platform::default(), &opts)
            .unwrap()
            .unwrap()
            .hull
    }

    fn opts() -> StripOptions {
        StripOptions {
            panels: 16,
            density: RHO,
            gravity: G,
            roll_damping: 0.0,
        }
    }

    /// A wave many hull lengths long lifts the hull with the water and tilts
    /// it with the slope: heave → 1, pitch → i k cos β (bow up where the
    /// surface rises ahead).
    #[test]
    fn a_long_wave_is_followed() {
        let l = 3.0;
        let hull = wigley(l);
        let mass = MassProperties::floating(&hull, RHO, 0.25 * l);
        for heading in [PI, 0.0] {
            let k = 2.0 * PI / (40.0 * l);
            let wave = Wave {
                omega: (k * G).sqrt(),
                heading,
                speed: 0.0,
            };
            let r = response(&hull, &mass, &wave, &opts()).unwrap();
            let slope = C64::new(0.0, k * heading.cos());
            assert!(
                (r.heave - C64::ONE).abs() < 0.02,
                "β {heading}: heave {:?}",
                r.heave
            );
            assert!(
                (r.pitch - slope).abs() < 0.03 * k,
                "β {heading}: pitch {:?} vs {slope:?}",
                r.pitch
            );
        }
    }

    /// At zero speed the added-mass and damping matrices are symmetric, and
    /// the forward-speed terms are the only source of asymmetry.
    #[test]
    fn zero_speed_coefficients_are_symmetric() {
        let l = 3.0;
        let hull = wigley(l);
        let mass = MassProperties::floating(&hull, RHO, 0.25 * l);
        let wave = Wave {
            omega: (2.0 * PI / l * G).sqrt(),
            heading: PI,
            speed: 0.0,
        };
        let c = response(&hull, &mass, &wave, &opts()).unwrap().coefficients;
        let (a, b) = (c.added_mass, c.damping);
        assert!((a[0][1] - a[1][0]).abs() < 1e-9 * a[0][0].abs() * l);
        assert!((b[0][1] - b[1][0]).abs() < 1e-9 * b[0][0].abs().max(1e-12) * l);
        // A fore-aft symmetric hull about its centre has no coupling at all.
        assert!(a[0][1].abs() < 1e-3 * a[0][0] * l, "A35 {}", a[0][1]);
        let with_speed = response(&hull, &mass, &Wave { speed: 1.0, ..wave }, &opts())
            .unwrap()
            .coefficients;
        assert!(
            (with_speed.damping[0][1] - with_speed.damping[1][0]).abs()
                > 1e-3 * with_speed.damping[0][0] * l
        );
    }

    /// Gerritsma–Beukelman added resistance in head seas at Fn 0.3: positive
    /// throughout, peaking with the motions (λ/L ≈ 1–1.2) and vanishing in
    /// long waves, where the hull follows the water. The peak here,
    /// `R_aw/(ρgζ²B²/L)` ≈ 44 at λ/L = 1 on this parabolic Wigley, is of the
    /// order strip theory gives Journée's Wigley III at Fn 0.3 (≈ 49 at λ/L
    /// 1.05, against a measured ≈ 20 at 1.25; see `crate::validation`). The
    /// magnitude bound is a sanity guard only.
    #[test]
    fn added_resistance_peaks_with_the_motions() {
        let l = 3.0;
        let hull = wigley(l);
        let mass = MassProperties::floating(&hull, RHO, 0.25 * l);
        let u = 0.3 * (G * l).sqrt();
        let mut rows = Vec::new();
        for i in 0..10 {
            let lam = l * (0.8 + 0.2 * i as f64);
            let k = 2.0 * PI / lam;
            let wave = Wave {
                omega: (k * G).sqrt(),
                heading: PI,
                speed: u,
            };
            let r = added_resistance(&hull, &mass, &wave, &opts()).unwrap();
            let sigma = r.coefficient(RHO, G, 0.1 * l, l);
            eprintln!(
                "λ/L {:.2}: heave {:.3} pitch {:.3} σ_aw {:.3}",
                lam / l,
                r.response.heave_rao(),
                r.response.pitch_rao(),
                sigma
            );
            rows.push((lam / l, sigma));
        }
        assert!(rows.iter().all(|r| r.1 > 0.0), "{rows:?}");
        let peak = rows
            .iter()
            .cloned()
            .fold((0.0, 0.0), |a, r| if r.1 > a.1 { r } else { a });
        assert!((1.0..=1.6).contains(&peak.0), "peak at {peak:?}");
        assert!(peak.1 > 1.0 && peak.1 < 80.0, "peak {peak:?}");
    }

    /// Two identical hulls far apart in head seas carry twice the forces
    /// and twice the mass of one: the platform moves exactly like the
    /// single hull. In beam-ish seas the hulls meet the wave at different
    /// phases, so the platform heaves less than a lone hull would.
    #[test]
    fn a_catamaran_of_twins_heaves_like_one_hull() {
        let l = 3.0;
        let hull = wigley(l);
        let one = MassProperties::floating(&hull, RHO, 0.25 * l);
        let two = MassProperties {
            mass: 2.0 * one.mass,
            ..one
        };
        let span = 1.2;
        let cat = [
            (
                &hull,
                Placement {
                    x: 0.0,
                    y: 0.5 * span,
                },
            ),
            (
                &hull,
                Placement {
                    x: 0.0,
                    y: -0.5 * span,
                },
            ),
        ];
        let u = 0.3 * (G * l).sqrt();
        let k = 2.0 * PI / (1.2 * l);
        let head = Wave {
            omega: (k * G).sqrt(),
            heading: PI,
            speed: u,
        };
        let solo = response(&hull, &one, &head, &opts()).unwrap();
        let pair = response_fleet(&cat, &two, &head, &opts()).unwrap();
        assert!(
            (solo.heave - pair.heave).abs() < 1e-9,
            "{:?} vs {:?}",
            solo.heave,
            pair.heave
        );
        assert!((solo.pitch - pair.pitch).abs() < 1e-9 * k);
        // Bow-quartering: the hulls 1.2 m apart see the wave out of phase.
        let oblique = Wave {
            heading: 0.75 * PI,
            ..head
        };
        let solo_ob = response(&hull, &one, &oblique, &opts()).unwrap();
        let pair_ob = response_fleet(&cat, &two, &oblique, &opts()).unwrap();
        assert!(
            pair_ob.heave_rao() < solo_ob.heave_rao(),
            "{} vs {}",
            pair_ob.heave_rao(),
            solo_ob.heave_rao()
        );
        let raw_solo = added_resistance(&hull, &one, &head, &opts())
            .unwrap()
            .per_amplitude_sq;
        let raw_pair = added_resistance_fleet(&cat, &two, &head, &opts())
            .unwrap()
            .per_amplitude_sq;
        assert!((raw_pair - 2.0 * raw_solo).abs() < 1e-9 * raw_solo.abs());
    }

    /// The Kochin normalisation: a hull heaving at unit velocity (no speed,
    /// no waves diffracted) radiates the power `½B` with
    /// `B = (ωρk/4π) ∫ |H(k,θ)|² dθ`. Where the waves are short against the
    /// hull but long against its beam (`kL ≫ 1 ≫ kB`) the sections radiate
    /// independently and sideways, and that is strip theory's `∫ b dx`; on a
    /// hull as wide as the waves are long its 2-D sources also radiate fore
    /// and aft, which strip theory never sees, so the check needs a very
    /// slender hull (L/B = 50).
    #[test]
    fn the_kochin_function_carries_the_strip_damping() {
        let l = 3.0;
        let surfaces = iges::wigley_surfaces(l, 0.02 * l, 0.0125 * l).unwrap();
        let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
        let so = SectionalOptions {
            stations: 61,
            ..SectionalOptions::default()
        };
        let hull = source
            .situate_sectional(0, 0.0, &HullPose::default(), &Platform::default(), &so)
            .unwrap()
            .unwrap()
            .hull;
        let mass = MassProperties::floating(&hull, RHO, 0.25 * l);
        let members = [(&hull, Placement::default())];
        for lam in [0.25 * l, 0.4 * l] {
            let k = 2.0 * PI / lam;
            let wave = Wave {
                omega: (k * G).sqrt(),
                heading: PI,
                speed: 0.0,
            };
            let (resp, fleet) = solve(&members, &mass, &wave, &opts()).unwrap();
            // Unit heave velocity: η₃ = i/ω (e^{−iωt}), no pitch, and the
            // heave sources alone.
            let heave_only = Response {
                heave: C64::new(0.0, 1.0 / resp.omega_e),
                pitch: C64::ZERO,
                ..resp
            };
            let radiating: Vec<Vec<Strip>> = fleet
                .into_iter()
                .map(|st| {
                    st.into_iter()
                        .map(|s| Strip {
                            x: s.x,
                            sol: s.sol.map(|mut sol| {
                                sol.clear_diffraction();
                                sol
                            }),
                            lat: None,
                        })
                        .collect()
                })
                .collect();
            let hsq = |t: f64| kochin(&members, &radiating, &heave_only, &wave, k, t).abs_sq();
            let b3d = resp.omega_e * RHO * k / (4.0 * PI) * integrate(0.0, 2.0 * PI, hsq);
            let b2d = resp.coefficients.damping[0][0];
            assert!(
                (b3d - b2d).abs() < 0.1 * b2d,
                "λ/L {}: B from Kochin {b3d} vs ∫b dx {b2d}",
                lam / l
            );
        }
    }

    /// A beam-sea wave many hull lengths long carries the hull with the
    /// water: heave → 1, sway → 1 (the surface particles' orbit), roll →
    /// the wave slope `k` (for a roll natural frequency well above the
    /// wave's), and a fore-aft symmetric hull does not yaw.
    #[test]
    fn a_long_beam_wave_is_followed() {
        let l = 3.0;
        let hull = wigley(l);
        let mut mass = MassProperties::floating(&hull, RHO, 0.25 * l);
        mass.bg = -0.05; // G below B: stiff in roll
        let k = 2.0 * PI / (60.0 * l);
        let wave = Wave {
            omega: (k * G).sqrt(),
            heading: 0.5 * PI,
            speed: 0.0,
        };
        let r = response(&hull, &mass, &wave, &opts()).unwrap();
        assert!(
            (r.heave_rao() - 1.0).abs() < 0.03,
            "heave {}",
            r.heave_rao()
        );
        assert!((r.sway_rao() - 1.0).abs() < 0.05, "sway {}", r.sway_rao());
        assert!((r.roll_rao() - 1.0).abs() < 0.1, "roll/k {}", r.roll_rao());
        assert!(r.yaw_rao() < 1e-3, "yaw/k {}", r.yaw_rao());
    }

    /// At zero speed the five-mode added mass and damping are symmetric
    /// (reciprocity), sway, roll and yaw included.
    #[test]
    fn zero_speed_five_mode_coefficients_are_reciprocal() {
        let l = 3.0;
        let hull = wigley(l);
        let mass = MassProperties::floating(&hull, RHO, 0.25 * l);
        let k = 2.0 * PI / (1.2 * l);
        let wave = Wave {
            omega: (k * G).sqrt(),
            heading: 0.6 * PI,
            speed: 0.0,
        };
        let f = response(&hull, &mass, &wave, &opts())
            .unwrap()
            .coefficients
            .full;
        let scale = f.added_mass[0][0].abs().max(f.added_mass[1][1].abs());
        for a in 0..5 {
            for b in 0..5 {
                let arm = if a > 1 || b > 1 { l } else { 1.0 };
                assert!(
                    (f.added_mass[a][b] - f.added_mass[b][a]).abs() < 1e-6 * scale * arm * arm,
                    "A{a}{b} {} vs A{b}{a} {}",
                    f.added_mass[a][b],
                    f.added_mass[b][a]
                );
                let dscale = f.damping[a][a].abs().max(f.damping[b][b].abs()).max(1e-9);
                assert!(
                    (f.damping[a][b] - f.damping[b][a]).abs() < 1e-2 * dscale,
                    "B{a}{b} {} vs B{b}{a} {} (diag {dscale})",
                    f.damping[a][b],
                    f.damping[b][a]
                );
            }
        }
        assert!(f.added_mass[0][0] > 0.0 && f.added_mass[2][2] > 0.0 && f.added_mass[4][4] > 0.0);
    }

    /// Mirror images move as mirror images: a proa with its ama to port
    /// rolls and yaws opposite to one with it to starboard (heave and pitch
    /// the same), and a head-sea wave rolls it at all — the asymmetric
    /// platform couples roll to heave. A symmetric catamaran in head seas
    /// neither rolls nor sways.
    #[test]
    fn mirrored_platforms_move_as_mirror_images() {
        let main = wigley(3.0);
        let ama = wigley(1.8);
        let mass = |ms: &[(&SectionalHull, Placement)]| {
            let vol: f64 = ms.iter().map(|(h, _)| h.displaced_volume()).sum();
            MassProperties {
                mass: RHO * vol,
                lcg: 0.0,
                radius_of_gyration: 0.75,
                bg: 0.0,
                roll_radius_of_gyration: 0.6,
                yaw_radius_of_gyration: 0.8,
            }
        };
        let port = [
            (&main, Placement::default()),
            (&ama, Placement { x: 0.0, y: 1.2 }),
        ];
        let star = [
            (&main, Placement::default()),
            (&ama, Placement { x: 0.0, y: -1.2 }),
        ];
        let k = 2.0 * PI / 3.6;
        let head = Wave {
            omega: (k * G).sqrt(),
            heading: PI,
            speed: 1.0,
        };
        let a = response_fleet(&port, &mass(&port), &head, &opts()).unwrap();
        let b = response_fleet(&star, &mass(&star), &head, &opts()).unwrap();
        assert!(
            a.roll.abs() > 1e-3 * k,
            "the proa should roll: {:?}",
            a.roll
        );
        assert!(
            (a.roll + b.roll).abs() < 1e-6 * a.roll.abs(),
            "{:?} vs {:?}",
            a.roll,
            b.roll
        );
        assert!((a.heave - b.heave).abs() < 1e-6 * a.heave.abs());
        // Oblique seas from either side on a symmetric catamaran.
        let cat = [
            (&main, Placement { x: 0.0, y: 0.8 }),
            (&main, Placement { x: 0.0, y: -0.8 }),
        ];
        let m = mass(&cat);
        let from_port = Wave {
            heading: 0.75 * PI,
            ..head
        };
        let from_star = Wave {
            heading: 1.25 * PI,
            ..head
        };
        let (p, q) = (
            response_fleet(&cat, &m, &from_port, &opts()).unwrap(),
            response_fleet(&cat, &m, &from_star, &opts()).unwrap(),
        );
        assert!((p.heave - q.heave).abs() < 1e-6 * p.heave.abs());
        for (x, y) in [(p.sway, q.sway), (p.roll, q.roll), (p.yaw, q.yaw)] {
            assert!((x + y).abs() < 1e-6 * x.abs().max(1e-12), "{x:?} vs {y:?}");
        }
        let straight = response_fleet(&cat, &m, &head, &opts()).unwrap();
        assert!(
            straight.sway.abs() < 1e-12
                && straight.roll.abs() < 1e-12
                && straight.yaw.abs() < 1e-12
        );
    }

    /// Head seas. At rest the heave resonance lies in short waves that
    /// barely excite it, so heave rises monotonically toward 1 with wave
    /// length; at Fn 0.3 the encounter frequency meets the natural frequency
    /// where the waves are about a hull length long, and heave and pitch
    /// resonate well above 1 there — the familiar Wigley shape. The bounds
    /// are qualitative.
    #[test]
    fn wigley_head_seas_have_the_expected_shape() {
        let l = 3.0;
        let hull = wigley(l);
        let mass = MassProperties::floating(&hull, RHO, 0.25 * l);
        let sweep = |fnum: f64| -> Vec<(f64, f64, f64)> {
            let u = fnum * (G * l).sqrt();
            (0..12)
                .map(|i| {
                    let lam = l * (0.6 + 0.2 * i as f64);
                    let k = 2.0 * PI / lam;
                    let wave = Wave {
                        omega: (k * G).sqrt(),
                        heading: PI,
                        speed: u,
                    };
                    let r = response(&hull, &mass, &wave, &opts()).unwrap();
                    (lam / l, r.heave_rao(), r.pitch_rao())
                })
                .collect()
        };
        let still = sweep(0.0);
        for w in still.windows(2).skip(1) {
            assert!(w[1].1 > w[0].1, "Fn 0 heave not rising: {w:?}");
        }
        assert!(still.last().unwrap().1 < 1.0);
        assert!(
            (still.last().unwrap().2 - 1.0).abs() < 0.1,
            "Fn 0 long-wave pitch {:?}",
            still.last()
        );
        let fast = sweep(0.3);
        let heave = fast
            .iter()
            .cloned()
            .fold((0.0, 0.0, 0.0), |a, r| if r.1 > a.1 { r } else { a });
        let pitch = fast
            .iter()
            .cloned()
            .fold((0.0, 0.0, 0.0), |a, r| if r.2 > a.2 { r } else { a });
        assert!(
            heave.1 > 1.2 && heave.1 < 2.0 && (0.9..=1.6).contains(&heave.0),
            "Fn 0.3 heave peak {heave:?}"
        );
        assert!(
            pitch.2 > 1.2 && pitch.2 < 2.2 && (0.9..=1.8).contains(&pitch.0),
            "Fn 0.3 pitch peak {pitch:?}"
        );
    }
}

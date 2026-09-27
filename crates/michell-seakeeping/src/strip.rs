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
use crate::section2d::{HeaveSolution, Section};
use michell_geometry::parallel::map_indexed;
use michell_geometry::{Error, Result, SectionalHull, C64};

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
    /// subtracts `ρg∇·BG` from the pitch stiffness (0 leaves only the
    /// waterplane's `ρg I_L`, adequate for a slender hull).
    pub bg: f64,
}

impl MassProperties {
    /// A hull floating freely at its cut attitude: mass `ρ∇`, centre of
    /// gravity over the centre of buoyancy, radius of gyration `k_yy`.
    pub fn floating(hull: &SectionalHull, density: f64, radius_of_gyration: f64) -> Self {
        MassProperties {
            mass: density * hull.displaced_volume(),
            lcg: hull.lcb_x(),
            radius_of_gyration,
            bg: 0.0,
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
}

impl Default for StripOptions {
    fn default() -> Self {
        StripOptions {
            panels: 20,
            density: 1025.0,
            gravity: michell_geometry::STANDARD_GRAVITY,
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
}

/// One station's 2-D solution (or none, for a dry or degenerate station).
struct Strip {
    x: f64,
    sol: Option<HeaveSolution>,
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
    let (rho, g) = (opts.density, opts.gravity);
    let w0 = wave.omega;
    let k = w0 * w0 / g;
    let cb = wave.heading.cos();
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
    // Stations and their 2-D solutions at the encounter frequency.
    let curves: Vec<(f64, Vec<(f64, f64)>)> = hull.curves().map(|(x, c)| (x, c.to_vec())).collect();
    let strips: Vec<Strip> = map_indexed(
        curves.len(),
        || (),
        |_, i| {
            let (x, c) = &curves[i];
            let sol = Section::from_curve(c, opts.panels)
                .filter(|s| s.half_beam() > 0.0 && s.draft() > 0.0)
                .and_then(|s| s.heave(we, g, rho));
            Strip { x: x - xg, sol }
        },
    );
    let real = |v: f64| C64::new(v, 0.0);
    let a0 = trapz(&strips, |s| real(s.a())).re;
    let b0 = trapz(&strips, |s| real(s.b())).re;
    let xa1 = trapz(&strips, |s| real(s.x * s.a())).re;
    let xb1 = trapz(&strips, |s| real(s.x * s.b())).re;
    let xa2 = trapz(&strips, |s| real(s.x * s.x * s.a())).re;
    let xb2 = trapz(&strips, |s| real(s.x * s.x * s.b())).re;
    // The transom: the aft station's section, if it has one.
    let aft = strips.first().filter(|s| s.sol.is_some());
    let (x_a, a_a, b_a) = aft.map_or((0.0, 0.0, 0.0), |s| (s.x, s.a(), s.b()));
    let (uw, uw2) = (u / (we * we), u * u / (we * we));
    // STF convention (e^{iωt}, pitch bow down).
    let a33 = a0 - uw * b_a;
    let b33 = b0 + u * a_a;
    let a35 = -xa1 - uw * b0 + uw * x_a * b_a;
    let b35 = -xb1 + u * a0 - u * x_a * a_a;
    let a53 = -xa1 + uw * b0 + uw * x_a * b_a;
    let b53 = -xb1 - u * a0 - u * x_a * a_a;
    let a55 = xa2 + uw2 * a0 - uw * x_a * x_a * b_a + uw2 * x_a * a_a;
    let b55 = xb2 + uw2 * b0 + u * x_a * x_a * a_a + uw2 * x_a * b_a;
    // Restoring about G (this crate's convention; C35 flips with pitch).
    let c = restoring(hull, rho, g, xg);
    let c55 = c.c55 - rho * g * hull.displaced_volume() * mass.bg;
    let (c33, c35) = (c.c33, -c.c35);
    let m = [mass.mass, mass.mass * mass.radius_of_gyration.powi(2)];
    // Exciting force. Froude–Krylov: whole-hull closed form, converted.
    let fk = froude_krylov(hull, rho, g, k, wave.heading, xg);
    let (fk3, fk5) = (fk.heave.conj(), -fk.pitch.conj());
    // Diffraction per station (this crate's convention, then conjugated).
    let hd = |s: &Strip| -> C64 {
        s.sol.as_ref().map_or(C64::ZERO, |sol| {
            let f = sol.diffraction(k, wave.heading, g, rho).scale(we / w0);
            (f * C64::cis(k * cb * s.x)).conj()
        })
    };
    let h0 = trapz(&strips, &hd);
    let h1 = trapz(&strips, |s| hd(s).scale(s.x));
    let h_a = aft.map_or(C64::ZERO, &hd);
    let u_iw = C64::new(0.0, -u / we); // U/(iω)
    let fd3 = h0 + u_iw * h_a;
    let fd5 = -h1 - u_iw * h0 - u_iw * h_a.scale(x_a);
    // Solve [−ω²(M + A) + iωB + C] η = F (STF convention).
    let iw = C64::new(0.0, we);
    let w2 = we * we;
    let z = |mm: f64, aa: f64, bb: f64, cc: f64| real(-w2 * (mm + aa) + cc) + iw.scale(bb);
    let z33 = z(m[0], a33, b33, c33);
    let z35 = z(0.0, a35, b35, c35);
    let z53 = z(0.0, a53, b53, c35);
    let z55 = z(m[1], a55, b55, c55);
    let (f3, f5) = (fk3 + fd3, fk5 + fd5);
    let det = z33 * z55 - z35 * z53;
    if det.abs() == 0.0 {
        return Err(Error::InvalidInput("singular heave–pitch system".into()));
    }
    let eta3 = (f3 * z55 - z35 * f5) / det;
    let eta5 = (z33 * f5 - z53 * f3) / det;
    // To this crate's convention: conjugate (time), negate pitch (bow up).
    let coefficients = Coefficients {
        added_mass: [[a33, -a35], [-a53, a55]],
        damping: [[b33, -b35], [-b53, b55]],
        restoring: [[c33, -c35], [-c35, c55]],
        mass: m,
        froude_krylov: [fk3.conj(), -fk5.conj()],
        diffraction: [fd3.conj(), -fd5.conj()],
    };
    Ok(Response {
        wave: *wave,
        omega_e: we,
        k,
        heave: eta3.conj(),
        pitch: -eta5.conj(),
        coefficients,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use michell_geometry::iges::{self, HullPose, Platform, SectionalOptions};
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

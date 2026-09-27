//! The **Froude–Krylov** heave force and pitch moment: the undisturbed
//! incident wave's pressure integrated over the hull.
//!
//! The linear dynamic pressure of a deep-water wave of unit amplitude is
//! `p = ρ g e^{−kz} e^{i k_x (x − x_ref)}` (`z` down, `k_x = k cos β`).
//! Closing the wetted surface with its waterplane and applying Gauss's
//! theorem, the upward force per unit length is the waterplane's pressure
//! less the volume's pressure gradient:
//!
//! ```text
//! f₃(x) = ρ g e^{i k_x (x − x_ref)} [ b(x) − 2k Z(x; k) ],
//! Z(x; k) = ∫ f(x, z) e^{−kz} dz,     b(x) = 2 f(x, 0)
//! ```
//!
//! and `F₃ = ∫ f₃ dx`, `M₅ = ∫ (x − x_ref) f₃ dx`. Both are the sectional
//! hull's own closed-form transforms along `x`
//! ([`SectionalHull::x_transform`], [`SectionalHull::waterline_transform`]),
//! so they carry no quadrature error beyond the stations' depth integrals.
//!
//! In oblique seas the pressure also varies across each section, as
//! `cos(κy)`, `κ = k sin β`. Green's theorem turns the section's share into
//! a line integral along its curve,
//!
//! ```text
//! f₃(x) = ρ g e^{i k_x (x − x_ref)} [ 2 sin(κ f₀)/κ − 2k ∫ (sin(κy)/κ) e^{−kz} dz ]
//! ```
//!
//! (`f₀` the waterline half-beam), which reduces to the uniform form as
//! `κ → 0`. It is added station by station as a correction to the closed
//! form — the exact minus the uniform section force — so the head- and
//! following-sea result is the closed form itself and the oblique one joins
//! it continuously. As `k → 0` the force tends to the hydrostatic `C₃₃` and
//! `C₃₅` of a unit rise of the water ([`crate::restoring`]).

use michell_geometry::sectional::SectionalContracted;
use michell_geometry::{SectionalHull, C64};

/// Complex heave force and pitch moment per unit wave amplitude, with phase
/// relative to a crest at `x_ref`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WaveLoad {
    /// Upward force [N/m of wave amplitude].
    pub heave: C64,
    /// Bow-up moment about `x_ref` [N·m/m].
    pub pitch: C64,
}

/// The Froude–Krylov load on `hull` from a unit-amplitude wave of
/// wavenumber `k` at `heading` (`π` head seas), about `x_ref` in the hull's
/// own x, in a fluid of `density` under `gravity`.
pub fn froude_krylov(
    hull: &SectionalHull,
    density: f64,
    gravity: f64,
    k: f64,
    heading: f64,
    x_ref: f64,
) -> WaveLoad {
    let kx = k * heading.cos();
    let mut zc = SectionalContracted::default();
    hull.contract(k, &mut zc);
    // Transforms about x_c, phase-referenced to x_c.
    let (z0, z1) = hull.x_transform(&zc, kx);
    let (w0, w1) = hull.waterline_transform(kx);
    let f0 = (w0 - z0.scale(k)).scale(2.0);
    let f1 = (w1 - z1.scale(k)).scale(2.0);
    // Re-reference to x_ref: e^{ik_x(x − x_ref)} = e^{ik_x(x − x_c)} e^{ik_x d},
    // and (x − x_ref) = (x − x_c) + d, with d = x_c − x_ref.
    let d = hull.x_center() - x_ref;
    let phase = C64::cis(kx * d);
    let rg = density * gravity;
    let (ch, cp) = transverse_correction(hull, k, heading, x_ref);
    WaveLoad {
        heave: (phase * f0 + ch).scale(rg),
        pitch: (phase * (f1 + f0.scale(d)) + cp).scale(rg),
    }
}

/// The oblique-sea correction to the uniform-pressure force and moment
/// (per `ρg`): each station's exact section force less its uniform one,
/// integrated along the stations by the trapezoidal rule. Zero when
/// `sin β = 0`.
fn transverse_correction(hull: &SectionalHull, k: f64, heading: f64, x_ref: f64) -> (C64, C64) {
    let kappa = k * heading.sin();
    if kappa.abs() < 1e-12 * k.max(1e-300) {
        return (C64::ZERO, C64::ZERO);
    }
    let kx = k * heading.cos();
    let sinc = |y: f64| (kappa * y).sin() / kappa;
    let stations: Vec<(f64, f64)> = hull
        .curves()
        .map(|(x, curve)| {
            let f0 = curve.first().map_or(0.0, |p| p.0);
            // ∫ (g(y) − y) e^{−kz} dz along the curve, g = sin(κy)/κ.
            let line: f64 = curve
                .windows(2)
                .map(|w| {
                    let (a, b) = (w[0], w[1]);
                    let ga = (sinc(a.0) - a.0) * (-k * a.1).exp();
                    let gb = (sinc(b.0) - b.0) * (-k * b.1).exp();
                    0.5 * (ga + gb) * (b.1 - a.1)
                })
                .sum();
            (x, 2.0 * (sinc(f0) - f0) - 2.0 * k * line)
        })
        .collect();
    let (mut heave, mut pitch) = (C64::ZERO, C64::ZERO);
    for w in stations.windows(2) {
        let dx = w[1].0 - w[0].0;
        for &(x, c) in &[w[0], w[1]] {
            let f = C64::cis(kx * (x - x_ref)).scale(0.5 * dx * c);
            heave = heave + f;
            pitch = pitch + f.scale(x - x_ref);
        }
    }
    (heave, pitch)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::restoring::restoring;
    use michell_geometry::iges::{self, HullPose, Platform, SectionalOptions};
    use michell_geometry::quadrature::gauss_legendre;
    use std::f64::consts::PI;

    const RHO: f64 = 1025.0;
    const G: f64 = 9.81;
    const L: f64 = 10.0;
    const B: f64 = 1.0;
    const T: f64 = 0.625;

    fn wigley() -> SectionalHull {
        let surfaces = iges::wigley_surfaces(L, B, T).unwrap();
        let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
        source
            .situate_sectional(
                0,
                0.0,
                &HullPose::default(),
                &Platform::default(),
                &SectionalOptions::default(),
            )
            .unwrap()
            .expect("the hull is wet")
            .hull
    }

    /// The load by brute-force quadrature of the Wigley half-breadth
    /// `f = (B/2)(1 − ξ²)(1 − (z/T)²)` itself, `ξ = 2(x − x_c)/L`.
    fn brute_force(x_c: f64, k: f64, heading: f64, x_ref: f64) -> WaveLoad {
        let kx = k * heading.cos();
        let (gx, gw) = gauss_legendre(40);
        let panels = 40;
        let half = |xi: f64, z: f64| 0.5 * B * (1.0 - xi * xi) * (1.0 - (z / T).powi(2));
        // ∫ cos(κy) over the section's width at depth z: 2 sin(κ f)/κ.
        let kappa = k * heading.sin();
        let f = |xi: f64, z: f64| {
            let h = half(xi, z);
            if kappa.abs() < 1e-12 { h } else { (kappa * h).sin() / kappa }
        };
        let (mut heave, mut pitch) = (C64::ZERO, C64::ZERO);
        for p in 0..panels {
            let (a, b) = (
                -L / 2.0 + L * p as f64 / panels as f64,
                -L / 2.0 + L * (p + 1) as f64 / panels as f64,
            );
            for (&t, &w) in gx.iter().zip(&gw) {
                let s = 0.5 * (a + b) + 0.5 * (b - a) * t;
                let wx = 0.5 * (b - a) * w;
                let xi = 2.0 * s / L;
                // Z(x; k) by the same rule in depth.
                let z_int: f64 = gx
                    .iter()
                    .zip(&gw)
                    .map(|(&tz, &wz)| {
                        let z = 0.5 * T * (1.0 + tz);
                        0.5 * T * wz * f(xi, z) * (-k * z).exp()
                    })
                    .sum();
                let x = x_c + s;
                let f3 = C64::cis(kx * (x - x_ref))
                    .scale(RHO * G * (2.0 * f(xi, 0.0) - 2.0 * k * z_int));
                heave = heave + f3.scale(wx);
                pitch = pitch + f3.scale(wx * (x - x_ref));
            }
        }
        WaveLoad { heave, pitch }
    }

    fn close(a: C64, b: C64, scale: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * scale
    }

    #[test]
    fn matches_brute_force_quadrature_of_the_wigley_hull() {
        let hull = wigley();
        let x_c = hull.x_center();
        let x_ref = x_c - 0.3; // off-centre pivot: exercises the re-referencing
        // Head and following seas: the closed form, to quadrature precision.
        // Oblique and beam seas: the station-by-station transverse
        // correction, to the trapezoidal rule's precision along the stations.
        for (k, heading, tol) in [
            (2.0 * PI / L, PI, 1e-6),
            (2.0 * PI / (0.5 * L), PI, 1e-6),
            (2.0 * PI / L, 0.0, 1e-6),
            (1.3, 2.4, 1e-3),
            (4.0, 0.5 * PI, 1e-3),
            (2.0 * PI / (0.4 * L), 2.0, 1e-3),
        ] {
            let got = froude_krylov(&hull, RHO, G, k, heading, x_ref);
            let want = brute_force(x_c, k, heading, x_ref);
            let scale = RHO * G * hull.waterplane_area();
            assert!(close(got.heave, want.heave, scale, tol), "k {k} β {heading}: heave {:?} vs {:?}", got.heave, want.heave);
            assert!(close(got.pitch, want.pitch, scale * L, tol), "k {k} β {heading}: pitch {:?} vs {:?}", got.pitch, want.pitch);
        }
    }

    #[test]
    fn a_long_wave_is_a_hydrostatic_rise() {
        let hull = wigley();
        let x_ref = hull.x_center() + 0.7;
        let c = restoring(&hull, RHO, G, x_ref);
        let fk = froude_krylov(&hull, RHO, G, 1e-7, PI, x_ref);
        assert!(
            (fk.heave.re - c.c33).abs() < 1e-5 * c.c33,
            "{} vs {}",
            fk.heave.re,
            c.c33
        );
        assert!(
            (fk.pitch.re - c.c35).abs() < 1e-5 * c.c33 * L,
            "{} vs {}",
            fk.pitch.re,
            c.c35
        );
        assert!(fk.heave.im.abs() < 1e-5 * c.c33);
    }
}

//! Against experiment: Journée's four Wigley hulls in head waves (DUT-SHL
//! Report 0909, 1992; data freely distributed by the author). The values
//! below are the report's own, for Wigley III (C_m 0.667, L/B 10) and, for
//! the zero-speed motions, also IV (L/B 5).
//!
//! What strip theory is held to here is what it gets right on these hulls:
//! the heave coefficients, the wave loads, and heave away from resonance;
//! pitch at zero speed loosely. Pitch at speed is not asserted — the pitch
//! added inertia comes out ~30% low and the pitch damping grows with U²
//! where the tank shows none, which moves the pitch resonance to shorter
//! waves; Journée (MARIND 2001, "Discrepancies in Hydrodynamic Coefficients
//! of Wigley Hull Forms") found the same with his own Frank and Ursell
//! codes. The `journee_wigley` example prints the full comparison.

use crate::strip::{added_resistance, response, MassProperties, StripOptions, Wave};
use michell_geometry::sectional::{DepthQuadrature, SectionNodes, SectionalHull};
use std::f64::consts::PI;

const L: f64 = 3.0;
const D: f64 = 0.1875;
const RHO: f64 = 1000.0;
const G: f64 = 9.81;

/// Journée's Wigley `model` (1–4): `η = (1 − ζ²)(1 − ξ²)(1 + 0.2ξ²) +
/// α ζ²(1 − ζ⁸)(1 − ξ²)⁴`, α = 1 for I and II, 0 for III and IV; B = 0.3 m
/// (I, III) or 0.6 m (II, IV). Returns the hull and its beam.
pub(crate) fn journee_wigley(model: usize, spans: usize) -> (SectionalHull, f64) {
    let b = if model % 2 == 1 { 0.3 } else { 0.6 };
    let alpha = if model <= 2 { 1.0 } else { 0.0 };
    let half = move |x: f64, z: f64| {
        let xi = 2.0 * x / L;
        let ze = z / D;
        let eta = (1.0 - ze * ze) * (1.0 - xi * xi) * (1.0 + 0.2 * xi * xi)
            + alpha * ze * ze * (1.0 - ze.powi(8)) * (1.0 - xi * xi).powi(4);
        0.5 * b * eta.max(0.0)
    };
    let p = 3usize;
    let mut knots = vec![-0.5 * L; p + 1];
    for i in 1..spans {
        knots.push(-0.5 * L + L * i as f64 / spans as f64);
    }
    knots.extend(vec![0.5 * L; p + 1]);
    let n = knots.len() - p - 1;
    let xs: Vec<f64> = (0..n)
        .map(|i| knots[i + 1..=i + p].iter().sum::<f64>() / p as f64)
        .collect();
    let q = DepthQuadrature::default();
    let sections = xs
        .iter()
        .map(|&x| SectionNodes::from_depth_function(D, &[], |z| half(x, z), &q))
        .collect();
    (SectionalHull::new(p, knots, &xs, sections).unwrap(), b)
}

fn setup(model: usize) -> (SectionalHull, f64, MassProperties, StripOptions) {
    let (hull, b) = journee_wigley(model, 30);
    let mass = MassProperties {
        mass: RHO * hull.displaced_volume(),
        lcg: 0.0,
        radius_of_gyration: 0.75,
        bg: 0.0,
        roll_radius_of_gyration: 0.35 * b,
        yaw_radius_of_gyration: 0.75,
    };
    let opts = StripOptions {
        panels: 16,
        density: RHO,
        gravity: G,
    };
    (hull, b, mass, opts)
}

fn head_wave(lam_over_l: f64, fnum: f64) -> Wave {
    let k = 2.0 * PI / (lam_over_l * L);
    Wave {
        omega: (k * G).sqrt(),
        heading: PI,
        speed: fnum * (G * L).sqrt(),
    }
}

fn within(got: f64, want: f64, rel: f64) -> bool {
    (got - want).abs() <= rel * want.abs()
}

#[test]
fn the_hull_has_the_reported_displacement() {
    for (model, vol) in [(1, 0.0946), (2, 0.1892), (3, 0.0780), (4, 0.1560)] {
        let (hull, _) = journee_wigley(model, 30);
        assert!(
            within(hull.displaced_volume(), vol, 5e-3),
            "Wigley {model}: ∇ {}",
            hull.displaced_volume()
        );
    }
}

/// Table 2-III, Fn 0.3: forced-heave added mass `A33/(ρ∇)` and damping
/// `B33/(ρ∇)·√(L/g)` against `ω√(L/g)`.
#[test]
fn heave_coefficients_match_wigley_iii() {
    let (hull, _, mass, opts) = setup(3);
    let vol = hull.displaced_volume();
    let u = 0.3 * (G * L).sqrt();
    for &(wn, a_exp, b_exp) in &[
        (2.22, 0.707, 1.963),
        (2.78, 0.537, 1.994),
        (3.34, 0.465, 1.854),
        (3.88, 0.425, 1.706),
        (4.45, 0.421, 1.518),
    ] {
        let we = wn / (L / G).sqrt();
        // The head wave met at this encounter frequency: ω_e = ω₀ + ω₀²U/g.
        let a = u / G;
        let w0 = (-1.0 + (1.0 + 4.0 * a * we).sqrt()) / (2.0 * a);
        let wave = Wave {
            omega: w0,
            heading: PI,
            speed: u,
        };
        let c = response(&hull, &mass, &wave, &opts).unwrap().coefficients;
        let a_c = c.added_mass[0][0] / (RHO * vol);
        let b_c = c.damping[0][0] / (RHO * vol) * (L / G).sqrt();
        assert!(
            within(a_c, a_exp, 0.15),
            "ω' {wn}: A33'' {a_c:.3} vs {a_exp}"
        );
        assert!(
            within(b_c, b_exp, 0.15),
            "ω' {wn}: B33'' {b_c:.3} vs {b_exp}"
        );
    }
}

/// Table 6-III, Fn 0.3: wave force `|X3|/(C33 ζ)` and moment
/// `|X5|/(k C55 ζ)` on the restrained hull.
#[test]
fn wave_loads_match_wigley_iii() {
    let (hull, _, mass, opts) = setup(3);
    for &(lam, x3, x5) in &[
        (1.25, 0.345, 0.444),
        (1.5, 0.454, 0.546),
        (1.75, 0.515, 0.600),
        (2.0, 0.577, 0.658),
    ] {
        let wave = head_wave(lam, 0.3);
        let k = wave.omega * wave.omega / G;
        let c = response(&hull, &mass, &wave, &opts).unwrap().coefficients;
        let f3 = (c.froude_krylov[0] + c.diffraction[0]).abs() / c.restoring[0][0];
        let f5 = (c.froude_krylov[1] + c.diffraction[1]).abs() / (k * c.restoring[1][1]);
        assert!(within(f3, x3, 0.10), "λ/L {lam}: X3'' {f3:.3} vs {x3}");
        assert!(within(f5, x5, 0.15), "λ/L {lam}: X5'' {f5:.3} vs {x5}");
    }
}

/// Table 8-III: heave RAO in head waves at Fn 0.2 and 0.3, away from the
/// pitch resonance.
#[test]
fn heave_matches_wigley_iii_in_head_waves() {
    let (hull, _, mass, opts) = setup(3);
    for &(fnum, lam, z) in &[
        (0.2, 1.25, 0.71),
        (0.2, 1.5, 0.81),
        (0.2, 2.0, 0.90),
        (0.3, 1.25, 1.28),
        (0.3, 1.384, 1.11),
        (0.3, 1.5, 1.01),
        (0.3, 2.0, 0.99),
    ] {
        let r = response(&hull, &mass, &head_wave(lam, fnum), &opts).unwrap();
        assert!(
            within(r.heave_rao(), z, 0.15),
            "Fn {fnum} λ/L {lam}: heave {:.3} vs {z}",
            r.heave_rao()
        );
    }
}

/// Table 8-I and 8-III: the peak added resistance at Fn 0.3 and 0.4, where
/// Gerritsma–Beukelman overshoots it two- to six-fold. Maruo's far-field
/// estimate lands within a factor of two of the measured peak (taken over
/// the resonance band, since the computed pitch resonance falls at
/// slightly shorter waves).
#[test]
fn far_field_added_resistance_peaks_near_the_measured_ones() {
    for (model, fnum, measured) in [
        (1, 0.3, 27.6),
        (3, 0.3, 20.2),
        (1, 0.4, 14.9),
        (3, 0.4, 27.7),
    ] {
        let (hull, b, mass, opts) = setup(model);
        let peak = [1.0, 1.05, 1.108, 1.25, 1.384, 1.5]
            .iter()
            .map(|&lam| {
                crate::strip::added_resistance_maruo(&hull, &mass, &head_wave(lam, fnum), &opts)
                    .unwrap()
                    .coefficient(RHO, G, b, L)
            })
            .fold(0.0f64, f64::max);
        assert!(
            peak > 0.5 * measured && peak < 2.0 * measured,
            "Wigley {model} Fn {fnum}: Maruo peak {peak:.1} vs measured {measured}"
        );
    }
}

/// Table 10-III and 10-IV: heave and pitch at zero speed, pitch in the
/// report's normalisation `θ_a/(2π ζ_a/L)`. Pitch comes out 10–20% low.
#[test]
fn zero_speed_motions_match_wigley_iii_and_iv() {
    for (model, points) in [
        (3, [(1.0, 0.29, 0.57), (1.5, 0.57, 0.57), (2.0, 0.74, 0.50)]),
        (4, [(1.0, 0.30, 0.57), (1.5, 0.65, 0.69), (2.0, 0.79, 0.54)]),
    ] {
        let (hull, b, mass, opts) = setup(model);
        for (lam, z, th) in points {
            let r = added_resistance(&hull, &mass, &head_wave(lam, 0.0), &opts).unwrap();
            let pitch = r.response.pitch_rao() / lam;
            assert!(
                within(r.response.heave_rao(), z, 0.12),
                "Wigley {model} λ/L {lam}: heave {:.3} vs {z}",
                r.response.heave_rao()
            );
            assert!(
                within(pitch, th, 0.25),
                "Wigley {model} λ/L {lam}: pitch {pitch:.3} vs {th}"
            );
            assert!(r.coefficient(RHO, G, b, L) >= 0.0);
        }
    }
}

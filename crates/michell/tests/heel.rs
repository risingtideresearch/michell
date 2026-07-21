//! Heeled-hull wave resistance via the tilted-centreplane (complex-κ) kernel.
//!
//! Heel keeps the sources on the ship's own tilted centreplane, turning the
//! vertical decay complex: `κ = νλ²cosφ + i·νλ√(λ²−1)·sinφ`. The upright kernel
//! is untouched, so `heel = 0` must reproduce `wave_resistance` exactly, the
//! result must be even in the heel angle, and (at fixed displacement) heel
//! raises the wave resistance — growing like `sin²φ` for small angles.

use michell::{heel_wave_resistance, hulls, wave_resistance, Conditions, WaveOptions};

const DEG: f64 = std::f64::consts::PI / 180.0;

fn hull_and_cond(u: f64) -> (michell::Hull, Conditions) {
    (hulls::wigley(10.0, 1.0, 0.625).unwrap(), Conditions::seawater(u))
}

/// Zero heel reproduces the upright Michell resistance to full precision — the
/// complex kernel collapses onto the real one.
#[test]
fn zero_heel_reproduces_upright() {
    for u in [2.0, 3.0, 4.5] {
        let (hull, cond) = hull_and_cond(u);
        let upright = wave_resistance(&hull, &cond).unwrap().resistance;
        let heeled0 = heel_wave_resistance(&hull, &cond, 0.0, &WaveOptions::default())
            .unwrap()
            .resistance;
        let rel = (heeled0 - upright).abs() / upright;
        assert!(rel < 1e-10, "u={u}: heel-0 {heeled0} vs upright {upright} (rel {rel:e})");
    }
}

/// Resistance is even in the heel angle: port and starboard heel are mirror
/// images of the same physical flow.
#[test]
fn resistance_is_even_in_heel() {
    let (hull, cond) = hull_and_cond(3.0);
    let opts = WaveOptions::default();
    for deg in [5.0, 12.0, 25.0] {
        let plus = heel_wave_resistance(&hull, &cond, deg * DEG, &opts).unwrap().resistance;
        let minus = heel_wave_resistance(&hull, &cond, -deg * DEG, &opts).unwrap().resistance;
        assert!((plus - minus).abs() < 1e-9 * plus, "±{deg}°: {plus} vs {minus}");
    }
}

/// Heel raises the wave resistance (the tilt reduces the effective depth-decay
/// and radiates an asymmetric pattern), and the excess grows like sin²φ for
/// small angles.
#[test]
fn heel_adds_resistance_quadratically() {
    let (hull, cond) = hull_and_cond(3.0);
    let opts = WaveOptions::default();
    let r0 = heel_wave_resistance(&hull, &cond, 0.0, &opts).unwrap().resistance;
    let small = 4.0 * DEG;
    let big = 8.0 * DEG;
    let d_small = heel_wave_resistance(&hull, &cond, small, &opts).unwrap().resistance - r0;
    let d_big = heel_wave_resistance(&hull, &cond, big, &opts).unwrap().resistance - r0;
    assert!(d_small > 0.0, "heel must add resistance: {d_small}");
    // sin²φ scaling: ratio ≈ sin²(8°)/sin²(4°).
    let expect = (big.sin() / small.sin()).powi(2);
    let ratio = d_big / d_small;
    assert!(
        (ratio - expect).abs() < 0.08 * expect,
        "sin²φ scaling: ratio {ratio} vs {expect}"
    );
}

/// Beyond the linear regime the effect keeps growing monotonically with heel.
#[test]
fn resistance_increases_monotonically_with_heel() {
    let (hull, cond) = hull_and_cond(3.5);
    let opts = WaveOptions::default();
    let mut prev = -1.0;
    for deg in [0.0, 5.0, 10.0, 20.0, 30.0] {
        let r = heel_wave_resistance(&hull, &cond, deg * DEG, &opts).unwrap().resistance;
        assert!(r > prev, "heel {deg}°: {r} did not exceed previous {prev}");
        prev = r;
    }
}

/// A capsizing heel (|φ| ≥ 90°) is rejected.
#[test]
fn extreme_heel_is_rejected() {
    let (hull, cond) = hull_and_cond(3.0);
    assert!(heel_wave_resistance(&hull, &cond, 1.6, &WaveOptions::default()).is_err());
}

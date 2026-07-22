//! Heeled-hull wave resistance via the tilted-centreplane (complex-κ) kernel.
//!
//! Heel keeps the sources on the ship's own tilted centreplane, turning the
//! vertical decay complex: `κ = νλ²cosφ + i·νλ√(λ²−1)·sinφ`. The upright kernel
//! is untouched, so `heel = 0` must reproduce `wave_resistance` exactly, the
//! result must be even in the heel angle, and (at fixed displacement) heel
//! raises the wave resistance — growing like `sin²φ` for small angles.

use michell::{
    heel_wave_resistance, hulls, multihull_heel_wave_resistance, multihull_wave_resistance_with,
    wave_resistance, Conditions, Placement, WaveOptions,
};

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

// --- Multihull heel -------------------------------------------------------

/// A fleet of one heeled hull is exactly the single-hull path (the placement
/// phase is unity), so `heel_wave_resistance` and the fleet function must agree
/// to the bit.
#[test]
fn fleet_of_one_matches_single_hull() {
    let (hull, cond) = hull_and_cond(3.0);
    let opts = WaveOptions::default();
    for deg in [0.0, 7.0, 20.0] {
        let phi = deg * DEG;
        let single = heel_wave_resistance(&hull, &cond, phi, &opts).unwrap().resistance;
        let fleet = multihull_heel_wave_resistance(&[(&hull, Placement::default())], &cond, phi, &opts)
            .unwrap()
            .resistance;
        assert_eq!(single, fleet, "{deg}°: single {single} vs fleet-of-one {fleet}");
    }
}

/// Zero heel on a catamaran reproduces the upright multihull resistance (with
/// its interference) to full precision.
#[test]
fn zero_heel_catamaran_reproduces_upright() {
    let (hull, cond) = hull_and_cond(3.0);
    let opts = WaveOptions::default();
    let members = [
        (&hull, Placement { x: 0.0, y: 3.0 }),
        (&hull, Placement { x: 0.0, y: -3.0 }),
    ];
    let upright = multihull_wave_resistance_with(&members, &cond, &opts).unwrap().resistance;
    let heeled0 = multihull_heel_wave_resistance(&members, &cond, 0.0, &opts).unwrap().resistance;
    let rel = (heeled0 - upright).abs() / upright;
    assert!(rel < 1e-10, "heel-0 {heeled0} vs upright {upright} (rel {rel:e})");
}

/// A catamaran that is mirror-symmetric about its mean centreplane heels evenly:
/// reflecting y → −y maps the fleet to itself and φ → −φ, so R(+φ) = R(−φ).
/// (An asymmetric or staggered arrangement need not be even — hence the
/// per-half-system treatment.)
#[test]
fn symmetric_catamaran_heel_is_even() {
    let (hull, cond) = hull_and_cond(3.0);
    let opts = WaveOptions::default();
    let members = [
        (&hull, Placement { x: 0.0, y: 2.5 }),
        (&hull, Placement { x: 0.0, y: -2.5 }),
    ];
    for deg in [6.0, 15.0, 28.0] {
        let plus = multihull_heel_wave_resistance(&members, &cond, deg * DEG, &opts)
            .unwrap()
            .resistance;
        let minus = multihull_heel_wave_resistance(&members, &cond, -deg * DEG, &opts)
            .unwrap()
            .resistance;
        assert!((plus - minus).abs() < 1e-9 * plus, "±{deg}°: {plus} vs {minus}");
    }
}

/// Heel raises a catamaran's wave resistance too — the per-demihull tilt effect
/// survives the interference.
#[test]
fn heel_adds_resistance_to_catamaran() {
    let (hull, cond) = hull_and_cond(3.0);
    let opts = WaveOptions::default();
    let members = [
        (&hull, Placement { x: 0.0, y: 2.5 }),
        (&hull, Placement { x: 0.0, y: -2.5 }),
    ];
    let r0 = multihull_heel_wave_resistance(&members, &cond, 0.0, &opts).unwrap().resistance;
    let r20 = multihull_heel_wave_resistance(&members, &cond, 20.0 * DEG, &opts)
        .unwrap()
        .resistance;
    assert!(r20 > r0, "heel must add resistance to the catamaran: {r20} vs {r0}");
}

/// Far-apart demihulls stop interfering: the combined heeled resistance
/// approaches the sum of the members' standalone heeled resistances.
#[test]
fn far_apart_heeled_fleet_is_additive() {
    let (hull, cond) = hull_and_cond(3.0);
    let opts = WaveOptions::default();
    let phi = 15.0 * DEG;
    let solo = heel_wave_resistance(&hull, &cond, phi, &opts).unwrap().resistance;
    // 40 hull-lengths apart transversely: the interference phase oscillates
    // fast enough that the cross term integrates away.
    let members = [
        (&hull, Placement { x: 0.0, y: 200.0 }),
        (&hull, Placement { x: 0.0, y: -200.0 }),
    ];
    let combined = multihull_heel_wave_resistance(&members, &cond, phi, &opts).unwrap().resistance;
    let ratio = combined / (2.0 * solo);
    assert!((ratio - 1.0).abs() < 0.05, "far-apart additivity: combined/2·solo = {ratio}");
}

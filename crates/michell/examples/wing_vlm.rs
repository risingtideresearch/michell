//! 3D horseshoe vortex-lattice engine — the core of the level-2 asymmetric-hull
//! (centreplane lifting-surface) solve.
//!
//! Prints the lift-curve slope of a flat rectangular wing versus aspect ratio,
//! bracketed by its two exact limits: `2π` as `AR → ∞` (2-D) and `(π/2)·AR` as
//! `AR → 0` (slender-wing / R.T. Jones), with the Prandtl lifting-line estimate
//! `2π/(1 + 2/AR)` alongside.
//!
//! Run with: cargo run --release --example wing_vlm

use michell::lifting3d::{rectangular_wing, solve};
use std::f64::consts::PI;

const DEG: f64 = PI / 180.0;

fn slope(ar: f64) -> f64 {
    let a = 2.0 * DEG;
    solve(&rectangular_wing(ar, 6, 60), a).cl / a.sin()
}

fn main() {
    println!("Flat rectangular wing — lift-curve slope dC_L/dα [1/rad]");
    println!("2-D limit 2π = {:.4}\n", 2.0 * PI);
    println!(
        "{:>6} {:>12} {:>16} {:>16}",
        "AR", "VLM slope", "LLT 2π/(1+2/AR)", "Jones (π/2)AR"
    );
    for ar in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0] {
        let vlm = slope(ar);
        let llt = 2.0 * PI / (1.0 + 2.0 / ar);
        let jones = 0.5 * PI * ar;
        println!("{ar:>6.1} {vlm:>12.4} {llt:>16.4} {jones:>16.4}");
    }
    println!(
        "\nAs AR grows the VLM slope climbs toward 2π; as AR → 0 it follows the\n\
         slender-wing line (π/2)·AR. The Prandtl estimate brackets it in between."
    );
}

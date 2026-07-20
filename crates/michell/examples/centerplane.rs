//! Hull centreplane lifting solve (double-body model) — applies the 3D vortex
//! lattice to the vertical centreplane of an asymmetric hull and reports the
//! solved doublet density that replaces the strip closure μ = 2U f_a.
//!
//! Run with: cargo run --release --example centerplane

use michell::centerplane::{drift_slope, solve_centerplane};
use michell::lifting3d::{rectangular_wing, solve};

const DEG: f64 = std::f64::consts::PI / 180.0;

fn main() {
    // Effective-aspect-ratio doubling: a draft-T chord-L plate behaves as a
    // span-2T wing because the free surface acts as a rigid wall.
    let alpha = 3.0 * DEG;
    println!("Free-surface image doubles the effective aspect ratio:");
    println!("  geometric AR = T/L, effective AR = 2T/L\n");
    println!(
        "{:>6} {:>8} {:>14} {:>16}",
        "T/L", "eff AR", "plate C_Y slope", "AR-wing C_L slope"
    );
    for t in [1.0, 2.0, 3.0, 4.0] {
        let cp = solve_centerplane(1.0, t, drift_slope(alpha.tan()), 5, 15);
        let wing = solve(&rectangular_wing(2.0 * t, 5, 30), alpha);
        println!(
            "{t:>6.1} {:>8.1} {:>14.4} {:>16.4}",
            2.0 * t,
            cp.side_force.abs() / alpha.sin(),
            wing.cl.abs() / alpha.sin()
        );
    }

    // Solved doublet density on the physical centreplane (drift/yaw case).
    let (nx, nz) = (8, 10);
    let cp = solve_centerplane(1.0, 2.0, drift_slope(0.05), nx, nz);
    println!(
        "\nSolved doublet density μ(x,z)/U at the trailing edge, by depth\n\
         (root at the waterline, relieved toward the keel tip):\n"
    );
    println!("{:>8} {:>12}", "z/T", "μ_TE/U");
    for strip in 0..nz {
        let d = cp.doublet[strip * nx + (nx - 1)];
        println!("{:>8.3} {:>12.5}", d.z / 2.0, d.mu);
    }
    println!("\nThis μ is what feeds the free-wave amplitude G(λ) = ∬ μ e^(−κz) e^(iνλx),");
    println!("replacing the prescribed strip closure in the asymmetric wave resistance.");
}

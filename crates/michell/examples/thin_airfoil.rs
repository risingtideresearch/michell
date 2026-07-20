//! 2D thin-airfoil lifting solve — the proving ground for the level-2
//! asymmetric-hull (centreplane-dipole) lifting solve.
//!
//! Cross-checks the vortex-lattice solver against the closed-form Glauert
//! result for a flat plate (C_L = 2πα) and a parabolic-camber mean line
//! (C_L = 2π(α + 2m), C_{m,c/4} = −πm).
//!
//! Run with: cargo run --release --example thin_airfoil

use michell::lifting::{glauert, parabolic_camber_slope, solve_vortex_lattice};

const DEG: f64 = std::f64::consts::PI / 180.0;

fn main() {
    println!("Flat plate — expect C_L = 2πα, x_cp = 0.25c\n");
    println!(
        "{:>6} {:>12} {:>12} {:>10} {:>10}",
        "α[°]", "C_L (VLM)", "C_L (Glau)", "2πα", "x_cp/c"
    );
    for deg in [0.0, 2.0, 4.0, 6.0, 8.0] {
        let a = deg * DEG;
        let vl = solve_vortex_lattice(a, |_| 0.0, 40);
        let gl = glauert(a, |_| 0.0, 6, 64);
        println!(
            "{deg:>6.1} {:>12.5} {:>12.5} {:>10.5} {:>10.4}",
            vl.cl,
            gl.cl,
            2.0 * std::f64::consts::PI * a,
            vl.center_of_pressure
        );
    }

    let m = 0.05;
    println!("\nParabolic camber, max-camber ratio m = {m}");
    println!("— expect C_L = 2π(α + 2m), C_m,c/4 = −πm = {:.5}\n", -std::f64::consts::PI * m);
    println!(
        "{:>6} {:>12} {:>12} {:>12} {:>12}",
        "α[°]", "C_L (VLM)", "C_L (Glau)", "C_m (VLM)", "C_m (Glau)"
    );
    let slope = parabolic_camber_slope(m);
    for deg in [0.0, 2.0, 4.0, 6.0] {
        let a = deg * DEG;
        let vl = solve_vortex_lattice(a, &slope, 120);
        let gl = glauert(a, &slope, 6, 128);
        println!(
            "{deg:>6.1} {:>12.5} {:>12.5} {:>12.5} {:>12.5}",
            vl.cl, gl.cl, vl.cm_quarter, gl.cm_quarter
        );
    }
    println!("\nZero-lift angle for this camber: α₀ = −2m = {:.3}°", -2.0 * m / DEG);
}

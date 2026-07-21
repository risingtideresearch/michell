//! Asymmetric-hull wave resistance: the camber (dipole) system from the
//! *solved* centreplane lifting distribution vs the prescribed strip closure.
//!
//! Builds a hull with a full-beam starboard side and a half-beam port side,
//! then compares, across speed, the symmetric (source-only) resistance, the
//! strip-closure asymmetric resistance, and the lifting-solve asymmetric
//! resistance.
//!
//! Run with: cargo run --release --example asymmetric_lifting

use michell::{
    asymmetric_wave_resistance_lifting, hulls, wave_resistance, BSplineSurface, Conditions, Hull,
    LiftingGrid, STANDARD_GRAVITY,
};

fn scaled(s: &BSplineSurface, k: f64) -> BSplineSurface {
    BSplineSurface::new(
        s.degree_x(),
        s.degree_z(),
        s.knots_x().to_vec(),
        s.knots_z().to_vec(),
        s.control().iter().map(|c| k * c).collect(),
    )
    .unwrap()
}

fn main() {
    let l = 10.0;
    let base = hulls::wigley(l, 1.0, 0.625).unwrap();
    let starboard = scaled(base.surface(), 1.0);
    let port = scaled(base.surface(), 0.5); // half-beam port ⇒ f_a = 0.25·base
    let asym = Hull::new_asymmetric(scaled(&port, 1.0), scaled(&starboard, 1.0)).unwrap();
    let mean = Hull::new(scaled(asym.surface(), 1.0)).unwrap(); // f_sym

    println!("Asymmetric Wigley (starboard 1.0, port 0.5) — wave resistance [N]");
    println!("dipole from the SOLVED centreplane lifting distribution vs the strip closure\n");
    println!(
        "{:>5} {:>7} {:>11} {:>11} {:>11} {:>11} {:>11}",
        "Fn", "U[m/s]", "R_source", "R_strip", "R_lift", "dip_strip", "dip_lift"
    );
    let grid = LiftingGrid { nx: 48, nz: 16 };
    let mut fr = 0.20;
    while fr <= 0.451 {
        let u = fr * (STANDARD_GRAVITY * l).sqrt();
        let cond = Conditions::seawater(u);
        let rs = wave_resistance(&mean, &cond).unwrap().resistance;
        let rstrip = wave_resistance(&asym, &cond).unwrap().resistance;
        let rlift = asymmetric_wave_resistance_lifting(&asym, &cond, &Default::default(), grid)
            .unwrap()
            .resistance;
        println!(
            "{fr:>5.2} {u:>7.3} {rs:>11.2} {rstrip:>11.2} {rlift:>11.2} {:>11.2} {:>11.2}",
            rstrip - rs,
            rlift - rs
        );
        fr += 0.05;
    }
    println!(
        "\ndip_* is the camber contribution (asymmetric − symmetric). The lifting\n\
         dipole is a genuine lifting-surface solution; the strip closure μ = 2·f_a\n\
         is the crude prescribed estimate it replaces."
    );
}

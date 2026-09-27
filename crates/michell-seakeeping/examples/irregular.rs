//! Convergence of the infinite-frequency heave added mass of a semicircle
//! toward `ρπR²/2`, and of the finite-frequency one toward it.
use michell_seakeeping::section2d::Section;
use std::f64::consts::PI;

fn main() {
    let (r, rho, g) = (1.0, 1000.0, 9.81);
    let exact = rho * PI * r * r / 2.0;
    for n in [16, 32, 64, 128, 256] {
        let a = Section::semicircle(r, n).added_mass_infinite(rho).unwrap();
        println!("N {n:4}: a∞/exact − 1 = {:+.2e}", a / exact - 1.0);
    }
    let sec = Section::semicircle(r, 64);
    for nur in [2.0, 5.0, 10.0, 20.0, 40.0] {
        let s = sec.heave((nur / r * g).sqrt(), g, rho).unwrap();
        println!(
            "νR {nur:4}: a/exact {:.4}  interpolated {}",
            s.added_mass / exact,
            s.interpolated
        );
    }
}

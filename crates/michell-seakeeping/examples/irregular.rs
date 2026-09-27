//! Added mass and damping of a Wigley midship section swept across its
//! first irregular frequency, with the energy check `b = ρω|C|²`.
use michell_seakeeping::section2d::Section;

fn main() {
    let (b, t, g, rho) = (0.15, 0.1875, 9.81, 1000.0);
    let curve: Vec<(f64, f64)> = (0..=64)
        .map(|i| {
            let z = t * i as f64 / 64.0;
            (b * (1.0 - (z / t).powi(2)), z)
        })
        .collect();
    let sec = Section::from_curve(&curve, 24).unwrap();
    let m = rho * 2.0 * b * t * 2.0 / 3.0;
    for i in 0..45 {
        let nu = 4.0 + 0.5 * i as f64;
        let w = (nu * g).sqrt();
        let s = sec.heave(w, g, rho).unwrap();
        let far = rho * w * s.far.abs_sq();
        println!(
            "ν {nu:5.2}  a/m {:8.4}  b/(mω) {:8.4}  energy {:+.3}",
            s.added_mass / m,
            s.damping / (m * w),
            s.damping / far - 1.0
        );
    }
}

//! The interior lid against the plain Frank method: on a semicircle the two
//! must converge to the same coefficients; across a Wigley II midship
//! section's first irregular frequency (ν ≈ 7/m) only the lid stays smooth.
use michell_geometry::C64;
use michell_seakeeping::section2d::Section;
use std::f64::consts::PI;

fn main() {
    let (g, rho, r) = (9.81, 1000.0, 1.0);
    let m = rho * PI * r * r / 2.0;
    for nur in [0.3, 1.0, 2.0] {
        let w = (nur / r * g).sqrt();
        for n in [16, 64, 256] {
            let a = Section::semicircle(r, n).heave(w, g, rho).unwrap();
            let b = Section::semicircle(r, n)
                .without_lid()
                .heave(w, g, rho)
                .unwrap();
            println!(
                "semicircle νR {nur} N {n:3}: lid a/m {:.4} b/mω {:.4} | no lid a/m {:.4} b/mω {:.4}",
                a.added_mass / m,
                a.damping / (m * w),
                b.added_mass / m,
                b.damping / (m * w)
            );
        }
    }
    let (b, t) = (0.3, 0.1875);
    let curve: Vec<(f64, f64)> = (0..=128)
        .map(|i| {
            let z = t * i as f64 / 128.0;
            (b * (1.0 - (z / t).powi(10)), z)
        })
        .collect();
    let k = 2.79; // the incident wave of λ/L 0.75 on a 3 m hull
    for i in 0..13 {
        let nu = 6.4 + 0.1 * i as f64;
        let w = (nu * g).sqrt();
        let row = |sec: Section| {
            let s = sec.heave_with_diffraction(w, k, PI, g, rho).unwrap();
            (
                s.energy_error(),
                s.source_spectrum(nu, 1.0, C64::ONE, C64::ZERO).abs(),
                s.source_spectrum(nu, 1.0, C64::ZERO, C64::ONE).abs(),
            )
        };
        let lid = row(Section::from_curve(&curve, 20).unwrap());
        let bare = row(Section::from_curve(&curve, 20).unwrap().without_lid());
        println!(
            "Wigley II midship ν {nu:4.2}: lid err {:.3} |S_rad| {:.4} |S_diff| {:.4} | no lid err {:.3} |S_rad| {:.4} |S_diff| {:.4}",
            lid.0, lid.1, lid.2, bare.0, bare.1, bare.2
        );
    }
}

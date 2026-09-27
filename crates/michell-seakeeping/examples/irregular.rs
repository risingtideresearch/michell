//! Diagnostics: semicircle sections with differently graded panels, with
//! the energy check `b = ρω|C|²`.
use michell_seakeeping::section2d::Section;
use std::f64::consts::PI;

fn semi(r: f64, n: usize, grade: impl Fn(f64) -> f64) -> Section {
    Section::from_nodes(
        (0..=n)
            .map(|i| {
                let th = 0.5 * PI * grade(i as f64 / n as f64);
                [r * th.cos(), -r * th.sin()]
            })
            .collect(),
    )
}

fn main() {
    let (g, r) = (9.81, 0.15);
    let nu = 4.0f64;
    let w = (nu * g).sqrt();
    let cases: Vec<(&str, Section)> = vec![
        ("uniform", semi(r, 24, |t| t)),
        ("both ends", semi(r, 24, |t| 0.5 * (1.0 - (PI * t).cos()))),
        ("waterline", semi(r, 24, |t| 1.0 - (0.5 * PI * t).cos())),
        ("keel", semi(r, 24, |t| (0.5 * PI * t).sin())),
        ("mild both", semi(r, 24, |t| 0.5 * t + 0.25 * (1.0 - (PI * t).cos()))),
    ];
    let curve: Vec<(f64, f64)> = (0..=128)
        .map(|i| {
            let th = 0.5 * PI * i as f64 / 128.0;
            (r * th.cos(), r * th.sin())
        })
        .collect();
    let fc = Section::from_curve(&curve, 24).unwrap();
    let exact = semi(r, 24, |t| 0.5 * (1.0 - (PI * t).cos()));
    for (a, b) in fc.nodes().iter().zip(exact.nodes()) {
        println!("curve {:+.6} {:+.6}   exact {:+.6} {:+.6}", a[0], a[1], b[0], b[1]);
    }
    let cases: Vec<(&str, Section)> = cases.into_iter().chain([("from_curve", fc)]).collect();
    for (name, sec) in cases {
        let s = sec.heave(w, g, 1000.0).unwrap();
        println!("{name:10} a {:9.4} b {:9.4}  ρω|C|² {:9.4}", s.added_mass, s.damping, 1000.0 * w * s.far.abs_sq());
    }
}

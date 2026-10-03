//! Cost of strip-theory responses on a 3 m Wigley across wave lengths at
//! Fn 0.3 (the short waves meet high encounter frequencies).
use hullgeom::iges::{self, HullPose, Platform, SectionalOptions};
use seakeeping::section2d::Section;
use seakeeping::strip::{response, MassProperties, StripOptions, Wave};
use std::f64::consts::PI;

fn main() {
    let (l, g, rho) = (3.0, 9.81, 1000.0);
    let surfaces = iges::wigley_surfaces(l, 0.1 * l, 0.0625 * l).unwrap();
    let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
    let hull = source
        .situate_sectional(
            0,
            0.0,
            &HullPose::default(),
            &Platform::default(),
            &SectionalOptions::default(),
        )
        .unwrap()
        .unwrap()
        .hull;
    println!("{} stations", hull.stations());
    let mass = MassProperties::floating(&hull, rho, 0.25 * l);
    let opts = StripOptions {
        panels: 20,
        density: rho,
        gravity: g,
        roll_damping: 0.0,
    };
    let u = 0.3 * (g * l).sqrt();
    for lam in [1.2, 0.75, 0.5] {
        let k = 2.0 * PI / (lam * l);
        let t = std::time::Instant::now();
        let r = response(
            &hull,
            &mass,
            &Wave {
                omega: (k * g).sqrt(),
                heading: PI,
                speed: u,
            },
            &opts,
        )
        .unwrap();
        println!(
            "λ/L {lam}: ωe {:.1}, {:.0} ms",
            r.omega_e,
            t.elapsed().as_secs_f64() * 1e3
        );
    }
    // One midship section at high frequency: plain and bridged solves.
    let (_, c) = hull.curves().nth(hull.stations() / 2).unwrap();
    let sec = Section::from_curve(c, 20).unwrap();
    for we in [8.0, 16.0, 24.0] {
        let t = std::time::Instant::now();
        let s = sec.heave(we, g, rho).unwrap();
        println!(
            "section ωe {we}: {:.2} ms, interpolated {}, energy err {:.3}",
            t.elapsed().as_secs_f64() * 1e3,
            s.interpolated,
            s.energy_error()
        );
    }
}

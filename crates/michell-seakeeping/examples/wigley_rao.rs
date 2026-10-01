//! Heave and pitch RAOs of a 3 m parabolic Wigley in head seas, with the
//! exciting-force breakdown, for inspection.
use michell_geometry::iges::{self, HullPose, Platform, SectionalOptions};
use michell_seakeeping::strip::{response, MassProperties, StripOptions, Wave};
use std::f64::consts::PI;

fn main() {
    let (l, g, rho) = (3.0, 9.81, 1000.0);
    let surfaces = iges::wigley_surfaces(l, 0.1 * l, 0.0625 * l).unwrap();
    let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).unwrap();
    let so = SectionalOptions {
        stations: 41,
        ..Default::default()
    };
    let hull = source
        .situate_sectional(0, 0.0, &HullPose::default(), &Platform::default(), &so)
        .unwrap()
        .unwrap()
        .hull;
    let mass = MassProperties::floating(&hull, rho, 0.25 * l);
    let opts = StripOptions {
        panels: 16,
        density: rho,
        gravity: g,
        roll_damping: 0.0,
    };
    for fnum in [0.0, 0.2, 0.3] {
        let u = fnum * (g * l).sqrt();
        println!("Fn {fnum}");
        for i in 0..12 {
            let lam = l * (0.6 + 0.2 * i as f64);
            let k = 2.0 * PI / lam;
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
            let c = r.coefficients;
            let wn = (c.restoring[0][0] / (c.mass[0] + c.added_mass[0][0])).sqrt();
            let zeta = c.damping[0][0] / (2.0 * (c.mass[0] + c.added_mass[0][0]) * wn);
            println!(
                "  λ/L {:.2} ωe {:.2} ωn {:.2} ζ {:.2}  heave {:.3} pitch {:.3}  |FK3| {:.1} |FD3| {:.1} |F3| {:.1} C33 {:.1}",
                lam / l,
                r.omega_e,
                wn,
                zeta,
                r.heave_rao(),
                r.pitch_rao(),
                c.froude_krylov[0].abs(),
                c.diffraction[0].abs(),
                (c.froude_krylov[0] + c.diffraction[0]).abs(),
                c.restoring[0][0]
            );
        }
    }
}

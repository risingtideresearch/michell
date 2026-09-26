//! Resistance curve for the standard Wigley hull (L/B = 10, B/T = 1.6),
//! cut into sections from its exact CAD surfaces.
//!
//! Run with: cargo run --release --example wigley_curve

use michell::iges::{self, HullPose, Platform, SectionalOptions};
use michell::{sectional, Conditions, STANDARD_GRAVITY};

fn main() {
    let (l, b, t) = (10.0, 1.0, 0.625);
    let surfaces = iges::wigley_surfaces(l, b, t).expect("valid hull");
    let source = iges::source_fleet_from_surfaces(surfaces.to_vec(), 1.0, 0.0).expect("one hull");
    let cut = source
        .situate_sectional(
            0,
            0.0,
            &HullPose::default(),
            &Platform::default(),
            &SectionalOptions::default(),
        )
        .expect("sections")
        .expect("wet");
    let hull = &cut.hull;
    println!(
        "Wigley hull: L = {l} m, B = {b} m, T = {t} m, S = {:.3} m^2, V = {:.3} m^3\n",
        hull.wetted_surface(),
        hull.displaced_volume()
    );
    println!(
        "{:>5} {:>7} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "Fn", "U[m/s]", "Rw[N]", "Rv[N]", "Rt[N]", "Cw*1e3", "Ct*1e3"
    );
    let members = [(hull, cut.placement)];
    let mut fr = 0.10;
    while fr <= 0.501 {
        let u = fr * (STANDARD_GRAVITY * l).sqrt();
        let cond = Conditions::freshwater(u);
        let r = sectional::multihull_resistance(
            &members,
            &cond,
            &Default::default(),
            &Default::default(),
        )
        .expect("resistance");
        println!(
            "{fr:>5.2} {u:>7.3} {:>10.2} {:>10.2} {:>10.2} {:>10.4} {:>10.4}",
            r.wave.resistance,
            r.viscous_total,
            r.total,
            r.cw * 1e3,
            r.ct * 1e3
        );
        fr += 0.02;
    }
}

//! Resistance curve for the standard Wigley hull (L/B = 10, B/T = 1.6).
//!
//! Run with: cargo run --release --example wigley_curve

use michell::{hulls, Conditions, STANDARD_GRAVITY};

fn main() {
    let (l, b, t) = (10.0, 1.0, 0.625);
    let hull = hulls::wigley(l, b, t).expect("valid hull");
    println!(
        "Wigley hull: L = {l} m, B = {b} m, T = {t} m, S = {:.3} m^2, V = {:.3} m^3\n",
        hull.wetted_surface(),
        hull.displaced_volume()
    );
    println!(
        "{:>5} {:>7} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "Fn", "U[m/s]", "Rw[N]", "Rv[N]", "Rt[N]", "Cw*1e3", "Ct*1e3"
    );
    let mut fr = 0.10;
    while fr <= 0.501 {
        let u = fr * (STANDARD_GRAVITY * l).sqrt();
        let cond = Conditions::freshwater(u);
        let r = michell::resistance(&hull, &cond).expect("resistance");
        println!(
            "{fr:>5.2} {u:>7.3} {:>10.2} {:>10.2} {:>10.2} {:>10.4} {:>10.4}",
            r.wave.resistance,
            r.viscous.resistance,
            r.total,
            r.cw * 1e3,
            r.ct * 1e3
        );
        fr += 0.02;
    }
}

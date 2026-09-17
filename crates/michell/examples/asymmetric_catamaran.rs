//! Asymmetric-demihull catamaran — the motivating case for solved-dipole
//! multihull (Kaklis & Papanikolaou). Two mirror-image demihulls, each
//! asymmetric about its own centreplane, carry solved centreplane-lifting
//! dipoles whose wave systems interfere across the catamaran.
//!
//! Compares, across speed, the symmetric-mean catamaran with the
//! asymmetric-demihull catamaran (both at the same demihull separation).
//!
//! Run with: cargo run --release --example asymmetric_catamaran

use michell::{
    multihull_wave_resistance_lifting, BSplineSurface, Conditions, Hull, LiftingGrid, Placement,
    STANDARD_GRAVITY,
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
    let base = michell::hulls::wigley(l, 1.0, 0.625).unwrap();
    // Each demihull: flat-ish inboard (port 0.5), full outboard (stbd 1.0).
    let stb = scaled(base.surface(), 1.0);
    let port = scaled(base.surface(), 0.5);
    let h = Hull::new_asymmetric(scaled(&port, 1.0), scaled(&stb, 1.0)).unwrap();
    let h_m = Hull::new_asymmetric(scaled(&stb, 1.0), scaled(&port, 1.0)).unwrap();
    let mean = Hull::new(scaled(h.surface(), 1.0)).unwrap();

    let sep = 3.0; // demihull centreplane separation [m]
    let grid = LiftingGrid { nx: 48, nz: 16 };
    let opts = Default::default();

    println!("Asymmetric-demihull catamaran vs symmetric-mean catamaran");
    println!("(demihull separation {sep} m, camber facing outboard)\n");
    println!(
        "{:>5} {:>7} {:>12} {:>12} {:>9}",
        "Fn", "U[m/s]", "R_sym [N]", "R_asym [N]", "Δ%"
    );
    let mut fr = 0.20;
    while fr <= 0.451 {
        let u = fr * (STANDARD_GRAVITY * l).sqrt();
        let c = Conditions::seawater(u);
        let sym = [
            (
                &mean,
                Placement {
                    x: 0.0,
                    y: sep / 2.0,
                },
            ),
            (
                &mean,
                Placement {
                    x: 0.0,
                    y: -sep / 2.0,
                },
            ),
        ];
        let asym = [
            (
                &h,
                Placement {
                    x: 0.0,
                    y: sep / 2.0,
                },
            ),
            (
                &h_m,
                Placement {
                    x: 0.0,
                    y: -sep / 2.0,
                },
            ),
        ];
        let rs = multihull_wave_resistance_lifting(&sym, &c, &opts, grid)
            .unwrap()
            .resistance;
        let ra = multihull_wave_resistance_lifting(&asym, &c, &opts, grid)
            .unwrap()
            .resistance;
        println!(
            "{fr:>5.2} {u:>7.3} {rs:>12.2} {ra:>12.2} {:>8.1}%",
            100.0 * (ra - rs) / rs
        );
        fr += 0.05;
    }
    println!(
        "\nR_sym is the thickness-only (source) catamaran; R_asym adds the solved\n\
         centreplane-dipole systems of the two asymmetric demihulls, which\n\
         interfere across the catamaran. Demihull asymmetry is a design knob for\n\
         tuning that interference."
    );
}

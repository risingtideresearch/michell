//! Heeled-hull wave resistance — the tilted-centreplane (complex-κ) kernel.
//!
//! Tabulates the added wave resistance of a Wigley hull versus heel angle at a
//! few Froude numbers. Heel keeps the sources on the ship's tilted centreplane,
//! turning the vertical decay complex; the upright kernel is untouched, so 0°
//! reproduces the standard Michell result exactly.
//!
//! Run with: cargo run --release --example heel

use michell::{heel_wave_resistance, hulls, Conditions, WaveOptions, STANDARD_GRAVITY};

const DEG: f64 = std::f64::consts::PI / 180.0;

fn main() {
    let l = 10.0;
    let hull = hulls::wigley(l, 1.0, 0.625).unwrap();
    let opts = WaveOptions::default();
    let angles = [0.0, 5.0, 10.0, 15.0, 20.0, 30.0];

    println!("Wigley hull — wave resistance [N] vs heel angle (added % over upright)\n");
    print!("{:>5} {:>7}", "Fn", "U[m/s]");
    for a in angles {
        print!("{:>13}", format!("{a:.0}deg"));
    }
    println!();

    for fr in [0.25, 0.32, 0.40] {
        let u = fr * (STANDARD_GRAVITY * l).sqrt();
        let cond = Conditions::seawater(u);
        let r0 = heel_wave_resistance(&hull, &cond, 0.0, &opts)
            .unwrap()
            .resistance;
        print!("{fr:>5.2} {u:>7.3}");
        for a in angles {
            let r = heel_wave_resistance(&hull, &cond, a * DEG, &opts)
                .unwrap()
                .resistance;
            let pct = 100.0 * (r - r0) / r0;
            print!("{:>13}", format!("{r:.1} (+{pct:.1}%)"));
        }
        println!();
    }
    println!(
        "\nHeel raises wave resistance: the tilt reduces the effective depth-decay\n\
         (e^(−νλ²z·cosφ)) and radiates a port/starboard-asymmetric wave pattern.\n\
         Models the tilted thickness distribution's wave-making; the lifting\n\
         side-force of a heeled-and-yawed hull is a separate effect."
    );
}

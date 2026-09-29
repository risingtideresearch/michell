//! Vugts' (1970) horizontal cylinders in beam waves: this crate's 2-D sway,
//! heave and roll coefficients and wave loads, non-dimensionalised as in
//! Journée's SEAWAY validation report (1213): `ω' = ω√(B/2g)`,
//! `a'₂₂ = a₂₂/(ρA)`, `a'₄₄ = a₄₄/(ρAB²)`, `a'₂₄ = a₂₄/(ρAB)` (damping the
//! same times `√(B/2g)`), `F'₂ = F₂/(ρgkA)`, `F'₃ = F₃/(ρgB)`,
//! `F'₄ = F₄/(ρgkB³/12)`, roll about G at `OG` above the waterline. Prints
//! one tab-separated row per case and frequency; `python/tools/vugts/`
//! digitises the report's figures and compares.
//!
//! Against the report's own conventions: its heave-force phase is this
//! one's; its sway-force phase is this one's plus 180° (the opposite sign of
//! sway), and its roll-moment phase is 90° minus this one's (taken against
//! the wave slope rather than the elevation). For the two sections with G
//! above the waterline (B/d = 4, 8) its plotted sway–roll couplings match
//! this crate's about O, with the sign reversed, rather than about G —
//! while its roll moment and roll inertia there match this crate's about G.
use michell_seakeeping::section2d::Section;
use std::f64::consts::PI;

fn main() {
    let (g, rho, b) = (9.81, 1000.0, 1.0);
    // (name, section, area, OG)
    let cases: Vec<(&str, Section, f64, f64)> = vec![
        (
            "circle",
            Section::semicircle(0.5 * b, 64),
            PI * 0.125 * b * b,
            0.0,
        ),
        (
            "rect2",
            Section::rectangle(0.5 * b, 0.5 * b, 64),
            0.5 * b * b,
            0.0,
        ),
        (
            "rect4",
            Section::rectangle(0.5 * b, 0.25 * b, 64),
            0.25 * b * b,
            0.25 * b,
        ),
        (
            "rect8",
            Section::rectangle(0.5 * b, 0.125 * b, 64),
            0.125 * b * b,
            3.0 * 0.125 * b,
        ),
    ];
    println!(
        "case\tw'\ta22\tb22\ta24\tb24\ta42\tb42\ta44\tb44\ta33\tb33\tF2\tph2\tF3\tph3\tF4\tph4"
    );
    for (name, sec, area, og) in &cases {
        let (area, og) = (*area, *og);
        for i in 1..=60 {
            let wn = 0.025 * i as f64;
            let w = wn / (b / (2.0 * g)).sqrt();
            let k = w * w / g;
            let heading = 0.5 * PI;
            let (Some(lat), Some(hv)) = (
                sec.lateral(w, g, rho, Some((k, heading))),
                sec.heave_with_diffraction(w, k, heading, g, rho),
            ) else {
                continue;
            };
            let (a, bb) = (lat.added_mass, lat.damping);
            // About G at height og: ψ₄ᴳ = ψ₄ + og ψ₂.
            let cog = og;
            let a24 = a[0][1] + cog * a[0][0];
            let a42 = a[1][0] + cog * a[0][0];
            let a44 = a[1][1] + og * (a[0][1] + a[1][0]) + og * og * a[0][0];
            let b24 = bb[0][1] + cog * bb[0][0];
            let b42 = bb[1][0] + cog * bb[0][0];
            let b44 = bb[1][1] + og * (bb[0][1] + bb[1][0]) + og * og * bb[0][0];
            let s = (b / (2.0 * g)).sqrt();
            let fl = lat.froude_krylov(k, heading, g, rho);
            let dl = lat.diffraction(k, heading, g, rho);
            let f2 = fl[0] + dl[0];
            let f4 = fl[1] + dl[1] + f2.scale(og);
            let f3 = hv.froude_krylov(k, heading, g, rho) + hv.diffraction(k, heading, g, rho);
            // Phase ε of F = Fa cos(ωt + ε) against ζ = cos(ωt) at the centre:
            // with e^{−iωt} amplitudes, ε = −arg F.
            let ph = |f: michell_geometry::C64| -(f.im.atan2(f.re)).to_degrees();
            println!(
                "{name}\t{wn:.3}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.1}\t{:.4}\t{:.1}\t{:.4}\t{:.1}",
                a[0][0] / (rho * area),
                bb[0][0] / (rho * area) * s,
                a24 / (rho * area * b),
                b24 / (rho * area * b) * s,
                a42 / (rho * area * b),
                b42 / (rho * area * b) * s,
                a44 / (rho * area * b * b),
                b44 / (rho * area * b * b) * s,
                hv.added_mass / (rho * area),
                hv.damping / (rho * area) * s,
                f2.abs() / (rho * g * k * area),
                ph(f2),
                f3.abs() / (rho * g * b),
                ph(f3),
                f4.abs() / (rho * g * k * b * b * b / 12.0),
                ph(f4)
            );
        }
    }
}

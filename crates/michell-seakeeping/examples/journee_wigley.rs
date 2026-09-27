//! Comparison with Journée's experiments on four Wigley hulls in head waves
//! (DUT-SHL Report 0909, 1992): forced-heave coefficients (Table 2) and
//! heave, pitch and added resistance in regular waves (Table 8). Pitch is
//! compared in the report's normalisation, `θ_a / (2π ζ_a / L)`.
//!
//! The report's ASCII data files (`0909-DUT-92.zip`, freely distributed by
//! the author; not in this repository) are read from the directory given as
//! the first argument:
//!
//! ```text
//! cargo run --release -p michell-seakeeping --example journee_wigley -- /path/to/wigley
//! ```
//!
//! Hulls (report §2): `η = (1 − ζ²)(1 − ξ²)(1 + 0.2ξ²) + α ζ²(1 − ζ⁸)(1 − ξ²)⁴`,
//! `ξ = 2x/L`, `η = 2y/B`, `ζ = z/d`, α = 1 for I and II (C_m 0.909), 0 for
//! III and IV (C_m 0.667); L = 3 m, d = 0.1875 m, B = 0.3 m (I, III) or
//! 0.6 m (II, IV), k_yy = 0.75 m, fresh water.

use michell_geometry::sectional::{DepthQuadrature, SectionNodes, SectionalHull};
use michell_seakeeping::strip::{
    added_resistance, added_resistance_maruo, MassProperties, StripOptions, Wave,
};
use std::f64::consts::PI;

const L: f64 = 3.0;
const D: f64 = 0.1875;
const RHO: f64 = 1000.0;
const G: f64 = 9.81;

/// Wigley `model` (1–4) as a sectional hull: its exact half-breadth,
/// sampled at the Greville stations of a clamped cubic B-spline in x.
fn wigley(model: usize) -> (SectionalHull, f64) {
    let b = if model % 2 == 1 { 0.3 } else { 0.6 };
    let alpha = if model <= 2 { 1.0 } else { 0.0 };
    let half = move |x: f64, z: f64| {
        let xi = 2.0 * x / L;
        let ze = z / D;
        let eta = (1.0 - ze * ze) * (1.0 - xi * xi) * (1.0 + 0.2 * xi * xi)
            + alpha * ze * ze * (1.0 - ze.powi(8)) * (1.0 - xi * xi).powi(4);
        0.5 * b * eta.max(0.0)
    };
    let (p, spans) = (3usize, 60usize);
    let mut knots = vec![-0.5 * L; p + 1];
    for i in 1..spans {
        knots.push(-0.5 * L + L * i as f64 / spans as f64);
    }
    knots.extend(vec![0.5 * L; p + 1]);
    let n = knots.len() - p - 1;
    let xs: Vec<f64> = (0..n)
        .map(|i| knots[i + 1..=i + p].iter().sum::<f64>() / p as f64)
        .collect();
    let q = DepthQuadrature::default();
    let sections = xs
        .iter()
        .map(|&x| SectionNodes::from_depth_function(D, &[], |z| half(x, z), &q))
        .collect();
    (SectionalHull::new(p, knots, &xs, sections).unwrap(), b)
}

/// Rows of a data file; `-1111` marks a missing value.
fn read(dir: &str, name: &str) -> Option<Vec<Vec<Option<f64>>>> {
    let text = std::fs::read_to_string(format!("{dir}/{name}")).ok()?;
    Some(
        text.lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                l.split_whitespace()
                    .map(|v| v.parse::<f64>().ok().filter(|&x| x > -1000.0))
                    .collect()
            })
            .collect(),
    )
}

fn main() {
    let dir = std::env::args()
        .nth(1)
        .expect("usage: journee_wigley <dir with the .exp files>");
    let opts = StripOptions {
        panels: 20,
        density: RHO,
        gravity: G,
    };
    let names = ["I", "II", "III", "IV"];
    for model in 1..=4 {
        let (hull, b) = wigley(model);
        let vol = hull.displaced_volume();
        let mass = MassProperties {
            mass: RHO * vol,
            lcg: 0.0,
            radius_of_gyration: 0.75,
            bg: 0.0,
            roll_radius_of_gyration: 0.35 * b,
            yaw_radius_of_gyration: 0.75,
        };
        let name = names[model - 1];
        println!("\n=== Wigley {name}: B {b} m, ∇ {vol:.4} m³");
        // Table 10: zero forward speed (III and IV): no speed terms at all.
        if let Some(rows) = read(&dir, &format!("10-{name}-00.exp")) {
            println!("Fn 0    head waves:  λ/L   heave exp/calc   pitch exp/calc   Raw'' exp / GB / Maruo");
            for r in rows {
                let (Some(lam), Some(z), Some(th)) = (r[0], r[4], r[6]) else {
                    continue;
                };
                let k = 2.0 * PI / (lam * L);
                let wave = Wave {
                    omega: (k * G).sqrt(),
                    heading: PI,
                    speed: 0.0,
                };
                let r_c = added_resistance(&hull, &mass, &wave, &opts).unwrap();
                let raw = r
                    .get(8)
                    .copied()
                    .flatten()
                    .map_or("  -  ".into(), |v| format!("{v:5.2}"));
                let maruo = added_resistance_maruo(&hull, &mass, &wave, &opts).unwrap();
                println!(
                    "            {lam:5.3}   {z:5.2} / {:5.2}    {th:5.2} / {:5.2}    {raw} / {:6.2} / {:6.2}",
                    r_c.response.heave_rao(),
                    r_c.response.pitch_rao() / lam,
                    r_c.coefficient(RHO, G, b, L),
                    maruo.coefficient(RHO, G, b, L)
                );
            }
        }
        for fnum in [0.2, 0.3, 0.4] {
            let u = fnum * (G * L).sqrt();
            let tag = format!("{:02}", (fnum * 100.0f64).round() as usize);
            // Table 2: forced heave, A33'' = A33/(ρ∇), B33'' = B33/(ρ∇)·√(L/g).
            if let Some(rows) = read(&dir, &format!("02-{name}-{tag}.exp")) {
                println!("Fn {fnum}  forced heave:   ω'   A33'' exp/calc   B33'' exp/calc");
                for r in rows {
                    let (Some(wn), Some(a), Some(bb)) = (r[0], r[1], r[2]) else {
                        continue;
                    };
                    let we = wn / (L / G).sqrt();
                    // Heave-only coefficients at ω_e: from the strip solver
                    // with any wave of that encounter frequency; zero-pitch
                    // A33, B33 are wave-independent. Pick the head-sea wave.
                    let w0 = solve_w0(we, u);
                    let resp = added_resistance(
                        &hull,
                        &mass,
                        &Wave {
                            omega: w0,
                            heading: PI,
                            speed: u,
                        },
                        &opts,
                    )
                    .unwrap()
                    .response;
                    let c = resp.coefficients;
                    let a_c = c.added_mass[0][0] / (RHO * vol);
                    let b_c = c.damping[0][0] / (RHO * vol) * (L / G).sqrt();
                    println!("            {wn:5.2}   {a:6.3} / {a_c:6.3}    {bb:6.3} / {b_c:6.3}");
                }
            }
            // Table 3: forced pitch, A55'' = A55/(ρ∇L²), B55'' = B55/(ρ∇L²)·√(L/g),
            // A35'' = A35/(ρ∇L), B35'' = B35/(ρ∇L)·√(L/g).
            if let Some(rows) = read(&dir, &format!("03-{name}-{tag}.exp")) {
                println!("Fn {fnum}  forced pitch:   ω'   A55'' exp/calc    B55'' exp/calc    A35'' exp/calc    B35'' exp/calc");
                for r in rows {
                    let (Some(wn), Some(a55), Some(b55), Some(a35), Some(b35)) =
                        (r[0], r[1], r[2], r[3], r[4])
                    else {
                        continue;
                    };
                    let we = wn / (L / G).sqrt();
                    let w0 = solve_w0(we, u);
                    let c = added_resistance(
                        &hull,
                        &mass,
                        &Wave {
                            omega: w0,
                            heading: PI,
                            speed: u,
                        },
                        &opts,
                    )
                    .unwrap()
                    .response
                    .coefficients;
                    let s = (L / G).sqrt();
                    println!(
                        "            {wn:5.2}   {a55:7.4} / {:7.4}   {b55:7.4} / {:7.4}   {a35:7.4} / {:7.4}   {b35:7.4} / {:7.4}",
                        c.added_mass[1][1] / (RHO * vol * L * L),
                        c.damping[1][1] / (RHO * vol * L * L) * s,
                        c.added_mass[0][1] / (RHO * vol * L),
                        c.damping[0][1] / (RHO * vol * L) * s,
                    );
                }
            }
            // Table 6: wave loads on the restrained model, X3'' = |X3|/(C33 ζ),
            // X5'' = |X5|/(k C55 ζ).
            if let Some(rows) = read(&dir, &format!("06-{name}-{tag}.exp")) {
                println!("Fn {fnum}  wave loads:  λ/L   X3'' exp/calc   X5'' exp/calc");
                for r in rows {
                    let Some(lam) = r[0] else { continue };
                    let Some(g5) = r[1..].chunks(5).find(|g| g.len() == 5 && g[1].is_some()) else {
                        continue;
                    };
                    let k = 2.0 * PI / (lam * L);
                    let wave = Wave {
                        omega: (k * G).sqrt(),
                        heading: PI,
                        speed: u,
                    };
                    let c = added_resistance(&hull, &mass, &wave, &opts)
                        .unwrap()
                        .response
                        .coefficients;
                    let f3 = (c.froude_krylov[0] + c.diffraction[0]).abs() / c.restoring[0][0];
                    let f5 =
                        (c.froude_krylov[1] + c.diffraction[1]).abs() / (k * c.restoring[1][1]);
                    let f = |v: Option<f64>| v.map_or("  -  ".into(), |x| format!("{x:5.3}"));
                    println!(
                        "            {lam:5.3}   {} / {f3:5.3}    {} / {f5:5.3}",
                        f(g5[1]),
                        f(g5[3])
                    );
                }
            }
            // Table 8: motions and added resistance in regular head waves.
            if let Some(rows) = read(&dir, &format!("08-{name}-{tag}.exp")) {
                println!("Fn {fnum}  head waves:  λ/L   heave exp/calc   pitch exp/calc   Raw'' exp / GB / Maruo");
                for r in rows {
                    let Some(lam) = r[0] else { continue };
                    // First amplitude group with data.
                    let group = r[1..].chunks(6).find(|g| g.len() == 6 && g[1].is_some());
                    let Some(g6) = group else { continue };
                    let k = 2.0 * PI / (lam * L);
                    let wave = Wave {
                        omega: (k * G).sqrt(),
                        heading: PI,
                        speed: u,
                    };
                    let r_c = added_resistance(&hull, &mass, &wave, &opts).unwrap();
                    let raw_c = r_c.coefficient(RHO, G, b, L);
                    let raw_m = added_resistance_maruo(&hull, &mass, &wave, &opts)
                        .unwrap()
                        .coefficient(RHO, G, b, L);
                    let f = |v: Option<f64>| v.map_or("  -  ".into(), |x| format!("{x:5.2}"));
                    println!(
                        "            {lam:5.3}   {} / {:5.2}    {} / {:5.2}    {} / {:6.2} / {:6.2}",
                        f(g6[1]),
                        r_c.response.heave_rao(),
                        f(g6[3]),
                        r_c.response.pitch_rao() / lam,
                        f(g6[5]),
                        raw_c,
                        raw_m
                    );
                }
            }
        }
    }
}

/// The absolute frequency whose head-sea encounter frequency at `u` is `we`:
/// `we = w0 + w0² U/g`.
fn solve_w0(we: f64, u: f64) -> f64 {
    if u == 0.0 {
        return we;
    }
    let a = u / G;
    (-1.0 + (1.0 + 4.0 * a * we).sqrt()) / (2.0 * a)
}

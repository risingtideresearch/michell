//! `michell seakeeping`: heave, pitch and added resistance in regular
//! waves by strip theory (`michell-seakeeping`), and optionally their
//! statistics in an irregular sea.

use super::{fleet, parse_args, parse_range, KNOT};
use michell::sectional::dynamic_load_closure;
use michell::squat::SquatOptions;
use michell_geometry::float::{solve_equilibrium_sectional_dynamic, FleetState, LoadCase};
use michell_geometry::iges::HullPose;
use michell_geometry::source::SourceHull;
use michell_geometry::{Placement, SectionalHull};
use michell_seakeeping::sea::{sea_response_fleet, Spectrum};
use michell_seakeeping::strip::{added_resistance_both, MassProperties, StripOptions, Wave};
use std::f64::consts::PI;

pub(crate) const USAGE: &str =
    "usage: michell seakeeping <hull>[@x=DX,y=Y]... (--speed U | --froude F) \
[--heading DEG] [--lambda A:B:STEP] [--kyy FRAC] [--mass KG] [--lcg X] [--panels N] \
[--sea hs=H,tp=T[,gamma=G]] [--dynamic] [--csv]";

pub(crate) fn cmd_seakeeping(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err(USAGE.into());
    }
    let csv = p.switch("--csv");
    // In --csv mode stdout carries only the table; notes go to stderr.
    let note = |line: String| {
        if csv {
            eprintln!("{line}");
        } else {
            println!("{line}");
        }
    };
    let settings = p.load_settings()?;
    let fleet = fleet::load(&p.positional, &settings)?;
    let static_members = fleet.hulls();
    let members = &static_members;
    let cond = p.conditions(0.0)?;
    let (rho, g) = (cond.fluid.density, cond.gravity);
    let l_ref = members
        .iter()
        .map(|(h, _)| h.length())
        .fold(0.0f64, f64::max);
    let speed = match (p.f64_flag("speed")?, p.f64_flag("froude")?) {
        (Some(_), Some(_)) => return Err("give either --speed or --froude, not both".into()),
        (Some(u), None) => {
            if p.switch("--knots") {
                u * KNOT
            } else {
                u
            }
        }
        (None, Some(f)) => f * (g * l_ref).sqrt(),
        (None, None) => return Err(format!("select a speed with --speed or --froude\n{USAGE}")),
    };
    let heading = p.f64_flag("heading")?.unwrap_or(180.0).to_radians();
    let lambdas = match p.flag("lambda") {
        Some(s) => parse_range(s)?,
        None => parse_range("0.5:3:0.125")?,
    };
    let volume: f64 = members.iter().map(|(h, _)| h.displaced_volume()).sum();
    if !(volume > 0.0) {
        return Err("the hulls are dry at this waterline".into());
    }
    let lcb = members
        .iter()
        .map(|(h, pl)| h.displaced_volume() * (h.lcb_x() + pl.x))
        .sum::<f64>()
        / volume;
    let mass = MassProperties {
        mass: p.f64_flag("mass")?.unwrap_or(rho * volume),
        lcg: p.f64_flag("lcg")?.unwrap_or(lcb),
        radius_of_gyration: p.f64_flag("kyy")?.unwrap_or(0.25) * l_ref,
        bg: 0.0,
    };
    // With --dynamic, float the platform at its dynamic attitude at this
    // speed (thin-ship sinkage and trim, `michell::squat`) and take the
    // motions about that attitude rather than the loaded waterline.
    let dynamic_state: Option<FleetState<SectionalHull>> = if p.switch("--dynamic") {
        let posed: Vec<SourceHull> = fleet
            .members
            .iter()
            .map(|m| SourceHull {
                source: fleet.source(m),
                index: m.index,
                waterline_z: fleet.files[m.file].waterline_z,
                pose: HullPose {
                    dx: m.shift.x,
                    dy: m.shift.y,
                    ..HullPose::default()
                },
            })
            .collect();
        let cond = p.conditions(speed)?;
        let squat = SquatOptions::default();
        let eq = solve_equilibrium_sectional_dynamic(
            &posed,
            &LoadCase {
                mass: mass.mass,
                lcg: Some(mass.lcg),
            },
            rho,
            g,
            &settings.sectional(0.0),
            dynamic_load_closure(&cond, mass.lcg, &squat),
            None,
        )
        .map_err(|e| format!("dynamic equilibrium: {e}"))?;
        note(format!(
            "dynamic attitude: sinkage {:.1} mm, trim {:.3}° bow up (dynamic lift {:.1}% of weight)",
            1000.0 * eq.sinkage,
            eq.trim.to_degrees(),
            100.0 * eq.lift_fraction
        ));
        Some(eq.fleet)
    } else {
        None
    };
    let dynamic_members: Vec<(&SectionalHull, Placement)> = dynamic_state
        .as_ref()
        .map(|st| st.members.iter().map(|(h, pl)| (h, *pl)).collect())
        .unwrap_or_default();
    let members = if dynamic_state.is_some() {
        &dynamic_members
    } else {
        members
    };
    let opts = StripOptions {
        panels: p.f64_flag("panels")?.map_or(20, |n| n.max(4.0) as usize),
        density: rho,
        gravity: g,
    };
    let beam = members
        .iter()
        .map(|(h, _)| {
            let (a, b) = h.x_range();
            (0..=200)
                .map(|i| 2.0 * h.waterline_half_beam(a + (b - a) * i as f64 / 200.0))
                .fold(0.0f64, f64::max)
        })
        .fold(0.0f64, f64::max);
    note(format!(
        "platform: {} hull(s), L {:.3} m, B(hull) {:.3} m, mass {:.1} kg, LCG {:.3} m, k_yy {:.3} m; \
         U {:.3} m/s (Fn {:.3}), heading {:.0}°",
        members.len(),
        l_ref,
        beam,
        mass.mass,
        mass.lcg,
        mass.radius_of_gyration,
        speed,
        speed / (g * l_ref).sqrt(),
        heading.to_degrees()
    ));
    if csv {
        println!("lambda_over_L,omega,omega_e,heave,heave_phase_deg,pitch_over_k,pitch_phase_deg,raw_gb_per_zeta2,sigma_aw_gb,raw_maruo_per_zeta2,sigma_aw_maruo");
    } else {
        println!(
            "{:>7} {:>8} {:>8} {:>8} {:>9} {:>8} {:>9} {:>11} {:>8} {:>11} {:>8}",
            "λ/L",
            "ω[r/s]",
            "ωe[r/s]",
            "heave",
            "ph3[deg]",
            "pitch",
            "ph5[deg]",
            "Raw/ζ² GB",
            "σ_aw GB",
            "Raw/ζ² far",
            "σ_aw far"
        );
    }
    for &lam in &lambdas {
        let k = 2.0 * PI / (lam * l_ref);
        let wave = Wave {
            omega: (k * g).sqrt(),
            heading,
            speed,
        };
        match added_resistance_both(&members, &mass, &wave, &opts) {
            Ok((resp, gb, far)) => {
                let sigma = |raw: f64| raw / (rho * g * beam * beam / l_ref);
                let row = [
                    lam,
                    wave.omega,
                    resp.omega_e,
                    resp.heave_rao(),
                    resp.heave.im.atan2(resp.heave.re).to_degrees(),
                    resp.pitch_rao(),
                    resp.pitch.im.atan2(resp.pitch.re).to_degrees(),
                    gb,
                    sigma(gb),
                    far,
                    sigma(far),
                ];
                if csv {
                    let cells: Vec<String> = row.iter().map(|v| format!("{v:.6}")).collect();
                    println!("{}", cells.join(","));
                } else {
                    println!(
                        "{:7.3} {:8.3} {:8.3} {:8.3} {:9.1} {:8.3} {:9.1} {:11.2} {:8.3} {:11.2} {:8.3}",
                        row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10]
                    );
                }
            }
            Err(e) => note(format!("{lam:7.3} {:8.3}  — {e}", wave.omega)),
        }
    }
    note(
        "heave per unit wave amplitude; pitch as |η5|/(kζ); phases relative to a crest at the LCG; \
         added resistance by radiated energy (GB, Gerritsma–Beukelman) and far-field momentum \
         (far, Maruo) — on Journée's Wigley hulls GB is nearer the tank at Fn 0.2, the far field \
         at Fn 0.3–0.4; neither ranks hulls reliably"
            .into(),
    );
    if let Some(sea) = p.flag("sea") {
        let spectrum = parse_sea(sea)?;
        let (_, x1) = members
            .iter()
            .map(|(h, pl)| (h.x_range().0 + pl.x, h.x_range().1 + pl.x))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |a, r| {
                (a.0.min(r.0), a.1.max(r.1))
            });
        let stations = [x1, mass.lcg];
        let r = sea_response_fleet(
            &members, &mass, &spectrum, heading, speed, &stations, 41, &opts,
        )
        .map_err(|e| format!("{e}"))?;
        note(format!("irregular sea {spectrum:?}:"));
        note(format!("  significant heave amplitude   {:.3} m", r.heave));
        note(format!(
            "  significant pitch amplitude   {:.2}°",
            r.pitch.to_degrees()
        ));
        note(format!(
            "  significant vertical accel.   {:.2} m/s² at the bow, {:.2} m/s² at the LCG",
            r.accelerations[0], r.accelerations[1]
        ));
        note(format!(
            "  mean added resistance         {:.1} N (GB), {:.1} N (far field)",
            r.added_resistance, r.added_resistance_far_field
        ));
        if r.skipped_energy > 0.0 {
            note(format!(
                "  ({:.0}% of the sea's energy was overtaken and left out)",
                100.0 * r.skipped_energy
            ));
        }
    }
    Ok(())
}

/// `hs=H,tp=T[,gamma=G]`: Bretschneider, or JONSWAP when `gamma` is given.
fn parse_sea(s: &str) -> Result<Spectrum, String> {
    let (mut hs, mut tp, mut gamma) = (None, None, None);
    for part in s.split(',') {
        let (k, v) = part
            .split_once('=')
            .ok_or_else(|| format!("--sea: expected key=value, got {part:?}"))?;
        let v: f64 = v
            .trim()
            .parse()
            .map_err(|_| format!("--sea: cannot parse {v:?}"))?;
        match k.trim() {
            "hs" => hs = Some(v),
            "tp" => tp = Some(v),
            "gamma" => gamma = Some(v),
            other => return Err(format!("--sea: unknown key {other:?} (hs, tp, gamma)")),
        }
    }
    let (hs, tp) = (hs.ok_or("--sea needs hs=")?, tp.ok_or("--sea needs tp=")?);
    Ok(match gamma {
        Some(gamma) => Spectrum::Jonswap { hs, tp, gamma },
        None => Spectrum::Bretschneider { hs, tp },
    })
}

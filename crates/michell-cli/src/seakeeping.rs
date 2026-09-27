//! `michell seakeeping`: heave, pitch and added resistance in regular
//! waves by strip theory (`michell-seakeeping`), and optionally their
//! statistics in an irregular sea.

use super::{fleet, parse_args, parse_range, KNOT};
use michell_seakeeping::sea::{sea_response_fleet, Spectrum};
use michell_seakeeping::strip::{added_resistance_fleet, MassProperties, StripOptions, Wave};
use std::f64::consts::PI;

pub(crate) const USAGE: &str = "usage: michell seakeeping <hull>[@x=DX,y=Y]... (--speed U | --froude F) \
[--heading DEG] [--lambda A:B:STEP] [--kyy FRAC] [--mass KG] [--lcg X] [--panels N] \
[--sea hs=H,tp=T[,gamma=G]]";

pub(crate) fn cmd_seakeeping(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err(USAGE.into());
    }
    let fleet = fleet::load(&p.positional, &p.load_settings()?)?;
    let members = fleet.hulls();
    let cond = p.conditions(0.0)?;
    let (rho, g) = (cond.fluid.density, cond.gravity);
    let l_ref = members.iter().map(|(h, _)| h.length()).fold(0.0f64, f64::max);
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
    println!(
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
    );
    println!(
        "{:>7} {:>8} {:>8} {:>8} {:>9} {:>8} {:>9} {:>11} {:>9}",
        "λ/L", "ω[r/s]", "ωe[r/s]", "heave", "ph3[deg]", "pitch", "ph5[deg]", "Raw/ζ²[N/m²]", "σ_aw"
    );
    for &lam in &lambdas {
        let k = 2.0 * PI / (lam * l_ref);
        let wave = Wave {
            omega: (k * g).sqrt(),
            heading,
            speed,
        };
        match added_resistance_fleet(&members, &mass, &wave, &opts) {
            Ok(r) => {
                let resp = r.response;
                println!(
                    "{:7.3} {:8.3} {:8.3} {:8.3} {:9.1} {:8.3} {:9.1} {:11.2} {:9.3}",
                    lam,
                    wave.omega,
                    resp.omega_e,
                    resp.heave_rao(),
                    resp.heave.im.atan2(resp.heave.re).to_degrees(),
                    resp.pitch_rao(),
                    resp.pitch.im.atan2(resp.pitch.re).to_degrees(),
                    r.per_amplitude_sq,
                    r.coefficient(rho, g, beam, l_ref)
                );
            }
            Err(e) => println!("{lam:7.3} {:8.3}  — {e}", wave.omega),
        }
    }
    println!("heave per unit wave amplitude; pitch as |η5|/(kζ); phases relative to a crest at the LCG");
    if let Some(sea) = p.flag("sea") {
        let spectrum = parse_sea(sea)?;
        let (_, x1) = members
            .iter()
            .map(|(h, pl)| (h.x_range().0 + pl.x, h.x_range().1 + pl.x))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |a, r| (a.0.min(r.0), a.1.max(r.1)));
        let stations = [x1, mass.lcg];
        let r = sea_response_fleet(&members, &mass, &spectrum, heading, speed, &stations, 41, &opts)
            .map_err(|e| format!("{e}"))?;
        println!("irregular sea {spectrum:?}:");
        println!("  significant heave amplitude   {:.3} m", r.heave);
        println!("  significant pitch amplitude   {:.2}°", r.pitch.to_degrees());
        println!("  significant vertical accel.   {:.2} m/s² at the bow, {:.2} m/s² at the LCG", r.accelerations[0], r.accelerations[1]);
        println!("  mean added resistance         {:.1} N", r.added_resistance);
        if r.skipped_energy > 0.0 {
            println!("  ({:.0}% of the sea's energy was overtaken and left out)", 100.0 * r.skipped_energy);
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
        let v: f64 = v.trim().parse().map_err(|_| format!("--sea: cannot parse {v:?}"))?;
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

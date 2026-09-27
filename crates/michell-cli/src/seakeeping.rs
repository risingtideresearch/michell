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
use michell_seakeeping::platform::{self, Loading};
use michell_seakeeping::sea::{sea_response_fleet, Spectrum};
use michell_seakeeping::strip::StripOptions;

pub(crate) const USAGE: &str =
    "usage: michell seakeeping <hull>[@x=DX,y=Y]... (--speed U | --froude F) \
[--heading DEG] [--lambda A:B:STEP] [--kyy FRAC] [--mass KG] [--lcg X] [--panels N] \
[--sea hs=H,tp=T[,gamma=G]] [--vcg Z] [--kxx K] [--kzz K] [--roll-damping ZETA] [--dynamic] [--csv]";

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
    // Vertical centre of gravity: --vcg metres above the waterline (default
    // on it), giving BG against the platform's centre of buoyancy.
    let vcg = p.f64_flag("vcg")?.unwrap_or(0.0);
    let mass = platform::mass_properties(
        members,
        rho,
        &Loading {
            mass: p.f64_flag("mass")?,
            lcg: p.f64_flag("lcg")?,
            vcg: Some(vcg),
            k_yy: p.f64_flag("kyy")?.map(|f| f * l_ref),
            k_xx: p.f64_flag("kxx")?,
            k_zz: p.f64_flag("kzz")?,
        },
    );
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
        roll_damping: p.f64_flag("roll-damping")?.unwrap_or(0.0),
    };
    let beam = platform::hull_beam(members);
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
    let roll = platform::roll_stability(members, &mass, &opts);
    let gm = roll.gm;
    note(format!(
        "roll: VCG {vcg:.3} m above the waterline, BG {:.3} m, GM_T {gm:.3} m, k_xx {:.3} m, k_zz {:.3} m; \
         natural roll period {} (without added inertia)",
        mass.bg,
        mass.roll_radius_of_gyration,
        mass.yaw_radius_of_gyration,
        match roll.period_dry {
            Some(t) => format!("{t:.2} s"),
            None => "— unstable (GM_T ≤ 0)".into(),
        }
    ));
    if let Some(t) = roll.period {
        note(format!(
            "natural roll period with added inertia {t:.2} s (at rest)"
        ));
    }
    if csv {
        println!("lambda_over_L,omega,omega_e,heave,heave_phase_deg,pitch_over_k,pitch_phase_deg,sway,roll_over_k,yaw_over_k,raw_gb_per_zeta2,sigma_aw_gb,raw_maruo_per_zeta2,sigma_aw_maruo");
    } else {
        println!(
            "{:>7} {:>8} {:>8} {:>8} {:>9} {:>8} {:>9} {:>7} {:>7} {:>7} {:>11} {:>8} {:>11} {:>8}",
            "λ/L",
            "ω[r/s]",
            "ωe[r/s]",
            "heave",
            "ph3[deg]",
            "pitch",
            "ph5[deg]",
            "sway",
            "roll",
            "yaw",
            "Raw/ζ² GB",
            "σ_aw GB",
            "Raw/ζ² far",
            "σ_aw far"
        );
    }
    let points = platform::rao_sweep(members, &mass, heading, speed, &lambdas, &opts, &mut |_| {
        true
    });
    for (&lam, point) in lambdas.iter().zip(points) {
        match point {
            Ok(pt) => {
                let sigma = |raw: f64| raw / (rho * g * beam * beam / l_ref);
                let row = [
                    lam,
                    pt.omega,
                    pt.omega_e,
                    pt.heave.abs(),
                    pt.heave.im.atan2(pt.heave.re).to_degrees(),
                    pt.pitch.abs() / pt.k,
                    pt.pitch.im.atan2(pt.pitch.re).to_degrees(),
                    pt.sway.abs(),
                    pt.roll.abs() / pt.k,
                    pt.yaw.abs() / pt.k,
                    pt.added_resistance,
                    sigma(pt.added_resistance),
                    pt.added_resistance_far_field,
                    sigma(pt.added_resistance_far_field),
                ];
                if csv {
                    let cells: Vec<String> = row.iter().map(|v| format!("{v:.6}")).collect();
                    println!("{}", cells.join(","));
                } else {
                    println!(
                        "{:7.3} {:8.3} {:8.3} {:8.3} {:9.1} {:8.3} {:9.1} {:7.3} {:7.3} {:7.3} {:11.2} {:8.3} {:11.2} {:8.3}",
                        row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13]
                    );
                }
            }
            Err(e) => note(format!("{lam:7.3}  — {e}")),
        }
    }
    note(
        "heave and sway per unit wave amplitude; pitch, roll and yaw as |η|/(kζ); phases relative \
         to a crest at the LCG; roll is potential-flow damped only (plus --roll-damping), so its \
         resonance is overstated on a monohull; \
         added resistance by radiated energy (GB, Gerritsma–Beukelman) and far-field momentum \
         (far, Maruo; head and following seas only — oblique, it lacks the antisymmetric \
         diffraction) — on Journée's Wigley hulls GB is nearer the tank at Fn 0.2, the far field \
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
pub fn parse_sea(s: &str) -> Result<Spectrum, String> {
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

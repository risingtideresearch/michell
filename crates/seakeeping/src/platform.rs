//! A platform's seakeeping in one call: its mass properties from a few
//! loading figures (with the usual defaults), its roll stability and natural
//! period, and its responses and added resistance over a range of
//! wavelengths — what the command line and the web page both report.

use crate::restoring::{buoyancy_depth, transverse_metacentric_height};
use crate::strip::{added_resistance_both, response_fleet, MassProperties, StripOptions, Wave};
use hullgeom::{Placement, SectionalHull, C64};
use std::f64::consts::PI;

/// How the platform is loaded. `None` takes the default noted on each.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Loading {
    /// Mass [kg]; default `ρ∇`.
    pub mass: Option<f64>,
    /// Longitudinal centre of gravity [m, fleet x]; default the LCB.
    pub lcg: Option<f64>,
    /// Height of the centre of gravity above the waterline [m]; default 0.
    pub vcg: Option<f64>,
    /// Pitch radius of gyration [m]; default 0.25 of the longest hull.
    pub k_yy: Option<f64>,
    /// Roll radius of gyration [m]; default each hull's mass at its
    /// centreplane plus 0.35 of its beam.
    pub k_xx: Option<f64>,
    /// Yaw radius of gyration [m]; default the hull spread plus `k_yy`.
    pub k_zz: Option<f64>,
}

/// The longest hull's length [m].
pub fn length(members: &[(&SectionalHull, Placement)]) -> f64 {
    members.iter().map(|(h, _)| h.length()).fold(0.0, f64::max)
}

/// The widest hull's greatest waterline beam [m].
pub fn hull_beam(members: &[(&SectionalHull, Placement)]) -> f64 {
    members
        .iter()
        .map(|(h, _)| {
            let (a, b) = h.x_range();
            (0..=200)
                .map(|i| 2.0 * h.waterline_half_beam(a + (b - a) * i as f64 / 200.0))
                .fold(0.0f64, f64::max)
        })
        .fold(0.0, f64::max)
}

/// The platform's mass properties from its loading.
pub fn mass_properties(
    members: &[(&SectionalHull, Placement)],
    density: f64,
    loading: &Loading,
) -> MassProperties {
    let volume: f64 = members.iter().map(|(h, _)| h.displaced_volume()).sum();
    let volume = volume.max(f64::MIN_POSITIVE);
    let lcb = members
        .iter()
        .map(|(h, pl)| h.displaced_volume() * (h.lcb_x() + pl.x))
        .sum::<f64>()
        / volume;
    // Hull offsets set a multihull's roll and yaw inertia: mass at each
    // hull's centreplane (weighted by displacement), plus its own spread.
    let spread_sq = members
        .iter()
        .map(|(h, pl)| h.displaced_volume() * pl.y * pl.y)
        .sum::<f64>()
        / volume;
    let k_yy = loading.k_yy.unwrap_or(0.25 * length(members));
    MassProperties {
        mass: loading.mass.unwrap_or(density * volume),
        lcg: loading.lcg.unwrap_or(lcb),
        radius_of_gyration: k_yy,
        bg: loading.vcg.unwrap_or(0.0) + buoyancy_depth(members),
        roll_radius_of_gyration: loading
            .k_xx
            .unwrap_or((spread_sq + (0.35 * hull_beam(members)).powi(2)).sqrt()),
        yaw_radius_of_gyration: loading.k_zz.unwrap_or((spread_sq + k_yy * k_yy).sqrt()),
    }
}

/// Roll stability: the transverse metacentric height and the natural roll
/// period, without and with the added inertia (at rest; `None` when
/// unstable, or when the section solves fail).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RollStability {
    /// `GM_T` [m].
    pub gm: f64,
    /// `2π k_xx/√(g GM_T)` [s].
    pub period_dry: Option<f64>,
    /// With the roll added inertia at that frequency, iterated [s].
    pub period: Option<f64>,
}

pub fn roll_stability(
    members: &[(&SectionalHull, Placement)],
    mass: &MassProperties,
    opts: &StripOptions,
) -> RollStability {
    let gm = transverse_metacentric_height(members, mass.bg);
    if !(gm > 0.0) {
        return RollStability {
            gm,
            period_dry: None,
            period: None,
        };
    }
    let volume: f64 = members.iter().map(|(h, _)| h.displaced_volume()).sum();
    let c44 = opts.density * opts.gravity * volume * gm;
    let i44 = mass.mass * mass.roll_radius_of_gyration.powi(2);
    let mut w = (c44 / i44).sqrt();
    let period_dry = Some(2.0 * PI / w);
    for _ in 0..4 {
        let beam_wave = Wave {
            omega: w,
            heading: 0.5 * PI,
            speed: 0.0,
        };
        match response_fleet(members, mass, &beam_wave, opts) {
            Ok(r) => {
                let a44 = r.coefficients.full.added_mass[2][2];
                w = (c44 / (i44 + a44.max(0.0))).sqrt();
            }
            Err(_) => {
                return RollStability {
                    gm,
                    period_dry,
                    period: None,
                }
            }
        }
    }
    RollStability {
        gm,
        period_dry,
        period: Some(2.0 * PI / w),
    }
}

/// The response at one wavelength, per unit wave amplitude; phases relative
/// to a crest at the centre of gravity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RaoPoint {
    /// Wavelength over the longest hull's length.
    pub lambda_over_l: f64,
    /// Wave and encounter frequencies [rad/s].
    pub omega: f64,
    pub omega_e: f64,
    /// Wavenumber [1/m].
    pub k: f64,
    /// Heave and sway [m/m]; pitch (bow up), roll and yaw [rad/m].
    pub heave: C64,
    pub pitch: C64,
    pub sway: C64,
    pub roll: C64,
    pub yaw: C64,
    /// Mean added resistance per wave amplitude squared [N/m²]:
    /// Gerritsma–Beukelman radiated energy, and Maruo's far field.
    pub added_resistance: f64,
    pub added_resistance_far_field: f64,
}

/// The responses at each wavelength `λ/L` (L the longest hull), at heading
/// `heading` [rad, π head seas] and speed `speed` [m/s]. Each point is
/// solved independently; a failed one reports its error. `progress(i)` is
/// called after each; returning `false` stops the sweep there.
pub fn rao_sweep(
    members: &[(&SectionalHull, Placement)],
    mass: &MassProperties,
    heading: f64,
    speed: f64,
    lambdas: &[f64],
    opts: &StripOptions,
    progress: &mut dyn FnMut(usize) -> bool,
) -> Vec<Result<RaoPoint, String>> {
    let l = length(members);
    let mut out = Vec::with_capacity(lambdas.len());
    for (i, &lam) in lambdas.iter().enumerate() {
        let k = 2.0 * PI / (lam * l);
        let wave = Wave {
            omega: (k * opts.gravity).sqrt(),
            heading,
            speed,
        };
        out.push(
            added_resistance_both(members, mass, &wave, opts)
                .map(|(r, gb, far)| RaoPoint {
                    lambda_over_l: lam,
                    omega: wave.omega,
                    omega_e: r.omega_e,
                    k,
                    heave: r.heave,
                    pitch: r.pitch,
                    sway: r.sway,
                    roll: r.roll,
                    yaw: r.yaw,
                    added_resistance: gb,
                    added_resistance_far_field: far,
                })
                .map_err(|e| e.to_string()),
        );
        if !progress(i) {
            break;
        }
    }
    out
}

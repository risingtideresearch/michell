//! A case's statics — its float at rest, its roll stability and
//! its GZ curve — and a study in waves: the platform's responses from one
//! heading, about the attitude its calm-water study floated it at.
//!
//! Heights: a case's VCG is above the hull's **design** waterline
//! (a point on the hull); at an attitude that sinks the hull by `s`, the
//! centre of gravity is `vcg − s` above the water.

use crate::params::{default_lambdas, CaseParams, StudyParams};
use crate::{setup, LoftRequest, Progress, Setup, Tracker, CANCELLED};
use michell_geometry::iges::Platform;
use michell_geometry::stability::{
    gz_curve, FullSection, GzCurve, GzOptions, StabilityHull, StabilityLoad, Wave,
};
use michell_geometry::{Conditions, Placement, SectionalHull, STANDARD_GRAVITY};
use michell_seakeeping::platform::{self as sk, Loading};
use michell_seakeeping::strip::StripOptions;
use serde_json::{json, Value};

/// The platform's mass properties at an attitude sinking it by `sinkage`.
fn loading(
    s: &Setup,
    c: &CaseParams,
    members: &[(&SectionalHull, Placement)],
    sinkage: f64,
) -> michell_seakeeping::strip::MassProperties {
    let l = sk::length(members);
    sk::mass_properties(
        members,
        michell_geometry::Fluid::SEAWATER_15C.density,
        &Loading {
            mass: Some(s.mass),
            lcg: Some(s.lcg),
            vcg: Some(c.vcg.unwrap_or(0.0) - sinkage),
            k_yy: c.kyy.map(|f| f * l),
            k_xx: c.kxx,
            k_zz: c.kzz,
        },
    )
}

fn strip_options(c: &CaseParams) -> StripOptions {
    let fluid = michell_geometry::Fluid::SEAWATER_15C;
    StripOptions {
        density: fluid.density,
        gravity: STANDARD_GRAVITY,
        roll_damping: c.roll_damping,
        ..StripOptions::default()
    }
}

/// Each hull's whole sections, keel to sheer, at the design pose: the hull
/// cut with the water raised past its top, so every station is "under".
fn whole_sections(s: &Setup) -> Result<Vec<StabilityHull>, String> {
    let mut out = Vec::new();
    for (i, pose) in &s.layout {
        let src = s.cut.file.source.as_ref();
        let wl = s.cut.file.waterline_z;
        let (v, _) = src
            .posed_tessellation(*i, wl, pose, &Platform::default())
            .map_err(|e| e.to_string())?;
        let top = v.iter().map(|p| p[2]).fold(f64::NEG_INFINITY, f64::max);
        if !top.is_finite() {
            return Err("the hull has no geometry above its keel".into());
        }
        // A little past the top, so the sheer is a section's own top.
        let over = top + 1e-4 * s.l_ref.max(1e-3);
        let raised = Platform {
            sinkage: over,
            trim: 0.0,
            pivot_x: 0.0,
        };
        let imp = src
            .situate_sectional(*i, wl, pose, &raised, &s.cut.opts)
            .map_err(|e| e.to_string())?
            .ok_or("the hull could not be cut whole")?;
        out.push(StabilityHull {
            sections: imp
                .hull
                .curves()
                .map(|(x, curve)| FullSection {
                    x: x + imp.placement.x,
                    half: curve.iter().map(|&(h, d)| (h, over - d)).collect(),
                })
                .collect(),
            y: imp.placement.y,
        });
    }
    Ok(out)
}

/// A case's float at rest, `(sinkage [m], trim [rad])`, and nothing else:
/// what a study held at rest needs, without the GZ curve.
pub fn at_rest(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
) -> Result<(f64, f64), String> {
    use michell_geometry::float::{solve_equilibrium_sectional, LoadCase};
    use michell_geometry::source::SourceHull;
    let s = setup(name, bytes, cut, c)?;
    let rho = michell_geometry::Fluid::SEAWATER_15C.density;
    let sources: Vec<SourceHull> = s
        .layout
        .iter()
        .map(|&(index, pose)| SourceHull {
            source: s.cut.file.source.as_ref(),
            index,
            waterline_z: s.cut.file.waterline_z,
            pose,
        })
        .collect();
    let eq = solve_equilibrium_sectional(
        &sources,
        &LoadCase {
            mass: s.mass,
            lcg: Some(s.lcg),
        },
        rho,
        &s.cut.opts,
    )
    .map_err(|e| format!("the float at rest: {e}"))?;
    Ok((eq.sinkage, eq.trim))
}

/// A case's statics: its float at rest, its hydrostatics there,
/// its roll stability (GM_T, the natural roll period) and its GZ curve.
pub fn statics(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
) -> Result<Value, String> {
    use michell_geometry::float::{solve_equilibrium_sectional, LoadCase};
    use michell_geometry::source::SourceHull;
    let t0 = std::time::Instant::now();
    let s = setup(name, bytes, cut, c)?;
    let rho = michell_geometry::Fluid::SEAWATER_15C.density;
    let sources: Vec<SourceHull> = s
        .layout
        .iter()
        .map(|&(index, pose)| SourceHull {
            source: s.cut.file.source.as_ref(),
            index,
            waterline_z: s.cut.file.waterline_z,
            pose,
        })
        .collect();
    let eq = solve_equilibrium_sectional(
        &sources,
        &LoadCase {
            mass: s.mass,
            lcg: Some(s.lcg),
        },
        rho,
        &s.cut.opts,
    )
    .map_err(|e| format!("the float at rest: {e}"))?;
    let members: Vec<(&SectionalHull, Placement)> =
        eq.fleet.members.iter().map(|(h, p)| (h, *p)).collect();
    if members.is_empty() {
        return Err("the platform is dry at rest".into());
    }
    let hydro = json!({
        "volume": eq.volume,
        "displacement": rho * eq.volume,
        "lcb": eq.lcb,
        "waterplane_area": eq.waterplane_area,
        "wetted_surface": members.iter().map(|(h, _)| h.wetted_surface()).sum::<f64>(),
        "draft": members.iter().map(|(h, _)| h.draft()).fold(0.0, f64::max),
        "length": sk::length(&members),
        "beam": sk::hull_beam(&members),
    });

    let props = loading(&s, c, &members, eq.sinkage);
    let roll = sk::roll_stability(&members, &props, &strip_options(c));

    let vcg = c.vcg.unwrap_or(0.0);
    let gz = match whole_sections(&s).and_then(|hulls| {
        gz_curve(
            &hulls,
            &StabilityLoad {
                mass: s.mass,
                lcg: s.lcg,
                vcg,
            },
            &GzOptions {
                density: rho,
                ..GzOptions::default()
            },
        )
        .map_err(|e| e.to_string())
    }) {
        Ok(g) => json!({
            "heel_deg": g.points.iter().map(|p| p.heel.to_degrees()).collect::<Vec<_>>(),
            "gz": g.points.iter().map(|p| p.gz).collect::<Vec<_>>(),
            "sinkage": g.points.iter().map(|p| p.sinkage).collect::<Vec<_>>(),
            "trim_deg": g.points.iter().map(|p| p.trim.to_degrees()).collect::<Vec<_>>(),
            "gm": g.gm,
            "max_gz": g.max_gz,
            "heel_at_max_deg": g.heel_at_max.to_degrees(),
            "vanishing_deg": g.vanishing.map(f64::to_degrees),
            "area_30": g.area_30,
            "area_40": g.area_40,
            "area_total": g.area_total,
            "free_trim": true,
            // For drawing the platform at each heel: the centre of buoyancy
            // there and of gravity, in the hull's own axes (x fore, y port,
            // z up from the design waterline). The attitude maps them to
            // the water's: heel φ about x (lifting +y), trim θ bow up about
            // x = 0, then down by the sinkage —
            //   Z = x sinθ + cosθ (y sinφ + z cosφ) − s.
            "cb": g.points.iter().map(|p| p.cb).collect::<Vec<_>>(),
            "trim_rad": g.points.iter().map(|p| p.trim).collect::<Vec<_>>(),
            "g": [s.lcg, 0.0, vcg],
        }),
        Err(e) => json!({ "error": e }),
    };

    let at_rest = s.platform(eq.sinkage, eq.trim);
    Ok(json!({
        "mass": s.mass,
        "lcg": s.lcg,
        "vcg": vcg,
        "vcg_given": c.vcg.is_some(),
        "at_rest": {
            "sinkage": eq.sinkage,
            "trim_rad": eq.trim,
            "trim_deg": eq.trim.to_degrees(),
        },
        "hydrostatics": hydro,
        "roll": {
            "gm_t": roll.gm,
            "period_dry": roll.period_dry,
            "period": roll.period,
            "k_xx": props.roll_radius_of_gyration,
            "k_yy": props.radius_of_gyration,
            "k_zz": props.yaw_radius_of_gyration,
            "bg": props.bg,
            "roll_damping": c.roll_damping,
        },
        "gz": gz,
        "meshes": s.meshes(&at_rest)?,
        // The hull at its design pose, in its own axes (the GZ curve's).
        "body_meshes": s.meshes(&Platform::default())?,
        "seconds": t0.elapsed().as_secs_f64(),
    }))
}

/// A case's GZ curves in a regular wave of length `length` and height
/// `height` held still around it (quasi-static: see
/// [`michell_geometry::stability::Wave`]), beam, stern quartering and bow
/// quartering, with the crest at eight places along a wavelength and
/// under each hull; for each heading the worst of them — the lowest peak
/// arm — and in beam seas the crest under the windward and under the lee
/// hull. Beside each, the calm-water curve at the same heels, and for each
/// the energy to heel from where it rests to the peak arm, `W ∫ GZ dφ`;
/// the worst crest of a heading is the one needing least.
pub fn wave_gz(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
    length: f64,
    height: f64,
) -> Result<Value, String> {
    use std::f64::consts::PI;
    if !(length > 0.0 && height > 0.0 && length.is_finite() && height.is_finite()) {
        return Err("the wave needs a positive length and height".into());
    }
    if height / length > 1.0 / 7.0 {
        return Err(format!(
            "H/λ {:.3}: steeper than a wave can stand (1/7)",
            height / length
        ));
    }
    let t0 = std::time::Instant::now();
    let s = setup(name, bytes, cut, c)?;
    let hulls = whole_sections(&s)?;
    let rho = michell_geometry::Fluid::SEAWATER_15C.density;
    let load = StabilityLoad {
        mass: s.mass,
        lcg: s.lcg,
        vcg: c.vcg.unwrap_or(0.0),
    };
    let heels: Vec<f64> = (0..=40).map(|i| (3.0 * i as f64).to_radians()).collect();
    let weight = s.mass * STANDARD_GRAVITY;
    // What a curve is summed up by: its peak arm and heel, the heel it
    // rests at, and the energy to heel from there to the peak. On a wave
    // whose slope heels it (a negative arm upright), it rests where the arm
    // crosses zero on the way up; one that pushes it back (a positive arm
    // upright) rests at a negative heel, off the curve, so upright stands
    // in for it — which understates the energy, erring toward the smaller
    // wave.
    let summary = |c: &GzCurve| {
        let pts = &c.points;
        let (mut k, mut best) = (0usize, f64::NEG_INFINITY);
        for (i, q) in pts.iter().enumerate() {
            if q.gz > best {
                best = q.gz;
                k = i;
            }
        }
        let (rest, from) = if pts[0].gz < 0.0 {
            match (0..k).find(|&i| pts[i].gz < 0.0 && pts[i + 1].gz >= 0.0) {
                Some(i) => {
                    let f = -pts[i].gz / (pts[i + 1].gz - pts[i].gz);
                    (pts[i].heel + f * (pts[i + 1].heel - pts[i].heel), i)
                }
                None => (pts[k].heel, k),
            }
        } else {
            (0.0, 0)
        };
        // ∫ GZ dφ from the resting heel to the peak, trapezoidal, the
        // first interval from the crossing (where the arm is zero).
        let mut area = 0.0;
        if from < k {
            let q = &pts[from + 1];
            area += 0.5 * (q.heel - rest) * q.gz.max(0.0);
            area += pts[from + 1..=k]
                .windows(2)
                .map(|w| 0.5 * (w[1].heel - w[0].heel) * (w[0].gz + w[1].gz))
                .sum::<f64>();
        }
        json!({
            "gz": pts.iter().map(|q| q.gz).collect::<Vec<_>>(),
            "max_gz": best,
            "heel_at_max_deg": pts[k].heel.to_degrees(),
            "rest_heel_deg": rest.to_degrees(),
            // Pushed back past upright: it rests at a negative heel, and
            // heeling that way is the mirrored crest heeling this way.
            "rests_negative": pts[0].gz > 1e-4,
            // To draw it at each heel: its attitude and centre of buoyancy
            // (the hull's axes; see the statics' `gz` for the map).
            "sinkage": pts.iter().map(|q| q.sinkage).collect::<Vec<_>>(),
            "trim_rad": pts.iter().map(|q| q.trim).collect::<Vec<_>>(),
            "cb": pts.iter().map(|q| q.cb).collect::<Vec<_>>(),
            "energy_to_peak": weight * area,
        })
    };
    let opts = |wave| GzOptions {
        heels: heels.clone(),
        density: rho,
        wave,
        ..GzOptions::default()
    };
    let calm = gz_curve(&hulls, &load, &opts(None)).map_err(|e| e.to_string())?;
    // The windward hull is the port one (a positive heel lifts +y).
    let half = 0.5 * c.span.unwrap_or(0.0);
    let headings = [
        (90.0, "beam"),
        (45.0, "stern quartering"),
        (135.0, "bow quartering"),
    ];
    // (heading, phase, label) for every curve, computed across the cores.
    let mut jobs: Vec<(usize, f64, String)> = Vec::new();
    for (h, &(deg, _)) in headings.iter().enumerate() {
        let mu = f64::to_radians(deg);
        for j in 0..8 {
            jobs.push((h, 2.0 * PI * j as f64 / 8.0, format!("crest at {}/8 λ", j)));
        }
        if half > 0.0 {
            jobs.push((
                h,
                Wave::crest_at(length, mu, s.lcg, half),
                "crest under the windward hull".into(),
            ));
            jobs.push((
                h,
                Wave::crest_at(length, mu, s.lcg, -half),
                "crest under the lee hull".into(),
            ));
        } else {
            jobs.push((
                h,
                Wave::crest_at(length, mu, s.lcg, 0.0),
                "crest amidships".into(),
            ));
            jobs.push((
                h,
                Wave::crest_at(length, mu, s.lcg, 0.0) + PI,
                "trough amidships".into(),
            ));
        }
    }
    let curves = michell_geometry::parallel::map_indexed(
        jobs.len(),
        || (),
        |_, i| {
            let (h, phase, _) = &jobs[i];
            let wave = Wave {
                length,
                height,
                heading: headings[*h].0.to_radians(),
                phase: *phase,
            };
            gz_curve(&hulls, &load, &opts(Some(wave))).map_err(|e| e.to_string())
        },
    );
    let mut out = Vec::new();
    for (h, &(deg, what)) in headings.iter().enumerate() {
        let mine: Vec<(usize, Value)> = jobs
            .iter()
            .enumerate()
            .filter(|(_, j)| j.0 == h)
            .map(|(i, j)| {
                let v = match &curves[i] {
                    Ok(c) => {
                        let mut v = summary(c);
                        v["label"] = json!(j.2);
                        v["phase"] = json!(j.1);
                        v
                    }
                    Err(e) => json!({ "label": j.2, "phase": j.1, "error": e }),
                };
                (i, v)
            })
            .collect();
        // The worst crest: the least energy to reach the peak arm (what a
        // breaking wave must supply), of those it rests at a positive heel
        // on (the others are their mirrored crests, heeling the other way).
        let worst = mine
            .iter()
            .filter(|(_, v)| v["error"].is_null() && v["rests_negative"] != true)
            .min_by(|a, b| {
                a.1["energy_to_peak"]
                    .as_f64()
                    .unwrap_or(f64::INFINITY)
                    .total_cmp(&b.1["energy_to_peak"].as_f64().unwrap_or(f64::INFINITY))
            })
            .map(|(_, v)| v.clone());
        out.push(json!({
            "heading_deg": deg,
            "name": what,
            "worst": worst,
            "curves": mine.into_iter().map(|(_, v)| v).collect::<Vec<_>>(),
        }));
    }
    let l = hulls
        .iter()
        .flat_map(|h| h.sections.iter().map(|q| q.x))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), x| {
            (a.min(x), b.max(x))
        });
    Ok(json!({
        "wave": { "length": length, "height": height },
        "heel_deg": heels.iter().map(|h| h.to_degrees()).collect::<Vec<_>>(),
        "calm": summary(&calm),
        "headings": out,
        "mass": s.mass,
        "weight": weight,
        "vcg": load.vcg,
        "vcg_given": c.vcg.is_some(),
        // The length a beam-on breaking crest strikes along: the hull's.
        "struck_length": l.1 - l.0,
        "seconds": t0.elapsed().as_secs_f64(),
    }))
}

/// A study in waves: the platform held at `hold` (its calm-water study's
/// attitude) and its responses from the study's heading over its
/// wavelengths, with the irregular sea's statistics when it has one.
pub fn waves_with_progress(
    name: &str,
    bytes: Vec<u8>,
    cut: &LoftRequest,
    c: &CaseParams,
    study: &StudyParams,
    hold: (f64, f64),
    report: &mut dyn FnMut(&Progress) -> bool,
) -> Result<Value, String> {
    use michell_seakeeping::sea::sea_response_fleet;
    let w = study.waves.as_ref().ok_or("not a study in waves")?;
    let lambdas = w.lambdas.clone().unwrap_or_else(default_lambdas);
    let t0 = std::time::Instant::now();
    let n = lambdas.len();
    let mut track = Tracker {
        stages: vec![
            ("Cutting sections", 0.02),
            ("Roll stability", 0.05),
            ("Regular waves", 0.93 * n as f64 / 21.0),
            ("Irregular sea", if w.sea.is_some() { 0.93 } else { 0.0 }),
        ],
        report,
    };
    track.at(0, 0.0, name.to_string())?;
    let s = setup(name, bytes, cut, c)?;
    let (sinkage, trim) = hold;
    let platform = s.platform(sinkage, trim);
    let owned = s.situate(&platform)?;
    let members: Vec<(&SectionalHull, Placement)> = owned.iter().map(|(h, p)| (h, *p)).collect();
    let speed = study.froude * (STANDARD_GRAVITY * s.l_ref).sqrt();
    let cond = Conditions::seawater(speed);
    let opts = strip_options(c);
    let props = loading(&s, c, &members, sinkage);
    let l = sk::length(&members);
    track.at(1, 0.0, "roll stability at this attitude".into())?;
    let roll = sk::roll_stability(&members, &props, &opts);
    let deg = w.heading;
    // The far field's sections carry only the symmetric part of the
    // diffraction problem: head and following seas only.
    let far = |v: f64| (deg.to_radians().sin().abs() < 0.1).then_some(v);
    let c64 = |z: michell_geometry::C64| json!([z.re, z.im]);
    let mut stopped = false;
    let points = sk::rao_sweep(
        &members,
        &props,
        deg.to_radians(),
        speed,
        &lambdas,
        &opts,
        &mut |i| {
            stopped = track
                .at(
                    2,
                    (i + 1) as f64 / n as f64,
                    format!("heading {deg:.0}° · λ/L {:.2}", lambdas[i]),
                )
                .is_err();
            !stopped
        },
    );
    if stopped {
        return Err(CANCELLED.into());
    }
    let points: Vec<Value> = points
        .into_iter()
        .zip(&lambdas)
        .map(|(p, lam)| match p {
            Ok(p) => json!({
                "lambda": p.lambda_over_l,
                "omega": p.omega,
                "omega_e": p.omega_e,
                "k": p.k,
                "heave": c64(p.heave),
                "pitch": c64(p.pitch),
                "sway": c64(p.sway),
                "roll": c64(p.roll),
                "yaw": c64(p.yaw),
                "raw_gb": p.added_resistance,
                "raw_far": far(p.added_resistance_far_field),
            }),
            Err(e) => json!({ "lambda": lam, "error": e }),
        })
        .collect();
    let x_bow = members
        .iter()
        .map(|(h, pl)| h.x_range().1 + pl.x)
        .fold(f64::NEG_INFINITY, f64::max);
    let sea = match &w.sea {
        Some(sea) => {
            track.at(3, 0.0, format!("heading {deg:.0}° · irregular sea"))?;
            match sea_response_fleet(
                &members,
                &props,
                &sea.spectrum(),
                deg.to_radians(),
                speed,
                &[x_bow, s.lcg],
                41,
                &opts,
            ) {
                Ok(r) => json!({
                    "heave": r.heave,
                    "pitch_deg": r.pitch.to_degrees(),
                    "accel_bow": r.accelerations[0],
                    "accel_lcg": r.accelerations[1],
                    "raw_gb": r.added_resistance,
                    "raw_far": far(r.added_resistance_far_field),
                    "skipped_energy": r.skipped_energy,
                }),
                Err(e) => json!({ "error": e.to_string() }),
            }
        }
        None => Value::Null,
    };
    track.at(3, 1.0, "done".into())?;
    let meshes = s.meshes(&platform)?;
    Ok(json!({
        "froude": study.froude,
        "speed": cond.speed,
        "seconds": t0.elapsed().as_secs_f64(),
        "attitude": { "sinkage": sinkage, "trim_rad": trim, "trim_deg": trim.to_degrees() },
        "hulls": meshes.into_iter().map(|m| json!({ "mesh": m })).collect::<Vec<_>>(),
        "seakeeping": {
            "length": l,
            "beam": sk::hull_beam(&members),
            "mass": props.mass,
            "lcg": props.lcg,
            "bg": props.bg,
            // G: at the LCG, on the centreplane (y = 0), this far above the
            // water at this attitude — the point the motions are about.
            "vcg": c.vcg.unwrap_or(0.0) - sinkage,
            "k_yy": props.radius_of_gyration,
            "k_xx": props.roll_radius_of_gyration,
            "k_zz": props.yaw_radius_of_gyration,
            "roll_damping": c.roll_damping,
            "gm_t": roll.gm,
            "roll_period_dry": roll.period_dry,
            "roll_period": roll.period,
            "sea": w.sea,
            "headings": [{ "heading": deg, "points": points, "sea": sea }],
        },
    }))
}

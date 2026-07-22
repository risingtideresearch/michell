//! JSON sweep manifests: a study definition referencing full-band `.hull`
//! bodies, with every varying quantity (speed, load, waterline, hull poses)
//! expressed as a sweep axis. See the README for the schema.

use crate::formats::{body_options, load_body, parse_pair, LoadSettings};
use crate::json::{parse as parse_json, Json};
use michell::body::{Body, BodyOptions};
use michell::float::{solve_equilibrium_heeled, FleetState, HeeledEquilibrium, LoadCase};
use michell::iges::{HullPose, Platform};
use michell::inclined::InclinedGrid;
use michell::{Conditions, Hull, Placement, WaveOptions, STANDARD_GRAVITY};

const KNOT: f64 = 1852.0 / 3600.0;

#[derive(Clone, Copy, PartialEq)]
enum PoseParam {
    Dx,
    Dy,
    Dz,
    Spread,
    TrimDeg,
}

enum Target {
    Weight,
    Lcg,
    Vcg,
    Waterline,
    Pose(Vec<usize>, PoseParam),
}

#[derive(Clone, Copy)]
enum SpeedUnit {
    Ms,
    Knots,
    Froude,
}

struct Axis {
    label: String,
    values: Vec<f64>,
    target: Target,
}

struct MHull {
    id: String,
    body: Body,
    base: HullPose,
}

/// Configuration for the heel roll-up metrics (active in vcg/equilibrium mode).
/// Heel is no longer a sweep axis; instead each row carries the resistance rise
/// at a few fixed angles plus GZ-curve summaries.
struct HeelCfg {
    /// Heel angles (deg) at which to report the fractional resistance rise.
    res_angles: Vec<f64>,
    /// GZ-curve scan step (deg).
    gz_step: f64,
    /// GZ-curve scan cap (deg) — the search for the angle of vanishing
    /// stability stops here (also bounded below 90° by the inclined solver).
    gz_max: f64,
}

impl Default for HeelCfg {
    fn default() -> Self {
        HeelCfg {
            res_angles: vec![5.0, 10.0],
            gz_step: 2.5,
            gz_max: 90.0,
        }
    }
}

/// Summaries rolled up from a hull's GZ (righting-arm) curve.
#[derive(Clone, Copy)]
struct GzStats {
    /// Heel angle of peak righting moment [deg].
    peak_deg: f64,
    /// Peak righting moment [N·m].
    rm_peak: f64,
    /// Area under the GZ curve to the vanishing angle [m·rad].
    area: f64,
    /// Angle of vanishing stability (first GZ zero-crossing) [deg].
    vanish_deg: f64,
}

/// Reduce a sampled GZ curve to its summary metrics. `samples` are
/// `(heel_rad, gz_m)` pairs, ascending in heel and starting at `(0, 0)`, at a
/// uniform step. `weight_times_g` [N] converts GZ to righting moment.
///
/// Returns the stats and a `capped` flag that is `true` when the curve never
/// crossed zero within the samples (so `vanish_deg` is the scan cap, not a real
/// angle of vanishing stability). The area is integrated (trapezoid, in
/// radians) only up to the interpolated zero-crossing.
fn gz_curve_stats(samples: &[(f64, f64)], weight_times_g: f64) -> (GzStats, bool) {
    // Locate the first positive→negative crossing and its interpolated angle.
    let mut vanish_rad = samples.last().map(|s| s.0).unwrap_or(0.0);
    let mut cross_idx = None; // index i such that the crossing is in [i, i+1]
    let mut capped = true;
    for i in 0..samples.len().saturating_sub(1) {
        let (x0, y0) = samples[i];
        let (x1, y1) = samples[i + 1];
        if y0 >= 0.0 && y1 < 0.0 {
            let t = y0 / (y0 - y1); // y0 - y1 > 0 here
            vanish_rad = x0 + t * (x1 - x0);
            cross_idx = Some(i);
            capped = false;
            break;
        }
    }

    // Area under the curve up to the vanishing angle (trapezoid, m·rad).
    let last = cross_idx.map(|i| i + 1).unwrap_or(samples.len());
    let mut area = 0.0;
    for w in samples[..last].windows(2) {
        area += 0.5 * (w[0].1 + w[1].1) * (w[1].0 - w[0].0);
    }
    if let Some(i) = cross_idx {
        // Final wedge from the last positive sample to (vanish_rad, 0).
        let (x0, y0) = samples[i];
        area += 0.5 * y0 * (vanish_rad - x0);
    }

    // Peak righting moment: argmax GZ over the pre-vanishing samples, with a
    // parabolic refinement when the max is interior and the neighbours dip.
    let scan = &samples[..last.max(1)];
    let mut pk = 0usize;
    for (i, &(_, y)) in scan.iter().enumerate() {
        if y > scan[pk].1 {
            pk = i;
        }
    }
    let (mut peak_rad, mut peak_gz) = scan[pk];
    if pk > 0 && pk + 1 < scan.len() {
        let (xm, ym) = scan[pk - 1];
        let (x0, y0) = scan[pk];
        let (xp, yp) = scan[pk + 1];
        let denom = ym - 2.0 * y0 + yp;
        let h = 0.5 * ((x0 - xm) + (xp - x0)); // uniform step
        if denom < 0.0 {
            let delta = 0.5 * (ym - yp) / denom; // vertex offset in units of h
            if delta.abs() <= 1.0 {
                peak_rad = x0 + delta * h;
                peak_gz = y0 - 0.25 * (ym - yp) * delta;
            }
        }
    }

    (
        GzStats {
            peak_deg: peak_rad.to_degrees(),
            rm_peak: weight_times_g * peak_gz,
            area,
            vanish_deg: vanish_rad.to_degrees(),
        },
        capped,
    )
}

/// Format an angle for a column name: `5` not `5.0`, but `7.5` kept.
fn fmt_num(v: f64) -> String {
    if v.fract().abs() < 1e-9 {
        format!("{}", v as i64)
    } else {
        format!("{v}")
    }
}

pub fn run(manifest_path: &str) -> Result<(), String> {
    let text = std::fs::read_to_string(manifest_path)
        .map_err(|e| format!("cannot read {manifest_path}: {e}"))?;
    let doc = parse_json(&text).map_err(|e| format!("{manifest_path}: {e}"))?;
    let dir = std::path::Path::new(manifest_path)
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_default();

    // Options.
    let mut settings = LoadSettings::default();
    let mut wave_opts = WaveOptions::default();
    let mut form_factor = 0.0;
    let mut gravity = STANDARD_GRAVITY;
    let mut rho_override = None;
    let mut nu_override = None;
    let mut heel_cfg = HeelCfg::default();
    if let Some(o) = doc.get("options") {
        if let Some(v) = o.get("samples").and_then(Json::as_str) {
            settings.samples = parse_pair(v)?;
        }
        if let Some(v) = o.get("fit_degree").and_then(Json::as_str) {
            let (px, pz) = parse_pair(v)?;
            settings.fit.degree_x = px;
            settings.fit.degree_z = pz;
            settings.fit_explicit = true;
        }
        if let Some(v) = o.get("fit_control").and_then(Json::as_str) {
            let (nx, nz) = parse_pair(v)?;
            settings.fit.n_ctrl_x = nx;
            settings.fit.n_ctrl_z = nz;
            settings.fit_explicit = true;
        }
        if let Some(v) = o.get("rel_tol").and_then(Json::as_f64) {
            wave_opts.rel_tol = v;
        }
        if let Some(v) = o.get("form_factor").and_then(Json::as_f64) {
            form_factor = v;
        }
        if let Some(v) = o.get("gravity").and_then(Json::as_f64) {
            gravity = v;
        }
        rho_override = o.get("rho").and_then(Json::as_f64);
        nu_override = o.get("nu").and_then(Json::as_f64);
        if let Some(h) = o.get("heel") {
            if let Some(a) = h.get("resistance_angles").and_then(Json::as_arr) {
                heel_cfg.res_angles = a
                    .iter()
                    .map(|v| {
                        v.as_f64()
                            .ok_or("options.heel.resistance_angles must be numbers")
                    })
                    .collect::<Result<_, _>>()?;
            }
            if let Some(v) = h.get("gz_step").and_then(Json::as_f64) {
                heel_cfg.gz_step = v;
            }
            if let Some(v) = h.get("gz_max").and_then(Json::as_f64) {
                heel_cfg.gz_max = v;
            }
        }
    }
    for &a in &heel_cfg.res_angles {
        if !a.is_finite() || a.abs() >= 90.0 {
            return Err(format!(
                "options.heel.resistance_angles: {a} must be finite and within ±90 degrees"
            ));
        }
    }
    if heel_cfg.gz_step <= 0.0 || !heel_cfg.gz_step.is_finite() {
        return Err("options.heel.gz_step must be a positive number".into());
    }
    if heel_cfg.gz_max <= 0.0 || !heel_cfg.gz_max.is_finite() {
        return Err("options.heel.gz_max must be a positive number".into());
    }
    let bopts: BodyOptions = body_options(&settings);
    let fluid_name = doc
        .get("fluid")
        .and_then(Json::as_str)
        .unwrap_or("seawater");
    let make_cond = |speed: f64| -> Result<Conditions, String> {
        let mut c = match fluid_name {
            "seawater" => Conditions::seawater(speed),
            "freshwater" => Conditions::freshwater(speed),
            other => return Err(format!("fluid {other:?}: expected seawater or freshwater")),
        };
        if let Some(r) = rho_override {
            c.fluid.density = r;
        }
        if let Some(n) = nu_override {
            c.fluid.kinematic_viscosity = n;
        }
        c.gravity = gravity;
        Ok(c)
    };
    let density = make_cond(1.0)?.fluid.density;

    // Hulls.
    let hull_specs = doc
        .get("hulls")
        .and_then(Json::as_arr)
        .ok_or("manifest needs a \"hulls\" array")?;
    let mut hulls: Vec<MHull> = Vec::new();
    for h in hull_specs {
        let id = h
            .get("id")
            .and_then(Json::as_str)
            .ok_or("every hull needs an \"id\"")?
            .to_string();
        if hulls.iter().any(|x| x.id == id) {
            return Err(format!("duplicate hull id {id:?}"));
        }
        let file = h
            .get("file")
            .and_then(Json::as_str)
            .ok_or_else(|| format!("hull {id:?} needs a \"file\""))?;
        let path = dir.join(file);
        let body = load_body(path.to_str().unwrap_or(file))?;
        let mut base = HullPose::default();
        if let Some(pz) = h.get("pose") {
            base.dx = pz.get("dx").and_then(Json::as_f64).unwrap_or(0.0);
            base.dy = pz.get("dy").and_then(Json::as_f64).unwrap_or(0.0);
            base.dz = pz.get("dz").and_then(Json::as_f64).unwrap_or(0.0);
            base.trim = pz
                .get("trim")
                .and_then(Json::as_f64)
                .unwrap_or(0.0)
                .to_radians();
        }
        eprintln!(
            "hull {id}: {file} (centerplane {:.4}, base y {:.4})",
            body.centerplane(),
            body.centerplane() + base.dy
        );
        hulls.push(MHull { id, body, base });
    }
    if hulls.is_empty() {
        return Err("manifest has no hulls".into());
    }

    // Axes.
    let sweep = doc
        .get("sweep")
        .and_then(Json::as_arr)
        .ok_or("manifest needs a \"sweep\" array")?;
    let mut speed_axis: Option<(Vec<f64>, SpeedUnit)> = None;
    let mut axes: Vec<Axis> = Vec::new();
    for entry in sweep {
        let values = axis_values(entry)?;
        let targets: Vec<String> = match entry.get("target") {
            Some(Json::Str(s)) => vec![s.clone()],
            Some(Json::Arr(a)) => a
                .iter()
                .map(|t| {
                    t.as_str()
                        .map(str::to_string)
                        .ok_or_else(|| "target list entries must be strings".to_string())
                })
                .collect::<Result<_, _>>()?,
            _ => return Err("every sweep entry needs a \"target\"".into()),
        };
        let param = entry.get("param").and_then(Json::as_str);
        if targets.len() == 1 && param.is_none() {
            match targets[0].as_str() {
                "speed" => {
                    if speed_axis.is_some() {
                        return Err("only one speed axis is allowed".into());
                    }
                    let unit = match entry.get("unit").and_then(Json::as_str).unwrap_or("ms") {
                        "ms" | "m/s" => SpeedUnit::Ms,
                        "knots" | "kn" => SpeedUnit::Knots,
                        "froude" | "fn" => SpeedUnit::Froude,
                        other => return Err(format!("unknown speed unit {other:?}")),
                    };
                    speed_axis = Some((values, unit));
                    continue;
                }
                "weight" => {
                    axes.push(Axis {
                        label: "weight".into(),
                        values,
                        target: Target::Weight,
                    });
                    continue;
                }
                "lcg" => {
                    axes.push(Axis {
                        label: "lcg".into(),
                        values,
                        target: Target::Lcg,
                    });
                    continue;
                }
                "vcg" => {
                    axes.push(Axis {
                        label: "vcg".into(),
                        values,
                        target: Target::Vcg,
                    });
                    continue;
                }
                "waterline" => {
                    axes.push(Axis {
                        label: "waterline".into(),
                        values,
                        target: Target::Waterline,
                    });
                    continue;
                }
                _ => {}
            }
        }
        // Hull pose axis.
        let param =
            param.ok_or_else(|| format!("sweep entry targeting {targets:?} needs a \"param\""))?;
        let pp = match param {
            "dx" => PoseParam::Dx,
            "dy" => PoseParam::Dy,
            "dz" => PoseParam::Dz,
            "spread" => PoseParam::Spread,
            "trim" => PoseParam::TrimDeg,
            other => return Err(format!("unknown pose param {other:?}")),
        };
        let idxs: Vec<usize> = targets
            .iter()
            .map(|t| {
                hulls
                    .iter()
                    .position(|h| &h.id == t)
                    .ok_or_else(|| format!("sweep target {t:?} is not a hull id"))
            })
            .collect::<Result<_, _>>()?;
        if pp == PoseParam::Spread {
            for &i in &idxs {
                let y = hulls[i].body.centerplane() + hulls[i].base.dy;
                if y.abs() < 1e-9 {
                    return Err(format!(
                        "spread targets hull {:?} which sits on the centerline",
                        hulls[i].id
                    ));
                }
            }
        }
        axes.push(Axis {
            label: format!("{}:{param}", targets.join("+")),
            values,
            target: Target::Pose(idxs, pp),
        });
    }
    let Some((speed_values, speed_unit)) = speed_axis else {
        return Err("the sweep needs a speed axis (target \"speed\")".into());
    };
    let float_mode = axes.iter().any(|a| matches!(a.target, Target::Weight));
    if axes.iter().any(|a| matches!(a.target, Target::Lcg)) && !float_mode {
        return Err("an lcg axis requires a weight axis".into());
    }
    if float_mode && axes.iter().any(|a| matches!(a.target, Target::Waterline)) {
        return Err(
            "a waterline axis cannot be combined with a weight axis (the \
                    waterline is solved)"
                .into(),
        );
    }
    let vcg_mode = axes.iter().any(|a| matches!(a.target, Target::Vcg));
    if vcg_mode && !float_mode {
        return Err("a vcg axis requires a weight axis (gz is computed at a \
                    solved equilibrium)"
            .into());
    }

    // Base situate: reference length for Froude numbers.
    let mut l_ref = 0.0f64;
    for h in &hulls {
        if let Some(sb) = h
            .body
            .situate(0.0, &h.base, &Platform::default(), &bopts)
            .map_err(|e| format!("hull {}: {e}", h.id))?
        {
            l_ref = l_ref.max(sb.hull.length());
        }
    }
    if l_ref <= 0.0 {
        return Err("no hull is wetted at the design waterline".into());
    }
    let speeds: Vec<f64> = speed_values
        .iter()
        .map(|&v| match speed_unit {
            SpeedUnit::Ms => v,
            SpeedUnit::Knots => v * KNOT,
            SpeedUnit::Froude => v * (gravity * l_ref).sqrt(),
        })
        .collect();

    let points: usize = axes
        .iter()
        .map(|a| a.values.len())
        .product::<usize>()
        .max(1);
    if points * speeds.len() > 100_000 {
        return Err(format!(
            "sweep would produce {} rows; narrow the axes",
            points * speeds.len()
        ));
    }
    if let Some(name) = doc.get("name").and_then(Json::as_str) {
        eprintln!("study: {name}");
    }
    eprintln!(
        "sweep: {points} point(s) x {} speed(s){}",
        speeds.len(),
        if float_mode { ", equilibrium mode" } else { "" }
    );

    // Output setup.
    let (format, out_file) = match doc.get("output") {
        Some(o) => (
            o.get("format")
                .and_then(Json::as_str)
                .unwrap_or("csv")
                .to_string(),
            o.get("file").and_then(Json::as_str).map(str::to_string),
        ),
        None => ("csv".to_string(), None),
    };
    if format != "csv" && format != "json" {
        return Err(format!("output format {format:?}: expected csv or json"));
    }
    // Heel roll-up columns (vcg/equilibrium mode only): per-point GZ summaries
    // and a per-(point×speed) resistance rise at each configured angle.
    let gz_cols: Vec<String> = if vcg_mode {
        ["gz_peak_deg", "rm_peak", "gz_area", "gz_vanish_deg"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    } else {
        Vec::new()
    };
    let rise_cols: Vec<String> = if vcg_mode {
        heel_cfg
            .res_angles
            .iter()
            .map(|a| format!("rt_rise_{}deg", fmt_num(*a)))
            .collect()
    } else {
        Vec::new()
    };
    let header: Vec<String> = axes
        .iter()
        .map(|a| a.label.clone())
        .chain(
            ["sinkage", "trim_deg", "volume", "lcb"]
                .iter()
                .map(|s| s.to_string()),
        )
        .chain(gz_cols.iter().cloned())
        .chain(
            [
                "dry",
                "band_exceeded",
                "speed",
                "froude",
                "rw",
                "rv",
                "rt",
                "pe",
                "interference",
                "cw",
                "ct",
            ]
            .iter()
            .map(|s| s.to_string()),
        )
        .chain(rise_cols.iter().cloned())
        .collect();
    let mut out = String::new();
    if format == "json" {
        out.push('[');
    } else {
        out.push_str(&header.join(","));
        out.push('\n');
    }
    let mut first_row = true;

    let bodies: Vec<&Body> = hulls.iter().map(|h| &h.body).collect();
    // Section-integration resolution for the heeled inclined-waterplane
    // hydrostatics (volume balance, trim, and GZ).
    let incl_grid = InclinedGrid::default();
    let mut idx = vec![0usize; axes.len()];
    for point in 0..points {
        let vals: Vec<f64> = axes.iter().zip(&idx).map(|(a, &i)| a.values[i]).collect();
        let mut poses: Vec<HullPose> = hulls.iter().map(|h| h.base).collect();
        let mut waterline = 0.0f64;
        let mut weight = None;
        let mut lcg = None;
        let mut vcg = None;
        for (a, &v) in axes.iter().zip(&vals) {
            match &a.target {
                Target::Waterline => waterline = v,
                Target::Weight => weight = Some(v),
                Target::Lcg => lcg = Some(v),
                Target::Vcg => vcg = Some(v),
                Target::Pose(idxs, pp) => {
                    for &hi in idxs {
                        let pose = &mut poses[hi];
                        match pp {
                            PoseParam::Dx => pose.dx = hulls[hi].base.dx + v,
                            PoseParam::Dy => pose.dy = hulls[hi].base.dy + v,
                            PoseParam::Dz => pose.dz = hulls[hi].base.dz + v,
                            PoseParam::Spread => {
                                let side = hulls[hi].body.centerplane() + hulls[hi].base.dy;
                                pose.dy = hulls[hi].base.dy + if side < 0.0 { -v } else { v };
                            }
                            PoseParam::TrimDeg => pose.trim = hulls[hi].base.trim + v.to_radians(),
                        }
                    }
                }
            }
        }

        // Upright equilibrium. Heel is no longer a sweep axis — it is rolled up
        // into the per-row metrics below — so this solve (and the resistance
        // columns) are always upright. The weight branch solves flotation; the
        // waterline branch situates at a fixed cut.
        let (state, sinkage, trim_deg, volume, lcb) = if let Some(mass) = weight {
            let eq = solve_equilibrium_heeled(
                &bodies,
                0.0,
                &poses,
                &LoadCase { mass, lcg },
                density,
                0.0,
                vcg.unwrap_or(0.0),
                &bopts,
                incl_grid,
            )
            .map_err(|e| format!("point {}: {e}", point + 1))?;
            (
                eq.fleet,
                eq.sinkage,
                eq.trim.to_degrees(),
                eq.volume,
                eq.lcb,
            )
        } else {
            let mut members = Vec::new();
            let mut dry = 0usize;
            let mut band_exceeded = 0usize;
            let mut volume = 0.0;
            let mut moment = 0.0;
            for (h, pose) in hulls.iter().zip(&poses) {
                match h
                    .body
                    .situate(waterline, pose, &Platform::default(), &bopts)
                    .map_err(|e| format!("point {}: {e}", point + 1))?
                {
                    Some(sb) => {
                        volume += sb.hull.displaced_volume();
                        moment += (sb.hull.lcb_x() + sb.placement.x) * sb.hull.displaced_volume();
                        band_exceeded += sb.band_exceeded;
                        members.push((sb.hull, sb.placement));
                    }
                    None => dry += 1,
                }
            }
            let lcb = if volume > 0.0 { moment / volume } else { 0.0 };
            (
                FleetState {
                    members,
                    dry,
                    band_exceeded,
                },
                0.0,
                0.0,
                volume,
                lcb,
            )
        };

        // Heel roll-up metrics (vcg/equilibrium mode only). GZ is speed-
        // independent, so its curve — and the heeled equilibria used for the
        // resistance rise — are solved once per point here. vcg_mode implies a
        // weight axis, so `weight` is always present in this branch.
        let (gz_stats, heeled): (Option<GzStats>, Vec<Option<HeeledEquilibrium>>) = if vcg_mode {
            let mass = weight.expect("vcg mode implies a weight axis");
            let vcg_v = vcg.unwrap_or(0.0);
            let solve_at = |deg: f64| {
                solve_equilibrium_heeled(
                    &bodies,
                    0.0,
                    &poses,
                    &LoadCase { mass, lcg },
                    density,
                    deg.to_radians(),
                    vcg_v,
                    &bopts,
                    incl_grid,
                )
            };
            // Scan the GZ curve from upright until it crosses zero (the angle of
            // vanishing stability) or hits the cap. The inclined solver is only
            // valid below 90°, so the cap is bounded there.
            let cap = heel_cfg.gz_max.min(89.5);
            let mut samples = vec![(0.0f64, 0.0f64)];
            let mut deg = heel_cfg.gz_step;
            while deg <= cap + 1e-9 {
                match solve_at(deg) {
                    Ok(eq) => {
                        let vanished = eq.gz < 0.0;
                        samples.push((deg.to_radians(), eq.gz));
                        if vanished {
                            break;
                        }
                    }
                    Err(e) => {
                        eprintln!("point {}: GZ scan stopped at {deg}°: {e}", point + 1);
                        break;
                    }
                }
                deg += heel_cfg.gz_step;
            }
            let (stats, capped) = gz_curve_stats(&samples, mass * gravity);
            if capped {
                eprintln!(
                    "point {}: GZ still positive at {cap}° — gz_vanish_deg/gz_area are capped there",
                    point + 1
                );
            }
            // Heeled equilibria at the resistance angles (flotation only; the
            // wave rise per speed is computed from each fleet in the loop).
            let heeled = heel_cfg
                .res_angles
                .iter()
                .map(|&a| match solve_at(a) {
                    Ok(eq) => Some(eq),
                    Err(e) => {
                        eprintln!("point {}: heeled solve at {a}° failed: {e}", point + 1);
                        None
                    }
                })
                .collect();
            (Some(stats), heeled)
        } else {
            (None, Vec::new())
        };

        let members: Vec<(&Hull, Placement)> = state.members.iter().map(|(h, p)| (h, *p)).collect();

        for &u in &speeds {
            let cond = make_cond(u)?;
            let (rw, rv, rt, pe, iff, cw, ct) = if members.is_empty() {
                (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
            } else {
                let r =
                    michell::multihull_resistance_with(&members, &cond, &wave_opts, form_factor)
                        .map_err(|e| format!("point {} U={u}: {e}", point + 1))?;
                (
                    r.wave.resistance,
                    r.viscous_total,
                    r.total,
                    r.effective_power,
                    r.interference,
                    r.cw,
                    r.ct,
                )
            };
            // Fractional total-resistance rise at each heel angle, relative to
            // the upright total at this speed. Uses the asymmetric heel wave
            // kernel on the heeled-flotation fleet. NaN if the fleet is dry, the
            // upright total is non-positive, or the heeled solve failed.
            let mut rises: Vec<f64> = Vec::with_capacity(heeled.len());
            for (he, &angle) in heeled.iter().zip(&heel_cfg.res_angles) {
                let rise = match he {
                    Some(eq) if rt > 0.0 => {
                        let hm: Vec<(&Hull, Placement)> =
                            eq.fleet.members.iter().map(|(h, p)| (h, *p)).collect();
                        if hm.is_empty() {
                            f64::NAN
                        } else {
                            let rh = michell::multihull_resistance_heeled(
                                &hm,
                                &cond,
                                &wave_opts,
                                form_factor,
                                angle.to_radians(),
                            )
                            .map_err(|e| format!("point {} U={u} heel {angle}°: {e}", point + 1))?;
                            rh.total / rt - 1.0
                        }
                    }
                    _ => f64::NAN,
                };
                rises.push(rise);
            }
            let froude = u / (gravity * l_ref).sqrt();
            let nums: Vec<f64> = vals
                .iter()
                .cloned()
                .chain([sinkage, trim_deg, volume, lcb])
                .chain(
                    gz_stats
                        .map(|s| [s.peak_deg, s.rm_peak, s.area, s.vanish_deg])
                        .into_iter()
                        .flatten(),
                )
                .chain([
                    state.dry as f64,
                    state.band_exceeded as f64,
                    u,
                    froude,
                    rw,
                    rv,
                    rt,
                    pe,
                    iff,
                    cw,
                    ct,
                ])
                .chain(rises)
                .collect();
            if format == "json" {
                if !first_row {
                    out.push(',');
                }
                first_row = false;
                out.push('{');
                for (i, (k, v)) in header.iter().zip(&nums).enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    out.push_str(&format!("{k:?}:{v}"));
                }
                out.push('}');
            } else {
                let row: Vec<String> = nums.iter().map(|v| format!("{v}")).collect();
                out.push_str(&row.join(","));
                out.push('\n');
            }
        }
        eprintln!("point {}/{points} done", point + 1);
        for (i, a) in axes.iter().enumerate().rev() {
            idx[i] += 1;
            if idx[i] < a.values.len() {
                break;
            }
            idx[i] = 0;
        }
    }
    if format == "json" {
        out.push(']');
    }
    match out_file {
        Some(f) => {
            let path = dir.join(&f);
            std::fs::write(&path, out).map_err(|e| format!("cannot write {f}: {e}"))?;
            eprintln!("wrote {}", path.display());
        }
        None => print!("{out}"),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "expected {b}, got {a}");
    }

    #[test]
    fn fmt_num_trims_integers() {
        assert_eq!(fmt_num(5.0), "5");
        assert_eq!(fmt_num(10.0), "10");
        assert_eq!(fmt_num(7.5), "7.5");
    }

    #[test]
    fn gz_stats_triangle_curve() {
        // Symmetric triangle peaking at 0.5 rad, back to 0 at 1.0 rad, then
        // negative. Vanishing angle = 1.0 rad; area to there = 0.5 m·rad.
        let s = [
            (0.0, 0.0),
            (0.25, 0.5),
            (0.5, 1.0),
            (0.75, 0.5),
            (1.0, 0.0),
            (1.25, -0.5),
        ];
        let (st, capped) = gz_curve_stats(&s, 2.0); // W·g = 2 N
        assert!(!capped);
        approx(st.vanish_deg, 1.0_f64.to_degrees(), 1e-6);
        approx(st.area, 0.5, 1e-9);
        approx(st.peak_deg, 0.5_f64.to_degrees(), 1e-6); // symmetric → sample apex
        approx(st.rm_peak, 2.0 * 1.0, 1e-9);
    }

    #[test]
    fn gz_stats_interpolates_zero_crossing() {
        // Crossing lands between samples: gz 0.2 → -0.2 over 0.4→0.6 rad ⇒ 0.5 rad.
        let s = [(0.0, 0.0), (0.2, 0.4), (0.4, 0.2), (0.6, -0.2)];
        let (st, capped) = gz_curve_stats(&s, 1.0);
        assert!(!capped);
        approx(st.vanish_deg, 0.5_f64.to_degrees(), 1e-9);
    }

    #[test]
    fn gz_stats_capped_when_never_vanishes() {
        // Monotonic-ish, always positive within the samples.
        let s = [(0.0, 0.0), (0.3, 0.3), (0.6, 0.5), (0.9, 0.6)];
        let (st, capped) = gz_curve_stats(&s, 1.0);
        assert!(capped);
        approx(st.vanish_deg, 0.9_f64.to_degrees(), 1e-9); // last sample = cap
        assert!(st.area > 0.0);
    }

    #[test]
    fn gz_stats_parabolic_peak_refines_off_sample() {
        // Peak between samples: parabola apex should land near 0.55 rad, above
        // the nearest sample value.
        let s = [(0.0, 0.0), (0.4, 0.8), (0.6, 0.9), (0.8, 0.7), (1.0, -0.1)];
        let (st, _) = gz_curve_stats(&s, 1.0);
        let peak_rad = st.peak_deg.to_radians();
        assert!(peak_rad > 0.4 && peak_rad < 0.8, "peak_rad = {peak_rad}");
        assert!(st.rm_peak >= 0.9, "rm_peak = {}", st.rm_peak);
    }
}

/// Resolve an axis's values: `value`, `values`, or `range: [a, b]` with an
/// optional `step` (default: a fifth of the span).
fn axis_values(entry: &Json) -> Result<Vec<f64>, String> {
    if let Some(v) = entry.get("value").and_then(Json::as_f64) {
        return Ok(vec![v]);
    }
    if let Some(a) = entry.get("values").and_then(Json::as_arr) {
        let vals: Option<Vec<f64>> = a.iter().map(Json::as_f64).collect();
        return vals.ok_or_else(|| "\"values\" must be an array of numbers".to_string());
    }
    if let Some(r) = entry.get("range").and_then(Json::as_arr) {
        let [a, b] = r else {
            return Err("\"range\" must be [start, stop]".into());
        };
        let (a, b) = (
            a.as_f64().ok_or("range start must be a number")?,
            b.as_f64().ok_or("range stop must be a number")?,
        );
        if b < a {
            return Err(format!("range [{a}, {b}]: stop must be >= start"));
        }
        if b == a {
            return Ok(vec![a]);
        }
        let step = match entry.get("step").and_then(Json::as_f64) {
            Some(s) if s > 0.0 => s,
            Some(s) => return Err(format!("step must be positive, got {s}")),
            None => (b - a) / 5.0,
        };
        let n = ((b - a) / step + 1e-9).floor() as usize;
        return Ok((0..=n).map(|i| a + i as f64 * step).collect());
    }
    Err("sweep entry needs one of: value, values, range".into())
}

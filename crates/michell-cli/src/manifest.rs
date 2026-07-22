//! JSON sweep manifests: a study definition referencing full-band `.hull`
//! bodies, with every varying quantity (speed, load, waterline, hull poses)
//! expressed as a sweep axis. See the README for the schema.

use crate::formats::{body_options, load_body, parse_pair, LoadSettings};
use crate::json::{parse as parse_json, Json};
use michell::body::{Body, BodyOptions};
use michell::float::{fleet_cg, solve_equilibrium_heeled, FleetState, HullLoad, LoadCase};
use michell::inclined::InclinedGrid;
use michell::iges::{HullPose, Platform};
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

/// Per-hull load parameters, swept independently per hull (offsets from the
/// hull's base load, exactly as pose params offset the base pose). The fleet CG
/// is always derived from these — never set directly.
#[derive(Clone, Copy, PartialEq)]
enum LoadParam {
    Mass,
    Lcg,
    Vcg,
}

enum Target {
    HeelDeg,
    Waterline,
    Pose(Vec<usize>, PoseParam),
    Load(Vec<usize>, LoadParam),
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
    load: HullLoad,
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
        // Per-hull load: mass and local CG. The fleet CG is derived by summing
        // these across hulls (carried through each pose), never set directly.
        // An unspecified longitudinal CG defaults to the hull's midship (like
        // the pose pivot), so a load with no `lcg` trims to ~zero, not to x=0.
        let (x0, x1) = body.surface().x_domain();
        let midship = 0.5 * (x0 + x1);
        let mut load = HullLoad {
            mass: 0.0,
            lcg: midship,
            vcg: 0.0,
        };
        if let Some(ld) = h.get("load") {
            load.mass = ld.get("mass").and_then(Json::as_f64).unwrap_or(0.0);
            load.lcg = ld.get("lcg").and_then(Json::as_f64).unwrap_or(midship);
            load.vcg = ld.get("vcg").and_then(Json::as_f64).unwrap_or(0.0);
            if load.mass < 0.0 {
                return Err(format!("hull {id:?}: load mass must be >= 0"));
            }
        }
        eprintln!(
            "hull {id}: {file} (centerplane {:.4}, base y {:.4}, mass {:.1} kg)",
            body.centerplane(),
            body.centerplane() + base.dy,
            load.mass,
        );
        hulls.push(MHull {
            id,
            body,
            base,
            load,
        });
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
                "heel" => {
                    axes.push(Axis {
                        label: "heel".into(),
                        values,
                        target: Target::HeelDeg,
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
        // Per-hull axis: either a pose param (geometry) or a load param (mass /
        // CG). Both offset the hull's base value and target hull ids.
        let param = param.ok_or_else(|| {
            format!("sweep entry targeting {targets:?} needs a \"param\"")
        })?;
        let idxs: Vec<usize> = targets
            .iter()
            .map(|t| {
                hulls
                    .iter()
                    .position(|h| &h.id == t)
                    .ok_or_else(|| format!("sweep target {t:?} is not a hull id"))
            })
            .collect::<Result<_, _>>()?;
        let target = match param {
            "dx" => Target::Pose(idxs, PoseParam::Dx),
            "dy" => Target::Pose(idxs, PoseParam::Dy),
            "dz" => Target::Pose(idxs, PoseParam::Dz),
            "trim" => Target::Pose(idxs, PoseParam::TrimDeg),
            "spread" => {
                for &i in &idxs {
                    let y = hulls[i].body.centerplane() + hulls[i].base.dy;
                    if y.abs() < 1e-9 {
                        return Err(format!(
                            "spread targets hull {:?} which sits on the centerline",
                            hulls[i].id
                        ));
                    }
                }
                Target::Pose(idxs, PoseParam::Spread)
            }
            "mass" => Target::Load(idxs, LoadParam::Mass),
            "lcg" => Target::Load(idxs, LoadParam::Lcg),
            "vcg" => Target::Load(idxs, LoadParam::Vcg),
            other => return Err(format!("unknown param {other:?}")),
        };
        axes.push(Axis {
            label: format!("{}:{param}", targets.join("+")),
            values,
            target,
        });
    }
    let Some((speed_values, speed_unit)) = speed_axis else {
        return Err("the sweep needs a speed axis (target \"speed\")".into());
    };
    // The fleet floats (solves sinkage/pitch to its weight) whenever it carries
    // mass — a base hull mass or a swept mass axis. Its CG is always derived
    // from the per-hull loads, so there is no global weight/lcg/vcg axis.
    let base_mass: f64 = hulls.iter().map(|h| h.load.mass).sum();
    let mass_axis = axes
        .iter()
        .any(|a| matches!(a.target, Target::Load(_, LoadParam::Mass)));
    let float_mode = base_mass > 0.0 || mass_axis;
    let has_cg_axis = axes
        .iter()
        .any(|a| matches!(a.target, Target::Load(_, LoadParam::Lcg | LoadParam::Vcg)));
    if float_mode && axes.iter().any(|a| matches!(a.target, Target::Waterline)) {
        return Err("a waterline axis cannot be combined with hull mass (the \
                    waterline is solved from the load)"
            .into());
    }
    if has_cg_axis && !float_mode {
        return Err("an lcg/vcg axis needs the fleet to carry mass; give a hull a \
                    \"load\": { \"mass\": ... } or sweep a mass axis"
            .into());
    }
    if axes.iter().any(|a| matches!(a.target, Target::HeelDeg)) && !float_mode {
        return Err("a heel axis needs the fleet to carry mass so the righting \
                    arm is computed at a solved equilibrium; give a hull a \
                    \"load\": { \"mass\": ... }"
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

    let points: usize = axes.iter().map(|a| a.values.len()).product::<usize>().max(1);
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
    let header: Vec<String> = axes
        .iter()
        .map(|a| a.label.clone())
        .chain(["sinkage", "trim_deg", "volume", "lcb"].iter().map(|s| s.to_string()))
        .chain(
            if float_mode {
                &["mass", "lcg", "vcg", "tcg", "gz", "rm"][..]
            } else {
                &[]
            }
            .iter()
            .map(|s| s.to_string()),
        )
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
        let mut loads: Vec<HullLoad> = hulls.iter().map(|h| h.load).collect();
        let mut waterline = 0.0f64;
        let mut heel = 0.0f64;
        for (a, &v) in axes.iter().zip(&vals) {
            match &a.target {
                Target::Waterline => waterline = v,
                Target::HeelDeg => heel = v.to_radians(),
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
                            PoseParam::TrimDeg => {
                                pose.trim = hulls[hi].base.trim + v.to_radians()
                            }
                        }
                    }
                }
                Target::Load(idxs, lp) => {
                    for &hi in idxs {
                        let load = &mut loads[hi];
                        match lp {
                            LoadParam::Mass => load.mass = hulls[hi].load.mass + v,
                            LoadParam::Lcg => load.lcg = hulls[hi].load.lcg + v,
                            LoadParam::Vcg => load.vcg = hulls[hi].load.vcg + v,
                        }
                    }
                }
            }
        }
        // The fleet CG is always derived by summing the per-hull loads carried
        // through their poses — so it tracks dx/dy/dz and the swept masses.
        let cg = fleet_cg(&bodies, &loads, &poses);

        // In float mode the fleet solves its sinkage/pitch to the derived
        // weight, and heel (if any) enters as a true inclined-waterplane cut
        // inside the solve; the righting arm uses the derived vcg/tcg.
        let (state, sinkage, trim_deg, volume, lcb, gz_solved) = if float_mode {
            if !(cg.mass > 0.0) {
                return Err(format!(
                    "point {}: fleet carries no mass (all hull masses are zero)",
                    point + 1
                ));
            }
            let eq = solve_equilibrium_heeled(
                &bodies,
                0.0,
                &poses,
                &LoadCase {
                    mass: cg.mass,
                    lcg: Some(cg.lcg),
                },
                density,
                heel,
                cg.vcg,
                cg.tcg,
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
                eq.gz,
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
                        moment +=
                            (sb.hull.lcb_x() + sb.placement.x) * sb.hull.displaced_volume();
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
                0.0,
            )
        };
        // Derived-CG columns, reported in float mode: the mass-weighted fleet
        // CG (mass, lcg, vcg, tcg), the inclined-cut righting arm gz, and the
        // righting moment rm = m·g·gz.
        let cg_cols = if float_mode {
            Some([
                cg.mass,
                cg.lcg,
                cg.vcg,
                cg.tcg,
                gz_solved,
                cg.mass * gravity * gz_solved,
            ])
        } else {
            None
        };

        let members: Vec<(&Hull, Placement)> =
            state.members.iter().map(|(h, p)| (h, *p)).collect();

        for &u in &speeds {
            let cond = make_cond(u)?;
            let (rw, rv, rt, pe, iff, cw, ct) = if members.is_empty() {
                (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
            } else {
                // The heel axis repositions each demihull (transverse offset +
                // immersion via `heel_poses`); the heel wave kernel adds the
                // remaining rotation-about-own-axis effect, so the resistance
                // column is consistent with the heeled GZ state.
                let r = if heel != 0.0 {
                    michell::multihull_resistance_heeled(
                        &members,
                        &cond,
                        &wave_opts,
                        form_factor,
                        heel,
                    )
                } else {
                    michell::multihull_resistance_with(&members, &cond, &wave_opts, form_factor)
                }
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
            let froude = u / (gravity * l_ref).sqrt();
            let nums: Vec<f64> = vals
                .iter()
                .cloned()
                .chain([sinkage, trim_deg, volume, lcb])
                .chain(cg_cols.into_iter().flatten())
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

//! JSON sweep manifests: a study definition referencing full-band `.hull`
//! bodies, with every varying quantity (speed, load, waterline, hull poses)
//! expressed as a sweep axis. See the README for the schema.

use crate::formats::{body_options, load_body, parse_pair, LoadSettings};
use crate::json::{parse as parse_json, Json};
use michell::body::{Body, BodyOptions};
use michell::float::{solve_equilibrium_heeled, FleetState, LoadCase};
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

enum Target {
    Weight,
    Lcg,
    Vcg,
    HeelDeg,
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
        // Hull pose axis.
        let param = param.ok_or_else(|| {
            format!("sweep entry targeting {targets:?} needs a \"param\"")
        })?;
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
        return Err("a waterline axis cannot be combined with a weight axis (the \
                    waterline is solved)"
            .into());
    }
    let vcg_mode = axes.iter().any(|a| matches!(a.target, Target::Vcg));
    if vcg_mode && !float_mode {
        return Err("a vcg axis requires a weight axis (gz is computed at a \
                    solved equilibrium)"
            .into());
    }
    if axes.iter().any(|a| matches!(a.target, Target::HeelDeg)) && !vcg_mode {
        return Err("a heel axis requires a vcg axis so gz is well defined \
                    (vcg is metres above the design floatplane; use \
                    { \"target\": \"vcg\", \"value\": 0 } to put G on it)"
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
    if format != "csv" && format != "json" && format != "pdf" {
        return Err(format!("output format {format:?}: expected csv, json, or pdf"));
    }
    let want_pdf = format == "pdf";
    if want_pdf && out_file.is_none() {
        return Err("output format \"pdf\" needs an output \"file\"".into());
    }
    // GZ-curve heel sweep for the PDF detail pages (weight mode only).
    let (gz_heel_max, gz_heel_step) = match doc.get("output") {
        Some(o) => (
            o.get("gz_heel_max").and_then(Json::as_f64).unwrap_or(50.0),
            o.get("gz_heel_step").and_then(Json::as_f64).unwrap_or(5.0),
        ),
        None => (50.0, 5.0),
    };
    let header: Vec<String> = axes
        .iter()
        .map(|a| a.label.clone())
        .chain(["sinkage", "trim_deg", "volume", "lcb"].iter().map(|s| s.to_string()))
        .chain(if vcg_mode { &["gz", "rm"][..] } else { &[] }.iter().map(|s| s.to_string()))
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
    // PDF report accumulators (one entry per output row).
    let mut report_rows: Vec<crate::report::Row> = Vec::new();
    let mut index_rows: Vec<Vec<String>> = Vec::new();
    let index_header: Vec<String> = std::iter::once("row".to_string())
        .chain(header.iter().cloned())
        .collect();
    let mut idx = vec![0usize; axes.len()];
    for point in 0..points {
        let vals: Vec<f64> = axes.iter().zip(&idx).map(|(a, &i)| a.values[i]).collect();
        let mut poses: Vec<HullPose> = hulls.iter().map(|h| h.base).collect();
        let mut waterline = 0.0f64;
        let mut weight = None;
        let mut lcg = None;
        let mut vcg = None;
        let mut heel = 0.0f64;
        for (a, &v) in axes.iter().zip(&vals) {
            match &a.target {
                Target::Waterline => waterline = v,
                Target::Weight => weight = Some(v),
                Target::Lcg => lcg = Some(v),
                Target::Vcg => vcg = Some(v),
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
            }
        }

        // Heel enters the hydrostatics as a true inclined-waterplane rotation
        // inside `solve_equilibrium_heeled` (a heel axis requires a vcg axis
        // requires a weight axis, so only the weight branch below can heel).
        let (state, sinkage, trim_deg, volume, lcb, gz_solved) = if let Some(mass) = weight {
            let eq = solve_equilibrium_heeled(
                &bodies,
                0.0,
                &poses,
                &LoadCase { mass, lcg },
                density,
                heel,
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
        // Righting arm (from the inclined cut) and moment. Reported only with a
        // vcg axis, which the constraints tie to a weight axis, so `gz_solved`
        // is the solved-equilibrium value here.
        let gz_rm = vcg.map(|_| (gz_solved, weight.unwrap_or(0.0) * gravity * gz_solved));

        let members: Vec<(&Hull, Placement)> =
            state.members.iter().map(|(h, p)| (h, *p)).collect();

        // Per-point report geometry: fleet outlines and the righting-arm curve
        // (both independent of speed). Only built for the PDF output.
        let point_hulls: Vec<crate::report::HullGeom> = if want_pdf {
            state.members.iter().map(|(h, p)| hull_geom(h, *p)).collect()
        } else {
            Vec::new()
        };
        let point_gz: Option<Vec<(f64, f64)>> = if want_pdf {
            // The plotted curve re-solves the equilibrium at every heel, so run
            // it at reduced loft/section resolution — plenty for a plot, and a
            // few times cheaper. (The table's operating-point `gz` above keeps
            // the full-resolution solve.)
            let mut gz_opts = bopts;
            gz_opts.stations = (bopts.stations / 2).clamp(31, bopts.stations.max(31));
            gz_opts.waterlines = (bopts.waterlines / 2).clamp(11, bopts.waterlines.max(11));
            let gz_grid = InclinedGrid {
                stations: 81,
                band: 41,
            };
            weight.map(|mass| {
                gz_curve(
                    &bodies,
                    &poses,
                    &LoadCase { mass, lcg },
                    density,
                    vcg.unwrap_or(0.0),
                    &gz_opts,
                    gz_grid,
                    gz_heel_max,
                    gz_heel_step,
                )
            })
        } else {
            None
        };
        let l_ref_row = l_ref;

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
                .chain(gz_rm.map(|(gz, rm)| [gz, rm]).into_iter().flatten())
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
            } else if !want_pdf {
                let row: Vec<String> = nums.iter().map(|v| format!("{v}")).collect();
                out.push_str(&row.join(","));
                out.push('\n');
            }

            if want_pdf {
                let n = report_rows.len() + 1;
                // Value table and index row from the header/number pairing.
                let values: Vec<(String, String)> = index_header[1..]
                    .iter()
                    .cloned()
                    .zip(nums.iter().map(|v| fmt_value(*v)))
                    .collect();
                let mut index_row = vec![format!("{n}")];
                index_row.extend(nums.iter().map(|v| fmt_value(*v)));

                let wave = if members.is_empty() {
                    None
                } else {
                    wave_plan(&members, &cond, l_ref_row)
                };
                let axis_desc: Vec<String> = axes
                    .iter()
                    .zip(&vals)
                    .map(|(a, v)| format!("{}={}", a.label, fmt_value(*v)))
                    .collect();
                report_rows.push(crate::report::Row {
                    label: format!("Row {n}"),
                    subtitle: format!(
                        "{}   |   U = {:.3} m/s (Fn {:.3})",
                        axis_desc.join(", "),
                        u,
                        froude
                    ),
                    values,
                    hulls: point_hulls.clone(),
                    gz_curve: point_gz.clone(),
                    heel_deg: heel.to_degrees(),
                    wave,
                });
                index_rows.push(index_row);
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
    if want_pdf {
        let title = doc
            .get("name")
            .and_then(Json::as_str)
            .unwrap_or("Sweep report");
        let bytes = crate::report::build(title, &index_header, &index_rows, &report_rows);
        let f = out_file.expect("pdf output requires a file (checked above)");
        let path = dir.join(&f);
        std::fs::write(&path, bytes).map_err(|e| format!("cannot write {f}: {e}"))?;
        eprintln!("wrote {} ({} page(s))", path.display(), report_rows.len() + 1);
        return Ok(());
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

/// Reduce a situated hull to fleet-frame outlines for the report views
/// (metres; z up, waterline at z = 0). The situated `Hull` carries only the
/// immersed band, so these are the underwater silhouette / sections.
fn hull_geom(hull: &Hull, place: Placement) -> crate::report::HullGeom {
    let surf = hull.surface();
    let (x0, x1) = surf.x_domain();
    let (_, draft) = surf.z_domain();
    let nx = 80usize;
    let nz = 24usize;
    let xf = |i: usize| x0 + (x1 - x0) * i as f64 / nx as f64;

    let mut maxb = 0.0f64;
    for i in 0..=nx {
        maxb = maxb.max(surf.eval(xf(i), 0.0));
    }
    let eps = (1e-3 * maxb).max(1e-6);

    // Profile silhouette: waterline top edge, keel depth beneath.
    let mut top = Vec::with_capacity(nx + 1);
    let mut bot = Vec::with_capacity(nx + 1);
    for i in 0..=nx {
        let x = xf(i);
        let mut keel = 0.0f64;
        for j in 0..=nz {
            let zd = draft * j as f64 / nz as f64;
            if surf.eval(x, zd) > eps {
                keel = zd;
            }
        }
        top.push((x + place.x, 0.0));
        bot.push((x + place.x, -keel));
    }
    let mut profile = top;
    for &p in bot.iter().rev() {
        profile.push(p);
    }

    // Waterplane outline (starboard forward, port aft).
    let mut waterplane = Vec::with_capacity(2 * (nx + 1));
    for i in 0..=nx {
        let x = xf(i);
        waterplane.push((x + place.x, place.y + surf.eval(x, 0.0)));
    }
    for i in (0..=nx).rev() {
        let x = xf(i);
        waterplane.push((x + place.x, place.y - surf.eval(x, 0.0)));
    }

    // Body-plan half sections at interior stations.
    let nsec = 13usize;
    let mut sections = Vec::with_capacity(nsec);
    for s in 0..nsec {
        let x = x0 + (x1 - x0) * (s as f64 + 0.5) / nsec as f64;
        let mut sec = Vec::with_capacity(2 * (nz + 1));
        for j in 0..=nz {
            let zd = draft * j as f64 / nz as f64;
            sec.push((place.y + surf.eval(x, zd).max(0.0), -zd));
        }
        for j in (0..=nz).rev() {
            let zd = draft * j as f64 / nz as f64;
            sec.push((place.y - surf.eval(x, zd).max(0.0), -zd));
        }
        sections.push(sec);
    }

    crate::report::HullGeom {
        profile,
        sections,
        waterplane,
    }
}

/// Righting-arm curve: re-solve the heeled equilibrium (holding displacement,
/// and LCG if given) at each heel from 0 to `heel_max` degrees, reading the
/// inclined-cut `gz`. Stops early if a solve fails (e.g. past deck immersion).
#[allow(clippy::too_many_arguments)]
fn gz_curve(
    bodies: &[&Body],
    poses: &[HullPose],
    load: &LoadCase,
    density: f64,
    vcg: f64,
    opts: &BodyOptions,
    grid: InclinedGrid,
    heel_max: f64,
    step: f64,
) -> Vec<(f64, f64)> {
    let step = if step > 0.0 { step } else { 5.0 };
    let mut curve = Vec::new();
    let mut a = 0.0f64;
    while a <= heel_max + 1e-9 {
        let phi = a.to_radians();
        if phi.abs() >= std::f64::consts::FRAC_PI_2 {
            break;
        }
        match solve_equilibrium_heeled(bodies, 0.0, poses, load, density, phi, vcg, opts, grid) {
            Ok(e) => curve.push((a, e.gz)),
            Err(_) => break,
        }
        a += step;
    }
    curve
}

/// Colour-mapped plan-view wave field around the fleet at this speed.
fn wave_plan(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    l_ref: f64,
) -> Option<crate::report::WavePlan> {
    let (mut x_lo, mut x_hi, mut y_abs) = (f64::MAX, f64::MIN, 0.0f64);
    for (h, p) in members {
        let (a, b) = h.surface().x_domain();
        x_lo = x_lo.min(a + p.x);
        x_hi = x_hi.max(b + p.x);
        y_abs = y_abs.max(p.y.abs());
    }
    if x_lo > x_hi {
        return None;
    }
    let x1 = x_hi + 0.30 * l_ref;
    let x0 = x_lo - 2.0 * l_ref;
    let yh = (0.34 * (x1 - x0)).max(y_abs + 0.6 * l_ref);
    let (y0, y1) = (-yh, yh);
    let nx = 220usize;
    let ny = (((nx as f64) * (y1 - y0) / (x1 - x0)).round() as usize).clamp(48, 300);
    let mut spec = michell::FreeWaveSpectrum::new(members, cond).ok()?;
    let grid = spec.elevation_grid(x0, x1, y0, y1, nx, ny).ok()?;
    Some(crate::report::WavePlan::from_grid(&grid))
}

/// Compact human-readable number for the report tables.
fn fmt_value(v: f64) -> String {
    if !v.is_finite() {
        return "-".to_string();
    }
    if v == 0.0 {
        return "0".to_string();
    }
    let a = v.abs();
    if !(1e-3..1e6).contains(&a) {
        return format!("{v:.3e}");
    }
    if (v.round() - v).abs() < 1e-9 && a < 1e5 {
        return format!("{}", v.round() as i64);
    }
    let s = format!("{v:.4}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    s.to_string()
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

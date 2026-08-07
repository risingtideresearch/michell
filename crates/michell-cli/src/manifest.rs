//! JSON sweep manifests: a study definition referencing full-band `.hull`
//! bodies, with every varying quantity (speed, load, waterline, hull poses)
//! expressed as a sweep axis. See the README for the schema.

use crate::archive::{
    Archive, Rows, SpecSample, KIND_HULLFILE, KIND_MANIFEST, KIND_META, KIND_ROWS,
};
use crate::formats::{body_options, load_body, parse_pair, LoadSettings};
use crate::json::{parse as parse_json, Json};
use michell::body::{Body, BodyOptions};
use michell::float::{
    fleet_cg, solve_equilibrium_heeled, FleetState, HeeledEquilibrium, HullLoad, LoadCase,
    PointLoad,
};
use michell::iges::{HullPose, Platform};
use michell::inclined::InclinedGrid;
use michell::{Conditions, FreeWaveSpectrum, Hull, Placement, WaveOptions, STANDARD_GRAVITY};

const KNOT: f64 = 1852.0 / 3600.0;

#[derive(Clone, Copy, PartialEq)]
enum PoseParam {
    Dx,
    Dy,
    Dz,
    Spread,
    TrimDeg,
    Scale,
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

/// Point-load parameters, swept independently per point (offsets from the
/// point's base value). `Dx`/`Dy`/`Dz` move the mass relative to the hull
/// centerpoint; `Mass` scales it.
#[derive(Clone, Copy, PartialEq)]
enum PointParam {
    Mass,
    Dx,
    Dy,
    Dz,
}

enum Target {
    Waterline,
    Pose(Vec<usize>, PoseParam),
    Load(Vec<usize>, LoadParam),
    /// `(hull index, point index)` pairs and the swept point parameter.
    Point(Vec<(usize, usize)>, PointParam),
}

/// What a sweep-target id resolves to.
enum TargetKind {
    Hull(usize),
    Point(usize, usize),
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
    /// Ids of the point loads in `load.points`, in the same order.
    point_ids: Vec<String>,
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

/// Per-point heel roll-up: the GZ summary, the heeled equilibria solved at the
/// resistance angles, and the raw GZ curve `(heel_rad, gz_m)` samples (retained
/// for the binary archive; empty outside equilibrium mode).
type HeelRollup = (
    Option<GzStats>,
    Vec<Option<HeeledEquilibrium>>,
    Vec<(f64, f64)>,
);

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

/// Run a JSON sweep manifest. Progress and informational lines flow through
/// `report` (the CLI echoes them to stderr); the returned string is the CSV/JSON
/// the CLI prints to stdout, or empty when the manifest names an `output.file`
/// (written here directly).
pub fn run(manifest_path: &str, report: &mut crate::Reporter) -> Result<String, String> {
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
    // Raw bytes of each distinct referenced hull file, kept for the binary
    // archive so a study is fully self-contained.
    let mut hull_files: Vec<(String, Vec<u8>)> = Vec::new();
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
        if !hull_files.iter().any(|(f, _)| f == file) {
            let raw = std::fs::read(&path)
                .map_err(|e| format!("cannot re-read hull file {file}: {e}"))?;
            hull_files.push((file.to_string(), raw));
        }
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
            if let Some(s) = pz.get("scale").and_then(Json::as_f64) {
                if !(s > 0.0 && s.is_finite()) {
                    return Err(format!(
                        "hull {id:?} pose scale must be positive and finite; got {s}"
                    ));
                }
                base.scale = s;
            }
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
            points: Vec::new(),
        };
        if let Some(ld) = h.get("load") {
            load.mass = ld.get("mass").and_then(Json::as_f64).unwrap_or(0.0);
            load.lcg = ld.get("lcg").and_then(Json::as_f64).unwrap_or(midship);
            load.vcg = ld.get("vcg").and_then(Json::as_f64).unwrap_or(0.0);
            if load.mass < 0.0 {
                return Err(format!("hull {id:?}: load mass must be >= 0"));
            }
        }
        // Extra point masses mounted on the hull, offset from its centerpoint.
        let mut point_ids: Vec<String> = Vec::new();
        if let Some(pts) = h.get("points") {
            let arr = pts
                .as_arr()
                .ok_or_else(|| format!("hull {id:?}: \"points\" must be an array"))?;
            for p in arr {
                let pid = p
                    .get("id")
                    .and_then(Json::as_str)
                    .ok_or_else(|| format!("hull {id:?}: every point load needs an \"id\""))?
                    .to_string();
                let mass = p.get("mass").and_then(Json::as_f64).unwrap_or(0.0);
                if mass < 0.0 {
                    return Err(format!("hull {id:?} point {pid:?}: mass must be >= 0"));
                }
                point_ids.push(pid);
                load.points.push(PointLoad {
                    mass,
                    dx: p.get("dx").and_then(Json::as_f64).unwrap_or(0.0),
                    dy: p.get("dy").and_then(Json::as_f64).unwrap_or(0.0),
                    dz: p.get("dz").and_then(Json::as_f64).unwrap_or(0.0),
                });
            }
        }
        let point_mass: f64 = load.points.iter().map(|p| p.mass).sum();
        report(
            &format!(
                "hull {id}: {file} (centerplane {:.4}, base y {:.4}, mass {:.1} kg + {} point(s) {:.1} kg)",
                body.centerplane(),
                body.centerplane() + base.dy,
                load.mass,
                load.points.len(),
                point_mass,
            ),
            None,
        );
        hulls.push(MHull {
            id,
            body,
            base,
            load,
            point_ids,
        });
    }
    if hulls.is_empty() {
        return Err("manifest has no hulls".into());
    }
    // Ids (hull and point) must be globally unique so a sweep target resolves
    // to exactly one thing.
    {
        let mut seen = std::collections::HashSet::new();
        for h in &hulls {
            for id in std::iter::once(&h.id).chain(h.point_ids.iter()) {
                if !seen.insert(id.as_str()) {
                    return Err(format!(
                        "duplicate id {id:?} (hull and point ids must be unique)"
                    ));
                }
            }
        }
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
        // Per-hull or per-point axis. Resolve every target id, then dispatch by
        // whether they are hulls (pose / load params) or point loads
        // (mass/dx/dy/dz). A list must be all one kind so the param is
        // unambiguous.
        let param =
            param.ok_or_else(|| format!("sweep entry targeting {targets:?} needs a \"param\""))?;
        let kinds: Vec<TargetKind> = targets
            .iter()
            .map(|t| {
                if let Some(hi) = hulls.iter().position(|h| &h.id == t) {
                    return Ok(TargetKind::Hull(hi));
                }
                for (hi, h) in hulls.iter().enumerate() {
                    if let Some(pi) = h.point_ids.iter().position(|p| p == t) {
                        return Ok(TargetKind::Point(hi, pi));
                    }
                }
                Err(format!("sweep target {t:?} is not a hull or point-load id"))
            })
            .collect::<Result<_, _>>()?;
        let all_points = kinds.iter().all(|k| matches!(k, TargetKind::Point(..)));
        let all_hulls = kinds.iter().all(|k| matches!(k, TargetKind::Hull(_)));
        let target = if all_points {
            let pairs: Vec<(usize, usize)> = kinds
                .iter()
                .map(|k| match k {
                    TargetKind::Point(h, p) => (*h, *p),
                    TargetKind::Hull(_) => unreachable!(),
                })
                .collect();
            match param {
                "mass" => Target::Point(pairs, PointParam::Mass),
                "dx" => Target::Point(pairs, PointParam::Dx),
                "dy" => Target::Point(pairs, PointParam::Dy),
                "dz" => Target::Point(pairs, PointParam::Dz),
                other => {
                    return Err(format!(
                        "unknown point-load param {other:?} (use mass, dx, dy, dz)"
                    ))
                }
            }
        } else if all_hulls {
            let idxs: Vec<usize> = kinds
                .iter()
                .map(|k| match k {
                    TargetKind::Hull(h) => *h,
                    TargetKind::Point(..) => unreachable!(),
                })
                .collect();
            match param {
                "dx" => Target::Pose(idxs, PoseParam::Dx),
                "dy" => Target::Pose(idxs, PoseParam::Dy),
                "dz" => Target::Pose(idxs, PoseParam::Dz),
                "trim" => Target::Pose(idxs, PoseParam::TrimDeg),
                "scale" => {
                    if values.iter().any(|&v| !(v > 0.0 && v.is_finite())) {
                        return Err(format!(
                            "sweep entry targeting {targets:?}: scale factors must be \
                             positive and finite"
                        ));
                    }
                    Target::Pose(idxs, PoseParam::Scale)
                }
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
            }
        } else {
            return Err(format!(
                "sweep entry targeting {targets:?} mixes hull and point-load ids"
            ));
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
    // mass — a base hull/point mass or a swept mass axis. Its CG is always
    // derived from the per-hull loads, so there is no global weight/lcg/vcg
    // axis, and the heel roll-up rides on any floated sweep.
    let base_mass: f64 = hulls
        .iter()
        .map(|h| h.load.mass + h.load.points.iter().map(|p| p.mass).sum::<f64>())
        .sum();
    let mass_axis = axes.iter().any(|a| {
        matches!(
            a.target,
            Target::Load(_, LoadParam::Mass) | Target::Point(_, PointParam::Mass)
        )
    });
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
        return Err(
            "an lcg/vcg axis needs the fleet to carry mass; give a hull a \
                    \"load\": { \"mass\": ... } or sweep a mass axis"
                .into(),
        );
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
        report(&format!("study: {name}"), None);
    }
    report(
        &format!(
            "sweep: {points} point(s) x {} speed(s){}",
            speeds.len(),
            if float_mode { ", equilibrium mode" } else { "" }
        ),
        Some((0, points)),
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
    if format != "csv" && format != "json" && format != "binary" {
        return Err(format!(
            "output format {format:?}: expected csv, json, or binary"
        ));
    }
    let binary = format == "binary";
    // Spectrum sampling for the binary archive (uniform over the significant
    // angular range, matching `michell spectrum`). Ignored for csv/json.
    let spec_points: usize = doc
        .get("output")
        .and_then(|o| o.get("spectrum"))
        .and_then(|s| s.get("points"))
        .and_then(Json::as_f64)
        .map(|v| (v as usize).max(9))
        .unwrap_or(721);
    // Derived fleet-CG columns (equilibrium mode): the mass-weighted CG that
    // the solve and roll-up ride on, so it is visible as the hulls/loads move.
    let cg_cols: Vec<String> = if float_mode {
        ["mass", "lcg", "vcg", "tcg"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    } else {
        Vec::new()
    };
    // Heel roll-up columns (equilibrium mode only): per-point GZ summaries
    // and a per-(point×speed) resistance rise at each configured angle.
    let gz_cols: Vec<String> = if float_mode {
        ["gz_peak_deg", "rm_peak", "gz_area", "gz_vanish_deg"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    } else {
        Vec::new()
    };
    let rise_cols: Vec<String> = if float_mode {
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
        .chain(cg_cols.iter().cloned())
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
    let mut rows = binary.then(|| Rows::new(axes.len(), header.len() - axes.len()));

    let bodies: Vec<&Body> = hulls.iter().map(|h| &h.body).collect();
    // Section-integration resolution for the heeled inclined-waterplane
    // hydrostatics (volume balance, trim, and GZ).
    let incl_grid = InclinedGrid::default();
    let mut idx = vec![0usize; axes.len()];
    for point in 0..points {
        let vals: Vec<f64> = axes.iter().zip(&idx).map(|(a, &i)| a.values[i]).collect();
        let mut poses: Vec<HullPose> = hulls.iter().map(|h| h.base).collect();
        let mut loads: Vec<HullLoad> = hulls.iter().map(|h| h.load.clone()).collect();
        let mut waterline = 0.0f64;
        for (a, &v) in axes.iter().zip(&vals) {
            match &a.target {
                Target::Waterline => waterline = v,
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
                            PoseParam::Scale => pose.scale = hulls[hi].base.scale * v,
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
                Target::Point(pairs, pp) => {
                    for &(hi, pi) in pairs {
                        let base = &hulls[hi].load.points[pi];
                        let p = &mut loads[hi].points[pi];
                        match pp {
                            PointParam::Mass => p.mass = base.mass + v,
                            PointParam::Dx => p.dx = base.dx + v,
                            PointParam::Dy => p.dy = base.dy + v,
                            PointParam::Dz => p.dz = base.dz + v,
                        }
                    }
                }
            }
        }
        // The fleet CG is always derived by summing the per-hull loads (and
        // point loads) carried through their poses — so it tracks dx/dy/dz and
        // the swept masses.
        let cg = fleet_cg(&bodies, &loads, &poses);

        // Upright equilibrium. Heel is no longer a sweep axis — it is rolled up
        // into the per-row metrics below — so this solve (and the resistance
        // columns) are always upright. In float mode the fleet solves flotation
        // to the derived weight; otherwise it situates at a fixed cut. `gz0` is
        // the upright righting arm (nonzero only for a laterally asymmetric CG).
        let (state, sinkage, trim_deg, volume, lcb, gz0) = if float_mode {
            if cg.mass <= 0.0 {
                return Err(format!(
                    "point {}: fleet carries no mass (all hull and point masses are zero)",
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
                0.0,
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
                0.0,
            )
        };

        // Heel roll-up metrics (equilibrium mode only). GZ is speed-independent,
        // so its curve — and the heeled equilibria used for the resistance rise
        // — are solved once per point here, riding on the derived fleet CG.
        let (gz_stats, heeled, gz_samples): HeelRollup = if float_mode {
            let load = LoadCase {
                mass: cg.mass,
                lcg: Some(cg.lcg),
            };
            let solve_at = |deg: f64| {
                solve_equilibrium_heeled(
                    &bodies,
                    0.0,
                    &poses,
                    &load,
                    density,
                    deg.to_radians(),
                    cg.vcg,
                    cg.tcg,
                    &bopts,
                    incl_grid,
                )
            };
            // Scan the GZ curve from the upright arm (`gz0`; zero for a symmetric
            // CG, nonzero when the load is laterally offset) until it crosses
            // zero (the angle of vanishing stability) or hits the cap. The
            // inclined solver is only valid below 90°, so the cap is bounded.
            let cap = heel_cfg.gz_max.min(89.5);
            let mut samples = vec![(0.0f64, gz0)];
            let mut deg = heel_cfg.gz_step;
            let mut stopped_early = false;
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
                        report(
                            &format!("point {}: GZ scan stopped at {deg}°: {e}", point + 1),
                            None,
                        );
                        stopped_early = true;
                        break;
                    }
                }
                deg += heel_cfg.gz_step;
            }
            let (stats, capped) = gz_curve_stats(&samples, cg.mass * gravity);
            if capped {
                // `capped` means GZ never crossed zero within the samples, so
                // `gz_vanish_deg`/`gz_area` only reach the last scanned angle —
                // whether the scan hit the cap or stopped early on a
                // non-converged solve. Report the real last angle either way.
                let last_deg = samples.last().map(|s| s.0.to_degrees()).unwrap_or(0.0);
                let why = if stopped_early {
                    "the heeled solve stopped converging"
                } else {
                    "GZ was still positive at the scan cap"
                };
                report(
                    &format!(
                        "point {}: {why} — gz_vanish_deg/gz_area truncated at {last_deg:.1}°",
                        point + 1
                    ),
                    None,
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
                        report(
                            &format!("point {}: heeled solve at {a}° failed: {e}", point + 1),
                            None,
                        );
                        None
                    }
                })
                .collect();
            (Some(stats), heeled, samples)
        } else {
            (None, Vec::new(), Vec::new())
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
                    float_mode
                        .then_some([cg.mass, cg.lcg, cg.vcg, cg.tcg])
                        .into_iter()
                        .flatten(),
                )
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
            if let Some(rows) = rows.as_mut() {
                let (nu, twl, samples) = sample_spectrum(&members, &cond, spec_points)
                    .map_err(|e| format!("point {} U={u} spectrum: {e}", point + 1))?;
                rows.push(&vals, &nums[axes.len()..], &gz_samples, nu, twl, &samples);
            } else if format == "json" {
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
        report(
            &format!("point {}/{points} done", point + 1),
            Some((point + 1, points)),
        );
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

    if let Some(rows) = rows {
        // Assemble the self-contained binary archive.
        let mut ar = Archive::default();
        let manifest_name = std::path::Path::new(manifest_path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("manifest.json");
        ar.add(KIND_MANIFEST, manifest_name, text.into_bytes());
        for (file, raw) in hull_files {
            ar.add(KIND_HULLFILE, &file, raw);
        }
        ar.add(
            KIND_META,
            "meta.json",
            build_meta(
                doc.get("name").and_then(Json::as_str),
                fluid_name,
                gravity,
                density,
                l_ref,
                &speeds,
                &axes,
                &header[axes.len()..],
                spec_points,
                &heel_cfg,
                float_mode,
            )
            .into_bytes(),
        );
        ar.add(KIND_ROWS, "rows", rows.into_blob());
        let bytes = ar.into_bytes();

        // Binary must go to a file, never stdout. Default to the manifest stem.
        let f = out_file.unwrap_or_else(|| {
            let stem = std::path::Path::new(manifest_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("sweep");
            format!("{stem}.msw")
        });
        let path = dir.join(&f);
        std::fs::write(&path, &bytes).map_err(|e| format!("cannot write {f}: {e}"))?;
        report(
            &format!("wrote {} ({} bytes)", path.display(), bytes.len()),
            None,
        );
        return Ok(String::new());
    }

    match out_file {
        Some(f) => {
            let path = dir.join(&f);
            std::fs::write(&path, &out).map_err(|e| format!("cannot write {f}: {e}"))?;
            report(&format!("wrote {}", path.display()), None);
            Ok(String::new())
        }
        None => Ok(out),
    }
}

/// Build the `META` blob: a JSON object naming the row columns and recording
/// study-level scalars, so a reader can interpret `ROWS` without hard-coding
/// the schema.
#[allow(clippy::too_many_arguments)]
fn build_meta(
    name: Option<&str>,
    fluid: &str,
    gravity: f64,
    density: f64,
    l_ref: f64,
    speeds: &[f64],
    axes: &[Axis],
    metric_labels: &[String],
    spectrum_points: usize,
    heel: &HeelCfg,
    float_mode: bool,
) -> String {
    let arr = |xs: &[String]| -> String {
        let items: Vec<String> = xs.iter().map(|s| json_str(s)).collect();
        format!("[{}]", items.join(","))
    };
    let nums = |xs: &[f64]| -> String {
        let items: Vec<String> = xs.iter().map(|v| format!("{v}")).collect();
        format!("[{}]", items.join(","))
    };
    let axis_labels: Vec<String> = axes.iter().map(|a| a.label.clone()).collect();
    let mut s = String::from("{");
    s.push_str("\"format\":\"michell-sweep v1\"");
    if let Some(n) = name {
        s.push_str(&format!(",\"name\":{}", json_str(n)));
    }
    s.push_str(&format!(",\"fluid\":{}", json_str(fluid)));
    s.push_str(&format!(",\"gravity\":{gravity}"));
    s.push_str(&format!(",\"density\":{density}"));
    s.push_str(&format!(",\"l_ref\":{l_ref}"));
    s.push_str(&format!(",\"float_mode\":{float_mode}"));
    s.push_str(&format!(",\"speeds_ms\":{}", nums(speeds)));
    s.push_str(&format!(",\"axis_labels\":{}", arr(&axis_labels)));
    s.push_str(&format!(",\"metric_labels\":{}", arr(metric_labels)));
    s.push_str(&format!(",\"spectrum_points\":{spectrum_points}"));
    s.push_str(&format!(
        ",\"gz_scan\":{{\"step_deg\":{},\"max_deg\":{}}}",
        heel.gz_step, heel.gz_max
    ));
    s.push_str(
        ",\"row_layout\":\"per row: f64[axis_labels.len] params, \
         f64[metric_labels.len] metrics, u32 gz_n, (f64 heel_rad, f64 gz_m)*gz_n, \
         f64 wavenumber, f64 transverse_wavelength, u32 spec_n, \
         (f64 theta, f64 amp_re, f64 amp_im, f64 drw_dtheta)*spec_n\"",
    );
    s.push('}');
    s
}

/// Minimal JSON string escaping (quotes and backslashes; control chars are not
/// expected in labels).
fn json_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            _ => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Sample the free-wave spectrum uniformly over its significant angular range
/// (the same detection `michell spectrum` uses), returning the wavenumber, the
/// transverse wavelength, and `n` samples. Empty when the fleet is dry.
fn sample_spectrum(
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    n: usize,
) -> Result<(f64, f64, Vec<SpecSample>), String> {
    if members.is_empty() {
        return Ok((0.0, 0.0, Vec::new()));
    }
    let mut spec = FreeWaveSpectrum::new(members, cond).map_err(|e| format!("{e}"))?;
    let nu = spec.wavenumber();
    let twl = spec.transverse_wavelength();

    // Significant range: where dRw/dθ still matters (peak-relative threshold).
    let lim = 89.5f64.to_radians();
    let scan = 4096;
    let mut peak = 0.0f64;
    for i in 0..=scan {
        let theta = -lim + 2.0 * lim * i as f64 / scan as f64;
        peak = peak.max(spec.resistance_density(theta));
    }
    let mut theta_max = 0.0f64;
    for i in 0..=scan {
        let theta = -lim + 2.0 * lim * i as f64 / scan as f64;
        if spec.resistance_density(theta) > 1e-5 * peak {
            theta_max = theta_max.max(theta.abs());
        }
    }
    let theta_max = (theta_max * 1.05).min(lim).max(1e-3);

    let h = 2.0 * theta_max / (n - 1) as f64;
    let mut samples = Vec::with_capacity(n);
    for i in 0..n {
        let theta = -theta_max + h * i as f64;
        let a = spec.amplitude(theta);
        let d = spec.resistance_density(theta);
        samples.push(SpecSample {
            theta,
            amp_re: a.re,
            amp_im: a.im,
            drw_dtheta: d,
        });
    }
    Ok((nu, twl, samples))
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

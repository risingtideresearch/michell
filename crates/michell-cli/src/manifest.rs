//! JSON sweep manifests: a study definition referencing full-band `.hull`
//! bodies, with every varying quantity (speed, load, waterline, hull poses)
//! expressed as a sweep axis. See the README for the schema.

use crate::formats::{body_options, load_body, parse_pair, LoadSettings};
use crate::json::{parse as parse_json, Json};
use michell::body::{Body, BodyOptions};
use michell::float::{
    fleet_cg, solve_equilibrium_bodies, solve_equilibrium_bodies_dynamic, FleetState, HullLoad,
    LoadCase, PointLoad,
};
use michell::iges::{HullPose, Platform};
use michell::squat::{dynamic_load_closure, SquatOptions};
use michell::{Conditions, Hull, Placement, WaveOptions, STANDARD_GRAVITY};

/// One knot in m/s.
pub(crate) const KNOT: f64 = 1852.0 / 3600.0;

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

pub(crate) struct Axis {
    pub(crate) label: String,
    pub(crate) values: Vec<f64>,
    target: Target,
}

pub(crate) struct MHull {
    id: String,
    pub(crate) body: Body,
    pub(crate) base: HullPose,
    load: HullLoad,
    /// Ids of the point loads in `load.points`, in the same order.
    point_ids: Vec<String>,
}

/// Fluid/gravity settings resolved from `options`/`fluid`, bundled so a
/// per-speed `Conditions` can be rebuilt at any call site. A closure would do
/// inside one function, but the parsed manifest has to cross a function
/// boundary now that `report` shares it, and a struct travels where a
/// borrowing closure does not.
pub(crate) struct FluidCfg {
    fluid_name: String,
    rho_override: Option<f64>,
    nu_override: Option<f64>,
    pub(crate) gravity: f64,
}

impl FluidCfg {
    pub(crate) fn make_cond(&self, speed: f64) -> Result<Conditions, String> {
        let mut c = match self.fluid_name.as_str() {
            "seawater" => Conditions::seawater(speed),
            "freshwater" => Conditions::freshwater(speed),
            other => return Err(format!("fluid {other:?}: expected seawater or freshwater")),
        };
        if let Some(r) = self.rho_override {
            c.fluid.density = r;
        }
        if let Some(n) = self.nu_override {
            c.fluid.kinematic_viscosity = n;
        }
        c.gravity = self.gravity;
        Ok(c)
    }
}

/// A manifest parsed and validated, ready to be swept: every option resolved,
/// every hull loaded, every axis converted to SI. Split out of `run` so the
/// image-producing `report` command can consume the same schema without a
/// second parser drifting away from this one.
pub(crate) struct ParsedManifest {
    pub(crate) doc: Json,
    pub(crate) dir: std::path::PathBuf,
    pub(crate) hulls: Vec<MHull>,
    pub(crate) axes: Vec<Axis>,
    pub(crate) speeds: Vec<f64>,
    pub(crate) points: usize,
    pub(crate) float_mode: bool,
    pub(crate) dynamic_mode: bool,
    pub(crate) squat_opts: SquatOptions,
    pub(crate) wave_opts: WaveOptions,
    pub(crate) viscous: michell::ViscousOptions,
    pub(crate) gravity: f64,
    pub(crate) density: f64,
    pub(crate) bopts: BodyOptions,
    pub(crate) l_ref: f64,
    pub(crate) fluid: FluidCfg,
}

/// The per-point state an axis combination resolves to: every hull's pose and
/// load, plus the fixed waterline cut if the manifest sweeps one. Shared with
/// the `report` command so a sweep and its PDF agree about what a grid point
/// means.
pub(crate) struct PointState {
    pub(crate) poses: Vec<HullPose>,
    pub(crate) loads: Vec<HullLoad>,
    pub(crate) waterline: f64,
}

pub(crate) fn point_state(hulls: &[MHull], axes: &[Axis], vals: &[f64]) -> PointState {
    let mut poses: Vec<HullPose> = hulls.iter().map(|h| h.base).collect();
    let mut loads: Vec<HullLoad> = hulls.iter().map(|h| h.load.clone()).collect();
    let mut waterline = 0.0f64;
    for (a, &v) in axes.iter().zip(vals) {
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
    PointState {
        poses,
        loads,
        waterline,
    }
}

/// Run a JSON sweep manifest. Progress and informational lines flow through
/// `report` (the CLI echoes them to stderr); the returned string is the CSV/JSON
/// the CLI prints to stdout, or empty when the manifest names an `output.file`
/// (written here directly).
pub fn run(manifest_path: &str, report: &mut crate::Reporter) -> Result<String, String> {
    let ParsedManifest {
        doc,
        dir,
        hulls,
        axes,
        speeds,
        points,
        float_mode,
        dynamic_mode,
        squat_opts,
        wave_opts,
        viscous,
        gravity,
        density,
        bopts,
        l_ref,
        fluid,
    } = parse_manifest(manifest_path, report)?;

    if let Some(name) = doc.get("name").and_then(Json::as_str) {
        report(&format!("study: {name}"), None);
    }
    report(
        &format!(
            "sweep: {points} point(s) x {} speed(s){}{}",
            speeds.len(),
            if float_mode { ", equilibrium mode" } else { "" },
            if dynamic_mode { ", dynamic" } else { "" }
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
    if format != "csv" && format != "json" {
        return Err(format!("output format {format:?}: expected csv or json"));
    }
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
    // Dynamic-mode columns: the near-field vertical force and its share of
    // the weight, at the speed-dependent equilibrium each row is solved at.
    let dyn_cols: Vec<String> = if dynamic_mode {
        ["fz", "lift_pct"].iter().map(|s| s.to_string()).collect()
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
        .chain(dyn_cols.iter().cloned())
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
    // Grid order: the last axis varies fastest (odometer order).
    let mut strides = vec![1usize; axes.len()];
    for i in (0..axes.len().saturating_sub(1)).rev() {
        strides[i] = strides[i + 1] * axes[i + 1].values.len();
    }

    /// One output row.
    struct RowOut {
        nums: Vec<f64>,
    }
    /// One evaluated grid point: its rows, one per speed. The axis values
    /// are already the leading columns of each row.
    struct PointOut {
        rows: Vec<RowOut>,
    }

    // Evaluate one grid point. Everything it reads — axes, bodies, options,
    // speeds — is immutable, so points are independent and evaluate in any
    // order, on any thread. Informational lines go through `say`, which the
    // main thread relays to `report` as they arrive.
    let eval_point = |point: usize, _say: &dyn Fn(String)| -> Result<PointOut, String> {
        let vals: Vec<f64> = axes
            .iter()
            .enumerate()
            .map(|(i, a)| a.values[(point / strides[i]) % a.values.len()])
            .collect();
        let PointState {
            poses,
            loads,
            waterline,
        } = point_state(&hulls, &axes, &vals);
        // The fleet CG is always derived by summing the per-hull loads (and
        // point loads) carried through their poses — so it tracks dx/dy/dz and
        // the swept masses.
        let cg = fleet_cg(&bodies, &loads, &poses);

        // Upright equilibrium: the fleet solves flotation to the derived
        // weight in float mode, or situates at a fixed cut otherwise.
        let (state, sinkage, trim_deg, volume, lcb) = if float_mode {
            if cg.mass <= 0.0 {
                return Err(format!(
                    "point {}: fleet carries no mass (all hull and point masses are zero)",
                    point + 1
                ));
            }
            let eq = solve_equilibrium_bodies(
                &bodies,
                0.0,
                &poses,
                &LoadCase {
                    mass: cg.mass,
                    lcg: Some(cg.lcg),
                },
                density,
                &bopts,
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
            let mut band_overshoot = 0.0f64;
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
                        band_overshoot = band_overshoot.max(sb.band_overshoot);
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
                    band_overshoot,
                },
                0.0,
                0.0,
                volume,
                lcb,
            )
        };

        let members: Vec<(&Hull, Placement)> = state.members.iter().map(|(h, p)| (h, *p)).collect();

        let mut rows_out: Vec<RowOut> = Vec::with_capacity(speeds.len());
        // Dynamic mode re-solves the equilibrium per speed (the near-field
        // pressure grows with U); each speed warm-starts from the last, since
        // sinkage and trim vary smoothly along the curve.
        let mut warm: Option<(f64, f64)> = None;
        for &u in &speeds {
            let cond = fluid.make_cond(u)?;
            let dyn_eq = if dynamic_mode {
                let eq = solve_equilibrium_bodies_dynamic(
                    &bodies,
                    0.0,
                    &poses,
                    &LoadCase {
                        mass: cg.mass,
                        lcg: Some(cg.lcg),
                    },
                    density,
                    gravity,
                    &bopts,
                    dynamic_load_closure(&cond, cg.lcg, &squat_opts),
                    warm,
                )
                .map_err(|e| format!("point {} U={u}: {e}", point + 1))?;
                warm = Some((eq.sinkage, eq.trim));
                Some(eq)
            } else {
                None
            };
            // The attitude (and wetted fleet) each row's resistance is taken
            // at: the dynamic equilibrium in dynamic mode, else the upright one.
            let (members_u, sinkage_u, trim_u, volume_u, lcb_u, dry_u, band_u) = match &dyn_eq {
                Some(eq) => (
                    eq.fleet
                        .members
                        .iter()
                        .map(|(h, p)| (h, *p))
                        .collect::<Vec<(&Hull, Placement)>>(),
                    eq.sinkage,
                    eq.trim.to_degrees(),
                    eq.volume,
                    eq.lcb,
                    eq.fleet.dry,
                    eq.fleet.band_exceeded,
                ),
                None => (
                    members.clone(),
                    sinkage,
                    trim_deg,
                    volume,
                    lcb,
                    state.dry,
                    state.band_exceeded,
                ),
            };
            let (rw, rv, rt, pe, iff, cw, ct) = if members_u.is_empty() {
                (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
            } else {
                let r = michell::multihull_resistance_with(&members_u, &cond, &wave_opts, &viscous)
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
                .chain([sinkage_u, trim_u, volume_u, lcb_u])
                .chain(
                    float_mode
                        .then_some([cg.mass, cg.lcg, cg.vcg, cg.tcg])
                        .into_iter()
                        .flatten(),
                )
                .chain(
                    dyn_eq
                        .as_ref()
                        .map(|eq| [eq.dynamic.force_up, 100.0 * eq.lift_fraction])
                        .into_iter()
                        .flatten(),
                )
                .chain([
                    dry_u as f64,
                    band_u as f64,
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
            rows_out.push(RowOut { nums });
        }
        Ok(PointOut { rows: rows_out })
    };

    // Fan the points out across the cores: workers pull the next point off a
    // shared counter (so slow and fast points balance), and each worker's
    // library calls get an equal share of the remaining cores — a study with
    // fewer points than cores still fills the machine through the solver's
    // own θ-node parallelism. Progress and warnings are relayed live over a
    // channel; rows are gathered per point and written in grid order
    // afterwards, so the output is independent of scheduling.
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::mpsc;
    enum Msg {
        Line(String),
        Done(usize, Result<PointOut, String>),
    }
    let cores = michell::parallel::available();
    let workers = cores.min(points).max(1);
    let inner_threads = (cores / workers).max(1);
    let next = AtomicUsize::new(0);
    let (tx, rx) = mpsc::channel::<Msg>();
    let eval_point = &eval_point;
    let mut results: Vec<(usize, Result<PointOut, String>)> = Vec::with_capacity(points);
    std::thread::scope(|scope| {
        for _ in 0..workers {
            let tx = tx.clone();
            let next = &next;
            scope.spawn(move || {
                michell::parallel::set_threads(inner_threads);
                let say = |line: String| {
                    let _ = tx.send(Msg::Line(line));
                };
                loop {
                    let point = next.fetch_add(1, Ordering::Relaxed);
                    if point >= points {
                        break;
                    }
                    let r = eval_point(point, &say);
                    if r.is_err() {
                        // Stop handing out work; the error surfaces below.
                        next.fetch_max(points, Ordering::Relaxed);
                    }
                    if tx.send(Msg::Done(point, r)).is_err() {
                        break;
                    }
                }
            });
        }
        drop(tx);
        // The channel closes once every worker has finished.
        let mut done = 0usize;
        for msg in rx {
            match msg {
                Msg::Line(line) => report(&line, None),
                Msg::Done(point, r) => {
                    done += 1;
                    report(
                        &format!("point {}/{points} done", point + 1),
                        Some((done, points)),
                    );
                    results.push((point, r));
                }
            }
        }
    });
    results.sort_by_key(|(point, _)| *point);
    for (_, r) in results {
        let PointOut { rows: point_rows } = r?;
        for row in point_rows {
            let nums = row.nums;
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
    }
    if format == "json" {
        out.push(']');
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

/// Parse and validate a manifest without sweeping it. `run` and the `report`
/// command both start here, so the schema has exactly one reader.
pub(crate) fn parse_manifest(
    manifest_path: &str,
    report: &mut crate::Reporter,
) -> Result<ParsedManifest, String> {
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
    let mut viscous = michell::ViscousOptions::default();
    let mut gravity = STANDARD_GRAVITY;
    let mut rho_override = None;
    let mut nu_override = None;
    let mut dynamic_mode = false;
    let mut squat_opts = SquatOptions::default();
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
            viscous.form_factor = v;
        }
        if let Some(v) = o.get("roughness") {
            let spec = v.as_str().ok_or_else(|| {
                "options.roughness: expected a string — off | cf=DELTA_CF | ks=HEIGHT".to_string()
            })?;
            viscous.roughness = crate::formats::parse_roughness(spec)?;
        }
        if let Some(v) = o.get("transom") {
            let spec = v.as_str().ok_or_else(|| {
                "options.transom: expected a string — off | ballistic[=COEFF] | hollow=METRES"
                    .to_string()
            })?;
            wave_opts.transom = crate::parse_transom(Some(spec))?;
        }
        if let Some(v) = o.get("gravity").and_then(Json::as_f64) {
            gravity = v;
        }
        rho_override = o.get("rho").and_then(Json::as_f64);
        nu_override = o.get("nu").and_then(Json::as_f64);
        if let Some(v) = o.get("dynamic") {
            dynamic_mode = match v {
                Json::Bool(b) => *b,
                _ => return Err("options.dynamic: expected true or false".into()),
            };
        }
        if let Some(v) = o.get("squat_tol").and_then(Json::as_f64) {
            squat_opts.rel_tol = v;
        }
    }
    squat_opts.wave.transom = wave_opts.transom;
    let bopts: BodyOptions = body_options(&settings);
    let fluid_name = doc
        .get("fluid")
        .and_then(Json::as_str)
        .unwrap_or("seawater");
    let fluid = FluidCfg {
        fluid_name: fluid_name.to_string(),
        rho_override,
        nu_override,
        gravity,
    };
    let density = fluid.make_cond(1.0)?.fluid.density;

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
    if dynamic_mode && !float_mode {
        return Err(
            "options.dynamic needs the fleet to carry mass (dynamic sinkage/trim \
                    is solved at the platform's equilibrium); give a hull a \
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

    Ok(ParsedManifest {
        doc,
        dir,
        hulls,
        axes,
        speeds,
        points,
        float_mode,
        dynamic_mode,
        squat_opts,
        wave_opts,
        viscous,
        gravity,
        density,
        bopts,
        l_ref,
        fluid,
    })
}

//! `michell-web` — a browser front end for the `michell` tools.
//!
//! Upload a hull and see it the way the physics sees it: cut into sections.
//! IGES and STL hulls are cut straight from their patches or triangles —
//! the same loader the CLI uses ([`michell_cli::fleet`]). The page is sent each
//! station's section curve (what the depth integral integrates), the CAD
//! ray hits it was interpolated from, the
//! depth-integral curve the kernel interpolates along x, the hydrostatics,
//! and the transom with what the closure needs to draw its virtual appendage
//! at any speed.

use michell::{Conditions, Placement, TransomClosure, WaveOptions, STANDARD_GRAVITY};
use michell_cli::fleet::{open_source_bytes, Kind, LoadSettings};
use michell_geometry::iges::{HullPose, Platform, SectionalImport};
use michell_geometry::SectionalHull;
use serde_json::{json, Value};

pub mod case;
pub mod store;
pub mod worker;

/// Largest upload accepted [bytes]. Big enough for a finely tessellated STL.
pub const MAX_UPLOAD: usize = 128 << 20;

/// Import options the page can set.
#[derive(Default, Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LoftRequest {
    /// Design waterline height in the file's frame [m] (IGES).
    pub waterline: Option<f64>,
    /// Centreplane override [m] (IGES).
    pub centerplane: Option<f64>,
    /// Stations along the hull (IGES).
    pub stations: Option<usize>,
    /// Rays across each section (IGES).
    pub rays: Option<usize>,
    /// Scale to metres (STL, which carries no units): mm, m, in, ... or a
    /// number.
    pub units: Option<f64>,
}

impl LoftRequest {
    /// Parse from `key=value` query pairs; unknown keys are ignored.
    pub fn from_query(pairs: &[(String, String)]) -> Result<LoftRequest, String> {
        let mut r = LoftRequest::default();
        let num = |k: &str, v: &str| {
            v.trim()
                .parse::<f64>()
                .map_err(|_| format!("{k}: expected a number, got {v:?}"))
        };
        let count = |k: &str, v: &str, min: usize| match v.trim().parse::<usize>() {
            Ok(n) if n >= min => Ok(n),
            _ => Err(format!(
                "{k}: expected a count of at least {min}, got {v:?}"
            )),
        };
        for (k, v) in pairs {
            if v.trim().is_empty() {
                continue;
            }
            match k.as_str() {
                "waterline" => r.waterline = Some(num(k, v)?),
                "centerplane" => r.centerplane = Some(num(k, v)?),
                "stations" => r.stations = Some(count(k, v, 8)?),
                "rays" => r.rays = Some(count(k, v, 5)?),
                "units" => r.units = Some(michell_cli::parse_units(v.trim())?),
                _ => {}
            }
        }
        Ok(r)
    }
}

/// Cut an uploaded hull into sections and describe it as the JSON the page
/// draws.
pub fn loft(name: &str, bytes: Vec<u8>, req: &LoftRequest) -> Result<Value, String> {
    let t0 = std::time::Instant::now();
    let cut = cut(name, bytes, req)?;
    let hulls: Vec<Value> = cut
        .hulls
        .iter()
        .map(|h| sectioned_json(h, cut.kind))
        .collect();
    Ok(json!({
        "name": name,
        "seconds": t0.elapsed().as_secs_f64(),
        "notes": cut.notes,
        "hulls": hulls,
    }))
}

/// The whole of an upload's geometry for display, topsides and all, in the
/// file's own frame (z up, metres): each hull's tessellation as base64
/// little-endian arrays, and the geometry's z range — what the waterline
/// can be set within.
pub fn geometry(name: &str, bytes: Vec<u8>, units: Option<f64>) -> Result<Value, String> {
    let settings = LoadSettings {
        units,
        ..LoadSettings::default()
    };
    let file = michell_cli::fleet::open_geometry(name, bytes, &settings)?;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut meshes = Vec::new();
    for i in 0..file.source.len() {
        let (v, t) = file
            .source
            .posed_tessellation(i, 0.0, &HullPose::default(), &Platform::default())
            .map_err(|e| e.to_string())?;
        for p in &v {
            lo = lo.min(p[2]);
            hi = hi.max(p[2]);
        }
        let flat: Vec<f32> = v.iter().flat_map(|p| p.map(|c| c as f32)).collect();
        let idx: Vec<u32> = t.iter().flatten().copied().collect();
        meshes.push(json!({
            "vertices": b64(bytemuck_f32(&flat)),
            "triangles": b64(bytemuck_u32(&idx)),
        }));
    }
    Ok(json!({ "z_range": [lo, hi], "meshes": meshes }))
}

/// An upload cut into sections at its design pose, with the source kept so
/// its hulls can be re-cut at another attitude.
struct Cut {
    kind: Kind,
    file: michell_cli::fleet::SourceFile,
    /// Each cut hull's index in the file.
    index: Vec<usize>,
    hulls: Vec<SectionalImport>,
    opts: michell_geometry::iges::SectionalOptions,
    notes: Vec<String>,
}

fn cut(name: &str, bytes: Vec<u8>, req: &LoftRequest) -> Result<Cut, String> {
    let d = LoadSettings::default();
    let settings = LoadSettings {
        waterline_z: req.waterline.unwrap_or(0.0),
        centerplane: req.centerplane,
        stations: req.stations.unwrap_or(d.stations),
        rays: req.rays.unwrap_or(d.rays),
        units: req.units,
    };
    let file = open_source_bytes(name, bytes, &settings)?;
    let opts = settings.sectional(file.waterline_z);
    let mut notes = Vec::new();
    let mut hulls = Vec::new();
    let mut index = Vec::new();
    for i in 0..file.source.len() {
        match file.source.situate_sectional(
            i,
            file.waterline_z,
            &HullPose::default(),
            &Platform::default(),
            &opts,
        ) {
            Ok(Some(h)) => {
                hulls.push(h);
                index.push(i);
            }
            Ok(None) => notes.push(format!("hull {} is dry at this waterline", i + 1)),
            Err(e) => notes.push(format!("hull {} not sectioned: {e}", i + 1)),
        }
    }
    if hulls.is_empty() {
        return Err(notes
            .first()
            .cloned()
            .unwrap_or_else(|| "no hull found".into()));
    }
    Ok(Cut {
        kind: file.kind,
        file,
        index,
        hulls,
        opts,
        notes,
    })
}

/// The solver's dynamic-load model with a progress report on every full
/// evaluation (the solver's cheap slope probes pass straight through), so a
/// cancelled request also stops the Newton loop.
struct Counted<'a, R> {
    inner: michell::sectional::SectionalDynamic<'a>,
    report: R,
}

impl<R> michell_geometry::float::DynamicModel<SectionalHull> for Counted<'_, R>
where
    R: FnMut(&michell_geometry::float::DynamicLoad, f64) -> michell::Result<()>,
{
    fn load(
        &mut self,
        fleet: &michell_geometry::float::FleetState<SectionalHull>,
    ) -> michell::Result<michell_geometry::float::DynamicLoad> {
        let d = self.inner.load(fleet)?;
        let vol: f64 = fleet
            .members
            .iter()
            .map(|(h, _)| h.displaced_volume())
            .sum();
        (self.report)(&d, vol)?;
        Ok(d)
    }

    fn probe(
        &mut self,
        fleet: &michell_geometry::float::FleetState<SectionalHull>,
    ) -> Option<michell::Result<michell_geometry::float::DynamicLoad>> {
        self.inner.probe(fleet)
    }
}

/// How a case's mass is carried (see [`FlowRequest::mass_by`]).
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default, serde::Serialize, serde::Deserialize)]
pub enum MassBy {
    /// The hull floats deeper or shallower: the waterline moves.
    #[default]
    #[serde(rename = "sinking")]
    Sinking,
    /// The hull is scaled uniformly, `k = (m/m₀)^(1/3)`, about its design
    /// waterline: the waterline stays.
    #[serde(rename = "scale")]
    ScaleXyz,
    /// Beam and draft are scaled, `k = (m/m₀)^(1/2)`, the length kept: the
    /// waterline stays.
    #[serde(rename = "scale_yz")]
    ScaleYz,
}

/// The flow request: a speed and the transom closure, on top of the cut.
pub struct FlowRequest {
    pub cut: LoftRequest,
    /// Length Froude number on the longest hull.
    pub froude: f64,
    pub closure: TransomClosure,
    /// Free-surface grid columns (rows follow the aspect).
    pub grid: usize,
    /// Float at the dynamic equilibrium (sinkage and trim at speed) rather
    /// than the design attitude.
    pub dynamic: bool,
    /// Load [kg] and LCG [m, fleet x]; default: the design displacement and
    /// LCB, so the hull floats at its design waterline at rest.
    pub mass: Option<f64>,
    pub lcg: Option<f64>,
    /// How a given mass is carried: by the hull sinking to it, or by scaling
    /// the hull (in all three axes, or in beam and draft only) until it
    /// displaces the mass at its design waterline.
    pub mass_by: MassBy,
    /// Double the (single) hull into a catamaran with this centre span:
    /// the distance between the demihulls' centreplanes [m].
    pub span: Option<f64>,
    /// A previous solution `(sinkage [m], trim [rad])` to start the
    /// equilibrium from — the last speed's, as the slider moves.
    pub warm: Option<(f64, f64)>,
    /// Hold the platform at this `(sinkage [m], trim [rad])` instead of
    /// solving for it (a fast span sweep's shared attitude).
    pub hold: Option<(f64, f64)>,
    /// Also compute the platform's seakeeping at this attitude.
    pub seakeeping: Option<SeakeepingRequest>,
}

/// Seakeeping in regular waves (and optionally one irregular sea) for a
/// flow case, at the case's load and attitude.
#[derive(Debug, Clone, PartialEq)]
pub struct SeakeepingRequest {
    /// Headings [deg]: 180 head seas, 90 beam, 0 following.
    pub headings: Vec<f64>,
    /// Wavelengths over the longest hull's length.
    pub lambdas: Vec<f64>,
    /// Centre of gravity above the waterline [m] (default 0).
    pub vcg: Option<f64>,
    /// Radii of gyration: pitch as a fraction of L, roll and yaw in metres.
    pub kyy: Option<f64>,
    pub kxx: Option<f64>,
    pub kzz: Option<f64>,
    /// Roll damping as a fraction of critical.
    pub roll_damping: f64,
    pub sea: Option<michell_seakeeping::sea::Spectrum>,
}

impl SeakeepingRequest {
    /// `sk=1` turns it on; `sk_heading` (a list or range, degrees),
    /// `sk_lambda` (a range), `sk_vcg`, `sk_kyy` (of L), `sk_kxx`,
    /// `sk_kzz`, `sk_damping`, `sk_sea` (`hs=H,tp=T[,gamma=G]`).
    pub fn from_query(pairs: &[(String, String)]) -> Result<Option<SeakeepingRequest>, String> {
        let get = |k: &str| {
            pairs
                .iter()
                .find(|(q, _)| q == k)
                .map(|(_, v)| v.trim())
                .filter(|v| !v.is_empty())
        };
        if !matches!(get("sk"), Some("1" | "true" | "on")) {
            return Ok(None);
        }
        let num = |k: &str| -> Result<Option<f64>, String> {
            get(k)
                .map(|v| {
                    v.parse::<f64>()
                        .map_err(|_| format!("{k}: expected a number, got {v:?}"))
                })
                .transpose()
        };
        let list = |k: &str, default: &str| -> Result<Vec<f64>, String> {
            let mut out = Vec::new();
            for part in get(k).unwrap_or(default).split(',') {
                out.extend(michell_cli::parse_range(part.trim()).map_err(|e| format!("{k}: {e}"))?);
            }
            Ok(out)
        };
        let headings = list("sk_heading", "180")?;
        let lambdas = list("sk_lambda", "0.5:3:0.125")?;
        if lambdas.is_empty() || lambdas.len() > 200 || lambdas.iter().any(|&l| !(l > 0.0)) {
            return Err("sk_lambda: expected 1 to 200 positive wavelengths (λ/L)".into());
        }
        if headings.is_empty() || headings.len() > 12 {
            return Err("sk_heading: expected 1 to 12 headings per case".into());
        }
        Ok(Some(SeakeepingRequest {
            headings,
            lambdas,
            vcg: num("sk_vcg")?,
            kyy: num("sk_kyy")?,
            kxx: num("sk_kxx")?,
            kzz: num("sk_kzz")?,
            roll_damping: num("sk_damping")?.unwrap_or(0.0).max(0.0),
            sea: get("sk_sea").map(michell_cli::parse_sea).transpose()?,
        }))
    }
}

impl FlowRequest {
    /// `froude`, `closure` (`ballistic` | `fixed` | `off`) with `param` (the
    /// ballistic coefficient or the fixed hollow length), `grid`, and the
    /// cut's own keys.
    pub fn from_query(pairs: &[(String, String)]) -> Result<FlowRequest, String> {
        let get = |k: &str| pairs.iter().find(|(q, _)| q == k).map(|(_, v)| v.trim());
        let num = |k: &str| -> Result<Option<f64>, String> {
            get(k)
                .filter(|v| !v.is_empty())
                .map(|v| {
                    v.parse::<f64>()
                        .map_err(|_| format!("{k}: expected a number, got {v:?}"))
                })
                .transpose()
        };
        let froude = num("froude")?.ok_or("froude is required")?;
        if !(froude > 0.0 && froude < 5.0) {
            return Err(format!("froude {froude}: expected 0 < Fn < 5"));
        }
        let param = num("param")?;
        let closure = match get("closure").unwrap_or("ballistic") {
            "off" | "none" => TransomClosure::None,
            "fixed" => TransomClosure::Fixed {
                length: param.unwrap_or(0.3).max(0.0),
            },
            "ballistic" => match param {
                Some(c) => TransomClosure::Ballistic { coeff: c.max(0.0) },
                None => TransomClosure::default(),
            },
            other => {
                return Err(format!(
                    "closure {other:?}: expected ballistic, fixed or off"
                ))
            }
        };
        let grid = num("grid")?.map_or(640, |g| (g as usize).clamp(40, 1200));
        let dynamic = match get("attitude").unwrap_or("dynamic") {
            "dynamic" => true,
            "design" => false,
            other => return Err(format!("attitude {other:?}: expected dynamic or design")),
        };
        let mass = num("mass")?;
        if mass.is_some_and(|m| !(m > 0.0)) {
            return Err("mass must be positive".into());
        }
        let warm = match get("warm").filter(|v| !v.is_empty()) {
            None => None,
            Some(v) => {
                let (a, b) = v.split_once(',').ok_or("warm: expected sinkage,trim")?;
                let p = |t: &str| {
                    t.trim()
                        .parse::<f64>()
                        .map_err(|_| format!("warm: bad number {t:?}"))
                };
                let (s, t) = (p(a)?, p(b)?);
                (s.is_finite() && t.is_finite() && t.abs() < 0.3).then_some((s, t))
            }
        };
        Ok(FlowRequest {
            cut: LoftRequest::from_query(pairs)?,
            froude,
            closure,
            grid,
            dynamic,
            mass,
            lcg: num("lcg")?,
            hold: None,
            mass_by: match get("mass_by").unwrap_or("sinking") {
                "sinking" => MassBy::Sinking,
                "scale" => MassBy::ScaleXyz,
                "scale_yz" => MassBy::ScaleYz,
                other => {
                    return Err(format!(
                        "mass_by {other:?}: expected sinking, scale or scale_yz"
                    ))
                }
            },
            span: match num("span")? {
                Some(s) if !(s > 0.0 && s.is_finite()) => {
                    return Err(format!("span {s}: expected a positive centre span"))
                }
                s => s,
            },
            warm,
            seakeeping: SeakeepingRequest::from_query(pairs)?,
        })
    }
}

/// Where a [`flow_with_progress`] computation has got to.
#[derive(Debug, Clone)]
pub struct Progress {
    /// The stage under way, e.g. "Floating at speed".
    pub stage: &'static str,
    /// Stage number (1-based) and how many there are.
    pub step: usize,
    pub steps: usize,
    /// What the stage is doing, e.g. the force evaluation count.
    pub detail: String,
    /// The whole computation's estimated fraction done, `0..=1`. Stages are
    /// weighted by their typical cost; the equilibrium's iteration count is
    /// not known in advance, so its share approaches its end asymptotically.
    pub fraction: f64,
}

/// The stages of a flow and their typical shares of its time: floating the
/// hull at speed dominates when it is solved; otherwise the free surface
/// and the dynamic force (then computed at the forces stage) do.
fn flow_stages(dynamic: bool, seakeeping: Option<&SeakeepingRequest>) -> [(&'static str, f64); 6] {
    // Measured on e12 (Fn 0.3, dynamic: 0.04 / 8.4 / 0.6 / 3.6 / 7.2 s):
    // the forces stage also encodes the hull meshes into the answer.
    // Seakeeping: ~7 s per heading for 21 head-sea wavelengths, as much
    // again for an irregular sea.
    let sk = seakeeping.map_or(0.0, |s| {
        let sea = if s.sea.is_some() { 1.0 } else { 0.0 };
        0.4 * s.headings.len() as f64 * (s.lambdas.len() as f64 / 21.0 + sea)
    });
    let [a, b, c, d, e] = if dynamic {
        [
            ("Cutting sections", 0.01),
            ("Floating at speed", 0.5),
            ("Hull pressure", 0.04),
            ("Free surface", 0.2),
            ("Forces", 0.25),
        ]
    } else {
        [
            ("Cutting sections", 0.02),
            ("Floating at speed", 0.0),
            ("Hull pressure", 0.08),
            ("Free surface", 0.35),
            ("Forces", 0.55),
        ]
    };
    [a, b, c, d, e, ("Seakeeping", sk)]
}

/// Reports progress through the caller's callback, and turns its refusal
/// (the client has gone) into an error that stops the computation.
struct Tracker<'a> {
    stages: [(&'static str, f64); 6],
    report: &'a mut dyn FnMut(&Progress) -> bool,
}

/// The error a flow stops with when its progress callback asks it to.
pub const CANCELLED: &str = "cancelled";

impl Tracker<'_> {
    /// Stage `i` (0-based) is `within` (0..1) of the way through.
    fn at(&mut self, i: usize, within: f64, detail: String) -> Result<(), String> {
        let total: f64 = self.stages.iter().map(|s| s.1).sum();
        let before: f64 = self.stages[..i].iter().map(|s| s.1).sum();
        let p = Progress {
            stage: self.stages[i].0,
            step: i + 1,
            steps: self.stages.len(),
            detail,
            fraction: ((before + self.stages[i].1 * within.clamp(0.0, 1.0)) / total).min(1.0),
        };
        if (self.report)(&p) {
            Ok(())
        } else {
            Err(CANCELLED.into())
        }
    }
}

/// The steady flow at one speed: each hull's near-field pressure, the free
/// surface around the fleet (local field and waves), and the forces.
pub fn flow(name: &str, bytes: Vec<u8>, req: &FlowRequest) -> Result<Value, String> {
    flow_with_progress(name, bytes, req, &mut |_| true)
}

/// [`flow`], reporting each stage (and each force evaluation of the
/// equilibrium solve) through `report`; returning `false` from it stops the
/// computation with the error [`CANCELLED`].
pub fn flow_with_progress(
    name: &str,
    bytes: Vec<u8>,
    req: &FlowRequest,
    report: &mut dyn FnMut(&Progress) -> bool,
) -> Result<Value, String> {
    flow_inner(name, bytes, req, report, None)
}

/// One demihull's free surface at a held attitude, computed once for a fast
/// span sweep and summed at each span's offsets.
struct SweepField {
    /// The largest span the half-field must cover.
    max_span: f64,
    /// `(dy, rows, grid)`: the half-field (y ≥ 0) on rows `dy` apart.
    field: Option<(f64, usize, michell::WaveGrid)>,
}

fn flow_inner(
    name: &str,
    bytes: Vec<u8>,
    req: &FlowRequest,
    report: &mut dyn FnMut(&Progress) -> bool,
    mut shared: Option<&mut SweepField>,
) -> Result<Value, String> {
    use michell::nearfield::{free_surface, hull_pressure, NearFieldOptions};
    use michell_geometry::float::{solve_equilibrium_sectional_dynamic, LoadCase};
    use michell_geometry::source::SourceHull;
    let t0 = std::time::Instant::now();
    let mut track = Tracker {
        stages: flow_stages(req.dynamic, req.seakeeping.as_ref()),
        report,
    };
    track.at(0, 0.0, name.to_string())?;
    let cut = cut(name, bytes, &req.cut)?;
    track.at(
        0,
        1.0,
        format!(
            "{} hull{}",
            cut.hulls.len(),
            if cut.hulls.len() == 1 { "" } else { "s" }
        ),
    )?;
    // The fleet: the file's hulls where they are, or its one hull doubled
    // into a catamaran — each copy as `(hull index, pose)`, the pose's `dy`
    // putting the copies' centreplanes at ±span/2.
    let layout: Vec<(usize, HullPose)> = match req.span {
        None => cut
            .index
            .iter()
            .map(|&i| (i, HullPose::default()))
            .collect(),
        Some(span) => {
            if cut.hulls.len() != 1 {
                return Err(format!(
                    "a catamaran doubles a single hull; this file holds {}",
                    cut.hulls.len()
                ));
            }
            let (yc, beam) = (
                cut.hulls[0].placement.y,
                michell_cli::fleet::max_beam(&cut.hulls[0].hull),
            );
            if span <= beam {
                return Err(format!(
                    "span {span} m: the demihulls overlap (beam {beam:.3} m)"
                ));
            }
            [0.5 * span, -0.5 * span]
                .map(|y| {
                    (
                        cut.index[0],
                        HullPose {
                            dy: y - yc,
                            ..HullPose::default()
                        },
                    )
                })
                .to_vec()
        }
    };
    let design_owned: Vec<(SectionalHull, Placement)> = layout
        .iter()
        .map(|(i, pose)| {
            let h = &cut.hulls[cut.index.iter().position(|j| j == i).expect("a cut hull")];
            (
                h.hull.clone(),
                Placement {
                    x: h.placement.x,
                    y: h.placement.y + pose.dy,
                },
            )
        })
        .collect();
    let design: Vec<(&SectionalHull, Placement)> =
        design_owned.iter().map(|(h, p)| (h, *p)).collect();
    let l_ref = design
        .iter()
        .map(|(h, _)| h.length())
        .fold(0.0f64, f64::max);
    let cond = Conditions::seawater(req.froude * (STANDARD_GRAVITY * l_ref).sqrt());
    let rho = cond.fluid.density;
    let wave = WaveOptions {
        transom: req.closure,
        ..WaveOptions::default()
    };
    let squat = michell::squat::SquatOptions {
        wave,
        ..Default::default()
    };
    // The load: by default the design displacement at its LCB, so at rest
    // the fleet floats exactly at its design waterline.
    let vol: f64 = design.iter().map(|(h, _)| h.displaced_volume()).sum();
    let mass = req.mass.unwrap_or(rho * vol);
    // A mass carried by scaling: every hull scaled alike until the fleet
    // displaces it at the design waterline, and re-cut there — the design
    // data (and the LCB the LCG defaults to) are the scaled hull's.
    let mut layout = layout;
    let design_owned = match (req.mass, req.mass_by) {
        (Some(m), by) if by != MassBy::Sinking => {
            let ratio = m / (rho * vol);
            for (_, pose) in layout.iter_mut() {
                match by {
                    MassBy::ScaleXyz => pose.scale = ratio.cbrt(),
                    _ => pose.scale_yz = ratio.sqrt(),
                }
            }
            let opts = &cut.opts;
            layout
                .iter()
                .map(|(i, pose)| {
                    cut.file
                        .source
                        .situate_sectional(
                            *i,
                            cut.file.waterline_z,
                            pose,
                            &Platform::default(),
                            opts,
                        )
                        .map_err(|e| e.to_string())?
                        .map(|h| (h.hull, h.placement))
                        .ok_or_else(|| "the scaled hull is dry".to_string())
                })
                .collect::<Result<Vec<_>, String>>()?
        }
        _ => design_owned,
    };
    let design: Vec<(&SectionalHull, Placement)> =
        design_owned.iter().map(|(h, p)| (h, *p)).collect();
    let vol: f64 = design.iter().map(|(h, _)| h.displaced_volume()).sum();
    let lcb = design
        .iter()
        .map(|(h, pl)| h.displaced_volume() * (h.lcb_x() + pl.x))
        .sum::<f64>()
        / vol.max(f64::MIN_POSITIVE);
    let lcg = req.lcg.unwrap_or(lcb);

    // The attitude: the dynamic equilibrium at this speed, or the design one.
    let mut platform = Platform::default();
    let mut solved = None;
    let at_attitude: Vec<(SectionalHull, Placement)> = if let Some((sinkage, trim)) = req.hold {
        // Held: cut the fleet at the given attitude, no solve.
        platform = Platform {
            sinkage,
            trim,
            pivot_x: lcg,
        };
        track.at(1, 0.5, "at the held attitude".into())?;
        let mut out: Vec<(SectionalHull, Placement)> = Vec::new();
        for (k, (i, pose)) in layout.iter().enumerate() {
            // A copy of an earlier layout entry, shifted sideways, is that cut.
            if let Some(j) = layout[..k].iter().position(|(o, q)| {
                o == i && HullPose { dy: 0.0, ..*q } == HullPose { dy: 0.0, ..*pose }
            }) {
                let (h, p) = out[j].clone();
                out.push((
                    h,
                    Placement {
                        y: p.y + pose.dy - layout[j].1.dy,
                        ..p
                    },
                ));
                continue;
            }
            let h = cut
                .file
                .source
                .situate_sectional(*i, cut.file.waterline_z, pose, &platform, &cut.opts)
                .map_err(|e| e.to_string())?
                .ok_or("the hull is dry at the held attitude")?;
            out.push((h.hull, h.placement));
        }
        out
    } else if req.dynamic {
        let sources: Vec<SourceHull> = layout
            .iter()
            .map(|&(index, pose)| SourceHull {
                source: cut.file.source.as_ref(),
                index,
                waterline_z: cut.file.waterline_z,
                pose,
            })
            .collect();
        track.at(1, 0.0, "starting the Newton solve".into())?;
        // Count the solve's force evaluations (each is most of an
        // iteration's cost) and report the latest lift; the count is open
        // ended, so the stage's share fills asymptotically.
        let inner = michell::sectional::dynamic_load_closure(&cond, lcg, &squat);
        let mut evaluations = 0usize;
        let weight = mass * cond.gravity;
        let track_ref = &mut track;
        let report = move |d: &michell_geometry::float::DynamicLoad, vol: f64| {
            evaluations += 1;
            track_ref
                .at(
                    1,
                    1.0 - (-(evaluations as f64) / 10.0).exp(),
                    format!(
                        "force evaluation {evaluations} · lift {:.1}% of weight · \
                         displacing {:.1}% of the load",
                        100.0 * d.force_up / weight,
                        100.0 * rho * vol / mass,
                    ),
                )
                .map_err(michell::Error::InvalidInput)
        };
        let counted = Counted { inner, report };
        let eq = solve_equilibrium_sectional_dynamic(
            &sources,
            &LoadCase {
                mass,
                lcg: Some(lcg),
            },
            rho,
            cond.gravity,
            &cut.opts,
            counted,
            req.warm,
        )
        .map_err(|e| {
            if e.to_string().contains(CANCELLED) {
                return CANCELLED.to_string();
            }
            let hint = if e.to_string().contains("waterplane") {
                " — if the geometry ends at the waterline (no topsides), it cannot sink or \
                 trim; choose the design-waterline attitude"
            } else {
                ""
            };
            format!("dynamic equilibrium at Fn {}: {e}{hint}", req.froude)
        })?;
        platform = Platform {
            sinkage: eq.sinkage,
            trim: eq.trim,
            pivot_x: lcg,
        };
        solved = Some((
            eq.sinkage,
            eq.trim,
            eq.iterations,
            eq.dynamic,
            eq.lift_fraction,
        ));
        eq.fleet.members
    } else {
        design_owned.clone()
    };
    let members: Vec<(&SectionalHull, Placement)> =
        at_attitude.iter().map(|(h, p)| (h, *p)).collect();
    let t_attitude = t0.elapsed().as_secs_f64();

    let nf = NearFieldOptions {
        closure: req.closure,
        ..NearFieldOptions::default()
    };
    track.at(2, 0.0, format!("{} rows deep per hull", nf.depth_rows))?;
    let pressures = hull_pressure(&members, &cond, &nf).map_err(|e| e.to_string())?;

    // The free surface: from ahead of the bows to ~1.5 lengths astern.
    let (mut xa, mut xb, mut yh) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
    for (h, pl) in &members {
        let (a, b) = h.x_range();
        xa = xa.min(a + pl.x);
        xb = xb.max(b + pl.x);
        yh = yh.max(pl.y.abs() + 0.5 * michell_cli::fleet::max_beam(h));
    }
    let (x0, x1) = (xa - 1.5 * l_ref, xb + 0.4 * l_ref);
    let yh = (yh + 0.45 * l_ref).max(0.3 * (x1 - x0));
    let nx = req.grid;
    // An odd row count, so y = 0 is a row and a symmetric fleet mirrors.
    let nh = ((nx as f64) * yh / (x1 - x0)).round().clamp(8.0, 600.0) as usize + 1;
    let ny = 2 * nh - 1;
    track.at(
        3,
        0.0,
        format!("{nx} × {ny} grid over {:.0} × {:.0} m", x1 - x0, 2.0 * yh),
    )?;
    // A fleet symmetric about y = 0 — one hull on the centreline, or
    // mirrored pairs — has a symmetric field: compute the half with y ≥ 0.
    let symmetric = members.iter().all(|(h, pl)| {
        let tol = 1e-6 * l_ref;
        pl.y.abs() < tol
            || members.iter().any(|(o, q)| {
                (q.y + pl.y).abs() < tol
                    && (q.x - pl.x).abs() < tol
                    && (o.displaced_volume() - h.displaced_volume()).abs()
                        <= 1e-6 * h.displaced_volume()
                    && (o.length() - h.length()).abs() < tol
            })
    });
    let g = if let (Some(span), Some(sf)) = (req.span, shared.as_deref_mut()) {
        // A fast span sweep: one demihull half-field at the held attitude,
        // wide enough for the largest span, summed at this span's offsets
        // with linear interpolation between its rows (display only).
        let dy = (x1 - x0) / (nx - 1) as f64;
        if sf.field.is_none() {
            let reach = yh.max(0.0) + 0.5 * sf.max_span + 0.6 * l_ref;
            let rows = (reach / dy).ceil() as usize + 2;
            let (h0, p0) = members[0];
            let one = [(h0, Placement { x: p0.x, y: 0.0 })];
            let half = free_surface(
                &one,
                &cond,
                &nf,
                x0,
                x1,
                0.0,
                (rows - 1) as f64 * dy,
                nx,
                rows,
            )
            .map_err(|e| e.to_string())?;
            sf.field = Some((dy, rows, half));
        }
        let (dy1, rows, half) = sf.field.as_ref().expect("just computed");
        let at = |y: f64, i: usize| -> f64 {
            let u = (y.abs() / dy1).min((rows - 1) as f64);
            let j = (u.floor() as usize).min(rows - 2);
            let t = u - j as f64;
            (1.0 - t) * half.zeta[j * nx + i] + t * half.zeta[(j + 1) * nx + i]
        };
        let nh = (yh / dy).ceil() as usize + 1;
        let (yh, ny) = ((nh - 1) as f64 * dy, 2 * nh - 1);
        let mut zeta = Vec::with_capacity(nx * ny);
        for iy in 0..ny {
            let y = (iy.abs_diff(nh - 1)) as f64 * dy;
            for i in 0..nx {
                zeta.push(at(y - 0.5 * span, i) + at(y + 0.5 * span, i));
            }
        }
        michell::WaveGrid {
            y0: -yh,
            y1: yh,
            ny,
            zeta,
            ..half.clone()
        }
    } else if let Some(span) = req.span {
        // A catamaran's field is its demihull's, twice, shifted by ±span/2
        // (exact in thin-ship theory): one hull's field on a y-grid whose
        // spacing divides span/2, summed at offset rows. Both demihulls are
        // the same cut at the same attitude, and each is symmetric about its
        // own centreplane, so the sum is symmetric too — half is computed.
        let dxg = (x1 - x0) / (nx - 1) as f64;
        let m = ((0.5 * span) / dxg).round().max(1.0) as usize;
        let dy = 0.5 * span / m as f64;
        let nh = (yh / dy).ceil() as usize + 1;
        let (yh, ny) = ((nh - 1) as f64 * dy, 2 * nh - 1);
        let rows = nh + m;
        let (h0, p0) = members[0];
        let one = [(h0, Placement { x: p0.x, y: 0.0 })];
        let half = free_surface(
            &one,
            &cond,
            &nf,
            x0,
            x1,
            0.0,
            (rows - 1) as f64 * dy,
            nx,
            rows,
        )
        .map_err(|e| e.to_string())?;
        let row = |j: usize| &half.zeta[j * nx..(j + 1) * nx];
        let mut zeta = Vec::with_capacity(nx * ny);
        for iy in 0..ny {
            let k = iy.abs_diff(nh - 1);
            let (a, b) = (row(k.abs_diff(m)), row(k + m));
            zeta.extend(a.iter().zip(b).map(|(a, b)| a + b));
        }
        michell::WaveGrid {
            y0: -yh,
            y1: yh,
            ny,
            zeta,
            ..half
        }
    } else if symmetric {
        let half = free_surface(&members, &cond, &nf, x0, x1, 0.0, yh, nx, nh)
            .map_err(|e| e.to_string())?;
        let mut zeta = Vec::with_capacity(nx * ny);
        for iy in 0..ny {
            let k = iy.abs_diff(nh - 1);
            zeta.extend_from_slice(&half.zeta[k * nx..(k + 1) * nx]);
        }
        michell::WaveGrid {
            y0: -yh,
            ny,
            zeta,
            ..half
        }
    } else {
        free_surface(&members, &cond, &nf, x0, x1, -yh, yh, nx, ny).map_err(|e| e.to_string())?
    };
    track.at(
        4,
        0.0,
        if req.dynamic {
            "resistance, then the hull meshes".into()
        } else {
            "resistance, dynamic force, then the hull meshes".into()
        },
    )?;
    let t_field = t0.elapsed().as_secs_f64() - t_attitude;

    let res = michell::sectional::multihull_resistance(
        &members,
        &cond,
        &wave,
        &michell::ViscousOptions::default(),
    )
    .map_err(|e| e.to_string())?;

    // Forces: the solved balance, or at the design attitude the force and
    // the first-order sinkage and trim it implies.
    let forces = match (solved, req.hold) {
        (None, Some((sinkage, trim))) => {
            // Held, not balanced: the dynamic load at this attitude shows
            // how far from its own equilibrium this span sits. An indication,
            // so from the quadrature's first passes (~0.1%), not its last.
            let quick = michell::squat::SquatOptions {
                max_refinements: 1,
                ..squat
            };
            let d = michell::sectional::multihull_dynamic_force(&members, &cond, lcg, &quick)
                .map_err(|e| e.to_string())?;
            json!({
                "fz": d.force_up,
                "moment": d.moment_bow_up,
                // Of the weight, as the solved cases report it.
                "lift_fraction": d.force_up / (mass * cond.gravity),
                "sinkage": sinkage,
                "trim_deg": trim.to_degrees(),
                "trim_rad": trim,
                "solved": false,
                "held": true,
            })
        }
        (s, _) => match s {
            Some((sinkage, trim, iterations, d, lift)) => json!({
                "fz": d.force_up,
                "moment": d.moment_bow_up,
                "lift_fraction": lift,
                "sinkage": sinkage,
                "trim_deg": trim.to_degrees(),
                "trim_rad": trim,
                "iterations": iterations,
                "solved": true,
            }),
            None => {
                let (mut aw, mut mw, mut iw) = (0.0, 0.0, 0.0);
                for (h, pl) in &members {
                    let a = h.waterplane_area();
                    aw += a;
                    mw += h.waterplane_moment() + pl.x * a;
                    iw += h.waterplane_second_moment()
                        + 2.0 * pl.x * h.waterplane_moment()
                        + pl.x * pl.x * a;
                }
                let lcf = mw / aw.max(f64::MIN_POSITIVE);
                let i_l = iw - mw * mw / aw.max(f64::MIN_POSITIVE);
                let d = michell::sectional::multihull_dynamic_force(&members, &cond, lcf, &squat)
                    .map_err(|e| e.to_string())?;
                json!({
                    "fz": d.force_up,
                    "moment": d.moment_bow_up,
                    "lift_fraction": d.lift_fraction,
                    "sinkage": -d.force_up / (rho * cond.gravity * aw),
                    "trim_deg": (d.moment_bow_up / (rho * cond.gravity * i_l)).to_degrees(),
                    "solved": false,
                })
            }
        },
    };

    // The whole hull at the attitude, topsides included, for display: the
    // source's tessellation, x forward, y across, z up from the water.
    let mut meshes = Vec::new();
    for (i, pose) in &layout {
        let (v, t) = cut
            .file
            .source
            .posed_tessellation(*i, cut.file.waterline_z, pose, &platform)
            .map_err(|e| e.to_string())?;
        let flat: Vec<f32> = v.iter().flat_map(|p| p.map(|c| c as f32)).collect();
        let idx: Vec<u32> = t.iter().flatten().copied().collect();
        meshes.push(json!({
            "vertices": b64(bytemuck_f32(&flat)),
            "triangles": b64(bytemuck_u32(&idx)),
        }));
    }

    let hulls: Vec<Value> = pressures
        .iter()
        .zip(meshes)
        .map(|(p, mesh)| {
            json!({
                "x": round(&p.x),
                "depth": round(&p.depth),
                "half_beam": round(&p.half_beam),
                "cp": round(&p.cp),
                "y": p.y,
                "force_up": p.force_up,
                "mesh": mesh,
            })
        })
        .collect();
    let mut forces = forces;
    for (k, v) in [
        ("rw", res.wave.resistance),
        ("rv", res.viscous_total),
        ("rt", res.total),
        ("pe", res.effective_power),
        ("cw", res.cw),
        ("ct", res.ct),
        ("interference", res.interference),
        ("mass", mass),
        ("lcg", lcg),
        (
            "displaced_volume",
            members.iter().map(|(h, _)| h.displaced_volume()).sum(),
        ),
        (
            "length",
            members.iter().map(|(h, _)| h.length()).fold(0.0, f64::max),
        ),
    ] {
        forces[k] = json!(v);
    }
    let seakeeping = match &req.seakeeping {
        Some(sk) => Some(seakeeping_json(
            &members, mass, lcg, cond.speed, sk, &mut track,
        )?),
        None => None,
    };
    Ok(json!({
        "seakeeping": seakeeping,
        "froude": req.froude,
        "speed": cond.speed,
        "transverse_wavelength": 2.0 * std::f64::consts::PI * cond.speed * cond.speed / cond.gravity,
        "seconds": t0.elapsed().as_secs_f64(),
        "timing": { "attitude": t_attitude, "field": t_field },
        "hulls": hulls,
        "surface": {
            "x0": g.x0, "x1": g.x1, "y0": g.y0, "y1": g.y1, "nx": g.nx, "ny": g.ny,
            "zeta": b64(bytemuck_f32(&g.zeta.iter().map(|&z| z as f32).collect::<Vec<_>>())),
        },
        "forces": forces,
    }))
}

/// The seakeeping block of a flow case: roll stability, and per heading the
/// responses over the wavelengths and the irregular sea's statistics.
fn seakeeping_json(
    members: &[(&SectionalHull, Placement)],
    mass: f64,
    lcg: f64,
    speed: f64,
    sk: &SeakeepingRequest,
    track: &mut Tracker,
) -> Result<Value, String> {
    use michell_seakeeping::platform::{self, Loading};
    use michell_seakeeping::sea::sea_response_fleet;
    use michell_seakeeping::strip::StripOptions;
    let cond = Conditions::seawater(speed);
    let opts = StripOptions {
        density: cond.fluid.density,
        gravity: cond.gravity,
        roll_damping: sk.roll_damping,
        ..StripOptions::default()
    };
    let l = platform::length(members);
    let props = platform::mass_properties(
        members,
        opts.density,
        &Loading {
            mass: Some(mass),
            lcg: Some(lcg),
            vcg: sk.vcg,
            k_yy: sk.kyy.map(|f| f * l),
            k_xx: sk.kxx,
            k_zz: sk.kzz,
        },
    );
    track.at(5, 0.0, "roll stability".into())?;
    let roll = platform::roll_stability(members, &props, &opts);
    // Progress: each heading's sweep, then its sea, in equal shares.
    let n = sk.lambdas.len();
    let per = n + if sk.sea.is_some() { n } else { 0 };
    let total = (sk.headings.len() * per).max(1) as f64;
    let x_bow = members
        .iter()
        .map(|(h, pl)| h.x_range().1 + pl.x)
        .fold(f64::NEG_INFINITY, f64::max);
    let c = |z: michell_geometry::C64| json!([z.re, z.im]);
    let mut headings = Vec::new();
    for (j, &deg) in sk.headings.iter().enumerate() {
        let done = (j * per) as f64;
        // The far field's sections carry only the symmetric part of the
        // diffraction problem: head and following seas only.
        let far = |v: f64| (deg.to_radians().sin().abs() < 0.1).then_some(v);
        let mut stopped = false;
        let points = platform::rao_sweep(
            members,
            &props,
            deg.to_radians(),
            speed,
            &sk.lambdas,
            &opts,
            &mut |i| {
                stopped = track
                    .at(
                        5,
                        (done + i as f64 + 1.0) / total,
                        format!("heading {deg:.0}° · λ/L {:.2}", sk.lambdas[i]),
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
            .zip(&sk.lambdas)
            .map(|(p, lam)| match p {
                Ok(p) => json!({
                    "lambda": p.lambda_over_l,
                    "omega": p.omega,
                    "omega_e": p.omega_e,
                    "k": p.k,
                    "heave": c(p.heave),
                    "pitch": c(p.pitch),
                    "sway": c(p.sway),
                    "roll": c(p.roll),
                    "yaw": c(p.yaw),
                    "raw_gb": p.added_resistance,
                    "raw_far": far(p.added_resistance_far_field),
                }),
                Err(e) => json!({ "lambda": lam, "error": e }),
            })
            .collect();
        let sea = match &sk.sea {
            Some(spectrum) => {
                track.at(
                    5,
                    (done + n as f64) / total,
                    format!("heading {deg:.0}° · irregular sea"),
                )?;
                match sea_response_fleet(
                    members,
                    &props,
                    spectrum,
                    deg.to_radians(),
                    speed,
                    &[x_bow, lcg],
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
        headings.push(json!({ "heading": deg, "points": points, "sea": sea }));
    }
    track.at(5, 1.0, "done".into())?;
    Ok(json!({
        "length": l,
        "beam": platform::hull_beam(members),
        "mass": props.mass,
        "lcg": props.lcg,
        "bg": props.bg,
        // G: at the LCG, on the centreplane (y = 0), this far above the
        // waterline — the point the motions are about.
        "vcg": sk.vcg.unwrap_or(0.0),
        "k_yy": props.radius_of_gyration,
        "k_xx": props.roll_radius_of_gyration,
        "k_zz": props.yaw_radius_of_gyration,
        "roll_damping": sk.roll_damping,
        "gm_t": roll.gm,
        "roll_period_dry": roll.period_dry,
        "roll_period": roll.period,
        "sea": sk.sea.map(|s| format!("{s:?}")),
        "headings": headings,
    }))
}

fn bytemuck_f32(v: &[f32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

fn bytemuck_u32(v: &[u32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

/// Standard base64, for little-endian arrays the page decodes into typed
/// arrays (a third the size of JSON numbers, and no parsing).
fn b64(bytes: Vec<u8>) -> String {
    const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for c in bytes.chunks(3) {
        let n = (c[0] as u32) << 16
            | (*c.get(1).unwrap_or(&0) as u32) << 8
            | *c.get(2).unwrap_or(&0) as u32;
        for k in 0..4 {
            if k <= c.len() {
                out.push(A[(n >> (18 - 6 * k) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

/// A fast catamaran span sweep: the pair floated once, at the middle span,
/// and every span evaluated at that held attitude — resistance, pressure,
/// the free surface (summed from one demihull field) and the dynamic load
/// the attitude leaves unbalanced. Each span's result goes to `emit` (with
/// its index in `spans`) as it is done; returning `false` from `emit` or
/// `report` stops the sweep.
pub fn span_sweep_with_progress(
    name: &str,
    bytes: Vec<u8>,
    req: &FlowRequest,
    spans: &[f64],
    emit: &mut dyn FnMut(usize, Value) -> bool,
    report: &mut dyn FnMut(&Progress) -> bool,
) -> Result<(), String> {
    if spans.is_empty() || spans.iter().any(|s| !(*s > 0.0)) {
        return Err("spans: expected positive centre spans".into());
    }
    let mut order: Vec<usize> = (0..spans.len()).collect();
    order.sort_by(|&a, &b| spans[a].total_cmp(&spans[b]));
    let mid = order[order.len() / 2];
    // The middle span first (solved), then the rest at its attitude.
    order.retain(|&i| i != mid);
    order.insert(0, mid);
    let mut shared = SweepField {
        max_span: spans.iter().cloned().fold(0.0, f64::max),
        field: None,
    };
    let n = order.len();
    let mut hold = None;
    for (k, &i) in order.iter().enumerate() {
        let sub = FlowRequest {
            span: Some(spans[i]),
            hold,
            cut: LoftRequest { ..req.cut },
            seakeeping: req.seakeeping.clone(),
            ..*req
        };
        // The solved span is most of the work; weight it as such.
        let (lo, w) = if k == 0 {
            (0.0, 0.6)
        } else {
            (
                0.6 + 0.4 * (k - 1) as f64 / (n - 1) as f64,
                0.4 / (n - 1) as f64,
            )
        };
        let mut sub_report = |p: &Progress| {
            report(&Progress {
                detail: format!("span {:.2} m ({}/{n}) · {}", spans[i], k + 1, p.detail),
                fraction: lo + w * p.fraction,
                ..p.clone()
            })
        };
        let v = flow_inner(
            name,
            bytes.clone(),
            &sub,
            &mut sub_report,
            Some(&mut shared),
        )?;
        if hold.is_none() {
            // The solved attitude; at the design attitude, the design one.
            let f = &v["forces"];
            hold = Some(if f["solved"].as_bool() == Some(true) {
                (
                    f["sinkage"].as_f64().unwrap_or(0.0),
                    f["trim_rad"].as_f64().unwrap_or(0.0),
                )
            } else {
                (0.0, 0.0)
            });
        }
        if !emit(i, v) {
            return Err(CANCELLED.into());
        }
    }
    Ok(())
}

/// Values rounded to 6 significant figures: a third of the JSON.
fn round(v: &[f64]) -> Vec<Value> {
    v.iter()
        .map(|&x| {
            if x == 0.0 || !x.is_finite() {
                json!(0.0)
            } else {
                let e = 5 - x.abs().log10().floor() as i32;
                let k = 10f64.powi(e);
                json!((x * k).round() / k)
            }
        })
        .collect()
}

fn sectioned_json(imp: &SectionalImport, kind: Kind) -> Value {
    let r = &imp.report;
    let sides = if r.two_sided {
        format!("both sides averaged about y = {:.4} m", r.centerplane)
    } else {
        format!("one side about y = {:.4} m", r.centerplane)
    };
    let what = match kind {
        Kind::Iges => format!("IGES, {} patches", r.patches),
        Kind::Stl => format!("STL, {} triangles", r.patches),
    };
    let mut lines = vec![format!(
        "source: {what}, cut at {} stations over x {:.4}..{:.4} m ({sides})",
        r.stations, r.x_range.0, r.x_range.1
    )];
    if r.dropped_stations > 0 {
        lines.push(format!(
            "{} interior stations could not be sectioned (bridged by the interpolant)",
            r.dropped_stations
        ));
    }
    if r.max_asymmetry > 1e-3 * r.draft.max(1e-9) {
        lines.push(format!(
            "port and starboard differ by up to {:.3e} m (averaged)",
            r.max_asymmetry
        ));
    }
    hull_json(&imp.hull, imp.placement, Some(&imp.sections), lines)
}

/// One hull, drawn as the physics uses it (see the module docs).
fn hull_json(
    hull: &SectionalHull,
    placement: Placement,
    rays: Option<&Vec<(f64, Vec<(f64, f64)>)>>,
    lines: Vec<String>,
) -> Value {
    // Each station as the curve the quadrature integrates (not a polyline
    // through its quadrature nodes, which are graded toward the waterline).
    let stations: Vec<(f64, Vec<(f64, f64)>)> =
        hull.curves().map(|(x, o)| (x, o.to_vec())).collect();
    // A see-through surface between stations, for orientation only (the
    // kernel interpolates each station's depth integral along x, not a
    // surface): rows joined at equal fractions of girth, so neighbouring
    // sections with different node spacing still meet cleanly.
    let rows = 48usize;
    let (mut mx, mut mz, mut my) = (Vec::new(), Vec::new(), Vec::new());
    for (x, o) in &stations {
        for (y, z) in by_girth(o, rows) {
            mx.push(*x);
            mz.push(z);
            my.push(y);
        }
    }
    let keel: Vec<(f64, f64)> = stations
        .iter()
        .map(|(x, o)| (*x, o.last().map_or(0.0, |p| p.1)))
        .collect();
    let beam = 2.0
        * stations
            .iter()
            .flat_map(|(_, o)| o.first().map(|p| p.0))
            .fold(0.0, f64::max);
    // The depth-integral curves at κ = 0 (sectional area) and at a short
    // wave's decay rate: λ = 2 at Fn 0.15, κ = νλ², ν = g/U² = 1/(0.0225 L).
    let kappas = [0.0, 4.0 / (0.0225 * hull.length())];
    let curves: Vec<Value> = kappas
        .iter()
        .map(|&kappa| {
            let (st, c) = hull.depth_integral_curve(kappa, 8);
            json!({ "kappa": kappa, "stations": st, "curve": c })
        })
        .collect();
    // The transom and its section — the aft end station's outline, whose
    // depth integral is the closure's depth factor — so the page can draw
    // the virtual appendage at whatever speed and closure it is asked about.
    let transom = hull.transom().map(|t| {
        let aft = stations
            .iter()
            .min_by(|a, b| (a.0 - t.x).abs().total_cmp(&(b.0 - t.x).abs()));
        json!({
            "x": t.x,
            "depth": t.depth,
            "half_beam": t.half_beam,
            "area": t.area,
            "area_ratio": t.area / hull.max_section_area().max(1e-300),
            "outline": aft.map(|s| by_girth(&s.1, rows)).unwrap_or_default(),
            // Z_T at each curve's κ: the kernel extends Z over the hollow as
            // Z_T · φ(s).
            "z_t": kappas
                .iter()
                .map(|&k| hull.depth_integral_curve(k, 1).0.first().map_or(0.0, |p| p.1))
                .collect::<Vec<_>>(),
        })
    });
    json!({
        "placement": { "x": placement.x, "y": placement.y },
        "diagnostics": lines,
        "length": hull.length(),
        "beam": beam,
        "draft": hull.draft(),
        "displaced_volume": hull.displaced_volume(),
        "wetted_surface": hull.wetted_surface(),
        "lcb_x": hull.lcb_x(),
        "waterplane_area": hull.waterplane_area(),
        "stations": stations,
        "rays": rays,
        "mesh": { "nx": stations.len(), "nz": rows + 1, "x": mx, "z": mz, "y": my },
        "keel": keel,
        "area": curves,
        "transom": transom,
    })
}

/// A section outline resampled at `rows + 1` points evenly spaced in girth
/// (a point for an empty section).
fn by_girth(o: &[(f64, f64)], rows: usize) -> Vec<(f64, f64)> {
    if o.len() < 2 {
        return vec![(0.0, 0.0); rows + 1];
    }
    let mut cum = vec![0.0];
    for w in o.windows(2) {
        let d = ((w[1].0 - w[0].0).powi(2) + (w[1].1 - w[0].1).powi(2)).sqrt();
        cum.push(cum.last().unwrap() + d);
    }
    let total = *cum.last().unwrap();
    (0..=rows)
        .map(|k| {
            let g = total * k as f64 / rows as f64;
            let j = cum.partition_point(|&c| c < g).clamp(1, o.len() - 1);
            let t = if cum[j] > cum[j - 1] {
                (g - cum[j - 1]) / (cum[j] - cum[j - 1])
            } else {
                0.0
            };
            (
                o[j - 1].0 + t * (o[j].0 - o[j - 1].0),
                o[j - 1].1 + t * (o[j].1 - o[j - 1].1),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The exact Wigley as IGES, shown by sections.
    #[test]
    fn a_wigley_is_shown_in_sections() {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell_geometry::iges::write(&surfaces, "wigley").unwrap();
        let v = loft("w.igs", text.into_bytes(), &LoftRequest::default()).unwrap();
        let h = &v["hulls"][0];
        let vol = h["displaced_volume"].as_f64().unwrap();
        // The Wigley's exact volume, 4/9 L B T.
        let exact = 4.0 / 9.0 * 10.0 * 0.625;
        assert!((vol - exact).abs() < 1e-6 * exact, "{vol}");
        assert!(h["transom"].is_null());
        assert!(!h["stations"].as_array().unwrap().is_empty());
    }

    /// The flow at one speed: pressure per hull, the surface grid asked for,
    /// and a dynamic lift that sinks the hull.
    #[test]
    fn a_wigley_flow_has_pressure_waves_and_lift() {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell_geometry::iges::write(&surfaces, "wigley").unwrap();
        // This Wigley ends at its waterline (no topsides to sink into), so
        // at the design attitude.
        let pairs: Vec<(String, String)> =
            [("froude", "0.35"), ("grid", "60"), ("attitude", "design")]
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect();
        let req = FlowRequest::from_query(&pairs).unwrap();
        let v = flow("w.igs", text.into_bytes(), &req).unwrap();
        let h = &v["hulls"][0];
        let (nx, nz) = (
            h["x"].as_array().unwrap().len(),
            h["depth"].as_array().unwrap().len(),
        );
        assert_eq!(h["cp"].as_array().unwrap().len(), nx * nz);
        assert_eq!(v["surface"]["nx"], 60);
        let f = &v["forces"];
        assert!(
            f["fz"].as_f64().unwrap() < 0.0 && f["sinkage"].as_f64().unwrap() > 0.0,
            "{f}"
        );
        assert!(f["rw"].as_f64().unwrap() > 0.0);
    }

    fn wigley_flow_request() -> (Vec<u8>, FlowRequest) {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell_geometry::iges::write(&surfaces, "wigley").unwrap();
        let pairs: Vec<(String, String)> =
            [("froude", "0.35"), ("grid", "40"), ("attitude", "design")]
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect();
        (text.into_bytes(), FlowRequest::from_query(&pairs).unwrap())
    }

    /// Progress walks every stage in order, its fraction never goes back,
    /// and a refusal from the callback stops the flow.
    #[test]
    fn flow_progress_advances_and_can_be_cancelled() {
        let (bytes, req) = wigley_flow_request();
        let mut seen: Vec<(usize, f64)> = Vec::new();
        flow_with_progress("w.igs", bytes.clone(), &req, &mut |p| {
            seen.push((p.step, p.fraction));
            true
        })
        .unwrap();
        assert!(seen.windows(2).all(|w| w[1].1 >= w[0].1), "{seen:?}");
        assert_eq!(seen.first().map(|s| s.0), Some(1));
        assert_eq!(seen.last().map(|s| s.0), Some(5));
        let mut calls = 0;
        let e = flow_with_progress("w.igs", bytes, &req, &mut |_| {
            calls += 1;
            calls < 3
        })
        .unwrap_err();
        assert_eq!(e, CANCELLED);
        assert_eq!(calls, 3);
    }

    /// A catamaran's free surface, summed from one demihull's field shifted
    /// by ±span/2, is the two-hull field computed directly.
    #[test]
    fn a_catamaran_field_is_its_demihulls_summed() {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell_geometry::iges::write(&surfaces, "wigley").unwrap();
        let pairs: Vec<(String, String)> = [
            ("froude", "0.35"),
            ("grid", "60"),
            ("attitude", "design"),
            ("span", "3.1"),
        ]
        .iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();
        let v = flow(
            "w.igs",
            text.clone().into_bytes(),
            &FlowRequest::from_query(&pairs).unwrap(),
        )
        .unwrap();
        let sfc = &v["surface"];
        let get = |k: &str| sfc[k].as_f64().unwrap();
        let (nx, ny) = (get("nx") as usize, get("ny") as usize);
        let bytes = unb64(sfc["zeta"].as_str().unwrap());
        let summed: Vec<f64> = bytes
            .chunks(4)
            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]) as f64)
            .collect();
        // The same two hulls, directly, on the same grid.
        let hull = michell_geometry::iges::source_fleet(&text, 0.0)
            .unwrap()
            .situate_sectional(
                0,
                0.0,
                &HullPose::default(),
                &Platform::default(),
                &Default::default(),
            )
            .unwrap()
            .unwrap()
            .hull;
        let fleet = [
            (&hull, Placement { x: 0.0, y: 1.55 }),
            (&hull, Placement { x: 0.0, y: -1.55 }),
        ];
        let cond = Conditions::seawater(0.35 * (STANDARD_GRAVITY * hull.length()).sqrt());
        let direct = michell::nearfield::free_surface(
            &fleet,
            &cond,
            &Default::default(),
            get("x0"),
            get("x1"),
            get("y0"),
            get("y1"),
            nx,
            ny,
        )
        .unwrap();
        let peak = direct.zeta.iter().fold(0.0f64, |m, z| m.max(z.abs()));
        let worst = summed
            .iter()
            .zip(&direct.zeta)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(
            worst < 1e-4 * peak,
            "summed vs direct: {worst} of peak {peak}"
        );
        assert_eq!(v["hulls"].as_array().unwrap().len(), 2);
    }

    fn unb64(s: &str) -> Vec<u8> {
        const A: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let val = |c: u8| A.iter().position(|&a| a == c).unwrap() as u32;
        let mut out = Vec::new();
        for q in s.as_bytes().chunks(4) {
            let n = q.iter().take_while(|&&c| c != b'=').count();
            let v = q[..n].iter().fold(0u32, |acc, &c| acc << 6 | val(c)) << (6 * (4 - n));
            out.extend(&[(v >> 16) as u8, (v >> 8) as u8, v as u8][..n - 1]);
        }
        out
    }

    /// A mass carried by scaling keeps the design waterline: the scaled
    /// Wigley displaces exactly that mass there, its length scaled by the
    /// cube root of the mass ratio (x, y, z) or kept (beam and draft only).
    #[test]
    fn a_mass_carried_by_scaling_keeps_the_waterline() {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell_geometry::iges::write(&surfaces, "wigley").unwrap();
        let rho = Conditions::seawater(1.0).fluid.density;
        let v0 = 4.0 / 9.0 * 10.0 * 0.625;
        let mass = 1.2 * rho * v0;
        for (by, length) in [("scale", 10.0 * 1.2f64.cbrt()), ("scale_yz", 10.0)] {
            let m = format!("{mass}");
            let pairs: Vec<(String, String)> = [
                ("froude", "0.3"),
                ("grid", "40"),
                ("attitude", "design"),
                ("mass", m.as_str()),
                ("mass_by", by),
            ]
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
            let v = flow(
                "w.igs",
                text.clone().into_bytes(),
                &FlowRequest::from_query(&pairs).unwrap(),
            )
            .unwrap();
            let f = &v["forces"];
            let (vol, len) = (
                f["displaced_volume"].as_f64().unwrap(),
                f["length"].as_f64().unwrap(),
            );
            assert!(
                (vol - 1.2 * v0).abs() < 1e-5 * v0,
                "{by}: volume {vol} vs {}",
                1.2 * v0
            );
            assert!(
                (len - length).abs() < 1e-6 * length,
                "{by}: length {len} vs {length}"
            );
        }
    }

    #[test]
    fn seakeeping_keys_parse() {
        let pairs = |kv: &[(&str, &str)]| -> Vec<(String, String)> {
            kv.iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        };
        assert_eq!(
            SeakeepingRequest::from_query(&pairs(&[("sk_vcg", "1")])),
            Ok(None)
        );
        let r = SeakeepingRequest::from_query(&pairs(&[("sk", "1")]))
            .unwrap()
            .unwrap();
        assert_eq!(r.headings, vec![180.0]);
        assert_eq!(r.lambdas.len(), 21);
        assert!(r.sea.is_none() && r.vcg.is_none());
        let r = SeakeepingRequest::from_query(&pairs(&[
            ("sk", "1"),
            ("sk_heading", "180, 90:150:30"),
            ("sk_lambda", "1,2"),
            ("sk_sea", "hs=0.5,tp=4,gamma=3.3"),
            ("sk_damping", "0.05"),
        ]))
        .unwrap()
        .unwrap();
        assert_eq!(r.headings, vec![180.0, 90.0, 120.0, 150.0]);
        assert_eq!(r.lambdas, vec![1.0, 2.0]);
        assert_eq!(r.roll_damping, 0.05);
        assert!(matches!(
            r.sea,
            Some(michell_seakeeping::sea::Spectrum::Jonswap { .. })
        ));
        for bad in [("sk_lambda", "0"), ("sk_sea", "hs=1"), ("sk_vcg", "x")] {
            assert!(
                SeakeepingRequest::from_query(&pairs(&[("sk", "1"), bad])).is_err(),
                "{bad:?}"
            );
        }
    }

    /// A flow case with seakeeping carries its responses: in long head
    /// waves a Wigley hull heaves and pitches with the water, and the far
    /// field is reported only where it applies.
    #[test]
    fn a_flow_case_carries_its_seakeeping() {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell_geometry::iges::write(&surfaces, "wigley").unwrap();
        let pairs: Vec<(String, String)> = [
            ("froude", "0.2"),
            ("grid", "40"),
            ("attitude", "design"),
            ("sk", "1"),
            ("sk_heading", "180,120"),
            ("sk_lambda", "6"),
            // A deep, narrow Wigley hull rolls over with G at the waterline.
            ("sk_vcg", "-0.3"),
        ]
        .iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();
        let v = flow(
            "w.igs",
            text.into_bytes(),
            &FlowRequest::from_query(&pairs).unwrap(),
        )
        .unwrap();
        let sk = &v["seakeeping"];
        assert!(sk["gm_t"].as_f64().unwrap() > 0.0, "{}", sk);
        let heads = sk["headings"].as_array().unwrap();
        assert_eq!(heads.len(), 2);
        let p = &heads[0]["points"][0];
        let abs = |z: &Value| z[0].as_f64().unwrap().hypot(z[1].as_f64().unwrap());
        let k = p["k"].as_f64().unwrap();
        assert!(
            (abs(&p["heave"]) - 1.0).abs() < 0.1,
            "heave {}",
            abs(&p["heave"])
        );
        assert!(
            (abs(&p["pitch"]) / k - 1.0).abs() < 0.15,
            "pitch {}",
            abs(&p["pitch"]) / k
        );
        assert!(p["raw_far"].is_f64());
        assert!(heads[1]["points"][0]["raw_far"].is_null());
    }

    /// A span sweep at the design attitude is, span by span, the single
    /// catamaran flows: the held path and the shared field change nothing
    /// but the cost.
    #[test]
    fn a_span_sweep_is_its_spans_flows() {
        let surfaces = michell_geometry::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell_geometry::iges::write(&surfaces, "wigley").unwrap();
        let q = |extra: &[(&str, &str)]| -> FlowRequest {
            let mut pairs: Vec<(String, String)> =
                [("froude", "0.35"), ("grid", "60"), ("attitude", "design")]
                    .iter()
                    .map(|(k, v)| (k.to_string(), v.to_string()))
                    .collect();
            pairs.extend(extra.iter().map(|(k, v)| (k.to_string(), v.to_string())));
            FlowRequest::from_query(&pairs).unwrap()
        };
        let spans = [2.0, 3.0, 4.5];
        let mut got = vec![None; spans.len()];
        span_sweep_with_progress(
            "w.igs",
            text.clone().into_bytes(),
            &q(&[]),
            &spans,
            &mut |i, v| {
                got[i] = Some(v);
                true
            },
            &mut |_| true,
        )
        .unwrap();
        for (i, s) in spans.iter().enumerate() {
            let one = flow(
                "w.igs",
                text.clone().into_bytes(),
                &q(&[("span", &s.to_string())]),
            )
            .unwrap();
            let v = got[i].as_ref().expect("every span reported");
            let (a, b) = (
                v["forces"]["rt"].as_f64().unwrap(),
                one["forces"]["rt"].as_f64().unwrap(),
            );
            assert!(
                (a - b).abs() < 1e-9 * b,
                "span {s}: sweep {a} vs single {b}"
            );
        }
    }

    #[test]
    fn an_stl_needs_its_units() {
        let e = loft("x.stl", vec![0; 200], &LoftRequest::default()).unwrap_err();
        assert!(e.contains("units"), "{e}");
    }
}

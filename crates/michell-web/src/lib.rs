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

use michell::iges::{HullPose, Platform, SectionalImport};
use michell::sectional::SectionalHull;
use michell::{Conditions, Placement, TransomClosure, WaveOptions, STANDARD_GRAVITY};
use michell_cli::fleet::{open_source_bytes, Kind, LoadSettings};
use serde_json::{json, Value};

/// Largest upload accepted [bytes]. Big enough for a finely tessellated STL.
pub const MAX_UPLOAD: usize = 128 << 20;

/// Import options the page can set.
#[derive(Default)]
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
    let hulls: Vec<Value> = cut.hulls.iter().map(|h| sectioned_json(h, cut.kind)).collect();
    Ok(json!({
        "name": name,
        "seconds": t0.elapsed().as_secs_f64(),
        "notes": cut.notes,
        "hulls": hulls,
    }))
}

/// An upload cut into sections at its design pose, with the source kept so
/// its hulls can be re-cut at another attitude.
struct Cut {
    kind: Kind,
    file: michell_cli::fleet::SourceFile,
    /// Each cut hull's index in the file.
    index: Vec<usize>,
    hulls: Vec<SectionalImport>,
    opts: michell::iges::SectionalOptions,
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
    /// A previous solution `(sinkage [m], trim [rad])` to start the
    /// equilibrium from — the last speed's, as the slider moves.
    pub warm: Option<(f64, f64)>,
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
                .map(|v| v.parse::<f64>().map_err(|_| format!("{k}: expected a number, got {v:?}")))
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
            other => return Err(format!("closure {other:?}: expected ballistic, fixed or off")),
        };
        let grid = num("grid")?.map_or(320, |g| (g as usize).clamp(40, 800));
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
                let p = |t: &str| t.trim().parse::<f64>().map_err(|_| format!("warm: bad number {t:?}"));
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
            warm,
        })
    }
}

/// The steady flow at one speed: each hull's near-field pressure, the free
/// surface around the fleet (local field and waves), and the forces.
pub fn flow(name: &str, bytes: Vec<u8>, req: &FlowRequest) -> Result<Value, String> {
    use michell::float::{solve_equilibrium_sectional_dynamic, LoadCase};
    use michell::nearfield::{free_surface, hull_pressure, NearFieldOptions};
    use michell::source::SourceHull;
    let t0 = std::time::Instant::now();
    let cut = cut(name, bytes, &req.cut)?;
    let design: Vec<(&SectionalHull, Placement)> =
        cut.hulls.iter().map(|h| (&h.hull, h.placement)).collect();
    let l_ref = design.iter().map(|(h, _)| h.length()).fold(0.0f64, f64::max);
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
    let lcb = design
        .iter()
        .map(|(h, pl)| h.displaced_volume() * (h.lcb_x() + pl.x))
        .sum::<f64>()
        / vol.max(f64::MIN_POSITIVE);
    let mass = req.mass.unwrap_or(rho * vol);
    let lcg = req.lcg.unwrap_or(lcb);

    // The attitude: the dynamic equilibrium at this speed, or the design one.
    let mut platform = Platform::default();
    let mut solved = None;
    let at_attitude: Vec<(SectionalHull, Placement)> = if req.dynamic {
        let sources: Vec<SourceHull> = cut
            .index
            .iter()
            .map(|&index| SourceHull {
                source: cut.file.source.as_ref(),
                index,
                waterline_z: cut.file.waterline_z,
                pose: HullPose::default(),
            })
            .collect();
        let eq = solve_equilibrium_sectional_dynamic(
            &sources,
            &LoadCase {
                mass,
                lcg: Some(lcg),
            },
            rho,
            cond.gravity,
            &cut.opts,
            michell::sectional::dynamic_load_closure(&cond, lcg, &squat),
            req.warm,
        )
        .map_err(|e| {
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
        solved = Some((eq.sinkage, eq.trim, eq.iterations, eq.dynamic, eq.lift_fraction));
        eq.fleet.members
    } else {
        cut.hulls.iter().map(|h| (h.hull.clone(), h.placement)).collect()
    };
    let members: Vec<(&SectionalHull, Placement)> =
        at_attitude.iter().map(|(h, p)| (h, *p)).collect();
    let t_attitude = t0.elapsed().as_secs_f64();

    let nf = NearFieldOptions {
        closure: req.closure,
        ..NearFieldOptions::default()
    };
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
    let ny = ((nx as f64) * 2.0 * yh / (x1 - x0)).round().clamp(16.0, 400.0) as usize;
    let g = free_surface(&members, &cond, &nf, x0, x1, -yh, yh, nx, ny).map_err(|e| e.to_string())?;
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
    let forces = match solved {
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
    };

    // The whole hull at the attitude, topsides included, for display: the
    // source's tessellation, x forward, y across, z up from the water.
    let mut meshes = Vec::new();
    for &i in &cut.index {
        let (v, t) = cut
            .file
            .source
            .posed_tessellation(i, cut.file.waterline_z, &HullPose::default(), &platform)
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
    ] {
        forces[k] = json!(v);
    }
    Ok(json!({
        "froude": req.froude,
        "speed": cond.speed,
        "transverse_wavelength": 2.0 * std::f64::consts::PI * cond.speed * cond.speed / cond.gravity,
        "seconds": t0.elapsed().as_secs_f64(),
        "timing": { "attitude": t_attitude, "field": t_field },
        "hulls": hulls,
        "surface": {
            "x0": g.x0, "x1": g.x1, "y0": g.y0, "y1": g.y1, "nx": g.nx, "ny": g.ny,
            "zeta": round(&g.zeta),
        },
        "forces": forces,
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
        let n = (c[0] as u32) << 16 | (*c.get(1).unwrap_or(&0) as u32) << 8 | *c.get(2).unwrap_or(&0) as u32;
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
        let surfaces = michell::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell::iges::write(&surfaces, "wigley").unwrap();
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
        let surfaces = michell::iges::wigley_surfaces(10.0, 1.0, 0.625).unwrap();
        let text = michell::iges::write(&surfaces, "wigley").unwrap();
        // This Wigley ends at its waterline (no topsides to sink into), so
        // at the design attitude.
        let pairs: Vec<(String, String)> = [("froude", "0.35"), ("grid", "60"), ("attitude", "design")]
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        let req = FlowRequest::from_query(&pairs).unwrap();
        let v = flow("w.igs", text.into_bytes(), &req).unwrap();
        let h = &v["hulls"][0];
        let (nx, nz) = (h["x"].as_array().unwrap().len(), h["depth"].as_array().unwrap().len());
        assert_eq!(h["cp"].as_array().unwrap().len(), nx * nz);
        assert_eq!(v["surface"]["nx"], 60);
        let f = &v["forces"];
        assert!(f["fz"].as_f64().unwrap() < 0.0 && f["sinkage"].as_f64().unwrap() > 0.0, "{f}");
        assert!(f["rw"].as_f64().unwrap() > 0.0);
    }

    #[test]
    fn an_stl_needs_its_units() {
        let e = loft("x.stl", vec![0; 200], &LoftRequest::default()).unwrap_err();
        assert!(e.contains("units"), "{e}");
    }
}

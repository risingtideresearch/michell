//! `michell-web` — a browser front end for the `michell` tools.
//!
//! The first step is the loft viewer: upload a hull file in any input format
//! the CLI reads (native `.hull`, offsets table, `*.grid.json`, IGES, STL),
//! and see what it lofted to. The server does the work with the CLI's own
//! loaders ([`michell_cli::load_hulls_from_bytes`]), so a hull looks here
//! exactly as `michell info` would describe it, and ships the browser a
//! sampled surface, the control net, and the sample grid with per-sample
//! residuals so the loft's fidelity can be judged by eye.

use michell::fit::FitOptions;
use michell::{Hull, SampleGrid};
use michell_cli::{describe_source, load_hulls_from_bytes, LoadSettings, LoadedHull, Source};
use serde_json::{json, Value};

mod variants;

/// Largest upload accepted [bytes]. Big enough for a finely tessellated STL.
pub const MAX_UPLOAD: usize = 128 << 20;

/// Loft options the page can set, mirroring the CLI flags of the same names.
#[derive(Default)]
pub struct LoftRequest {
    /// Design waterline height in the file's frame [m] (`--waterline`).
    pub waterline: Option<f64>,
    /// Centerplane override [m] (`--centerplane`).
    pub centerplane: Option<f64>,
    /// Unit scale for STL (`--units`), e.g. "mm".
    pub units: Option<String>,
    /// Control-net size (`--fit-control NX,NZ`).
    pub fit_control: Option<(usize, usize)>,
    /// Fairing weight (`--fit-fairing`).
    pub fairing: Option<f64>,
    /// Control net of the experimental keel-following loft.
    pub keel_net: Option<(usize, usize)>,
    /// Keel segments per knot span for the experimental lofts.
    pub keel_subdiv: Option<usize>,
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
        for (k, v) in pairs {
            if v.trim().is_empty() {
                continue;
            }
            match k.as_str() {
                "waterline" => r.waterline = Some(num(k, v)?),
                "centerplane" => r.centerplane = Some(num(k, v)?),
                "units" => r.units = Some(v.clone()),
                "fairing" => r.fairing = Some(num(k, v)?),
                "ksub" => {
                    r.keel_subdiv = Some(
                        v.trim()
                            .parse()
                            .map_err(|_| format!("{k}: expected a count, got {v:?}"))?,
                    )
                }
                "knx" | "kns" => {
                    let n = v
                        .trim()
                        .parse::<usize>()
                        .map_err(|_| format!("{k}: expected a count, got {v:?}"))?;
                    let (nx, ns) = r.keel_net.get_or_insert((20, 12));
                    if k == "knx" {
                        *nx = n;
                    } else {
                        *ns = n;
                    }
                }
                "nx" | "nz" => {
                    let n = v
                        .trim()
                        .parse::<usize>()
                        .map_err(|_| format!("{k}: expected a count, got {v:?}"))?;
                    let (nx, nz) = r.fit_control.get_or_insert((0, 0));
                    if k == "nx" {
                        *nx = n;
                    } else {
                        *nz = n;
                    }
                }
                _ => {}
            }
        }
        if let Some((nx, nz)) = r.fit_control {
            if nx == 0 || nz == 0 {
                return Err("control net: give both nx and nz".into());
            }
        }
        Ok(r)
    }

    fn settings(&self) -> Result<LoadSettings, String> {
        let mut s = LoadSettings::default();
        if let Some(w) = self.waterline {
            s.waterline_z = w;
        }
        s.centerplane = self.centerplane;
        if let Some(u) = &self.units {
            s.units = Some(michell_cli::parse_units(u)?);
        }
        if let Some((nx, nz)) = self.fit_control {
            s.fit = FitOptions {
                n_ctrl_x: nx,
                n_ctrl_z: nz,
                ..s.fit
            };
            s.fit_explicit = true;
        }
        if let Some(f) = self.fairing {
            s.fit.fairing = f;
        }
        Ok(s)
    }
}

/// Loft an uploaded file and describe the result as the JSON the page draws.
pub fn loft(name: &str, bytes: Vec<u8>, req: &LoftRequest) -> Result<Value, String> {
    let t0 = std::time::Instant::now();
    // CAD files are also imported by sections, straight from the patches.
    let sectional = is_iges(name, &bytes).then(|| {
        let text = String::from_utf8_lossy(&bytes).into_owned();
        let opts = michell::iges::SectionalOptions {
            waterline_z: req.waterline.unwrap_or(0.0),
            centerplane: req.centerplane,
            ..Default::default()
        };
        let t = std::time::Instant::now();
        (
            michell::iges::import_sectional(&text, &opts),
            t.elapsed().as_secs_f64(),
        )
    });
    let hulls = load_hulls_from_bytes(name, bytes, &req.settings()?)?;
    let elapsed = t0.elapsed().as_secs_f64();
    let vopts = variants::VariantOptions {
        fairing: req.fairing.unwrap_or(0.0),
        keel_net: req.keel_net.unwrap_or((20, 12)),
        keel_subdiv: req.keel_subdiv.unwrap_or(1),
    };
    let mut out: Vec<Value> = hulls.iter().map(|h| hull_json(h, &vopts)).collect();
    if let Some((result, secs)) = sectional {
        match result {
            Ok(fleet) => {
                // Pair each sectional hull with the lofted one at the same
                // centreplane.
                for sh in &fleet.hulls {
                    let nearest = hulls
                        .iter()
                        .enumerate()
                        .min_by(|(_, a), (_, b)| {
                            (a.placement.y - sh.placement.y)
                                .abs()
                                .total_cmp(&(b.placement.y - sh.placement.y).abs())
                        })
                        .map(|(i, _)| i);
                    if let Some(v) =
                        nearest.and_then(|i| out[i]["compare"]["variants"].as_array_mut())
                    {
                        v.push(variants::sectional_variant(sh, secs));
                    }
                }
            }
            Err(e) => eprintln!("sectional import of {name}: {e}"),
        }
    }
    Ok(json!({
        "name": name,
        "seconds": elapsed,
        "hulls": out,
    }))
}

/// IGES by extension or by the section letter in column 73.
fn is_iges(name: &str, bytes: &[u8]) -> bool {
    let lower = name.to_ascii_lowercase();
    if lower.ends_with(".igs") || lower.ends_with(".iges") {
        return true;
    }
    let first = bytes.split(|&b| b == b'\n').next().unwrap_or(&[]);
    first.len() >= 73 && matches!(first[72], b'S' | b'G')
}

fn hull_json(l: &LoadedHull, vopts: &variants::VariantOptions) -> Value {
    let h = &l.hull;
    let fit = match &l.source {
        Source::Native => None,
        Source::Body(r) | Source::Offsets(r) | Source::Grid(r) => Some(r),
        Source::Iges(r) | Source::Stl(r) => Some(&r.fit),
    };
    json!({
        "placement": { "x": l.placement.x, "y": l.placement.y },
        "source": source_kind(&l.source),
        "diagnostics": describe_source(&l.source),
        "under_resolved": fit.is_some_and(|r| r.under_resolved()),
        "relative_rms": fit.map(|r| r.relative_rms()),
        "length": h.length(),
        "beam": michell_cli::max_beam(h),
        "draft": h.draft(),
        "wetted_surface": h.wetted_surface(),
        "displaced_volume": h.displaced_volume(),
        "lcb_x": h.lcb_x(),
        "waterplane_area": h.waterplane_area(),
        "transom": h.transom().map(|t| json!({
            "x": t.x, "depth": t.depth, "half_beam": t.half_beam, "area": t.area,
        })),
        "surface": surface_json(h),
        "control": control_json(h),
        "grid": l.grid.as_ref().map(|g| grid_json(g, h)),
        // Experimental alternative lofts of the same samples, for comparison.
        "compare": l.grid.as_ref().map(|g| variants::variants(g, h, vopts)),
    })
}

fn source_kind(s: &Source) -> &'static str {
    match s {
        Source::Native => "control net",
        Source::Body(_) => "body",
        Source::Offsets(_) => "offsets",
        Source::Grid(_) => "sample grid",
        Source::Iges(_) => "IGES",
        Source::Stl(_) => "STL",
    }
}

/// The half-beam surface sampled for display: stations at every knot and
/// uniformly between (so chines and knuckles land on a mesh line), uniform
/// waterlines likewise augmented with the z knots.
fn surface_json(h: &Hull) -> Value {
    let s = h.surface();
    let xs = display_axis(s.x_domain(), s.knots_x(), 240);
    let zs = display_axis(s.z_domain(), s.knots_z(), 48);
    let mut y = Vec::with_capacity(xs.len() * zs.len());
    for &x in &xs {
        for &z in &zs {
            y.push(s.eval(x, z));
        }
    }
    json!({ "x": xs, "z": zs, "y": y })
}

fn display_axis((a, b): (f64, f64), knots: &[f64], n: usize) -> Vec<f64> {
    let mut v: Vec<f64> = (0..=n)
        .map(|i| a + (b - a) * i as f64 / n as f64)
        .chain(knots.iter().copied())
        .filter(|t| (a..=b).contains(t))
        .collect();
    v.sort_by(f64::total_cmp);
    let tol = 1e-9 * (b - a).abs().max(1e-12);
    v.dedup_by(|p, q| (*p - *q).abs() < tol);
    v
}

/// Control points at their Greville abscissae: where each control value
/// "sits" on the surface, which is how a control net is conventionally drawn.
fn control_json(h: &Hull) -> Value {
    let s = h.surface();
    let greville = |knots: &[f64], p: usize, n: usize| -> Vec<f64> {
        (0..n)
            .map(|i| knots[i + 1..=i + p].iter().sum::<f64>() / p.max(1) as f64)
            .collect()
    };
    let (nx, nz) = (s.n_ctrl_x(), s.n_ctrl_z());
    json!({
        "degree": [s.degree_x(), s.degree_z()],
        "x": greville(s.knots_x(), s.degree_x(), nx),
        "z": greville(s.knots_z(), s.degree_z(), nz),
        "y": s.control(),
    })
}

/// The samples the loft was fitted to, with the loft's residual at each
/// (`null` where the sample was excluded, e.g. a failed CAD inversion).
fn grid_json(g: &SampleGrid, h: &Hull) -> Value {
    let s = h.surface();
    let (st, wl, hb) = (g.stations(), g.waterlines(), g.half_beams());
    let w = g.weights();
    let mut residual = Vec::with_capacity(hb.len());
    for (i, &x) in st.iter().enumerate() {
        for (j, &z) in wl.iter().enumerate() {
            let k = g.idx(i, j);
            let used = w.is_none_or(|w| w[k] > 0.0) && hb[k].is_finite();
            residual.push(if used {
                Some(s.eval(x, z) - hb[k])
            } else {
                None
            });
        }
    }
    json!({
        "x": st,
        "z": wl,
        "y": hb.iter().map(|v| v.is_finite().then_some(*v)).collect::<Vec<_>>(),
        "residual": residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const OFFSETS: &str = "michell-offsets v1
waterlines 0 0.125 0.25 0.375 0.5 0.625
station -5.0   0    0    0    0    0    0
station -4.0   0.36 0.35 0.32 0.27 0.19 0
station -3.0   0.64 0.62 0.57 0.48 0.34 0
station -2.0   0.84 0.81 0.75 0.63 0.44 0
station -1.0   0.96 0.93 0.86 0.72 0.51 0
station  0.0   1.00 0.97 0.89 0.75 0.53 0
station  1.0   0.96 0.93 0.86 0.72 0.51 0
station  2.0   0.84 0.81 0.75 0.63 0.44 0
station  3.0   0.64 0.62 0.57 0.48 0.34 0
station  4.0   0.36 0.35 0.32 0.27 0.19 0
station  5.0   0    0    0    0    0    0
";

    #[test]
    fn lofts_an_offsets_table_with_its_samples() {
        let v = loft("t.offsets", OFFSETS.into(), &LoftRequest::default()).unwrap();
        let hulls = v["hulls"].as_array().unwrap();
        assert_eq!(hulls.len(), 1);
        let h = &hulls[0];
        assert_eq!(h["source"], "offsets");
        assert!((h["length"].as_f64().unwrap() - 10.0).abs() < 1e-9);
        let surf = &h["surface"];
        let (nx, nz) = (
            surf["x"].as_array().unwrap().len(),
            surf["z"].as_array().unwrap().len(),
        );
        assert_eq!(surf["y"].as_array().unwrap().len(), nx * nz);
        let grid = &h["grid"];
        assert_eq!(grid["residual"].as_array().unwrap().len(), 11 * 6);
        let c = &h["control"];
        assert_eq!(
            c["y"].as_array().unwrap().len(),
            c["x"].as_array().unwrap().len() * c["z"].as_array().unwrap().len()
        );
    }

    #[test]
    fn rejects_an_unrecognised_file() {
        assert!(loft("x.txt", b"hello".to_vec(), &LoftRequest::default()).is_err());
    }

    #[test]
    fn stl_without_units_says_so() {
        let e = loft("x.stl", vec![0; 200], &LoftRequest::default()).unwrap_err();
        assert!(e.contains("units"), "{e}");
    }
}

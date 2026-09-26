//! `michell report`: run a dynamic-squat speed sweep from a JSON manifest and
//! render a self-contained PDF — an index page plus, per sweep row, a
//! plan-view wake image and a profile view of sinkage/trim (with a wake
//! elevation cut along the fleet centreline). One CLI invocation, no
//! external tooling and no per-run agent orchestration: the same manifest
//! schema `michell sweep` reads, via [`crate::manifest::parse_manifest`].

use crate::json::Json;
use crate::manifest::{parse_manifest, point_state, Axis, MHull, PointState, KNOT};
use crate::pdf::{Document, Page};
use crate::png;
use crate::manifest::source_hulls;
use michell::float::{fleet_cg, solve_equilibrium_sectional_dynamic, LoadCase};
use michell::iges::{HullPose, Platform};
use michell::sectional::{dynamic_load_closure, multihull_resistance, SectionalHull};
use michell::source::SourceHull;
use michell::{Conditions, FreeWaveSpectrum, Placement, TransomClosure};

const PAGE_W: f64 = 792.0;
const PAGE_H: f64 = 612.0;
const MARGIN: f64 = 36.0;

const BLACK: [f64; 3] = [0.0, 0.0, 0.0];
const GRAY: [f64; 3] = [0.45, 0.45, 0.45];
const BLUE: [f64; 3] = [0.10, 0.35, 0.75];
const ORANGE: [f64; 3] = [0.85, 0.45, 0.05];
/// Hull silhouette fills in the profile view: the part above the actual
/// water surface, and the submerged part the wave integral actually sees.
const HULL_DRY: [f64; 3] = [0.84, 0.83, 0.81];
const HULL_WET: [f64; 3] = [0.60, 0.65, 0.72];
/// `tan` of the Kelvin wedge half-angle, `asin(1/3)` = 19.4712 degrees.
const KELVIN_TAN: f64 = 0.353_553_390_593_273_76;
/// The profile's keel line is the contour where a station's half-beam falls
/// to this fraction of that station's own design-waterline half-beam,
/// floored at `KEEL_MIN_BEAM`. It is a contour and not literally the keel
/// because this kind of loft has no sharp keel edge to find: under the hull
/// the fitted half-beam decays into a few millimetres of ringing that
/// wanders on down to the bottom of the band.
/// Overall propulsive coefficient assumed when turning effective power into
/// brake power on the index. It is one stand-in for the whole chain - hull
/// efficiency, relative rotative efficiency, open-water propeller efficiency
/// and shaft losses - and NOT a propeller calculation. 0.55 is a reasonable
/// working figure for a small craft; change it here if the installation is
/// known. The index states the value next to the column so nobody reads BP
/// as a prediction.
const PROPULSIVE_EFFICIENCY: f64 = 0.55;

pub fn run(manifest_path: &str, out_path: &str, cache_path: Option<&str>) -> Result<(), String> {
    // The report renders its own progress to stderr, which is what the CLI
    // would have done with these lines anyway.
    let mut say = |line: &str, _frac: Option<(usize, usize)>| eprintln!("{line}");
    let pm = parse_manifest(manifest_path, &mut say)?;
    if !pm.dynamic_mode {
        return Err(
            "michell report currently supports dynamic-mode manifests only \
             (options.dynamic: true — a weight axis with per-speed \
             sinkage/trim); for a non-dynamic sweep use `michell \
             sweep` for the CSV/JSON table instead"
                .into(),
        );
    }

    // The Newton solve is what makes this command slow (minutes per row);
    // everything downstream of a (sinkage, trim) pair — re-lofting the wetted
    // hulls, the resistance breakdown, the images — is cheap. So the cache
    // holds only the solved scalars, keyed positionally (row order is
    // deterministic for a given manifest), and is validated against each
    // row's own axis values and speed before being trusted.
    let cached_rows: Vec<CachedRow> = match cache_path {
        Some(p) if std::path::Path::new(p).exists() => load_cache(p)?,
        _ => Vec::new(),
    };
    let mut cache_out: Vec<CachedRow> = Vec::new();

    let midships: Vec<f64> = pm.hulls.iter().map(|h| h.midship).collect();
    let total_rows = pm.points * pm.speeds.len();
    let mut doc = Document::new();
    let mut detail_pages: Vec<Page> = Vec::new();
    let mut rows: Vec<RowSummary> = Vec::new();

    let mut idx = vec![0usize; pm.axes.len()];
    for point in 0..pm.points {
        let vals: Vec<f64> = pm
            .axes
            .iter()
            .zip(&idx)
            .map(|(a, &i)| a.values[i])
            .collect();
        let PointState { poses, loads, .. } = point_state(&pm.hulls, &pm.axes, &vals);
        // Mass and LCG are derived from the hull and point loads carried
        // through their poses, the same way the sweep derives them, so a
        // report and a sweep of the same manifest agree on what floats.
        let cg = fleet_cg(&midships, &loads, &poses);
        let posed = source_hulls(&pm.hulls, &poses);
        if cg.mass <= 0.0 {
            return Err(format!(
                "point {}: fleet carries no mass (all hull and point masses are zero)",
                point + 1
            ));
        }
        let mass = cg.mass;
        let lcg = Some(cg.lcg);
        let pivot_x = cg.lcg;
        let point_label = axis_label(&pm.axes, &vals);

        let mut warm: Option<(f64, f64)> = None;
        for &u in &pm.speeds {
            let row_started = std::time::Instant::now();
            let froude = u / (pm.gravity * pm.l_ref).sqrt();
            let prefix = if point_label.is_empty() {
                String::new()
            } else {
                format!("{point_label}, ")
            };
            let row_no = rows.len() + 1;
            let cond = pm.fluid.make_cond(u)?;

            let cached = cached_rows
                .get(row_no - 1)
                .filter(|c| c.point == point + 1 && (c.speed - u).abs() <= 1e-6 * u.max(1.0))
                .filter(|c| axes_match(&c.axis, &pm.axes, &vals));

            let (sinkage, trim, rt, pe) = if let Some(c) = cached {
                eprintln!(
                    "row {row_no}/{total_rows}: {prefix}U = {u:.3} m/s ({:.2} kn, \
                     Fn {froude:.3}) from cache (sinkage {:.4} m, trim {:.3} deg)",
                    u / KNOT,
                    c.sinkage,
                    c.trim_rad.to_degrees()
                );
                (c.sinkage, c.trim_rad, c.rt, c.pe)
            } else {
                eprintln!(
                    "row {row_no}/{total_rows}: solving {prefix}U = {u:.3} m/s \
                     ({:.2} kn, Fn {froude:.3})...",
                    u / KNOT
                );
                let closure = dynamic_load_closure(&cond, pivot_x, &pm.squat_opts);
                let dyn_eq = solve_equilibrium_sectional_dynamic(
                    &posed,
                    &LoadCase { mass, lcg },
                    pm.density,
                    pm.gravity,
                    &pm.sopts,
                    closure,
                    warm,
                )
                .map_err(|e| format!("point {} U={u}: {e}", point + 1))?;
                eprintln!(
                    "row {row_no}/{total_rows} done in {:.1}s: sinkage {:.4} m, trim {:.3} deg",
                    row_started.elapsed().as_secs_f64(),
                    dyn_eq.sinkage,
                    dyn_eq.trim.to_degrees()
                );
                (dyn_eq.sinkage, dyn_eq.trim, None, None)
            };
            warm = Some((sinkage, trim));

            let platform = Platform {
                sinkage,
                trim,
                pivot_x,
            };
            let members_owned = situate_at(&posed, &platform, &pm.sopts)?;
            let members: Vec<(&SectionalHull, Placement)> =
                members_owned.iter().map(|(h, p)| (h, *p)).collect();
            let resistance = if rt.is_some() {
                None
            } else if members.is_empty() {
                None
            } else {
                Some(
                    multihull_resistance(&members, &cond, &pm.wave_opts, &pm.viscous)
                        .map_err(|e| format!("point {} U={u}: {e}", point + 1))?,
                )
            };
            let rt = rt.or_else(|| resistance.as_ref().map(|r| r.total));
            let pe = pe.or_else(|| resistance.as_ref().map(|r| r.effective_power));

            let page = build_detail_page(
                &mut doc,
                row_no,
                &point_label,
                u,
                froude,
                &cond,
                &members,
                pm.l_ref,
                &pm.hulls,
                &posed,
                &platform,
            )?;
            rows.push(RowSummary {
                row_no,
                point_label: point_label.clone(),
                speed: u,
                froude,
                sinkage,
                trim_deg: trim.to_degrees(),
                rt,
                pe,
            });
            detail_pages.push(page);
            cache_out.push(CachedRow {
                point: point + 1,
                axis: pm
                    .axes
                    .iter()
                    .zip(&vals)
                    .map(|(a, &v)| (a.label.clone(), v))
                    .collect(),
                speed: u,
                froude,
                sinkage,
                trim_rad: trim,
                rt,
                pe,
            });
            // Persist after every row, not once at the end. These solves run
            // for minutes each; a sweep that is interrupted, killed, or fails
            // on a later row used to throw away every row it had already
            // paid for, which defeats the point of having a cache.
            if let Some(p) = cache_path {
                write_cache(p, &cache_out)?;
            }
        }

        for (i, a) in pm.axes.iter().enumerate().rev() {
            idx[i] += 1;
            if idx[i] < a.values.len() {
                break;
            }
            idx[i] = 0;
        }
    }

    if let Some(p) = cache_path {
        write_cache(p, &cache_out)?;
        eprintln!("wrote {p} ({} row(s) cached)", cache_out.len());
    }

    let index = build_index_page(manifest_path, &rows);
    doc.push_page(index);
    for p in detail_pages {
        doc.push_page(p);
    }
    let bytes = doc.finish();
    std::fs::write(out_path, &bytes).map_err(|e| format!("cannot write {out_path}: {e}"))?;
    eprintln!(
        "wrote {out_path}: {} page(s) ({} sweep row(s) + index)",
        rows.len() + 1,
        rows.len()
    );
    Ok(())
}

struct RowSummary {
    row_no: usize,
    point_label: String,
    speed: f64,
    froude: f64,
    sinkage: f64,
    trim_deg: f64,
    rt: Option<f64>,
    pe: Option<f64>,
}

/// A compact repr of a point's non-speed axis values (e.g. a spacing or lcg
/// sweep alongside speed), or an empty string when speed is the only axis —
/// the common case this command was built for.
fn axis_label(axes: &[Axis], vals: &[f64]) -> String {
    axes.iter()
        .zip(vals)
        .map(|(a, v)| format!("{}={v:.3}", a.label))
        .collect::<Vec<_>>()
        .join(", ")
}

/// One row's solved state, as read from or written to `--cache`. Deliberately
/// minimal: just the two Newton unknowns plus enough of the row's own inputs
/// (point index, axis values, speed) to sanity-check that a cached entry
/// still belongs to the row it's about to replace solving for. `rt`/`pe` are
/// cached too so a cache-only run never needs to touch the resistance
/// integral either, but their absence just means "recompute them" — they
/// are not load-bearing for validation.
struct CachedRow {
    point: usize,
    axis: Vec<(String, f64)>,
    speed: f64,
    froude: f64,
    sinkage: f64,
    trim_rad: f64,
    rt: Option<f64>,
    pe: Option<f64>,
}

fn axes_match(cached: &[(String, f64)], axes: &[Axis], vals: &[f64]) -> bool {
    cached.len() == axes.len()
        && axes.iter().zip(vals).all(|(a, &v)| {
            cached
                .iter()
                .any(|(k, cv)| k == &a.label && (cv - v).abs() <= 1e-9 * v.abs().max(1.0))
        })
}

fn load_cache(path: &str) -> Result<Vec<CachedRow>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let doc = crate::json::parse(&text).map_err(|e| format!("{path}: {e}"))?;
    let rows = doc
        .get("rows")
        .and_then(Json::as_arr)
        .ok_or_else(|| format!("{path}: expected an object with a \"rows\" array"))?;
    let mut out = Vec::with_capacity(rows.len());
    for r in rows {
        let field = |name: &str| {
            r.get(name)
                .and_then(Json::as_f64)
                .ok_or_else(|| format!("{path}: row missing numeric \"{name}\""))
        };
        let mut axis = Vec::new();
        if let Some(Json::Obj(kv)) = r.get("axis") {
            for (k, v) in kv {
                if let Some(f) = v.as_f64() {
                    axis.push((k.clone(), f));
                }
            }
        }
        out.push(CachedRow {
            point: field("point")? as usize,
            axis,
            speed: field("speed")?,
            froude: field("froude").unwrap_or(0.0),
            sinkage: field("sinkage")?,
            trim_rad: field("trim_rad")?,
            rt: field("rt").ok(),
            pe: field("pe").ok(),
        });
    }
    Ok(out)
}

fn write_cache(path: &str, rows: &[CachedRow]) -> Result<(), String> {
    let mut out = String::from("{\n  \"rows\": [\n");
    for (i, r) in rows.iter().enumerate() {
        if i > 0 {
            out.push_str(",\n");
        }
        out.push_str(&format!("    {{\"point\":{},\"axis\":{{", r.point));
        for (j, (k, v)) in r.axis.iter().enumerate() {
            if j > 0 {
                out.push(',');
            }
            out.push_str(&format!("{k:?}:{v}"));
        }
        out.push_str(&format!(
            "}},\"speed\":{},\"froude\":{},\"sinkage\":{},\"trim_rad\":{}",
            r.speed, r.froude, r.sinkage, r.trim_rad
        ));
        if let Some(rt) = r.rt {
            out.push_str(&format!(",\"rt\":{rt}"));
        }
        if let Some(pe) = r.pe {
            out.push_str(&format!(",\"pe\":{pe}"));
        }
        out.push('}');
    }
    out.push_str("\n  ]\n}\n");
    std::fs::write(path, out).map_err(|e| format!("cannot write {path}: {e}"))
}

/// Cut the fleet directly at a known (sinkage, trim) — the cheap,
/// non-iterative half of what the Newton solver does on every call. Used
/// both for a fresh solve's result and to rebuild a cached row's geometry
/// without re-solving it.
fn situate_at(
    hulls: &[SourceHull],
    platform: &Platform,
    opts: &michell::iges::SectionalOptions,
) -> Result<Vec<(SectionalHull, Placement)>, String> {
    let mut members = Vec::new();
    for h in hulls {
        let o = michell::iges::SectionalOptions {
            waterline_z: h.waterline_z,
            ..*opts
        };
        if let Some(sh) = h
            .source
            .situate_sectional(h.index, h.waterline_z, &h.pose, platform, &o)
            .map_err(|e| format!("{e}"))?
        {
            members.push((sh.hull, sh.placement));
        }
    }
    Ok(members)
}

#[allow(clippy::too_many_arguments)]
fn build_detail_page(
    doc: &mut Document,
    row_no: usize,
    point_label: &str,
    speed: f64,
    froude: f64,
    cond: &Conditions,
    members: &[(&SectionalHull, Placement)],
    l_ref: f64,
    hulls: &[MHull],
    posed: &[SourceHull],
    platform: &Platform,
) -> Result<Page, String> {
    let mut page = Page::new(PAGE_W, PAGE_H);
    let knots = speed / KNOT;
    let title = if point_label.is_empty() {
        format!("Row {row_no}: U = {speed:.3} m/s ({knots:.2} kn, Fn {froude:.3})")
    } else {
        format!("Row {row_no}: {point_label}, U = {speed:.3} m/s ({knots:.2} kn, Fn {froude:.3})")
    };
    page.text(MARGIN, PAGE_H - 26.0, 14.0, BLACK, &title);
    page.text_right(PAGE_W - MARGIN, PAGE_H - 26.0, 9.0, GRAY, "index");
    page.link_to_page(
        [
            PAGE_W - MARGIN - 40.0,
            PAGE_H - 34.0,
            PAGE_W - MARGIN,
            PAGE_H - 20.0,
        ],
        0,
    );

    // The plan panel takes the taller share. Its wake is drawn to true
    // proportions now (see draw_plan_view), so panel height buys wake
    // coverage astern rather than just pixels.
    let plan_rect = [MARGIN, 240.0, PAGE_W - 2.0 * MARGIN, 322.0];
    let profile_rect = [MARGIN, 66.0, PAGE_W - 2.0 * MARGIN, 140.0];

    let plan_caption = draw_plan_view(doc, &mut page, plan_rect, members, cond, l_ref)?;
    draw_caption(&mut page, plan_rect[0], plan_rect[1] - 12.0, &plan_caption);
    page.text(
        plan_rect[0],
        plan_rect[1] + plan_rect[3] + 4.0,
        9.0,
        BLACK,
        "Plan view over wave field (ship advances toward +x)",
    );

    let profile_hull_y = hulls
        .iter()
        .zip(posed)
        .map(|(h, p)| h.centerplane + p.pose.dy)
        .find(|y| y.abs() > 1e-9)
        .unwrap_or(0.0);
    let profile_caption = draw_profile_view(
        &mut page,
        profile_rect,
        posed,
        platform,
        members,
        cond,
        profile_hull_y,
    )?;
    draw_caption(
        &mut page,
        profile_rect[0],
        profile_rect[1] - 12.0,
        &profile_caption,
    );
    page.text(
        profile_rect[0],
        profile_rect[1] + profile_rect[3] + 4.0,
        9.0,
        BLACK,
        "Profile: hull at the solved sinkage and trim, with a wave elevation cut",
    );

    Ok(page)
}

/// Draw a panel caption, one line per `\n`, running downward from `top`.
fn draw_caption(page: &mut Page, x: f64, top: f64, text: &str) {
    for (i, line) in text.split('\n').enumerate() {
        page.text(x, top - 9.5 * i as f64, 8.0, GRAY, line);
    }
}

/// Render the plan-view Kelvin wake heatmap for this row directly into the
/// page as an embedded raster image (the same field `michell wake` plots,
/// simplified to fixed sizing/region since this is a canned report panel).
fn draw_plan_view(
    doc: &mut Document,
    page: &mut Page,
    rect: [f64; 4],
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    l_ref: f64,
) -> Result<String, String> {
    if members.is_empty() {
        page.rect_stroke(rect, GRAY, 0.75);
        page.text(
            rect[0] + 8.0,
            rect[1] + rect[3] / 2.0,
            10.0,
            GRAY,
            "fleet is dry",
        );
        return Ok("no wetted hulls at this row".to_string());
    }

    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_abs = 0.0f64;
    for (h, p) in members {
        let (h0, h1) = h.x_range();
        x_lo = x_lo.min(h0 + p.x);
        x_hi = x_hi.max(h1 + p.x);
        y_abs = y_abs.max(p.y.abs());
    }
    // Frame the wake to TRUE PROPORTIONS. The panel is far wider than it is
    // tall, so asking for a fixed astern reach and letting the raster stretch
    // to fill the rect squashed the transverse axis about 2.6x here: the
    // Kelvin wedge drew at roughly 8 degrees instead of its 19.47, and the
    // hulls came out as needles. Derive the astern reach from the panel's own
    // aspect instead, so the wedge exactly fills the panel height:
    //   isotropic  =>  yh = a * Lx,  a = rect_h / (2 rect_w)
    //   wedge      =>  yh = KELVIN_TAN * d,  Lx = d + fore_aft
    //   hence      =>  d = a * fore_aft / (KELVIN_TAN - a)
    let ahead = 0.35 * l_ref;
    let fore_aft = (x_hi - x_lo) + ahead;
    let a = 0.5 * rect[3] / rect[2];
    let mut lx = if a < KELVIN_TAN {
        fore_aft * KELVIN_TAN / (KELVIN_TAN - a)
    } else {
        // Panel tall enough that the wedge never leaves its sides; fall back
        // to a fixed reach rather than dividing by something non-positive.
        fore_aft + 3.0 * l_ref
    };
    let mut yh = a * lx;
    // Never clip the fleet itself, even at the cost of some wake.
    let need = 1.15 * y_abs + 0.12 * l_ref;
    if yh < need {
        yh = need;
        lx = yh / a;
    }
    let x1 = x_hi + ahead;
    let x0 = x1 - lx;
    let (y0, y1) = (-yh, yh);

    let nx = 760usize;
    let ny = ((nx as f64) * (y1 - y0) / (x1 - x0))
        .round()
        .clamp(64.0, 900.0) as usize;

    // Plain thin-ship field, not the resistance/squat integrals' transom
    // virtual-appendage closure: that closure fixes up an integrated force
    // and is not a model of the actual (breaking, unsteady) near-transom
    // sea surface, so drawing it here would show fabricated structure
    // behind a wet transom that doesn't correspond to anything real.
    let mut spec = FreeWaveSpectrum::new_sectional(members, cond, TransomClosure::None)
        .map_err(|e| format!("{e}"))?;
    let grid = spec
        .elevation_grid(x0, x1, y0, y1, nx, ny)
        .map_err(|e| format!("{e}"))?;

    let mut vmax = 0.0f64;
    {
        let mut abs: Vec<f64> = grid.zeta.iter().map(|v| v.abs()).collect();
        abs.sort_by(|a, b| a.total_cmp(b));
        if !abs.is_empty() {
            vmax = abs[((abs.len() - 1) as f64 * 0.995) as usize].max(1e-12);
        }
    }

    // Fade the field from the stern forward: the free-wave spectrum is a
    // downstream representation, so over and ahead of the ship it is not the
    // real surface. This used to be a hard switch at x_lo, which painted a
    // full-height vertical seam straight down the panel at exactly the
    // transom - a fake transom wave, drawn by the renderer rather than the
    // physics. Ramp it smoothly across the stern instead.
    const FADE_TOWARD: [u8; 3] = [0xf0, 0xef, 0xec];
    const FADE_FRACTION: f64 = 0.55;
    let ramp = 0.45 * l_ref;
    let mut rgb = vec![0u8; 3 * nx * ny];
    for iy in 0..ny {
        let row = ny - 1 - iy; // row 0 = top = +y edge, matching Document::image's convention
        for ix in 0..nx {
            let t = grid.get(ix, iy) / vmax;
            let mut c = png::diverging(t);
            let u = ((grid.x(ix) - (x_lo - ramp)) / (2.0 * ramp)).clamp(0.0, 1.0);
            let smooth = u * u * (3.0 - 2.0 * u); // smoothstep, C1 at both ends
            if smooth > 0.0 {
                c = png::fade(c, FADE_TOWARD, FADE_FRACTION * smooth);
            }
            rgb[3 * (row * nx + ix)..3 * (row * nx + ix) + 3].copy_from_slice(&c);
        }
    }
    const HULL_GRAY: [u8; 3] = [0x52, 0x51, 0x4e];
    for (h, p) in members {
        let (h0, h1) = h.x_range();
        for ix in 0..nx {
            let x = grid.x(ix) - p.x;
            if x < h0 || x > h1 {
                continue;
            }
            let half_beam = h.waterline_half_beam(x);
            for iy in 0..ny {
                if (grid.y(iy) - p.y).abs() <= half_beam {
                    let row = ny - 1 - iy;
                    rgb[3 * (row * nx + ix)..3 * (row * nx + ix) + 3].copy_from_slice(&HULL_GRAY);
                }
            }
        }
    }

    let id = doc.image(nx, ny, rgb);
    page.rect_stroke(rect, GRAY, 0.75);
    page.image(id, rect);

    Ok(format!(
        "true proportions, {:.1} hull lengths of wake astern; zeta +-{vmax:.4} m; \
         x {x0:.1}..{x1:.1} m, y {y0:.1}..{y1:.1} m{}",
        (x_lo - x0) / l_ref,
        if grid.resolution_limited {
            " (grid-resolution limited)"
        } else {
            ""
        }
    ))
}

/// A hull's centreplane profile at the solved attitude, in the water frame
/// (x, depth below the actual free surface — down-positive): the keel line,
/// the top of the geometry, the design-waterline mark, and the (top, keel)
/// pairs that bound the hull's silhouette.
struct HullProfile {
    keel: Vec<[f64; 2]>,
    /// The top of the source geometry (sheer, or deck edge).
    top: Vec<[f64; 2]>,
    design_wl: Vec<[f64; 2]>,
    /// `[top, keel]` pairs at the stations that carry any geometry.
    silhouette: Vec<[[f64; 2]; 2]>,
}

/// The profile from the hull's tessellation at its pose and the platform
/// state: per x bin, the highest and lowest vertex. The design waterline is
/// the posed image of the design floatplane; the pose and platform move the
/// hull rigidly (with uniform scale), so that image is the affine map fitted
/// between the unposed and posed vertices.
fn hull_profile_at(h: &SourceHull, platform: &Platform) -> Result<HullProfile, String> {
    const STATIONS: usize = 120;
    let (at_rest, _) = h
        .source
        .posed_tessellation(h.index, h.waterline_z, &HullPose::default(), &Platform::default())
        .map_err(|e| format!("{e}"))?;
    let (verts, _) = h
        .source
        .posed_tessellation(h.index, h.waterline_z, &h.pose, platform)
        .map_err(|e| format!("{e}"))?;
    let (mut x0, mut x1) = (f64::INFINITY, f64::NEG_INFINITY);
    for v in &verts {
        x0 = x0.min(v[0]);
        x1 = x1.max(v[0]);
    }
    let mut hi = vec![f64::NEG_INFINITY; STATIONS];
    let mut lo = vec![f64::INFINITY; STATIONS];
    let bin = |x: f64| (((x - x0) / (x1 - x0) * STATIONS as f64) as usize).min(STATIONS - 1);
    for v in &verts {
        let b = bin(v[0]);
        hi[b] = hi[b].max(v[2]);
        lo[b] = lo[b].min(v[2]);
    }
    let mut keel = Vec::new();
    let mut top = Vec::new();
    let mut silhouette = Vec::new();
    for b in 0..STATIONS {
        if hi[b] < lo[b] {
            continue;
        }
        let x = x0 + (x1 - x0) * (b as f64 + 0.5) / STATIONS as f64;
        let (t, k) = ([x, -hi[b]], [x, -lo[b]]);
        top.push(t);
        keel.push(k);
        silhouette.push([t, k]);
    }
    // The design floatplane at rest is z = 0; fit x', z' = A·(x, z, 1) over
    // the vertex correspondence (exact for the affine pose) and map it.
    let fit = affine_xz(&at_rest, &verts);
    let (ax0, ax1) = at_rest
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v[0]), b.max(v[0])));
    let design_wl = (0..=40)
        .map(|i| {
            let x = ax0 + (ax1 - ax0) * i as f64 / 40.0;
            let (xw, zw) = fit(x, 0.0);
            [xw, -zw]
        })
        .collect();
    Ok(HullProfile {
        keel,
        top,
        design_wl,
        silhouette,
    })
}

/// The affine map in the (x, z) plane taking `from` to `to` (paired
/// vertices), by least squares.
fn affine_xz(from: &[[f64; 3]], to: &[[f64; 3]]) -> impl Fn(f64, f64) -> (f64, f64) {
    // Normal equations for [x z 1] against each target coordinate.
    let mut m = [[0.0f64; 3]; 3];
    let mut rx = [0.0f64; 3];
    let mut rz = [0.0f64; 3];
    for (a, b) in from.iter().zip(to) {
        let r = [a[0], a[2], 1.0];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += r[i] * r[j];
            }
            rx[i] += r[i] * b[0];
            rz[i] += r[i] * b[2];
        }
    }
    let solve = |rhs: [f64; 3]| -> [f64; 3] {
        let det = |m: &[[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let d = det(&m);
        let mut out = [0.0; 3];
        for (k, o) in out.iter_mut().enumerate() {
            let mut mk = m;
            for i in 0..3 {
                mk[i][k] = rhs[i];
            }
            *o = if d != 0.0 { det(&mk) / d } else { 0.0 };
        }
        out
    };
    let (cx, cz) = (solve(rx), solve(rz));
    move |x, z| (cx[0] * x + cx[1] * z + cx[2], cz[0] * x + cz[1] * z + cz[2])
}

#[allow(clippy::too_many_arguments)]
fn draw_profile_view(
    page: &mut Page,
    rect: [f64; 4],
    posed: &[SourceHull],
    platform: &Platform,
    members: &[(&SectionalHull, Placement)],
    cond: &Conditions,
    wave_cut_y: f64,
) -> Result<String, String> {
    let profiles: Vec<HullProfile> = posed
        .iter()
        .map(|h| hull_profile_at(h, platform))
        .collect::<Result<_, _>>()?;

    // Wave elevation cut along the fleet centreline at the reference hull's
    // transverse offset ("the wake, in profile"); z is up here, so it adds
    // to the flat sea surface as a negative depth.
    let mut wave_cut: Vec<[f64; 2]> = Vec::new();
    if !members.is_empty() {
        let mut spec = FreeWaveSpectrum::new_sectional(members, cond, TransomClosure::None)
            .map_err(|e| format!("{e}"))?;
        let (x0, x1) = data_x_range(&profiles);
        const N: usize = 160;
        for i in 0..N {
            let x = x0 + (x1 - x0) * i as f64 / (N - 1) as f64;
            if let Ok(eta) = spec.elevation_at(x, wave_cut_y) {
                wave_cut.push([x, -eta]);
            }
        }
    }

    // Data bounds (down-positive z; we flip to page y-up at draw time). The
    // top curve spans the hull's full length regardless of local wetness,
    // so it sets the x-range.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut z_lo = 0.0f64; // includes the still-water line
    let mut z_hi = 0.0f64;
    for p in &profiles {
        for pt in p.keel.iter().chain(&p.top).chain(&p.design_wl) {
            x_lo = x_lo.min(pt[0]);
            x_hi = x_hi.max(pt[0]);
            z_lo = z_lo.min(pt[1]);
            z_hi = z_hi.max(pt[1]);
        }
    }
    for pt in &wave_cut {
        z_lo = z_lo.min(pt[1]);
        z_hi = z_hi.max(pt[1]);
    }
    if !(x_hi > x_lo) {
        x_lo = -1.0;
        x_hi = 1.0;
    }
    let x_pad = 0.03 * (x_hi - x_lo);
    x_lo -= x_pad;
    x_hi += x_pad;
    let z_pad = 0.15 * (z_hi - z_lo).max(0.02);
    z_lo -= z_pad;
    z_hi += z_pad;

    // A hull's draft is a small fraction of its length (this one is roughly
    // 30:1), so true-proportion scaling would draw everything of interest
    // into a sliver a few points tall. Scale x and z independently instead,
    // filling most of the panel's height with the z data, and report the
    // resulting exaggeration so the distortion is stated, not hidden.
    let sx = rect[2] / (x_hi - x_lo);
    let sz = 0.88 * rect[3] / (z_hi - z_lo);
    let z_mid = 0.5 * (z_hi + z_lo);
    let y_mid = rect[1] + 0.5 * rect[3];
    let to_page =
        |x: f64, z: f64| -> [f64; 2] { [rect[0] + (x - x_lo) * sx, y_mid - (z - z_mid) * sz] };

    page.rect_stroke(rect, GRAY, 0.75);

    // Hull silhouette first, so every line below reads on top of it. Filling
    // the body's centreplane section in two tones, split at the ACTUAL water
    // surface, is what makes this panel show a boat: the darker area is the
    // hull the wave integral sees, and the gap between the flat water line
    // and the orange design waterline is the sinkage and trim this row
    // solved for. Two stroked curves alone showed neither.
    for p in &profiles {
        let wet: Vec<[[f64; 2]; 2]> = p
            .silhouette
            .iter()
            .filter(|sl| sl[1][1] > 0.0)
            .map(|sl| [[sl[0][0], sl[0][1].max(0.0)], sl[1]])
            .collect();
        let dry: Vec<[[f64; 2]; 2]> = p
            .silhouette
            .iter()
            .filter(|sl| sl[0][1] < 0.0)
            .map(|sl| [sl[0], [sl[1][0], sl[1][1].min(0.0)]])
            .collect();
        for (band, fill) in [(dry, HULL_DRY), (wet, HULL_WET)] {
            if band.len() < 2 {
                continue;
            }
            let mut poly: Vec<[f64; 2]> =
                band.iter().map(|sl| to_page(sl[0][0], sl[0][1])).collect();
            poly.extend(band.iter().rev().map(|sl| to_page(sl[1][0], sl[1][1])));
            page.filled_polygon(&poly, fill);
        }
    }

    // Still water: exactly flat by construction (sinkage/trim are baked
    // into the hull curves below), always drawn across the full panel width.
    page.polyline(&[to_page(x_lo, 0.0), to_page(x_hi, 0.0)], BLUE, 0.5);

    for p in &profiles {
        let top_pts: Vec<[f64; 2]> = p.top.iter().map(|pt| to_page(pt[0], pt[1])).collect();
        page.polyline(&top_pts, GRAY, 0.6);
        let keel_pts: Vec<[f64; 2]> = p.keel.iter().map(|pt| to_page(pt[0], pt[1])).collect();
        page.polyline(&keel_pts, BLACK, 1.2);
        let wl_pts: Vec<[f64; 2]> = p.design_wl.iter().map(|pt| to_page(pt[0], pt[1])).collect();
        page.polyline(&wl_pts, ORANGE, 0.9);
    }

    if wave_cut.len() > 1 {
        let cut_pts: Vec<[f64; 2]> = wave_cut.iter().map(|pt| to_page(pt[0], pt[1])).collect();
        page.polyline(&cut_pts, BLUE, 0.9);
    }

    let caption = format!(
        "shaded: hull profile, darker = submerged; black: keel; gray: top of the \
         geometry; orange: design waterline\n\
         blue: still water and the wake cut at y = {wave_cut_y:.2} m; vertical \
         exaggeration {:.1}x",
        sz / sx
    );
    Ok(caption)
}

fn data_x_range(profiles: &[HullProfile]) -> (f64, f64) {
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    for p in profiles {
        for pt in p.top.iter() {
            x_lo = x_lo.min(pt[0]);
            x_hi = x_hi.max(pt[0]);
        }
    }
    if x_hi > x_lo {
        (x_lo, x_hi)
    } else {
        (-1.0, 1.0)
    }
}

fn build_index_page(manifest_path: &str, rows: &[RowSummary]) -> Page {
    let mut page = Page::new(PAGE_W, PAGE_H);
    page.text(MARGIN, PAGE_H - 30.0, 16.0, BLACK, "Speed sweep report");
    page.text(MARGIN, PAGE_H - 46.0, 9.0, GRAY, manifest_path);

    let cols = [
        "row",
        "point",
        "U [m/s]",
        "U [kn]",
        "Fn",
        "sinkage [m]",
        "trim [deg]",
        "Rt [kN]",
        "Pe [kW]",
        "BP [kW]",
    ];
    let col_x = [
        MARGIN,
        MARGIN + 34.0,
        MARGIN + 210.0,
        MARGIN + 268.0,
        MARGIN + 322.0,
        MARGIN + 372.0,
        MARGIN + 450.0,
        MARGIN + 528.0,
        MARGIN + 592.0,
        MARGIN + 656.0,
    ];
    let mut y = PAGE_H - 74.0;
    for (c, &x) in cols.iter().zip(&col_x) {
        page.text(x, y, 9.0, GRAY, c);
    }
    y -= 4.0;
    page.polyline(&[[MARGIN, y], [PAGE_W - MARGIN, y]], GRAY, 0.5);

    let row_h = 14.0;
    for row in rows {
        y -= row_h;
        if y < MARGIN + 20.0 {
            // A long sweep can overflow one index page; keep going on the
            // same page rather than paginating the index — links still work
            // (each entry points at its own detail page), it just runs off
            // the printable area for very large sweeps.
        }
        let target = row.row_no; // detail pages start at document index 1
        page.text(col_x[0], y, 9.0, BLUE, &format!("{}", row.row_no));
        page.link_to_page([col_x[0] - 2.0, y - 3.0, col_x[1] - 4.0, y + 9.0], target);
        page.text(col_x[1], y, 9.0, BLACK, &row.point_label);
        page.text(col_x[2], y, 9.0, BLACK, &format!("{:.3}", row.speed));
        page.text(col_x[3], y, 9.0, BLACK, &format!("{:.2}", row.speed / KNOT));
        page.text(col_x[4], y, 9.0, BLACK, &format!("{:.3}", row.froude));
        page.text(col_x[5], y, 9.0, BLACK, &format!("{:.4}", row.sinkage));
        page.text(col_x[6], y, 9.0, BLACK, &format!("{:.3}", row.trim_deg));
        if let Some(rt) = row.rt {
            page.text(col_x[7], y, 9.0, BLACK, &format!("{:.3}", rt / 1000.0));
        }
        if let Some(pe) = row.pe {
            page.text(col_x[8], y, 9.0, BLACK, &format!("{:.2}", pe / 1000.0));
            page.text(
                col_x[9],
                y,
                9.0,
                BLACK,
                &format!("{:.2}", pe / 1000.0 / PROPULSIVE_EFFICIENCY),
            );
        }
    }

    // Brake power is the one column here that is not computed from the
    // geometry, so state the assumption under the table rather than leaving
    // a reader to infer it from a column heading.
    y -= 22.0;
    draw_caption(
        &mut page,
        MARGIN,
        y,
        &format!(
            "Pe is effective power, Rt x U. BP = Pe / {PROPULSIVE_EFFICIENCY:.2}, an assumed \
             overall propulsive coefficient.\n\
             That one number stands in for hull, relative-rotative, open-water propeller and \
             shaft efficiency together - a placeholder for a real propeller, not a prediction \
             from one."
        ),
    );
    page
}

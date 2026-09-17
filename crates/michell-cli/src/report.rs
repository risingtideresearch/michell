//! `michell report`: run a dynamic-squat speed sweep from a JSON manifest and
//! render a self-contained PDF — an index page plus, per sweep row, a
//! plan-view wake image and a profile view of sinkage/trim (with a wake
//! elevation cut along the fleet centreline). One CLI invocation, no
//! external tooling and no per-run agent orchestration: the same manifest
//! schema `michell sweep` reads, via [`crate::manifest::parse_manifest`].

use crate::json::Json;
use crate::manifest::{parse_manifest, point_state, Axis, MHull, PointState};
use crate::pdf::{Document, Page};
use crate::png;
use michell::body::{Body, BodyOptions};
use michell::float::{solve_equilibrium_bodies_dynamic, LoadCase};
use michell::iges::{HullPose, Platform};
use michell::squat::dynamic_load_closure;
use michell::{
    multihull_resistance_with, Conditions, FreeWaveSpectrum, Hull, Placement, TransomClosure,
};

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
const KEEL_BEAM_FRACTION: f64 = 0.10;
const KEEL_MIN_BEAM: f64 = 2.0e-3;

pub fn run(manifest_path: &str, out_path: &str, cache_path: Option<&str>) -> Result<(), String> {
    let pm = parse_manifest(manifest_path)?;
    if !pm.dynamic_mode {
        return Err(
            "michell report currently supports dynamic-mode manifests only \
             (options.dynamic: true — a weight axis with per-speed \
             sinkage/trim); for a non-dynamic or heeled sweep use `michell \
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

    let bodies: Vec<&Body> = pm.hulls.iter().map(|h| &h.body).collect();
    let total_rows = pm.points * pm.speeds.len();
    let mut doc = Document::new();
    let mut detail_pages: Vec<Page> = Vec::new();
    let mut rows: Vec<RowSummary> = Vec::new();

    let mut idx = vec![0usize; pm.axes.len()];
    for point in 0..pm.points {
        let vals: Vec<f64> = pm.axes.iter().zip(&idx).map(|(a, &i)| a.values[i]).collect();
        let PointState {
            poses,
            weight,
            lcg,
            ..
        } = point_state(&pm.hulls, &pm.axes, &vals);
        let mass = weight.ok_or("dynamic report requires a weight axis (checked above)")?;
        let pivot_x = lcg.unwrap_or(0.0);
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
                    "row {row_no}/{total_rows}: {prefix}U = {u:.3} m/s (Fn {froude:.3}) \
                     from cache (sinkage {:.4} m, trim {:.3} deg)",
                    c.sinkage,
                    c.trim_rad.to_degrees()
                );
                (c.sinkage, c.trim_rad, c.rt, c.pe)
            } else {
                eprintln!(
                    "row {row_no}/{total_rows}: solving {prefix}U = {u:.3} m/s (Fn {froude:.3})..."
                );
                let closure = dynamic_load_closure(&cond, pivot_x, &pm.squat_opts);
                let dyn_eq = solve_equilibrium_bodies_dynamic(
                    &bodies,
                    0.0,
                    &poses,
                    &LoadCase { mass, lcg },
                    pm.density,
                    pm.gravity,
                    &pm.bopts,
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
            let (members_owned, band_exceeded) =
                situate_at(&bodies, &poses, &platform, &pm.bopts)?;
            if band_exceeded > 0 {
                eprintln!(
                    "row {row_no}/{total_rows}: WARNING {band_exceeded} wetted sample(s) rose \
                     above the top of the lofted band. That geometry is not in the .hull file \
                     and was taken as zero half-beam, so this row understates the immersed \
                     hull. Re-loft with a taller --band."
                );
            }
            let members: Vec<(&Hull, Placement)> =
                members_owned.iter().map(|(h, p)| (h, *p)).collect();
            let resistance = if rt.is_some() {
                None
            } else if members.is_empty() {
                None
            } else {
                Some(
                    multihull_resistance_with(&members, &cond, &pm.wave_opts, pm.form_factor)
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
                &poses,
                &platform,
                band_exceeded,
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
                band_exceeded,
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
    /// Wetted samples that fell above the lofted band at this row's pose.
    /// Non-zero means the hull the forces were computed on is missing its
    /// immersed upper stern, so `rt`/`pe` in the index understate the row.
    band_exceeded: usize,
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

/// Re-loft the fleet directly at a known (sinkage, trim) — the cheap,
/// non-iterative half of what the Newton solver does on every call. Used
/// both for a fresh solve's result and to rebuild a cached row's geometry
/// without re-solving it.
fn situate_at(
    bodies: &[&Body],
    poses: &[HullPose],
    platform: &Platform,
    opts: &BodyOptions,
) -> Result<(Vec<(Hull, Placement)>, usize), String> {
    let mut members = Vec::new();
    let mut band_exceeded = 0usize;
    for (body, pose) in bodies.iter().zip(poses) {
        if let Some(sb) = body
            .situate(0.0, pose, platform, opts)
            .map_err(|e| format!("{e}"))?
        {
            band_exceeded += sb.band_exceeded;
            members.push((sb.hull, sb.placement));
        }
    }
    Ok((members, band_exceeded))
}

#[allow(clippy::too_many_arguments)]
fn build_detail_page(
    doc: &mut Document,
    row_no: usize,
    point_label: &str,
    speed: f64,
    froude: f64,
    cond: &Conditions,
    members: &[(&Hull, Placement)],
    l_ref: f64,
    hulls: &[MHull],
    poses: &[HullPose],
    platform: &Platform,
    band_exceeded: usize,
) -> Result<Page, String> {
    let mut page = Page::new(PAGE_W, PAGE_H);
    let title = if point_label.is_empty() {
        format!("Row {row_no}: U = {speed:.3} m/s (Fn {froude:.3})")
    } else {
        format!("Row {row_no}: {point_label}, U = {speed:.3} m/s (Fn {froude:.3})")
    };
    page.text(MARGIN, PAGE_H - 26.0, 14.0, BLACK, &title);
    page.text_right(PAGE_W - MARGIN, PAGE_H - 26.0, 9.0, GRAY, "index");
    page.link_to_page([PAGE_W - MARGIN - 40.0, PAGE_H - 34.0, PAGE_W - MARGIN, PAGE_H - 20.0], 0);

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
        .zip(poses)
        .map(|(h, p)| h.body.centerplane() + p.dy)
        .find(|y| y.abs() > 1e-9)
        .unwrap_or(0.0);
    let profile_caption = draw_profile_view(
        &mut page,
        profile_rect,
        hulls,
        poses,
        platform,
        members,
        cond,
        profile_hull_y,
        band_exceeded,
    )?;
    draw_caption(&mut page, profile_rect[0], profile_rect[1] - 12.0, &profile_caption);
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
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    l_ref: f64,
) -> Result<String, String> {
    if members.is_empty() {
        page.rect_stroke(rect, GRAY, 0.75);
        page.text(rect[0] + 8.0, rect[1] + rect[3] / 2.0, 10.0, GRAY, "fleet is dry");
        return Ok("no wetted hulls at this row".to_string());
    }

    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_abs = 0.0f64;
    for (h, p) in members {
        let (h0, h1) = h.surface().x_domain();
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
    let ny = ((nx as f64) * (y1 - y0) / (x1 - x0)).round().clamp(64.0, 900.0) as usize;

    // Plain thin-ship field, not the resistance/squat integrals' transom
    // virtual-appendage closure: that closure fixes up an integrated force
    // and is not a model of the actual (breaking, unsteady) near-transom
    // sea surface, so drawing it here would show fabricated structure
    // behind a wet transom that doesn't correspond to anything real.
    let mut spec = FreeWaveSpectrum::new_with_transom(members, cond, TransomClosure::None)
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
        let (h0, h1) = h.surface().x_domain();
        for ix in 0..nx {
            let x = grid.x(ix) - p.x;
            if x < h0 || x > h1 {
                continue;
            }
            let half_beam = h.surface().eval(x, 0.0);
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

/// Positive angle raises the +x side, down-positive z — mirrors
/// `michell::body`'s private `FrameMap::body_to_water`. That struct isn't
/// exported (situating never needs the inverse outside the crate), so the
/// handful of trig lines are duplicated here rather than plumbing a new
/// public API through the library for a one-off report view.
fn rot_down(x: f64, zd: f64, px: f64, pzd: f64, sin: f64, cos: f64) -> (f64, f64) {
    let (dx, dz) = (x - px, zd - pzd);
    (px + dx * cos + dz * sin, pzd + dz * cos - dx * sin)
}

/// Map a point in a hull's own body frame (`xb`, `zb` down from the band top)
/// to the water frame (`xw`, depth below the *actual* free surface — always
/// 0 there by construction, since the platform's sinkage/trim is baked in).
fn body_to_water(
    xb: f64,
    zb: f64,
    body: &Body,
    pose: &HullPose,
    pose_pivot_x: f64,
    platform: &Platform,
) -> (f64, f64) {
    let mut x = xb;
    let mut z = zb - body.waterline();
    if pose.trim != 0.0 {
        let (s, c) = pose.trim.sin_cos();
        (x, z) = rot_down(x, z, pose_pivot_x, 0.0, s, c);
    }
    x += pose.dx;
    z += pose.dz;
    let zw = -platform.sinkage;
    if platform.trim != 0.0 {
        let (s, c) = platform.trim.sin_cos();
        (x, z) = rot_down(x, z, platform.pivot_x, zw, s, c);
    }
    (x, z - zw)
}

/// A hull's centreplane profile at the solved attitude: the keel/rocker
/// line, the top of the lofted band, the design-waterline mark, and the
/// (band-top, keel) pairs that bound the hull's silhouette.
struct HullProfile {
    keel: Vec<[f64; 2]>,
    /// The top of the body's modelled band (`z_b = 0`). This is where the
    /// loft was cut, NOT a sheer line: `michell loft --band` deliberately
    /// stops the fit a short way above the design waterline, because a deck
    /// is a cliff for a height-field loft. The source CAD may well carry
    /// full topsides above it; the `.hull` body simply does not. Drawn, and
    /// labelled as the band top, so a pose that immerses past it is visible
    /// rather than silently missing.
    band_top: Vec<[f64; 2]>,
    design_wl: Vec<[f64; 2]>,
    /// `[band top, keel]` pairs at the stations that carry any beam.
    silhouette: Vec<[[f64; 2]; 2]>,
    /// Band top's height above the design waterline [m], for the caption.
    band_above_wl: f64,
}

fn hull_profile_at(body: &Body, pose: &HullPose, platform: &Platform) -> HullProfile {
    const STATIONS: usize = 120;
    const SCAN: usize = 160;
    const SMOOTH_PASSES: usize = 3;
    let (x0, x1) = body.surface().x_domain();
    let (_, depth) = body.surface().z_domain();
    let pivot_x = pose.pivot_x.unwrap_or(0.5 * (x0 + x1));
    let xs: Vec<f64> = (0..STATIONS)
        .map(|i| x0 + (x1 - x0) * i as f64 / (STATIONS - 1) as f64)
        .collect();
    let z_of = |j: usize| depth * j as f64 / (SCAN - 1) as f64;

    // Raw keel depth per station, in the body's own frame.
    let mut raw: Vec<Option<f64>> = Vec::with_capacity(STATIONS);
    let mut fs = vec![0.0f64; SCAN];
    for &xb in &xs {
        let f_wl = body.surface().eval(xb, body.waterline()).max(0.0);
        let eps = (KEEL_BEAM_FRACTION * f_wl).max(KEEL_MIN_BEAM);
        for (j, f) in fs.iter_mut().enumerate() {
            *f = body.surface().eval(xb, z_of(j));
        }
        // Walk DOWN from the band top and stop at the first crossing back
        // below eps. Taking the DEEPEST crossing instead - the obvious
        // reading of "where does the hull stop" - lands in the ringing tail
        // this loft leaves under the hull, a few millimetres of half-beam
        // wandering all the way to the band bottom. The detected depth then
        // hops between ripple lobes from one station to the next and the
        // keel draws as a sawtooth.
        let mut entered = false;
        let mut keel = None;
        for j in 0..SCAN {
            if !entered {
                entered = fs[j] > eps;
            } else if fs[j] <= eps {
                let t = ((fs[j - 1] - eps) / (fs[j - 1] - fs[j])).clamp(0.0, 1.0);
                keel = Some(z_of(j - 1) + t * (z_of(j) - z_of(j - 1)));
                break;
            }
        }
        raw.push(match (entered, keel) {
            (true, None) => Some(depth), // beam all the way to the band bottom
            (_, k) => k,
        });
    }

    // Light binomial smoothing of the detected depths. A hull's underside is
    // smooth; what survives the crossing rule above is detection jitter of a
    // millimetre or two, which this panel's vertical exaggeration magnifies
    // into a visible sawtooth. The window shrinks at the ends so the transom
    // and the stem stay exactly where they were found rather than being
    // rounded off by the filter.
    let mut zs: Vec<f64> = raw.iter().filter_map(|v| *v).collect();
    let n = zs.len();
    for _ in 0..SMOOTH_PASSES {
        let o = zs.clone();
        for i in 0..n {
            zs[i] = match i.min(n - 1 - i).min(2) {
                0 => o[i],
                1 => 0.25 * (o[i - 1] + 2.0 * o[i] + o[i + 1]),
                _ => {
                    (o[i - 2] + 4.0 * o[i - 1] + 6.0 * o[i] + 4.0 * o[i + 1] + o[i + 2]) / 16.0
                }
            };
        }
    }

    let mut keel = Vec::with_capacity(STATIONS);
    let mut band_top = Vec::with_capacity(STATIONS);
    let mut design_wl = Vec::with_capacity(STATIONS);
    let mut silhouette = Vec::with_capacity(STATIONS);
    let mut k = 0usize;
    for (i, &xb) in xs.iter().enumerate() {
        let top = body_to_water(xb, 0.0, body, pose, pivot_x, platform);
        band_top.push([top.0, top.1]);
        if raw[i].is_some() {
            let (xw, zw) = body_to_water(xb, zs[k], body, pose, pivot_x, platform);
            k += 1;
            keel.push([xw, zw]);
            silhouette.push([[top.0, top.1], [xw, zw]]);
        }
        let (xw, zw) = body_to_water(xb, body.waterline(), body, pose, pivot_x, platform);
        design_wl.push([xw, zw]);
    }
    HullProfile {
        keel,
        band_top,
        design_wl,
        silhouette,
        band_above_wl: body.waterline(),
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_profile_view(
    page: &mut Page,
    rect: [f64; 4],
    hulls: &[MHull],
    poses: &[HullPose],
    platform: &Platform,
    members: &[(&Hull, Placement)],
    cond: &Conditions,
    wave_cut_y: f64,
    band_exceeded: usize,
) -> Result<String, String> {
    let profiles: Vec<HullProfile> = hulls
        .iter()
        .zip(poses)
        .map(|(h, p)| hull_profile_at(&h.body, p, platform))
        .collect();

    // Wave elevation cut along the fleet centreline at the reference hull's
    // transverse offset ("the wake, in profile"); z is up here, so it adds
    // to the flat sea surface as a negative depth.
    let mut wave_cut: Vec<[f64; 2]> = Vec::new();
    if !members.is_empty() {
        let mut spec = FreeWaveSpectrum::new_with_transom(members, cond, TransomClosure::None)
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
    // band-top curve spans the hull's full modelled length regardless of
    // local wetness, so it sets the x-range.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut z_lo = 0.0f64; // includes the still-water line
    let mut z_hi = 0.0f64;
    for p in &profiles {
        for pt in p.keel.iter().chain(&p.band_top).chain(&p.design_wl) {
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
    let to_page = |x: f64, z: f64| -> [f64; 2] {
        [rect[0] + (x - x_lo) * sx, y_mid - (z - z_mid) * sz]
    };

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
        let top_pts: Vec<[f64; 2]> = p.band_top.iter().map(|pt| to_page(pt[0], pt[1])).collect();
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

    let band_h = profiles.first().map(|p| p.band_above_wl).unwrap_or(0.0);
    let mut caption = format!(
        "shaded: hull centreplane section, darker = submerged; black: keel, the {:.0}% \
         half-beam contour lightly smoothed (this loft has no sharp keel edge); \
         orange: design waterline\n\
         blue: still water and the wake cut at y = {wave_cut_y:.2} m; gray: top of the \
         lofted band, {band_h:.2} m above the design waterline - topsides above it are \
         not in the .hull file; vertical exaggeration {:.1}x",
        100.0 * KEEL_BEAM_FRACTION,
        sz / sx
    );
    if band_exceeded > 0 {
        caption.push_str(&format!(
            "\nWARNING: the water rose above the band top at {band_exceeded} wetted \
             sample(s) - that hull is missing from the file and was counted as zero \
             beam, so the immersed hull shown (and this row's forces) are understated"
        ));
    }
    Ok(caption)
}

fn data_x_range(profiles: &[HullProfile]) -> (f64, f64) {
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    for p in profiles {
        for pt in p.band_top.iter() {
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

    let cols = ["row", "point", "U [m/s]", "Fn", "sinkage [m]", "trim [deg]", "Rt [N]", "Pe [W]"];
    let col_x = [
        MARGIN,
        MARGIN + 40.0,
        MARGIN + 220.0,
        MARGIN + 300.0,
        MARGIN + 360.0,
        MARGIN + 450.0,
        MARGIN + 540.0,
        MARGIN + 620.0,
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
        page.text(col_x[3], y, 9.0, BLACK, &format!("{:.3}", row.froude));
        page.text(col_x[4], y, 9.0, BLACK, &format!("{:.4}", row.sinkage));
        page.text(col_x[5], y, 9.0, BLACK, &format!("{:.3}", row.trim_deg));
        if let Some(rt) = row.rt {
            page.text(col_x[6], y, 9.0, BLACK, &format!("{rt:.1}"));
        }
        if let Some(pe) = row.pe {
            page.text(col_x[7], y, 9.0, BLACK, &format!("{pe:.1}"));
        }
        if row.band_exceeded > 0 {
            page.text(col_x[7] + 56.0, y, 9.0, ORANGE, "*");
        }
    }

    // A row whose pose lifted water above the lofted band was integrated over
    // a hull that is missing its immersed upper stern, so its Rt and Pe are
    // understated. Say so next to the numbers, not only on the detail page.
    let flagged = rows.iter().filter(|r| r.band_exceeded > 0).count();
    if flagged > 0 {
        y -= 24.0;
        page.text(
            MARGIN,
            y,
            8.0,
            ORANGE,
            &format!(
                "* {flagged} of {} row(s): the solved attitude immersed the hull above the \
                 top of the lofted band.",
                rows.len()
            ),
        );
        page.text(
            MARGIN,
            y - 10.0,
            8.0,
            GRAY,
            "That geometry is not in the .hull file and was taken as zero half-beam, so Rt \
             and Pe are understated. Re-loft with a taller --band to close the gap.",
        );
    }
    page
}

//! `michell report`: run a dynamic-squat speed sweep from a JSON manifest and
//! render a self-contained PDF — an index page plus, per sweep row, a
//! plan-view wake image and a profile view of sinkage/trim (with a wake
//! elevation cut along the fleet centreline). One CLI invocation, no
//! external tooling and no per-run agent orchestration: the same manifest
//! schema `michell sweep` reads, via [`crate::manifest::parse_manifest`].

use crate::manifest::{parse_manifest, point_state, Axis, MHull, PointState};
use crate::pdf::{Document, Page};
use crate::png;
use michell::body::Body;
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

pub fn run(manifest_path: &str, out_path: &str) -> Result<(), String> {
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
            eprintln!(
                "row {}/{total_rows}: solving {prefix}U = {u:.3} m/s (Fn {froude:.3})...",
                rows.len() + 1
            );
            let cond = pm.fluid.make_cond(u)?;
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
            warm = Some((dyn_eq.sinkage, dyn_eq.trim));

            let members: Vec<(&Hull, Placement)> =
                dyn_eq.fleet.members.iter().map(|(h, p)| (h, *p)).collect();
            let resistance = if members.is_empty() {
                None
            } else {
                Some(
                    multihull_resistance_with(&members, &cond, &pm.wave_opts, pm.form_factor)
                        .map_err(|e| format!("point {} U={u}: {e}", point + 1))?,
                )
            };

            let row_no = rows.len() + 1;
            let platform = Platform {
                sinkage: dyn_eq.sinkage,
                trim: dyn_eq.trim,
                pivot_x,
            };
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
            )?;
            rows.push(RowSummary {
                row_no,
                point_label: point_label.clone(),
                speed: u,
                froude,
                sinkage: dyn_eq.sinkage,
                trim_deg: dyn_eq.trim.to_degrees(),
                rt: resistance.as_ref().map(|r| r.total),
                pe: resistance.as_ref().map(|r| r.effective_power),
            });
            detail_pages.push(page);
            eprintln!(
                "row {row_no}/{total_rows} done in {:.1}s: sinkage {:.4} m, trim {:.3} deg",
                row_started.elapsed().as_secs_f64(),
                dyn_eq.sinkage,
                dyn_eq.trim.to_degrees()
            );
        }

        for (i, a) in pm.axes.iter().enumerate().rev() {
            idx[i] += 1;
            if idx[i] < a.values.len() {
                break;
            }
            idx[i] = 0;
        }
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

    let plan_rect = [MARGIN, PAGE_H - 300.0, PAGE_W - 2.0 * MARGIN, 232.0];
    let profile_rect = [MARGIN, 56.0, PAGE_W - 2.0 * MARGIN, 200.0];

    let plan_caption = draw_plan_view(doc, &mut page, plan_rect, members, cond, l_ref)?;
    page.text(plan_rect[0], plan_rect[1] - 12.0, 8.0, GRAY, &plan_caption);
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
    )?;
    page.text(profile_rect[0], profile_rect[1] - 12.0, 8.0, GRAY, &profile_caption);
    page.text(
        profile_rect[0],
        profile_rect[1] + profile_rect[3] + 4.0,
        9.0,
        BLACK,
        "Profile: sinkage/trim (keel and design waterline) with a wave elevation cut",
    );

    Ok(page)
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
    let x1 = x_hi + 0.35 * l_ref;
    let x0 = x_lo - 3.0 * l_ref;
    let yh = (0.42 * (x1 - x0)).max(y_abs + 0.8 * l_ref);
    let (y0, y1) = (-yh, yh);

    let nx = 520usize;
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

    const FADE_TOWARD: [u8; 3] = [0xf0, 0xef, 0xec];
    const FADE_FRACTION: f64 = 0.55;
    let mut rgb = vec![0u8; 3 * nx * ny];
    for iy in 0..ny {
        let row = ny - 1 - iy; // row 0 = top = +y edge, matching Document::image's convention
        for ix in 0..nx {
            let t = grid.get(ix, iy) / vmax;
            let mut c = png::diverging(t);
            if grid.x(ix) > x_lo {
                c = png::fade(c, FADE_TOWARD, FADE_FRACTION);
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
        "zeta +-{vmax:.4} m; x {x0:.1}..{x1:.1} m, y {y0:.1}..{y1:.1} m{}",
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

/// A hull's keel/rocker line and design-waterline mark, mapped into the
/// water frame at the solved sinkage/trim.
struct HullProfile {
    keel: Vec<[f64; 2]>,
    /// The top of the body's modelled band (`z_b = 0`). Not a real sheer or
    /// deck line: the source geometry for this crate's hulls typically stops
    /// a few centimetres above the design waterline (it was exported as the
    /// wetted hull only, with a small margin for the sinkage/trim range of
    /// interest, not as a full topsides model). Kept only to bound the
    /// x-range for the wave-elevation cut; deliberately not drawn.
    deck: Vec<[f64; 2]>,
    design_wl: Vec<[f64; 2]>,
}

fn hull_profile_at(body: &Body, pose: &HullPose, platform: &Platform) -> HullProfile {
    const STATIONS: usize = 60;
    const SCAN: usize = 48;
    let (x0, x1) = body.surface().x_domain();
    let (_, depth) = body.surface().z_domain();
    let pivot_x = pose.pivot_x.unwrap_or(0.5 * (x0 + x1));
    let eps = 1e-4 * body.surface().eval(0.5 * (x0 + x1), body.waterline()).max(1e-6);

    let mut keel = Vec::with_capacity(STATIONS);
    let mut deck = Vec::with_capacity(STATIONS);
    let mut design_wl = Vec::with_capacity(STATIONS);
    for i in 0..STATIONS {
        let xb = x0 + (x1 - x0) * i as f64 / (STATIONS - 1) as f64;
        let mut bottom = None;
        for j in 0..SCAN {
            let zb = depth * j as f64 / (SCAN - 1) as f64;
            if body.surface().eval(xb, zb) > eps {
                bottom = Some(zb);
            }
        }
        if let Some(zb) = bottom {
            let (xw, zw) = body_to_water(xb, zb, body, pose, pivot_x, platform);
            keel.push([xw, zw]);
        }
        let (xw, zw) = body_to_water(xb, 0.0, body, pose, pivot_x, platform);
        deck.push([xw, zw]);
        let (xw, zw) = body_to_water(xb, body.waterline(), body, pose, pivot_x, platform);
        design_wl.push([xw, zw]);
    }
    HullProfile {
        keel,
        deck,
        design_wl,
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
    // "deck" curve isn't drawn (see HullProfile's doc comment: it is not a
    // real sheer line) but still contributes to the x-range, since it spans
    // the hull's full modelled length regardless of local wetness.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut z_lo = 0.0f64; // includes the still-water line
    let mut z_hi = 0.0f64;
    for p in &profiles {
        for pt in p.keel.iter().chain(&p.deck).chain(&p.design_wl) {
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
    // Still water: exactly flat by construction (sinkage/trim are baked
    // into the hull curves below), always drawn across the full panel width.
    page.polyline(
        &[to_page(x_lo, 0.0), to_page(x_hi, 0.0)],
        BLUE,
        0.5,
    );

    for p in &profiles {
        let keel_pts: Vec<[f64; 2]> = p.keel.iter().map(|pt| to_page(pt[0], pt[1])).collect();
        page.polyline(&keel_pts, BLACK, 1.2);
        let wl_pts: Vec<[f64; 2]> = p.design_wl.iter().map(|pt| to_page(pt[0], pt[1])).collect();
        page.polyline(&wl_pts, ORANGE, 0.9);
    }

    if wave_cut.len() > 1 {
        let cut_pts: Vec<[f64; 2]> = wave_cut.iter().map(|pt| to_page(pt[0], pt[1])).collect();
        page.polyline(&cut_pts, BLUE, 0.9);
    }

    Ok(format!(
        "black: keel/rocker (modelled draft only - the source geometry \
         has no topsides above the waterline); orange: design waterline; \
         blue: still water and wake elevation cut at y = {wave_cut_y:.2} m; \
         vertical exaggeration {:.1}x",
        sz / sx
    ))
}

fn data_x_range(profiles: &[HullProfile]) -> (f64, f64) {
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    for p in profiles {
        for pt in p.deck.iter() {
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
    }
    page
}

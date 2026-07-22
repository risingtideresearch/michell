//! Sweep PDF report: an index page whose table rows hyperlink to a page per
//! sweep row, each showing the fleet in profile and body-plan (front) views
//! against the waterline, the righting-arm (GZ) curve, a plan view over the
//! wave-amplitude field, and a table of the row's values.
//!
//! Geometry arrives already reduced to fleet-frame polylines (z up, waterline
//! at z = 0) so this module never touches the hull maths; it only lays ink on
//! the page via [`crate::pdf`].

use crate::pdf::{text_width, Pdf};
use crate::png::diverging;
use michell::WaveGrid;

// US Letter, landscape.
const PAGE_W: f64 = 792.0;
const PAGE_H: f64 = 612.0;
const MARGIN: f64 = 36.0;

const INK: [f64; 3] = [0.12, 0.12, 0.14];
const GRIDLINE: [f64; 3] = [0.82, 0.82, 0.85];
const PANEL_EDGE: [f64; 3] = [0.55, 0.55, 0.60];
const WATERLINE: [f64; 3] = [0.10, 0.42, 0.72];
const HULL_FILL: [f64; 3] = [0.80, 0.85, 0.90];
const HULL_EDGE: [f64; 3] = [0.30, 0.34, 0.40];
const GZ_LINE: [f64; 3] = [0.12, 0.45, 0.72];
const MARKER: [f64; 3] = [0.80, 0.18, 0.18];
const LINK_BLUE: [f64; 3] = [0.10, 0.30, 0.65];

/// A hull reduced to fleet-frame outlines (metres; z up, waterline at 0).
#[derive(Clone)]
pub struct HullGeom {
    /// Closed underwater profile silhouette in the x–z plane.
    pub profile: Vec<(f64, f64)>,
    /// Body-plan half-sections, each a `(y, z)` polyline.
    pub sections: Vec<Vec<(f64, f64)>>,
    /// Closed waterplane outline in the x–y plane.
    pub waterplane: Vec<(f64, f64)>,
}

/// A colour-mapped wave-elevation raster plus its fleet-frame extent.
pub struct WavePlan {
    pub rgb: Vec<u8>,
    pub nx: usize,
    pub ny: usize,
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub vmax: f64,
    pub resolution_limited: bool,
}

impl WavePlan {
    /// Colour-map a wave grid with the diverging heatmap (row 0 = +y edge),
    /// saturating at the 99.5th percentile of |ζ|.
    pub fn from_grid(grid: &WaveGrid) -> WavePlan {
        let mut abs: Vec<f64> = grid.zeta.iter().map(|v| v.abs()).collect();
        abs.sort_by(|a, b| a.total_cmp(b));
        let vmax = if abs.is_empty() {
            1e-12
        } else {
            abs[((abs.len() - 1) as f64 * 0.995) as usize].max(1e-12)
        };
        let (nx, ny) = (grid.nx, grid.ny);
        let mut rgb = vec![0u8; 3 * nx * ny];
        for iy in 0..ny {
            let row = ny - 1 - iy; // image row 0 = +y edge (top)
            for ix in 0..nx {
                let c = diverging(grid.get(ix, iy) / vmax);
                let o = 3 * (row * nx + ix);
                rgb[o..o + 3].copy_from_slice(&c);
            }
        }
        WavePlan {
            rgb,
            nx,
            ny,
            x0: grid.x0,
            x1: grid.x1,
            y0: grid.y0,
            y1: grid.y1,
            vmax,
            resolution_limited: grid.resolution_limited,
        }
    }
}

/// Everything one sweep row contributes to the report.
pub struct Row {
    /// Short page heading (e.g. "Row 3").
    pub label: String,
    /// Longer subtitle describing the axis coordinates.
    pub subtitle: String,
    /// Value table: `(name, formatted)` pairs.
    pub values: Vec<(String, String)>,
    pub hulls: Vec<HullGeom>,
    /// GZ curve samples `(heel_deg, gz_m)`, if a righting arm is defined.
    pub gz_curve: Option<Vec<(f64, f64)>>,
    /// Operating heel highlighted on the GZ curve [deg].
    pub heel_deg: f64,
    pub wave: Option<WavePlan>,
}

/// Build the whole report document.
pub fn build(title: &str, header: &[String], table: &[Vec<String>], rows: &[Row]) -> Vec<u8> {
    let mut pdf = Pdf::new();

    // The index paginates; reserve its pages so row-page indices are known.
    let rows_per_page = 30usize;
    let n_index = table.len().div_ceil(rows_per_page).max(1);
    let index_pages: Vec<usize> = (0..n_index).map(|_| pdf.add_page(PAGE_W, PAGE_H)).collect();

    // Per-row pages.
    let mut row_pages = Vec::with_capacity(rows.len());
    for row in rows {
        let p = pdf.add_page(PAGE_W, PAGE_H);
        draw_row_page(&mut pdf, p, index_pages[0], row);
        row_pages.push(p);
    }

    // Fill the index pages now that targets are known.
    for (k, &page) in index_pages.iter().enumerate() {
        let lo = k * rows_per_page;
        let hi = (lo + rows_per_page).min(table.len());
        draw_index_page(
            &mut pdf,
            page,
            title,
            header,
            &table[lo..hi],
            &row_pages[lo..hi],
            k + 1,
            n_index,
        );
    }

    pdf.finish()
}

#[allow(clippy::too_many_arguments)]
fn draw_index_page(
    pdf: &mut Pdf,
    page: usize,
    title: &str,
    header: &[String],
    rows: &[Vec<String>],
    targets: &[usize],
    page_no: usize,
    n_pages: usize,
) {
    let ncol = header.len();
    let x_left = MARGIN;
    let x_right = PAGE_W - MARGIN;
    let usable = x_right - x_left;

    // Natural column widths from the widest cell, then scale to fit.
    let fs = 7.0;
    let pad = 5.0;
    let mut widths = vec![0.0f64; ncol];
    for (c, h) in header.iter().enumerate() {
        widths[c] = text_width(h, fs, true);
    }
    for r in rows {
        for (c, cell) in r.iter().enumerate().take(ncol) {
            widths[c] = widths[c].max(text_width(cell, fs, false));
        }
    }
    let natural: f64 = widths.iter().map(|w| w + 2.0 * pad).sum();
    let scale = (usable / natural).min(1.0);
    let widths: Vec<f64> = widths.iter().map(|w| (w + 2.0 * pad) * scale).collect();
    let total_w: f64 = widths.iter().sum();

    let pg = pdf.page(page);
    if page_no == 1 {
        pg.text(x_left, PAGE_H - 34.0, 15.0, INK, true, title);
        pg.text(
            x_left,
            PAGE_H - 50.0,
            8.5,
            [0.4, 0.4, 0.45],
            false,
            "Each row links to its detail page: profile & body-plan views, GZ curve, and plan-view wave field.",
        );
    }
    if n_pages > 1 {
        pg.text_right(
            x_right,
            PAGE_H - 34.0,
            9.0,
            [0.4, 0.4, 0.45],
            false,
            &format!("index {page_no}/{n_pages}"),
        );
    }

    let row_h = 13.0;
    let top = PAGE_H - 64.0;

    // Header band.
    pg.rect(x_left, top - row_h, total_w, row_h, 0.0, None, Some([0.90, 0.92, 0.96]));
    let mut x = x_left;
    for (c, h) in header.iter().enumerate() {
        pg.text(x + pad, top - row_h + 4.0, fs, INK, true, h);
        x += widths[c];
    }
    pg.line(x_left, top - row_h, x_left + total_w, top - row_h, 0.6, PANEL_EDGE);

    // Data rows.
    for (i, r) in rows.iter().enumerate() {
        let y = top - row_h * (i + 2) as f64;
        if i % 2 == 1 {
            pg.rect(x_left, y, total_w, row_h, 0.0, None, Some([0.965, 0.965, 0.98]));
        }
        let mut x = x_left;
        for (c, w) in widths.iter().enumerate() {
            let cell = r.get(c).map(String::as_str).unwrap_or("");
            let color = if c == 0 { LINK_BLUE } else { INK };
            let bold = c == 0;
            pg.text(x + pad, y + 4.0, fs, color, bold, cell);
            x += w;
        }
        // The whole row is a hyperlink to its detail page.
        pg.link([x_left, y, x_left + total_w, y + row_h], targets[i]);
    }
    let bottom = top - row_h * (rows.len() + 1) as f64;
    pg.rect(x_left, bottom, total_w, top - bottom, 0.6, Some(PANEL_EDGE), None);
}

fn draw_row_page(pdf: &mut Pdf, page: usize, index_page: usize, row: &Row) {
    // Header with a back-link to the index.
    {
        let pg = pdf.page(page);
        pg.text(MARGIN, PAGE_H - 30.0, 14.0, INK, true, &row.label);
        pg.text(MARGIN, PAGE_H - 44.0, 8.5, [0.4, 0.4, 0.45], false, &row.subtitle);
        let back = "\u{2190} index";
        let w = text_width(back, 9.0, false);
        pg.text(PAGE_W - MARGIN - w, PAGE_H - 30.0, 9.0, LINK_BLUE, false, "< index");
        pg.link(
            [PAGE_W - MARGIN - w - 4.0, PAGE_H - 34.0, PAGE_W - MARGIN, PAGE_H - 22.0],
            index_page,
        );
    }

    // Layout: a value-table column on the right; a 2x2 grid of panels on the
    // left (profile over plan; front over GZ) so the two longitudinal views
    // share a horizontal axis.
    let table_w = 172.0;
    let content_top = PAGE_H - 58.0;
    let content_bot = MARGIN;
    let gap = 14.0;
    let plots_x0 = MARGIN;
    let plots_x1 = PAGE_W - MARGIN - table_w - gap;
    let col_w = (plots_x1 - plots_x0 - gap) / 2.0;
    let row_h = (content_top - content_bot - gap) / 2.0;

    let profile = (plots_x0, content_top - row_h, col_w, row_h);
    let plan = (plots_x0, content_bot, col_w, row_h);
    let front = (plots_x0 + col_w + gap, content_top - row_h, col_w, row_h);
    let gz = (plots_x0 + col_w + gap, content_bot, col_w, row_h);

    draw_profile(pdf, page, profile, row);
    draw_front(pdf, page, front, row);
    draw_plan(pdf, page, plan, row);
    draw_gz(pdf, page, gz, row);
    draw_value_table(pdf, page, PAGE_W - MARGIN - table_w, content_bot, table_w, content_top - content_bot, row);
}

/// Draw a titled panel; return the inner drawing rectangle `(x, y, w, h)`.
fn panel(pdf: &mut Pdf, page: usize, rect: (f64, f64, f64, f64), title: &str) -> (f64, f64, f64, f64) {
    let (x, y, w, h) = rect;
    let pg = pdf.page(page);
    pg.rect(x, y, w, h, 0.7, Some(PANEL_EDGE), Some([1.0, 1.0, 1.0]));
    pg.rect(x, y + h - 16.0, w, 16.0, 0.0, None, Some([0.93, 0.94, 0.97]));
    pg.text(x + 6.0, y + h - 12.0, 8.5, INK, true, title);
    pg.line(x, y + h - 16.0, x + w, y + h - 16.0, 0.6, PANEL_EDGE);
    (x + 8.0, y + 8.0, w - 16.0, h - 16.0 - 12.0)
}

/// Fit a data-space bounding box into an inner rect with equal aspect ratio,
/// returning a closure mapping (data x, data y) -> (page x, page y).
fn fit_equal(
    inner: (f64, f64, f64, f64),
    bbox: (f64, f64, f64, f64),
) -> impl Fn(f64, f64) -> (f64, f64) {
    let (ix, iy, iw, ih) = inner;
    let (bx0, bx1, by0, by1) = bbox;
    let dw = (bx1 - bx0).max(1e-9);
    let dh = (by1 - by0).max(1e-9);
    let s = (iw / dw).min(ih / dh);
    // Centre the drawing in the panel.
    let ox = ix + 0.5 * (iw - s * dw);
    let oy = iy + 0.5 * (ih - s * dh);
    move |x: f64, y: f64| (ox + (x - bx0) * s, oy + (y - by0) * s)
}

fn hull_bbox_xz(row: &Row) -> (f64, f64, f64, f64) {
    let (mut x0, mut x1, mut z0, mut z1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for hg in &row.hulls {
        for &(x, z) in &hg.profile {
            x0 = x0.min(x);
            x1 = x1.max(x);
            z0 = z0.min(z);
            z1 = z1.max(z);
        }
    }
    if x0 > x1 {
        return (-1.0, 1.0, -1.0, 0.1);
    }
    z1 = z1.max(0.0); // always show the waterline
    // A little vertical breathing room.
    let pad = 0.06 * (z1 - z0).max(1e-6);
    (x0, x1, z0 - pad, z1 + pad)
}

fn draw_profile(pdf: &mut Pdf, page: usize, rect: (f64, f64, f64, f64), row: &Row) {
    let inner = panel(pdf, page, rect, "Profile (looking to port)");
    if row.hulls.is_empty() {
        return;
    }
    let (x0, x1, z0, z1) = hull_bbox_xz(row);
    let map = fit_equal(inner, (x0, x1, z0, z1));
    let pg = pdf.page(page);
    for hg in &row.hulls {
        let pts: Vec<(f64, f64)> = hg.profile.iter().map(|&(x, z)| map(x, z)).collect();
        pg.fill_poly(&pts, HULL_FILL);
        pg.polyline(&pts, 0.8, HULL_EDGE, true);
    }
    // Waterline at z = 0.
    let (wx0, wy) = map(x0, 0.0);
    let (wx1, _) = map(x1, 0.0);
    pg.line(wx0, wy, wx1, wy, 1.0, WATERLINE);
    axis_note(pg, inner, "x (m)", "z (m)");
}

fn draw_front(pdf: &mut Pdf, page: usize, rect: (f64, f64, f64, f64), row: &Row) {
    let inner = panel(pdf, page, rect, "Body plan (looking forward)");
    if row.hulls.is_empty() {
        return;
    }
    let (mut y0, mut y1, mut z0, mut z1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for hg in &row.hulls {
        for sec in &hg.sections {
            for &(y, z) in sec {
                y0 = y0.min(y);
                y1 = y1.max(y);
                z0 = z0.min(z);
                z1 = z1.max(z);
            }
        }
    }
    if y0 > y1 {
        return;
    }
    z1 = z1.max(0.0);
    let padz = 0.06 * (z1 - z0).max(1e-6);
    let pady = 0.06 * (y1 - y0).max(1e-6);
    let map = fit_equal(inner, (y0 - pady, y1 + pady, z0 - padz, z1 + padz));
    let pg = pdf.page(page);
    for hg in &row.hulls {
        for sec in &hg.sections {
            let pts: Vec<(f64, f64)> = sec.iter().map(|&(y, z)| map(y, z)).collect();
            pg.polyline(&pts, 0.6, HULL_EDGE, false);
        }
    }
    let (wx0, wy) = map(y0 - pady, 0.0);
    let (wx1, _) = map(y1 + pady, 0.0);
    pg.line(wx0, wy, wx1, wy, 1.0, WATERLINE);
    axis_note(pg, inner, "y (m)", "z (m)");
}

fn draw_plan(pdf: &mut Pdf, page: usize, rect: (f64, f64, f64, f64), row: &Row) {
    let inner = panel(pdf, page, rect, "Plan view over wave field");
    let (ix, iy, _iw, ih) = inner;
    let Some(wave) = &row.wave else {
        let pg = pdf.page(page);
        pg.text(ix + 6.0, iy + ih - 14.0, 8.0, [0.5, 0.5, 0.55], false, "no wave field");
        return;
    };
    // Fit the wave extent into the panel with equal aspect; image fills that.
    let map = fit_equal(inner, (wave.x0, wave.x1, wave.y0, wave.y1));
    let (px0, py0) = map(wave.x0, wave.y0);
    let (px1, py1) = map(wave.x1, wave.y1);
    let id = pdf.add_image(wave.nx, wave.ny, wave.rgb.clone());
    let pg = pdf.page(page);
    pg.image(id, px0, py0, px1 - px0, py1 - py0);
    // Overlay hull waterplanes in the same mapping.
    for hg in &row.hulls {
        let pts: Vec<(f64, f64)> = hg.waterplane.iter().map(|&(x, y)| map(x, y)).collect();
        pg.polyline(&pts, 0.7, INK, true);
    }
    pg.rect(px0, py0, px1 - px0, py1 - py0, 0.6, Some(PANEL_EDGE), None);
    let note = if wave.resolution_limited {
        format!("zeta max {:.3} m (grid-limited); ship advances toward +x", wave.vmax)
    } else {
        format!("zeta max {:.3} m; ship advances toward +x", wave.vmax)
    };
    pg.text(ix, iy - 1.0, 6.5, [0.5, 0.5, 0.55], false, &note);
}

fn draw_gz(pdf: &mut Pdf, page: usize, rect: (f64, f64, f64, f64), row: &Row) {
    let inner = panel(pdf, page, rect, "Righting arm GZ");
    let (ix, iy, iw, ih) = inner;
    let Some(curve) = &row.gz_curve else {
        let pg = pdf.page(page);
        pg.text(
            ix + 6.0,
            iy + ih - 14.0,
            8.0,
            [0.5, 0.5, 0.55],
            false,
            "GZ needs a weight axis",
        );
        return;
    };
    if curve.len() < 2 {
        return;
    }
    let phi_max = curve.iter().map(|c| c.0).fold(0.0f64, f64::max).max(1e-6);
    let mut gz_lo = 0.0f64;
    let mut gz_hi = 0.0f64;
    for &(_, gz) in curve {
        gz_lo = gz_lo.min(gz);
        gz_hi = gz_hi.max(gz);
    }
    let span = (gz_hi - gz_lo).max(1e-4);
    gz_hi += 0.12 * span;
    gz_lo -= 0.12 * span;
    // Plot area with room for labels.
    let pl = (ix + 26.0, iy + 16.0, iw - 32.0, ih - 24.0);
    let map = |phi: f64, gz: f64| -> (f64, f64) {
        (
            pl.0 + phi / phi_max * pl.2,
            pl.1 + (gz - gz_lo) / (gz_hi - gz_lo) * pl.3,
        )
    };
    let pg = pdf.page(page);
    // Zero line and axes.
    let (zx0, zy) = map(0.0, 0.0);
    let (zx1, _) = map(phi_max, 0.0);
    pg.line(zx0, zy, zx1, zy, 0.5, GRIDLINE);
    pg.line(pl.0, pl.1, pl.0, pl.1 + pl.3, 0.6, PANEL_EDGE);
    pg.line(pl.0, pl.1, pl.0 + pl.2, pl.1, 0.6, PANEL_EDGE);
    // Curve.
    let pts: Vec<(f64, f64)> = curve.iter().map(|&(p, g)| map(p, g)).collect();
    pg.polyline(&pts, 1.2, GZ_LINE, false);
    // Operating-heel marker.
    if row.heel_deg.abs() > 1e-6 && row.heel_deg <= phi_max + 1e-6 {
        // Interpolate GZ at the operating heel.
        let phi = row.heel_deg;
        let mut gz = curve[0].1;
        for w in curve.windows(2) {
            if phi >= w[0].0 && phi <= w[1].0 {
                let t = (phi - w[0].0) / (w[1].0 - w[0].0).max(1e-9);
                gz = w[0].1 + t * (w[1].1 - w[0].1);
                break;
            }
        }
        let (mx, my) = map(phi, gz);
        pg.rect(mx - 1.8, my - 1.8, 3.6, 3.6, 0.0, None, Some(MARKER));
    }
    // Tick labels.
    pg.text(pl.0, pl.1 - 9.0, 6.5, INK, false, "0");
    pg.text_right(pl.0 + pl.2, pl.1 - 9.0, 6.5, INK, false, &format!("{phi_max:.0}\u{00b0}"));
    pg.text_center(pl.0 + 0.5 * pl.2, pl.1 - 9.0, 6.5, [0.45, 0.45, 0.5], false, "heel");
    pg.text(ix, pl.1 + pl.3 - 3.0, 6.5, INK, false, &format!("{gz_hi:.3} m"));
    pg.text(ix, pl.1 - 3.0, 6.5, INK, false, &format!("{gz_lo:.3} m"));
}

fn draw_value_table(
    pdf: &mut Pdf,
    page: usize,
    x: f64,
    y: f64,
    w: f64,
    h: f64,
    row: &Row,
) {
    let pg = pdf.page(page);
    pg.rect(x, y, w, h, 0.7, Some(PANEL_EDGE), Some([1.0, 1.0, 1.0]));
    pg.rect(x, y + h - 16.0, w, 16.0, 0.0, None, Some([0.93, 0.94, 0.97]));
    pg.text(x + 6.0, y + h - 12.0, 8.5, INK, true, "Values");
    pg.line(x, y + h - 16.0, x + w, y + h - 16.0, 0.6, PANEL_EDGE);
    let fs = 7.5;
    let line_h = ((h - 24.0) / row.values.len().max(1) as f64).min(13.0);
    let mut cy = y + h - 22.0;
    for (i, (k, v)) in row.values.iter().enumerate() {
        cy -= line_h;
        if i % 2 == 1 {
            pg.rect(x + 1.0, cy - 1.5, w - 2.0, line_h, 0.0, None, Some([0.965, 0.965, 0.98]));
        }
        pg.text(x + 6.0, cy, fs, [0.4, 0.4, 0.45], false, k);
        pg.text_right(x + w - 6.0, cy, fs, INK, false, v);
    }
}

fn axis_note(pg: &mut crate::pdf::Page, inner: (f64, f64, f64, f64), xl: &str, yl: &str) {
    let (ix, iy, iw, ih) = inner;
    pg.text_right(ix + iw, iy - 1.0, 6.5, [0.5, 0.5, 0.55], false, xl);
    pg.text(ix, iy + ih - 6.0, 6.5, [0.5, 0.5, 0.55], false, yl);
}

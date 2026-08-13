//! An egui view of a `.msw` sweep archive. Decodes the self-contained binary
//! format (via `michell_cli::archive`) and renders it: the study metadata, a
//! scatter plot of any column against any other across all rows, a selectable
//! rows table, and — for the selected row — the righting-arm (GZ) curve and a
//! Kelvin-wake heatmap reconstructed from the stored free-wave spectrum (the
//! same wave field the `michell wake` CLI draws, with the same colormap). The
//! scatter and GZ plots are hand-drawn with an `egui::Painter`; the wake is a
//! reconstructed elevation grid uploaded as a texture.
//!
//! It is used two ways: embedded in the editor (a sweep that writes `binary`
//! output pops this open on its result) and as the standalone `michell-viewer`
//! binary.

use egui::{Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, TextureHandle, Vec2};
use michell_cli::archive::{self, Row, SweepArchive};
use michell_cli::png;
use std::path::{Path, PathBuf};

const RAD: f64 = 180.0 / std::f64::consts::PI;
const PI: f64 = std::f64::consts::PI;

// Series palette (kept legible on both light and dark backgrounds).
const C_PRIMARY: Color32 = Color32::from_rgb(70, 150, 240);
const C_ACCENT: Color32 = Color32::from_rgb(90, 190, 110);
const C_WARN: Color32 = Color32::from_rgb(230, 150, 50);

/// A loaded archive plus the view's selection state.
pub struct Viewer {
    pub path: Option<PathBuf>,
    archive: Result<SweepArchive, String>,
    /// Index into `rows` of the selected row.
    selected: usize,
    /// Column indices (into the combined axis+metric column list) for the XY
    /// plot's horizontal and vertical axes.
    x_col: usize,
    y_col: usize,
    /// The reconstructed wake texture for the selected row, rebuilt lazily when
    /// the selection changes (the reconstruction is too heavy for every frame).
    wake: Option<WakeTex>,
}

impl Viewer {
    /// An empty viewer with a message (the standalone binary's initial state).
    pub fn empty(message: impl Into<String>) -> Self {
        Viewer {
            path: None,
            archive: Err(message.into()),
            selected: 0,
            x_col: 0,
            y_col: 0,
            wake: None,
        }
    }

    /// Decode raw `.msw` bytes.
    pub fn from_bytes(bytes: &[u8], path: Option<PathBuf>) -> Self {
        let archive = archive::read(bytes);
        let mut v = Viewer {
            path,
            archive,
            selected: 0,
            x_col: 0,
            y_col: 0,
            wake: None,
        };
        v.pick_default_columns();
        v
    }

    /// Read and decode a `.msw` file from disk.
    pub fn load(path: &Path) -> Self {
        match std::fs::read(path) {
            Ok(bytes) => Self::from_bytes(&bytes, Some(path.to_path_buf())),
            Err(e) => Viewer::empty(format!("cannot read {}: {e}", path.display())),
        }
    }

    /// Default the XY plot to a natural pairing: speed on X, total resistance
    /// (else the first metric) on Y.
    fn pick_default_columns(&mut self) {
        let Ok(a) = &self.archive else { return };
        let cols = column_names(&a.meta.axis_labels, &a.meta.metric_labels);
        let find = |want: &str| cols.iter().position(|c| c == want);
        self.x_col = find("speed").unwrap_or(0);
        self.y_col = find("rt")
            .or_else(|| find("rw"))
            .unwrap_or_else(|| a.meta.axis_labels.len().min(cols.len().saturating_sub(1)));
    }

    /// The decode/load error, if the archive failed to open.
    pub fn error(&self) -> Option<&str> {
        self.archive.as_ref().err().map(String::as_str)
    }

    /// Short one-line label for a window title / status.
    pub fn title(&self) -> String {
        let file = self
            .path
            .as_ref()
            .and_then(|p| p.file_name())
            .map(|n| n.to_string_lossy().into_owned());
        match (&self.archive, file) {
            (Ok(a), Some(f)) => match &a.meta.name {
                Some(n) => format!("{f} — {n}"),
                None => f,
            },
            (Ok(a), None) => a.meta.name.clone().unwrap_or_else(|| "sweep archive".into()),
            (Err(_), Some(f)) => f,
            (Err(_), None) => "sweep archive".into(),
        }
    }

    /// Render the whole viewer into `ui`.
    pub fn ui(&mut self, ui: &mut egui::Ui) {
        let a = match &self.archive {
            Ok(a) => a,
            Err(e) => {
                ui.colored_label(ui.visuals().error_fg_color, e);
                ui.label("Open a .msw sweep archive to view its results.");
                return;
            }
        };
        if self.selected >= a.rows.len() {
            self.selected = 0;
        }

        summary(ui, a);
        ui.separator();

        let cols = column_names(&a.meta.axis_labels, &a.meta.metric_labels);
        let n_axes = a.meta.axis_labels.len();

        // --- XY plot across all rows ---
        egui::CollapsingHeader::new("Plot")
            .default_open(true)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Y");
                    col_combo(ui, "y-col", &mut self.y_col, &cols);
                    ui.label("vs  X");
                    col_combo(ui, "x-col", &mut self.x_col, &cols);
                });
                let x_col = self.x_col.min(cols.len().saturating_sub(1));
                let y_col = self.y_col.min(cols.len().saturating_sub(1));
                let pts: Vec<[f64; 2]> = a
                    .rows
                    .iter()
                    .map(|r| [col_value(r, x_col, n_axes), col_value(r, y_col, n_axes)])
                    .collect();
                let highlight = a
                    .rows
                    .get(self.selected)
                    .map(|r| [col_value(r, x_col, n_axes), col_value(r, y_col, n_axes)]);
                draw_plot(
                    ui,
                    "xy",
                    220.0,
                    cols.get(x_col).map(String::as_str).unwrap_or("x"),
                    cols.get(y_col).map(String::as_str).unwrap_or("y"),
                    &[Series {
                        name: "",
                        color: C_PRIMARY,
                        pts,
                        connect: false,
                        markers: true,
                    }],
                    highlight,
                );
            });

        ui.separator();

        // --- rows table ---
        ui.label(format!("Rows ({})", a.rows.len()));
        rows_table(ui, a, &cols, n_axes, &mut self.selected);

        ui.separator();

        // --- selected-row detail: GZ curve + spectrum ---
        if let Some(row) = a.rows.get(self.selected) {
            ui.heading(format!("Row {} detail", self.selected + 1));
            row_params_line(ui, a, row);

            if row.gz_curve.is_empty() {
                ui.label(
                    egui::RichText::new(
                        "No GZ curve for this row (needs a loaded fleet: mass + vcg).",
                    )
                    .weak(),
                );
            } else {
                let pts: Vec<[f64; 2]> = row
                    .gz_curve
                    .iter()
                    .map(|&(heel, gz)| [heel * RAD, gz])
                    .collect();
                draw_plot(
                    ui,
                    "gz",
                    200.0,
                    "heel [deg]",
                    "GZ [m]",
                    &[Series {
                        name: "GZ",
                        color: C_ACCENT,
                        pts,
                        connect: true,
                        markers: false,
                    }],
                    None,
                );
            }

            if row.spectrum.is_empty() {
                ui.label(
                    egui::RichText::new("No free-wave spectrum for this row (dry fleet).").weak(),
                );
            } else {
                ui.label("Kelvin wake — surface elevation reconstructed from the stored spectrum");
                // Rebuild the texture when the selection changes.
                let stale = self.wake.as_ref().map(|w| w.row) != Some(self.selected);
                if stale {
                    self.wake = reconstruct_wake(row, a.meta.l_ref)
                        .map(|w| WakeTex::upload(ui.ctx(), self.selected, w));
                }
                if let Some(w) = &self.wake {
                    draw_wake(ui, w);
                } else {
                    ui.colored_label(
                        ui.visuals().error_fg_color,
                        "wake could not be reconstructed (degenerate spectrum)",
                    );
                }
            }
        }
    }
}

/// The study-metadata header block.
fn summary(ui: &mut egui::Ui, a: &SweepArchive) {
    let m = &a.meta;
    if let Some(name) = &m.name {
        ui.heading(name);
    }
    egui::Grid::new("summary")
        .num_columns(4)
        .spacing([18.0, 4.0])
        .show(ui, |ui| {
            ui.label("rows");
            ui.strong(a.rows.len().to_string());
            ui.label("fluid");
            ui.strong(&m.fluid);
            ui.end_row();

            ui.label("gravity");
            ui.strong(format!("{:.4} m/s²", m.gravity));
            ui.label("density");
            ui.strong(format!("{:.1} kg/m³", m.density));
            ui.end_row();

            ui.label("L_ref");
            ui.strong(format!("{:.3} m", m.l_ref));
            ui.label("mode");
            ui.strong(if m.float_mode { "equilibrium" } else { "fixed" });
            ui.end_row();

            ui.label("speeds");
            ui.strong(fmt_list(&m.speeds_ms, " m/s"));
            ui.label("spectrum pts");
            ui.strong(m.spectrum_points.to_string());
            ui.end_row();
        });
    if !m.axis_labels.is_empty() {
        ui.label(
            egui::RichText::new(format!("axes: {}", m.axis_labels.join(", ")))
                .weak()
                .small(),
        );
    }
    if !a.hull_files.is_empty() {
        let names: Vec<&str> = a.hull_files.iter().map(|(n, _)| n.as_str()).collect();
        ui.label(
            egui::RichText::new(format!("bundled: {}", names.join(", ")))
                .weak()
                .small(),
        );
    }
}

/// The scrollable, selectable rows table.
fn rows_table(
    ui: &mut egui::Ui,
    a: &SweepArchive,
    cols: &[String],
    n_axes: usize,
    selected: &mut usize,
) {
    egui::ScrollArea::both()
        .id_salt("rows-table")
        .max_height(220.0)
        .auto_shrink([false, false])
        .show(ui, |ui| {
            egui::Grid::new("rows-grid")
                .striped(true)
                .num_columns(cols.len() + 1)
                .spacing([12.0, 2.0])
                .show(ui, |ui| {
                    ui.strong("#");
                    for c in cols {
                        ui.strong(c.as_str());
                    }
                    ui.end_row();
                    for (i, row) in a.rows.iter().enumerate() {
                        if ui
                            .selectable_label(*selected == i, format!("{}", i + 1))
                            .clicked()
                        {
                            *selected = i;
                        }
                        for (col, _) in cols.iter().enumerate() {
                            ui.label(fmt_sig(col_value(row, col, n_axes)));
                        }
                        ui.end_row();
                    }
                });
        });
}

/// A one-line readout of the selected row's swept-parameter values.
fn row_params_line(ui: &mut egui::Ui, a: &SweepArchive, row: &Row) {
    if a.meta.axis_labels.is_empty() {
        return;
    }
    let parts: Vec<String> = a
        .meta
        .axis_labels
        .iter()
        .zip(&row.params)
        .map(|(k, v)| format!("{k} = {}", fmt_sig(*v)))
        .collect();
    ui.label(egui::RichText::new(parts.join("   ")).weak());
}

fn col_combo(ui: &mut egui::Ui, id: &str, sel: &mut usize, cols: &[String]) {
    let current = cols.get(*sel).map(String::as_str).unwrap_or("—");
    egui::ComboBox::from_id_salt(id)
        .selected_text(current)
        .show_ui(ui, |ui| {
            for (i, c) in cols.iter().enumerate() {
                ui.selectable_value(sel, i, c.as_str());
            }
        });
}

/// Combined column names: swept axes first, then metrics (matching row layout).
fn column_names(axis_labels: &[String], metric_labels: &[String]) -> Vec<String> {
    axis_labels
        .iter()
        .chain(metric_labels.iter())
        .cloned()
        .collect()
}

/// A row's value for combined-column `col`: axes index `params`, the rest index
/// `metrics`. Out-of-range (a malformed row) reads as NaN.
fn col_value(row: &Row, col: usize, n_axes: usize) -> f64 {
    if col < n_axes {
        row.params.get(col).copied().unwrap_or(f64::NAN)
    } else {
        row.metrics.get(col - n_axes).copied().unwrap_or(f64::NAN)
    }
}

// ---------------------------------------------------------------------------
// Hand-drawn plotting
// ---------------------------------------------------------------------------

struct Series<'a> {
    name: &'a str,
    color: Color32,
    pts: Vec<[f64; 2]>,
    /// Draw a polyline through the points (in given order).
    connect: bool,
    /// Draw a marker at each point.
    markers: bool,
}

/// Draw a labelled XY plot: a framed panel with gridlines, numeric ticks, axis
/// titles, one or more series, and an optional highlighted point. Bounds are
/// auto-fit to the finite data across all series.
fn draw_plot(
    ui: &mut egui::Ui,
    id: &str,
    height: f32,
    x_title: &str,
    y_title: &str,
    series: &[Series],
    highlight: Option<[f64; 2]>,
) {
    ui.push_id(id, |ui| {
        draw_plot_inner(ui, height, x_title, y_title, series, highlight)
    });
}

fn draw_plot_inner(
    ui: &mut egui::Ui,
    height: f32,
    x_title: &str,
    y_title: &str,
    series: &[Series],
    highlight: Option<[f64; 2]>,
) {
    let width = ui.available_width().max(220.0);
    let (resp, painter) = ui.allocate_painter(Vec2::new(width, height), Sense::hover());
    let outer = resp.rect;
    let bg = ui.visuals().extreme_bg_color;
    painter.rect_filled(outer, 3.0, bg);

    // Data bounds over finite points.
    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for s in series {
        for p in &s.pts {
            if p[0].is_finite() && p[1].is_finite() {
                xmin = xmin.min(p[0]);
                xmax = xmax.max(p[0]);
                ymin = ymin.min(p[1]);
                ymax = ymax.max(p[1]);
            }
        }
    }
    let weak = ui.visuals().weak_text_color();
    if !xmin.is_finite() || !ymin.is_finite() {
        painter.text(
            outer.center(),
            Align2::CENTER_CENTER,
            "no data",
            FontId::proportional(12.0),
            weak,
        );
        return;
    }
    // Pad a flat axis so a constant series still gets a visible band.
    if (xmax - xmin).abs() < 1e-12 {
        xmin -= 1.0;
        xmax += 1.0;
    }
    if (ymax - ymin).abs() < 1e-12 {
        let pad = ymin.abs().max(1.0) * 0.05;
        ymin -= pad;
        ymax += pad;
    } else {
        let pad = (ymax - ymin) * 0.06;
        ymin -= pad;
        ymax += pad;
    }

    // Plot area, leaving margins for tick labels and axis titles.
    let plot = Rect::from_min_max(
        Pos2::new(outer.left() + 52.0, outer.top() + 8.0),
        Pos2::new(outer.right() - 10.0, outer.bottom() - 30.0),
    );
    let x_of = |x: f64| plot.left() + ((x - xmin) / (xmax - xmin)) as f32 * plot.width();
    let y_of = |y: f64| plot.bottom() - ((y - ymin) / (ymax - ymin)) as f32 * plot.height();

    let grid = weak.gamma_multiply(0.25);
    let tick_font = FontId::proportional(9.5);
    // Gridlines + ticks.
    for k in 0..=4 {
        let fx = xmin + (xmax - xmin) * k as f64 / 4.0;
        let px = x_of(fx);
        painter.line_segment(
            [Pos2::new(px, plot.top()), Pos2::new(px, plot.bottom())],
            Stroke::new(1.0_f32, grid),
        );
        painter.text(
            Pos2::new(px, plot.bottom() + 3.0),
            Align2::CENTER_TOP,
            fmt_sig(fx),
            tick_font.clone(),
            weak,
        );
    }
    for k in 0..=4 {
        let fy = ymin + (ymax - ymin) * k as f64 / 4.0;
        let py = y_of(fy);
        painter.line_segment(
            [Pos2::new(plot.left(), py), Pos2::new(plot.right(), py)],
            Stroke::new(1.0_f32, grid),
        );
        painter.text(
            Pos2::new(plot.left() - 4.0, py),
            Align2::RIGHT_CENTER,
            fmt_sig(fy),
            tick_font.clone(),
            weak,
        );
    }
    // Frame.
    painter.rect_stroke(plot, 0.0, Stroke::new(1.0_f32, weak.gamma_multiply(0.6)));

    // A zero line, when the y-range straddles it.
    if ymin < 0.0 && ymax > 0.0 {
        let py = y_of(0.0);
        painter.line_segment(
            [Pos2::new(plot.left(), py), Pos2::new(plot.right(), py)],
            Stroke::new(1.0_f32, weak.gamma_multiply(0.5)),
        );
    }

    // Axis titles.
    painter.text(
        Pos2::new(plot.center().x, outer.bottom() - 2.0),
        Align2::CENTER_BOTTOM,
        x_title,
        FontId::proportional(11.0),
        weak,
    );
    painter.text(
        Pos2::new(outer.left() + 2.0, plot.top() - 2.0),
        Align2::LEFT_BOTTOM,
        y_title,
        FontId::proportional(11.0),
        weak,
    );

    // Series.
    for s in series {
        let screen: Vec<Pos2> = s
            .pts
            .iter()
            .filter(|p| p[0].is_finite() && p[1].is_finite())
            .map(|p| Pos2::new(x_of(p[0]), y_of(p[1])))
            .collect();
        if s.connect && screen.len() >= 2 {
            painter.add(egui::Shape::line(screen.clone(), Stroke::new(1.8_f32, s.color)));
        }
        if s.markers {
            for p in &screen {
                painter.circle_filled(*p, 2.5, s.color);
            }
        }
    }

    // Highlight the selected point.
    if let Some(h) = highlight {
        if h[0].is_finite() && h[1].is_finite() {
            let p = Pos2::new(x_of(h[0]), y_of(h[1]));
            painter.circle_stroke(p, 5.0, Stroke::new(2.0_f32, C_WARN));
        }
    }

    // Legend (only when named series would otherwise be ambiguous).
    let named: Vec<&Series> = series.iter().filter(|s| !s.name.is_empty()).collect();
    if named.len() > 1 {
        let mut y = plot.top() + 4.0;
        for s in named {
            painter.line_segment(
                [
                    Pos2::new(plot.right() - 60.0, y + 5.0),
                    Pos2::new(plot.right() - 46.0, y + 5.0),
                ],
                Stroke::new(2.0_f32, s.color),
            );
            painter.text(
                Pos2::new(plot.right() - 42.0, y),
                Align2::LEFT_TOP,
                s.name,
                tick_font.clone(),
                ui.visuals().text_color(),
            );
            y += 14.0;
        }
    }
}

// ---------------------------------------------------------------------------
// Kelvin-wake reconstruction
//
// The free-wave elevation is the same superposition the CLI's `wake` command
// integrates: ζ(x,y) = Re ∫ A(θ) e^{i(kx·x + ky·y)} dθ, with kx = ν secθ and
// ky = ν secθ tanθ (θ over the stored range, both signs). The stored samples
// are A(θ) at uniform θ already scaled to metres, so a trapezoidal sum
// reproduces the field. x is measured from the fleet centroid (the writer's
// phase reference cancels), so the source sits at x ≈ 0 and the wake trails
// toward −x.
// ---------------------------------------------------------------------------

/// A reconstructed elevation field over a rectangle, in the fleet frame.
struct Wake {
    nx: usize,
    ny: usize,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    /// Row-major `zeta[iy * nx + ix]` [m], iy from y0 to y1, ix from x0 to x1.
    zeta: Vec<f64>,
    /// Colour-scale saturation elevation (99.5th percentile of |ζ|) [m].
    zmax: f64,
    transverse_wavelength: f64,
    /// Number of spectrum samples that contributed (after the resolution cut).
    used_samples: usize,
    /// Longitudinal fade band: full colour for x ≤ `x_solid`, fully neutral for
    /// x ≥ `x_gone` — the free-wave field is only physical astern of the hull.
    x_solid: f64,
    x_gone: f64,
}

/// Quintic smoothstep on [0, 1] (zero slope at both ends).
fn smoother(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

/// Reconstruct the wake for one row. `None` when the spectrum is degenerate or
/// the wavenumber is unusable.
fn reconstruct_wake(row: &Row, l_ref: f64) -> Option<Wake> {
    let nu = row.wavenumber;
    if !(nu.is_finite() && nu > 0.0) || row.spectrum.len() < 3 {
        return None;
    }
    let twl = if row.transverse_wavelength > 0.0 {
        row.transverse_wavelength
    } else {
        2.0 * PI / nu
    };
    // Window: a few reference lengths of wake astern, a little ahead of the
    // source, and wide enough for the ±19.5° Kelvin wedge.
    let base = if l_ref.is_finite() && l_ref > 0.0 {
        l_ref
    } else {
        twl.max(1.0)
    };
    let (x0, x1) = (-3.2 * base, 0.6 * base);
    let yext = 1.5 * base;
    let (y0, y1) = (-yext, yext);

    // A generous internal grid so divergent (short-wave) arms resolve without
    // aliasing; the texture is scaled to the panel, so this sets quality, not
    // size. Recomputed once per selection, not per frame.
    let nx = 440usize;
    let hx = (x1 - x0) / (nx - 1) as f64;
    // Square pixels in metres: derive ny from the same spacing, then clamp.
    let ny = (((y1 - y0) / hx).round() as usize).clamp(60, 360);
    let hy = (y1 - y0) / (ny - 1) as f64;

    // Uniform θ spacing → trapezoidal weight. Bail if it is not increasing.
    let dtheta = row.spectrum[1].theta - row.spectrum[0].theta;
    if !(dtheta.is_finite() && dtheta > 0.0) {
        return None;
    }
    let theta_max = row.spectrum.last().unwrap().theta.abs().max(1e-6);

    // Pre-scale each sample: weight = dθ × resolution taper (drop components
    // whose wave the grid cannot resolve, cosine-tapering the marginal band,
    // as the CLI's elevation grid does), then fold into the amplitude.
    struct Comp {
        kx: f64,
        ky: f64,
        ar: f64,
        ai: f64,
    }
    let mut comps: Vec<Comp> = Vec::with_capacity(row.spectrum.len());
    for s in &row.spectrum {
        let c = s.theta.cos();
        if !(c.is_finite() && c.abs() > 1e-6) {
            continue;
        }
        let sec = 1.0 / c;
        let kx = nu * sec;
        let ky = nu * sec * s.theta.tan();
        // Drop components the grid cannot resolve (aliasing), fully keeping
        // only those with ≳4 px per wavelength and smoothly tapering the rest.
        let res = (PI / (kx.abs() * hx + 1e-12)).min(PI / (ky.abs() * hy + 1e-12));
        let taper = smoother((res - 1.5) / (2.8 - 1.5));
        // Also soften the very edge of the sampled range to curb ringing.
        let edge = {
            let f = s.theta.abs() / theta_max;
            if f > 0.8 {
                0.5 * (1.0 + (PI * (f - 0.8) / 0.2).cos())
            } else {
                1.0
            }
        };
        let w = dtheta * taper * edge;
        if w <= 0.0 {
            continue;
        }
        comps.push(Comp {
            kx,
            ky,
            ar: s.amp_re * w,
            ai: s.amp_im * w,
        });
    }
    if comps.is_empty() {
        return None;
    }

    // Accumulate ζ. For each component, march across x by complex recurrence
    // (one cis() per row, a complex multiply per pixel).
    let mut zeta = vec![0.0f64; nx * ny];
    for comp in &comps {
        let (sr, si) = (comp.kx * hx).sin_cos();
        let (step_im, step_re) = (sr, si); // e^{i kx hx}
        for iy in 0..ny {
            let y = y0 + iy as f64 * hy;
            let ph = comp.kx * x0 + comp.ky * y;
            let (sph, cph) = ph.sin_cos();
            let mut cur_re = comp.ar * cph - comp.ai * sph;
            let mut cur_im = comp.ar * sph + comp.ai * cph;
            let base = iy * nx;
            for ix in 0..nx {
                zeta[base + ix] += cur_re;
                let nr = cur_re * step_re - cur_im * step_im;
                let ni = cur_re * step_im + cur_im * step_re;
                cur_re = nr;
                cur_im = ni;
            }
        }
    }

    let zmax = percentile_abs(&zeta, 0.995).max(1e-9);
    Some(Wake {
        nx,
        ny,
        x0,
        x1,
        y0,
        y1,
        zeta,
        zmax,
        transverse_wavelength: twl,
        used_samples: comps.len(),
        x_solid: -0.5 * base,
        x_gone: 0.15 * base,
    })
}

/// The `p`-quantile (0..1) of `|values|`.
fn percentile_abs(values: &[f64], p: f64) -> f64 {
    let mut abs: Vec<f64> = values.iter().map(|v| v.abs()).collect();
    if abs.is_empty() {
        return 0.0;
    }
    abs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let idx = ((abs.len() - 1) as f64 * p.clamp(0.0, 1.0)).round() as usize;
    abs[idx]
}

/// A wake elevation field uploaded to the GPU, with the scalars needed to
/// label it. Rebuilt only when the selected row changes.
struct WakeTex {
    row: usize,
    tex: TextureHandle,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    zmax: f64,
    transverse_wavelength: f64,
    used_samples: usize,
}

impl WakeTex {
    fn upload(ctx: &egui::Context, row: usize, w: Wake) -> WakeTex {
        // Row 0 of the image is the top; put +y at the top, x0 at the left.
        let hx = (w.x1 - w.x0) / (w.nx - 1) as f64;
        let neutral = png::diverging(0.0); // undisturbed-water colour
        let mut pixels = Vec::with_capacity(w.nx * w.ny);
        for r in 0..w.ny {
            let iy = w.ny - 1 - r;
            for ix in 0..w.nx {
                let t = (w.zeta[iy * w.nx + ix] / w.zmax).clamp(-1.0, 1.0);
                // Fade toward neutral ahead of the hull — free waves only trail
                // astern, so the forward field is not physical.
                let x = w.x0 + ix as f64 * hx;
                let vis = 1.0 - smoother((x - w.x_solid) / (w.x_gone - w.x_solid));
                let [cr, cg, cb] = png::fade(png::diverging(t), neutral, 1.0 - vis);
                pixels.push(Color32::from_rgb(cr, cg, cb));
            }
        }
        let image = egui::ColorImage {
            size: [w.nx, w.ny],
            pixels,
        };
        let tex = ctx.load_texture(
            format!("wake-{row}"),
            image,
            egui::TextureOptions::LINEAR,
        );
        WakeTex {
            row,
            tex,
            x0: w.x0,
            x1: w.x1,
            y0: w.y0,
            y1: w.y1,
            zmax: w.zmax,
            transverse_wavelength: w.transverse_wavelength,
            used_samples: w.used_samples,
        }
    }
}

/// Draw the wake texture with axes, a Kelvin-wedge overlay, the source marker,
/// and a colour bar.
fn draw_wake(ui: &mut egui::Ui, w: &WakeTex) {
    let aspect = ((w.x1 - w.x0) / (w.y1 - w.y0)).abs().max(0.1) as f32;
    let (ml, mr, mt, mb) = (46.0f32, 64.0f32, 8.0f32, 24.0f32);
    let avail = ui.available_width().max(280.0);
    let mut img_w = avail - ml - mr;
    let mut img_h = img_w / aspect;
    let max_h = 360.0;
    if img_h > max_h {
        img_h = max_h;
        img_w = img_h * aspect;
    }
    let total = Vec2::new(img_w + ml + mr, img_h + mt + mb);
    let (resp, painter) = ui.allocate_painter(total, Sense::hover());
    let rect = resp.rect;
    let img_rect = Rect::from_min_size(
        Pos2::new(rect.left() + ml, rect.top() + mt),
        Vec2::new(img_w, img_h),
    );
    painter.image(
        w.tex.id(),
        img_rect,
        Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
        Color32::WHITE,
    );
    let weak = ui.visuals().weak_text_color();
    painter.rect_stroke(img_rect, 0.0, Stroke::new(1.0_f32, weak.gamma_multiply(0.6)));

    // Screen mappings: x0..x1 → left..right; y1..y0 → top..bottom (+y up).
    let x_of = |x: f64| img_rect.left() + ((x - w.x0) / (w.x1 - w.x0)) as f32 * img_rect.width();
    let y_of = |y: f64| img_rect.top() + ((w.y1 - y) / (w.y1 - w.y0)) as f32 * img_rect.height();
    let tick_font = FontId::proportional(9.5);

    // Kelvin wedge (±19.47°) from the source, trailing aft (−x).
    if w.x0 < 0.0 {
        let half = 19.471_f64.to_radians().tan();
        let wedge = Color32::from_rgba_unmultiplied(255, 255, 255, 120);
        for sign in [-1.0f64, 1.0] {
            let x_end = w.x0;
            let y_end = (sign * (-x_end) * half).clamp(w.y0, w.y1);
            painter.line_segment(
                [
                    Pos2::new(x_of(0.0), y_of(0.0)),
                    Pos2::new(x_of(x_end), y_of(y_end)),
                ],
                Stroke::new(1.0_f32, wedge),
            );
        }
    }
    // Source marker at (0, 0).
    if (w.x0..=w.x1).contains(&0.0) && (w.y0..=w.y1).contains(&0.0) {
        painter.circle_stroke(
            Pos2::new(x_of(0.0), y_of(0.0)),
            3.0,
            Stroke::new(1.5_f32, Color32::from_rgb(30, 30, 30)),
        );
    }

    // Axis ticks.
    for k in 0..=4 {
        let fx = w.x0 + (w.x1 - w.x0) * k as f64 / 4.0;
        painter.text(
            Pos2::new(x_of(fx), img_rect.bottom() + 2.0),
            Align2::CENTER_TOP,
            fmt_sig(fx),
            tick_font.clone(),
            weak,
        );
    }
    for k in 0..=4 {
        let fy = w.y0 + (w.y1 - w.y0) * k as f64 / 4.0;
        painter.text(
            Pos2::new(img_rect.left() - 4.0, y_of(fy)),
            Align2::RIGHT_CENTER,
            fmt_sig(fy),
            tick_font.clone(),
            weak,
        );
    }
    painter.text(
        Pos2::new(img_rect.center().x, rect.bottom() - 1.0),
        Align2::CENTER_BOTTOM,
        "x  (aft ← 0 → fwd) [m]",
        FontId::proportional(10.5),
        weak,
    );
    painter.text(
        Pos2::new(rect.left() + 1.0, img_rect.top() - 1.0),
        Align2::LEFT_BOTTOM,
        "y [m]",
        FontId::proportional(10.5),
        weak,
    );

    // Colour bar on the right: +zmax (crest) at top → −zmax (trough) at bottom.
    let bar = Rect::from_min_size(
        Pos2::new(img_rect.right() + 12.0, img_rect.top()),
        Vec2::new(12.0, img_rect.height()),
    );
    let bands = 48;
    for b in 0..bands {
        let t = 1.0 - 2.0 * b as f64 / (bands - 1) as f64; // +1 at top
        let [cr, cg, cb] = png::diverging(t);
        let y = bar.top() + b as f32 / bands as f32 * bar.height();
        let h = bar.height() / bands as f32 + 1.0;
        painter.rect_filled(
            Rect::from_min_size(Pos2::new(bar.left(), y), Vec2::new(bar.width(), h)),
            0.0,
            Color32::from_rgb(cr, cg, cb),
        );
    }
    painter.rect_stroke(bar, 0.0, Stroke::new(1.0_f32, weak.gamma_multiply(0.6)));
    painter.text(
        Pos2::new(bar.right() + 3.0, bar.top()),
        Align2::LEFT_TOP,
        format!("+{} m", fmt_sig(w.zmax)),
        tick_font.clone(),
        weak,
    );
    painter.text(
        Pos2::new(bar.right() + 3.0, bar.center().y),
        Align2::LEFT_CENTER,
        "0",
        tick_font.clone(),
        weak,
    );
    painter.text(
        Pos2::new(bar.right() + 3.0, bar.bottom()),
        Align2::LEFT_BOTTOM,
        format!("−{} m", fmt_sig(w.zmax)),
        tick_font.clone(),
        weak,
    );

    ui.label(
        egui::RichText::new(format!(
            "λ_transverse = {:.2} m · {} spectral components · colour saturates at the 99.5th \
             percentile of |ζ| (blue trough → red crest)",
            w.transverse_wavelength, w.used_samples
        ))
        .weak()
        .small(),
    );
}

// ---------------------------------------------------------------------------
// Number formatting
// ---------------------------------------------------------------------------

/// Compact ~4-significant-figure rendering: fixed notation in a sane range,
/// scientific for very large/small magnitudes, and trailing zeros trimmed.
fn fmt_sig(v: f64) -> String {
    if !v.is_finite() {
        return "—".into();
    }
    if v == 0.0 {
        return "0".into();
    }
    let a = v.abs();
    let s = if !(1e-3..1e6).contains(&a) {
        format!("{v:.3e}")
    } else {
        let decimals = if a >= 100.0 {
            1
        } else if a >= 1.0 {
            3
        } else {
            5
        };
        let mut s = format!("{v:.*}", decimals);
        if s.contains('.') {
            while s.ends_with('0') {
                s.pop();
            }
            if s.ends_with('.') {
                s.pop();
            }
        }
        s
    };
    s
}

fn fmt_list(xs: &[f64], suffix: &str) -> String {
    if xs.is_empty() {
        return "—".into();
    }
    let shown: Vec<String> = xs.iter().take(6).map(|v| fmt_sig(*v)).collect();
    let more = if xs.len() > 6 {
        format!(" … (+{})", xs.len() - 6)
    } else {
        String::new()
    };
    format!("{}{suffix}{more}", shown.join(", "))
}

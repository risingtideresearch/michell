//! An egui view of a `.msw` sweep archive. Decodes the self-contained binary
//! format (via `michell_cli::archive`) and renders it: the study metadata, an
//! XY plot of any column against any other across all rows, a selectable rows
//! table, and — for the selected row — the righting-arm (GZ) curve and the
//! free-wave spectrum. Plots are hand-drawn with an `egui::Painter`, matching
//! the import-preview style and keeping the crate dependency-light.
//!
//! It is used two ways: embedded in the editor (a sweep that writes `binary`
//! output pops this open on its result) and as the standalone `michell-viewer`
//! binary.

use egui::{Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};
use michell_cli::archive::{self, Row, SweepArchive};
use std::path::{Path, PathBuf};

const RAD: f64 = 180.0 / std::f64::consts::PI;

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
                // Connect the line only when X is strictly increasing (a clean
                // 1-D sweep); otherwise a multi-axis grid would draw zig-zags,
                // so show markers alone.
                let monotone = pts.windows(2).all(|w| w[1][0] > w[0][0]);
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
                        connect: monotone,
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
                let amp: Vec<[f64; 2]> = row
                    .spectrum
                    .iter()
                    .map(|s| [s.theta * RAD, (s.amp_re * s.amp_re + s.amp_im * s.amp_im).sqrt()])
                    .collect();
                draw_plot(
                    ui,
                    "spec-amp",
                    180.0,
                    "theta [deg]",
                    "|A| [m/rad]",
                    &[Series {
                        name: "|A|",
                        color: C_PRIMARY,
                        pts: amp,
                        connect: true,
                        markers: false,
                    }],
                    None,
                );
                let drw: Vec<[f64; 2]> = row
                    .spectrum
                    .iter()
                    .map(|s| [s.theta * RAD, s.drw_dtheta])
                    .collect();
                draw_plot(
                    ui,
                    "spec-drw",
                    180.0,
                    "theta [deg]",
                    "dRw/dtheta [N/rad]",
                    &[Series {
                        name: "dRw/dθ",
                        color: C_WARN,
                        pts: drw,
                        connect: true,
                        markers: false,
                    }],
                    None,
                );
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

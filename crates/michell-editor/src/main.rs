//! A small egui app for authoring `michell` sweep-study manifests. It edits a
//! typed model with a form, shows a live JSON preview and validation, reads and
//! writes the `.json` files the CLI consumes, and can drive the CLI itself:
//! import CAD/mesh geometry (lofting it to full-band `.hull` bodies) and run the
//! sweep, both with live progress.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use michell_editor::model::{HullSpec, Manifest};
use michell_editor::preview::Preview;
use michell_editor::runner::{self, Job, JobKind};
use michell_editor::validate::Level;
use michell_editor::{jsonio, ui, validate};
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::time::Duration;

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1000.0, 760.0])
            .with_min_inner_size([640.0, 480.0])
            .with_title("michell — sweep manifest editor"),
        ..Default::default()
    };
    eframe::run_native(
        "michell-editor",
        options,
        Box::new(|_cc| Ok(Box::new(EditorApp::default()))),
    )
}

/// Units picker for STL import (STL files carry no units).
#[derive(Clone, Copy, PartialEq, Eq)]
enum StlUnits {
    Mm,
    Cm,
    M,
    In,
    Ft,
}

impl StlUnits {
    const ALL: [StlUnits; 5] = [
        StlUnits::Mm,
        StlUnits::Cm,
        StlUnits::M,
        StlUnits::In,
        StlUnits::Ft,
    ];
    fn as_str(self) -> &'static str {
        match self {
            StlUnits::Mm => "mm",
            StlUnits::Cm => "cm",
            StlUnits::M => "m",
            StlUnits::In => "in",
            StlUnits::Ft => "ft",
        }
    }

    /// Metres per file unit.
    fn scale(self) -> f32 {
        match self {
            StlUnits::Mm => 0.001,
            StlUnits::Cm => 0.01,
            StlUnits::M => 1.0,
            StlUnits::In => 0.0254,
            StlUnits::Ft => 0.3048,
        }
    }
}

/// Parameters for an in-progress "import & loft" the user is filling in.
struct ImportDialog {
    source: PathBuf,
    is_stl: bool,
    waterline: String,
    band: String,
    units: StlUnits,
    prefix: String,
    error: Option<String>,
    /// Transverse silhouette of the source, for placing the waterline/band.
    preview: Result<Preview, String>,
    /// Which handle is being dragged in the preview (0 = DWL, 1 = band top).
    dragging: Option<u8>,
}

struct EditorApp {
    manifest: Manifest,
    path: Option<PathBuf>,
    dirty: bool,
    status: String,
    preview: Result<String, String>,
    show_preview: bool,
    michell_bin: PathBuf,
    /// A running loft/sweep subprocess, if any (one at a time).
    job: Option<Job>,
    import: Option<ImportDialog>,
}

impl Default for EditorApp {
    fn default() -> Self {
        let manifest = Manifest::default();
        let preview = jsonio::to_string(&manifest);
        EditorApp {
            manifest,
            path: None,
            dirty: false,
            status: "new manifest".into(),
            preview,
            show_preview: true,
            michell_bin: runner::default_michell_bin(),
            job: None,
            import: None,
        }
    }
}

impl EditorApp {
    fn refresh_preview(&mut self) {
        self.preview = jsonio::to_string(&self.manifest);
    }

    fn manifest_dir(&self) -> Option<PathBuf> {
        self.path
            .as_ref()
            .and_then(|p| p.parent())
            .map(Path::to_path_buf)
    }

    fn title_file(&self) -> String {
        self.path
            .as_ref()
            .and_then(|p| p.file_name())
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "untitled".into())
    }

    fn error_count(&self) -> usize {
        validate::validate(&self.manifest)
            .iter()
            .filter(|i| i.level == Level::Error)
            .count()
    }

    fn busy(&self) -> bool {
        self.job.is_some()
    }

    // --- file ops ---

    fn new_manifest(&mut self) {
        self.manifest = Manifest::default();
        self.path = None;
        self.dirty = false;
        self.status = "new manifest".into();
        self.refresh_preview();
    }

    fn open(&mut self) {
        let Some(path) = rfd::FileDialog::new()
            .add_filter("manifest", &["json"])
            .pick_file()
        else {
            return;
        };
        match std::fs::read_to_string(&path)
            .map_err(|e| e.to_string())
            .and_then(|t| jsonio::from_str(&t))
        {
            Ok(m) => {
                self.manifest = m;
                self.status = format!("opened {}", path.display());
                self.path = Some(path);
                self.dirty = false;
                self.refresh_preview();
            }
            Err(e) => self.status = format!("open failed: {e}"),
        }
    }

    fn save(&mut self) {
        match self.path.clone() {
            Some(p) => self.save_to(p),
            None => self.save_as(),
        }
    }

    fn save_as(&mut self) {
        let mut dialog = rfd::FileDialog::new().add_filter("manifest", &["json"]);
        if let Some(dir) = self.manifest_dir() {
            dialog = dialog.set_directory(dir);
        }
        dialog = dialog.set_file_name(self.suggested_name());
        if let Some(path) = dialog.save_file() {
            self.save_to(path);
        }
    }

    fn suggested_name(&self) -> String {
        if self.manifest.name.trim().is_empty() {
            return "study.json".into();
        }
        let slug: String = self
            .manifest
            .name
            .trim()
            .chars()
            .map(|c| {
                if c.is_alphanumeric() {
                    c.to_ascii_lowercase()
                } else {
                    '-'
                }
            })
            .collect();
        format!("{}.json", slug.trim_matches('-'))
    }

    fn save_to(&mut self, path: PathBuf) {
        match jsonio::to_string(&self.manifest) {
            Ok(text) => match std::fs::write(&path, text) {
                Ok(()) => {
                    self.status = format!("saved {}", path.display());
                    self.path = Some(path);
                    self.dirty = false;
                }
                Err(e) => self.status = format!("save failed: {e}"),
            },
            Err(e) => self.status = format!("cannot save — {e}"),
        }
    }

    // --- running the CLI ---

    /// Prompt for a CAD/mesh file, then open the loft-parameters dialog.
    fn begin_import(&mut self) {
        if self.busy() {
            return;
        }
        let Some(src) = rfd::FileDialog::new()
            .add_filter("CAD / mesh", &["igs", "iges", "stl"])
            .pick_file()
        else {
            return;
        };
        let is_stl = src
            .extension()
            .map(|e| e.eq_ignore_ascii_case("stl"))
            .unwrap_or(false);
        let prefix = src
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "hull".into());
        let preview = Preview::load(&src, is_stl);
        self.import = Some(ImportDialog {
            source: src,
            is_stl,
            waterline: "0".into(),
            band: String::new(),
            units: StlUnits::Mm,
            prefix,
            error: None,
            preview,
            dragging: None,
        });
    }

    /// Validate the import dialog and spawn `michell loft`.
    fn start_loft(&mut self) {
        let Some(d) = &mut self.import else { return };
        if d.waterline.trim().parse::<f64>().is_err() {
            d.error = Some("design waterline must be a number".into());
            return;
        }
        if !d.band.trim().is_empty() && d.band.trim().parse::<f64>().is_err() {
            d.error = Some("band must be a number (or blank)".into());
            return;
        }
        if d.prefix.trim().is_empty() {
            d.error = Some("output prefix is required".into());
            return;
        }

        let d = self.import.take().expect("checked above");
        // Lofted bodies land next to the manifest (so their paths stay relative)
        // when it is saved; otherwise beside the source file.
        let cwd = self
            .manifest_dir()
            .or_else(|| d.source.parent().map(Path::to_path_buf))
            .unwrap_or_else(|| PathBuf::from("."));

        let mut args: Vec<OsString> = vec![
            "loft".into(),
            d.source.clone().into_os_string(),
            "--waterline".into(),
            d.waterline.trim().into(),
        ];
        if !d.band.trim().is_empty() {
            args.push("--band".into());
            args.push(d.band.trim().into());
        }
        if d.is_stl {
            args.push("--units".into());
            args.push(d.units.as_str().into());
        }
        args.push("-o".into());
        args.push(d.prefix.trim().into());

        self.status = format!("lofting {}…", d.source.display());
        self.job = Some(runner::spawn(
            &self.michell_bin,
            args,
            &cwd,
            JobKind::Loft,
            "Lofting geometry".into(),
        ));
    }

    /// Save (if needed) and spawn `michell sweep`.
    fn run_sweep(&mut self) {
        if self.busy() {
            return;
        }
        if self.error_count() > 0 {
            self.status = "fix the validation errors before running".into();
            return;
        }
        if self.dirty || self.path.is_none() {
            self.save();
        }
        let Some(path) = self.path.clone() else {
            self.status = "save the manifest before running".into();
            return;
        };
        if self.dirty {
            return; // save failed or was cancelled
        }
        let cwd = self.manifest_dir().unwrap_or_else(|| PathBuf::from("."));
        let args: Vec<OsString> = vec!["sweep".into(), path.into_os_string()];
        self.status = "running sweep…".into();
        self.job = Some(runner::spawn(
            &self.michell_bin,
            args,
            &cwd,
            JobKind::Sweep,
            "Running sweep".into(),
        ));
    }

    /// Handle a finished job: add lofted hulls, or report the sweep result.
    fn finish_job(&mut self, job: Job) {
        let outcome = job.outcome.expect("finished job has an outcome");
        match (job.kind, outcome) {
            (JobKind::Loft, Ok(success)) => {
                let produced = parse_lofted(&success.stdout);
                if produced.is_empty() {
                    self.status = "loft finished but produced no .hull files".into();
                    return;
                }
                let dir = self.manifest_dir();
                let mut ids: Vec<String> =
                    self.manifest.hulls.iter().map(|h| h.id.clone()).collect();
                for name in &produced {
                    let abs = success.cwd.join(name);
                    let file = ui::rel_to(dir.as_deref(), &abs);
                    let base = Path::new(name)
                        .file_stem()
                        .map(|s| s.to_string_lossy().into_owned())
                        .unwrap_or_else(|| name.clone());
                    let id = unique_id(&base, &mut ids);
                    self.manifest.hulls.push(HullSpec {
                        id,
                        file,
                        ..Default::default()
                    });
                }
                self.dirty = true;
                self.refresh_preview();
                self.status = format!("lofted and added {} hull(s)", produced.len());
            }
            (JobKind::Sweep, Ok(success)) => {
                self.status = match &self.manifest.output.file {
                    f if f.trim().is_empty() => {
                        let rows = success.stdout.lines().count().saturating_sub(1);
                        format!("sweep done — {rows} row(s) printed to stdout (no output file set)")
                    }
                    f => {
                        let where_ = self
                            .manifest_dir()
                            .map(|d| d.join(f))
                            .unwrap_or_else(|| PathBuf::from(f));
                        format!("sweep done — wrote {}", where_.display())
                    }
                };
            }
            (JobKind::Loft, Err(e)) => self.status = format!("loft failed: {e}"),
            (JobKind::Sweep, Err(e)) => self.status = format!("sweep failed: {e}"),
        }
    }
}

/// Pull the produced `.hull` filenames out of `michell loft`'s stdout table
/// (the first column of each non-header row is the file).
fn parse_lofted(stdout: &str) -> Vec<String> {
    stdout
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter(|t| t.ends_with(".hull"))
        .map(str::to_string)
        .collect()
}

fn unique_id(base: &str, existing: &mut Vec<String>) -> String {
    let base = if base.is_empty() { "hull" } else { base };
    let mut id = base.to_string();
    let mut n = 2;
    while existing.iter().any(|e| e == &id) {
        id = format!("{base}_{n}");
        n += 1;
    }
    existing.push(id.clone());
    id
}

impl eframe::App for EditorApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Poll a running job; process it when it finishes.
        if let Some(job) = self.job.as_mut() {
            job.poll();
            if job.is_running() {
                ctx.request_repaint_after(Duration::from_millis(80));
            } else {
                let job = self.job.take().expect("just checked");
                self.finish_job(job);
            }
        }

        // Keyboard shortcuts.
        ctx.input_mut(|i| {
            use egui::{Key, Modifiers};
            if i.consume_key(Modifiers::COMMAND, Key::S) {
                self.save();
            }
            if i.consume_key(Modifiers::COMMAND, Key::O) {
                self.open();
            }
            if i.consume_key(Modifiers::COMMAND, Key::N) {
                self.new_manifest();
            }
        });

        let mut want_import = false;
        let mut want_run = false;
        let mut want_locate = false;

        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("New").clicked() {
                        self.new_manifest();
                        ui.close_menu();
                    }
                    if ui.button("Open…").clicked() {
                        self.open();
                        ui.close_menu();
                    }
                    if ui.button("Save").clicked() {
                        self.save();
                        ui.close_menu();
                    }
                    if ui.button("Save As…").clicked() {
                        self.save_as();
                        ui.close_menu();
                    }
                });
                ui.menu_button("View", |ui| {
                    ui.checkbox(&mut self.show_preview, "JSON preview");
                });
                ui.menu_button("Tools", |ui| {
                    if ui.button("Locate michell binary…").clicked() {
                        want_locate = true;
                        ui.close_menu();
                    }
                    ui.label(
                        egui::RichText::new(format!("using: {}", self.michell_bin.display()))
                            .weak()
                            .small(),
                    );
                });

                ui.separator();
                ui.add_enabled_ui(!self.busy(), |ui| {
                    if ui
                        .button("▶ Run sweep")
                        .on_hover_text("Save and run `michell sweep` on this manifest")
                        .clicked()
                    {
                        want_run = true;
                    }
                });

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let title = format!(
                        "{}{}",
                        self.title_file(),
                        if self.dirty { " *" } else { "" }
                    );
                    ui.label(title);
                });
            });
        });

        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            let issues = validate::validate(&self.manifest);
            ui.horizontal(|ui| {
                let errors = issues.iter().filter(|i| i.level == Level::Error).count();
                if errors == 0 {
                    ui.colored_label(egui::Color32::from_rgb(80, 170, 90), "✓ valid");
                } else {
                    ui.colored_label(ui.visuals().error_fg_color, format!("✗ {errors} error(s)"));
                }
                ui.separator();
                ui.label(&self.status);
            });
            if !issues.is_empty() {
                egui::ScrollArea::vertical()
                    .id_salt("issues")
                    .max_height(96.0)
                    .show(ui, |ui| {
                        for issue in &issues {
                            let color = match issue.level {
                                Level::Error => ui.visuals().error_fg_color,
                                Level::Warning => egui::Color32::from_rgb(200, 150, 40),
                            };
                            ui.colored_label(color, format!("• {}", issue.msg));
                        }
                    });
            }
        });

        if self.show_preview {
            egui::SidePanel::right("preview")
                .default_width(340.0)
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        ui.heading("JSON");
                        if ui.button("copy").clicked() {
                            if let Ok(text) = &self.preview {
                                ui.ctx().copy_text(text.clone());
                            }
                        }
                    });
                    ui.separator();
                    egui::ScrollArea::both().show(ui, |ui| match &self.preview {
                        Ok(text) => {
                            ui.add(
                                egui::TextEdit::multiline(&mut text.as_str())
                                    .code_editor()
                                    .desired_width(f32::INFINITY),
                            );
                        }
                        Err(e) => {
                            ui.colored_label(ui.visuals().error_fg_color, e);
                        }
                    });
                });
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                let dir = self.manifest_dir();
                let resp = ui::manifest_form(ui, &mut self.manifest, dir.as_deref());
                if resp.changed {
                    self.dirty = true;
                    self.refresh_preview();
                }
                if resp.import_clicked {
                    want_import = true;
                }
            });
        });

        self.import_window(ctx);
        self.job_window(ctx);

        // Act on intents gathered above, now that the panel borrows are gone.
        if want_locate {
            if let Some(p) = rfd::FileDialog::new().pick_file() {
                self.michell_bin = p;
                self.status = format!("michell binary: {}", self.michell_bin.display());
            }
        }
        if want_import {
            self.begin_import();
        }
        if want_run {
            self.run_sweep();
        }
    }
}

impl EditorApp {
    fn import_window(&mut self, ctx: &egui::Context) {
        let mut do_loft = false;
        let mut cancel = false;
        if let Some(d) = &mut self.import {
            egui::Window::new("Import & loft geometry")
                .collapsible(false)
                .resizable(false)
                .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
                .show(ctx, |ui| {
                    ui.label(format!("source: {}", d.source.display()));
                    ui.add_space(6.0);
                    egui::Grid::new("import-grid")
                        .num_columns(2)
                        .spacing([12.0, 6.0])
                        .show(ui, |ui| {
                            ui.label("design waterline z");
                            ui.text_edit_singleline(&mut d.waterline);
                            ui.end_row();
                            ui.label("band [m]");
                            ui.horizontal(|ui| {
                                ui.text_edit_singleline(&mut d.band);
                                ui.label(egui::RichText::new("(blank = draft/2)").weak());
                            });
                            ui.end_row();
                            if d.is_stl {
                                ui.label("units");
                                egui::ComboBox::from_id_salt("stl-units")
                                    .selected_text(d.units.as_str())
                                    .show_ui(ui, |ui| {
                                        for u in StlUnits::ALL {
                                            ui.selectable_value(&mut d.units, u, u.as_str());
                                        }
                                    });
                                ui.end_row();
                            }
                            ui.label("output prefix");
                            ui.text_edit_singleline(&mut d.prefix);
                            ui.end_row();
                        });
                    ui.add_space(8.0);
                    ui.label(
                        egui::RichText::new(
                            "forward view — drag the lines to set the waterline and band",
                        )
                        .weak()
                        .small(),
                    );
                    draw_geometry_view(ui, d);
                    if let Some(e) = &d.error {
                        ui.colored_label(ui.visuals().error_fg_color, e);
                    }
                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        if ui.button("Loft").clicked() {
                            do_loft = true;
                        }
                        if ui.button("Cancel").clicked() {
                            cancel = true;
                        }
                    });
                });
        }
        if cancel {
            self.import = None;
        }
        if do_loft {
            self.start_loft();
        }
    }

    fn job_window(&mut self, ctx: &egui::Context) {
        let Some(job) = &self.job else { return };
        egui::Window::new(&job.title)
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                ui.set_min_width(420.0);
                match job.fraction() {
                    Some(f) => {
                        ui.add(egui::ProgressBar::new(f).show_percentage());
                    }
                    None => {
                        ui.add(egui::ProgressBar::new(0.0).animate(true));
                        ui.label("working…");
                    }
                }
                ui.add_space(6.0);
                egui::ScrollArea::vertical()
                    .id_salt("job-log")
                    .max_height(180.0)
                    .stick_to_bottom(true)
                    .show(ui, |ui| {
                        for line in &job.log {
                            ui.monospace(line);
                        }
                    });
            });
    }
}

/// A forward (body-plan) view of the import source: the geometry's transverse
/// silhouette with draggable design-waterline and band lines. Dragging writes
/// back into the dialog's `waterline`/`band` fields (and typing moves the
/// lines, since they are read from those fields each frame).
fn draw_geometry_view(ui: &mut egui::Ui, d: &mut ImportDialog) {
    use egui::{Align2, Color32, FontId, Pos2, Sense, Stroke, Vec2};

    let preview = match &d.preview {
        Ok(p) => p,
        Err(e) => {
            ui.weak(format!("(no preview: {e})"));
            return;
        }
    };

    let scale = if d.is_stl { d.units.scale() } else { 1.0 };
    let dwl = d.waterline.trim().parse::<f32>().unwrap_or(0.0);
    let keel = preview.z_min * scale;
    let band_blank = d.band.trim().is_empty();
    // Blank band → the loft's default of half the design draft; show it so the
    // user sees (and can grab) the effective band.
    let band = if band_blank {
        (0.5 * (dwl - keel)).max(0.0)
    } else {
        d.band.trim().parse::<f32>().unwrap_or(0.0).max(0.0)
    };
    let band_top = dwl + band;

    // Fixed view bounds from the geometry alone — never the dragged lines — so
    // dragging the waterline/band cannot rescale or compress the view.
    let ymid = 0.5 * (preview.y_min + preview.y_max) * scale;
    let zmid = 0.5 * (preview.z_min + preview.z_max) * scale;
    let yspan = ((preview.y_max - preview.y_min) * scale).max(1e-6);
    let zspan = ((preview.z_max - preview.z_min) * scale).max(1e-6);

    let (resp, painter) = ui.allocate_painter(Vec2::new(400.0, 240.0), Sense::click_and_drag());
    let rect = resp.rect;
    painter.rect_filled(rect, 2.0, ui.visuals().extreme_bg_color);

    // One pixels-per-metre scale for both axes (true proportions), with an 8%
    // margin, centred in the rect. y → x, z → screen y (inverted).
    let s_px = (rect.width() / yspan).min(rect.height() / zspan) * 0.92;
    let (cx, cy) = (rect.center().x, rect.center().y);
    let x_of = |y: f32| cx + (y - ymid) * s_px;
    let y_of = |z: f32| cy - (z - zmid) * s_px;
    let z_of = |py: f32| zmid - (py - cy) / s_px;
    // Keep a dragged line within the visible height.
    let z_half = rect.height() * 0.5 / s_px;
    let (z_lo, z_hi) = (zmid - z_half, zmid + z_half);

    let weak = ui.visuals().weak_text_color();
    // Midship section outline: the raw model's segments where the plane x=x_mid
    // cut it (all hulls at that station, unclustered).
    let hull_col = Color32::from_rgb(120, 165, 205);
    for seg in &preview.segments {
        painter.line_segment(
            [
                Pos2::new(x_of(seg[0][0] * scale), y_of(seg[0][1] * scale)),
                Pos2::new(x_of(seg[1][0] * scale), y_of(seg[1][1] * scale)),
            ],
            Stroke::new(1.3_f32, hull_col),
        );
    }
    let x0 = x_of(0.0);
    if (rect.left()..=rect.right()).contains(&x0) {
        painter.line_segment(
            [Pos2::new(x0, rect.top()), Pos2::new(x0, rect.bottom())],
            Stroke::new(1.0_f32, weak.gamma_multiply(0.4)),
        );
    }

    let dwl_col = Color32::from_rgb(70, 150, 240);
    let band_col = Color32::from_rgb(90, 190, 110);
    painter.line_segment(
        [
            Pos2::new(rect.left(), y_of(dwl)),
            Pos2::new(rect.right(), y_of(dwl)),
        ],
        Stroke::new(2.0_f32, dwl_col),
    );
    let band_w: f32 = if band_blank { 1.0 } else { 2.0 };
    let band_a: f32 = if band_blank { 0.6 } else { 1.0 };
    painter.line_segment(
        [
            Pos2::new(rect.left(), y_of(band_top)),
            Pos2::new(rect.right(), y_of(band_top)),
        ],
        Stroke::new(band_w, band_col.gamma_multiply(band_a)),
    );
    painter.text(
        Pos2::new(rect.left() + 4.0, y_of(dwl) - 2.0),
        Align2::LEFT_BOTTOM,
        format!("DWL  z = {dwl:.3}"),
        FontId::proportional(11.0),
        dwl_col,
    );
    painter.text(
        Pos2::new(rect.left() + 4.0, y_of(band_top) - 2.0),
        Align2::LEFT_BOTTOM,
        format!(
            "band top  z = {band_top:.3}   (band {band:.3}{})",
            if band_blank { ", default" } else { "" }
        ),
        FontId::proportional(11.0),
        band_col,
    );

    // Grab the nearer handle on press, then track the pointer.
    if resp.drag_started() {
        if let Some(pos) = resp.interact_pointer_pos() {
            let near_band = (pos.y - y_of(band_top)).abs();
            let near_dwl = (pos.y - y_of(dwl)).abs();
            d.dragging = Some(if near_band < near_dwl { 1 } else { 0 });
        }
    }
    if resp.dragged() {
        if let Some(pos) = resp.interact_pointer_pos() {
            let z = z_of(pos.y).clamp(z_lo, z_hi);
            match d.dragging {
                Some(1) => d.band = fmt3((z - dwl).max(0.0)),
                _ => d.waterline = fmt3(z),
            }
        }
    }
    if resp.drag_stopped() {
        d.dragging = None;
    }
}

/// Format a coordinate for a text field: up to 3 decimals, trailing zeros
/// trimmed.
fn fmt3(v: f32) -> String {
    let mut s = format!("{v:.3}");
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    s
}

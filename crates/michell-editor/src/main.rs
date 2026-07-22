//! A small egui app for authoring `michell` sweep-study manifests. It edits a
//! typed model with a form, shows a live JSON preview and validation, and
//! reads/writes the same `.json` files the CLI consumes:
//!
//! ```text
//! michell sweep study.json
//! ```

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use michell_editor::model::Manifest;
use michell_editor::validate::Level;
use michell_editor::{jsonio, ui, validate};
use std::path::PathBuf;

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

struct EditorApp {
    manifest: Manifest,
    path: Option<PathBuf>,
    dirty: bool,
    status: String,
    /// Cached pretty-printed JSON (or the serialization error), refreshed
    /// whenever the model changes.
    preview: Result<String, String>,
    show_preview: bool,
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
        }
    }
}

impl EditorApp {
    fn refresh_preview(&mut self) {
        self.preview = jsonio::to_string(&self.manifest);
    }

    fn title_file(&self) -> String {
        self.path
            .as_ref()
            .and_then(|p| p.file_name())
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "untitled".into())
    }

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
        if let Some(name) = self.suggested_name() {
            dialog = dialog.set_file_name(name);
        }
        if let Some(path) = dialog.save_file() {
            self.save_to(path);
        }
    }

    fn suggested_name(&self) -> Option<String> {
        if !self.manifest.name.trim().is_empty() {
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
            Some(format!("{}.json", slug.trim_matches('-')))
        } else {
            Some("study.json".into())
        }
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
}

impl eframe::App for EditorApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
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
                ui.separator();
                let title = format!(
                    "{}{}",
                    self.title_file(),
                    if self.dirty { " *" } else { "" }
                );
                ui.label(title);
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
                let changed = ui::manifest_form(ui, &mut self.manifest);
                if changed {
                    self.dirty = true;
                    self.refresh_preview();
                }
            });
        });
    }
}

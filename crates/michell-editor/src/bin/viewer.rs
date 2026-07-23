//! `michell-viewer` — a standalone egui viewer for `.msw` sweep archives.
//!
//! Open it on a file (`michell-viewer study.msw`) or launch it bare and pick a
//! file from the File menu. The rendering is the same `michell_editor::viewer`
//! widget the editor pops open when a binary sweep finishes.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use michell_editor::viewer::Viewer;
use std::path::PathBuf;

fn main() -> eframe::Result<()> {
    // An optional path argument opens straight into a file.
    let initial = std::env::args().nth(1).map(PathBuf::from);
    let viewer = match initial {
        Some(p) => Viewer::load(&p),
        None => Viewer::empty("Open a .msw sweep archive from the File menu."),
    };

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([820.0, 720.0])
            .with_min_inner_size([560.0, 420.0])
            .with_title("michell — sweep archive viewer"),
        ..Default::default()
    };
    eframe::run_native(
        "michell-viewer",
        options,
        Box::new(|_cc| Ok(Box::new(ViewerApp { viewer }))),
    )
}

struct ViewerApp {
    viewer: Viewer,
}

impl ViewerApp {
    fn open(&mut self) {
        let mut dialog = rfd::FileDialog::new().add_filter("sweep archive", &["msw"]);
        if let Some(dir) = self
            .viewer
            .path
            .as_ref()
            .and_then(|p| p.parent())
            .map(|p| p.to_path_buf())
        {
            dialog = dialog.set_directory(dir);
        }
        if let Some(path) = dialog.pick_file() {
            self.viewer = Viewer::load(&path);
        }
    }
}

impl eframe::App for ViewerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.input_mut(|i| {
            if i.consume_key(egui::Modifiers::COMMAND, egui::Key::O) {
                self.open();
            }
        });

        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open…").clicked() {
                        self.open();
                        ui.close_menu();
                    }
                });
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(self.viewer.title());
                });
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| self.viewer.ui(ui));
        });
    }
}

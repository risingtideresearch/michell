//! The manifest form. `manifest_form` renders the whole editing surface into a
//! panel and mutates the model in place; it returns what happened this frame
//! (whether anything changed, and whether the user asked to import geometry) so
//! the app can mark the document dirty, refresh the preview, and start a loft.

use crate::model::*;
use egui::{ComboBox, DragValue, Grid, RichText, Ui};
use std::path::Path;

/// What the form did this frame.
#[derive(Default)]
pub struct FormResponse {
    pub changed: bool,
    pub import_clicked: bool,
}

/// `dir` is the directory the manifest is (or will be) saved in — file paths
/// picked via Browse are stored relative to it when possible, matching the
/// "hull files load relative to the manifest" contract.
pub fn manifest_form(ui: &mut Ui, m: &mut Manifest, dir: Option<&Path>) -> FormResponse {
    let mut r = FormResponse::default();

    ui.heading("Study");
    Grid::new("study")
        .num_columns(2)
        .spacing([12.0, 6.0])
        .show(ui, |ui| {
            ui.label("Name");
            r.changed |= ui.text_edit_singleline(&mut m.name).changed();
            ui.end_row();

            ui.label("Fluid");
            ComboBox::from_id_salt("fluid")
                .selected_text(m.fluid.as_str())
                .show_ui(ui, |ui| {
                    for f in [Fluid::Seawater, Fluid::Freshwater] {
                        r.changed |= ui.selectable_value(&mut m.fluid, f, f.as_str()).changed();
                    }
                });
            ui.end_row();
        });

    ui.add_space(12.0);
    r.changed |= hulls_section(ui, m, dir, &mut r.import_clicked);

    ui.add_space(12.0);
    let ids: Vec<String> = m.hulls.iter().map(|h| h.id.clone()).collect();
    r.changed |= sweep_section(ui, &mut m.axes, &ids);

    ui.add_space(12.0);
    r.changed |= output_section(ui, &mut m.output, dir);

    ui.add_space(12.0);
    r.changed |= options_section(ui, &mut m.options);

    r
}

fn hulls_section(ui: &mut Ui, m: &mut Manifest, dir: Option<&Path>, import: &mut bool) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.heading("Hulls");
        if ui.button("+ add hull").clicked() {
            m.hulls.push(HullSpec::new());
            changed = true;
        }
        if ui
            .button("Import IGS/STL…")
            .on_hover_text("Loft a CAD/mesh file into full-band .hull bodies and add them")
            .clicked()
        {
            *import = true;
        }
    });

    let mut remove: Option<usize> = None;
    for (i, h) in m.hulls.iter_mut().enumerate() {
        let title = if h.id.trim().is_empty() {
            format!("hull #{}", i + 1)
        } else {
            h.id.clone()
        };
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.strong(title);
                if ui.button("remove").clicked() {
                    remove = Some(i);
                }
            });
            Grid::new(("hull", i))
                .num_columns(2)
                .spacing([12.0, 4.0])
                .show(ui, |ui| {
                    ui.label("id");
                    changed |= ui.text_edit_singleline(&mut h.id).changed();
                    ui.end_row();
                    ui.label("file");
                    ui.horizontal(|ui| {
                        changed |= ui.text_edit_singleline(&mut h.file).changed();
                        if ui.button("Browse…").clicked() {
                            if let Some(p) = pick_file_relative(dir, "hull body", &["hull"]) {
                                h.file = p;
                                changed = true;
                            }
                        }
                    });
                    ui.end_row();
                });
            changed |= ui.checkbox(&mut h.pose.enabled, "base pose").changed();
            if h.pose.enabled {
                Grid::new(("pose", i))
                    .num_columns(4)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        changed |= labelled_drag(ui, "dx", &mut h.pose.dx, 0.01);
                        changed |= labelled_drag(ui, "dy", &mut h.pose.dy, 0.01);
                        ui.end_row();
                        changed |= labelled_drag(ui, "dz", &mut h.pose.dz, 0.01);
                        changed |= labelled_drag(ui, "trim°", &mut h.pose.trim_deg, 0.1);
                        ui.end_row();
                    });
            }
        });
    }
    if let Some(i) = remove {
        m.hulls.remove(i);
        changed = true;
    }
    changed
}

fn sweep_section(ui: &mut Ui, axes: &mut Vec<AxisSpec>, ids: &[String]) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.heading("Sweep axes");
        if ui.button("+ add axis").clicked() {
            axes.push(AxisSpec::new(AxisKind::Weight));
            changed = true;
        }
    });

    let mut remove: Option<usize> = None;
    let mut move_up: Option<usize> = None;
    for (i, a) in axes.iter_mut().enumerate() {
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.horizontal(|ui| {
                ComboBox::from_id_salt(("axis-kind", i))
                    .selected_text(a.kind.label())
                    .show_ui(ui, |ui| {
                        for k in AxisKind::ALL {
                            changed |= ui.selectable_value(&mut a.kind, k, k.label()).changed();
                        }
                    });
                if a.kind == AxisKind::Speed {
                    ComboBox::from_id_salt(("axis-unit", i))
                        .selected_text(a.unit.as_str())
                        .show_ui(ui, |ui| {
                            for u in SpeedUnit::ALL {
                                changed |=
                                    ui.selectable_value(&mut a.unit, u, u.as_str()).changed();
                            }
                        });
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("remove").clicked() {
                        remove = Some(i);
                    }
                    if i > 0 && ui.button("↑").clicked() {
                        move_up = Some(i);
                    }
                });
            });

            if a.kind == AxisKind::Pose {
                changed |= pose_targets(ui, i, a, ids);
            }

            changed |= value_editor(ui, i, &mut a.values);
        });
    }
    if let Some(i) = remove {
        axes.remove(i);
        changed = true;
    }
    if let Some(i) = move_up {
        axes.swap(i - 1, i);
        changed = true;
    }
    changed
}

fn pose_targets(ui: &mut Ui, i: usize, a: &mut AxisSpec, ids: &[String]) -> bool {
    let mut changed = false;
    ui.horizontal_wrapped(|ui| {
        ui.label("param");
        ComboBox::from_id_salt(("pose-param", i))
            .selected_text(a.param.as_str())
            .show_ui(ui, |ui| {
                for p in PoseParam::ALL {
                    changed |= ui.selectable_value(&mut a.param, p, p.as_str()).changed();
                }
            });
    });
    ui.horizontal_wrapped(|ui| {
        ui.label("hulls");
        if ids.is_empty() {
            ui.label(RichText::new("(add a hull first)").italics().weak());
        }
        for id in ids {
            if id.trim().is_empty() {
                continue;
            }
            let mut on = a.targets.iter().any(|t| t == id);
            if ui.selectable_label(on, id).clicked() {
                on = !on;
                if on {
                    a.targets.push(id.clone());
                } else {
                    a.targets.retain(|t| t != id);
                }
                changed = true;
            }
        }
    });
    // Targets that no longer match any hull (e.g. a renamed hull) still show so
    // the user can see and clear them.
    let orphans: Vec<String> = a
        .targets
        .iter()
        .filter(|t| !ids.iter().any(|id| id == *t))
        .cloned()
        .collect();
    if !orphans.is_empty() {
        ui.horizontal_wrapped(|ui| {
            ui.colored_label(ui.visuals().error_fg_color, "unknown targets:");
            for t in orphans {
                if ui.selectable_label(true, format!("{t} ✕")).clicked() {
                    a.targets.retain(|x| x != &t);
                    changed = true;
                }
            }
        });
    }
    changed
}

fn value_editor(ui: &mut Ui, i: usize, v: &mut ValueSpec) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.label("values");
        for (mode, name) in [
            (ValueMode::Scalar, "single"),
            (ValueMode::List, "list"),
            (ValueMode::Range, "range"),
        ] {
            changed |= ui.selectable_value(&mut v.mode, mode, name).changed();
        }
    });
    ui.horizontal_wrapped(|ui| match v.mode {
        ValueMode::Scalar => {
            ui.label("value");
            changed |= text(ui, ("scalar", i), &mut v.scalar, 100.0);
        }
        ValueMode::List => {
            ui.label("values (comma-separated)");
            changed |= text(ui, ("list", i), &mut v.list, 240.0);
        }
        ValueMode::Range => {
            ui.label("start");
            changed |= text(ui, ("rs", i), &mut v.range_start, 70.0);
            ui.label("stop");
            changed |= text(ui, ("re", i), &mut v.range_stop, 70.0);
            ui.label("step");
            changed |= text(ui, ("rstep", i), &mut v.step, 70.0);
            ui.label(RichText::new("(blank = span/5)").weak());
        }
    });
    changed
}

fn output_section(ui: &mut Ui, o: &mut Output, dir: Option<&Path>) -> bool {
    let mut changed = false;
    ui.heading("Output");
    Grid::new("output")
        .num_columns(2)
        .spacing([12.0, 6.0])
        .show(ui, |ui| {
            ui.label("format");
            ComboBox::from_id_salt("out-format")
                .selected_text(o.format.as_str())
                .show_ui(ui, |ui| {
                    for f in [OutputFormat::Csv, OutputFormat::Json] {
                        changed |= ui.selectable_value(&mut o.format, f, f.as_str()).changed();
                    }
                });
            ui.end_row();
            ui.label("file");
            ui.horizontal(|ui| {
                changed |= ui.text_edit_singleline(&mut o.file).changed();
                if ui.button("Browse…").clicked() {
                    let ext = o.format.as_str();
                    if let Some(p) = save_file_relative(dir, ext, &format!("study.{ext}")) {
                        o.file = p;
                        changed = true;
                    }
                }
            });
            ui.end_row();
        });
    ui.label(RichText::new("blank file → results print to stdout").weak());
    changed
}

fn options_section(ui: &mut Ui, o: &mut Options) -> bool {
    let mut changed = false;
    egui::CollapsingHeader::new("Options (advanced)").show(ui, |ui| {
        Grid::new("options")
            .num_columns(2)
            .spacing([12.0, 4.0])
            .show(ui, |ui| {
                changed |= opt_row(ui, "samples (NxM)", &mut o.samples);
                changed |= opt_row(ui, "fit_degree (NxM)", &mut o.fit_degree);
                changed |= opt_row(ui, "fit_control (NxM)", &mut o.fit_control);
                changed |= opt_row(ui, "rel_tol", &mut o.rel_tol);
                changed |= opt_row(ui, "form_factor", &mut o.form_factor);
                changed |= opt_row(ui, "gravity", &mut o.gravity);
                changed |= opt_row(ui, "rho (fluid density)", &mut o.rho);
                changed |= opt_row(ui, "nu (kinematic visc.)", &mut o.nu);
            });
    });
    changed
}

fn opt_row(ui: &mut Ui, label: &str, f: &mut OptField) -> bool {
    let mut changed = ui.checkbox(&mut f.enabled, label).changed();
    ui.add_enabled_ui(f.enabled, |ui| {
        changed |= ui.text_edit_singleline(&mut f.text).changed();
    });
    ui.end_row();
    changed
}

// --- small helpers ---

fn labelled_drag(ui: &mut Ui, label: &str, v: &mut f64, speed: f64) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.label(label);
        changed = ui.add(DragValue::new(v).speed(speed)).changed();
    });
    changed
}

fn text(ui: &mut Ui, id: impl std::hash::Hash, s: &mut String, width: f32) -> bool {
    ui.add_sized([width, 20.0], egui::TextEdit::singleline(s).id_salt(id))
        .changed()
}

/// Express `target` relative to the manifest directory when it sits inside it,
/// so the stored path stays portable; otherwise keep it absolute.
pub fn rel_to(dir: Option<&Path>, target: &Path) -> String {
    if let Some(d) = dir {
        if let Ok(stripped) = target.strip_prefix(d) {
            return stripped.to_string_lossy().into_owned();
        }
    }
    target.to_string_lossy().into_owned()
}

fn pick_file_relative(dir: Option<&Path>, name: &str, exts: &[&str]) -> Option<String> {
    let mut d = rfd::FileDialog::new().add_filter(name, exts);
    if let Some(base) = dir {
        d = d.set_directory(base);
    }
    d.pick_file().map(|p| rel_to(dir, &p))
}

fn save_file_relative(dir: Option<&Path>, ext: &str, suggested: &str) -> Option<String> {
    let mut d = rfd::FileDialog::new()
        .add_filter(ext, &[ext])
        .set_file_name(suggested);
    if let Some(base) = dir {
        d = d.set_directory(base);
    }
    d.save_file().map(|p| rel_to(dir, &p))
}

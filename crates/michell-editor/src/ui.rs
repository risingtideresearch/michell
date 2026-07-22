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
    let hull_ids: Vec<String> = m.hulls.iter().map(|h| h.id.clone()).collect();
    let point_ids = m.point_ids();
    r.changed |= sweep_section(ui, &mut m.axes, &hull_ids, &point_ids);

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

            changed |= pose_editor(ui, i, &mut h.pose);
            changed |= load_editor(ui, i, &mut h.load);
            changed |= points_editor(ui, i, &mut h.points);
        });
    }
    if let Some(i) = remove {
        m.hulls.remove(i);
        changed = true;
    }
    changed
}

fn pose_editor(ui: &mut Ui, i: usize, pose: &mut Pose) -> bool {
    let mut changed = ui.checkbox(&mut pose.enabled, "base pose").changed();
    if pose.enabled {
        Grid::new(("pose", i))
            .num_columns(4)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                changed |= labelled_drag(ui, "dx", &mut pose.dx, 0.01);
                changed |= labelled_drag(ui, "dy", &mut pose.dy, 0.01);
                ui.end_row();
                changed |= labelled_drag(ui, "dz", &mut pose.dz, 0.01);
                changed |= labelled_drag(ui, "trim°", &mut pose.trim_deg, 0.1);
                ui.end_row();
                changed |= labelled_drag(ui, "scale", &mut pose.scale, 0.01);
                ui.end_row();
            });
    }
    changed
}

fn load_editor(ui: &mut Ui, i: usize, load: &mut Load) -> bool {
    let mut changed = ui
        .checkbox(&mut load.enabled, "load (mass + CG)")
        .on_hover_text("Total mass and centre of gravity in the hull's own frame")
        .changed();
    if load.enabled {
        Grid::new(("load", i))
            .num_columns(4)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                changed |= labelled_drag(ui, "mass [kg]", &mut load.mass, 1.0);
                changed |= labelled_drag(ui, "vcg [m]", &mut load.vcg, 0.01);
                ui.end_row();
                changed |= ui.checkbox(&mut load.lcg_set, "lcg [m]").changed();
                ui.add_enabled_ui(load.lcg_set, |ui| {
                    changed |= ui.add(DragValue::new(&mut load.lcg).speed(0.01)).changed();
                });
                ui.label(RichText::new("(blank = midship)").weak());
                ui.end_row();
            });
    }
    changed
}

fn points_editor(ui: &mut Ui, i: usize, points: &mut Vec<PointLoad>) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.label(format!("point loads ({})", points.len()));
        if ui.button("+ point").clicked() {
            points.push(PointLoad::default());
            changed = true;
        }
    });
    let mut remove: Option<usize> = None;
    for (j, p) in points.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            ui.label("id");
            changed |= ui
                .add_sized(
                    [90.0, 20.0],
                    egui::TextEdit::singleline(&mut p.id).id_salt(("pid", i, j)),
                )
                .changed();
            changed |= labelled_drag(ui, "mass", &mut p.mass, 1.0);
            changed |= labelled_drag(ui, "dx", &mut p.dx, 0.01);
            changed |= labelled_drag(ui, "dy", &mut p.dy, 0.01);
            changed |= labelled_drag(ui, "dz", &mut p.dz, 0.01);
            if ui.button("✕").clicked() {
                remove = Some(j);
            }
        });
    }
    if let Some(j) = remove {
        points.remove(j);
        changed = true;
    }
    changed
}

fn sweep_section(
    ui: &mut Ui,
    axes: &mut Vec<AxisSpec>,
    hull_ids: &[String],
    point_ids: &[String],
) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.heading("Sweep axes");
        if ui.button("+ add axis").clicked() {
            axes.push(AxisSpec::new(AxisKind::Hull));
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

            match a.kind {
                AxisKind::Hull => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("param");
                        ComboBox::from_id_salt(("hull-param", i))
                            .selected_text(a.hull_param.as_str())
                            .show_ui(ui, |ui| {
                                for p in HullParam::ALL {
                                    changed |= ui
                                        .selectable_value(&mut a.hull_param, p, p.as_str())
                                        .changed();
                                }
                            });
                    });
                    changed |=
                        target_picker(ui, i, &mut a.targets, hull_ids, "hulls", "add a hull first");
                }
                AxisKind::Point => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("param");
                        ComboBox::from_id_salt(("point-param", i))
                            .selected_text(a.point_param.as_str())
                            .show_ui(ui, |ui| {
                                for p in PointParam::ALL {
                                    changed |= ui
                                        .selectable_value(&mut a.point_param, p, p.as_str())
                                        .changed();
                                }
                            });
                    });
                    changed |= target_picker(
                        ui,
                        i,
                        &mut a.targets,
                        point_ids,
                        "points",
                        "add a point load to a hull first",
                    );
                }
                _ => {}
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

/// Toggle-button set for selecting one or more target ids (a coupled axis).
/// Also surfaces targets that match no id, so a rename can be cleared.
fn target_picker(
    ui: &mut Ui,
    i: usize,
    targets: &mut Vec<String>,
    ids: &[String],
    label: &str,
    empty_hint: &str,
) -> bool {
    let mut changed = false;
    ui.horizontal_wrapped(|ui| {
        ui.label(label);
        let selectable: Vec<&String> = ids.iter().filter(|s| !s.trim().is_empty()).collect();
        if selectable.is_empty() {
            ui.label(RichText::new(format!("({empty_hint})")).italics().weak());
        }
        for id in selectable {
            let on = targets.iter().any(|t| t == id);
            if ui.selectable_label(on, id).clicked() {
                if on {
                    targets.retain(|t| t != id);
                } else {
                    targets.push(id.clone());
                }
                changed = true;
            }
        }
    });
    let orphans: Vec<String> = targets
        .iter()
        .filter(|t| !ids.iter().any(|id| id == *t))
        .cloned()
        .collect();
    if !orphans.is_empty() {
        ui.horizontal_wrapped(|ui| {
            ui.colored_label(ui.visuals().error_fg_color, "unknown targets:");
            for t in orphans {
                if ui
                    .selectable_label(true, format!("{t} ✕"))
                    .on_hover_text("click to remove")
                    .clicked()
                {
                    targets.retain(|x| x != &t);
                    changed = true;
                }
            }
        });
    }
    let _ = i;
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

        ui.add_space(6.0);
        changed |= ui
            .checkbox(&mut o.heel.enabled, "heel roll-up metrics")
            .on_hover_text("GZ-curve summaries + resistance rise (equilibrium mode only)")
            .changed();
        if o.heel.enabled {
            Grid::new("heel")
                .num_columns(2)
                .spacing([12.0, 4.0])
                .show(ui, |ui| {
                    ui.label("resistance_angles [deg]");
                    changed |= ui
                        .text_edit_singleline(&mut o.heel.resistance_angles)
                        .changed();
                    ui.end_row();
                    ui.label("gz_step [deg]");
                    changed |= ui.text_edit_singleline(&mut o.heel.gz_step).changed();
                    ui.end_row();
                    ui.label("gz_max [deg]");
                    changed |= ui.text_edit_singleline(&mut o.heel.gz_max).changed();
                    ui.end_row();
                });
        }
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

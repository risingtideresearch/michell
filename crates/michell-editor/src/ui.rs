//! The manifest form. `manifest_form` renders the whole editing surface into a
//! panel and mutates the model in place; it returns `true` if anything changed
//! this frame so the app can mark the document dirty and refresh the preview.

use crate::model::*;
use egui::{ComboBox, DragValue, Grid, RichText, Ui};

pub fn manifest_form(ui: &mut Ui, m: &mut Manifest) -> bool {
    let mut changed = false;

    ui.heading("Study");
    Grid::new("study")
        .num_columns(2)
        .spacing([12.0, 6.0])
        .show(ui, |ui| {
            ui.label("Name");
            changed |= ui.text_edit_singleline(&mut m.name).changed();
            ui.end_row();

            ui.label("Fluid");
            ComboBox::from_id_salt("fluid")
                .selected_text(m.fluid.as_str())
                .show_ui(ui, |ui| {
                    for f in [Fluid::Seawater, Fluid::Freshwater] {
                        changed |= ui.selectable_value(&mut m.fluid, f, f.as_str()).changed();
                    }
                });
            ui.end_row();
        });

    ui.add_space(12.0);
    changed |= hulls_section(ui, m);

    ui.add_space(12.0);
    let ids: Vec<String> = m.hulls.iter().map(|h| h.id.clone()).collect();
    changed |= sweep_section(ui, &mut m.axes, &ids);

    ui.add_space(12.0);
    changed |= output_section(ui, &mut m.output);

    ui.add_space(12.0);
    changed |= options_section(ui, &mut m.options);

    changed
}

fn hulls_section(ui: &mut Ui, m: &mut Manifest) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.heading("Hulls");
        if ui.button("+ add hull").clicked() {
            m.hulls.push(HullSpec::new());
            changed = true;
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
                    changed |= ui.text_edit_singleline(&mut h.file).changed();
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

fn output_section(ui: &mut Ui, o: &mut Output) -> bool {
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
            changed |= ui.text_edit_singleline(&mut o.file).changed();
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

//! Editing model, JSON I/O, validation, and egui form for `michell` sweep
//! manifests. The binary (`main.rs`) is a thin eframe shell over these; the
//! library surface is what the tests exercise.

pub mod jsonio;
pub mod model;
pub mod preview;
pub mod runner;
pub mod ui;
pub mod validate;

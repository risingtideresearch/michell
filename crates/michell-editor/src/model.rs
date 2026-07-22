//! The in-memory editing model for a sweep-study manifest. It mirrors the
//! schema consumed by `michell sweep <study.json>` (see the CLI's
//! `manifest.rs`), but is shaped for editing rather than for a direct
//! (de)serialize — axes are a tagged union, and every "one of value / values /
//! range" choice is an explicit mode so the form can round-trip a partially
//! filled entry without losing what the user typed.

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Fluid {
    Seawater,
    Freshwater,
}

impl Fluid {
    pub fn as_str(self) -> &'static str {
        match self {
            Fluid::Seawater => "seawater",
            Fluid::Freshwater => "freshwater",
        }
    }
}

/// A hull's base pose. The manifest omits `pose` entirely when it is all
/// zeros, so `enabled` tracks whether to emit the object.
#[derive(Clone, Default)]
pub struct Pose {
    pub enabled: bool,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub trim_deg: f64,
}

#[derive(Clone, Default)]
pub struct HullSpec {
    pub id: String,
    pub file: String,
    pub pose: Pose,
}

impl HullSpec {
    pub fn new() -> Self {
        Self::default()
    }
}

/// What an axis drives. `Pose` carries its own target/param; the scalar kinds
/// map straight onto the reserved manifest targets.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AxisKind {
    Speed,
    Weight,
    Lcg,
    Vcg,
    Heel,
    Waterline,
    Pose,
}

impl AxisKind {
    pub const ALL: [AxisKind; 7] = [
        AxisKind::Speed,
        AxisKind::Weight,
        AxisKind::Lcg,
        AxisKind::Vcg,
        AxisKind::Heel,
        AxisKind::Waterline,
        AxisKind::Pose,
    ];

    pub fn label(self) -> &'static str {
        match self {
            AxisKind::Speed => "speed",
            AxisKind::Weight => "weight",
            AxisKind::Lcg => "lcg",
            AxisKind::Vcg => "vcg",
            AxisKind::Heel => "heel",
            AxisKind::Waterline => "waterline",
            AxisKind::Pose => "hull pose",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SpeedUnit {
    Ms,
    Knots,
    Froude,
}

impl SpeedUnit {
    pub const ALL: [SpeedUnit; 3] = [SpeedUnit::Ms, SpeedUnit::Knots, SpeedUnit::Froude];

    /// The token written to the manifest.
    pub fn as_str(self) -> &'static str {
        match self {
            SpeedUnit::Ms => "ms",
            SpeedUnit::Knots => "knots",
            SpeedUnit::Froude => "froude",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PoseParam {
    Dx,
    Dy,
    Dz,
    Spread,
    Trim,
}

impl PoseParam {
    pub const ALL: [PoseParam; 5] = [
        PoseParam::Dx,
        PoseParam::Dy,
        PoseParam::Dz,
        PoseParam::Spread,
        PoseParam::Trim,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            PoseParam::Dx => "dx",
            PoseParam::Dy => "dy",
            PoseParam::Dz => "dz",
            PoseParam::Spread => "spread",
            PoseParam::Trim => "trim",
        }
    }
}

/// How an axis's values are given: a single point, an explicit list, or a
/// start/stop range with an optional step.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ValueMode {
    Scalar,
    List,
    Range,
}

/// The value spec is held as text so half-typed numbers survive a frame and
/// parse errors can be surfaced next to the field rather than swallowed.
#[derive(Clone)]
pub struct ValueSpec {
    pub mode: ValueMode,
    pub scalar: String,
    pub list: String,
    pub range_start: String,
    pub range_stop: String,
    pub step: String,
}

impl Default for ValueSpec {
    fn default() -> Self {
        ValueSpec {
            mode: ValueMode::Range,
            scalar: "0".into(),
            list: String::new(),
            range_start: "0".into(),
            range_stop: "1".into(),
            step: String::new(),
        }
    }
}

#[derive(Clone)]
pub struct AxisSpec {
    pub kind: AxisKind,
    pub unit: SpeedUnit,
    /// Hull ids driven by a pose axis (a multi-hull list is a coupled axis).
    pub targets: Vec<String>,
    pub param: PoseParam,
    pub values: ValueSpec,
}

impl AxisSpec {
    pub fn new(kind: AxisKind) -> Self {
        AxisSpec {
            kind,
            unit: SpeedUnit::Knots,
            targets: Vec::new(),
            param: PoseParam::Spread,
            values: ValueSpec::default(),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum OutputFormat {
    Csv,
    Json,
}

impl OutputFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            OutputFormat::Csv => "csv",
            OutputFormat::Json => "json",
        }
    }
}

#[derive(Clone)]
pub struct Output {
    pub format: OutputFormat,
    /// Empty = write to stdout (the manifest omits `file`).
    pub file: String,
}

impl Default for Output {
    fn default() -> Self {
        Output {
            format: OutputFormat::Csv,
            file: String::new(),
        }
    }
}

/// One optional numeric option, with a checkbox gating whether it is written.
#[derive(Clone)]
pub struct OptField {
    pub enabled: bool,
    pub text: String,
}

impl OptField {
    fn off(default: &str) -> Self {
        OptField {
            enabled: false,
            text: default.into(),
        }
    }
}

/// The `options` object. `samples`/`fit_degree`/`fit_control` are `"NxM"`
/// strings; the rest are scalars.
#[derive(Clone)]
pub struct Options {
    pub samples: OptField,
    pub fit_degree: OptField,
    pub fit_control: OptField,
    pub rel_tol: OptField,
    pub form_factor: OptField,
    pub gravity: OptField,
    pub rho: OptField,
    pub nu: OptField,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            samples: OptField::off("80x64"),
            fit_degree: OptField::off("3x3"),
            fit_control: OptField::off("12x14"),
            rel_tol: OptField::off("1e-5"),
            form_factor: OptField::off("0.05"),
            gravity: OptField::off("9.80665"),
            rho: OptField::off("1025"),
            nu: OptField::off("1.19e-6"),
        }
    }
}

impl Options {
    /// True when nothing is enabled, so `to_json` can omit the whole object.
    pub fn is_empty(&self) -> bool {
        !(self.samples.enabled
            || self.fit_degree.enabled
            || self.fit_control.enabled
            || self.rel_tol.enabled
            || self.form_factor.enabled
            || self.gravity.enabled
            || self.rho.enabled
            || self.nu.enabled)
    }
}

#[derive(Clone)]
pub struct Manifest {
    pub name: String,
    pub fluid: Fluid,
    pub hulls: Vec<HullSpec>,
    pub axes: Vec<AxisSpec>,
    pub output: Output,
    pub options: Options,
}

impl Default for Manifest {
    fn default() -> Self {
        Manifest {
            name: String::new(),
            fluid: Fluid::Seawater,
            hulls: Vec::new(),
            axes: vec![AxisSpec::new(AxisKind::Speed)],
            output: Output::default(),
            options: Options::default(),
        }
    }
}

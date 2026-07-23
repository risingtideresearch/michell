//! The in-memory editing model for a sweep-study manifest. It mirrors the
//! schema consumed by `michell sweep <study.json>` (see the CLI's
//! `manifest.rs`), but is shaped for editing rather than for a direct
//! (de)serialize.
//!
//! The schema, in brief: each hull carries a `pose`, a `load` (mass + centre of
//! gravity), and any number of `points` (discrete masses). A sweep axis is
//! either the reserved `speed`/`waterline`, or it targets one or more hull ids
//! (pose/load params) or point-load ids (mass/offset params) and offsets their
//! base value. Heel is not an axis — it is a roll-up configured under
//! `options.heel`.

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

/// A hull's base pose. The manifest omits `pose` entirely when it is at rest, so
/// `enabled` tracks whether to emit the object. `scale` multiplies the built
/// size (default 1.0); the others are offsets.
#[derive(Clone)]
pub struct Pose {
    pub enabled: bool,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
    pub trim_deg: f64,
    pub scale: f64,
}

impl Default for Pose {
    fn default() -> Self {
        Pose {
            enabled: false,
            dx: 0.0,
            dy: 0.0,
            dz: 0.0,
            trim_deg: 0.0,
            scale: 1.0,
        }
    }
}

/// A hull's load: total mass and centre of gravity in the hull's own frame.
/// `lcg` defaults to the hull's midship (which only the CLI can compute), so it
/// is written only when `lcg_set`.
#[derive(Clone, Default)]
pub struct Load {
    pub enabled: bool,
    pub mass: f64,
    pub lcg_set: bool,
    pub lcg: f64,
    pub vcg: f64,
}

/// A discrete point mass mounted on a hull, offset from its centerpoint
/// (`dx` forward, `dy` to +y, `dz` down).
#[derive(Clone, Default)]
pub struct PointLoad {
    pub id: String,
    pub mass: f64,
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
}

#[derive(Clone, Default)]
pub struct HullSpec {
    pub id: String,
    pub file: String,
    pub pose: Pose,
    pub load: Load,
    pub points: Vec<PointLoad>,
}

/// Computed hydrostatics of a loaded `.hull` body, from `michell info`.
#[derive(Clone, Copy)]
pub struct HullInfo {
    pub length: f64,
    pub beam: f64,
    pub draft: f64,
    pub wetted_surface: f64,
    pub displaced_volume: f64,
}

impl HullSpec {
    pub fn new() -> Self {
        Self::default()
    }
}

/// What an axis drives. `Speed`/`Waterline` are reserved global axes; `Hull`
/// and `Point` target ids and carry their own param.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AxisKind {
    Speed,
    Waterline,
    Hull,
    Point,
}

impl AxisKind {
    pub const ALL: [AxisKind; 4] = [
        AxisKind::Speed,
        AxisKind::Waterline,
        AxisKind::Hull,
        AxisKind::Point,
    ];

    pub fn label(self) -> &'static str {
        match self {
            AxisKind::Speed => "speed",
            AxisKind::Waterline => "waterline",
            AxisKind::Hull => "hull param",
            AxisKind::Point => "point-load param",
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

    pub fn as_str(self) -> &'static str {
        match self {
            SpeedUnit::Ms => "ms",
            SpeedUnit::Knots => "knots",
            SpeedUnit::Froude => "froude",
        }
    }
}

/// A parameter swept on a hull target. Pose params offset (or, for `Scale`,
/// multiply) the base pose; load params offset the base load.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum HullParam {
    Dx,
    Dy,
    Dz,
    Trim,
    Spread,
    Scale,
    Mass,
    Lcg,
    Vcg,
}

impl HullParam {
    pub const ALL: [HullParam; 9] = [
        HullParam::Dx,
        HullParam::Dy,
        HullParam::Dz,
        HullParam::Trim,
        HullParam::Spread,
        HullParam::Scale,
        HullParam::Mass,
        HullParam::Lcg,
        HullParam::Vcg,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            HullParam::Dx => "dx",
            HullParam::Dy => "dy",
            HullParam::Dz => "dz",
            HullParam::Trim => "trim",
            HullParam::Spread => "spread",
            HullParam::Scale => "scale",
            HullParam::Mass => "mass",
            HullParam::Lcg => "lcg",
            HullParam::Vcg => "vcg",
        }
    }
}

/// A parameter swept on a point-load target.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PointParam {
    Mass,
    Dx,
    Dy,
    Dz,
}

impl PointParam {
    pub const ALL: [PointParam; 4] = [
        PointParam::Mass,
        PointParam::Dx,
        PointParam::Dy,
        PointParam::Dz,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            PointParam::Mass => "mass",
            PointParam::Dx => "dx",
            PointParam::Dy => "dy",
            PointParam::Dz => "dz",
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
    /// Target ids: hull ids for a `Hull` axis, point-load ids for a `Point`
    /// axis (a multi-id list is a coupled axis). Unused for speed/waterline.
    pub targets: Vec<String>,
    pub hull_param: HullParam,
    pub point_param: PointParam,
    pub values: ValueSpec,
}

impl AxisSpec {
    pub fn new(kind: AxisKind) -> Self {
        AxisSpec {
            kind,
            unit: SpeedUnit::Knots,
            targets: Vec::new(),
            hull_param: HullParam::Mass,
            point_param: PointParam::Mass,
            values: ValueSpec::default(),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Csv,
    Json,
    /// Self-contained binary sweep archive (`.msw`): manifest, hull files, and
    /// per-row metrics + GZ curve + spectrum in one file, viewable in the app.
    Binary,
}

impl OutputFormat {
    pub const ALL: [OutputFormat; 3] =
        [OutputFormat::Csv, OutputFormat::Json, OutputFormat::Binary];

    /// The `output.format` token written to the manifest.
    pub fn as_str(self) -> &'static str {
        match self {
            OutputFormat::Csv => "csv",
            OutputFormat::Json => "json",
            OutputFormat::Binary => "binary",
        }
    }

    /// The file extension the format produces.
    pub fn ext(self) -> &'static str {
        match self {
            OutputFormat::Csv => "csv",
            OutputFormat::Json => "json",
            OutputFormat::Binary => "msw",
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
            // The binary archive is the default: it captures everything a study
            // produces (metrics, GZ curves, spectra) in one file the app can
            // open and plot.
            format: OutputFormat::Binary,
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

/// The heel roll-up config (`options.heel`). Written only when `enabled`; each
/// field is emitted only when non-empty, else the CLI default applies.
#[derive(Clone)]
pub struct Heel {
    pub enabled: bool,
    /// Comma-separated heel angles (deg) for the resistance-rise columns.
    pub resistance_angles: String,
    pub gz_step: String,
    pub gz_max: String,
}

impl Default for Heel {
    fn default() -> Self {
        Heel {
            enabled: false,
            resistance_angles: "5, 10".into(),
            gz_step: "2.5".into(),
            gz_max: "90".into(),
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
    pub heel: Heel,
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
            heel: Heel::default(),
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
            || self.nu.enabled
            || self.heel.enabled)
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

impl Manifest {
    /// All point-load ids across every hull (for point-axis target pickers and
    /// uniqueness checks).
    pub fn point_ids(&self) -> Vec<String> {
        self.hulls
            .iter()
            .flat_map(|h| h.points.iter().map(|p| p.id.clone()))
            .collect()
    }
}

//! Loading hull specs into a working fleet of sectional hulls.
//!
//! Every input becomes source geometry that can be re-posed and re-cut
//! ([`michell::source::HullSource`]), and each of its hulls is cut into
//! sections at the design waterline:
//!
//! * `*.igs` / `*.iges`: NURBS patches, clustered into hulls;
//! * `*.stl`: a triangle mesh (binary or ASCII; `--units` gives its scale),
//!   clustered likewise.

use crate::formats::looks_binary_stl;
use michell::iges::{
    source_fleet, HullPose, Platform, SectionalOptions, SectionalReport,
};
use michell::sectional::SectionalHull;
use michell::source::{HullSource, SourceHull};
use michell::Placement;
use std::collections::HashMap;

/// What a file held.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Iges,
    Stl,
}

/// Import settings shared by every command that reads a hull.
#[derive(Clone, Copy)]
pub struct LoadSettings {
    /// CAD height of the design waterline.
    pub waterline_z: f64,
    pub centerplane: Option<f64>,
    pub stations: usize,
    pub rays: usize,
    /// Scale to metres for STL, which carries no units.
    pub units: Option<f64>,
}

impl Default for LoadSettings {
    fn default() -> Self {
        let d = SectionalOptions::default();
        LoadSettings {
            waterline_z: 0.0,
            centerplane: None,
            stations: d.stations,
            rays: d.rays,
            units: None,
        }
    }
}

impl LoadSettings {
    /// Sectioning options for a file whose design waterline is at
    /// `waterline_z`.
    pub fn sectional(&self, waterline_z: f64) -> SectionalOptions {
        SectionalOptions {
            waterline_z,
            centerplane: self.centerplane,
            stations: self.stations,
            rays: self.rays,
            ..Default::default()
        }
    }
}

/// A parsed file: its hulls as re-cuttable geometry.
pub struct SourceFile {
    pub path: String,
    pub kind: Kind,
    pub source: Box<dyn HullSource>,
    /// CAD height of the design waterline in this file's frame.
    pub waterline_z: f64,
}

/// One hull of the working fleet, cut at its design pose.
pub struct Member {
    pub path: String,
    /// Index into [`Fleet::files`].
    pub file: usize,
    /// Hull index within that file's source.
    pub index: usize,
    pub hull: SectionalHull,
    /// Where the fleet places the hull (detected centreplane plus the spec's
    /// shifts).
    pub placement: Placement,
    /// The spec's shift from the placement the hull was cut at.
    pub shift: Placement,
    pub report: SectionalReport,
}

pub struct Fleet {
    pub files: Vec<SourceFile>,
    pub members: Vec<Member>,
}

impl Fleet {
    /// `(hull, placement)` pairs, as the physics takes them.
    pub fn hulls(&self) -> Vec<(&SectionalHull, Placement)> {
        self.members.iter().map(|m| (&m.hull, m.placement)).collect()
    }

    pub fn source(&self, m: &Member) -> &dyn HullSource {
        self.files[m.file].source.as_ref()
    }

    /// Every member as source geometry at the default design pose.
    pub fn source_hulls(&self) -> Vec<SourceHull<'_>> {
        self.members
            .iter()
            .map(|m| SourceHull {
                source: self.source(m),
                index: m.index,
                waterline_z: self.files[m.file].waterline_z,
                pose: HullPose::default(),
            })
            .collect()
    }

    /// Reference length for Froude numbers: the longest hull.
    pub fn l_ref(&self) -> f64 {
        self.members
            .iter()
            .map(|m| m.hull.length())
            .fold(0.0f64, f64::max)
    }
}

/// Placement request from a hull spec suffix.
#[derive(Default, Clone, Copy)]
pub struct SpecPlacement {
    /// Absolute centerplane position (single-hull files only).
    pub y_abs: Option<f64>,
    /// Transverse shift applied to the file's detected placements.
    pub dy: f64,
    /// Longitudinal shift added to the file's x coordinates.
    pub dx: f64,
}

/// Parse `path` or `path@key=V,...` (keys: `y` absolute centerplane,
/// `dy` transverse shift, `x`/`dx` longitudinal shift).
pub fn parse_hull_spec(spec: &str) -> Result<(String, SpecPlacement), String> {
    let Some((path, rest)) = spec.split_once('@') else {
        return Ok((spec.to_string(), SpecPlacement::default()));
    };
    let mut place = SpecPlacement::default();
    for part in rest.split(',') {
        let (k, v) = part
            .split_once('=')
            .ok_or_else(|| format!("bad placement {rest:?}: expected key=value pairs"))?;
        let val: f64 = v
            .trim()
            .parse()
            .map_err(|_| format!("bad placement value {v:?} in {spec:?}"))?;
        match k.trim() {
            "y" => place.y_abs = Some(val),
            "dy" => place.dy = val,
            "x" | "dx" => place.dx = val,
            other => return Err(format!("unknown placement key {other:?} (use y, dy, x/dx)")),
        }
    }
    if place.y_abs.is_some() && place.dy != 0.0 {
        return Err(format!(
            "{spec:?}: give either y (absolute) or dy (shift), not both"
        ));
    }
    Ok((path.to_string(), place))
}

/// Parse one file into source geometry.
pub fn open_source(path: &str, settings: &LoadSettings) -> Result<SourceFile, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    open_source_bytes(path, bytes, settings)
}

/// [`open_source`] on contents already in memory; `path` names the file and
/// sniffs its format by extension.
pub fn open_source_bytes(
    path: &str,
    bytes: Vec<u8>,
    settings: &LoadSettings,
) -> Result<SourceFile, String> {
    let lower = path.to_ascii_lowercase();
    // STL: by extension or binary layout (binary STL is not UTF-8).
    if lower.ends_with(".stl") || looks_binary_stl(&bytes) || std::str::from_utf8(&bytes).is_err()
    {
        let scale = settings.units.ok_or_else(|| {
            format!(
                "{path}: STL files carry no units; pass --units mm|cm|m|in|ft \
                 (or a scale to metres)"
            )
        })?;
        let src = michell::stl::mesh_fleet(&bytes, scale, settings.waterline_z)
            .map_err(|e| format!("{path}: STL import failed: {e}"))?;
        return Ok(SourceFile {
            path: path.into(),
            kind: Kind::Stl,
            source: Box::new(src),
            waterline_z: settings.waterline_z,
        });
    }
    let text = String::from_utf8(bytes).expect("checked utf8");
    let first = text
        .lines()
        .find(|l| !l.trim().is_empty())
        .unwrap_or("")
        .trim_end();
    if first.starts_with("michell-hull")
        || first.starts_with("michell-offsets")
        || text.trim_start().starts_with('{')
    {
        return Err(format!(
            "{path}: .hull control nets, offsets tables and sample grids are no longer \
             read; import the CAD geometry (IGES or STL) instead"
        ));
    }
    let looks_iges = lower.ends_with(".igs")
        || lower.ends_with(".iges")
        || first.len() >= 73 && matches!(first.as_bytes()[72], b'S' | b'G');
    if looks_iges {
        let src = source_fleet(&text, settings.waterline_z)
            .map_err(|e| format!("{path}: IGES import failed: {e}"))?;
        return Ok(SourceFile {
            path: path.into(),
            kind: Kind::Iges,
            source: Box::new(src),
            waterline_z: settings.waterline_z,
        });
    }
    Err(format!(
        "cannot determine the format of {path}: expected an IGES or STL file"
    ))
}

/// Load a fleet of hull specs, reading each unique file once. A file
/// containing several hulls contributes all of them, each at its detected
/// placement plus the spec's offset. Hulls that cannot be sectioned (a sliver
/// of appendage geometry, say) are reported and skipped, as are dry ones; a
/// fleet with none left is an error.
pub fn load(specs: &[String], settings: &LoadSettings) -> Result<Fleet, String> {
    let mut files: Vec<SourceFile> = Vec::new();
    let mut cut: Vec<Vec<(usize, michell::iges::SectionalImport)>> = Vec::new();
    let mut by_path: HashMap<String, usize> = HashMap::new();
    let mut members = Vec::new();
    for spec in specs {
        let (path, sp) = parse_hull_spec(spec)?;
        let fi = match by_path.get(&path) {
            Some(&i) => i,
            None => {
                let f = open_source(&path, settings)?;
                let opts = settings.sectional(f.waterline_z);
                let mut hulls = Vec::new();
                for i in 0..f.source.len() {
                    match f.source.situate_sectional(
                        i,
                        f.waterline_z,
                        &HullPose::default(),
                        &Platform::default(),
                        &opts,
                    ) {
                        Ok(Some(h)) => hulls.push((i, h)),
                        Ok(None) => {}
                        Err(e) => eprintln!("{path}: hull {} not sectioned, skipped: {e}", i + 1),
                    }
                }
                if hulls.is_empty() {
                    return Err(format!("{path}: no hull is wetted at this waterline"));
                }
                files.push(f);
                cut.push(hulls);
                by_path.insert(path.clone(), files.len() - 1);
                files.len() - 1
            }
        };
        let hulls = &cut[fi];
        if sp.y_abs.is_some() && hulls.len() > 1 {
            return Err(format!(
                "{spec:?}: absolute y placement is ambiguous for a file with {} hulls; \
                 use dy=SHIFT instead",
                hulls.len()
            ));
        }
        for (index, h) in hulls {
            let placement = Placement {
                x: h.placement.x + sp.dx,
                y: sp.y_abs.unwrap_or(h.placement.y + sp.dy),
            };
            members.push(Member {
                path: path.clone(),
                file: fi,
                index: *index,
                hull: h.hull.clone(),
                placement,
                shift: Placement {
                    x: placement.x - h.placement.x,
                    y: placement.y - h.placement.y,
                },
                report: h.report.clone(),
            });
        }
    }
    Ok(Fleet { files, members })
}

/// Largest beam [m] of a sectional hull: twice its widest section.
pub fn max_beam(hull: &SectionalHull) -> f64 {
    2.0 * hull
        .curves()
        .flat_map(|(_, c)| c.iter().map(|p| p.0))
        .fold(0.0f64, f64::max)
}

/// What the reports print about where a hull came from and how it was cut.
pub fn describe(fleet: &Fleet, m: &Member) -> Vec<String> {
    let r = &m.report;
    let sides = if r.two_sided {
        format!("both sides averaged about y = {:.4} m", r.centerplane)
    } else if r.mirrored {
        format!("port half mirrored about y = {:.4} m", r.centerplane)
    } else {
        format!("one side about y = {:.4} m", r.centerplane)
    };
    let what = match fleet.files[m.file].kind {
        Kind::Iges => format!("IGES ({} patches, units scale {})", r.patches, r.units_scale),
        Kind::Stl => format!("STL ({} triangles, units scale {})", r.patches, r.units_scale),
    };
    let mut lines = vec![format!(
        "source: {what}, {} stations over x {:.4}..{:.4} m, {sides}",
        r.stations, r.x_range.0, r.x_range.1
    )];
    if r.dropped_stations > 0 || r.max_asymmetry > 1e-3 * r.draft.max(1e-9) {
        lines.push(format!(
            "sections: {} interior stations dropped, port/starboard differ by up to {:.3e} m",
            r.dropped_stations, r.max_asymmetry
        ));
    }
    if let Some(t) = &r.transom {
        lines.push(format!(
            "transom: immersed at x {:.4} m, {:.1}% of max section, equivalent \
             depth {:.4} m, waterline half-beam {:.4} m (closed by --transom)",
            t.x,
            100.0 * t.area / m.hull.max_section_area().max(f64::MIN_POSITIVE),
            t.depth,
            t.half_beam
        ));
    }
    lines
}

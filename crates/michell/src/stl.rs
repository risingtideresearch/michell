//! STL mesh import: triangle soup as a sampling source for the same
//! clip → sample → loft pipeline as IGES.
//!
//! The mesh is never differentiated — half-beams are extracted by casting
//! transverse rays against the triangles (exact, no Newton iteration, and
//! empty results *are* the footprint), then lofted by least squares like
//! every other source. Loft quality is therefore tied to the export's chord
//! tolerance: fine CAD tessellations (≤ ~1 mm sag) are indistinguishable
//! from IGES; heavily decimated meshes add geometry noise that wave
//! resistance is sensitive to.
//!
//! The **sectional** import ([`MeshFleet::situate_sectional`],
//! [`import_sectional`]) is the IGES one run on the mesh itself: the same
//! stations, ray fan, outermost-fold reach, end-cap rejection, hull ends,
//! centreplane and transom detection, with each ray's hits taken straight
//! from the station plane's cut through the facets — the mesh is the exact
//! geometry here, so there is nothing to polish onto. Section accuracy is
//! the export's chord error.
//!
//! Conventions match the IGES path: x longitudinal, z **up**, y transverse;
//! one shell (or a whole multihull) per file. STL has no units field, so the
//! caller must supply the scale to metres.

use crate::error::{Error, Result};
use crate::fit::fit_grid;
use crate::grid::SampleGrid;
use crate::iges::{
    collect_sectional, rotate_xz, sectional_mesh, HullPose, ImportOptions, ImportReport,
    ImportedHull, Platform, PoseMap, PosedMesh, SectionalFleet, SectionalImport, SectionalOptions,
    SectionalState, SituatedFleet, Tessellation,
};
use crate::michell::Placement;

/// One triangle, vertices in CAD coordinates (metres after scaling).
pub type Tri = [[f64; 3]; 3];

/// Parse STL bytes (binary or ASCII, auto-detected), scaling coordinates by
/// `units_scale` (STL carries no units; e.g. 0.001 for a millimetre export).
pub fn parse_stl(bytes: &[u8], units_scale: f64) -> Result<Vec<Tri>> {
    if !(units_scale.is_finite() && units_scale > 0.0) {
        return Err(Error::InvalidConditions(
            "units scale must be finite and positive".into(),
        ));
    }
    // Binary detection: exact size match beats the unreliable "solid" prefix.
    if bytes.len() >= 84 {
        let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        if bytes.len() == 84 + 50 * n {
            return parse_binary(bytes, n, units_scale);
        }
    }
    let text = std::str::from_utf8(bytes)
        .map_err(|_| Error::Parse("not a binary STL, and not valid UTF-8 ASCII STL".into()))?;
    if !text.trim_start().starts_with("solid") {
        return Err(Error::Parse(
            "not an STL file (no binary size match, no `solid` header)".into(),
        ));
    }
    parse_ascii(text, units_scale)
}

fn parse_binary(bytes: &[u8], n: usize, scale: f64) -> Result<Vec<Tri>> {
    let mut tris = Vec::with_capacity(n);
    for i in 0..n {
        let at = 84 + 50 * i + 12; // skip the normal
        let mut tri = [[0.0f64; 3]; 3];
        for (v, tv) in tri.iter_mut().enumerate() {
            for (c, coord) in tv.iter_mut().enumerate() {
                let o = at + 12 * v + 4 * c;
                let f = f32::from_le_bytes([bytes[o], bytes[o + 1], bytes[o + 2], bytes[o + 3]]);
                if !f.is_finite() {
                    return Err(Error::Parse(format!("non-finite vertex in triangle {i}")));
                }
                *coord = f as f64 * scale;
            }
        }
        tris.push(tri);
    }
    Ok(tris)
}

fn parse_ascii(text: &str, scale: f64) -> Result<Vec<Tri>> {
    let mut tris = Vec::new();
    let mut verts: Vec<[f64; 3]> = Vec::with_capacity(3);
    for (ln, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if let Some(rest) = line.strip_prefix("vertex") {
            let vals: Vec<f64> = rest
                .split_whitespace()
                .map(|t| t.parse::<f64>())
                .collect::<std::result::Result<_, _>>()
                .map_err(|_| Error::Parse(format!("line {}: bad vertex", ln + 1)))?;
            if vals.len() != 3 || vals.iter().any(|v| !v.is_finite()) {
                return Err(Error::Parse(format!(
                    "line {}: vertex needs 3 finite coordinates",
                    ln + 1
                )));
            }
            verts.push([vals[0] * scale, vals[1] * scale, vals[2] * scale]);
            if verts.len() == 3 {
                tris.push([verts[0], verts[1], verts[2]]);
                verts.clear();
            }
        } else if line.starts_with("outer loop") {
            verts.clear();
        }
    }
    if !verts.is_empty() {
        return Err(Error::Parse("dangling vertices at end of ASCII STL".into()));
    }
    Ok(tris)
}

/// Parsed and clustered mesh, kept in CAD coordinates so hulls can be
/// re-situated repeatedly — the mesh analogue of [`crate::iges::SourceFleet`].
pub struct MeshFleet {
    units_scale: f64,
    /// Triangles of each detected hull, CAD frame, metres, sorted by y.
    hulls: Vec<Vec<Tri>>,
    /// Each hull as an indexed mesh (shared vertices merged), for posing and
    /// sectioning.
    meshes: Vec<Tessellation>,
    /// Each hull's vertex x and y mids: the default trim and scale pivots.
    mids: Vec<(f64, f64)>,
}

/// Parse an STL file and cluster its triangles into hulls at a reference
/// waterline (CAD z, up). Triangles connect through shared vertices; the
/// resulting components merge by wetted-bbox proximity, exactly like IGES
/// patches, so dry structure cannot bridge two hulls.
pub fn mesh_fleet(bytes: &[u8], units_scale: f64, reference_waterline: f64) -> Result<MeshFleet> {
    let tris = parse_stl(bytes, units_scale)?;
    if tris.is_empty() {
        return Err(Error::Parse("the STL contains no triangles".into()));
    }

    // Union triangles sharing (quantized) vertices.
    let mut parent: Vec<usize> = (0..tris.len()).collect();
    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let mut seen: std::collections::HashMap<(i64, i64, i64), usize> =
        std::collections::HashMap::new();
    let q = 1e8; // ~10 nm quantization: exported shared vertices coincide
    for (ti, tri) in tris.iter().enumerate() {
        for v in tri {
            let key = (
                (v[0] * q).round() as i64,
                (v[1] * q).round() as i64,
                (v[2] * q).round() as i64,
            );
            match seen.get(&key) {
                Some(&other) => {
                    let (a, b) = (find(&mut parent, ti), find(&mut parent, other));
                    parent[a] = b;
                }
                None => {
                    seen.insert(key, ti);
                }
            }
        }
    }
    // Component id -> triangle indices.
    let mut comp_of: Vec<usize> = (0..tris.len()).map(|i| find(&mut parent, i)).collect();
    let mut roots: Vec<usize> = comp_of.clone();
    roots.sort_unstable();
    roots.dedup();
    for c in comp_of.iter_mut() {
        *c = roots.binary_search(c).expect("root present");
    }
    let ncomp = roots.len();

    // Wetted bbox per component at the reference waterline (z' = wl - z).
    let mut wet_box: Vec<Option<[f64; 6]>> = vec![None; ncomp];
    let mut all_box: Vec<[f64; 6]> = vec![
        [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ];
        ncomp
    ];
    let mut global_wet: Option<[f64; 6]> = None;
    fn grow(b: &mut [f64; 6], x: f64, y: f64, zd: f64) {
        b[0] = b[0].min(x);
        b[1] = b[1].max(x);
        b[2] = b[2].min(y);
        b[3] = b[3].max(y);
        b[4] = b[4].min(zd);
        b[5] = b[5].max(zd);
    }
    for (tri, &c) in tris.iter().zip(&comp_of) {
        let wet = tri.iter().any(|v| reference_waterline - v[2] >= -1e-12);
        for v in tri {
            let (x, y, zd) = (v[0], v[1], reference_waterline - v[2]);
            grow(&mut all_box[c], x, y, zd);
            if wet {
                grow(wet_box[c].get_or_insert([x, x, y, y, zd, zd]), x, y, zd);
                grow(global_wet.get_or_insert([x, x, y, y, zd, zd]), x, y, zd);
            }
        }
    }
    let Some(g) = global_wet else {
        return Err(Error::InvalidGeometry(
            "the mesh lies entirely above the specified waterline".into(),
        ));
    };
    let scale = (g[1] - g[0]).max(g[3] - g[2]).max(g[5] - g[4]);
    if scale <= 0.0 || !scale.is_finite() {
        return Err(Error::InvalidGeometry(
            "the wetted part of the mesh is degenerate".into(),
        ));
    }
    let eps = 0.01 * scale;

    // Merge wetted components by bbox proximity (like IGES patch clusters).
    let mut cluster_of: Vec<Option<usize>> = vec![None; ncomp];
    let mut clusters: Vec<[f64; 6]> = Vec::new();
    for c in 0..ncomp {
        let Some(wb) = wet_box[c] else { continue };
        let touch = |a: &[f64; 6], b: &[f64; 6]| {
            (0..3).all(|k| a[2 * k] - eps <= b[2 * k + 1] && b[2 * k] - eps <= a[2 * k + 1])
        };
        let mut joined = None;
        for (ci, cb) in clusters.iter_mut().enumerate() {
            if touch(&wb, cb) {
                for k in 0..3 {
                    cb[2 * k] = cb[2 * k].min(wb[2 * k]);
                    cb[2 * k + 1] = cb[2 * k + 1].max(wb[2 * k + 1]);
                }
                joined = Some(ci);
                break;
            }
        }
        cluster_of[c] = Some(match joined {
            Some(ci) => ci,
            None => {
                clusters.push(wb);
                clusters.len() - 1
            }
        });
    }
    // Note: single-pass merging can in principle leave two clusters that a
    // later component would bridge; a second pass catches the common cases.
    for _ in 0..2 {
        let mut merged = false;
        let mut i = 0;
        while i < clusters.len() {
            let mut j = i + 1;
            while j < clusters.len() {
                let touch = (0..3).all(|k| {
                    clusters[i][2 * k] - eps <= clusters[j][2 * k + 1]
                        && clusters[j][2 * k] - eps <= clusters[i][2 * k + 1]
                });
                if touch {
                    let cj = clusters[j];
                    for (k, ci) in clusters[i].iter_mut().enumerate() {
                        *ci = if k % 2 == 0 {
                            ci.min(cj[k])
                        } else {
                            ci.max(cj[k])
                        };
                    }
                    for co in cluster_of.iter_mut().flatten() {
                        if *co == j {
                            *co = i;
                        } else if *co > j {
                            *co -= 1;
                        }
                    }
                    clusters.remove(j);
                    merged = true;
                } else {
                    j += 1;
                }
            }
            i += 1;
        }
        if !merged {
            break;
        }
    }
    // Dry components attach to the nearest cluster.
    for c in 0..ncomp {
        if cluster_of[c].is_some() {
            continue;
        }
        let pb = all_box[c];
        let dist = |a: &[f64; 6], b: &[f64; 6]| -> f64 {
            (0..3)
                .map(|k| {
                    let gap = (a[2 * k] - b[2 * k + 1])
                        .max(b[2 * k] - a[2 * k + 1])
                        .max(0.0);
                    gap * gap
                })
                .sum()
        };
        cluster_of[c] = clusters
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| dist(&pb, a).total_cmp(&dist(&pb, b)))
            .map(|(i, _)| i);
    }

    // Gather triangles per cluster, ordered by wetted-y midpoint.
    let mut order: Vec<usize> = (0..clusters.len()).collect();
    order.sort_by(|&a, &b| {
        let mid = |b: &[f64; 6]| (b[2] + b[3]) / 2.0;
        mid(&clusters[a]).total_cmp(&mid(&clusters[b]))
    });
    let rank: Vec<usize> = {
        let mut r = vec![0; order.len()];
        for (pos, &ci) in order.iter().enumerate() {
            r[ci] = pos;
        }
        r
    };
    let mut hulls: Vec<Vec<Tri>> = vec![Vec::new(); clusters.len()];
    for (tri, &c) in tris.iter().zip(&comp_of) {
        if let Some(ci) = cluster_of[c] {
            hulls[rank[ci]].push(*tri);
        }
    }
    let meshes: Vec<Tessellation> = hulls.iter().map(|h| indexed(h)).collect();
    let mids = meshes
        .iter()
        .map(|m| {
            let mid = |k: usize| {
                let (lo, hi) = m
                    .verts
                    .iter()
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                        (lo.min(v[k]), hi.max(v[k]))
                    });
                0.5 * (lo + hi)
            };
            (mid(0), mid(1))
        })
        .collect();
    Ok(MeshFleet {
        units_scale,
        hulls,
        meshes,
        mids,
    })
}

/// Triangle soup as an indexed mesh: coincident vertices (to ~10 nm, as in
/// the clustering) merged, so a pose moves each once.
fn indexed(tris: &[Tri]) -> Tessellation {
    let q = 1e8;
    let mut seen: std::collections::HashMap<[i64; 3], u32> =
        std::collections::HashMap::with_capacity(tris.len());
    let mut verts = Vec::with_capacity(tris.len() / 2 + 3);
    let tris = tris
        .iter()
        .map(|t| {
            t.map(|v| {
                let key = v.map(|c| (c * q).round() as i64);
                *seen.entry(key).or_insert_with(|| {
                    verts.push(v);
                    (verts.len() - 1) as u32
                })
            })
        })
        .collect();
    Tessellation::from_mesh(verts, tris)
}

/// Import every hull of an STL file by sections, at the fixed waterline in
/// `opts` (CAD z, up), coordinates scaled by `units_scale` — the mesh
/// counterpart of [`crate::iges::import_sectional`]. Hulls cluster as in
/// [`mesh_fleet`].
pub fn import_sectional(
    bytes: &[u8],
    units_scale: f64,
    opts: &SectionalOptions,
) -> Result<SectionalFleet> {
    let fleet = mesh_fleet(bytes, units_scale, opts.waterline_z)?;
    collect_sectional(fleet.len(), |i| {
        fleet.situate_sectional(
            i,
            opts.waterline_z,
            &HullPose::default(),
            &Platform::default(),
            opts,
        )
    })
}

impl MeshFleet {
    pub fn len(&self) -> usize {
        self.hulls.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hulls.is_empty()
    }

    /// A hull's vertex x mid: the default pivot of its design trim.
    pub fn x_mid(&self, idx: usize) -> f64 {
        self.mids[idx].0
    }

    pub fn units_scale(&self) -> f64 {
        self.units_scale
    }

    /// Highest z (CAD frame, up) of a hull's vertices.
    pub fn hull_z_top(&self, idx: usize) -> f64 {
        self.hulls[idx]
            .iter()
            .flat_map(|t| t.iter())
            .fold(f64::NEG_INFINITY, |m, v| m.max(v[2]))
    }

    /// Lowest z (CAD frame, up) of a hull's vertices.
    pub fn hull_z_bottom(&self, idx: usize) -> f64 {
        self.hulls[idx]
            .iter()
            .flat_map(|t| t.iter())
            .fold(f64::INFINITY, |m, v| m.min(v[2]))
    }

    fn check_idx(&self, idx: usize) -> Result<()> {
        if idx >= self.hulls.len() {
            return Err(Error::InvalidInput(format!(
                "hull index {idx} out of range ({} hulls)",
                self.hulls.len()
            )));
        }
        Ok(())
    }

    /// The pose as a map on CAD-frame points: the IGES importer's map, with
    /// the mesh's x and y mids as the default trim and scale pivots.
    fn pose_map(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> PoseMap {
        let (xm, ym) = self.mids[idx];
        PoseMap::with_mids(waterline_z, pose, platform, || xm, || ym)
    }

    /// Situate one hull (pose and platform state as for
    /// [`crate::iges::SourceFleet::situate_sectional`]) and build it by
    /// sections cut from the mesh; `Ok(None)` when it is dry.
    pub fn situate_sectional(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &SectionalOptions,
    ) -> Result<Option<SectionalImport>> {
        self.situate_sectional_warm(
            idx,
            waterline_z,
            pose,
            platform,
            opts,
            &mut SectionalState::default(),
        )
    }

    /// [`MeshFleet::situate_sectional`], warm-started from (and updating)
    /// `state` — the cheap way to re-pose (see [`SectionalState`]).
    pub fn situate_sectional_warm(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &SectionalOptions,
        state: &mut SectionalState,
    ) -> Result<Option<SectionalImport>> {
        self.check_idx(idx)?;
        if opts.stations < 8 || opts.rays < 5 {
            return Err(Error::InvalidInput(
                "need at least 8 stations and 5 rays per section".into(),
            ));
        }
        let wl = waterline_z + platform.sinkage;
        let map = self.pose_map(idx, waterline_z, pose, platform);
        let mesh = PosedMesh::new(&self.meshes[idx], |mut p| {
            map.apply(&mut p);
            [p[0], p[1], wl - p[2]]
        });
        sectional_mesh(&mesh, opts, self.units_scale, state)
    }

    /// A hull's mesh at a pose, in the water frame: `x` forward, `y`
    /// transverse, `z` **up** from the effective waterline — the frame of
    /// [`crate::iges::SourceFleet::posed_tessellation`].
    pub fn posed_tessellation(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>)> {
        self.check_idx(idx)?;
        let map = self.pose_map(idx, waterline_z, pose, platform);
        let wl = waterline_z + platform.sinkage;
        let tess = &self.meshes[idx];
        let verts = tess
            .verts
            .iter()
            .map(|&p| {
                let mut q = p;
                map.apply(&mut q);
                [q[0], q[1], q[2] - wl]
            })
            .collect();
        Ok((verts, tess.tris.clone()))
    }

    /// Situate the fleet (same contract as [`crate::iges::SourceFleet::situate`]).
    pub fn situate(
        &self,
        waterline_z: f64,
        poses: &[HullPose],
        platform: &Platform,
        opts: &ImportOptions,
    ) -> Result<SituatedFleet> {
        if poses.len() != self.hulls.len() {
            return Err(Error::InvalidInput(format!(
                "{} poses supplied for {} hulls",
                poses.len(),
                self.hulls.len()
            )));
        }
        if self.hulls.len() > 1 && opts.centerplane.is_some() {
            return Err(Error::InvalidConditions(format!(
                "a centerplane override is ambiguous: {} separate hulls detected",
                self.hulls.len()
            )));
        }
        let mut members = Vec::new();
        let mut dry = Vec::new();
        for (hi, pose) in poses.iter().enumerate() {
            match self.situate_hull(hi, waterline_z, pose, platform, opts, &mut |_| {})? {
                Some(m) => members.push(m),
                None => dry.push(hi),
            }
        }
        Ok(SituatedFleet { members, dry })
    }

    /// Situate a single hull; `Ok(None)` when it is dry.
    pub fn situate_one(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &ImportOptions,
    ) -> Result<Option<ImportedHull>> {
        self.situate_one_progress(idx, waterline_z, pose, platform, opts, &mut |_| {})
    }

    /// Like [`MeshFleet::situate_one`], but reports loft-sampling progress as a
    /// fraction in `0.0..=1.0` (one call per station) through `progress`, so a
    /// front-end can show a bar. The final surface fit is not subdivided.
    pub fn situate_one_progress(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &ImportOptions,
        progress: &mut dyn FnMut(f32),
    ) -> Result<Option<ImportedHull>> {
        if idx >= self.hulls.len() {
            return Err(Error::InvalidInput(format!(
                "hull index {idx} out of range ({} hulls)",
                self.hulls.len()
            )));
        }
        self.situate_hull(idx, waterline_z, pose, platform, opts, progress)
    }

    fn situate_hull(
        &self,
        hi: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &ImportOptions,
        progress: &mut dyn FnMut(f32),
    ) -> Result<Option<ImportedHull>> {
        if opts.stations < 8 || opts.waterlines < 6 {
            return Err(Error::InvalidInput(
                "need at least 8 stations and 6 waterlines to sample".into(),
            ));
        }
        let wl = waterline_z + platform.sinkage;
        // Transform vertices (CAD frame), then convert to hull frame
        // (x, y, z' = wl - z).
        let src = &self.hulls[hi];
        let px = pose.pivot_x.unwrap_or_else(|| {
            let (lo, hi) = src
                .iter()
                .flat_map(|t| t.iter())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                    (lo.min(v[0]), hi.max(v[0]))
                });
            0.5 * (lo + hi)
        });
        let (pose_sin, pose_cos) = pose.trim.sin_cos();
        let (plat_sin, plat_cos) = platform.trim.sin_cos();
        let mut tris: Vec<Tri> = Vec::with_capacity(src.len());
        for t in src {
            let mut tri = *t;
            for v in tri.iter_mut() {
                if pose.trim != 0.0 {
                    rotate_xz(v, px, waterline_z, pose_cos, pose_sin);
                }
                v[0] += pose.dx;
                v[1] += pose.dy;
                v[2] -= pose.dz;
                if platform.trim != 0.0 {
                    rotate_xz(v, platform.pivot_x, wl, plat_cos, plat_sin);
                }
                v[2] = wl - v[2]; // hull frame: depth below water
            }
            tris.push(tri);
        }

        // Wetted statistics (mesh extremes are at vertices).
        let mut draft = 0.0f64;
        let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
        let (mut y_lo, mut y_hi) = (f64::INFINITY, f64::NEG_INFINITY);
        let mut y_sum = 0.0f64;
        let mut wet_count = 0usize;
        for t in &tris {
            for v in t {
                if v[2] >= -1e-12 {
                    draft = draft.max(v[2]);
                    x_min = x_min.min(v[0]);
                    x_max = x_max.max(v[0]);
                    y_lo = y_lo.min(v[1]);
                    y_hi = y_hi.max(v[1]);
                    y_sum += v[1];
                    wet_count += 1;
                }
            }
        }
        if wet_count == 0 {
            return Ok(None);
        }
        if !(draft > 0.0 && x_max > x_min) {
            return Err(Error::InvalidGeometry(
                "a hull's wetted mesh is degenerate (zero draft or length)".into(),
            ));
        }
        let length = x_max - x_min;
        let scale = length.max(draft).max(y_hi - y_lo);

        // Bin triangles over (x, z') for ray casting.
        let index = TriIndex::build(&tris, x_min, x_max, draft);

        // Centerplane: probe mid-depth wetted vertices for intersection counts.
        let mut probe_counts = Vec::new();
        let mut probe_mids = Vec::new();
        {
            let targets: Vec<(f64, f64)> = tris
                .iter()
                .flat_map(|t| t.iter())
                .filter(|v| v[2] >= 0.3 * draft && v[2] <= 0.7 * draft)
                .map(|v| (v[0], v[2]))
                .collect();
            let step = (targets.len() / 48).max(1);
            for t in targets.iter().step_by(step) {
                let ys = index.intersections(&tris, t.0, t.1, scale);
                if !ys.is_empty() {
                    probe_counts.push(ys.len());
                    if ys.len() >= 2 {
                        probe_mids.push((ys[0] + ys[ys.len() - 1]) / 2.0);
                    }
                }
            }
        }
        probe_counts.sort_unstable();
        let two_sided = !probe_counts.is_empty()
            && probe_counts[probe_counts.len() / 2] >= 2
            && !probe_mids.is_empty();
        let y_c = opts.centerplane.unwrap_or(if two_sided {
            probe_mids.iter().sum::<f64>() / probe_mids.len() as f64
        } else {
            0.0
        });
        let mirrored = !two_sided && y_sum / wet_count as f64 <= y_c;
        if !two_sided {
            let nearest = tris
                .iter()
                .flat_map(|t| t.iter())
                .filter(|v| v[2] >= -1e-12)
                .fold(f64::INFINITY, |m, v| m.min((v[1] - y_c).abs()));
            if nearest > 0.2 * (y_hi - y_lo).max(1e-12) {
                return Err(Error::InvalidGeometry(format!(
                    "the mesh is one-sided but never approaches the centerplane \
                     y = {y_c}; if this is an offset hull, supply the centerplane \
                     position explicitly"
                )));
            }
        }

        // Sample grid; empty ray results are the footprint.
        let ns = opts.stations;
        let nw = opts.waterlines;
        let stations: Vec<f64> = (0..ns)
            .map(|i| {
                let c = (std::f64::consts::PI * i as f64 / (ns - 1) as f64).cos();
                x_min + length * (1.0 - c) / 2.0
            })
            .collect();
        let waterlines: Vec<f64> = (0..nw)
            .map(|j| draft * j as f64 / (nw - 1) as f64)
            .collect();
        let mut grid = vec![0.0f64; ns * nw];
        let mut ambiguous = 0usize;
        let mut max_asym = 0.0f64;
        for (i, &xi) in stations.iter().enumerate() {
            for (j, &zj) in waterlines.iter().enumerate() {
                let ys = index.intersections(&tris, xi, zj, scale);
                if ys.is_empty() {
                    continue;
                }
                if ys.len() > 2 {
                    ambiguous += 1;
                }
                let folded = ys.iter().fold(0.0f64, |m, y| m.max((y - y_c).abs()));
                if two_sided && ys.len() >= 2 {
                    let stb = (ys[ys.len() - 1] - y_c).abs();
                    let prt = (ys[0] - y_c).abs();
                    max_asym = max_asym.max((stb - prt).abs());
                }
                grid[i * nw + j] = folded;
            }
            progress((i + 1) as f32 / ns as f32);
        }

        // Value-only grid: a tessellated mesh carries no usable slopes.
        let sample_grid = SampleGrid::new(stations, waterlines, grid)?;
        let (hull, fit_report) = fit_grid(&sample_grid, &opts.fit)?;
        Ok(Some(ImportedHull {
            placement: Placement { x: 0.0, y: y_c },
            hull,
            report: ImportReport {
                units_scale: self.units_scale,
                patches: src.len(),
                two_sided,
                centerplane: y_c,
                mirrored,
                draft,
                x_range: (x_min, x_max),
                max_asymmetry: max_asym,
                ambiguous_samples: ambiguous,
                failed_inversions: 0,
                derivative_gaps: 0,
                fit: fit_report,
            },
            grid: sample_grid,
        }))
    }
}

/// Uniform (x, z') bin grid over the wetted region for ray casting.
struct TriIndex {
    x0: f64,
    z0: f64,
    dx: f64,
    dz: f64,
    nx: usize,
    nz: usize,
    bins: Vec<Vec<u32>>,
}

impl TriIndex {
    fn build(tris: &[Tri], x_min: f64, x_max: f64, draft: f64) -> TriIndex {
        let (nx, nz) = (96usize, 48usize);
        let dx = ((x_max - x_min) / nx as f64).max(1e-12);
        let dz = (draft / nz as f64).max(1e-12);
        let mut idx = TriIndex {
            x0: x_min,
            z0: 0.0,
            dx,
            dz,
            nx,
            nz,
            bins: vec![Vec::new(); nx * nz],
        };
        for (ti, t) in tris.iter().enumerate() {
            let (mut xl, mut xh) = (f64::INFINITY, f64::NEG_INFINITY);
            let (mut zl, mut zh) = (f64::INFINITY, f64::NEG_INFINITY);
            for v in t {
                xl = xl.min(v[0]);
                xh = xh.max(v[0]);
                zl = zl.min(v[2]);
                zh = zh.max(v[2]);
            }
            if zh < -1e-12 {
                continue; // entirely above water
            }
            let (i0, i1) = (idx.clamp_x(xl), idx.clamp_x(xh));
            let (j0, j1) = (idx.clamp_z(zl), idx.clamp_z(zh));
            for i in i0..=i1 {
                for j in j0..=j1 {
                    idx.bins[i * idx.nz + j].push(ti as u32);
                }
            }
        }
        idx
    }

    fn clamp_x(&self, x: f64) -> usize {
        (((x - self.x0) / self.dx).floor().max(0.0) as usize).min(self.nx - 1)
    }

    fn clamp_z(&self, z: f64) -> usize {
        (((z - self.z0) / self.dz).floor().max(0.0) as usize).min(self.nz - 1)
    }

    /// All y where the transverse line through (x, z') crosses the mesh,
    /// deduped and sorted.
    fn intersections(&self, tris: &[Tri], x: f64, z: f64, scale: f64) -> Vec<f64> {
        let mut ys: Vec<f64> = Vec::new();
        let bin = &self.bins[self.clamp_x(x) * self.nz + self.clamp_z(z)];
        for &ti in bin {
            let t = &tris[ti as usize];
            let (p0, p1, p2) = (t[0], t[1], t[2]);
            let (e1x, e1z) = (p1[0] - p0[0], p1[2] - p0[2]);
            let (e2x, e2z) = (p2[0] - p0[0], p2[2] - p0[2]);
            let det = e1x * e2z - e2x * e1z;
            if det.abs() < 1e-14 * scale * scale {
                continue; // edge-on to the ray; neighbours cover it
            }
            let (rx, rz) = (x - p0[0], z - p0[2]);
            let u = (rx * e2z - e2x * rz) / det;
            let v = (e1x * rz - rx * e1z) / det;
            if u >= -1e-9 && v >= -1e-9 && u + v <= 1.0 + 1e-9 {
                ys.push(p0[1] + u * (p1[1] - p0[1]) + v * (p2[1] - p0[1]));
            }
        }
        ys.sort_by(f64::total_cmp);
        ys.dedup_by(|a, b| (*a - *b).abs() <= 1e-6 * scale);
        ys
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quad_tris(p: [[f64; 3]; 4]) -> Vec<Tri> {
        vec![[p[0], p[1], p[2]], [p[0], p[2], p[3]]]
    }

    #[test]
    fn binary_roundtrip() {
        // A single triangle, millimetre units.
        let tri = [[0.0f32, 0.0, 0.0], [1000.0, 0.0, 0.0], [0.0, 1000.0, 500.0]];
        let mut bytes = vec![0u8; 80];
        bytes.extend_from_slice(&1u32.to_le_bytes());
        bytes.extend_from_slice(&[0u8; 12]); // normal
        for v in tri {
            for c in v {
                bytes.extend_from_slice(&c.to_le_bytes());
            }
        }
        bytes.extend_from_slice(&0u16.to_le_bytes());
        let tris = parse_stl(&bytes, 0.001).unwrap();
        assert_eq!(tris.len(), 1);
        assert!((tris[0][1][0] - 1.0).abs() < 1e-12);
        assert!((tris[0][2][2] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn ascii_parse() {
        let text = "solid t\n facet normal 0 0 0\n  outer loop\n   vertex 0 0 0\n   \
                    vertex 1 0 0\n   vertex 0 1 0\n  endloop\n endfacet\nendsolid t\n";
        let tris = parse_stl(text.as_bytes(), 1.0).unwrap();
        assert_eq!(tris.len(), 1);
        assert_eq!(tris[0][1][0], 1.0);
    }

    #[test]
    fn ray_casting_finds_both_walls() {
        // Two vertical walls at y = ±1 spanning x 0..10, z' handled via a
        // fleet situated at waterline 2 (walls span z_cad 0..2).
        let mut tris = quad_tris([
            [0.0, 1.0, 0.0],
            [10.0, 1.0, 0.0],
            [10.0, 1.0, 2.0],
            [0.0, 1.0, 2.0],
        ]);
        tris.extend(quad_tris([
            [0.0, -1.0, 0.0],
            [10.0, -1.0, 0.0],
            [10.0, -1.0, 2.0],
            [0.0, -1.0, 2.0],
        ]));
        // Hull frame: z' = 2 - z.
        let hf: Vec<Tri> = tris
            .iter()
            .map(|t| {
                let mut t = *t;
                for v in t.iter_mut() {
                    v[2] = 2.0 - v[2];
                }
                t
            })
            .collect();
        let index = TriIndex::build(&hf, 0.0, 10.0, 2.0);
        let ys = index.intersections(&hf, 5.0, 1.0, 10.0);
        assert_eq!(ys.len(), 2, "{ys:?}");
        assert!((ys[0] + 1.0).abs() < 1e-12 && (ys[1] - 1.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod sectional_tests {
    use super::*;
    use crate::conditions::Conditions;
    use crate::michell::{TransomClosure, WaveOptions};
    use crate::sectional::{wave_resistance, DepthQuadrature, SectionalHull};

    /// Binary STL of a triangle list (f32, as STL stores it).
    pub(super) fn stl_bytes(tris: &[Tri]) -> Vec<u8> {
        let mut bytes = vec![0u8; 80];
        bytes.extend_from_slice(&(tris.len() as u32).to_le_bytes());
        for t in tris {
            bytes.extend_from_slice(&[0u8; 12]);
            for v in t {
                for c in v {
                    bytes.extend_from_slice(&(*c as f32).to_le_bytes());
                }
            }
            bytes.extend_from_slice(&0u16.to_le_bytes());
        }
        bytes
    }

    /// Every hull of an IGES fleet, tessellated and posed into the water
    /// frame (z up from `wl`), as one triangle list.
    pub(super) fn fleet_tris(fleet: &crate::iges::SourceFleet, wl: f64) -> Vec<Tri> {
        let mut tris: Vec<Tri> = Vec::new();
        for i in 0..fleet.len() {
            let (v, t) = fleet
                .posed_tessellation(i, wl, &HullPose::default(), &Platform::default())
                .unwrap();
            tris.extend(t.iter().map(|k| k.map(|j| v[j as usize])));
        }
        tris
    }

    /// The Wigley hull `y = ±(B/2)(1 − (2x/L)²)(1 − (z/T)²)` below the
    /// waterline (CAD z = 0), carried on wall-sided to `top` above it, as a
    /// two-sided mesh: `nx` columns along the length, `nz` rows down the
    /// draft.
    fn wigley_tris(l: f64, b: f64, t: f64, top: f64, nx: usize, nz: usize) -> Vec<Tri> {
        let half = |x: f64, depth: f64| {
            let xi = 2.0 * x / l;
            let zeta = (depth / t).clamp(0.0, 1.0);
            0.5 * b * (1.0 - xi * xi) * (1.0 - zeta * zeta)
        };
        // Depth rows: `top` above the water, then the waterline to the keel.
        let depths: Vec<f64> = std::iter::once(-top)
            .chain((0..=nz).map(|j| t * j as f64 / nz as f64))
            .collect();
        let xs: Vec<f64> = (0..=nx)
            .map(|i| -0.5 * l + l * i as f64 / nx as f64)
            .collect();
        let mut tris = Vec::new();
        for side in [1.0, -1.0] {
            let p = |i: usize, j: usize| {
                let (x, d) = (xs[i], depths[j]);
                [x, side * half(x, d), -d]
            };
            for i in 0..nx {
                for j in 0..depths.len() - 1 {
                    tris.push([p(i, j), p(i + 1, j), p(i + 1, j + 1)]);
                    tris.push([p(i, j), p(i + 1, j + 1), p(i, j + 1)]);
                }
            }
        }
        tris
    }

    fn untransomed() -> WaveOptions {
        WaveOptions {
            transom: TransomClosure::None,
            ..WaveOptions::default()
        }
    }

    fn rel(a: f64, b: f64) -> f64 {
        (a - b).abs() / b.abs()
    }

    /// A finely tessellated Wigley, imported by sections straight from its
    /// facets, against the exact hull's sections: volume and R_w agree to
    /// the tessellation's chord error.
    #[test]
    fn stl_sections_of_a_wigley_match_the_exact_hull() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let hull = crate::hulls::wigley(l, b, t).unwrap();
        let exact = SectionalHull::from_hull(&hull, &DepthQuadrature::default()).unwrap();
        let bytes = stl_bytes(&wigley_tris(l, b, t, 0.2, 400, 64));
        let fleet = import_sectional(&bytes, 1.0, &SectionalOptions::default()).unwrap();
        assert_eq!(fleet.hulls.len(), 1, "{:?}", fleet.failed);
        let imp = &fleet.hulls[0];
        eprintln!("{:?}", imp.report);
        assert!(imp.report.two_sided && imp.report.centerplane.abs() < 1e-9);
        assert!(imp.report.transom.is_none());
        assert!((imp.report.x_range.0 + 0.5 * l).abs() < 1e-6);
        assert!((imp.report.x_range.1 - 0.5 * l).abs() < 1e-6);
        let dv = rel(imp.hull.displaced_volume(), hull.displaced_volume());
        eprintln!("volume {dv:.2e}");
        assert!(dv < 2e-4, "volume {dv:.2e}");
        let wave = untransomed();
        for fn_ in [0.25, 0.3, 0.35, 0.5] {
            let cond = Conditions::seawater(fn_ * (9.81 * l).sqrt());
            let a = wave_resistance(&exact, &cond, &wave).unwrap().resistance;
            let s = wave_resistance(&imp.hull, &cond, &wave).unwrap().resistance;
            let e = rel(s, a);
            eprintln!("Fn {fn_}: Rw exact {a:.6e}, stl {s:.6e} ({e:.2e})");
            assert!(e < 5e-4, "Fn {fn_}: {e:.2e}");
        }
    }

    /// The same Wigley posed (design dz and trim, platform sinkage and
    /// trim) from its STL and from its IGES patches: the pose moves both
    /// alike, and the sections agree to the chord error.
    #[test]
    fn a_posed_stl_wigley_matches_the_posed_iges_one() {
        let (l, b, t) = (10.0, 1.0, 0.625);
        let hull = crate::hulls::wigley(l, b, t).unwrap();
        let surfs = crate::iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
        let text = crate::iges::write(&surfs, "wigley").unwrap();
        let iges = crate::iges::source_fleet(&text, 0.0).unwrap();
        let bytes = stl_bytes(&wigley_tris(l, b, t, 0.2, 400, 64));
        let mesh = mesh_fleet(&bytes, 1.0, 0.0).unwrap();
        // Raised overall, so the IGES shell's top rim (at the design
        // waterline) stays dry: the two wetted shapes are the same.
        let pose = HullPose {
            dz: -0.1,
            trim: 0.4f64.to_radians(),
            ..Default::default()
        };
        let plat = Platform {
            sinkage: 0.02,
            trim: -0.2f64.to_radians(),
            pivot_x: 1.0,
        };
        let opts = SectionalOptions::default();
        let a = iges
            .situate_sectional(0, 0.0, &pose, &plat, &opts)
            .unwrap()
            .unwrap();
        let s = mesh
            .situate_sectional(0, 0.0, &pose, &plat, &opts)
            .unwrap()
            .unwrap();
        eprintln!("iges {:?}\nstl  {:?}", a.report, s.report);
        let dv = rel(s.hull.displaced_volume(), a.hull.displaced_volume());
        eprintln!(
            "volume {:.6} vs {:.6} ({dv:.2e}), lcb {:.5} vs {:.5}",
            s.hull.displaced_volume(),
            a.hull.displaced_volume(),
            s.hull.lcb_x(),
            a.hull.lcb_x()
        );
        assert!(dv < 2e-4, "volume {dv:.2e}");
        assert!((s.hull.lcb_x() - a.hull.lcb_x()).abs() < 1e-4 * l);
        assert!((s.report.draft - a.report.draft).abs() < 1e-6);
        // The pose really moved the hull (not a vacuous comparison).
        assert!(a.report.draft < t - 0.05);
        let wave = untransomed();
        for fn_ in [0.3, 0.5] {
            let cond = Conditions::seawater(fn_ * (9.81 * l).sqrt());
            let ra = wave_resistance(&a.hull, &cond, &wave).unwrap().resistance;
            let rs = wave_resistance(&s.hull, &cond, &wave).unwrap().resistance;
            let e = rel(rs, ra);
            eprintln!("Fn {fn_}: Rw iges {ra:.6e}, stl {rs:.6e} ({e:.2e})");
            assert!(e < 5e-4, "Fn {fn_}: {e:.2e}");
        }
        // The posed tessellations share their frame: the same keel depth.
        let (vi, _) = iges.posed_tessellation(0, 0.0, &pose, &plat).unwrap();
        let (vs, _) = mesh.posed_tessellation(0, 0.0, &pose, &plat).unwrap();
        let low = |v: &[[f64; 3]]| v.iter().fold(f64::INFINITY, |m, p| m.min(p[2]));
        assert!(
            (low(&vi) - low(&vs)).abs() < 1e-5,
            "{} {}",
            low(&vi),
            low(&vs)
        );
    }

    /// Real CAD through STL: e12's IGES tessellation written as binary STL
    /// (every hull, in the water frame) and imported by sections, against
    /// the IGES sectional import. They differ by the tessellation's chord
    /// error.
    #[test]
    fn e12_through_stl_matches_its_iges_import() {
        let Some(text) = crate::cad_fixture("e12.igs") else {
            return;
        };
        let wl = -0.95;
        let iges = crate::iges::source_fleet(&text, wl).unwrap();
        let tris = fleet_tris(&iges, wl);
        let bytes = stl_bytes(&tris);
        let so = SectionalOptions {
            waterline_z: wl,
            ..Default::default()
        };
        let a = crate::iges::import_sectional(&text, &so).unwrap();
        let s = import_sectional(
            &bytes,
            1.0,
            &SectionalOptions {
                waterline_z: 0.0,
                ..so
            },
        )
        .unwrap();
        eprintln!(
            "{} triangles; iges {} hulls ({} failed), stl {} hulls ({} failed)",
            tris.len(),
            a.hulls.len(),
            a.failed.len(),
            s.hulls.len(),
            s.failed.len()
        );
        assert_eq!(
            a.hulls.len(),
            s.hulls.len(),
            "{:?} / {:?}",
            a.failed,
            s.failed
        );
        let l = a.hulls.iter().map(|h| h.hull.length()).fold(0.0, f64::max);
        let cond = Conditions::seawater(0.3 * (9.81 * l).sqrt());
        let wave = WaveOptions::default();
        for (ha, hs) in a.hulls.iter().zip(&s.hulls) {
            let dv = rel(hs.hull.displaced_volume(), ha.hull.displaced_volume());
            let ra = wave_resistance(&ha.hull, &cond, &wave).unwrap().resistance;
            let rs = wave_resistance(&hs.hull, &cond, &wave).unwrap().resistance;
            let er = rel(rs, ra);
            let ta = ha.report.transom.as_ref().map_or(0.0, |t| t.area);
            let ts = hs.report.transom.as_ref().map_or(0.0, |t| t.area);
            eprintln!(
                "y {:.4}/{:.4}: volume {:.5}/{:.5} ({dv:.2e}), x {:.4}..{:.4} / {:.4}..{:.4}, \
                 transom {ta:.5}/{ts:.5}, Rw(Fn 0.3) {ra:.5e}/{rs:.5e} ({er:.2e}), \
                 ambiguous {}/{}",
                ha.report.centerplane,
                hs.report.centerplane,
                ha.hull.displaced_volume(),
                hs.hull.displaced_volume(),
                ha.report.x_range.0,
                ha.report.x_range.1,
                hs.report.x_range.0,
                hs.report.x_range.1,
                ha.report.ambiguous_rays,
                hs.report.ambiguous_rays,
            );
            assert_eq!(ha.report.transom.is_some(), hs.report.transom.is_some());
            assert!(dv < 1e-3, "volume {dv:.2e}");
            assert!(er < 3e-3, "Rw {er:.2e}");
        }
    }
}

#[cfg(test)]
mod pose_timing {
    use super::*;

    /// How long a re-pose of an STL takes (pose the vertices, bucket, frame
    /// and ends, sections), on e12's tessellation as STL, warm and cold.
    /// A report, not a check.
    #[test]
    #[ignore = "timing report"]
    fn stl_sectional_repose_cost() {
        let Some(text) = crate::cad_fixture("e12.igs") else {
            return;
        };
        let wl = -0.95;
        let iges = crate::iges::source_fleet(&text, wl).unwrap();
        let tris = super::sectional_tests::fleet_tris(&iges, wl);
        let bytes = super::sectional_tests::stl_bytes(&tris);
        let t = std::time::Instant::now();
        let fleet = mesh_fleet(&bytes, 1.0, 0.0).unwrap();
        eprintln!(
            "load (parse, cluster, index) {:.1} ms, {} triangles",
            1e3 * t.elapsed().as_secs_f64(),
            tris.len()
        );
        let idx = (0..fleet.len())
            .max_by_key(|&i| fleet.hulls[i].len())
            .unwrap();
        eprintln!("hull {idx}: {} triangles", fleet.hulls[idx].len());
        let opts = SectionalOptions::default();
        let mut state = SectionalState::default();
        for (sink, trim_deg) in [
            (0.0, 0.0),
            (0.005, 0.0),
            (0.01, 0.0),
            (0.01, 0.1),
            (0.01, 0.2),
            (0.012, 0.3),
        ] {
            let plat = Platform {
                sinkage: sink,
                trim: f64::to_radians(trim_deg),
                pivot_x: 0.0,
            };
            let t = std::time::Instant::now();
            let warm = fleet
                .situate_sectional_warm(idx, 0.0, &HullPose::default(), &plat, &opts, &mut state)
                .unwrap()
                .unwrap();
            let t_warm = t.elapsed().as_secs_f64();
            let t = std::time::Instant::now();
            let cold = fleet
                .situate_sectional(idx, 0.0, &HullPose::default(), &plat, &opts)
                .unwrap()
                .unwrap();
            let t_cold = t.elapsed().as_secs_f64();
            let dv = (warm.hull.displaced_volume() / cold.hull.displaced_volume() - 1.0).abs();
            eprintln!(
                "sinkage {sink} trim {trim_deg}°: warm {:.1} ms, cold {:.1} ms, \
                 volume {:.5} (warm vs cold {dv:.1e})",
                1e3 * t_warm,
                1e3 * t_cold,
                warm.hull.displaced_volume()
            );
        }
    }
}

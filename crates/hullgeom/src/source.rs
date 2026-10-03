//! Geometry that sectional hulls are cut from: a parsed, clustered file kept
//! in its own (CAD) frame so each hull can be re-posed and re-cut as often as
//! a sweep or an equilibrium solve needs.
//!
//! [`crate::iges::SourceFleet`] (IGES patches) and [`crate::stl::MeshFleet`]
//! (a triangle mesh) are the two; either plugs into the solver and the CLI through [`HullSource`].

use crate::error::Result;
use crate::iges::{
    HullPose, Platform, SectionalImport, SectionalOptions, SectionalState, SourceFleet,
};
use crate::stl::MeshFleet;

/// A file's hulls, re-cuttable into sections at any pose.
pub trait HullSource: Send + Sync {
    /// Number of hulls in the file.
    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The hull's x mid in the file frame: the default pivot of its design
    /// trim ([`HullPose::pivot_x`]).
    fn x_mid(&self, idx: usize) -> f64;

    /// Cut hull `idx` into sections at the design pose and platform state,
    /// warm-started from (and updating) `state`; `Ok(None)` when it is dry.
    fn situate_sectional_warm(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &SectionalOptions,
        state: &mut SectionalState,
    ) -> Result<Option<SectionalImport>>;

    /// [`HullSource::situate_sectional_warm`] from a cold start.
    fn situate_sectional(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &SectionalOptions,
    ) -> Result<Option<SectionalImport>> {
        let mut state = SectionalState::default();
        self.situate_sectional_warm(idx, waterline_z, pose, platform, opts, &mut state)
    }

    /// The hull's surface patches at a pose, in the CAD frame (z up, metres)
    /// with the water back at `waterline_z` — for exporting a studied
    /// configuration. `None` for a source without patches (a mesh).
    fn posed_surfaces(
        &self,
        _idx: usize,
        _waterline_z: f64,
        _pose: &HullPose,
        _platform: &Platform,
    ) -> Option<Result<Vec<crate::iges::NurbsSurface3>>> {
        None
    }

    /// The whole hull (above water too) as triangles at a pose, for display:
    /// `x` forward, `y` transverse, `z` up from the effective waterline.
    fn posed_tessellation(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>)>;
}

impl HullSource for SourceFleet {
    fn len(&self) -> usize {
        SourceFleet::len(self)
    }

    fn posed_surfaces(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> Option<Result<Vec<crate::iges::NurbsSurface3>>> {
        Some(SourceFleet::posed_surfaces(
            self,
            idx,
            waterline_z,
            pose,
            platform,
        ))
    }

    fn x_mid(&self, idx: usize) -> f64 {
        SourceFleet::x_mid(self, idx)
    }

    fn situate_sectional_warm(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &SectionalOptions,
        state: &mut SectionalState,
    ) -> Result<Option<SectionalImport>> {
        SourceFleet::situate_sectional_warm(self, idx, waterline_z, pose, platform, opts, state)
    }

    fn posed_tessellation(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>)> {
        SourceFleet::posed_tessellation(self, idx, waterline_z, pose, platform)
    }
}

impl HullSource for MeshFleet {
    fn len(&self) -> usize {
        MeshFleet::len(self)
    }

    fn x_mid(&self, idx: usize) -> f64 {
        MeshFleet::x_mid(self, idx)
    }

    fn situate_sectional_warm(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
        opts: &SectionalOptions,
        state: &mut SectionalState,
    ) -> Result<Option<SectionalImport>> {
        MeshFleet::situate_sectional_warm(self, idx, waterline_z, pose, platform, opts, state)
    }

    fn posed_tessellation(
        &self,
        idx: usize,
        waterline_z: f64,
        pose: &HullPose,
        platform: &Platform,
    ) -> Result<(Vec<[f64; 3]>, Vec<[u32; 3]>)> {
        MeshFleet::posed_tessellation(self, idx, waterline_z, pose, platform)
    }
}

/// One hull of a fleet under study: which file, which hull in it, where its
/// design waterline is, and how it is mounted.
#[derive(Clone, Copy)]
pub struct SourceHull<'a> {
    pub source: &'a dyn HullSource,
    pub index: usize,
    /// CAD height of the design waterline in the source's frame.
    pub waterline_z: f64,
    pub pose: HullPose,
}

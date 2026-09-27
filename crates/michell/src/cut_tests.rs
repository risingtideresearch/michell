//! Geometry checked through the wave kernel: hulls cut from IGES, STL and
//! B-spline surfaces must give the same resistance and squat as the exact
//! oracle, posed or not. These live here rather than in `michell-geometry`
//! because the resistance is the sharpest measure of a cut's fidelity.

mod iges_ends {
    use michell_geometry::iges::*;

    fn e12() -> Option<SourceFleet> {
        let text = crate::cad_fixture("e12.igs")?;
        Some(source_fleet(&text, -0.95).unwrap())
    }

    /// Steeply trimmed at speed (Fn 0.68 floats e12 about 1.5° bow up), the
    /// leaning transom leaves only a wedge of hull aft of its full section,
    /// and the bow's plumb stem a thin wide sliver forward. Whether an end
    /// station caught those or nothing once decided whether the hull ended
    /// on a transom — the lift jumped 6× between neighbouring trims, and
    /// the dynamic equilibrium never converged. The ends are now pulled in
    /// to the full sections, and the end steps are always in the force, so
    /// the lift moves smoothly with trim.
    #[test]
    fn a_steeply_trimmed_cut_keeps_its_ends_and_a_smooth_lift() {
        let Some(src) = e12() else {
            return;
        };
        let opts = SectionalOptions {
            waterline_z: -0.95,
            ..Default::default()
        };
        let design = src
            .situate_sectional(0, -0.95, &HullPose::default(), &Platform::default(), &opts)
            .unwrap()
            .unwrap()
            .hull;
        let cond =
            crate::Conditions::seawater(0.68 * (9.81f64 * design.length()).sqrt());
        let mut last: Option<f64> = None;
        for k in 0..=5 {
            let trim = (1.0 + 0.2 * k as f64).to_radians();
            let pf = Platform {
                sinkage: 0.0,
                trim,
                pivot_x: design.lcb_x(),
            };
            let h = src
                .situate_sectional(0, -0.95, &HullPose::default(), &pf, &opts)
                .unwrap()
                .unwrap();
            assert!(h.hull.transom().is_some(), "trim {k}: transom lost");
            let f = crate::sectional::multihull_dynamic_force(
                &[(&h.hull, h.placement)],
                &cond,
                design.lcb_x(),
                &Default::default(),
            )
            .unwrap()
            .force_up;
            if let Some(prev) = last {
                assert!((f / prev - 1.0).abs() < 0.1, "lift {prev} -> {f} at step {k}");
            }
            last = Some(f);
        }
    }

    /// e12's transom is a flat face its side skins run on past. Trimmed bow
    /// up, the face tilts across the aft station plane; read as section
    /// boundary it collapsed the end section to a sliver, closing the hull
    /// over one span (the near-field force doubled from trim 0 to +0.02°).
    /// Faces are end caps, not shell: the aft station keeps the full
    /// transom section at every trim.
    #[test]
    fn a_trimmed_transom_keeps_its_section() {
        let Some(fleet) = e12() else {
            return;
        };
        let idx = (0..fleet.len())
            .max_by_key(|&i| fleet.patch_count(i))
            .unwrap();
        let so = SectionalOptions {
            waterline_z: -0.95,
            stations: 61,
            rays: 17,
            ..Default::default()
        };
        let mut last_vol = 0.0;
        for trim_deg in [-0.1f64, 0.0, 0.02, 0.05, 0.1] {
            let plat = Platform {
                sinkage: 0.012,
                trim: trim_deg.to_radians(),
                pivot_x: 4.5,
            };
            let h = fleet
                .situate_sectional(idx, -0.95, &HullPose::default(), &plat, &so)
                .unwrap()
                .unwrap();
            let (st, _) = h.hull.depth_integral_curve(0.0, 1);
            assert!(
                st[0].1 > 0.9 * st[1].1,
                "trim {trim_deg}°: end section {:.3e} against its neighbour's {:.3e}",
                st[0].1,
                st[1].1
            );
            let vol = h.hull.displaced_volume();
            assert!(vol > last_vol, "bow-up trim sinks this stern deeper: {vol}");
            last_vol = vol;
        }
    }}

mod stl_sections {
    use crate::hull::FromHull;
    use crate::michell::{TransomClosure, WaveOptions};
    use crate::sectional::wave_resistance;
    use crate::Conditions;
    use michell_geometry::iges::{HullPose, Platform, SectionalOptions};
    use michell_geometry::sectional::{DepthQuadrature, SectionalHull};
    use michell_geometry::stl::*;

    /// Binary STL of a triangle list (f32, as STL stores it).
    fn stl_bytes(tris: &[Tri]) -> Vec<u8> {
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
    fn fleet_tris(fleet: &michell_geometry::iges::SourceFleet, wl: f64) -> Vec<Tri> {
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
        let surfs = michell_geometry::iges::halfbreadth_surfaces(hull.surface(), 0.0, 0.0);
        let text = michell_geometry::iges::write(&surfs, "wigley").unwrap();
        let iges = michell_geometry::iges::source_fleet(&text, 0.0).unwrap();
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
        let iges = michell_geometry::iges::source_fleet(&text, wl).unwrap();
        let tris = fleet_tris(&iges, wl);
        let bytes = stl_bytes(&tris);
        let so = SectionalOptions {
            waterline_z: wl,
            ..Default::default()
        };
        let a = michell_geometry::iges::import_sectional(&text, &so).unwrap();
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
    }}

mod source_recut {
    use crate::hull::FromHull;
    use crate::sectional::wave_resistance;
    use crate::{Conditions, WaveOptions};
    use michell_geometry::iges::{HullPose, Platform, SectionalOptions, SourceFleet};
    use michell_geometry::sectional::SectionalHull;

    /// A B-spline hull re-cut from its exact surfaces agrees with the
    /// direct conversion, and re-poses: 5 cm deeper displaces more.
    #[test]
    fn a_spline_hull_recuts_from_its_surfaces() {
        let hull = crate::hulls::wigley(10.0, 1.0, 0.625).unwrap();
        let direct = SectionalHull::from_hull(&hull, &Default::default()).unwrap();
        let src = SourceFleet::from_halfbreadth(hull.surface(), 0.0, 0.0).unwrap();
        let opts = SectionalOptions::default();
        let cut = src
            .situate_sectional(0, 0.0, &HullPose::default(), &Platform::default(), &opts)
            .unwrap()
            .unwrap();
        let cond = Conditions::seawater(0.35 * (9.81f64 * 10.0).sqrt());
        let w = WaveOptions::default();
        let (a, b) = (
            wave_resistance(&direct, &cond, &w).unwrap().resistance,
            wave_resistance(&cut.hull, &cond, &w).unwrap().resistance,
        );
        let vol = (direct.displaced_volume(), cut.hull.displaced_volume());
        eprintln!(
            "Rw direct {a} cut {b}; volume {vol:?}; y {}",
            cut.placement.y
        );
        assert!((a - b).abs() < 1e-3 * a, "Rw {a} vs {b}");
        assert!((vol.0 - vol.1).abs() < 1e-5 * vol.0, "{vol:?}");
        let deeper = HullPose {
            dz: 0.05,
            ..Default::default()
        };
        let d = src
            .situate_sectional(0, 0.0, &deeper, &Platform::default(), &opts)
            .unwrap()
            .unwrap();
        assert!(d.hull.displaced_volume() > vol.0);
    }}

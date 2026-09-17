//! A transverse **section cut** through the raw, unclustered model at midship,
//! shown in the import dialog so the design waterline and band can be placed
//! against the real section shape. Clustering into hulls happens later, in the
//! loft — this view is deliberately the whole model at one station, so a
//! trimaran shows its three section curves side by side.
//!
//! The cut is a plane `x = x_mid` intersected with the geometry: STL triangles
//! are sliced exactly; IGES surfaces are sampled on a (u, v) grid into quads and
//! those are sliced. It yields line segments in the transverse `(y, z)` plane —
//! the section outline(s). Coordinates are the file's own (CAD frame, z up); the
//! caller applies any STL units scale at draw time.

use michell::iges::NurbsSurface3;

/// Section segments and their bounds, in file units (pre-scale).
pub struct Preview {
    /// Section outline as line segments: `[[y0, z0], [y1, z1]]`.
    pub segments: Vec<[[f32; 2]; 2]>,
    pub y_min: f32,
    pub y_max: f32,
    pub z_min: f32,
    pub z_max: f32,
}

/// (u, v) samples per IGES surface when slicing.
const GRID: usize = 48;

impl Preview {
    /// Load a midship section from a CAD/mesh file. `is_stl` selects the parser.
    pub fn load(path: &std::path::Path, is_stl: bool) -> Result<Preview, String> {
        let bytes =
            std::fs::read(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        let tris: Vec<[[f64; 3]; 3]> = if is_stl {
            michell::stl::parse_stl(&bytes, 1.0)
                .map_err(|e| e.to_string())?
                .to_vec()
        } else {
            let text =
                std::str::from_utf8(&bytes).map_err(|_| "IGES file is not UTF-8".to_string())?;
            let file = michell::iges::parse(text).map_err(|e| e.to_string())?;
            tessellate(&file.surfaces)
        };
        if tris.is_empty() {
            return Err("no geometry found in the file".into());
        }
        Self::section(&tris)
    }

    /// Cut `tris` with the plane `x = x_mid` and collect the (y, z) segments.
    fn section(tris: &[[[f64; 3]; 3]]) -> Result<Preview, String> {
        let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for t in tris {
            for v in t {
                x_min = x_min.min(v[0]);
                x_max = x_max.max(v[0]);
            }
        }
        let x_mid = 0.5 * (x_min + x_max);

        let mut segments = Vec::new();
        let (mut y_min, mut y_max, mut z_min, mut z_max) = (
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::INFINITY,
            f32::NEG_INFINITY,
        );
        for t in tris {
            if let Some(seg) = slice_tri(t, x_mid) {
                for p in &seg {
                    y_min = y_min.min(p[0]);
                    y_max = y_max.max(p[0]);
                    z_min = z_min.min(p[1]);
                    z_max = z_max.max(p[1]);
                }
                segments.push(seg);
            }
        }
        if segments.is_empty() {
            return Err("the midship section is empty (no geometry at x_mid)".into());
        }
        Ok(Preview {
            segments,
            y_min,
            y_max,
            z_min,
            z_max,
        })
    }
}

/// Sample every IGES surface into a grid of quads (two triangles each), so the
/// section slicer can treat IGES and STL identically.
fn tessellate(surfaces: &[NurbsSurface3]) -> Vec<[[f64; 3]; 3]> {
    let mut tris = Vec::new();
    for s in surfaces {
        let (u0, u1) = s.u_domain();
        let (v0, v1) = s.v_domain();
        if !(u1 > u0 && v1 > v0) {
            continue;
        }
        // Sample points on a (GRID+1)^2 lattice over the surface's domain.
        let mut grid = Vec::with_capacity((GRID + 1) * (GRID + 1));
        for iu in 0..=GRID {
            let u = u0 + (u1 - u0) * iu as f64 / GRID as f64;
            for iv in 0..=GRID {
                let v = v0 + (v1 - v0) * iv as f64 / GRID as f64;
                grid.push(s.point(u, v));
            }
        }
        let at = |iu: usize, iv: usize| grid[iu * (GRID + 1) + iv];
        for iu in 0..GRID {
            for iv in 0..GRID {
                let (a, b, c, d) = (
                    at(iu, iv),
                    at(iu + 1, iv),
                    at(iu + 1, iv + 1),
                    at(iu, iv + 1),
                );
                tris.push([a, b, c]);
                tris.push([a, c, d]);
            }
        }
    }
    tris
}

/// Intersect a triangle with the plane `x = x0`; returns the `(y, z)` segment
/// where the triangle crosses it, if any. The `>= 0` sign test counts a shared
/// vertex once, so adjacent triangles don't double up.
fn slice_tri(tri: &[[f64; 3]; 3], x0: f64) -> Option<[[f32; 2]; 2]> {
    let mut pts = [[0.0f32; 2]; 2];
    let mut n = 0;
    for e in 0..3 {
        let a = tri[e];
        let b = tri[(e + 1) % 3];
        let (da, db) = (a[0] - x0, b[0] - x0);
        if (da < 0.0) != (db < 0.0) {
            let t = da / (da - db);
            let y = (a[1] + t * (b[1] - a[1])) as f32;
            let z = (a[2] + t * (b[2] - a[2])) as f32;
            if n < 2 {
                pts[n] = [y, z];
            }
            n += 1;
        }
    }
    if n == 2 {
        Some(pts)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slices_a_triangle_across_the_plane() {
        // Triangle spanning x = -1..1; the plane x = 0 cuts two of its edges.
        let tri = [[-1.0, 0.0, 0.0], [1.0, 2.0, 0.0], [1.0, 0.0, 3.0]];
        let seg = slice_tri(&tri, 0.0).expect("crosses x=0");
        // Both crossing points lie at x = 0: edge 0 (→ y=1,z=0) and edge 2 (→ y=0,z=1.5).
        let ys: Vec<f32> = seg.iter().map(|p| p[0]).collect();
        assert!(ys.contains(&1.0));
        assert!(seg.iter().any(|p| (p[1] - 1.5).abs() < 1e-6));
    }

    #[test]
    fn triangle_entirely_on_one_side_yields_nothing() {
        let tri = [[1.0, 0.0, 0.0], [2.0, 1.0, 0.0], [3.0, 0.0, 1.0]];
        assert!(slice_tri(&tri, 0.0).is_none());
    }

    #[test]
    fn section_of_a_prism_has_segments_and_bounds() {
        // A little triangular prism spanning x = -1..1; midship is x = 0.
        let tris = [
            [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0]],
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 2.0], [-1.0, 0.0, 2.0]],
        ];
        let p = Preview::section(&tris).unwrap();
        assert!(!p.segments.is_empty());
        assert!(p.z_max >= 1.9);
    }

    #[test]
    fn tiny_ascii_stl_sections() {
        let stl = "solid t\nfacet normal 0 0 0\nouter loop\n\
                   vertex -1 0 0\nvertex 1 2 0\nvertex 1 0 3\n\
                   endloop\nendfacet\nendsolid t\n";
        let tris = michell::stl::parse_stl(stl.as_bytes(), 1.0).unwrap();
        let p = Preview::section(&tris).unwrap();
        assert_eq!(p.segments.len(), 1);
    }
}

//! A cheap transverse (body-plan) silhouette of an IGES/STL file, used by the
//! import dialog so the design waterline and band can be placed against the
//! geometry rather than typed blind.
//!
//! It reads the raw geometry only — IGES control nets, STL vertices — projects
//! every point onto the (y, z) plane, and reduces it to a transverse *envelope*:
//! the min/max beam at each height. Drawn as a filled outline that reads as the
//! hull's forward section, far clearer than a raw point scatter. No lofting or
//! sampling happens here, so it runs the moment a file is picked. Coordinates
//! are the file's own (CAD frame, z up); the caller applies any STL units scale
//! at draw time, matching the frame the waterline is entered in.

/// The transverse envelope and its bounds, in file units (pre-scale).
pub struct Preview {
    /// Per-height slices, ascending in z: `[z, y_lo, y_hi]` — the widest port
    /// and starboard reach of the geometry at that height.
    pub slices: Vec<[f32; 3]>,
    pub y_min: f32,
    pub y_max: f32,
    pub z_min: f32,
    pub z_max: f32,
}

/// Number of height bands the envelope is sampled into.
const BINS: usize = 120;

impl Preview {
    /// Load a preview from a CAD/mesh file. `is_stl` selects the parser.
    pub fn load(path: &std::path::Path, is_stl: bool) -> Result<Preview, String> {
        let bytes =
            std::fs::read(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        let raw: Vec<[f64; 3]> = if is_stl {
            // Parse in file units (scale 1.0); the dialog applies the chosen
            // units scale when drawing.
            let tris = michell::stl::parse_stl(&bytes, 1.0).map_err(|e| e.to_string())?;
            tris.iter().flat_map(|t| t.iter().copied()).collect()
        } else {
            let text =
                std::str::from_utf8(&bytes).map_err(|_| "IGES file is not UTF-8".to_string())?;
            let file = michell::iges::parse(text).map_err(|e| e.to_string())?;
            file.surfaces
                .iter()
                .flat_map(|s| s.ctrl.iter().copied())
                .collect()
        };
        Self::from_raw(&raw)
    }

    /// Project `(x, y, z)` points onto the transverse `(y, z)` plane, then bin
    /// by height into a min/max-beam envelope.
    fn from_raw(raw: &[[f64; 3]]) -> Result<Preview, String> {
        if raw.is_empty() {
            return Err("no geometry found in the file".into());
        }
        let (mut y_min, mut y_max, mut z_min, mut z_max) = (
            f32::INFINITY,
            f32::NEG_INFINITY,
            f32::INFINITY,
            f32::NEG_INFINITY,
        );
        for p in raw {
            let (y, z) = (p[1] as f32, p[2] as f32);
            y_min = y_min.min(y);
            y_max = y_max.max(y);
            z_min = z_min.min(z);
            z_max = z_max.max(z);
        }

        // Bin points by height; track the port/starboard reach in each band.
        let span = (z_max - z_min).max(1e-6);
        let mut lo = vec![f32::INFINITY; BINS];
        let mut hi = vec![f32::NEG_INFINITY; BINS];
        for p in raw {
            let (y, z) = (p[1] as f32, p[2] as f32);
            let t = ((z - z_min) / span * BINS as f32).floor() as usize;
            let b = t.min(BINS - 1);
            lo[b] = lo[b].min(y);
            hi[b] = hi[b].max(y);
        }
        let slices: Vec<[f32; 3]> = (0..BINS)
            .filter(|&b| hi[b] >= lo[b])
            .map(|b| {
                let z = z_min + (b as f32 + 0.5) / BINS as f32 * span;
                [z, lo[b], hi[b]]
            })
            .collect();

        Ok(Preview {
            slices,
            y_min,
            y_max,
            z_min,
            z_max,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn envelope_covers_the_y_range() {
        // A little box of points spanning y in [-2, 1.5], z in [-1, 3].
        let raw = [
            [10.0, -2.0, 0.0],
            [-5.0, 1.5, 3.0],
            [0.0, 0.0, -1.0],
            [3.0, -1.0, 1.0],
        ];
        let p = Preview::from_raw(&raw).unwrap();
        assert_eq!(p.y_min, -2.0);
        assert_eq!(p.y_max, 1.5);
        assert_eq!(p.z_min, -1.0);
        assert_eq!(p.z_max, 3.0);
        assert!(!p.slices.is_empty());
        // The envelope's overall reach matches the point bounds.
        let lo = p.slices.iter().map(|s| s[1]).fold(f32::INFINITY, f32::min);
        let hi = p
            .slices
            .iter()
            .map(|s| s[2])
            .fold(f32::NEG_INFINITY, f32::max);
        assert_eq!(lo, -2.0);
        assert_eq!(hi, 1.5);
        // Slices are ascending in z.
        assert!(p.slices.windows(2).all(|w| w[0][0] <= w[1][0]));
    }

    #[test]
    fn empty_geometry_errors() {
        assert!(Preview::from_raw(&[]).is_err());
    }

    #[test]
    fn parses_a_tiny_ascii_stl() {
        let stl = "solid t\nfacet normal 0 0 0\nouter loop\n\
                   vertex 0 0 0\nvertex 1 2 0\nvertex 0 0 4\n\
                   endloop\nendfacet\nendsolid t\n";
        let tris = michell::stl::parse_stl(stl.as_bytes(), 1.0).unwrap();
        let raw: Vec<[f64; 3]> = tris.iter().flat_map(|t| t.iter().copied()).collect();
        let p = Preview::from_raw(&raw).unwrap();
        assert_eq!(p.z_max, 4.0);
        assert_eq!(p.y_max, 2.0);
        assert!(!p.slices.is_empty());
    }
}

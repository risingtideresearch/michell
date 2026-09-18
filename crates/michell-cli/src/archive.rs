//! Self-describing binary sweep archive (`.msw`).
//!
//! A single-file container bundling everything needed to reproduce and
//! post-process a study: the manifest, the referenced hull files, and, for
//! every output row, the swept parameter values, the scalar metrics, the full
//! righting-arm (GZ) curve, and the free-wave spectrum A(θ). Zero-dependency —
//! the format is a hand-rolled, little-endian, length-prefixed blob stream, so
//! a reader needs nothing but `struct.unpack`/`DataView`.
//!
//! ```text
//! Header
//!   magic    4 bytes  "MSWP"
//!   version  u32      = 1
//! Blob stream (repeats to EOF)
//!   kind     u32      1=MANIFEST 2=HULLFILE 3=META 4=ROWS
//!   name_len u32
//!   name     name_len bytes, UTF-8
//!   data_len u64
//!   data     data_len bytes
//! ```
//!
//! The `META` blob is a JSON object naming the axis and metric columns (so the
//! `ROWS` blob is interpretable without hard-coding the schema) plus study
//! scalars. The `ROWS` blob is:
//!
//! ```text
//!   n_rows    u32
//!   n_axes    u32
//!   n_metrics u32
//!   per row:
//!     f64 * n_axes                              swept parameter values
//!     f64 * n_metrics                           scalar metrics (see META)
//!     gz_n      u32
//!     (f64 heel_rad, f64 gz_m) * gz_n           righting-arm curve
//!     f64 spec_wavenumber                       ν = g/U²
//!     f64 spec_transverse_wavelength
//!     spec_n    u32
//!     (f64 theta, f64 amp_re, f64 amp_im, f64 drw_dtheta) * spec_n
//! ```
//!
//! All multi-byte integers and floats are little-endian.

pub const KIND_MANIFEST: u32 = 1;
pub const KIND_HULLFILE: u32 = 2;
pub const KIND_META: u32 = 3;
pub const KIND_ROWS: u32 = 4;

const MAGIC: &[u8; 4] = b"MSWP";
const VERSION: u32 = 1;

/// One sample of the free-wave spectrum at a propagation angle θ.
pub struct SpecSample {
    pub theta: f64,
    pub amp_re: f64,
    pub amp_im: f64,
    pub drw_dtheta: f64,
}

/// The `ROWS` blob under construction. Rows are appended in output order.
pub struct Rows {
    n_axes: u32,
    n_metrics: u32,
    n_rows: u32,
    body: Vec<u8>,
}

impl Rows {
    pub fn new(n_axes: usize, n_metrics: usize) -> Self {
        Rows {
            n_axes: n_axes as u32,
            n_metrics: n_metrics as u32,
            n_rows: 0,
            body: Vec::new(),
        }
    }

    /// Append a row. `gz_curve` is `(heel_rad, gz_m)` pairs (empty when GZ is
    /// undefined, i.e. no weight+vcg loading); `wavenumber` /
    /// `transverse_wavelength` describe the spectrum whose samples follow.
    pub fn push(
        &mut self,
        params: &[f64],
        metrics: &[f64],
        gz_curve: &[(f64, f64)],
        wavenumber: f64,
        transverse_wavelength: f64,
        spectrum: &[SpecSample],
    ) {
        debug_assert_eq!(params.len() as u32, self.n_axes);
        debug_assert_eq!(metrics.len() as u32, self.n_metrics);
        for &v in params {
            put_f64(&mut self.body, v);
        }
        for &v in metrics {
            put_f64(&mut self.body, v);
        }
        put_u32(&mut self.body, gz_curve.len() as u32);
        for &(heel, gz) in gz_curve {
            put_f64(&mut self.body, heel);
            put_f64(&mut self.body, gz);
        }
        put_f64(&mut self.body, wavenumber);
        put_f64(&mut self.body, transverse_wavelength);
        put_u32(&mut self.body, spectrum.len() as u32);
        for s in spectrum {
            put_f64(&mut self.body, s.theta);
            put_f64(&mut self.body, s.amp_re);
            put_f64(&mut self.body, s.amp_im);
            put_f64(&mut self.body, s.drw_dtheta);
        }
        self.n_rows += 1;
    }

    /// Serialize to the `ROWS` blob bytes (header + accumulated rows).
    pub fn into_blob(self) -> Vec<u8> {
        let mut out = Vec::with_capacity(12 + self.body.len());
        put_u32(&mut out, self.n_rows);
        put_u32(&mut out, self.n_axes);
        put_u32(&mut out, self.n_metrics);
        out.extend_from_slice(&self.body);
        out
    }
}

/// The archive container: an ordered list of named, typed blobs.
#[derive(Default)]
pub struct Archive {
    blobs: Vec<(u32, String, Vec<u8>)>,
}

impl Archive {
    pub fn add(&mut self, kind: u32, name: &str, data: Vec<u8>) {
        self.blobs.push((kind, name.to_string(), data));
    }

    pub fn into_bytes(self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        put_u32(&mut out, VERSION);
        for (kind, name, data) in self.blobs {
            put_u32(&mut out, kind);
            put_u32(&mut out, name.len() as u32);
            out.extend_from_slice(name.as_bytes());
            put_u64(&mut out, data.len() as u64);
            out.extend_from_slice(&data);
        }
        out
    }
}

fn put_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_f64(buf: &mut Vec<u8>, v: f64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

/// A parsed archive — used by the round-trip test and available to readers.
#[cfg(test)]
pub struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

#[cfg(test)]
impl<'a> Reader<'a> {
    pub fn new(data: &'a [u8]) -> Result<Self, String> {
        if data.len() < 8 || &data[0..4] != MAGIC {
            return Err("not an MSWP archive".into());
        }
        let version = u32::from_le_bytes(data[4..8].try_into().unwrap());
        if version != VERSION {
            return Err(format!("unsupported archive version {version}"));
        }
        Ok(Reader { data, pos: 8 })
    }

    /// Next blob as `(kind, name, data)`, or None at EOF.
    pub fn next_blob(&mut self) -> Option<(u32, String, &'a [u8])> {
        if self.pos >= self.data.len() {
            return None;
        }
        let kind = u32::from_le_bytes(self.data[self.pos..self.pos + 4].try_into().unwrap());
        self.pos += 4;
        let nl = u32::from_le_bytes(self.data[self.pos..self.pos + 4].try_into().unwrap()) as usize;
        self.pos += 4;
        let name = String::from_utf8(self.data[self.pos..self.pos + nl].to_vec()).unwrap();
        self.pos += nl;
        let dl = u64::from_le_bytes(self.data[self.pos..self.pos + 8].try_into().unwrap()) as usize;
        self.pos += 8;
        let data = &self.data[self.pos..self.pos + dl];
        self.pos += dl;
        Some((kind, name, data))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Read a little-endian f64 from `data` at `*pos`, advancing it.
    fn get_f64(data: &[u8], pos: &mut usize) -> f64 {
        let v = f64::from_le_bytes(data[*pos..*pos + 8].try_into().unwrap());
        *pos += 8;
        v
    }
    fn get_u32(data: &[u8], pos: &mut usize) -> u32 {
        let v = u32::from_le_bytes(data[*pos..*pos + 4].try_into().unwrap());
        *pos += 4;
        v
    }

    #[test]
    fn round_trips_blobs_and_rows() {
        let mut rows = Rows::new(2, 3);
        rows.push(
            &[1.0, 2.0],
            &[10.0, 20.0, 30.0],
            &[(0.0, 0.0), (0.1, 0.5)],
            0.5,
            12.5,
            &[
                SpecSample {
                    theta: -0.2,
                    amp_re: 1.0,
                    amp_im: -1.0,
                    drw_dtheta: 3.0,
                },
                SpecSample {
                    theta: 0.2,
                    amp_re: 2.0,
                    amp_im: 0.5,
                    drw_dtheta: 4.0,
                },
            ],
        );
        rows.push(&[3.0, 4.0], &[11.0, 21.0, 31.0], &[], 0.6, 10.0, &[]);

        let mut ar = Archive::default();
        ar.add(KIND_MANIFEST, "study.json", b"{\"a\":1}".to_vec());
        ar.add(KIND_HULLFILE, "h.hull", b"michell-hull v1\n".to_vec());
        ar.add(KIND_META, "meta.json", b"{}".to_vec());
        ar.add(KIND_ROWS, "rows", rows.into_blob());
        let bytes = ar.into_bytes();

        let mut r = Reader::new(&bytes).unwrap();
        let (k, name, data) = r.next_blob().unwrap();
        assert_eq!(k, KIND_MANIFEST);
        assert_eq!(name, "study.json");
        assert_eq!(data, b"{\"a\":1}");
        let (k, name, _) = r.next_blob().unwrap();
        assert_eq!((k, name.as_str()), (KIND_HULLFILE, "h.hull"));
        let (k, _, _) = r.next_blob().unwrap();
        assert_eq!(k, KIND_META);
        let (k, _, rowbytes) = r.next_blob().unwrap();
        assert_eq!(k, KIND_ROWS);
        assert!(r.next_blob().is_none());

        // Decode the ROWS blob.
        let mut p = 0usize;
        let n_rows = get_u32(rowbytes, &mut p);
        let n_axes = get_u32(rowbytes, &mut p);
        let n_metrics = get_u32(rowbytes, &mut p);
        assert_eq!((n_rows, n_axes, n_metrics), (2, 2, 3));

        // Row 0.
        assert_eq!(get_f64(rowbytes, &mut p), 1.0);
        assert_eq!(get_f64(rowbytes, &mut p), 2.0);
        for want in [10.0, 20.0, 30.0] {
            assert_eq!(get_f64(rowbytes, &mut p), want);
        }
        let gz_n = get_u32(rowbytes, &mut p);
        assert_eq!(gz_n, 2);
        for _ in 0..gz_n {
            get_f64(rowbytes, &mut p);
            get_f64(rowbytes, &mut p);
        }
        assert_eq!(get_f64(rowbytes, &mut p), 0.5); // wavenumber
        assert_eq!(get_f64(rowbytes, &mut p), 12.5); // transverse wavelength
        let spec_n = get_u32(rowbytes, &mut p);
        assert_eq!(spec_n, 2);
        for _ in 0..spec_n {
            for _ in 0..4 {
                get_f64(rowbytes, &mut p);
            }
        }

        // Row 1: empty gz curve and spectrum.
        assert_eq!(get_f64(rowbytes, &mut p), 3.0);
        assert_eq!(get_f64(rowbytes, &mut p), 4.0);
        for want in [11.0, 21.0, 31.0] {
            assert_eq!(get_f64(rowbytes, &mut p), want);
        }
        assert_eq!(get_u32(rowbytes, &mut p), 0); // gz_n
        get_f64(rowbytes, &mut p);
        get_f64(rowbytes, &mut p);
        assert_eq!(get_u32(rowbytes, &mut p), 0); // spec_n
        assert_eq!(p, rowbytes.len(), "consumed the whole ROWS blob");
    }

    #[test]
    fn rejects_bad_magic() {
        assert!(Reader::new(b"NOPE\0\0\0\0").is_err());
    }
}

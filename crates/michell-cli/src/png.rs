//! Minimal zero-dependency PNG writer (8-bit RGB, stored deflate blocks) and
//! the diverging colormap used by the wake heatmap.

/// Encode an RGB8 image (`rgb` is row-major, `3 * width * height` bytes,
/// row 0 at the top) as a PNG file image.
pub fn encode_rgb(width: usize, height: usize, rgb: &[u8]) -> Vec<u8> {
    assert_eq!(rgb.len(), 3 * width * height, "pixel buffer size");
    // Raw scanlines, each prefixed by filter byte 0 (None).
    let stride = 3 * width;
    let mut raw = Vec::with_capacity(height * (stride + 1));
    for row in 0..height {
        raw.push(0);
        raw.extend_from_slice(&rgb[row * stride..(row + 1) * stride]);
    }
    // zlib stream: header + stored (uncompressed) deflate blocks + adler32.
    let mut z = vec![0x78, 0x01];
    let mut off = 0;
    while off < raw.len() {
        let n = (raw.len() - off).min(65535);
        z.push(u8::from(off + n == raw.len()));
        z.extend_from_slice(&(n as u16).to_le_bytes());
        z.extend_from_slice(&(!(n as u16)).to_le_bytes());
        z.extend_from_slice(&raw[off..off + n]);
        off += n;
    }
    z.extend_from_slice(&adler32(&raw).to_be_bytes());

    let mut png = vec![0x89, b'P', b'N', b'G', 0x0d, 0x0a, 0x1a, 0x0a];
    let mut ihdr = Vec::with_capacity(13);
    ihdr.extend_from_slice(&(width as u32).to_be_bytes());
    ihdr.extend_from_slice(&(height as u32).to_be_bytes());
    // bit depth 8, color type 2 (truecolor), default compression/filter/interlace
    ihdr.extend_from_slice(&[8, 2, 0, 0, 0]);
    chunk(&mut png, b"IHDR", &ihdr);
    chunk(&mut png, b"IDAT", &z);
    chunk(&mut png, b"IEND", &[]);
    png
}

fn chunk(out: &mut Vec<u8>, kind: &[u8; 4], data: &[u8]) {
    out.extend_from_slice(&(data.len() as u32).to_be_bytes());
    out.extend_from_slice(kind);
    out.extend_from_slice(data);
    let mut crc = Crc32::new();
    crc.update(kind);
    crc.update(data);
    out.extend_from_slice(&crc.finish().to_be_bytes());
}

struct Crc32 {
    table: [u32; 256],
    value: u32,
}

impl Crc32 {
    fn new() -> Crc32 {
        let mut table = [0u32; 256];
        for (n, slot) in table.iter_mut().enumerate() {
            let mut c = n as u32;
            for _ in 0..8 {
                c = if c & 1 != 0 {
                    0xedb8_8320 ^ (c >> 1)
                } else {
                    c >> 1
                };
            }
            *slot = c;
        }
        Crc32 {
            table,
            value: 0xffff_ffff,
        }
    }

    fn update(&mut self, data: &[u8]) {
        for &b in data {
            self.value = self.table[((self.value ^ b as u32) & 0xff) as usize] ^ (self.value >> 8);
        }
    }

    fn finish(&self) -> u32 {
        self.value ^ 0xffff_ffff
    }
}

fn adler32(data: &[u8]) -> u32 {
    const MOD: u32 = 65521;
    let (mut a, mut b) = (1u32, 0u32);
    for chunk in data.chunks(5552) {
        for &byte in chunk {
            a += byte as u32;
            b += a;
        }
        a %= MOD;
        b %= MOD;
    }
    (b << 16) | a
}

/// Diverging blue–gray–red colormap on t ∈ [−1, 1] (trough → crest, gray at
/// undisturbed water), interpolated in linear-light RGB between fixed
/// anchors (cool and warm poles around a neutral midpoint).
pub fn diverging(t: f64) -> [u8; 3] {
    const ANCHORS: &[(f64, [u8; 3])] = &[
        (-1.0, [0x0d, 0x36, 0x6b]),
        (-0.75, [0x1c, 0x5c, 0xab]),
        (-0.5, [0x2a, 0x78, 0xd6]),
        (-0.3, [0x55, 0x98, 0xe7]),
        (-0.12, [0x9e, 0xc5, 0xf4]),
        (0.0, [0xf0, 0xef, 0xec]),
        (0.12, [0xf5, 0xb8, 0xab]),
        (0.3, [0xee, 0x8f, 0x77]),
        (0.5, [0xe3, 0x49, 0x48]),
        (0.75, [0xb0, 0x2a, 0x2a]),
        (1.0, [0x6b, 0x14, 0x14]),
    ];
    let t = t.clamp(-1.0, 1.0);
    let t = if t.is_nan() { 0.0 } else { t };
    let mut i = 0;
    while i + 2 < ANCHORS.len() && t > ANCHORS[i + 1].0 {
        i += 1;
    }
    let (t0, c0) = ANCHORS[i];
    let (t1, c1) = ANCHORS[i + 1];
    let s = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
    let mut out = [0u8; 3];
    for k in 0..3 {
        let a = srgb_to_linear(c0[k]);
        let b = srgb_to_linear(c1[k]);
        out[k] = linear_to_srgb(a + (b - a) * s);
    }
    out
}

/// Blend `c` a fraction `f` toward `toward` in linear-light RGB.
pub fn fade(c: [u8; 3], toward: [u8; 3], f: f64) -> [u8; 3] {
    let f = f.clamp(0.0, 1.0);
    let mut out = [0u8; 3];
    for k in 0..3 {
        let a = srgb_to_linear(c[k]);
        let b = srgb_to_linear(toward[k]);
        out[k] = linear_to_srgb(a + (b - a) * f);
    }
    out
}

fn srgb_to_linear(v: u8) -> f64 {
    let x = v as f64 / 255.0;
    if x <= 0.04045 {
        x / 12.92
    } else {
        ((x + 0.055) / 1.055).powf(2.4)
    }
}

fn linear_to_srgb(x: f64) -> u8 {
    let x = x.clamp(0.0, 1.0);
    let v = if x <= 0.003_130_8 {
        12.92 * x
    } else {
        1.055 * x.powf(1.0 / 2.4) - 0.055
    };
    (v * 255.0).round() as u8
}

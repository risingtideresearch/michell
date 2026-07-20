//! Zero-dependency 3D software renderer for the `render` subcommand: a
//! z-buffered Gouraud rasterizer over triangle meshes built from the wave
//! heightfield and the hull B-spline surfaces.
//!
//! Coordinates are the fleet frame with z **up** (the hull surface's
//! z-downward convention is flipped when meshing): x longitudinal (+x is
//! the direction of advance), y transverse, waterline at z = 0.

pub type V3 = [f64; 3];

pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

pub fn normalize(a: V3) -> V3 {
    let n = norm(a);
    if n > 0.0 {
        scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Blend `a` toward `b` by `t` (componentwise; colors are linear-light RGB).
pub fn lerp3(a: V3, b: V3, t: f64) -> V3 {
    add(scale(a, 1.0 - t), scale(b, t))
}

/// One mesh vertex: world position, base color (linear-light RGB in
/// [0, 1]), and Phong material (specular strength and exponent).
#[derive(Clone, Copy)]
pub struct Vertex {
    pub pos: V3,
    pub base: V3,
    pub ks: f64,
    pub spec_p: f64,
}

/// Triangle soup with shared vertices.
#[derive(Default)]
pub struct Scene {
    pub verts: Vec<Vertex>,
    pub tris: Vec<[usize; 3]>,
}

impl Scene {
    pub fn push_vert(&mut self, v: Vertex) -> usize {
        self.verts.push(v);
        self.verts.len() - 1
    }

    pub fn push_quad(&mut self, a: usize, b: usize, c: usize, d: usize) {
        self.tris.push([a, b, c]);
        self.tris.push([a, c, d]);
    }

    /// Area-weighted vertex normals (orientation-free: the shader flips
    /// each normal toward the eye, so winding need not be consistent).
    fn vertex_normals(&self) -> Vec<V3> {
        let mut acc = vec![[0.0; 3]; self.verts.len()];
        for t in &self.tris {
            let [a, b, c] = *t;
            let n = cross(
                sub(self.verts[b].pos, self.verts[a].pos),
                sub(self.verts[c].pos, self.verts[a].pos),
            );
            for &i in t {
                // Sign-align before accumulating so mirrored patches with
                // opposite winding do not cancel at shared seam vertices.
                let s = if dot(acc[i], n) < 0.0 { -1.0 } else { 1.0 };
                acc[i] = add(acc[i], scale(n, s));
            }
        }
        acc.into_iter().map(normalize).collect()
    }
}

pub struct Camera {
    pub eye: V3,
    pub target: V3,
    /// Vertical field of view [degrees].
    pub fov_deg: f64,
}

pub struct Light {
    /// Direction *toward* the sun (unit).
    pub dir: V3,
    pub ambient: f64,
    pub diffuse: f64,
}

/// Sky gradient colors, linear-light RGB (zenith at image top).
const SKY_TOP: V3 = [0.42, 0.60, 0.82];
const SKY_HORIZON: V3 = [0.87, 0.91, 0.95];

/// Render the scene to an RGB8 buffer (row 0 at the top), 2×2 supersampled.
pub fn render(scene: &Scene, cam: &Camera, light: &Light, width: usize, height: usize) -> Vec<u8> {
    const SS: usize = 2;
    let w = width * SS;
    let h = height * SS;

    // Camera basis: forward, right, up.
    let fwd = normalize(sub(cam.target, cam.eye));
    let right = normalize(cross(fwd, [0.0, 0.0, 1.0]));
    let up = cross(right, fwd);
    let fl = (h as f64 / 2.0) / (cam.fov_deg.to_radians() / 2.0).tan();
    let znear = 1e-3 * norm(sub(cam.target, cam.eye)).max(1.0);

    // Shade every vertex once (Gouraud), keeping view-space position.
    let normals = scene.vertex_normals();
    struct SV {
        view: V3,
        color: V3,
    }
    let shaded: Vec<SV> = scene
        .verts
        .iter()
        .zip(&normals)
        .map(|(v, &n0)| {
            let to_eye = normalize(sub(cam.eye, v.pos));
            // Two-sided: light the face the viewer sees.
            let n = if dot(n0, to_eye) < 0.0 { scale(n0, -1.0) } else { n0 };
            let lambert = dot(n, light.dir).max(0.0);
            let mut c = scale(v.base, light.ambient + light.diffuse * lambert);
            if v.ks > 0.0 {
                let refl = sub(scale(n, 2.0 * dot(n, light.dir)), light.dir);
                let spec = dot(refl, to_eye).max(0.0).powf(v.spec_p);
                c = add(c, scale([1.0, 1.0, 1.0], v.ks * spec));
            }
            let rel = sub(v.pos, cam.eye);
            SV {
                view: [dot(rel, right), dot(rel, up), dot(rel, fwd)],
                color: c,
            }
        })
        .collect();

    let mut zbuf = vec![f64::INFINITY; w * h];
    let mut color = vec![[0.0f64; 3]; w * h];
    let mut covered = vec![false; w * h];

    // Clip each triangle to z >= znear in view space, then rasterize.
    let mut poly: Vec<(V3, V3)> = Vec::with_capacity(4);
    for t in &scene.tris {
        poly.clear();
        for k in 0..3 {
            let a = &shaded[t[k]];
            let b = &shaded[t[(k + 1) % 3]];
            if a.view[2] >= znear {
                poly.push((a.view, a.color));
            }
            if (a.view[2] >= znear) != (b.view[2] >= znear) {
                let s = (znear - a.view[2]) / (b.view[2] - a.view[2]);
                poly.push((lerp3(a.view, b.view, s), lerp3(a.color, b.color, s)));
            }
        }
        if poly.len() < 3 {
            continue;
        }
        // Project the (convex) clipped polygon and fan-triangulate.
        let proj: Vec<(f64, f64, f64, V3)> = poly
            .iter()
            .map(|&(v, c)| {
                (
                    w as f64 / 2.0 + fl * v[0] / v[2],
                    h as f64 / 2.0 - fl * v[1] / v[2],
                    v[2],
                    c,
                )
            })
            .collect();
        for k in 1..proj.len() - 1 {
            raster_tri(
                &[proj[0], proj[k], proj[k + 1]],
                w,
                h,
                &mut zbuf,
                &mut color,
                &mut covered,
            );
        }
    }

    // Downsample, filling uncovered samples with the sky gradient.
    let mut rgb = vec![0u8; 3 * width * height];
    for py in 0..height {
        for px in 0..width {
            let mut acc = [0.0; 3];
            for sy in 0..SS {
                for sx in 0..SS {
                    let i = (py * SS + sy) * w + px * SS + sx;
                    let c = if covered[i] {
                        color[i]
                    } else {
                        let f = (py * SS + sy) as f64 / h as f64;
                        lerp3(SKY_TOP, SKY_HORIZON, (f * 1.6).min(1.0))
                    };
                    acc = add(acc, c);
                }
            }
            let c = scale(acc, 1.0 / (SS * SS) as f64);
            let o = 3 * (py * width + px);
            for k in 0..3 {
                rgb[o + k] = linear_to_srgb(c[k]);
            }
        }
    }
    rgb
}

/// Rasterize one screen-space triangle with affine attribute interpolation.
#[allow(clippy::too_many_arguments)]
fn raster_tri(
    v: &[(f64, f64, f64, V3); 3],
    w: usize,
    h: usize,
    zbuf: &mut [f64],
    color: &mut [V3],
    covered: &mut [bool],
) {
    let area = (v[1].0 - v[0].0) * (v[2].1 - v[0].1) - (v[1].1 - v[0].1) * (v[2].0 - v[0].0);
    if area.abs() < 1e-12 {
        return;
    }
    let x0 = v.iter().map(|p| p.0).fold(f64::INFINITY, f64::min).floor().max(0.0) as usize;
    let x1 = (v.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max).ceil() as isize)
        .clamp(0, w as isize - 1) as usize;
    let y0 = v.iter().map(|p| p.1).fold(f64::INFINITY, f64::min).floor().max(0.0) as usize;
    let y1 = (v.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max).ceil() as isize)
        .clamp(0, h as isize - 1) as usize;
    if x0 > x1 || y0 > y1 {
        return;
    }
    let inv = 1.0 / area;
    for py in y0..=y1 {
        let y = py as f64 + 0.5;
        for px in x0..=x1 {
            let x = px as f64 + 0.5;
            let w0 = ((v[1].0 - x) * (v[2].1 - y) - (v[1].1 - y) * (v[2].0 - x)) * inv;
            let w1 = ((v[2].0 - x) * (v[0].1 - y) - (v[2].1 - y) * (v[0].0 - x)) * inv;
            let w2 = 1.0 - w0 - w1;
            if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 {
                continue;
            }
            let z = w0 * v[0].2 + w1 * v[1].2 + w2 * v[2].2;
            let i = py * w + px;
            if z < zbuf[i] {
                zbuf[i] = z;
                covered[i] = true;
                color[i] =
                    std::array::from_fn(|k| w0 * v[0].3[k] + w1 * v[1].3[k] + w2 * v[2].3[k]);
            }
        }
    }
}

pub fn srgb_to_linear_f(v: u8) -> f64 {
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

/// Linear-light color from a `0xRRGGBB` literal.
pub fn hex(c: u32) -> V3 {
    [
        srgb_to_linear_f((c >> 16) as u8),
        srgb_to_linear_f((c >> 8) as u8),
        srgb_to_linear_f(c as u8),
    ]
}

// ---------------------------------------------------------------------------
// Scene building
// ---------------------------------------------------------------------------

const WATER_DEEP: u32 = 0x1c4f86;
const WATER_LIGHT: u32 = 0xa9cfee;
/// Pale tint blended over the water where the free-wave field is not
/// physical (forward of the aft-most stern), matching the 2D heatmap.
const WATER_PALE: u32 = 0xdfe6ec;
const FADE_3D: f64 = 0.3;
const WATER_KS: f64 = 0.09;
const WATER_SPEC_P: f64 = 36.0;

/// Fade parameters for the non-physical part of the free-wave field: pale
/// from `x_phys` (the aft-most stern) forward, feathered over `feather`
/// metres so the boundary reads as a gradient, not a seam.
#[derive(Clone, Copy)]
pub struct PhysFade {
    pub x_phys: f64,
    pub feather: f64,
}

fn water_vertex(x: f64, y: f64, zeta: f64, z_scale: f64, vmax: f64, fade: PhysFade) -> Vertex {
    let t = (zeta / vmax).clamp(-1.0, 1.0);
    // Mild gamma keeps small-amplitude structure visible.
    let t = t.signum() * t.abs().powf(0.65);
    let mut base = lerp3(hex(WATER_DEEP), hex(WATER_LIGHT), 0.5 * (t + 1.0));
    let s = ((x - fade.x_phys) / fade.feather.max(1e-9)).clamp(0.0, 1.0);
    if s > 0.0 {
        base = lerp3(base, hex(WATER_PALE), FADE_3D * s * s * (3.0 - 2.0 * s));
    }
    Vertex {
        pos: [x, y, zeta * z_scale],
        base,
        ks: WATER_KS,
        spec_p: WATER_SPEC_P,
    }
}

/// Wave heightfield over the grid, tinted by elevation (saturating at
/// `vmax`) and faded ahead of `x_phys` (the aft-most stern).
pub fn add_water(
    scene: &mut Scene,
    grid: &michell::WaveGrid,
    z_scale: f64,
    vmax: f64,
    fade: PhysFade,
) {
    let base_idx = scene.verts.len();
    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            // Taper the elevation to zero over the outermost few percent of
            // the grid so the chopped wake dissolves into the flat apron
            // instead of ending in a torn edge.
            let ex = edge_taper(ix, grid.nx);
            let ey = edge_taper(iy, grid.ny);
            let v = water_vertex(
                grid.x(ix),
                grid.y(iy),
                grid.get(ix, iy) * ex * ey,
                z_scale,
                vmax,
                fade,
            );
            scene.push_vert(v);
        }
    }
    for iy in 0..grid.ny - 1 {
        for ix in 0..grid.nx - 1 {
            let a = base_idx + iy * grid.nx + ix;
            let b = a + 1;
            let c = a + grid.nx + 1;
            let d = a + grid.nx;
            scene.push_quad(a, b, c, d);
        }
    }
}

/// Smoothstep from 0 at the grid boundary to 1 past a 4% margin.
fn edge_taper(i: usize, n: usize) -> f64 {
    let m = (0.04 * n as f64).max(2.0);
    let d = (i.min(n - 1 - i) as f64 / m).min(1.0);
    d * d * (3.0 - 2.0 * d)
}

/// Flat still-water apron from the grid edge out to the horizon so the
/// heightfield does not end in visible sky. Inner ring matches the grid
/// boundary exactly (no cracks); outer rings sit at z = 0.
pub fn add_skirt(
    scene: &mut Scene,
    grid: &michell::WaveGrid,
    z_scale: f64,
    vmax: f64,
    fade: PhysFade,
) {
    // Perimeter of the grid, counter-clockwise, sampled coarsely.
    let step = (grid.nx.max(grid.ny) / 192).max(1);
    let mut loop_ij: Vec<(usize, usize)> = Vec::new();
    let mut ix = 0;
    while ix < grid.nx {
        loop_ij.push((ix, 0));
        ix += step;
    }
    let mut iy = 0;
    while iy < grid.ny {
        loop_ij.push((grid.nx - 1, iy));
        iy += step;
    }
    let mut ix = grid.nx as isize - 1;
    while ix >= 0 {
        loop_ij.push((ix as usize, grid.ny - 1));
        ix -= step as isize;
    }
    let mut iy = grid.ny as isize - 1;
    while iy >= 0 {
        loop_ij.push((0, iy as usize));
        iy -= step as isize;
    }
    loop_ij.dedup();
    if loop_ij.last() == loop_ij.first() {
        loop_ij.pop();
    }
    let n = loop_ij.len();

    let cx = 0.5 * (grid.x0 + grid.x1);
    let cy = 0.5 * (grid.y0 + grid.y1);
    const RINGS: usize = 16;
    const REACH: f64 = 25.0;
    let mut ring_start = Vec::with_capacity(RINGS + 1);
    for r in 0..=RINGS {
        ring_start.push(scene.verts.len());
        // Geometric expansion factor: 1 at the grid edge, REACH at the rim.
        let f = REACH.powf(r as f64 / RINGS as f64);
        for &(i, j) in &loop_ij {
            let (x, y) = (grid.x(i), grid.y(j));
            let (xs, ys) = (cx + (x - cx) * f, cy + (y - cy) * f);
            // The water mesh tapers to ζ = 0 at its boundary, so the whole
            // apron is flat and ring 0 meets it without a crack.
            let mut v = water_vertex(xs, ys, 0.0, z_scale, vmax, fade);
            // Matte out distant water: specular interpolated across the
            // ever-larger outer quads reads as false glint edges.
            v.ks *= (1.0 - r as f64 / RINGS as f64).powi(2);
            scene.push_vert(v);
        }
    }
    for r in 0..RINGS {
        for k in 0..n {
            let a = ring_start[r] + k;
            let b = ring_start[r] + (k + 1) % n;
            let c = ring_start[r + 1] + (k + 1) % n;
            let d = ring_start[r + 1] + k;
            scene.push_quad(a, b, c, d);
        }
    }
}

const HULL_WETTED: u32 = 0x7a4034;
const HULL_TOPSIDE: u32 = 0xe9e5db;
const HULL_DECK: u32 = 0xd8d2c6;
const HULL_KS: f64 = 0.06;
const HULL_SPEC_P: f64 = 16.0;

/// One hull: wetted surface (both sides, keel to waterline), wall-sided
/// topsides extruded to `freeboard` above the waterline, and a deck cap.
/// `px`/`py` place the hull in the fleet frame; z is flipped to point up.
pub fn add_hull(
    scene: &mut Scene,
    surf: &michell::BSplineSurface,
    px: f64,
    py: f64,
    freeboard: f64,
) {
    let (hx0, hx1) = surf.x_domain();
    let (_, draft) = surf.z_domain();
    const NX: usize = 96;
    const NZ: usize = 16;
    const NFB: usize = 3;
    const NDK: usize = 6;
    let xs: Vec<f64> = (0..=NX)
        .map(|i| hx0 + (hx1 - hx0) * i as f64 / NX as f64)
        .collect();
    let half: Vec<f64> = xs.iter().map(|&x| surf.eval(x, 0.0).max(0.0)).collect();

    let mat = |base: u32| Vertex {
        pos: [0.0; 3],
        base: hex(base),
        ks: HULL_KS,
        spec_p: HULL_SPEC_P,
    };

    for side in [1.0f64, -1.0] {
        // Wetted: z downward 0..draft, world z = -zd.
        let start = scene.verts.len();
        for j in 0..=NZ {
            let zd = draft * j as f64 / NZ as f64;
            for (i, &x) in xs.iter().enumerate() {
                let hb = if j == 0 {
                    half[i]
                } else {
                    surf.eval(x, zd).max(0.0)
                };
                let mut v = mat(HULL_WETTED);
                v.pos = [x + px, py + side * hb, -zd];
                scene.push_vert(v);
            }
        }
        grid_quads(scene, start, NX + 1, NZ + 1);

        // Topsides: waterline section extruded up to the freeboard.
        let start = scene.verts.len();
        for j in 0..=NFB {
            let z = freeboard * j as f64 / NFB as f64;
            for (i, _) in xs.iter().enumerate() {
                let mut v = mat(HULL_TOPSIDE);
                v.pos = [xs[i] + px, py + side * half[i], z];
                scene.push_vert(v);
            }
        }
        grid_quads(scene, start, NX + 1, NFB + 1);
    }

    // Deck: strip across the sheer line at z = freeboard.
    let start = scene.verts.len();
    for j in 0..=NDK {
        let f = 2.0 * j as f64 / NDK as f64 - 1.0;
        for (i, _) in xs.iter().enumerate() {
            let mut v = mat(HULL_DECK);
            v.pos = [xs[i] + px, py + f * half[i], freeboard];
            scene.push_vert(v);
        }
    }
    grid_quads(scene, start, NX + 1, NDK + 1);
}

/// Quads over a `cols` × `rows` vertex lattice appended at `start`
/// (row-major, column index fastest).
fn grid_quads(scene: &mut Scene, start: usize, cols: usize, rows: usize) {
    for j in 0..rows - 1 {
        for i in 0..cols - 1 {
            let a = start + j * cols + i;
            let b = a + 1;
            let c = a + cols + 1;
            let d = a + cols;
            scene.push_quad(a, b, c, d);
        }
    }
}

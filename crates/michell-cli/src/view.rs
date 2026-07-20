//! `michell view` — an interactive local viewer for fleet wave fields.
//!
//! A zero-dependency HTTP server (in the hand-rolled spirit of the rest of the
//! crate — see [`crate::png`]) serves a single-page browser front end. The
//! physics runs native: each hull's *standalone* free-wave field is computed
//! once per (speed, displacement) on a local grid, and the browser composites
//! the fleet field by translating and summing those per-hull grids as hulls
//! are dragged around. That superposition is exact in thin-ship theory (the
//! fleet amplitude is a sum of per-hull amplitudes, each carrying only a
//! placement phase; see [`michell::FreeWaveSpectrum`]), so dragging is a pure
//! client-side recomposite at interactive rates with no physics re-run.
//!
//! Only speed and displacement changes force a native recompute (ν changes, or
//! the wetted hull is re-floated): both are parallelised across hulls. These
//! are explicit, one-shot actions in the UI (a speed dropdown; a displacement
//! slider that re-floats on release) with a spinner while they run, since each
//! blocks the single-threaded server for a fraction of a second up to ~2 s.
//!
//! Endpoints (all GET; this is a single-user local tool):
//!   /                         the page
//!   /api/state                fleet + view + colour scale as JSON
//!   /api/field?hull=I         hull I's local ζ grid as little-endian f32
//!   /api/speed?u=U            recompute all fields at speed U, return state
//!   /api/displacement?mass=M  re-float the assembly at total mass M
//!   /api/resistance?p=x,y;... resistance for the given world placements

use crate::formats::{body_options, load_body};
use crate::{load_fleet, parse_args, Member};
use michell::body::{Body, BodyOptions};
use michell::float::{solve_equilibrium_bodies, LoadCase};
use michell::iges::{source_fleet, HullPose, ImportOptions, Platform};
use michell::{Conditions, FreeWaveSpectrum, Hull, Placement};
use std::io::{BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

const PAGE: &str = include_str!("view.html");

/// A hull's standalone wave field on a local grid (hull-fixed coordinates,
/// row-major over y then x). The client shifts it by the hull's world
/// placement and sums across the fleet.
struct Field {
    lx0: f64,
    lx1: f64,
    ly0: f64,
    ly1: f64,
    nx: usize,
    ny: usize,
    zeta: Vec<f32>,
}

/// One fleet member the viewer can manipulate.
struct ViewHull {
    name: String,
    hull: Hull,
    /// Full-band body, retained so displacement can re-float; `None` for
    /// sources that cannot be re-floated (wetted-only nets, IGES/STL fleets).
    body: Option<Body>,
    /// Design placement in world coordinates; drag offsets are added by the
    /// client on top of this and never mutate it.
    home: Placement,
    /// Beam profile at the waterline: (local x, half-beam) for the glyph.
    beam: Vec<(f64, f64)>,
    x0: f64,
    x1: f64,
    field: Field,
    /// True when the last re-float lifted this hull entirely out of the water
    /// (a light load can raise slender outriggers clear): it then makes no
    /// wake and contributes no resistance.
    dry: bool,
}

struct ViewState {
    hulls: Vec<ViewHull>,
    cond: Conditions,
    l_ref: f64,
    view: [f64; 4],
    margin: f64,
    base_px: usize,
    design_mass: f64,
    mass: f64,
    /// Colour saturation elevation [m]: fixed while dragging, refreshed on
    /// speed/displacement change.
    vmax: f64,
    body_opts: BodyOptions,
    /// Placement-independent resistance pieces, refreshed with the fields on
    /// every speed/displacement change so the live readout only has to redo
    /// the (placement-dependent) combined wave integral.
    viscous_total: f64,
    wetted_surface: f64,
    /// Σ of each hull's standalone wave resistance, by the SAME fixed-grid
    /// estimate the combined figure uses, so their ratio (the interference
    /// factor) is consistent and → 1 as hulls separate.
    solo_wave_total: f64,
}

pub fn cmd_view(args: &[String]) -> Result<(), String> {
    let p = parse_args(args)?;
    if p.positional.is_empty() {
        return Err("usage: michell view <hull>... --speed U [--port N] [options]".into());
    }
    let settings = p.load_settings()?;
    let members = load_fleet(&p.positional, &settings)?;

    let l_ref = members
        .iter()
        .map(|m| m.hull.length())
        .fold(0.0f64, f64::max);
    let u = crate::single_speed(&p, l_ref)?;
    let cond = p.conditions(u)?;
    let port: u16 = match p.flag("port") {
        Some(s) => s
            .parse()
            .map_err(|_| format!("--port: cannot parse {s:?}"))?,
        None => 8737,
    };
    let base_px = match p.flag("size") {
        Some(s) => s
            .parse::<usize>()
            .map_err(|_| format!("--size: cannot parse {s:?}"))?
            .clamp(120, 1600),
        None => 420,
    };

    let body_opts = body_options(&settings);

    // Retain a Body per member where the source is a single-body `.hull`
    // file, so displacement can re-float it. A file contributing several
    // hulls (an IGES fleet) is not re-floatable member-by-member here.
    let mut path_counts = std::collections::HashMap::new();
    for m in &members {
        *path_counts.entry(m.path.clone()).or_insert(0usize) += 1;
    }

    // A single IGES file contributes several hulls that float as one platform.
    // Decompose it once into full-band bodies so displacement can re-float via
    // the fast body path (re-floating the raw IGES source directly would
    // re-loft the NURBS patches every Newton step — tens of seconds).
    let iges_bodies = bodies_from_iges(&p.positional, &settings);

    let mut hulls: Vec<ViewHull> = Vec::new();
    let mut design_mass = 0.0;
    for (idx, m) in members.iter().enumerate() {
        let body = match &iges_bodies {
            Some(bodies) if bodies.len() == members.len() => Some(bodies[idx].clone()),
            _ if path_counts[&m.path] == 1 => load_body(&m.path).ok(),
            _ => None,
        };
        design_mass += m.hull.displaced_volume() * cond.fluid.density;
        hulls.push(build_view_hull(m, body));
    }

    // Disambiguate duplicate names (a multi-hull IGES gives every hull the
    // file stem) so the fleet reads "full #1 · full #2 · …".
    let mut name_counts = std::collections::HashMap::new();
    for vh in &hulls {
        *name_counts.entry(vh.name.clone()).or_insert(0) += 1;
    }
    if name_counts.values().any(|&c| c > 1) {
        for (i, vh) in hulls.iter_mut().enumerate() {
            vh.name = format!("{} #{}", vh.name, i + 1);
        }
    }

    // World view: as the `wake` default region, over the whole fleet.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_abs = 0.0f64;
    for h in &hulls {
        x_lo = x_lo.min(h.x0 + h.home.x);
        x_hi = x_hi.max(h.x1 + h.home.x);
        y_abs = y_abs.max(h.home.y.abs());
    }
    let x1 = x_hi + 0.35 * l_ref;
    let x0 = x_lo - 3.0 * l_ref;
    let yh = (0.42 * (x1 - x0)).max(y_abs + 0.8 * l_ref);
    let view = [x0, x1, -yh, yh];
    // Drag margin: how far a hull may roam from home in each direction. The
    // per-hull local grids are sized to cover the whole view for any
    // placement within this margin; keep it modest so recompute stays snappy.
    let margin = (0.22 * (x1 - x0)).max(1.3 * l_ref);

    let mut state = ViewState {
        hulls,
        cond,
        l_ref,
        view,
        margin,
        base_px,
        design_mass,
        mass: design_mass,
        vmax: 0.0,
        body_opts,
        viscous_total: 0.0,
        wetted_surface: 0.0,
        solo_wave_total: 0.0,
    };
    recompute_fields(&mut state)?;

    let addr = format!("127.0.0.1:{port}");
    let listener = TcpListener::bind(&addr)
        .map_err(|e| format!("cannot bind {addr}: {e} (try a different --port)"))?;
    let url = format!("http://{addr}/");
    println!(
        "michell view: {} hull(s) at U = {:.3} m/s (Fn {:.3})\n  open {url}\n  \
         drag hulls to move them; sliders change speed and displacement; Ctrl-C to stop",
        state.hulls.len(),
        state.cond.speed,
        state.cond.froude_number(state.l_ref),
    );
    let _ = std::process::Command::new("open").arg(&url).spawn();

    // One thread per connection. Browsers open speculative "preconnect"
    // sockets that may send no request; a single-threaded accept loop would
    // block on reading one of those and stall every real request behind it.
    // Shared state is behind the mutex, so recomputes still serialise safely.
    let state = Arc::new(Mutex::new(state));
    for stream in listener.incoming() {
        match stream {
            Ok(s) => {
                let st = Arc::clone(&state);
                std::thread::spawn(move || {
                    if let Err(e) = handle(s, &st) {
                        eprintln!("view: connection error: {e}");
                    }
                });
            }
            Err(e) => eprintln!("view: accept error: {e}"),
        }
    }
    Ok(())
}

fn build_view_hull(m: &Member, body: Option<Body>) -> ViewHull {
    let (x0, x1) = m.hull.surface().x_domain();
    let n = 96usize;
    let beam = (0..n)
        .map(|i| {
            let x = x0 + (x1 - x0) * i as f64 / (n - 1) as f64;
            (x, m.hull.surface().eval(x, 0.0).max(0.0))
        })
        .collect();
    ViewHull {
        name: hull_name(&m.path),
        hull: m.hull.clone(),
        body,
        home: m.placement,
        beam,
        x0,
        x1,
        field: Field {
            lx0: 0.0,
            lx1: 0.0,
            ly0: 0.0,
            ly1: 0.0,
            nx: 0,
            ny: 0,
            zeta: Vec::new(),
        },
        dry: false,
    }
}

/// Decompose a single IGES file into per-hull full-band bodies (as `michell
/// loft` does), so displacement can re-float them via the fast body path.
/// `None` unless the input is exactly one IGES file with no `@` suffix.
fn bodies_from_iges(
    positional: &[String],
    settings: &crate::formats::LoadSettings,
) -> Option<Vec<Body>> {
    if positional.len() != 1 {
        return None;
    }
    let path = &positional[0];
    if path.contains('@') {
        return None;
    }
    let lower = path.to_ascii_lowercase();
    if !(lower.ends_with(".igs") || lower.ends_with(".iges")) {
        return None;
    }
    let text = std::fs::read_to_string(path).ok()?;
    let design_wl = settings.waterline_z;
    let opts = ImportOptions {
        waterline_z: design_wl,
        stations: settings.samples.0,
        waterlines: settings.samples.1,
        fit: if settings.fit_explicit {
            settings.fit
        } else {
            ImportOptions::default().fit
        },
        centerplane: settings.centerplane,
    };
    let src = source_fleet(&text, design_wl).ok()?;
    let mut bodies = Vec::with_capacity(src.len());
    for idx in 0..src.len() {
        let bottom = src.hull_z_bottom(idx);
        let top = src.hull_z_top(idx);
        let draft_est = design_wl - bottom;
        if draft_est <= 0.0 {
            return None;
        }
        // Band from keel to half a draft above the design waterline (as loft).
        let band_top = (design_wl + 0.5 * draft_est).min(top);
        let mut hull_opts = opts;
        if hull_opts.centerplane.is_none() {
            // Detect the centerplane at the design waterline (coarse).
            let mut d = opts;
            d.stations = 61;
            d.waterlines = 17;
            d.fit.n_ctrl_x = d.fit.n_ctrl_x.min(10);
            d.fit.n_ctrl_z = d.fit.n_ctrl_z.min(7);
            hull_opts.centerplane = src
                .situate_one(idx, design_wl, &HullPose::default(), &Platform::default(), &d)
                .ok()?
                .map(|m| m.report.centerplane);
        }
        let m = src
            .situate_one(
                idx,
                band_top,
                &HullPose::default(),
                &Platform::default(),
                &hull_opts,
            )
            .ok()??;
        let wl_depth = band_top - design_wl;
        bodies.push(Body::new(m.hull.surface().clone(), wl_depth, m.report.centerplane).ok()?);
    }
    Some(bodies)
}

fn hull_name(path: &str) -> String {
    std::path::Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(path)
        .to_string()
}

/// Recompute every hull's standalone field at the current speed, in parallel
/// across hulls, then refresh the shared colour scale from the composited
/// fleet field at the home placements.
fn recompute_fields(state: &mut ViewState) -> Result<(), String> {
    let cond = state.cond;
    let l_ref = state.l_ref;
    let [vx0, vx1, vy0, vy1] = state.view;
    let margin = state.margin;
    let base_px = state.base_px;

    // Common pixel spacing: the view resolution, but never coarser than ~8
    // points per transverse wavelength (which shrinks with speed) so the grid
    // always resolves the pattern.
    let nu = cond.gravity / (cond.speed * cond.speed);
    let lambda = 2.0 * std::f64::consts::PI / nu;
    let dx = ((vx1 - vx0) / (base_px as f64 - 1.0)).min(lambda / 8.0);
    let dy = dx;

    // Geometry of each hull's local grid (world = local + placement, the
    // placement roaming within home ± margin, so local spans the view minus
    // that range).
    let geoms: Vec<Field> = state
        .hulls
        .iter()
        .map(|vh| {
            if vh.dry {
                // Out of the water: no field (sampling an empty grid is zero).
                return Field { lx0: 0.0, lx1: 0.0, ly0: 0.0, ly1: 0.0, nx: 0, ny: 0, zeta: Vec::new() };
            }
            let lx0 = vx0 - vh.home.x - margin;
            let lx1 = vx1 - vh.home.x + margin;
            let ly0 = vy0 - vh.home.y - margin;
            let ly1 = vy1 - vh.home.y + margin;
            let nx = (((lx1 - lx0) / dx).round() as usize + 1).clamp(2, 1400);
            let ny = (((ly1 - ly0) / dy).round() as usize + 1).clamp(2, 1400);
            Field { lx0, lx1, ly0, ly1, nx, ny, zeta: vec![0.0; nx * ny] }
        })
        .collect();

    // Fan out over (hull × row-band) so every core is busy even for a
    // single-hull fleet. Bands share the row spacing exactly, so stitching is
    // seamless. Each band returns its rows; we copy them into place after.
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    let nh = geoms.len();
    let bands_per_hull = ((2 * cores) / nh.max(1)).max(1);
    let hulls = &state.hulls;

    type BandTask = (usize, usize, usize); // (hull, row0, row1)
    let mut tasks: Vec<BandTask> = Vec::new();
    for (i, g) in geoms.iter().enumerate() {
        let k = bands_per_hull.min(g.ny.div_ceil(24)).max(1);
        for b in 0..k {
            let r0 = b * g.ny / k;
            let r1 = (b + 1) * g.ny / k;
            if r1 > r0 {
                tasks.push((i, r0, r1));
            }
        }
    }

    let cond_ref = &cond;
    // (hull index, first row, that band's rows).
    type BandResult = Result<(usize, usize, Vec<f32>), String>;
    let results: Vec<BandResult> = std::thread::scope(|scope| {
        let handles: Vec<_> = tasks
            .iter()
            .map(|&(i, r0, r1)| {
                let g = &geoms[i];
                let hull = &hulls[i].hull;
                scope.spawn(move || {
                    let by0 = grid_row_coord(g.ly0, g.ly1, g.ny, r0);
                    let by1 = grid_row_coord(g.ly0, g.ly1, g.ny, r1 - 1);
                    let members = [(hull, Placement { x: 0.0, y: 0.0 })];
                    let mut spec =
                        FreeWaveSpectrum::new(&members, cond_ref).map_err(|e| format!("{e}"))?;
                    let grid = spec
                        .elevation_grid(g.lx0, g.lx1, by0, by1, g.nx, r1 - r0)
                        .map_err(|e| format!("{e}"))?;
                    Ok((i, r0, grid.zeta.iter().map(|&v| v as f32).collect()))
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    let mut fields = geoms;
    for r in results {
        let (i, r0, rows) = r?;
        let nx = fields[i].nx;
        fields[i].zeta[r0 * nx..r0 * nx + rows.len()].copy_from_slice(&rows);
    }
    for (vh, f) in state.hulls.iter_mut().zip(fields) {
        vh.field = f;
    }

    // Placement-independent resistance pieces (viscous, wetted surface, and
    // each hull's standalone wave resistance) — recomputed here so the live
    // readout only redoes the combined wave integral.
    let (sec_cap, n_theta) = resistance_grid(&cond, l_ref);
    let mut viscous_total = 0.0;
    let mut wetted_surface = 0.0;
    let mut solo_wave_total = 0.0;
    for vh in state.hulls.iter().filter(|h| !h.dry) {
        viscous_total += michell::viscous_resistance(&vh.hull, &cond)
            .map(|v| v.resistance)
            .unwrap_or(0.0);
        wetted_surface += vh.hull.wetted_surface();
        let mut solo = FreeWaveSpectrum::new(&[(&vh.hull, Placement { x: 0.0, y: 0.0 })], &cond)
            .map_err(|e| format!("{e}"))?;
        solo_wave_total += fast_wave_resistance(&mut solo, sec_cap, n_theta);
    }
    state.viscous_total = viscous_total;
    state.wetted_surface = wetted_surface;
    state.solo_wave_total = solo_wave_total;

    state.vmax = colour_scale(state);
    Ok(())
}

/// Fixed-grid estimate of wave resistance R_w = ∫ (dR_w/dθ) dθ over
/// (−θ_max, θ_max). Simpson on a uniform θ grid: unlike the library's adaptive
/// integrator its cost is independent of hull separation, so the live readout
/// stays responsive when hulls are dragged far apart (where the interference
/// integrand oscillates fast and the exact integrator slows to ~a second).
///
/// The diverging-wave tail (θ → ±π/2, λ = sec θ → ∞) decays only slowly, and
/// carries more of the resistance at high speed (small ν), so `sec_cap` — how
/// far toward π/2 to integrate — is raised with speed rather than fixed at the
/// renderer's cap of 15. This matches the CLI `resistance` figure to display
/// precision across the speed range; that command remains the reference.
fn fast_wave_resistance(spec: &mut FreeWaveSpectrum, sec_cap: f64, n: usize) -> f64 {
    let theta_max = (1.0 / sec_cap).acos();
    let n = n & !1; // force even for Simpson
    let h = 2.0 * theta_max / n as f64;
    let mut sum = spec.resistance_density(-theta_max) + spec.resistance_density(theta_max);
    for i in 1..n {
        let th = -theta_max + h * i as f64;
        let w = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += w * spec.resistance_density(th);
    }
    sum * h / 3.0
}

/// Integration cap (largest sec θ = λ) and node count for the wave-resistance
/// estimate at the current speed. Higher speed ⇒ slower spectral decay ⇒
/// integrate further toward π/2, with proportionally more nodes to keep the
/// diverging-wave sliver resolved.
fn resistance_grid(cond: &Conditions, l_ref: f64) -> (f64, usize) {
    // The spectral tail reaches to larger λ = sec θ at higher Froude number
    // (the exp(−νλ²z) decay is slower when ν = g/U² is small), so scale the
    // integration cap with Fn. Calibrated against the CLI `resistance` figure.
    let fn_ = cond.froude_number(l_ref.max(1e-6));
    let sec_cap = (15.0 + 60.0 * fn_).clamp(15.0, 90.0);
    // The dR_w/dθ envelope is smooth, so accuracy is set by the cap, not node
    // density; a modest count keeps the call fast (tens of ms).
    let n = ((sec_cap * 120.0) as usize).clamp(4000, 12000);
    (sec_cap, n)
}

/// y coordinate of row `i` of an inclusive `n`-point grid over [a, b] — the
/// same rule `WaveGrid` uses, so row-bands land on identical sample points.
fn grid_row_coord(a: f64, b: f64, n: usize, i: usize) -> f64 {
    if n <= 1 {
        a
    } else {
        a + (b - a) * i as f64 / (n - 1) as f64
    }
}

/// 99.5th-percentile |ζ| of the fleet field composited at the home
/// placements — the same saturation rule the `wake` PNG uses.
fn colour_scale(state: &ViewState) -> f64 {
    let [vx0, vx1, vy0, vy1] = state.view;
    let nx = 240usize;
    let ny = 200usize;
    let mut abs = Vec::with_capacity(nx * ny);
    for iy in 0..ny {
        let y = vy0 + (vy1 - vy0) * iy as f64 / (ny - 1) as f64;
        for ix in 0..nx {
            let x = vx0 + (vx1 - vx0) * ix as f64 / (nx - 1) as f64;
            let mut z = 0.0;
            for vh in &state.hulls {
                z += sample_field(&vh.field, x - vh.home.x, y - vh.home.y);
            }
            abs.push(z.abs());
        }
    }
    abs.sort_by(|a, b| a.total_cmp(b));
    abs[((abs.len() - 1) as f64 * 0.995) as usize].max(1e-12)
}

/// Bilinear sample of a field at local (x, y); zero outside the grid.
fn sample_field(f: &Field, x: f64, y: f64) -> f64 {
    if f.nx < 2 || f.ny < 2 || x < f.lx0 || x > f.lx1 || y < f.ly0 || y > f.ly1 {
        return 0.0;
    }
    let fx = (x - f.lx0) / (f.lx1 - f.lx0) * (f.nx - 1) as f64;
    let fy = (y - f.ly0) / (f.ly1 - f.ly0) * (f.ny - 1) as f64;
    let ix = (fx.floor() as usize).min(f.nx - 2);
    let iy = (fy.floor() as usize).min(f.ny - 2);
    let tx = fx - ix as f64;
    let ty = fy - iy as f64;
    let g = |ax: usize, ay: usize| f.zeta[ay * f.nx + ax] as f64;
    let a = g(ix, iy) * (1.0 - tx) + g(ix + 1, iy) * tx;
    let b = g(ix, iy + 1) * (1.0 - tx) + g(ix + 1, iy + 1) * tx;
    a * (1.0 - ty) + b * ty
}

/// Re-float the whole assembly (bodies only) at a new total mass and rebuild
/// the wetted hulls, then recompute fields. Members without a retained body
/// keep their design geometry.
fn set_displacement(state: &mut ViewState, mass: f64) -> Result<(), String> {
    let opts = state.body_opts;
    let density = state.cond.fluid.density;

    // Solve the assembly's common flotation for this total mass.
    let (sinkage, trim) = {
        let bodies: Vec<&Body> = state.hulls.iter().filter_map(|h| h.body.as_ref()).collect();
        if bodies.is_empty() {
            return Err("no re-floatable bodies in this fleet (displacement is fixed)".into());
        }
        let poses = vec![HullPose::default(); bodies.len()];
        let eq = solve_equilibrium_bodies(
            &bodies,
            0.0,
            &poses,
            &LoadCase { mass, lcg: None },
            density,
            &opts,
        )
        .map_err(|e| format!("{e}"))?;
        (eq.sinkage, eq.trim)
    };

    // Re-loft each body at the solved platform state. A body that comes back
    // dry (lifted clear at a light load) is flagged, not an error.
    let platform = Platform {
        sinkage,
        trim,
        pivot_x: 0.0,
    };
    for vh in state.hulls.iter_mut() {
        let situated = match &vh.body {
            Some(body) => body
                .situate(0.0, &HullPose::default(), &platform, &opts)
                .map_err(|e| format!("{e}"))?,
            None => continue, // not re-floatable: leave as loaded
        };
        match situated {
            Some(sb) => {
                set_hull_geometry(vh, sb.hull);
                vh.dry = false;
            }
            None => vh.dry = true,
        }
    }
    state.mass = mass;
    recompute_fields(state)
}

/// Replace a view hull's wetted geometry (and its waterline beam profile for
/// the glyph) after a re-float.
fn set_hull_geometry(vh: &mut ViewHull, hull: Hull) {
    let (x0, x1) = hull.surface().x_domain();
    let n = vh.beam.len().max(2);
    vh.beam = (0..n)
        .map(|i| {
            let x = x0 + (x1 - x0) * i as f64 / (n - 1) as f64;
            (x, hull.surface().eval(x, 0.0).max(0.0))
        })
        .collect();
    vh.x0 = x0;
    vh.x1 = x1;
    vh.hull = hull;
}

// ---------------------------------------------------------------------------
// HTTP
// ---------------------------------------------------------------------------

fn handle(stream: TcpStream, state: &Mutex<ViewState>) -> std::io::Result<()> {
    let _ = stream.set_nodelay(true);
    // Idle timeout: a keep-alive connection (or an idle speculative socket)
    // that sends no next request is closed after this, freeing the thread.
    let _ = stream.set_read_timeout(Some(Duration::from_secs(20)));
    let mut reader = BufReader::new(stream.try_clone()?);
    let mut stream = stream;

    // HTTP/1.1 keep-alive: serve requests on this connection until the client
    // closes it or goes idle. Forcing a new connection per request (Connection:
    // close) made browsers stall on their per-host connection limit — the
    // dominant latency for the real UI, invisible to curl.
    loop {
        let accepted = Instant::now();
        let mut request_line = String::new();
        match reader.read_line(&mut request_line) {
            Ok(0) => return Ok(()),  // client closed the connection
            Ok(_) => {}
            Err(_) => return Ok(()), // idle timeout or read error: drop it
        }
        if request_line.trim().is_empty() {
            continue; // tolerate a stray blank line between requests
        }
        // Consume headers to the blank line; note an explicit close request.
        let mut client_close = false;
        loop {
            let mut line = String::new();
            let n = reader.read_line(&mut line)?;
            if n == 0 || line == "\r\n" || line == "\n" {
                break;
            }
            let l = line.to_ascii_lowercase();
            if l.starts_with("connection:") && l.contains("close") {
                client_close = true;
            }
        }
        let read_ms = accepted.elapsed().as_millis();
        let (method, target) = {
            let mut parts = request_line.split_whitespace();
            (
                parts.next().unwrap_or("").to_string(),
                parts.next().unwrap_or("/").to_string(),
            )
        };
        let (path, query) = match target.split_once('?') {
            Some((p, q)) => (p, q),
            None => (target.as_str(), ""),
        };

        let work = Instant::now();
        let routed = route(path, query, state);
        let route_ms = work.elapsed().as_millis();

        let keep_alive = !client_close;
        let r = match routed {
            Ok(resp) => write_response(&mut stream, 200, resp.ctype, &resp.body, keep_alive),
            Err(msg) => {
                write_response(&mut stream, 400, "text/plain; charset=utf-8", msg.as_bytes(), keep_alive)
            }
        };
        if path.starts_with("/api/") {
            eprintln!(
                "view: {method} {target} — read {read_ms}ms, route {route_ms}ms, total {}ms",
                accepted.elapsed().as_millis()
            );
        }
        r?;
        if client_close {
            return Ok(());
        }
    }
}

struct Resp {
    ctype: &'static str,
    body: Vec<u8>,
}

fn route(path: &str, query: &str, state: &Mutex<ViewState>) -> Result<Resp, String> {
    match path {
        "/" | "/index.html" => Ok(Resp {
            ctype: "text/html; charset=utf-8",
            body: PAGE.as_bytes().to_vec(),
        }),
        "/api/state" => {
            let s = state.lock().unwrap();
            Ok(Resp {
                ctype: "application/json",
                body: state_json(&s).into_bytes(),
            })
        }
        "/api/field" => {
            let i: usize = param(query, "hull")
                .ok_or("field: missing hull index")?
                .parse()
                .map_err(|_| "field: bad hull index".to_string())?;
            let s = state.lock().unwrap();
            let f = &s.hulls.get(i).ok_or("field: hull index out of range")?.field;
            let mut body = Vec::with_capacity(f.zeta.len() * 4);
            for &v in &f.zeta {
                body.extend_from_slice(&v.to_le_bytes());
            }
            Ok(Resp {
                ctype: "application/octet-stream",
                body,
            })
        }
        "/api/speed" => {
            let u: f64 = param(query, "u")
                .ok_or("speed: missing u")?
                .parse()
                .map_err(|_| "speed: bad value".to_string())?;
            if !(u.is_finite() && u > 0.0) {
                return Err("speed must be positive".into());
            }
            let mut s = state.lock().unwrap();
            s.cond.speed = u;
            recompute_fields(&mut s)?;
            Ok(Resp {
                ctype: "application/json",
                body: state_json(&s).into_bytes(),
            })
        }
        "/api/displacement" => {
            let mass: f64 = param(query, "mass")
                .ok_or("displacement: missing mass")?
                .parse()
                .map_err(|_| "displacement: bad value".to_string())?;
            if !(mass.is_finite() && mass > 0.0) {
                return Err("mass must be positive".into());
            }
            let mut s = state.lock().unwrap();
            set_displacement(&mut s, mass)?;
            Ok(Resp {
                ctype: "application/json",
                body: state_json(&s).into_bytes(),
            })
        }
        "/api/resistance" => {
            let s = state.lock().unwrap();
            let places = parse_placements(query, s.hulls.len())?;
            Ok(Resp {
                ctype: "application/json",
                body: resistance_json(&s, &places)?.into_bytes(),
            })
        }
        _ => Err(format!("no such path {path:?}")),
    }
}

fn parse_placements(query: &str, n: usize) -> Result<Vec<Placement>, String> {
    let raw = param(query, "p").ok_or("resistance: missing placements")?;
    let mut out = Vec::new();
    for pair in raw.split(';').filter(|s| !s.is_empty()) {
        let (xs, ys) = pair
            .split_once(',')
            .ok_or("resistance: placement must be x,y")?;
        let x: f64 = xs.parse().map_err(|_| "resistance: bad x".to_string())?;
        let y: f64 = ys.parse().map_err(|_| "resistance: bad y".to_string())?;
        out.push(Placement { x, y });
    }
    if out.len() != n {
        return Err(format!(
            "resistance: expected {n} placements, got {}",
            out.len()
        ));
    }
    Ok(out)
}

fn resistance_json(s: &ViewState, places: &[Placement]) -> Result<String, String> {
    let members: Vec<(&Hull, Placement)> = s
        .hulls
        .iter()
        .zip(places)
        .filter(|(vh, _)| !vh.dry)
        .map(|(vh, &p)| (&vh.hull, p))
        .collect();
    if members.is_empty() {
        return Ok("{\"total\":null,\"wave\":null,\"viscous\":null,\"interference\":null,\
                   \"cw\":null,\"ct\":null,\"effective_power\":null,\"froude\":null}"
            .to_string());
    }
    // Only the combined wave resistance depends on placement; integrate it on
    // the fixed grid (fast regardless of separation). Viscous, wetted surface,
    // and the solo-wave total are precomputed for the current speed/mass.
    let mut spec = FreeWaveSpectrum::new(&members, &s.cond).map_err(|e| format!("{e}"))?;
    let (sec_cap, n_theta) = resistance_grid(&s.cond, s.l_ref);
    let wave = fast_wave_resistance(&mut spec, sec_cap, n_theta);
    let total = wave + s.viscous_total;
    let u = s.cond.speed;
    let q = 0.5 * s.cond.fluid.density * u * u * s.wetted_surface;
    let interference = if s.solo_wave_total > f64::MIN_POSITIVE {
        wave / s.solo_wave_total
    } else {
        1.0
    };
    Ok(format!(
        "{{\"total\":{},\"wave\":{},\"viscous\":{},\"interference\":{},\
         \"cw\":{},\"ct\":{},\"effective_power\":{},\"froude\":{}}}",
        jnum(total),
        jnum(wave),
        jnum(s.viscous_total),
        jnum(interference),
        jnum(if q > 0.0 { wave / q } else { f64::NAN }),
        jnum(if q > 0.0 { total / q } else { f64::NAN }),
        jnum(total * u),
        jnum(s.cond.froude_number(s.l_ref)),
    ))
}

/// A finite f64 as a JSON number, else JSON `null` (JSON has no NaN/Inf, and
/// emitting a bare `NaN` would make the whole response unparseable — which
/// would silently freeze the live readout rather than degrade one field).
fn jnum(x: f64) -> String {
    if x.is_finite() {
        format!("{x}")
    } else {
        "null".to_string()
    }
}

fn state_json(s: &ViewState) -> String {
    let [x0, x1, y0, y1] = s.view;
    let nu = s.cond.gravity / (s.cond.speed * s.cond.speed);
    let mut out = String::new();
    out.push_str(&format!(
        "{{\"view\":{{\"x0\":{x0},\"x1\":{x1},\"y0\":{y0},\"y1\":{y1}}},\
         \"speed\":{},\"froude\":{},\"transverseWavelength\":{},\"vmax\":{},\
         \"mass\":{},\"designMass\":{},\"lRef\":{},\"margin\":{},\
         \"fadeToward\":[240,239,236],\"fadeFraction\":0.55,\"hulls\":[",
        s.cond.speed,
        s.cond.froude_number(s.l_ref),
        2.0 * std::f64::consts::PI / nu,
        s.vmax,
        s.mass,
        s.design_mass,
        s.l_ref,
        s.margin,
    ));
    for (i, vh) in s.hulls.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let f = &vh.field;
        out.push_str(&format!(
            "{{\"id\":{i},\"name\":\"{}\",\"homeX\":{},\"homeY\":{},\"x0\":{},\"x1\":{},\
             \"canFloat\":{},\"dry\":{},\"field\":{{\"lx0\":{},\"lx1\":{},\"ly0\":{},\"ly1\":{},\
             \"nx\":{},\"ny\":{}}},\"beam\":[",
            json_escape(&vh.name),
            vh.home.x,
            vh.home.y,
            vh.x0,
            vh.x1,
            vh.body.is_some(),
            vh.dry,
            f.lx0,
            f.lx1,
            f.ly0,
            f.ly1,
            f.nx,
            f.ny,
        ));
        for (j, (x, b)) in vh.beam.iter().enumerate() {
            if j > 0 {
                out.push(',');
            }
            out.push_str(&format!("[{x},{b}]"));
        }
        out.push_str("]}");
    }
    out.push_str("]}");
    out
}

fn json_escape(s: &str) -> String {
    let mut o = String::new();
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o
}

fn param(query: &str, key: &str) -> Option<String> {
    query.split('&').find_map(|kv| {
        let (k, v) = kv.split_once('=')?;
        if k == key {
            Some(url_decode(v))
        } else {
            None
        }
    })
}

fn url_decode(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'%' if i + 2 < b.len() => {
                let h = hex(b[i + 1]).zip(hex(b[i + 2]));
                if let Some((hi, lo)) = h {
                    out.push(hi * 16 + lo);
                    i += 3;
                    continue;
                }
                out.push(b[i]);
                i += 1;
            }
            b'+' => {
                out.push(b' ');
                i += 1;
            }
            c => {
                out.push(c);
                i += 1;
            }
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn hex(c: u8) -> Option<u8> {
    match c {
        b'0'..=b'9' => Some(c - b'0'),
        b'a'..=b'f' => Some(c - b'a' + 10),
        b'A'..=b'F' => Some(c - b'A' + 10),
        _ => None,
    }
}

fn write_response(
    stream: &mut TcpStream,
    status: u16,
    ctype: &str,
    body: &[u8],
    keep_alive: bool,
) -> std::io::Result<()> {
    let reason = if status == 200 { "OK" } else { "Bad Request" };
    let conn = if keep_alive { "keep-alive" } else { "close" };
    let header = format!(
        "HTTP/1.1 {status} {reason}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\n\
         Cache-Control: no-store\r\nConnection: {conn}\r\n\r\n",
        body.len()
    );
    stream.write_all(header.as_bytes())?;
    stream.write_all(body)?;
    stream.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use michell::hulls;

    /// The viewer composites the fleet field client-side by translating and
    /// summing each hull's *standalone* field. That must reproduce the fleet
    /// field computed directly with all hulls in one spectrum (thin-ship
    /// superposition), and it must do so with this crate's local-coordinate
    /// sampling convention — the correctness foundation of the whole viewer.
    #[test]
    fn superposition_matches_direct_fleet_field() {
        let hull = hulls::wigley(12.0, 1.2, 0.75).unwrap();
        let cond = Conditions::seawater(4.0);
        let places = [
            Placement { x: 0.0, y: 2.5 },
            Placement { x: -3.0, y: -2.5 },
        ];

        // Direct: both hulls in one spectrum.
        let members: Vec<(&Hull, Placement)> = places.iter().map(|&p| (&hull, p)).collect();
        let mut direct = FreeWaveSpectrum::new(&members, &cond).unwrap();

        // Superposition: each hull alone at the origin, evaluated at the
        // placement-shifted point and summed — exactly what the browser does.
        let mut solos: Vec<FreeWaveSpectrum> = places
            .iter()
            .map(|_| FreeWaveSpectrum::new(&[(&hull, Placement { x: 0.0, y: 0.0 })], &cond).unwrap())
            .collect();

        let mut max_err = 0.0f64;
        let mut max_mag = 0.0f64;
        for &x in &[-30.0, -20.0, -12.0, -6.0] {
            for &y in &[-6.0, -1.0, 0.0, 3.0, 7.0] {
                let truth = direct.elevation_at(x, y).unwrap();
                let mut sum = 0.0;
                for (solo, p) in solos.iter_mut().zip(&places) {
                    sum += solo.elevation_at(x - p.x, y - p.y).unwrap();
                }
                max_err = max_err.max((truth - sum).abs());
                max_mag = max_mag.max(truth.abs());
            }
        }
        assert!(
            max_err <= 1e-6 * max_mag.max(1e-9),
            "superposition drifted from the direct fleet field: \
             max_err={max_err:e}, field_mag={max_mag:e}"
        );
    }

    #[test]
    fn bilinear_sample_reads_grid_corners_and_center() {
        // A 3x3 field with known values; check exact hits and a midpoint.
        let f = Field {
            lx0: 0.0,
            lx1: 2.0,
            ly0: 0.0,
            ly1: 2.0,
            nx: 3,
            ny: 3,
            zeta: vec![0.0, 1.0, 2.0, 10.0, 11.0, 12.0, 20.0, 21.0, 22.0],
        };
        assert!((sample_field(&f, 0.0, 0.0) - 0.0).abs() < 1e-12);
        assert!((sample_field(&f, 2.0, 2.0) - 22.0).abs() < 1e-12);
        assert!((sample_field(&f, 1.0, 1.0) - 11.0).abs() < 1e-12);
        // Between (0,0)=0 and (1,0)=1 at x=0.5 → 0.5.
        assert!((sample_field(&f, 0.5, 0.0) - 0.5).abs() < 1e-12);
        // Outside the grid is zero.
        assert_eq!(sample_field(&f, -0.1, 0.0), 0.0);
        assert_eq!(sample_field(&f, 0.0, 2.1), 0.0);
    }
}

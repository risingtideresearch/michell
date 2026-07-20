# `michell` — thin-ship wave resistance from B-spline hulls

A pure-Rust, zero-dependency library computing **wave resistance by Michell's
integral** and **viscous resistance by the ITTC-57 line** for hulls described
as clamped, polynomial (non-rational) tensor-product **B-spline half-breadth
surfaces**. Companion to [`resistance`](../resistance) (Holtrop–Mennen /
Savitsky blend), which deliberately excludes shape-sensitive wave methods —
this crate is that missing shape diagnostic.

```rust
use michell::{hulls, Conditions};

let hull = hulls::wigley(10.0, 1.0, 0.625)?;   // L, B, T [m]
let cond = Conditions::seawater(3.0);          // U [m/s]
let r = michell::resistance(&hull, &cond)?;
println!("Rw = {:.1} N, Rv = {:.1} N, Cw = {:.4e}",
         r.wave.resistance, r.viscous.resistance, r.cw);
```

## Geometry contract

The hull is port/starboard symmetric, given by its half-beam `y = f(x, z) ≥ 0`
as a B-spline surface:

- `x` runs along the hull (arbitrary origin), **`z` runs vertically downward**
  from the undisturbed waterline; domain `z ∈ [0, T]`. SI units throughout.
- Knot vectors must be **clamped** (end knots repeated degree+1 times);
  interior knots may repeat up to the degree — a multiplicity-`degree` knot in
  `z` is how you represent a **chine** exactly.
- **Polynomial only** (all NURBS weights = 1). This is a deliberate contract:
  polynomial spans are what allow the Michell inner integrals to be evaluated
  in closed form (below). CAD hull surfaces are almost always unit-weight.
- The control net must be non-negative (sufficient for `f ≥ 0` by the
  convex-hull property).

## Theory and numerics

With `ν = g/U²` and half-beam `f` (Tuck 1989; Dambrine–Pierre–Rousseaux 2016):

```
R_w = (4 ρ g²)/(π U²) ∫₁^∞ (I² + J²) λ²/√(λ² − 1) dλ
I(λ) + i J(λ) = ∬ (∂f/∂x) · exp(−ν λ² z) · exp(i ν λ x) dx dz
```

The numerical strategy exploits the spline structure end to end:

1. On every knot span, `∂f/∂x` is an exact local polynomial (extracted once
   per hull via corner Taylor data). The inner integrals then reduce to
   closed-form moments `∫ tᵃ e^{ikt} dt` and `∫ tᵇ e^{−κt} dt` (series for
   small argument to avoid cancellation, stable recurrences otherwise) — **I
   and J carry no quadrature error at all**.
2. The outer integral is transformed by `λ = sec θ`, removing the `λ = 1`
   singularity. It is integrated by 16-point Gauss–Legendre panels sized to
   the local oscillation rate, truncated only when a full multi-period window
   of accumulated phase contributes negligibly, and refined (panel halving)
   until a requested relative tolerance is met, with the achieved estimate
   reported in the result.

Viscous resistance: ITTC-57 `C_F = 0.075/(log₁₀Re − 2)²` with optional form
factor, on the thin-ship wetted surface `S = 2∬√(1 + fx² + fz²) dx dz`.

**Multihulls**: thin-ship far-field amplitudes superpose, so hull `j` placed
at longitudinal offset `Δx_j` and transverse position `y_j` contributes
`F_j(λ)·exp(iν(λΔx_j + λ√(λ²−1) y_j))` and the resistance uses `|Σ_j …|²` —
wave interference is exact within the theory (for a catamaran this reduces to
the classical `4cos²(½νsλ√(λ²−1))` factor). Viscous resistance sums per
member. `multihull_resistance` also reports the interference factor
(combined R_w / Σ standalone R_w). Demihulls must each be symmetric about
their own centerplane.

Performance: a full 21-speed Wigley resistance curve at the default 1e-5
tolerance runs in ~20 ms (release build).

## Validation

`cargo test` checks, among others:

- B-spline evaluation, derivatives, and per-span polynomial extraction against
  closed forms (single- and multi-span, with chine knots);
- moment integrals against brute-force quadrature across both branches;
- `I(λ), J(λ)` against the **analytic Wigley-hull expressions** to ~1e-10;
- total `R_w` against an independent 2-million-point Simpson reference;
- Wigley displaced volume against the exact `4BLT/9`, wetted surface against a
  dense reference, and the Cw(Fn) curve for the expected humps and hollows.

## Assumptions and limitations

- Michell linearisation: slender hull (`|∂f/∂x| ≪ 1`), no sinkage/trim, no
  breaking, deep water, infinite domain, monohull.
- Half-breadth should close at both ends; **transom sterns are not yet
  modelled** (planned: virtual closure).
- Viscous model is a flat-plate correlation; supply your own form factor.

## Input front-ends

Every front-end reduces to the same intermediate representation
(`grid::SampleGrid`): a station × waterline grid of half-beam samples,
optionally augmented with `∂f/∂x`/`∂f/∂z` channels (`NaN` = unknown at that
sample) and per-sample weights (`0` excludes a sample — e.g. a failed CAD
inversion, which is *unknown* geometry rather than zero beam). The grid is
lofted to the spline by `fit::fit_grid`: weighted tensor-product least
squares over every channel present (quantile knot placement, per-channel
residuals reported so you can judge fit quality). Derivative observations
are scaled by the local sample spacing so slopes and values are
commensurate, and a slope that predicts more change across one sample cell
than the half-beam anywhere on its stencil is skipped — it describes
geometry (a bilge wall, the keel fold) that no loft at that sampling can
resolve, and fitting it would only distort the values.
`FitOptions::derivative_weight` tunes or disables the channels.

- **Offsets**: `fit::fit_offsets(stations, waterlines, half_beams, opts)` —
  the human-authorable path: a station × waterline table of half-beams
  (a value-only grid).
- **STL**: `stl::mesh_fleet(bytes, units_scale, waterline)` — binary or ASCII
  triangle meshes. Half-beams are extracted by transverse **ray casting**
  (exact, no Newton iteration; empty results are the footprint), hulls
  cluster by shared vertices + wetted proximity, and everything downstream
  (folding, bodies, sweeps) is shared with IGES. STL has no units field, so
  `--units mm|cm|m|in|ft` is required. Quality tracks the export's chord
  tolerance: fine CAD tessellations match IGES; decimated meshes add
  geometry noise that wave resistance is sensitive to.
- **IGES**: `iges::import_hull(text, opts)` — reads **one or many untrimmed**
  NURBS patches (entity 128, unit weights; 124 transforms and unit conversion
  handled; naturally-bounded 143/141 wrappers, as produced by SubD → NURBS
  exports, are tolerated). Patches are mapped to the hull frame
  (`waterline_z` picks the DWL), the hull **centerplane is auto-detected** —
  a full both-sided shell is folded about the midplane of its shell
  intersections, a half hull measures from y = 0, and `centerplane` overrides
  either — then the wetted region is sampled by per-patch Newton inversion
  (outermost fold wins) and lofted; the surface slopes `∂y/∂x`, `∂y/∂z` come
  for free from the converged Newton Jacobian (implicit function theorem)
  and join the loft as derivative observations. Genuinely trimmed surfaces
  (142/144) and rational weights are rejected with specific messages; fold
  asymmetry, ambiguous samples, failed inversions, slope gaps, and loft
  residuals (with the location of the worst one) are reported so a bad
  import is visible. Handles both the
  "export selected surface" workflow and multi-patch SubD hull exports,
  including hulls modelled off-centre (e.g. an ama in position).

## CLI

The `michell-cli` crate builds a `michell` binary (`cargo build --release`,
binary at `target/release/michell`):

```text
michell wigley -o wigley.hull                       # reference hull
michell info wigley.hull                            # geometry & diagnostics
michell resistance wigley.hull --froude 0.2:0.5:0.05
michell resistance hull.igs --waterline 2.6 --speeds 4:9:0.5 --knots --json
michell resistance vaka.hull ama.igs@y=1.9 ama.igs@y=-1.9 --speeds 3:8:0.5
michell loft table.offsets -o hull.hull             # offsets -> control net
michell place ama.igs@dy=1.7,dz=0.05 ama.igs@dy=-1.7 -o boat.igs  # posed CAD geometry
michell spectrum wigley.hull --speed 3              # free-wave spectrum (CSV)
michell wake boat-*.hull --speed 8 --knots -o wake.png   # Kelvin wake heatmap
```

**Wave field** (`spectrum`, `wake`; one speed via `--speed` or `--froude`):
the far-field wave pattern is reconstructed from the same exactly-evaluated
amplitude function `F = I + iJ` the resistance uses. `spectrum` tabulates the
free-wave spectrum by propagation angle θ — elevation amplitude density |A(θ)|
[m/rad], phase, and the angular resistance density dR_w/dθ, whose integral
reproduces R_w (cross-checked on stderr) — showing where the wave energy goes
(transverse θ ≈ 0 vs diverging θ → ±90°) and which angles a multihull's
interference cancels. `wake` evaluates the Kelvin pattern

```text
ζ(x, y) = Re ∫ A(θ) e^{iν secθ (x + y tanθ)} dθ,   A = −(2ν/π) sec³θ conj(F)
```

on a grid and renders a PNG heatmap (blue trough / red crest, hull
waterplanes in gray), or emits CSV/JSON for other tooling. Conventions: the
ship advances toward +x, so the wake trails toward −x; the reconstruction is
the far-field free-wave part of the linear solution, physical astern of each
hull (not on or ahead of it). Grids too coarse for the shortest diverging
waves are smoothly band-limited and flagged. The magnitude of A is pinned by
the deep-water free-wave resistance identity R_w = ½πρU² ∫|A|²cos³θ dθ; the
phase follows Tuck, Scullen & Lazauskas mapped to these conventions.

**Sweeps** (`michell sweep study.json`): long-form CSV/JSON over the
Cartesian product of axes — every varying quantity (speed, weight, lcg,
waterline, hull poses) is an axis, with fixed values as single-valued axes.
A `weight` axis switches to **equilibrium mode**: each point's platform
sinkage (and pitch, with `lcg`) is solved by a Newton iteration whose
Jacobian comes from the waterplane properties, so counterfactuals like
"what if the boat were heavier / the CG further forward" are swept at
physically consistent attitudes. Every record carries the solved state,
displacement, LCB, dry-hull count, and the resistance breakdown. All poses
are hydrostatic (no speed-dependent squat).

```json
{
  "name": "ama placement study",
  "fluid": "seawater",
  "hulls": [
    { "id": "vaka",  "file": "boat-center.hull" },
    { "id": "ama_s", "file": "boat-starboard.hull" },
    { "id": "ama_p", "file": "boat-port.hull", "pose": { "trim": 0.5 } }
  ],
  "sweep": [
    { "target": "speed", "unit": "knots", "range": [4, 10], "step": 0.5 },
    { "target": "weight", "range": [1800, 2600], "step": 200 },
    { "target": "lcg", "value": -5.8 },
    { "target": ["ama_s", "ama_p"], "param": "spread", "range": [1.5, 2.5] },
    { "target": "ama_s", "param": "trim", "values": [-2, 0, 2] }
  ],
  "output": { "format": "csv", "file": "study.csv" },
  "options": { "rel_tol": 1e-5, "form_factor": 0.05 }
}
```

Axis values: `range: [start, stop]` with optional `step` (default: a fifth
of the span), `values: [...]`, or scalar `value`. Speed axes take `unit`
(`ms` | `knots` | `froude`). Pose params: `dx`, `dy`, `dz` (+down),
`spread` (outboard, sign follows each hull's side), `trim` (degrees,
+ raises the +x end); a target list moves several hulls as one coupled
axis. Hull files load relative to the manifest. A flag-based sweep over raw
IGES (`--axis`, `--float`) remains for one-liners.

**GZ curves**: a `heel` axis (degrees, + puts the +y side down) heels the
platform rigidly about the centerline at the design floatplane and re-solves
the equilibrium at every angle, so the displaced volume is held while
buoyancy transfers between hulls — the windward hull flying shows up in the
`dry` column, and resistance is computed on the heeled fleet. It requires a
`weight` axis and a `vcg` axis (centre of gravity in metres above the design
floatplane — itself sweepable for KG studies); with `vcg` present every row
carries `gz` (righting arm, m; positive rights the boat) and `rm` (righting
moment, N·m). The inter-hull buoyancy transfer — the dominant multihull
mechanism — is exact; each hull's *own* heel cannot be represented by a
symmetric half-breadth surface and enters metacentrically, as
`sin φ·(I_T/∇ − KB)` per hull, so a single slender monohull reduces to
`GZ = GM_T·sin φ` and hard-chine form stability at large heel is
underestimated.

```json
  "sweep": [
    { "target": "speed", "unit": "knots", "value": 8 },
    { "target": "weight", "value": 2200 },
    { "target": "vcg", "value": 1.1 },
    { "target": "heel", "range": [-15, 15], "step": 1 }
  ]
```

**Bodies**: sweep manifests reference **full-band** `.hull` files — the
half-breadth spline over the hull's band from keel to above the design
waterline, written by `michell loft`:

```text
michell loft boat.igs --waterline 0.42 -o boat
  -> boat-port.hull  boat-center.hull  boat-starboard.hull
```

Each body records `waterline` (design WL depth below its band top) and
`centerplane` (detected transverse position, its default placement); hulls
are named by role when the layout is recognizable. The band reaches
`--band` metres above the design WL (default: half the design draft) — decks
are deliberately excluded, since a deck is a cliff for a height-field loft;
poses that rise past the band are reported per record (`band_exceeded`).
Re-situating a body needs no Newton inversion (pitch rotation of a height
field is an exact domain reparametrization), so equilibrium points run
seconds-fast; the loft itself defaults to a dense net (28x32 at 241x97)
because wave resistance is sensitive to loft resolution near the keel
rocker and the body is fit once, reused thousands of times.

**Export back to CAD** (`michell place`): once a sweep has found a good
configuration — amas at a chosen `dx`/`dy`/`dz`, a solved sinkage and trim —
`place` writes that posed geometry as a new IGES file to import back into
CAD:

```text
michell place ama.igs@dy=3.9,dz=-0.008 ama.igs@dy=-0.1,dz=-0.008 \
  --waterline 0 --sinkage -0.008 -o boat.igs
```

Each spec's pose (`dx`/`x`, `dy`, `y` absolute centerplane, `dz` immersion,
`trim` degrees about `pivot`) applies rigidly to every hull in that file, and
`--sinkage`/`--platform-trim`/`--pivot-x` apply the whole-platform state a
`weight`/`lcg` equilibrium row reports. IGES inputs pass their surfaces
through **exactly** (untrimmed 128 patches, unit weights, metres; bounded
143/141 bases are written as full surfaces with the parameter range
restricted to the bounded box). `.hull` control nets and full-band bodies
convert **exactly** too: the half-breadth graph `y = ±f(x, z)` is a B-spline
surface whose control net sits on the graph's coordinate lines (linear
precision at the Greville abscissae), emitted as the mirrored port/starboard
pair. Platform sinkage is re-expressed as the hulls moving down, so the water
surface stays at `--waterline` in the output frame. The written file
round-trips through the importer: re-importing the reconstructed boat above
reproduces the study's resistance to within the loft tolerance.

Hydrostatics on every hull: displaced volume, LCB, KB, waterplane area and
moments (longitudinal and transverse), LCF — exact spline integrals.

Multihulls: list several hulls, each with an optional placement suffix —
`@y=Y` places a (single-hull file's) centerplane absolutely, `@dy=S` shifts
transversely, `@x=DX`/`@dx=DX` shifts longitudinally. An IGES file holding a
whole multihull imports as a fleet automatically: hulls are detected by
clustering wetted patches, each at its detected centerplane, and dry
structure (beams, decks) is dropped. The `IF` column / `interference` JSON
field reports combined R_w over the sum of standalone R_w.

All input kinds are accepted everywhere (sniffed by header/extension):
the canonical `.hull` control net, a `michell-offsets v1` station × waterline
table (lofted on load), a `*.grid.json` sample grid, and IGES (sampled +
lofted; `--waterline` sets the DWL). `--json` gives machine-readable output;
`--fluid`, `--rho`, `--nu`, `--form-factor`, `--rel-tol` control the physics.
Conversion diagnostics (loft residuals, mirroring, failed inversions) are
always printed so a bad import can't pass silently.

`--dump-grid PATH` writes the sample grid a load produced — stations,
waterlines, half-beams, slope channels, weights — as `*.grid.json`
(multihull files get `-0`, `-1`, ... suffixes), so you can inspect or diff
exactly what the importer sampled, and re-loft it later without the source
CAD file. `--fit-deriv-weight W` scales the slope observations (0 = fit
values only).

**`.hull` format** (canonical, SI, `#` comments): `michell-hull v1`,
`degree-x/z`, `knots-x/z`, then one `row` of control values per x index.
**Offsets format**: `michell-offsets v1`, a `waterlines` line (z downward
from DWL, starting 0), then `station <x> <half-beams...>` lines.
**Grid format**: JSON with `"michell": "sample-grid"`, `stations`,
`waterlines`, row-major `half_beams` (waterline index fastest), optional
`dfdx`/`dfdz` (JSON `null` = unknown at that sample), optional `weights`,
optional `centerplane`.

## API sketch

- `BSplineSurface::new(degree_x, degree_z, knots_x, knots_z, control)` —
  validated surface (control row-major, z fastest).
- `Hull::new(surface)` — validates hull semantics, precomputes span
  polynomials, wetted surface, displaced volume.
- `Conditions::seawater(u)` / `::freshwater(u)` / custom `Fluid`.
- `wave_resistance[_with]`, `viscous_resistance[_with]`,
  `resistance[_with]` → forces, effective power P_E = R_t·U, coefficients,
  quadrature diagnostics.
- `multihull_resistance[_with]`, `multihull_wave_resistance[_with]`,
  `Placement` — fleets with exact wave interference.
- `iges::import_fleet` — every hull in a file, with detected placements;
  `iges::import_hull` — exactly one (errors on multihull files).
- `iges::write` — serialize surfaces to an IGES 5.3 file (untrimmed 128
  patches, metres); `iges::halfbreadth_surfaces` — the exact mirrored
  surface pair of a half-breadth spline; `iges::apply_pose` /
  `SourceFleet::posed_surfaces` — pose CAD geometry for re-export.
- `iges::source_fleet` + `SourceFleet::situate[_one](waterline, poses,
  platform)` — re-situate hulls repeatedly (immersion, mount trim, position).
- `body::Body` — full-band half-breadth spline; `situate` via exact
  reparametrization (no Newton), fast enough for solver inner loops.
- `float::solve_equilibrium[_bodies|_with]` — hydrostatic sinkage/pitch
  balance for a mass + LCG load case, over IGES fleets, body assemblies, or
  any custom situate closure.
- `float::heel_poses` + `float::righting_arm` — rigid platform heel and the
  GZ of the re-solved fleet (exact inter-hull transfer, metacentric per-hull
  term).
- `inner_integrals(hull, cond, λ)` — free-wave amplitude functions.
- `FreeWaveSpectrum::new(&members, &cond)` — far-field spectrum of a fleet:
  `amplitude(θ)`, `resistance_density(θ)` (dR_w/dθ), and Kelvin-wake
  reconstruction via `elevation_at(x, y)` / `elevation_grid(...)`.
- `hulls::wigley(l, b, t)` — exact reference hull.

## Roadmap

1. Parallel sweep evaluation (each equilibrium point is independent).
2. STEP reader feeding the same sample-and-loft pipeline; OBJ via the mesh
   path.
3. Transom closure, Python bindings.
4. Longitudinal wave cuts against published Wigley measurements; wake
   animation over a speed range.

## References

- E. O. Tuck, *The wave resistance formula of J.H. Michell (1898) and its
  significance to recent research in ship hydrodynamics*, J. Austral. Math.
  Soc. B 30 (1989).
- E. O. Tuck, D. C. Scullen & L. Lazauskas, *Ship-wave patterns in the
  spirit of Michell*, IUTAM Symposium (2001); *Wave patterns and minimum
  wave resistance for high-speed vessels*, 24th Symp. Naval Hydrodynamics
  (2002) — far-field free-wave spectrum and wave-pattern evaluation.
- J. Dambrine, M. Pierre, G. Rousseaux, *A theoretical and numerical
  determination of optimal ship forms based on Michell's wave resistance*,
  ESAIM: COCV (2016), arXiv:1410.2800.
- ITTC — Recommended Procedures: *1957 ITTC Performance Prediction Method*.

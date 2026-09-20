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

The hull is, by default, port/starboard symmetric, given by its half-beam
`y = f(x, z) ≥ 0` as a B-spline surface (asymmetric hulls — port ≠ starboard —
are supported too; see [Asymmetric hulls](#asymmetric-hulls) below):

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

Viscous resistance follows the ITTC-78 shape on the thin-ship wetted surface
`S = 2∬√(1 + fx² + fz²) dx dz`:

```text
C_V = (1 + k)·C_F(Re) + ΔC_F,      C_F = 0.075/(log₁₀Re − 2)²
```

The **form factor** `k` multiplies flat-plate friction — it is a property of
the *shape* (streamline curvature speeds the flow over most of the hull, and
the stern boundary layer costs a viscous pressure defect). The **roughness
allowance** `ΔC_F` is added *outside* it, because surface finish is a property
of the *skin*. They are kept apart deliberately: they come from different
places (one from geometry, one from the paint), and on a small hull the
roughness term can be the larger of the two, so folding it into `k` hides it.
Both default to zero, which is the bare flat plate.

`ΔC_F` can be given directly (`--roughness cf=4e-4`, which is also where
ITTC-78's correlation allowance `C_A` goes) or estimated from an equivalent
sand-grain height (`--roughness ks=150um`). The estimate is the textbook Moody
construction — the excess of Prandtl–Schlichting fully-rough friction over the
smooth line, floored at zero — so it returns **exactly** zero while the surface
is hydraulically smooth at that Reynolds number, and grows with speed as the
smooth line falls away beneath the Re-independent rough one. It reports the
roughness Reynolds number `k_s⁺ = k_s·u_τ/ν` with it, because that is what says
whether the finish matters at all (`≲ 5` smooth, `≳ 70` fully rough) and how
far to trust the number: the transitional band between is *bridged by the
crossover*, not fitted, so it reads high in there. Published ship correlations
are no help at this size — on an 8 m hull at `Re ≈ 9×10⁶`, Bowden–Davison
returns about 60% of `C_F` and Townsin returns a negative number.

**Multihulls**: thin-ship far-field amplitudes superpose, so hull `j` placed
at longitudinal offset `Δx_j` and transverse position `y_j` contributes
`F_j(λ)·exp(iν(λΔx_j + λ√(λ²−1) y_j))` and the resistance uses `|Σ_j …|²` —
wave interference is exact within the theory (for a catamaran this reduces to
the classical `4cos²(½νsλ√(λ²−1))` factor). Viscous resistance sums per
member. `multihull_resistance` also reports the interference factor
(combined R_w / Σ standalone R_w). Demihulls may individually be asymmetric
(see below); their dipole systems superpose with the source systems.

Performance: a full 21-speed Wigley resistance curve at the default 1e-5
tolerance runs in ~20 ms (release build). The outer quadrature and the
near-field (sinkage/trim) integrals fan their independent nodes out across
the machine's cores with `std::thread::scope` (still no dependencies); the
reduction order is the serial one, so the answer is bit-for-bit independent
of the thread count. The worker budget is per thread (`michell::parallel`),
so a caller that already runs jobs in parallel can hand each job a share of
the cores — the manifest sweep does — and `MICHELL_THREADS=1` in the
environment disables threading altogether.

### Asymmetric hulls

A hull whose two sides differ — starboard `y = +f₊(x, z)`, port `y = −f₋(x, z)`
— is built with `Hull::new_asymmetric(port, starboard)`. It is split into a
**symmetric thickness** part `f_sym = (f₊ + f₋)/2` and an **antisymmetric
camber** part `f_a = (f₊ − f₋)/2` on a shared parametrisation (the two surfaces
must share degrees and knot vectors; only their control nets differ):

```
R_w = R_source(f_sym) + R_dipole(f_a)
```

The thickness part is the classical Michell **source** system. The camber part
adds a centreplane **y-dipole** (normal-doublet) system: the same inner
integral over `∂f_a/∂x`, weighted by the transverse wavenumber factor
`√(λ²−1)`. Because the source amplitude is even in the wave angle θ and the
dipole amplitude is odd, their cross term integrates to zero over the Kelvin
fan — the two resistances add with no interference, and a symmetric hull
(`f_a ≡ 0`) reproduces classical Michell to full floating-point precision.

### Transom sterns

Michell's integral is over a **closed** body — `∬(∂f/∂x)…` assumes the
half-breadth reaches zero at both ends. A transom leaves a step there, and
taking that step at face value models a hull that shuts instantaneously, which
radiates far too much. A real ventilated transom instead lets the flow leave
the edge cleanly and close in a hollow some way downstream.

A **virtual appendage** supplies that hollow: running aft from the transom over
a length `L_v`, the half-beam decays as `f_v(x, z) = f_T(z)·φ(s)` with
`s = (x_T − x)/L_v` and the smoothstep `φ(s) = (1 − s)²(1 + 2s)`, which is flat
at both ends so `f` stays continuous at the transom and closes tangentially.
Substituting `x = x_T − s·L_v` separates the amplitude completely:

```text
F_v = −e^{iνλx_T} · ∫ f_T(z) e^{−κz} dz · ∫₀¹ φ'(s) e^{−iνλ L_v s} ds
```

The transom section `f_T` is a piecewise polynomial on the hull's own z-spans,
so the first factor is the *same* per-span z-moments the hull uses and the
second is a three-term oscillatory moment — the closure is exact too, at the
cost of one dot product per λ, and `L_v → 0` reproduces the bare step
analytically.

The hollow length defaults to the **ballistic** estimate `L_v = √2·U·√(d_T/g)`
— water leaving the transom horizontally at `U` falls the transom depth `d_T`
under gravity — and is tunable (`--transom ballistic=C`, or `hollow=METRES`
for a fixed length; `off` restores classical Michell). `Hull::transom()`
reports the transom's immersed area, waterline beam, and equivalent-rectangle
depth, and `michell info` prints `A_T/A_X` so a wet transom is never silent.
**For a hull that closes aft the whole mechanism is inert**, bit for bit.

Two caveats. The hollow length is a modelling choice, not a derived quantity,
so transom results inherit that uncertainty — sweep `ballistic=C` to see how
much it matters. And a lofted half-breadth cannot hold the transom's sharp
lower edge: the fit smears it downward and overstates `A_T` (by ~38% on a test
case), which feeds straight into `f_T`.

### Dynamic sinkage and trim

Every hydrostatic result above holds the hull at its still-water attitude.
`michell::squat` adds the missing piece: the steady near-field pressure's
**vertical force and pitch moment**, so a sweep can float the platform at its
*dynamic* attitude at speed rather than its at-rest one.

Write the Kelvin source in 2-D Fourier form. With the stream toward −x and `z`
down, the free-surface condition fixes the image amplitude
`A(k_x, k) = (k_x² + νk)/(k_x² − νk)` — rigid-wall (−1) for long modes,
free-surface (+1) for short ones, with the steady-wave dispersion curve
`k = k_x²/ν` the pole between them. Integrating the linearised pressure over
the hull and reducing by parts along `x` gives the force and moment as
wavenumber integrals of the **same per-span transforms** the wave integral
already evaluates exactly — one extra transform of `f` itself alongside
`∂f/∂x`, contracted against the existing z-moments (`InnerIntegral::contract_z`
/ `transforms_at`), so a hull with hundreds of spans stays affordable. Two
structural facts do the rest: the radiation condition's half-residue is odd in
`k_x` and so drops out of the force (even integrand) but survives in the
moment — **sinkage is a local-field effect, trim is mostly a wave effect** —
and the unbounded-fluid part of the source (`−1/r`) contributes no force at
all (d'Alembert), only a Munk-type moment for a fore-aft asymmetric hull.

```text
F_up  = −(ρU²/2π²) PV∬ d²k [ A·Re(q̄q) − ((A−1)/k)·Re(w̄q) ]
M_res = −(4ρU²ν/π) ∫₁^∞ dλ λ/√(λ²−1) · Im[ νλ² r̄q − r̄_wl q ]   on k = k_x²/ν
```

Validated three ways: an independent from-scratch derivation panel checked the
formulation (Rankine/wave-term split, radiation condition, low-Froude sign)
against Havelock (1939), Yeung (1972), and the Wigley sinkage/trim literature;
an independent NumPy oracle implementing the same integrals with adaptive
quadrature agrees with this crate's evaluation to **<0.05% on force** and
**~0.1% on moment** across Fn 0.05–0.45; and on the Wigley hull the sign,
Fn²-scaling (`s/L/Fn² ≈ 0.021`–`0.032`, bracketing Havelock's rigid-wall
ellipsoid limit and the Neumann–Michell/experimental range), and the trim
sign-reversal near Fn ≈ 0.34–0.36 all reproduce the published pattern.

`squat::multihull_dynamic_force` gives the force/moment for a fleet directly
(demihull interaction included, through the same placement phase the wave
superposition uses); `squat::dynamic_load_closure` adapts it to
`float::solve_equilibrium_dynamic_with` / `solve_equilibrium_bodies_dynamic`,
which balance it against buoyancy in the same Newton loop as the hydrostatic
solver (`DynamicLoad`/`DynamicEquilibrium` — bit-for-bit the hydrostatic
solver when the dynamic load is zero). `michell squat <hull>...` reports the
force, moment, lift fraction, and first-order equivalent sinkage/trim
directly; `options.dynamic: true` in a sweep manifest re-solves equilibrium at
**every speed** (attitude now depends on `U`), warm-started from the previous
speed's solution.

Two things to know before using it. It does not yet compose with the transom
closure's own
appendage moment contribution beyond what the hollow's transforms already
carry — both are scoped out for now, and a manifest sweep rejects the
combination rather than silently ignoring it. And thin-ship theory overstates
sinkage/trim by the same 20–40% it overstates wave resistance by, at the same
Fn 0.3–0.4 range, against surface-panel (Neumann–Michell) linear theory; the
linearisation itself expires once the dynamic force is a real share of the
weight (`DynamicForce::lift_fraction`/`DynamicEquilibrium::lift_fraction`
report exactly that, so the boundary is visible rather than silent). The
near-field quadrature is also markedly more expensive than the wave integral
— tens of thousands of transform evaluations per force/moment call on a
finely-lofted hull — so a dynamic sweep costs real wall-clock time per point;
budget accordingly, especially before a large speed × load grid.

## Validation

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

- Michell linearisation: slender hull (`|∂f/∂x| ≪ 1`), no breaking, deep
  water, infinite domain, monohull. Sinkage and trim are hydrostatic by
  default (the still-water attitude); `michell::squat` / `options.dynamic`
  (see below) solve the speed-dependent attitude instead, at real
  computational cost and only up to the same linearisation.
- Half-breadth should close at the bow. A **transom stern** is closed by a
  virtual appendage (`--transom`, `TransomClosure`); its hollow length is a
  modelling choice, so transom-sterned results carry that uncertainty.
- Viscous model is a flat-plate correlation plus two allowances you supply:
  a form factor `k` and a roughness `ΔC_F` (see above). Neither is derived
  from the hull — `k` from a regression such as Holtrop–Mennen or a
  double-body solve, `ΔC_F` from the finish. Do **not** back `k` out of a
  measured `C_T` using this crate's `C_W`: thin-ship theory overstates `C_W`,
  and a fit would quietly absorb that error into `k`.

## Input front-ends

Every front-end reduces to the same intermediate representation
(`grid::SampleGrid`): a station × waterline grid of half-beam samples,
optionally augmented with `∂f/∂x`/`∂f/∂z` channels (`NaN` = unknown at that
sample) and per-sample weights (`0` excludes a sample — e.g. a failed CAD
inversion, which is *unknown* geometry rather than zero beam). The grid is
lofted to the spline by `fit::fit_grid`: weighted tensor-product least
squares over every channel present (quantile knot placement, banded
Cholesky, per-channel residuals reported so you can judge fit quality). Derivative observations
are scaled by the local sample spacing so slopes and values are
commensurate, and a slope that predicts more change across one sample cell
than the half-beam anywhere on its stencil is skipped — it describes
geometry (a bilge wall, the keel fold) that no loft at that sampling can
resolve, and fitting it would only distort the values.
`FitOptions::derivative_weight` tunes or disables the channels.

**Loft resolution matters, and is cheap.** The normal equations of a
tensor-product loft are *banded* — a sample's basis row reaches only
`degree` control points in each direction, so `A[i][j]` vanishes beyond
`|i − j| > degree_x·n_ctrl_z + degree_z` — and factoring the band rather than
the full matrix turns the solve from `O(n³)` into `O(n·b²)`. A 160×30 net
lofts in ~0.14 s instead of ~14 s, so there is no longer a reason to run a
control net too coarse to hold the hull. That matters more than it sounds:
a least-squares fit that cannot reach its samples removes exactly the
short-scale content of `∂f/∂x` that feeds the **diverging** end of the
free-wave spectrum, and the first thing it biases is `R_w` at low Froude
number. On a real CAD import the previous 20×12 default overstated `R_w` by
**more than 3×** at Fn 0.15 and 2.5× at Fn 0.25, converging only around
80×24 — which is now the import default (with a 301×61 sample grid).
`FitReport::under_resolved` flags a fit whose RMS residual is still more
than 2% of the hull's own half-beam scale, and the CLI prints that as a
warning rather than letting it pass into a resistance curve.

A finer net does cost time downstream — every inner integral is linear in the
**span** count — so the kernel stops walking z-spans once `e^{−κ z₀}` falls
below `1e-20`. At large `λ` the decay `κ = νλ²` confines the integrand to a
sliver under the waterline, which is exactly where the outer quadrature spends
most of its evaluations, and the dropped terms are four orders below double
epsilon: on a real import this halves a resistance sweep with bit-identical
output. Net effect of resolving the geometry properly: a 31-speed sweep on an
8 m IGES hull goes from 9 s to 56 s, and stops being wrong by 3x. The residual is
a proxy for what actually matters (error in `∂f/∂x`, not in `f`), so treat it
as a floor on the problem, not a measure of it; where a grid carries observed
slopes, `FitReport::fx_residual` is the sharper signal.

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
michell view boat-*.hull --speed 8                       # interactive fleet viewer
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
waterplanes in gray), or emits CSV/JSON for other tooling. Conventions: the ship advances toward +x, so the wake trails toward
−x; the reconstruction is the far-field free-wave part of the linear
solution, physical astern of each hull (not on or ahead of it) — the PNG
fades the water where that caveat bites. Grids too coarse for the
shortest diverging waves are smoothly band-limited and flagged. The magnitude of A is pinned by
the deep-water free-wave resistance identity R_w = ½πρU² ∫|A|²cos³θ dθ; the
phase follows Tuck, Scullen & Lazauskas mapped to these conventions.

**Interactive viewer** (`michell view <hull>... --speed U`): serves a local
web page (default `http://127.0.0.1:8737`, `--port` to change) where the fleet
wave field is shown live and hulls can be **dragged** to reposition them
relative to one another. This exploits the superposition structure directly:
the fleet field is the exact sum of each hull's *standalone* field translated
to its placement (thin-ship amplitudes superpose, each carrying only a
placement phase), so each hull's field is computed once on a local grid, native
and parallel, and the browser recomposites the fleet by translate-and-sum as
you drag — no physics re-run, so dragging is instant. Re-running the solver is
reserved for the explicit controls (each an on-release action with a spinner):
**speed** rebuilds the per-hull fields (ν changes; a fraction of a second,
parallelised across hulls and row-bands), while **displacement** and
**ama immersion** re-float the assembly (full-band bodies; an IGES fleet is
decomposed to bodies once on load) and re-loft the wetted hulls. Same physicality caveat as `wake`: the field is faded ahead
of the aft-most stern. The server is dependency-free, in the spirit of the rest
of the crate. Build with `--release`; a debug build runs the integrals ~40×
slower and the viewer warns about it.

**Sweeps** (`michell sweep study.json`): long-form CSV/JSON over the
Cartesian product of axes — every varying quantity (speed, hull loads,
waterline, hull poses) is an axis, with fixed values as single-valued axes.
Each hull carries its own **load** — `mass` [kg] and a centre of gravity
(`lcg` longitudinal, `vcg` metres above the design floatplane) in the hull's
own frame — plus any number of discrete **point loads** (`points`: batteries,
crew, ballast), each a `mass` at an offset from the hull's centerpoint
(`dx` forward, `dy` to +y, `dz` **down**). The **fleet CG is never set
directly**: it is always the mass-weighted sum of every hull load and point
load, carried through each hull's pose, so mounting a hull moves its weight
with it (`dx`/`dz` translate each CG, design `trim` rotates it) and the
platform CG tracks the geometry automatically. An unspecified `lcg` defaults
to the hull's midship.

When the fleet carries mass it runs in **equilibrium mode**: each point's
platform sinkage (and pitch, from the derived LCG) is solved by a Newton
iteration whose Jacobian comes from the waterplane properties, so
counterfactuals like "what if this hull were heavier / its CG further forward"
are swept at physically consistent attitudes. Every such record carries the
solved state, displacement, LCB, the **derived** fleet CG (`mass`, `lcg`,
`vcg`, `tcg`), dry-hull count, and the resistance breakdown.
All poses are hydrostatic (no speed-dependent squat).

```json
{
  "name": "ama placement study",
  "fluid": "seawater",
  "hulls": [
    { "id": "vaka",  "file": "boat-center.hull",
      "load": { "mass": 1800, "lcg": -5.8, "vcg": 1.1 },
      "points": [
        { "id": "battery", "mass": 200, "dx": 1.0, "dz": 0.6 },
        { "id": "crew",    "mass": 160, "dx": -2.0, "dz": -0.4 }
      ] },
    { "id": "ama_s", "file": "boat-starboard.hull",
      "load": { "mass": 120, "vcg": 0.4 } },
    { "id": "ama_p", "file": "boat-port.hull", "pose": { "trim": 0.5 },
      "load": { "mass": 120, "vcg": 0.4 } }
  ],
  "sweep": [
    { "target": "speed", "unit": "knots", "range": [4, 10], "step": 0.5 },
    { "target": "vaka", "param": "mass", "range": [0, 800], "step": 200 },
    { "target": "battery", "param": "dz", "values": [0.0, 0.6, 1.2] },
    { "target": ["ama_s", "ama_p"], "param": "spread", "range": [1.5, 2.5] },
    { "target": "ama_s", "param": "trim", "values": [-2, 0, 2] },
    { "target": "vaka", "param": "scale", "values": [0.9, 1.0, 1.1] }
  ],
  "output": { "format": "csv", "file": "study.csv" },
  "options": { "rel_tol": 1e-5, "form_factor": 0.05, "transom": "ballistic=1.4" }
}
```

`options.transom` picks the transom-stern closure — `off`,
`ballistic[=COEFF]` (default, `COEFF = √2`), or `hollow=METRES` — matching the
`--transom` flag; it is inert on a hull whose half-breadth closes aft.

`options.dynamic: true` switches from hydrostatic to **dynamic** sinkage/trim
(see [Dynamic sinkage and trim](#dynamic-sinkage-and-trim) above): equilibrium
is re-solved at every speed rather than once per point, adding `fz` (dynamic
force, N) and `lift_pct` (as a fraction of the weight, %) columns; `sinkage`/
`trim_deg`/`volume`/`lcb` become the dynamic-attitude values. Requires a
`weight` axis. Real cost: the
near-field quadrature is far more expensive than the wave integral it reuses
parts of, and equilibrium calls it every Newton iteration at every speed —
budget minutes, not seconds, per sweep point, and prefer a handful of speed
values over a fine grid until you know how much resolution you need.
`options.squat_tol` overrides the closure's own quadrature tolerance
(default `1e-4`; loosening it trades sweep speed for a rougher force/moment).

Two situations hand the dynamic Newton solver a state that was converged
*somewhere else* rather than validated against what it is about to evaluate:
the coarse-to-fine handoff (a deliberately coarsened control net for early
iterations, full resolution for the polish — a coarsened loft cannot resolve
fine stern detail, a transom or a chine, the way the full-resolution one
does, so the same `(sinkage, trim)` can loft to a visibly different fleet),
and a sweep's speed-to-speed warm start (seeded from a *different* speed's
converged solution, whose dynamic force can be a different scale entirely).
Either can perturb the dynamic force enough that an undamped first Newton
step overshoots correcting for what is really a model or operating-point
shift, not a residual to chase — so that first step is damped whenever the
dynamic load is genuinely nonzero (inert, bit for bit, when it is zero or
absent); ordinary within-phase oscillation detection recovers full speed
within a couple more iterations regardless.

That damping alone isn't the whole story once a hull's transom is
substantial (a large `lift_pct`, order 15–20% of the weight, is the warning
sign): the near-field force can have genuine local curvature there — its own
hollow length depends on the current transom depth, which depends on
attitude — that the Newton core's analytic Jacobian (purely hydrostatic
waterplane properties) has no way to see, since it treats the dynamic load as
a *constant* added to the residual. Confirmed on the motivating case by
comparing loose- and tight-quadrature force evaluations across the operating
range: they agreed to <1%, ruling out quadrature noise, while the force
itself showed a real sign change in local slope — a genuine Jacobian
mismatch, not roughness. So on every phase but the initial, cold, from-scratch
one (i.e. exactly where the handoff damping above applies), the solver also
takes two extra finite-difference evaluations per iteration — perturbing
sinkage, then trim — and folds the dynamic load's own local sensitivity into
the Newton Jacobian alongside the hydrostatic terms. Deliberately **not**
applied during the initial coarse phase: that phase already converges
reliably on the hydrostatic Jacobian alone, and probing a finite difference
from a wild, far-from-solution starting guess turned out to be actively
harmful there (found by testing it unscoped: it sent a cold solve to a
multi-metre "sinkage" and a trim past the 20° abort limit). Real cost: this
roughly triples the per-iteration evaluation count on top of the near-field
quadrature's own expense, so a `dynamic: true` sweep on a transom-heavy hull
is priced in minutes per point, not seconds — but it is what took a case that
previously failed outright (a 3-speed sweep erroring at the last point,
residual 24× tolerance) to a clean converged solve at every speed. A hull
that closes cleanly aft (no `Transom` reported by `michell info`) never
exercises any of this — the Jacobian addition is `None`, bit for bit,
whenever the dynamic load has no local sensitivity to add.

Axis values: `range: [start, stop]` with optional `step` (default: a fifth
of the span), `values: [...]`, or scalar `value`. Speed axes take `unit`
(`ms` | `knots` | `froude`). Axis targets name hull ids or point-load ids
(all ids are unique) and **offset the base value** (a `scale` axis instead
*multiplies* the base) — for a **hull**, pose: `dx`, `dy`, `dz` (+down),
`spread` (outboard, sign follows each hull's side), `trim` (degrees, + raises
the +x end), `scale` (uniform size factor, `> 0`; grows or shrinks the hull in
place — length, beam, and draft all scale together — about its design waterline
and centre, so `1.0` leaves it unchanged and displacement goes as the cube),
and load: `mass`, `lcg`, `vcg`; for a **point load**, `mass`, `dx`, `dy`, `dz`
(relative to the hull centerpoint). A target list moves several targets as one
coupled axis (e.g. sweep both amas' `mass` together, or two symmetric ballast
points), but must be all hulls or all points. A `scale` in a hull's base `pose`
sets its built size. Hull files load relative to the manifest. A flag-based
sweep over raw IGES (`--axis`, `--float`) remains for one-liners.

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
`--sinkage`/`--platform-trim`/`--pivot-x` apply the whole-platform state an
equilibrium (floated-load) row reports. IGES inputs pass their surfaces
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
- `inner_integrals(hull, cond, λ)` — free-wave amplitude functions.
- `FreeWaveSpectrum::new(&members, &cond)` — far-field spectrum of a fleet:
  `amplitude(θ)`, `resistance_density(θ)` (dR_w/dθ), and Kelvin-wake
  reconstruction via `elevation_at(x, y)` / `elevation_grid(...)`.
- `hulls::wigley(l, b, t)` — exact reference hull.

## Roadmap

1. STEP reader feeding the same sample-and-loft pipeline; OBJ via the mesh
   path.
2. Python bindings.
3. Longitudinal wave cuts against published Wigley measurements; wake
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
- P. D. Kaklis & A. Papanikolaou et al., *Hydrodynamic optimization of
  fast-displacement catamarans*, 21st Symposium on Naval Hydrodynamics (1997),
  Appendix A — thin-ship theory for asymmetric demihulls via centreplane
  source + normal-dipole distributions (the asymmetric-hull extension here).
- ITTC — Recommended Procedures: *1957 ITTC Performance Prediction Method*.

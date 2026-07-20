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

Both front-ends reduce to the same deterministic back-end
(`fit::fit_offsets`): a grid of half-beam samples lofted to the spline by
separable tensor-product least squares (quantile knot placement, residuals
reported so you can judge fit quality).

- **Offsets**: `fit::fit_offsets(stations, waterlines, half_beams, opts)` —
  the human-authorable path: a station × waterline table of half-beams.
- **IGES**: `iges::import_hull(text, opts)` — reads **one or many untrimmed**
  NURBS patches (entity 128, unit weights; 124 transforms and unit conversion
  handled; naturally-bounded 143/141 wrappers, as produced by SubD → NURBS
  exports, are tolerated). Patches are mapped to the hull frame
  (`waterline_z` picks the DWL), the hull **centerplane is auto-detected** —
  a full both-sided shell is folded about the midplane of its shell
  intersections, a half hull measures from y = 0, and `centerplane` overrides
  either — then the wetted region is sampled by per-patch Newton inversion
  (outermost fold wins) and lofted. Genuinely trimmed surfaces (142/144) and
  rational weights are rejected with specific messages; fold asymmetry,
  ambiguous samples, failed inversions, and loft residuals (with the location
  of the worst one) are reported so a bad import is visible. Handles both the
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
```

**Sweeps** (`michell sweep`, IGES inputs): long-form CSV/JSON over the
Cartesian product of design and load axes. `--axis waterline=…` sweeps the
raw waterline; `--axis stem:dz|dx|dy|spread|trim=…` sweeps a hull's mount
pose (trim in degrees rotates the control nets exactly — affine); and
`--float weight=… [--float lcg=…]` switches to **equilibrium mode**: for
each load the platform's sinkage (and pitch) are solved by a Newton
iteration whose Jacobian comes from the waterplane properties, so
counterfactuals like "what if the boat were heavier / the CG further
forward" are directly sweepable at physically consistent attitudes. Every
record carries the solved state, displacement, LCB, and the resistance
breakdown. Hulls that fly dry at a pose contribute zero and are counted,
not errored. All poses are hydrostatic (no speed-dependent squat).

Hydrostatics on every hull: displaced volume, LCB, waterplane area and
moments, LCF — exact spline integrals.

Multihulls: list several hulls, each with an optional placement suffix —
`@y=Y` places a (single-hull file's) centerplane absolutely, `@dy=S` shifts
transversely, `@x=DX`/`@dx=DX` shifts longitudinally. An IGES file holding a
whole multihull imports as a fleet automatically: hulls are detected by
clustering wetted patches, each at its detected centerplane, and dry
structure (beams, decks) is dropped. The `IF` column / `interference` JSON
field reports combined R_w over the sum of standalone R_w.

All three input kinds are accepted everywhere (sniffed by header/extension):
the canonical `.hull` control net, a `michell-offsets v1` station × waterline
table (lofted on load), and IGES (sampled + lofted; `--waterline` sets the
DWL). `--json` gives machine-readable output; `--fluid`, `--rho`, `--nu`,
`--form-factor`, `--rel-tol` control the physics. Conversion diagnostics
(loft residuals, mirroring, failed inversions) are always printed so a bad
import can't pass silently.

**`.hull` format** (canonical, SI, `#` comments): `michell-hull v1`,
`degree-x/z`, `knots-x/z`, then one `row` of control values per x index.
**Offsets format**: `michell-offsets v1`, a `waterlines` line (z downward
from DWL, starting 0), then `station <x> <half-beams...>` lines.

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
- `iges::source_fleet` + `SourceFleet::situate(waterline, poses, platform)` —
  re-situate hulls repeatedly (immersion, mount trim, position) for sweeps.
- `float::solve_equilibrium(fleet, waterline, poses, load, ρ)` — hydrostatic
  sinkage/pitch balance for a mass + LCG load case.
- `inner_integrals(hull, cond, λ)` — free-wave amplitude functions.
- `hulls::wigley(l, b, t)` — exact reference hull.

## Roadmap

1. Parallel sweep evaluation (each equilibrium point is independent).
2. STEP reader feeding the same sample-and-loft pipeline.
3. Transom closure, wave spectrum output, Python bindings.

## References

- E. O. Tuck, *The wave resistance formula of J.H. Michell (1898) and its
  significance to recent research in ship hydrodynamics*, J. Austral. Math.
  Soc. B 30 (1989).
- J. Dambrine, M. Pierre, G. Rousseaux, *A theoretical and numerical
  determination of optimal ship forms based on Michell's wave resistance*,
  ESAIM: COCV (2016), arXiv:1410.2800.
- ITTC — Recommended Procedures: *1957 ITTC Performance Prediction Method*.

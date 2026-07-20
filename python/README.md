# pymichell

A small, **pedagogical** NumPy package that reproduces the core calculation
of the Rust [`michell`](../README.md) crate: from a hull described as a
B-spline half-breadth **control net**, compute the thin-ship **wave
amplitude grid** — the far-field Kelvin wave pattern — and the wave
resistance that pattern carries away.

It is deliberately *not* a port of the whole tool. But it reproduces the
mathematically interesting core faithfully: like the crate, the inner
amplitude integral `I + iJ` is evaluated **exactly, span by span** — on each
knot span `∂f/∂x` is a local polynomial and the two 1-D integrals collapse
onto closed-form moments — so `I` and `J` match the crate to machine
precision (~1e-16). Only the smooth *outer* integrals over propagation angle
(resistance, wake reconstruction) use plain trapezoidal quadrature, chosen
for legibility over the crate's adaptive oscillation-aware panels.

## The calculation

With `ν = g/U²` and the half-beam `f(x, z)` (`z` positive downward from the
waterline), everything flows from one integral over the hull surface:

1. **Amplitude function** — the Michell inner integral at `λ = sec θ`:

   ```
   F(λ) = I(λ) + i J(λ) = ∬ (∂f/∂x) · e^(−ν λ² z) · e^(i ν λ (x − x_c)) dx dz
   ```

   Evaluated **exactly**: on each knot span `∂f/∂x` is a local polynomial,
   extracted from corner Taylor data off the control net
   (`BSplineSurface.corner_partials`); the two 1-D integrals then reduce to
   the closed-form moments `∫ Xᵃ e^(ikX) dX` and `∫ Zᵇ e^(−κZ) dZ`
   (`moments.py`) and are summed span by span (`WaveField.inner_integral`).

2. **Free-wave spectrum** — the complex amplitude density per propagation
   angle:

   ```
   A(θ) = −(2ν/π) sec³θ · conj(F(sec θ))     [m/rad]
   ```

3. **Wave resistance** — computed two algebraically identical ways as a
   mutual check:

   ```
   R_w = (4 ρ g²)/(π U²) ∫₀^{π/2} |F(sec θ)|² sec³θ dθ      (Michell)
       = ½ π ρ U²        ∫_{−π/2}^{π/2} |A(θ)|² cos³θ dθ    (spectrum)
   ```

4. **The wave amplitude grid** — the Kelvin wake elevation, a superposition
   of one plane wave per angle:

   ```
   ζ(x, y) = Re ∫_{−π/2}^{π/2} A(θ) · e^(i ν sec θ ((x − x_c) + y tan θ)) dθ
   ```

## Install

```bash
pip install -e .            # numpy only
pip install -e ".[plot]"    # + matplotlib for the demo heatmap
```

## Use

```python
import numpy as np
from pymichell import wigley, load_hull, wave_field

surface = wigley(10.0, 1.0, 0.625)       # L, B, T [m]  — or load_hull("boat.hull")
field = wave_field(surface, speed=3.0)     # U [m/s], seawater

field.wave_resistance()                    # -> ~145.7 N
field.amplitude(np.radians(-20.0))         # -> A(θ), complex [m/rad]

x = np.linspace(surface.x_center - 50, surface.x_center + 10, 320)
y = np.linspace(-20, 20, 240)
zeta = field.elevation_grid(x, y)          # (len(y), len(x)) elevation grid [m]
```

Or run the demo:

```bash
python examples/demo.py                    # Wigley hull -> wake.png
python examples/demo.py ../ama.hull 8      # a real .hull control net at 8 m/s
```

## Layout

| file | contents |
|------|----------|
| `pymichell/bspline.py` | `.hull` parser, Cox–de Boor B-spline surface, knot spans, corner Taylor data |
| `pymichell/moments.py` | closed-form moment integrals `∫ tᵃ e^(ikt) dt`, `∫ tᵇ e^(−κt) dt` |
| `pymichell/wave.py`    | exact inner integral `F`, amplitude `A(θ)`, wake grid `ζ`, `R_w` |
| `tests/test_wigley.py` | validation against the crate's Wigley reference values |
| `tests/test_moments.py`| moments vs. direct quadrature (both branches) |
| `examples/demo.py`     | compute + plot a wake |

## Validation

`tests/test_wigley.py` checks against numbers produced by the Rust
`michell` binary for the Wigley hull at U = 3 m/s (seawater): the transverse
wavelength, the amplitude `A(θ)` at several angles, both resistance forms
agreeing with each other, and `R_w ≈ 145.76 N`. Run:

```bash
pytest            # or:  python tests/test_wigley.py
```

## Hull files: wetted surfaces vs. full-band bodies

Two flavours of `.hull` file exist, and both load:

* **Wetted-surface hulls** — `z = 0` is the waterline and the whole `z`
  domain is submerged (e.g. `michell wigley` output, or a lofted offsets
  table). Here `pymichell` reproduces the crate's core wave path *exactly*
  (Wigley `R_w` matches to ~0.01%).

* **Full-band bodies** — carry a `waterline` line: they model the hull up
  past the waterline, and the submerged part is `z ∈ [waterline, keel]`.
  `pymichell` runs Michell **directly on that control net**, restricting the
  depth integral below the waterline. The crate's CLI instead *re-samples
  and re-lofts* the wetted band into a new spline before integrating, so the
  two differ by ~1% on such files — a re-lofting difference, not a
  difference in the wave calculation. (That sample-and-loft trimming
  pipeline is the "full functionality" this package intentionally omits.)

## Scope

Single symmetric hull; deep water; no sinkage/trim; the free-wave
reconstruction is physical only *astern* of the hull. Multihull
interference, adaptive oscillation-aware outer quadrature, body re-lofting,
and the CAD/STL/IGES front-ends all live in the Rust crate. (The exact
closed-form *inner* integral, by contrast, is reproduced here in full.)

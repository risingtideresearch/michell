# Michell wave-resistance deep dive

Date: 2026-08-01
Branch: `story-wave-resistance-frontier`
Starting revision: `abe2b1f`
Remote operations: none

## Executive result

The project now has a published-value validation anchor, five physical/numerical
property checks, a repeatable benchmark harness, two reproduced numerical or
validation bugs with fixes, one measured quadrature optimization, and an exact
gradient of Michell wave resistance with respect to a symmetric B-spline control
net.

The gradient is the Phase-4 advance. It agrees with centered finite differences
for every Wigley control and is 14.42 times faster for the nine-control test case
(1.507 ms versus 21.734 ms median). The underlying quadratic dependence of
Michell resistance on a discretized hull is known in the literature; the new
result here is an engineering implementation that reverse-accumulates through
this code's exact per-span moment machinery. It is not claimed as novel.

## Claim classification

The required honesty classes are used throughout this report:

1. **Known result reproduced.** An external or analytic result reproduced here.
2. **Engineering improvement to this codebase.** A new test, fix, API, or
   measured optimization, without a priority claim.
3. **Matches published state of the art.** The implementation reaches a method
   or capability demonstrated in the cited literature.
4. **Possibly novel.** A result for which a documented search did not find prior
   art.

| Claim | Class | Evidence |
|---|---:|---|
| Standard Wigley `10^3 Cw` at `Fn=0.35` | 1 | 1.250831 computed versus 1.2486 published; 0.179% difference |
| Froude similarity, zero camber for symmetric sides, and classical catamaran `4 cos²` interference | 1 | Phase-0 property tests |
| Stable high-degree moment switchover | 2 | Failing regressions at degrees 19–20, followed by the endpoint-series fix |
| Per-side asymmetric half-breadth validation | 2 | Failing negative-port regression, followed by constructor validation |
| Adaptive outer-loop trigonometric reuse | 2 | Interleaved before/after benchmark at unchanged validation accuracy |
| Exact polynomial-times-kernel span integration | 3 | This code and Dambrine–Pierre–Rousseaux both evaluate polynomial basis integrals exactly |
| Exact B-spline control gradient | 1 and 2 | Known quadratic-form result reproduced; reverse implementation and public API are new here |

There are **no class-4 claims**. In particular, searches for “analytic gradient
Michell wave resistance”, “Michell integral shape derivative”, “B-spline hull
optimization Michell”, and “adjoint thin-ship wave resistance” found prior
quadratic-form and sensitivity work. The gradient is therefore not described as
a breakthrough.

## Phase 0 — ground truth and validation

### Documentation and initial state

`README.md`, `docs/michell-calculation.tex`, and
`docs/michell-calculation-asymmetric.tex` were read before numerical changes.
The checkout started clean and detached at `abe2b1f`; work moved to the local
branch named above.

### Baseline commands

| Command | Initial result |
|---|---|
| `cargo build --workspace` | pass, 15.66 s after fetching missing crates |
| `cargo build --workspace --release` | pass, 28.34 s |
| `cargo test --workspace` | all 176 pre-existing Rust tests pass |
| `uv run --with pytest --with numpy pytest -q` in `python/` | 11 passed, 0.60 s |

Plain `python3 -m pytest -q` could not run because `pytest` was not installed in
the ambient interpreter. That is an environment failure, not a test pass; the
isolated `uv` invocation above is the Python baseline.

### Published Wigley anchor

Doctors and Beck's Table 1 gives the classical thin-ship result
`10^3 Cw = 1.2486` for the standard Wigley hull with `B/L=0.1`,
`T/L=0.0625`, and `Fn=0.35` ([DOI record](https://doi.org/10.5957/jsr.1987.31.1.1)).
The implementation gives `10^3 Cw = 1.250831`, a relative difference of 0.179%.
The regression tolerance is 0.5% to avoid fitting differences in printed
precision, physical constants, or wetted-area convention.

This is a class-1 reproduction, not evidence that Michell theory is physically
accurate for arbitrary hulls.

### Harness coverage

The committed Phase-0 harness checks:

- the published Wigley value;
- cubic resistance scaling for geometrically similar hulls at fixed Froude
  number;
- identical port/starboard surfaces reducing exactly to classical Michell;
- an identical-demihull catamaran against an independent dense Simpson integral
  carrying the `4 cos²(ky s/2)` factor;
- monotonic self-convergence as requested tolerance tightens;
- reported versus actual error at `Fn = 0.08, 0.20, 0.35, 0.50`.

The first error-estimator regression was deliberately committed red. Diagnosis
showed that the alleged reference value had exceeded the four-million-evaluation
per-pass safety cap: its `max_lambda` retreated from about 32 to 11.8 as nominal
refinement increased. The production result was stable through the finest
uncapped pass. The reference was corrected to six forced refinements and now
asserts that it did not truncate earlier than the result under test. This was a
harness bug, not a numerical-core bug.

### Baseline benchmark

The dependency-free benchmark runs 30 warmed samples by default and supports
`MICHELL_BENCH_SAMPLES`. The initial measurements were:

| Case | Median | Best | Diagnostics |
|---|---:|---:|---|
| 21 speeds, `Fn=0.10…0.50` | 25.806 ms | 25.176 ms | checksum `2.553250998079e5` |
| Low Froude, `Fn=0.05` | 9.088 ms | 8.893 ms | 235,584 evaluations, `lambda_max=28.259`, estimated relative error `1.077e-9` |

Absolute wall-clock values vary with thermals and scheduling; performance claims
below use interleaved binaries from committed revisions rather than comparing
isolated runs.

## Phase 1 — numerical audit and bugs

Audited files: `michell.rs`, `moments.rs`, `quadrature.rs`, `bspline.rs`,
`inclined.rs`, `lifting.rs`, and `lifting3d.rs`, plus the asymmetric split in
`hull.rs`.

### Confirmed bug 1: high-degree moment recurrence

The fixed switch at `|kh|=4` handed all degrees to an upward recurrence. That
recurrence amplifies roundoff once degree becomes large relative to phase or
decay. Legal high-degree splines showed a discontinuity across the switch,
first exceeding the regression tolerance around oscillatory degree 19 and
exponential degree 20.

The failing tests were committed in `624f256` before the fix. The fix in
`5772f43` uses the complementary endpoint expansions

```text
integral(u^a exp(i x u), 0..1)
  = exp(i x) sum_m (-i x)^m a!/(a+m+1)!

integral(u^b exp(-x u), 0..1)
  = exp(-x) sum_m x^m b!/(b+m+1)!
```

when degree overtakes kernel magnitude. Their terms contract precisely where
the upward recurrence becomes unstable. The original origin series and fast
recurrence remain in their stable regimes. Real, oscillatory, and complex-decay
switch tests pass.

Classification: class 2.

### Confirmed bug 2: asymmetric side validation

`Hull::new_asymmetric` validated only the symmetric mean surface. A negative
physical port half-breadth could therefore be accepted when a positive
starboard surface masked it in the mean, contradicting the constructor's
geometry contract.

The failing regression was committed in `31f44aa`; `daa3f92` validates each
physical control net independently before decomposition.

Classification: class 2.

### Reproduced edge cases without defects

- Interior knots with multiplicity equal to degree: accepted as a `C0` chine;
  zero-length repeated spans are filtered; corner partials and one-sided values
  remain finite.
- Zero-length knot spans: excluded from all exact-moment accumulation.
- Odd/even Gauss–Legendre rules: polynomial exactness tests cover orders through
  24.
- Symmetric thickness/camber split: identical sides produce bit-for-bit zero
  camber contribution.
- Lifting solvers: existing analytic limits, Kutta behavior, symmetry, and grid
  convergence remained green.
- Inclined hydrostatics: existing box, equilibrium, heel-symmetry, and grid
  tests remained green.

### Suspected or unverified; no speculative fix

- The outer tail stop uses one phase-width window and a hard evaluation cap.
  The Phase-0 cases and low-Froude benchmark did not reproduce an incorrect
  supported-tolerance result, but this is not a rigorous tail bound for every
  possible multi-span or multihull interference envelope. Backlog: require
  multiple decreasing windows and expose an explicit “safety cap reached” flag.
- `solve_dense` uses a debug-only singularity assertion. Invalid or degenerate
  user-supplied 3-D panels could yield non-finite release results. No failure was
  reproduced for the library's builders; input conditioning belongs in a
  separate hardening change.
- The asymmetric lifting/dipole closure is documented as approximate. Its
  physical normalization was not independently validated against experimental
  asymmetric-hull data, so no stronger claim is made.

## Phase 2 — profile and performance

A 10-second macOS sampled profile of a 1,000-sample release benchmark placed the
dominant cost in `integrate_outer`, `InnerIntegral::accumulate`, moment
evaluation, and their trigonometric/exponential kernels. Gauss-node generation
and allocation were negligible.

The retained optimization computes the fixed GL16 rule once per resistance
solve, reuses `sin_cos` at panel boundaries, and avoids repeated endpoint
cosines. Interleaved 100-sample runs against the immediately preceding committed
binary gave:

| Case | Before median | After median | Change |
|---|---:|---:|---:|
| 21-speed Wigley sweep | 25.681 ms | 25.296 ms | 1.5% faster |
| `Fn=0.05` Wigley | 9.024 ms | 8.876 ms | 1.6% faster |

The resistance checksums were identical to 12 printed digits and the Phase-0
harness passed. A second proposed zero-placement phase shortcut was discarded:
interleaved results were 24.133 versus 24.057 ms for the sweep and 8.469 versus
8.480 ms at low Froude, i.e. noise rather than a defensible improvement.

Classification: class 2.

## Phase 3 — source-verified frontier map

### Tuck, Michell numerics, and Michlet

Tuck's 1987 report uses a special exponential-aware trapezoidal rule vertically,
Filon quadrature longitudinally, and fixed-interval Simpson quadrature in wave
angle, typically 40 intervals. It explicitly warns that the longitudinal
integral becomes highly oscillatory toward the diverging-wave endpoint
([actual report PDF](https://cembercikutuphanesi.com.tr/catamaran/design/t8701.pdf)).
Tuck's later historical/numerical assessment documents Michell's formula and the
then-current computational context
([actual paper PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6D0B69CE2AE6BDC1D06BA675F1C4DEDD/S0334270000006329a.pdf/wave_resistance_formula_of_jh_michell_1898_and_its_significance_to_recent_research_in_ship_hydrodynamics.pdf)).

Lazauskas's Michlet work extends practical thin-ship calculation to multihulls,
wave patterns, and viscous/boundary-layer corrections; his thesis emphasizes
viscous smoothing of Michell humps and hollows
([Adelaide thesis](https://digital.library.adelaide.edu.au/bitstreams/d054ebb4-f5a1-41ae-8bbf-678421fcfa80/download)).
The inspected Michlet 8 manual documents a user-selected even angular count of
10–4096 and recommends at least 160 for monohulls, increasing with hull count
([manual copy](https://www.scribd.com/document/47263938/MICHLET-users-manual-thin-ship-theory)).

Compared with those implementations, this code's exact polynomial B-spline
inner integrals eliminate longitudinal and vertical sampling error, and its
outer rule is adaptive with an error estimate. That is an engineering advantage
over the documented fixed grids, not a claim to dominate every modern Michell
implementation.

### Exact discretizations and hull optimization

Dambrine, Pierre, and Rousseaux write the discretized Michell resistance as the
quadratic form `F^T M_w F`, compute polynomial basis integrals exactly, and solve
a regularized hull-optimization problem. They also show that minimizing wave
resistance alone is ill-posed and requires regularization/constraints
([actual paper PDF](https://www.numdam.org/item/10.1051/cocv/2014067.pdf)).

This code matches the exact-polynomial-inner-integral aspect for arbitrary
polynomial tensor B-splines rather than their bilinear tent functions (class 3).
The Phase-4 gradient is the matrix-free differential of the same known
quadratic structure (class 1 reproduced, class 2 implementation).

Recent work continues Michell-based optimization rather than rendering it
obsolete: robust optimization under uncertain cruise speed appeared in 2024
([DOI](https://doi.org/10.1002/mma.9693)), and a 2024 geometric-operator method
targets inexpensive sensitivity/dimensionality reduction using slender-body
wave resistance ([preprint](https://arxiv.org/abs/2403.06990)).

### Neumann–Michell theory

Noblesse, Huang, and Yang's Neumann–Michell theory revises the linear
free-surface boundary-integral formulation and removes the classical
Neumann–Kelvin waterline integral by a consistent linearization and integration
by parts ([DOI](https://doi.org/10.1007/s10665-012-9568-7)). It is a different,
higher-fidelity boundary-integral model, not a quadrature tweak to Michell's
centerplane source formula.

Adopting it would require a hull-surface potential solve and new Green-function
machinery. It could improve non-slender applicability, but it would sacrifice
much of the present solver's closed-form speed and simplicity. It belongs as a
separate solver sharing geometry and validation infrastructure.

### Oscillatory quadrature

Relevant mature families are:

- Filon/Filon–Clenshaw–Curtis, which interpolates the slowly varying amplitude
  and integrates oscillations analytically; stable rules can be `O(N log N)`
  and improve as frequency grows
  ([Domínguez, Graham, Smyshlyaev](https://doi.org/10.1093/imanum/drq036)).
- Levin collocation, which converts the oscillatory integral into an auxiliary
  differential problem and can have decreasing relative error with frequency
  ([actual author-hosted paper](https://www.math.tau.ac.il/~levin/OscIntAnal.pdf)).
- Numerical steepest descent/analytic continuation, whose accuracy can improve
  with frequency for analytic integrands
  ([Huybrechs and Vandewalle](https://doi.org/10.1137/050636814)).
- Double-exponential Fourier formulas for infinite oscillatory tails
  ([Ooura](https://ems.press/journals/prims/articles/2320)).

These methods are directly attractive for the complex inner amplitude. This
code already removes that difficulty more strongly by integrating each
piecewise-polynomial inner span in closed form. Applying the same methods to the
outer resistance is less direct: after `lambda=sec(theta)`, the integrand is a
nonnegative squared modulus containing many cross-frequency interference terms,
not one analytic amplitude times `exp(i omega g)`. Complex steepest descent is
also obstructed by complex conjugation unless the square is first expanded into
analytic cross terms. A Filon/Levin outer method remains promising, especially
at low Froude number, but needs a careful phase decomposition and independent
tail proof.

### Finite depth and restricted water

Sretensky-type theory replaces the unrestricted angular continuum by channel
modes; Dambrine et al. note that their functional framework also covers a
truncated Sretensky summation. Finite width and finite depth additionally change
the dispersion relation and introduce subcritical/critical/supercritical
regimes. A primary symposium treatment gives the channel formula and implicit
finite-depth dispersion relation
([National Academies chapter](https://nap.nationalacademies.org/read/9771/chapter/5)).

This code is deep-water and laterally unbounded. A finite-depth kernel is a
high-value known extension, but validation must cover the critical-depth
singularity and the disappearance/change of transverse modes; simply replacing
the exponential decay is insufficient.

### Physical corrections newer than classical Michell

Boundary-layer displacement/tangency corrections can improve resistance
predictions for less-slender hulls. Bašić, Blagojević, and Andrun report improved
results over original Michell for five hull families by including boundary-layer
effects and a tangency correction
([article](https://doi.org/10.1016/j.oceaneng.2020.107079)). Tsubogo proposes a
depth-gradient modification that shifts and damps humps and hollows
([DOI](https://doi.org/10.2534/jjasnaoe.19.19)). These address model-form error,
whereas the work in Phases 0–2 addresses numerical error.

### What this code already does well

- Exact inner integrals for every span of an arbitrary-degree polynomial
  tensor B-spline; no station/waterline sampling error in Michell amplitude.
- Stable origin, recurrence, and endpoint regimes for real, oscillatory, and
  complex-decay moments.
- Adaptive endpoint-regularized outer integration with diagnostics.
- Coherent multihull phase superposition, heel via complex vertical decay, and
  experimental asymmetric/dipole paths.
- A matrix-free exact control gradient whose cost does not grow by one primal
  solve per design variable.

The first, and now the gradient capability, match known published practice in
their respective areas. The combined arbitrary-degree B-spline implementation
is useful, but no exhaustive search established priority for that exact
combination, so no novelty claim is made.

## Phase 4 — exact control-net gradient

### Candidate ranking

| Rank | Candidate | Leverage | Risk | Decision |
|---:|---|---|---|---|
| 1 | Exact `dR_w/dP` for B-spline controls | Enables deterministic hull optimization; natural quadratic structure | Low–medium | Implemented |
| 2 | Phase-decomposed Filon/Levin outer integral | Could make very low Froude cheap and tail-safe | High; squared-modulus phase decomposition and tail proof | Backlog |
| 3 | Finite-depth/channel kernel | Large practical scope increase | High; critical modes and new validation corpus | Backlog |
| 4 | Boundary-layer tangency correction | Improves physical prediction beyond slender hulls | Medium–high; needs empirical/CFD calibration | Backlog |
| 5 | Neumann–Michell surface solver | Higher-fidelity linear theory | Very high; a separate BIE solver | Backlog |

### Method

For fixed knots, degrees, speed, and fluid conditions, every inner amplitude is
linear in the control net. On the final accepted outer-quadrature pass the code
accumulates

```text
d|F|²/dc_k = 2 Re(conj(F) dF/dc_k)
```

for the local span-polynomial coefficients `c_k`, reusing exact x/z moments.
It then applies the transpose of the exact control-to-corner-derivative map to
obtain `dR_w/dP_ij`. The public API returns the primal `WaveResistance`, the
row-major control gradient, and the number of reverse-pass inner evaluations.

The initial API intentionally accepts only symmetric hulls. An asymmetric hull
has two physical nets and multiple dipole closures; returning a derivative of
only its stored symmetric mean would be misleading.

### Validation and claimed advantage

| Check | Result |
|---|---|
| Every one of 9 controls versus centered finite differences | pass; scaled error below `2e-6` |
| Constant half-breadth shift (null direction of `df/dx`) | pass |
| Returned primal versus ordinary resistance API | bit-for-bit equal in test |
| Analytic work versus control count | one primal convergence plus one reverse pass |
| Full Phase-0 harness | pass |
| Full Rust workspace | 190 tests pass |

Release benchmark, 30 samples, default tolerance:

| Gradient method | Median | Best |
|---|---:|---:|
| Exact reverse, 9 controls | 1.507 ms | 1.483 ms |
| Centered finite differences | 21.734 ms | 21.617 ms |

Median speedup: **14.42×**. Aggregate absolute-gradient checksums agree to about
`2.4e-11` relative (`1.734433365029e4` versus `1.734433364988e4`, including
the benchmark's 30-run accumulation).

Classification: the quadratic derivative is class 1; this matrix-free B-spline
implementation and API are class 2. No class-4 claim.

## Commit map

| Commit | Logical change |
|---|---|
| `05fd665` | Phase-0 validation and benchmark harness (with red capped-reference test) |
| `c0e1299` | Correct and guard the capped validation reference |
| `624f256` | Failing high-degree moment switchover regressions |
| `5772f43` | Stable endpoint expansions for high-degree moments |
| `31f44aa` | Failing asymmetric negative-half-breadth regression |
| `daa3f92` | Validate both physical asymmetric sides |
| `7ec095b` | Chine/repeated-span audit coverage |
| `6eb45ad` | Profiled outer-quadrature optimization |
| `95569fb` | Exact control gradient, validation, and benchmark |

## Ranked backlog

1. Add explicit convergence status (`converged`, `tail_cap`, `evaluation_cap`)
   and make invalid reference runs impossible to mistake for converged values.
2. Replace the single quiet-window tail heuristic with a tested envelope or
   multi-window remainder bound; build adversarial multihull interference cases.
3. Develop an outer-integral phase decomposition, then compare Filon–Clenshaw–
   Curtis, Levin, and numerical steepest descent at `Fn <= 0.05` against a
   high-precision independent reference.
4. Extend exact gradients to port and starboard nets, multihull placements,
   speed, and constrained objectives (volume, wetted area, fairness).
5. Add a small constrained optimizer example with non-negativity, closure,
   displacement, and curvature regularization; do not optimize wave resistance
   alone because the published problem is ill-posed.
6. Implement finite-depth/infinite-width first, then finite-width channel modes;
   validate both sides of critical depth Froude number.
7. Add boundary-layer displacement/tangency corrections behind an explicit
   model option and validate against the five-hull literature set.
8. Harden lifting solvers with public dimension/geometry validation and
   condition estimates instead of debug-only singularity assertions.
9. Treat Neumann–Michell as a separate higher-fidelity solver sharing the same
   B-spline geometry, benchmarks, and published Wigley cases.

## Final validation protocol

The release handoff should run and record:

```text
cargo build --workspace
cargo build --workspace --release
cargo test --workspace
cargo bench -p michell --bench wigley
cd python && uv run --with pytest --with numpy pytest -q
```

No database-backed tests exist in this repository. No remote push was made.

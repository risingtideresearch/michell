# Michell wave-resistance deep dive

Date: 2026-08-02
Branch: `story-kernel-hardening`
Kernel-hardening starting revision: `30d8f53`
Design-tool upgrade starting revision: `64925d4`
Original deep-dive starting revision: `abe2b1f`
Remote operations: none

## Executive result

The project now has a published-value validation anchor, physical and numerical
property checks, a repeatable benchmark harness, reproduced numerical bugs with
red-before-fix commits, a measured general-regime quadrature optimization, an
exact design-variable adjoint, attributable wave signatures, a new low-Froude
solver, and a six-design ranking-stability gate.

The design-tool upgrade from `64925d4` completed all four requested workstreams.
Results now identify both their numerical method and termination outcome;
multihull placement, symmetric/asymmetric control-net, displacement, LCB, and
wetted-area gradients are available; angular waves and low-Froude endpoint
pairs carry signed attribution; and design rankings are checked across three
tolerances, two numerical routes, and exact knot insertion. The new ranking
gate exposed one additional diagnostic bug: at `Fn=0.05` the marcher reported
`5.123e-9` while its actual error was `5.781e-7`. A power-law tail extrapolation
fixed that low-Froude failure. Independent follow-up then found that its first
1.25 safety factor was still 2.97 times optimistic at `Fn=0.35`. A red
λ-to-4000 analytic-Wigley regression anchored the design-Froude regime, but a
second independent extension found undercoverage again above `Fn≈0.5`.
Assertion-red commit `cc25d18` reproduces it. The decisive fix separates true
x/y oscillation phase from vertical exponential decay when forming the quiet
window and adds, rather than maximizes, refinement and tail diagnostics. The
12-point `Fn=0.02…1.00` sweep now has 3.87×–4.68× coverage.

The second Phase-4 advance is the larger one. It rewrites each polynomial
B-spline span exactly as endpoint waves, discards only depth-damped endpoints
under an explicit absolute bound, expands the squared amplitude into pairwise
kernels, and evaluates those kernels on a Gaussian-decaying steepest-descent
contour. At `Fn=0.02` the result needs 1,152 kernel evaluations instead of
511,504 inner-amplitude evaluations and is **713.48 times faster** (0.029 ms
versus 20.691 ms median in the final run). Against an independent real-axis
reference its relative difference is `3.56e-11`; the reference's finite tail,
rather than the new solver, limits that comparison. Cost is effectively
independent of the oscillation frequency in the tested low-Froude range.

The constituent mathematics has a substantial analytic lineage: Birkhoff and
Kotik's kernel separation, Michelsen's 1960 polynomial reduction and 1972 JSR
sequel, the Sendagorta--Grases tabulation program, low-speed endpoint
asymptotics, Bickley--Naylor functions, and numerical steepest descent. The
validated B-spline implementation, cancellation accounting, omission bound,
contour evaluator, error-gated dispatch, and real-axis fallback are class-2
engineering improvements. This report makes no priority claim for their
combination.

The first Phase-4 advance, the exact gradient, remains useful: it is 17.59 times
faster than centered finite differences for the nine-control benchmark. Its
quadratic structure is known; the matrix-free B-spline reverse pass is an
engineering implementation, not a novelty claim.

## Kernel hardening after coordinated adversarial review

Classification: the failure reproduction and degree-elevation identities are
class 1; the refusal gate, finite-value checks, compensated endpoint
accumulation, and revised diagnostics are class 2. No novelty claim is made for
this work.

Three independent reviews showed that the former “arbitrary-degree” contract
was false. `Hull::fx_coeff` converted corner derivatives to local Taylor
coefficients with factorial scaling. At high degree that map is ill-conditioned,
and both the endpoint reducer and general marcher consumed the same corrupted
array. Their agreement was therefore correlated, not independent validation.
The original reviewer programs and hashes are preserved in
`probes/adversarial/`; commit `b2cb256` ports all required cases as red tests.

| Red case on `30d8f53` | Accepted result | Independent evidence |
|---|---:|---:|
| R1 degree `(48,2)`, `Fn=0.05` | `2.240158007412e11 N`, reported `1.271e-11` | closed-form amplitude plus independent GL: `9.840371404040e-3 N` |
| R1 degree sweep | first requested-`1e-5` violation at degree 24: `1.459e-5` actual, `2.917e-12` reported | analytic Bernstein scaling |
| R2 degree `(48,2)`, `Fn=0.10` | `5.624385122311e13 N` | rigorous variation/envelope upper bound `4.9899908e3 N` |
| R2 degree `(192,2)` | `Ok(NaN)`, `EvalCap` | finiteness contract |
| R3 degree `(24,16)`, `Fn=0.05` | `1.305597e-5` actual, `7.417e-14` reported, `Converged` | independent high-precision reference `6.3825284969346485e-3 N` |

The fix deliberately chooses the explicit validated-envelope option. Replacing
only the derivative extraction would leave high-order cancellation in the
downstream local power basis and endpoint expansion. `Hull` now supports degrees
up to 16 in each direction, subject to the reconstruction gate, and returns
`Error::Unsupported` above that cap.
Within the cap, construction independently reconstructs `fx` from the local
coefficients at a tensor Gauss grid and refuses a normalized residual above
`1e-8`. Every public resistance route rejects non-finite outputs rather than
returning `Ok(NaN)`.

Error reporting now charges both solvers a `2e-8` resistance-level coefficient
floor. That empirical floor is about 18 times the independent R3 degree
`(16,12)` discrepancy (`1.113940e-9`); it is not a formal floating-point proof.
The endpoint combined map uses compensated complex accumulation and propagates
a standard gamma-style roundoff bound through every endpoint pair. Public
rustdoc states explicitly that marcher/reducer agreement is not independent
evidence because both share `fx_coeff`.

Exact Bézier degree elevation of one Wigley geometry validates degrees 2–16 in
x alone, z alone, and along the tensor diagonal. Degree 1 remains covered by
the existing analytic linear/wedge tests; a quadratic Wigley cannot be
represented at degree 1. At the worst supported tensor corner `(16,16)`, the
relative discrepancy was `2.514206e-10` for the endpoint route at `Fn=0.05`
and `4.896228e-10` for the general marcher at `Fn=0.35`; reported estimates
were `2.000378e-8` and `6.100313e-8`. Degree 17 in either direction refuses.
The unchanged exact knot-insertion ranking gate also passes.

The low-degree resistance values did not change. The 30-sample release sweep
checksum remains exactly `2.553251156101e5`. Current medians were 12.726 ms for
the 21-speed sweep, 0.029 ms for the `Fn=0.05` endpoint case, 20.852 ms for the
direct `Fn=0.02` marcher, and 0.030 ms for its endpoint reduction. The new
diagnostic floor changes printed error estimates, as intended, but not those
resistance checksums.

## Claim classification

The active honesty classes used throughout this report are:

1. **Known result reproduced.** An external or analytic result reproduced here.
2. **Engineering improvement to this codebase.** A new test, fix, API, or
   measured optimization, without a priority claim.
3. **Matches published state of the art.** The implementation reaches a method
   or capability demonstrated in the cited literature.

The former priority category has been withdrawn after equation-level review of
the Michelsen and Birkhoff--Kotik lineage. No current result is assigned a
priority-seeking classification.

| Claim | Class | Evidence |
|---|---:|---|
| Standard Wigley `10^3 Cw` at `Fn=0.35` | 1 | 1.250831 computed versus 1.2486 published; 0.179% difference |
| Froude similarity, zero camber for symmetric sides, and classical catamaran `4 cos²` interference | 1 | Phase-0 property tests |
| Stable high-degree moment switchover | 2 | Failing regressions at degrees 19–20, followed by the endpoint-series fix |
| Validated spline-degree envelope and loud refusal above degree 16 | 2 | Five adversarial red cases, independent high-precision references, exact degree-elevation sweep, and finite-result guards |
| Per-side asymmetric half-breadth validation | 2 | Failing negative-port regression, followed by constructor validation |
| Adaptive outer-loop trigonometric reuse | 2 | Interleaved before/after benchmark at unchanged validation accuracy |
| Exact polynomial-times-kernel span integration | 3 | This code and Dambrine–Pierre–Rousseaux both evaluate polynomial basis integrals exactly |
| Exact B-spline control gradient | 1 and 2 | Known quadratic-form result reproduced; reverse implementation and public API are new here |
| Explicit method/outcome status and `.msw` v2 diagnostics | 2 | Cap regressions, CLI/archive compatibility tests, and v1 reader coverage |
| Power-law marcher-tail diagnostic | 1 and 2 | Known algebraic endpoint decay; independent λ-to-4000 Wigley regressions cover reproduced low- and design-Froude failures |
| Multihull placement and asymmetric per-side gradients | 1 and 2 | Analytic phase/chain-rule derivatives; every requested component is centered-FD checked |
| Displacement, LCB, and wetted-area control derivatives | 1 and 2 | Standard spline/calculus derivatives; exact polynomial and differentiated-quadrature implementation |
| Member/pair wave and endpoint attribution | 1 and 2 | Standard quadratic interference expansion; new typed attribution API and summation regressions |
| Six-design ranking-stability gate | 2 | Ordering checks across tolerances, solver routes, and exact knot insertion |
| Low-speed dominance by waterline bow/stern data | 1 | Keller–Ahluwalia, Wehausen/Kotik, and Gotman endpoint results reproduced computationally |
| Steepest-descent evaluation of oscillatory ship-wave integrals | 1 | Motygin and the general numerical-steepest-descent literature |
| Bounded low-Froude solver and automatic fallback | 2 | New implementation, independent reference tests, and routing regression |
| Historical analytic reduction lineage: Birkhoff--Kotik, Michelsen, and Sendagorta--Grases | 1 and 3 | Known kernel/basis separation and polynomial or orthogonal-basis reductions, reproduced or matched in implementation scope |
| B-spline endpoint/Bickley/steepest-descent solver, error accounting, and fallback | 2 | Validated implementation within the degree-16 and endpoint-separation envelope; no method-priority claim |

The gradient searches (“analytic gradient Michell wave resistance”, “Michell
integral shape derivative”, “B-spline hull optimization Michell”, and “adjoint
thin-ship wave resistance”) found prior quadratic-form and sensitivity work.
The low-Froude review found the older analytic-reduction lineage summarized in
the manuscript's equation-level comparison table. Those findings support the
class-1/class-3 historical classification and the class-2 implementation
classification above.

## Design-tool-grade upgrade from `64925d4`

### Workstream 1 — honest convergence status: complete

`WaveResistance` now carries explicit `method: WaveMethod` and
`outcome: WaveOutcome`. Methods are `GeneralMarcher` and `EndpointReduction`;
outcomes are `Converged`, `TailCap`, `EvalCap`, and the additional honest state
`RefinementCap` for a completed tail march whose requested panel tolerance was
not met. The old `max_lambda == infinity` route sentinel remains a useful
physical diagnostic but is no longer an API discriminator.

The status propagates through single- and multihull resistance, upright and
heeled paths, spectrum/wake diagnostics, CLI JSON, and sweep archives. The
`.msw` container is explicitly version 2 for the added scalar fields; row
framing is unchanged and the reader retains version-1 compatibility. New CLI
and archive fields include method, outcome, estimated relative error,
`max_lambda`, and evaluation count.

The cap-status contract tests in `6b8c42d` are deliberately compile-red: they
name `WaveOutcome`, which did not yet exist. This establishes the missing
observable API at the type level, but is not an assertion-red demonstration of
the old runtime behavior. Commit `ca2d648` added the statuses and propagation;
the same scenarios then assert the distinct `TailCap` and `EvalCap` outcomes.

The ranking work later exposed a subtler diagnostic failure. At `Fn=0.05`:

| General marcher diagnostic | Refinement only | First tail estimator | Corrected estimator | Independent actual error |
|---|---:|---:|---:|---:|
| Relative error | `5.123e-9` | about `6.356e-7` | about `2.034e-6` | `5.781e-7` |

The old value only measured panel refinement. The tail diagnostic adds a
safety-factor extrapolation of the terminating phase window. The extrapolation
uses the endpoint result `F = O(lambda^-3)`, hence a transformed resistance
density `O(lambda^-5)` and tail proportional to one quarter of the local
density times `lambda`.

Independent review found that the initial 1.25 factor covered this low-Froude
case but not the finite-λ transition at design Froude numbers. Red commit
`30bc6ef` reproduces the worst reported case without production moments or
outer quadrature: analytic Wigley amplitude, phase-resolved GL16 panels,
compensated summation, and `lambda_max=4000`. Commit `8a7ac15` raises the factor
to 4.0 and extends the reference sweep:

| Fn | Initial estimate | Actual relative error | Corrected estimate | Corrected/actual |
|---:|---:|---:|---:|---:|
| 0.12 | `1.747e-7` | `2.128e-7` | `5.589e-7` | 2.63× |
| 0.20 | `9.040e-8` | `1.452e-7` | `2.893e-7` | 1.99× |
| 0.35 | `6.129e-8` | `1.822e-7` | `1.961e-7` | 1.08× |

That factor-only correction moved rather than eliminated the boundary.
Independent extension found only 0.642× coverage at `Fn=0.70` and 0.346× at
`Fn=1.00`; assertion-red commit `cc25d18` adds both cases. The proposed
`max(refinement, tail) → refinement + tail` change was necessary but not
sufficient: at `Fn=0.70` the two terms summed to about `1.36e-7`, still below
`1.82e-7` actual error.

The structural cause was subtler. Panel sizing legitimately included the
vertical exponential-decay rate, but the same combined rate advanced the
supposed oscillation window. At high Froude number, already-vanished submerged
terms could therefore advance the window through `8π` while its actual
longitudinal phase covered only a narrow `cos²` trough. Commit `a2aa894` uses
only physical longitudinal/transverse phase for the stopping window and sums
the two distinct error diagnostics. The λ-to-4000 compensated reference now
covers the experimental program's full range:

| Fn | Reported estimate | Actual relative error | Reported/actual |
|---:|---:|---:|---:|
| 0.02 | `2.497e-6` | `6.334e-7` | 3.94× |
| 0.03 | `1.300e-6` | `3.322e-7` | 3.91× |
| 0.05 | `5.774e-7` | `1.493e-7` | 3.87× |
| 0.08 | `2.696e-7` | `6.918e-8` | 3.90× |
| 0.12 | `1.437e-7` | `3.540e-8` | 4.06× |
| 0.20 | `6.897e-8` | `1.726e-8` | 4.00× |
| 0.35 | `4.100e-8` | `9.626e-9` | 4.26× |
| 0.40 | `2.866e-8` | `6.970e-9` | 4.11× |
| 0.45 | `2.599e-8` | `5.552e-9` | 4.68× |
| 0.50 | `2.286e-8` | `5.385e-9` | 4.24× |
| 0.70 | `1.740e-8` | `4.114e-9` | 4.23× |
| 1.00 | `1.372e-8` | `3.124e-9` | 4.39× |

The estimator remains a heuristic, not a rigorous bound or a claim about all
hulls. When its fixed tail floor exceeds the requested tolerance, panel
refinement stops once its own change is smaller and the result reports
`RefinementCap` rather than burning the entire allowance or claiming success.
The public `WaveOptions`, `WaveOutcome`, and `WaveResistance` documentation now
states prominently that `rel_tol` is a target and that callers must inspect
both `outcome` and `est_rel_error`.

Classification: endpoint decay is class 1; status plumbing, the estimator,
archive versioning, and regressions are class 2. No novelty claim.

### Workstream 2 — design-variable gradients: complete

The exact reverse pass now handles a fleet in one final outer pass. For each
member it differentiates both placement phases,

```text
exp(i nu (lambda dx_j +/- lambda sqrt(lambda^2-1) y_j)),
```

and the member's local source/camber amplitude. The asymmetric chain rule maps
the mean/camber adjoints back to independent physical controls:

```text
d/d port      = 1/2 (d/d source - d/d camber)
d/d starboard = 1/2 (d/d source + d/d camber).
```

This exactly covers the default asymmetric strip closure. It deliberately does
not claim to differentiate the optional solved-lifting closure. Both primal and
reverse paths always use `GeneralMarcher`; a regression proves that at
`Fn=0.05` the ordinary primal selects `EndpointReduction` while the gradient's
embedded primal reports `GeneralMarcher`.

Constraint derivatives are exposed separately on `Hull`. Displaced volume and
its longitudinal first moment integrate B-spline bases exactly to floating-
point roundoff; LCB uses the quotient rule. Wetted area differentiates through
the identical 24-point-per-span Gauss–Legendre rule used for the reported area.
Symmetric nets include both physical sides; asymmetric results return separate
port/starboard arrays. A zero-volume hull returns an explicit error because its
LCB derivative is undefined.

Validation results:

| Derivative family | Coverage | Maximum allowed scaled FD error |
|---|---|---:|
| Symmetric wave controls | every control | `2e-6` |
| Interfering multihull controls | selected controls on every member | `1e-4` |
| Multihull `x` and `y` placement | every member/component | `1e-4` |
| Asymmetric wave controls | every port and starboard control | `3e-6` |
| Displaced volume | every symmetric and per-side asymmetric control | `2e-8` |
| LCB and wetted area | every symmetric and per-side asymmetric control | `2e-7` |

The final 30-sample release benchmark retains the constant-cost advantage:
exact reverse `0.820 ms` median versus `14.431 ms` for centered finite
differences, a `17.59x` speedup for nine controls. The aggregate gradient
checksums agree to about `8.1e-12` relative.

Classification: analytic phase derivatives, spline-basis integrals, quotient
rules, and reverse differentiation of a quadratic form are class 1. The fleet,
asymmetric, constraint APIs and their implementation are class 2. No class-4
claim.

### Workstream 3 — attributable wave signatures: complete

`FreeWaveSpectrum::signature(theta)` returns complex amplitudes in input-member
order, their total, the total amplitude squared, and an upper-triangular ledger
of signed self/pair contributions. Diagonal terms are `|A_j|^2`; off-diagonal
terms are `2 Re(A_j conj(A_k))` and may be negative at favourable-interference
angles. Each term also carries its signed resistance density. The spectrum now
includes the same asymmetric strip-camber contribution as the default
resistance path rather than silently returning mean-thickness waves only.

`LowFroudeResistance::endpoint_pairs` exposes every retained upper-triangular
waterline term pair. Each entry identifies bow, stern, or interior knot;
coordinates; lambda power; complex coefficient; signed resistance and share;
quadrature estimate; and kernel work. This is attribution of the reduced
retained result, not of the bounded omitted submerged terms.

Regressions prove that:

- member amplitudes sum to the fleet amplitude;
- all self/pair terms sum to the total angular integrand and resistance density;
- an identical catamaran reproduces the classical `4 cos^2` factor;
- asymmetric spectrum integration reproduces default asymmetric resistance;
- endpoint-pair resistances, shares, and error estimates sum to their reported
  totals; and
- a full-multiplicity interior chine is labelled as an interior-knot source.

Classification: quadratic pair expansion and endpoint-wave decomposition are
class 1; the typed attribution APIs and asymmetric spectrum correction are
class 2. No priority claim is made.

### Workstream 4 — design-ranking stability: complete

The gate uses six exact degree-4 polynomial Wigley-family variants: base,
narrow/wide beam, finer/fuller longitudinal distribution, and a forward-LCB
perturbation. It verifies:

1. identical resistance ordering at `rel_tol = 1e-4`, `1e-6`, and `1e-8`;
2. identical ordering at `Fn=0.05` between dispatched endpoint reduction and
   forced general marching, with every pairwise margin difference inside the
   sum of the four reported absolute error estimates; and
3. identical geometry and ordering after exact degree-preserving knot insertion
   at mid-length and half-draft.

The first version was intentionally committed red as `5ed4f61`. Ordering did
not flip, but the base-versus-narrow-beam margin differed by `1.821e-10 N` while
the old combined estimate allowed only `1.069e-10 N`. The tail diagnostic fix
in `fedd88a` raises the combined allowance to the physically relevant tail
scale and closes the test without relaxing its assertion.

At tight general-marcher requests, variants honestly report `RefinementCap`
when the tail floor exceeds `rel_tol`. The phase-honest window reduces that
floor enough for the `Fn=0.05`, `rel_tol=1e-6` route comparison to converge,
while the `Fn=0.30`, `rel_tol=1e-8` cases remain honestly cap-limited. Rankings
are identical and every pairwise margin difference remains inside the sum of
the four reported absolute error estimates. The route-comparison test now
again requires `Converged`; the broader tolerance sweep accepts honest
`RefinementCap` results but rejects `TailCap` and `EvalCap` and never relabels
an unmet tolerance as convergence.

Classification: the gate, exact test-only knot insertion, and six-design corpus
are class 2. Stable ranking is demonstrated for this corpus, not generalized to
all hulls or operating points.

### New public API surface

| API | Purpose |
|---|---|
| `WaveMethod`, `WaveOutcome`, fields on `WaveResistance` | explicit route and termination status |
| `ControlNetGradient`, `PlacementGradient` | symmetric/asymmetric control and rigid-placement derivatives |
| `MemberWaveResistanceGradient`, `MultihullWaveResistanceGradient` | per-member derivative results |
| `multihull_wave_resistance_gradient[_with]` | exact fleet reverse pass |
| `ConstraintGradient`, `HullConstraintGradients` | volume, LCB, and wetted-area derivatives |
| `Hull::constraint_gradients` | constraint derivative entry point |
| `WaveSignature`, `WaveInterferenceContribution` | angular complex amplitude and pair ledger |
| `FreeWaveSpectrum::signature` | per-angle attribution entry point |
| `EndpointKind`, `EndpointWave`, `EndpointPairContribution` | low-Froude endpoint descriptors and pair results |
| `LowFroudeResistance::endpoint_pairs` | retained endpoint-pair attribution |

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

### Confirmed bug 3: low-Froude tail error estimate

The general outer marcher stops after one phase-width window is small relative
to the accumulated integral. At `Fn=0.05` it reported `5.123e-9` relative error
while its actual error against the independent analytic-Wigley reference was
`5.781e-7`, an underestimate by a factor of about 113. A single nearly cancelling
window is not a bound on the accumulated algebraic sequence of later windows.

The failing regression was committed in `406dcc3`. Commit `5a397f7` first fixed
the public single-hull result by attempting the independently bounded endpoint
solver and accepting it only when its combined bound/estimate meets the
requested tolerance; unsupported geometry or an insufficient bound falls back
to the general marcher. The route-selection regression in `db1c764` proves both
sides of this gate: default tolerance rejects the reduction at `Fn=0.08` and
accepts it at `Fn=0.05`.

The design-ranking gate later proved that the forced marcher's own diagnostic
was still too small for pairwise design margins. Red commit `5ed4f61` captures
that failure; `fedd88a` adds the algebraic-tail estimate described in Workstream
1. Independent review then found the first estimator optimistic at design
Froude numbers; assertion-red commit `30bc6ef` reproduces the 2.97× `Fn=0.35`
shortfall, and `8a7ac15` calibrates and checks the finite-λ safety factor across
`Fn=0.12`, `0.20`, and `0.35`. A further independent extension found the moved
boundary above `Fn≈0.5`; red commit `cc25d18` captures it. Commit `a2aa894`
separates stopping-window phase from vertical decay and sums tail/refinement
diagnostics. The resulting 12-point `Fn=0.02…1.00` regression covers every
measured actual error. It remains heuristic for general multihulls.

Classification: class 2.

### Suspected or unverified; no speculative fix

- The improved algebraic tail diagnostic is reproduced for the Wigley family
  from `Fn=0.02` through `Fn=1.00`, but adversarial multi-span and multihull
  beating envelopes have not been characterized. It is still a heuristic, and
  the bounded low-Froude route does not yet support those configurations.
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

The phase-honest stopping fix initially exposed a new cost already located in
that profiled hot loop: the vertical-decay rate kept forcing shrinking panels
after the associated submerged term was exponentially absent. Commit
`b79d46b` caps only this panel-sizing contribution once its exponent exceeds
eight; it does not discard the term. Same-session 30-sample medians were:

| Case | Uncapped phase-honest | Capped | Change |
|---|---:|---:|---:|
| 21-speed Wigley sweep | 72.090 ms | 12.663 ms | 82.4% faster |
| Forced marcher, `Fn=0.02` | 55.888 ms | 20.691 ms | 63.0% faster |
| Exact nine-control gradient | 6.924 ms | 0.820 ms | 88.2% faster |

Work at `Fn=0.02` fell from 1,378,528 to 511,504 inner evaluations. This is not
a bit-identical comparison because adaptive panel locations change: the sweep
checksums differ by `3.43e-10` relative. Accuracy did not regress—the actual
λ-to-4000 Wigley error decreased at every one of the 12 Froude numbers, the
minimum reported/actual margin remained 3.87×, and the Phase-0 harness passed.

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
wave patterns, and viscous/boundary-layer corrections. The actual thesis also
records the numerical history: Tuck used a Filon-like exact treatment for
piecewise-linear longitudinal data, Lazauskas extended that idea to exact
piecewise quadratics, and a 512-node outer angular rule was within 0.1% in his
tests *except at very low Froude number*. That exception is precisely the
regime attacked in Phase 4
([actual Adelaide thesis](https://digital.library.adelaide.edu.au/server/api/core/bitstreams/d054ebb4-f5a1-41ae-8bbf-678421fcfa80/content)).
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

### Low-speed endpoint asymptotics

The core physical asymptotic insight is old and unusually specific. Keller and
Ahluwalia's 1976 small-Froude analysis says that resistance depends on four
waterline quantities: hull and profile slopes at bow and stern, producing
longitudinal and transverse waves at each end
([fetched bibliographic record and abstract](https://trid.trb.org/View/66344)).
Wehausen's survey reports the Kotik expansion in endpoint slopes and notes the
additional contribution of corners
([actual ONR chapter](https://www.onr.navy.mil/media/document/wave-resistance-thin-ships)).
Gotman's fetched 2002 paper integrates polynomial longitudinal forms by parts,
separates nonoscillatory bow/stern self terms from their interference, and writes
the result in endpoint derivatives
([actual paper PDF](https://shipdesign.ru/Gotman/Study_of_Michells_Integral.pdf)).

These sources invalidate any broad claim that “endpoint reduction of Michell's
integral” is new. The code's exact finite per-span reduction reproduces and
generalizes that mechanism to tensor B-splines with interior knots; its new part
is the numerical treatment of the resulting *outer* endpoint-pair kernels and
the quantitative rule for safely omitting submerged endpoints.

There is a second, important boundary on interpretation. Modern exponential-
asymptotic work finds beyond-all-orders waves and Stokes switching in nonlinear
low-Froude free-surface flow, including smooth bodies whose relevant singularity
can lie at infinity
([Trinh and Chapman](https://arxiv.org/abs/1403.7182),
[Johnson-Llambias and Trinh](https://arxiv.org/abs/2402.03764)). That is not a
better quadrature for this linear Michell integral. It warns that the present
piecewise-polynomial endpoint dominance is a statement about this model and
geometry representation, not a universal low-speed law for the full nonlinear
ship-wave problem.

### Oscillatory quadrature and adjacent fields

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

Motygin has already applied both Levin collocation and steepest-descent plus
Clenshaw--Curtis to the oscillatory Kelvin Green-function integral in linear
ship-wave theory
([actual paper](https://arxiv.org/abs/1411.0321)). Thus “numerical steepest
descent in ship-wave hydrodynamics” is firmly class 1, although that Green-
function integral is not Michell resistance.

The obstacle in Michell resistance is the nonanalytic `|F(lambda)|^2`. The
decisive move was borrowed conceptually from high-frequency scattering: expose
known oscillatory phases and keep their amplitudes nonoscillatory. Expanding the
exact endpoint sum before squaring removes complex conjugation and produces
analytic pair integrals. Hybrid numerical-asymptotic boundary elements use the
same phase/amplitude separation to obtain frequency-independent cost
([Gibbs et al.](https://arxiv.org/abs/1912.09916)); automated contour methods now
handle coalescing endpoints and stationary points for general polynomial phases
([PathFinder](https://arxiv.org/abs/2307.07261)). These are analogies and future
machinery, not direct prior implementations of the solver here.

After `lambda=cosh(t)`, each endpoint-pair kernel is
`Ki_s(-i omega)`, the analytic continuation of a Bickley--Naylor function.
The standard integral definition and recurrence are documented by DLMF
([definition](https://dlmf.nist.gov/10.43.E3)). A June 2026 preprint expresses
all integer-order Bickley functions using a four-generator module of modified
Bessel and Struve functions
([actual preprint](https://arxiv.org/abs/2606.26415)). That new representation
may eventually replace contour quadrature with special-function evaluation,
but stability on a purely imaginary argument has not been established here.

The implemented contour starts from `lambda=1+t^2` and rotates
`t=exp(i*pi/4)y/sqrt(omega)`, turning the oscillation into Gaussian decay. Its
fixed 24/48-node work is therefore independent of `omega`. This realizes the
outer phase decomposition that the first research pass left as backlog.

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

- Exact inner integrals for every span of a polynomial tensor B-spline through
  the validated degree-16-per-direction envelope; no station/waterline sampling
  error in Michell amplitude inside that envelope.
- Stable origin, recurrence, and endpoint regimes for real, oscillatory, and
  complex-decay moments.
- Adaptive endpoint-regularized outer integration with diagnostics.
- Frequency-independent low-Froude endpoint-pair integration with an analytic
  bound on every omitted submerged term and error-gated automatic fallback.
- Coherent multihull phase superposition, heel via complex vertical decay, and
  experimental asymmetric/dipole paths.
- A matrix-free exact control gradient whose cost does not grow by one primal
  solve per design variable.

The exact inner integration and gradient match known published practice in their
respective areas. The low-Froude components reproduce known mathematics, while
their validated B-spline realization and routing are engineering improvements.

## Phase 4 — two validated advances

### Candidate ranking

| Rank | Candidate | Leverage | Risk | Decision |
|---:|---|---|---|---|
| 1 | Phase-decomposed low-Froude outer integral | Makes the numerically hostile regime cheap and bounded | High; squared modulus, infinite tail, branch endpoint | Implemented and validated |
| 2 | Exact `dR_w/dP` for B-spline controls | Enables deterministic hull optimization; natural quadratic structure | Low–medium | Implemented and validated |
| 3 | Finite-depth/channel kernel | Large practical scope increase | High; critical modes and new validation corpus | Backlog |
| 4 | Boundary-layer tangency correction | Improves physical prediction beyond slender hulls | Medium–high; needs empirical/CFD calibration | Backlog |
| 5 | Neumann–Michell surface solver | Higher-fidelity linear theory | Very high; a separate BIE solver | Backlog |

### Advance A: exact control-net gradient

#### Method

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

The initial API at `64925d4` intentionally accepted only symmetric hulls. The
design-tool upgrade now returns separate physical port/starboard gradients for
the default strip closure and still refuses to imply coverage of the optional
solved-lifting closure.

#### Validation and claimed advantage

| Check | Result |
|---|---|
| Every one of 9 controls versus centered finite differences | pass; scaled error below `2e-6` |
| Constant half-breadth shift (null direction of `df/dx`) | pass |
| Returned primal versus ordinary resistance API | bit-for-bit equal in test |
| Analytic work versus control count | one primal convergence plus one reverse pass |
| Full Phase-0 harness | pass |
| Full Rust workspace | 215 test cases pass, including one doctest |

Final release benchmark, 30 samples, default tolerance:

| Gradient method | Median | Best |
|---|---:|---:|
| Exact reverse, 9 controls | 0.820 ms | 0.808 ms |
| Centered finite differences | 14.431 ms | 14.156 ms |

Median speedup: **17.59×**. Aggregate absolute-gradient checksums agree to about
`8.1e-12` relative (`1.734433595800e4` versus `1.734433595786e4`, including
the benchmark's 30-run accumulation).

Classification: the quadratic derivative is class 1; this matrix-free B-spline
implementation and API are class 2. No class-4 claim.

### Advance B: bounded endpoint/steepest-descent solver at low Froude

#### Derivation

On each polynomial x span, repeated integration by parts is finite and exact:

```text
integral p(x) exp(i nu lambda x) dx
  = sum over span endpoints of
      c_(endpoint,n) exp(i nu lambda x_endpoint) / lambda^n.
```

The exact z moments attach `exp(-nu lambda^2 z_endpoint)` to those terms. At low
Froude number (`nu` large), every endpoint below the waterline is exponentially
small. Retaining the `z=0` terms turns the full inner amplitude into a finite sum
of endpoint waves. Expanding its squared modulus gives pairwise kernels

```text
K_s(omega) = integral_1^infinity
  exp(i omega lambda) / (lambda^s sqrt(lambda^2 - 1)) d lambda.
```

`lambda=cosh(t)` identifies `K_s(omega)=Ki_s(-i omega)`. For nonzero frequency,
`lambda=1+t^2` removes the square-root endpoint and the exact contour rotation
`t=exp(i*pi/4)y/sqrt(|omega|)` makes the phase Gaussian-decaying. Zero-frequency
self terms use the beta-function recurrence exactly. Coarse/fine contour rules
estimate quadrature error; triangle inequality plus analytic `K_s(0)` bounds
every pair involving an omitted submerged endpoint.

The current public specialization deliberately requires an upright symmetric
single hull, both spline degrees in the validated range 1--16, and every
distinct active endpoint frequency at least 25. For `m` equal longitudinal
spans with active adjacent endpoints, that last gate gives
`m <= 1/(25 Fn^2)`. An illustrative 360-pair retained map uses 25,920 ordinary
contour nodes, versus 1,152 for Wigley's 16 pairs. It does not silently
approximate asymmetric, heeled, multihull, closely spaced, or moderate-frequency
cases. `wave_resistance_with` accepts it only when the combined estimate meets
the caller's tolerance and otherwise uses the general solver.

The existing harness contains Wigley and one synthetic chine geometry, not a
curated representative multi-span design corpus. An acceptance-rate study
would therefore measure an arbitrary sampling choice rather than practical
coverage; the narrow scope and exact span arithmetic are reported instead.

#### Independent validation

The reference does not call production moments or outer quadrature. It uses the
closed-form Wigley inner amplitude, 16-point Gauss panels on the positive real
axis with at most `pi/2` bow-phase advance per panel, and `lambda_max=500`.

| Fn | Reference resistance (N) | Reduced relative difference | Reduced estimate | General marcher difference | Work: marcher / reduced |
|---:|---:|---:|---:|---:|---:|
| 0.08 | `1.945940872497e-1` | `1.172e-5` | `3.814e-5` | `6.917e-8` | 56,032 / 1,152 |
| 0.05 | `1.057847519834e-2` | `1.058e-11` | `3.728e-12` | `1.493e-7` | 119,936 / 1,152 |
| 0.03 | `5.113923599287e-4` | `2.328e-11` | `3.462e-12` | `3.322e-7` | 267,856 / 1,152 |
| 0.02 | `4.417234459536e-5` | `3.558e-11` | `7.620e-13` | `6.333e-7` | 511,504 / 1,152 |

At `Fn=0.08` the estimate correctly refuses default-tolerance dispatch. At
`Fn<=0.05` the observed `1e-11`-scale discrepancies exceed the solver's analytic
omission bound because the independent reference has a finite real-axis tail;
tests therefore include a separately justified `2e-10` reference floor. The
contour part is checked independently against dense real-axis integration for
orders `s=4,7,10,50,98,128` and frequencies `25,100,400`. The public degree
envelope reaches `s=98`; `s=128` provides margin. This extension first exposed
a `5.13e-6` relative error at `s=98`, `omega=25` in the old 48-point result.
The red regression is commit `a4d011a`; commit `69f9ff4` selects a 48/96 rule
when `s/|omega| >= 2`, reducing the same independent-reference discrepancy
below `1e-9` while retaining the 24/48 rule for Wigley. The exact endpoint
decomposition is also checked against production exact moments for both Wigley
and a multi-span full-multiplicity chine hull.

The total reported estimate is not an interval proof: the omitted-term part is
analytically bounded, while the selected 24/48- or 48/96-node contour difference
is an empirical quadrature estimate. “Bounded” in this report refers to the
discarded submerged physics terms, not to a formally certified floating-point
result.

#### Claimed advantage

Final 30-sample release benchmark:

| Case | Median | Best | Work/diagnostics |
|---|---:|---:|---|
| General marcher, `Fn=0.02` | 20.691 ms | 20.561 ms | 511,504 inner evaluations |
| Endpoint/NSD, `Fn=0.02` | 0.029 ms | 0.028 ms | 1,152 kernel evaluations |
| Default API, `Fn=0.05` | 0.029 ms | 0.028 ms | estimate `3.728e-12` |
| 21 speeds, `Fn=0.10…0.50` | 12.663 ms | 12.463 ms | checksum `2.553251156101e5` |

The direct low-Froude speedup is **713.48×** at `Fn=0.02`. The 21-speed
production sweep remains on the general method where appropriate; its changed
checksum is the corrected positive tail, independently checked above.

Classification: endpoint integration by parts, endpoint low-speed dominance,
Bickley functions, and numerical steepest descent are class 1. The Rust solver,
error-bound routing, tests, and diagnostics are class 2.

#### Literature due diligence and reclassification

Searches were run across general web indexing, arXiv, DOI/publisher pages,
TRID, the Adelaide repository, DLMF, and references inside the fetched Tuck,
Lazauskas, Gotman, Motygin, and oscillatory-quadrature papers. Queries included:

- `"Michell integral" "Bickley-Naylor"`
- `"ship wave resistance" "Bickley function"`
- `"Ki_n" "wave resistance" thin ship`
- `"Michell wave resistance" "numerical steepest descent"`
- `"Michell integral" endpoint asymptotic numerical method`
- `"endpoint decomposition" Michell integral B-spline`
- `thin-ship steepest descent wave resistance`
- `Michell wave resistance Fresnel integral bow stern interference low Froude`

Equation-level follow-up established a broader analytic lineage than the first
pass recognized. Birkhoff and Kotik separate hull data from a reusable kernel;
Michelsen's 1960 dissertation reduces polynomial hull functions to tabulatable
special-function expressions; its verified 1972 JSR record describes a finite
Gegenbauer double sum; and the verified Sendagorta--Grases record describes
rapidly convergent, shape-separated Michell/Havelock series for design use.
Gotman supplies endpoint-derivative structure, Motygin supplies ship-wave
steepest descent, and Keller--Ahluwalia supplies low-speed endpoint dominance.
The full Michelsen 1972 and Sendagorta--Grases papers remain interlibrary-loan
due-diligence items, but their records already justify withdrawing the earlier
priority-seeking classification. The method lineage is class 1/class 3; the
validated B-spline endpoint implementation, cancellation accounting, omission
bound, contour evaluation, error-gated dispatch, and fallback are class 2.

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
| `3ceb6ec` | Initial frontier report and results |
| `73452a1` | Isolated formatter-only cleanup |
| `12336ec` | Exact endpoint/Bickley/steepest-descent low-Froude solver |
| `406dcc3` | Failing low-Froude tail-estimate regression |
| `5a397f7` | Tolerance-bounded low-Froude dispatch and bug fix |
| `db1c764` | Acceptance/fallback routing regression |
| `e6568d4` | Remove rustdoc ambiguity from the new public API |
| `6b8c42d` | Failing tail/evaluation-cap status regressions |
| `ca2d648` | Explicit method/outcome status, CLI propagation, and `.msw` v2 |
| `45d8c36` | Multihull placement and asymmetric control-net gradients |
| `2ff748f` | Displacement, LCB, and wetted-area gradients |
| `3e921b6` | Per-member angular wave signatures and interference attribution |
| `e8f952e` | Low-Froude endpoint-pair attribution |
| `5ed4f61` | Failing six-design ranking-margin gate |
| `fedd88a` | Algebraic marcher-tail diagnostic and tail-limited refinement stop |
| `282f0f0` | Isolated formatter-only cleanup |
| `8c338e0` | Remove rustdoc ambiguity from new unit annotations |
| `30bc6ef` | Failing design-Froude error-estimate coverage regression |
| `8a7ac15` | Calibrate the finite-λ tail estimate and document cap semantics |
| `cc25d18` | Failing high-Froude error-estimate coverage regressions |
| `a2aa894` | Separate quiet-window phase from decay and sum error sources |
| `b79d46b` | Cap vanished depth-envelope panel rates with measured speedups |
| `b2cb256` | Failing adversarial high-degree regression suite and preserved reviewer probes |
| `0ddda93` | Degree-16 support cap and non-finite-result refusal |
| `84e2885` | Shared coefficient validation and endpoint-cancellation error accounting |
| `7b824f5` | Exact degree-elevation validation of the supported envelope |

## Ranked backlog

1. Generalize the endpoint expansion to multihulls. The transverse placement
   phase introduces `lambda*sqrt(lambda^2-1)` and new saddles; use Motygin and
   PathFinder-style contour topology rather than assuming the present rotation.
2. Retain submerged endpoints rather than bound-and-drop them. Their mixed
   linear/quadratic complex phase has moving saddles; derive a uniform contour
   for shallow endpoints and phase coalescence so the method reaches moderate
   Froude number smoothly.
3. Replace Bickley contour quadrature with a cache/recurrence or the 2026
   Bessel--Struve module after validating stability for imaginary arguments.
   This may reduce 1,152 kernel evaluations further and amortize speed sweeps.
4. Make the entire result rigorous with directed rounding or ball arithmetic,
   a proven Gaussian-contour quadrature remainder, and a certified independent
   reference. The present submerged-term bound is rigorous; contour error is not.
5. Differentiate through the endpoint solver so low-Froude exact gradients get
   both advances at once. Extend the solved-lifting closure only after deriving
   and validating its own adjoint; do not silently reuse the strip derivative.
6. Stress the general tail diagnostic with adversarial multispan and multihull
   beating envelopes. Replace the power-law heuristic with a certified or
   envelope-aware bound if practical.
7. Add speed, draft/waterline, fairness, and curvature derivatives. A production
   knot-insertion API would also make adaptive design parametrisations easier;
   current knot insertion exists only as an exact ranking regression helper.
8. Add a small constrained optimizer example with non-negativity, closure,
   displacement, and curvature regularization; do not optimize wave resistance
   alone because the published problem is ill-posed. This remains deliberately
   out of the current diff.
9. Implement finite-depth/infinite-width first, then finite-width channel modes.
   Near critical depth the saddle/mode structure changes, so reuse the uniform-
   asymptotic contour work rather than treating this as a kernel substitution.
10. Add boundary-layer displacement/tangency corrections behind an explicit
   model option and validate against the five-hull literature set.
11. Harden lifting solvers with public dimension/geometry validation and
   condition estimates instead of debug-only singularity assertions.
12. Treat Neumann–Michell as a separate higher-fidelity solver sharing the same
   B-spline geometry, benchmarks, and published Wigley cases.

## Final validation

| Command | Result |
|---|---|
| `cargo build --workspace --all-targets` | pass |
| `cargo test --workspace` | pass: 226 test cases including one doctest; 0 failed |
| `uv run --with pytest --with numpy pytest -q` in `python/` | pass: 11 passed in 0.10 s |
| `cargo test --release -p michell --test high_degree_hardening -- --nocapture` before K2 | expected red: all five independent regressions failed |
| `cargo test --release -p michell --test degree_envelope -- --nocapture` | pass: degree 2–16 endpoint/marcher sweeps covered; degree 17 refused |
| `cargo test -p michell --test ranking_stability ordering_is_invariant_under_exact_knot_insertion -- --exact` | pass: unchanged knot-insertion gate |
| `MICHELL_BENCH_SAMPLES=30 cargo bench -p michell --bench wigley` | pass; sweep checksum exactly `2.553251156101e5`; exact-gradient/FD checksums agree to `8.1e-12` relative |
| `cargo doc -p michell --no-deps` | generated successfully; pre-existing broken-link/unit-bracket warnings remain |
| allowed-lint Clippy command reproduced below | pass |
| `cargo fmt --all --check` | pass |
| `git diff --check` | pass |

The structural stopping fix deliberately changes resistance values at the old
`1e-7`-scale truncation-error level; the checksum change is therefore expected,
not a performance-regression artifact. The new values are checked directly
against the independent λ-to-4000 reference rather than assumed equivalent to
the old checksum.

The exact passing Clippy command was:

```text
cargo clippy -p michell --all-targets -- -D warnings \
  -A clippy::neg-cmp-op-on-partial-ord \
  -A clippy::redundant-closure \
  -A clippy::needless-borrows-for-generic-args \
  -A clippy::unnecessary-cast
```

No database-backed tests exist in this repository. No test was skipped. The
pre-existing untracked `python/uv.lock` was preserved outside this worktree and
excluded from every commit. No remote push was made.

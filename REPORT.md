# Michell wave-resistance deep dive

Date: 2026-08-01  
Branch: `story-wave-resistance-frontier`  
Starting revision: `abe2b1f`  
Remote operations: none

## Executive result

The project now has a published-value validation anchor, physical and numerical
property checks, a repeatable benchmark harness, three reproduced bugs with
fixes, a measured general-regime quadrature optimization, an exact control-net
gradient, and a new low-Froude solver.

The second Phase-4 advance is the larger one. It rewrites each polynomial
B-spline span exactly as endpoint waves, discards only depth-damped endpoints
under an explicit absolute bound, expands the squared amplitude into pairwise
kernels, and evaluates those kernels on a Gaussian-decaying steepest-descent
contour. At `Fn=0.02` the result needs 1,152 kernel evaluations instead of
877,424 inner-amplitude evaluations and is **971.79 times faster** (0.029 ms
versus 28.141 ms median). Against an independent real-axis reference its
relative difference is `3.56e-11`; the reference's finite tail, rather than the
new solver, limits that comparison. Cost is effectively independent of the
oscillation frequency in the tested low-Froude range.

The constituent mathematics is not new: endpoint low-speed asymptotics,
Bickley--Naylor functions, and numerical steepest descent all have substantial
literatures. The implementation is class 2. The *combination* of an exact
arbitrary-degree B-spline endpoint reduction, analytically continued Bickley
kernels, frequency-independent contour quadrature, and a computable omitted-
endpoint bound in Michell resistance is classified only as **possibly novel**
(class 4), after the documented searches below found no prior instance. This is
not a proof of priority and is not called a breakthrough.

The first Phase-4 advance, the exact gradient, remains useful: it is 14.76 times
faster than centered finite differences for the nine-control benchmark. Its
quadratic structure is known; the matrix-free B-spline reverse pass is an
engineering implementation, not a novelty claim.

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
| Low-speed dominance by waterline bow/stern data | 1 | Keller–Ahluwalia, Wehausen/Kotik, and Gotman endpoint results reproduced computationally |
| Steepest-descent evaluation of oscillatory ship-wave integrals | 1 | Motygin and the general numerical-steepest-descent literature |
| Bounded low-Froude solver and automatic fallback | 2 | New implementation, independent reference tests, and routing regression |
| Exact B-spline endpoint/Bickley/steepest-descent Michell reduction | 4, qualified | No prior instance found by the searches recorded below; priority is unproved |

The gradient searches (“analytic gradient Michell wave resistance”, “Michell
integral shape derivative”, “B-spline hull optimization Michell”, and “adjoint
thin-ship wave resistance”) found prior quadratic-form and sensitivity work, so
the gradient is not class 4. The low-Froude search was broader and found close
precursors, discussed explicitly below rather than hidden behind a novelty
label.

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

The failing regression was committed in `406dcc3`. Commit `5a397f7` fixes the
public result by first attempting the independently bounded endpoint solver and
accepting it only when its combined bound/estimate meets the requested tolerance;
unsupported geometry or an insufficient bound falls back to the general marcher.
The route-selection regression in `db1c764` proves both sides of this gate:
default tolerance rejects the reduction at `Fn=0.08` and accepts it at `Fn=0.05`.

Classification: class 2. The general multihull marcher's diagnostic remains a
heuristic, so callers that specifically request that path still need the explicit
convergence-status work below.

### Suspected or unverified; no speculative fix

- The outer tail stop's failure is reproduced for the low-Froude Wigley case,
  but adversarial multi-span and multihull interference envelopes have not been
  characterized. The bounded low-Froude route does not yet support them.
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

- Exact inner integrals for every span of an arbitrary-degree polynomial
  tensor B-spline; no station/waterline sampling error in Michell amplitude.
- Stable origin, recurrence, and endpoint regimes for real, oscillatory, and
  complex-decay moments.
- Adaptive endpoint-regularized outer integration with diagnostics.
- Frequency-independent low-Froude endpoint-pair integration with an analytic
  bound on every omitted submerged term and conservative automatic fallback.
- Coherent multihull phase superposition, heel via complex vertical decay, and
  experimental asymmetric/dipole paths.
- A matrix-free exact control gradient whose cost does not grow by one primal
  solve per design variable.

The exact inner integration and gradient match known published practice in their
respective areas. The low-Froude components also reproduce known mathematics;
only their specific composition is the qualified class-4 claim. No exhaustive
search can establish priority, and no patent search was performed.

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

The initial API intentionally accepts only symmetric hulls. An asymmetric hull
has two physical nets and multiple dipole closures; returning a derivative of
only its stored symmetric mean would be misleading.

#### Validation and claimed advantage

| Check | Result |
|---|---|
| Every one of 9 controls versus centered finite differences | pass; scaled error below `2e-6` |
| Constant half-breadth shift (null direction of `df/dx`) | pass |
| Returned primal versus ordinary resistance API | bit-for-bit equal in test |
| Analytic work versus control count | one primal convergence plus one reverse pass |
| Full Phase-0 harness | pass |
| Full Rust workspace | 196 test cases pass, including one doctest |

Release benchmark, 30 samples, default tolerance (latest run):

| Gradient method | Median | Best |
|---|---:|---:|
| Exact reverse, 9 controls | 1.426 ms | 1.377 ms |
| Centered finite differences | 21.046 ms | 20.761 ms |

Median speedup: **14.76×**. Aggregate absolute-gradient checksums agree to about
`2.4e-11` relative (`1.734433365029e4` versus `1.734433364988e4`, including
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
single hull, x degree at least one, and sufficiently separated endpoint phases.
It does not silently approximate asymmetric, heeled, multihull, or moderate-
frequency cases. `wave_resistance_with` accepts it only when the combined
estimate meets the caller's tolerance and otherwise uses the general solver.

#### Independent validation

The reference does not call production moments or outer quadrature. It uses the
closed-form Wigley inner amplitude, 16-point Gauss panels on the positive real
axis with at most `pi/2` bow-phase advance per panel, and `lambda_max=500`.

| Fn | Reference resistance (N) | Reduced relative difference | Reduced estimate | General marcher difference | Work: marcher / reduced |
|---:|---:|---:|---:|---:|---:|
| 0.08 | `1.945940872497e-1` | `1.172e-5` | `3.814e-5` | `3.125e-7` | 114,432 / 1,152 |
| 0.05 | `1.057847519834e-2` | `1.058e-11` | `3.728e-12` | `5.781e-7` | 235,808 / 1,152 |
| 0.03 | `5.113923599287e-4` | `2.328e-11` | `3.462e-12` | `1.063e-6` | 501,904 / 1,152 |
| 0.02 | `4.417234459536e-5` | `3.558e-11` | `7.620e-13` | `2.024e-6` | 877,424 / 1,152 |

At `Fn=0.08` the estimate correctly refuses default-tolerance dispatch. At
`Fn<=0.05` the observed `1e-11`-scale discrepancies exceed the solver's analytic
omission bound because the independent reference has a finite real-axis tail;
tests therefore include a separately justified `2e-10` reference floor. The
contour part is checked independently against dense real-axis integration for
orders `s=4,7,10` and frequencies `25,100,400`. The exact endpoint decomposition
is also checked against production exact moments for both Wigley and a
multi-span full-multiplicity chine hull.

The total reported estimate is not an interval proof: the omitted-term part is
analytically bounded, while the 24/48-node contour difference is an empirical
quadrature estimate. “Bounded” in this report refers to the discarded submerged
physics terms, not to a formally certified floating-point result.

#### Claimed advantage

Latest 30-sample release benchmark:

| Case | Median | Best | Work/diagnostics |
|---|---:|---:|---|
| General marcher, `Fn=0.02` | 28.141 ms | 27.794 ms | 877,424 inner evaluations |
| Endpoint/NSD, `Fn=0.02` | 0.029 ms | 0.028 ms | 1,152 kernel evaluations |
| Default API, `Fn=0.05` | 0.029 ms | 0.028 ms | estimate `3.728e-12` |
| 21 speeds, `Fn=0.10…0.50` | 21.897 ms | 21.584 ms | unchanged checksum |

The direct low-Froude speedup is **971.79×** at `Fn=0.02`. Relative to the
pre-dispatch `Fn=0.05` baseline of 7.494 ms from the same development run, the
default API is about **268×** faster. The 21-speed production sweep remains on
the general method where appropriate and retains its numerical checksum.

Classification: endpoint integration by parts, endpoint low-speed dominance,
Bickley functions, and numerical steepest descent are class 1. The Rust solver,
error-bound routing, tests, and diagnostics are class 2. Their specific combined
Michell/B-spline construction is the qualified class-4 claim described next.

#### Novelty falsification log

Searches were run across general web indexing, arXiv, DOI/publisher pages,
TRID, the Adelaide repository, DLMF, and references inside the fetched Tuck,
Lazauskas, Gotman, Motygin, and oscillatory-quadrature papers. Exact queries
that did **not** find this combination included:

- `"Michell integral" "Bickley-Naylor"`
- `"ship wave resistance" "Bickley function"`
- `"Ki_n" "wave resistance" thin ship`
- `"Michell wave resistance" "numerical steepest descent"`
- `"Michell integral" endpoint asymptotic numerical method`
- `"endpoint decomposition" Michell integral B-spline`
- `thin-ship steepest descent wave resistance`
- `Michell wave resistance Fresnel integral bow stern interference low Froude`

The search *did* find substantial near-prior art: Gotman's endpoint-derivative
series for polynomial, separable hulls; de Sendagorta and Grases's 1988 abstract
describing rapidly convergent Michell/Havelock series and tabulatable velocity
functions ([fetched record](https://trid.trb.org/View/397494)); Motygin's
steepest-descent Kelvin Green function; Keller--Ahluwalia low-speed endpoints;
and the Bickley literature. The full de Sendagorta--Grases article was not
available for inspection, so it is a material uncertainty, explicitly not
silence that proves novelty. No source found all four elements: arbitrary-degree
B-spline span endpoints, Bickley analytic continuation of the pair kernel,
fixed-cost contour evaluation, and a submerged-endpoint omission bound with
automatic tolerance routing.

Accordingly, class 4 means only “possibly novel implementation-level
combination after a serious but non-exhaustive search.” Establishing priority
would require a professional database and patent search plus expert review.

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
5. Add explicit convergence status (`converged`, `tail_cap`, `evaluation_cap`)
   to the legacy marcher and make capped runs impossible to mistake for
   converged values. Build adversarial multihull interference cases.
6. Differentiate through the endpoint solver so low-Froude exact gradients get
   both advances at once; then extend gradients to port and starboard nets,
   multihull placements, speed, and constrained objectives (volume, wetted
   area, fairness).
7. Add a small constrained optimizer example with non-negativity, closure,
   displacement, and curvature regularization; do not optimize wave resistance
   alone because the published problem is ill-posed.
8. Implement finite-depth/infinite-width first, then finite-width channel modes.
   Near critical depth the saddle/mode structure changes, so reuse the uniform-
   asymptotic contour work rather than treating this as a kernel substitution.
9. Add boundary-layer displacement/tangency corrections behind an explicit
   model option and validate against the five-hull literature set.
10. Harden lifting solvers with public dimension/geometry validation and
   condition estimates instead of debug-only singularity assertions.
11. Treat Neumann–Michell as a separate higher-fidelity solver sharing the same
   B-spline geometry, benchmarks, and published Wigley cases.

## Final validation

| Command | Result |
|---|---|
| `cargo build --workspace` | pass |
| `cargo build --workspace --release` | pass |
| `cargo test --workspace` | pass: 196 test cases including one doctest; 0 failed |
| `uv run --with pytest --with numpy pytest -q` in `python/` | pass: 11 passed in 0.10 s |
| `MICHELL_BENCH_SAMPLES=30 cargo bench -p michell --bench wigley` | pass; final numbers recorded above |
| `cargo doc -p michell --no-deps` | generated successfully; pre-existing broken-link warnings remain outside the new API |
| strict `cargo clippy -p michell --all-targets -- -D warnings` | does not pass: seven pre-existing current-Clippy findings in `body.rs`, `lifting3d.rs`, `lifting.rs`, and `tests/inclined.rs` |
| same Clippy command with the four named pre-existing lints allowed | pass |
| `cargo fmt --all --check` | does not pass: current rustfmt disagrees with already committed formatting in several files; no formatter rewrite was applied |
| `git diff --check` | pass |

The exact passing Clippy command was:

```text
cargo clippy -p michell --all-targets -- -D warnings \
  -A clippy::neg-cmp-op-on-partial-ord \
  -A clippy::redundant-closure \
  -A clippy::needless-borrows-for-generic-args \
  -A clippy::unnecessary-cast
```

No database-backed tests exist in this repository. No test was skipped. No
remote push was made.

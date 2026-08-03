# Pre-registered Insel-Wigley acceptance criteria

Committed before any `michell` prediction for this study is computed. These
rules will be applied unchanged to the frozen digitized data. If a rule proves
ill-posed, the finding will be reported as unscored; it will not be replaced by
a more favourable post-hoc rule.

## Scope and invariants

The primary target is the fixed-attitude C2-FX series. The monohull is scored
first, followed independently by the four catamarans at `S/L = 0.2, 0.3, 0.4,
0.5`. Every computation uses the physical 1.800 m model in 15-degree freshwater
with `rho = 1000 kg/m^3`, `g = 9.80665 m/s^2`, and the thesis geometry recorded
in `METHOD.md`. No full-scale extrapolation is allowed.

The isolated-hull coefficient uses `S_W = 0.482 m^2`. A catamaran coefficient
uses total static wetted area `2 S_W = 0.964 m^2`; this convention makes a
non-interacting pair have the same coefficient as the isolated demihull. The
two hull placements are `y = -S/2` and `y = +S/2`, where `S` is demihull
centreline spacing.

The harness must call `multihull_wave_resistance_with` for every catamaran,
request `rel_tol = 1e-6`, and fail unless every outcome is `Converged`. Every
raw row must record resistance, coefficient, library interference, method,
outcome, `est_rel_error`, and evaluation count. Failure of any invariant is a
harness failure, not an experimental disagreement.

## Deterministic experimental curve

The source markers are unsnapped and include replicates. Comparisons that need
a continuous experimental curve use this fixed smoother:

```text
E(q) = sum_i w_i(q) v_i / sum_i w_i(q)
w_i(q) = exp[-0.5 ((q - Fn_i) / 0.0125)^2]
```

Only markers within `|q - Fn_i| <= 0.0375` contribute, and a value is emitted
only when at least two markers contribute. No extrapolation is permitted. The
scoring grid is `Fn = 0.25, 0.275, ..., 0.80`; interference is scored only
through `Fn = 0.55`, where both numerator and denominator are resolved. A
`0.001` grid is used only to locate extrema. Digitization uncertainty is not
used to move a central value, but is shown on plots and included in the
interpretation.

Predictions are linearly interpolated from the committed fine-grid output. No
fit parameter, horizontal shift, vertical scale, or empirical correction may
be applied to a prediction.

## 1. Wave-pattern coefficient

At every valid experimental scoring point, define

```text
d(q) = |C_W,michell(q) - C_WP,experiment(q)|
       / max(|C_WP,experiment(q)|, 0.0005).
```

The floor prevents small low-speed wave ordinates from dominating a percentage
score. Each configuration receives a pointwise grade:

| grade | median `d` | 80th percentile `d` |
|---|---:|---:|
| agreement | <= 0.25 | <= 0.45 |
| partial | <= 0.40 | <= 0.75 |
| disagreement | otherwise | otherwise |

Shape is scored over the principal-wave range. The experimental and predicted
principal hump are the global maxima on `0.38 <= Fn <= 0.58`, evaluated on the
`0.001` grid. Relative amplitude error uses the experimental hump value as its
denominator, with the same `0.0005` floor as `d`. Hump position and amplitude
receive separate grades:

| grade | `|delta Fn_hump|` | relative amplitude error |
|---|---:|---:|
| agreement | <= 0.020 | <= 0.25 |
| partial | <= 0.040 | <= 0.45 |
| disagreement | otherwise | otherwise |

The coefficient verdict for a configuration is agreement only if all three
grades agree; it is disagreement if any grade disagrees; otherwise it is
partial. The five fixed configurations are reported separately so a good
catamaran score cannot hide a bad monohull baseline.

## 2. Wave-interference ratio

Insel's decomposition defines `tau` from resistance after a fitted viscous
multiplier `beta`. The marker archive does not supply `beta`, and Figures
155--158 plot several derived curves without raw markers. Therefore the study
will not manufacture a resistance-decomposition `tau` from `C_T`.

The directly observable experimental comparison is instead

```text
tau_WP(q, S/L) = C_WP,cat(q, S/L) / C_WP,mono(q).
```

On the documented total-wetted-area coefficient convention, unity means no
wave-pattern interference and this ratio is dimensionally aligned with the
library's `interference` output. It will always be labelled `tau_WP`, never
silently relabelled as Insel's resistance-decomposition `tau`. Points are
scored only when both smoothed coefficients exist and the monohull denominator
is at least `0.0005`.

Three properties are scored for each separation on `0.25 <= Fn <= 0.55`:

1. **Level:** median and 80th-percentile absolute errors
   `|interference_michell - tau_WP|`. Agreement requires <= 0.25 and <= 0.45;
   partial requires <= 0.40 and <= 0.70; otherwise disagreement.
2. **Sign:** values within `|tau_WP - 1| <= 0.10` are experimentally neutral
   and excluded from sign scoring. A prediction is constructive above 1.05,
   destructive below 0.95, and neutral otherwise. Agreement requires at least
   80% matching signs; partial requires at least 60%; otherwise disagreement.
3. **Hump/hollow:** the dominant constructive hump is the maximum on
   `0.38 <= Fn <= 0.52`; its preceding hollow is the minimum on
   `0.30 <= Fn <= 0.42`. Each feature is scored only if its experimental
   prominence from unity is at least 0.10. Location agreement is
   `|delta Fn| <= 0.020` (partial <= 0.040); amplitude agreement is absolute
   ratio error <= 0.25 (partial <= 0.50). A missing predicted feature is a
   disagreement.

The interference verdict is agreement when all applicable properties agree,
disagreement when any applicable property disagrees, and partial otherwise.
All component scores remain visible.

## 3. Trend with separation

At each `0.025` scoring point where all four experimental ratios exist, rank
the four separations by `tau_WP` and by predicted interference. The pointwise
Spearman correlation is computed with average ranks for ties. The trend grade
is agreement if the median correlation is at least 0.70, partial if at least
0.40, and disagreement otherwise.

A second envelope check computes, for each separation,

```text
A(S/L) = RMS[tau(q, S/L) - 1],  0.40 <= Fn <= 0.55.
```

The predicted and experimental four-value `A` rankings are scored by Spearman
correlation with the same thresholds. This tests decay of interference with
separation without assuming that the instantaneous phase ordering is
monotone.

## 4. Trust-envelope rule

A fixed-attitude `(Fn, S/L)` cell on the `0.025` grid is provisionally
design-grade only when all of the following hold:

- monohull normalized coefficient error `d <= 0.35`;
- catamaran normalized coefficient error `d <= 0.35`;
- absolute interference error is <= 0.35;
- if the experiment is non-neutral, predicted and experimental interference
  signs agree; and
- the solve is `Converged` with a finite recorded error estimate.

The final trust envelope is the union of passing cells. A quoted continuous
interval must contain at least three consecutive grid points with no skipped
or failing cell. Free-attitude runs can never enter this envelope.

## 5. Secondary free-attitude comparison

Free C2 values are plotted and summarized but not given a fixed-attitude grade.
Before prediction, the expected qualitative result is increasing divergence
above approximately `Fn = 0.4`, because the present calculation neither solves
running trim and sinkage nor changes immersed geometry. The comparison will
report median coefficient discrepancy below and above `Fn = 0.4`; it will not
be used to tune the fixed model or its score.

## Expected discrepancies recorded before comparison

These are hypotheses, not exemptions from the numerical scores:

- `C_WP` can lie systematically below true wave resistance because of viscous
  attenuation between hull and probes, finite longitudinal record, finite
  harmonic/matrix resolution, and unresolved small wave heights.
- The supplied brief's “transverse-cut truncation” explanation is inapplicable:
  Insel used a multiple-longitudinal-cut matrix method. It will not be cited as
  a post-hoc explanation.
- Linear thin-ship theory may exaggerate hump/hollow amplitudes for this
  `B/L = 0.1` hull and cannot represent nonlinear or breaking waves.
- The library contains no viscous-interference model (`beta`). This limits a
  comparison to Insel's resistance decomposition but does not excuse an error
  in the directly comparable wave-pattern ratio.
- Free-to-trim discrepancies are expected to grow once running attitude changes
  materially, approximately above `Fn = 0.4`.
- The highest speeds carry documented acceleration, shallow-water, breaking-
  wave, and measurement cautions. They remain plotted; they are not silently
  deleted from the coefficient score through `Fn = 0.8`.

Report 72 is methodological precedent for different NPL hulls, not an exact C2
comparison. No numerical “Report 72 match” will be claimed for C2. This source
mismatch is an open item rather than permission to invent a published baseline.

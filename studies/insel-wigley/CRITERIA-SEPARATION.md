# Preregistered separation-definition probe

Committed after the complete attributed source family and before any
separation-grid prediction or fit is produced. This probe does not alter the
frozen kernel or the agreement thresholds in `CRITERIA-THEORY.md`.

## Prediction grid and objective

For the exact C2 Wigley demihull, the validated library solver evaluates
centreline separation

```text
S_eff/L = 0.080, 0.085, ..., 0.550
```

at `0.005` spacing. Each separation uses `Fn = 0.15, 0.155, ..., 1.00`, the
existing model-scale geometry, fluid properties, and `rel_tol = 1e-6`. Every
pair and standalone solve must converge.

For source Figure `j`, let `D_j` be its admitted digitized C2 points over the
source curve's full archived range. At candidate separation `s`, linearly
interpolate the library prediction in `Fn` and compute

```text
J_j(s) = median over D_j of |tau_source - tau_library(s)|.
```

The best fit is the smallest grid separation attaining `min J_j`; choosing
the smaller value makes ties deterministic. Its fit-uncertainty interval is
the contiguous grid component containing the best fit for which

```text
J_j(s) <= min J_j + 0.010.
```

The `0.010 tau` allowance is half the base source ordinate uncertainty and is
fixed before fitting. The grid contributes an additional `+/-0.0025` only
when reporting interval endpoints; it does not change membership.

## Accuracy at a tested separation

At each label separation, best-fit separation, and transformation-predicted
separation, apply the unchanged `CRITERIA-THEORY.md` tests:

- at least 50 of the 61 anchors `Fn = 0.20, 0.21, ..., 0.80` must be covered;
- median absolute error `<= 0.05`;
- 90th-percentile absolute error `<= 0.10`;
- principal-hump position error `<= 0.010 Fn` on
  `0.35 <= Fn <= 0.55`; and
- library/source principal-peak ratio in `[0.90, 1.10]`.

A panel meets the agreement threshold only when every component passes.
Interpolation, source uncertainty, and the fit objective do not relax these
limits.

## Transformations fixed before fitting

The probe tests exactly three shared transformations of the printed label:

```text
gap confusion:          S_eff/L = S_label/L - 0.10
half-spacing confusion: S_eff/L = (S_label/L) / 2
null definition:        S_eff/L = S_label/L
```

Here `B/L = 0.10` for C2. A transformation is consistent with a panel only if
its predicted separation lies within the reported fit-uncertainty interval
and its prediction passes every agreement component above. Values between
the `0.005` grid nodes may be linearly interpolated in separation solely for
evaluating the three preregistered transformations.

## Outcome

The separation-definition hypothesis is **CONFIRMED** if one transformation
is consistent with all four figures. It is **PARTIAL** if one transformation
is consistent with exactly three figures; the fourth panel and every failed
component are reported. It is **REFUTED** if no transformation is consistent
with at least three figures, if the best-fit intervals demand mutually
incompatible mappings, or if no shared transformation meets the agreement
thresholds.

The study reports all best fits, intervals, and component scores. It does not
fit a free affine transformation or invent another mapping after viewing the
results. Any descriptive post-hoc pattern is labelled as such and cannot
change the registered outcome.

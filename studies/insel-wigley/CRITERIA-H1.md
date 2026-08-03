# Pre-registered line-attribution hypothesis criteria

Committed after the independent legend audit and before any full-family
digitization is overlaid on a prediction or scored. These rules supplement but
do not alter `CRITERIA-THEORY.md`; the original C2 source, scores, and plots
remain frozen.

## Inputs and blind digitization

H1 concerns the four plotted model curves in Insel Figures 359--362, the
theoretical wave-resistance interference ratios at `S/L = 0.2, 0.3, 0.4,
0.5` (printed 358--359, PDF 368--369). Every newly digitized pass records both:

- the visible style followed: solid, dashed, dash-dot, or dotted; and
- the curve's relative vertical and phase position at source-read checkpoints.

Two source-only passes use independent calibration and path descriptions. They
must not inspect a library curve, score, comparison overlay, or the other pass.
Reconciliation uses the existing `DIGITIZATION.md` admission bounds: no more
than `0.015` disagreement in `Fn` and `0.08` in `tau`; unresolved crossings or
style loss are omitted rather than inferred. The first full-family figure is
Figure 359 because its amplitude separation is visually the greatest. One
source-data commit is made per figure.

## Curve-to-library match

Each candidate source curve is independently compared with the frozen library
C2 prediction under the unchanged `CRITERIA-THEORY.md` definitions:

- common scoring points are `Fn = 0.20, 0.21, ..., 0.80`, with at least 50 of
  61 points and complete principal-hump coverage required;
- pointwise agreement requires median absolute error `<= 0.05` and 90th-
  percentile absolute error `<= 0.10`;
- principal-hump position agreement on `0.35 <= Fn <= 0.55` requires
  `|delta Fn_hump| <= 0.010`; and
- principal-hump amplitude agreement requires a library/source peak ratio in
  `[0.90, 1.10]`.

A library match requires agreement on all three components. No horizontal
shift, amplitude scale, fitted phase, smoothing substitution, or empirical
correction is allowed.

## H1 decision rule

**H1 CONFIRMED** requires all of the following:

1. the library C2 prediction matches one specific plotted family member at
   every separation under the criteria above;
2. the independent Step 1 legend audit identifies that matched member as C2;
3. the originally digitized line is identified by the same legends as a
   different model; and
4. that different model's `L/B` position is consistent with both the observed
   interference-amplitude excess and feature-position offset.

**H1 REFUTED** if either:

- Step 1 unambiguously confirms the original solid-line assignment as C2; or
- no one plotted family member matches the library across all four separations.

The legend condition is decisive: a numerically closer non-C2 curve cannot be
relabeled as C2. All candidate scores remain reported even when the legend has
already refuted H1.

## Family-consistency checks

For each model, define its principal constructive amplitude as
`A = max(tau) - 1` and its hump location as the maximizing `Fn` on the same
`0.35 <= Fn <= 0.55` window and `0.001` grid used by
`CRITERIA-THEORY.md`. Propagate the archived digitization uncertainty at the
maximum; two values whose uncertainty intervals overlap count as tied rather
than reversed.

Insel states that higher `L/B` gives smaller wave-interference amplitude,
especially at `S/L = 0.2` and `0.3`, and that smaller `L/B` moves humps and
hollows to higher Froude number (printed 124, PDF 134). Applying the requested
four-curve check gives the registered amplitude order

```text
C3 (L/B 7) >= C4 (L/B 9) >= C2 (L/B 10) >= C5 (L/B 11).
```

At `S/L = 0.5`, Insel specifically describes the `L/B = 9` and 11 amplitudes
as practically the same; their uncertainty intervals must overlap or their
central amplitudes must differ by at most `0.05`.

The separately registered hull-form check requires the Wigley C2 principal
hump to have the lowest `Fn` of the four family members at each separation;
uncertainty overlap counts as a tie. Amplitude and hump orderings are reported
per figure before any causal interpretation.

Any resolved ordering reversal is a **source anomaly**. Per the requested stop
condition, the study records the anomaly and halts for review rather than
discarding a curve, changing the window, or interpreting around it.

## Branch after H1

If H1 is confirmed and no family anomaly occurs, the original theory scores are
preserved but marked superseded, the true C2 curves are rescored unchanged, and
the methods lesson is added to `DIGITIZATION.md` before closeout.

If H1 is refuted and no family anomaly occurs, the harness-only H2 archaeology
is limited to the four registered single deviations: literal equation (4.29),
doubled interference cross term, missing or doubled off-axis mode multiplicity,
and cosine versus cosine-squared pair coupling. A candidate must fit the printed
C2 curve family across all four separations simultaneously; a one-spacing fit
is not explanatory. Any fit is reported only as a plausible, unverifiable
reconstruction of a historical implementation deviation.

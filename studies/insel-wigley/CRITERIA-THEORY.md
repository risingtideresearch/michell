# Pre-registered Insel theory-to-theory criteria

Committed after the Chapter 4 source audit and before any C2 theory curve is
digitized, before any theory overlay is drawn, and before the comparison is
scored. These rules supplement but do not alter `CRITERIA.md` or any existing
experimental result.

## Scope and expectation

The source targets are the solid C2 curves in Insel Figures 359--362: the
theoretical catamaran-to-monohull wave-resistance interference ratio at
`S/L = 0.2, 0.3, 0.4, 0.5` (printed 358--359, PDF 368--369). Figures 355--358
are excluded because they reconstruct interference from measured monohull wave
patterns rather than from Insel's source-theory calculation.

Both calculations are linear thin-ship theories on the same fixed C2 Wigley
geometry and use the same centreline separations. They should agree much more
closely with each other than either agrees with experiment. The known
formulation differences are fixed in advance:

- Insel uses a finite-width, finite-depth canal Green function and a discrete
  transverse-mode sum; the library uses unbounded, infinite-depth Michell
  theory.
- Insel replaces rectangular centreplane elements by point sources. The thesis
  does not give the C2 element counts. The library integrates its B-spline
  representation analytically inside the outer integral.
- Fixed C2 uses no running-trim or sinkage correction, and its pointed Wigley
  stern has no transom contribution.

Insel estimated tank-wall effects below 1% and shallow-water effects below 4%
for the catamaran away from the depth-critical band (printed 63, PDF 73).
Allowing for those formulation differences plus reading a thick scanned line,
the criteria below are intentionally looser than solver tolerance but far
tighter than the experimental criteria.

## Blind digitization and deterministic curve

Two source-only passes will independently trace the C2 solid line on a 360 dpi
render. Neither pass may inspect an overlay, a score, or the other pass. Each
pass records raw calibration endpoints and raw pixel coordinates. Both use
affine axis calibration from the inner plot border: `0.1 <= Fn <= 1.0` and
`0 <= tau <= 2.5`.

The passes sample the visible C2 curve at common `Fn = 0.005` anchors wherever
the line is unambiguous. At an anchor present in both passes, the reconciled
ordinate is their mean. A point is admitted only if the pass difference is at
most `0.015` in `Fn` and `0.08` in `tau`; otherwise the source is re-inspected
without reference to predictions. If the C2 line cannot be separated from
another line, the point is omitted rather than inferred. Base digitization
uncertainty is `0.003` in `Fn` and `0.02` in `tau`; the archived uncertainty is
the larger of the base value and half the pass difference.

The reconciled source and library prediction are linearly interpolated without
extrapolation onto `Fn = 0.20, 0.21, ..., 0.80`. A configuration is scorable
only if at least 50 of those 61 points survive and the whole principal-hump
window below is covered. No horizontal shift, amplitude scale, fitted phase,
or empirical correction is permitted.

## 1. Pointwise agreement

At each common scoring point define

```text
e(q) = |tau_michell(q) - tau_insel(q)|.
```

The pointwise grade is:

| grade | median `e` | 90th percentile `e` |
|---|---:|---:|
| agreement | <= 0.05 | <= 0.10 |
| partial | <= 0.10 | <= 0.20 |
| disagreement | otherwise | otherwise |

Absolute error is used because unity has the physical meaning of no
interference and avoids unstable percentages near Insel's deep hollows. RMS
error, maximum error, and median signed error are also reported as diagnostics
but do not replace the pre-registered grade.

## 2. Principal constructive hump

For each separation, the principal constructive hump is the global maximum on
`0.35 <= Fn <= 0.55`, evaluated from each linearly interpolated curve on a
`0.001` grid. This window contains the broad design-speed hump in all four
source panels while excluding the rapid low-speed oscillations.

Hump position is graded:

| grade | `|delta Fn_hump|` |
|---|---:|
| agreement | <= 0.010 |
| partial | <= 0.020 |
| disagreement | > 0.020 |

Hump amplitude uses the ratio

```text
A = tau_michell(Fn_hump,michell) / tau_insel(Fn_hump,insel).
```

and is graded:

| grade | amplitude ratio `A` |
|---|---:|
| agreement | 0.90 <= A <= 1.10 |
| partial | 0.80 <= A <= 1.20 |
| disagreement | otherwise |

The `0.010` position allowance is roughly three times the base horizontal
reading uncertainty. The 10% amplitude allowance exceeds the source's stated
ordinary tank-effect estimates and the line-reading uncertainty, while still
requiring substantially closer agreement than the experiment study.

## 3. Configuration and study verdicts

A separation agrees only when its pointwise, hump-position, and hump-amplitude
grades all agree. It disagrees when any component disagrees; otherwise it is
partial. The study agrees only if all four separations agree, materially
disagrees if any separation disagrees, and is partial otherwise.

Any material disagreement triggers the requested stop condition. Before
stopping, the study may check calibration, curve identity, archived pass
reconciliation, normalization, separation convention, and prediction
reproducibility. It may not widen a threshold, shift or scale a curve, invoke
the experimental phase correction, or continue to the favourable closeout.

## 4. Critical-Froude consistency check

This check uses only the library prediction and does not affect the overlay
grade. On a committed `0.005` prediction grid spanning `0.20 <= Fn <= 0.95`,
define the critical Froude number for each separation as the smallest grid
value `Fn_c` for which

```text
|tau_michell(Fn) - 1| < 0.05
```

at that point and every subsequent grid point through `Fn = 0.95`. If none
exists, report `not reached`; isolated later excursions cannot be discarded.
Linear interpolation or smoothing is not used for this threshold. The result
is compared, without a pass/fail grade, with Insel's stated progression from
about `Fn = 0.55` at `S/L = 0.5` to about `Fn = 0.8` at `S/L = 0.2` (printed
131, PDF 141). The intermediate separations are reported but no undocumented
source values are invented for them.

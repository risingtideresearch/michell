# Pre-registered finite-canal attribution criteria

Committed after the independent canal implementation and gates G1/G2, but
before generating or inspecting any prediction at the physical tank dimensions
`W = 3.7 m`, `H = 1.85 m`. These rules supplement, and do not revise,
`CRITERIA.md` or `CRITERIA-THEORY.md`.

## Frozen inputs and comparison grid

The three theory curves are:

- `tau_unbounded`: the already committed library result;
- `tau_canal`: the independent continuous-Wigley finite-canal modal reference;
- `tau_Insel`: the already committed two-pass digitization of the solid C2
  curves in Insel Figures 359--362 (printed 358--359, PDF 368--369).

All use the fixed C2 geometry, `S/L = 0.2, 0.3, 0.4, 0.5`, and
`tau = R_pair/(2 R_mono)`. Canal values will be generated at the committed
`0.005` Froude-number anchors over `0.20 <= Fn <= 0.95`. Scores use every
source anchor in that interval that survived the blind digitization; prediction
curves are linearly interpolated without extrapolation if a source ordinate is
not exactly on a machine-identical anchor. No horizontal shift, smoothing,
amplitude scale, fitted phase, or empirical correction is permitted.

## Principal hump and gap closure

For each separation, the principal hump location `q_h` is the global maximum
of the linearly interpolated Insel curve on `0.35 <= Fn <= 0.55`, evaluated on
a `0.001` grid. This source-defined location is fixed before consulting either
prediction and prevents each model from selecting a different favourable
ordinate. Define the three hump ordinates at that common location as
`tau_X^h = tau_X(q_h)`.

The gap-closure fraction is

```text
f = (tau_canal^h - tau_unbounded^h)
    / (tau_Insel^h - tau_unbounded^h).
```

`f = 0` means the canal calculation leaves the unbounded discrepancy
unchanged; `f = 1` means it closes the signed discrepancy at the source hump.
The value is reported without clipping, so overshoot and motion in the wrong
direction remain visible. If the denominator has magnitude below `0.02`, the
fraction is declared unscorable and the study stops for review rather than
substituting another location or threshold.

The hump-amplitude relative error is evaluated at the same source-defined
location:

```text
h = |tau_canal^h - tau_Insel^h| / |tau_Insel^h|.
```

Each curve's own maximum and its location in the same window are reported as
diagnostics, but they do not replace `f` or `h`.

## Whole-curve residual

For every admitted source anchor `q` over `0.20 <= Fn <= 0.95`, define

```text
e_canal(q) = |tau_canal(q) - tau_Insel(q)|.
```

The median, 90th percentile, RMS, maximum, and median signed residual are
reported per separation. Only the median absolute residual enters the outcome
rule; the other statistics remain visible to expose localized disagreement.

## Outcome rule and stop condition

Outcomes are evaluated in this order:

1. **REPRODUCED**: at every `S/L`, median `e_canal <= 0.02` and `h <= 0.05`.
   The disagreement is fully attributed to water geometry, and both
   implementations are validated within the registered comparison resolution.
2. **UNEXPLAINED**: if REPRODUCED fails and `f < 0.5` at either `S/L = 0.2` or
   `0.3`. The canal hypothesis fails at close spacing. Append the result to
   `FINDINGS.md` and stop without `CLOSEOUT.md`; this outcome requires review
   because it reopens a possible implementation discrepancy.
3. **PARTIAL**: if REPRODUCED fails but `f >= 0.5` at both `S/L = 0.2` and
   `0.3`. The canal is the dominant identified cause. Any residual may be
   bounded against the observed score and associated with the known difference
   between the continuous analytic Wigley amplitude and Insel's undocumented
   point-source mesh and with registered curve-digitization uncertainty. It
   will not be assigned uniquely to either source without evidence.

No threshold will be widened after the tank-size curves are inspected.

## Experimental overlay and critical Froude diagnostics

Where fixed-attitude experimental `tau_WP` is available, it is plotted as a
fourth curve under the existing smoothing and labelling rules in
`CRITERIA.md`. It is not used to grade the theory-to-theory attribution and
does not revise the committed design trust envelope.

For each separation, the canal critical Froude number is the smallest committed
`0.005` grid value for which `|tau_canal - 1| < 0.05` at that point and every
subsequent point through `Fn = 0.95`; if none exists, report `not reached`.
The unbounded value is recomputed by the identical rule from
`CRITERIA-THEORY.md`. These are diagnostics against Insel's stated trend, not
additional pass/fail gates.

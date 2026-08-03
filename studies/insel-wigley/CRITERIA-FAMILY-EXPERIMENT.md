# Preregistered cross-model family-expectation test

This rule is committed before opening or rendering Insel Figures 347--350.
It tests whether the Figure 359 family-order stop reflects a source anomaly or
an invalid cross-hull extrapolation of an NPL-family trend. It does not alter
`CRITERIA-H1.md` or relabel any theoretical curve.

## Inputs and attribution

The evidence consists of the experimental cross-model interference curves in
Figures 347--350 (printed 352--353, PDF 362--363) and the accompanying Chapter
8 discussion. Before tracing, the audit records each caption and the mapping
from printed line style to C2--C5. Unreadable styles, crossings, or axes remain
unassigned.

Two source-only passes use independent calibrations and path descriptions.
Each row records the model, printed line style, and relative curve position.
The existing reconciliation limits apply: the passes may differ by at most
`0.015` in `Fn` and `0.08` in `tau`. Unresolved anchors are omitted.

## Quantities

For each model and separation, define the principal constructive-interference
amplitude on `0.35 <= Fn <= 0.55` as

```text
A_model = max(tau_model) - 1.
```

Propagate the reconciled ordinate uncertainty at the maximum. Define the NPL
envelope and the measured C2 advantage as

```text
A_NPL = max(A_C3, A_C4, A_C5)
D_exp = A_C2 - A_NPL
Q_exp = A_C2 / A_NPL,
```

provided both amplitudes are positive. An excess is resolved only when the
lower uncertainty bound of `A_C2` exceeds the upper bound of `A_NPL`.

Figure 359 gives the already archived theoretical close-spacing values
`A_C2 = 0.979` and `A_NPL = 0.582`, hence `D_theory = 0.397` and
`Q_theory = 1.68`. A measured excess is **qualitatively comparable** when it
is resolved and either:

- `D_exp >= 0.20`, half the Figure 359 theoretical excess rounded outward; or
- `Q_exp >= 1.34`, halfway from unity to the Figure 359 theoretical ratio.

A measured excess is **material** when it is resolved and either
`D_exp >= 0.10` or `Q_exp >= 1.15`. These lower bounds separate a real
cross-form effect from digitization-scale overlap without requiring the tank
curve to reproduce the theoretical magnitude.

## Decision rule

The family-order anomaly stop **DISSOLVES** if, at `S/L = 0.2` and `0.3`, C2
is never resolved below the NPL envelope, at least one separation has a
qualitatively comparable C2 excess, and the other has at least a material C2
excess. The theoretical ordering then agrees with Insel's own measured
cross-model ordering, even if its magnitude differs.

Figure 359 is **ANOMALOUS-VS-OWN-EXPERIMENT** if neither close separation has
a material measured C2 excess, or if an NPL model is resolved above C2 at
either close separation. In that outcome, the theoretical C2 separation from
the NPL family lacks support in the paired experimental panels.

Any remaining combination is **MIXED**: the cross-hull extrapolation remains
unsupported, but the experiment does not establish a contradiction. The
study reports the registered outcome and proceeds to the separation probe in
all cases.

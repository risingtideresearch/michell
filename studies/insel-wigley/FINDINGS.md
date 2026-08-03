# Insel-Wigley external validation findings

Criteria were frozen in commit `3f81c2a` before the first prediction. Raw
scores are under `data/analysis/`; raw numerical results are
`data/predictions/predictions.csv`. Claim classes used below are:

- **K** — known result reproduced;
- **E** — engineering result for this codebase;
- **P** — matches published practice;
- **N?** — possibly novel.

No claim in this study is class N?. This is the first external experimental
comparison archived for this codebase, but that is an engineering milestone,
not a scientific breakthrough.

## Verdict

**E — The frozen kernel passes the numerical contract but does not earn a
continuous design-grade trust envelope against these experiments.** All 578
prediction rows, including every pair and standalone solve, converged at the
requested `rel_tol = 1e-6`. Nevertheless, only seven isolated preregistered
trust cells pass and no `S/L` has three consecutive passing `Fn` points.

**E — Absolute wave-pattern coefficient agreement is the limiting result.**
The fixed monohull is partial; catamarans are disagreement at `S/L = 0.2` and
0.3 and partial at 0.4 and 0.5. Michell gets the principal monohull hump
position and amplitude inside the agreement thresholds, but generally predicts
the catamaran hump too early and too high.

**E — The interference ratio is substantially more useful than the absolute
coefficient at the wider separations.** The preregistered verdict progresses
from disagreement at `S/L = 0.2`, to partial at 0.3 and 0.4, to agreement at
0.5. That result supports interference-guided comparison at `S/L = 0.5` for
this Wigley geometry; it does not by itself make the absolute resistance
design-grade.

## Numerical and normalization checks

**K — The exact C2 product-parabola geometry was reproduced.** The library's
single-span biquadratic constructor represents `L = 1.8 m`, `B = 0.18 m`, and
`T = 0.1125 m` exactly. Catamaran centreplanes are at `+/- S/2`, with
`S/L = 0.2, 0.3, 0.4, 0.5`.

**E — The study is numerically resolved relative to the experimental
discrepancy.** All rows use the general marcher and report `Converged`. The
largest pair `est_rel_error` is `9.75e-7`; the largest standalone estimate is
`9.20e-8`. The largest pair solve used 375,200 inner evaluations. These errors
are many orders below the measured/model coefficient differences and cannot
explain the verdict.

**K — All coefficients are at tank-model scale.** Calculations use
`rho = 1000 kg/m^3`, `nu = 1.141e-6 m^2/s`, and `g = 9.80665 m/s^2`, with
`S_W = 0.482 m^2` for the monohull and `0.964 m^2` for a pair. No full-scale
conversion appears anywhere in the harness or outputs.

## Fixed-attitude coefficient score

The pointwise statistic is the normalized absolute discrepancy `d` in
`CRITERIA.md`. Hump amplitude error is relative to the experimental hump;
negative `delta Fn` means the prediction peaks early.

| configuration | median `d` | 80th `d` | `delta Fn` hump | hump amplitude error | verdict |
|---|---:|---:|---:|---:|---|
| monohull | 0.315 | 0.643 | +0.002 | +23.9% | partial |
| `S/L = 0.2` | 0.643 | 1.196 | -0.033 | +60.9% | disagreement |
| `S/L = 0.3` | 0.368 | 1.083 | -0.033 | +50.8% | disagreement |
| `S/L = 0.4` | 0.315 | 0.592 | -0.024 | +30.3% | partial |
| `S/L = 0.5` | 0.375 | 0.667 | -0.016 | +29.1% | partial |

**E — The main absolute-coefficient failure is an amplitude bias, not numerical
noise.** The predicted catamaran principal humps are 29--61% above measured
`C_WP`. The bias is largest for the closest spacings and decreases as the
demihulls separate. The plots also show the expected high prediction relative
to `C_WP` over much of the post-hump range.

**K — A prediction above `C_WP` is physically plausible but is not automatically
agreement.** `C_WP` can miss wave energy through probe-field attenuation,
finite longitudinal record, finite harmonic/matrix resolution, and weak-wave
resolution. Insel did not use a transverse-cut integral, so transverse-cut
truncation is not offered as an explanation. The preregistered scores give no
credit for the expected direction of this systematic difference.

## Fixed-attitude wave-pattern interference

The experimental observable is explicitly `tau_WP = C_WP,cat/C_WP,mono`, not
Insel's resistance-decomposition `tau`, which additionally requires the fitted
viscous multiplier `beta`.

| `S/L` | median abs. error | 80th abs. error | sign agreement | hump: `delta Fn`, abs. error | hollow: `delta Fn`, abs. error | verdict |
|---:|---:|---:|---:|---:|---:|---|
| 0.2 | 0.320 | 0.406 | 7/9 (77.8%) | -0.074, 0.321 | -0.035, 0.021 | disagreement |
| 0.3 | 0.241 | 0.431 | 9/10 (90.0%) | -0.031, 0.331 | -0.020, 0.167 | partial |
| 0.4 | 0.066 | 0.258 | 9/11 (81.8%) | -0.026, 0.024 | -0.018, 0.057 | partial |
| 0.5 | 0.035 | 0.147 | 6/6 (100%) | -0.014, 0.054 | -0.017, 0.087 | agreement |

**E — The closest-spacing problem is primarily phase.** At `S/L = 0.2` the
dominant constructive hump is 0.074 in Froude number early, which fails even
the partial threshold. The preceding hollow amplitude is accurate but its
position remains only partial. At `S/L = 0.5`, level, sign, hump, and hollow all
meet the agreement criteria.

**E — Separation trends are captured incompletely.** The median pointwise
Spearman correlation for ranking the four separations is 0.40 (partial). The
rank correlation of the RMS interference envelope over `0.40 <= Fn <= 0.55`
is 0.80 (agreement). Thus the code better captures how overall interference
strength decays with spacing than the instantaneous ordering of phase-shifted
curves.

## Trust envelope

The passing fixed-attitude cells are:

| `S/L` | passing `Fn` cells |
|---:|---|
| 0.2 | 0.325, 0.350 |
| 0.3 | 0.350 |
| 0.4 | 0.325, 0.500 |
| 0.5 | 0.350, 0.500 |

**E — No continuous interval is design-grade under the preregistered rule.**
Every list above has fewer than three consecutive `0.025` grid points. The
defensible narrow statement is that `S/L = 0.5` interference shape and level
agree on this hull over the scored range; absolute `C_W` and the joint
trust-cell rule remain weaker.

This is a validated negative result, not a kernel defect. The numerical solver
converges honestly; the remaining gap is between linear fixed-geometry theory
and measured wave-pattern resistance plus its experimental transfer function.

## Free-attitude secondary comparison

The conservative blind pass leaves no valid smoothed free-attitude point at or
below `Fn = 0.4` for four configurations; `S/L = 0.3` has only one. Therefore
the preregistered low/high running-attitude contrast is not scoreable without
relaxing the digitization rule.

Above `Fn = 0.4`, the descriptive median `d` values are 0.031 (monohull), 0.835
(`S/L = 0.2`), 0.070 (0.3), 0.098 (0.4), and 0.166 (0.5).

**E — The preregistered expectation of universally increasing free-attitude
divergence above `Fn ≈ 0.4` is not supported.** It is strongly visible at
`S/L = 0.2` but not in the other archived CWP series. This result is descriptive
and cannot enter the fixed-attitude trust envelope: the calculation did not
apply measured trim/sinkage or re-immerse the hull.

## Published precedent

**P — The model-scale fixed/free comparison policy follows published
practice.** Ship Science Report 72 inserted measured trim and sinkage into its
theoretical NPL hulls, regenerated the source panels, and treated the remaining
linear-theory discrepancy qualitatively. This study is stricter for its primary
score by using fixed C2-FX before looking at free C2.

There is no same-geometry Report 72 number to place beside ours. Report 72
compares different round-bilge, transom-stern NPL forms from Report 71, not
Insel's C2 Wigley hull. Its qualitative conclusion—useful slender-hull trends,
with larger discrepancies outside the favourable regime—is context only.
Claiming a numerical C2 match to Report 72 would be a source error.

## What remains untested or blocked

- **Resistance-decomposition `tau`:** blocked by the absence of raw `beta`
  values paired with the marker data. `tau_WP` is not silently substituted.
- **Static wetted-area convention:** `0.482/0.964 m^2` is the documented best
  reading but the thesis does not state the catamaran factor in one explicit
  sentence. Raw dimensional resistance and area remain in every output row.
- **Hull-form generality:** only one exact Wigley form is tested. The NPL series
  is the required next study before generalizing the envelope.
- **Viscous interference:** absent from Michell and not tested by `tau_WP`.
- **Running attitude:** free data are secondary and do not test a re-immersed
  hull calculation.
- **Experimental transfer:** the relation between true radiated wave resistance
  and measured `C_WP` is not independently calibrated.
- **Tank effects:** the high-speed end retains Insel's acceleration,
  shallow-water, wave-breaking, and weak-signal cautions.

## Claim classification summary

| claim | class | status |
|---|---|---|
| exact C2 geometry, placement, and model-scale normalization reproduced | K | validated |
| every study solve converges at the requested tolerance | E | validated |
| monohull coefficient is partial; close catamarans disagree | E | validated negative |
| `S/L = 0.5` `tau_WP` interference meets all criteria | E | validated |
| a continuous design-grade envelope exists | E | rejected |
| procedure matches Report 72's model-scale published practice | P | contextual match |
| numerical C2 result matches Report 72 | P | not claimable; source mismatch |
| any result here is possibly novel | N? | no claim made |

## Reproduction

```sh
python3 studies/insel-wigley/data/digitized/validate_digitization.py
cargo run --release --manifest-path studies/insel-wigley/harness/Cargo.toml
uv run --project python --extra plot python \
  studies/insel-wigley/plots/analyze.py
```

The last command regenerates all score CSVs and 14 comparison plots from the
committed digitized and prediction CSVs. It reads no PDF and no solver state.

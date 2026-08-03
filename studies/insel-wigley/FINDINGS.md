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

## Interpretation against the primary source

The preregistered scores above remain unchanged, but their physical
interpretation closes differently when read against Insel's own conclusions.
Insel states that both theoretical calculations and predictions reconstructed
from the monohull wave pattern can carry a Froude-number phase shift. His
recommended procedure is empirical: locate the experimental principal hump in
Figures 375--377 and shift the theoretical wave-resistance curve to it
(printed 125, PDF 135). The measured theoretical lead in this study,
`delta Fn = -0.016` to `-0.074`, therefore quantifies a known `phase shift` on
the exact C2 comparison rather than exposing an undocumented solver behavior.

The separation dependence is equally explicit in the primary source. Insel
warns that predictions using theoretical `tau` at lower speeds, especially at
small separation, should be `treated with caution` (printed 128, PDF 138).
The disagreement at `S/L = 0.2` is thus inside the regime he singles out. His
mechanistic interpretation is that neighbour-induced asymmetric flow changes
the wave phase but is absent from the symmetric thin-ship approximation; he
says that effect must be `corrected by empirical methods` (printed 130,
PDF 140). He also reports inter-hull bow/stern wave breaking and says it is
likely to make wave-pattern analysis underpredict wave resistance over some
speed ranges (printed 130, PDF 140). That measurement bias acts in the same
direction as part of the present 29--61% theory-over-`C_WP` hump-amplitude gap
and should be strongest where inter-hull interaction is strongest. The
available record does not calibrate the bias, so this study does not assign a
numerical fraction of the gap to it.

Insel's preferred design observable also matches the score hierarchy here: he
concludes that using interference factors instead of direct wave-resistance
prediction gives `much better accuracy` (printed 131, PDF 141). This is the
published counterpart of the present progression from weak absolute-`C_W`
scores to agreement for `tau_WP` at `S/L = 0.5`. The same pattern predates the
thesis. Insel's review reports that Yokoo and Tasaki observed a Froude-number
phase shift between calculated and measured interaction humps and hollows
(printed 17, PDF 27), while Eggers obtained satisfactory theory/experiment
agreement only for `S/L > 0.4` (printed 16--17, PDF 26--27). The improvement
with separation in this study therefore reproduces the historical shape of
linear-theory validation rather than creating a new failure mode.

There is also an implementation-level anchor independent of the experimental
wave-pattern transfer. At `Fn = 0.35`, the committed monohull row in
`predictions.csv` gives `C_W = 1.248133344362e-3`; Doctors and Beck's published
classical thin-ship value for the same `B/L = 0.1`, `T/L = 0.0625` Wigley
geometry is `1.2486e-3`. The relative difference is 0.037%, or 0.04% rounded.
That agreement, together with converged study rows, is direct evidence that
the implementation reproduces the classical C2 calculation.

The classification consequently separates the phenomenon from its
measurement. The close-spacing absolute-`C_W` limitation and hump phase shift
are class K/P: known behavior reproduced and published practice matched. The
class-E contribution is the preregistered numerical score and explicit trust
envelope, including its rejection of any continuous interval. No result in
this study indicts the frozen numerical kernel.

### Claim-classification addendum

| claim | class | status |
|---|---|---|
| Froude-number phase shift in linear catamaran theory | K | known result reproduced and quantified on C2 |
| degraded close-spacing agreement and improvement beyond `S/L = 0.4` | K/P | historical result reproduced; published caution matched |
| interference factors outperform direct wave-resistance prediction | P | Insel's stated design practice matched |
| Doctors--Beck C2 thin-ship value at `Fn = 0.35` | K | reproduced within 0.04% |
| preregistered cell-by-cell trust envelope | E | quantitative codebase contribution; no continuous interval admitted |
| frozen numerical kernel is defective | E | rejected by the evidence in this study |

## Preregistered theory-to-theory comparison

This follow-up adds a stronger test than the experiment overlay: the frozen
library prediction is compared directly with Insel's calculated C2
interference curves in Figures 359--362 (printed pages 358--359, PDF pages
368--369). `CRITERIA-THEORY.md` was committed before either source curve was
digitized or any prediction overlay was produced. The source archive retains
two source-only readings, rejects 36 over-tolerance or line-ambiguous anchors,
and leaves 54--59 of the 61 exact scoring anchors per panel.

The preregistered study verdict is **material disagreement**:

| `S/L` | median `|delta tau|` | 90th percentile | pointwise | `|delta Fn_hump|` | hump position | peak ratio, library/Insel | hump amplitude | configuration |
|---|---:|---:|---|---:|---|---:|---|---|
| 0.2 | 0.0820 | 0.2993 | disagreement | 0.005 | agreement | 0.844 | partial | disagreement |
| 0.3 | 0.0413 | 0.1450 | partial | 0.005 | agreement | 0.917 | agreement | partial |
| 0.4 | 0.0369 | 0.1270 | partial | 0.010 | agreement | 0.933 | agreement | partial |
| 0.5 | 0.0203 | 0.0714 | agreement | 0.010 | agreement | 0.964 | agreement | agreement |

The result is not a repeat of the experimental phase-shift finding. All four
principal-hump positions agree within `0.010` in `Fn`; the discrepancy is
predominantly an amplitude excess in Insel's finite-canal calculation that
grows as separation decreases. At `S/L = 0.2`, the unambiguous broad source
peak is `tau = 1.980` at `Fn = 0.450`, versus `1.672` at `Fn = 0.445` from the
library. That `0.309` peak difference is much larger than the archived source
reading uncertainty. The pointwise 90th-percentile failure therefore does not
depend on the visibly difficult high-speed line crossings in the scan.

### Post-result audit permitted by the preregistration

- The plotted source is the solid C2 curve, not the dashed C3--C5 curves or
  Figures 355--358's wave-pattern reconstruction. Plot borders and affine axes
  were checked against the rendered source; both passes and both calibrations
  remain archived.
- `S` is demihull-centreline spacing, and the placements remain `+/- S/2`.
  The ratio is `R_w,pair / (2 R_w,solo)`, equivalent to Insel's catamaran-
  coefficient/monohull-coefficient normalization when the catamaran reference
  area is twice the demihull area. Unity therefore has the same noninteracting
  meaning in both calculations.
- Every one of the 604 new pair solves and 151 cached standalone solves
  converged at the frozen settings. Regenerating the original 578-row harness
  output leaves its SHA-256 unchanged at
  `5addfb836c80e5385de6b7b92819addfc529d22c140a433d50a6d2221d9fefff`.
- The Doctors--Beck monohull anchor still agrees within 0.04%. The discrepancy
  is separation-dependent, rather than a uniform monohull scale error.

The known formulation difference is a serious candidate, not a resolution:
Insel used a finite-width, finite-depth modal canal Green function and an
undocumented point-source mesh, while the library uses exact hull moments with
the unbounded deep-water Michell kernel. The widening disagreement toward
small `S/L` is qualitatively compatible with a transverse-boundary effect, but
this study has not isolated wall, depth, or mesh contributions. It would be
incorrect to use the source's ordinary tank-effect estimates to explain away a
16% principal-peak difference without reproducing the finite-canal
calculation.

Accordingly, the earlier implementation evidence remains valid but the last
sentence of the interpretation addendum is no longer sufficient to close the
multihull kernel question. This direct theory comparison is a genuine red flag
for close-spacing equivalence, not proof of a numerical bug in either solver.
Per the preregistered stop rule, no favourable `CLOSEOUT.md` is written. The
next review must reproduce Insel's finite-canal calculation or otherwise
separate Green-function physics from numerical implementation before the
experimental result is promoted as kernel closure.

### Theory-comparison claim classification

| claim | class | status |
|---|---|---|
| unbounded and finite-canal C2 hump positions agree within `0.010 Fn` | E | validated at all four separations |
| theory-to-theory amplitude and pointwise equivalence | E | rejected; material disagreement at `S/L = 0.2` |
| wide-spacing (`S/L = 0.5`) theory equivalence | E | validated under all preregistered components |
| finite-canal physics explains the close-spacing gap | E | suspected, not isolated or verified |
| frozen numerical kernel is defective | E | unresolved by the direct multihull comparison; not established |

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

At this checkpoint the known formulation difference was a serious candidate,
not a resolution:
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
Per the preregistered stop rule, no favourable `CLOSEOUT.md` was written at
this checkpoint. The finite-canal follow-up below performs the required next
test.

### Theory-comparison claim classification

| claim | class | status |
|---|---|---|
| unbounded and finite-canal C2 hump positions agree within `0.010 Fn` | E | validated at all four separations |
| theory-to-theory amplitude and pointwise equivalence | E | rejected; material disagreement at `S/L = 0.2` |
| wide-spacing (`S/L = 0.5`) theory equivalence | E | validated under all preregistered components |
| finite-canal physics explains the close-spacing gap | E | suspected at this checkpoint; tested below |
| frozen numerical kernel is defective | E | unresolved by the direct multihull comparison; not established |

## Critical-Froude consistency check

The preregistered diagnostic defines `Fn_c` as the first point on the committed
`0.005` grid after which `|tau - 1| < 0.05` remains true through `Fn = 0.95`.
It is intentionally not interpolated, smoothed, or graded. Insel describes the
separation-dependent progression as approximately `Fn = 0.55` at `S/L = 0.5`
to `Fn = 0.8` at `S/L = 0.2` (printed 131, PDF 141); he does not give numerical
intermediate thresholds there.

| `S/L` | library `Fn_c` | Insel statement |
|---|---:|---:|
| 0.2 | not reached by 0.95 | about 0.8 |
| 0.3 | not reached by 0.95 | not stated |
| 0.4 | 0.585 | not stated |
| 0.5 | 0.550 | about 0.55 |

The wide-spacing endpoint reproduces Insel's statement exactly on this grid.
The close-spacing endpoint does not: although the library curve approaches the
5% band near `Fn = 0.8`, it does not stay there, ending at `tau = 0.9189` at
`Fn = 0.95`. The `S/L = 0.3` curve likewise ends just outside at `tau =
0.9434`. This mixed result is consistent with, but does not identify the cause
of, the separation-dependent theory-to-theory discrepancy. It strengthened the
reason to stop before closeout and to isolate the finite-canal Green function
in the follow-up below.

## Finite-canal attribution follow-up

`CRITERIA-CANAL.md` was committed before any physical-tank result was computed
or inspected. The independent reference follows Insel's discrete transverse
modes, finite-depth dispersion and vertical profile, centered-catamaran
eigenfunction factors, and modal resistance sum (Chapter 4, printed 40--54,
PDF 50--64). It integrates the continuous C2 product parabola analytically;
Insel's source-panel counts remain undocumented.

The implementation exposed one source inconsistency before the tank run.
Equation (4.25), printed page 47 (PDF 57), contains the dimensional factor
`K_0 + K_n cos^2(theta_n)` in the wave-elevation coefficient, while equation
(4.29), printed page 48 (PDF 58), omits it. Literal omission is dimensionally
inconsistent and failed the preregistered wide/deep limit by the exact
asymptotic factor `4 K_0^2`. The reference retains the factor from the governing
equation (4.25); the rationale and algebra are recorded in `METHOD.md`. No
alternate tank curve using the failed transcription was computed.

Both pre-run gates pass. G1 checks `Fn = 0.25, 0.35, 0.50` at every separation
against an independent analytic-Wigley quadrature and the frozen library for
canals scaled by 5, 20, and 80, with the 80-times endpoint inside 0.5% for
monohull resistance and 0.02 for `tau`. G2 checks the physical tank dimensions
without exposing comparison values and shows that 128 versus 256 modes changes
both observables by less than 0.002. The production calculation uses successive
mode doubling with tighter `5e-6` relative-resistance and `2e-6` absolute-
interference criteria; all 755 output rows converged.

### Preregistered outcome: UNEXPLAINED

The physical canal reference remains nearly coincident with the unbounded
library curve at every principal hump. It therefore does not close the
published-curve amplitude gap:

| `S/L` | median `|tau_canal - tau_Insel|` | canal hump error | gap closure `f` |
|---|---:|---:|---:|
| 0.2 | 0.0835 | 15.53% | 0.004 |
| 0.3 | 0.0510 | 8.43% | 0.010 |
| 0.4 | 0.0425 | 7.31% | 0.019 |
| 0.5 | 0.0171 | 3.83% | 0.059 |

The registered close-spacing gate requires `f >= 0.5` at both `S/L = 0.2`
and 0.3. Values of 0.004 and 0.010 fail by two orders of magnitude. This is not
a marginal threshold decision: at the `S/L = 0.2` source hump (`Fn = 0.450`),
the source is `tau = 1.9801`, the unbounded reference is `1.6714`, and the canal
reference is `1.6725`.

The critical-Froude diagnostic changes only modestly:

| `S/L` | unbounded `Fn_c` | canal `Fn_c` | Insel statement |
|---|---:|---:|---:|
| 0.2 | not reached | not reached | about 0.8 |
| 0.3 | not reached | 0.640 | not stated |
| 0.4 | 0.585 | 0.595 | not stated |
| 0.5 | 0.550 | 0.565 | about 0.55 |

The finite canal affects the high-speed 5% settling rule at intermediate
spacing, but it does not produce the close-spacing principal-hump amplitude in
Insel's published curves. The original hypothesis is rejected under the
preregistered definition. The remaining discrepancy cannot be assigned to the
undocumented point-source mesh or digitization as a bounded residual because
the PARTIAL gate was not reached. It reopens the possibility of an
implementation difference in Insel's unpublished calculation, this
transcription, or the published curve identification; the present evidence
does not choose among them.

The experimental trust envelope is unchanged. Both the unbounded and canal
references still overpredict close-spacing interference against measured
`tau_WP`, so the earlier model-form conclusion remains in force. Per the
registered UNEXPLAINED stop condition, this study halts here and does not
create `CLOSEOUT.md`.

### Canal-follow-up claim classification

| claim | class | status |
|---|---|---|
| independent canal solver recovers the deep unbounded limit | E | validated by G1 against two references |
| physical-tank modal sum is truncated below comparison tolerance | E | validated by G2 and tighter production convergence |
| finite-canal water geometry explains the close-spacing amplitude gap | E | rejected; `f = 0.004` and `0.010` at close spacing |
| finite-canal and unbounded references are equivalent at the principal hump | E | matched to within 0.33% in `tau` at all four separations |
| either frozen kernel or Insel's published calculation is defective | E | unresolved; no side identified |

## H1 line-attribution follow-up: source-anomaly stop

The independent legend audit refutes H1 without numerical matching. In every
one of Figures 359--362, the legend assigns the solid line to `C2 (WIGLEY
HULL)`, the dashed line to C3, the dash-dot line to C4, and the dotted line to
C5 (Figures 359--360: printed 358, PDF 368; Figures 361--362: printed 359, PDF
369). This exactly matches the original two-pass attribution. Figures 355--358
are separately captioned as predictions from monohull wave-pattern analysis,
so no hybrid/pure panel swap was found (printed 356--357, PDF 366--367).

The preregistered family check then encountered its own stop condition in the
first and most discriminating panel. Two independent source-only traces of all
four Figure 359 curves admitted 503 of 510 common anchors. On the registered
`0.35 <= Fn <= 0.55` window, `A = max(tau) - 1` is:

| model | legend identity | peak `Fn` | `A` | ordinate uncertainty |
|---|---|---:|---:|---:|
| C2 | Wigley hull, solid | 0.450 | 0.979 | 0.020 |
| C3 | RBH `L/B = 7`, dashed | 0.525 | 0.582 | 0.020 |
| C4 | RBH `L/B = 9`, dash-dot | 0.525 | 0.582 | 0.020 |
| C5 | RBH `L/B = 11`, dotted | 0.500 | 0.552 | 0.020 |

The C3/C4 intervals overlap and count as tied. C2's hump occurs below the
earliest RBH hump and passes the registered feature-position ordering. But the
registered amplitude relation `C4 >= C2` is reversed by `0.397`; even the
nearest uncertainty bounds remain separated by `0.357`. This is not a line-
crossing or threshold-edge result. It conflicts with the thesis statement
that higher `L/B` gives smaller interference amplitude at `S/L = 0.2` and
0.3 (printed 124, PDF 134), when that statement is applied to the legend
identities exactly as preregistered.

This is recorded as a **source anomaly**, not as evidence that a plotted curve
should be relabelled. The machine-readable result is in
`data/digitized/h1_figure_359_family_checks.csv`; both raw passes, the seven
omitted anchors, reconciled family, and peak metrics remain archived beside
it. Per `CRITERIA-H1.md`, work stops here for review. No library overlay was
used to reinterpret the family, Figures 360--362 were not redigitized, no H2
implementation variants were run, prior theory scores and predictions remain
unchanged, and no `CLOSEOUT.md` is created.

### H1-follow-up claim classification

| claim | class | status |
|---|---|---|
| original solid-line attribution is wrong | E | rejected by unambiguous legends |
| Figure 359 satisfies the registered L/B amplitude ordering | E | rejected; source-anomaly stop triggered |
| Wigley C2 has the lowest broad-hump `Fn` in Figure 359 | K | reproduced within archived digitization uncertainty |
| a different plotted family member explains the library discrepancy | E | not tested; overlays prohibited after the source-anomaly stop |
| a specific historical implementation deviation explains the curves | E | not tested; H2 was not entered |

## Experimental re-adjudication of the family expectation

The earlier family rule overreached. Insel states that smaller `L/B` moves
humps and hollows to higher Froude number and that higher `L/B` reduces the
interference amplitude at `S/L = 0.2` and `0.3` (printed 124, PDF 134). His
three round-bilge NPL forms establish that trend within one hull-form family.
The text does not extend the amplitude ordering across the deep, parabolic,
transom-free Wigley C2 and the shallow, round-bilge, transom NPL forms. The
registered `C3 >= C4 >= C2 >= C5` order therefore included an untested
cross-form extrapolation.

`CRITERIA-FAMILY-EXPERIMENT.md` replaced that extrapolation with a direct
test before Figures 347--350 were opened. The independent legend audit maps
solid to C2, dashed to C3, dash-dot to C4, and dotted to C5 in all four panels;
the separate dash-dot-dot curve is C2-FIXED (Figures 347--348: printed 352,
PDF 362; Figures 349--350: printed 353, PDF 363). Two source-only passes of the
principal-hump window admitted 652 of 656 anchors.

The close-spacing experimental amplitudes are:

| `S/L` | C2 `A` | NPL envelope | NPL `A` | C2 minus NPL | result |
|---:|---:|---|---:|---:|---|
| 0.2 | 0.775 | C3 | 0.883 | -0.109 | NPL resolved above C2 |
| 0.3 | 0.560 | C3 | 0.697 | -0.137 | NPL resolved above C2 |

Each peak ordinate carries `0.020` digitization uncertainty. The intervals do
not overlap. Neither close separation has a material C2 excess; both place C3
above C2. The preregistered result is therefore
**ANOMALOUS-VS-OWN-EXPERIMENT**. Figure 359's theoretical C2 amplitude does
not merely violate an extrapolated NPL trend: it reverses the cross-model
ordering in Insel's paired experimental comparison. Chapter 8 says that
Figures 347--350 show pronounced `L/B` dependence at small separation and
that amplitude depends on both hull form and separation (printed 124, PDF
134); it does not explain this reversal.

### What centred differences do and do not exclude

Insel computes element source density from centred station differences of
corner offsets, then replaces each continuous element with a point source
(equation 4.53, printed 56, PDF 66). For a quadratic waterline
`q(s) = a s^2 + b s + c`, the longitudinal difference is exact:

```text
[q(s+h/2) - q(s-h/2)] / h = 2as + b = q'(s).
```

Because each Wigley waterline is quadratic in `x`, station spacing cannot bias
that longitudinal derivative. The stronger claim that this makes every C2
source or resistance exact at any mesh is false. Equation 4.53 also averages
vertical corner offsets, which is trapezoidal rather than exact for the
Wigley's quadratic depth factor, and it applies the continuous element at one
point, which is not exact for the oscillatory wave phase. The deduction
excludes longitudinal source-gradient resolution as the cause; it does not
exclude vertical averaging or point-source quadrature without a mesh study.

### Family-experiment claim classification

| claim | class | status |
|---|---|---|
| Insel claimed one amplitude order across Wigley and NPL forms | K | rejected; the text states an NPL `L/B` trend and separate hull-form dependence |
| Figure 359 agrees with the measured close-spacing cross-model ordering | E | rejected at both close separations |
| centred station differences recover the Wigley longitudinal derivative | K | exact by the quadratic identity above |
| every mesh-resolution mechanism is thereby excluded | E | rejected; vertical averaging and point-source phase quadrature remain |

## Separation-definition probe

`CRITERIA-SEPARATION.md` was committed before the harness evaluated any new
separation. The exact C2 solver then produced 16,245 converged pair rows on
`S/L = 0.080:0.005:0.550` and `Fn = 0.150:0.005:1.000`. One `Fn = 1.0` row hit
the original eight-refinement harness cap; the completed run raised only that
harness allowance to ten and retained `rel_tol = 1e-6`.

The registered median-error fits are:

| figure | label `S/L` | best `S_eff/L` | fit interval | best median error | best p90 | best-fit agreement |
|---:|---:|---:|---:|---:|---:|---|
| 359 | 0.200 | 0.195 | 0.1875--0.2025 | 0.0913 | 0.2952 | no |
| 360 | 0.300 | 0.300 | 0.2775--0.3325 | 0.0413 | 0.1431 | no |
| 361 | 0.400 | 0.345 | 0.3275--0.3725 | 0.0215 | 0.1512 | no |
| 362 | 0.500 | 0.460 | 0.4225--0.5500 | 0.0185 | 0.0967 | yes |

The fitted values do not follow one definition change. `S_label - B` and
`S_label/2` are consistent with zero of four panels. The label itself passes
the interval and full-agreement requirements only for Figure 362. The null
mapping therefore scores one of four; no mapping reaches the three-panel
PARTIAL threshold. The preregistered outcome is **REFUTED**.

The result also rejects a looser spacing explanation. Allowing each figure its
own free separation cannot bring Figures 359--361 inside every unchanged
agreement threshold: their 90th-percentile errors remain 0.295, 0.143, and
0.151. The close-spacing amplitude excess is therefore not a gap-versus-
centreline or half-spacing transcription hidden in the labels. Insel defines
`S` as demihull-centreline separation in the theoretical derivation (printed
52, PDF 62), consistent with the null interpretation tested here.

### Separation-probe claim classification

| claim | class | status |
|---|---|---|
| printed `S/L` means clear gap rather than centreline separation | E | rejected by the registered fit and primary definition |
| printed `S/L` is twice the implemented half-spacing | E | rejected; zero panels fit the half-spacing mapping |
| one shared separation transformation explains all C2 curves | E | rejected; best fits are mutually inconsistent |
| free per-panel separation restores theory agreement | E | rejected for Figures 359--361; Figure 362 agrees |

## Time-boxed historical-variant probe

Because the separation hypothesis was refuted, the final registered probe
tested five literal or historically plausible transcription variants against
the same attributed C2 traces. The geometry, face-value centreline separations,
canal dimensions, scoring range, and every `CRITERIA-THEORY.md` threshold were
held fixed. A hit required pointwise, hump-position, and hump-amplitude
agreement in all four panels.

| variant | `S/L = 0.2` | 0.3 | 0.4 | 0.5 | four-panel result |
|---|---|---|---|---|---|
| resolved equations (baseline) | disagreement | partial | partial | agreement | miss |
| literal printed equation 4.29 | disagreement | partial | partial | agreement | miss |
| doubled interference cross term | disagreement | disagreement | disagreement | disagreement | miss |
| half nonzero-mode multiplicity | disagreement | partial | disagreement | disagreement | miss |
| doubled nonzero-mode multiplicity | disagreement | disagreement | disagreement | disagreement | miss |
| `cos` rather than `cos^2` pair coupling | disagreement | disagreement | disagreement | disagreement | miss |

The baseline Figure 359 errors are median `0.1038`, 90th percentile `0.2958`,
hump-position error `0.005`, and hump-amplitude ratio `0.845`. Literal equation
4.29 gives the same displayed Figure 359 metrics and only a negligible change
in the Figure 360 median (`0.0505` versus `0.0496`); it does not recover the
printed close-spacing amplitude. Doubling the cross term raises every hump too
far (`A = 1.185--1.253`). Halving the nonzero modes improves the Figure 360
amplitude but fails the other panels. Doubling those modes and replacing
`cos^2` by `cos` fail all four.

No variant is a plausible match to the published four-curve family under the
registered rule. Consequently there is no apparent hit whose historical use
needs to be labelled plausible-but-unverifiable. These misses do not establish
which unpublished implementation produced the thesis curves; they only remove
the five specified transcription mechanisms from the live explanations.

### Historical-variant claim classification

| claim | class | status |
|---|---|---|
| literal equation 4.29 reproduces the printed C2 family | E | rejected; Figure 359 still disagrees |
| an interference cross-term factor-of-two explains the family | E | rejected on all four panels |
| off-axis modal multiplicity explains the family | E | rejected for both half and double variants |
| `cos`/`cos^2` pair coupling explains the family | E | rejected on all four panels |
| one of these misses identifies an Insel implementation error | E | not claimable; unpublished implementation remains unknown |

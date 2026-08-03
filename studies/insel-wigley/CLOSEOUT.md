# Insel-Wigley validation closeout

## Final status

The direct theory-to-theory comparison closes as an **unresolved archival
annex**. The unresolved item is the amplitude of the printed 1990 C2 curves in
Figures 359--362, especially Figure 359 at `S/L = 0.2` (printed 358--359,
PDF 368--369). It is not an unresolved question about whether the current
library evaluates classical Michell theory correctly.

At `S/L = 0.2`, the printed curve peaks at `tau = 1.980` and `Fn = 0.450`; the
unbounded library peaks at `tau = 1.672` and `Fn = 0.445`. The registered hump
position agrees, but the amplitude ratio is `0.844` and the pointwise 90th-
percentile error is `0.299`. No tested water-geometry, attribution, separation,
or equation-transcription mechanism recovers all four printed curves. The
source code and mesh used to produce them are unavailable, so assigning the
remaining difference to Insel, the plot, or one untested implementation detail
would exceed the evidence.

The experimental result is unchanged: the solver is numerically resolved, but
the fixed-attitude comparison admits no continuous design-grade trust interval.
The useful positive result remains the wave-pattern interference ratio at
`S/L = 0.5`; the absolute wave coefficient and close spacings are materially
weaker.

The claim classes used here are **K**, known result reproduced; **E**,
engineering result for this codebase; and **P**, published practice matched.
No claim is presented as possibly novel.

## Complete hypothesis ledger

| hypothesis | outcome | decisive evidence | archive pointer |
|---|---|---|---|
| finite canal water geometry causes the close-spacing excess | **REFUTED** | registered gap closure is `0.004`, `0.010`, `0.019`, and `0.059` at `S/L = 0.2--0.5`; the two close-spacing values fail the `0.5` gate by about two orders of magnitude | `CRITERIA-CANAL.md`; `data/analysis/canal_scores.csv`; `FINDINGS.md` |
| the original trace followed the wrong family member | **REFUTED** | all four legends identify solid as C2, dashed as C3, dash-dot as C4, and dotted as C5 | `LEGEND-AUDIT-H1.md`; Figures 359--362, printed 358--359, PDF 368--369 |
| the registered NPL amplitude trend can be applied across NPL and Wigley forms | **INVALID AS REGISTERED** | Insel's text establishes an `L/B` trend within the NPL family, not a universal cross-form order; direct experiment places C3 above C2 by `0.109` at `S/L = 0.2` and `0.137` at `S/L = 0.3` | `CRITERIA-FAMILY-EXPERIMENT.md`; `data/digitized/final_experimental_family_checks.csv`; Figures 347--350, printed 352--353, PDF 362--363; discussion printed 124, PDF 134 |
| a gap, half-spacing, or centreline-spacing convention explains the four curves | **REFUTED** | gap-minus-beam and half-spacing fit zero panels; the face-value label fits only Figure 362; even free per-panel fits fail the full agreement gate for Figures 359--361 | `CRITERIA-SEPARATION.md`; `data/analysis/separation_fit_scores.csv`; `data/analysis/separation_transformation_scores.csv` |
| a tested historical equation or modal variant reproduces the family | **REFUTED** | literal equation 4.29, doubled cross term, half/double nonzero-mode multiplicity, and `cos` coupling all miss the unchanged gates on at least one panel; none fits all four | `data/analysis/historical_variant_scores.csv`; `FINDINGS.md` |
| C2 mesh resolution explains the excess | **LONGITUDINAL PART EXCLUDED; WHOLE-MESH CLAIM OPEN** | centred station differences recover the derivative of each quadratic Wigley waterline exactly; vertical corner averaging and replacing an element by a point source are not exact for the quadratic depth dependence and oscillatory phase | equation 4.53, printed 56, PDF 66; `FINDINGS.md` |

The last entry deliberately narrows the requested mesh conclusion. For a
quadratic waterline `q(s) = as^2 + bs + c`,

```text
[q(s+h/2) - q(s-h/2)] / h = 2as + b = q'(s).
```

Thus longitudinal station spacing cannot bias Insel's C2 source gradient.
Equation 4.53 also averages offsets between vertical corners and then places
the continuous element at one point (printed 56, PDF 66). Those operations are
not made exact by the identity. It would be inaccurate to say that every mesh
mechanism has been eliminated without the unpublished mesh or a dedicated
vertical/phase convergence study.

### Separation fits

| source figure | label `S/L` | fitted `S_eff/L` | registered fit interval | best median error | best 90th percentile | full agreement |
|---:|---:|---:|---:|---:|---:|---|
| 359 | 0.200 | 0.195 | 0.1875--0.2025 | 0.0913 | 0.2952 | no |
| 360 | 0.300 | 0.300 | 0.2775--0.3325 | 0.0413 | 0.1431 | no |
| 361 | 0.400 | 0.345 | 0.3275--0.3725 | 0.0215 | 0.1512 | no |
| 362 | 0.500 | 0.460 | 0.4225--0.5500 | 0.0185 | 0.0967 | yes |

These are fits to the solid C2 curves in Figures 359--362 (printed 358--359,
PDF 368--369). Their inconsistent movement with label separation, combined
with the unchanged 90th-percentile failures, rules out a shared convention
change rather than merely failing to identify one.

## Results independent of the archival puzzle

### Exact model and normalization

**K — The C2 geometry is independently reconstructed.** Insel identifies C2
as a Wigley model with parabolic waterlines and sections (printed 66--67,
PDF 76--77). Table 1 gives `L = 1.800 m`, `L/B = 10`, `B/T = 1.6`,
`C_B = 0.444`, and `C_P = C_M = 0.667` (printed 173, PDF 183). These values
give `B = 0.180 m`, `T = 0.1125 m`, and the ordinary product-parabola form,
which the library represents exactly as a single biquadratic surface.

**K/E — The comparison remains at tank-model scale and its normalization is
auditable.** The harness uses the physical model, freshwater at 15 degrees C,
and the recorded demihull wetted surface `0.482 m^2`; every row retains
dimensional resistance and reference area. Insel defines the coefficient with
`rho S_W U^2/2` (printed 64, PDF 74), while Table 1 supplies the C2 wetted area
(printed 173, PDF 183). No full-scale conversion or undocumented scale factor
is used.

### Classical published anchor and independent implementations

**K — The monohull reproduces the Doctors--Beck classical value.** Doctors and
Beck Table 1 reports `C_W = 1.2486e-3` at `Fn = 0.35` for the same
`B/L = 0.1`, `T/L = 0.0625` Wigley geometry. The model-scale study row gives
`1.248133344362e-3`, a relative difference of `0.037%`. This tests the
classical monohull result independently of Insel's multihull plots
([Doctors--Beck DOI record](https://doi.org/10.5957/jsr.1987.31.1.1)).

**E — Two in-house formulations agree where their domains meet.** The main
implementation evaluates the unbounded Michell integral using exact B-spline
moments. The separate study harness solves Insel's finite-width, finite-depth
modal equations and integrates the continuous Wigley source analytically.
Gate G1 takes the modal calculation toward wide, deep water and checks it
against both the library and a third analytic-Wigley quadrature at
`Fn = 0.25, 0.35, 0.50` and all four separations. Its 80-times-water endpoint
is within `0.5%` in monohull resistance and `0.02` in interference. Gate G2
shows that doubling the physical-tank modal sum from 128 to 256 modes changes
both observables by less than `0.002`. These are implementation gates, not
fits to the disputed curves.

### Equation 4.29 source erratum

**E — The finite-canal transcription resolves a dimensional omission in the
printed thesis.** Equation 4.25 contains
`K_0 + K_n cos^2(theta_n)` in the wave-elevation coefficient (printed 47,
PDF 57); the abbreviated `tau_n` in equation 4.29 omits it (printed 48,
PDF 58). The omission removes one inverse-length factor and is dimensionally
inconsistent. In the wide/deep limit,
`K_n cos^2(theta_n) = K_0`, so the missing amplitude factor is `2 K_0` and the
resistance misses exactly `4 K_0^2`. A literal implementation fails G1 by that
factor. Retaining the factor from the preceding governing equation restores
the independently required unbounded limit. This correction is an engineering
resolution of an internal source inconsistency, not a claim about Insel's
unavailable program.

### Southampton canal confinement

**K/E/P — The physical canal produces negligible principal-hump confinement
for this C2 comparison.** At `W = 3.7 m`, `H = 1.85 m`—the Southampton tank
dimensions (printed 62, PDF 72)—the canal changes the source-hump interference
from the unbounded result by only `0.0011`, `0.0014`, `0.0022`, and `0.0032`
at `S/L = 0.2--0.5`, at most about `0.25%` of the unbounded ordinate. That
supports Insel's estimate that C2 wall interference was below `1%` (printed
63, PDF 73). It does not erase the separate shallow-water caution near the
depth-critical band on the same page.

### Experimental trust envelope

**E — The preregistered experimental envelope is unchanged.** All original
prediction rows converged, but no configuration has three consecutive passing
`0.025 Fn` cells. The isolated passing fixed-attitude cells are:

| `S/L` | passing `Fn` cells |
|---:|---|
| 0.2 | 0.325, 0.350 |
| 0.3 | 0.350 |
| 0.4 | 0.325, 0.500 |
| 0.5 | 0.350, 0.500 |

The fixed C2-FX experimental markers come from Figures 135--144 (printed
244--248, PDF 254--258). The later archival probes change neither those data
nor the registered score. Absolute catamaran hump coefficients remain
29--61% above measured `C_WP`; interference progresses from disagreement at
`S/L = 0.2` to agreement at `S/L = 0.5`. This is a bounded validation result,
not a universal accuracy claim.

## Why the printed theory curves remain an annex

The evidence now separates three questions that were initially entangled:

1. The library reproduces classical monohull Michell theory and agrees with an
   independently implemented finite-canal formulation in the wide/deep limit.
2. Canal confinement at the actual tank dimensions is too small to create the
   disputed close-spacing amplitude.
3. The solid curves really are labelled C2, but their cross-model amplitude
   order at close spacing is opposite to Insel's experimental C2/NPL panels.
   The theoretical family is in Figures 359--362 (printed 358--359,
   PDF 368--369); the experimental family is in Figures 347--350 (printed
   352--353, PDF 362--363).

The separation grid and historical variants remove several simple archival
explanations. They do not reconstruct the unpublished 1990 executable. The
remaining live possibilities include untested details of vertical element
averaging or point-source phase quadrature, another undocumented implementation
choice, or a plotting/archival error. None is verified. The correct status is
therefore **unresolved printed-curve provenance, independently validated
library calculation**.

## Methods lessons

### Attribution is a protocol step

Curve identity should be committed from captions and legends before tracing.
This study originally recorded the right solid-C2 identity, but the later
attribution-only audit made that fact independently testable and prevented a
numerically convenient relabelling. The relevant legends are in Figures
359--362 (printed 358--359, PDF 368--369) and Figures 347--350 (printed
352--353, PDF 362--363).

### Family trends stay within their demonstrated family

Insel's reported lower-`L/B`, larger-interference trend concerns the NPL
round-bilge family; the same discussion separately says that amplitude depends
on hull form and separation (printed 124, PDF 134). Registering
`C3 >= C4 >= C2 >= C5` converted that NPL trend into an unsupported cross-form
rule. Future criteria should name the population over which an ordering has
actually been demonstrated.

### Stop conditions worked

Three preregistered stops fired and each prevented a favourable narrative from
outrunning the evidence:

1. `CRITERIA-THEORY.md` stopped the first closeout when Figure 359 materially
   disagreed.
2. `CRITERIA-CANAL.md` returned `UNEXPLAINED` when physical water geometry
   closed only `0.4%` and `1.0%` of the two close-spacing gaps.
3. `CRITERIA-H1.md` stopped on the Figure 359 family-order anomaly before any
   implementation variant could be selected to match it.

Each stop was later resumed by an explicit, separately preregistered question.
That sequence produced narrower conclusions rather than changing a threshold
after seeing a result.

## Design-use guidance

- Use the solver as a validated implementation of linear, fixed-geometry
  Michell theory, not as a stand-alone certificate of tank or full-scale total
  resistance.
- For this C2 hull, interference shape and level at `S/L = 0.5` are the
  strongest experimentally supported use. Do not describe any continuous
  `Fn` interval as design-grade under the registered joint-cell rule.
- Treat close spacings, especially `S/L = 0.2`, as trend guidance. Insel also
  cautions against low-speed small-separation prediction (printed 128,
  PDF 138) and attributes part of the phase difference to neighbour-induced
  asymmetric flow absent from symmetric thin-ship theory (printed 130,
  PDF 140).
- If experimental calibration is available, align the principal-hump phase as
  an explicit empirical correction rather than silently tuning the kernel;
  Insel describes that procedure at printed 125 (PDF 135).
- Prefer interference factors over direct absolute wave-resistance prediction
  when the design question permits it, consistent with Insel's conclusion at
  printed 131 (PDF 141). Preserve dimensional resistance and the coefficient
  area convention alongside any ratio.
- Do not transfer this envelope to NPL, transom, non-Wigley, running-attitude,
  or viscous-interference problems without new validation. The present solver
  does not model `beta`, measured trim/sinkage, wave breaking, or the
  wave-pattern-analysis transfer function.

## Reproducibility pointers

- Source interpretation and equations: `METHOD.md`
- Digitization protocol and attribution records: `DIGITIZATION.md`,
  `LEGEND-AUDIT-H1.md`, `LEGEND-AUDIT-FAMILY-EXPERIMENT.md`
- Frozen decisions: `CRITERIA.md`, `CRITERIA-THEORY.md`,
  `CRITERIA-CANAL.md`, `CRITERIA-FAMILY-EXPERIMENT.md`,
  `CRITERIA-SEPARATION.md`
- Complete numerical findings and claim classifications: `FINDINGS.md`
- Machine-readable scores: `data/analysis/`
- Deterministic predictions: `data/predictions/`
- Harness and commands: `harness/README.md`

The preregistration/evidence commit order is preserved in local history:
family rule `ab75d85` before legend audit `e9a3e3a` and evidence `017cb90`;
separation rule `f878a9a` before results `03019d0`; historical variants
`4f14e20`; and this closeout follows all evidence. Nothing in these rounds
modifies the frozen numerical kernel.

## Final validation

| command or check | result |
|---|---|
| `git diff --exit-code 30d8f53 -- crates python` | pass; no frozen-kernel diff |
| `cargo fmt --all --check` | pass |
| `cargo test --workspace` | pass; 215 tests, including doc tests |
| `cargo test --manifest-path studies/insel-wigley/harness/Cargo.toml` | pass; G1 and G2 both pass |
| `uv run --with pytest --with numpy pytest -q` from `python/` | pass; 11 tests |
| regenerate digitized archives, five prediction grids, all score tables, and all plots, then `git diff --exit-code HEAD -- studies/insel-wigley/data studies/insel-wigley/plots` | pass; byte-for-byte identical |

The four theory-overlay PNGs require the original Matplotlib cache to reproduce
their raster bytes in the sandboxed desktop environment. Regeneration with
that cache produced no diff. This is a rendering-environment detail; the
source, prediction, and score CSVs reproduced identically in either setting.

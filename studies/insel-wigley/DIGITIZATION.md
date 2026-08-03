# Insel C2 plot digitization record

The Insel thesis contains no numerical resistance tables. The values in
`data/digitized/` are readings of printed experimental markers, not transcribed
table cells and not samples of the fitted or theoretical curves. The source
correction and figure inventory are in `METHOD.md`.

## Blinding and admission rule

Three source-only readings were made before any `michell` prediction existed:

1. Pass A used image-component candidates from 360 dpi rendered pages, followed
   by a labelled visual overlay check. It retained raw pixel centres and did not
   snap Froude numbers to a presumed test grid.
2. Pass B was an independent review of the same rendered figures. It did not
   inspect pass A, any repository digitized CSV, or solver output. Its rule was
   deliberately conservative: retain only visibly separable marker glyphs.
3. Pass C was a second independent review of the free-attitude figures. It did
   not inspect either earlier pass or solver output. This supplemental pass is
   an adversarial audit of marker-family ambiguity; it is not used to admit
   final values.

Pass B is the canonical independent pass. A pass-A/pass-B pair is matched only
when the raw centres are within 18 horizontal and 24 vertical rendered pixels
on the same figure and observable. Matched centres are averaged. A pass-B-only
glyph is admitted with the base uncertainty after a direct source and curve-
family recheck. A pass-A-only candidate is excluded: it is never rescued from
the printed curve by numerical plausibility. The complete decisions are in
`data/digitized/passes/mismatches.csv`.

This rule admitted 874 markers: 796 two-pass matches and 78 conservative
pass-B-only glyphs. It excluded 525 pass-A-only candidates. The admitted counts
by source configuration are:

| attitude/configuration | two-pass | pass-B only | total |
|---|---:|---:|---:|
| fixed monohull | 49 | 3 | 52 |
| fixed `S/L = 0.2` | 41 | 2 | 43 |
| fixed `S/L = 0.3` | 68 | 6 | 74 |
| fixed `S/L = 0.4` | 82 | 3 | 85 |
| fixed `S/L = 0.5` | 61 | 3 | 64 |
| free monohull | 109 | 2 | 111 |
| free `S/L = 0.2` | 68 | 23 | 91 |
| free `S/L = 0.3` | 87 | 14 | 101 |
| free `S/L = 0.4` | 113 | 15 | 128 |
| free `S/L = 0.5` | 118 | 7 | 125 |

No unresolved marker is represented by a guessed value. Omission is preferred
when a fitted curve crosses a marker, several replicates merge, or a printed
symbol cannot be assigned to an observable without tracing the curve.

## Calibration and uncertainty

The canonical axis endpoints are preserved in
`data/digitized/passes/pass_b_calibration.tsv`. Every panel uses an affine map
from its rendered plot-border intersections. The coefficient panels span
`0.1 <= Fn <= 1.0`, `0.003 <= C_T <= 0.013`, and
`0 <= C_WP <= 0.010`. The shared running-attitude panels span -3 to 6 degrees
for trim and -0.1 to 0.2 for sinkage/draught. The inner plot border, rather than
the surrounding page box, is used; this distinction matters on Figure 174.

The base digitization uncertainties are `Fn = 0.002`, coefficient `0.00008`,
trim `0.05 degree`, and sinkage/draught `0.002`. For a two-pass match, the
reported uncertainty is the larger of the base value and half the separation
between the two calibrated readings. These are digitization uncertainties, not
Insel's experimental uncertainty. The final CSVs retain both raw centres,
source figure and page, pass count, and adjudication text so each value can be
reconstructed.

Repeated experimental markers are preserved. `replicate_id` is an archive row
identifier within an observable, not a claim that Insel assigned a run number.

## Independent-pass disagreement

Pass C deliberately probed the hard free-attitude panels more aggressively.
Against pass B it produced 779 candidates versus 556; 389 paired within the
same pixel tolerance, while 167 pass-B and 390 pass-C candidates remained
unpaired. The observable-level audit is
`data/digitized/passes/free_pass_b_c_audit.csv`.

Direct overlays showed why this is not harmless sampling variation: pass-C-only
groups included portions of fitted curves and assignments to the wrong marker
family, most visibly in the `S/L = 0.2` wave-pattern and attitude panels. The
canonical data therefore remain the conservative pass B, including where pass
C would add points that look numerically smooth. Pass C remains archived as a
negative audit result and a warning that automated component detection alone
is not reliable on these scans.

## Visual curve check

`plots/plot_digitized.py` plots only the final CSVs; it does not read the PDFs,
raw passes, mismatch log, or predictions. All ten plots were compared at full
rendered resolution with their source marker families:

- fixed C2-FX: Figures 135--144, printed pages 244--248 (PDF 254--258);
- free C2: Figures 161--179, printed pages 258--267 (PDF 268--277).

The final series follow the printed marker trajectories and preserve visible
replicate scatter. Sparse ranges correspond to genuinely ambiguous or merged
print, especially low-`Fn` wave-pattern markers; no fitted line was used to fill
those gaps. The PNGs under `plots/digitized-*.png` are the checked outputs.

## Reproduction and validation

Fetch the primary PDFs, render the relevant thesis pages as images, regenerate
pass A, reconcile it with the archived blind pass B, validate the archive, and
rebuild the plots:

```sh
cd studies/insel-wigley
./fetch.sh
mkdir -p tmp/pdfs/insel-wigley/thesis-c2-highres
pdftoppm -jpeg -r 360 -f 254 -l 277 \
  data/raw/insel-1990-thesis.pdf \
  tmp/pdfs/insel-wigley/thesis-c2-highres/page
cd ../..
uv run --project python --extra plot python \
  studies/insel-wigley/data/digitized/extract_fixed_pass_a.py
uv run --project python --extra plot python \
  studies/insel-wigley/data/digitized/extract_free_pass_a.py
python3 studies/insel-wigley/data/digitized/reconcile.py
python3 studies/insel-wigley/data/digitized/validate_digitization.py
uv run --project python --extra plot python \
  studies/insel-wigley/plots/plot_digitized.py
```

The extraction scripts also write labelled candidate overlays beneath the
ignored `tmp/` tree for manual inspection. Neither OCR nor PDF text extraction
is used. Reconciliation is deterministic from the committed raw pixel passes;
plot regeneration is deterministic from the committed final CSVs.

## C2 theoretical-interference curves

The later theory comparison uses the solid C2 curves in Figures 359--362
(printed pages 358--359, PDF pages 368--369). Those curves are distinct from
the experimental markers above and were digitized only after
`CRITERIA-THEORY.md` was committed. The source was rendered at 360 dpi. Each
panel was calibrated independently from its inner plot border with affine axes
`0.1 <= Fn <= 1.0` and `0 <= tau <= 2.5`.

Two source-only passes sampled common `Fn = 0.005` anchors. Pass A used local
pixel-ink centring around source-read guide ordinates. Pass B used a
continuity-constrained whole-curve path, independent plot-border readings, and
separately source-read checkpoints before local centring. Neither extraction
script reads the other pass, a solver prediction, or a comparison overlay.
The temporary labelled overlays were inspected only against the scanned source
to check C2 curve identity, including line-style crossings.

`reconcile_theory.py` calibrates each pass separately and admits an anchor only
when the readings differ by no more than `0.015` in `Fn` and `0.08` in `tau`,
as preregistered. It averages the two calibrated ordinates and records the
larger of the preregistered base uncertainty (`0.003` in `Fn`, `0.02` in
`tau`) and half the pass difference. A source-only reinspection left 36
over-tolerance or line-ambiguous anchors unresolved; all 36 are omitted in
`passes/theory_mismatches.csv` rather than assigned from curve smoothness.

| figure | `S/L` | raw anchors per pass | admitted | omitted | admitted 0.01 scoring anchors, `0.20 <= Fn <= 0.80` |
|---|---:|---:|---:|---:|---:|
| 359 | 0.2 | 131 | 121 | 10 | 54 |
| 360 | 0.3 | 161 | 153 | 8 | 57 |
| 361 | 0.4 | 161 | 153 | 8 | 59 |
| 362 | 0.5 | 161 | 151 | 10 | 57 |

Thus every panel exceeds the preregistered minimum of 50 of 61 scoring anchors
and covers the full principal-hump window. The reconciled curve is archived in
`data/digitized/theory_interference.csv`; it retains both raw pixel centres,
both calibrated values, pass differences, uncertainty, and source locations.

Reproduce the theory digitization from the fetched thesis without extracting
PDF text:

```sh
mkdir -p tmp/pdfs/insel-followup/theory-highres
pdftoppm -jpeg -r 360 -f 368 -l 369 \
  studies/insel-wigley/data/raw/insel-1990-thesis.pdf \
  tmp/pdfs/insel-followup/theory-highres/page
uv run --project python --extra plot python \
  studies/insel-wigley/data/digitized/extract_theory_pass_a.py
uv run --project python --extra plot python \
  studies/insel-wigley/data/digitized/extract_theory_pass_b.py
python3 studies/insel-wigley/data/digitized/reconcile_theory.py
python3 studies/insel-wigley/data/digitized/validate_digitization.py
```

## Figure 359 full-family attribution audit

The H1 follow-up separately traces all four legend-identified curves in Figure
359 (printed 358, PDF 368). This source-only audit is additive: it does not
modify `theory_interference.csv` or any prior pass. Both new passes record the
model, hull description, printed line style, and a relative-position
description at source-read checkpoints. Pass A uses a compact local ink
search; pass B uses an independently read calibration and checkpoints with a
wider horizontal score to bridge dash and dot gaps. Neither reads solver data
or the other pass.

The unchanged admission limits accept 503 of 510 anchors and omit seven. The
reconciled family, mismatch log, peak metrics, and preregistered family checks
are `theory_family_359.csv`, `passes/h1_figure_359_mismatches.csv`,
`h1_figure_359_family_metrics.csv`, and
`h1_figure_359_family_checks.csv`. The resolved amplitude-order reversal
triggers `CRITERIA-H1.md`'s source-anomaly stop, so no later figure or solver
overlay is part of this audit.

Reproduce the archived source trace and checks with:

```sh
mkdir -p tmp/pdfs/insel-h1-family/theory-highres
pdftoppm -jpeg -r 360 -f 368 -l 369 \
  studies/insel-wigley/data/raw/insel-1990-thesis.pdf \
  tmp/pdfs/insel-h1-family/theory-highres/page
uv run --project python --extra plot python \
  studies/insel-wigley/data/digitized/extract_family_359_pass_a.py
uv run --project python --extra plot python \
  studies/insel-wigley/data/digitized/extract_family_359_pass_b.py
python3 studies/insel-wigley/data/digitized/reconcile_family_359.py
```

## Experimental cross-model family audit

Figures 347--350 (printed 352--353, PDF 362--363) were rendered at 480 dpi
after `CRITERIA-FAMILY-EXPERIMENT.md` and the independent legend audit were
committed. The audit traces C2--C5 on `0.35 <= Fn <= 0.55`; the distinct
C2-FIXED curve is identified but is not an input to the registered free-model
cross-form test.

Both passes record model identity, line style, and relative position. Pass A
uses compact local ink centring. Pass B uses separate plot-border readings,
source checkpoints, and a wider horizontal response for broken styles. The
unchanged `0.015 Fn` and `0.08 tau` bounds admit 652 of 656 anchors and omit
four. The final curve, mismatches, peak metrics, and registered result are:

- `experimental_family_347_350.csv`;
- `passes/final_experimental_family_mismatches.csv`;
- `final_experimental_family_metrics.csv`; and
- `final_experimental_family_checks.csv`.

Reproduce the extraction and reconciliation with:

```sh
mkdir -p tmp/pdfs/insel-final/experimental-figures
pdftoppm -jpeg -r 480 -f 362 -l 363 \
  studies/insel-wigley/data/raw/insel-1990-thesis.pdf \
  tmp/pdfs/insel-final/experimental-figures/page
uv run --project python --extra plot python \
  studies/insel-wigley/data/digitized/extract_experimental_family_pass_a.py
uv run --project python --extra plot python \
  studies/insel-wigley/data/digitized/extract_experimental_family_pass_b.py
python3 studies/insel-wigley/data/digitized/reconcile_experimental_family.py
```

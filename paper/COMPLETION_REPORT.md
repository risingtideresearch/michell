# Paper A completion report

Date: 2 August 2026  
Target: *Journal of Ship Research* (JSR)  
Local branch: `story-paper-a-submission`  
Frozen numerical kernel: `30d8f53ba8fb48c3c77d3f2a0d95bddde74381e8`  
Remote operations: none

## Outcome

Paper A has been taken from technical draft to a locally reproducible JSR
submission package. The manuscript is a named, six-page, two-column paper by
Rob Story and Avi Bryant for the Rising Tide Research Foundation. It follows
the public SNAME journal template's visible constraints, carries a transparent
AI-use disclosure, links the public repository, and reserves a Zenodo DOI
slot. No journal was contacted, no submission was made, and no branch or tag
was pushed.

The paper's leading result is deliberately bounded: the real-axis marcher's
historical error diagnostic understated measured error by about 113-fold in a
reproduced case; the frozen correction is conservative over the tested Wigley
range. The endpoint--Bickley construction then gives an approximately
880-fold median improvement at `Fn=0.02` in the frozen 30-sample run, while
matching an independently coded real-axis reference to `1.526e-13`
relatively. The implementation is an engineering improvement. Only the
specific combination of known ideas is classified as possibly novel.

## Reproducible build and validation

Final toolchain and host:

- macOS Darwin 25.5.0, Apple M5 MacBook Air, 10 cores, 24 GB memory.
- `rustc 1.96.0 (ac68faa20 2026-05-25)`.
- `cargo 1.96.0 (30a34c682 2026-05-25)`.
- Tectonic 0.17.0.

Final validation commands:

```sh
cargo test --workspace
cd python
UV_CACHE_DIR=/private/tmp/michell-uv-cache \
  uv run --no-project --with pytest --with numpy python -m pytest -q
cd ..
MICHELL_BENCH_SAMPLES=30 cargo bench -p michell --bench wigley
cd paper
tectonic main.tex --outdir ../output/pdf --keep-logs --keep-intermediates
```

Results:

- Rust: 215 passed, zero failed, zero ignored; this total includes one
  documentation test. The longest CLI integration group completed in 92.57 s.
- Python: 11 passed in 0.08 s.
- Final verification benchmark: checksum
  `2.553251156101e5` for the 21-speed sweep and identical direct/endpoint
  resistance checksums to the frozen precision. The current host-state run
  measured 12.788 ms for the 21-speed median, 22.293 ms for the `Fn=0.02`
  marcher, and 0.030 ms for endpoint/NSD, a 736.96-fold median speedup. These
  timings are a verification rerun, not replacements for the frozen,
  same-process numbers reported in the paper.
- Manuscript evidence remains the committed 30-sample frozen run: 26.289 ms
  for the 21-speed median, 39.233 ms for the `Fn=0.02` marcher, 0.045 ms for
  endpoint/NSD, and 879.17-fold measured speedup (reported as approximately
  880-fold). Deterministic checksums and work counts match the final rerun.
- LaTeX: Tectonic completed with zero undefined references, zero undefined
  citations, and zero overfull boxes. Font-loader and underfull-box messages
  are non-fatal and were checked visually rather than suppressed.
- PDF: six letter-size pages, rendered and inspected page by page. Final file:
  `output/pdf/endpoint-bickley-michell-jsr.pdf`.

One initial Python invocation from the repository root was invalid because the
package is rooted in `python/`; rerunning the same suite from that directory is
the result reported above. No database-backed tests exist in this repository.

## Numerical-claim audit

Every manuscript measurement is mapped to its generator in
`paper/results/number-audit.md`. The final audit made one substantive wording
correction: the draft's unsupported “3--100” improvement range was replaced
with the four-case evidence, 3.2--4.5-fold lower measured error and 42--51%
fewer evaluations. The hard stale-number search found no occurrence of the
superseded 971.79-fold speedup, `3.558e-11` reduced error, pre-fix evaluation
counts, or cutoff `lambda=500` in `paper/main.tex`.

The paper's most important numerical anchors are:

| Quantity | Final manuscript value | Evidence |
|---|---:|---|
| Published Wigley anchor at `Fn=0.35` | `10^3 Cw = 1.247922`, 0.0543% from 1.2486 | Doctors and Beck (1987), Table 1; committed Phase-0 test |
| Historical error understatement at `Fn=0.05` | about 113-fold | `5.781e-7 / 5.123e-9` at pre-fix revision `64925d4` |
| Corrected marcher coverage at `Fn=0.05` | 3.87-fold | `5.774e-7 / 1.493e-7` |
| Endpoint/reference difference at `Fn=0.02` | `1.526e-13` | independent analytic-Wigley real-axis reference |
| Reference cutoff check, 4,000 to 8,000 | `1.227e-15` | same independent reference |
| Endpoint work at `Fn=0.02` | 1,152 nodes | 16 pairs times 72 contour nodes |
| Marcher work at `Fn=0.02` | 511,504 evaluations | benchmark diagnostics |
| Frozen median runtime at `Fn=0.02` | 39.233 ms to 0.045 ms | 30 warmed, same-process optimized samples |
| Frozen speedup | 879.17-fold; about 880-fold in prose | same benchmark |
| 21-speed checksum | `2.553251156101e5` | frozen and final verification runs |

## Claim-class audit

The manuscript follows the four-class honesty contract from `REPORT.md`.

| Paper claim | Class | Final disposition |
|---|---:|---|
| Michell theory, low-speed endpoint dominance, Bickley kernels, and numerical steepest descent | 1 — known result reproduced | Attributed to primary or authoritative sources; no novelty wording |
| Published Wigley coefficient and classical physical properties | 1 — known result reproduced | Reproduced by committed validation tests |
| Corrected marcher termination diagnostic | 2 — engineering improvement | Described as a reproduced historical defect and codebase correction, not a general theorem |
| Exact spline-span endpoint implementation and conservative hybrid dispatch | 2 — engineering improvement | Claimed only for the supported upright symmetric monohull path |
| Exact polynomial inner integration comparable with published practice | 3 — matches published state of the art | No breakthrough language and no head-to-head Michlet performance claim |
| Arbitrary-degree B-spline endpoint pairs plus Bickley continuation, frequency-scaled Gaussian contour, and submerged-endpoint gate | 4 — possibly novel | Preserved only with explicit negative-search, full-text, expert-review, and no-patent-search caveats |

The class-4 search covered publisher and DOI records, arXiv, TRID, the
University of Adelaide repository, DLMF, and backward references in the
verified sources. Searches included endpoint/asymptotic Michell resistance,
Bickley/Bickley--Naylor ship-wave kernels, polynomial and B-spline Michell
integrals, steepest-descent ship-wave quadrature, and combinations of those
terms. The closest unresolved source is de Sendagorta and Grases (1988), DOI
`10.5957/jsr.1988.32.1.19`: its record and abstract were verified, but the full
text was unavailable. No patent search was performed. The manuscript therefore
uses “possibly new combination,” never “breakthrough” or an unqualified claim
of priority.

## Submission and archive status

Prepared locally:

- JSR-formatted manuscript and review PDF.
- JSR cover-letter draft and submission checklist.
- Venue decision memo based on official SNAME and APNUM materials.
- `.zenodo.json`, DOI placeholder, and archive checklist.
- Exact frozen measurement output and numerical audit.
- Transparent AI-use acknowledgment.

Still requires human action before submission:

1. Choose the corresponding author and add contact details.
2. Obtain both authors' and any required foundation approval of the exact PDF.
3. Obtain and inspect the full de Sendagorta--Grases paper, then use the
   prepared novelty-paragraph variant if its contents change the comparison.
4. Confirm JSR's current blinding and generative-AI rules in the live
   ScholarOne form or with SNAME publications staff.
5. Confirm funding, conflicts, author contributions, permissions, membership
   information, and reviewer nominations requested by the live form.
6. Add a repository-level license if the foundation intends the package-wide
   MIT declaration currently present only in package manifests.
7. Push and release only after authorization, mint the Zenodo DOI, replace
   `10.5281/zenodo.REPLACE-ME`, rebuild, and inspect the archived PDF.
8. Consider a patent search or professional novelty review before making any
   stronger priority statement.

Skipped by design: live ScholarOne submission, journal contact, DOI minting,
remote push, patent search, and access-controlled full-text review. None is
reported as passed.

## Ranked post-submission backlog

1. Extend the endpoint kernels to multihull transverse phases using uniform or
   automated steepest-descent contours; this has the highest practical leverage.
2. Add interval-certified contour quadrature and floating-point error bounds so
   the full estimate, not only omission and Gaussian tails, is rigorous.
3. Differentiate the endpoint representation for cheap low-Froude design
   gradients and connect it to the existing exact general-route adjoint.
4. Retain shallow submerged endpoints through mixed linear/quadratic contours,
   widening the acceptance region above the present waterline-only reduction.
5. Evaluate stable imaginary-argument Bickley recurrences or the recent
   Bessel--Struve representation as a faster kernel backend.
6. Develop finite-depth/Sretensky endpoint kernels and test the changed
   dispersion singularities.
7. Compare systematically with Filon and Levin methods on broader polynomial
   hull families, including adverse cancellation and multiple chine knots.
8. Add a true external implementation comparison if a reproducible Michlet or
   equivalent executable and licensing path become available.

## Local history and handoff

The paper branch consists of one logical commit per revision stage: draft,
authorship, frozen measurements, number audit, literature/reproducibility
passes, venue selection, claim tightening, JSR formatting, and submission
logistics. The local annotated tag `paper-a-jsr-v1` is to identify the final
review candidate. It must not be pushed until the authors authorize release.

The final branch should be clean except for ignored LaTeX intermediates. No
remote operation is part of this completion.

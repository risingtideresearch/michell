# Paper A reconstruction completion report

Date: 3 August 2026
Target: *Journal of Ship Research* (JSR)
Local branch: `story-paper-a-r2`
Frozen measurement revision: `1a0deac42494c5775098a3601fd8677e9d23c02a`
Release tag: immutable local tag `paper-a-jsr-v3` on the final artifact commit
Remote operations: none

## Outcome

Paper A is reconstructed as “Error-Gated Endpoint Evaluation of Michell
Resistance for Low-Froude B-Spline Hulls.” The former method-priority claim is
withdrawn without a hedged remnant. The paper now presents a largely forgotten
analytic lineage—Birkhoff and Kotik, Michelsen's 1960 dissertation and 1972 JSR
sequel, and the verified Sendagorta--Grases record—revived as a modern,
validated, error-bounded implementation for B-spline hulls.

The historical reduction lineage is class 1/class 3. The degree-bounded
B-spline endpoint implementation, coefficient and cancellation accounting,
submerged-pair omission bound, contour evaluator, error gate, and real-axis
fallback are class 2. The paper makes no claim of method priority.

The final numerical claim is deliberately protocol-specific. In the frozen
30-sample paired run, the Fn=0.02 marcher and endpoint medians were 21.541 ms
and 0.030 ms, 718-fold apart and reported as about 700-fold. Repeated protocols
span roughly 608--865x. The manuscript reports the dispersion and this range
instead of treating a three-significant-figure ratio as portable.

## Old-to-new audit

| Topic | Reviewed v1 state | Final v2 state |
|---|---|---|
| Title | Previous Endpoint--Bickley reduction title | “Error-Gated Endpoint Evaluation…” |
| Historical framing | Near-prior-art discussion supporting a priority category | Equation-level Birkhoff--Kotik/Michelsen/Sendagorta lineage; no priority category |
| Spline scope | Unbounded-degree language | Both degrees explicitly validated only for 1--16; degree 17 refuses |
| Geometry scope | Upright symmetric monohull | Additionally restricted to coarse, well-separated active endpoints |
| Span consequence | Unstated | `m <= 1/(25 Fn^2)` for equal active spans; one short spacing forces fallback |
| Multi-span cost | Unstated | 25,920-node illustrative arithmetic versus Wigley's 1,152 nodes |
| Contour conditions | Rotation stated without full conditions | Principal branch, poles, branch points, sector, closing-arc decay, and `s>=4` stated |
| Kernel-order evidence | `s={4,7,10}` | Public maximum `s=98` plus margin at `s=128` |
| Contour rule | Fixed 24/48 nodes | 24/48 normally; 48/96 for stiff `s/abs(omega)>=2` kernels |
| Benchmark protocol | One unrelated warm call, fixed case order, median/best | Per-case warmup, alternating paired comparisons, median/IQR/best |
| Frozen runtime claim | Approximately 880-fold from 39.233/0.045 ms | About 700-fold from 21.541/0.030 ms, with protocol range disclosed |
| Rust suite | 215 tests in the reviewed report | 229 tests, including one doctest |
| Independent reference | One analytic Wigley marcher | Independently coded `lambda=1+t^2` and `lambda=sec(theta)` endpoint regularizations, each refined in order and cutoff |
| Dispatch estimate | Components described qualitatively | Equations for weighting, aggregation, construction, omission, contour difference, guarded denominator, floor, refusal, and combination |
| Marcher status | Accuracy rows could imply requested tolerance was met | Every row labels `RefinementCap`, work, estimate, and observed error |
| Archive | v1 tag had been moved during concurrent review | Immutable-tag rule preserved; v2 remains fixed and the revision receives v3 |

## P4 representativeness choice

The narrow option was chosen. The existing harness contains Wigley and one
synthetic full-multiplicity-chine geometry, not a curated corpus of realistic
multi-span designs. A cheap acceptance-rate survey would therefore measure an
arbitrary sampling distribution. Rather than present that as representative,
the abstract, algorithm section, limitations, and conclusion now scope the
contribution to degrees 1--16 and coarse, well-separated endpoint maps. The
paper gives the exact equal-span gate and labels the 25,920-node example as
arithmetic rather than a measured benchmark.

## Additional defect found and fixed

Extending the contour test to the public degree envelope exposed a genuine
coverage defect. At `s=98`, `omega=25`, the old 48-point contour value differed
from an independent resolved-real-axis value by `5.13e-6` relatively. Commit
`a4d011a` preserves the failing test. Commit `69f9ff4` adds a lazy 48/96 rule
for `s/abs(omega)>=2`; the permanent 18-case matrix covers
`s={4,7,10,50,98,128}` and `omega={25,100,400}` below the scaled `1e-9`
threshold. Ordinary Wigley cases retain 72 nodes per nonzero pair.

The addendum found a second genuine defect in the endpoint estimator. The
absolute-to-relative conversion used `B/abs(R_e)`, which is optimistic when
the requested error is relative to the unknown exact resistance. Commit
`195c633` preserves a failing regression; `101ac9b` now uses
`B/(abs(R_e)-B)` and refuses dispatch when `B >= abs(R_e)`. The manuscript
publishes every estimator component and decomposes all four accuracy rows.

Independent-reference due diligence is also stronger. Commit `ee306bb` adds
separately implemented `lambda=1+t^2` and `lambda=sec(theta)` endpoint-aware
references. Their 16-node results agree within `4.7e-15`; 16/24-node changes
are at most `6.2e-15`; and both cutoff-doubling changes at Fn=0.02 are about
`1e-15`. At Fn 0.03 and 0.02 the physical omission terms are `6.153e-32` and
`5.508e-70`, so the manuscript correctly describes those rows as quadrature
agreement rather than omission-physics validation.

Commit `fc39e8c` quantifies the fixed 48-node contour boundary. It is accurate
through approximately `s=50` and degrades above that range, while the 24/48
difference covers the meaningful 48/192 discrepancy in scans to `s=400`.
This is the quantitative basis for the adaptive 48/96 rule and public
`s<=98` envelope.

## Claim-class audit

| Claim | Class | Disposition |
|---|---:|---|
| Michell theory and published Wigley coefficient | 1 | Known results reproduced by committed tests |
| Birkhoff--Kotik kernel separation and Michelsen polynomial/Gegenbauer reduction lineage | 1/3 | Attributed through inspected equations or verified records |
| Sendagorta--Grases shape/velocity separation and design-tabulation program | 3 | Described only to the extent established by its verified abstract |
| Low-speed endpoint dominance, Bickley kernels, and steepest descent | 1 | Known ingredients, explicitly attributed |
| Validated degree envelope, endpoint construction, cancellation accounting, omission bound, contour evaluation, dispatch, and fallback | 2 | Engineering improvements to this codebase |
| High-order contour rule and its regression matrix | 2 | Red-before-fix engineering correction |
| Guarded endpoint relative-bound conversion | 2 | Red-before-fix engineering correction; conditioned on `B < abs(R_e)` |
| Two endpoint-aware real-axis references and convergence study | 1/2 | Tuck's regularization reproduced; independent validation engineering added here |

The downgrade is justified by direct inspection of Michelsen's 1960 equations
3.3--3.4, 3.21--3.22, and 3.38--3.40; Wehausen's equations 39--45 describing
the Birkhoff--Kotik lineage; and verified records for Michelsen (1972) and de
Sendagorta and Grases (1988). Full method-by-method performance benchmarking is
future work.

## Validation

Host and toolchain:

- Darwin 25.5.0, Apple M5 MacBook Air, arm64.
- `rustc 1.96.0 (ac68faa20 2026-05-25)`.
- `cargo 1.96.0 (30a34c682 2026-05-25)`.
- Tectonic 0.17.0.

Commands and results:

```sh
cargo test --workspace
```

Result: 229 passed, zero failed, zero ignored, including one doctest. The
longest CLI integration group completed in 91.57 s.

```sh
cd python
UV_CACHE_DIR=/private/tmp/michell-uv-cache \
  uv run --no-project --with pytest --with numpy python -m pytest -q
```

Result: 11 passed in 0.10 s. The first sandboxed invocation could not resolve
PyPI; the authorized network retry obtained the declared environment and is
the result reported here.

```sh
MICHELL_BENCH_SAMPLES=30 ./paper/reproduce_measurements.sh
```

Result: both focused validation tests passed; sweep checksum
`2.553251156101e5`; direct and endpoint 30-run checksums
`1.325169498584e-3` and `1.325170337908e-3`; frozen timing and IQR values are
recorded in `paper/results/frozen-paper-a-r2-2026-08-03.txt`. Fast cases used
256-call timed batches; the Fn=0.02 medians were 21.541 ms and 0.030 ms.

```sh
cd paper
tectonic main.tex --outdir ../output/pdf --keep-logs --keep-intermediates
```

Result: zero undefined references, zero undefined citations, and zero overfull
boxes. The nine letter-size pages were rendered at 120 dpi and inspected page
by page. Equations, tables, links, and margins are intact; the accuracy/work
figure uses black solid/dashed lines and distinct circle/square markers, so it
does not depend on color discrimination.

The exported-tag validation extracts `paper-a-jsr-v3` without `.git`, verifies
the tar listing contains no `.git` entry, and runs
`MICHELL_BENCH_SAMPLES=3 paper/reproduce_measurements.sh` from that exact tree.
The checksum is stored beside the ignored local tarball under `output/archive/`.

## Source and submission status

Still requires Rob or the authors:

1. Keep the interlibrary-loan requests open for Michelsen (1972) and de
   Sendagorta and Grases (1988). Their verified records suffice for the current
   wording; update the table with positive equation-level facts when received.
2. Choose the corresponding author and complete contact details.
3. Obtain both authors' and any required foundation approval of the exact PDF.
4. Confirm live JSR blinding and generative-AI requirements, funding,
   conflicts, contributions, permissions, membership information, and reviewer
   nominations.
5. Push only after authorization, create the release, mint the Zenodo DOI,
   replace `10.5281/zenodo.REPLACE-ME`, rebuild, inspect, and create a new tag.
   Do not move `paper-a-jsr-v2` or `paper-a-jsr-v3`.
6. Decide whether to add a repository-level license before public archiving.

The DOI is explicitly a placeholder; no text asserts that the archive already
exists. No journal was contacted, no submission was made, no DOI was minted,
and nothing was pushed.

## Local commit map

| Commit | Logical change |
|---|---|
| `649c351` | Reclassify claims and restore the historical analytic lineage |
| `a4d011a` | Add failing high-order contour-kernel coverage |
| `69f9ff4` | Resolve stiff high-order endpoint kernels |
| `a54acb2` | Bound endpoint scope and state contour conditions |
| `cd22f44` | Pair benchmark comparisons and report dispersion |
| `2f4d2fa` | Make exported archives self-identifying and reproducible |
| `195c633` | Preserve the failing endpoint relative-bound conversion test |
| `101ac9b` | Guard the denominator and correct the relative bound |
| `ee306bb` | Add independent endpoint-reference convergence |
| `fc39e8c` | Measure fixed contour-rule degree dependence |
| `701711a` | Batch sub-millisecond timing cases |
| `1c49e5a` | Publish the estimator and self-contained validation tables |
| `7ce69bc` | Extend and qualify the analytic lineage |
| `1a0deac` | State the measured contour envelope |
| `3bbc636` | Freeze reference and timing evidence |
| final artifact commit | Freeze evidence, reports, and the visually checked PDF |

The final branch is expected to be clean after the ignored local archive is
created. The v1 retagging incident is recorded here to make the standing rule
unambiguous: release tags are immutable; every revised state receives a new
tag.

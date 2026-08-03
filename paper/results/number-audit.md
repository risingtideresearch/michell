# Paper A numerical audit

Unless noted otherwise, new values come from
`./paper/reproduce_measurements.sh` on frozen kernel base `30d8f53`; captured
output is in `frozen-kernel-2026-08-02.txt`. Timings used the optimized profile,
30 warmed samples, and default `rel_tol=1e-5` on the same development host.

| Manuscript quantity | Draft at `64925d4` | Frozen kernel | Generator |
|---|---:|---:|---|
| Marcher rel. diff., Fn 0.08 | 3.125e-7 | 6.917e-8 | `cargo test -p michell --test low_froude endpoint_reduction_converges_as_froude_number_falls -- --nocapture` |
| Marcher rel. diff., Fn 0.05 | 5.781e-7 | 1.493e-7 | same |
| Marcher rel. diff., Fn 0.03 | 1.063e-6 | 3.322e-7 | same |
| Marcher rel. diff., Fn 0.02 | 2.024e-6 | 6.333e-7 | same |
| Marcher estimate, Fn 0.05 | 5.123e-9 (pre-fix) | 5.774e-7 | same |
| Marcher evaluations, Fn 0.08 | 114,432 | 56,032 | same |
| Marcher evaluations, Fn 0.05 | 235,808 | 119,936 | same |
| Marcher evaluations, Fn 0.03 | 501,904 | 267,856 | same |
| Marcher evaluations, Fn 0.02 | 877,424 | 511,504 | same |
| 21-speed median / best (ms) | 21.897 / 21.584 | 26.289 / 22.843 | `MICHELL_BENCH_SAMPLES=30 cargo bench -p michell --bench wigley` |
| Default API Fn 0.05 median / best (ms) | 0.029 / 0.028 | 0.050 / 0.050 | same |
| Marcher Fn 0.02 median / best (ms) | 28.141 / 27.794 | 39.233 / 35.950 | same |
| Endpoint Fn 0.02 median / best (ms) | 0.029 / 0.028 | 0.045 / 0.044 | same |
| Fn 0.02 median speedup | 971.79x | 879.17x measured; approximately 880-fold reported | same |
| 21-speed checksum | 2.553250998079e5 | 2.553251156101e5 | same |
| Rust test count | 196 | 215 | `cargo test --workspace` |
| Independent reference cutoff | 500 | 4,000 (8,000 cutoff check at Fn 0.02) | `cargo test -p michell --test low_froude endpoint_reduction_converges_as_froude_number_falls -- --nocapture` |
| Reduced rel. diff., Fn 0.05 | 1.058e-11 | 2.070e-13 | same |
| Reduced rel. diff., Fn 0.03 | 2.328e-11 | 1.034e-12 | same |
| Reduced rel. diff., Fn 0.02 | 3.558e-11 | 1.526e-13 | same |
| Fn 0.02 reference cutoff check | not measured | 1.227e-15 relative (4,000 to 8,000) | same |
| Published Wigley anchor, Fn 0.35 | 1.2486e-3 published | 1.2479219624e-3 reproduced; 5.430383e-4 relative difference | `cargo test -p michell --test phase0_validation wigley_matches_published_thin_ship_value -- --nocapture` |

The checksum change is intentional: the corrected marcher retains a positive
tail that the old stopping rule truncated. Timing changes are host-state
sensitive; deterministic resistance values, evaluation counts, and checksums
are the primary reproducibility data.

## Derived and structural claims

| Manuscript claim | Value | Derivation or generator |
|---|---:|---|
| Historical error understatement | 112.84-fold; reported as about 113-fold | `5.781e-7 / 5.123e-9` using the pre-fix values above |
| Frozen/pre-fix marcher error improvement over the four reported low-Froude cases | 3.20--4.52-fold; reported as 3.2--4.5 | row-wise ratios of the draft and frozen marcher relative differences above |
| Frozen evaluation-count reduction over the same cases | 41.7--51.0%; reported as 42--51% | row-wise reductions from the draft to frozen evaluation counts above |
| Gaussian-tail bound at the frequency gate | below `4e-30` | `exp(-64) / (8 sqrt(2*25))` from the analytic bound in `paper/main.tex` |
| Wigley frequency gate | `Fn <= 0.2` | `omega = 1/Fn^2` and `omega >= 25` |
| Wigley endpoint/contour work | 1,152 scalar nodes | 16 nonzero-frequency pairs times `(24 + 48)` Gauss--Legendre nodes |
| Accepted reduced/reference range, Fn 0.05--0.02 | `1.526e-13`--`1.034e-12` | frozen low-Froude test rows above |
| Fn 0.05 corrected marcher-estimate coverage | 3.87-fold | `5.774e-7 / 1.493e-7` |
| Published Wigley difference | 0.05430383% | `100 * 5.430383e-4` |
| Fn 0.02 cutoff check | `1.227e-15` | independent analytic-Wigley reference with cutoffs 4,000 and 8,000 |
| Kernel test matrix | 9 cases, scaled discrepancy below `1e-9` | `s in {4,7,10}` crossed with `abs(omega) in {25,100,400}` in the low-Froude unit test |
| Rust/Python test totals | 215 / 11 | `cargo test --workspace` (214 ordinary plus one documentation test) and Python `pytest` |

The runtime table remains the frozen, same-process 30-sample measurement used
to prepare the manuscript. A later verification run may differ in absolute
time because of host state; matching checksums and work counts, not timing
identity, are the accuracy invariant.

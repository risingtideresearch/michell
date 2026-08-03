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

The checksum change is intentional: the corrected marcher retains a positive
tail that the old stopping rule truncated. Timing changes are host-state
sensitive; deterministic resistance values, evaluation counts, and checksums
are the primary reproducibility data.

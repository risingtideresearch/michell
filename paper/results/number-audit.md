# Paper A numerical audit

Frozen measurement revision: `1a0deac42494c5775098a3601fd8677e9d23c02a`
Captured output: `frozen-paper-a-r2-2026-08-03.txt`

The frozen run used `MICHELL_BENCH_SAMPLES=30`, the optimized Rust benchmark
profile, no custom `RUSTFLAGS`, and default `rel_tol=1e-5` on a 10-core Apple
M5 MacBook Air with 24 GB RAM, Darwin 25.5.0, Rust/Cargo 1.96.0, and LLVM
22.1.2. Each benchmark case received one unmeasured warm call. Calls below
0.1 ms were timed in batches of 256 and divided by 256. Direct/reduced and
adjoint/finite-difference comparisons alternated order on each sample. The
table reports the upper-middle median and empirical
quartiles selected by integer indices `n/4` and `3n/4` after sorting 30 samples.

## Frozen numerical values

| Manuscript quantity | Frozen value | Generator |
|---|---:|---|
| Marcher rel. diff., Fn 0.08 / 0.05 / 0.03 / 0.02 | `6.918e-8` / `1.493e-7` / `3.322e-7` / `6.334e-7` | `cargo test -p michell --test low_froude endpoint_reduction_converges_as_froude_number_falls -- --nocapture` |
| Marcher estimate, Fn 0.08 / 0.05 / 0.03 / 0.02 | `2.896e-7` / `5.974e-7` / `1.320e-6` / `2.517e-6` | same |
| Marcher status, all four rows | `RefinementCap` | same |
| Marcher evaluations, Fn 0.08 / 0.05 / 0.03 / 0.02 | 56,032 / 119,936 / 267,856 / 511,504 | same |
| Reduced rel. diff., Fn 0.05 / 0.03 / 0.02 | `1.963e-13` / `1.045e-12` / `1.829e-13` | same |
| Endpoint-map difference, Fn 0.08 / 0.05 / 0.03 / 0.02 | `1.284e-15` / `3.280e-16` / `4.664e-15` / `1.534e-15` | same |
| Reference cutoff check, 4,000 to 8,000 at Fn 0.02 | `1.381e-15` (`1+t^2`); `1.227e-15` (`sec theta`) | same |
| 21-speed median; IQR; best | 13.059 ms; 12.952--13.126 ms; 12.686 ms | `MICHELL_BENCH_SAMPLES=30 cargo bench -p michell --bench wigley` |
| Default API Fn 0.05 median; IQR; best | 0.030 ms; 0.030--0.030 ms; 0.030 ms | same |
| Marcher Fn 0.02 median; IQR; best | 21.541 ms; 21.311--21.666 ms; 20.621 ms | same |
| Endpoint Fn 0.02 median; IQR; best | 0.030 ms; 0.030--0.030 ms; 0.030 ms | same |
| Fn 0.02 frozen median speedup | `718.06x`; about 700-fold in prose | same |
| 21-speed checksum | `2.553251156101e5` | same |
| Published Wigley anchor, Fn 0.35 | `1.2479219624e-3`, 0.05430383% from published `1.2486e-3` | Phase-0 published-value test |
| Full suites | 229 Rust tests, including one doctest; 11 Python tests | `cargo test --workspace`; Python `pytest` |

Repeated paired benchmark runs span roughly 608--865x, including 702--718x
under the final batched protocol. The manuscript reports the paired frozen
result as about 700-fold and states the measured protocol sensitivity; it does
not present a three-significant-figure point ratio as portable.

## Derived and structural claims

| Claim | Value | Derivation or generator |
|---|---:|---|
| Historical error understatement | about 113-fold | `5.781e-7 / 5.123e-9` at pre-fix revision `64925d4` |
| Corrected marcher-estimate coverage at Fn 0.05 | 4.00-fold | `5.974e-7 / 1.493e-7` |
| Gaussian-tail bound at the frequency gate | below `4e-30` | `exp(-64) / (8 sqrt(2*25))` |
| Wigley frequency gate | `Fn <= 0.2` | `omega = 1/Fn^2`, `omega >= 25` |
| Equal-span frequency gate | `m <= 1/(25 Fn^2)` | adjacent spacing `L/m` gives `omega_adj=1/(m Fn^2)` |
| Illustrative multi-span work | 25,920 nodes | 360 ordinary pairs times 72 nodes; arithmetic example, not runtime evidence |
| Wigley endpoint work | 1,152 nodes | 32 nonzero-frequency ordered pairs reduce by conjugacy to 16 kernel evaluations, each using 72 nodes |
| Kernel order maximum | `s=98` | `n_max=p+2q+2=50`, then `s_max=2n_max-2`, for `p,q<=16` |
| Kernel test matrix | 18 cases through `s=128` | `s in {4,7,10,50,98,128}` crossed with `abs(omega) in {25,100,400}` |
| Fixed-rule degree scan | accurate through `s=50`; maximum discrepancy `3.157e-1` above; minimum 24/48 coverage `2.222x` | scan through `s=400`, `omega in {25,100,400}`, 192-node same-contour reference with `1e-10` comparison floor |

The reduced estimator sums omission, coarse/fine contour, and endpoint-map
accumulation in absolute units, refuses when their sum `B >= abs(R_e)`, and
otherwise converts with `B/(abs(R_e)-B)`. The independent test deliberately
distinguishes this guarded conversion from the optimistic `B/abs(R_e)` form.
The per-row relative decomposition is:

| Fn | Omission | Contour | Endpoint rounding | Coefficient/safety floor |
|---:|---:|---:|---:|---:|
| 0.08 | `3.815e-5` | `1.355e-12` | `9.153e-15` | `2.000e-8` |
| 0.05 | `3.656e-12` | `7.222e-14` | `9.702e-15` | `2.000e-8` |
| 0.03 | `6.153e-32` | `3.462e-12` | `9.292e-15` | `2.000e-8` |
| 0.02 | `5.508e-70` | `7.620e-13` | `9.428e-15` | `2.000e-8` |

At Fn 0.03 and 0.02 the physical omitted contribution is effectively zero;
the reduced/reference column therefore measures quadrature agreement, not the
physics of the omission approximation.

The high-order extension was red before the fix: the old 48-point result at
`s=98`, `omega=25` differed from the independent real-axis reference by
`5.13e-6` relatively. Commit `a4d011a` records the failing coverage; `69f9ff4`
selects a 48/96 rule when `s/abs(omega) >= 2`, while ordinary Wigley kernels
retain the 24/48 rule and 1,152-node count.

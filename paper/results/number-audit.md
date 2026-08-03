# Paper A numerical audit

Frozen measurement revision: `2f4d2faacd301537e4e28ccc6349f9722fc9fe94`
Captured output: `frozen-paper-a-r2-2026-08-03.txt`

The frozen run used `MICHELL_BENCH_SAMPLES=30`, the optimized Rust benchmark
profile, and default `rel_tol=1e-5` on an Apple M5 MacBook Air running Darwin
25.5.0 with Rust/Cargo 1.96.0. Each benchmark case received one unmeasured warm
call. Direct/reduced and adjoint/finite-difference comparisons alternated order
on each sample. The table reports the upper-middle median and empirical
quartiles selected by integer indices `n/4` and `3n/4` after sorting 30 samples.

## Frozen numerical values

| Manuscript quantity | Frozen value | Generator |
|---|---:|---|
| Marcher rel. diff., Fn 0.08 / 0.05 / 0.03 / 0.02 | `6.918e-8` / `1.493e-7` / `3.322e-7` / `6.334e-7` | `cargo test -p michell --test low_froude endpoint_reduction_converges_as_froude_number_falls -- --nocapture` |
| Marcher estimate, Fn 0.05 | `5.974e-7` | same |
| Marcher evaluations, Fn 0.08 / 0.05 / 0.03 / 0.02 | 56,032 / 119,936 / 267,856 / 511,504 | same |
| Reduced rel. diff., Fn 0.05 / 0.03 / 0.02 | `2.070e-13` / `1.034e-12` / `1.526e-13` | same |
| Reference cutoff check, 4,000 to 8,000 at Fn 0.02 | `1.227e-15` relative | same |
| 21-speed median; IQR; best | 13.035 ms; 12.996--13.064 ms; 12.945 ms | `MICHELL_BENCH_SAMPLES=30 cargo bench -p michell --bench wigley` |
| Default API Fn 0.05 median; IQR; best | 0.030 ms; 0.030--0.030 ms; 0.029 ms | same |
| Marcher Fn 0.02 median; IQR; best | 21.506 ms; 21.464--21.557 ms; 21.412 ms | same |
| Endpoint Fn 0.02 median; IQR; best | 0.031 ms; 0.030--0.031 ms; 0.030 ms | same |
| Fn 0.02 frozen median speedup | `702.23x`; about 700-fold in prose | same |
| 21-speed checksum | `2.553251156101e5` | same |
| Published Wigley anchor, Fn 0.35 | `1.2479219624e-3`, 0.05430383% from published `1.2486e-3` | Phase-0 published-value test |
| Full suites | 227 Rust tests, including one doctest; 11 Python tests | `cargo test --workspace`; Python `pytest` |

Two immediately repeated paired benchmark runs gave 702.23x and 706.51x.
An earlier cold-host paired run gave 607.94x. The former fixed-order reviewer
protocol gave approximately 750--865x. The manuscript reports the paired frozen
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
| Wigley endpoint work | 1,152 nodes | 16 nonzero-frequency pairs times 72 nodes |
| Kernel order maximum | `s=98` | `n_max=p+2q+2=50`, then `s_max=2n_max-2`, for `p,q<=16` |
| Kernel test matrix | 18 cases through `s=128` | `s in {4,7,10,50,98,128}` crossed with `abs(omega) in {25,100,400}` |

The high-order extension was red before the fix: the old 48-point result at
`s=98`, `omega=25` differed from the independent real-axis reference by
`5.13e-6` relatively. Commit `a4d011a` records the failing coverage; `69f9ff4`
selects a 48/96 rule when `s/abs(omega) >= 2`, while ordinary Wigley kernels
retain the 24/48 rule and 1,152-node count.

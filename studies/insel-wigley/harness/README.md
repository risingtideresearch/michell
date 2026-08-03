# Insel-Wigley prediction harness

This isolated binary evaluates the exact biquadratic C2 Wigley reconstruction
at every digitized fixed/free `C_WP` Froude number and on the pre-registered
`Fn = 0.25:0.005:0.55` grid. It does not modify or join the repository
workspace.

```sh
cargo run --release --manifest-path studies/insel-wigley/harness/Cargo.toml
```

The deterministic output is
`studies/insel-wigley/data/predictions/predictions.csv`. Catamaran rows are
computed with `multihull_wave_resistance_with` at centreplane placements
`y = +/- S/2`. The recorded interference is the library's documented ratio,
`R_w,pair / (2 R_w,solo)`, using a standalone solve made through the same
multihull path and options.

Every solve requests `rel_tol = 1e-6` and the binary aborts on a non-converged
outcome. Both pair and standalone diagnostics are retained in every row.

The preregistered theory comparison and critical-Froude check use a separate
`0.20 <= Fn <= 0.95` grid at `0.005` spacing for all four catamaran
separations. It is generated without rewriting the frozen experimental-grid
output:

```sh
cargo run --release --manifest-path studies/insel-wigley/harness/Cargo.toml -- --theory
```

That mode writes `data/predictions/theory_predictions.csv` and retains the
same geometry, fluid properties, solver options, convergence assertions, and
standalone-path normalization.

The independent finite-width, finite-depth canal reference and its
wide/deep-limit and mode-doubling gates are exercised with:

```sh
cargo test --manifest-path studies/insel-wigley/harness/Cargo.toml
```

Its transcription, the resolved equation (4.25)/(4.29) inconsistency, and the
explicit modal truncation rule are recorded in `../METHOD.md`. The tests do not
write prediction data or expose the physical-tank comparison.

After `CRITERIA-CANAL.md` is committed, the physical `W = 3.7 m`, `H = 1.85 m`
canal grid and its preregistered scores, plots, and critical-Froude diagnostics
are regenerated with:

```sh
cargo run --release --manifest-path studies/insel-wigley/harness/Cargo.toml -- --canal
python/.venv/bin/python studies/insel-wigley/plots/analyze_canal.py
python/.venv/bin/python studies/insel-wigley/plots/critical_froude.py
```

The first command writes `data/predictions/canal_predictions.csv`; the analysis
commands write `data/analysis/canal_scores.csv`, update
`data/analysis/critical_froude.csv`, and regenerate the eight `canal-*.png`
overlays. The Python environment is the repository's existing plotting extra.

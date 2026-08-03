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

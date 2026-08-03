# Adversarial high-degree evidence

These files are exact copies of the temporary programs supplied by the three
independent kernel reviewers. They are preserved here so the provenance and
independence of the regression references remain auditable after the original
`/private/tmp` worktrees disappear.

| Reviewer | Preserved file | SHA-256 |
| --- | --- | --- |
| R1 | `r1/analytic_ref.rs` | `9e8354506e7cf2f93b3c00afb67ec0eae081cb39d64bc7107b52a077c36f810c` |
| R1 | `r1/degree_threshold.rs` | `f6d51742d913fa836716edcb424eae8a9f52241b89689e3a924b7dd8bc69fe51` |
| R1 | `r1/scope_probes.rs` | `7066075cb9d3e4e80acb37f788f1b17d1d222876d4de87895837928c7902b35e` |
| R2 | `r2/high_degree.rs` | `f11132a6f219481e75293f6fdeb82ddc9c292c4ce7ad25bf673b08c84a1a576e` |
| R3 | `r3/reference.py` | `79ca513e49d6cb2bac5fa9d75d212da9b7ae714c8bba6a4e7b5b73745a3fe61a` |

R1 derives a closed-form inner amplitude for its single-span Bernstein family
and integrates the outer Michell integral with independently implemented,
phase-resolved 32-point Gauss-Legendre quadrature. Its degree sweep uses the
analytic Bernstein scaling of the same geometry. R2 supplies a variation-based
upper bound independent of either production solver. R3 recursively evaluates
the B-spline basis with `mpmath`, reconstructs span polynomials by independent
high-precision interpolation, performs its own integration by parts, and uses
independently generated quadrature for the endpoint kernels.

The committed Rust regression tests port the required reviewer cases. Run them
with:

```sh
cargo test -p michell --test high_degree_hardening
```

The original R3 oracle can be rerun (after building the reviewer's compatible
stress-probe binary) with:

```sh
uv run --with mpmath --with numpy python probes/adversarial/r3/reference.py \
  --binary target/release/stress-probe --case degree24x16_deep_fn005
```

The independently measured R3 reference is
`0.0063825284969346485 N`; the unhardened production path differed by
`1.3055972e-5` relative while reporting `7.4172953e-14` and `Converged`.

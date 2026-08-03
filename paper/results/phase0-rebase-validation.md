# Phase 0 rebase validation

- Submission branch: `story-paper-a-submission`
- Frozen kernel base: `30d8f53`
- Replayed paper commits: `2ea3b06`, `692c4ec`
- Command: `cargo test --workspace`
- Result: 215 passed, 0 failed, 0 ignored, including one documentation test.
- Command: `tectonic main.tex --outdir ../output/pdf --keep-logs --keep-intermediates`
  from `paper/`
- Result: pass with zero undefined references, zero undefined citations, and
  zero overfull boxes in `output/pdf/main.log`.

The endpoint solver invariance is covered directly by:

- `low_froude::tests::endpoint_decomposition_reproduces_exact_inner_amplitude`
- `low_froude::tests::endpoint_decomposition_handles_multiple_spans_and_a_chine`
- `low_froude::tests::steepest_descent_kernel_matches_resolved_real_axis_reference`
- `endpoint_reduction_converges_as_froude_number_falls`
- `low_froude_reported_error_covers_actual_error`
- `default_solver_accepts_only_a_reduction_within_tolerance`

All six passed in the full workspace run.

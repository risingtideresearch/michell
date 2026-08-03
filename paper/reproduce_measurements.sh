#!/bin/sh
set -eu

samples="${MICHELL_BENCH_SAMPLES:-30}"

if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git rev-parse HEAD
elif [ -f ARCHIVE_REVISION ]; then
    sed -n '1p' ARCHIVE_REVISION
else
    printf '%s\n' 'revision unavailable (non-Git source tree)'
fi
rustc --version
cargo --version
uname -a

cargo test -p michell --test low_froude \
  endpoint_reduction_converges_as_froude_number_falls -- --nocapture
cargo test -p michell --test low_froude \
  default_solver_accepts_only_a_reduction_within_tolerance -- --nocapture
MICHELL_BENCH_SAMPLES="$samples" cargo bench -p michell --bench wigley

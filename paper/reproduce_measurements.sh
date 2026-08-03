#!/bin/sh
set -eu

samples="${MICHELL_BENCH_SAMPLES:-30}"

git rev-parse HEAD
rustc --version
cargo --version
uname -a

cargo test -p michell --test low_froude \
  endpoint_reduction_converges_as_froude_number_falls -- --nocapture
cargo test -p michell --test low_froude \
  default_solver_accepts_only_a_reduction_within_tolerance -- --nocapture
MICHELL_BENCH_SAMPLES="$samples" cargo bench -p michell --bench wigley

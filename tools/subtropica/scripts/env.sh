#!/bin/bash
# SubTropica environment template — EDIT the paths for your install.
#   source tools/subtropica/scripts/env.sh
# The engine (HyperFLINT CLI port) is built separately; see GUIDE.md.
export SUBTROPICA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Run-dir root for all driver entry points (default: ./subtropica_runs under cwd)
export SUBTROPICA_RUNS=${SUBTROPICA_RUNS:-"$PWD/subtropica_runs"}
# Engine wrapper + MZV data table (REQUIRED: point these at your engine build)
export SUBTROPICA_HF_BIN=${SUBTROPICA_HF_BIN:-hyperflint.sh}
export SUBTROPICA_MZV_DATA=${SUBTROPICA_MZV_DATA:-mzv_reductions.json}
# MZV table pin (R7): set to `sha256sum "$SUBTROPICA_MZV_DATA"` of YOUR table.
# When set, mzv_provenance_gate refuses on mismatch — any table change
# re-triggers the gate against the new hash. Empty = unpinned (recorded).
export SUBTROPICA_MZV_SHA256=${SUBTROPICA_MZV_SHA256:-}
# Required memory cap: every engine invocation runs under ulimit -v 32505856.
# Period atoms: keep HF_PERIOD_TUPLES=0 (named MZV symbols, not opaque ids).
export HF_PERIOD_TUPLES=0

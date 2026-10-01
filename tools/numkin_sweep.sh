#!/bin/bash
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
# Flat-path alias of tools/numkin/numkin_sweep.sh: an exec wrapper, so argv, exit
# code and signal disposition are preserved (no subshell lingers).
_T="$(dirname "$0")/numkin/numkin_sweep.sh"
if [ ! -f "$_T" ]; then
  echo "numkin_sweep.sh wrapper: MISSING target '$_T' — expected the real tool at tools/numkin/numkin_sweep.sh" >&2
  exit 127
fi
exec bash "$_T" "$@"

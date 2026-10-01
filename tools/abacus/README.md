# ABACUS — certified counting on abelian fourfolds

The canonical body lives under `build/abcount/`; the harnesses resolve their sha-pinned files
relative to the package root, and the `abcount_*` names are the import
surface.

## Layout (why build/abcount/ is nested here)

The harnesses resolve their pinned files relative to the package root,
`PKG_ROOT = dirname(dirname(HERE))`, with paths like
`build/abcount/abcount_s0.py`. This directory reproduces that geometry so
the harnesses run UNMODIFIED:

    tools/abacus/                      <- the package root (PKG_ROOT)
      M3T_REGISTRATION.md              (sha-pinned spec doc)
      M3T_VALIDATION_BATTERY.md        (sha-pinned battery spec)
      S0PRE_RECEIPT.md                 (sha-pinned input-contract pins)
      SHA_PINS_M3T.txt                 (sha pin table)
      build/abcount/                   <- the CANONICAL tool body
        abcount_s0..s3.py, hecke_desk.py, run_s0..s3.py,
        theta_g4_bridge.jl, period_g2_bridge.jl,
        battery/SEEDS.json (BATTERY_RESULT.json is written here on each
        battery run), manuals/abcount.md
      abcount_s0..s3.py, hecke_desk.py <- COMPAT SHIMS (header-marked)
      run_s0..s3.py                    <- COMPAT SHIMS (header-marked)

## Call forms (key-by-name: every abcount form shim-forwards)

- Canonical battery/stage run (from the repo root):
  `python3 tools/abacus/build/abcount/run_s3.py`
- Compat CLI (same thing):
  `python3 tools/abacus/run_s3.py`
- Compat import (with tools/abacus on sys.path): `import abcount_s0` etc. —
  the top-level shims load the canonical body under the same module names.

## Acceptance rule

Outputs are EVIDENCE only behind the shipped battery (B1-B4 + B4s + B5 + N1-N5); ANY tool change re-runs the FULL battery
(run_s3.py) before its outputs are used.

Manual + math/proof authority: build/abcount/manuals/abcount.md (beside
the code).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

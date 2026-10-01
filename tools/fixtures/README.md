# fixtures — shared test data

Small data files consumed by other packages' batteries:

- `amflow_smoke/` — job specs and reference output for the amflow-cpp smoke
  test (used by `tools/amflow-kit`).
- `gatekeeper/amflow_ball_nonzero.json` — an amflow-cpp `solve_integrals`
  ball-string dump used by `tools/gatekeeper` (a nonzero ball-form midpoint
  must not be flagged as a suspicious zero).

Nothing here is executable; the consuming packages document how each file is
read.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

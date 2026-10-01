# reference/ — vendored upstream sources (READ-ONLY)

Vendored for spec citation stability. All files MIT-licensed
(github.com/SubTropica/SubTropica) under two copyright lines: `SubTropica.wl`
and `LICENSE` are Copyright (c) 2025-2026 Sebastian Mizera, Mathieu Giroux,
Giulio Salvatori (`reference/LICENSE`); the HyperFLINT files
(`hyperflint_cli_main.cpp`, `hyperflint_bridge_handlers.cpp`,
`HyperFLINT_README.md`) are Copyright (c) 2026 Sebastian Mizera and SubTropica
contributors (`reference/LICENSE-HyperFLINT`, upstream `HyperFLINT/LICENSE`).
Verify with `sha256sum -c SHA256SUMS` in this directory. Do NOT edit —
re-vendor from upstream and update SHA256SUMS + this file instead.

## Upstream pins — TWO SHAs, note the skew

| Artifact | Upstream SHA | What |
|---|---|---|
| `.wl` source clone | `adac2f722be64337aa095b2b2e7266628b03289b` | upstream commit (git log verified). All vendored files here are byte copies at THIS commit. |
| engine binary the port was validated against (not distributed here) | `ead8c6e` | upstream commit of the `hyperflint` v1.2.8 release binary (`bin/VERSION`; binary sha256 `5b9893238520fab5667d6865cce33c9d5f0d7568533d371dcf6bddf5b8dd34a8`) |

**SKEW**: the upstream release binary this port was validated against (not
distributed here) is built from git `ead8c6e`; the `.wl`/C++ source
pin is `adac2f72` (the v1.2.8 release commit — adjacent, NOT identical).
Consequences:
- Spec citations (CONTRACTS.md, DESIGN.md) refer to the vendored copies here,
  which are at `adac2f72`. The binary's behavior may differ in `ead8c6e`-only
  commits; any behavior/spec mismatch found during implementation must be
  resolved by LIVE probe of the binary, and the probe result recorded in
  CONTRACTS.md.
- The version gate is therefore stamp + behavioral canary, never stamp alone
  (see DESIGN.md correction vi).
- Any future upstream pull (either artifact) re-runs the FULL fixture suite
  (`test/runtests.jl`) before the new pin is accepted.

## File map (vendored name → upstream path @ adac2f72)

| Vendored file | Upstream path |
|---|---|
| `SubTropica.wl` | `SubTropica.wl` (35,606 LOC front-end; all `SubTropica.wl:NNNN` citations in CONTRACTS.md/DESIGN.md are line numbers in THIS file) |
| `hyperflint_cli_main.cpp` | `HyperFLINT/bridge/cli/main.cpp` (CLI op dispatch, ~69 ops; per-op request schemas inline as `// Request:` comments) |
| `hyperflint_bridge_handlers.cpp` | `HyperFLINT/src/bridge/handlers.cpp` (transport-neutral op handlers; e.g. `verify_order` FIELD of find_lr_orders at :741) |
| `HyperFLINT_README.md` | `HyperFLINT/README.md` |
| `LICENSE` | `LICENSE` (MIT; SubTropica) |
| `LICENSE-HyperFLINT` | `HyperFLINT/LICENSE` (MIT; governs the three HyperFLINT files above) |

## License pointer

MIT. Full texts in `reference/LICENSE` (SubTropica) and
`reference/LICENSE-HyperFLINT` (HyperFLINT). Redistribution of these copies in
this repository is within the MIT grant; keep this notice and both LICENSE
files next to the vendored sources.

## Modifications

The two C++ files (`hyperflint_bridge_handlers.cpp`, `hyperflint_cli_main.cpp`)
are modified copies: engineering review-note comments inside them were
neutrally reworded (no code bytes changed) and `SHA256SUMS` pins the modified
files. They still carry upstream's own dated development comments. All other
files are byte-identical to upstream at the commit above.

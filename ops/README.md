# ops/ — operations packages

Infrastructure for *running* BootLoops on a shared machine, kept apart from the
scientific instruments under [`../tools/`](../tools/README.md). Nothing under
`tools/` imports from here; nothing here computes a scientific result. Each
package carries a `GUIDE.md` (purpose, when to reach for it, gates, hazards), a
`README.md`, and its own self-test, in the same format as the `tools/` packages.

`../run_selftests.py` walks `tools/` only; an ops package's battery is run by
hand with the command in the table (from the repository root). The verification
classes mean the same as in `tools/README.md`.

| package | what | class | battery (from repo root) |
|---|---|---|---|
| `turnstile/` | Turnstile: admission control for long jobs on a shared machine — priority token, RAM and CPU-width ledgers (Linux-only; named SKIP on darwin) | selftest | `bash ops/turnstile/selftest.sh` |

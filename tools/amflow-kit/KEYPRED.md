# amflow-kit member manual: keypred — offline IBP-cache key predictor (yaml-byte canonicalized)

`amflow_kit.keypred` (`tools/amflow-kit/amflow_kit/keypred.py`), a member of
amflow-kit, the house kit around the AMFlow.cpp fork (overview: `GUIDE.md`).
Tool page: https://bootloops.ai/tools/amflow.html (the kit rides the AMFlow entry; the
cache-key context is on the Kira page, https://bootloops.ai/tools/kira-stack.html)

KIND: module + command line (stdlib only)

PURPOSE: Offline ibp-cache key predictor that hashes the CANONICAL jobs.yaml form — the
form amflow-cpp's own ibp_cache_key hashes — so cache injections/replica predictions land
on the keys the live run actually looks up. Kills the class where a raw-byte prediction
from one d-cap byte convention misses live keys computed from the other.

USE-WHEN:
- Predicting the cache key for an injection BEFORE the run that will look it up starts.
- Forensics on a key MISS: `key --raw` reproduces the defective raw-byte convention; diffing
  `key` vs `key --raw` on a dir proves/disproves the byte-convention mechanism in one command.
- Keying a hand-built mirror dir: `predict` repairs CRLF/missing-final-newline onto the
  kira_yaml.cpp writer convention first.
- Predicting a numeric-d Reduce leg's key: `--numd=P/Q` salts with the exact live string
  `AMFLOW_DIFFEQ_NUMD=P/Q` (bare `--numd` reads the env var, as the live run does);
  `--salt=S` is the generic verbatim door. The salt is folded after the input files under
  its own `\0@salt\0=` frame, mirroring the C++ extra-salt fold.

NOT-FOR: mirror SHAPE drift — a staged dir whose semantic content differs from what amflow
writes (wrong sectors/r box, run_firefly vs triangular mode, different target set) still
misses; that class is handled at run time instead (read the live MISS key from
the log, re-point the entry). `--explain` prints the canonical bytes so shape drift is
inspectable.

INVOKE (from `tools/amflow-kit`): `python3 -m amflow_kit.keypred key DIR [--raw] [--explain]
[--salt=S | --numd[=P/Q]]`; `python3 -m amflow_kit.keypred predict DIR [--explain] [--salt=S |
--numd[=P/Q]]`; `python3 -m amflow_kit.keypred --selftest` (<1 s); `python3 amflow_kit/keypred.py
ARGS` is the same command. Library: `from amflow_kit import keypred`;
`keypred.ibp_cache_key(DIR, raw=False, pin=False, explain=None, salt="")`,
`keypred.canonicalize_jobs_yaml(text)`, `keypred.pin_jobs_yaml(text)`, `keypred.numd_salt("P/Q")`,
`keypred.validate_numd(v)`.

INPUTS: a kira input dir (config/{integralfamilies,kinematics}.yaml, jobs.yaml,
jobs_kira2math.yaml, preferred, target — missing files skipped exactly as the C++ does).

OUTPUTS: 16-hex FNV-1a key on stdout; `--explain` adds pin repairs + the canonical
jobs.yaml bytes as hashed to stderr.

GATES: `--selftest` 16/16 — 5 fixture positive controls (`tests/fixtures/keypred/` = synthetic toy
input dirs in the kira format, one per byte convention, planted so the key relations
are known in advance: the d-cap pair and the n_thread dir share ONE canonical key
while their raw keys differ — the miss mechanism — and the two Masters-mode dirs
differ only in the result-determining `select_mandatory_recursively` d:, so their
canonical keys must differ) + 4 mutation predicts (crlf, no-final-newline,
dcap-stripped, thread-line all recover the base key) + 6 salt legs (pinned salted
key; empty salt byte-identical to unsalted; salted != unsalted — the numeric-d/
symbolic collision class; distinct numd values key distinctly; the salted trio still
shares one canonical key; salted predict on a crlf mirror recovers the pin) + the
numd format gate (valid P and P/Q accepted, malformed/overflow/Q<=0 refused). The
pinned hex constants are regression pins generated from these fixtures by the
replica itself; the FNV-1a core is anchored to its published public test vector.
The convention derives from ibp_cache.cpp / kira_yaml.cpp / reduce.cpp /
numd_subst.cpp, which ship in the AMFlow.cpp fork (the sibling repository amflow-cpp).
Run after any edit and after any amflow-cpp ibp_cache.cpp/kira_yaml.cpp/
reduce.cpp/numd_subst.cpp change.

FOOTGUNS:
- The canonicalizer must track ibp_cache.cpp BYTE-FOR-BYTE — never "improve" it (e.g.
  whitespace normalization) inside `canonicalize_jobs_yaml`; mirror-convention repair
  belongs ONLY in the `pin_jobs_yaml` layer. If the C++ canonicalizer changes, this
  replica and its fixtures must be re-validated in the same cycle.
- The d: inside select_mandatory_recursively (Masters mode, after select_integrals) is
  RESULT-determining and is KEPT by both the C++ and this replica — do not strip it when
  hand-reasoning about keys (AMFLOW_PREHEAT_DOT changes that byte AND the result, hence
  legitimately the key).
- The cache key also hashes kinematics.yaml (per-point) — a warm cache still cold-MISSes
  at a new kinematic point unless AMFLOW_SYMBOLIC_IBP=1 (see the amflow wrapper env rows in
  the AMFlow.cpp fork, the sibling repository amflow-cpp).
- Raw mode exists for forensics ONLY; never build an injection off a `--raw` key.
- Numeric-d salting applies to Reduce-mode FireFly legs ONLY. Masters-mode legs keep the
  UNSALTED key (their result is d-independent, so numeric-d and symbolic runs deliberately
  share preheat entries) — never salt a Masters prediction. The salt must be the FULL
  `AMFLOW_DIFFEQ_NUMD=<value>` string; `--numd` builds exactly that, `--salt=` folds its
  value verbatim (a bare value without the env-assignment prefix predicts a key no live
  run computes).

RECEIPTS: the canonicalization convention and its upstream pins, below.

RELATED: the sibling repository amflow-cpp (the amflow-cpp fork whose ibp-cache convention this
replicates); tools/kira-stack (the kira extension scripts); the kit's other members (GUIDE.md).

## The canonicalization convention

Convention derived from (upstream pins, the AMFlow.cpp fork in the sibling repository amflow-cpp):

- `src/ibp/ibp_cache.cpp` — kInputFiles order, `\0@name\0=` framing, FNV-1a 64,
  `canonicalize_jobs_yaml` (seed-line `, d: <int>` strip before
  `select_integrals`; `n_thread`/`n_threads`/`threads` line drop;
  `select_mandatory_recursively` `d:` KEPT).
- `src/ibp/kira_yaml.cpp:599-611` — `AMFLOW_REDUCE_DCAP` appends `, d: N` to
  the Reduce seed line AT WRITE TIME (Reduce mode only); the key is computed
  from the on-disk post-write bytes, i.e. the canonicalizer at hash time is
  what makes both byte conventions key identically.
- `src/ibp/ibp_cache.cpp` extra-salt fold — a non-empty `extra_salt` is hashed
  AFTER the input files under its own 8-byte `\0@salt\0=` frame; empty salt
  folds nothing, so unsalted keys are byte-identical to the salt-free
  convention.
- `src/ibp/reduce.cpp` — numeric-d Reduce legs (`AMFLOW_DIFFEQ_NUMD` set,
  FireFly mode) pass `AMFLOW_DIFFEQ_NUMD=<value>` as the salt: the d -> P/Q
  substitution rewrites tmp/SYSTEM files that are NOT hashed inputs, so
  without the salt the numeric-d and symbolic entries would collide in both
  directions. Masters-mode legs stay unsalted (d-independent result; preheat
  entries deliberately shared).
- `src/ibp/numd_subst.cpp` validate_numd — the value is `P` or `P/Q`, P a
  signed and Q an unsigned decimal integer, each within signed 64-bit, Q > 0.

Mechanism: offline predictions that hash RAW `jobs.yaml` bytes of one d-cap
convention miss live keys computed from the other; the canonical form predicts
the live key from EITHER d-cap byte convention.

Scope: byte-convention misses only. Mirror SHAPE drift (wrong sectors/r box,
wrong reduce mode, wrong target set) still misses and stays covered by
re-reading the live run's inputs.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

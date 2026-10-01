# ffcapital — guide (THE FireFly ff_save salvage package)

Tool page: see the tools index, https://bootloops.ai/tools/ (ffcapital row
in the "Reduction and IBP" group).

KIND: package (members: harvest_ffsave_eta.py, harvest_ffsave_2var.py,
recon_symbolic_d.py)

PURPOSE: Exact F(d,x) out of dead/capped FireFly state, fail-closed. Every member
decodes kira+FireFly `ff_save/` reconstruction state (RatReconst::save_state) with
ZERO solver compute and banks only what passes a typed gate battery — a fn that closed
before the node died carries its complete rational function in `states/`, and
`validation.gz` black-box probes give an independent per-fn oracle. The 1-var, 2-var, and cross-slice symbolic-d
salvage tools are MEMBERS of this one package, not standalone tools.

USE-WHEN (route by save shape):
- 1-var (numeric-d slice; farm node killed at cap) → **harvest_ffsave_eta** member.
- 2-var symbolic (d,eta) save dead mid-CRT, slot names wanted → **harvest_ffsave_2var**
  member.
- numeric-d slice FARM harvested, symbolic-d dependence wanted back →
  **recon_symbolic_d** member.

NOT-FOR: advancing/resuming the solve (that is FireFly's job — see the amflow wrapper
rows in the AMFlow.cpp fork, the sibling repository amflow-cpp, for resumable
embedded saves); >2-var states (extend the monomial parser first); per-entry degree censuses (ffsave_degree_census.py, same state format — a
GPL-3.0-or-later tool derived from Kira/FireFly sources, shipped in the sibling repository
kira under `bootloops-tools/`; `--census` here takes its JSON output, and
`tools/kira-stack/kira_gpl_tools.py` locates it: `../kira/bootloops-tools`
or env `BOOTLOOPS_KIRA_TOOLS`).

INVOKE (public; each member is a CLI script):
- `harvest_ffsave_eta.py census --ffsave <ff_save-dir> --out census.jsonl [--workers 16]`
- `harvest_ffsave_eta.py extract --ffsave <dir> --out funcs.jsonl (--fns 0,17,26683 |
  --all-done) [--workers 16] [--no-validate]`
- `recon_symbolic_d.py control|index|probe3|run --outdir <dir> [--coverage ...]
  [--wave1 ...] [--bank-states ...] [--ff-helper <ReconstHelper.cpp>] [--p2-file ...]
  [--config <json>] [--base <root>] [--skip-banked <prior>]` — flag > --config key >
  derived-from---base; unset+needed = loud abort (`CONFIG: <NAME> unset — ... no
  silent defaults`). `control` = the mandatory synthetic-truth battery, run FIRST and
  after ANY pipeline patch.
- harvest_ffsave_2var: decode + validate per-save, ASSIGN each fn to its
  (target, master) slot by fingerprint + constraint discovery, certify by exact
  cross-multiplication vs ≥2 exact rational-d slice reductions; typed quarantine;
  >`--kill-rate` hard-failure rate = CONTAMINATED_STOP, no bank.
All members need a USER-SUPPLIED `ff_save/` state (and, for recon_symbolic_d,
FireFly's ReconstHelper.cpp as the prime-table source) — no public fixtures ship.

GATES (fail-closed, per member):
- eta: prime auto-detect UNANIMOUS and unique over ≤8 sample done fns (assert; a
  corruption INSIDE the sample aborts loud instead of voiding); validation.gz entry
  count == state-file count (assert); per-fn black-box match mod the detected prime;
  rc=2 + named fn list on any validation failure (voids the fn, never the node).
- recon_symbolic_d (ALL must pass to bank; any failure voids the fn → RETRY):
  (control FIRST) synthetic truth with planted d-vanishing coefficient + planted
  cancellation slice, exact recovery + negative controls each FAIL loudly; (fit)
  overdetermined exact nullspace (need_fit=nd+dd+5); (internal) G_aN·H_aD == 1 exactly
  + deg_d(P)≤nd, deg_d(Q)≤dd; (heldout) EVERY validation slice by cross-multiplication,
  single λ; (loo) every fit slice: identical refit without it; (oracle) never-in-fit
  multi-prime bank.
- Package-level battery (receipts not shipped): synthetic-ff_save truth case 10/10 exact
  recovery + 2 fail-closed negatives; parked-node replay vs its banked harvest
  receipt; control battery byte-identical to the banked control_result.json;
  cancellation-slice regression (naive coeff-vector REJECTS / cross-mult ACCEPTS /
  corrupted slice still REJECTED).

FOOTGUNS (package-wide):
- Harvest AT-REST only: verify no kira/FireFly pid has the save as cwd (pgrep +
  /proc/<pid>/cwd) BEFORE and re-verify AFTER any remote pull; live nodes rotate state
  files. Re-pulling a DEEPENED node needs `rsync --delete` (stale twins = typed
  abort/quarantine).
- A poisoned slice can PASS its own black-box validation (internally consistent, wrong
  for fitting) — node health is certified by the downstream cross-slice battery, never
  by extraction alone.
- Slice-local num/den cancellation yields a REDUCED but correct function: consumers
  compare by CROSS-MULTIPLICATION, never coefficient-vector diff.
- FireFly's per-slice normalization pins SOME coefficient to 1/1 (per-fn stable, NOT
  stable across runs) — fit only scalar-free anchor ratios; cancellation slices
  auto-EXCLUDED from fits (joint (degN,degD) drop ≥1 below generic), validation-only
  via cross-mult; *_HELDOUT slices + quarantined nodes forced out of every fit set
  (never_fit law); re-run `control` after ANY pipeline patch.
- recon_symbolic_d slice-name SCOPE is named, not silently generalized: slice names
  follow `node_eps_<a>_<b>[_HELDOUT]` (d=4−2a/b); the control battery's planted
  vanishing/cancellation d-points are picked FROM the run's own grid
  (`control_plants`) — the reference grid keeps its banked d=118/29 plant so the
  byte-identical control_result.json law still holds; a grid too small for the
  synthetic coverage law (< nd+dd+9 = 16 slices) or without ≥2 distinct d-values
  is refused loudly by name, never fitted thin. Coverage law (need_fit/need_total)
  comes from the upstream degree census, not from this tool.

BATTERY (harvest legs data-gated; control leg data-FREE): the ff_save-decoding
gates are fail-closed and need a user-supplied ff_save/ state — no public
fixtures. Smoke: `--help` on both harvest members (rc=0);
`recon_symbolic_d.py control` refuses loudly with the named unset CONFIG
(fail-closed refusal fires as designed). The synthetic-truth control itself
needs only a slice-name grid plus the FireFly prime table, so
`python3 selftest.py` runs it FOR REAL (generated prime-table fixture;
reference grid full-battery PASS, generic grid without d=118/29 PASS with an
on-grid plant, too-small grid refused by name; python-flint absent = named
skip), then exits nonzero with the named data gate for the harvest members. MEASURED (production runs): 26,719
fns extracted+validated in 28 s on one node; ~1.4M validated extractions across 55
slices, 0 failures; 2-var 5,125/5,125 union-done banked, 0 failures, 78 s @24 workers;
symbolic-d 24,924 fns banked (93.3%), 0 oracle failures.

RELATED: ffsave_degree_census (GPL, kira/bootloops-tools — degree census,
same state format);
numkin (builds the freeze-one-scale slice farms these salvage); ratfit (1-var
exact-fit-plus-gates analogue); the AMFlow.cpp fork
the sibling repository amflow-cpp (resume vs salvage).

CREDIT: decodes the save-state format of FireFly by Jonas Klappert and Fabian Lange
[FF1] and Klappert, Sven Yannick Klein and Lange [FF2] (the layout written by
RatReconst.cpp:2537-2716 and Reconstructor.hpp:2157ff, read as the file-format
specification; FireFly is GPL-3.0-or-later, obtained upstream, not bundled, and never
linked, imported or executed here; no FireFly code is included — see THIRD_PARTY.md
B2), as written by Kira [Kira2, Kira3]. We are grateful to the FireFly authors for a state
format regular enough to salvage. Bracketed keys resolve in REFERENCES.md at the
repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

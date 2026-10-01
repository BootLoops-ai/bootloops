# Ansatzer — guide

Ansatzer — `tools/ansatzer` — https://bootloops.ai/tools/ansatzer.html

KIND: package (`ansatzer.py`, `alphabets/` validated fixtures,
`test_ansatzer.py`). Per-orbit path: `--per-orbit`.

PURPOSE: predicts the residual symbol-ansatz dimension per weight (and per master, if a
connection is supplied) after the standard big structural cuts — BEFORE oracle-farming.
Residual = number of unknowns the value-fit must determine ≈ number of high-precision
oracle points needed (× safety 1.3). Know the fit cost before farming.

USE-WHEN:
- Residual fit size before farming oracle points (alphabet ± connection → residual-count
  table).
- Deciding whether a value-fit is affordable at all (rule of thumb: a value-fit is
  practical below ~100-1000 residual functions).
- Per-orbit residual K_w, finer than the global count → `--per-orbit` (recursive
  symbol-alphabet builder over the subsector inclusion lattice, mod-p convergence check
  at 12-24 generic points; hexabox: K_w4=1 sec255, K_w4=6 sub-orbits).

NOT-FOR: does not evaluate anything; coproduct cut is SKIPPED without a connection
(flagged in output).

INVOKE:
`python ansatzer.py --alphabet A.json [--connection C.json] --wmax 4
[--method auto|exact_Q|modp] [--no-steinmann] --out report.json
[--compare B.json]`
- Per-orbit: `python ansatzer.py --per-orbit --alphabet A.json [--wmax 4]
  [--npts 96] [--out r.json]` → JSON `{Kw_LE, Kw_full_noLE, converged, wallclock_s}`
  (adaptive 6→npts mod-p; same input adapter as the global path; channel tags
  not applied on this path).

INPUTS:
- alphabet.json (collapse-schema; landau-alphabet output auto-adapted):
  `{"name","vars","letters":{name:expr},"first_entry":[...],"last_entry":[...],
  "parity":{...},"parity_target":[...],"channel":{...}}`. Algebraic letters must be
  supplied on a RATIONALIZED chart so dlogs are rational in vars (e.g.
  O_A3=((1+y)/(1-y))² on t=4/(1+y²) reproduces radical-aware p4 counts exactly —
  `alphabets/p4_J25.json`). `channel` (tag list or bare tag string per letter)
  activates the Steinmann cut: adjacent entries with disjoint tag sets are
  forbidden (same-channel repeats allowed; untagged letters unconstrained);
  `--no-steinmann` disables it.
- conn.json (optional): `{"masters":[...],"A":{letter:[[i,j,"q"],...]}}` — only the
  SUPPORT (which (i,j) nonzero) is used.

OUTPUTS: report.json; columns
`weight | raw | integrability | first | last | steinmann | parity | coproduct |
RESIDUAL | pts_needed`. exact_Q mode: each column = subspace dimension after that
cut. modp mode: first/last/steinmann/parity/coproduct are word-prefilter counts;
integrability=RESIDUAL is the converged mod-p nullspace dim on surviving words.

GATES: honesty flags in output — `coproduct: SKIPPED (no connection supplied)`;
`steinmann: applied | vacuous (no channel tags) | SKIPPED (disabled)`;
`rank_method: exact_Q | modp | trivial` per weight; `bound: ≥N` = integrability timed
out at that weight ⇒ residual is an UPPER bound only. Validated fixtures in
`alphabets/`: trivial_2var (hand-check), p4_J23 (tower 3,7,15,31 → residual 2,4,8,16),
p4_J25 (exact match on rationalized chart), synth_7letter (synthetic 7-letter,
invented letters, w1 cuts and w2 Steinmann word count hand-checkable),
5pt_2mass_nonplanar (19-letter), p4_J23_conn (per-master coproduct demo). Fixture names
are physics family tags (synth_* = synthetic). Battery:
`python3 test_ansatzer.py` (T1-T4, T6 Steinmann, T7 per-orbit adapter
self-contained; T5 adapter leg SKIPs unless a landau-alphabet `alphabet.json`
is supplied via `COLLAPSE_T5_ALPHABET`).

FOOTGUNS:
- Algebraic letters on the raw chart give wrong counts — rationalize the chart first.
- A `bound: ≥N` residual is not a farm budget; re-run exact or per-orbit before pricing
  a farm.

RELATED: upstream — landau-alphabet (alphabet input; see the Landau alphabet extractor
tool page); downstream — lockpick (which points to farm), oracle farms.

CREDIT: the cuts are the standard symbol-bootstrap constraints: symbols and
integrability after Goncharov, Spradlin, Vergu and Volovich [GSVV] and Dixon, Drummond
and Henn [DDH]; the first-entry condition of Gaiotto, Maldacena, Sever and Vieira
[GMSV]; Steinmann relations as used by Caron-Huot, Dixon, McLeod and von Hippel [CDMvH]
and their extended form [CDDvHMP]; the coproduct cut in Duhr's formulation [Duhr12] of
Goncharov's coproduct [Gon]; and the sequential-discontinuity constraints of Bourjaily,
Hannesdottir, McLeod, Schwartz and Vergu [BHMSV] and Hannesdottir, McLeod, Schwartz and
Vergu [LB2, LB3]. The residual-count bookkeeping and the mod-p engine are ours.
Bracketed keys resolve in REFERENCES.md at the repository root.

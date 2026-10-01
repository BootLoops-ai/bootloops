# Dogtag — integral-family identity checks (isomorphism up to loop-momentum relabeling, catalog match, planarity, cut signature, drawn-graph compare; loomcheck applicability screen)

Tool page: https://bootloops.ai/tools/dogtag.html (listed in the tools index)

NAME: a dog tag identifies its wearer, and Dogtag is the tag that establishes an
integral family's identity — which family the definition file actually is, up to a
relabeling of the loop momenta, which catalog topology, planar or crossed in the physical
leg order, which cut signature, the same graph as its drawing or not — before any
reduction or numerical evaluation is spent on it. The package directory is `tools/dogtag`;
the script keeps its historical file name `topology_audit.py` (module `topology_audit`),
so every command below reads `python3 topology_audit.py ...` from this directory.
(Formerly topology-audit.)

KIND: script (sympy + networkx). Single file `topology_audit.py` + README + `tests/`
(the batteries and the vendored fixtures) + the member subfolder `loomcheck/` (the
Yangian/loom/fishnet applicability screen, section "Member: loomcheck" below);
catalog-driven, fully self-contained. networkx
is a declared dependency: without it the realizer / planarity / cut-signature /
isomorphism legs of the self-test and the pytest battery SKIP by name (`networkx
unavailable (declared dependency; pip install networkx)`) and the battery exits 2 (a skip
is an unknown, never a pass). This is the one statement of that contract; the battery
section below refers back to it.

## PURPOSE

Feynman-family provenance/topology auditor — catches the NP-dbox-class mislabel (a
"nonplanar crossed box" that is actually the planar Smirnov double box under k2 -> -k2).
Computes (a) the loop-relabeling isomorphism against a canonical catalog — signed
permutations k_i -> ±k_sigma(i) first, then the affine class k -> U k + c·p (U unimodular,
the shift c solved, never enumerated) that describes a second routing of the same graph —
with the line masses, the propagator multiplicities (dots) and the external-leg
virtuality classes carried as labels, and the realized-graph isomorphism (networkx VF2 on
the vertex-realized labeled multigraph) as the identity test; (b) PHYSICAL planarity via
the external-legs-closed-in-cyclic-order test (abstract-graph planarity cannot
distinguish planar vs crossed box — the closure cycle does), the closure order taken from
the record's declared cyclic order or the canonical index order by name, flagged
"not discriminating (N < 4 legs)" where the closure cannot separate anything, with the
legs-joined-at-infinity criterion reported beside it; (c) per-channel (s,t,u) unitarity
min-cut signature as a label-free fingerprint, plus a routing-independent canonical hash
of the realized labeled graph beside the routed one; (d) the compare mode: a family of
record against its drawn graph gives exactly one of `IDENTITY-PASS` / `FAMILY-MISMATCH` /
`NON-GRAPH` / `NOT-CHECKABLE` with its evidence (the identity census's verdict strings).

## USE-WHEN

- Ingesting ANY new family file (integralfamilies.yaml with its kinematics.yaml,
  AmflowFamily *.jl, an AMFlow-port JSON config, a pySecDec graph script) — one-line
  provenance check before reductions are launched on a mislabeled graph.
- A family's claimed planarity/provenance disagrees with its analytic behavior (e.g.
  Im=0 where a nonplanar signature was expected).
- A figure's graph must be checked against the family file that claims to be it: the
  compare mode (`--drawn`), on an edge list, a graphspec JSON or a routed yaml whose
  header comments state the legs, their p^2, the cyclic order and the multiplicities.

## NOT-FOR

The planarity closure test discriminates only with >= 4 external legs; with 2 or 3 the
field reads "not discriminating (N < 4 legs)" by name and the legs-joined-at-infinity
criterion is reported beside it — never a bare planar=True. A family whose leg set is not
declared (no kinematics, no leg keys) is read with the LEG-SET-INFERRED warning and no
planarity verdict. Very large families can exhaust the realizer's budget and degrade to
planar=None with the cause named (honest refusal, never a false flag). Catalog-match
warnings only cover families in its canonical-representative catalog (9 entries, the
README table).

## INVOKE

`python3 tools/dogtag/topology_audit.py FAMILY_FILE [--name NAME] [--leg-perms]
[--kinematics K.yaml] [--integral N] [--json]` = the family audit; `FAMILY_FILE --drawn
DRAWN_GRAPH [--name NAME] [--kinematics K.yaml] [--kinematics-drawn KD.yaml] [--integral N]
[--verdict] [--json]` = the compare mode (the family of record against its drawn graph:
`VERDICT: IDENTITY-PASS | FAMILY-MISMATCH | NON-GRAPH | NOT-CHECKABLE` first, then the leg
map, the identified masses, the relabeling, the cyclic-order sub-finding, legs per vertex
and V, E of both sides, planarity in the drawn order; `--verdict` = that and no audit;
`--name NAME` = the family of that name in a multi-family file is the side compared -- the
compare battery's route -- and a name the file does not carry is refused by name, rc 2, the
file's family names listed, no verdict); `--selftest`
(alias `--self-test`) = the validation suite of 17 legs, `SELF_TEST_LEGS` +
`EXTRA_SELF_TEST_LEGS` in the file ([1] NP-dbox; [2] catalog planar Smirnov dbox; [3] catalog
crossed box entries; [4] C3 dbox; [5] 3-loop banana; [6] catalog integrity; [7]
signed-relabeling positive control; [8] reader battery; [realizer] vertex realizer;
[catalog] catalog; [kinematics] kinematics; [legperms] leg permutations; [labels] labels;
[planarity] planarity; [isomorphism] isomorphism beyond signed permutations; [compare]
compare mode; [loomcheck] the loomcheck member's structural battery); the run prints the
executed / skipped / failed counts and the registered sub-leg totals. The member's own
entry points: `python3 -m loomcheck --selftest | --basis FILE ...` from this package
directory (section "Member: loomcheck").

`--integral N` (AMFlow-port JSON): which integral of the record defines the family is the
record's call — pass `--integral N` (the compare battery's row-21 entry uses index 2); the
reader's default is the integral with the most positive indices (first among ties). When
the record lists several integrals, the compare's evidence names the index it read with
its index vector and lists the others, in the printed evidence and in the JSON verdict
block (`integral_index`, `n_integrals`, `integral_indices`, `other_integrals` on each
side); a one-integral record prints nothing new. Row 21 by default: `integral 0 of 3: [1, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]; others: 1 [1, 1, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 2 [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]  (reader default: the integral with the most positive indices (first among ties))`;
with `--integral 2`: `integral 2 of 3: [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]; others: 0 [1, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1 [1, 1, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]  (--integral N (the record's call))`.

INPUTS: integralfamilies.yaml (+ kinematics.yaml, or the yaml's own external_momenta /
momentum_conservation / leg_virtualities keys) or AmflowFamily .jl, an AMFlow-port JSON
config (family block + integrals[].indices + numeric_values), a pySecDec graph script, a
drawn-graph edge list / graphspec JSON / routed yaml with the census's header comments.
OUTPUTS: fingerprint (routed canonical hash + the routing-independent realized hash, loop /
line counts, mass and multiplicity multisets, leg classes, cut signature, planarity with
its closure order and the discriminating flag) + warnings ("graph-isomorphic to <family>
under <relabeling>" with the labels carried, "LABEL/PROVENANCE MISMATCH: labeled nonplanar
but planar", "LEG-SET-INFERRED"); with --drawn the compare verdict and its evidence (JSON
key "verdict").

## GATES

--self-test before trusting on a new install: `SELF-TEST: ALL PASS (rc 0)` is the only
passing line — printed only at 0 skipped AND 0 expected-fail; rc 2 on any SKIP (networkx
absent, a fixture absent: named), rc 1 on any FAIL.

## FOOTGUNS

- q² == (−q)² makes sign-flipped relabelings invisible to propagator-set eyeballing —
  that is exactly the class this automates; don't "verify by inspection" instead.
- Leg INDEX order is a convention, not a cyclic order: the census's rule is that a leg
  map which does not preserve the drawn cyclic order is a named sub-finding
  (`CYCLIC-ORDER MISMATCH`) that never changes the verdict.
- `--integral N` is the record's call (see INVOKE); the default reading of a record
  listing several integrals can be the dotted one, and the evidence says which.
- N < 4 legs: the closure test is not discriminating; read the flag, not the boolean.
- The printed leg map is VF2's pick among the valid maps of a symmetric graph (the double
  box has four); the identity is the isomorphism, the map one witness.

## Member: loomcheck

The subfolder `loomcheck/` (one file of logic, `loomcheck/loomcheck.py`, importable as the
module `loomcheck` with this package directory on `sys.path`) is a second, independent
graph screen that shares the package's networkx dependency and its self-test.

PURPOSE: a mechanical applicability screen for the Yangian / loom / fishnet integrability
theorems on position-space conformal integrals. A graph is given as a dictionary of
x_ij^2 powers (a numerator factor is a negative power), or parsed from a Mathematica basis
file in the arXiv:2607.11645 ancillary format. Three published conditions are checked per
graph: (1) the conformal weight at every internal vertex, sum of powers = D; (2) planarity
of the FULL graph with the numerator edges counted as edges (every Yangian theorem needs
it: 1708.00007, 2304.04654, 2505.05550); (3) the level-one-momentum face condition of
2505.05550, sum of powers around each k-gon face = (k-2)D/2. `yangian_class` is true only
when all three hold.

USE-WHEN: deciding whether integrability (Yangian differential equations, separation of
variables) is an available route for a conformal integral before spending design time on
it; re-screening the default targets if a successor theorem weakens the conditions.

NOT-FOR (honest limits): SCREEN-ONLY — a FAIL proves exclusion from the cited theorems
as published; a PASS only fails to exclude and constructs nothing. The face check counts
the planar embedding's outer face like any other face: conservative for the nine default
targets (margin of at least 7 violations per graph), but a near-threshold verdict (0-2
violations) or a graph with tree-like appendages needs outer-face exclusion added before
it is trusted, and no fishnet positive control ships for the same reason. The basis-file
parser assumes the 2607.11645 conventions (single-digit vertex labels, entry =
numerator/(denominators), internal vertices labeled 5 and up); for other sources pass the
power dictionary to `screen_graph()` directly.

INVOKE (cwd = this package directory): `python3 -m loomcheck --basis 4Loop_int_basis.m
[--targets 144,163,...] [--D 4] [--json RECEIPT.json]` (the basis file is the
arXiv:2607.11645 ancillary; `LOOMCHECK_BASIS` sets the default; default targets = entries
144, 163, 167, 168, 169, 212, 214, 233, 260), or the script form `python3
loomcheck/loomcheck.py ...` from any cwd. Battery: `python3 -m loomcheck --selftest`, exit
0/1, ten structural known-answer checks that need no data files (K5 with a numerator edge
is non-planar, the four-point star's internal weight with unit and squared powers,
quadrilateral versus pentagon face sums). API: `from loomcheck import screen_graph;
screen_graph({(1,5):1, ...}, D=4, internal=[5,6])` returns `conformal_internal_weight`,
`planar_full_graph`, `numerator_edges`, `faces`, `phat_face_violations`, `yangian_class`.

GATES: the structural battery is the validation story, one hand-checkable graph per
condition; it runs three ways with the same checks — `python3 -m loomcheck --selftest`,
the package self-test's `[loomcheck]` leg (which adds the conjunction invariant on the
same graphs, 14 checks), and `tests/test_loomcheck.py` (those plus the parser on a
synthetic basis file and the CLI surface). The nine-target verdict (0/9 in any proven
class: four non-planar, the planar five fail the face condition) reproduces in one command
from the public ancillary: `python3 -m loomcheck --basis 4Loop_int_basis.m --json
RECEIPT.json`.

## Self-test and battery

- `python3 topology_audit.py --selftest` (cwd = this package or a scratch cwd, env -i)
  runs every registered leg and prints each as PASS / SKIP (named) / FAIL with a count
  line `SELF-TEST: N executed / M skipped (...) of K legs` and the registered sub-leg
  totals. Contract: `SELF-TEST: ALL PASS` is printed only with 0 skipped AND 0 open
  expected-fail sub-legs; the exit code is 2 whenever any leg was skipped (a skip is an
  unknown), 1 on any FAIL, else 0. The networkx SKIP contract is the KIND paragraph at
  the top of this file.
- Expected output from the package cwd:
  `SELF-TEST: 17 executed / 0 skipped (17 PASS, 0 FAIL) of 17 legs` /
  `SELF-TEST: registered sub-legs N executed / 0 skipped / 0 expected-fail` /
  `SELF-TEST: ALL PASS (rc 0)`; the `[compare]` leg's own line reads
  `[compare] battery: 53 cases, 53 PASS, 0 FAIL, 0 SKIP  ->  PASS` and the `[loomcheck]`
  leg's `[loomcheck] battery: 14 checks, 14 PASS, 0 FAIL  ->  PASS`. The pytest suite
  (`python3 -m pytest tests`) covers the same legs case by case plus the CLI surface
  (`tests/test_loomcheck.py` is the member's file; `python3 tests/test_loomcheck.py` runs
  it without pytest).
  The fixture legs run on vendored copies (`DOGTAG_NPDBOX_FAMILY` /
  `DOGTAG_C3_FAMILY` / `DOGTAG_BANANA_FAMILY` override them with your own files; the
  older names `TOPOLOGY_AUDIT_NPDBOX_FAMILY` / `TOPOLOGY_AUDIT_C3_FAMILY` /
  `TOPOLOGY_AUDIT_BANANA_FAMILY` are still read when the `DOGTAG_*` one is unset).
- Indicative timing on a workstation: the self-test about 2 minutes, the full pytest suite
  about 5 minutes, the compare battery alone (`python3 tests/census_battery.py`) about
  20 seconds; a single family audit or compare is about 1 second.
- Fixtures under `tests/fixtures/<row>/` are sha-pinned copies of the identity-census
  corpus files (30 published families and the graphs the paper draws for them; the row
  numbers are the corpus's). `PROVENANCE.json` per row pins each copy's sha256 and size
  (every battery refuses on drift) beside the original file's digest; the two differ
  where a comment or an unread configuration line was rewritten in the copy, never in a
  field a reader consumes. Fixture directory names are corpus identifiers: `cured_*` is a
  corrected family, `retired_*` the family it replaced. `tests/fixtures/CENSUS_BATTERY.json`
  is the compare battery's manifest (fixtures and expected verdict per entry, rollups,
  verdict counts).
- Maintenance: after an intentional edit of a fixture, `python3 tests/repin_fixtures.py`
  rewrites its pins in PROVENANCE.json, in the manifest and in the PINS / FIX tables of
  the test files and the tool (`--check` reports stale pins without writing); then re-run
  the self-test and the pytest suite.

CREDIT: realized-graph isomorphism uses the VF2 algorithm of Cordella, Foggia, Sansone
and Vento [VF2] through networkx (Hagberg, Schult, Swart); equivalent-family detection
is the problem addressed by Pak's algorithm [Pak] and by Vitaly Magerya's feynson, to
which this auditor is complementary (it reports the explicit unimodular map and physical
planarity). Catalog families are named for their authors: Smirnov's planar double box
[Smirnov99], Tausk's non-planar double box [Tausk], the Ussyukina–Davydychev ladders
[UD] (three-loop routing after Broadhurst and Davydychev [BD10]). Bracketed keys resolve
in REFERENCES.md at the repository root.

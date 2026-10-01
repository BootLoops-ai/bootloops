# Dogtag — integral-family identity checks

Dogtag establishes an integral family's identity from its definition file before
compute is spent on it (the package directory is `tools/dogtag`; the script keeps its
historical file name `topology_audit.py`, formerly the topology-audit package). It is a
one-line check that catches the **NP-dbox-style mislabel**: a family labeled
"nonplanar crossed box" that is actually the **planar Smirnov double box**
relabeled (under `k2 -> -k2` its propagator *set* maps onto Smirnov's; `Im = 0` at `u > 0` confirms the planar signature). Dogtag makes catching this class of mislabel automatic.

## What it computes

Given a family file (`integralfamilies.yaml` **or** an `AmflowFamily` `*.jl`):

1. **Loop-relabeling isomorphism** — searches every signed permutation of the
   loop momenta `k_i -> +-k_{sigma(i)}` (and, with `--leg-perms`, external-leg
   permutations). Because `q^2 == (-q)^2`, two "different" families that are the
   *same graph* under such a relabeling have identical canonical propagator
   sets. Matches against a catalog of canonical representatives and reports the
   relabeling (e.g. `k2 -> -k2`).
2. **Physical planarity** — reconstructs the Feynman graph from the
   physical-propagator momenta by momentum-conservation vertex realization, then
   tests planarity **with the external legs closed in their canonical cyclic
   order** (`p1->p2->...->pN`). This is the physical question (a planar box and a
   crossed box are the *same abstract graph* — only the leg order differs), so
   abstract-graph planarity alone does **not** discriminate them; the closure
   cycle does. `planar=True` <=> drawable with legs in canonical order (planar
   box); `planar=False` <=> forces a crossing (crossed box).
3. **Cut signature** — per Mandelstam channel `(s,t,u)`, the minimal number of
   internal propagators a unitarity cut must sever (graph min-cut between the
   two external clusters). The `(s,t,u)` multiset is a label-free fingerprint.

Output: a **fingerprint** (canonical hash, loop/prop counts, cut signature,
planarity) plus **warnings** — `graph-isomorphic to <known family> under
<relabeling>` and `LABEL/PROVENANCE MISMATCH: labeled nonplanar but planar`.

## Usage

```bash
python3 topology_audit.py FAMILY_FILE [--name NAME] [--leg-perms] [--json]
python3 topology_audit.py FAMILY_FILE --drawn DRAWN_GRAPH [--name NAME] [--kinematics K.yaml] \
        [--kinematics-drawn KD.yaml] [--integral N] [--verdict] [--json]
                                               # compare mode: the verdict line
                                               # first, then its evidence; --name
                                               # picks the family of a multi-family
                                               # file (unknown name: refused, rc 2)
python3 topology_audit.py --self-test          # validation suite: ALL PASS (rc 0) only
                                               # at 0 skipped AND 0 expected-fail;
                                               # rc 1 any FAIL, rc 2 any SKIP (named)
```

**Compare mode** (`--drawn`): the family of record against a drawn graph (an
edge list, a graphspec JSON, or a routed yaml whose header comments state the
legs, their p^2, the cyclic order and the multiplicities) emits exactly one of
four verdict strings, `VERDICT: IDENTITY-PASS` /
`FAMILY-MISMATCH` / `NON-GRAPH` / `NOT-CHECKABLE`, with the evidence beneath
it: the leg map, the identified mass labels (a drawing's `m2` against a
record's resolved `1`), the relabeling (a signed loop permutation, or the
affine loop redefinition when the two propagator sets differ by a routing),
the cyclic-order sub-finding (`CYCLIC-ORDER MISMATCH` when both sides declare
an order the leg map does not carry up to rotation and reflection -- reported,
never changing the verdict), the legs per vertex and V, E
of both sides, and the planarity of both graphs with the legs closed in the
drawn order. `--verdict` prints the compare alone; `--json` carries it under
`"verdict"` (and the family's audit under `"audit"` unless `--verdict`).
IDENTITY-PASS is the realized-graph isomorphism with masses, multiplicities
and leg virtuality classes; FAMILY-MISMATCH names the first failing level
(bare / masses / legs); NON-GRAPH is the realizer's verdict on either side;
NOT-CHECKABLE names the side that cannot be read, whose legs do not close, or
that declares no virtuality class for a leg (three or more legs; with two the
class level is vacuous and said so).

`FAMILY_FILE` may be a Kira `integralfamilies.yaml` (uses `top_level_sectors`
to pick the physical propagators) or an `AmflowFamily` `*.jl` (uses the
`*_INDICES = [...]` vector). ISPs / numerator dot-products are excluded from the
graph automatically.

## The NP-dbox catch (validated)

```
[1] NP-dbox  family.jl
   planar=True   cut_signature={'s':2,'t':3,'u':4}
   >> ISO to 'planar_smirnov_dbox' under [k2 -> -k2]
   !! LABEL/PROVENANCE MISMATCH: file is labeled nonplanar/crossed but the
      reconstructed graph is PLANAR in the canonical external-leg order.
      This is the NP-dbox-style mislabel.
```

`--self-test` validates: the NP-dbox **flags** (planar + iso to planar_smirnov
under `k2->-k2`), and the planar Smirnov dbox, the genuine crossed Tausk box,
the C3 rung-mass dbox, and the FF 3-loop banana **do not** false-flag. Its
`[compare]` leg runs the compare battery (`tests/census_battery.py`) on the
identity-census corpus: 30 published Feynman-integral families, each paired
with the graph the paper draws for it, several rows carrying more than one
family (a corrected family beside the one it replaced, or two variants). The
manifest `tests/fixtures/CENSUS_BATTERY.json` names, per family-level entry,
the sha-pinned fixture files (family, kinematics, drawn graph) and the
census's verdict string, and the compare mode must reproduce every string
exactly; the row-level rollups ("see families: X / Y") are joined from the
family verdicts; the dual cases (row 24's corrected family against two
readings of its drawing, row 24's double pentagon against the figure and
against the family's own edge list) carry their own expected strings; a
multi-family file is read by the family name the census judged (row 18); and
planted controls must flip or be refused by name (a mass moved, a leg class
moved, a seven-line top sector, undeclared leg classes, the cyclic-order
sub-finding, the other family of a multi-family file). The manifest also
carries the census's verdict counts and counting rule, which the battery
re-derives from the entries. The self-test also runs a **catalog integrity
sweep** (every catalog entry parses to its declared
propagator count, matches itself under the identity relabeling, reads its
declared planarity, and cross-matches no other entry — in particular the
planar/crossed dbox pair stays distinct) and a **signed-relabeling positive
control** (the one-loop box written with `k1 -> -k1` must be caught as ISO to
`massless_box_1l` under the sign flip, with no fixture file needed), and the
**loomcheck member's structural battery** as the `[loomcheck]` leg (below).

## Member: loomcheck (Yangian/loom/fishnet applicability screen)

The subfolder `loomcheck/` is a second graph screen shipped in this package: given a
position-space conformal integral as a dictionary of x_ij^2 powers (numerators as
negative powers), or an entry of the arXiv:2607.11645 ancillary basis file, it checks
the three published applicability conditions of the Yangian / loom / fishnet theorems
(internal-vertex conformal weight = D; planarity of the full graph with numerator
edges counted; the 2505.05550 face condition) and reports `yangian_class`.
SCREEN-ONLY: a FAIL proves exclusion from the cited theorems as published, a PASS
constructs nothing; the face check counts the outer face (conservative). From this
directory:

```bash
python3 -m loomcheck --selftest                       # 10 structural checks, rc 0/1
python3 -m loomcheck --basis 4Loop_int_basis.m [--targets 144,163,...] [--D 4] [--json OUT]
python3 loomcheck/loomcheck.py --selftest             # the script form, any cwd
```

`python3 topology_audit.py --self-test` runs the member's battery as its `[loomcheck]`
leg; `tests/test_loomcheck.py` is its pytest file. Details and limits: GUIDE.md,
"Member: loomcheck", and `loomcheck/README.md`.

## Scope and honest limits

- Designed for the **2-loop 4-point box family** (where the mislabel lives) and
  related 1-/2-loop boxes/triangles/ladders. The momentum-conservation vertex
  realizer is exact for these.
- **Needs >= 4 external legs on the boundary** for the planarity discriminator;
  with 2-3 legs a closure cycle is trivially planar, so the verdict is reported
  as `INCONCLUSIVE` rather than a (false) mismatch (e.g. a 3-point-routed family yaml).
- For large/many-loop families (e.g. a 3-loop banana with a single current) the realizer may not complete; it then returns `planar=None`
  ("realization-failed") and degrades gracefully — it never false-flags.
- The cut `t<->u` symmetry is **routing-dependent** and only an informational
  hint; **planarity is the rigorous discriminator**.

## Extending the catalog

Edit `CATALOG_SOURCES` in `topology_audit.py` (loops, exts, ext_subs,
propagators, note, and an optional `planar` key declaring the expected
leg-order planarity verdict — the self-test's catalog-integrity leg checks
every entry against it). Each new entry is automatically used in the
isomorphism search. Nine canonical representatives ship (the table below is generated from
`CATALOG_SOURCES`):

| entry | loops | note |
|---|---|---|
| `planar_smirnov_dbox` | 2 | massless planar double box (Smirnov hep-ph/9905323) |
| `crossed_tausk_dbox` | 2 | genuine nonplanar crossed double box, the (2,1,1) theta graph (Tausk hep-ph/9909506) |
| `massless_box_1l` | 1 | massless on-shell one-loop box |
| `massless_triangle_1l` | 1 | massless one-loop triangle (three-point) |
| `massless_bubble_1l` | 1 | massless one-loop bubble (two-point) |
| `sunrise_2l` | 2 | massless two-loop sunrise (two-point) |
| `ladder_vertex_2l` | 2 | massless planar two-loop ladder vertex (three-point) |
| `planar_triple_box_3l` | 3 | massless planar triple box (three-rung ladder) |
| `ud_ladder_3pt_offshell_3l` | 3 | massless three-loop three-point off-shell ladder, Usyukina-Davydychev Phi^(3) (legs p1 and p3 off shell, p2 on shell) |

Two-point entries use `exts: [p1, p2]` with `p2 = -p1` (the realizer needs
every external leg present so the vertex momentum sums close to zero).

Pure `sympy` + `networkx`. No amflow / Kira / heavy compute. networkx is a declared
dependency: without it the realizer / planarity / cut-signature / isomorphism legs of
`--self-test` and of the pytest battery SKIP by name (`networkx unavailable (declared
dependency; pip install networkx)`) and the battery exits 2 (a skip is an unknown, never
a pass); the tool itself still runs, printing `None` with that reason for every value
that needs the realized graph.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

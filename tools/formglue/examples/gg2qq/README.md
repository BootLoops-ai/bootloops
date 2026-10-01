# gg -> q qbar in FORM (formglue worked example)

The worked example of the FORM tool page (https://bootloops.ai/tools/form.html):
the tree-level squared amplitude for g g -> q qbar in massless QCD, computed
entirely in FORM and verified symbolically against the textbook result, plus
the one-loop diagram enumeration and numerator algebra that feed the
`formglue.form_route` momentum router. Everything here is plain FORM and
Python; FORM itself is an external engine (GPL-3.0, github.com/form-dev/form)
that you obtain upstream -- no FORM source or binary is included.

Files:

- `gg2qq.frm` -- tree-level |M|^2: three diagrams, d-dimensional Dirac
  traces, exact SU(N) color (inline Fierz reduction, T_R = 1/2), physical
  axial-gauge gluon polarization sums. Output: the spin/polarization/color
  summed |M|^2/g^4 as a Laurent polynomial in s, t, u with d and N symbolic,
  written to `gg2qq.out`. Run: `form gg2qq.frm`. Mutation knobs for the
  self-test: `form -D MUT=1 gg2qq.frm` (flipped three-gluon vertex sign),
  `form -D MUT=2 gg2qq.frm` (broken color Fierz identity).
- `check_gg2qq.py` -- runs FORM (bare `form` on PATH, or set `FORM_BIN`),
  parses the output, averages over initial states, sets d = 4, N = 3,
  u = -s-t, and checks exact symbolic equality with the Ellis-Stirling-Webber
  reference (t^2+u^2)/(6tu) - (3/8)(t^2+u^2)/s^2, then reruns both mutations
  and requires the comparison to fail. Needs python3 + sympy. Exit 0 only if
  all gates pass; exit 2 (SKIP) if FORM or sympy is absent.
- `enum_1loop.frm` -- one-loop diagram enumeration for the same process with
  FORM 5's integrated diagram generator (QCD model with ghosts, one flavor):
  3 tree diagrams (control), 109 connected one-loop diagrams, 30 after
  removing tadpoles, snails and external self-energies, 7 one-particle
  irreducible. Needs FORM >= 5 (`diagrams_`, `OnShell_`, `NoTadpole_`).
- `box_numerator.frm` -- full numerator algebra for one representative
  one-loop box interfered with the conjugate tree amplitude: d-dimensional
  trace, color reduction, projection onto the loop scalar products
  {k^2, k.p1, k.p2, k.p3}. The output numerator is the integrand ready for
  tensor/IBP reduction; the loop integration itself is not performed here.
- `gg2qq_topologies.json` -- the parsed topology data of the enumeration:
  3 tree + 30 one-loop entries extracted from the generator's `node_`/`edge_`
  output by `form_route` (per-edge endpoints, field, momentum, arrow
  direction). It ships alongside so the renderer runs without FORM installed.
- `render_gg2qq_diagrams.py` -- renders the two diagram figures
  (`form-gg2qq-tree.png`, `form-gg2qq-1loop.png`) from
  `gg2qq_topologies.json` with deterministic layout rules; line styles are
  keyed off the field labels (curly gluons, solid arrowed quarks, dashed
  arrowed ghosts) and nothing is drawn by hand. Needs python3 + matplotlib
  only (no FORM); output is byte-reproducible for a fixed matplotlib version.

How to run (from this directory, or from a scratch copy of it):

```sh
form gg2qq.frm                     # writes gg2qq.out
python3 check_gg2qq.py             # PASS x3, "All gates passed." (rc 0)
form box_numerator.frm             # TERMS after trace 4427 -> integrand-ready 1225
form enum_1loop.frm                # FORM 5 only: COUNT lines 3 / 109 / 30 / 7
python3 render_gg2qq_diagrams.py   # the two PNG figures
```

To route momenta on a fresh enumeration, feed the `diagrams_` output through
`formglue.form_route` (see the package GUIDE.md, INVOKE); the router's own
tests use fixtures generated from the same reference-manual QCD model.

`gg2qq.frm` and `box_numerator.frm` run on FORM 4.x as well as FORM 5
(`gg2qq.frm` output verified byte-identical between 4.3 and 5.0.1);
`enum_1loop.frm` needs the FORM 5 diagram generator. These files are not
collected by the package battery (`pytest tools/formglue`), which stays
engine-free by default.

Credit: FORM is by Jos Vermaseren and the form-dev community (see the package
GUIDE.md CREDIT paragraph and REFERENCES.md at the repository root); the
reference cross section is the standard result as given in Ellis, Stirling &
Webber, "QCD and Collider Physics".

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

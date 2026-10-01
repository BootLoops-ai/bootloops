# ABACUS — certified counting on abelian fourfolds

Tool page: https://bootloops.ai/tools/abacus.html

KIND: package (python + Arb/acb_theta via python-flint + a Julia Eichler bridge; PARI
only on planted-truth cross-paths).

PURPOSE: certified point counts / Frobenius charpolys / L-factors for abelian
FOURFOLDS presented as a period lattice with polarization and (optionally)
Z[sqrt(-5)]-action: hermitian/period front door -> tau in H_4 -> certified theta
receipts + Weil-box integer isolation + exact rho/Rosati receipts.

USE-WHEN:

- Counting an abelian fourfold (KS route) — nothing else in the kit counts abelian
  fourfolds.
- Need certified genus-4 (or genus-2 block) theta constants with a PROVED desk tail
  bound dominating acb_theta's claimed radius.
- Need exact O_K = Z[sqrt(-5)] action receipts (rho^2=-5, Rosati, sign-pin s) on a
  period lattice.

NOT-FOR: outputs are EVIDENCE only behind the shipped battery (B1-B4 + B4s + B5 +
N1-N5) — any tool change re-runs the FULL battery first; not a generic curve counter
(use PARI ellap/hyperellcharpoly directly for curves); sign conversions happen ONLY at
abcount_s0.convert_sign_convention (never add a second site).

## Layout and invocation

The canonical, sha-pinned body lives at `build/abcount/` inside this directory; the
`abcount_s0..s3.py`, `hecke_desk.py`, and `run_s0..s3.py` files at the package root
are header-marked compat shims that forward to it. The harnesses resolve their
sha-pins relative to their own file location, so this tree runs unmodified.

- CLI: `python3 tools/abacus/run_s0.py` and `run_s3.py` (`run_s1.py`/`run_s2.py`
  refuse against the shipped SEEDS.json — see GATES) — routes T (theta receipts),
  D (integer assembly), C (counting truth side) with declared preconditions per route.
- Import form: `import abcount_s0` etc. with `tools/abacus/` on `sys.path`.

INPUTS: period matrix Pi (balls), integer polarization E, optional integer 8x8 rho;
`build/abcount/battery/SEEDS.json` pins (dyadic pilot tau, primes, walls);
sign_convention mandatory on every artifact.

OUTPUTS: per-stage RUN_OUTPUT.json receipts; battery BATTERY_RESULT.json; L-factors
as exact integer polys; named-predicate refusals (FAIL-WEIL-BOX-EMPTY,
UNDECIDED-PRECISION+needed-dps, FAIL-SIGN-PIN, ...).

Dependencies: `python-flint` (Arb/acb_theta; MIT over LGPL FLINT) and `cypari2`
(GPL-2.0-or-later; PARI ellap / hyperellcharpoly on the planted-truth cross-paths;
`run_s0.py` and the battery import it at run time, it is not bundled), both
pip-installable; the PARI/GP binary `gp` on PATH (the S3
battery's B4/B4s/B5 truth-side legs drive it directly and refuse by name when
it is missing). Two harness paths reference software
outside this package: `abcount_s1.py` reads the Eichler.jl project root from
`ABACUS_EICHLER_PROJECT` for the theta/period bridges (the same package ships in
this repository under `upgrades/Eichler.jl/` — point the variable at your
install; the bridges refuse loudly when it is unset), and `run_s1/run_s3` add
the in-repo Lockpick (PSLQ/recognition) toolkit path for the j-recognition
control. The bridges run whatever `julia` is on PATH; `ABACUS_JULIA`
overrides the executable (an absolute path or a name on PATH), and the
julia-call sites refuse by name when neither resolves. The julia bridges
need the Eichler project's pinned packages
installed once (`julia --project=<project> -e 'import Pkg; Pkg.instantiate()'`);
`selftest.py` performs that step itself before entering the battery — the
first run may download and precompile for a few minutes — and refuses with a
named error when it cannot.

GATES: C7 dominance gate (desk tail bound must dominate acb_theta radius at matched
truncation; ANY violation = STOP); Weil-box FULL ENUMERATION (never per-coefficient
rounding); negative controls N1-N5 must FAIL by name; battery = acceptance authority.
Battery: `run_s0.py` → S0-PASS 7/7; `run_s3.py` → BATTERY-PASS(11/11, B5
RUNS+PASS), 22/22 checks, ~47 s / 0.11 GB from a scratch cwd. `run_s3.py` is
the acceptance harness; `run_s1.py`/`run_s2.py` pin other SEEDS revisions
and refuse against the shipped SEEDS.json. One command:
`python3 tools/abacus/selftest.py` runs both stages from a scratch working
directory (and defaults ABACUS_EICHLER_PROJECT to the in-repo
`upgrades/Eichler.jl` when unset).

FOOTGUNS:

- Negating rho flips the sign pin s and every emitted b_P (a_P invariant) — the B4s
  detector exists precisely for this; check s before consuming b_P.
- Odd-characteristic thetas at z=0 are exactly 0; acb_theta's tiny radius there is a
  rounding artifact, not a truncation tail — use the exact-zero containment gate, not
  dominance.
- The module/harness filenames keep the abcount_* stems (file and stem are what
  the verification harness sha-pins); ABACUS is the package name — do not
  "clean up" the stems: renaming a pinned file makes verification fail.

Math/proof authority ships in-tree:
`build/abcount/manuals/abcount.md` (C7 theta-tail + C3/C2b Weil-box theorems, proved).

CREDIT: the certified theta evaluations are FLINT's acb_theta module by Jean Kieffer
(Kieffer, arXiv:2203.02000; Elkies & Kieffer, arXiv:2505.22382) on top of Arb (F.
Johansson 2017) via python-flint; planted-truth cross-paths use PARI/GP (The PARI Group,
Univ. Bordeaux) through cypari2, with curve labels from Cremona's tables / the LMFDB;
the Weil-box isolation rests on the classical Weil bounds. We thank J. Kieffer and the
FLINT developers — genus-4 certified theta constants exist as a library call because of
them.
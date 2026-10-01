"""GALOIS — motivic Galois coaction engine + exact MZV relation ring.

The coaction-cut engine, with the exact MZV relation ring as its
relation-ring MEMBER (subpackage mzvring) — sector21's reduction oracle.

DISAMBIGUATION: this is the MOTIVIC Galois coaction (cosmic Galois,
F. Brown) — not differential-Galois-group computations.

Members:
  core     — Goncharov-coproduct Galois derivations D3/D5/D7 on MZVs (zpoly,
             t0deep IKZ tables), on H-words, and on Brown svHPLs (calL) via
             the formal sv T-tables; sv-space decomposition; conjugation;
             ansatz-element conversion (the I[z,w,0]/f-token conventions).
  vspace   — allowed-space machinery: parses in-tree SOLVED conformal-
             integral results into pure sv block vectors, saturates under
             D3/D5/D7, exposes weight+parity fibers with exact annihilator
             functionals (the coaction-stability rows).
  sector21 — D_m (m odd 3..19) + sv projection on depth<=3 odd-index MZVs at
             odd weights 11..21, reduction oracle = the mzvring member's
             reference ring pkls (in place of the w<=8 t0deep tables).
  svmap21  — ALGEBRAIC single-valued map zeta_sv[this package's] on depth<=3
             odd-index MZVs (w<=21): Brown's sv projection on the f-alphabet
             (arXiv:1309.5309 eq (7.3)), this package's orient-B transport
             (61/61-pinned), exact phi^-1 solve over the ring basis;
             `selftest --receipt PATH` writes a validation receipt on PASS.
  mzvring  — the relation-ring member (subpackage): exact double-shuffle
             canonicalization of {zeta(n), mzv(a,b), mzv(a,b,c)} monomials
             over Q, depth<=3, NO PSLQ; survivors reproduce the BK depth<=3
             dims. `from galois import mzvring; R = mzvring.load_ring()`.

NAME-COLLISION FENCE: PyPI ships an unrelated `galois` package (GF/finite
fields). Anything with
this repo's tools/ dir ahead of site-packages on sys.path gets
THIS package under `import galois`. If you need both, import the
PyPI one first or by explicit path.

DATA PLACEMENT (tools link, never copy): the reference ring
pickles live at mzvring.RING_BANK; the reference sector
tables at sector21.TABLE_BANK; the t0deep w<=8 tables at
core.T0DEEP — all pointed at via env vars or explicit paths. Nothing
multi-MB lives in tools/.

Gate receipts are not shipped; the selftests regenerate them against
your own banks.
"""

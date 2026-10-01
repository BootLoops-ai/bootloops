#!/usr/bin/env python3
"""
gvderive.py (FRAME/GV-DERIVATION, periods wing) — derive + CERTIFY a
GV working set from the card toric layer instead of trusting a repo table.
GENERALIZED from one-card derivation scripts (not included in this release);
the extraction engine is DELEGATED to pipeline/geff_series.gv_extract — no
extraction logic duplicated here.

Capabilities:
  crosswalk_repo   exact repo->card-frame charge map E^T q_card = B q_repo
                   (a supplied map is verified, held-out)
  derive_gv        box extraction from toric data ONLY (sigma/eps pinned)
  gate_working_set h1 held-out exactness, h2 integrality, h3 level-0
                   conifold-ray law, h4 margin shells -> margin-CERTIFIED
                   complete level, h5 new-class census
  ladder_from_gv   exact racetrack ladders S_L, sum n (q.M)^2 per level
  completeness_band  the box-to-box shift band, ALWAYS labeled

GV-BAND LABELING LAW (honesty line): completeness is certified
only THROUGH the margin level (min shell level - 1); beyond it the W0 error
is a LABELED band = box-to-box shift x5 — a measured law, NOT a theorem.
Any ball quoting this module MUST carry the band as a separate labeled term
(never folded silently into a certified radius).
FLOAT PIN SCOPING (honesty line): vacuum values taken from a float64 source
are point pins promoted to rationals; results are CONDITIONAL on them with
a ~1e-13 rel sensitivity floor.  M2 AMENDMENT: uniqueness backing this
dictionary is zero-kernel exact identities (see conifold.py header), which
is weaker than ads-5-81-class rank saturation — state it wherever quoted.
Battery: none registered in this release (the worked example this module
was validated on is not included in the package).
"""
import os, sys
from fractions import Fraction as Fr

_PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
if _PIPE not in sys.path:
    sys.path.insert(0, _PIPE)
import geff_series as GS                                    # noqa: E402


def _solve_int(A, b, h):
    """exact solve A x = b (h x h), assert integer solution."""
    M = [[Fr(A[i][j]) for j in range(h)] + [Fr(b[i])] for i in range(h)]
    for c in range(h):
        pr = next(r for r in range(c, h) if M[r][c] != 0)
        M[c], M[pr] = M[pr], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(h):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][t] - f * M[c][t] for t in range(h + 1)]
    x = [M[i][h] for i in range(h)]
    assert all(v.denominator == 1 for v in x), f"non-integer solve {x}"
    return [int(v) for v in x]


def crosswalk_repo(card, spec, nu):
    """Repo mirror-GV table -> card Mori frame: q_card solves
    E^T q_card = B q_repo (E = card frames.E_rows_coni, B = spec basis
    transformation).  Gates: every image nonnegative; card pins are a
    value-consistent subset.  Returns [(q_card, n, level)]."""
    h = card.h
    E = card.raw["frames"]["E_rows_coni"]
    B = spec["data"]["basis transformation"]
    ET = [[E[j][i] for j in range(h)] for i in range(h)]
    out = []
    for row in spec["data"]["mirror GVs"]:
        qc, val = row[:h], row[h]
        qb = [sum(B[i][j] * qc[j] for j in range(h)) for i in range(h)]
        qcard = _solve_int(ET, qb, h)
        lv = sum(qcard[a] * nu[a] for a in range(h))
        out.append((tuple(qcard), int(val), lv))
    assert all(all(x >= 0 for x in q) for q, _v, _l in out), \
        "crosswalk: negative charge in card frame"
    for q, v in card.gv_pinned.items():
        match = [vv for qq, vv, _l in out if qq == tuple(q)]
        assert match in ([int(v)], []), f"crosswalk pin mismatch {q}"
    return out


def derive_gv(card, box):
    """Toric-data-only GV extraction on the box (h2 integrality asserted).
    NO repo curve is an input.  Returns {q: n} (nonzero classes)."""
    nGV, _G = GS.gv_extract(card, tuple(box), tuple(card.sigma), card.eps)
    nz = {q: v for q, v in nGV.items() if v != 0}
    assert all(Fr(v).denominator == 1 for v in nz.values()), "h2 FAIL"
    return {q: int(v) for q, v in nz.items()}


def gate_working_set(card, nz, box, nu, q_cf, held_out=(), extras=()):
    """h1/h3/h4/h5 gates on a derived set; AssertionError = gate FAIL.
    held_out: [(q, n), ...] classes the extraction must REPRODUCE (they are
    never inputs); extras likewise.  Returns receipts incl. the margin-
    certified complete level (h4: min shell level - 1)."""
    h = card.h
    lv = lambda q: sum(q[a] * nu[a] for a in range(h))
    # h1 held-out exactness (inside the box only — scope reported)
    n_in, bad = 0, []
    for q, v in list(held_out) + list(extras):
        q = tuple(q)
        if all(q[a] <= box[a] for a in range(h)):
            n_in += 1
            if nz.get(q, 0) != int(v):
                bad.append((q, int(v), nz.get(q, 0)))
    assert not bad, f"h1 FAIL held-out: {bad[:5]}"
    for q, v in card.gv_pinned.items():
        if all(q[a] <= box[a] for a in range(h)):
            assert nz.get(tuple(q), 0) == int(v), f"h1 FAIL pin {q}"
    # h3 level-0 = conifold ray only (k=1), no negative levels
    lev0 = sorted((q, v) for q, v in nz.items() if lv(q) == 0)
    assert [q for q, _v in lev0] == [tuple(q_cf)], f"h3 FAIL level-0: {lev0}"
    assert not [q for q in nz if lv(q) < 0], "h3 FAIL negative level"
    # h4 margin shells -> completeness horizon
    shell = sorted((q, v) for q, v in nz.items()
                   if any(q[a] >= box[a] for a in range(h)))
    complete_level = (min(lv(q) for q, _v in shell) - 1 if shell
                      else max(map(lv, nz)))
    # h5 census: new classes per level beyond the held-out+extras set
    known = {tuple(q) for q, _v in list(held_out) + list(extras)}
    new_by_level = {}
    for q in nz:
        if q not in known:
            new_by_level.setdefault(lv(q), []).append(q)
    return dict(n_classes=len(nz), held_out_in_box=n_in,
                complete_level=complete_level,
                shell_classes=[[list(q), v] for q, v in shell],
                new_by_level={L: sorted(v) for L, v in new_by_level.items()},
                label="completeness CERTIFIED through complete_level only "
                      "(margin-shell law); beyond it use completeness_band")


def ladder_from_gv(card, nz, nu):
    """Exact ladders: S[L] = sum n (q.M); S2[L] = sum n (q.M)^2 (ints)."""
    h = card.h
    S, S2 = {}, {}
    for q, v in nz.items():
        L = sum(q[a] * nu[a] for a in range(h))
        qM = sum(int(q[a]) * int(card.M[a]) for a in range(h))
        S[L] = S.get(L, 0) + v * qM
        S2[L] = S2.get(L, 0) + v * qM * qM
    return S, S2


def completeness_band(shift_rel, margin=5):
    """GV-BAND LABELING LAW: band = |box-to-box W0 shift| x margin.
    Returns a dict that MUST travel with any ball quoting the table."""
    band = abs(float(shift_rel)) * margin
    return dict(band_rel=band, shift_rel=float(shift_rel), margin=margin,
                label="LABELED GV-completeness band (box-shift x%d law; "
                      "measured, NOT a certified bound — never fold into "
                      "a certified radius)" % margin)

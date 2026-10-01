#!/usr/bin/env python3
"""mplll_graded.py — weight-graded block reduction (--method graded).

Bases are transcendentality-graded; the true relation respects the grading, so
the (K+1)-lattice is block-triangular in weight.  Per-block QR-BKZ at dim
~2K_w << 2K makes β=K_w exact-SVP cheap AND localises span-deficiency to a
single weight instead of poisoning the whole fit.

v2: per-weight lattices at BOTH precisions are built up-front and reduced in a
single reduce_rows_batch call (one julia/fplll dispatch instead of 2W).
Sequential residual peel r ← r − Σ c_i f_i is preserved: the top-weight block
sees I directly; lower blocks whose pre-built (target=I) lattice fails to lock
against the peeled residual get a single-shot rebuild (rare fallback)."""
import sys, os, math, re
import mpmath as mp
from . import mplll                                       # g3 package-relative
from .mplll_lattice import reduce_rows_batch, reduce_rows


def infer_weight(name):
    """Heuristic transcendental weight from a basis-function name."""
    n = name.strip()
    if n in ('1',): return 0
    w = 0
    for m in re.finditer(r'Li(\d)', n): w += int(m.group(1))
    for m in re.finditer(r'L[stu]\d?0?\^?(\d)?', n): w += int(m.group(1) or 1)
    for m in re.finditer(r'Lu(?:4|16)?\^?(\d)?', n): w += int(m.group(1) or 1)
    if 'pi2' in n or 'pi^2' in n: w += 2
    if 'zeta' in n or 'z3' in n: w += 3
    for m in re.finditer(r'G(\d+)', n): w += len(m.group(1))
    if n.startswith(('m1', 'dm1', 'd2m1', 'I0_', 'I4_', 'I16_', 'I2_')): w = 5
    if 'Eb' in n or 'Cl2' in n: w = 3
    return max(w, 0)


def _lock(rel, rel2, r, Fw, fit, d):
    if rel is None or rel2 is None: return False, None
    from .pslq_gate import canonicalize
    if canonicalize(rel) != canonicalize(rel2): return False, None
    isd, _ = mplll.heldout_digits(rel, r, Fw, fit, d)
    return isd > max(5, d * 0.4), isd


def graded_fit(I, F, names, weights, fit, ho, d, hmax=10 ** 30, backend=None):
    """Weight-graded QR-BKZ.  weights: list[int] len K, or None → infer_weight+WARN."""
    K, N = len(F), len(I)
    if weights is None:
        weights = [infer_weight(n) for n in names]
        sys.stderr.write("[mplll_graded] WARN: inferring weights heuristically; "
                         "pass --weights for production\n")
    ws = sorted(set(weights), reverse=True)
    blocks = {w: [i for i, wi in enumerate(weights) if wi == w] for w in ws}
    d2 = max(10, d - 20)
    # Pre-build ALL per-weight lattices at (d, d2) with target=I; batch-reduce once.
    latt, tags = [], []
    for w in ws:
        Fw = [F[i] for i in blocks[w]]
        for dd in (d, d2):
            rows, _ = mplll._qrbkz_lattice(I, Fw, fit, dd)
            latt.append(rows); tags.append((w, dd))
    Rs, bknd, wall = reduce_rows_batch(
        latt, backend, bkz_beta=min(max(len(b) for b in blocks.values()) + 1, 40))
    Rmap = {tags[i]: Rs[i] for i in range(len(tags))}
    # Sequential peel with fallback rebuild if pre-built lattice doesn't lock.
    with mp.workdps(d + 40):
        r = [mplll._mpc(I[j]) for j in range(N)]
    per_w, coeffs = {}, [mp.mpf(0)] * K
    for w in ws:
        idx = blocks[w]; Fw = [F[i] for i in idx]; Kw = len(idx)
        rel, _ = mplll._extract(Rmap[(w, d)], Kw, hmax)
        rel2, _ = mplll._extract(Rmap[(w, d2)], Kw, hmax)
        locked, isd = _lock(rel, rel2, r, Fw, fit, d)
        with mp.workdps(d + 40):    # never compare residuals at ambient dps
            differs = any(abs(r[j] - mplll._mpc(I[j])) > mp.mpf(10) ** (-d)
                          for j in fit)
        if not locked and differs:
            # residual differs from I → rebuild against peeled r (fallback)
            mats = [mplll._qrbkz_lattice(r, Fw, fit, dd)[0] for dd in (d, d2)]
            R2s, _, w2 = reduce_rows_batch(mats, backend, bkz_beta=min(Kw + 1, 40))
            wall += w2
            rel, _ = mplll._extract(R2s[0], Kw, hmax)
            rel2, _ = mplll._extract(R2s[1], Kw, hmax)
            locked, isd = _lock(rel, rel2, r, Fw, fit, d)
        rec = {"K_w": Kw, "locked": locked, "insample_d": isd,
               "relation": rel, "floor_d": None if locked else isd}
        if locked:
            with mp.workdps(d + 40):
                cw0 = mp.mpf(rel[0])
                for k, i in enumerate(idx): coeffs[i] += mp.mpf(rel[k + 1]) / cw0
                for j in range(N):
                    # heldout_digits semantics: block contributes -Σ rel[k+1]·f_k/c0,
                    # so the peel ADDS Σ rel[k+1]·f_k/c0 (a '-=' here would
                    # double the residual and falsely fail lower weights)
                    r[j] += mp.fsum(mp.mpf(rel[k + 1]) * mplll._mpc(Fw[k][j])
                                    for k in range(Kw)) / cw0
            rec["heldout_d"], _ = mplll.heldout_digits(rel, I, [F[i] for i in idx], ho, d)
        per_w[w] = rec
    from fractions import Fraction
    fr = [Fraction(mp.nstr(c, 25)).limit_denominator(10 ** 15) for c in coeffs]
    lcm = 1
    for f in fr: lcm = lcm * f.denominator // math.gcd(lcm, f.denominator)
    rel = tuple([lcm] + [int(f * lcm) for f in fr])
    hd, _ = mplll.heldout_digits(rel, I, F, ho, d)
    isd, _ = mplll.heldout_digits(rel, I, F, fit, d)
    all_locked = all(v["locked"] for v in per_w.values())
    return {"method": "graded", "K": K, "dps": d, "per_weight": per_w,
            "relation": rel, "heldout_min_d": hd, "insample_min_d": isd,
            "all_locked": all_locked, "lll_wall_s": round(wall, 2), "backend": bknd,
            "capacity": mplll.capacity(rel, K, len(fit), d),
            "status": "HIT" if hd >= mplll.GATE_D and all_locked else
                      ("PARTIAL" if any(v["locked"] for v in per_w.values()) else "NULL")}

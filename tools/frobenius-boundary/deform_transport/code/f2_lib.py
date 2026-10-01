#!/usr/bin/env python3
"""f2_lib.py — exact-arithmetic Legendre-curve helpers (vendored engine for
the deform_transport member).

Everything here is exact integer / F_p arithmetic; no floats in any
verdict-bearing quantity.  Three helpers:

  1. chi-table Legendre a_p engine — a_p = -sum_u chi(u(u-1)(u-lam)), with
     the quadratic character precomputed as a table per prime (pure caching).
  2. F_p^2 char-sum a_{p^2} engine: for the Legendre curve with lam in F_p
     read over F_p^2 = F_p[s]/(s^2 - d), d a non-residue,
     a_{p^2} = -sum_{u in F_p^2} chi_2(f(u)),  chi_2(t) = chi_p(Norm(t)),
     Norm(a+bs) = a^2 - d b^2.  At fiber-inert primes the fiber has no
     F_p-point but the curve (lam = x0^2 in F_p) does have an F_p^2 story,
     and supersingularity (p | a_{p^2}) is readable there.  Gate per
     (fiber, prime): the curve is defined over F_p, so the exact identity
     a_{p^2} = a_p^2 - 2p must hold; the char-sum route is independent of
     the a_p route, so the identity is a genuine cross-check.
  3. Deviant-residue lookup — per-prime deviant-residue sets loaded from
     exhaustive-scan receipt files.  The receipts are NOT distributed with
     the repo: point DEFORM_F1_RECEIPTS at a directory holding
     FANF1_RECEIPT.json / FANF1_EXT_RECEIPT.json to use load_f1_deviants
     (it refuses loudly otherwise).
"""
import json, os


def _f1_dir():
    d = os.environ.get("DEFORM_F1_RECEIPTS")
    if not d or not os.path.isdir(d):
        raise RuntimeError(
            "f2_lib.load_f1_deviants: exhaustive-scan receipt files are not "
            "distributed with the repo; set DEFORM_F1_RECEIPTS to a directory "
            "holding FANF1_RECEIPT.json (+ optionally FANF1_EXT_RECEIPT.json)")
    return d


def chi_table(p):
    """chi[t] = quadratic character of t mod p (0 for t=0)."""
    tab = [0] * p
    # squares via direct enumeration (no pow calls): mark QRs
    for t in range(1, p):
        tab[t] = -1
    for t in range(1, (p + 1) // 2):
        tab[t * t % p] = 1
    return tab


def legendre_ap_fast(lam, p, chi):
    """a_p of y^2 = u(u-1)(u-lam) over F_p; chi = chi_table(p)."""
    s = 0
    for u in range(p):
        s += chi[u * (u - 1) % p * ((u - lam) % p) % p]
    return -s


def nonresidue(p, chi):
    for d in range(2, p):
        if chi[d % p] == -1:
            return d
    raise RuntimeError("no non-residue found (p=2?)")


def ap2_charsum(lam, p, chi, d=None):
    """a_{p^2} of y^2 = u(u-1)(u-lam), lam in F_p, read over F_p^2 = F_p[s]/(s^2-d)
    by the norm-character sum.  O(p^2)."""
    if d is None:
        d = nonresidue(p, chi)
    lam %= p
    s = 0
    for a in range(p):
        aa1 = a * (a - 1) % p
        al = (a - lam) % p
        for b in range(p):
            db2 = d * b * b % p
            re1 = (aa1 + db2) % p          # (a+bs)(a-1+bs) real
            im1 = b * (2 * a - 1) % p      # imag
            re2 = (re1 * al + d * im1 * b) % p
            im2 = (re1 * b + im1 * al) % p
            s += chi[(re2 * re2 - d * im2 * im2) % p]
    return -s


def load_f1_deviants():
    """{prime:int -> set of deviant residues} + per-prime pole counts, from
    exhaustive-scan receipt files under DEFORM_F1_RECEIPTS (refuses if unset)."""
    dev = {}
    meta = {}
    d = _f1_dir()
    main = json.load(open(os.path.join(d, "FANF1_RECEIPT.json")))
    ext_path = os.path.join(d, "FANF1_EXT_RECEIPT.json")
    ext = json.load(open(ext_path)) if os.path.exists(ext_path) else {"scan_ext": {}}
    for src, key in ((main, "scan"), (ext, "scan_ext")):
        for pstr, row in src[key].items():
            p = int(pstr)
            dev[p] = {int(dd["x0"]) for dd in row.get("deviants", [])}
            meta[p] = {"n_scanned": row.get("n_scanned"),
                       "majority_profile": row.get("majority_profile")}
    return dev, meta


def rational_heights_menu(H):
    """Reduced a/b with 1 <= a,b <= H, gcd(a,b)=1, (a,b) != (1,1), both signs.
    Deterministic order: height max(a,b) asc, then b, then a, then sign (+ first).
    Returns list of (num, den) with num carrying the sign."""
    from math import gcd
    pts = []
    for h in range(1, H + 1):
        row = []
        for b in range(1, h + 1):
            for a in range(1, h + 1):
                if max(a, b) != h or gcd(a, b) != 1 or (a, b) == (1, 1):
                    continue
                row.append((a, b))
        for a, b in sorted(row, key=lambda t: (t[1], t[0])):
            pts.append((a, b))
            pts.append((-a, b))
    return pts

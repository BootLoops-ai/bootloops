#!/usr/bin/env python3
"""gen_refined.py — REFINED S-unit law:
the fibre at z has small conductor <=> z is S-integral on P1 minus the
singular divisor: for ONE small prime set S, ALL of
   lead(f_z), const(f_z)  (0/infinity punctures: z and 1/z S-supported), and
   N(m_i(z)) for EVERY singular factor m_i of the card
have support inside S.  m_i = primitive ascending integer minpoly of each
finite nonzero singular value (union: carded singular_locus + leadpoly
factors; origin factor dropped).  For rational z=p/q, N(m(z)) = m(z); for
quadratic z, N = field norm of m(z) (equivalently Res(m,f_z)/lc(f_z)^deg m).
Verified by hand on all six knowns before coding: T1 {2,7}, T2 {2},
T3 {3,11}, T4 {2,7}, X118 {2,3,11}, X290 {2,3,7}.

S families (generation nets): S subset of {2,3,5,7,11,13} with |S| <= 3,
plus S = {2,p} for prime p <= 29.  ACCEPTANCE is family-free: compute
S_min(z) = union of the supports above and accept iff S_min fits some
allowed family.  S_min is recorded per candidate (and per gate point).

NET TRUNCATIONS (stated with every run; acceptance itself is exact):
  rational net: exponent |a_p| <= 12 per prime, Weil height <= 40
    (transport-price law); no de Weger reduction is applied.
  quadratic net: per squarefree D <= 100, pool = group elements gamma
    (fund. unit u_D, prime-ideal generators pi_p over S u ram(D), rational
    p in S; gen_sunit machinery) with h(gamma) <= 20, <= 600/cell,
    best-first — PLUS singular-anchored shifts z = s_i + gamma for every
    rational s_i (and quadratic s_i whose field is Q(sqrt D)); then the
    full family-free acceptance check; h(z) <= 40.
  Completeness caveats are recorded in the SIZES2 receipt meta.

Paths: output dir env-overridable (TERRIER_SUNIT_OUT, default
HERE/candidates2) via _outdir(); atlas/known-orbit paths as in gen_sunit
(TERRIER_SUNIT_ATLAS, TERRIER_KNOWN_ORBIT).
"""
import json, math, os, sys, time
from fractions import Fraction as Fr
from math import gcd, log

import gen_sunit as G

HERE = os.path.dirname(os.path.abspath(__file__))


def _outdir():
    return os.environ.get("TERRIER_SUNIT_OUT",
                          os.path.join(HERE, "candidates2"))


BIG6 = (2, 3, 5, 7, 11, 13)
EXTRA = (17, 19, 23, 29)
ALLOWED_P = BIG6 + EXTRA
A_MAX = 12          # rational exponent bound (net cap)
H_RAT = 40.0        # rational height cap
H_POOL = 20.0       # quadratic gamma-pool height cap (net cap)
POOL_CAP = 600      # rows per (Sgen, D) cell (net cap)
H_QUAD = 40.0       # final quadratic candidate height cap
BLOWUP = 10 ** 5    # per-op list alarm threshold (design default)


BOX = frozenset(BIG6)   # FROZEN S-box {2,3,5,7,11,13}
SCOPE21 = ["AESZ34", "AESZ22", "AESZ118", "AESZ17", "AESZ290", "AESZ100",
           "AESZt101", "h5", "h10", "h2222", "h33", "h223", "h24", "h8",
           "h6", "h212", "h44", "h34", "h46", "h66", "h26"]


def allowed_smin(sm):
    """FROZEN box law: S_min within {2,3,5,7,11,13}.
    (Supersedes the earlier family law |S|<=3 or {2,p<=29}: the EXTRA
    {2,p>13} families are OUT OF BOX and not accepted.)"""
    return sm <= BOX


def strip_p(n):
    """divide out ALLOWED_P primes -> (leftover, support set)."""
    n, sup = abs(n), set()
    for p in ALLOWED_P:
        if n % p == 0:
            sup.add(p)
            while n % p == 0:
                n //= p
    return n, sup


def sing_factors(card):
    """primitive ascending integer factor tuples of the singular divisor
    (carded singular_locus minpolys UNION leadpoly irreducible factors,
    origin dropped, content stripped, dedup up to sign)."""
    import sympy as sp
    z = sp.symbols("z")
    out = set()
    P = sum(a * z ** i for i, a in enumerate(card["coeffs"][-1]))
    for f, _m in sp.factor_list(P)[1]:
        cs = [int(c) for c in reversed(sp.Poly(f, z).all_coeffs())]
        if cs[0] == 0:
            continue
        g = 0
        for c in cs:
            g = gcd(g, abs(c))
        t = tuple(c // g for c in cs)
        out.add(t if t[-1] > 0 else tuple(-c for c in t))
    for e in card.get("singular_locus", []):
        cs = [int(c) for c in e["minpoly"]]
        if cs[0] == 0:
            continue
        g = 0
        for c in cs:
            g = gcd(g, abs(c))
        t = tuple(c // g for c in cs)
        out.add(t if t[-1] > 0 else tuple(-c for c in t))
    return sorted(out)

# --------------------------------------------------------- rational tier
RAT_SUPP_MAX = 4    # rational net: max #distinct primes in supp(z) (NET
#   TRUNCATION — the frozen box allows supp up to all 6; 5- and 6-prime
#   support monomials are NOT enumerated (24^5*6 + 24^6 ~ 2.4e8, infeasible).
#   PRECISE completeness statement: rational tier complete for
#   z = +-prod p^a_p, supp(z) <= 4 primes of the box, 1 <= |a_p| <= 12,
#   Weil height <= 40; acceptance (all singular distances in-box) is exact.


def rational_net():
    """Exact-support monomial stream over the FROZEN box: yields
    (num, den, h), one per monomial (no overlap: all exponents nonzero
    per chosen prime subset), supp <= RAT_SUPP_MAX, h <= H_RAT."""
    from itertools import combinations, product
    lim = math.exp(H_RAT)
    rng = [a for a in range(-A_MAX, A_MAX + 1) if a]
    yield (1, 1, 0.0)
    for k in range(1, RAT_SUPP_MAX + 1):
        for ps in combinations(BIG6, k):
            for exps in product(rng, repeat=k):
                num = den = 1
                for p, a in zip(ps, exps):
                    if a > 0:
                        num *= p ** a
                    else:
                        den *= p ** -a
                if max(num, den) <= lim:
                    yield (num, den, log(max(num, den)))


def check_rational(num, den, sgn, factors):
    """family-free acceptance for z = sgn*num/den.  -> S_min or None."""
    p, q = sgn * num, den
    lo, sup = strip_p(num)
    if lo != 1:
        return None
    lo, s2 = strip_p(den)
    if lo != 1:
        return None
    sup |= s2
    for m in factors:
        d = len(m) - 1
        r, qq = 0, 1
        for i in range(d, -1, -1):        # sum c_i p^i q^(d-i)
            r += m[i] * p ** i * qq
            qq *= q
        if r == 0:
            return None                   # z IS a singular point
        lo, s3 = strip_p(r)
        if lo != 1:
            return None
        sup |= s3
    return sup if allowed_smin(sup) else None

# -------------------------------------------------------- quadratic tier
def q_add(x, y, D):
    a1, b1, c1 = x
    a2, b2, c2 = y
    return G.q_reduce(a1 * c2 + a2 * c1, b1 * c2 + b2 * c1, c1 * c2)


def key_to_triple(key, D):
    """minpoly key 'A,B,C' -> one root triple in Q(sqrt D), or None."""
    A, B, C = (int(x) for x in key.split(","))
    disc = B * B - 4 * A * C
    if disc <= 0:
        return None
    e2, e = divmod_sq(disc, D)
    if e is None:
        return None
    return G.q_reduce(-B, e, 2 * A)


def divmod_sq(disc, D):
    """disc = e^2 * D exactly? -> (D, e) else (None, None)."""
    if disc % D:
        return None, None
    e = G.is_sq(disc // D)
    return (D, e) if e >= 0 else (None, None)


def families():
    """quadratic GENERATION nets: S subsets of the frozen box, |S| <= 3
    (truncation of the generation net only; acceptance is box-exact).
    EXTRA {2,p>13} families dropped at freeze (out of box)."""
    from itertools import combinations
    return [frozenset(c) for k in (1, 2, 3)
            for c in combinations(BIG6, k)]


def pool_for_D(D, cache):
    """union of gamma cells over all S families -> {triple} (both signs),
    h(gamma) <= H_POOL, <= POOL_CAP rows/cell (net truncation)."""
    G.H_CAP, G.QUAD_CAP, G.STATE_CEIL = H_POOL, POOL_CAP, 15000
    pool = {}
    for S in families():
        ck = ("cell", S, D)
        if ck not in cache:
            cache[ck] = G.quad_cell(sorted(S), D, cache)
        rows, _meta = cache[ck]
        for h, key, _prov in rows[:POOL_CAP]:
            if key not in pool:
                t = key_to_triple(key, D)
                if t is not None:
                    pool[key] = t
    return list(pool.values())


def check_quad(z, D, sqD, factors):
    """family-free acceptance for quadratic z (triple).  -> (S_min, h,
    unit, key) or None."""
    a, b, c = z
    if b == 0:
        return None
    h = G.q_height(z, D, sqD)
    if h > H_QUAD:
        return None
    A, C = c * c, a * a - D * b * b
    lo, sup = strip_p(A)
    if lo != 1 or C == 0:
        return None
    lo, s2 = strip_p(C)
    if lo != 1:
        return None
    # raw (unreduced) lead/const supports — conservative (never under-counts)
    sup = sup | s2
    nz = Fr(C, A)
    for m in factors:
        d = len(m) - 1
        val = (m[d], 0, 1)
        for i in range(d - 1, -1, -1):
            val = G.q_mul(val, z, D)
            val = q_add(val, (m[i], 0, 1), D)
        nv = G.q_norm(val, D)
        if nv == 0:
            return None          # z is a singular point
        lo, s3 = strip_p(nv.numerator)
        if lo != 1:
            return None
        lo, s4 = strip_p(nv.denominator)
        if lo != 1:
            return None
        sup |= s3 | s4
    if not allowed_smin(sup):
        return None
    return sup, h, abs(nz) == 1, G.q_minpoly_key(z, D)

# ---------------------------------------------------- assembly and gate
def assemble_op(op, card, cache):
    factors = sing_factors(card)
    rat_sing = [(-m[0], m[1]) for m in factors if len(m) == 2]
    rows = []
    for num, den, h in rational_net():
        for sgn in (1, -1):
            sm = check_rational(num, den, sgn, factors)
            if sm is not None:
                f = Fr(sgn * num, den)
                rows.append((not (num == den == 1), len(sm), h,
                             f"rat:{f}", {"z": str(f), "kind": "rational",
                             "S_min": sorted(sm), "h": round(h, 4),
                             "unit": num == den == 1}))
    seen = {r[3] for r in rows}
    for D in G.SQFREE_D:
        sqD = math.sqrt(D)
        if ("pool", D) not in cache:
            cache[("pool", D)] = pool_for_D(D, cache)
        pool = cache[("pool", D)]
        anchors = [(0, 1)] + [(p, q) for p, q in rat_sing]
        cands = {t for t in pool}
        for p, q in rat_sing:
            s = (p, 0, q) if q > 0 else (-p, 0, -q)
            for t in pool:
                cands.add(q_add(s, t, D))
        for m in factors:              # quadratic sing anchors in this field
            if len(m) == 3:
                dd, e = divmod_sq(m[1] * m[1] - 4 * m[2] * m[0], D)
                if dd:
                    for sg in (1, -1):
                        s = G.q_reduce(-m[1], sg * e, 2 * m[2])
                        for t in pool:
                            cands.add(q_add(s, t, D))
        for z in cands:
            r = check_quad(z, D, sqD, factors)
            if r is None:
                continue
            sm, h, unit, key = r
            mk = f"alg:{key}"
            if mk in seen:
                continue
            seen.add(mk)
            a, b, c = z
            rows.append((not unit, len(sm), h, mk,
                         {"z": f"({a}{'+' if b >= 0 else ''}{b}*sqrt{D})/{c}",
                          "D": D, "kind": "quadratic",
                          "S_min": sorted(sm), "h": round(h, 4),
                          "unit": unit}))
    rows.sort(key=lambda r: r[:4])
    known = cache["known"].get(op, set())
    fp = os.path.join(_outdir(), f"{op}.jsonl")
    nk = 0
    with open(fp, "w") as f:
        for _u, _ns, _h, mk, prov in rows:
            row = {"op": op, "match_key": mk, **prov, "method": "sunit",
                   "law": "refined-conductor-boxS", "status": "GENERATED"}
            if mk in known:
                row["dedup"] = "KNOWN-ORBIT-REDISCOVERY"
                nk += 1
            f.write(json.dumps(row) + "\n")
    st = {"rows": len(rows), "rational": sum(1 for r in rows if
          r[4]["kind"] == "rational"), "known_orbit_strips": nk,
          "factors": [list(m) for m in factors],
          "blowup_alarm": len(rows) > BLOWUP}
    st["run_rows_after_dedup"] = st["rows"] - nk
    return st

GATE_POINTS = [("T1", "AESZ34", "rat:-1/7"), ("T2pair", "AESZ34",
               "alg:1,-66,1"), ("T3", "AESZ22", "rat:-1"),
               ("T4", "AESZ17", "rat:-1"), ("X118", "AESZ118", "rat:-1/32"),
               ("X290", "AESZ290", "rat:1/729"),
               ("AESZ4", "h33", "rat:-1/5832"),      # BCM control points
               ("AESZ11", "h34", "rat:-1/432")]


def run_gate():
    rows, npass = [], 0
    for name, op, mk in GATE_POINTS:
        fp = os.path.join(_outdir(), f"{op}.jsonl")
        hit = None
        for ln in open(fp):
            r = json.loads(ln)
            if r["match_key"] == mk:
                hit = r
                break
        row = {"gate": name, "op": op, "target": mk,
               "law": "refined-conductor-boxS",
               "status": "PASS" if hit else "FAIL"}
        if hit:
            row["S_min"] = hit["S_min"]
            npass += 1
        rows.append(row)
    with open(os.path.join(_outdir(), "GATE_REFINED.jsonl"), "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    print(f"GATE(refined): {npass}/{len(GATE_POINTS)} PASS")
    for r in rows:
        print(" ", r["gate"], r["status"], r.get("S_min"))
    return npass == len(GATE_POINTS)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "gate"
    shard = int(os.environ.get("TERRIER_SHARD", "0"))
    nshard = int(os.environ.get("TERRIER_NSHARD", "1"))
    os.makedirs(_outdir(), exist_ok=True)
    cards, _al = G.load_cards()
    cache = {"known": G.load_known()}
    if mode == "gate":
        ops = [o for o in GATE_OPS if o in cards]
    elif mode == "census":
        ops = [o for o in SCOPE21 if o in cards]
    else:
        ops = sorted(cards)
    ops = [o for i, o in enumerate(ops) if i % nshard == shard]
    sizes = {}
    for op in ops:
        t0 = time.time()
        sizes[op] = assemble_op(op, cards[op], cache)
        print(f"{op}: rows={sizes[op]['rows']} "
              f"(rat={sizes[op]['rational']}) strips="
              f"{sizes[op]['known_orbit_strips']} "
              f"alarm={sizes[op]['blowup_alarm']} "
              f"({time.time() - t0:.0f}s)", flush=True)
    tag = mode + (f".s{shard}" if nshard > 1 else "")
    json.dump({"law": "refined-conductor-boxS",
               "pilot_truncations": {"A_MAX": A_MAX, "H_RAT": H_RAT,
               "RAT_SUPP_MAX": RAT_SUPP_MAX, "H_POOL": H_POOL,
               "POOL_CAP": POOL_CAP, "H_QUAD": H_QUAD,
               "quad_gen_families": "S in box, |S|<=3 (net only)"},
               "ops": sizes}, open(os.path.join(
                   _outdir(), f"SIZES2_{tag}.json"), "w"), indent=1)
    if mode in ("gate", "census") and nshard == 1:
        run_gate()


GATE_OPS = ["AESZ34", "AESZ22", "AESZ17", "AESZ118", "AESZ290",
            "h33", "h34"]

if __name__ == "__main__":
    main()

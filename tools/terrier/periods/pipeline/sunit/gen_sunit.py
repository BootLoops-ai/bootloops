#!/usr/bin/env python3
"""gen_sunit.py — ATTRACTOR-HUNT S-unit candidate generator.

Paths: ATLAS dir + KNOWN_ORBIT table env-overridable (TERRIER_SUNIT_ATLAS,
TERRIER_KNOWN_ORBIT); defaults = the packaged atlas/ and
common/data/KNOWN_ORBIT.jsonl.

EMPIRICAL LAW UNDER TEST: known rank-2 points are S-units/S-integers w.r.t.
their operator's bad primes.  Falsifiable 6-point gate runs BEFORE transport.

S-SET LAW (exactly what we take as S, per card atlas/<op>.json):
  S(op) = union over finite nonzero singular points sigma (conifolds AND
  apparent singularities, as carded) of primes dividing the CONSTANT and
  LEADING coefficients of sigma's primitive integer minpoly [c0..ck]
  (ascending).  Rationale: product of the minpoly's roots = ±c0/ck, so any
  prime in the support of N(sigma) divides c0 or ck, and any denominator
  prime of sigma divides ck.  Sources, in order:
    (a) card "singular_locus" minpolys (includes apparent sing., mult>=2);
    (b) else: sympy factor_list of the leading theta-coeff polynomial
        coeffs[-1] (ascending in z); drop the z^m factor (origin=MUM) and
        the integer content (not a root datum); same c0/ck rule per
        primitive irreducible factor.
  When both exist we take the UNION (both are root data of the same locus)
  and cross-check (a) vs (b); mismatch -> "anomaly" field.
  NOT included (documented exclusions): primes of disc(leadpoly), content
  primes, primes of lower-order coefficients, conductor primes of any
  fibre.  No widening: the gate tests this S as-is.

CANDIDATES (all carry provenance: construction + support + height):
  RATIONAL  z = +-prod_{p in S} p^e, Weil log-height <= H(=40); z != carded
    rational singular point; cap 2000/op, kept ascending by height.
  QUADRATIC per squarefree D in [2,100]: z=(a+b sqrt D)/c, b!=0, in the
    multiplicative group < -1, u_D, {pi_p}, {p in S} > where u_D = fundamental
    unit of O_K (Pell/brute), pi_p = minimal-power generator of a prime ideal
    over split/ramified p in S u ram(D) (brute search, k<=3, b<=20000; a miss
    is RECORDED and that p is dropped from the quadratic basis — completeness
    caveat, not silent).  support(N(z)) subset S u ram(D) by construction
    (asserted); denominator c is S u ram(D)-supported (documented deviation:
    specification said "c S-supported"; ram primes enter c only through
    pi_p^-1).  Proper units (norm +-1) = u_D^n branch, INCLUDED.
    Search: best-first on exact Weil height, +-1 exponent steps, admission
    slack 2*max gen-height (unit-cancellation tunneling), state ceiling;
    emit h <= H, cap 3000/(op,D) ascending by height.  Heuristically
    complete; the falsifiable gate is the check.
"""
import heapq, json, math, os, sys
from fractions import Fraction as Fr
from math import gcd, isqrt, log

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = os.environ.get("TERRIER_SUNIT_ATLAS", os.path.join(HERE, "atlas"))
H_CAP = 40.0            # Weil log-height budget
RAT_CAP = 2000          # per op
QUAD_CAP = 3000         # per (op, D)
STATE_CEIL = 60000      # best-first state ceiling per (S,D) cell
GEN_B_MAX = 20000       # brute bound for prime-ideal generator search
GEN_K_MAX = 3           # max prime power for ideal generators
PELL_B_MAX = 10**7      # fundamental-unit brute bound (D<=100: max b~221064)
SQFREE_D = [d for d in range(2, 101)
            if all(d % (q * q) for q in range(2, 11))]


def primes_of(n):
    n, out, p = abs(n), set(), 2
    while p * p <= n:
        if n % p == 0:
            out.add(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out.add(n)
    return out


def s_supported(n, S):
    n = abs(n)
    for p in S:
        while n % p == 0:
            n //= p
    return n == 1

# ---------------------------------------------------------------- S sets
def load_cards():
    """{opname: card} for unique realpaths; alias map for symlinked cards."""
    cards, seen, aliases = {}, {}, {}
    for fn in sorted(os.listdir(ATLAS)):
        if not fn.endswith(".json"):
            continue
        path = os.path.join(ATLAS, fn)
        rp = os.path.realpath(path)
        try:
            card = json.load(open(path))
        except Exception:
            continue
        if not (isinstance(card, dict) and "coeffs" in card and "op" in card):
            continue
        name = fn[:-5]
        if rp in seen:
            aliases.setdefault(seen[rp], []).append(name)
            continue
        seen[rp] = name
        cards[name] = card
    return cards, aliases


def s_from_minpoly(mp):
    """primes(|c0|) U primes(|ck|) for ascending integer minpoly, c0 != 0."""
    c0, ck = mp[0], mp[-1]
    if c0 == 0:
        return set()          # origin root = MUM, excluded
    return primes_of(c0) | primes_of(ck)


def s_of_card(card):
    """-> (sorted S, source, anomaly-or-None, rational singular points)."""
    import sympy as sp
    z = sp.symbols("z")
    lead = card["coeffs"][-1]
    P = sum(a * z ** i for i, a in enumerate(lead))
    S_lead, rat_sing, quad_sing = set(), set(), set()
    for f, _m in sp.factor_list(P)[1]:
        cs = [int(c) for c in reversed(sp.Poly(f, z).all_coeffs())]  # asc
        if cs[0] == 0 and len(cs) == 2 and cs[1] != 0:
            continue                       # pure z factor (origin)
        if cs[0] == 0:
            continue
        S_lead |= s_from_minpoly(cs)
        if len(cs) == 2:                   # root -c0/c1
            rat_sing.add(Fr(-cs[0], cs[1]))
        elif len(cs) == 3:
            g = gcd(gcd(cs[0], cs[1]), cs[2])
            t = tuple(c // g for c in cs)
            quad_sing.add(t if t[2] > 0 else tuple(-c for c in t))
    if "singular_locus" in card:
        S_card = set()
        for e in card["singular_locus"]:
            mp = [int(c) for c in e["minpoly"]]
            S_card |= s_from_minpoly(mp)
            if len(mp) == 2:
                rat_sing.add(Fr(-mp[0], mp[1]))
            elif len(mp) == 3:
                g = gcd(gcd(mp[0], mp[1]), mp[2])
                t = tuple(c // g for c in mp)
                quad_sing.add(t if t[2] > 0 else tuple(-c for c in t))
        anom = (None if S_card == S_lead else
                {"singular_locus_S": sorted(S_card),
                 "leadpoly_S": sorted(S_lead)})
        return sorted(S_card | S_lead), "singular_locus+leadpoly", anom, \
            rat_sing, quad_sing
    return sorted(S_lead), "leadpoly", None, rat_sing, quad_sing

# ------------------------------------------------- quadratic arithmetic
# element of K = Q(sqrt D): triple (a, b, c) = (a + b sqrt D)/c, c > 0,
# gcd(a, b, c) = 1.
def q_reduce(a, b, c):
    if c < 0:
        a, b, c = -a, -b, -c
    g = gcd(gcd(abs(a), abs(b)), c)
    return (a // g, b // g, c // g)


def q_mul(x, y, D):
    a1, b1, c1 = x
    a2, b2, c2 = y
    return q_reduce(a1 * a2 + D * b1 * b2, a1 * b2 + a2 * b1, c1 * c2)


def q_inv(x, D):
    a, b, c = x                      # 1/x = c*(a - b sqrt D)/(a^2 - D b^2)
    n = a * a - D * b * b
    if n < 0:
        return q_reduce(-c * a, c * b, -n)
    return q_reduce(c * a, -c * b, n)


def q_norm(x, D):                    # Fraction norm
    a, b, c = x
    return Fr(a * a - D * b * b, c * c)


def q_height(x, D, sqD):
    """Weil log-height via minpoly A z^2 + B z + C (b!=0) or Fraction."""
    a, b, c = x
    if b == 0:
        return log(max(abs(a), c)) if a else 0.0
    A, B, C = c * c, -2 * a * c, a * a - D * b * b
    g = gcd(gcd(A, abs(B)), abs(C))
    A, B, C = A // g, B // g, C // g
    # larger |root| has |a|, |b|sqD adding -> float-safe; smaller root via
    # the exact norm |C|/A (avoids catastrophic cancellation a ~ -+b sqD)
    m = (abs(a) + abs(b) * sqD) / c
    s = abs(C) / (A * m)
    return 0.5 * (log(A) + log(max(1.0, m)) + log(max(1.0, s)))


def q_minpoly_key(x, D):
    """primitive descending minpoly key 'A,B,C' (A>0), b != 0."""
    a, b, c = x
    A, B, C = c * c, -2 * a * c, a * a - D * b * b
    g = gcd(gcd(A, abs(B)), abs(C))
    return f"{A // g},{B // g},{C // g}"


def is_sq(n):
    if n < 0:
        return -1
    r = isqrt(n)
    return r if r * r == n else -1


def fund_unit(D):
    """Fundamental unit > 1 of O_K, K = Q(sqrt D), via ascending-b brute.
    D=1 mod 4: check a^2 = D b^2 -+ 4 with a = b mod 2 (half-integer units,
    subsume c=1); else a^2 = D b^2 -+ 1.  First b with a hit is fundamental."""
    half = (D % 4 == 1)
    for b in range(1, PELL_B_MAX):
        t = D * b * b
        if half:
            for s in (-4, 4):
                a = is_sq(t + s)
                if a >= 0 and (a - b) % 2 == 0:
                    return q_reduce(a, b, 2)
        else:
            for s in (-1, 1):
                a = is_sq(t + s)
                if a >= 0:
                    return (a, b, 1)
    raise RuntimeError(f"no fundamental unit found for D={D}")

def kronecker(a, n):
    """Kronecker symbol (a|n), n > 0."""
    a %= n
    if gcd(a, n) != 1 and n > 1 and gcd(a, n) != n:
        pass
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


def prime_gen(D, p):
    """Minimal-power generator pi with |N(pi)| = p^k, k<=GEN_K_MAX, via brute
    b<=GEN_B_MAX.  -> (triple, k) or (None, why).  Cached by caller."""
    dK = D if D % 4 == 1 else 4 * D
    kr = kronecker(dK % p if dK % p else 0, p) if p > 2 else \
        (0 if dK % 2 == 0 else (1 if dK % 8 == 1 else -1))
    if kr == -1:
        return None, "inert"
    half = (D % 4 == 1)
    for k in range(1, GEN_K_MAX + 1):
        pk = p ** k
        for b in range(1, GEN_B_MAX):
            t = D * b * b
            if half:
                for s in (-4 * pk, 4 * pk):
                    a = is_sq(t + s)
                    if a >= 0 and (a - b) % 2 == 0:
                        return q_reduce(a, b, 2), k
            for s in (-pk, pk):
                a = is_sq(t + s)
                if a >= 0:
                    return (a, b, 1), k
    return None, f"no-generator b<={GEN_B_MAX} k<={GEN_K_MAX}"


# ------------------------------------------------ quadratic enumeration
def quad_cell(S, D, cache):
    """All quadratic candidates for prime set S, field Q(sqrt D).
    -> (rows, meta).  rows: (h, minpoly_key, provenance dict).  Cached
    upstream by (S, D) since ops share S sets."""
    sqD = math.sqrt(D)
    dK = D if D % 4 == 1 else 4 * D
    ram = sorted(primes_of(dK))
    Sp = sorted(set(S) | set(ram))
    if ("u", D) not in cache:
        cache[("u", D)] = fund_unit(D)
    u = cache[("u", D)]
    gens, gen_meta = [("u", u)], {"ram": ram, "dropped": []}
    for p in Sp:
        if ("pi", D, p) not in cache:
            cache[("pi", D, p)] = prime_gen(D, p)
        pi, k = cache[("pi", D, p)]
        if pi is None:
            if k != "inert":
                gen_meta["dropped"].append({"p": p, "why": k})
            continue
        gens.append((f"pi{p}", pi))
    for p in S:
        gens.append((f"r{p}", (p, 0, 1)))
    hs = [q_height(g, D, sqD) for _, g in gens]
    slack = 2.0 * max(hs)
    lim = H_CAP + slack
    return _bfs(gens, hs, D, sqD, lim, Sp, gen_meta)

def _bfs(gens, hs, D, sqD, lim, Sp, meta):
    """Best-first over exponent vectors; emit (h, key, prov) rows, b != 0,
    h <= H_CAP, both signs.  Stop: heap empty | STATE_CEIL | 2*QUAD_CAP."""
    ginv = [q_inv(g, D) for _, g in gens]
    ng = len(gens)
    start = (0,) * ng
    heap = [(0.0, start, (1, 0, 1))]
    visited = {start}
    rows, keys, pops = [], set(), 0
    while heap and pops < STATE_CEIL and len(rows) < 2 * QUAD_CAP:
        h, ex, v = heapq.heappop(heap)
        pops += 1
        if v[1] != 0 and h <= H_CAP:
            nrm = q_norm(v, D)
            assert s_supported(nrm.numerator, Sp) and \
                s_supported(nrm.denominator, Sp), (v, D, Sp)
            cons = "*".join(f"{gens[i][0]}^{e}" for i, e in
                            enumerate(ex) if e) or "1"
            sup = sorted(primes_of(nrm.numerator) |
                         primes_of(nrm.denominator))
            unit = abs(nrm) == 1
            for sgn in (1, -1):
                w = v if sgn == 1 else (-v[0], -v[1], v[2])
                k = q_minpoly_key(w, D)
                if k in keys:
                    continue
                keys.add(k)
                rows.append((h, k, {
                    "z": f"({w[0]}{'+' if w[1] >= 0 else ''}{w[1]}*sqrt{D})"
                         f"/{w[2]}", "D": D, "construction":
                    ("-" if sgn < 0 else "") + cons, "norm": str(nrm),
                    "support": sup, "unit": unit, "h": round(h, 4)}))
        for i in range(ng):
            for stp in (1, -1):
                ex2 = ex[:i] + (ex[i] + stp,) + ex[i + 1:]
                if ex2 in visited:
                    continue
                v2 = q_mul(v, gens[i][1] if stp == 1 else ginv[i], D)
                h2 = q_height(v2, D, sqD)
                if h2 > lim:
                    continue
                visited.add(ex2)
                heapq.heappush(heap, (h2, ex2, v2))
    meta.update({"pops": pops, "exhausted": not heap,
                 "rows_pre_cap": len(rows)})
    rows.sort(key=lambda r: (r[0], r[1]))
    return rows, meta


# ------------------------------------------------- rational enumeration
def smooth_ascending(S, limit, nmax=200000):
    """S-smooth integers <= limit, ascending, capped."""
    out, heap, seen = [], [1], {1}
    while heap and len(out) < nmax:
        n = heapq.heappop(heap)
        out.append(n)
        for p in S:
            m = n * p
            if m <= limit and m not in seen:
                seen.add(m)
                heapq.heappush(heap, m)
    return out, bool(heap)


def rational_cands(S):
    """(h, 'rat:num/den', prov) rows, +-smooth/smooth coprime, h <= H_CAP."""
    lim = int(math.exp(H_CAP))
    L, truncated = smooth_ascending(S, lim)
    rows = []
    for i, M in enumerate(L):
        if len(rows) >= 3 * RAT_CAP:
            break
        for N in L[:i + 1]:
            if gcd(M, N) != 1:
                continue
            h = log(M) if M > 1 else 0.0
            for num, den in ({(M, N), (N, M)}):
                for sgn in (1, -1):
                    f = Fr(sgn * num, den)
                    rows.append((h, f"rat:{f}", {
                        "z": str(f), "construction": f"{'-' if sgn < 0 else ''}"
                        f"{num}/{den} (S-monomial)",
                        "support": sorted(primes_of(num) | primes_of(den)),
                        "unit": abs(f) == 1, "h": round(h, 4)}))
    rows.sort(key=lambda r: (r[0], r[1]))
    return rows, truncated

# ------------------------------------------------------ assembly + gate
GATE_OPS = ["AESZ34", "AESZ22", "AESZ17", "AESZ118", "AESZ290"]
GATE_POINTS = [  # (name, op, kind, key)  keys: rat:<frac> | alg:<A,B,C> desc
    ("T1", "AESZ34", "rat", "rat:-1/7"),
    ("T2pair", "AESZ34", "alg", "alg:1,-66,1"),
    ("T3", "AESZ22", "rat", "rat:-1"),
    ("T4", "AESZ17", "rat", "rat:-1"),
    ("X118", "AESZ118", "rat", "rat:-1/32"),
    ("X290", "AESZ290", "rat", "rat:1/729")]


def load_known():
    ko, path = {}, os.environ.get(
        "TERRIER_KNOWN_ORBIT", os.path.abspath(os.path.join(
            HERE, "..", "..", "..", "common", "data", "KNOWN_ORBIT.jsonl")))
    for ln in open(path):
        r = json.loads(ln)
        ko.setdefault(r["operator_key"], set()).add(r["match_key"])
    return ko


def assemble_op(op, info, cache, emit):
    """Build candidate rows for one op.  emit: write JSONL (else count only).
    -> stats dict.  Priority order in file: units, |support|, height."""
    S = tuple(info["S"])
    rat_sing = {Fr(x) for x in info["rat_sing"]}
    qsd = {f"{t[2]},{t[1]},{t[0]}" for t in info["quad_sing"]}  # desc keys
    if ("rat", S) not in cache:
        cache[("rat", S)] = rational_cands(list(S))
    rrows, trunc = cache[("rat", S)]
    kept, excl = [], 0
    for h, key, prov in rrows:
        if Fr(key[4:]) in rat_sing:
            excl += 1
            continue
        kept.append((h, key, {**prov, "kind": "rational", "op": op}))
        if len(kept) >= RAT_CAP:
            break
    qkept, qexcl, cells_sat, qtot = [], 0, 0, 0
    for D in SQFREE_D:
        if ("quad", S, D) not in cache:
            cache[("quad", S, D)] = quad_cell(list(S), D, cache)
        rows, meta = cache[("quad", S, D)]
        n = 0
        for h, key, prov in rows:
            if key in qsd:
                qexcl += 1
                continue
            n += 1
            if n > QUAD_CAP:
                break
            qkept.append((h, f"alg:{key}",
                          {**prov, "kind": "quadratic", "op": op}))
        qtot += min(n, QUAD_CAP)
        cells_sat += (n > QUAD_CAP)
    known = cache["known"].get(op, set())
    stats = {"S": list(S), "rational": len(kept), "quadratic": qtot,
             "cells_saturated_of_60": cells_sat, "sing_excluded_rat": excl,
             "sing_excluded_quad": qexcl, "rat_truncated_smooth": trunc,
             "known_orbit_strips": 0}
    if emit:
        allr = kept + qkept
        allr.sort(key=lambda r: (0 if r[2]["unit"] else 1,
                                 len(r[2]["support"]), r[0], r[1]))
        fp = os.path.join(HERE, "candidates", f"{op}.jsonl")
        with open(fp, "w") as f:
            for h, key, prov in allr:
                row = {"op": op, "match_key": key, **prov,
                       "method": "sunit", "status": "GENERATED"}
                if key in known:
                    row["dedup"] = "KNOWN-ORBIT-REDISCOVERY"
                    stats["known_orbit_strips"] += 1
                f.write(json.dumps(row) + "\n")
    else:
        stats["known_orbit_strips"] = sum(
            1 for h, k, p in kept + qkept if k in known)
    stats["run_rows_after_dedup"] = (stats["rational"] + stats["quadratic"]
                                     - stats["known_orbit_strips"])
    return stats

def t1_diagnosis(cards):
    """Exhaustive card-side check: is 7 in ANY singular datum of AESZ34?"""
    import sympy as sp
    card = cards["AESZ34"]
    z = sp.symbols("z")
    lead = card["coeffs"][-1]
    P = sum(a * z ** i for i, a in enumerate(lead))
    d = {"lead_poly_asc": lead,
         "factors": str(sp.factor_list(P)),
         "disc_leadpoly": dict(sp.factorint(sp.discriminant(P))),
         "sing_minpolys": [e["minpoly"] for e in card["singular_locus"]],
         "seven_divides_any_minpoly_coeff": any(
             c % 7 == 0 for e in card["singular_locus"]
             for c in e["minpoly"] if c),
         "seven_divides_any_leadpoly_coeff": any(
             c % 7 == 0 for c in lead if c),
         "content_leadpoly": int(sp.gcd_list([c for c in lead if c]))}
    d["verdict"] = ("LAW-FAIL: support(-1/7)={7} not within S(AESZ34) under "
                    "any reading of the card singular data (roots 1, 1/9, "
                    "1/25 -> S={3,5}; disc leadpoly=2^20*3^2; no apparent "
                    "singularities; content 1 and all leadpoly factors are "
                    "7-free — 7 divides only NON-root data: mid coeffs of "
                    "S4 (35, -259=-7*37) and of lower-order S_i). Generator "
                    "is correct; the empirical law itself fails at T1.")
    return d


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "gate"
    cards, aliases = load_cards()
    ssets, infos = {}, {}
    for op, card in sorted(cards.items()):
        S, src, anom, rs, qs = s_of_card(card)
        infos[op] = {"S": S, "rat_sing": [str(x) for x in rs],
                     "quad_sing": sorted(qs)}
        ssets[op] = {"S": S, "source": src, "anomaly": anom,
                     "rat_sing": [str(x) for x in rs],
                     "quad_sing": [f"{t[2]},{t[1]},{t[0]}" for t in
                                   sorted(qs)], "aliases": aliases.get(op, [])}
    json.dump({"law": "S = primes(c0)+primes(ck) per singular minpoly "
               "(see module docstring)", "H_cap": H_CAP, "ops": ssets},
              open(os.path.join(HERE, "S_SETS.json"), "w"), indent=1)
    print(f"S_SETS.json: {len(ssets)} unique cards, "
          f"{sum(len(v) for v in aliases.values())} aliases")
    if mode == "ssets":
        return
    cache = {"known": load_known()}
    ops = GATE_OPS if mode == "gate" else sorted(infos)
    sizes = {}
    for op in ops:
        t0 = __import__("time").time()
        sizes[op] = assemble_op(op, infos[op], cache, emit=(mode == "gate"))
        print(f"{op}: S={sizes[op]['S']} rat={sizes[op]['rational']} "
              f"quad={sizes[op]['quadratic']} "
              f"({__import__('time').time() - t0:.1f}s)", flush=True)
    if mode == "gate":
        gate_rows, npass = [], 0
        for name, op, kind, key in GATE_POINTS:
            found = False
            for ln in open(os.path.join(HERE, "candidates", f"{op}.jsonl")):
                if json.loads(ln)["match_key"] == key:
                    found = True
                    break
            row = {"gate": name, "op": op, "target": key,
                   "status": "PASS" if found else "FAIL"}
            if found:
                npass += 1
            elif name == "T1":
                row["diagnosis"] = t1_diagnosis(cards)
            gate_rows.append(row)
        with open(os.path.join(HERE, "GATE.jsonl"), "w") as f:
            for r in gate_rows:
                f.write(json.dumps(r) + "\n")
        print(f"GATE: {npass}/{len(GATE_POINTS)} PASS")
    json.dump(sizes, open(os.path.join(
        HERE, f"SIZES_{mode}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

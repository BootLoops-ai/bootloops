#!/usr/bin/env python3
# bench_rung1_gauss.py -- the certified Gauss-Jacobi rule builder, reference
# sweep, and exact-rational N=8 surrogate ("bench rung 1" = the surrogate-scale
# tier of the engine).
#
# (i) Certified Gauss-Jacobi rule builder (python-flint/arb): interval-Newton
#     polish of scipy seeds; weight u^(c0-1) on [0,1] at DYADIC p0 = 709/2048
#     (c0 = 2*lam*p0*(1-p0), lam = 10); n = (637, 478) for u1,u2 and the
#     y3-parity-reduced v = u3^2 rule n = 160 with weight v^(c0/2-1).
#     Acceptance: all node balls pairwise disjoint, all weight balls strictly
#     positive, sum(weights) reproduces the exact moment 1/(alpha+1) to ball
#     width; PLUS the kill round's endpoint-cluster node check.
# (ii) ONE full certified production sweep: N=318 real pattern counts at p0
#     through the 637 x 478 x 160 tensor rule (streamed, never stored) ->
#     a single ball Z_s(p0) = c0^3 * (1/2) * sum w1 w2 w3 F(u1,u2,v3),
#     F = prod_x P_x^{n_x} over the 14 observed patterns (RAW positive product,
#     per-node value evaluation A_x + B_x*v3 -- NO coefficient expansion).
#     Acceptance: relwidth <= 2^-150 at prec 192; magnitude consistent with the
#     raw-product reference class (float-shadow diagnostic); the N<=8
#     exact-rational cross-check reruns inside the certified harness to
#     ball-width agreement (+ planted-count control).
# KILL: any weight ball straddling 0; sweep relwidth worse than 2^-100;
#     endpoint cancellation blowing 192 bits.
#
# Write law: writes ONLY the RUNG1_GAUSS_SWEEP.json receipt + scratch dumps.
import json, time, sys, math, argparse, resource, os, tempfile
from fractions import Fraction
from flint import arb, ctx

import os as _os
BASE = _os.environ.get("POSQ_STAGE0", _os.path.dirname(_os.path.abspath(__file__)))
# Receipt/scratch roots are env-configured.
# RESULTS = optional dir of original production inputs/results (not shipped).
RESULTS = _os.environ.get("POSQ_RESULTS_DIR", "")
SCRATCH = _os.environ.get("POSQ_SCRATCH") or _os.path.join(
    tempfile.gettempdir(), "posq_bench_rungs")
_os.makedirs(SCRATCH, exist_ok=True)
# The production pattern-count data is NOT shipped with this package (only
# the N=8 surrogate paths, which need no data file, are exercised by the
# battery). Point POSQ_STAGE0 at a directory containing QUARTET_COLLAPSE.json
# to enable the production sweep commands; without it they refuse by name.
_QC_PATH = f"{BASE}/QUARTET_COLLAPSE.json"
if _os.path.exists(_QC_PATH):
    QC = json.load(open(_QC_PATH))
    PATS = [(tuple(int(ch) for ch in p), int(n))
            for p, n in QC["pattern_counts"].items()]
    assert sum(n for _, n in PATS) == 318 and len(PATS) == 14
else:
    QC, PATS = None, None


def _require_pats():
    if PATS is None:
        raise SystemExit(
            "REFUSE: POSQ-PRODUCTION-DATA-ABSENT -- the production "
            "pattern-count file QUARTET_COLLAPSE.json is not shipped; point "
            "POSQ_STAGE0 at a directory containing it to run production "
            "sweep commands (the surrogate/battery paths do not need it)")
    return PATS

LAM = 10
P0_NUM, P0_DEN = 709, 2048                    # dyadic p0 = 0.34619140625
C0 = Fraction(2 * LAM * P0_NUM * (P0_DEN - P0_NUM), P0_DEN * P0_DEN)
ALPHA1 = C0 - 1                               # weight u^(c0-1), dims 1,2
ALPHA3 = C0 / 2 - 1                           # weight v^(c0/2-1), parity dim
N1, N2, N3 = 637, 478, 160                    # exact for degrees (1272, 954, 318)
PREC_BUILD = 2048
PREC_SWEEP = 192

def fr2arb(q):
    return arb(q.numerator) / arb(q.denominator)

def man_exp_pair(x):
    m = x.mid(); r = x.rad()
    mm, me = m.man_exp(); rm, re = r.man_exp()
    return [str(mm), str(me), str(rm), str(re)]

def from_man_exp(q):
    mm, me, rm, re = int(q[0]), int(q[1]), int(q[2]), int(q[3])
    return arb(mm) * arb(2) ** me + (arb(rm) * arb(2) ** re) * arb(0, 1)

def lb_of(x):
    l = x.mid() - x.rad()
    return l - l.rad()

def ub_of(x):
    u = x.mid() + x.rad()
    return u + u.rad()

# ---------------- certified Gauss-Jacobi rule on [0,1], weight u^alpha ----------
def recur_coeffs(n, alpha):
    """monic three-term recurrence for weight u^alpha on [0,1] (shifted Jacobi
    a=0, b=alpha; Gautschi r_jacobi mapped by u=(x+1)/2): pi_{k+1} =
    (x - A[k]) pi_k - B[k] pi_{k-1}; B[0] = mu0 = 1/(alpha+1).
    All coefficients exact rationals -> arb balls at working prec."""
    al = fr2arb(alpha)
    A = []; B = []
    for k in range(n):
        if k == 0:
            ak = al / (al + 2)
        else:
            ak = al * al / ((2 * k + al) * (2 * k + al + 2))
        A.append((ak + 1) / 2)
    B.append(1 / (al + 1))
    for k in range(1, n + 1):
        if k == 1:
            bk = 4 * (al + 1) / ((al + 2) ** 2 * (al + 3))
        else:
            bk = 4 * k * k * (k + al) ** 2 / ((2 * k + al) ** 2
                                              * (2 * k + al + 1) * (2 * k + al - 1))
        B.append(bk / 4)
    return A, B

def poly_dp(x, n, A, sb, p0t):
    """(ptilde_n(x), ptilde_n'(x)) via the ORTHONORMAL recurrence
    b_{k+1} p_{k+1} = (x - a_k) p_k - b_k p_{k-1}, b_k = sqrt(B_k):
    values stay O(poly(n)) -- the monic recurrence underflows (4^-n decay
    below the 256-bit radius floor, measured). Same roots as monic pi_n."""
    pm = arb(0); dm = arb(0)
    p = p0t; d = arb(0)
    for k in range(n):
        inv = sb[k + 1]
        pn = ((x - A[k]) * p - (sb[k] * pm if k > 0 else arb(0))) / inv
        dn = (p + (x - A[k]) * d - (sb[k] * dm if k > 0 else arb(0))) / inv
        pm, dm, p, d = p, d, pn, dn
    return p, d

def build_rule(n, alpha, tag):
    """certified nodes/weights; returns dict with balls + receipts."""
    from scipy.special import roots_jacobi
    t0 = time.time()
    ctx.prec = PREC_BUILD
    A, B = recur_coeffs(n, alpha)
    sb = [B[k].sqrt() for k in range(n + 1)]
    p0t = 1 / sb[0]
    xs, _ = roots_jacobi(n, 0.0, float(alpha))
    seeds = [(x + 1) / 2 for x in xs]
    nodes = []
    RCERT = arb(2) ** (-1000)     # interval-Newton input radius (endpoint-cluster
    # amplification measured ~2^830 -- input radius must sit far below it)
    for i, s in enumerate(seeds):
        # stage A: POINT Newton polish at prec 2048 (never certifies -- seeds only)
        x = arb(s)
        for it in range(10):
            p, d = poly_dp(x, n, A, sb, p0t)
            if not (abs(d) > 0):
                raise RuntimeError(f"rule {tag}: polish derivative ~0 at node {i}")
            step = p / d
            x = (x - step).mid()
            if bool(abs(step) < arb(2) ** (-1150)):
                break
        # stage B: one interval-Newton step certifies existence + uniqueness
        ok = False
        for attempt in range(3):
            X = x + arb(0, 1) * (RCERT * arb(2) ** (60 * attempt))
            mid = X.mid()
            pm, _ = poly_dp(arb(mid), n, A, sb, p0t)
            _, dX = poly_dp(X, n, A, sb, p0t)
            if dX.contains(arb(0)):
                continue
            Nx = mid - pm / dX
            if (lb_of(Nx) > lb_of(X)) and (ub_of(Nx) < ub_of(X)):
                nodes.append(Nx)
                ok = True
                break
        if not ok:
            raise RuntimeError(f"rule {tag}: Newton certification failed node {i}")
    # disjointness + range
    nodes.sort(key=lambda b: float(b.mid()))
    min_gap = None
    for i in range(n - 1):
        gap = lb_of(nodes[i + 1]) - ub_of(nodes[i])
        if not (gap > 0):
            raise RuntimeError(f"rule {tag}: node balls {i},{i+1} not disjoint")
        if min_gap is None or float(gap) < min_gap:
            min_gap = float(gap)
    if not (lb_of(nodes[0]) > 0) or not (ub_of(nodes[-1]) < 1):
        raise RuntimeError(f"rule {tag}: node outside (0,1)")
    # weights: w_i = 1 / sum_{k<n} ptilde_k(x_i)^2 (orthonormal recurrence)
    weights = []
    for X in nodes:
        s = p0t * p0t
        pk = p0t; pkm = arb(0)
        for k in range(n - 1):
            pn = ((X - A[k]) * pk - (sb[k] * pkm if k > 0 else arb(0))) / sb[k + 1]
            s = s + pn * pn
            pkm, pk = pk, pn
        w = 1 / s
        if not (w > 0):
            raise RuntimeError(f"rule {tag}: weight ball not strictly positive")
        weights.append(w)
    # receipts
    mu0 = fr2arb(Fraction(1, 1) / (alpha + 1))
    sw = arb(0)
    for w in weights:
        sw = sw + w
    dsw = sw - mu0
    mono = {}
    for m in (1, 2, n, 2 * n - 2, 2 * n - 1):
        sm = arb(0)
        for X, w in zip(nodes, weights):
            sm = sm + w * X ** m
        ex = fr2arb(Fraction(1, 1) / (alpha + m + 1))
        d = sm - ex
        mono[str(m)] = {"absdiff_ub": float(ub_of(abs(d))),
                        "contains_exact": bool(d.contains(arb(0)) or abs(d) < sm.rad() * 4)}
    rec = {
        "tag": tag, "n": n, "alpha": str(alpha), "prec_build": PREC_BUILD,
        "min_node_gap": min_gap,
        "max_node_rad": max(float(x.rad()) for x in nodes),
        "min_weight_lb": float(lb_of(min(weights, key=lambda w: float(w.mid())))),
        "max_weight_relrad": max(float(w.rad() / w.mid()) for w in weights),
        "all_weights_positive": True, "nodes_pairwise_disjoint": True,
        "sum_w_minus_moment_ub": float(ub_of(abs(dsw))),
        "sum_w_contains_moment": bool(dsw.contains(arb(0))),
        "monomial_exactness": mono,
        "endpoint_cluster": {"u_min": float(nodes[0].mid()),
                             "one_minus_u_max": float(1 - nodes[-1].mid()),
                             "n^-2": 1.0 / n ** 2},
        "build_s": round(time.time() - t0, 2),
    }
    return {"nodes": nodes, "weights": weights, "receipts": rec}

def rules_cmd():
    t0 = time.time()
    out = {}
    for (n, al, tag) in ((N1, ALPHA1, "u1"), (N2, ALPHA1, "u2"), (N3, ALPHA3, "v3")):
        r = build_rule(n, al, tag)
        out[tag] = r
        print(json.dumps(r["receipts"]), flush=True)
        ok = (r["receipts"]["sum_w_contains_moment"]
              and all(v["contains_exact"] for v in r["receipts"]["monomial_exactness"].values()))
        if not ok:
            raise RuntimeError(f"rule {tag}: exactness receipts FAILED")
    dump = {tag: {"receipts": out[tag]["receipts"],
                  "nodes": [man_exp_pair(x) for x in out[tag]["nodes"]],
                  "weights": [man_exp_pair(w) for w in out[tag]["weights"]]}
            for tag in out}
    json.dump(dump, open(f"{SCRATCH}/rung1_rules.json", "w"))
    print(f"[rules] built+certified in {time.time()-t0:.1f}s", flush=True)
    return 0

def load_rules():
    d = json.load(open(f"{SCRATCH}/rung1_rules.json"))
    prec_save = ctx.prec
    ctx.prec = PREC_BUILD + 64        # exact reconstruction of 2048-bit dyadics
    out = {}
    for tag in d:
        out[tag] = {"nodes": [from_man_exp(q) for q in d[tag]["nodes"]],
                    "weights": [from_man_exp(q) for q in d[tag]["weights"]],
                    "receipts": d[tag]["receipts"]}
    ctx.prec = prec_save
    return out

# ---------------- per-node evaluator: P_x = A_x + B_x * v3 ----------------
def ab_fiber(u1, u2, pi0, pi1, pats):
    """per-fiber A_x, B_x lists for pattern list pats (bit tuples)."""
    one = arb(1)
    m1 = one - u1; m2 = one - u2
    u12 = u1 * u2; m12 = one - u12
    P1 = ((pi0 + pi1 * u1, pi1 * m1), (pi0 * m1, pi1 + pi0 * u1))
    P2 = ((pi0 + pi1 * u2, pi1 * m2), (pi0 * m2, pi1 + pi0 * u2))
    P12 = ((pi0 + pi1 * u12, pi1 * m12), (pi0 * m12, pi1 + pi0 * u12))
    g = [[None, None], [None, None]]
    gg = [[[None, None], [None, None]], [[None, None], [None, None]]]
    for w in (0, 1):
        for a in (0, 1):
            for b in (0, 1):
                gg[w][a][b] = (P2[w][0] * P1[0][a] * P1[0][b]
                               + P2[w][1] * P1[1][a] * P1[1][b])
    pi = (pi0, pi1)
    A = []; Bv = []
    for (a, b, c, d) in pats:
        base0 = pi0 * gg[0][a][b] * P12[0][c]
        base1 = pi1 * gg[1][a][b] * P12[1][c]
        A.append((base0 + base1) * pi[d])
        Bv.append((base0 * ((1 if d == 0 else 0) - pi[d])
                   + base1 * ((1 if d == 1 else 0) - pi[d])) * u12)
    return A, Bv

def anchors_cmd():
    """controls: (a) A/B evaluator vs the reference integrand.py pattern_prob at a
    random rational point (16 patterns); (b) T1/T2 50-digit uncorrected
    anchors through the A/B path; (c) planted count mutation FIRES."""
    _require_pats()
    ctx.prec = 256
    sys.path.insert(0, BASE)
    import integrand as ig
    ok_all = True
    lines = [f"== RUNG1 evaluator gate ({time.strftime('%Y-%m-%d %H:%M:%S %Z')}) =="]
    def chk(name, ok, got=""):
        nonlocal ok_all
        lines.append(f"{'PASS' if ok else 'FAIL'}  {name}  {got}")
        ok_all = ok_all and ok
    allpats = [tuple(int(ch) for ch in f"{i:04b}") for i in range(16)]
    # (a) structural check at rational point
    pi1 = arb(2) / 5; pi0 = 1 - pi1
    s1, s2, s3 = arb(1) / 2, arb(1) / 4, arb(1) / 8
    beta = 1 / (2 * pi0 * pi1)
    u1 = (-beta * s1).exp(); u2 = (-beta * s2).exp(); u3 = (-beta * s3).exp()
    A, Bv = ab_fiber(u1, u2, pi0, pi1, allpats)
    v3 = u3 * u3
    worst = arb(0)
    for i, x in enumerate(allpats):
        direct = ig.pattern_prob(pi1, s1, s2, s3, x)
        mine = A[i] + Bv[i] * v3
        d = abs(direct - mine)
        if not (d < worst):
            worst = d
    chk("A/B evaluator == reference pattern_prob (16 pats) < 1e-60",
        bool(worst < arb("1e-60")), worst.str(3))
    # (b) T1/T2 uncorrected anchors through the A/B path
    for (name, a, b, cc, p), (corr50, unc50) in ig.ANCHORS.items():
        s1, s2, s3 = ig.frac(a), ig.frac(b), ig.frac(cc)
        pi1 = ig.frac(p); pi0 = 1 - pi1
        beta = 1 / (2 * pi0 * pi1)
        u1 = (-beta * s1).exp(); u2 = (-beta * s2).exp(); u3 = (-beta * s3).exp()
        A, Bv = ab_fiber(u1, u2, pi0, pi1, [x for x, _ in PATS])
        v3 = u3 * u3
        ll = arb(0)
        for i, (_, n) in enumerate(PATS):
            ll = ll + n * (A[i] + Bv[i] * v3).log()
        hit = bool(abs(ll - arb(unc50)) < arb("1e-49"))
        chk(f"{name} uncorrected 50-digit anchor via A/B path", hit, ll.str(45))
    # (c) planted count mutation fires
    pats_bad = [((x, n + 1) if x == (0, 0, 0, 1) else (x, n)) for x, n in PATS]
    ll_bad = arb(0)
    for i, (x, n) in enumerate(pats_bad):
        ll_bad = ll_bad + n * (A[i] + Bv[i] * v3).log()
    chk("planted count mutation FIRES (T2 anchor moves >> width)",
        bool(abs(ll_bad - ll) > arb(1)), f"shift {float((ll_bad-ll).mid()):.4f}")
    out = "\n".join(lines + [f"RUNG1 EVALUATOR GATE: {'PASS' if ok_all else 'FAIL'}"])
    print(out, flush=True)
    with open(f"{SCRATCH}/rung1_gate.txt", "a") as f:
        f.write(out + "\n")
    return 0 if ok_all else 1

# ---------------- generic certified tensor sweep (streamed) ----------------
def sweep_range(rules, pats, i0, i1, prec):
    """partial sum over u1 nodes [i0,i1): sum_i w1 sum_j w2 sum_k w3 F."""
    ctx.prec = prec
    pi1 = arb(P0_NUM) / P0_DEN; pi0 = 1 - pi1
    x1, w1 = rules["u1"]["nodes"], rules["u1"]["weights"]
    x2, w2 = rules["u2"]["nodes"], rules["u2"]["weights"]
    x3, w3 = rules["v3"]["nodes"], rules["v3"]["weights"]
    pat_bits = [x for x, _ in pats]
    exps = [n for _, n in pats]
    npat = len(pats)
    tot = arb(0)
    nodes_done = 0
    t0 = time.time(); c0 = resource.getrusage(resource.RUSAGE_SELF)
    for i in range(i0, i1):
        u1 = x1[i]
        row = arb(0)
        for j in range(len(x2)):
            A, Bv = ab_fiber(u1, x2[j], pi0, pi1, pat_bits)
            fib = arb(0)
            for k in range(len(x3)):
                v = x3[k]
                F = None
                for t in range(npat):
                    P = A[t] + Bv[t] * v
                    n = exps[t]
                    Pn = P if n == 1 else P ** n
                    F = Pn if F is None else F * Pn
                fib = fib + w3[k] * F
            row = row + w2[j] * fib
            nodes_done += len(x3)
        tot = tot + w1[i] * row
    c1 = resource.getrusage(resource.RUSAGE_SELF)
    cpu = c1.ru_utime + c1.ru_stime - c0.ru_utime - c0.ru_stime
    return tot, nodes_done, cpu, time.time() - t0

def chunk_cmd(chunk, out_path):
    _require_pats()
    rules = load_rules()
    ci, cn = (int(t) for t in chunk.split("/"))
    i0 = ci * N1 // cn; i1 = (ci + 1) * N1 // cn
    tot, nd, cpu, wall = sweep_range(rules, PATS, i0, i1, PREC_SWEEP)
    part = {"chunk": chunk, "i0": i0, "i1": i1, "nodes": nd,
            "sum": man_exp_pair(tot), "cpu_s": round(cpu, 1),
            "wall_s": round(wall, 1),
            "us_per_node": round(cpu * 1e6 / nd, 3),
            "recorded": time.strftime("%Y-%m-%d %H:%M:%S %Z")}
    json.dump(part, open(out_path, "w"))
    print(f"[sweep {chunk}] nodes={nd} cpu={cpu:.0f}s {part['us_per_node']}us/node",
          flush=True)
    return 0

# ---------------- N=8 exact-rational surrogate (inside the same harness) -----
class TP:
    """tiny trivariate polynomial, dict (i,j,k) -> Fraction."""
    def __init__(self, d=None):
        self.d = dict(d or {})
    @staticmethod
    def const(c):
        return TP({(0, 0, 0): Fraction(c)})
    def __add__(self, o):
        d = dict(self.d)
        for m, c in o.d.items():
            d[m] = d.get(m, Fraction(0)) + c
        return TP(d)
    def __mul__(self, o):
        d = {}
        for m1, c1 in self.d.items():
            for m2, c2 in o.d.items():
                m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])
                d[m] = d.get(m, Fraction(0)) + c1 * c2
        return TP(d)
    def scale(self, c):
        return TP({m: v * Fraction(c) for m, v in self.d.items()})
    def powi(self, n):
        r = TP.const(1)
        b = self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
        return r

def surrogate_cmd():
    """N=8 exact-rational cross-check + planted control, run through the SAME
    rule builder + sweep code at (17, 13, 5)."""
    t0 = time.time()
    pi1 = Fraction(P0_NUM, P0_DEN); pi0 = 1 - pi1
    pi = (pi0, pi1)
    U1 = TP({(1, 0, 0): Fraction(1)}); U2 = TP({(0, 1, 0): Fraction(1)})
    U3 = TP({(0, 0, 1): Fraction(1)})
    U12 = U1 * U2; U123 = U12 * U3
    def trans(um):
        # P[i][j] = pi_j + (delta_ij - pi_j) * um
        return [[TP.const(pi[j]) + um.scale(int(i == j) - pi[j]) for j in (0, 1)]
                for i in (0, 1)]
    P1 = trans(U1); P2 = trans(U2); P3 = trans(U3)
    P12 = trans(U12); P123 = trans(U123)
    def pat_poly(x):
        a, b, c, d = x
        tot = TP()
        for r in (0, 1):
            for w in (0, 1):
                for v in (0, 1):
                    tot = tot + (P3[r][w] * P2[w][v] * P1[v][a] * P1[v][b]
                                 * P12[w][c] * P123[r][d]).scale(pi[r])
        return tot
    counts_s = [((0, 0, 1, 1), 6), ((1, 1, 1, 0), 1), ((1, 0, 1, 0), 1)]  # N=8
    counts_plant = [((0, 0, 1, 1), 5), ((1, 1, 1, 0), 2), ((1, 0, 1, 0), 1)]
    # y3-parity exact assert on ALL 16 patterns
    parity_ok = True
    for i in range(16):
        pp = pat_poly(tuple(int(ch) for ch in f"{i:04b}"))
        for (a, b, k), cf in pp.d.items():
            if k % 2 == 1 and cf != 0:
                parity_ok = False
    assert parity_ok, "y3-parity FAILED in exact surrogate"
    F = TP.const(1)
    for x, n in counts_s:
        F = F * pat_poly(x).powi(n)
    degs = (max(m[0] for m in F.d), max(m[1] for m in F.d), max(m[2] for m in F.d))
    assert degs == (32, 24, 16), f"surrogate degrees {degs}"
    I_exact = Fraction(0)
    for (i, j, k), cf in F.d.items():
        if cf:
            I_exact += cf / ((C0 + i) * (C0 + j) * (C0 + k))
    Z_exact = C0 ** 3 * I_exact
    # certified harness: same builder + same sweep code, v-rule at alpha3, n=5
    rules_s = {"u1": build_rule(17, ALPHA1, "s_u1"),
               "u2": build_rule(13, ALPHA1, "s_u2"),
               "v3": build_rule(5, ALPHA3, "s_v3")}
    tot, nd, cpu, wall = sweep_range(rules_s, counts_s, 0, 17, PREC_SWEEP)
    ctx.prec = 256
    Zball = fr2arb(C0) ** 3 * tot / 2
    Zex = fr2arb(Z_exact)
    agree = bool(Zball.contains(Zex))
    reld = float(ub_of(abs(Zball - Zex) / Zex))
    tot_p, _, _, _ = sweep_range(rules_s, counts_plant, 0, 17, PREC_SWEEP)
    Zplant = fr2arb(C0) ** 3 * tot_p / 2
    fired = bool(abs(Zplant - Zex) > (Zex * arb("1e-6")))
    out = {
        "surrogate": "N=8 exact-rational cross-check (counts 0011:6,1110:1,1010:1)",
        "rules": "(17,13,5) via the SAME certified builder; sweep via the SAME code",
        "y3_parity_exact": parity_ok, "degrees": degs,
        "Z_exact_str": f"{Z_exact.numerator}/{Z_exact.denominator}",
        "ln_Z_exact": float(Zex.log().mid()),
        "harness_ball_contains_exact": agree,
        "rel_diff_ub": reld,
        "harness_relwidth_log2": float((Zball.rad() / Zball.mid()).log()
                                       / arb(2).log()),
        "planted_count_control_fires": fired,
        "cpu_s": round(cpu + time.time() - t0, 1),
    }
    json.dump(out, open(f"{SCRATCH}/rung1_surrogate.json", "w"), indent=1)
    print(json.dumps(out, indent=1), flush=True)
    if not (agree and fired):
        raise RuntimeError("surrogate cross-check FAILED")
    return 0

# ---------------- float shadow (diagnostic only, NEVER evidence) -------------
def shadow_cmd():
    _require_pats()
    import numpy as np
    from scipy.special import roots_jacobi
    t0 = time.time()
    p0 = P0_NUM / P0_DEN; q0 = 1 - p0
    al1 = float(ALPHA1); al3 = float(ALPHA3)
    def rule_f(n, al):
        x, w = roots_jacobi(n, 0.0, al)
        return (x + 1) / 2, w / 2 ** (al + 1)
    x1, w1 = rule_f(N1, al1); x2, w2 = rule_f(N2, al1); x3, w3 = rule_f(N3, al3)
    lnw2 = np.log(w2); lnw3 = np.log(w3)
    pi = (q0, p0)
    M = np.full(N1, -np.inf)
    nneg = 0
    for i in range(N1):
        u1 = x1[i]
        u2 = x2                    # vector
        u12 = u1 * u2
        P1 = ((q0 + p0 * u1, p0 * (1 - u1)), (q0 * (1 - u1), p0 + q0 * u1))
        P2 = ((q0 + p0 * u2, p0 * (1 - u2)), (q0 * (1 - u2), p0 + q0 * u2))
        P12 = ((q0 + p0 * u12, p0 * (1 - u12)), (q0 * (1 - u12), p0 + q0 * u12))
        gg = {}
        for w in (0, 1):
            for a in (0, 1):
                for b in (0, 1):
                    gg[(w, a, b)] = (P2[w][0] * P1[0][a] * P1[0][b]
                                     + P2[w][1] * P1[1][a] * P1[1][b])
        lnF = np.zeros((N2, N3))
        for (a, b, c, d), n in PATS:
            b0 = q0 * gg[(0, a, b)] * P12[0][c]
            b1 = p0 * gg[(1, a, b)] * P12[1][c]
            A = (b0 + b1) * pi[d]
            B = (b0 * ((1 if d == 0 else 0) - pi[d])
                 + b1 * ((1 if d == 1 else 0) - pi[d])) * u12
            P = A[:, None] + B[:, None] * x3[None, :]
            bad = P <= 0
            nneg += int(bad.sum())
            P = np.where(bad, 1e-300, P)
            lnF += n * np.log(P)
        z = lnw2[:, None] + lnw3[None, :] + lnF
        zm = z.max()
        M[i] = zm + np.log(np.exp(z - zm).sum())
    lw1M = np.log(w1) + M
    zm = lw1M.max()
    lnT = zm + np.log(np.exp(lw1M - zm).sum())
    lnZ = 3 * math.log(float(C0)) + math.log(0.5) + lnT
    out = {"float_shadow_lnZs": lnZ, "neg_P_nodes_masked": nneg,
           "wall_s": round(time.time() - t0, 1),
           "note": "float64 log-domain shadow of the same tensor rule -- "
                   "DIAGNOSTIC ONLY, never evidence"}
    json.dump(out, open(f"{SCRATCH}/rung1_shadow.json", "w"), indent=1)
    print(json.dumps(out), flush=True)
    return 0

# ---------------- endpoint-cluster node check --------------------------------
def endpoint_cmd():
    _require_pats()
    rules = load_rules()
    ctx.prec = PREC_SWEEP
    pi1 = arb(P0_NUM) / P0_DEN; pi0 = 1 - pi1
    combos = []
    for i in (0, N1 - 1):
        for j in (0, N2 - 1):
            for k in (0, N3 - 1):
                combos.append((i, j, k))
    rows = []
    worst = -1e9
    for (i, j, k) in combos:
        u1 = rules["u1"]["nodes"][i]; u2 = rules["u2"]["nodes"][j]
        v = rules["v3"]["nodes"][k]
        A, Bv = ab_fiber(u1, u2, pi0, pi1, [x for x, _ in PATS])
        F = None; pw = -1e9
        for t, (_, n) in enumerate(PATS):
            P = A[t] + Bv[t] * v
            if not (P > 0):
                raise RuntimeError(f"P ball not positive at corner {(i,j,k)}")
            rw = float((P.rad() / P.mid()).log() / arb(2).log())
            pw = max(pw, rw)
            Pn = P if n == 1 else P ** n
            F = Pn if F is None else F * Pn
        Fw = float((F.rad() / F.mid()).log() / arb(2).log()) if F > 0 else None
        worst = max(worst, Fw if Fw is not None else 1e9)
        rows.append({"corner_ijk": [i, j, k], "worst_P_relwidth_log2": round(pw, 1),
                     "F_relwidth_log2": (round(Fw, 1) if Fw is not None else None),
                     "F_positive": bool(F > 0),
                     "lnF_mid": float(F.log().mid()) if F > 0 else None})
    out = {"endpoint_cluster_check": rows, "worst_F_relwidth_log2": round(worst, 1),
           "pass_192bit": bool(worst < -100)}
    json.dump(out, open(f"{SCRATCH}/rung1_endpoint.json", "w"), indent=1)
    print(json.dumps(out, indent=1), flush=True)
    return 0

# ---------------- final assembly ----------------
def assemble_cmd(chunk_paths, out_path):
    ctx.prec = PREC_SWEEP
    parts = [json.load(open(q)) for q in chunk_paths]
    ivs = sorted((p["i0"], p["i1"]) for p in parts)
    assert ivs[0][0] == 0 and ivs[-1][1] == N1
    for a, b in zip(ivs[:-1], ivs[1:]):
        assert a[1] == b[0], "chunk coverage gap"
    tot = arb(0)
    for p in parts:
        tot = tot + from_man_exp(p["sum"])
    Zs = fr2arb(C0) ** 3 * tot / 2
    assert Zs > 0, "sweep ball not positive"
    lnZ = Zs.log()
    relw = Zs.rad() / Zs.mid()
    relw_log2 = float(relw.log() / arb(2).log())
    rules = json.load(open(f"{SCRATCH}/rung1_rules.json"))
    surro = json.load(open(f"{SCRATCH}/rung1_surrogate.json"))
    shadow = json.load(open(f"{SCRATCH}/rung1_shadow.json"))
    endp = json.load(open(f"{SCRATCH}/rung1_endpoint.json"))
    cpu_sweep = sum(p["cpu_s"] for p in parts)
    nodes = sum(p["nodes"] for p in parts)
    us_node = cpu_sweep * 1e6 / nodes
    ln_lo = lb_of(lnZ); ln_hi = ub_of(lnZ)
    accept_relw = relw_log2 <= -150
    kill_relw = relw_log2 > -100
    endp_ok = endp["pass_192bit"]
    mag_ok = abs(float(lnZ.mid()) - shadow["float_shadow_lnZs"]) < 0.5
    verdict = "PASS" if (accept_relw and endp_ok and surro["harness_ball_contains_exact"]
                         and mag_ok) else ("KILL" if (kill_relw or not endp_ok) else "MARGINAL")
    row = {
        "artifact": "RUNG 1 -- GAUSS-COLLAPSE certified production sweep "
                    "(endpoint-cluster node checks included)",
        "object": "Z_s(p0) = c0^3 * int_{[0,1]^3} u^(c0-1) F(u;p0) du, F = raw "
                  "positive product of 318 pattern probabilities (14 observed "
                  "patterns), pinned caterpillar (((A,B),C),D); exact Gauss-"
                  "Jacobi tensor rule (637 x 478 x 160), y3-parity reduced "
                  "v=u3^2 third dim; per-node VALUE evaluation A_x + B_x*v3, "
                  "no coefficient expansion, streamed positive sum",
        "p0_dyadic": f"{P0_NUM}/{P0_DEN}", "lambda": LAM,
        "c0": f"{C0.numerator}/{C0.denominator}",
        "prec_sweep": PREC_SWEEP,
        "rule_certification": {tag: rules[tag]["receipts"] for tag in rules},
        "sweep_ball": {
            "Zs_lo": ln_lo.exp().str(40), "Zs_hi": ln_hi.exp().str(40),
            "lnZs_enclosure": [float(ln_lo) - 1e-12, float(ln_hi) + 1e-12],
            "lnZs_str": lnZ.str(45),
            "relwidth_log2": relw_log2,
            "nodes_evaluated": nodes,
        },
        "acceptance": {
            "relwidth_le_2^-150": accept_relw,
            "kill_relwidth_gt_2^-100": kill_relw,
            "weight_balls_all_positive": all(rules[t]["receipts"]["all_weights_positive"]
                                             for t in rules),
            "endpoint_cluster_192bit_ok": endp_ok,
            "endpoint_worst_F_relwidth_log2": endp["worst_F_relwidth_log2"],
            "exact_rational_crosscheck_contains": surro["harness_ball_contains_exact"],
            "exact_rational_rel_diff_ub": surro["rel_diff_ub"],
            "planted_count_control_fires": surro["planted_count_control_fires"],
            "float_shadow_lnZs": shadow["float_shadow_lnZs"],
            "float_shadow_agree_0.5nat": mag_ok,
        },
        "verdict": verdict,
        "cost": {
            "sweep_cpu_s": round(cpu_sweep, 1),
            "sweep_cpu_h": round(cpu_sweep / 3600, 3),
            "us_per_node_measured": round(us_node, 3),
            "chunks": [{"chunk": p["chunk"], "cpu_s": p["cpu_s"],
                        "us_per_node": p["us_per_node"]} for p in parts],
            "projected_per_sweep_production_cost": {
                "this_object_48.7M_nodes_cpu_h": round(us_node * 48717760 / 3.6e9, 3),
                "note": "pure python-flint prec 192 on this box; proposal claimed "
                        "8.9 us/node (0.12 core-h) -- honest measured number here",
            },
        },
        "endpoint_cluster_rows": endp["endpoint_cluster_check"],
        "surrogate": surro,
        "recorded": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
    }
    json.dump(row, open(out_path, "w"), indent=1)
    print(json.dumps({k: row[k] for k in ("sweep_ball", "verdict")}, indent=1),
          flush=True)
    print(json.dumps(row["cost"], indent=1), flush=True)
    return 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["rules", "anchors", "surrogate", "shadow",
                                    "endpoint", "chunk", "assemble", "pilot"])
    ap.add_argument("--chunk", type=str, default="0/1")
    ap.add_argument("--chunks", type=str, nargs="*", default=None)
    ap.add_argument("--out", type=str, default=None)
    a = ap.parse_args()
    if a.cmd == "rules":
        sys.exit(rules_cmd())
    if a.cmd == "anchors":
        sys.exit(anchors_cmd())
    if a.cmd == "surrogate":
        sys.exit(surrogate_cmd())
    if a.cmd == "shadow":
        sys.exit(shadow_cmd())
    if a.cmd == "endpoint":
        sys.exit(endpoint_cmd())
    if a.cmd == "pilot":
        _require_pats()
        rules = load_rules()
        tot, nd, cpu, wall = sweep_range(rules, PATS, 300, 302, PREC_SWEEP)
        print(json.dumps({"nodes": nd, "cpu_s": round(cpu, 2),
                          "us_per_node": round(cpu * 1e6 / nd, 3),
                          "projected_full_sweep_min": round(cpu / nd * 48717760 / 60, 1)}),
              flush=True)
        sys.exit(0)
    if a.cmd == "chunk":
        out = a.out or f"{SCRATCH}/rung1_part_{a.chunk.replace('/', '_')}.json"
        sys.exit(chunk_cmd(a.chunk, out))
    sys.exit(assemble_cmd(a.chunks, a.out or f"{RESULTS}/RUNG1_GAUSS_SWEEP.json"))

if __name__ == "__main__":
    main()

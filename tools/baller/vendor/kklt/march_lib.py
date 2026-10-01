# baller transport engine — certified Taylor step kernel.
#!/usr/bin/env python3
"""march_lib.py — arb-ball Taylor endpoint transport for the exact 16x16
Pfaffian system  den(s) * s * v'(s) = Anum(s) v(s)  (ball law end-to-end).
Rigorous per-step tail: ||A|| <= CA on |t|<=r (arb interval eval) =>
||v|| <= N0 e^{CA r} on the disk, Cauchy => ||d_k|| <= Nr/r^k, geometric
tail (|h|/r)^{K+1}/(1-|h|/r).  No float64 reduction anywhere."""
from flint import acb, arb, acb_poly, ctx


def taylor_at(pol, z0, K):
    """first K+1 Taylor coeffs of pol at z0 via synthetic division (acb)."""
    out = []
    lin = acb_poly([-z0, acb(1)])
    q = pol
    for _ in range(K + 1):
        q, r = divmod(q, lin)
        out.append(r[0] if r.length() > 0 else acb(0))
        if q.length() == 0:
            out += [acb(0)] * (K - len(out) + 1)
            break
    return out


def mat_norm_ub(rows):
    """upper bound on max-row-sum norm (arb upper)."""
    best = arb(0)
    for row in rows:
        s = arb(0)
        for e in row:
            s += abs(e)
        best = arb(max(float(best.upper()), 0)).union(s) if False else \
            (s if s > best else (best if best > s else best.union(s)))
    return best.upper()


def vec_norm_ub(v):
    s = arb(0)
    for e in v:
        s += abs(e)
    return s.upper()


def _ubf(x):
    """float() rounds to NEAREST — an arf upper bound can
    round DOWN below the certified value; convert radius-position floats
    outward (nextafter toward +inf)."""
    import math
    return math.nextafter(float(x), math.inf)


def step(NUMP, DENP, RK, z0, v0, h, r, K):
    """one Taylor step z0 -> z0+h; returns new ball vector (list of acb).
    NUMP: dict (i,j)->acb_poly; DENP: acb_poly (den(s)*s).  Requires
    |h| < r and DENP nonzero on |t-z0| <= r (checked)."""
    # input contract ENFORCED — with |h| >= r the
    # geometric tail goes NEGATIVE and arb(0, negative) silently ABSOLUTIZES
    # to a too-small radius (wrong certified enclosure, no error); negative r
    # underestimated the Cauchy majorant the same silent way. Refuse unless
    # provably in contract (flint comparisons are True only when provable).
    if len(v0) != RK:
        raise ValueError(f"march step: v0 length {len(v0)} != RK {RK}")
    if not (arb(r) > 0):
        raise ValueError(f"march step: r must be provably > 0 (got {r!r})")
    if not (abs(acb(h)) < arb(r)):
        raise ValueError(
            f"march step: requires |h| provably < r (got |h|={abs(acb(h))}, r={r!r})")
    ball = z0 + acb(arb(0, _ubf(abs(arb(r)).upper())),
                    arb(0, _ubf(abs(arb(r)).upper())))
    dv = DENP(ball)
    if dv.contains(0):     # explicit check, not assert — must survive -O
        raise ValueError("march step: den*s vanishes inside step disk")
    # CA = sup ||A|| on the disk (row-sum norm of Anum(ball)/den(ball))
    rows = []
    for i in range(RK):
        row = []
        for j in range(RK):
            p = NUMP.get((i, j))
            row.append((p(ball) / dv) if p is not None else acb(0))
        rows.append(row)
    CA = mat_norm_ub([[abs(e) for e in row] for row in rows])
    g = taylor_at(DENP, z0, K)
    An = {}
    for (i, j), p in NUMP.items():
        An[(i, j)] = taylor_at(p, z0, K)
    d = [list(v0)]
    for k in range(K):
        rhs = [acb(0)] * RK
        for i in range(k + 1):           # An_i . d[k-i]
            dk = d[k - i]
            for a in range(RK):
                for b in range(RK):
                    e = An.get((a, b))
                    if e is not None:
                        rhs[a] += e[i] * dk[b]
        for i in range(1, k + 2):        # - g_i (k+1-i) d[k+1-i]
            if i <= k + 1 and (k + 1 - i) <= k:
                co = g[i] * (k + 1 - i)
                dk = d[k + 1 - i]
                for a in range(RK):
                    rhs[a] -= co * dk[a]
        inv = 1 / (g[0] * (k + 1))
        d.append([x * inv for x in rhs])
    # evaluate at h + rigorous tail ball
    N0 = vec_norm_ub(v0)
    Nr = arb(N0) * (arb(CA) * arb(r)).exp()
    q = abs(acb(h)) / arb(r)   # h may be acb (complex step); arb(acb)
                               # raises TypeError under python-flint 0.8.0. acb(h)
                               # accepts float/arb/acb; abs(acb)->arb.
    tail = (Nr * q ** (K + 1) / (1 - q)).upper()
    out = []
    tb = acb(arb(0, _ubf(tail)), arb(0, _ubf(tail)))   # outward float on the tail ball
    for a in range(RK):
        acc = acb(0)
        for k in range(K, -1, -1):
            acc = acc * h + d[k][a]
        out.append(acc + tb)
    return out

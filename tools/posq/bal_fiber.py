# bal_fiber.py -- balanced-shape per-fiber evaluator, mirroring
# bench_rung1_gauss.ab_fiber: P_y = A_y + B_y * v for all 16 wedge-slot
# patterns y = (x1,x2,x3,x4) (old1, old2, young1, young2).
#   e1 = w1*w2 (old leaves), e2 = w1 (young leaves)
#   S1[y12] = sum_a pi_a C1_a,  S2[y34] = sum_b pi_b C2_b,
#   T[y12,y34] = sum_a pi_a C1_a C2_a   (root-edge collapse, stationarity)
#   A_y = S1*S2,  B_y = w2 * (T - S1*S2)          [v = w3^2]
from flint import arb

def bal_fiber(w1, w2, pi0, pi1, pats):
    one = arb(1)
    e1 = w1 * w2
    e2 = w1
    m1 = one - e1
    m2 = one - e2
    # T1[a][x] = pi_x + (delta_ax - pi_x) * e1 ; T2[b][x] likewise with e2
    T1 = ((pi0 + pi1 * e1, pi1 * m1), (pi0 * m1, pi1 + pi0 * e1))
    T2 = ((pi0 + pi1 * e2, pi1 * m2), (pi0 * m2, pi1 + pi0 * e2))
    C1 = {}
    C2 = {}
    for xa in (0, 1):
        for xb in (0, 1):
            C1[(xa, xb)] = (T1[0][xa] * T1[0][xb], T1[1][xa] * T1[1][xb])
            C2[(xa, xb)] = (T2[0][xa] * T2[0][xb], T2[1][xa] * T2[1][xb])
    S1 = {k: pi0 * v[0] + pi1 * v[1] for k, v in C1.items()}
    S2 = {k: pi0 * v[0] + pi1 * v[1] for k, v in C2.items()}
    A = []
    Bv = []
    for (x1, x2, x3, x4) in pats:
        s1 = S1[(x1, x2)]
        s2 = S2[(x3, x4)]
        a = s1 * s2
        t = pi0 * C1[(x1, x2)][0] * C2[(x3, x4)][0] \
            + pi1 * C1[(x1, x2)][1] * C2[(x3, x4)][1]
        A.append(a)
        Bv.append(w2 * (t - a))
    return A, Bv

def sweep_range_bal(rules, expvecs, pnum, pden, i0, i1, prec):
    """multi-case balanced partial sum: for each exponent 16-vector in
    expvecs, sum_i w1 sum_j w2 sum_k w3 prod_y P_y^vec[y]. Returns list of
    balls (one per vector), nodes_done, cpu, wall."""
    import time, resource
    from flint import ctx
    ctx.prec = prec
    pi1 = arb(pnum) / pden
    pi0 = 1 - pi1
    x1, w1 = rules["u1"]["nodes"], rules["u1"]["weights"]
    x2, w2 = rules["u2"]["nodes"], rules["u2"]["weights"]
    x3, w3 = rules["v3"]["nodes"], rules["v3"]["weights"]
    allpats = [tuple(int(ch) for ch in f"{y:04b}") for y in range(16)]
    used = sorted({y for v in expvecs for y in range(16) if v[y]})
    upats = [allpats[y] for y in used]
    plans = [[(ui, vec[y]) for ui, y in enumerate(used) if vec[y]]
             for vec in expvecs]
    nvec = len(expvecs)
    tots = [arb(0) for _ in range(nvec)]
    nodes_done = 0
    t0 = time.time()
    c0 = resource.getrusage(resource.RUSAGE_SELF)
    for i in range(i0, i1):
        u1 = x1[i]
        rows = [arb(0) for _ in range(nvec)]
        for j in range(len(x2)):
            A, Bv = bal_fiber(u1, x2[j], pi0, pi1, upats)
            fibs = [arb(0) for _ in range(nvec)]
            for k in range(len(x3)):
                v = x3[k]
                V = [A[t] + Bv[t] * v for t in range(len(used))]
                for s in range(nvec):
                    F = None
                    for ui, n in plans[s]:
                        Pn = V[ui] if n == 1 else V[ui] ** n
                        F = Pn if F is None else F * Pn
                    fibs[s] += w3[k] * F
                nodes_done += 1
            for s in range(nvec):
                rows[s] += w2[j] * fibs[s]
        for s in range(nvec):
            tots[s] += w1[i] * rows[s]
    c1 = resource.getrusage(resource.RUSAGE_SELF)
    cpu = c1.ru_utime + c1.ru_stime - c0.ru_utime - c0.ru_stime
    return tots, nodes_done, cpu, time.time() - t0

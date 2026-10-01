# baller engine — validated ball-Taylor Fuchsian monodromy transport.
# Validated (ball-arithmetic) monodromy transport of an order-5 operator L5 at
# rates (2,3,5). Rigorous Taylor stepping of the ODE sum_k b_k(x) f^(k)(x) = 0
# with EXACT integer coefficient polynomials b_k (from the pinned BC_cache.pkl
# beside this file), from the dyadic basepoint x_b = float(0.01).
#
# KNOWN GAP (documented, not hidden): the loop GEOMETRY (loop_paths) is NOT
# vendored — loop_waypoints() raises a typed error; pass explicit waypoint
# lists to transport_loop(..., waypoints=...) instead. __main__ runs a
# synthetic self-check (null-homotopic loop => identity monodromy).
#
# Rigor chain per step (center z0 dyadic double, step h = z1-z0 exact in acb):
#   * shifted coefficients beta_k[j] of b_k(z0+t): exact dyadic (integer polys,
#     double centers, prec 380 >> needed bits).
#   * disc radius rho chosen s.t. LB := |b5(z0)| - sum_{j>=1}|beta5_j| rho^j > 0
#     (series lower bound; rigorous |b5| >= LB on D(z0,rho)); also rho <= 0.7 *
#     (distance to nearest singular-point ball).
#   * Taylor coefficients c_0..c_N of all 5 fundamental columns from the EXACT
#     recurrence (ball arithmetic; initial data = ball jets; coefficients carry
#     no truncation error).
#   * Tail bound: on D(z0,rho), ||A||_inf <= M with A the companion matrix,
#     row bound (sum_m UB|b_m|)/LB|b5|, UB|b_m| = sum_j |beta_mj| rho^j. For
#     column j, sup_D ||Y_j||_inf <= u_j exp(M rho) (Gronwall along radii; disc
#     convex), so by Cauchy the d-th component's coefficients obey
#     |c^[d]_n| <= G_j / rho^n, and the truncation tail at |h| is
#     <= G_j q^(N-d+1)/(1-q), q = |h|/rho <= 0.51. Added as ball radius.
# The final 5x5 ball matrix rigorously encloses the true monodromy matrix in the
# jet frame (f, f', f'', f''', f'''') at the dyadic basepoint.
import sys, os, json, time, math
from flint import acb, arb, acb_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
# Output root for saved loop matrices: default the CURRENT working directory —
# never this module's directory (unpinned files under vendor/ trip
# baller.verify()).
OUTDIR = os.environ.get("VALMONO_OUT", os.getcwd())

def get_BC():
    import pickle
    cache = os.path.join(HERE, 'BC_cache.pkl')
    if not os.path.exists(cache):
        raise FileNotFoundError(
            'BC_cache.pkl (the pinned exact ODE coefficient cache, provenance '
            'in PROVENANCE_BC_cache.txt) is missing beside valmono.py — '
            'refusing; the derivation script that regenerates it is not '
            'vendored')
    return pickle.load(open(cache, 'rb'))

BC = get_BC()
DEG = [len(c) - 1 for c in BC]
MAXDEG = max(DEG)

def ff(n, k):
    r = 1
    for i in range(k):
        r *= (n - i)
    return r

def singular_balls(prec):
    ctx.prec = prec
    s2, s3, s5 = arb(2).sqrt(), arb(3).sqrt(), arb(5).sqrt()
    pts = [acb(0)]
    for e1 in (1, -1):
        for e2 in (1, -1):
            w = s2 + e1*s3 + e2*s5
            pts.append(acb(1)/(w*w))
    disc = (acb(400) - acb(4*216*3)).sqrt()   # apparent pair, 216x^2+20x+3
    pts.append((acb(-20) + disc)/acb(432))
    pts.append((acb(-20) - disc)/acb(432))
    return pts

def shifted_beta(z0c, prec):
    ctx.prec = prec
    z0 = acb(z0c.real, z0c.imag)
    zp = [acb(1)]
    for i in range(MAXDEG):
        zp.append(zp[-1]*z0)
    beta = []
    for k in range(6):
        d = DEG[k]
        bk = []
        for j in range(d + 1):
            s = acb(0)
            for i in range(j, d + 1):
                s += acb(BC[k][i]*math.comb(i, j))*zp[i - j]
            bk.append(s)
        beta.append(bk)
    return beta

def choose_rho(beta, dclear, prec):
    """Largest rho = 0.7*dclear/2^t (t<=14) with series LB(|b5|) >= 0.4*|b5(z0)|."""
    ctx.prec = prec
    b50a = beta[5][0].abs_lower()
    if not (b50a > arb(0)):
        return None, None
    rho = 0.7*dclear
    for _ in range(15):
        r = arb(rho)
        sub = arb(0)
        rp = arb(1)
        for j in range(1, DEG[5] + 1):
            rp = rp*r
            sub += beta[5][j].abs_upper()*rp
        lb = arb(b50a) - sub
        if lb.abs_lower() > arb(0) and lb > arb(b50a)*arb(0.4):
            return rho, lb
        rho *= 0.5
    return None, None

def step(beta, z0c, z1c, Ymat, rho, prec, nmax=6000):
    """One validated step z0c -> z1c (python complex doubles; h formed EXACTLY
    in acb — double subtraction would round and detach the evaluated point from
    the tracked path). beta from shifted_beta(z0c). Returns (Ynew, err)."""
    ctx.prec = prec
    h = acb(z1c.real, z1c.imag) - acb(z0c.real, z0c.imag)   # exact dyadic
    r = arb(rho)
    # |b_m| upper bounds and |b5| lower bound on the disc, from series coeffs
    ub = []
    for k in range(6):
        s = arb(0)
        rp = arb(1)
        for j in range(DEG[k] + 1):
            s += beta[k][j].abs_upper()*rp
            rp = rp*r
        ub.append(s)
    sub = arb(0)
    rp = arb(1)
    for j in range(1, DEG[5] + 1):
        rp = rp*r
        sub += beta[5][j].abs_upper()*rp
    b5lo = arb(beta[5][0].abs_lower()) - sub
    if not (b5lo.abs_lower() > arb(0)):
        return None, 'b5 lower bound not positive'
    Mb = arb(1)
    for m in range(5):
        Mb += ub[m]/b5lo
    q_up = arb(abs(h).abs_upper())/r
    qf = min(max(float(q_up.abs_upper().str(10, radius=False).replace('[','').replace(']','') if False else 0.51), 1e-9), 0.75)
    qf = 0.52  # h chosen <= 0.51*rho by caller; fixed budget
    N = int(prec*1.05/(-math.log2(qf))) + 40
    if N > nmax:
        return None, f'N={N} exceeds nmax'
    rows = [None]*(N + 6)
    fact = 1
    for d in range(5):
        if d:
            fact *= d
        rr_ = acb_mat(1, 5)
        for j in range(5):
            rr_[0, j] = Ymat[d, j]/acb(fact)
        rows[d] = rr_
    terms = []
    for k in range(6):
        for j in range(DEG[k] + 1):
            if not (k == 5 and j == 0):
                terms.append((k, j, beta[k][j]))
    b50 = beta[5][0]
    zero = acb_mat(1, 5)
    for m in range(N + 1):
        acc = acb_mat(1, 5)
        for (k, j, bkj) in terms:
            n2 = m - j + k
            if 0 <= n2 <= m + 4:
                f2 = ff(n2, k)
                if f2:
                    acc += (bkj*acb(f2))*rows[n2]
        rows[m + 5] = acc*(acb(-1)/(b50*acb(ff(m + 5, 5))))
    out = acb_mat(5, 5)
    for d in range(5):
        S = acb_mat(1, 5)
        for n in range(N + 5, d - 1, -1):
            S = S*h + acb(ff(n, d))*rows[n]
        for j in range(5):
            out[d, j] = S[0, j]
    # Fuchsian-scaled Gronwall: Yhat_k = rho^k f^(k); ||Ahat||_inf * rho <=
    # max(1, sum_m ub_m rho^(5-m)/b5lo) =: E (O(1) along the path).
    E = arb(1)
    S5 = arb(0)
    for m in range(5):
        S5 += ub[m]*r**(5 - m)/b5lo
    if S5 > E:
        E = S5
    expf = E.exp()
    rpow = [arb(1)]
    for _ in range(5):
        rpow.append(rpow[-1]*r)
    for j in range(5):
        uj = arb(0)
        for d in range(5):
            a = Ymat[d, j].abs_upper()*rpow[d]
            if a > uj:
                uj = a
        Gj = arb(uj)*expf
        for d in range(5):
            # |c^[d]_n| <= Gj / (rho^d rho^n)  =>  tail at |h| over n > N-d
            tail = (Gj/rpow[d])*q_up**(N + 6 - d)/(arb(1) - q_up)
            tu = tail.abs_upper()
            out[d, j] += acb(arb(0, tu), arb(0, tu))
    return out, None

def loop_waypoints():
    try:
        from monodromy_transport2 import loop_paths
    except ImportError as e:
        raise ImportError(
            'the loop-geometry module (loop_paths) is not vendored: pass an '
            'explicit waypoint list to transport_loop(..., waypoints=...) '
            'instead, or put a checkout providing monodromy_transport2 on '
            'sys.path') from e
    P = loop_paths(0.01)
    xb = 0.01
    return {li: [xb] + [complex(p) for p in P[li]] + [xb] for li in P}

def transport_loop(li, prec, log=print, waypoints=None):
    ctx.prec = prec
    W = waypoints if waypoints is not None else loop_waypoints()[li]
    sings = singular_balls(prec)
    sf = [complex(s.mid()) for s in sings]
    srad = [max(float(abs(complex(s.rad()))), 1e-30) for s in sings]
    Y = acb_mat(5, 5)
    for i in range(5):
        Y[i, i] = acb(1)
    t0 = time.time()
    nsteps = 0
    z = W[0]
    for w in W[1:]:
        if w == z:
            continue
        while z != w:
            d = min(abs(z - s) - rr for s, rr in zip(sf, srad))
            if d <= 0:
                raise RuntimeError(f'loop {li}: point {z} inside singular clearance')
            beta = shifted_beta(z, prec)
            rho, _lb = choose_rho(beta, d, prec)
            if rho is None:
                raise RuntimeError(f'loop {li}: no valid rho at {z}')
            rem = abs(w - z)
            if rem <= 0.5*rho:
                zn = w
            else:
                dirn = (w - z)/rem
                zn = z + dirn*(0.5*rho)
                zn = complex(float(zn.real), float(zn.imag))  # snap to double
                if abs(zn - z) > 0.51*rho:                    # paranoia after snap
                    zn = z + dirn*(0.45*rho)
                    zn = complex(float(zn.real), float(zn.imag))
            Ynew, err = step(beta, z, zn, Y, rho, prec)
            if Ynew is None:
                raise RuntimeError(f'loop {li} at {z}: {err}')
            Y = Ynew
            nsteps += 1
            z = zn
            if nsteps % 25 == 0:
                mr = max(float(abs(complex(Y[i, j].rad()))) for i in range(5) for j in range(5))
                log(f'  loop {li} step {nsteps} z={z:.4g} maxrad={mr:.3e} t={time.time()-t0:.1f}s')
    return Y, nsteps, time.time() - t0

def save_loop(li, Y, nsteps, dt, prec):
    maxrad = max(float(abs(complex(Y[i, j].rad()))) for i in range(5) for j in range(5))
    rows = [[[str(Y[i, j].real.mid()), str(Y[i, j].imag.mid()),
              str(Y[i, j].real.rad()), str(Y[i, j].imag.rad())] for j in range(5)] for i in range(5)]
    json.dump({'loop': li, 'prec': prec, 'steps': nsteps, 'secs': dt, 'maxrad': maxrad,
               'entries_mid_rad': rows},
              open(os.path.join(OUTDIR, f'valmono_L{li}_p{prec}.json'), 'w'))
    return maxrad

def selftest(prec=160):
    """Synthetic self-check: transport around a small null-homotopic loop
    (a diamond of radius 0.004 around the basepoint 0.01, enclosing no
    singular point), through the full shifted_beta/choose_rho/step machinery.
    The monodromy of a contractible loop is the identity: every entry ball
    must CONTAIN delta_ij with a small radius."""
    xb = 0.01
    r = 0.004
    W = [xb, xb + r, xb + 1j * r, xb - r, xb - 1j * r, xb + r, xb]
    Y, nsteps, dt = transport_loop(-1, prec, waypoints=W)
    maxrad = 0.0
    for i in range(5):
        for j in range(5):
            want = acb(1) if i == j else acb(0)
            d = Y[i, j] - want
            if not d.contains(acb(0)):
                raise AssertionError(
                    f'selftest FAIL: entry ({i},{j}) ball does not contain '
                    f'the identity value: {Y[i, j]}')
            maxrad = max(maxrad,
                         float(abs(complex(Y[i, j].rad()))))
    if maxrad > 1e-10:
        raise AssertionError(f'selftest FAIL: max entry radius {maxrad:.3e} '
                             f'> 1e-10 (enclosure blew up)')
    print(f'valmono selftest PASS: null loop -> identity, {nsteps} steps, '
          f'{dt:.1f}s, prec {prec}, max entry radius {maxrad:.3e}', flush=True)
    return maxrad


if __name__ == '__main__':
    prec = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    sys.exit(0 if selftest(prec) is not None else 1)

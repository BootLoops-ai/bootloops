#!/usr/bin/env python3
"""poly_de.py — exact polynomial DE  D(x)*dM/dx = G(x)*M  from a Kira table.

Takes the kira_targets .m path directly (no tag/points layout).  D = scalar
LCD of all Kira coef denominators at rational d; G fitted by DFT at scaled
roots of unity (cond=1), cross-checked at a fresh complex point.
"""
import os, re, time, pickle, hashlib
import mpmath as mp
mp.mp.dps = max(mp.mp.dps, 50)      # set FIRST
import gmpy2

from .de_core import assemble_A, evalA, g2mp

_A_CACHE = {}      # kira_path -> (A, N, miss)
_DEN_CACHE = {}    # kira_path -> sorted denom strings

# optional accelerators (both OFF unless env set; numerics unchanged):
#   FROBENIUS_POLYDE_CACHE=<dir>  — disk-cache the fitted poly-DE per
#       (kira file, d, prec, var, mass_slots, family); values round-trip via
#       prec_fit+5-digit strings (way below the ~1e-77 fit error).
#   FROBENIUS_SAMPLE_NCPU=<n>     — os.fork-parallel D*A sampling (safe inside
#       daemonic pool workers, unlike multiprocessing.Pool).


def _cache_file(kira_path, dd_rat, prec, var, mass_slots, family):
    cdir = os.environ.get('FROBENIUS_POLYDE_CACHE')
    if not cdir:
        return None
    os.makedirs(cdir, exist_ok=True)
    st = os.stat(kira_path)
    key = repr(('v1', os.path.realpath(kira_path), st.st_mtime, st.st_size,
                str(dd_rat), prec, var, tuple(mass_slots), family))
    return os.path.join(cdir, 'polyde_' +
                        hashlib.sha1(key.encode()).hexdigest()[:20] + '.pkl')


def _forkmap(fn, items, ncpu):
    """os.fork map (children inherit the parsed A; results as nstr strings)."""
    ncpu = min(ncpu, len(items))
    chunks = [list(range(i, len(items), ncpu)) for i in range(ncpu)]
    pipes, pids = [], []
    for ch in chunks:
        r, w = os.pipe()
        pid = os.fork()
        if pid == 0:
            os.close(r)
            try:
                out = {i: fn(items[i]) for i in ch}
                buf = pickle.dumps(out, protocol=4)
                with os.fdopen(w, 'wb') as f:
                    f.write(buf)
                os._exit(0)
            except BaseException:
                os._exit(1)
        os.close(w)
        pipes.append(r); pids.append(pid)
    res = {}
    for r, pid in zip(pipes, pids):
        with os.fdopen(r, 'rb') as f:
            buf = f.read()
        _, status = os.waitpid(pid, 0)
        if status != 0 or not buf:
            raise RuntimeError('poly_de._forkmap child failed')
        res.update(pickle.loads(buf))
    return [res[i] for i in range(len(items))]


def get_A(kira_path, basis, mass_slots, family=None):
    key = (kira_path, tuple(mass_slots))
    if key not in _A_CACHE:
        _A_CACHE[key] = assemble_A(kira_path, basis, mass_slots, family)
    return _A_CACHE[key]


def get_denom_strings(kira_path):
    if kira_path not in _DEN_CACHE:
        txt = open(kira_path).read()
        dens = set()
        for m in re.finditer(r'\)/\((.*?)\)\)\n', txt):
            dens.add(m.group(1))
        _DEN_CACHE[kira_path] = sorted(dens, key=len)
    return _DEN_CACHE[kira_path]


def get_denom_LCD_at_d(kira_path, dd_rat, var='m2'):
    """LCD of all Kira coef denoms as Poly(var) over QQ at d=dd_rat (gcd-based
    LCM; avoids factoring large irreducibles). Returns (coeffs ascending, deg).
    Fast path: Fraction-ring eval of each string (fast_lcd.parse_poly) instead
    of sympify+expand — ~600s CPU/eps -> seconds on the reference family (1690 strings,
    10-23KB Horner form).  Per-string sympy fallback on any parse failure."""
    from .fast_lcd import lcd_at_d
    return lcd_at_d(get_denom_strings(kira_path), dd_rat, var)


def build_poly_DE_full(kira_path, basis, dd_rat, prec, mass_slots,
                       var='m2', family=None, verbose=False):
    """Returns dict(Dc_rat, Dc, Gc, roots, degD, degG, fit_err_log10)."""
    import sympy as sp
    xs = sp.symbols(var)
    cf = _cache_file(kira_path, dd_rat, prec, var, mass_slots, family)
    if cf and os.path.exists(cf):
        d = pickle.load(open(cf, 'rb'))
        mp.mp.dps = prec + 60
        out = {'Dc_rat': [sp.Rational(c) for c in d['Dc_rat']],
               'Dc': [mp.mpf(c) for c in d['Dc']],
               'Gc': [mp.matrix([[mp.mpc(mp.mpf(re), mp.mpf(im))
                                  for re, im in row] for row in g])
                      for g in d['Gc']],
               'roots': [(mp.mpc(mp.mpf(re), mp.mpf(im)), 1)
                         for re, im in d['roots']],
               'N': d['N'], 'degD': d['degD'], 'degG': d['degG'],
               'fit_err_log10': d['fit_err_log10']}
        mp.mp.dps = prec
        if verbose:
            print(f'    [poly_de] cache hit {os.path.basename(cf)}', flush=True)
        return out
    mp.mp.dps = prec + 10
    Dc_rat, degD = get_denom_LCD_at_d(kira_path, dd_rat, var)
    dd = mp.mpf(int(dd_rat.p)) / int(dd_rat.q)
    # roots of square-free part (sympy nroots: huge-rational LCDs overflow np)
    Dpoly = sp.Poly(dict(enumerate(Dc_rat)), xs, domain='QQ')
    Dsf = sp.Poly(sp.cancel(Dpoly.as_expr() /
                            sp.gcd(Dpoly, Dpoly.diff()).as_expr()), xs, domain='QQ')
    rts = sp.nroots(Dsf, n=30, maxsteps=200)
    roots = [(mp.mpc(mp.mpf(str(sp.re(r))), mp.mpf(str(sp.im(r)))), 1) for r in rts]
    A, N, miss = get_A(kira_path, basis, mass_slots, family)
    if miss:
        raise RuntimeError(f'{kira_path}: {len(miss)} unreduced targets')
    # degG = degD + max positive var-power of A (2-point slope probe)
    rmax = max(abs(r) for r, _ in roots)
    mp.mp.dps = 50
    xa, xb = mp.mpf(10) * rmax, mp.mpf(100) * rmax
    Aa = evalA(A, N, xa, dd, var); Ab = evalA(A, N, xb, dd, var)
    qmax = 0
    for k in range(N):
        for l in range(N):
            a1, a2 = abs(Aa[k, l]), abs(Ab[k, l])
            if a1 > 1e-30 and a2 > 1e-30:
                q = int(round(float(mp.log10(a2 / a1))))
                qmax = max(qmax, q)
    degG = degD + qmax
    nS = degG + 1
    # sample D*A at |x|=R0 roots of unity (w=x/R0 makes the fit a pure DFT)
    R0 = mp.mpf('1.6') * rmax
    prec_fit = prec + int(degG * float(mp.log10(R0))) + 20
    mp.mp.dps = prec_fit
    Dc = [mp.mpf(int(c.p)) / int(c.q) for c in Dc_rat]
    pts = [R0 * mp.expjpi(2 * mp.mpf(k) / nS) for k in range(nS)]
    def evalD(x):
        v = mp.mpc(0)
        for c in reversed(Dc): v = v * x + c
        return v
    t0 = time.time()
    ncpu = int(os.environ.get('FROBENIUS_SAMPLE_NCPU', '1') or '1')
    if ncpu > 1:
        def _one(x):
            M = evalD(x) * evalA(A, N, x, dd, var)
            return [[(mp.nstr(M[i, j].real, prec_fit + 5),
                      mp.nstr(mp.im(M[i, j]), prec_fit + 5))
                     for j in range(N)] for i in range(N)]
        raws = _forkmap(_one, pts, ncpu)
        samples = [mp.matrix([[mp.mpc(mp.mpf(re), mp.mpf(im))
                               for re, im in row] for row in raw])
                   for raw in raws]
    else:
        samples = [evalD(x) * evalA(A, N, x, dd, var) for x in pts]
    if verbose:
        print(f'    [poly_de] sampled D*A at {nS} pts (|{var}|={float(R0):.2f}, '
              f'prec_fit={prec_fit}) in {time.time()-t0:.1f}s', flush=True)
    # inverse DFT entrywise in gmpy2; G_c = Ghat_c / R0^c
    gmpy2.get_context().precision = int(prec_fit * 3.33) + 30
    om = gmpy2.mpc(gmpy2.mpfr(mp.nstr(mp.cos(-2 * mp.pi / nS), prec_fit + 5)),
                   gmpy2.mpfr(mp.nstr(mp.sin(-2 * mp.pi / nS), prec_fit + 5)))
    R0g = gmpy2.mpfr(mp.nstr(R0, prec_fit + 5))
    sf = [[gmpy2.mpc(gmpy2.mpfr(mp.nstr(samples[r][i, j].real, prec_fit + 5)),
                     gmpy2.mpfr(mp.nstr(samples[r][i, j].imag, prec_fit + 5)))
           for i in range(N) for j in range(N)] for r in range(nS)]
    Gc = []
    R0p = gmpy2.mpfr(1)
    ocs = [om ** c for c in range(nS)]
    for c in range(nS):
        S = [gmpy2.mpc(0)] * (N * N); wc = gmpy2.mpc(1); oc = ocs[c]
        for r in range(nS):
            sr = sf[r]
            for ij in range(N * N): S[ij] += sr[ij] * wc
            wc *= oc
        sc = gmpy2.mpfr(1) / (nS * R0p)
        M = mp.matrix(N, N)
        for i in range(N):
            for j in range(N): M[i, j] = g2mp(S[i * N + j] * sc)
        Gc.append(M)
        R0p *= R0g
    mp.mp.dps = prec
    Dc = [mp.mpf(int(cR.p)) / int(cR.q) for cR in Dc_rat]
    g0n = max(abs(Gc[0][i, j]) for i in range(N) for j in range(N))
    while len(Gc) > 1 and max(abs(Gc[-1][i, j]) for i in range(N)
                              for j in range(N)) < g0n * mp.mpf(10) ** (-(prec - 10)):
        Gc.pop()
    degG = len(Gc) - 1
    # cross-check at a fresh complex point
    xt = mp.mpc(mp.mpf('1.7'), mp.mpf('2.31'))
    Gv = mp.matrix(N, N)
    for c in range(degG, -1, -1): Gv = Gv * xt + Gc[c]
    Dv = mp.mpc(0)
    for c in reversed(Dc): Dv = Dv * xt + c
    DA = Dv * evalA(A, N, xt, dd, var)
    err = max(abs(Gv[i, j] - DA[i, j]) for i in range(N) for j in range(N))
    nrm = max(abs(DA[i, j]) for i in range(N) for j in range(N))
    if verbose:
        print(f'    [poly_de] degD={degD} degG={degG} '
              f'fit-relerr={mp.nstr(err/nrm, 3)}', flush=True)
    out = {'Dc_rat': Dc_rat, 'Dc': Dc, 'Gc': Gc, 'roots': roots, 'N': N,
           'degD': degD, 'degG': degG,
           'fit_err_log10': float(mp.log10(err / nrm + mp.mpf(10) ** (-prec)))}
    if cf:
        nd = prec_fit + 5
        pickle.dump(
            {'Dc_rat': [str(c) for c in Dc_rat],
             'Dc': [mp.nstr(c, nd) for c in Dc],
             'Gc': [[[(mp.nstr(g[i, j].real, nd), mp.nstr(mp.im(g[i, j]), nd))
                      for j in range(N)] for i in range(N)] for g in Gc],
             'roots': [(mp.nstr(r.real, 40), mp.nstr(mp.im(r), 40))
                       for r, _ in roots],
             'N': N, 'degD': degD, 'degG': degG,
             'fit_err_log10': out['fit_err_log10']},
            open(cf, 'wb'), protocol=4)
    return out

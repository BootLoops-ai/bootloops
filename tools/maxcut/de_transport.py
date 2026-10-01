#!/usr/bin/env python3
r"""
maxcut_de_transport.py — 24×24 maximal-cut DE transport for a 3-loop
on-shell 4-point family (top sector, 24 top-sector masters of 41).

Pipeline (steps 1–6):
  1. assemble  A_s(s,t,d), A_t(s,t,d) from a source recipe + Kira reduction
  2. integrability  ∂_tA_s − ∂_sA_t + [A_s,A_t]  at N points, dps=DPS
  3. boundary at (s₀,t₀) via the corner-limit LS + route (a)/(b)
  4. transport per-ε (Bulirsch–Stoer stepper)
  5. ε-FFT extract ε⁰
  6. gate against an independently derived basis

STATUS: steps 1–2 exercised on production-shaped inputs; the steps 3–6
driver (boundary → per-ε transport → ε-FFT → gate) runs automatically
once step 2 PASSES, and its machinery is gated by the built-in synthetic
battery (`--selftest`: a commuting 24×24 toy pair with a closed-form
solution, gated on Laurent orders 0–3 plus a perturbed-basis negative
control). The default input files are NOT shipped in this repo; without
your own inputs the module fails LOUDLY at import/run. With a source
recipe whose s/t operators suffer the ∂_u contamination described below,
step 2 fails and the run exits 2 without executing steps 3–6.

INPUT FORMATS for steps 3–6 (JSON; numbers are strings — decimal or
`a/b` rational — and complex values are [re, im] string pairs):
  corner spec (MAXCUT_DE_CORNER / --corner), the boundary basepoint:
    {"point": {"s": "-1", "t": "-1"},
     "route": "a",                        per-master ε-series boundary
     "n_neg": 0,                          series starts at ε^{-n_neg}
     "series": [[c0, c1, ...] × 24]}      y0_i(ε) = Σ_k c_k·ε^{k-n_neg}
  or route (b), the corner-limit Γ-product times per-master weights:
    {"point": ..., "route": "b",
     "weights": [w × 24]}                 y0(ε) = mercedes_corner(ε)·w
  basis spec (MAXCUT_DE_BASIS / --basis), the independent gate values:
    {"point": {"s": ..., "t": ...},       transport target (s₁,t₁)
     "eps0": [v × 24],                    independently derived ε⁰ values
     "laurent": {"<order>": [v × 24]},    optional further gated orders
     "digits": 20}                        optional gate-threshold override
Transport path: s₀→s₁ at t=t₀, then t₀→t₁ at s=s₁ — the caller keeps it
clear of singular loci. Exit codes: 0 gate PASS · 2 integrability FAIL ·
3 gate FAIL.

OPERATOR CONTAMINATION (the failure class step 2 detects)
---------------------------------------------------------
A naive on-shell operator pair
   ∂_s = (1/s)·p1.∂_{p1},   ∂_t = (1/t)·p3.∂_{p3}
preserves p1²=p2²=p3²=0 **but NOT** p4²=(p1+p2+p3)²=0 (the 4th on-shell
condition):
   (1/s)p1.∂_{p1}[p4²] = (s+u)/s = −t/s,
   (1/t)p3.∂_{p3}[p4²] = (t+u)/t = −s/t.
Because the integrand DOES depend on p1.p3=u/2 *after loop integration*
(loop tensor reduction couples k.p1·k.p3 → p1.p3 etc., even though no
single propagator carries both p1 and p3), the contamination is nonzero:
   A_s^{naive} = ∂_s I − (t/s)·C,    A_t^{naive} = ∂_t I − (s/t)·C,
   C := ∂_{p4²}I |_{p4²=0}  (= "∂_u F" in the off-shell extension).
The combination s·A_s − t·A_t is contamination-free (both contaminations
cancel there, verified exact-rational) — a null check, so it HIDES the bug.

CORRECT operator pair (preserves all four pᵢ²=0; derived from
O_k := p_k.∂_{p_k}, k=1,2,3, all of which preserve every pᵢ²):
   ∂_s|_{on-shell} = (O₁+O₂−O₃)/(2s) − (O₁−O₂+O₃)/(2u),   u=−(s+t)
   ∂_t|_{on-shell} = (O₂+O₃−O₁)/(2t) − (O₁−O₂+O₃)/(2u).
O₁,O₃ recipes are s·A_s^{naive}, t·A_t^{naive} (already reduced); the O₂ =
p2.∂_{p2} recipe generates additional targets that must be present in the
Kira reduction file — a recipe target missing from the reduction raises
loudly in load_de24.

Diagnostic fact (measured on the reference inputs at dps=150, _RAT-safe
eval): under the contamination, the integrability residue ≈ |[A_s,A_t]| at
every point (curl ≪ commutator) — the contamination is rank-full, not a
single-row dropout. The maxcut reduction itself can still be correct (it
was, to 67–71 digits against independent 1-var-slice full-IBP reductions);
step 2 isolates the OPERATOR error from the reduction.
"""
import json, sys, os, re, time, argparse
import mpmath as mp
import sympy as sp

# ---------------------------------------------------------------------------
# The default input files below are not shipped in this repo: without your
# own files a production run fails LOUDLY at import/run — point the
# MAXCUT_DE_* env vars at your own reduction/recipe/corner/basis artifacts
# (formats in the module docstring). `--selftest` needs no inputs.
DEFAULTS = dict(
    src      = os.environ.get('MAXCUT_DE_SRC',      'source_recipe.json'),
    kira_red = os.environ.get('MAXCUT_DE_KIRA_RED', 'kira_sec_targets.m'),
    masters  = os.environ.get('MAXCUT_DE_MASTERS',  'masters'),
    corner   = os.environ.get('MAXCUT_DE_CORNER',   'TOPCUT_RESULT.json'),
    basis    = os.environ.get('MAXCUT_DE_BASIS',    'BASIS.json'),
    out      = os.environ.get('MAXCUT_DE_OUT',      'MAXCUT_DE_TRANSPORT.json'),
)
_FAMILY = os.environ.get('MAXCUT_DE_FAMILY')   # Kira family name in kira_red
# kira2math parser: the shipped frobenius-boundary copy, or your own
# kira_parse.py via MAXCUT_DE_PARSER_DIR (takes precedence when set).
_ktp = os.environ.get('MAXCUT_DE_PARSER_DIR')
_FB_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _cand in ([_ktp] if _ktp else []) + [
        os.path.join(_FB_ROOT, 'frobenius_boundary'),
        os.path.join(_FB_ROOT, 'frobenius-boundary', 'frobenius_boundary')]:
    if _cand and os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, _cand)
        break
from kira_parse import parse_kira2math   # noqa: E402  (loud when absent)

# ---- 1. assemble ----------------------------------------------------------
# PRECISION-CRITICAL: Kira coef strings carry ~3000 bare `int/int` literals
# per entry (Horner depth ~13×10).  Python eval() turns each into a 53-bit
# float BEFORE multiplying by mp.mpf vars → catastrophic precision loss
# (this trap once contaminated a "STILL FAILS at 50dps" verdict: the failure
# was REAL but the per-entry values were only ~12-d clean).
# Fix: regex-wrap every int that immediately precedes '/' so a/b → mpf(a)/b.
_RAT = re.compile(r'(?<![\w.])(\d+)(?=\s*/)')
_CCACHE = {}
def _fc(cs, dd, ss, tt):
    if cs == '1':
        return mp.mpf(1)
    co = _CCACHE.get(cs)
    if co is None:
        s2 = cs.replace('^', '**')
        s2 = _RAT.sub(r'_m(\1)', s2)
        co = compile(s2, '<coef>', 'eval')
        _CCACHE[cs] = co
    return eval(co, {'__builtins__': {}, '_m': mp.mpf}, {'d': dd, 's': ss, 't': tt})


def load_de24(src_path=DEFAULTS['src'], kira_path=DEFAULTS['kira_red']):
    """Return (M24, A24) where A24[v][i][j] = list[(pre,coef)] term-list."""
    src = json.load(open(src_path))
    M41 = [tuple(m) for m in src['masters']]
    SEC = [sum(1 << k for k in range(6) if m[k] > 0) for m in M41]
    TOP = [i for i in range(len(M41)) if SEC[i] == 63]
    M24 = [M41[i] for i in TOP]
    IDX24 = {m: i for i, m in enumerate(M24)}
    red = parse_kira2math(kira_path, _FAMILY)
    A24 = {'s': [[[] for _ in range(24)] for _ in range(24)],
           't': [[[] for _ in range(24)] for _ in range(24)]}
    miss = set()
    for ii, I in enumerate(TOP):
        for v in ('s', 't'):
            for pre, idx in src['rows_source'][str(I + 1)][v]:
                idx = tuple(idx)
                if not all(idx[p] >= 1 for p in range(6)):
                    continue   # sub-sector → 0 on the maxcut
                if idx in IDX24:
                    A24[v][ii][IDX24[idx]].append((pre, '1'))
                elif idx in red:
                    for rhs, cs in red[idx]:
                        A24[v][ii][IDX24[rhs]].append((pre, cs))
                else:
                    miss.add(idx)
    if miss:
        raise RuntimeError(f"{len(miss)} top-sector recipe targets unreduced: "
                           f"{sorted(miss)[:3]} …")
    return M24, A24, len(red)


def evalA(A24, v, dd, ss, tt):
    M = mp.zeros(24, 24)
    for i in range(24):
        for j in range(24):
            for pre, cs in A24[v][i][j]:
                M[i, j] += _fc(pre, dd, ss, tt) * _fc(cs, dd, ss, tt)
    return M


# ---- 2. integrability ----------------------------------------------------
def integrability(A24, points, dps=150, fd_h_exp=50):
    mp.mp.dps = dps
    h = mp.mpf(10) ** (-fd_h_exp)
    out = []
    for (d_, s_, t_) in points:
        d0, s0, t0 = mp.mpf(d_), mp.mpf(s_), mp.mpf(t_)
        As0 = evalA(A24, 's', d0, s0, t0)
        At0 = evalA(A24, 't', d0, s0, t0)
        dAs = (evalA(A24, 's', d0, s0, t0+h) - evalA(A24, 's', d0, s0, t0-h)) / (2*h)
        dAt = (evalA(A24, 't', d0, s0+h, t0) - evalA(A24, 't', d0, s0-h, t0)) / (2*h)
        comm = As0*At0 - At0*As0
        R = dAs - dAt + comm
        nrm = lambda M: max(abs(M[i, j]) for i in range(24) for j in range(24))
        res, c, a, b = nrm(R), nrm(comm), nrm(As0), nrm(At0)
        out.append({
            'd': d_, 's': s_, 't': t_,
            'residue': mp.nstr(res, 12),
            'curl_norm': mp.nstr(nrm(dAs - dAt), 12),
            'comm_norm': mp.nstr(c, 12),
            'As_norm': mp.nstr(a, 12), 'At_norm': mp.nstr(b, 12),
            'rel_digits': float(mp.nstr(-mp.log10(res / max(c, mp.mpf(1))), 4))
                          if res > 0 else 999.,
        })
    return out


# ---- 3. boundary (route a: per-master ε-series; route b: Γ-product) ------
def _num(v):
    """Real number from a spec value: mpf/int/float pass through; strings may
    be decimal or `a/b` rationals (mpf's parser takes both)."""
    return mp.mpf(v)


def _cnum(v):
    """Possibly-complex spec value: an [re, im] pair or a real per _num."""
    if isinstance(v, (list, tuple)):
        return mp.mpc(_num(v[0]), _num(v[1]))
    return _num(v)


def mercedes_corner(eps):
    """i·CONV·2^{2-ε}/(π² Γ(1-ε) Γ(3-2ε)),  CONV=i  — the (s,t)-CONST corner
    maxcut, validated to 44.3 digits against an independent top-cut evaluation."""
    return -mp.power(2, 2 - eps) / (mp.pi**2 * mp.gamma(1 - eps) * mp.gamma(3 - 2*eps))


def boundary_y0(corner, eps):
    """Step 3: boundary vector y0(ε) at the basepoint (s₀,t₀) from the corner
    spec (format in the module docstring).  Route (a): explicit per-master
    ε-series, y0_i(ε) = Σ_k series[i][k]·ε^{k-n_neg}.  Route (b): the
    corner-limit LS — mercedes_corner(ε) times a per-master weight vector."""
    if corner.get('route', 'a') == 'b':
        mc = mercedes_corner(eps)
        return [mc * _cnum(w) for w in corner['weights']]
    n_neg = int(corner.get('n_neg', 0))
    out = []
    for ser in corner['series']:
        acc = mp.mpc(0)
        for k, c in enumerate(ser):
            acc += _cnum(c) * eps**(k - n_neg)
        out.append(acc)
    return out


# ---- 4. transport (Bulirsch–Stoer stepper) --------------------------------
def transport_one_eps(A24, y0, var, x0, x1, fixed, eps, prec, n=24):
    """Integrate ∂_var Y = A_var(var,fixed,d=4-2ε)·Y from x0→x1 at fixed
    other-var.  BS / modified-midpoint + Richardson, tol=10^{-(prec-12)}."""
    mp.mp.dps = prec
    # mpmathify, not mpf: the ε-FFT feeds complex samples on the |ε|=R circle
    d_val = mp.mpf(4) - 2*mp.mpmathify(eps)
    if var == 's':
        Afun = lambda x: evalA(A24, 's', d_val, x, fixed)
    else:
        Afun = lambda x: evalA(A24, 't', d_val, fixed, x)
    def rhs(x, y):
        return Afun(x) * y
    tol = mp.mpf(10) ** (-(prec - 12))
    h = mp.mpf(x1 - x0) / 4
    x = mp.mpf(x0); y = mp.matrix(y0)
    def mmid(x, y, H, ns):
        hh = H / ns
        z0 = y; z1 = y + hh * rhs(x, y)
        for k in range(1, ns):
            z2 = z0 + 2*hh * rhs(x + k*hh, z1); z0, z1 = z1, z2
        return (z0 + z1 + hh * rhs(x + H, z1)) / 2
    while abs(x1 - x) > tol:
        H = (x1 - x) if abs(x1 - x) < abs(h) else h
        T = [mmid(x, y, H, 2)]; conv = False
        for kk in range(1, 12):
            ns = 2 * (kk + 1); T.append(mmid(x, y, H, ns))
            for j in range(kk, 0, -1):
                T[j-1] = T[j] + (T[j] - T[j-1]) / ((mp.mpf(ns)/(2*j))**2 - 1)
            err = max(abs(T[0][i] - T[1][i]) for i in range(n))
            if err < tol: conv = True; break
        if conv:
            y = T[0]; x += H; h = H * mp.mpf('1.5')
        else:
            h = H / 2
    return y


# ---- 5/6. ε-FFT + gate ---------------------------------------------------
def eps_fft(values, eps_samples, n_neg=0, n_pos=6):
    """values[k] at ε=eps_samples[k] on |ε|=R circle → Laurent coeffs ε^{-n_neg..n_pos-1}."""
    N = len(values); R = abs(eps_samples[0])
    out = []
    for p in range(-n_neg, n_pos):
        out.append(sum(values[k] * mp.expj(-2*mp.pi*k*p/N)
                       for k in range(N)) / (N * R**p))
    return out


def gate_laurent(coeffs, basis, n_neg, gate_digits):
    """Step 6: gate the extracted Laurent coefficients against independently
    derived basis values (format in the module docstring).  Matched digits per
    entry: −log10(|got−ref| / max(|ref|,1)); each gated order keeps its worst.
    PASS iff every gated order clears the digit threshold (basis['digits']
    overrides the CLI default)."""
    orders = {0: basis['eps0']} if 'eps0' in basis else {}
    for key, vec in basis.get('laurent', {}).items():
        orders[int(key)] = vec
    if not orders:
        raise RuntimeError("basis spec gates nothing: no 'eps0'/'laurent' key")
    n_pos = len(coeffs[0]) - n_neg
    need = float(basis.get('digits', gate_digits))
    rows, worst = [], 999.0
    for p in sorted(orders):
        if not (-n_neg <= p < n_pos):
            raise RuntimeError(f"gated order eps^{p} outside the FFT window "
                               f"eps^{-n_neg}..eps^{n_pos-1}")
        wd = 999.0
        for i, ref_raw in enumerate(orders[p]):
            ref = _cnum(ref_raw); got = coeffs[i][p + n_neg]
            err = abs(got - ref) / max(abs(ref), mp.mpf(1))
            wd = min(wd, float(-mp.log10(err)) if err > 0 else 999.0)
        rows.append({'order': p, 'min_digits': round(wd, 2)})
        worst = min(worst, wd)
    return {'orders': rows, 'digits_required': need,
            'min_digits': round(worst, 2), 'PASS': worst >= need}


# ---- driver: steps 3–6 ---------------------------------------------------
def transport_steps_3_6(A24, corner, basis, neps=16, radius='1/16', prec=60,
                        n_neg=0, n_pos=6, gate_digits=20, verbose=True):
    """Steps 3–6, run once integrability passes: boundary y0(ε) at the corner
    basepoint → per-ε BS transport s₀→s₁ (at t=t₀) then t₀→t₁ (at s=s₁) →
    ε-FFT on the |ε|=radius circle (neps samples, orders ε^{-n_neg}…
    ε^{n_pos-1}) → gate against the independent basis values.

    Returns (report_block, coeffs) where coeffs[i] is master i's Laurent
    coefficient list.  Aliasing: an order-p coefficient is clean only while
    |a_{p+neps}|·radius^neps sits below the gate tolerance — raise neps or
    shrink radius when the gate rides that floor."""
    if neps < n_neg + n_pos:
        raise RuntimeError(f"neps={neps} < n_neg+n_pos={n_neg + n_pos}: "
                           "FFT window undersampled")
    mp.mp.dps = prec
    R = _num(radius)
    s0, t0 = _num(corner['point']['s']), _num(corner['point']['t'])
    s1, t1 = _num(basis['point']['s']), _num(basis['point']['t'])
    eps_samples = [R * mp.expj(2*mp.pi*k/neps) for k in range(neps)]
    vals, tw = [], time.time()
    for k, eps in enumerate(eps_samples):
        y = boundary_y0(corner, eps)
        y = transport_one_eps(A24, y, 's', s0, s1, t0, eps, prec)
        y = transport_one_eps(A24, y, 't', t0, t1, s1, eps, prec)
        vals.append(y)
        if verbose:
            print(f"    eps sample {k + 1}/{neps} transported "
                  f"({time.time() - tw:.1f}s)", flush=True)
    coeffs = [eps_fft([vals[k][i] for k in range(neps)], eps_samples,
                      n_neg=n_neg, n_pos=n_pos) for i in range(24)]
    gate = gate_laurent(coeffs, basis, n_neg, gate_digits)
    rpt = {
        'step3_boundary': {'route': corner.get('route', 'a'),
                           'point': {'s': str(corner['point']['s']),
                                     't': str(corner['point']['t'])}},
        'step4_transport': {'neps': neps, 'radius': mp.nstr(R, 8),
                            'prec': prec,
                            'target': {'s': str(basis['point']['s']),
                                       't': str(basis['point']['t'])},
                            'wall_s': round(time.time() - tw, 1)},
        'step5_eps_fft': {'orders': f'eps^{-n_neg}..eps^{n_pos - 1}',
                          'eps0': [[mp.nstr(coeffs[i][n_neg].real, 20),
                                    mp.nstr(coeffs[i][n_neg].imag, 20)]
                                   for i in range(24)]},
        'step6_gate': gate,
    }
    return rpt, coeffs


# ---- selftest (battery leg: toy fixture, no external inputs) -------------
_TOY_C = [2, 1, -1, 0] * 6   # diagonal charges c_i (t-side uses b_i = 2c_i+1)


def _toy_A24():
    """Synthetic commuting 24×24 pair with a closed-form transport solution:
        A_s = ε·(M/s + I/(s+t)),   A_t = ε·(N/t + I/(s+t)),
        M = diag(c) + E₀₁,   N = 2M + 1   (⇒ [M,N] = 0 ⇒ integrable; the
    shared I/(s+t) term makes ∂_t A_s = ∂_s A_t ≠ 0, so step 2's curl legs do
    real work).  Along the driver's path (s₀→s₁ at t₀, then t₀→t₁ at s₁):
        Y = ((s₁+t₁)/(s₀+t₀))^ε · e^{ε·ln(t₁/t₀)·N} · e^{ε·ln(s₁/s₀)·M} · y0,
    and the E₀₁ coupling closes in 2×2 upper-triangular exponentials."""
    A = {'s': [[[] for _ in range(24)] for _ in range(24)],
         't': [[[] for _ in range(24)] for _ in range(24)]}
    for i, c in enumerate(_TOY_C):
        if c:
            A['s'][i][i].append(('1', f'({c})*(4-d)/(2*s)'))
        A['s'][i][i].append(('1', '(4-d)/(2*(s+t))'))
        A['t'][i][i].append(('1', f'({2 * c + 1})*(4-d)/(2*t)'))
        A['t'][i][i].append(('1', '(4-d)/(2*(s+t))'))
    A['s'][0][1] = [('1', '(4-d)/(2*s)')]    # E₀₁ in M
    A['t'][0][1] = [('1', '(4-d)/t')]        # 2·E₀₁ in N
    return A


def selftest():
    """Battery for the steps 3–6 driver (toy fixture, seconds-scale;
    portability, not physics).  Legs:
      T1 boundary route (a): ε-series eval vs a hand-summed polynomial
      T2 step 2 on the toy pair: integrability PASS with nonzero curl
      T3 steps 3–6 end-to-end: route (b) boundary → transport → ε-FFT,
         gated ≥10 digits on Laurent orders 0–3 vs closed-form references
      T4 negative control: one perturbed basis value must FAIL the gate
    Step 1 (assemble from a reduction file) is NOT covered — it needs
    recipe/reduction inputs that do not ship.  Returns 0 all-PASS, else 1."""
    ok = True

    def leg(name, cond, detail=''):
        nonlocal ok
        print(f"[selftest:{name}] {'PASS' if cond else 'FAIL'}  {detail}",
              flush=True)
        ok = ok and cond

    prec = 40
    mp.mp.dps = prec

    # T1: route (a) boundary series  y0(ε) = ε⁻¹ + 2 + 3ε at ε=1/4
    corner_a = {'point': {'s': '-1', 't': '-1'}, 'route': 'a', 'n_neg': 1,
                'series': [[['1', '0'], ['2', '0'], ['3', '0']]] * 24}
    e = mp.mpf(1) / 4
    got = boundary_y0(corner_a, e)[7]
    want = 1/e + 2 + 3*e
    leg('T1-boundary-a', abs(got - want) < mp.mpf(10)**(-30),
        f'|delta|={mp.nstr(abs(got - want), 3)}')

    A24 = _toy_A24()

    # T2: integrability on the toy pair (curl terms nonzero, residue ~ FD err)
    integ = integrability(A24, [('37/10', -3, -5), ('29/10', -9, -6)],
                          dps=80, fd_h_exp=25)
    leg('T2-integrability', all(p['rel_digits'] > 30 for p in integ),
        f"rel_digits={[p['rel_digits'] for p in integ]}")

    # T3: steps 3–6 end-to-end, gated vs closed-form Taylor references
    s0 = t0 = mp.mpf(-1); s1, t1 = mp.mpf('-6/5'), mp.mpf('-11/10')
    Ls, Lt = mp.log(s1/s0), mp.log(t1/t0)
    Lu = mp.log((s1 + t1)/(s0 + t0))

    def ref_fn(i):
        c = _TOY_C[i]

        def f(e):
            mc = mercedes_corner(e); cw = mp.exp(e*Lu)
            if i >= 2:
                return mc * cw * mp.exp(e*(c*Ls + (2*c + 1)*Lt))
            # E₀₁ block: M₀₁=[[2,1],[0,1]], N₀₁=[[5,2],[0,3]]
            xs, xt = e*Ls, e*Lt
            ys0 = (mp.exp(2*xs) + (mp.exp(2*xs) - mp.exp(xs))) * mc
            ys1 = mp.exp(xs) * mc
            yt0 = mp.exp(5*xt)*ys0 + (mp.exp(5*xt) - mp.exp(3*xt))*ys1
            yt1 = mp.exp(3*xt)*ys1
            return cw * (yt0 if i == 0 else yt1)
        return f

    with mp.workdps(100):
        refs = [mp.taylor(ref_fn(i), 0, 3) for i in range(24)]
    basis = {'point': {'s': '-6/5', 't': '-11/10'},
             'eps0': [[mp.nstr(refs[i][0], 40), '0'] for i in range(24)],
             'laurent': {str(p): [[mp.nstr(refs[i][p], 40), '0']
                                  for i in range(24)] for p in (1, 2, 3)},
             'digits': 10}
    corner_b = {'point': {'s': '-1', 't': '-1'}, 'route': 'b',
                'weights': [['1', '0']] * 24}
    r36, coeffs = transport_steps_3_6(A24, corner_b, basis, neps=8,
                                      radius='1/32', prec=prec,
                                      n_neg=0, n_pos=4, verbose=False)
    g = r36['step6_gate']
    leg('T3-steps3-6', g['PASS'],
        f"min_digits={g['min_digits']} (need {g['digits_required']:.0f}); "
        f"per-order={[(o['order'], o['min_digits']) for o in g['orders']]}; "
        f"wall={r36['step4_transport']['wall_s']}s")

    # T4: negative control — a perturbed basis value must FAIL the gate
    bad = json.loads(json.dumps(basis))
    with mp.workdps(100):
        bad['eps0'][7][0] = mp.nstr(refs[7][0] + mp.mpf('1e-6'), 40)
    g_bad = gate_laurent(coeffs, bad, 0, 10)
    leg('T4-negative-control', not g_bad['PASS'],
        f"perturbed min_digits={g_bad['min_digits']} (need "
        f"{g_bad['digits_required']:.0f}) -> gate "
        f"{'FAILS as required' if not g_bad['PASS'] else 'passed (BUG)'}")

    print(f"selftest: {'PASS' if ok else 'FAIL'} (4 legs)", flush=True)
    return 0 if ok else 1


# ---- main ---------------------------------------------------------------
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--dps', type=int, default=150)
    ap.add_argument('--out', default=DEFAULTS['out'])
    ap.add_argument('--corner', default=DEFAULTS['corner'],
                    help='step-3 boundary spec JSON (format: module docstring)')
    ap.add_argument('--basis', default=DEFAULTS['basis'],
                    help='step-6 independent basis spec JSON')
    ap.add_argument('--neps', type=int, default=16,
                    help='eps samples on the FFT circle (step 4)')
    ap.add_argument('--eps-radius', default='1/16',
                    help='FFT circle radius |eps| (decimal or a/b)')
    ap.add_argument('--prec', type=int, default=60,
                    help='transport working precision (digits)')
    ap.add_argument('--n-neg', type=int, default=0,
                    help='deepest extracted pole order eps^{-n_neg}')
    ap.add_argument('--n-pos', type=int, default=6,
                    help='extract orders up to eps^{n_pos-1}')
    ap.add_argument('--gate-digits', type=float, default=20,
                    help='matched digits every gated order must clear')
    ap.add_argument('--selftest', action='store_true',
                    help='toy-fixture battery for the steps 3-6 driver')
    a = ap.parse_args()

    if a.selftest:
        sys.exit(selftest())

    t0 = time.time()
    print(f"[1] assemble 24×24 from {DEFAULTS['kira_red']} ...", flush=True)
    M24, A24, nred = load_de24()
    print(f"    {nred} reductions; assemble OK", flush=True)

    print(f"[2] integrability at 3 (d,s,t) points, dps={a.dps} ...", flush=True)
    integ = integrability(A24, [('37/10', -3, -5),
                                ('41/10', -7, -3),
                                ('29/10', -9, -6)], dps=a.dps)
    for p in integ:
        print(f"    (d,s,t)=({p['d']},{p['s']},{p['t']}): "
              f"res={p['residue']}  comm={p['comm_norm']}  "
              f"rel_d={p['rel_digits']:.2f}", flush=True)
    PASS = all(p['rel_digits'] > 30 for p in integ)

    rpt = {
        'date': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'masters24': [list(m) for m in M24],
        'step1_assemble': {'n_reductions': nred, 'missing': 0},
        'step2_integrability': {
            'dps': a.dps, 'fd_h': '1e-50',
            'points': integ,
            'PASS': PASS,
        },
        'elapsed_s': round(time.time() - t0, 1),
    }
    if not PASS:
        rpt['VERDICT'] = ('INTEGRABILITY FAILS — see OPERATOR CONTAMINATION '
                          'in this file docstring and the written report. '
                          'Steps 3–6 NOT executed.')
        json.dump(rpt, open(a.out, 'w'), indent=1)
        print(json.dumps(rpt['step2_integrability'], indent=1))
        sys.exit(2)

    # ---- steps 3–6: boundary → transport → ε-FFT → gate ------------------
    for what, path in (('corner (step 3)', a.corner), ('basis (step 6)', a.basis)):
        if not os.path.isfile(path):
            raise SystemExit(
                f"integrability PASSES but the {what} spec is missing: {path} "
                "— point --corner/--basis (or MAXCUT_DE_CORNER/MAXCUT_DE_BASIS) "
                "at your own files; formats in the module docstring")
    corner = json.load(open(a.corner))
    basis = json.load(open(a.basis))
    print(f"[3] boundary route ({corner.get('route', 'a')}) at "
          f"(s0,t0)=({corner['point']['s']},{corner['point']['t']})", flush=True)
    print(f"[4] transport {a.neps} eps samples on |eps|={a.eps_radius}, "
          f"prec={a.prec}, target (s1,t1)=({basis['point']['s']},"
          f"{basis['point']['t']}) ...", flush=True)
    res36, _ = transport_steps_3_6(A24, corner, basis, neps=a.neps,
                                   radius=a.eps_radius, prec=a.prec,
                                   n_neg=a.n_neg, n_pos=a.n_pos,
                                   gate_digits=a.gate_digits)
    rpt.update(res36)
    g = res36['step6_gate']
    print(f"[5] eps-FFT extracted orders {res36['step5_eps_fft']['orders']}",
          flush=True)
    print(f"[6] gate vs {a.basis}: min_digits={g['min_digits']} "
          f"(need {g['digits_required']:.0f}) -> "
          f"{'PASS' if g['PASS'] else 'FAIL'}", flush=True)
    rpt['VERDICT'] = ('GATED — steps 1–6 complete, gate PASS' if g['PASS'] else
                      'GATE FAILS — transported eps-coefficients do not match '
                      'the independent basis; see step6_gate')
    rpt['elapsed_s'] = round(time.time() - t0, 1)
    json.dump(rpt, open(a.out, 'w'), indent=1)
    print(f"report -> {a.out}", flush=True)
    sys.exit(0 if g['PASS'] else 3)

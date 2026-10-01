#!/usr/bin/env python3
"""svmap21.py — ALGEBRAIC single-valued map zeta_sv on depth<=3 odd-index
MZVs (this package's conventions), Brown's sv projection computed exactly over the
GALOIS relation ring. ADDITIVE MODULE. galois core.py / vspace.py / sector21.py
untouched: this module imports galois.sector21 (Goncharov D_m machinery) and
galois.mzvring (word algebra) — the same module objects, never copies.

Route: F. Brown,
arXiv:1309.5309 sect 7.2 eq (7.3) — the f-alphabet model
  H ~= Q<f3,f5,f7,...> (x) Q[f2],  sv(f2) = 0,
  sv(w) = sum_{uv=w} u sh v~   [~ = word reversal; Brown's PREFIX-strip
                                model = orient A here]
transported to this package's f-model, which is TRAILING-strip (U1 gate-c
locks: rho(zeta(3,5)) = -5 f3f5, this package's D_m strips the trailing f_m), so
the reversal moves to the LEFT deconcatenation factor:
  sv_camp(w) = sum_{uv=w} u~ sh v   [orient B — this package's transport]
MIRROR-ORIENTATION EMPIRICAL PIN: orient B reproduces the 61-formula
SV_FORMULAS.json bank 61/61 EXACTLY; verbatim orient A matches only the
palindromic keys (5,5,5),(7,7,7) where the orientations provably coincide
(2/61). Both derived and independently verified before use — gate_regress61 re-enforces
this every selftest. The two agree at depth <= 2.

phi (f-alphabet decomposition of the ring): phi(pi2) = f2, phi(z_m) = f_m;
survivor generators m2/m3 by coeff_{u.f_m}(x) = coeff_u(D_m x) with D_m =
the sector21 Goncharov derivation. The ONE non-recursive convention: the
pure-f2^{w/2} coefficient of an even-weight survivor generator := 0 (the
f2-analogue of the pi_m splitting; validated through every stored
z2^j coefficient by the 61/61 gate).

zeta_sv[campaign](c) := phi^-1( sv_camp( phi( mzv(c) ) ) ), exact Fraction
linear solve over ALL weight-w ring canonical monomials. The ring is
depth<=3-truncated, so for w >= 12 the monomial count is SMALLER than the
Zagier dim d_w (deficit = depth>=4 generators, first zeta(6,4,1,1) at w=12):
the solve is RECTANGULAR with a full-column-rank assert (phi injective on
the ring) + exact per-target consistency assert — a target needing depth>=4
content HARD-FAILS (reported, never smoothed) — + a phi-roundtrip assert on
every solution.

DATA (tools link, never copy; every entry point takes an
explicit dir/path): SVALIGN default = <GALOIS_CAMPAIGN_BANK>/svalign
(a reference bank, not shipped) holding SV_FORMULAS.json (the
61-formula regression bank), SV_FORMULAS_EXT.json (the 22 ALGEBRAIC-DERIVED
blocked keys), and w_coeffs_700.json (the
INDEPENDENT W-associator sv numeric bank, ~700d, keyed Brown-order — the
numeric check source). Ring pkls: galois.mzvring.RING_BANK (ring21.pkl
default, via sector21.init_ring).

GATES (`selftest --receipt PATH` writes a receipt on full PASS; no receipt
is vendored — the gates re-run against your own bank):
  pin       — 4 literature phi locks (Brown arXiv:1102.1310: coeff(f3f5) in
              phi(m2(5,3)) = -5, coeff(f5f3) = 0; arXiv:2511.15883:
              coeff(f3f7) = -14, coeff(f5f5) = -6 in phi(m2(7,3)));
              Brown's printed (7.3) depth-3 example reproduced VERBATIM by
              orient A + its mirror by orient B; sv(f2) = 0.
  depth2    — all 27 all-odd depth-2 pairs w=8..18: A == B, pi2-free, pure
              odd-zeta products, exact lock zeta_sv[campaign](5,3) =
              -10 z3 z5, numerics vs the sv bank at 460d (measured: min
              480.7d, 27/27).
  regress61 — HARD GATE: the 61 reference SV_FORMULAS reproduced EXACTLY under
              orient B (61/61 or refuse to sign); orient A counts reported
              as the orientation adjudication.
  extcheck  — the 22 reference EXT formulas == in-process re-derivation, exact.

CLI:
  python3 svmap21.py selftest [--ring PKL] [--svalign DIR] [--dps N]
                              [--receipt J]     -> pin+depth2+regress61+ext
  python3 svmap21.py derive --comps "5,7,7;3,9,7" [--svalign DIR] [--dps N]
                              [--out J]         -> formulas (+numeric check)

Downstream membership tests (exact Fraction span-membership with witness
tracking, plant/wrong-slot sanity checks) build directly on
derive_formulas(); no harness for them ships here.
"""
import ast
import json
import os
import re
import sys
from fractions import Fraction as Fr

TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if TOOLS not in sys.path:
    sys.path.insert(0, TOOLS)

from galois import sector21 as S21                       # BY IDENTITY
from galois.sector21 import (sym_vec, dmul, dsum, dscale, to_name_form,  # noqa
                             Dm_comp, Dm_vec, sv_proj, parse_brown_name,
                             ms_for)
from galois.mzvring import znorm, mweight, shuffle as wshuffle  # noqa

CAMPAIGN = os.environ.get('GALOIS_CAMPAIGN_BANK')  # reference-table root (not shipped; point at your own)
SVALIGN = os.path.join(CAMPAIGN, 'svalign') if CAMPAIGN else None


def _bank(p, what):
    """fail-closed default-path guard (see sector21._bank)."""
    if p is None:
        raise RuntimeError(
            f"{what}: no explicit svalign dir given and GALOIS_CAMPAIGN_BANK "
            "is not set. The SV formula/numeric banks are reference "
            "data that do not ship with this repo — pass --svalign / "
            "svalign_dir explicitly or point GALOIS_CAMPAIGN_BANK at your own "
            "bank (subdir svalign/).")
    return p
CHECK_DPS = 460


def init_ring(path=None):
    """Bind the reduction oracle (forwards to sector21.init_ring) and clear
    the ring-keyed phi/monomial caches HERE. Use THIS when svmap21 is in
    play — calling sector21.init_ring directly leaves these caches stale."""
    R = S21.init_ring(path)
    _PHI.clear()
    _PHIM.clear()
    _MONS.clear()
    _COLS.clear()
    _BOUND['path'] = S21._STATE['path']
    return R


def _ring():
    R = S21._ring()
    if _BOUND['path'] != S21._STATE['path']:      # rebound behind our back
        _PHI.clear()
        _PHIM.clear()
        _MONS.clear()
        _COLS.clear()
        _BOUND['path'] = S21._STATE['path']
    return R


_BOUND = {'path': None}


# ---------------------------------------------------------- f-alphabet vectors
# fvec = {(j, word): Fraction}, word = tuple of odd ints >= 3, j = f2 exponent.

def fstrip(d):
    return {k: q for k, q in d.items() if q}


def fsum(*ds):
    out = {}
    for d in ds:
        for k, q in d.items():
            out[k] = out.get(k, Fr(0)) + q
    return fstrip(out)


def fscale(d, c):
    return {k: q * c for k, q in d.items()} if c else {}


def fmul(A, B):
    out = {}
    for (ja, ua), qa in A.items():
        for (jb, ub), qb in B.items():
            q = qa * qb
            for w, c in wshuffle(ua, ub).items():
                k = (ja + jb, w)
                out[k] = out.get(k, Fr(0)) + q * c
    return fstrip(out)


def fweight(k):
    return 2 * k[0] + sum(k[1])


# ------------------------------------------------------------------- phi map

_PHI = {}      # symbol -> fvec
_PHIM = {}     # monomial -> fvec


def phi_sym(s):
    """phi of a single ring generator symbol (pi2 / z odd / m2,m3 survivor)."""
    _ring()
    if s in _PHI:
        return _PHI[s]
    if s == ('pi2',):
        v = {(1, ()): Fr(1)}
    elif s[0] == 'z':
        assert s[1] % 2 == 1, s
        v = {(0, (s[1],)): Fr(1)}
    else:
        c = tuple(s[1:])
        w = sum(c)
        v = {}
        for m in range(3, w + 1, 2):
            img = Dm_comp(c, m)
            if not img:
                continue
            for (j, u), q in phi_vec(img).items():
                k = (j, u + (m,))
                v[k] = v.get(k, Fr(0)) + q
        v = fstrip(v)
        # pure-f2 part of an even-weight generator: 0 by the package's
        # splitting convention (never added).
    for k in v:
        assert fweight(k) == (2 if s == ('pi2',) else sum(s[1:])), (s, k)
    _PHI[s] = v
    return v


def phi_mono(mono):
    if mono in _PHIM:
        return _PHIM[mono]
    out = {(0, ()): Fr(1)}
    for s in mono:
        out = fmul(out, phi_sym(s))
    _PHIM[mono] = out
    return out


def phi_vec(vec):
    out = {}
    for mono, q in vec.items():
        for k, c in phi_mono(mono).items():
            out[k] = out.get(k, Fr(0)) + q * c
    return fstrip(out)


# --------------------------------------------------------- single-valued map

def sv_f(fv, orient='B'):
    """Brown (7.3) transported to outer-first words. orient 'B' (this
    package's transport: reversal on the LEFT deconcatenation factor — the
    61/61-pinned choice) or 'A' (Brown's verbatim prefix-model orientation,
    reversal on the RIGHT factor). sv(f2)=0: j>0 components die."""
    out = {}
    for (j, u), q in fv.items():
        if j:
            continue
        for k in range(len(u) + 1):
            a, b = u[:k], u[k:]
            sh = (wshuffle(tuple(reversed(a)), b) if orient == 'B'
                  else wshuffle(a, tuple(reversed(b))))
            for w, c in sh.items():
                key = (0, w)
                out[key] = out.get(key, Fr(0)) + q * c
    return fstrip(out)


# ------------------------------------------- weight-w ring monomial basis

def mono_key(symbols):
    v = {(): Fr(1)}
    for s in symbols:
        v = dmul(v, {(s,): Fr(1)})
    (k, c), = v.items()
    assert c == 1
    return k


def gens_upto(w):
    R = _ring()
    G = [(('pi2',), 2)]
    for n in range(3, w + 1, 2):
        G.append((('z', n), n))
    for wv in sorted(R.survivors):
        if wv > w:
            continue
        for s in sorted(R.survivors[wv]):
            G.append((s, wv))
    return G


_MONS = {}


def monomials(w):
    """all ring canonical basis monomials of weight w. NOTE: the ring is
    depth<=3-truncated, so for w >= 12 this is SMALLER than the Zagier
    dimension d_w (the deficit = depth>=4 generators, first one
    zeta(6,4,1,1) at w=12); the solver is therefore rectangular with
    per-target consistency checks."""
    _ring()
    if w in _MONS:
        return _MONS[w]
    G = gens_upto(w)
    out = []

    def rec(i, rem, cur):
        if rem == 0:
            out.append(mono_key(cur))
            return
        if i >= len(G):
            return
        rec(i + 1, rem, cur)
        s, wt = G[i]
        k = 1
        while k * wt <= rem:
            rec(i + 1, rem - k * wt, cur + [s] * k)
            k += 1
    rec(0, w, [])
    _MONS[w] = sorted(set(out))
    return _MONS[w]


# -------------------------------------------------------------- exact solver

_COLS = {}


def weight_cols(w):
    if w not in _COLS:
        _COLS[w] = [(m, phi_mono(m)) for m in monomials(w)]
    return _COLS[w]


def solve_many(w, targets):
    """solve phi(x_t) = t for each target fvec t; returns list of canonical
    vectors. Rectangular exact solve (the ring is depth<=3-truncated, so
    n_coords >= n_cols for w >= 12). Asserts: full COLUMN rank (phi
    injective on the ring), exact per-target consistency (a residual = the
    target needs depth>=4 content -> hard failure, reported), and a
    phi-roundtrip on every solution."""
    cols = weight_cols(w)
    coords = set()
    for _, fv in cols:
        coords.update(fv)
    for t in targets:
        for k in t:
            assert fweight(k) == w, (w, k)
        coords.update(t)
    coords = sorted(coords)
    n = len(cols)
    nt = len(targets)
    rows = []
    for k in coords:
        row = [fv.get(k, Fr(0)) for _, fv in cols]
        row += [t.get(k, Fr(0)) for t in targets]
        rows.append(row)
    # Gauss-Jordan on the first n columns
    piv_of_col = [None] * n
    r = 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        piv_of_col[c] = r
        r += 1
    assert r == n, f"phi matrix at w={w} rank {r} < {n} (not injective)"
    sols = []
    for ti in range(nt):
        # consistency: rows below the rank must carry zero rhs
        bad = [coords[i] for i in range(r, len(rows)) if rows[i][n + ti] != 0]
        assert not bad, (f"target {ti} at w={w} OUTSIDE the depth<=3 ring "
                         f"span; residual f-coords: {bad[:8]}")
        x = {}
        for c in range(n):
            q = rows[piv_of_col[c]][n + ti]
            if q:
                x[cols[c][0]] = q
        # exact roundtrip: phi(x) == target
        chk = {}
        for mono, q in x.items():
            for k, cc in phi_mono(mono).items():
                chk[k] = chk.get(k, Fr(0)) + q * cc
        assert fstrip(chk) == fstrip(targets[ti]), \
            f"phi roundtrip FAIL at w={w} target {ti}"
        sols.append(x)
    return sols


# ------------------------------------------------------------- sv reduction

def sv_reduce_many(comps, orient='B'):
    """comps: outer-first compositions (all same weight, depth 2 or 3).
    Returns {comp: canonical vector of zeta_sv[campaign](comp)}."""
    w = sum(comps[0])
    targets = []
    for c in comps:
        assert sum(c) == w
        s = (('m2',) if len(c) == 2 else ('m3',)) + tuple(c)
        targets.append(sv_f(phi_vec(sym_vec(s)), orient))
    sols = solve_many(w, targets)
    return dict(zip(comps, sols))


# --------------------------------------------------------- numeric machinery

_NCACHE = {}


def _val(s, dps):
    import mpmath as mp
    from formglue import form_oracle
    key = (s, dps)
    if key not in _NCACHE:
        if s == ('pi2',):
            _NCACHE[key] = mp.pi ** 2
        elif s[0] == 'z':
            _NCACHE[key] = mp.zeta(s[1])
        else:
            _NCACHE[key] = form_oracle.to_mpf(
                form_oracle.mzv(tuple(s[1:]), dps), dps)
    return _NCACHE[key]


def eval_vec(vec, dps):
    import mpmath as mp
    mp.mp.dps = dps + 20
    tot = mp.mpf(0)
    for mono, q in vec.items():
        t = mp.mpf(q.numerator) / q.denominator
        for s in mono:
            t *= _val(s, dps)
        tot += t
    return tot


def sv_bank(svalign_dir=None):
    """the W-associator sv numeric bank (~700d, keyed Brown-order)."""
    return json.load(open(os.path.join(svalign_dir or _bank(SVALIGN, '_load_bank'),
                                       'w_coeffs_700.json')))['sv']


def agree_digits(lhs, rhs):
    import mpmath as mp
    if lhs == rhs:
        return 9999.0
    return float(-mp.log10(abs(lhs - rhs) / max(abs(lhs), mp.mpf('1e-999'))))


# ------------------------------------------------------ Brown-name emission

def brown_name(mono):
    """ring monomial -> svalign coeffs_brown_basis name + the 6^j z2 factor.
    Returns (name, scale) with coefficient_out = coefficient_in * scale."""
    j = sum(1 for s in mono if s[0] == 'pi2')
    zs = sorted((s[1] for s in mono if s[0] == 'z'), reverse=True)
    zbs = [tuple(s[1:]) for s in mono if s[0] in ('m2', 'm3')]
    parts = []
    if j:
        parts.append('z2' if j == 1 else f'z2^{j}')
    parts += [f'z{n}' for n in zs]
    for idx in sorted(zbs):
        parts.append('zB(' + ', '.join(map(str, reversed(idx))) + ')')
    assert parts, "constant monomial in an sv formula"
    return '*'.join(parts), Fr(6) ** j


def coeffs_brown(vec, camp):
    """vec minus the leading 2*mzv3(camp), as {brown name: 'p/q' string}."""
    rest = dsum(vec, dscale(sym_vec(('m3',) + tuple(camp)), Fr(-2)))
    out = {}
    for mono, q in rest.items():
        name, sc = brown_name(mono)
        assert name not in out, (name, mono)
        out[name] = str(q * sc)
    # roundtrip: 2*mzv3 + parsed coeffs == vec  (sector21 parser, raw names)
    back = dscale(sym_vec(('m3',) + tuple(camp)), Fr(2))
    for name, qs in out.items():
        back = dsum(back, dscale(parse_brown_name(name.replace(' ', '')),
                                 Fr(qs)))
    assert back == vec, f"brown-name roundtrip FAIL for {camp}"
    return out


def formula_strings(camp, coeffs):
    brown = tuple(reversed(camp))
    bkey = ','.join(map(str, brown))

    def terms(conv):
        ts = []
        for name, qs in sorted(coeffs.items()):
            q = Fr(qs)
            s = str(q) if q < 0 else '+ ' + str(q)
            s = s.replace('-', '- ', 1) if q < 0 else s
            if conv == 'campaign':
                nm = name.replace('zB(', 'MZVTMP(')
                nm = re.sub(r'MZVTMP\(([\d, ]+)\)',
                            lambda m: 'mzv' + str(len(m.group(1).split(','))) +
                            '(' + ','.join(reversed(
                                [x.strip() for x in m.group(1).split(',')]))
                            + ')', nm)
                nm = re.sub(r'\bz(\d+)(\^(\d+))?\b',
                            lambda m: f"zeta{m.group(1)}" +
                            (f"^{m.group(3)}" if m.group(3) else ''), nm)
            else:
                nm = name
            ts.append(f"{s}*{nm}")
        return ' '.join(ts)
    fb = f"zeta_sv({bkey}) = 2*zB({bkey}) " + terms('brown')
    ck = ','.join(map(str, camp))
    fc = (f"zeta_sv[campaign]({ck}) = 2*mzv3({ck}) " + terms('campaign'))
    return fb, fc


# ------------------------------------------------------------------- gates

LIT_LOCKS = [
    # (symbol, word, expected coefficient, citation)
    (('m2', 5, 3), (3, 5), Fr(-5),
     'rho(zeta(3,5)) = -5 f3f5 + c f2^4, Brown arXiv:1102.1310'),
    (('m2', 5, 3), (5, 3), Fr(0), 'no word ending f3 in rho(zeta(3,5))'),
    (('m2', 7, 3), (3, 7), Fr(-14),
     'rho(zeta(3,7)) = -14 f3f7 - 6 f5f5 + q10 f10, arXiv:2511.15883'),
    (('m2', 7, 3), (5, 5), Fr(-6), 'same expansion, coeff of f5f5'),
]


def gate_pin():
    """literature phi locks + Brown (7.3) printed examples + sv(f2)=0 +
    phi built for every survivor generator."""
    R = _ring()
    locks = []
    for s, word, want, cite in LIT_LOCKS:
        got = phi_sym(s).get((0, word), Fr(0))
        ok = got == want
        locks.append({'symbol': str(s), 'word': str(word), 'want': str(want),
                      'got': str(got), 'pass': ok, 'cite': cite})
        assert ok, locks[-1]
    # sv examples lock (Brown sect 7.2 printed examples, orient-independent
    # at depth <=2; the depth-3 example fixes the orientation dictionary):
    exA = sv_f({(0, (3, 5, 7)): Fr(1)}, 'A')
    exB = sv_f({(0, (3, 5, 7)): Fr(1)}, 'B')
    brown_words = {(0, w): Fr(2) for w in
                   [(3, 5, 7), (3, 7, 5), (7, 3, 5), (7, 5, 3)]}
    mirror_words = {(0, w): Fr(2) for w in
                    [(3, 5, 7), (5, 3, 7), (5, 7, 3), (7, 5, 3)]}
    assert exA == brown_words, exA     # verbatim (7.3) reproduces Brown's list
    assert exB == mirror_words, exB    # orient-B mirror
    assert sv_f({(0, (3,)): Fr(1)}, 'B') == {(0, (3,)): Fr(2)}
    assert sv_f({(0, (3, 5)): Fr(1)}, 'B') == \
        {(0, (3, 5)): Fr(2), (0, (5, 3)): Fr(2)}
    assert sv_f({(1, ()): Fr(1)}, 'B') == {}          # sv(f2) = 0
    built = {}
    for wv in sorted(R.survivors):
        for s in sorted(R.survivors[wv]):
            v = phi_sym(s)
            built[str(s)] = {'weight': wv, 'n_fterms': len(v)}
    return {'pass': True, 'literature_locks': locks,
            'brown_7p3_example_reproduced_verbatim': True,
            'sv_f2_kill': True, 'n_phi_generators_built': len(built)}


def gate_depth2(check_dps=CHECK_DPS, svalign_dir=None):
    """all 27 all-odd depth-2 pairs w=8..18: orientation agreement, pi2-free
    pure odd-zeta products, exact (5,3) lock, numerics vs the sv bank."""
    import mpmath as mp
    svb = sv_bank(svalign_dir)
    rows = []
    for w in range(8, 20, 2):
        comps = [(a, w - a) for a in range(3, w - 2, 2)]
        red = sv_reduce_many(comps, 'B')
        redA = sv_reduce_many(comps, 'A')
        for c in comps:
            vec = red[c]
            assert vec == redA[c], f"orientation depth-2 disagreement {c}"
            assert vec == sv_proj(vec), f"pi2 term in depth-2 sv {c}"
            for mono in vec:
                assert all(s[0] == 'z' for s in mono), (c, mono)
            bkey = ','.join(map(str, reversed(c)))
            mp.mp.dps = check_dps + 20
            lhs = mp.mpf(svb[bkey])
            rhs = eval_vec(vec, check_dps)
            d = agree_digits(lhs, rhs)
            rows.append({'campaign': list(c), 'brown': bkey,
                         'agree_digits': round(d, 1)})
            assert d > check_dps - 30, rows[-1]
    # the classic exact lock: zeta_sv[campaign](5,3) = -10 z3 z5
    v53 = sv_reduce_many([(5, 3)], 'B')[(5, 3)]
    assert v53 == {mono_key([('z', 3), ('z', 5)]): Fr(-10)}, v53
    return {'pass': True, 'n_checked': len(rows), 'check_dps': check_dps,
            'min_agree_digits': min(r['agree_digits'] for r in rows),
            'exact_lock': 'zeta_sv[campaign](5,3) = -10*z3*z5'}


def _load_bank_vectors(path):
    """SV_FORMULAS(-EXT).json -> {camp: (weight, canonical vector)}."""
    SVF = json.load(open(path))
    out = {}
    for k, f in SVF['formulas'].items():
        brown = ast.literal_eval(k)
        camp = tuple(reversed(brown))
        vec = dscale(sym_vec(('m3',) + camp), Fr(2))
        for name, q in f['coeffs_brown_basis'].items():
            vec = dsum(vec, dscale(parse_brown_name(name.replace(' ', '')),
                                   Fr(q)))
        out[camp] = (f['weight'], vec)
    return out


def gate_regress61(svalign_dir=None):
    """HARD GATE: the 61 reference SV_FORMULAS reproduced EXACTLY under orient
    B; orient A counts reported (the mirror-orientation adjudication —
    A matches only the palindromic keys)."""
    bank = _load_bank_vectors(os.path.join(svalign_dir or _bank(SVALIGN, 'bank_vectors'),
                                           'SV_FORMULAS.json'))
    byw = {}
    for camp, (w, vec) in bank.items():
        byw.setdefault(w, []).append(camp)
    res = {'match_B': 0, 'mismatch_B': 0, 'match_A': 0, 'mismatch_A': 0}
    mism = []
    for w in sorted(byw):
        comps = byw[w]
        derB = sv_reduce_many(comps, 'B')
        derA = sv_reduce_many(comps, 'A')
        for camp in comps:
            ref = bank[camp][1]
            okB = derB[camp] == ref
            okA = derA[camp] == ref
            res['match_B' if okB else 'mismatch_B'] += 1
            res['match_A' if okA else 'mismatch_A'] += 1
            if not okB:
                mism.append(','.join(map(str, camp)))
    n = sum(len(v) for v in byw.values())
    gate = res['mismatch_B'] == 0 and res['match_B'] == n == 61
    assert gate, {'counts': res, 'n': n, 'mismatches_B': mism,
                  'HARD_GATE': 'orient B must be 61/61 exact — REFUSE'}
    return {'pass': True, 'counts': res,
            'orientation_adjudication': 'orient B 61/61 exact; orient A '
                                        f"{res['match_A']}/61 (palindromic "
                                        'keys only)'}


def gate_extcheck(svalign_dir=None):
    """the reference EXT formulas (the 22 ALGEBRAIC-DERIVED blocked keys) ==
    in-process re-derivation under orient B, exact."""
    p = os.path.join(svalign_dir or _bank(SVALIGN, 'ext_bank'), 'SV_FORMULAS_EXT.json')
    bank = _load_bank_vectors(p)
    byw = {}
    for camp, (w, vec) in bank.items():
        byw.setdefault(w, []).append(camp)
    n = 0
    for w, comps in sorted(byw.items()):
        der = sv_reduce_many(comps, 'B')
        for camp in comps:
            assert der[camp] == bank[camp][1], \
                f"EXT bank cross-check FAIL {camp}"
            n += 1
    return {'pass': True, 'n_checked': n, 'bank': p}


# ------------------------------------------------------------------- derive

def derive_formulas(comps, dps=None, svalign_dir=None):
    """derive zeta_sv[campaign] for outer-first compositions (depth 2 or 3;
    grouped per weight), orient B; optional numeric check vs the sv bank at
    dps. Returns {brown_key_str: entry} in SV_FORMULAS_EXT format."""
    import mpmath as mp
    byw = {}
    for c in comps:
        byw.setdefault(sum(c), []).append(tuple(c))
    svb = sv_bank(svalign_dir) if dps else None
    out = {}
    for w in sorted(byw):
        der = sv_reduce_many(byw[w], 'B')
        for camp in byw[w]:
            vec = der[camp]
            coeffs = coeffs_brown(vec, camp)
            fb, fc = formula_strings(camp, coeffs)
            brown = tuple(reversed(camp))
            entry = {'weight': w, 'status': 'ALGEBRAIC-DERIVED',
                     'formula_brown_conv': fb, 'formula_campaign_conv': fc,
                     'coeffs_brown_basis': coeffs,
                     'zeta_sv_campaign_reduced': to_name_form(vec)}
            if dps:
                bkey = ','.join(map(str, brown))
                assert bkey in svb, f"no sv bank value for {bkey}"
                mp.mp.dps = dps + 20
                d = agree_digits(mp.mpf(svb[bkey]), eval_vec(vec, dps))
                entry['numeric_check'] = {'dps': dps,
                                          'agree_digits': round(d, 1),
                                          'pass': d > dps - 30}
                assert d > dps - 30, (camp, d)
            out['(' + ', '.join(map(str, brown)) + ')'] = entry
    return out


# ------------------------------------------------------------------------ main

def selftest(receipt_path=None, check_dps=CHECK_DPS, svalign_dir=None):
    """full registration battery, hard order: pin -> depth2 -> regress61 ->
    extcheck. Any failure raises; a receipt is written only on full PASS."""
    import datetime
    import time
    R = {'date': datetime.datetime.now().isoformat(timespec='seconds'),
         'module': os.path.abspath(__file__), 'ring_pkl': None,
         'svalign': svalign_dir or _bank(SVALIGN, 'selftest'), 'gates': {}}
    _ring()
    R['ring_pkl'] = S21._STATE['path']
    for name, fn in [('pin', gate_pin),
                     ('depth2', lambda: gate_depth2(check_dps, svalign_dir)),
                     ('regress61', lambda: gate_regress61(svalign_dir)),
                     ('extcheck', lambda: gate_extcheck(svalign_dir))]:
        t0 = time.time()
        R['gates'][name] = fn()
        R['gates'][name]['wall_s'] = round(time.time() - t0, 1)
        print(f"gate {name} PASS ({R['gates'][name]['wall_s']}s)", flush=True)
    R['bank_receipts'] = os.path.join(
        CAMPAIGN or '(table bank, not shipped)', 'galois',
        'SV_COMPLETE_STATE.json')
    if receipt_path:
        json.dump(R, open(receipt_path, 'w'), indent=1)
        print(f"receipt -> {receipt_path}", flush=True)
    return R


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', nargs='?', default='selftest',
                    choices=['selftest', 'derive'])
    ap.add_argument('--ring', default=None, help='ring pkl path')
    ap.add_argument('--svalign', default=None, help='svalign bank dir')
    ap.add_argument('--dps', type=int, default=None,
                    help='numeric-check digits (selftest default 460; '
                         'derive default off)')
    ap.add_argument('--receipt', default=None, help='selftest receipt path')
    ap.add_argument('--comps', default=None,
                    help='derive: outer-first comps "5,7,7;3,9,7"')
    ap.add_argument('--out', default=None, help='derive: output json path')
    a = ap.parse_args()
    if a.ring:
        init_ring(a.ring)
    if a.cmd == 'selftest':
        selftest(a.receipt, a.dps or CHECK_DPS, a.svalign)
    elif a.cmd == 'derive':
        assert a.comps, 'derive needs --comps'
        comps = [tuple(int(x) for x in c.split(','))
                 for c in a.comps.split(';')]
        out = derive_formulas(comps, a.dps, a.svalign)
        for k, e in out.items():
            print(e['formula_campaign_conv'])
            if 'numeric_check' in e:
                print(f"  numeric: {e['numeric_check']['agree_digits']}d "
                      f"at {e['numeric_check']['dps']}d")
        if a.out:
            json.dump({'formulas': out}, open(a.out, 'w'), indent=1)
            print(f"-> {a.out}")


if __name__ == '__main__':
    main()

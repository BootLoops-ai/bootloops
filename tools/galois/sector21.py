#!/usr/bin/env python3
# galois member — weight-21 sector closure driver.
"""sector21.py — SECTOR EXTENSION of the core coaction engine: motivic Galois
derivations D_m (m odd, 3..19) on depth<=3 odd-index MZVs at odd weights
11..21, with the exact MZV relation ring (galois.mzvring, reference ring pkls)
as the reduction oracle in place of the w<=8 t0deep tables. sv projection
included.

ADDITIVE MODULE:
core.py / vspace.py are UNTOUCHED; this module
imports galois.core (gate/table route only) and galois.mzvring
(reduction oracle) — the same module objects, never copies.

DATA PLACEMENT: the reference sector tables live under the table-bank root
(GALOIS_CAMPAIGN_BANK env var; TABLE_BANK/SECTOR_TABLES.json); every entry
point accepts an explicit table/ring path. Nothing multi-MB is duplicated
into tools/.

CONVENTION LOCKS (each is
enforced by an assert or a gate below — re-run `selftest` after ANY upstream
convention change):
 1. MZV convention (outer-first; the code and table keys label this package's
    convention 'campaign', as opposed to 'brown' for Brown's convention):
    mzv(a1,..,ak) = sum_{m1>m2>..>mk>=1} prod mi^-ai,
    OUTER index first, a1>=2. Symbols ('z',n), ('m2',a,b), ('m3',a,b,c);
    ('pi2',) = pi^2; zeta(2k) -> rational*pi2^k at construction.
    form_oracle.mzv is pinned to the SAME outer-first convention.
 2. Word <-> composition: (s1..sk) <-> 0^{s1-1}1 ... 0^{sk-1}1 (0=x, 1=y);
    Remiddi-Vermaseren, leftmost letter outermost; zeta_sh(word) =
    outer-first mzv(word_comp(word)) for convergent words. [gate a + assert]
 3. I-convention/signs (galois.core, reused): H_u = (-1)^{#1s} I(0;rev u;z);
    I(0;B;1) = (-1)^{#1s(B)} zeta_sh(rev B); I(1;B;0) = (-1)^{|B|} I(0;rev B;1);
    I(x;B;x)=0; shuffle reg: zeta_sh(1-word)=0 AND zeta_sh((0))=0.
 4. Brown/svalign: Brown zeta is INNER-first — Brown-zeta(a,b,c) = outer-first
    mzv3(c,b,a), EXACT REVERSAL (parse_brown_name reverses). [gate d]
 5. pi_m (m odd 3..19) := coefficient of the single canonical monomial
    (('z',m),) in the ring canonical form. Coincides with the core engine's
    zpoly-generator pi_m for m<=7 (gate a); basis-DEPENDENT for m>=11
    (depth-3 survivors, e.g. m3(5,3,3) at w=11) — locked to the ring
    survivors basis. pi_m kills products automatically => D_m is a derivation.
 6. D_m = (id x pi_m) Delta (Goncharov), kept word = left factor; only single
    contiguous cuts of length m survive pi_m (equality with the core engine's
    all-subsequences Dm_iword confirmed by gate a).

CLI:  python3 sector21.py selftest [--ring PKL]   -> quick regression battery
          (gate a: w<=8 core-table agreement, 255 values + 765 D-images;
           gate c: literature derivation values; gate e-spot: 6 FORM numeric
           reductions at 200d) + optional --receipt PATH (writes the battery
           receipt JSON)
      python3 sector21.py gates|oracle           -> full U1 batteries (a-e)
      python3 sector21.py sector|sv [--tables J] -> (re)generate sector tables
"""
import json
import os
import re
import sys
from collections import Counter
from fractions import Fraction as Fr

TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if TOOLS not in sys.path:
    sys.path.insert(0, TOOLS)

from galois import mzvring as MZR         # the ring MEMBER, by identity
from galois.mzvring import (Fr, dadd, dsum, dscale, dmul, znorm, mweight)  # noqa

CAMPAIGN = os.environ.get('GALOIS_CAMPAIGN_BANK')  # reference-table root (not shipped; point at your own)
TABLE_BANK = os.path.join(CAMPAIGN, 'galois') if CAMPAIGN else None
SVALIGN = os.path.join(CAMPAIGN, 'svalign') if CAMPAIGN else None


def _bank(p, what):
    """fail-closed default-path guard: the reference tables do not
    ship with this repo."""
    if p is None:
        raise RuntimeError(
            f"{what}: no explicit path given and GALOIS_CAMPAIGN_BANK is not "
            "set. The SECTOR_TABLES/SV_FORMULAS reference tables "
            "do not ship with this repo — "
            "pass an explicit path (--tables/--svalign/load_tables(path)) or "
            "point GALOIS_CAMPAIGN_BANK at your own bank (subdirs galois/, "
            "svalign/).")
    return p

# ---------------------------------------------------------- ring (lazy bind)

_STATE = {'ring': None, 'surv': None, 'path': None}


def init_ring(path=None):
    """Bind the reduction oracle to a reference ring pkl.
    Default: mzvring.RING_BANK/ring21.pkl (covers all weights 3..21 incl.
    even 8..18 — every depth<=3 symbol this module emits; coverage is
    hard-asserted in sym_vec)."""
    if path is None:
        if not MZR.RING_BANK:
            raise RuntimeError(
                "init_ring: no explicit path and GALOIS_RING_BANK unset — the "
                "reference ring pkls are not shipped; build one with "
                "mzvring/build_ring.py or pass init_ring(path).")
        path = os.path.join(MZR.RING_BANK, 'ring21.pkl')
    R = MZR.load_ring(path)
    _STATE.update(ring=R, path=path,
                  surv={s for ss in R.survivors.values() for s in ss})
    _ZSH.clear()                          # zsh cache is ring-keyed
    return R


def _ring():
    if _STATE['ring'] is None:
        init_ring()
    return _STATE['ring']


def sym_vec(s):
    """canonical vector of a single symbol; hard-stop on coverage gaps."""
    R = _ring()
    if s[0] == 'z':
        if s[1] % 2 == 0:
            return znorm(s[1])
        return {(s,): Fr(1)}
    if s in R.canon:
        return R.canon[s]
    if s in _STATE['surv']:
        return {(s,): Fr(1)}
    raise KeyError(f"ring coverage gap: symbol {s} (weight {sum(s[1:])}) "
                   f"not canonicalized and not a survivor ({_STATE['path']})")


# ------------------------------------------------------------- words / comps

def comp_word(c):
    return MZR.comp_word(c)


def word_comp(w):
    return MZR.word_comp(w)


def nones(w):
    return sum(1 for a in w if a == 1)


def shuffle_words(u, v):
    return MZR.shuffle(tuple(u), tuple(v))


# --------------------------------------------- shuffle-regularized zeta_sh(u)
# zsh(u) -> canonical vector (outer-first) for ANY {0,1} H-word with <=3 ones.

_ZSH = {}


def zsh(u):
    u = tuple(u)
    if u in _ZSH:
        return _ZSH[u]
    if not u:
        v = {(): Fr(1)}
    elif u[0] == 1:
        # leading-1 regularization: 0 = zsh(1)*zsh(u[1:]) = sum over insertions
        v = _reg_insert(u, u[1:], 1)
    elif u[-1] == 0:
        # trailing-0 regularization: 0 = zsh(u[:-1])*zsh((0,))
        v = _reg_insert(u, u[:-1], 0)
    else:
        c = word_comp(u)
        assert c[0] >= 2, (u, c)          # convention lock 2 (convergence)
        if len(c) == 1:
            v = _ring().cvec(znorm(c[0]))
        elif len(c) == 2:
            v = sym_vec(('m2',) + c)
        elif len(c) == 3:
            v = sym_vec(('m3',) + c)
        else:
            raise ValueError(f"depth {len(c)} > 3: {c}")
    for m in v:
        assert mweight(m) == len(u), (u, m)   # weight-homogeneity lock
    _ZSH[u] = v
    return v


def _reg_insert(u, base, letter):
    """zsh(u) from 0 = sum_{w in (letter) sh base} mult(w) zsh(w)."""
    mult = Counter()
    for i in range(len(base) + 1):
        mult[base[:i] + (letter,) + base[i:]] += 1
    cu = mult.pop(u)
    acc = {}
    for w, cw in mult.items():
        acc = dsum(acc, dscale(zsh(w), Fr(cw)))
    return dscale(acc, Fr(-1, cu))


# --------------------------------------------------------------- arc values

def I01(B):
    """I(0;B;1) as canonical vector."""
    if not B:
        return {(): Fr(1)}
    return dscale(zsh(tuple(reversed(B))), Fr((-1) ** nones(B)))


def arc(x, B, y):
    """I(x;B;y), x,y in {0,1}."""
    if not B:
        return {(): Fr(1)}
    if x == y:
        return {}
    if x == 0:
        return I01(B)
    return dscale(I01(tuple(reversed(B))), Fr((-1) ** len(B)))


# ------------------------------------------------------------------- pi_m, D_m

def pi_m(vec, m):
    """coefficient of the single canonical monomial zeta(m), m odd >= 3
    (convention lock 5: ring-basis splitting)."""
    assert m % 2 == 1 and m >= 3, m
    return vec.get((('z', m),), Fr(0))


def Dm_word(u, m):
    """D_m zeta_sh(u) (H-word u, any regularizable word, <=3 ones)
    -> canonical vector at weight len(u)-m."""
    n = len(u)
    if n < m:
        return {}
    sign = Fr((-1) ** nones(u))
    A = tuple(reversed(u))
    out = {}
    for p in range(n - m + 1):
        x = A[p - 1] if p > 0 else 0
        y = A[p + m] if p + m < n else 1
        q = pi_m(arc(x, A[p:p + m], y), m)
        if q == 0:
            continue
        kept = A[:p] + A[p + m:]
        out = dsum(out, dscale(I01(kept), sign * q))
    for mono in out:
        assert mweight(mono) == n - m, (u, m, mono)
    return out


def Dm_comp(c, m):
    return Dm_word(comp_word(c), m)


# ------------------------------------------ D_m on canonical vectors (Leibniz)

def Dm_symbol(s, m):
    if s[0] == 'pi2':
        return {}
    if s[0] == 'z':
        return {(): Fr(1)} if s[1] == m else {}
    return Dm_comp(tuple(s[1:]), m)


def Dm_vec(vec, m):
    """derivation on a canonical vector (Leibniz on basis monomials)."""
    out = {}
    for mono, q in vec.items():
        for i, s in enumerate(mono):
            ds = Dm_symbol(s, m)
            if not ds:
                continue
            rest = {mono[:i] + mono[i + 1:]: q}
            out = dsum(out, dmul(rest, ds))
    return out


def sv_proj(vec):
    """sv kills f2: drop every monomial containing pi2."""
    return {mn: q for mn, q in vec.items()
            if not any(s[0] == 'pi2' for s in mn)}


# ------------------------------------------------------------- serialization
# atlas-style name form (byte-compatible with symbolic/close_exact.py).

def to_name_form(vec):
    terms = {}
    for mn, c in vec.items():
        j = sum(1 for s in mn if s[0] == 'pi2')
        parts = []
        if j:
            parts.append('z2' if j == 1 else f'z2^{j}')
            c = c * 6 ** j
        for s in sorted(s for s in mn if s[0] != 'pi2'):
            parts.append(f'z{s[1]}' if s[0] == 'z'
                         else f"mzv({','.join(map(str, s[1:]))})")
        terms['*'.join(parts) if parts else 'one'] = c
    if not terms:
        return {'den': 1, 'num': {}}
    from math import lcm, gcd
    den = 1
    for c in terms.values():
        den = lcm(den, c.denominator)
    nums = {n: int(c * den) for n, c in terms.items()}
    g = gcd(den, *[abs(v) for v in nums.values()])
    return {'den': den // g, 'num': {n: v // g for n, v in nums.items()}}


def veq(v1, v2):
    return dsum(v1, dscale(v2, Fr(-1))) == {}


# ------------------------------------------- core table-route conversion

def zpoly_to_campaign(zp):
    """t0deep zpoly {(ep,e3,e5,e7,e53): Fr} -> canonical vector (outer-first).
    p = i*pi:
    p^2 = -pi2. Asserts ep even (parity lock, gate d)."""
    out = {}
    for (ep, e3, e5, e7, e53), q in zp.items():
        assert ep % 2 == 0, f"parity lock violated: odd p-power {ep}"
        mono = ((('pi2',),) * (ep // 2) + (('z', 3),) * e3
                + (('z', 5),) * e5 + (('z', 7),) * e7
                + (('m2', 5, 3),) * e53)
        dadd(out, tuple(sorted(mono)), q * Fr((-1) ** (ep // 2)))
    return out


# ---------------------------------------------------------------- sector def

ODD_W = (11, 13, 15, 17, 19, 21)


def sector_words(w):
    """all-odd depth-3 compositions (a,b,c), a>=3 odd, b,c>=1 odd,
    a+b+c=w; plus the depth-1 (w,)."""
    out = [(w,)]
    for a in range(3, w - 1, 2):
        for b in range(1, w - a, 2):
            c = w - a - b
            if c >= 1 and c % 2 == 1:
                out.append((a, b, c))
    return out


def ms_for(w):
    return tuple(range(3, w - 1, 2))


def load_tables(path=None):
    """read the stored sector tables (plain 'D' + 'D_svproj' per word)."""
    return json.load(open(path or os.path.join(_bank(TABLE_BANK, 'load_tables'),
                                               'SECTOR_TABLES.json')))


# ------------------------------------------------------------------ gates a-e
# (full U1 battery; receipts land at <bank>/galois/GATES.json)

def gate_a():
    """this module vs the core t0deep table route: values + D3/D5/D7 for
    every w<=8 word with <=3 ones (255 words, 765 D-image checks)."""
    from galois import core as SG          # BY IDENTITY — untouched core
    t0 = SG.T0DEEP
    if t0 not in sys.path:
        sys.path.insert(0, t0)
    from mzv_formal import all_words
    words = [u for u in all_words(8) if nones(u) <= 3]
    nval = nimg = nimg_nz = 0
    for u in words:
        assert veq(zsh(u), zpoly_to_campaign(SG.ZSH[u])), f"value mismatch {u}"
        nval += 1
        for m in (3, 5, 7):
            ours = Dm_word(u, m)
            ref = zpoly_to_campaign(SG.Dm_word_value(u, m))
            assert veq(ours, ref), f"D{m} mismatch on {u}"
            nimg += 1
            if ours:
                nimg_nz += 1
    return {'pass': True, 'words': len(words), 'value_checks': nval,
            'D_checks': nimg, 'D_checks_nonzero': nimg_nz}


def gate_b(seed=20260809):
    """exact shuffle/stuffle/derivation identities in the new sector."""
    import random
    rng = random.Random(seed)
    checks = {'shuffle_value': 0, 'stuffle_value': 0, 'derivation_shuffle': 0}
    pairs = []
    for wt in (11, 13, 15, 17):
        for _ in range(3):
            a = rng.choice(range(2, wt - 3))
            u = comp_word((a,))
            rest = wt - a
            b = rng.choice(range(2, rest))
            cc = rest - b
            v = comp_word((b, cc) if cc >= 1 else (rest,))
            pairs.append((u, v, wt))
    for u, v, wt in pairs:
        lhs = dmul(zsh(u), zsh(v))
        rhs = {}
        for w, cnt in shuffle_words(u, v).items():
            rhs = dsum(rhs, dscale(zsh(w), Fr(cnt)))
        assert veq(lhs, rhs), f"shuffle FAIL {u} sh {v}"
        checks['shuffle_value'] += 1
    for (a, b, c) in [(3, 5, 3), (5, 5, 3), (3, 7, 3), (5, 7, 3), (3, 9, 3),
                      (7, 7, 3), (3, 11, 3), (5, 9, 5), (3, 13, 5), (9, 9, 3)]:
        lhs = dmul(_ring().cvec(znorm(a)), zsh(comp_word((b, c))))
        rhs = {}
        for t in [(a, b, c), (b, a, c), (b, c, a)]:
            rhs = dsum(rhs, zsh(comp_word(t)))
        rhs = dsum(rhs, zsh(comp_word((a + b, c))), zsh(comp_word((b, a + c))))
        assert veq(lhs, rhs), f"stuffle FAIL z{a}*mzv({b},{c})"
        checks['stuffle_value'] += 1
    trip = [((5, 3), (3,), 5), ((5, 3), (3,), 3), ((7, 3), (3,), 7),
            ((5, 3), (5,), 5), ((9, 3), (3,), 3), ((5, 3, 3), (3,), 5),
            ((7, 3, 3), (3,), 7), ((5, 5, 3), (5,), 9), ((3, 3), (5,), 11)]
    for c1, c2, m in trip:
        u, v = comp_word(c1), comp_word(c2)
        lhs = {}
        for w, cnt in shuffle_words(u, v).items():
            lhs = dsum(lhs, dscale(Dm_word(w, m), Fr(cnt)))
        rhs = dsum(dmul(Dm_word(u, m), zsh(v)), dmul(zsh(u), Dm_word(v, m)))
        assert veq(lhs, rhs), f"derivation FAIL D{m} on {c1} sh {c2}"
        checks['derivation_shuffle'] += 1
    checks['pass'] = True
    return checks


def gate_c():
    """literature-known derivation values (outer-first convention; citations)."""
    out = {'checks': []}

    def chk(name, got, want, cite):
        ok = veq(got, want)
        out['checks'].append({'name': name, 'pass': ok,
                              'got': to_name_form(got), 'cite': cite})
        assert ok, f"gate c FAIL {name}: {to_name_form(got)}"
    z3 = {(('z', 3),): Fr(1)}
    z5 = {(('z', 5),): Fr(1)}
    one = {(): Fr(1)}
    chk('D5 zeta(5,3) = -5 zeta3', Dm_comp((5, 3), 5), dscale(z3, Fr(-5)),
        'F. Brown arXiv:1102.1310 / 1102.1312 (zeta(3,5) ascending conv)')
    chk('D3 zeta(5,3) = 0', Dm_comp((5, 3), 3), {},
        'rho(zeta(3,5)) = -5 f3f5 + c f8, no word ending f3: Brown 1102.1310')
    chk('D7 zeta(7,3) = -14 zeta3', Dm_comp((7, 3), 7), dscale(z3, Fr(-14)),
        'rho(zeta(3,7)) = -14 f3f7 - 6 f5f5 + q10 f10, arXiv:2511.15883')
    chk('D5 zeta(7,3) = -6 zeta5', Dm_comp((7, 3), 5), dscale(z5, Fr(-6)),
        'same rho(zeta(3,7)) expansion, coeff of words ending f5')
    chk('D3 zeta(3) = 1', Dm_comp((3,), 3), one, 'normalization: pi_m(zm)=1')
    chk('D5 zeta(7) = 0 (odd zetas primitive)', Dm_comp((7,), 5), {},
        'zeta(odd) generators: D_m zeta(n)=0 for m<n, Brown 1102.1312')
    out['pass'] = True
    return out


def gate_d(svalign_dir=None):
    """parity / convention locks (G3 spirit)."""
    import mpmath as mp
    out = {}
    assert zsh((0, 0, 0, 0, 1, 0, 0, 1)) == {(('m2', 5, 3),): Fr(1)}
    out['word_convention'] = 'zsh(00001001) == m2(5,3) survivor exactly'
    from galois import core as SG
    for u, zp in SG.ZSH.items():
        zpoly_to_campaign(zp)   # raises on odd p-power
    out['parity_lock'] = 'all 510 t0deep zsh zpolys have even p-powers'
    sva = svalign_dir or _bank(SVALIGN, 'gate_d')
    cv = json.load(open(os.path.join(sva, 'conv_values_700.json')))['values']
    from formglue import form_oracle
    mp.mp.dps = 120
    locks = []
    for brown, camp in [('3,5', (5, 3)), ('3,7', (7, 3)),
                        ('3,3,5', (5, 3, 3)), ('3,5,7', (7, 5, 3))]:
        a = mp.mpf(cv[brown])
        b = form_oracle.to_mpf(form_oracle.mzv(camp, 100), 100)
        d = float(-mp.log10(abs(a - b) / abs(b)))
        locks.append({'zB': brown, 'campaign': list(camp),
                      'agree_digits': round(d, 1)})
        assert d > 90, (brown, camp, d)
    out['zB_reversal_lock'] = locks
    out['pass'] = True
    return out


SPOT6 = [(12, 3), (2, 11, 4), (9, 5, 3), (7, 7, 4), (11, 5, 2), (6, 5, 5)]


def gate_e_spot(picks=SPOT6, dps=200):
    """FORM numeric spot checks: LHS = form_oracle.mzv(comp) vs RHS = dps-digit
    numeric of the ring-canonical vector (survivor mzv's also via form_oracle,
    DIFFERENT indices — the relation itself is what's tested)."""
    import mpmath as mp
    from formglue import form_oracle
    mp.mp.dps = dps + 20
    cache = {}

    def val(s):
        if s not in cache:
            if s[0] == 'pi2':
                cache[s] = mp.pi ** 2
            elif s[0] == 'z':
                cache[s] = mp.zeta(s[1])
            else:
                cache[s] = form_oracle.to_mpf(
                    form_oracle.mzv(tuple(s[1:]), dps), dps)
        return cache[s]
    results = []
    for cmp_ in picks:
        vec = zsh(comp_word(cmp_))
        rhs = mp.mpf(0)
        for mono, q in vec.items():
            t = mp.mpf(q.numerator) / q.denominator
            for s in mono:
                t *= val(s)
            rhs += t
        lhs = form_oracle.to_mpf(form_oracle.mzv(cmp_, dps), dps)
        d = 999.0 if lhs == rhs else \
            float(-mp.log10(abs(lhs - rhs) / abs(lhs)))
        results.append({'mzv': list(cmp_), 'weight': sum(cmp_),
                        'terms': len(vec), 'agree_digits': round(d, 1)})
        assert d > dps - 25, (cmp_, d)
    return {'pass': True, 'dps': dps, 'n': len(results), 'samples': results}


# ------------------------------------------------------------------ sector run

def run_sector(out_path=None):
    out_path = out_path or os.path.join(_bank(TABLE_BANK, 'run_sector'), 'SECTOR_TABLES.json')
    doc = {'convention': 'outer-first mzv(a1,..)=sum_{m1>m2>..} prod mi^-ai; '
                         'basis = ring canonical (pi2-towers, z(odd), '
                         'survivor mzv monomials, products); '
                         'D_m = (id x pi_m)Goncharov, pi_m = coeff of zeta(m) '
                         'in ring basis; name form: z2 = zeta(2)',
           'sector': 'depth-3 all-odd-index words + depth-1 zeta(w), '
                     'odd w=11..21; atlas subset flag: all parts >= 3',
           'weights': {}}
    if os.path.exists(out_path):
        doc = json.load(open(out_path))
    for w in ODD_W:
        if str(w) in doc['weights']:
            continue
        wd = {}
        for c in sector_words(w):
            key = ','.join(map(str, c))
            entry = {'depth': len(c),
                     'atlas_sector': len(c) == 3 and min(c) >= 3, 'D': {},
                     'D_svproj': {}}
            for m in ms_for(w):
                img = Dm_comp(c, m)
                entry['D'][str(m)] = to_name_form(img)
                entry['D_svproj'][str(m)] = to_name_form(sv_proj(img))
            wd[key] = entry
        doc['weights'][str(w)] = wd
        json.dump(doc, open(out_path, 'w'), indent=1)
        print(f"w={w}: {len(wd)} words stored", flush=True)
    return doc


# ------------------------------------------------------- sv formula derivations

def parse_brown_name(name):
    """svalign coeff name (Brown order) -> canonical vector (outer-first;
    convention lock 4: EXACT REVERSAL of the index tuple)."""
    v = {(): Fr(1)}
    for tok in name.split('*'):
        mm = re.fullmatch(r'z(\d+)(?:\^(\d+))?', tok)
        if mm:
            n, p = int(mm.group(1)), int(mm.group(2) or 1)
            for _ in range(p):
                v = dmul(v, _ring().cvec(znorm(n)))
            continue
        mm = re.fullmatch(r'zB\(([\d, ]+)\)', tok)
        if mm:
            idx = tuple(int(x) for x in mm.group(1).split(','))[::-1]  # REVERSE
            s = (('m2',) if len(idx) == 2 else ('m3',)) + idx
            v = dmul(v, sym_vec(s))
            continue
        raise ValueError(f"unparsed sv token {tok!r}")
    return v


def run_sv(out_path=None, svalign_dir=None):
    out_path = out_path or os.path.join(_bank(TABLE_BANK, 'run_sv'), 'SECTOR_TABLES.json')
    doc = json.load(open(out_path))
    SVF = json.load(open(os.path.join(svalign_dir or _bank(SVALIGN, 'run_sv'),
                                      'SV_FORMULAS.json')))
    import ast
    done = 0
    for k, f in SVF['formulas'].items():
        brown = ast.literal_eval(k)
        camp = tuple(reversed(brown))
        w = f['weight']
        vec = dscale(sym_vec(('m3',) + camp), Fr(2))
        for name, q in f['coeffs_brown_basis'].items():
            vec = dsum(vec, dscale(parse_brown_name(name), Fr(q)))
        entry = {'brown_key': list(brown), 'status': f['status'],
                 'zeta_sv_campaign_reduced': to_name_form(vec),
                 'D': {}, 'D_svproj': {}}
        for m in ms_for(w):
            img = Dm_vec(vec, m)
            entry['D'][str(m)] = to_name_form(img)
            entry['D_svproj'][str(m)] = to_name_form(sv_proj(img))
        doc['weights'].setdefault(str(w), {})
        doc['weights'][str(w)].setdefault('_sv_formulas', {})[
            ','.join(map(str, camp))] = entry
        done += 1
    doc['sv_coverage'] = {'formulas': done,
                          'blocked': {k: v['weight']
                                      for k, v in SVF['blocked'].items()}}
    json.dump(doc, open(out_path, 'w'), indent=1)
    print(f"sv: {done} formulas stored", flush=True)


# ------------------------------------------------------------------------ main

def selftest(receipt_path=None):
    """quick regression battery (registration gate subset of U1 a-e)."""
    import datetime
    R = {'date': datetime.datetime.now().isoformat(timespec='seconds'),
         'module': os.path.abspath(__file__),
         'ring_pkl': None, 'gates': {}}
    _ring()
    R['ring_pkl'] = _STATE['path']
    R['gates']['a_w8_vs_core_tables'] = gate_a()
    print('gate a PASS', R['gates']['a_w8_vs_core_tables'], flush=True)
    R['gates']['c_literature_values'] = gate_c()
    print('gate c PASS', flush=True)
    R['gates']['e_spot_form_oracle'] = gate_e_spot()
    print('gate e-spot PASS', R['gates']['e_spot_form_oracle']['samples'],
          flush=True)
    R['full_battery'] = ((TABLE_BANK or '(table bank, not shipped)') + '/GATES.json') + \
        ' (U1 gates a-e + sv closure receipt)'
    if receipt_path:
        json.dump(R, open(receipt_path, 'w'), indent=1)
        print(f"receipt -> {receipt_path}", flush=True)
    return R


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', nargs='?', default='selftest',
                    choices=['selftest', 'gates', 'oracle', 'sector', 'sv'])
    ap.add_argument('--ring', default=None, help='ring pkl path')
    ap.add_argument('--tables', default=None, help='SECTOR_TABLES.json path')
    ap.add_argument('--svalign', default=None, help='svalign dir')
    ap.add_argument('--receipt', default=None, help='selftest receipt path')
    a = ap.parse_args()
    if a.ring:
        init_ring(a.ring)
    if a.cmd == 'selftest':
        selftest(a.receipt)
    elif a.cmd == 'gates':
        print('a', gate_a(), flush=True)
        print('b', gate_b(), flush=True)
        print('c PASS' if gate_c()['pass'] else 'c FAIL', flush=True)
        print('d', gate_d(a.svalign), flush=True)
    elif a.cmd == 'oracle':
        print('e-spot', gate_e_spot(), flush=True)
    elif a.cmd == 'sector':
        run_sector(a.tables)
    elif a.cmd == 'sv':
        run_sv(a.tables, a.svalign)


if __name__ == '__main__':
    main()

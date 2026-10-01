#!/usr/bin/env python3
"""pslq_gate.py — shared PSLQ closure harness (library + CLI).

Consolidates the proven ring15/consts protocol:
  * canonicalize(vec)        sign-of-first-nonzero + gcd division
                             (canonicalization removes sign/gcd false instability)
  * two_prec_stable(...)     identical CANONICALIZED vectors at two dps legs;
                             raises DigitsError if the target string carries fewer
                             digits than the high leg (which would silently
                             null the fit)
  * capacity(...)            refuses scans whose expected digit consumption
                             n*log10(height)+margin exceeds available digits
                             (the 'no relation = height failure' class)
  * graded_basket(...)       weight-graded combo generator
  * reverify(...)            independent high-dps residual check
  * controls(...)            mandatory positive control in the SAME pool/protocol;
                             scan refuses to report nulls unless controls pass
  * members_audit(...)       members-only PSLQ independence audit of the basket

PROTOCOL (mandatory):
  1. Run members_audit(members, names) BEFORE EVERY SCAN. A basket with an
     internal Q-relation is DEGENERATE: mp.pslq returns the lowest-height
     relation in the lattice — the internal one, with coefficient 0 on the
     target — and any NULL on that basket is VOID, not an exclusion. This
     false-NULL class has voided real closures (height 7.5e4, well inside
     the claimed exclusion range).
  2. r[0]==0 from mp.pslq is NOT a null. two_prec_stable now raises
     DegeneratePoolError naming the internal relation and its members;
     quotient the pool (drop one member in the relation's support) and rescan.
  3. Positive controls must live INSIDE the scan basket (target' = known combo
     of scan-basket members) with coefficients >= 10x BELOW maxcoeff.
     controls() now warns on controls planted at/near the cap (a measured
     false-alarm class). A control in a different, clean basket
     attests the protocol, not the scan.
  4. Dilog/trilog baskets are degeneracy minefields: Landen/Abel at rational
     args (1/4, 1/3, -1/3), Cl2/Im-Li duplicates at 6th roots of unity,
     Glaisher Sl_n(pi/3) reductions. Curate to a Q-basis FIRST.

Usage (library):
    from pslq_gate import two_prec_stable, controls, reverify
    pool = {'1': '1.000...', 'pi^2': '9.8696...'}        # >=200-digit strings
    vec = two_prec_stable(target_str, ['1','pi^2'], pool, dps_pair=(120,200),
                          maxcoeff=10**6)
Usage (CLI):
    pslq_gate.py --target out.json:'M0[1111]eps0' \
                 --pool pool.json --protocol default
    pslq_gate.py --selftest
pool.json = {"members": {name: value_string}, "weights": {name: int}(opt)}
"""
import argparse, hashlib, itertools, json, math, os, re, sys
import mpmath as mp

ARB = re.compile(r'\[?\s*([+\-]?[0-9][0-9.eE+\-]*)')
DEFAULT = {'dps_pair': (120, 200), 'maxcoeff': 10**6, 'maxsteps': 120000,
           'resid_slack': 7, 'reverify_dps': 210, 'reverify_tol': '1e-150',
           'capacity_margin': 20}

class DigitsError(ValueError): pass
class CapacityError(ValueError): pass

class DegeneratePoolError(ValueError):
    """mp.pslq returned a relation with target coefficient 0 => the basket
    itself carries an internal Q-relation. Every NULL on this basket is VOID
    (the internal relation shadows any target relation at any precision).
    Attributes: .relation {name: coeff}, .members [names], .residual (str)."""
    def __init__(self, relation, members, residual=None, dps=None):
        self.relation, self.members = dict(relation), list(members)
        self.residual, self.dps = residual, dps
        super().__init__(
            f'degenerate pool: internal relation {self.relation} among members '
            f'{self.members}' + (f' (resid {residual}' + (f' @ {dps}d)' if dps else ')')
                                 if residual else '') +
            ' — quotient the pool (drop one member in the relation support) and rescan;'
            ' NULLs on this basket are VOID')

def canonicalize(vec):
    """gcd-divide + flip so first nonzero entry is positive. None-safe."""
    if not vec: return None
    g = 0
    for v in vec: g = math.gcd(g, abs(int(v)))
    if g == 0: return None
    out = [int(v) // g for v in vec]
    for v in out:
        if v:
            if v < 0: out = [-q for q in out]
            break
    return tuple(out)

def ndigits(s):
    """Significant decimal digits carried by a numeric string."""
    s = str(s).strip().split('e')[0].split('E')[0]
    return len(s.replace('.', '').replace('-', '').replace('+', '').lstrip('0'))

def check_digits(target_str, dps_needed):
    nd = ndigits(target_str)
    if nd < dps_needed:
        raise DigitsError(f'target carries only {nd} digits < required {dps_needed} '
                          f'(silent-null guard); refusing to scan')
    return nd

def capacity(dps, n_members, height, margin=DEFAULT['capacity_margin']):
    """Expected digits consumed by an (n+1)-vector PSLQ at coeff height."""
    need = (n_members + 1) * math.log10(max(height, 2)) + margin
    return {'required': need, 'available': dps, 'ok': need <= dps}

def _pools_from_strings(members, names, dps_pair):
    pools = {}
    for leg in dps_pair:
        with mp.workdps(leg):
            pools[leg] = [mp.mpf(members[n]) for n in names]
    return pools

def two_prec_stable(target_str, names, members, dps_pair=DEFAULT['dps_pair'],
                    maxcoeff=DEFAULT['maxcoeff'], maxsteps=DEFAULT['maxsteps'],
                    resid_slack=DEFAULT['resid_slack'], margin=DEFAULT['capacity_margin']):
    """Canonicalized two-precision PSLQ. Returns canonical tuple or None (null).
    Raises DigitsError / CapacityError / DegeneratePoolError instead of
    silently nulling. Run members_audit(members, names) before scanning."""
    check_digits(target_str, max(dps_pair))
    for n in names: check_digits(members[n], max(dps_pair))
    cap = capacity(min(dps_pair), len(names), maxcoeff, margin)
    if not cap['ok']:
        raise CapacityError(f"capacity {cap['required']:.0f}d > available {cap['available']}d "
                            f"(n={len(names)}, height={maxcoeff:g}); shrink pool or raise dps")
    pools, vecs = _pools_from_strings(members, names, dps_pair), {}
    for leg in dps_pair:
        with mp.workdps(leg):
            x = [mp.mpf(target_str)] + pools[leg]
            r = mp.pslq(x, maxcoeff=maxcoeff, maxsteps=maxsteps)
            if not r: return None
            if r[0] == 0:   # internal relation among members shadows the target
                res = mp.nstr(abs(mp.fsum(c*v for c, v in zip(r, x))), 3)
                rel = {n: c for n, c in zip(names, r[1:]) if c}
                raise DegeneratePoolError(rel, [n for n in names if rel.get(n)],
                                          residual=res, dps=leg)
            res, scale = mp.fsum(c*v for c, v in zip(r, x)), max(abs(v) for v in x)
            if abs(res) > scale * mp.mpf(10)**(-(leg - 14 - resid_slack*len(x))): return None
            vecs[leg] = canonicalize(r)
    return vecs[dps_pair[0]] if vecs[dps_pair[0]] == vecs[dps_pair[1]] else None

def members_audit(members, names=None, dps=160, h=1e4, maxsteps=2000000,
                  resid_exp_frac=0.75):
    """Members-only independence audit — RUN BEFORE EVERY SCAN (<0.1 s typ.).
    PSLQs the basket alone (no target) at height h, quotients out each internal
    relation found (drops the last member in its support) and repeats.
    Returns [] iff the basket is Q-independent to height h at this dps;
    else a list of {'relation': {name: coeff}, 'residual': str, 'dps': int}.
    Any prior NULL on a basket with a nonempty audit is VOID (a basket can hide
    low-height internal relations, e.g. Abel–Landen relations among its members)."""
    names = list(names if names is not None else sorted(members))
    found = []
    with mp.workdps(dps):
        live = list(names)
        while len(live) >= 2:
            x = [mp.mpf(members[n]) for n in live]
            r = mp.pslq(x, maxcoeff=int(h), maxsteps=maxsteps)
            if not r: break
            res = abs(mp.fsum(c*v for c, v in zip(r, x)))
            scale = max(abs(v) for v in x)
            if res > scale * mp.mpf(10)**(-int(dps*resid_exp_frac)): break  # spurious
            rel = {n: c for n, c in zip(live, r) if c}
            found.append({'relation': rel, 'residual': mp.nstr(res, 3), 'dps': dps})
            drop = list(rel)[-1]                      # quotient: drop one member
            live = [n for n in live if n != drop]     # in the relation's support
    return found

def odd_zeta_pool(max_weight, dps=None):
    """Preset pool: all odd-zeta monomials ∏ζ(2k+1) of total weight ≤ max_weight.
    Returns (members: {name→value_str}, weights: {name→int}). This is the svMZV
    *product* subspace — CAVEAT: dim(svMZV_w) = #compositions of w into odd parts
    ≥3 (NOT #partitions); products alone don't span at w≥11 (depth-3 irreducibles
    enter)."""
    if dps: mp.mp.dps = dps
    members, weights = {'1': '1'}, {'1': 0}
    def parts(rem, mn):
        if rem == 0: yield (); return
        p = mn
        while p <= rem:
            for r in parts(rem-p, p): yield (p,)+r
            p += 2
    with mp.workdps(dps or mp.mp.dps):
        for w in range(3, max_weight+1):
            for P in parts(w, 3):
                name = '*'.join(f'z{p}' for p in P)
                members[name] = mp.nstr(mp.fprod(mp.zeta(p) for p in P),
                                        (dps or mp.mp.dps)-5, strip_zeros=False)
                weights[name] = w
    return members, weights

def graded_basket(weights, max_weight, sizes=(1, 2, 3)):
    """Yield name-tuples with total weight <= max_weight; skips all-weight-0 multis."""
    names = sorted(weights, key=lambda n: (weights[n], n))
    for size in sizes:
        for combo in itertools.combinations(names, size):
            if size > 1 and all(weights[n] == 0 for n in combo): continue
            if sum(weights[n] for n in combo) <= max_weight: yield combo

def reverify(hit, members, dps=DEFAULT['reverify_dps'], tol=DEFAULT['reverify_tol']):
    """Independent high-dps relative-residual check of {'target','members','vector'}."""
    with mp.workdps(dps):
        x = [mp.mpf(hit['target'])] + [mp.mpf(members[n]) for n in hit['members']]
        res = mp.fsum(c*v for c, v in zip(hit['vector'], x))
        rel = abs(res) / max(abs(v) for v in x)
        return {'relresid': mp.nstr(rel, 4), 'ok': rel < mp.mpf(tol)}

def _headroom(entry, vec, maxcoeff):
    """Planted-control coefficients must sit >=10x BELOW maxcoeff, else the
    control probes the cap, not the protocol (a false alarm otherwise)."""
    h = max(abs(int(c)) for c in vec)
    if 10 * h > maxcoeff:
        entry['headroom_warn'] = (f'control height {h} within 10x of maxcoeff '
                                  f'{maxcoeff:g}; plant controls >=10x below cap')
        sys.stderr.write(f"[pslq_gate] WARN {entry['control']}: "
                         f"{entry['headroom_warn']}\n")
    return entry

def controls(members, proto=DEFAULT, known_hits=()):
    """Positive controls in the SAME pool/protocol. Returns {'ok':bool,'detail':[...]}.
    NOTE: these attest ONLY the basket they run in — plant controls INSIDE the
    scan basket. Entries get 'headroom_warn' if coeffs sit within 10x of maxcoeff."""
    detail, names = [], sorted(members)[:2]
    hi = max(proto['dps_pair']) + 35
    with mp.workdps(hi):  # synthetic built ABOVE the high leg (the 235-dps fix)
        synth = (3*mp.mpf(members[names[0]]) - 7*mp.mpf(members[names[1]])) / 5
        sstr = mp.nstr(synth, hi - 7, strip_zeros=False)
    v = two_prec_stable(sstr, names, members, proto['dps_pair'], proto['maxcoeff'])
    detail.append(_headroom({'control': f'synthetic (3*{names[0]}-7*{names[1]})/5',
                             'got': v, 'want': (5, -3, 7), 'ok': v == (5, -3, 7)},
                            (5, -3, 7), proto['maxcoeff']))
    for kh in known_hits:
        got = two_prec_stable(kh['target'], kh['members'], members,
                              proto['dps_pair'], proto['maxcoeff'])
        detail.append(_headroom({'control': kh.get('name', 'known_hit'), 'got': got,
                                 'want': canonicalize(kh['vector']),
                                 'ok': got == canonicalize(kh['vector'])},
                                kh['vector'], proto['maxcoeff']))
    return {'ok': all(d['ok'] for d in detail), 'detail': detail}

def load_amflow_targets(fn, gamma_norm=True, build_dps=235, out_digits=228):
    """AMFlow out-json -> {'M<i>[<idx>]eps<n>': value-string}, e^{3*gamma*eps} norm."""
    d = json.load(open(fn))
    out = {}
    with mp.workdps(build_dps):
        g3 = [mp.mpf(3)**k * mp.euler**k / mp.factorial(k) for k in range(8)]
        for i, b in enumerate(d['result']):
            idx = ''.join(str(x) for x in b['integral']['indices'])
            cs = {c['order']: mp.mpf(ARB.match(str(c['value']['re'])).group(1))
                  for c in b['coefficients']}
            lo = min(cs)
            for n in sorted(cs):
                v = mp.fsum(cs[n-k]*g3[k] for k in range(n-lo+1)) if gamma_norm else cs[n]
                out[f'M{i}[{idx}]eps{n}'] = mp.nstr(v, out_digits, strip_zeros=False)
    return out

def load_target(spec):
    fn, key = spec.rsplit(':', 1)
    d = json.load(open(fn))
    return load_amflow_targets(fn)[key] if 'result' in d else str(d[key])

def _synthetic_amflow(path):
    """Write a synthetic AMFlow-format out-json with PLANTED integer relations.

    The raw per-order coefficients are built so that AFTER the e^{3*gamma*eps}
    normalization in load_amflow_targets the eps^0 values are exact rational
    combinations of the pool members:

      M0[1111]eps0 = (3 + 7*pi^2 - 5*zeta3)/11    -> vector (11, -3, -7, 5);
        split across orders -1 and 0 (cs[0] = t - 3*euler*cs[-1]), so a broken
        gamma normalization shifts the value and the planted refind FAILS —
        the normalization path is checked for real, not asserted.
      M1[2110]eps0 = (20011 + 911*pi^2)/73        -> vector (73, -20011, -911);
        single order (pass-through path).
    """
    with mp.workdps(250):
        pi2, z3, g = mp.pi**2, mp.zeta(3), mp.euler
        t_a = (3 + 7*pi2 - 5*z3) / 11
        c_m1 = (pi2 - 3) / 4                  # eps^-1 raw coefficient
        c_0 = t_a - 3*g*c_m1                  # eps^0 raw: normalization must undo this
        t_b = (mp.mpf(20011) + 911*pi2) / 73
        s = lambda v: mp.nstr(v, 240, strip_zeros=False)
        d = {'result': [
            {'integral': {'indices': [1, 1, 1, 1]},
             'coefficients': [{'order': -1, 'value': {'re': s(c_m1)}},
                              {'order': 0, 'value': {'re': s(c_0)}}]},
            {'integral': {'indices': [2, 1, 1, 0]},
             'coefficients': [{'order': 0, 'value': {'re': s(t_b)}}]}]}
    json.dump(d, open(path, 'w'))

def selftest():
    # member-integrity leg + synthetic planted-relation battery: the targets
    # are built on the fly in the AMFlow out-json format with known integer
    # relations (hand-provable answers), so the loader, the gamma
    # normalization, and the two-precision protocol are all checked for real.
    import tempfile
    mdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ring15')
    P, fails, total = dict(DEFAULT), 0, 0
    def rep(name, ok):
        nonlocal fails, total
        total += 1; fails += (not ok); print(('PASS' if ok else 'FAIL'), name, flush=True)
    # --- member-integrity leg: vendored shas vs the plain pin table ---
    pins = json.load(open(os.path.join(mdir, 'PINS.json')))
    sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
    rep(f"member-integrity: {len(pins)} vendored ring15 files match PINS.json shas",
        all(sha(os.path.join(mdir, fn)) == h for fn, h in pins.items()))
    with mp.workdps(235):
        members = {'1': '1.0' + '0'*226, 'pi^2': mp.nstr(mp.pi**2, 228, strip_zeros=False),
                   'zeta3': mp.nstr(mp.zeta(3), 228, strip_zeros=False)}
    fd, spath = tempfile.mkstemp(suffix='.json', prefix='pslq_gate_selftest_')
    os.close(fd)
    try:
        _synthetic_amflow(spath)
        targets = load_amflow_targets(spath)
    finally:
        os.unlink(spath)
    t_u0 = targets['M0[1111]eps0']
    v = two_prec_stable(t_u0, ['1', 'pi^2', 'zeta3'], members)
    rep('refind planted M0[1111]eps0 = [11,-3,-7,5] over [1,pi^2,zeta3] '
        '(gamma normalization exercised across two orders)',
        v == canonicalize([11, -3, -7, 5]))
    t_t4 = targets['M1[2110]eps0']
    v = two_prec_stable(t_t4, ['1', 'pi^2'], members)
    rep('refind planted M1[2110]eps0 = [73,-20011,-911] over [1,pi^2]',
        v == canonicalize([73, -20011, -911]))
    rep('canonicalize kills sign-flip instability ([9,-4,0] vs [-9,4,0])',
        canonicalize([9, -4, 0]) == canonicalize([-9, 4, 0]) == (9, -4, 0))
    try:
        two_prec_stable(t_u0, ['1'], members, dps_pair=(30, 200),
                        maxcoeff=10**8)  # 30d leg cannot hold (n+1)*8+20
        rep('capacity() rejects impossible scan (30d leg, height 1e8)', False)
    except CapacityError:
        rep('capacity() rejects impossible scan (30d leg, height 1e8)', True)
    rep('capacity() arithmetic (32 members @1e8 needs >30d)',
        not capacity(30, 32, 10**8)['ok'] and capacity(200, 3, 10**6)['ok'])
    try:
        two_prec_stable(t_u0[:120], ['1', 'pi^2'], members)  # ~115-digit string
        rep('digit-length guard raises on truncated (115d) target string', False)
    except DigitsError:
        rep('digit-length guard raises on truncated (115d) target string', True)
    c = controls(members, P, known_hits=[{'name': 'planted-hit-as-control', 'target': t_u0,
        'members': ['1', 'pi^2', 'zeta3'], 'vector': [11, -3, -7, 5]}])
    rep('controls() pass in same pool/protocol (synthetic + known hit)', c['ok'])
    r1 = reverify({'target': t_u0, 'members': ['1', 'pi^2', 'zeta3'],
                   'vector': [11, -3, -7, 5]}, members)
    r2 = reverify({'target': t_t4, 'members': ['1', 'pi^2'],
                   'vector': [73, -20011, -911]}, members)
    rep(f"reverify both hits @210d < 1e-150 ({r1['relresid']}, {r2['relresid']})",
        r1['ok'] and r2['ok'])
    rep('graded_basket grading (wmax=2 excludes zeta3-pairs of weight>2)',
        ('1', 'pi^2') in set(graded_basket({'1': 0, 'pi^2': 2, 'zeta3': 3}, 2))
        and all('zeta3' not in c for c in graded_basket({'1': 0, 'pi^2': 2, 'zeta3': 3}, 2)))
    # --- degenerate-pool cases ---
    with mp.workdps(235):
        dmem = {'pi^2': mp.pi**2, 'log2^2': mp.log(2)**2,
                'log2log3': mp.log(2)*mp.log(3), 'log3^2': mp.log(3)**2,
                'Li2(1/3)': mp.polylog(2, mp.mpf(1)/3),
                'Li2(1/4)': mp.polylog(2, mp.mpf(1)/4)}
        dmem = {k: mp.nstr(v, 228, strip_zeros=False) for k, v in dmem.items()}
    dnames = ['pi^2', 'log2^2', 'log2log3', 'log3^2', 'Li2(1/3)', 'Li2(1/4)']
    aud = members_audit(dmem, dnames)
    avec = canonicalize([aud[0]['relation'].get(n, 0) for n in dnames]) if aud else None
    rep('members_audit finds the height-12 MT internal relation '
        '(pi^2-12log2^2+12log2log3-6log3^2-12Li2(1/3)-6Li2(1/4)=0)',
        avec == (1, -12, 12, -6, -12, -6))
    rep('members_audit clean basket -> [] ([1,pi^2,zeta3])',
        members_audit(members, ['1', 'pi^2', 'zeta3']) == [])
    try:
        two_prec_stable(members['zeta3'], dnames, dmem)
        rep('degenerate basket raises DegeneratePoolError (not silent NULL)', False)
    except DegeneratePoolError as e:
        rep('degenerate basket raises DegeneratePoolError (not silent NULL)',
            set(e.members) == set(dnames) and e.relation.get('pi^2') in (1, -1))
    c2 = controls(members, {**P, 'maxcoeff': 10**5},
                  known_hits=[{'name': 'planted-at-cap', 'target': t_t4,
                               'members': ['1', 'pi^2'],
                               'vector': [73, -20011, -911]}])
    rep('controls() WARNS when control height sits within 10x of maxcoeff',
        c2['ok'] and any(d.get('headroom_warn') for d in c2['detail'])
        and not any(d.get('headroom_warn') for d in c['detail']))
    # --- graded-peel leg (mplll_graded.py sign convention): after a
    #     block locks, the sequential peel must strictly REDUCE the residual
    #     (the inverted '-=' DOUBLED it) so the next weight block locks.
    #     Planted (5,-3,7,0): weight-2 block {pi^2 x, log^2(2) x^2} carries the
    #     whole target, weight-1 block {log(3) x^3} must lock the peeled ~0
    #     residual — exactly the class the July kite receipt exposed. ---
    from . import mplll, mplll_graded
    gd = 60
    with mp.workdps(gd + 40):
        xs = [mp.mpf(p) / q for p, q in
              ((1, 7), (1, 5), (1, 3), (2, 5), (1, 2), (3, 5), (2, 3))]
        gfa = [mp.pi ** 2 * x for x in xs]
        gfb = [mp.log(2) ** 2 * x ** 2 for x in xs]
        gfc = [mp.log(3) * x ** 3 for x in xs]
        gI = [(3 * a - 7 * b) / 5 for a, b in zip(gfa, gfb)]
        gs = lambda v: [mp.nstr(t, gd + 30) for t in v]
        gI, gF = gs(gI), [gs(gfa), gs(gfb), gs(gfc)]
    gr = mplll_graded.graded_fit(gI, gF, ['pi2x', 'log2sq_x2', 'log3_x3'],
                                 [2, 2, 1], [0, 1, 2, 3, 4], [5, 6], d=gd)
    grel = canonicalize(gr['relation']) if gr['relation'] else None
    rep('graded peel: planted (5,-3,7,0) — w2 locks, peel strictly reduces, '
        'w1 locks the peeled residual (HIT, all_locked, w1 insample>40d)',
        gr['status'] == 'HIT' and gr['all_locked'] and grel == (5, -3, 7, 0)
        and gr['per_weight'][1]['locked']
        and (gr['per_weight'][1]['insample_d'] or 0) > 40)
    w2rel = gr['per_weight'][2]['relation']
    with mp.workdps(gd + 40):
        peeled = [abs(mplll._mpc(gI[j]) + mp.fsum(
                      mp.mpf(w2rel[k + 1]) * mplll._mpc(gF[k][j]) for k in range(2))
                      / mp.mpf(w2rel[0])) for j in range(len(xs))]
        gmax = max(abs(mplll._mpc(v)) for v in gI)
        gratio = float(max(peeled) / gmax)
    rep(f'graded peel direction matches heldout_digits semantics '
        f'(|r + sum(rel*F)/c0| / |I| = {gratio:.1e} < 1e-30)', gratio < 1e-30)
    print('SELFTEST', 'FAIL' if fails else 'PASS', f'({total-fails}/{total})', flush=True)
    return 1 if fails else 0

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--target', help='file.json:key (AMFlow out-json or flat dict)')
    ap.add_argument('--pool', help='pool json: {"members":{...},"weights":{...}}')
    ap.add_argument('--protocol', default='default')
    ap.add_argument('--max-weight', type=int, default=99)
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest: sys.exit(selftest())
    if not (a.target and a.pool): ap.error('--target and --pool required (or --selftest)')
    P, ts = DEFAULT, load_target(a.target)
    pj = json.load(open(a.pool))
    members = pj.get('members', pj)
    weights = pj.get('weights', {n: 1 for n in members})
    aud = members_audit(members, sorted(members))   # mandatory pre-scan audit
    if aud:
        print(json.dumps({'status': 'DEGENERATE_POOL', 'internal_relations': aud,
                          'note': 'quotient the pool and rerun; NULLs would be VOID'},
                         default=str, indent=1)); sys.exit(3)
    c = controls(members, P)
    if not c['ok']:
        print(json.dumps({'status': 'CONTROLS_FAILED', 'detail': c['detail']},
                         default=str, indent=1)); sys.exit(2)
    try:
        for combo in graded_basket(weights, a.max_weight):
            v = two_prec_stable(ts, list(combo), members, P['dps_pair'], P['maxcoeff'])
            if v:
                hit = {'target': ts, 'members': list(combo), 'vector': list(v)}
                hit['reverify'] = reverify(hit, members)
                print(json.dumps({'status': 'HIT', **{k: hit[k] for k in
                      ('members', 'vector', 'reverify')}}, default=str)); sys.exit(0)
    except DegeneratePoolError as e:   # height above audit's h — still not a NULL
        print(json.dumps({'status': 'DEGENERATE_POOL', 'internal_relation': e.relation,
                          'members': e.members, 'residual': e.residual},
                         default=str, indent=1)); sys.exit(3)
    print(json.dumps({'status': 'NULL', 'controls_ok': True, 'members_audit': 'clean'}))
    sys.exit(1)

if __name__ == '__main__':
    main()

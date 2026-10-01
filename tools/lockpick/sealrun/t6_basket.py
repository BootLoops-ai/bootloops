#!/usr/bin/env python3
# lockpick sealrun member — frozen-basket builder (double-routed members, members_audit, planted/negative gates; quotient BEFORE freezing — the DegeneratePoolError doctrine).
"""t6_basket.py — build the frozen 24-member basket consumed by t6_close.py.

Basket: weight-3 CM-newform (f4) periods/quasiperiods and critical L-values
(consumed from sha-pinned reference receipts), Dirichlet L-values, zeta values,
pi powers, log 2 — 24 members.

Discipline (the lockpick pslq_gate protocol):
  * every house-computed member double-routed (route A / route B agreement floored),
  * reference members consumed from sha-pinned receipts at their certified digit counts,
  * members_audit BEFORE any scan (degenerate-basket NULL is VOID),
  * planted positive control INSIDE the basket at >=10x headroom below maxcoeff,
  * matched-magnitude negative control (must NULL),
  * two-precision canonicalized scans; >=30 held-out digits at reverify.

Writes: BASKET_T6.json  {members:{name:str}, weights:{name:int}, provenance, audits}
FE trap honored: L(f4,3) is NEVER a member together with pi^2*L(f4,1)-class members
in a way that plants the functional-equation relation; we carry L1, L2 and
pi^2*Lf4_1 (=4*L3 exactly, checked as an internal gate, not a member pair).
Actually: since pi^2*Lf4_1 IS 4*L(f4,3), we include it as the L3-slot member and
never also include L3 itself.
"""
import json, hashlib, os, sys, time
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))  # package parent
from lockpick import pslq_gate  # noqa: E402

# Basket REBUILD input dir (the sha-pinned reference L-value/eta receipts).
# No default is shipped: the battery consumes the frozen BASKET_T6.json copy
# and does not rebuild; a rebuild requires the reference receipts.
R1 = os.environ.get('SEALRUN_R1')
PIN = {
    'LVALUES_RESULT.json': '2f0a94f1cb9f8df9b43fdb3fbadebded8b53b8551bc140dff5c0928711101c4e',
    'ETA_RESULT.json': '1beeed0b569a0b5037aeb5c185940a2589919c10a2eee71c25224d916dc6fbbd',
}

BUILD_DPS = 470       # working precision for house-computed members
SER_DIGITS = 452      # serialized digits for house members (< BUILD_DPS - guard)
CHECK_DPS_B = 60      # route-B independent spot-check floor (direct summation legs)


def sha256(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def floored_agree(a, b, dps):
    """floored agreed digits between two mpf at working dps."""
    with mp.workdps(dps):
        a, b = mp.mpf(a), mp.mpf(b)
        if a == b:
            return dps
        return int(mp.floor(-mp.log10(abs(a - b) / max(abs(a), abs(b)))))


def hurwitz_L(chi_vals, q, s):
    """L(chi, s) = q^-s sum_a chi(a) zeta(s, a/q); chi_vals[a] for a=1..q."""
    return mp.mpf(q) ** (-s) * mp.fsum(
        chi_vals[a - 1] * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q + 1))


def direct_L(chi_vals, q, s):
    """Route B (independent of the Hurwitz-zeta engine): accelerated block sum —
    g(m) = sum_a chi(a)/(qm+a)^s summed by mp.nsum's Richardson extrapolation.
    Runs at ambient workdps; measured full-window agreement at 160d."""
    g = lambda m: mp.fsum(chi_vals[a - 1] / (q * m + a) ** s
                          for a in range(1, q + 1) if chi_vals[a - 1])
    return mp.nsum(g, [0, mp.inf])


def main():
    t0 = time.time()
    if not R1:
        print('t6_basket: set SEALRUN_R1 to the directory holding the sha-pinned '
              'rebuild inputs (LVALUES_RESULT.json, ETA_RESULT.json). The battery '
              'consumes the frozen BASKET_T6.json copy and does not rebuild.')
        return 4
    out = {'producer': 't6_basket.py (lockpick sealrun member)',
           'stamp_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
           'build_dps': BUILD_DPS, 'serialized_digits': SER_DIGITS,
           'inputs': {}, 'route_checks': {}, 'members': {}, 'weights': {},
           'member_digits': {}}

    # --- pinned reference inputs ----------------------------------------------
    lv_p = os.path.join(R1, 'LVALUES_RESULT.json')
    et_p = os.path.join(R1, 'ETA_RESULT.json')
    out['inputs'][lv_p] = sha256(lv_p)
    out['inputs'][et_p] = sha256(et_p)
    assert out['inputs'][lv_p] == PIN['LVALUES_RESULT.json'], 'LVALUES sha drift'
    assert out['inputs'][et_p] == PIN['ETA_RESULT.json'], 'ETA sha drift'
    lv = json.load(open(lv_p))
    et = json.load(open(et_p))
    L1s = lv['values_430_digits']['L1']
    L2s = lv['values_430_digits']['L2']
    L3s = lv['values_430_digits']['L3']
    omp = et['values']['omega_plus']            # certified 461
    omm = et['values']['omega_minus_over_i']    # certified 460
    etm = et['values']['eta_minus_over_i']      # certified 460
    etp = et['values']['eta_plus']              # certified 459

    M, W, D = {}, {}, {}

    def put(name, val, w, digits):
        with mp.workdps(BUILD_DPS):
            M[name] = mp.nstr(mp.mpf(val), digits, strip_zeros=False)
        W[name] = w
        D[name] = digits

    with mp.workdps(BUILD_DPS):
        pi = mp.pi
        l2 = mp.log(2)
        z3 = mp.zeta(3)
        z5 = mp.zeta(5)
        cat = mp.catalan

        # route-B doubles for house constants (independent algorithms)
        rc = out['route_checks']
        rc['pi_vs_gamma_half_sq'] = floored_agree(pi, mp.gamma(mp.mpf(1) / 2) ** 2, BUILD_DPS)
        rc['log2_vs_2atanh_third'] = floored_agree(l2, 2 * mp.atanh(mp.mpf(1) / 3), BUILD_DPS)
        rc['zeta3_vs_altzeta'] = floored_agree(z3, mp.altzeta(3) / (1 - mp.mpf(2) ** (-2)), BUILD_DPS)
        rc['zeta5_vs_altzeta'] = floored_agree(z5, mp.altzeta(5) / (1 - mp.mpf(2) ** (-4)), BUILD_DPS)
        rc['catalan_vs_ImLi2_i'] = floored_agree(cat, mp.polylog(2, mp.mpc(0, 1)).imag, BUILD_DPS)

        # Dirichlet L-values (route A: Hurwitz; route B: 60d direct block sums)
        chi_m3 = [1, -1, 0]                                   # mod 3, odd
        chi_m4 = [1, 0, -1, 0]                                # mod 4, odd (Catalan control)
        chi_p5 = [1, -1, -1, 1, 0]                            # mod 5, even
        Lm3 = hurwitz_L(chi_m3, 3, 2)
        Lp5 = hurwitz_L(chi_p5, 5, 2)

    # exact chi_{-20}(a) = kronecker(-20, a) = (-1/a)*(5/a) for gcd(a,20)=1,
    # with 5 = 1 mod 4 so (5/a) = (a/5) by reciprocity; squares mod 5 = {1,4}.
    def chi_m20_val(a):
        if a % 2 == 0 or a % 5 == 0:
            return 0
        s = 1 if a % 4 == 1 else -1
        return s * (1 if a % 5 in (1, 4) else -1)

    chi_m20 = [chi_m20_val(a) for a in range(1, 21)]
    with mp.workdps(BUILD_DPS):
        Lm20 = hurwitz_L(chi_m20, 20, 2)

        # route B: independent accelerated block sums at 160d (measured full-window)
        with mp.workdps(160):
            rc2 = out['route_checks']
            rc2['L_m3_2_richardson160'] = floored_agree(
                mp.mpf(mp.nstr(Lm3, 158)), direct_L(chi_m3, 3, 2), 160)
            rc2['L_p5_2_richardson160'] = floored_agree(
                mp.mpf(mp.nstr(Lp5, 158)), direct_L(chi_p5, 5, 2), 160)
            rc2['L_m20_2_richardson160'] = floored_agree(
                mp.mpf(mp.nstr(Lm20, 158)), direct_L(chi_m20, 20, 2), 160)
            rc2['catalan_vs_richardson160'] = floored_agree(
                mp.mpf(mp.nstr(cat, 158)), direct_L(chi_m4, 4, 2), 160)

        # ---------------- the 24 members ----------------
        put('one', 1, 0, SER_DIGITS)
        put('pi', pi, 1, SER_DIGITS)
        put('pi^2', pi ** 2, 2, SER_DIGITS)
        put('pi^3', pi ** 3, 3, SER_DIGITS)
        put('pi^4', pi ** 4, 4, SER_DIGITS)
        put('log2', l2, 1, SER_DIGITS)
        put('log2^2', l2 ** 2, 2, SER_DIGITS)
        put('log2^3', l2 ** 3, 3, SER_DIGITS)
        put('pi*log2', pi * l2, 2, SER_DIGITS)
        put('pi^2*log2', pi ** 2 * l2, 3, SER_DIGITS)
        put('pi*log2^2', pi * l2 ** 2, 3, SER_DIGITS)
        put('zeta3', z3, 3, SER_DIGITS)
        put('pi*zeta3', pi * z3, 4, SER_DIGITS)
        put('zeta5', z5, 5, SER_DIGITS)
        put('catalan', cat, 2, SER_DIGITS)
        put('L_m3_2', Lm3, 2, SER_DIGITS)
        put('L_p5_2', Lp5, 2, SER_DIGITS)
        put('L_m20_2', Lm20, 2, SER_DIGITS)
        # stored f4 critical L-values (430d strings; two-route gated)
        M['Lf4_1'] = L1s; W['Lf4_1'] = 1; D['Lf4_1'] = 430
        M['Lf4_2'] = L2s; W['Lf4_2'] = 2; D['Lf4_2'] = 430
        # stored eichler f4 periods / quasiperiods (real realizations).
        # PERIOD/L-VALUE ENTANGLEMENT (measured):
        #   omega_plus       = 2 * pi^2 * L(f4,1)   (0 residual at 420d)
        #   omega_minus/i    = 4 * pi  * L(f4,2)    (0 residual at 420d)
        # so the periods ARE the pi^k*L(f4,s) mixing directions; the members
        # 'pi^2*Lf4_1' / 'pi*Lf4_2' / 'Lf4_3' must NEVER co-enter this basket
        # (guaranteed DegeneratePoolError, and FE makes L3 = pi^2 L1/4).
        M['omega_plus'] = omp; W['omega_plus'] = 3; D['omega_plus'] = 461
        M['omega_minus/i'] = omm; W['omega_minus/i'] = 3; D['omega_minus/i'] = 460
        M['eta_plus'] = etp; W['eta_plus'] = 3; D['eta_plus'] = 459
        M['eta_minus/i'] = etm; W['eta_minus/i'] = 3; D['eta_minus/i'] = 460

        # cross-gates: FE (L3 = pi^2 L1 / 4) + the measured period decompositions
        fe = floored_agree(pi ** 2 * mp.mpf(L1s) / 4, mp.mpf(L3s), 425)
        out['route_checks']['FE_gate_pi2L1_over4_vs_L3'] = fe
        assert fe >= 400, 'FE gate failed'
        out['route_checks']['omega_plus_vs_2pi2L1'] = floored_agree(
            mp.mpf(omp), 2 * pi ** 2 * mp.mpf(L1s), 420)
        out['route_checks']['omega_minus/i_vs_4piL2'] = floored_agree(
            mp.mpf(omm), 4 * pi * mp.mpf(L2s), 420)
        assert out['route_checks']['omega_plus_vs_2pi2L1'] >= 400
        assert out['route_checks']['omega_minus/i_vs_4piL2'] >= 400

    assert len(M) == 24, f'basket has {len(M)} members, want 24'
    out['members'], out['weights'], out['member_digits'] = M, W, D

    # --- members_audit (mandatory; degenerate-basket NULL is VOID) -----------
    aud = pslq_gate.members_audit(M, sorted(M), dps=300, h=1e5)
    out['members_audit'] = {'h': 1e5, 'dps': 300, 'found': aud}
    if aud:
        print('AUDIT FOUND INTERNAL RELATIONS:', json.dumps(aud, default=str))

    # --- planted positive control at 10x headroom ----------------------------
    # target' = (3*pi^2 - 7*catalan + 2*log2^3)/5  -> canonical (5,-3,7,-2)
    names_c = ['pi^2', 'catalan', 'log2^3']
    with mp.workdps(400):
        synth = (3 * mp.mpf(M['pi^2']) - 7 * mp.mpf(M['catalan'])
                 + 2 * mp.mpf(M['log2^3'])) / 5
        sstr = mp.nstr(synth, 390, strip_zeros=False)
    v = pslq_gate.two_prec_stable(sstr, sorted(M), M, dps_pair=(220, 265),
                                  maxcoeff=10 ** 6)
    want = {n: 0 for n in sorted(M)}
    want.update({'pi^2': -3, 'catalan': 7, 'log2^3': -2})
    wantvec = pslq_gate.canonicalize([5] + [want[n] for n in sorted(M)])
    out['planted_control'] = {'names': names_c, 'got': list(v) if v else None,
                              'want': list(wantvec),
                              'ok': v == wantvec,
                              'headroom_x': 10 ** 6 // 7}
    # --- negative control: matched-magnitude pseudo-random real -> must NULL -
    with mp.workdps(400):
        neg = mp.mpf(
            '0.' + hashlib.sha512(b'T6-negative-control-2026-08-29').hexdigest()
            .translate(str.maketrans('abcdef', '123456'))
            + hashlib.sha512(b'T6-negative-control-2026-08-29-b').hexdigest()
            .translate(str.maketrans('abcdef', '123456'))
            + hashlib.sha512(b'T6-negative-control-2026-08-29-c').hexdigest()
            .translate(str.maketrans('abcdef', '123456')))
        negs = mp.nstr(neg * mp.mpf(M['pi^2']), 390, strip_zeros=False)
    try:
        nv = pslq_gate.two_prec_stable(negs, sorted(M), M, dps_pair=(220, 265),
                                       maxcoeff=10 ** 6)
    except pslq_gate.DegeneratePoolError as e:
        nv = ('DEGENERATE', str(e))
    out['negative_control'] = {'got': None if nv is None else str(nv),
                               'ok': nv is None}

    out['wall_s'] = round(time.time() - t0, 2)
    json.dump(out, open(os.path.join(HERE, 'BASKET_T6.json'), 'w'),
              indent=1, default=str)
    ok = (not aud) and out['planted_control']['ok'] and out['negative_control']['ok']
    print('BASKET BUILD', 'PASS' if ok else 'FAIL',
          f"(audit {'clean' if not aud else 'RELATIONS'}, planted "
          f"{out['planted_control']['ok']}, negative {out['negative_control']['ok']}, "
          f"{out['wall_s']}s)")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

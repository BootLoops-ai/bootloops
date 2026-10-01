#!/usr/bin/env python3
"""run_cal_b.py — CAL-B (license gate 1): the rank-6 Sym^5-Legendre control
through the FULL matrix transport.  ALL SIX Frobenius eigenvalues (chi-twisted
alpha^i beta^{5-i}, the CAL-A-measured convention) must match EXACTLY at >= 3
ordinary fibers, p in {23, 29}, at a decisive certified depth.  Operator:
work/SYM5_OP.json (built + verified by build_sym5.py).  Twist law: (1-x) at
multiplicity 5 (the measured c6 = -64(1-x)^5).
Emits work/CAL_B_RECEIPT.json."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deform_kit as dk

t0 = time.time()
rec = {'producer': dk.producer('code/run_cal_b.py'),
       'spec': 'CAL-B: rank-6 exact-eigenvalue control; convention from work/CAL_A_RECEIPT.json',
       'sym5_operator': 'work/SYM5_OP.json (annihilates f0^5 to depth 418; theta-lead -64(1-x)^5)',
       'register': 'license gate 1; every verdict script-emitted',
       'primes': {}}
cala = json.load(open(os.path.join(dk.MEMBER_DIR, 'work', 'CAL_A_RECEIPT.json')))
rec['convention_consumed'] = cala['convention_measured']
assert cala['verdict'].startswith('CAL-A-PASS'), 'CAL-A must pass first'
op = dk.op_sym5()
K = 8
# The bar: ALL SIX eigenvalues match EXACTLY at the measured precision.
# MIN_CERT = 3: six coefficient-integers compared mod p^3 (chance collision
# < p^{-18} per fiber; printed).  The p=29 depth-3 phenomenon (d_star
# insensitive to twist power, measured at s=16 and s=32) is recorded as an
# instrument note; the p=23 column reaches cert >= 6 on the same recipe.
MIN_CERT = 3
overall = True
for p in (23, 29):
    row = {'K_target': K, 'min_cert_required': MIN_CERT}
    try:
        tr = dk.transport(op, p, K)
        r = tr['res']
        row['transport'] = {k: r[k] for k in ('p', 'op', 'NF', 'Ntot', 'jordan',
                                              'nilpotent', 'E_space_dim',
                                              'inv_residual_certified_depth')}
        row['march'] = r['march']
        row['pin'] = r['pin']
        row['twist_law'] = r['twist']
        E0, evec = dk.pinned_E(tr)
        F = dk.assemble_full_F(tr, E0)
        E0, F, baseline = dk.normalize_ray(tr, E0, F)
        row['ray_baseline'] = baseline
        gd = dk.gd_residual(tr, F)
        row['GD_residual_worst_deficit'] = gd
        row['GD_pass'] = bool((gd or 0) <= tr['keep'] - tr['keep_target'] + 64)
        gb = dk.gamma_branches(E0, p, wt=5, n=6, Kg=max(6, min((tr['d_star'] or 6) + 2, 16)))
        row['gamma'] = {'val_detE0': gb['val_detE0'], 'target_val': gb['target_val'],
                        'gamma_valuation': gb.get('gamma_valuation'),
                        'n_branches': len(gb['branches'])}
        chi = dk.f2.chi_table(p)
        chim1 = 1 if p % 4 == 1 else -1
        fibers = {}
        n_exact = 0
        for lam in (4, 9, 16, 6, 10):
            lam %= p
            if lam in (0, 1):
                continue
            ap = dk.f2.legendre_ap_fast(lam, p, chi)
            if ap % p == 0:
                continue
            if len(fibers) >= 4:
                break
            Fm, cert, tailv = dk.eval_F_at(tr, F, lam)
            mod = p ** cert
            cp = dk.berkowitz_charpoly(Fm, mod)
            c = dk.rp_coeffs(cp, p, cert)
            a = ap % mod
            for _ in range(cert.bit_length() + 4):
                a = (a - (a * a - ap * a + p) * pow((2 * a - ap) % mod, -1, mod)) % mod
            alpha, beta = a, (ap - a) % mod
            poly = [1]
            for i in range(6):
                lamb = chim1 * pow(alpha, i, mod) * pow(beta, 5 - i, mod) % mod
                new = [0] * (len(poly) + 1)
                for k2, ck in enumerate(poly):
                    new[k2] = (new[k2] + ck) % mod
                    new[k2 + 1] = (new[k2 + 1] - lamb * ck) % mod
                poly = new
            frow = {'lam': lam, 'ap_exact': ap, 'cert_depth': cert, 'tail_val': tailv,
                    'target': 'prod(1 - chi2(-1) alpha^i beta^(5-i) T), exact Hensel'}
            best = None
            for br in gb['branches']:
                g, v = br['gamma_unit'], br['gamma_val']
                ginv = pow(g, -1, mod)
                okall, depth = True, cert
                for k2 in range(1, 7):
                    t_ = poly[k2] * pow(ginv, k2, mod) % mod
                    if v >= 0:
                        pv = p ** (v * k2)
                        if t_ % pv:
                            okall = False
                            break
                        t_ //= pv
                    else:
                        t_ = t_ * pow(p, -v * k2, mod) % mod
                    diff = (c[k2] - t_) % mod
                    if diff:
                        okall = False
                        depth = min(depth, dk.valp(diff, p) or 0)
                if best is None or (okall, depth) > (best[2], best[3]):
                    best = (br['sign'], br['root_mod_p'], okall, depth)
            frow['best_branch'] = {'sign': best[0], 'root_mod_p': best[1],
                                   'ALL_SIX_EXACT_at_cert': bool(best[2]),
                                   'agreement_depth': best[3]}
            frow['decisive'] = bool(best[2] and cert >= MIN_CERT)
            frow['chance_collision_log10'] = round(-6 * cert * __import__('math').log10(p), 1)
            n_exact += int(frow['decisive'])
            fibers[str(lam)] = frow
        row['fibers'] = fibers
        row['n_decisive_exact_fibers'] = n_exact
        row['PASS'] = bool(n_exact >= 3 and row['GD_pass'] and
                           r['pin']['ray_unique_up_to_scalar'])
    except AssertionError as ex:
        row['WALL'] = str(ex.args)[:400]
        row['PASS'] = False
    overall = overall and row['PASS']
    rec['primes'][str(p)] = row
    print('p=%d PASS=%s (%.1fs)' % (p, row['PASS'], time.time() - t0), flush=True)

rec['verdict'] = ('CAL-B-PASS: ALL SIX Frobenius eigenvalues of the rank-6 Sym^5-Legendre control '
                  'reproduce EXACTLY through the matrix transport at >=3 ordinary fibers at each of '
                  'p in {23,29}, at certified depth >= %d — the control that refutes every '
                  'trace-grade route PASSES through this instrument' % MIN_CERT) if overall else \
                 'CAL-B-FAIL: see prime rows (K1 class: one debug cycle governs; license NOT issuable)'
rec['wall_s'] = round(time.time() - t0, 1)
dk.emit(rec, 'work/CAL_B_RECEIPT.json')
print('VERDICT:', rec['verdict'][:130])

#!/usr/bin/env python3
"""classify_mod.py — multi-eps classification of the exponent spectrum.

Port of analyze_frobenius.py + exclude_rational.py:
 * pair eigenvalue lists across eps by exhaustive 2-pt line fits (first vs
   last eps), each candidate line checked against ALL middle eps lists
   (the 4-eps verification pattern);
 * rationalize (a,b) of lambda = a + b*eps (limit_denominator);
 * classes: INT_EPS_INDEP (killed: vacuum theorem, when enabled),
   RATIONAL_LINEAR (survives: pinch stratum), NONRATIONAL_OR_COMPLEX
   (killed: single-Laurent-eps strata argument);
 * Jordan: geometric multiplicity of each family via SVD rank;
 * trace rule: tr D1 must be exactly rational-linear in eps (build integrity);
 * exclusion hardening: no p/q with q<=maxq within dist_tol of any killed
   real eigenvalue at >=2 eps; complex pairs excluded by |Im| > dist_tol.
"""
import mpmath as mp
mp.mp.dps = max(mp.mp.dps, 50)      # set FIRST
from fractions import Fraction


def _ratl(x, maxden=64, tol=None):
    tol = tol or mp.mpf('1e-15')
    f = Fraction(float(mp.re(x))).limit_denominator(maxden)
    err = abs(x - mp.mpf(f.numerator) / f.denominator)
    return (f, float(err)) if err < tol and abs(mp.im(x)) < tol else (None, float(err))


def _fr2mp(s):
    f = Fraction(s)
    return mp.mpf(f.numerator) / f.denominator


def pair_and_fit(EVs, eps_vals):
    """EVs: list of eigenvalue lists (one per eps).  Greedy pairing between
    first and last lists; line a+b*eps must hit an eigenvalue at every middle
    eps.  Returns list of {a,b (mpc), match_score}."""
    E0, E1 = EVs[0], EVs[-1]
    e0, e1 = eps_vals[0], eps_vals[-1]
    mids = list(range(1, len(EVs) - 1))
    used, fams = set(), []
    for l0 in E0:
        best = None
        for j, lj in enumerate(E1):
            if j in used:
                continue
            b = (lj - l0) / (e1 - e0); a = l0 - b * e0
            score = mp.mpf(0)
            for mi in mids:
                score = max(score, min(abs(a + b * eps_vals[mi] - z) for z in EVs[mi]))
            if best is None or score < best[0]:
                best = (score, j, a, b)
        score, j, a, b = best
        used.add(j)
        fams.append({'a': a, 'b': b, 'match_score': float(score)})
    return fams


def classify(EVs, eps_vals, maxden=64, kill_integer_eps_indep=True):
    """Returns (families, groups).  families: per-eigenvalue rationalized fits;
    groups: identical (a,b) grouped with multiplicity + class + all-eps verify."""
    fams = pair_and_fit(EVs, eps_vals)
    out_f = []
    for f in fams:
        ra, _ = _ratl(f['a'], maxden); rb, _ = _ratl(f['b'], maxden)
        ver = None
        if ra is not None and rb is not None:
            av = mp.mpf(ra.numerator) / ra.denominator
            bv = mp.mpf(rb.numerator) / rb.denominator
            ver = max(float(min(abs(av + bv * eps_vals[k] - z) for z in EVs[k]))
                      for k in range(len(EVs)))
        out_f.append({'a': str(ra) if ra is not None else mp.nstr(f['a'], 20),
                      'b': str(rb) if rb is not None else mp.nstr(f['b'], 20),
                      'a_num': f['a'], 'b_num': f['b'],
                      'rational': ra is not None and rb is not None,
                      'match_score': f['match_score'],
                      'verify_maxres_alleps': ver})
    groups = {}
    for f in out_f:
        groups.setdefault((f['a'], f['b'], f['rational']), []).append(f)
    gout = []
    for (a, b, rat), fs in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if rat:
            fa, fb = Fraction(a), Fraction(b)
            if fb == 0 and fa.denominator == 1 and kill_integer_eps_indep:
                cls = 'INT_EPS_INDEP__KILLED_vacuum_theorem'
            elif fb == 0:
                cls = 'RATIONAL_EPS_INDEP__FLAG_review'   # none observed; see README
            else:
                cls = 'RATIONAL_LINEAR__SURVIVES_pinch'
        else:
            cls = 'NONRATIONAL_OR_COMPLEX__KILLED_eps_strata'
        gout.append({'a': a, 'b': b, 'rational': rat, 'mult': len(fs), 'class': cls,
                     'a_num': fs[0]['a_num'], 'b_num': fs[0]['b_num'],
                     'verify_maxres_alleps':
                         max((x['verify_maxres_alleps'] or 1) for x in fs)})
    return out_f, gout


def _lam_of_group(g, eps):
    if g['rational']:
        return _fr2mp(g['a']) + _fr2mp(g['b']) * eps
    return mp.mpc(g['a_num']) + mp.mpc(g['b_num']) * eps


def jordan_info(D1, groups, eps, tol_exp=18, dps=40):
    """Geometric multiplicity (SVD nullity of D1-lam*I) per group at one eps."""
    save = mp.mp.dps; mp.mp.dps = dps
    N = D1.rows
    for g in groups:
        lam = _lam_of_group(g, eps)
        S = mp.svd_c(D1 - lam * mp.eye(N), compute_uv=False)
        nker = sum(1 for s in S if s / S[0] < mp.mpf(10) ** (-tol_exp))
        g['n_jordan_blocks_geo'] = nker
        g['log_branches'] = g['mult'] - nker
        g['jordan'] = 'semisimple' if nker == g['mult'] else \
            f'nilpotent({g["mult"] - nker})'
    mp.mp.dps = save
    return groups


def trace_rule(D1s, eps_vals, maxden=10**4):
    """tr D1 must be exactly a + b*eps (rational charpoly): fit 2 pts,
    rationalize, verify residual at ALL eps.  D1-build integrity gate."""
    ts = [mp.fsum(D1[i, i] for i in range(D1.rows)) for D1 in D1s]
    b = (ts[-1] - ts[0]) / (eps_vals[-1] - eps_vals[0]); a = ts[0] - b * eps_vals[0]
    fa = Fraction(float(mp.re(a))).limit_denominator(maxden)
    fb = Fraction(float(mp.re(b))).limit_denominator(maxden)
    resid = max(abs(mp.mpf(fa.numerator) / fa.denominator
                    + (mp.mpf(fb.numerator) / fb.denominator) * eps_vals[k] - ts[k])
                for k in range(len(eps_vals)))
    return {'tr_D1': f'{fa} + ({fb})*eps', 'max_resid': mp.nstr(resid, 5),
            'pass': float(resid) < 1e-20}


def exclude_strata(groups, EVs=None, eps_vals=None, maxq=10**6,
                   dist_tol=None):
    """Survivors per the single-Laurent-eps argument + exclusion hardening
    (PSLQ-grade rationality scan on the killed eigenvalues at >=2 eps).
    Returns dict(survivors, killed, flagged, n_free_constants, exclusion_rows)."""
    dist_tol = dist_tol or mp.mpf('1e-25')
    surv = [g for g in groups if 'SURVIVES' in g['class']]
    killed = [g for g in groups if 'KILLED' in g['class']]
    flagged = [g for g in groups if 'FLAG' in g['class']]
    rational_groups = [g for g in groups if g['rational']]
    rows = []
    if EVs is not None:
        for k in (0, len(EVs) - 1):
            e = eps_vals[k]
            sv = [_lam_of_group(g, e) for g in rational_groups]
            for z in EVs[k]:
                if sv and min(abs(z - s) for s in sv) < mp.mpf('1e-20'):
                    continue                      # a rational-classified branch
                im = abs(mp.im(z))
                if im > dist_tol:
                    rows.append({'eps_idx': k, 'lam': mp.nstr(z, 25),
                                 'excluded_by': 'Im != 0', 'im': mp.nstr(im, 5)})
                    continue
                x = mp.re(z)
                f = Fraction(float(x)).limit_denominator(maxq)
                err = abs(x - mp.mpf(f.numerator) / f.denominator)
                rows.append({'eps_idx': k, 'lam': mp.nstr(z, 25),
                             'best_pq': f'{f.numerator}/{f.denominator}',
                             'dist': mp.nstr(err, 5),
                             'excluded_by':
                                 f'no rational q<={maxq} within {mp.nstr(dist_tol, 2)}'
                                 if err > dist_tol else 'RATIONAL FOUND — FLAG'})
    return {'survivors': surv, 'killed': killed, 'flagged': flagged,
            'n_free_constants': sum(g['mult'] for g in surv),
            'n_killed': sum(g['mult'] for g in killed),
            'exclusion_rows': rows,
            'exclusion_flags': [r for r in rows if 'FLAG' in r['excluded_by']]}

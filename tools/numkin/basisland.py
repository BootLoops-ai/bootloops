#!/usr/bin/env python3
# numkin.basisland — basis-constrained exact function landing.
# Canonical home is tools/numkin/ (import as numkin.basisland; CLI: python3
# tools/numkin/basisland.py rank|solve).
"""basisland — basis-constrained exact function landing.

When an unknown function lives in a KNOWN finite-dimensional space
(span of exactly-evaluable basis functions), land it exactly from
#samples ~ dimension, INDEPENDENT of the number of variables:
mixed point-value + linear-functional constraints -> one
overdetermined exact linear solve (multi-prime modular rref + CRT +
rational reconstruction) with surplus-row zero-inconsistency and
mandatory exact verification as the certificate.

Reference measurement: a 14-kinematic-variable, dim-543 function pair
landed from 596 exact rows in 1.2 min wall, 596/596 exact-verified.

verbs:
  rank  IN.json OUT.json [nprimes=3]
        IN: {"rows": [{sym: "p/q", ...}, ...]}  (spanning-set
        evaluations; one row per sample point; absent sym = 0)
        OUT: {"rank": r, "pivot_symbols": [...], "n_rows": n}
        Pivot choice greedy in sorted-symbol order; the rank must
        agree across all nprimes primes (assert).
  solve IN.json OUT.json [nprimes=96]
        IN: {"cols": [pivot syms], "rows": [{"w": {sym: "p/q"},
             "rhs": ["p/q", ...]}, ...]}   (multi-RHS supported)
        OUT: {"x": [[...], ...] one exact vector per RHS,
              "n_rows", "n_primes", "exact_verified": [n, ...],
              "max_height_digits"}
        FAIL-CLOSED: rank deficit -> exit 2 with missing columns;
        ratrec overflow -> exit 3 (rerun with more primes); ANY
        exact-verification miss -> exit 4. Surplus rows are the
        certificate — never trim them.
"""
import json, sys
from fractions import Fraction as Fr
from math import isqrt
import numpy as np
from sympy import prevprime


def _primes(n):
    ps, q = [], (1 << 25)
    for _ in range(n):
        q = int(prevprime(q))
        ps.append(q)
    return ps


def _modval(fr, p):
    return fr.numerator % p * pow(fr.denominator % p, p - 2, p) % p


def rank_verb(inp, out, nprimes=3):
    rows = [{k: Fr(v) for k, v in r.items()}
            for r in json.load(open(inp))['rows']]
    syms = sorted({k for r in rows for k in r})
    results = []
    for p in _primes(nprimes):
        M = np.zeros((len(rows), len(syms)), dtype=np.int64)
        six = {s: j for j, s in enumerate(syms)}
        for i, r in enumerate(rows):
            for k, v in r.items():
                M[i, six[k]] = _modval(v, p)
        piv, r_ = [], 0
        for c in range(len(syms)):
            nz = np.nonzero(M[r_:, c])[0]
            if not len(nz):
                continue
            pr = r_ + int(nz[0])
            M[[r_, pr]] = M[[pr, r_]]
            M[r_] = M[r_] * pow(int(M[r_, c]), p - 2, p) % p
            col = M[:, c].copy()
            col[r_] = 0
            M = (M - np.outer(col, M[r_])) % p
            piv.append(syms[c])
            r_ += 1
        results.append(tuple(piv))
    assert len(set(results)) == 1, \
        f'rank/pivots disagree across primes: {[len(r) for r in results]}'
    json.dump({'rank': len(results[0]), 'pivot_symbols': list(results[0]),
               'n_rows': len(rows)}, open(out, 'w'))
    print(f'rank {len(results[0])} of {len(syms)} candidate symbols '
          f'({len(rows)} sample rows, {nprimes} primes agree)')


def solve_verb(inp, out, nprimes=96):
    D = json.load(open(inp))
    cols = D['cols']
    cix = {s: j for j, s in enumerate(cols)}
    rank = len(cols)
    rows = [({cix[k]: Fr(v) for k, v in r['w'].items() if k in cix},
             [Fr(x) for x in r['rhs']]) for r in D['rows']]
    nrhs = len(rows[0][1])
    assert all(len(r[1]) == nrhs for r in rows)
    NR = len(rows)
    assert NR >= rank, f'underdetermined: {NR} rows < {rank} cols'
    sols = []
    for p in _primes(nprimes):
        M = np.zeros((NR, rank + nrhs), dtype=np.int64)
        for i, (w, rhs) in enumerate(rows):
            for j, v in w.items():
                M[i, j] = _modval(v, p)
            for m, v in enumerate(rhs):
                M[i, rank + m] = _modval(v, p)
        r_, missing = 0, []
        for c in range(rank):
            nz = np.nonzero(M[r_:, c])[0]
            if not len(nz):
                missing.append(cols[c])
                continue
            pr = r_ + int(nz[0])
            M[[r_, pr]] = M[[pr, r_]]
            M[r_] = M[r_] * pow(int(M[r_, c]), p - 2, p) % p
            col = M[:, c].copy()
            col[r_] = 0
            M = (M - np.outer(col, M[r_])) % p
            r_ += 1
        if missing:
            print(f'RANK DEFICIT ({len(missing)}): {missing[:8]}')
            sys.exit(2)
        bad = int(np.count_nonzero(M[rank:, rank:]))
        assert bad == 0, f'prime {p}: {bad} surplus-row inconsistencies'
        sols.append((p, [[int(M[j, rank + m]) for j in range(rank)]
                         for m in range(nrhs)]))
    xs = []
    for m in range(nrhs):
        vec = []
        for j in range(rank):
            r, md = 0, 1
            for p, s in sols:
                rj = s[m][j]
                if md == 1:
                    r, md = rj, p
                else:
                    t_ = (rj - r) % p * pow(md % p, p - 2, p) % p
                    r, md = r + md * t_, md * p
            a0, a1, x0, x1 = md, r, 0, 1
            bound = isqrt(md >> 1)
            while a1 > bound:
                qq = a0 // a1
                a0, a1 = a1, a0 - qq * a1
                x0, x1 = x1, x0 - qq * x1
            if x1 == 0 or abs(x1) > bound:
                print(f'RATREC FAIL rhs {m} coord {j} — rerun with '
                      f'more primes (had {nprimes})')
                sys.exit(3)
            vec.append(Fr(a1 if x1 > 0 else -a1, abs(x1)))
        xs.append(vec)
    verified = []
    for m in range(nrhs):
        nok = 0
        for w, rhs in rows:
            if sum(v * xs[m][j] for j, v in w.items()) == rhs[m]:
                nok += 1
        verified.append(nok)
        if nok != NR:
            print(f'EXACT VERIFY MISS: rhs {m}: {nok}/{NR}')
            sys.exit(4)
    hmax = max(max(abs(x.numerator), x.denominator)
               for v in xs for x in v)
    json.dump({'x': [[str(x) for x in v] for v in xs], 'n_rows': NR,
               'n_primes': nprimes, 'exact_verified': verified,
               'max_height_digits': len(str(hmax))}, open(out, 'w'))
    print(f'LANDED: {nrhs} vector(s) of dim {rank}; exact-verified '
          f'{verified} of {NR} rows; max height ~1e{len(str(hmax))-1}')


if __name__ == '__main__':
    verb = sys.argv[1]
    if verb == 'rank':
        rank_verb(sys.argv[2], sys.argv[3],
                  int(sys.argv[4]) if len(sys.argv) > 4 else 3)
    elif verb == 'solve':
        solve_verb(sys.argv[2], sys.argv[3],
                   int(sys.argv[4]) if len(sys.argv) > 4 else 96)
    else:
        sys.exit(f'unknown verb {verb} (rank|solve)')

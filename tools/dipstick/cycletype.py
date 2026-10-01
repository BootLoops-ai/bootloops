#!/usr/bin/env python3
# DIPSTICK member — verb `cycletype`: Frobenius cycle-type statistics of the
# critical points; assumes squarefree reduction at the sampled primes (a
# non-squarefree prime shows as a repeated factor — discard that prime).
"""galoiscycle — Frobenius cycle-type statistics for the Galois/monodromy
group of a 0-dimensional critical system (validated on the X(3,6) system,
where the group comes out S26).
Usage: point MS_SRC at a critcount-style .ms file; the script re-runs
msolve at several primes, factors each eliminant with flint, and reports
factor-degree multisets (= Frobenius cycle types), the rational-factor
subset-sum test (transitivity), parity, and the Jordan-criterion summary."""
import os, re, subprocess, sys
import flint
from sympy import isprime, nextprime
MS_SRC = sys.argv[1] if len(sys.argv) > 1 else 'crit36.ms'
PRIMES = [1073741827, 1073741831, 1073741833, 1073741839, 1073741789,
          1073741783]
src = open(MS_SRC).read().split('\n')
pats = []
for p0 in PRIMES:
    p = p0 if isprime(p0) else int(nextprime(p0))
    fn, out = f'gc_{p}.ms', f'gc_{p}.out'
    open(fn, 'w').write(src[0] + '\n' + str(p) + '\n' + '\n'.join(src[2:]))
    subprocess.run([os.environ.get('DIPSTICK_MSOLVE', os.path.join(
                        os.path.dirname(os.path.abspath(__file__)),
                        'run_msolve_capped.sh')),
                    '-f', fn, '-o', out, '-t', '8'],
                   env=dict(os.environ, MSOLVE_CAP_GB='30'),
                   capture_output=True, text=True, timeout=1800)
    s = open(out).read().replace('\n', ' ')
    m = re.search(r'\[1, \[\[(\d+), \[([0-9, -]+)\]\]', s)
    if not m:
        print(f'p={p}: parse fail'); continue
    coeffs = [int(x) % p for x in m.group(2).split(',')]
    fac = flint.nmod_poly(coeffs, p).factor()
    degs = sorted(int(g.degree()) for g, mult in fac[1]
                  for _ in range(int(mult)))
    pats.append(degs)
    print(f'p={p}: cycle type {degs}', flush=True)
from itertools import combinations
def sums(pat):
    return {sum(pat[i] for i in c) for r in range(1, len(pat))
            for c in combinations(range(len(pat)), r)}
common = set.intersection(*[sums(p) for p in pats]) if pats else set()
n = sum(pats[0])
print(f'\nrational-factor candidate degrees: {sorted(common) or "NONE (irreducible/transitive)"}')
for p in pats:
    print(f'  {p}: {"ODD" if sum(l-1 for l in p) % 2 else "even"} permutation')
print(f'Jordan: a prime cycle q with n/2 < q <= n-3 forces A_n;'
      f' + any odd element forces S_n.  (n = {n})')

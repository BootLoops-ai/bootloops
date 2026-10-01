"""selftest.py — battery for nullcal + seam + dcm.

Runnable from any cwd on a copy. Legs:
  S0 fixture integrity: every file in FIXTURES_SHA256 verified by full
     sha256; `--planted fixture` mutates a temp copy and must be REFUSED BY
     NAME [selftest.fixture-sha] with exit code 3.
  S1 nullcal reference reproduction (Federalist 65-paper lists): counts and
     verdicts exact.
  S2 seam fast path == brute enumeration (pinned + 25 seeded), Fraction
     equality.
  S3 seam total-probability: sum over all length-4 sequences == 1 exactly.
  S4 named refusals fire: nullcal.EngineUnavailable outside the package;
     seam.ProfileError on a non-normalised profile; then an explicit exact
     brute engine drives calibrate_objects end to end.
  S5 dcm micro-gates: mixture and change-point on a tiny instance agree with
     direct enumeration over all token assignments (exact).
Exit 0 with '[SELFTEST] ALL PASS'.
"""
import hashlib
import json
import os
import random
import sys
from fractions import Fraction
from itertools import product
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dcm
import nullcal
import seam

FIXTURES_SHA256 = {
    'fixtures/fixtures_nullcal.json':
        'ccb98921766e7b5d94183ea2a131afe4215720102e1e47808e48128f4b348767',
    'fixtures/in_MW30.json':
        'd4cb02620c64bcbf8b7f649fe3c7ef272952e72fb932644ea8488ed6321fff86',
    'fixtures/outA_MW30.json':
        '65322aa50b6ad508bb8e6bddb079b6319e595f1966a47b4b43ea5c8e848e879c',
    'fixtures/our_no55_finals.json':
        '32f6806a42a8e8e736887368dbd206369fa09c349a3df8001410eb50e1502b10',
}


def sha256_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def verify_fixtures(base=HERE, planted=None):
    for rel, want in FIXTURES_SHA256.items():
        if want is None:
            continue
        p = os.path.join(base, rel)
        got = sha256_file(p) if planted is None or rel != planted else \
            hashlib.sha256(open(p, 'rb').read() + b'x').hexdigest()
        if got != want:
            raise SystemExit(
                f'REFUSED [selftest.fixture-sha]: {rel} sha256 {got[:16]} != '
                f'pinned {want[:16]} (exit 3)')
    return True


if '--planted' in sys.argv and 'fixture' in sys.argv:
    try:
        verify_fixtures(planted='fixtures/fixtures_nullcal.json')
    except SystemExit as e:
        print(e)
        sys.exit(3)
    print('PLANTED FIXTURE NOT CAUGHT')
    sys.exit(1)

verify_fixtures()
print('[S0 PASS] pinned fixtures verified by sha256')

# ---- S1
fx = json.load(open(os.path.join(HERE, 'fixtures', 'fixtures_nullcal.json')))
expect = {'fw70': (3, 'TAIL'), 'cc318': (6, 'INSIDE-NULL'),
          'top300': (9, 'INSIDE-NULL')}
for s, (n_exp, verdict_exp) in expect.items():
    r = nullcal.calibrate(fx['no55_claim'][s], fx['known_bfs'][s])
    assert r['n_known'] == 65 and r['n_at_or_above'] == n_exp \
        and r['verdict'] == verdict_exp, (s, r)
print('[S1 PASS] Federalist reference reproduced: ' + ', '.join(
    f'{n}/65 ({s}, {v})' for s, (n, v) in expect.items()))

# ---- S2
pA = [Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)]
pB = [Fraction(1, 6), Fraction(1, 3), Fraction(1, 2)]
zf, _ = seam.z_seam([0, 0, 1, 2, 2, 2, 1, 0], pA, pB)
assert zf == seam.z_seam_brute([0, 0, 1, 2, 2, 2, 1, 0], pA, pB)
rng = random.Random(20260910)
for _ in range(25):
    k = rng.randint(2, 5)
    raw = [rng.randint(1, 9) for _ in range(k)]
    qa = [Fraction(x, sum(raw)) for x in raw]
    raw2 = [rng.randint(1, 9) for _ in range(k)]
    qb = [Fraction(x, sum(raw2)) for x in raw2]
    sq = [rng.randrange(k) for _ in range(rng.randint(1, 12))]
    assert seam.z_seam(sq, qa, qb)[0] == seam.z_seam_brute(sq, qa, qb)
print('[S2 PASS] seam fast path == brute enumeration (pinned + 25 seeded)')

# ---- S3
qA = [Fraction(2, 3), Fraction(1, 3)]
qB = [Fraction(1, 4), Fraction(3, 4)]
tot = Fraction(0)
for sq in product(range(2), repeat=4):
    tot += seam.z_seam(list(sq), qA, qB)[0]
assert tot == 1, tot
print('[S3 PASS] sum over all length-4 sequences == 1 exactly')

# ---- S4
try:
    nullcal.blend_vs_pure_bf([2, 1])
    raise AssertionError('EngineUnavailable did not fire')
except nullcal.EngineUnavailable as e:
    assert 'nullcal.engine-unresolved' in str(e)
try:
    seam.z_seam([0, 1], [Fraction(1, 2), Fraction(1, 3)],
                [Fraction(1, 2), Fraction(1, 2)])
    raise AssertionError('ProfileError did not fire')
except seam.ProfileError as e:
    assert 'seam.profile-normalization' in str(e)


def Z_brute(U, g):
    """Exact exchangeable-mixture evidence by the collapsed identity
    (Dir(1) convention), independent of any vendored engine."""
    U = list(U)
    k, N = len(U), sum(U)
    if g == 1:
        z = Fraction(factorial(k - 1))
        for u in U:
            z *= factorial(u)
        return z / factorial(N + k - 1)
    # CT dict DP over per-channel splits into g parts
    from collections import defaultdict
    cur = defaultdict(int)
    cur[(0,) * (g - 1)] = 1
    for u in U:
        nxt = defaultdict(int)
        for e, ct in cur.items():
            for split in product(range(u + 1), repeat=g - 1):
                if sum(split) <= u:
                    ne = tuple(a + b for a, b in zip(e, split))
                    nxt[ne] += ct
        cur = nxt
    tot = Fraction(0)
    for e, ct in cur.items():
        mg = N - sum(e)
        if mg < 0:
            continue
        w = Fraction(1)
        for m in list(e) + [mg]:
            w *= Fraction(factorial(m), factorial(m + k - 1))
        tot += ct * w
    pref = Fraction(factorial(g - 1) * factorial(k - 1) ** g,
                    factorial(N + g - 1))
    for u in U:
        pref *= factorial(u)
    return pref * tot


r = nullcal.calibrate_objects([3, 1, 2], [[2, 2, 2], [4, 1, 1], [1, 3, 2]],
                              2, Z=Z_brute)
assert r['n_known'] == 3 and r['verdict'] in ('INSIDE-NULL', 'TAIL',
                                              'ABOVE-NULL')
print('[S4 PASS] named refusals fire; brute-engine calibrate_objects runs')

# ---- S5: dcm micro-gates vs direct enumeration
kn = 2500                      # kappa = 5/2
nH = [3, 1]
nM = [1, 3]
aH, DH = dcm.profile_scaled(nH, kn)
aM, DM = dcm.profile_scaled(nM, kn)
c = [2, 1]
r = dcm.dcm_mixture(c, aH, DH, aM, DM, kn)


def dcm_stretch_direct(xs, a, D):
    """P of an ordered stretch with per-token categories xs."""
    N = len(xs)
    num = Fraction(1)
    cnt = {}
    for v in xs:
        num *= Fraction(a[v] + cnt.get(v, 0) * D, D)
        cnt[v] = cnt.get(v, 0) + 1
    den = Fraction(1)
    for j in range(N):
        den *= Fraction(kn + 1000 * j, 1000)
    return num / den


# direct token-mixture: for the multiset {0,0,1}, sum over which tokens are
# H's (f-polynomial integrated: each assignment with K H-tokens weighs
# B(K+1, N-K+1) = K!(N-K)!/(N+1)!)
toks = [0, 0, 1]
N = len(toks)
Zdir = Fraction(0)
for mask in product((0, 1), repeat=N):
    K = sum(mask)
    hs = [t for t, m in zip(toks, mask) if m]
    ms = [t for t, m in zip(toks, mask) if not m]
    Zdir += dcm_stretch_direct(hs, aH, DH) * dcm_stretch_direct(ms, aM, DM) \
        * Fraction(factorial(K) * factorial(N - K), factorial(N + 1))
# dcm_mixture's Z counts category-splits (unordered within category), and
# sequence evidence is order-free given counts: divide the assignment sum by
# nothing — both are the ordered-sequence evidence. Compare directly.
assert r['Z'] == Zdir, (r['Z'], Zdir)
# change-point micro-gate: U[K] weights vs direct two-orientation enumeration
cp = dcm.dcm_seam(toks, aH, DH, aM, DM, kn)
UK = [Fraction(0)] * (N + 1)
for K in range(N + 1):
    UK[K] = (dcm_stretch_direct(toks[:K], aH, DH)
             * dcm_stretch_direct(toks[K:], aM, DM)
             + dcm_stretch_direct(toks[-K:] if K else [], aH, DH)
             * dcm_stretch_direct(toks[:N - K], aM, DM))
hmean = sum(Fraction(K) * u for K, u in enumerate(UK)) / (N * sum(UK))
assert cp['hshare_mean'] == hmean, (cp['hshare_mean'], hmean)
print('[S5 PASS] dcm mixture and change-point == direct enumeration (exact)')

print('[SELFTEST] ALL PASS')

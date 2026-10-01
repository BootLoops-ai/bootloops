#!/usr/bin/env python3
"""lambda_exact — exact Lambda-coalescent merger rates and expected branch-length spectrum.

Register: EXACT. Every rate, every expected length and every spectrum entry is
a fractions.Fraction; no float touches any returned value or any sign
decision. Floats appear only in display fields (identity check 6 and the CLI
printout). Standard library only.

What it computes
----------------
For a Lambda-coalescent (Pitman 1999; Sagitov 1999) with b active lineages,
each specific set of k of them merges at rate

    lambda_{b,k} = integral_0^1 x^{k-2} (1-x)^{b-k} Lambda(dx),   2 <= k <= b.

Rate families, all exactly rational at rational parameters:
  * kingman_rate():  Lambda = delta_0, lambda_{b,k} = 1{k=2}.
  * beta_rate(alpha):  Lambda = Beta(2-alpha, alpha) (Schweinsberg 2003),
        lambda_{b,k} = B(k-alpha, b-k+alpha) / B(2-alpha, alpha)
                     = [prod_{j=2}^{k-1} (j-alpha)] [prod_{j=0}^{b-k-1} (alpha+j)] / (b-1)!
    by Pochhammer telescoping of the Gamma ratios, hence rational in alpha.
    alpha = 2 gives exactly the Kingman rates; alpha = 1 gives exactly the
    Bolthausen-Sznitman coalescent, lambda_{b,k} = (k-2)! (b-k)! / (b-1)!.
  * dirac_rate(psi):  Lambda = delta_psi (the psi-coalescent of Eldon &
    Wakeley 2006), lambda_{b,k} = psi^(k-2) (1-psi)^(b-k). psi = 0 gives
    exactly Kingman; psi = 1 gives the star coalescent.
  * msprime_dirac_rate(psi, c):  the two-atom mixture
    Lambda = delta_0 + c psi^2 delta_psi, lambda_{b,k} = 1{k=2} + c psi^k (1-psi)^(b-k)
    (the parametrization of msprime's DiracCoalescent(psi, c) at ploidy = 1; provided so the
    exact recursion can be cross-checked against that simulator).

expected_lengths(n, lam) returns h[m][i], m = 1..n, i = 1..m-1, where h(m,i)
is the expected total branch length subtending exactly i of m sampled leaves
under the exchangeable Lambda-coalescent with rates lam(b, k), constant
population size, time in the units in which Lambda is the stated measure
(Kingman: each pair coalesces at rate 1). E[L_i] = h(n,i) is the expected
unfolded site-frequency spectrum up to the factor theta/2 (infinite sites,
mutations at rate theta/2 along every branch): E[xi_i] = (theta/2) h(n,i).
xi_hat(n, lam) returns the normalized spectrum h(n,i) / sum_j h(n,j),
i = 1..n-1, which is theta-free and sums to 1 exactly.

Recursion (first-transition decomposition; equivalent to the recursions of
Birkner, Blath & Eldon 2013 for the expected SFS under Lambda-coalescents):
  h(1,.) = 0
  h(m,i) = (m/L_m) 1{i=1}
         + sum_{k=2}^{m-1} P_{m,k} [ (j/m') h(m',j)|_{j=i-k+1, 1<=j<=m'-1}
                                   + ((m'-i)/m') h(m',i)|_{1<=i<=m'-1} ]
  with m' = m-k+1, L_m = sum_k C(m,k) lambda_{m,k}, P_{m,k} = C(m,k) lambda_{m,k}/L_m.
Justification: before the first event all m branches subtend 1 leaf
(duration 1/L_m); after a k-merger the remaining process is again the same
exchangeable Lambda-coalescent on m' lineages carrying leaf multiplicities
(k,1,...,1); by exchangeability the expected length of j-subtending branches
that contain the multiplicity-k lineage is (j/m') h(m',j).
Kingman closed form (Fu 1995): h(n,i) = 2/i exactly, so xi_hat_i = (1/i)/H_{n-1}.

Built-in identity checks (identity_checks(n), exact unless noted):
  check1  Kingman through the recursion == 2/i, normalized == (1/i)/H_{n-1}
  check2  beta_rate(2) == Kingman, rates entrywise and spectrum
  check3  dirac_rate(0) spectrum == Kingman spectrum
  check4  beta_rate(1) rates == Bolthausen-Sznitman (k-2)!(b-k)!/(b-1)!
  check5  dirac_rate(1) spectrum == star coalescent (all mass on i = 1)
  check6  (display float) l_inf distance of the Beta(199/100) spectrum to Kingman

Reference data (reference/lambda_coalescent/, n = N = 20):
  spectra_n20.json  exact normalized spectra, as num/den strings, on the grid
      Kingman; Beta(2-alpha, alpha) for alpha in {199/100, 19/10, 7/4, 3/2,
      5/4, 11/10, 101/100}; Dirac for psi in {1/100, 1/10, 1/4, 1/2, 3/4,
      9/10}; plus the identity-check block and each row's exact margin under
      the functional below. load_reference_spectra() reads it;
      reference_grid() rebuilds the same grid of rate functions.
  kingman_class_functional_n20.json  a fixed linear functional on the
      19-simplex, F(x) = <w, x> + w0 with 19 exact rationals w and the dyadic
      constant w0 = -4898883/2^26. The file carries F as a separating
      functional for the Kingman class at n = 20: F >= 0 on every expected
      normalized spectrum of the Kingman coalescent under a variable
      population size (a curve-positivity certificate of the kind
      popcorn.certificates.positivity decides; this module takes it as given
      and does not re-derive it), so an exact F(x) < 0 places x outside that
      class. load_witness() reads it and witness_margin(w, w0, x) evaluates F
      exactly. On the reference grid F > 0 for Kingman and Beta(199/100) and
      F < 0 for every other Beta row and every Dirac row.

Cost: pure-Python Fraction arithmetic, O(n^3) operations on rationals whose
size grows with n and with the rate family. Indicative single-core wall
times: any family at n = 20 in milliseconds; Beta(3/2) about 2 s at n = 100
and about 17 s at n = 200; Dirac(1/2) about 3 min at n = 200. Kingman is
O(n^2) and instantaneous. Peak memory stays in the tens of MB.

Scope: constant population size, unfolded spectrum, expected values only (no
second moments), simultaneous multiple mergers (Xi-coalescents) not covered.

CLI:  python3 lambda_exact.py [--n N] [--out FILE]
  runs the identity checks (exit 1 if any of checks 1-5 fails), evaluates the
  reference grid at sample size n, prints one line per row (with the exact
  sign of F when n = 20), and writes the results JSON to FILE if given
  (same schema as spectra_n20.json plus display floats). Writes nothing
  unless --out is passed.

References
  Pitman, J. (1999). Coalescents with multiple collisions. Ann. Probab. 27,
    1870-1902.
  Sagitov, S. (1999). The general coalescent with asynchronous mergers of
    ancestral lines. J. Appl. Probab. 36, 1116-1125.
  Bolthausen, E. & Sznitman, A.-S. (1998). On Ruelle's probability cascades
    and an abstract cavity method. Comm. Math. Phys. 197, 247-276.
  Schweinsberg, J. (2003). Coalescent processes obtained from supercritical
    Galton-Watson processes. Stochastic Process. Appl. 106, 107-139.
  Eldon, B. & Wakeley, J. (2006). Coalescent processes when the distribution
    of offspring number among individuals is highly skewed. Genetics 172,
    2621-2633.
  Birkner, M., Blath, J. & Eldon, B. (2013). Statistical properties of the
    site-frequency spectrum associated with Lambda-coalescents. Genetics 195,
    1037-1053.
  Fu, Y.-X. (1995). Statistical properties of segregating sites. Theor.
    Popul. Biol. 48, 172-197.
"""
import argparse
import json
import os
import sys
from fractions import Fraction as F
from math import comb, factorial

N = 20                                 # sample size of the shipped reference data
HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.join(HERE, 'reference', 'lambda_coalescent')
SPECTRA_REF = os.path.join(REFERENCE_DIR, 'spectra_n20.json')
WITNESS_REF = os.path.join(REFERENCE_DIR, 'kingman_class_functional_n20.json')

# ---------------------------------------------------------------- rate families

def beta_rate(alpha):
    """lambda_{b,k} for Beta(2-alpha, alpha), telescoped, rational in alpha."""
    alpha = F(alpha)
    def lam(b, k):
        num = F(1)
        for j in range(2, k):          # prod_{j=2}^{k-1} (j - alpha)
            num *= (j - alpha)
        for j in range(0, b - k):      # prod_{j=0}^{b-k-1} (alpha + j)
            num *= (alpha + j)
        return num / factorial(b - 1)
    return lam

def dirac_rate(psi):
    """lambda_{b,k} for pure Dirac Lambda = delta_psi (Eldon-Wakeley)."""
    psi = F(psi)
    def lam(b, k):
        return psi ** (k - 2) * (1 - psi) ** (b - k)
    return lam

def kingman_rate():
    def lam(b, k):
        return F(1) if k == 2 else F(0)
    return lam

def msprime_dirac_rate(psi, c):
    """msprime DiracCoalescent exact model: Lambda = delta_0 + c psi^2 delta_psi."""
    psi, c = F(psi), F(c)
    def lam(b, k):
        base = F(1) if k == 2 else F(0)
        return base + c * psi ** k * (1 - psi) ** (b - k)
    return lam

# ---------------------------------------------------------------- exact recursion

def expected_lengths(n, lam):
    """h[m][i] (Fraction), m=1..n, i=1..m-1; h[m][i] = E total length subtending i."""
    h = {1: {}}
    for m in range(2, n + 1):
        rates = {k: lam(m, k) for k in range(2, m + 1)}
        assert all(r >= 0 for r in rates.values()), f"negative rate at b={m}"
        L = sum(comb(m, k) * rates[k] for k in range(2, m + 1))
        assert L > 0, f"total rate not positive at b={m}"
        row = {i: F(0) for i in range(1, m)}
        row[1] = F(m) / L
        for k in range(2, m):          # k=m leaves 1 lineage: no further branches
            P = comb(m, k) * rates[k] / L
            if P == 0:
                continue
            mp = m - k + 1
            for i in range(1, m):
                acc = F(0)
                j = i - k + 1          # branch contains the merged (mult-k) lineage
                if 1 <= j <= mp - 1:
                    acc += F(j, mp) * h[mp][j]
                if 1 <= i <= mp - 1:   # branch avoids the merged lineage
                    acc += F(mp - i, mp) * h[mp][i]
                row[i] += P * acc
        h[m] = row
    return h

def xi_hat(n, lam):
    h = expected_lengths(n, lam)[n]
    tot = sum(h.values())
    return [h[i] / tot for i in range(1, n)]

# ---------------------------------------------------------------- identity checks

def identity_checks(n):
    """Built-in exact identities (checks 1-5, booleans) plus display float check 6."""
    checks = {}
    # check 1: Kingman rates through the recursion == closed form 2/i, and
    #          normalized == (1/i)/H_{n-1}. Exact rational equality demanded.
    h = expected_lengths(n, kingman_rate())[n]
    g1a = all(h[i] == F(2, i) for i in range(1, n))
    H = sum(F(1, j) for j in range(1, n))
    xk = xi_hat(n, kingman_rate())
    g1b = all(xk[i - 1] == F(1, i) / H for i in range(1, n))
    checks['check1_kingman_closed_form'] = bool(g1a and g1b)

    # check 2: Beta family at alpha=2 == Kingman rates, entrywise and spectrum.
    lb = beta_rate(2)
    g2a = all(lb(b, k) == (F(1) if k == 2 else F(0))
              for b in range(2, n + 1) for k in range(2, b + 1))
    g2b = xi_hat(n, beta_rate(2)) == xk
    checks['check2_beta_alpha2_is_kingman'] = bool(g2a and g2b)

    # check 3: Dirac at psi=0 == Kingman.
    g3 = xi_hat(n, dirac_rate(0)) == xk
    checks['check3_dirac_psi0_is_kingman'] = bool(g3)

    # check 4: Beta at alpha=1 rates == Bolthausen-Sznitman (k-2)!(b-k)!/(b-1)!.
    l1 = beta_rate(1)
    g4 = all(l1(b, k) == F(factorial(k - 2) * factorial(b - k), factorial(b - 1))
             for b in range(2, n + 1) for k in range(2, b + 1))
    checks['check4_beta_alpha1_is_bolthausen_sznitman'] = bool(g4)

    # check 5: Dirac at psi=1 == star coalescent (all mass on singletons).
    xs = xi_hat(n, dirac_rate(1))
    g5 = xs[0] == 1 and all(x == 0 for x in xs[1:])
    checks['check5_dirac_psi1_is_star'] = bool(g5)

    # check 6 (display-only float): alpha=199/100 close to Kingman in l_infty.
    xa = xi_hat(n, beta_rate(F(199, 100)))
    d = max(abs(float(a - b)) for a, b in zip(xa, xk))
    checks['check6_beta_alpha199_100_linf_to_kingman_float'] = d
    return checks

# ---------------------------------------------------------------- reference data

BETA_ALPHAS = ['199/100', '19/10', '7/4', '3/2', '5/4', '11/10', '101/100']
DIRAC_PSIS = ['1/100', '1/10', '1/4', '1/2', '3/4', '9/10']

def reference_grid():
    """[(family, parameter name, value or None, rate function)] of the reference grid."""
    return ([('kingman', 'const', None, kingman_rate())]
            + [('beta', 'alpha', F(a), beta_rate(F(a))) for a in BETA_ALPHAS]
            + [('dirac', 'psi', F(p), dirac_rate(F(p))) for p in DIRAC_PSIS])

def load_reference_spectra(path=None):
    """The shipped reference spectra (dict; exact entries as num/den strings)."""
    with open(path or SPECTRA_REF) as fh:
        return json.load(fh)

def load_witness(path=None):
    """(w, w0) of the shipped Kingman-class functional at n = N, as Fractions."""
    with open(path or WITNESS_REF) as fh:
        d = json.load(fh)
    w = [F(s) for s in d['w']]
    w0 = F(d['w0_dyadic'])
    assert w0 == F(-4898883, 1 << 26), "w0 pin mismatch"
    assert len(w) == N - 1
    return w, w0

def witness_margin(w, w0, x):
    """Affine functional on the simplex: F = <w, xi-hat> + w0 (exact)."""
    return sum(wi * xi for wi, xi in zip(w, x)) + w0

# ---------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="exact Lambda-coalescent spectra on the reference grid")
    ap.add_argument('--n', type=int, default=N, help=f"sample size (default {N})")
    ap.add_argument('--out', default=None,
                    help="write the results JSON to this file (default: print only)")
    a = ap.parse_args(argv)
    n = a.n

    out = {'n': n, 'checks': identity_checks(n), 'grid': []}
    hard = [k for k, v in out['checks'].items() if v is False]
    if hard:
        print('IDENTITY CHECK FAILURE:', hard)
        if a.out:
            with open(a.out, 'w') as fh:
                json.dump(out, fh, indent=1)
        return 1

    use_f = (n == N) and os.path.isfile(WITNESS_REF)
    if use_f:
        w, w0 = load_witness()
    for fam, pname, pval, lam in reference_grid():
        x = xi_hat(n, lam)
        row = {
            'family': fam, 'param': pname,
            'value': str(pval) if pval is not None else None,
            'xi_hat': [str(v) for v in x],
            'xi_hat_float': [float(v) for v in x],
        }
        tag = f"{fam} {pname}={pval}" if pval is not None else fam
        line = f"{tag:22s}  xi_hat_1 = {float(x[0]):.6f}"
        if use_f:
            m = witness_margin(w, w0, x)
            row['witness_margin_exact'] = str(m)
            row['witness_margin_float'] = float(m)
            row['witness_sign'] = ('>=0' if m >= 0 else '<0')
            line += f"  F = {float(m):+.6e}  sign {'>=0' if m >= 0 else '<0'}"
        out['grid'].append(row)
        print(line)

    print('checks:', out['checks'])
    if a.out:
        with open(a.out, 'w') as fh:
            json.dump(out, fh, indent=1)
        print('wrote', a.out)
    return 0

if __name__ == '__main__':
    sys.exit(main())

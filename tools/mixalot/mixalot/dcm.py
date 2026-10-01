"""dcm — exact Dirichlet-compound-multinomial burstiness arm.

Alternative to the iid multinomial idealisation: a stretch of n tokens by
author A has category probabilities theta_d ~ Dirichlet(kappa * theta_A) and
tokens iid theta_d; integrating theta_d out,

  P_A(x) = G(kappa)/G(kappa+n) prod_w G(alpha_w + x_w)/G(alpha_w),
  alpha_w = kappa * theta_A(w),  theta_A = (n_A + 1/2)/(N_A + (V+1)/2)

(add-half profiles over V listed categories + 1 residual bucket, so every
token counts). ONE kappa is shared by the components, fitted on the known-pure
objects and then treated as an exact rational (the reference convention:
float ML rounded to 3 decimals, as in the independent DCM reference computation
that the fixtures reproduce).

Token-mixture: each token is component H's with probability f ~ Beta(1,1); the
H tokens form one DCM stretch, the rest another. P(D|f) = sum_K A_K f^K
(1-f)^(N-K) with one exact integer polynomial per category (the same
linear-factor lattice class as the frozen-profile blend). Contiguous
(change-point) variant on the ordered token sequence included.

All evidences exact rationals. Reference fixtures: in_MW30.json holds the word
counts of Federalist No. 55 and of the Hamilton and Madison training sets over
a 30-word function-word list (derived from the public-domain texts);
outA_MW30.json holds the No. 55 outputs of an independent implementation of
the same model, which dcm.py reproduces to their printed digits.

Scaled-integer conventions (all integers):
  kappa = kn/1000;  D_A = 1000*(2*N_A + V + 1);  a_Aw = kn*(2*n_Aw + 1)
  so alpha_Aw = a_Aw / D_A exactly and sum_w alpha_Aw = kappa exactly.
"""
from fractions import Fraction
from math import comb


def profile_scaled(train_counts, kn):
    """(a, D): a_w = kn*(2 n_w + 1), D = 1000*(2 N + V+1); alpha_w = a_w/D."""
    V1 = len(train_counts)                      # V listed + 1 residual
    N = sum(train_counts)
    D = 1000 * (2 * N + V1)
    a = [kn * (2 * n + 1) for n in train_counts]
    assert sum(a) * 1000 == kn * D
    return a, D


def _rising(a, D, k):
    """Numerator of (a/D)^(k rising) * D^k = prod_{j<k} (a + j*D)."""
    out = 1
    for j in range(k):
        out *= a + j * D
    return out


def _kappa_rising(kn, K):
    """R[K] = prod_{j<K}(kn + 1000 j)  (= kappa^(K rising) * 1000^K)."""
    R = [1]
    for j in range(K):
        R.append(R[-1] * (kn + 1000 * j))
    return R


def dcm_pure(x, a, D, kn):
    """Exact sequence evidence of one DCM stretch: a Fraction."""
    N = sum(x)
    num = 1
    for xw, aw in zip(x, a):
        if xw:
            num *= _rising(aw, D, xw)
    Rk = _kappa_rising(kn, N)
    return Fraction(num * 1000 ** N, D ** N * Rk[N])


def dcm_mixture(c, aH, DH, aM, DM, kn):
    """Exact token-mixture evidence and posterior (f ~ Beta(1,1)).

    Returns dict: Z (Fraction), mean (Fraction), BF_mix_vs_H / _vs_M /
    BF_H_vs_M (Fractions), T (integer weights, T[K] prop. to
    A'_K B(K+1,N-K+1)), and float q05/q50/q95 of the f posterior.
    """
    V1 = len(c)
    N = sum(c)
    order = sorted(range(V1), key=lambda w: -c[w])
    poly = [1]
    for w in order:
        cw = c[w]
        if cw == 0:
            break
        wp = [comb(cw, k) * _rising(aH[w], DH, k) * _rising(aM[w], DM, cw - k)
              for k in range(cw + 1)]
        new = [0] * (len(poly) + cw)
        for i, pi in enumerate(poly):
            if pi:
                for j, wj in enumerate(wp):
                    if wj:
                        new[i + j] += pi * wj
        poly = new
    assert len(poly) == N + 1
    Rk = _kappa_rising(kn, N)
    RN = Rk[N]
    powDH = [1]
    powDM = [1]
    for _ in range(N):
        powDH.append(powDH[-1] * DH)
        powDM.append(powDM[-1] * DM)
    fact = [1]
    for k in range(1, N + 2):
        fact.append(fact[-1] * k)
    T = [poly[K] * powDH[N - K] * powDM[K] * (RN // Rk[K]) * (RN // Rk[N - K])
         * fact[K] * fact[N - K] for K in range(N + 1)]
    S = sum(T)
    # Z = S * 1000^N / (DH^N DM^N RN^2 (N+1)!)
    Z = Fraction(S * 1000 ** N, powDH[N] * powDM[N] * RN * RN * fact[N + 1])
    mean = Fraction(sum(T[K] * (K + 1) for K in range(N + 1)), (N + 2) * S)
    BF_mix_M = Fraction(S, poly[0] * powDH[N] * RN * fact[N + 1])
    BF_mix_H = Fraction(S, poly[N] * powDM[N] * RN * fact[N + 1])
    BF_H_M = Fraction(poly[N] * powDM[N], poly[0] * powDH[N])
    # float quantiles of the Beta-mixture posterior sum_K t_K Beta(K+1, N-K+1)
    from math import log10
    logs = [log10(t) if t else None for t in T]
    mx = max(l for l in logs if l is not None)
    wts = [(10.0 ** (l - mx)) if l is not None else 0.0 for l in logs]
    sw = sum(wts)
    wts = [x / sw for x in wts]
    q = _beta_mixture_quantiles(wts, N, (0.05, 0.5, 0.95))
    return {'Z': Z, 'mean': mean, 'BF_mix_vs_M': BF_mix_M,
            'BF_mix_vs_H': BF_mix_H, 'BF_H_vs_M': BF_H_M, 'T': T,
            'q05': q[0.05], 'q50': q[0.5], 'q95': q[0.95]}


def _beta_mixture_quantiles(wts, N, qs):
    from scipy.special import betainc
    import numpy as np
    j = np.arange(N + 1)
    w = np.asarray(wts)

    def cdf(x):
        return float(np.dot(w, betainc(j + 1, N - j + 1, x)))

    out = {}
    for q in qs:
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            if cdf(mid) < q:
                lo = mid
            else:
                hi = mid
        out[q] = (lo + hi) / 2
    return out


def _seam_T(tid, aH, DH, aM, DM, kn):
    """Weights T[K] for 'component H wrote the FIRST K tokens of tid'.
    Prefix pass: PH[K] = prod_{t<=K} (aH[v_t] + (#v_t among first t-1) * DH);
    suffix pass symmetric for M. Exact integers throughout."""
    N = len(tid)
    PH = [1] * (N + 1)
    cnt = {}
    for t in range(N):
        v = tid[t]
        PH[t + 1] = PH[t] * (aH[v] + cnt.get(v, 0) * DH)
        cnt[v] = cnt.get(v, 0) + 1
    SM = [1] * (N + 1)
    cnt = {}
    for t in range(N - 1, -1, -1):
        v = tid[t]
        SM[t] = SM[t + 1] * (aM[v] + cnt.get(v, 0) * DM)
        cnt[v] = cnt.get(v, 0) + 1
    Rk = _kappa_rising(kn, N)
    RN = Rk[N]
    powDH = [1]
    powDM = [1]
    for _ in range(N):
        powDH.append(powDH[-1] * DH)
        powDM.append(powDM[-1] * DM)
    return [PH[K] * SM[K] * powDH[N - K] * powDM[K] * (RN // Rk[K])
            * (RN // Rk[N - K]) for K in range(N + 1)]


def dcm_seam(tid, aH, DH, aM, DM, kn):
    """Exact contiguous (change-point) model, the reference convention
    (reference method A): H wrote the first K tokens (orientation F) or the last
    K tokens (orientation L), prior uniform over (orientation, K); U[K] =
    TF[K] + TL[K] is the posterior weight of Hamilton-COUNT K. Returns:
      hshare_mean = sum U[K] K / (N sum U)   (mean H token share)
      q05/q95     discrete K/N quantiles
      BF_seam_vs_M  mean of the 2(N-1) interior configs vs pure M
      BF_H_vs_M     U[N]/U[0]
      P_pureM/P_pureH  endpoint posterior mass
    """
    N = len(tid)
    TF = _seam_T(tid, aH, DH, aM, DM, kn)
    TL = _seam_T(list(reversed(tid)), aH, DH, aM, DM, kn)
    assert TF[0] == TL[0] and TF[N] == TL[N]
    U = [TF[K] + TL[K] for K in range(N + 1)]
    SU = sum(U)
    hshare_mean = Fraction(sum(U[K] * K for K in range(N + 1)), N * SU)

    def dq(qn, qd):
        acc = 0
        for K in range(N + 1):
            acc += U[K]
            if acc * qd >= qn * SU:
                return Fraction(K, N)

    return {
        'hshare_mean': hshare_mean,
        'q05': dq(1, 20), 'q95': dq(19, 20),
        'BF_seam_vs_M': Fraction(sum(U[1:N]), 2 * (N - 1)) / Fraction(U[0], 2),
        'BF_H_vs_M': Fraction(U[N], U[0]),
        'P_pureM': Fraction(U[0], SU), 'P_pureH': Fraction(U[N], SU),
        'U': U,
    }


def fit_kappa(doc_counts, doc_profiles, kn_grid):
    """Exact-comparison kappa fit on known-pure objects: for each candidate
    kn (kappa = kn/1000), total log evidence sum_d ln dcm_pure(x_d | kappa *
    theta_d); returns (best_kn, table). Comparisons via 60-dps logs of the
    exact rationals (differences are O(1) or larger; ties broken exact).
    doc_profiles: per-doc (train_counts_of_its_author_without_it) lists.
    """
    import mpmath as mp
    mp.mp.dps = 60
    table = []
    for kn in kn_grid:
        tot = mp.mpf(0)
        for x, tc in zip(doc_counts, doc_profiles):
            a, D = profile_scaled(tc, kn)
            z = dcm_pure(x, a, D, kn)
            tot += mp.log(mp.mpf(z.numerator)) - mp.log(mp.mpf(z.denominator))
        table.append((kn, tot))
    best = max(table, key=lambda t: t[1])
    return best[0], table

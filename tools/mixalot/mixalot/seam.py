"""seam — exact contiguous-block (prefix/suffix) evidence.

Model M_seam(A->B): an ordered sequence of channel draws has tokens 1..s from
frozen profile A and s+1..N from frozen profile B, with the seam s uniform on
{0..N} (the endpoints are the two pure models). Exact evidence with
single-denominator profile integers:

  Z_seam = 1/(N+1) * sum_s P_A(s) S_B(s) D_A^{N-s} D_B^{s} / (D_A^N D_B^N)

Unlike a token-level iid mixture (which is unidentifiable in the number of
components — it only ever sees the mean profile), the seam model uses ORDER
and is identifiable. Its BF against best-pure takes the same null calibration
as any mixed-vs-pure factor (nullcal.calibrate).

API: z_seam(seq, pA, pB) -> (Z exact Fraction, posterior summary dict).
seq = list of channel indices in order; pA, pB = profiles as Fractions.
"""
from fractions import Fraction
from math import gcd, log10


class ProfileError(ValueError):
    """Named refusal: a profile that is not a normalised probability
    vector. Callers and batteries match on '[seam.profile-normalization]'."""


def _profile_ints(p):
    D = 1
    for x in p:
        x = Fraction(x)
        D = D * x.denominator // gcd(D, x.denominator)
    nums = [Fraction(x).numerator * (D // Fraction(x).denominator) for x in p]
    if sum(nums) != D:
        raise ProfileError(
            'REFUSED [seam.profile-normalization]: profile does not sum to 1')
    return nums, D


def z_seam(seq, pA, pB):
    """Exact seam evidence + seam-posterior summary for order A->B."""
    nA, DA = _profile_ints(pA)
    nB, DB = _profile_ints(pB)
    N = len(seq)
    PA = [1] * (N + 1)
    for t in range(N):
        PA[t + 1] = PA[t] * nA[seq[t]]
    SB = [1] * (N + 1)
    for t in range(N - 1, -1, -1):
        SB[t] = SB[t + 1] * nB[seq[t]]
    powA = [1] * (N + 1)
    powB = [1] * (N + 1)
    for j in range(1, N + 1):
        powA[j] = powA[j - 1] * DA
        powB[j] = powB[j - 1] * DB
    terms = [PA[s] * SB[s] * powA[N - s] * powB[s] for s in range(N + 1)]
    tot = sum(terms)
    Z = Fraction(tot, (N + 1) * powA[N] * powB[N])
    logs = [log10(t) for t in terms]
    mx = max(logs)
    w = [10.0 ** (l - mx) for l in logs]
    sw = sum(w)
    w = [x / sw for x in w]
    mode = max(range(N + 1), key=lambda s: w[s])
    cdf, q = 0.0, {}
    for s in range(N + 1):
        cdf += w[s]
        for qq in (0.05, 0.5, 0.95):
            if qq not in q and cdf >= qq:
                q[qq] = s
    post = {'mode': mode, 'q05': q.get(0.05), 'q50': q.get(0.5),
            'q95': q.get(0.95), 'endpoint_mass': w[0] + w[N]}
    return Z, post


def z_seam_brute(seq, pA, pB):
    """Reference: direct enumeration of the N+1 seams with per-token
    Fractions. Exponentially slower constants; for gate-scale instances."""
    N = len(seq)
    tot = Fraction(0)
    for s in range(N + 1):
        term = Fraction(1)
        for t in range(N):
            term *= Fraction(pA[seq[t]]) if t < s else Fraction(pB[seq[t]])
        tot += term
    return tot / (N + 1)

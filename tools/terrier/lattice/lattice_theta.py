r"""
lattice_theta.py — Jacobi/lattice theta series + flatness factors + CVP for {Z^n, D_n, D_n^+, E8}.

Reusable standalone tool (does not inherit the lattice wing's frozen K3 conventions).

Functions
---------
theta2, theta3, theta4(q)            : Jacobi theta nullwerte (|q|<1).
theta_series(name, q)                : Theta_Lambda(q) = sum_{v in Lambda} q^{|v|^2}
                                       for name in {"Z8","D8","D8plus","E8","Leech","Zn","Dn"}.
flatness_eps(name, sigma, conv=...)  : eps_{Lambda_L}(sigma) for the GKP encoding
                                       Lambda_S = sqrt(2)*Lambda, Lambda_L = Lambda/sqrt(2),
                                       Lambda unit-covolume self-dual.
                                       conv="sqrt2":   Theta_{sqrt2*Lambda}(exp(-4*pi^2*sigma^2)) - 1
                                       conv="llbs"   : Theta_{Lambda}(exp(-4*pi^2*sigma^2)) - 1
                                                      (= sum_{v in Lambda_L^* \ 0} exp(-2*pi^2*sigma^2*|v|^2))
cvp_Zn, cvp_Dn, cvp_Dnplus, cvp_E8   : closest-vector (vectorized over leading axis).
e8_shell_radii(kmax)                 : sqrt(2k), k=1..kmax (E8 root-lattice shells).

Conventions: theta argument is the *nome* q (NOT tau). Theta_Lambda(q) := sum q^{|v|^2}.
All CVP functions take/return float64 arrays shaped (..., n).

README line: lattice_theta.py — Jacobi theta_{2,3,4}, Theta_{Z^n,D_n,D_n^+,E8,Leech}, GKP flatness eps, vectorized CVP.
"""
from __future__ import annotations
import numpy as np

# ---------------------------------------------------------------- Jacobi theta
def theta3(q: float, N: int = 200) -> float:
    """theta_3(q) = 1 + 2*sum_{n>=1} q^{n^2}. Auto-truncates."""
    s = 1.0
    for n in range(1, N + 1):
        t = q ** (n * n)
        s += 2.0 * t
        if abs(t) < 1e-18 * abs(s):
            break
    return s

def theta2(q: float, N: int = 200) -> float:
    """theta_2(q) = 2*sum_{n>=0} q^{(n+1/2)^2}."""
    s = 0.0
    for n in range(0, N + 1):
        t = q ** ((n + 0.5) ** 2)
        s += 2.0 * t
        if abs(t) < 1e-18 * max(abs(s), 1e-300):
            break
    return s

def theta4(q: float, N: int = 200) -> float:
    """theta_4(q) = 1 + 2*sum_{n>=1} (-1)^n q^{n^2}."""
    s = 1.0
    for n in range(1, N + 1):
        t = ((-1) ** n) * (q ** (n * n))
        s += 2.0 * t
        if abs(t) < 1e-18 * abs(s):
            break
    return s

# --------------------------------------------------------- lattice theta series
def _sigma3(n: int) -> int:
    s = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            s += d ** 3
            e = n // d
            if e != d:
                s += e ** 3
        d += 1
    return s

def theta_E8(q: float, K: int = 100) -> float:
    """Theta_{E8}(q) = E_4(q) = 1 + 240*sum_{n>=1} sigma_3(n) q^{2n}.
    (E8 minimal norm = 2, so exponents are even integers.)"""
    s = 1.0
    for n in range(1, K + 1):
        t = 240.0 * _sigma3(n) * (q ** (2 * n))
        s += t
        if abs(t) < 1e-18 * abs(s):
            break
    return s

def theta_E8_jacobi(q: float) -> float:
    """Theta_{E8}(q) = (theta_2^8 + theta_3^8 + theta_4^8)/2 — cross-check form."""
    return 0.5 * (theta2(q) ** 8 + theta3(q) ** 8 + theta4(q) ** 8)

def _tau11(n: int, _cache={0: 0, 1: 1, 2: -24, 3: 252, 4: -1472, 5: 4830,
                           6: -6048, 7: -16744, 8: 84480, 9: -113643, 10: -115920}) -> int:
    return _cache.get(n, None)

def theta_Leech(q: float, K: int = 30) -> float:
    """Theta_{Leech}(q) = E_4(q)^3 - 720*Delta(q) computed via Jacobi-form identity
       = (theta_2^8+theta_3^8+theta_4^8)^3/8 - 45/16*(theta_2 theta_3 theta_4)^8.
    Minimal norm 4; first terms 1 + 196560 q^4 + ... Exact for |q|<1."""
    t2, t3, t4 = theta2(q), theta3(q), theta4(q)
    E4 = 0.5 * (t2 ** 8 + t3 ** 8 + t4 ** 8)
    # Delta(tau) where q_nome = q here uses q^{|v|^2} convention; with q = exp(i*pi*tau),
    # Delta in this nome is (1/256)*(theta_2 theta_3 theta_4)^8.
    Delta = (t2 * t3 * t4) ** 8 / 256.0
    return E4 ** 3 - 720.0 * Delta

def theta_series(name: str, q: float, n: int = 8) -> float:
    """Theta_Lambda(q) for named lattice. `n` only used for "Zn"/"Dn"."""
    name = name.upper().replace("_", "").replace("+", "PLUS")
    if name in ("Z8", "ZN"):
        nn = 8 if name == "Z8" else n
        return theta3(q) ** nn
    if name in ("D8", "DN"):
        nn = 8 if name == "D8" else n
        return 0.5 * (theta3(q) ** nn + theta4(q) ** nn)
    if name in ("D8PLUS", "DNPLUS"):
        nn = 8 if name == "D8PLUS" else n
        # D_n^+ = D_n  union  ([1/2]^n + D_n);  Theta = 1/2(theta3^n+theta4^n) + 1/2 theta2^n
        return 0.5 * (theta3(q) ** nn + theta4(q) ** nn) + 0.5 * theta2(q) ** nn
    if name == "E8":
        return theta_E8_jacobi(q)
    if name == "LEECH":
        return theta_Leech(q)
    raise ValueError(f"unknown lattice name: {name!r}")

# -------------------------------------------------------------- flatness factor
def flatness_eps(name: str, sigma: float, conv: str = "sqrt2") -> float:
    """GKP flatness factor for encoding Lambda_S = sqrt2*Lambda, Lambda_L = Lambda/sqrt2.

    conv="sqrt2":   eps = Theta_{sqrt2*Lambda}(q) - 1 = Theta_Lambda(q^2) - 1,  q=exp(-4*pi^2*sigma^2).
    conv="llbs"   : eps = sum_{v in Lambda_L^* \\ 0} exp(-2*pi^2*sigma^2*|v|^2)
                       = Theta_{sqrt2*Lambda}(exp(-2*pi^2*sigma^2)) - 1
                       = Theta_Lambda(exp(-4*pi^2*sigma^2)) - 1.
    NOTE the two conventions differ by a factor of 2 in the exponent; both are
    monotone in the same Theta_Lambda so *ordering* across lattices is identical.
    """
    if conv == "sqrt2":
        q = np.exp(-4.0 * np.pi ** 2 * sigma ** 2)
        return theta_series(name, q ** 2) - 1.0
    elif conv == "llbs":
        q = np.exp(-4.0 * np.pi ** 2 * sigma ** 2)
        return theta_series(name, q) - 1.0
    else:
        raise ValueError("conv must be 'sqrt2' or 'llbs'")

# --------------------------------------------------------------------------- CVP
def cvp_Zn(x: np.ndarray) -> np.ndarray:
    return np.round(x)

def cvp_Dn(x: np.ndarray) -> np.ndarray:
    """Closest D_n point: round to Z^n; if coord-sum is odd, flip the coord with
    maximal rounding residual (Conway-Sloane Ch.20 Alg.)."""
    r = np.round(x)
    par = np.mod(np.sum(r, axis=-1), 2)  # (...,)
    # second-nearest integer per coord
    diff = x - r
    # cost of flipping coord i: 1 - 2|diff_i| (always >=0); choose argmax|diff|
    idx = np.argmax(np.abs(diff), axis=-1)
    flip = np.where(diff[np.arange(diff.shape[0]) if diff.ndim == 2 else ..., ] >= 0, 1.0, -1.0) \
        if False else None  # placeholder; do it explicitly below
    # explicit handling for 1D and 2D batch
    r2 = r.copy()
    if r.ndim == 1:
        if par % 2 == 1:
            i = int(idx)
            r2[i] += 1.0 if diff[i] >= 0 else -1.0
    else:
        bad = (par.astype(int) % 2 == 1)
        rows = np.nonzero(bad)[0]
        cols = idx[bad]
        steps = np.where(diff[rows, cols] >= 0, 1.0, -1.0)
        r2[rows, cols] += steps
    return r2

def cvp_Dnplus(x: np.ndarray) -> np.ndarray:
    """CVP for D_n^+ = D_n union (g + D_n), g = (1/2,...,1/2). NOTE: D_8^+ == E8."""
    n = x.shape[-1]
    g = 0.5 * np.ones(n)
    a = cvp_Dn(x)
    b = cvp_Dn(x - g) + g
    da = np.sum((x - a) ** 2, axis=-1)
    db = np.sum((x - b) ** 2, axis=-1)
    pick_b = db < da
    out = a.copy()
    if out.ndim == 1:
        return b if pick_b else a
    out[pick_b] = b[pick_b]
    return out

def cvp_E8(x: np.ndarray) -> np.ndarray:
    """CVP for E8 (== D_8^+)."""
    return cvp_Dnplus(x)

CVP = {"Z8": cvp_Zn, "D8": cvp_Dn, "D8plus": cvp_Dnplus, "E8": cvp_E8}

def e8_shell_radii(kmax: int = 6):
    return [np.sqrt(2.0 * k) for k in range(1, kmax + 1)]

# ----------------------------------------------------------------- self-test
if __name__ == "__main__":
    q = 0.05
    a = theta_E8(q)
    b = theta_E8_jacobi(q)
    c = theta_series("D8plus", q)
    print(f"E8 sigma3-form  = {a:.12g}")
    print(f"E8 jacobi-form  = {b:.12g}")
    print(f"D8+ jacobi-form = {c:.12g}   (should equal E8: D8+ == E8 in dim 8)")
    print(f"Leech(q=0.05)   = {theta_Leech(0.05):.12g}  (~ 1 + 196560*q^4 = {1+196560*0.05**4:.6g})")
    print(f"theta3^8(q)     = {theta_series('Z8', q):.12g}")
    # CVP sanity
    rng = np.random.default_rng(0)
    X = rng.normal(0, 0.2, size=(5, 8))
    print("cvp_E8 residual norms:", np.sqrt(np.sum((X - cvp_E8(X))**2, axis=-1)))

"""Floating-point check of ISD_CUTOFF_DERIVATION.md on 200 synthetic Hodge
frames: H symmetric positive definite, the |B| vacuum identity, the theorem
N_ISD <= |B|, tightness of the AM-GM step, the lam^2 gram_E corollary, the
sec 6 Gauss-reduction chain, and the b=1 witness family."""
import numpy as np
rng = np.random.default_rng(20260716)
S = np.zeros((4, 4)); S[0, 2] = S[1, 3] = 1; S[2, 0] = S[3, 1] = -1  # witness.sig(4)
def point():
    while True:
        Pi = rng.normal(size=4) + 1j * rng.normal(size=4)
        if (1j * np.conj(Pi) @ S @ Pi).real <= 0: Pi = np.conj(Pi)
        N0 = (1j * np.conj(Pi) @ S @ Pi).real
        for _ in range(50):
            D = rng.normal(size=4) + 1j * rng.normal(size=4)
            D = D - (Pi @ S @ D) / (Pi @ S @ np.conj(Pi)) * np.conj(Pi)  # Griffiths
            D = D - (np.conj(Pi) @ S @ D) / (np.conj(Pi) @ S @ Pi) * Pi  # type (4,2)=0
            N1 = (-1j * np.conj(D) @ S @ D).real
            if N1 > 0: return Pi, D, N0, N1
def gram(M, f, h):
    return (f @ M @ f) * (h @ M @ h) - (f @ M @ h) ** 2
ok = {k: 0 for k in ["sym", "pd", "vac", "B>0", "thm", "tight", "cor", "red"]}
T = 200
for trial in range(T):
    Pi, D, N0, N1 = point()
    SP, SD = S @ Pi, S @ D
    H = 2 * np.real(np.outer(SP, np.conj(SP))) / N0 + 2 * np.real(np.outer(SD, np.conj(SD))) / N1
    lam = np.linalg.eigvalsh(H)[0]
    ok["sym"] += np.allclose(H, H.T); ok["pd"] += lam > 1e-12
    a, b = rng.normal(size=2) + 1j * rng.normal(size=2)
    g = a * D + b * np.conj(Pi)                      # ISD: H^{2,1}+H^{0,3} charge vec
    r, t = rng.normal(), abs(rng.normal()) + 0.05
    h = -np.imag(g) / t; f = np.real(g) - (r / t) * np.imag(g)
    assert np.allclose(g, f - (r + 1j * t) * h)      # G3 = f - tau h reconstructs
    B = f @ S @ h
    x = f - r * h
    ok["vac"] += np.isclose(abs(B), (x @ H @ x + t * t * (h @ H @ h)) / (2 * t))  # |B| law
    ok["B>0"] += B < 0                               # file-cert orientation: B NEGATIVE
    NISD = np.sqrt(max(gram(H, f, h), 0.0))
    ok["thm"] += NISD <= abs(B) * (1 + 1e-7)         # sec3(vi), |B| form
    ts = np.geomspace(1e-3, 1e3, 4001); rs = np.linspace(r - 8, r + 8, 801)
    xs = f[None, :] - rs[:, None] * h[None, :]
    q = np.einsum("ij,jk,ik->i", xs, H, xs)
    mn = np.min(q[:, None] / (2 * ts[None, :]) + ts[None, :] / 2 * (h @ H @ h))
    ok["tight"] += mn <= NISD * 1.001 + 1e-9         # AM-GM bound is tight
    fi, hi = rng.integers(-9, 10, 4), rng.integers(-9, 10, 4)
    if gram(np.eye(4), fi, hi) > 0:
        ok["cor"] += gram(H, fi, hi) >= lam**2 * gram(np.eye(4), fi, hi) - 1e-9
    else: ok["cor"] += 1
    F, Hh = fi.astype(float), hi.astype(float)       # Gauss-reduce pair, sec 6 chain
    for _ in range(60):
        if F @ F < Hh @ Hh: F, Hh = Hh, -F
        m = round((F @ Hh) / (Hh @ Hh)) if Hh @ Hh > 0 else 0
        if m == 0: break
        F = F - m * Hh
    gE = gram(np.eye(4), F, Hh)
    ok["red"] += (gE <= 0) or (F @ F + Hh @ Hh <= (8 / 3) * gE + 1e-9)
print("trials", T, ok)
Pi, D, N0, N1 = point()
SP, SD = S @ Pi, S @ D
H = 2 * np.real(np.outer(SP, np.conj(SP))) / N0 + 2 * np.real(np.outer(SD, np.conj(SD))) / N1
lam = np.linalg.eigvalsh(H)[0]
print("witness family b=1 members (p,q,r,s); synthetic lam=%.6f" % lam)
for (p, q, r, s) in [(1, 0, 1, 0), (1, 0, 1, 10), (1, 0, 1, 100)]:
    f = np.array([p, q, 0, 0], float); h = np.array([0, 0, r, s], float)
    B = f @ S @ h; gE = gram(np.eye(4), f, h)
    NI = np.sqrt(gram(H, f, h))
    print(" B=%g gram_E=%d N_ISD=%.4f lam*sqrt(gE)=%.4f Tier1-excl(lam^2 gE>256): %s bound-ok:%s"
          % (B, gE, NI, lam * np.sqrt(gE), lam**2 * gE > 256, NI >= lam * np.sqrt(gE) - 1e-9))

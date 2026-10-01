"""Sign/orientation probe for ISD_CUTOFF_DERIVATION.md sec 2-3: 300 synthetic
ISD vacua per pairing order; reports sign(f^T Sigma h), the |B| identity and
N_ISD <= |B| in both orientations."""
import numpy as np
rng = np.random.default_rng(7)
S = np.zeros((4, 4)); S[0, 2] = S[1, 3] = 1; S[2, 0] = S[3, 1] = -1
def point(sgn):
    """sgn=+1: enforce file certificates N0=i Pib.S.Pi>0, N1=-i Db.S.D>0.
       sgn=-1: enforce the order-reversed convention (both flipped)."""
    while True:
        Pi = rng.normal(size=4) + 1j * rng.normal(size=4)
        if sgn * (1j * np.conj(Pi) @ S @ Pi).real <= 0: Pi = np.conj(Pi)
        N0 = sgn * (1j * np.conj(Pi) @ S @ Pi).real
        for _ in range(80):
            D = rng.normal(size=4) + 1j * rng.normal(size=4)
            D = D - (Pi @ S @ D) / (Pi @ S @ np.conj(Pi)) * np.conj(Pi)   # int Om^chi=0
            D = D - (np.conj(Pi) @ S @ D) / (np.conj(Pi) @ S @ Pi) * Pi  # int Ombar^chi=0
            assert abs(Pi @ S @ D) < 1e-10 and abs(np.conj(Pi) @ S @ D) < 1e-10
            N1 = sgn * (-1j * np.conj(D) @ S @ D).real
            if N1 > 0: return Pi, D, N0, N1
def run(sgn, T=300):
    bad_id = bad_thm = 0; Bs = []
    for _ in range(T):
        Pi, D, N0, N1 = point(sgn)
        SP, SD = S @ Pi, S @ D
        H = (2 * np.real(np.outer(SP, np.conj(SP))) / N0
             + 2 * np.real(np.outer(SD, np.conj(SD))) / N1)
        a, b = rng.normal(size=2) + 1j * rng.normal(size=2)
        g = a * D + b * np.conj(Pi)              # ISD charge vector
        r, t = rng.normal(), abs(rng.normal()) + 0.05
        h = -np.imag(g) / t; f = np.real(g) - (r / t) * np.imag(g)
        B = f @ S @ h; x = f - r * h
        rhs = (x @ H @ x + t * t * (h @ H @ h)) / (2 * t)
        Bs.append(int(np.sign(B)))
        if not np.isclose(abs(B), rhs): bad_id += 1     # |B| = ||G||^2/(2t)?
        NI = np.sqrt(max((f @ H @ f) * (h @ H @ h) - (f @ H @ h) ** 2, 0))
        if NI > abs(B) * (1 + 1e-9): bad_thm += 1       # N_ISD <= |B|?
    print("cert-sign %+d: sign(B) at ISD vacua: value %+d occurs %d/%d;"
          " |B|-identity fails %d; N_ISD<=|B| fails %d"
          % (sgn, Bs[0], Bs.count(Bs[0]), T, bad_id, bad_thm))
run(+1); run(-1)

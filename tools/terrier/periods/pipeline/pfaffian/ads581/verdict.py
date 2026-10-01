#!/usr/bin/env python3
"""Route-A verdict for ads-5-81-3213 (Pfaffian route, direct-evaluation form).

Because the exact 17x17 connection A(s) is pole-order-19 at s=0 in the B17
basis (apparent degeneration; no diagonal gauge exists), the
anchor (s=0) -> s_vac leg is done WITHOUT marching: |s_vac| =
1/180 sits deep inside a PROVEN convergence disk of the restricted rho-jet
Gamma towers (closed-form coefficient bound c(n) <= 2^m, see bound_level),
so every period is evaluated directly from its exact tower with a rigorous
tail ball.  The exact A(s) gates the analytic sector (PA-2, s^3000).

Stages:
  1  jet tower bank C_beta (|beta|<=3) to M_TOW via geff_series.jet_tower
  2  Phi_alpha towers (PL.build_phi), model_periods(MF=17)
  3  frame fit: v^3 Pi_i = sum c_{alpha,q} v^{3-|alpha|} z3^q Phi_alpha
     fit m<=14, HELD OUT m=15..17 exact (gate PA-4a); kernel dim reported
  4  certified evaluation at s_vac = -1/180 (branch log s = ln|s| + i pi,
     eps=+1), arb dps 60 and 90; tail balls from bound_level
  5  contraction (PT.contract conventions verbatim), W0 vs published
usage: verdict.py [dps]
"""
import json, os, sys, time
from fractions import Fraction as Fr
from math import factorial, comb

# reference-data root from the environment (no machine-local default;
# fail-closed). Layout: $TERRIER_REFERENCE_BANKS/cards/ads-5-81-3213_bank/
# with the flux frame, the s^3000 series and dmodule/ (pmodp, boxops,
# connection samples) -- not included in the package, see pfaffian/README.md
_CB = os.environ.get("TERRIER_REFERENCE_BANKS", "")
if not _CB:
    raise SystemExit("REFUSE: TERRIER_REFERENCE_BANKS is not set — point it "
                     "at the reference-data root holding "
                     "cards/ads-5-81-3213_bank/ (not included in the "
                     "package; see periods/pipeline/pfaffian/README.md)")
BANK = os.path.join(_CB, "cards", "ads-5-81-3213_bank")
HERE = os.path.dirname(os.path.abspath(__file__))
# periods/pipeline of this package (family.py, geff_series.py, pipe_lib.py, cards/)
PIPE = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, PIPE)
import family
import geff_series as GS
import pipe_lib as PL

MTOW = 120
FIT_LO, FIT_HI = 50, 55   # rank saturates at m<=50 (rank 80 == rank(120))
NU = (9, 10, 10, 10, 10)
SIGMA = (1, 1, 1, 1, 1)
PUBLISHED = 2.03778e-23

card = family.Card(json.load(open(os.path.join(
    PIPE, "cards", "ads-5-81-3213.json"))))
g2 = json.load(open(os.path.join(BANK, "gate_g2.json")))   # flux frame F, H
SVAC = card.s_star          # pinned vac point (30-digit decimal pin, F6)

# ---- stage 1: tower bank (cached)
tw_path = os.path.join(HERE, f"tower_{MTOW}.json")
if os.path.exists(tw_path):
    tw = json.load(open(tw_path))
    C = {tuple(map(int, k.split(","))):
         [{tuple(map(int, kk.split(","))): Fr(v) for kk, v in d.items()}
          for d in lst] for k, lst in tw["C"].items()}
else:
    t0 = time.time()
    C = GS.jet_tower(card, NU, SIGMA, MTOW)
    json.dump({"M": MTOW, "C": {",".join(map(str, b)):
               [{f"{k[0]},{k[1]}": str(v) for k, v in d.items()}
                for d in lst] for b, lst in C.items()}}, open(tw_path, "w"))
    print(f"[tower] jet_tower M={MTOW} built ({time.time()-t0:.0f}s)",
          flush=True)

# cross-link: C_0 coefficients == reference w0 series (exact)
import gzip
W0S = json.load(gzip.open(os.path.join(BANK, "w0_series_3000.json.gz")))
z5 = (0,) * 5
for m in range(MTOW + 1):
    have = C[z5][m].get((0, 0), Fr(0))
    assert have == Fr(int(W0S["a"][m])), f"C_0[{m}] != reference w0"
print("[gate] C_0 == reference w0 series through m=120  PASS", flush=True)

# ---- stage 2: configure PL manually (op-free; the matrix route owns the DE side)
PL.H = card.h; PL.NU = NU; PL.NPER = 12
PL.KAPPA = card.kappa; PL.A_MAT = card.a_mat; PL.B_VEC = card.b_vec
PL.XI_COEF = card.xi_coef
PL.GV = {tuple(map(int, k.split(","))): int(v) for k, v in
         json.load(open(os.path.join(HERE, "gv_full.json"))).items()}
for d, v in card.gv_pinned.items():
    assert PL.GV.get(d) == v, f"GV pin mismatch {d}"
PL.FFLUX = list(g2["F"]); PL.HFLUX = list(g2["H"])
PL.M_TOW = MTOW; PL.C_JET = C
PL.BETAS = sorted(C)
PL.PHI = {al: PL.build_phi(al) for al in PL.BETAS}
MP = PL.model_periods(FIT_HI)
print(f"[model] model_periods(MF={FIT_HI}) built; BETAS={len(PL.BETAS)}",
      flush=True)

# ---- stage 3: frame fit (restricted dictionary p = 3-|alpha|, q in {0,1})
unknowns = [(al, 3 - sum(al), q) for al in PL.BETAS for q in (0, 1)]
def col_tower(u, M):
    al, p, q = u
    out = {}
    for (lp, e2, e3), v in PL.PHI[al].items():
        key = (lp, p + 2 * e2, q + e3)
        dst = out.setdefault(key, [Fr(0)] * (M + 1))
        f = Fr(-1, 24) ** e2
        for m in range(M + 1):
            if v[m]:
                dst[m] += f * v[m]
    return out
cols = [col_tower(u, FIT_HI) for u in unknowns]

def gauss_rank(rows, bvec, ncol):
    A = [list(r) + [bvec[i]] for i, r in enumerate(rows)]
    n = len(A); rl = 0; piv = []
    for c in range(ncol):
        pr = next((i for i in range(rl, n) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[rl], A[pr] = A[pr], A[rl]
        pv = A[rl][c]
        A[rl] = [x / pv for x in A[rl]]
        for r in range(n):
            if r != rl and A[r][c] != 0:
                f = A[r][c]
                A[r] = [A[r][t] - f * A[rl][t] for t in range(ncol + 1)]
        piv.append(c); rl += 1
    for r in range(rl, n):
        assert A[r][ncol] == 0, "inconsistent fit system"
    sol = [Fr(0)] * ncol
    for i, c in enumerate(piv):
        sol[c] = A[i][ncol]
    return sol, len(piv)

coeffs, towers = [], []
kmin = 99
for i in range(PL.NPER):
    rhs = MP[i]
    keys = sorted(set(k for c in cols for k in c) | set(rhs))
    rows, bvec = [], []
    for key in keys:
        for m in range(FIT_LO + 1):
            row = [c.get(key, None) for c in cols]
            row = [(r[m] if r else Fr(0)) for r in row]
            b = rhs.get(key, None)
            b = b[m] if b else Fr(0)
            if any(row) or b:
                rows.append(row); bvec.append(b)
    sol, rank = gauss_rank(rows, bvec, len(unknowns))
    kdim = len(unknowns) - rank
    kmin = min(kmin, kdim)
    # held-out gate PA-4a: exact match m = 0..FIT_HI on every key
    for key in keys:
        for m in range(FIT_HI + 1):
            lhs = sum(sol[j] * cols[j].get(key, [Fr(0)] * (FIT_HI + 1))[m]
                      for j in range(len(unknowns)) if sol[j])
            b = rhs.get(key, [Fr(0)] * (FIT_HI + 1))[m]
            assert lhs == b, f"PA-4a FAIL period {i} key {key} m {m}"
    coeffs.append({unknowns[j]: sol[j]
                   for j in range(len(unknowns)) if sol[j]})
    # full tower to MTOW
    T = {}
    for j, u in enumerate(unknowns):
        if sol[j] == 0:
            continue
        ct = col_tower(u, MTOW)
        for key, v in ct.items():
            dst = T.setdefault(key, [Fr(0)] * (MTOW + 1))
            for m in range(MTOW + 1):
                if v[m]:
                    dst[m] += sol[j] * v[m]
    towers.append({k: v for k, v in T.items()
                   if any(x != 0 for x in v)})
    print(f"[fit] period {i}: rank {rank}/{len(unknowns)} "
          f"(kernel {kdim}), nnz coeffs {len(coeffs[-1])}, held-out PASS",
          flush=True)
json.dump({"unknown_format": "(alpha, p, q) -> str Fr",
           "coeffs": [{str(k): str(v) for k, v in cd.items()}
                      for cd in coeffs]},
          open(os.path.join(HERE, "f3_frame.json"), "w"))

# uniqueness note: kernel dim per period reported; kernel directions that
# affect the verdict are excluded by held-out if 0; else flagged in report.

# ---- stage 4/5: certified evaluation + contraction
from flint import arb, acb, fmpq, ctx

def fr2arb(x):
    return arb(fmpq(x.numerator, x.denominator))

def bound_level(m):
    """Rigorous bound on |C_beta[m]| coefficients (any |beta|<=3, any
    (e2,e3) sector, zeta-values NOT included -- they are evaluated exactly
    at contraction; here we bound the rational jet coefficients):
      c(n) = C(4n4,2n4)*multinom(2n4; a_1..a_8) <= 2^{10 n4} <= 2^m
      #points(level m) <= (m/9+1)(m/10+1)^3
      psi-jet factors (beta! included, |beta|<=3):
        B1 = 26(1+ln(1+4m)); B2 = B1^2+264; B3 = B1^3+3*B1*264+2880
    -> bound = #points * 2^m * (1+B1+B2+B3).  All inequalities elementary
    (H_N <= 1+ln(1+N), psi' <= zeta2, |psi''| <= 2zeta3+2, sum|l|_1 = 26,
    max|l|_1 = 4; Faa di Bruno with positive majorants)."""
    npts = (fr2arb(Fr(m, 9)) + 1) * (fr2arb(Fr(m, 10)) + 1) ** 3
    B1 = 26 * (1 + (arb(1 + 4 * m)).log())
    B2 = B1 ** 2 + 264
    B3 = B1 ** 3 + 3 * B1 * 264 + 2880
    return npts * arb(2) ** m * (1 + B1 + B2 + B3) * 8   # x8: sector count

def tail_ball(cmax_abs, s_abs, M):
    """sum_{m>M} bound_level(m) * prefmax * s^m, geometric domination:
    ratio bound r = 2*s_abs * ((M+2)/(M+1))^4 * (1+ln(5+4M+4))/(1+ln(4M+5))
    ... conservatively r = 2*s_abs*1.5 (valid for M >= 20)."""
    r = 2 * fr2arb(s_abs) * arb(1.5)
    assert r < 1
    t = bound_level(M + 1) * fr2arb(s_abs) ** (M + 1) / (1 - r)
    return (t * cmax_abs).upper()

def run(dps):
    prec = int(dps * 3.33) + 20
    ctx.prec = prec
    sA = fr2arb(abs(SVAC))
    L = acb(sA.log(), arb.pi())          # log s_vac = ln|s| + i*pi (eps=+1)
    v = acb(0, 2 * arb.pi())             # 2 pi i
    z3 = arb.zeta(arb(3))
    sv = fr2arb(SVAC)
    Pi = []
    for i in range(PL.NPER):
        val = acb(0)
        # prefactor magnitude for tail: sum_u |c_u| * build_phi prefs already
        # inside towers; bound tower coeff by (sum |c_u| * prefmax) * C-bound.
        cmax = arb(0)
        for u, cv in coeffs[i].items():
            # build_phi prefactor bound: sum_b C(al,b) nu^{|al-b|} al!-comb
            # <= (1+max nu)^{|al|} * al! <= 11^3 * 6 < 8000
            cmax += fr2arb(abs(cv)) * 8000
        tb = tail_ball(cmax, abs(SVAC), MTOW)
        for (lp, P, Q), cf in towers[i].items():
            s = arb(0)
            for m in range(MTOW, -1, -1):
                s = s * sv + fr2arb(cf[m])
            pm = arb(0, tb)
            term = acb(s + pm, pm) * L ** lp * v ** P * z3 ** Q
            val += term
        Pi.append(val)
    # contraction (PT.contract conventions)
    v3 = v ** 3
    Piv = [p / v3 for p in Pi]
    def sympl(a, b):
        n = 6
        return sum(a[i] * b[n + i] for i in range(n)) \
            - sum(a[n + i] * b[i] for i in range(n))
    Ff = [acb(x) for x in PL.FFLUX]
    Hf = [acb(x) for x in PL.HFLUX]
    A = sympl(Ff, Piv); B = sympl(Hf, Piv)
    # tau pin: Re tau = 12 (flux frame; tau_hat reproduces it to 22 digits),
    # Im tau = card racetrack pin.  Published W0 is UNDRESSED (paper 09064:
    # W0 = sqrt(2/pi) |W(z*,tau*)|), no e^K dressing.
    tau = acb(arb(12), fr2arb(card.tau_pin))
    sq2pi = (2 / arb.pi()).sqrt()
    W = sq2pi * (A - tau * B)
    PiC = [z.conjugate() for z in Piv]
    zz = sympl(PiC, Piv) * acb(0, -1)
    emKcs = zz.real
    emK = emKcs * 2 * tau.imag
    W0d = abs(W) / emK.sqrt()
    # diagnostics: anc 1/gs pin; self-consistent dilaton stationarity
    tau_anc = acb(arb(12), arb("19.83564187"))
    W_anc = sq2pi * abs(A - tau_anc * B)
    th = (A / B).conjugate()
    W_self = sq2pi * abs(A - th.conjugate() * B)
    return {"W0": abs(W), "W0_dressed": W0d, "absW": abs(W), "emK": emK,
            "A": A, "B": B, "W0_anc_tau": W_anc, "W0_self_tau": W_self,
            "tau_hat": (A / B), "imag_ok": zz.imag.contains(arb(0))}

out = {}
for dps in (60, 90):
    t0 = time.time()
    r = run(dps)
    out[dps] = r
    print(f"[dps {dps}] W0(undressed) = {r['W0']}  ({time.time()-t0:.0f}s)")
    print(f"        W0_anc_tau = {r['W0_anc_tau']}")
    print(f"        W0_self_tau = {r['W0_self_tau']}  "
          f"W0_dressed = {r['W0_dressed']}")
    print(f"        tau_hat = {r['tau_hat']} imag_ok={r['imag_ok']}")
w60, w90 = out[60]["W0"], out[90]["W0"]
mid60, mid90 = float(w60.mid()), float(w90.mid())
agree = abs(mid60 - mid90) <= 1e-15 * abs(mid90)
rel_rad = float(w90.rad()) / abs(mid90) if mid90 else 1.0
dev = abs(mid90 - PUBLISHED) / PUBLISHED
print(f"[PA-6] two-dps agree<=1e-15: {agree}; rel radius {rel_rad:.2e}")
print(f"[VERDICT] W0(90dps) = {mid90:.10e}  published {PUBLISHED:.5e}  "
      f"rel dev {dev:.3e}")
json.dump({"W0_60": str(w60), "W0_90": str(w90),
           "published": PUBLISHED, "rel_dev": dev,
           "two_dps_agree": bool(agree), "rel_radius": rel_rad},
          open(os.path.join(HERE, "verdict_raw.json"), "w"), indent=1)

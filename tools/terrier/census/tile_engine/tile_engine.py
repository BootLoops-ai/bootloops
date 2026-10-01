#!/usr/bin/env python3
"""tile_engine.py — TILE-CRITERION evaluator for continuous exclusion of flux
vacua over closed moduli tiles T of the six-parameter Hulek-Verrill fourfold
family HV4 (arXiv:2404.12422; 29-dimensional period vector). Evaluation
engine = the U1-S certified summation card
(periods/summation/) plus its frame gates, which are MANDATORY at every
gated point. FAIL-CLOSED throughout: a tile where any certified bound is
unavailable stays OPEN — never forced.

Criterion implemented, step by step (the labels OL-2 / OL-4 / OL-5 are the
step names used in the code below):
  OL-2 : Rayleigh box law. With H(t) the Hodge Gram on the tile and lam(T) a
         certified lower bound of its smallest eigenvalue over T, every
         vacuum on T with flux self-pairing Q(G) <= 60 (the tadpole bound)
         has n2(G) <= 60/lam(T).
  OL-5 : candidate box B(T) := {G != 0 : n2(G) <= 60/lam(T) AND Q(G) <= 60},
         Q the exact integer quadratic form G.Sigma.G. B(T) is NEVER
         enumerated here: the criterion is period-side — the engine issues
         (a) the B(T) BOUND as a formula in lam(T) (a number only when a
         certified lam(T) box floor exists), (b) the trivial period-side
         branch lam(T) > 60 => B(T) empty => EXCLUDED, (c) per-flux row
         certificates.
  OL-4 : a tile is EXCLUDED only by the conjunction of (i) certified
         clearance from the special loci (coincidence strata phi_i = phi_j
         and large-complex-structure sheets), (ii) immersivity: the 7x7
         sesquilinear Gram det of (Pi, d_K Pi) separated from 0, (iii) a
         certified lam(T) floor > 0, (iv) per-flux row kill
         inf_T |L_i.G| > 0 for every G in B(T). Layers (ii)-(iv) are
         certified at POINTS (tile center, tile corner); no tile-uniform
         (iii)/(iv) certificate is issued by this engine, so tile_verdict
         returns OPEN with named reasons (branch (b) needs a certified box
         floor). A missing certificate always means OPEN; a control vacuum
         whose rows all CONTAIN 0 additionally pins NOT-EXCLUDABLE-AT-POINT.

Reference data (not included in the package; root supplied with
TERRIER_TILE_BANK): the u1s card directory (summation pass, frame, Hodge
and floor modules with their sha-gated SIGMA; certified corner lam
brackets MU2_POINTS_s{0,1,2}.log) and the per-tile envelope table.
Requires python-flint."""
import sys, os, json, hashlib
from fractions import Fraction as Fr
_H = os.path.dirname(os.path.abspath(__file__))
# reference bank root comes from the environment (no default); unset =>
# loud refusal (fail-closed).
TILE_BANK = os.environ.get('TERRIER_TILE_BANK', '')
if not TILE_BANK:
    raise SystemExit("REFUSE: env TERRIER_TILE_BANK unset — must point at the "
                     "tile-engine reference bank root (u1s card + envelope "
                     "table; not included in the package)")
U1S = os.path.join(TILE_BANK, 'u1s')
sys.path.insert(0, U1S)
sys.path.insert(0, TILE_BANK)
from flint import arb, acb, acb_mat, ctx
ctx.dps = 60
import mu2_pass as MP
import mu2_frame as MF          # sha-gated frame load (fail-closed)
import mu2_hodge as MH
import mu_a_probe as MA         # SIGMA sha-gated in-module

N = 29
C0 = [11, 13, 15, 17, 19, 21]   # tile-center base coordinates, units 1/1024
OFF, HNUM = 4, 4                # tile-center offset, half-width h=4/1024=1/256
PIN_KG6 = "f2ac7902b739f864"    # pinned sha16 of the K_G6 extreme-tile set
LADDER = (50, 80, 120, 170, 230, 300)


def sha16(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 16), b''):
            h.update(b)
    return h.hexdigest()[:16]


# ---------- tile law (pinned K_G6 extreme set; fixed serialization) --------
def sig_of(k):
    return [1 if (k >> i) & 1 else -1 for i in range(6)]


def center_exact(sig):
    return [[C0[i] + OFF * sig[i], 1024, OFF, 1024] for i in range(6)]


def tileset_check():
    ser = json.dumps([[k, sig_of(k), center_exact(sig_of(k))]
                      for k in range(64)], separators=(",", ":")).encode()
    got = hashlib.sha256(ser).hexdigest()[:16]
    assert got == PIN_KG6, "K_G6 tile-set sha mismatch — HALT"
    return got


# ---------- OL-4 (i): locus clearance (exact Fractions, max-norm chart) -----
def locus_record(k):
    """Special-locus clearance of closed tile T_k (frozen special set +
    rho-limiting sheets). The clearance threshold eps0 for this geometry is
    not fixed (UNFROZEN) => NO certified (i) verdict is possible; record =
    certified distances + honest state."""
    v = [C0[i] + OFF * sig_of(k)[i] for i in range(6)]
    # coincidence strata phi_i = phi_j: tile-min dist = max(0,|dv|-2h)/2
    co = []
    for i in range(6):
        for j in range(i + 1, 6):
            d = abs(v[i] - v[j])
            co.append((max(0, d - 2 * HNUM), i + 1, j + 1))
    co.sort()
    co_lo = Fr(co[0][0], 2 * 1024)          # tile-uniform lower bound
    co_ctr = Fr(min(abs(v[i] - v[j]) for i in range(6)
                    for j in range(i + 1, 6)), 2 * 1024)   # center distance
    # LCS sheets phi_i = 0: tile-min |phi_i| >= (Re_min, Im_min=0) exactly
    lcs_lo = Fr(min(vi - HNUM for vi in v), 1024)
    return dict(tile=k, sig="".join("+" if s > 0 else "-" for s in sig_of(k)),
                coincidence_dist_lo=str(co_lo), coincidence_ctr=str(co_ctr),
                straddles_coincidence=bool(co_lo == 0),
                worst_pair=("%d=%d" % (co[0][1], co[0][2])),
                lcs_dist_lo=str(lcs_lo), eps0="UNFROZEN",
                i_verdict="FAIL-STRADDLE" if co_lo == 0 else "PENDING-eps0")


# ---------- flux dictionary (u1s frame convention) --------------------------
SIG = MA.SIG                     # sha-gated integer Sigma (mu_a_probe load)
PAIRS = [(i, j) for i in range(1, 7) for j in range(i + 1, 7)]


def flux_abc(a, b, c, n):
    """G(a,b,c;n) = (0, n_I a, -(n_I+n_J) b, n_I (c-a), 0), sum n_I = 0."""
    assert sum(n) == 0
    G = [0] + [n[i] * a for i in range(6)]
    G += [-(n[i - 1] + n[j - 1]) * b for (i, j) in PAIRS]
    G += [n[i] * (c - a) for i in range(6)] + [0]
    return G


def q_exact(G):
    return sum(G[i] * SIG[i][j] * G[j] for i in range(N) for j in range(N))


def n2_exact(G):
    return sum(g * g for g in G)


# ---------- OL-5 B(T) law — period-side, NEVER enumerated -------------------
def bT_law(lam_lo, lam_up):
    """N2_ISD(T) = 60/lam(T) is a FORMULA; a number exists only when a
    certified lam(T) BOX floor exists. Trivial period-side branch:
    lam_lo > 60 => B(T) has no nonzero integral point (n2 >= 1)."""
    if lam_lo is None or lam_up is None:
        # Fail-closed on EITHER missing bound: a certified lam_lo with an
        # absent lam_up would otherwise reach arb(60) / None below. The
        # census path always passes (None, None) here.
        return dict(gate="CLOSED (lam(T) box floor uncertified)",
                    N2_formula="60/lam(T)",
                    N2_value=None, bT_empty=None)
    empty = bool(arb(lam_lo) > 60)
    return dict(gate="OPEN", N2_formula="60/lam(T)",
                N2_interval=[float((arb(60) / lam_up).lower()),
                             float((arb(60) / lam_lo).upper())],
                bT_empty=empty)


# ---------- U1-S gated point layer (frame gates MANDATORY per point) -------
def phis_of(prat):
    """prat: 6 x [re_num, re_den, im_num, im_den] exact rationals -> acb."""
    return [acb(arb(a) / b, arb(c) / d) for (a, b, c, d) in prat]


def _pi_only(phis, s2, M, cb, twopii):
    """Pi alone (weight-0 pass) — S6-swap gate reference leg."""
    T = MP.run_passes(phis, M, [(0,) * 6], s2)
    tv = [p.log() / twopii for p in phis]
    A = MF.contract(T[(0,) * 6], cb, acb(1) / twopii)
    return MF.E_at(tv) * acb_mat([[x] for x in A])


def upoint(name, prat, tgt=8):
    """One gated U1-S summation point: card passes + numeric frame gates
    (transversality 28 multisets + S6-swap + T_1-shift) + immersivity Gram
    det (OL-4 ii) + certified lam bracket at the point (OL-4 iii input)."""
    import time
    t0 = time.time()
    phis = phis_of(prat)
    s2 = MP.s2_of(phis)
    M = MP.pick_M(s2, LADDER, tgt)
    assert M is not None, "HALT: no ladder M certifies tails (K4-class)"
    T = MP.run_passes(phis, M, MF.weights28(), s2)
    twopii = acb(0, 2 * arb.pi())
    tv = [p.log() / twopii for p in phis]
    cb = 80 * acb(3).zeta() / twopii ** 3
    Pi, dPi, ddPi, E = MF.jets29(T, tv, cb, twopii)
    # frame gate 1: transversality (identity (ii)), 28 multisets, teeth required
    tg, sc = MH.transv_gates(Pi, dPi, ddPi)
    gt = {k: MH.gate(v, arb(sc)) for k, v in tg.items()}
    g_tv = all(c for (c, p, r) in gt.values())
    g_teeth = all(p for (c, p, r) in gt.values())
    # frame gate 2: S6-swap (identity (iii)): MSWAP.Pi(t) == Pi(t2,t1,...)
    scp = arb(max(float(abs(x).abs_upper()) for x in Pi))
    Psw = _pi_only([phis[1], phis[0]] + phis[2:], s2, M, cb, twopii)
    RPi = MF.SWm * acb_mat([[x] for x in Pi])
    gs6 = [MH.gate(RPi[i, 0] - Psw[i, 0], scp) for i in range(N)]
    # frame gate 3: T_1-shift (identity (i)): E(t+e_1) A == T_1 Pi, same A-ball
    A0 = MF.contract(T[(0,) * 6], cb, acb(1) / twopii)
    lhs = MF.E_at([tv[0] + 1] + tv[1:]) * acb_mat([[x] for x in A0])
    T1m = acb_mat([[acb(MF.TsI[0][i][j]) for j in range(N)] for i in range(N)])
    rhs = T1m * acb_mat([[x] for x in Pi])
    gsh = [MH.gate(lhs[i, 0] - rhs[i, 0], scp) for i in range(N)]
    # OL-4 (ii) at the point: 7x7 sesquilinear Gram det, ball-separated
    F3 = [Pi] + dPi
    ScF3 = MA.sig_mul([MA.conj_vec(u) for u in F3])
    G7 = acb_mat([[MA.pair(F3[a], ScF3[b]) for b in range(7)]
                  for a in range(7)])
    det7 = G7.det()
    d_lo, d_up = float(abs(det7).lower()), float(abs(det7).upper())
    # OL-2 input: certified TRUE-period Hodge floor at the point
    H, (g_hf, g_n40, s_or) = MH.hodge_true(Pi, dPi, ddPi)
    alo, rup, lhat, mrad, (g_re, g_sy) = MA.certified_floor(H)
    cert = (g_tv and g_teeth and all(c for c, p, r in gs6) and
            all(p for c, p, r in gs6) and all(c for c, p, r in gsh) and
            all(p for c, p, r in gsh) and g_n40 and g_re and g_sy and
            d_lo > 0 and alo is not None)
    rec = dict(pt=name, M=M, s2=float(s2.abs_upper()),
               b4_transv=bool(g_tv and g_teeth),
               b4_transv_worst=max(r for (c, p, r) in gt.values()),
               b4_s6=bool(all(c and p for c, p, r in gs6)),
               b4_s6_worst=max(r for c, p, r in gs6),
               b4_shift=bool(all(c and p for c, p, r in gsh)),
               b4_shift_worst=max(r for c, p, r in gsh),
               imm_det_lo=d_lo, imm_det_up=d_up,
               imm_sep=bool(d_lo > 0), lam_lo=alo, lam_up=rup,
               lam_float=lhat, H_maxrad=mrad, n40=bool(g_n40), s=s_or,
               CERT=bool(cert), s_pt=round(time.time() - t0, 1))
    return rec, Pi, dPi


# ---------- OL-4 (iv) primitive: per-flux 7-row certificate (period-side) ---
def rows_vs_flux(Pi, dPi, G, tol=arb(1) / 10 ** 5):
    """W=0 rows L_0 = Sigma Pi, L_K = Sigma d_K Pi paired with
    integer flux G. POINT-layer verdict per OL-4(iv): KILLED if some row ball
    is separated from 0; VACUOUS(FAIL) if containment has no teeth (arb-ball
    law); else CONTAINS-0 (no kill possible against this G)."""
    SG = [sum(SIG[a][b] * G[b] for b in range(N)) for a in range(N)]
    rows = [sum(Pi[a] * SG[a] for a in range(N))]
    rows += [sum(dPi[K][a] * SG[a] for a in range(N)) for K in range(6)]
    sc = max(float(abs(x).abs_upper()) for x in Pi) * max(
        1, max(abs(g) for g in SG)) * N
    out, killed = [], []
    for i, r in enumerate(rows):
        lo = float(abs(r).lower())
        rad = max(float(r.real.rad()), float(r.imag.rad()))
        cont = r.real.abs_lower() == 0 and r.imag.abs_lower() == 0
        teeth = arb(rad) < tol * arb(sc)
        out.append(dict(row=i, abs_lo=lo, rad=rad, contain0=bool(cont)))
        if lo > 0:
            killed.append(i)
    all0 = all(o["contain0"] for o in out)
    teeth_ok = all(arb(o["rad"]) < tol * arb(sc) for o in out)
    verdict = ("POINT-KILL" if killed else
               ("CONTAINS-0" if (all0 and teeth_ok) else
                ("VACUOUS-FAIL" if all0 else "MIXED-OPEN")))
    marg = max((o["abs_lo"] for o in out), default=0.0) / sc
    return dict(Q=q_exact(G), n2=n2_exact(G), rows=out, killed_rows=killed,
                verdict=verdict, kill_margin_rel=marg, scale=sc)


# ---------- reference corner/entropy data per tile (OL-2 lam inputs) --------
def reference_corner_bracket(k):
    """Certified TRUE-period lam bracket at corner_k — the outer real corner
    of extreme tile k (same bit law), read from the reference point logs
    u1s/MU2_POINTS_s*.log under TERRIER_TILE_BANK (every row must be CERT)."""
    tag = '"corner_%02d"' % k
    for s in range(3):
        p = os.path.join(U1S, "MU2_POINTS_s%d.log" % s)
        for ln in open(p):
            if ln.startswith("{") and tag in ln:
                r = json.loads(ln)
                assert r["CERT"], "reference corner row not CERT — HALT"
                return dict(lam_lo=r["lam_lo"], lam_up=r["lam_up"],
                            M=r["M"], src="MU2_POINTS_s%d.log" % s)
    raise RuntimeError("corner_%02d missing from MU2 logs" % k)


def reference_envelope_row(k):
    """Reference envelope row for tile k (rho-limiting-sheet clearance,
    necessary-condition layer), read from the envelope table under
    TERRIER_TILE_BANK."""
    for ln in open(os.path.join(TILE_BANK, "CORNER_TABLE.log")):
        p = ln.split()
        if len(p) == 9 and not ln.startswith("#") and p[0] == str(k):
            return dict(budget_lo=float(p[6]), limiter=p[7], verdict=p[8])
    raise RuntimeError("reference envelope row %d missing" % k)


def entropy_hull(k):
    """s^2 hull of tile k (corner-monotone majorant scale of the card)."""
    h = arb(1) / 256
    ph = [acb(arb(C0[i] + OFF * sig_of(k)[i]) / 1024 + arb(0, h),
              arb(OFF) / 1024 + arb(0, h)) for i in range(6)]
    s2 = MP.s2_of(ph)
    return float(s2.abs_upper())


# ---------- OL-4 verdict composer (fail-closed conjunction) -----------------
def tile_verdict(loc, ctr_rec, corner, env, s2h, fluxes):
    """EXCLUDED only if (i) AND (ii on T) AND (iii lam(T) box floor > 0) AND
    ((iv) all of B(T) killed tile-wide OR B(T) empty). Any missing certified
    bound => OPEN with named reasons; a control vacuum flux with rows
    CONTAINS-0 additionally pins NOT-EXCLUDABLE-AT-POINT."""
    reasons = []
    if loc["i_verdict"] != "PASS":
        reasons.append("(i) " + loc["i_verdict"] +
                       " [eps0 unfrozen; coincidence dist_lo=" +
                       loc["coincidence_dist_lo"] + "]")
    reasons.append("(ii) point-certified only (imm_sep=%s at center); "
                   "no tile-uniform certificate" % ctr_rec["imm_sep"])
    lam_box = None                       # no certified box floor available
    reasons.append("(iii) lam(T) box floor UNCERTIFIED; "
                   "POINT layer: center [%s, %s], corner [%s, %s]"
                   % (ctr_rec["lam_lo"], ctr_rec["lam_up"],
                      corner["lam_lo"], corner["lam_up"]))
    bT = bT_law(lam_box, None)
    reasons.append("(iv) per-flux POINT certificates only; no tile-wide "
                   "inf_T row bound available")
    excl = False                         # fail-closed: no tile-wide (iv)
    vac_pinned = any(f["verdict"] == "CONTAINS-0" for f in fluxes)
    return dict(outcome="EXCLUDED" if excl else "OPEN",
                not_excludable_at_point=bool(vac_pinned),
                bT=bT, s2_hull=s2h, envelope_nec=env, reasons=reasons)

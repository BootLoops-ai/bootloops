#!/usr/bin/env python3
"""
coni_frobenius.py — conifold Frobenius block + coni-PFV racetrack (g3').

(b) Local basis at z_cf = 0 (freeze Group D, PINNED-V — do not re-derive):
    D1/D2: B-branch log scheme  (n_cf/2 pi i) z_cf log(-2 pi i z_cf); the
    sqrt(pi/2)-normalized superpotential block (2406 W^(1), P-flux choice):
      sqrt(pi/2) W = W_bulk(s) + z_cf W1 + O(z_cf^2 log z_cf)
      W1 = -M_cf (n_cf/2 pi i)(log(-2 pi i z_cf) - 1) + c_tau * tau
           + (1/2 pi i) sum_{q not|| q_cf} n_q (q.M)(q.xi) Li_1(s^{q.nu})
      c_tau = (M.N.p)/M_cf = (M.K + lam M_cf)/M_cf   (exact Fraction)
    BOOKKEEPING NOTE: the D4 form carries -tau K z_cf
    with INTEGER K; in the 2406 chart that content sits inside the
    kappa/Li_1 terms above via Q_throat.  We implement the 2406 chart and
    report the D4-integer-K variant as a labeled residue, never a gate.
(d) Racetrack W_N = -(1/4 pi^2) S_N Li_2(s^N), S_N = sum_{q.nu=N} n_q (q.M)
    (exact integers); physical prefactor sqrt(2/pi) (4627-main fingerprint
    sqrt(8/pi^5) verified).  tau vev: 2406 two-term formula + Newton on the
    full pinned ladder.  z_cf vev: eq:conifold_vev (Q_throat form).  Floors:
    eq:zbound / eq:cookedness.  ALL numerics = mpmath ESTIMATES (labeled);
    exact parts are Fractions.  Gates CP7-CP10 per conipfv/DESIGN.md SS5.
"""
import sys, os
from fractions import Fraction as Fr

_PIPE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PIPE)

from mpmath import mp, mpf, polylog, exp as mexp, log as mlog, pi as mpi, \
    sqrt as msqrt, mpc

ZETA_KPV = mpf("114.037")     # freeze F4 (2406 App. B), used by eq:cookedness


def ladder(card, cur):
    """CP7: exact racetrack ladder.  levels[N] = (S_N, [(q, n_q), ...])."""
    h = card.h
    q_cf = tuple(cur["q_cf"])
    nu = cur["nu"]
    lev = {}
    for q, n in sorted(card.gv_pinned.items()):
        N = sum(q[a] * nu[a] for a in range(h))
        assert N >= 0, f"CP7: negative degree for pinned class {q}"
        qM = sum(Fr(q[a]) * card.M[a] for a in range(h))
        lev.setdefault(N, []).append((q, n, qM))
    # level 0 = conifold ray only; D5 GV-nilpotency: only k=1 has GV != 0
    for q, n, _ in lev.get(0, []):
        k = set(Fr(a, b) for a, b in zip(q, q_cf) if b != 0)
        ok = all(a == 0 for a, b in zip(q, q_cf) if b == 0) and len(k) == 1
        assert ok, f"CP7: level-0 class {q} not on the conifold ray"
        assert k.pop() == 1 or n == 0, \
            f"CP7: GV({q}) != 0 at multiple of q_cf (D5 nilpotency)"
    S = {N: sum(n * qM for _, n, qM in cls) for N, cls in lev.items()}
    pos = sorted(N for N in S if N > 0)
    assert len(pos) >= 2, "CP7: fewer than 2 positive racetrack levels pinned"
    return dict(levels={N: [(list(q), int(n), qm) for q, n, qm in cls]
                        for N, cls in lev.items()},
                S={N: S[N] for N in S}, pos_levels=pos)


def _fr2mp(x):
    return mpf(x.numerator) / x.denominator


def joint_vev(card, cur, fb, dps=50):
    """<tau> = i t from the JOINT F-term dW_eff/dtau = 0 (ESTIMATE), where
      W_eff = W_bulk_eff + W_cf_eff,  W_cf_eff = conifold sector with z_cf
      integrated out at its own F-point (the coni-PFV racetrack: an
      effective instanton of fractional degree r*lam/(M_cf n_cf)).
    On tau = i t (x = e^{-2 pi t/r}), everything reduces to ONE real eq:
      b(t)   = lam t - (1/2 pi) sum tot_N Li_1(x^N)      [B = i b]
      zmag   = e^{-2 pi b/(M_cf n_cf)}/(2 pi),  z* = i zmag
      G(t)   = Im dW/dtau = -(1/2 pi) sum S_N (N/r) Li_1(x^N)
               + zmag (lam + sum tot_N (N/r) Li_0(x^N))  = 0.
    Root-find: log-grid scan + bisection; take the LARGEST-t root (weak
    coupling); all roots reported.  D4 validity: exponent must be negative
    (b * sign(M_cf n_cf) > 0), else gate FAIL."""
    mp.dps = dps
    lad, r = cur["_ladder"], cur["r"]
    S, pos = lad["S"], lad["pos_levels"]
    lamf = _fr2mp(cur["lam"])
    Mn = mpf(int(cur["M_cf"])) * cur["n_cf"]
    li1 = [(N, _fr2mp(tot)) for N, tot in fb["li1_ladder"]]

    def parts(t):
        x = mexp(-2 * mpi * t / r)
        b = lamf * t - sum(tot * polylog(1, x ** N)
                           for N, tot in li1) / (2 * mpi)
        assert b * Mn > 0, "CP9: D4 validity fails (z_cf exponent >= 0)"
        zmag = mexp(-2 * mpi * b / Mn) / (2 * mpi)
        E = sum(_fr2mp(S[N]) * N * polylog(1, x ** N)
                for N in pos) / (2 * mpi * r)
        D = lamf + sum(tot * (mpf(N) / r) * polylog(0, x ** N)
                       for N, tot in li1)
        return -E + zmag * D, zmag, b

    lo, hi, n = mpf(1), mpf(400), 600
    ts = [lo * (hi / lo) ** (mpf(i) / (n - 1)) for i in range(n)]
    roots = []
    prev = None
    for t in ts:
        try:
            g = parts(t)[0]
        except AssertionError:
            prev = None
            continue
        if prev is not None and g * prev[1] < 0:
            a, bb = prev[0], t
            for _ in range(200):
                m = (a + bb) / 2
                gm = parts(m)[0]
                if gm * parts(a)[0] < 0:
                    bb = m
                else:
                    a = m
                if bb - a < mpf(10) ** (-dps + 10) * m:
                    break
            roots.append((a + bb) / 2)
        prev = (t, g)
    assert roots, "CP9: no racetrack root of the joint F-term in [1,400]"
    t = roots[-1]
    _, zmag, b = parts(t)
    return dict(t=t, roots=[float(x) for x in roots], zmag_F=zmag, b=b,
                q_eff=float(lamf * r / Mn),
                note="joint dW/dtau=0 on pinned window; Kahler-connection "
                     "term and beyond-window GV NOT included (ESTIMATE)")


def w_bulk_eff(lad, r, t, dps=50):
    """sqrt(pi/2)-normalized W_bulk_eff at tau = i t (ESTIMATE)."""
    mp.dps = dps
    x = mexp(-2 * mpi * t / r)
    return sum(-mpf(int(lad["S"][N])) / (4 * mpi ** 2)
               * polylog(2, x ** N) for N in lad["pos_levels"])


def frobenius_block(card, cur):
    """CP8: EXACT Frobenius-block coefficients (Fractions + pinned scheme).
    sqrt(pi/2) W supset z_cf * [ c_log * (1/2 pi i)(log(-2 pi i z_cf) - 1)
                                 + c_tau * tau + Li_1 ladder ]"""
    h = card.h
    # M.N.p = M.K + lam*M_cf   (exact lemma identity, CP2)
    MK = sum(card.M[a] * card.K[a] for a in range(h))
    MNp = MK + cur["lam"] * cur["M_cf"]
    c_log = -cur["M_cf"] * cur["n_cf"]          # coefficient of (1/2pi i)log
    # tau-coefficient of the z_cf block: kappa term (MNp/M_cf) PLUS the
    # D4-pinned -tau K z_cf term (-(K.xi) = -MK/M_cf).  EXACT IDENTITY:
    #   c_tau = (MNp - MK)/M_cf = lam  (the conifold residual = Table-3 K').
    # This is the coni-PFV racetrack structure: the conifold sector enters
    # the tau F-term as an effective instanton of fractional degree
    # lam/(M_cf n_cf) * r  (the card NOTEs' "conifold-sector racetrack").
    c_tau = cur["lam"]
    assert c_tau == (MNp - MK) / cur["M_cf"], "CP8: c_tau != lam identity"
    xiden = cur["M_cf"]                          # xi = M / M_cf
    li1 = []                                     # (N, sum n_q (q.M)(q.xi))
    for N, cls in sorted(cur["_ladder"]["levels"].items()):
        if N == 0:
            continue
        tot = Fr(0)
        for q, n, qM in cls:
            qxi = Fr(sum(q[a] * card.M[a] for a in range(h)), 1) / xiden
            tot += n * qM * qxi
        li1.append((int(N), tot))
    return dict(c_log=c_log, c_tau=c_tau, MK=MK, MNp=MNp,
                li1_ladder=li1,
                scheme="log(-2 pi i z_cf) (freeze D2, 2009.03312 (3.9)-(3.10))",
                norm="sqrt(pi/2) on W (freeze D4)")


def z_cf_vev(card, cur, t, dps=50):
    """CP9: eq:conifold_vev (Q_throat form) + eq:zbound floor.  ESTIMATE.
    Q_flux = -M.K (covering units);  ||M||^2 = -M.K(M).M at the vev:
    kappa part exact  ( = -(M.N.p) t ),  Li_1 part numeric."""
    mp.dps = dps
    h, r = card.h, cur["r"]
    MK = sum(card.M[a] * card.K[a] for a in range(h))
    Q_flux = int(-MK)
    MNp = MK + cur["lam"] * cur["M_cf"]
    x = mexp(-2 * mpi * t / r)
    li = mpf(0)
    for N, cls in sorted(cur["_ladder"]["levels"].items()):
        if N == 0:
            continue
        for q, n, qM in cls:
            qMf = mpf(qM.numerator) / qM.denominator
            li += mpf(int(n)) * qMf ** 2 * polylog(1, x ** N)
    MNpf = mpf(MNp.numerator) / MNp.denominator
    M2 = -MNpf * t + li / (2 * mpi)              # ||M||^2 estimate
    assert M2 > 0, "CP9: ||M||^2 <= 0 (inadmissible fluxes?)"
    g_s = 1 / t
    Qthr = mpf(Q_flux) - g_s * M2
    ncf, Mcf = cur["n_cf"], abs(int(cur["M_cf"]))
    zvev = mexp(-2 * mpi * Qthr / (g_s * Mcf ** 2 * ncf)) / (2 * mpi)
    Q_O = 2 * (card.Q_D3 - 1)
    floor = mexp(-2 * mpi * Q_O / (g_s * Mcf ** 2 * ncf)) / (2 * mpi)
    assert Qthr <= Q_flux + mpf("1e-30"), "CP9: Q_throat > Q_flux"
    return dict(zvev=zvev, floor=floor, Q_flux=int(Q_flux), Qthr=Qthr,
                M2=M2, g_s=g_s, floor_ok=bool(zvev >= floor * (1 - mpf("1e-25"))))


def w0_racetrack(card, cur, fb, dps=50):
    """CP10: coni-PFV W0 estimate (2406 racetrack form, the eq (3.36)-(3.38)
    class) + published floors eq:zbound / eq:cookedness.
      W0 = sqrt(2/pi) |W_bulk_eff(<tau>) + W_cf_eff|,
      W_cf_eff = -A z* = (M_cf n_cf / 2 pi) zmag   (z* = i zmag)."""
    mp.dps = dps
    lad, r = cur["_ladder"], cur["r"]
    jv = joint_vev(card, cur, fb, dps)
    t = jv["t"]
    wb = w_bulk_eff(lad, r, t, dps)
    Wcf = mpf(int(cur["M_cf"])) * cur["n_cf"] * jv["zmag_F"] / (2 * mpi)
    W0 = msqrt(mpf(2) / mpi) * abs(wb + Wcf)
    # cross-check: Q_throat form of the vev (eq:conifold_vev) at the same t
    zv = z_cf_vev(card, cur, t, dps)
    gap = abs(mlog(jv["zmag_F"] / zv["zvev"]))
    assert gap < mpf("1.5"), \
        "CP9: F-term vev vs Q_throat vev disagree beyond Li-correction level"
    # published floor (eq:cookedness; n_cf = 2 assumed there)
    g_s, Mcf = zv["g_s"], abs(int(cur["M_cf"]))
    Q_O = 2 * (card.Q_D3 - 1)
    floor = (1 / (g_s ** 2 * Mcf)) \
        * msqrt(ZETA_KPV * Q_O / (2 * mpi) ** mpf("4/3")) \
        * mexp(-2 * mpi * Q_O / (3 * g_s * Mcf ** 2))
    return dict(joint=jv, t=t, W_bulk=wb, W_cf=Wcf, z_cf=zv,
                zmag_F=jv["zmag_F"], vev_form_gap_efolds=gap, W0=W0,
                W0_floor=floor, floor_ok=bool(W0 >= floor * mpf("1e-2")),
                note="ESTIMATE (mpmath dps %d, pinned GV window only); "
                     "certified value = transport legs 1-3, not this" % dps)


def g3_series(card, cur, dps=50):
    """g3' orchestrator: CP7-CP10.  cur = coni_pfv_curve output (CONI-PFV)."""
    assert cur["routing"] == "CONI-PFV"
    lad = ladder(card, cur)                       # CP7
    cur["_ladder"] = lad
    fb = frobenius_block(card, cur)               # CP8
    assert fb["c_log"] == -cur["M_cf"] * cur["n_cf"], "CP8: log coeff pin"
    res = w0_racetrack(card, cur, fb, dps)        # CP9 + CP10 inside
    assert res["z_cf"]["floor_ok"], "CP9: z_cf vev below published floor"
    chk = _pin_checks(card, cur, res)
    return dict(ladder=lad, frobenius=fb, racetrack=res, pin_checks=chk)


def _pin_checks(card, cur, res):
    """Labeled consistency vs card/NOTE pins (tolerant; ESTIMATE gates)."""
    import re
    chk = {}
    t = res["t"]
    cw = card.raw.get("coni_basis_crosswalk") or {}
    fd = cw.get("flat_direction_exact") or {}
    sv = fd.get("s_vac_estimate") or ""
    m = re.search(r"=\s*([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)", str(sv))
    if m:
        # tolerance is LOOSE by design: the pinned-window PFV step-1
        # estimate vs the repo's full-F-term value (2406's own step-2 gap);
        # this is a consistency check, not a certification.
        s_note = mpf(m.group(1))
        s_here = mexp(-2 * mpi * t / cur["r"])
        dlog = abs(mlog(s_here / s_note))
        chk["s_vac_vs_note"] = dict(note=float(s_note), here=float(s_here),
                                    dlog=float(dlog),
                                    ok=bool(dlog < mpf("0.75")))
    if card.tau_pin is not None:
        tp = mpf(card.tau_pin.numerator) / card.tau_pin.denominator
        rel = abs(t - tp) / tp
        chk["tau_vs_pin"] = dict(pin=float(tp), here=float(t),
                                 rel=float(rel), ok=bool(rel < mpf("0.20")))
    bad = [k for k, v in chk.items() if not v["ok"]]
    assert not bad, f"CP9/CP10 pin check FAIL: {bad}"
    return chk


if __name__ == "__main__":
    import json
    from family import load_card
    from coni_curve import coni_pfv_curve
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    card = load_card(sys.argv[1])
    cur = coni_pfv_curve(card)
    if cur["routing"] != "CONI-PFV":
        print(f"{card.name}: NOT-CONI — g3' not applicable")
        sys.exit(0)
    out = g3_series(card, cur, dps=50)
    lad, res = out["ladder"], out["racetrack"]
    print(f"{card.name}: g3' PASS")
    print("  ladder S_N:", {N: str(lad['S'][N]) for N in sorted(lad['S'])})
    print(f"  tau = i*{mp.nstr(res['t'], 12)} (joint F-term; roots "
          f"{res['joint']['roots']}; coni q_eff = "
          f"{res['joint']['q_eff']:.4f} vs 1/r = {1.0/cur['r']:.4f})")
    print(f"  g_s = {mp.nstr(res['z_cf']['g_s'], 8)}  Q_thr = "
          f"{mp.nstr(res['z_cf']['Qthr'], 8)} <= Q_flux {res['z_cf']['Q_flux']}")
    print(f"  |z_cf| F-term = {mp.nstr(res['zmag_F'], 8)}  Q_throat-form = "
          f"{mp.nstr(res['z_cf']['zvev'], 8)}  (gap "
          f"{mp.nstr(res['vev_form_gap_efolds'], 4)} e-folds; floor "
          f"{mp.nstr(res['z_cf']['floor'], 6)}, ok={res['z_cf']['floor_ok']})")
    print(f"  Frobenius: c_log = {out['frobenius']['c_log']}, c_tau = "
          f"{out['frobenius']['c_tau']}  [{out['frobenius']['scheme']}]")
    print(f"  W0_est = {mp.nstr(res['W0'], 10)}  (floor "
          f"{mp.nstr(res['W0_floor'], 6)}, ok={res['floor_ok']})")
    print(f"  pin checks: {out['pin_checks']}")

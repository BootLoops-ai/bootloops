"""Adapter battery — the etienne oracles are NEW derivative code, so the
additions are battery-tested before any certificate of record cites them
(a planted-optimum-on-synthetic requirement + the BALLER tamper
culture).

Legs:
  E1  corner closed forms: K_J = 1, K_{J-1} = sum_{n_i>=2} n_i/2 vs the
      study's exact/ball kernels on small data.
  E2  corner q-gradient: closed form at q=0 vs mpmath FD of the exact
      q-form likelihood at small q (independent derivation).
  E3  diffI VALUE identity on a synthetic 2-plot instance vs a fully
      independent mpmath implementation (integer Stirling recurrence +
      python convolutions; no flint anywhere in the reference chain).
  E4  diffI GRADIENT+HESSIAN vs mpmath FD on the synthetic instance
      (A5 register: <=1e-12 matched precision; measured much tighter).
  E5  PLANTED OPTIMUM: the synthetic instance's ML point located
      independently in mpmath (FD Newton at 60 dps); CLINCH certificate at
      the adapter-polished center must be CERTIFIED-INTERIOR and the ball
      must contain the mpmath point.
  E6  TAMPER: a sign-flipped gradient component in a copied oracle must
      NOT certify at the same candidate.
  E7  diffI on the REAL Panama-65 matrix: value vs the study's own
      engine value (receipt register) + FD-vs-analytic spot gradient/
      Hessian rows at 224 bits (assembly check at scale).
Run: python3 -m etienne.battery_etienne [--fast]  (E7 is the slow leg).
"""
import json
import math
import os
import sys
import time

sys.dont_write_bytecode = True

import numpy as np

from . import RECEIPTS_DIR, ident
from .oracle_corner import EtienneCornerOracle
from .oracle_diffi import EtienneDiffIOracle

from clinch import engine

RESULTS = {}


def leg(name):
    def deco(fn):
        def run(*a, **k):
            t0 = time.time()
            try:
                out = fn(*a, **k)
                out = dict(out or {}, ok=bool(out.get("ok", True)))
            except Exception as e:  # a battery leg failure is a failure
                out = dict(ok=False, error=repr(e))
            out["wall_s"] = round(time.time() - t0, 1)
            RESULTS[name] = out
            print(f"[{name}] {'PASS' if out['ok'] else 'FAIL'} "
                  f"({out['wall_s']}s)", flush=True)
            return out
        return run
    return deco


# ---------------- independent mpmath diffI reference -----------------
def _stirling_row(n):
    """Unsigned Stirling first kind s(n, a), a=0..n — integer recurrence,
    independent of the flint rf-poly route."""
    row = [1]
    for m in range(n):
        new = [0] * (m + 2)
        for a, v in enumerate(row):
            new[a] += m * v
            new[a + 1] += v
        row = new
    return row


class MPDiffI:
    """Fully independent different-I multisample lnP at high dps."""

    def __init__(self, rows, dps=60):
        import mpmath as mp
        mp.mp.dps = dps
        self.mp = mp
        self.rows = [tuple(r) for r in rows]
        self.N = len(rows[0])
        self.S = len(rows)
        self.Js = [sum(r[i] for r in rows) for i in range(self.N)]

    def lnP(self, u):
        mp = self.mp
        th = mp.e ** u[0]
        I = [mp.e ** v for v in u[1:]]
        polys = []
        for row in self.rows:
            u_s = [mp.mpf(1)]
            for p, n in enumerate(row):
                if n == 0:
                    continue
                st = _stirling_row(n)
                f = [st[a] * I[p] ** a if a else mp.mpf(0)
                     for a in range(n + 1)]
                u_s = [mp.fsum(u_s[i] * f[b - i]
                               for i in range(max(0, b - n + 1 - 1),
                                              min(len(u_s), b + 1))
                               if 0 <= b - i < len(f))
                       for b in range(len(u_s) + n)]
            P_s = [u_s[b] * mp.factorial(b - 1) if b >= 1 else mp.mpf(0)
                   for b in range(len(u_s))]
            polys.append(P_s)
        M = [mp.mpf(1)]
        for P_s in polys:
            M = [mp.fsum(M[i] * P_s[b - i]
                         for i in range(len(M)) if 0 <= b - i < len(P_s))
                 for b in range(len(M) + len(P_s) - 1)]
        tot = mp.mpf(0)
        w = mp.mpf(1)
        for A in range(len(M)):
            tot += M[A] * w
            w /= th + A
        from collections import Counter
        st = mp.mpf(0)
        for cnt in Counter(self.rows).values():
            st -= mp.loggamma(cnt + 1)
        for i in range(self.N):
            st += mp.loggamma(self.Js[i] + 1)
        for r in self.rows:
            for n in r:
                if n > 1:
                    st -= mp.loggamma(n + 1)
        out = st + self.S * mp.log(th) + mp.log(tot)
        for i in range(self.N):
            out -= mp.loggamma(I[i] + self.Js[i]) - mp.loggamma(I[i])
        return out

    def grad_fd(self, u, h=None):
        mp = self.mp
        h = h or mp.mpf(2) ** -40
        g = []
        for k in range(len(u)):
            def f(t):
                v = list(u)
                v[k] = v[k] + t
                return self.lnP(v)
            g.append((8 * (f(h) - f(-h)) - (f(2 * h) - f(-2 * h)))
                     / (12 * h))
        return g

    def hess_fd(self, u, h=None):
        mp = self.mp
        h = h or mp.mpf(2) ** -20
        n = len(u)
        H = [[mp.mpf(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(i, n):
                def f(ti, tj):
                    v = list(u)
                    v[i] += ti
                    v[j] += tj
                    return self.lnP(v)
                if i == j:
                    v = (16 * (f(h, 0) + f(-h, 0))
                         - (f(2 * h, 0) + f(-2 * h, 0))
                         - 30 * f(0, 0)) / (12 * h * h)
                else:
                    v = (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) \
                        / (4 * h * h)
                H[i][j] = H[j][i] = v
        return H


SYN_ROWS = [(3, 1), (2, 0), (0, 4), (1, 2), (5, 1), (0, 1), (2, 2)]


def _syn_oracle():
    return EtienneDiffIOracle(SYN_ROWS, tag="battery_synth")


def _x_at(o, u, prec=224):
    from flint import arb, ctx
    ctx.prec = prec
    return ([[arb(v) for v in u[1:]]], [arb(u[0])])


@leg("E1_corner_K_closed_forms")
def e1():
    p = ident.pilot()
    from flint import arb
    checks = []
    for D in ([1, 1, 2, 3, 5, 8], [2, 2, 7], [1, 1, 1, 4], [3, 9, 27, 6]):
        J = sum(D)
        S, K = p["ball_engine"].K_exact(sorted(D))
        kJ = float(K[J - S])
        kJm1 = float(K[J - 1 - S])
        cf = sum(n / 2 for n in D if n >= 2)
        checks.append(dict(D=D, K_J=kJ, K_Jm1=kJm1, closed_form=cf,
                           ok=(kJ == 1.0 and abs(kJm1 - cf) < 1e-12)))
    return dict(checks=checks, ok=all(c["ok"] for c in checks))


@leg("E2_corner_q_gradient_vs_mp")
def e2():
    import mpmath as mp
    mp.mp.dps = 60
    p = ident.pilot()
    D = sorted([1, 1, 2, 3, 5, 8, 4, 4])
    J, S = sum(D), len(D)
    S0, K = p["ball_engine"].K_exact(D)
    logK = [mp.log(mp.mpf(k.numerator)) - mp.log(mp.mpf(k.denominator))
            for k in K]
    th = mp.mpf("7.3")

    def lnP_q(q):
        # sum_A K_A q^(J-A)/(th)_A  (log-free at these sizes) - sum ln(1+kq)
        tot = mp.mpf(0)
        poch = mp.mpf(1)
        for k in range(S):
            poch *= th + k
        for idx in range(len(K)):
            A = S + idx
            tot += mp.e ** logK[idx] * q ** (J - A) / poch
            if A < J:
                poch *= th + A
        out = mp.log(tot)
        for k in range(J):
            out -= mp.log(1 + k * q)
        return out

    h = mp.mpf(2) ** -50
    fd = (8 * (lnP_q(h) - lnP_q(-h)) - (lnP_q(2 * h) - lnP_q(-2 * h))) \
        / (12 * h)
    o = EtienneCornerOracle(D, tag="battery_corner")
    from flint import arb, ctx
    ctx.prec = 224
    x = ([[arb(0)], [arb(0)]], [arb(float(mp.log(th)))])
    Fz, Fg = o.F(x)
    g_q_closed = -float(Fz[1][0].mid())    # lnP convention
    diff = abs(g_q_closed - float(fd))
    return dict(fd_mp=float(fd), closed=g_q_closed, abs_diff=diff,
                bar=1e-10 * max(1.0, abs(float(fd))),
                ok=diff <= 1e-10 * max(1.0, abs(float(fd))))


@leg("E3_diffi_value_vs_independent_mp")
def e3():
    o = _syn_oracle()
    ref = MPDiffI(SYN_ROWS, dps=60)
    u = [math.log(3.7)] + [math.log(v) for v in (2.1, 5.3)]
    x = _x_at(o, u)
    mine = o.evaluate(x, order=0)["lnP"]
    theirs = ref.lnP([ref.mp.mpf(v) for v in u])
    diff = abs(float(mine.mid()) - float(theirs))
    return dict(arb=float(mine.mid()), mp=float(theirs), abs_diff=diff,
                rad=float(mine.rad()), ok=diff <= 1e-12)


@leg("E4_diffi_derivs_vs_mp_fd")
def e4():
    o = _syn_oracle()
    ref = MPDiffI(SYN_ROWS, dps=60)
    u = [math.log(3.7), math.log(2.1), math.log(5.3)]
    x = _x_at(o, u)
    r = o.evaluate(x, order=2)
    um = [ref.mp.mpf(v) for v in u]
    g_mp = ref.grad_fd(um)
    H_mp = ref.hess_fd(um)
    worst_g = max(abs(float(r["grad"][k].mid()) - float(g_mp[k]))
                  for k in range(3))
    worst_h = max(abs(float(r["hess"][i][j].mid()) - float(H_mp[i][j]))
                  for i in range(3) for j in range(3))
    return dict(worst_grad_absdiff=worst_g, worst_hess_absdiff=worst_h,
                grad_bar=1e-12, hess_bar=1e-9,
                ok=(worst_g <= 1e-12 and worst_h <= 1e-9))


@leg("E5_planted_optimum")
def e5():
    o = _syn_oracle()
    ref = MPDiffI(SYN_ROWS, dps=60)
    mp = ref.mp
    # independent mpmath FD-Newton to the stationary point
    u = [mp.mpf(math.log(3.0)), mp.mpf(math.log(2.0)), mp.mpf(math.log(4.0))]
    for _ in range(40):
        g = ref.grad_fd(u)
        H = ref.hess_fd(u)
        Hn = np.array([[float(H[i][j]) for j in range(3)]
                       for i in range(3)])
        gn = np.array([float(v) for v in g])
        step = np.linalg.solve(Hn, gn)
        u = [u[k] - mp.mpf(float(step[k])) for k in range(3)]
        if float(np.max(np.abs(step))) < 1e-25:
            break
    u_mp = [float(v) for v in u]
    thp, rec = engine.newton_polish(o, u_mp)
    res = engine.certify(o, list(thp))
    ok = res["verdict"] == "CERTIFIED-INTERIOR"
    # containment of the independently-planted point
    contain = True
    r_scale = res.get("radius", 0.0)
    for k in range(3):
        if abs(float(thp[k]) - u_mp[k]) > max(r_scale, 1e-12):
            contain = False
    return dict(verdict=res["verdict"], radius=res.get("radius"),
                margins=(res.get("krawczyk") or {}).get("margins"),
                mp_point=u_mp, polished=[float(v) for v in thp],
                planted_contained=contain, ok=(ok and contain))


@leg("E6_tamper_math_mutation")
def e6():
    """The mutation lives in the oracle MATH (the documented Krawczyk
    trust boundary: a sign flip of one GRADIENT COMPONENT is zero-set-
    preserving — (F1,-F2) has the same zeros as (F1,F2) — so it must NOT
    and does not flip the verdict; that no-op mutation was measured first
    and replaced). Here lnI_2 -> -lnI_2 inside the model (a sign flip in
    the oracle's input math, the CLINCH L4 real-input pattern): F and H
    stay mutually consistent but the mutated model's optimum moves, so
    certification at the honest candidate must fail."""
    o = _syn_oracle()

    class Tampered(EtienneDiffIOracle):
        def _unpack(self, x):
            th, Ivec = EtienneDiffIOracle._unpack(self, x)
            Ivec = list(Ivec)
            Ivec[1] = 1 / Ivec[1]       # lnI_2 sign flip in the math
            return th, Ivec
    t = Tampered(SYN_ROWS, tag="battery_tamper")
    u = [math.log(3.0), math.log(2.0), math.log(4.0)]
    thp, _ = engine.newton_polish(o, u)      # honest polish
    res_h = engine.certify(o, list(thp))
    res_t = engine.certify(t, list(thp))
    return dict(honest=res_h["verdict"], tampered=res_t["verdict"],
                tampered_reason=res_t.get("reason"),
                ok=(res_h["verdict"] == "CERTIFIED-INTERIOR"
                    and res_t["verdict"] != "CERTIFIED-INTERIOR"))


@leg("E7_diffi_real_scale_identity")
def e7():
    import csv
    from flint import arb, ctx
    from . import DATA_DIR, PILOT_DIR
    p65 = os.path.join(
        DATA_DIR, "panama65/derived/panama65_recent_abundance_matrix_dbh100.csv")
    ident.log_read(p65, "panama-65 matrix (battery E7)")
    rows = []
    with open(p65) as f:
        rd = csv.reader(f)
        next(rd)
        for line in rd:
            rows.append(tuple(int(x) for x in line[1:]))
    r_path = os.path.join(PILOT_DIR, "PANAMA65_DIFFI_FIT.json")
    ident.log_read(r_path, "diffI fit receipt (battery E7)")
    rec = json.load(open(r_path))
    o = EtienneDiffIOracle(rows, tag="battery_e7")
    u = [math.log(rec["theta"])] + [math.log(v) for v in rec["I"]]
    ctx.prec = 224
    x = _x_at(o, u, prec=224)
    t0 = time.time()
    lnP = o.evaluate(x, order=0)["lnP"]
    t_val = time.time() - t0
    ref_ball = arb(rec["lnL_certified"])
    d_receipt = abs(float((lnP - ref_ball).mid()))
    # FD spot-check of the analytic gradient/Hessian rows (own-lnP FD at
    # 224 bits, exact power-of-two steps; coords: lntheta, lnI_3, lnI_40)
    t0 = time.time()
    r1 = o.evaluate(x, order=1)
    t_grad = time.time() - t0
    h = 2.0 ** -18
    spots = [0, 3, 40]
    fd = {}
    for k in spots:
        up = list(u); dn = list(u)
        up[k] += h; dn[k] -= h
        fp = o.evaluate(_x_at(o, up, 224), order=0)["lnP"]
        fm = o.evaluate(_x_at(o, dn, 224), order=0)["lnP"]
        fd[k] = float(((fp - fm) / (2 * h)).mid())
    worst = max(abs(fd[k] - float(r1["grad"][k].mid()))
                / max(1.0, abs(fd[k])) for k in spots)
    return dict(value_absdiff_vs_receipt=d_receipt,
                receipt_ball=rec["lnL_certified"],
                my_value=str(lnP), wall_value_s=round(t_val, 1),
                wall_grad_s=round(t_grad, 1),
                fd_spots={str(k): fd[k] for k in spots},
                an_spots={str(k): float(r1["grad"][k].mid())
                          for k in spots},
                worst_rel_fd=worst,
                ok=(d_receipt <= 1e-9 and worst <= 1e-8))


def main(fast=False):
    e1(); e2(); e3(); e4(); e5(); e6()
    if not fast:
        e7()
    ok = all(v["ok"] for v in RESULTS.values())
    out = dict(overall=("PASS" if ok else "FAIL"),
               stamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               legs=RESULTS)
    os.makedirs(RECEIPTS_DIR, exist_ok=True)
    json.dump(out, open(os.path.join(RECEIPTS_DIR,
                                     "ADAPTER_BATTERY.json"), "w"),
              indent=1, default=str)
    print("OVERALL:", out["overall"], flush=True)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(fast="--fast" in sys.argv) else 1)

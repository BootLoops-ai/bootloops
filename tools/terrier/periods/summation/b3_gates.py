#!/usr/bin/env python3
"""b3_gates.py — B3 battery gates G1-G4 (fail-closed; any miss => SystemExit).
G1 exact 0-diff: separable GF pass == independent direct-lattice exact pass
   (generic 6-var jet-exp path), both modes, complex + real points.
G2 harmonic-law vs mpmath Gamma (dps 60): order-5 residual scaling (the jet law
   is EXACT to deg 4 => residual ~ |rho|^5; two-scale ratio must be ~2^5).
G3 S6 equivariance: towers at permuted phi == permuted towers (exact).
G4 ball containment: arb/acb GF ball (prec 200) contains the exact GF value."""
import sys, time
from fractions import Fraction as Fr
import towers6 as tw
from towers6 import CQ, Fr as _F

C0 = (11, 13, 15, 17, 19, 21)
PT_C = [CQ(Fr(c + 4, 1024), Fr(4, 1024)) for c in C0]     # worst extreme tile ctr
PT_R = [CQ(Fr(c + 8, 1024)) for c in C0]                  # worst Stage-A corner


def gf_exact(phis, M, mode, Ht):
    g = tw.g_tables_exact(phis, M, mode, Ht)
    w = tw.contraction_weights_exact(M, mode, Ht)
    wts = [[(CQ(x) if x else None) for x in row] for row in w]
    gcq = [[{m: c for m, c in cell.items()} for cell in col] for col in g]
    return tw.gf_eval(gcq, wts, M, tw.rpow_expansions())


def main():
    t0 = time.time()
    Ht = tw.harmonic_tables(400)
    # --- G1: exact equality GF vs lattice ---
    for (pt, M, tag) in ((PT_C, 6, "cplx"), (PT_R, 8, "real")):
        for mode in ("rho", "theta"):
            A = gf_exact(pt, M, mode, Ht)
            B = tw.lattice_eval_exact(pt, M, mode, Ht)
            ka = {k for k, v in A.items() if v.re != 0 or v.im != 0}
            kb = {k for k, v in B.items() if v.re != 0 or v.im != 0}
            A = {k: A[k] for k in ka}
            assert ka == kb, f"G1 support mismatch {tag}/{mode}"
            for m in ka:
                assert A[m] == B[m], f"G1 FAIL {tag}/{mode} mono {m}"
            print(f"[G1] GF == lattice EXACT, {len(ka)} comps, {tag}/{mode} "
                  f"M={M} PASS ({time.time()-t0:.1f}s)", flush=True)
    # --- G2: harmonic jet law vs mpmath Gamma, order-5 residual scaling ---
    from mpmath import mp, mpf, gamma, mpc
    mp.dps = 60
    RP = tw.rpow_expansions()
    for k in ((0, 1, 2, 0, 3, 1), (7, 3, 11, 0, 5, 9), (40, 1, 2, 8, 0, 12)):
        n = sum(k)
        Q = {}
        for m in range(1, 5):
            cm = 2 * tw.s_m(m) * Ht[m][n]
            for mo, mult in RP[m].items():
                Q[mo] = Q.get(mo, Fr(0)) + cm * mult
            for i in range(6):
                Q[m * tw.PB[i]] = (Q.get(m * tw.PB[i], Fr(0))
                                   - 2 * tw.s_m(m) * Ht[m][k[i]])
        E, pw = {0: Fr(1)}, {0: Fr(1)}
        for j in range(1, 5):
            pw = tw.jetmul(pw, Q)
            pw = {mo: c / j for mo, c in pw.items()}
            for mo, c in pw.items():
                E[mo] = E.get(mo, Fr(0)) + c
        wdir = (1, -2, 1, 3, -1, 2)          # fixed rational direction
        res = []
        for eps in (mpf(1) / 1000, mpf(1) / 2000):
            rho = [eps * w for w in wdir]
            r = sum(rho)
            Rj = mpf(0)
            for mo, c in E.items():
                u = tw.unpack(mo)
                t = mpf(c.numerator) / c.denominator
                for i in range(6):
                    t *= rho[i] ** u[i]
                Rj += t
            Rg = (gamma(1 + n + r) ** 2 / gamma(1 + r) ** 2)
            for i in range(6):
                Rg *= gamma(1 + rho[i]) ** 2 / gamma(1 + k[i] + rho[i]) ** 2
            from math import factorial
            cd = Fr(factorial(n) ** 2)
            for x in k:
                cd /= factorial(x) ** 2
            Rg /= mpf(cd.numerator) / cd.denominator
            res.append(abs(Rg - Rj))
        ratio = res[0] / res[1]
        assert res[1] < mpf("1e-6"), f"G2 abs residual too big k={k}"
        assert 27 < ratio < 37, f"G2 order-5 ratio {ratio} k={k}"
        print(f"[G2] Gamma-vs-jet k={k}: res(1e-3)={float(res[0]):.3e} "
              f"ratio={float(ratio):.2f} (~32=2^5) PASS", flush=True)
    # --- G3: S6 equivariance (exact, M=6, complex point) ---
    sig = (2, 0, 1, 4, 5, 3)                 # image positions
    pt2 = [None] * 6
    for i in range(6):
        pt2[sig[i]] = PT_C[i]
    for mode in ("rho", "theta"):
        A = gf_exact(PT_C, 6, mode, Ht)
        B = gf_exact(pt2, 6, mode, Ht)
        for m, v in A.items():
            u = tw.unpack(m)
            m2 = sum(u[i] * tw.PB[sig[i]] for i in range(6))
            assert B[m2] == v, f"G3 FAIL {mode} mono {m}"
        print(f"[G3] S6 equivariance {mode} M=6 PASS", flush=True)
    # --- G4: ball containment, arb (real pt) + acb (complex pt), M=20 ---
    from flint import arb, acb, ctx
    ctx.prec = 200
    for (pt, tag) in ((PT_R, "arb-real"), (PT_C, "acb-cplx")):
        cv = (tw.conv_arb_maker(arb) if tag == "arb-real"
              else tw.conv_acb_maker(arb, acb))
        g = tw.g_tables_exact(pt, 20, "rho", Ht)
        w = tw.contraction_weights_exact(20, "rho", Ht)
        gr, wr = tw.ring_tables(g, w, cv)
        TB = tw.gf_eval(gr, wr, 20, tw.rpow_expansions())
        TE = gf_exact(pt, 20, "rho", Ht)
        for m, v in TE.items():
            assert TB[m].contains(cv(v)), f"G4 containment FAIL {tag} {m}"
        print(f"[G4] {tag} M=20 ball contains exact, {len(TE)} comps PASS",
              flush=True)
    print(f"[gates] ALL PASS wall {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()

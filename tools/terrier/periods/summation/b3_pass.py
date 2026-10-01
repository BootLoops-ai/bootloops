#!/usr/bin/env python3
"""b3_pass.py — B3 payload pass: certified tower balls (rho + theta modes,
210 comps each) at the 65 Stage-A points (mu_a convention: center + 64 real
box corners (c0 +- 8)/1024) + the worst complex extreme tile center (demo).
M=50 (B3 pricing line). flint balls prec 200 (dps ~60), tails by the
entropy-majorant law (towers6.tails_at), radii = ball.rad + tail (outward).
Timed-pilot: center-first M-ladder probe (M=10,20) + extrapolation printed BEFORE
the batch; per-point checkpoint lines (2h law); fail-closed on tail hypotheses.
Usage: b3_pass.py OUT.json [M] [deep]   (deep: adds M=150 PT_R + M=85 PT_C)"""
import json, sys, time
from fractions import Fraction as Fr
from flint import arb, acb, ctx
import towers6 as tw
from towers6 import CQ

ctx.prec = 200
C0 = (11, 13, 15, 17, 19, 21)
RPOW = tw.rpow_expansions()


def points():
    pts = [("center", [CQ(Fr(c, 1024)) for c in C0], True)]
    for kbit in range(64):
        phis = [CQ(Fr(C0[i] + (8 if (kbit >> i) & 1 else -8), 1024))
                for i in range(6)]
        pts.append((f"corner_{kbit:02d}", phis, True))
    pts.append(("tilectr_worst",
                [CQ(Fr(c + 4, 1024), Fr(4, 1024)) for c in C0], False))
    return pts


def eval_point(phis, real, M, Ht, mode):
    cv = tw.conv_arb_maker(arb) if real else tw.conv_acb_maker(arb, acb)
    g = tw.g_tables_exact(phis, M, mode, Ht)
    w = tw.contraction_weights_exact(M, mode, Ht)
    gr, wr = tw.ring_tables(g, w, cv)
    return tw.gf_eval(gr, wr, M, RPOW)


def s2_ball(phis):
    s = arb(0)
    for p in phis:
        q = tw.sq_modulus(p)
        s += (arb(q.numerator) / q.denominator).sqrt().sqrt()
    return s * s


def bank_point(tag, phis, real, M, Ht, out):
    t1 = time.time()
    s2 = s2_ball(phis)
    th_t, rh_t = tw.tails_at(s2, M, Ht, arb)
    rec = {"pt": tag, "M": M, "real": real,
           "phi": [[str(p.re), str(p.im)] for p in phis],
           "s2_str": str(s2), "s2": float(s2),
           "tail_theta": [str(x) if x is not None else "FAIL" for x in th_t],
           "tail_rho": [str(x) if x is not None else "FAIL" for x in rh_t]}
    for mode, tails in (("rho", rh_t), ("theta", th_t)):
        T = eval_point(phis, real, M, Ht, mode)
        comps, maxab = {}, [0.0] * 5
        for m, v in T.items():
            d = tw.DEG[m]
            tail = tails[d]
            assert tail is not None, f"tail hypothesis FAIL {tag} d={d}"
            if real:
                # true value real (real phi, real series); tail on re radius
                rr = (v.rad() + tail).upper()
                comps[str(m)] = [v.mid().str(40, radius=False), "0",
                                 rr.str(8, radius=False), "0"]
                ab = abs(float(v.mid()))
            else:
                rr = (v.real.rad() + tail).upper()
                ri = (v.imag.rad() + tail).upper()
                comps[str(m)] = [v.real.mid().str(40, radius=False),
                                 v.imag.mid().str(40, radius=False),
                                 rr.str(8, radius=False),
                                 ri.str(8, radius=False)]
                ab = abs(complex(float(v.real.mid()), float(v.imag.mid())))
            maxab[d] = max(maxab[d], ab)
        rec[mode] = {"comps": comps, "max_abs_by_deg": maxab,
                     "tail_by_deg": [float(x) for x in tails]}
    rec["wall_s"] = round(time.time() - t1, 2)
    out.append(rec)
    return rec


def main():
    outfn = sys.argv[1] if len(sys.argv) > 1 else "b3_payload.json"
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    deep = len(sys.argv) > 3 and sys.argv[3] == "deep"
    t0 = time.time()
    Ht = tw.harmonic_tables(max(M, 160) + 2)
    pts = points()
    # timed-pilot probe: center-first M-ladder, extrapolate before the batch
    probes = []
    for Mp in (10, 20):
        tp = time.time()
        eval_point(pts[0][1], True, Mp, Ht, "rho")
        probes.append(time.time() - tp)
    # cost model ~ (M+1)^2 (x-pairs) x const (jet dims saturate)
    est1 = probes[1] * ((M + 1) / 21) ** 2
    est = est1 * 2 * len(pts)
    print(f"[pilot-probe] center rho M=10:{probes[0]:.2f}s M=20:{probes[1]:.2f}s"
          f" -> est {est1:.1f}s/pt/mode, batch ~{est/60:.1f} min "
          f"({len(pts)} pts x 2 modes, M={M})", flush=True)
    assert est < 7200, "timed-pilot: batch extrapolates past 2h envelope, FREEZE"
    out = []
    for i, (tag, phis, real) in enumerate(pts):
        rec = bank_point(tag, phis, real, M, Ht, out)
        print(f"[ckpt] {i+1}/{len(pts)} {tag} wall {rec['wall_s']}s "
              f"elapsed {time.time()-t0:.0f}s s2={rec['s2']:.4f} "
              f"tail0_rho={rec['rho']['tail_by_deg'][0]:.2e} "
              f"tail4_rho={rec['rho']['tail_by_deg'][4]:.2e}", flush=True)
    if deep:
        for (tag, phis, real, Md) in (("corner_63", pts[64][1], True, 150),
                                      ("tilectr_worst", pts[65][1], False, 85)):
            rec = bank_point(tag + f"_M{Md}", phis, real, Md, Ht, out)
            print(f"[ckpt-deep] {tag} M={Md} wall {rec['wall_s']}s "
                  f"tail0_rho={rec['rho']['tail_by_deg'][0]:.2e} "
                  f"tail4_rho={rec['rho']['tail_by_deg'][4]:.2e}", flush=True)
    with open(outfn, "w") as f:
        json.dump({"prec_bits": 200, "M": M, "law": "towers6.py B3",
                   "points": out}, f)
    print(f"[done] {len(out)} point-records -> {outfn}; "
          f"wall {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()

"""S1 driver: both Brown-linear eliminations (stages 1-2) of a first-pair GROUP of (channel, sigma) on the dipole slice, in the chart/order of
that group (chart x_gauge = 1, integrate x_k, x_j of the pair, last variable x_l), with the GI engine fibr_gi (fibr.py arithmetic transported to
Q[x,t,j]/(j^2+1)) on the auto-extending Q(i) alphabet alphabet_gi.GIAlphabet.  Input = slice_common.build_group (exact gate cells, d_ab -> d(t)).
Outputs (under $SURD_WORK): ckpt/<tag>/stage{1,2}.pkl, receipt SLICE_S12_<tag>.json (censuses, t-degrees of the coefficients, letters, wall, RSS)."""
import os, sys, json, time, argparse, resource
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc
from fractions import Fraction as Fr
from flint import fmpq, fmpq_mpoly_ctx
import fibr_gi, alphabet_gi
WORK = slice_prov.WORK
def rss_mb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
def make_tag(channel, sigma, pair, gauge, cell=None, support_idx=None):
    return "%s_s%s_p%s_g%d%s%s" % (channel, sigma, pair, gauge, "" if cell is None else "_c%d" % cell, "" if support_idx is None else "_u%d" % support_idx)
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--channel", default="q_qbqgq"); ap.add_argument("--pair", default="34"); ap.add_argument("--gauge", type=int, default=0, help="0 = the default chart of the pair (slice_common.CHART)")
    ap.add_argument("--sigma", default="1234"); ap.add_argument("--cell", type=int, default=None); ap.add_argument("--support-idx", dest="support_idx", type=int, default=None, help="denominator-support subclass index (g_gggg sub-split)")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    gauge = a.gauge or sc.CHART[a.pair]; sigma = tuple(int(c) for c in a.sigma)
    tag = make_tag(a.channel, a.sigma, a.pair, gauge, a.cell, a.support_idx)
    rpath = os.path.join(WORK, "SLICE_S12_%s.json" % tag); cdir = slice_prov.wdir("ckpt", tag)
    if os.path.exists(rpath) and os.path.exists(os.path.join(cdir, "stage2.pkl")) and not a.force:
        print("exists (resume): %s" % rpath); return
    slice_prov.log_pid("slice_run:" + tag)
    logf = open(os.path.join(slice_prov.wdir("logs"), "slice_run_%s.log" % tag), "a")
    def log(s):
        l = "%s %s" % (time.strftime("%H:%M:%S"), s); print(l, flush=True); logf.write(l + "\n"); logf.flush()
    t00 = time.time()
    support = None
    if a.support_idx is not None:
        import piece_select, oracle_cells as OC
        g0 = gauge; classes = piece_select.support_classes(OC.build_cells(a.channel, g0, verbose=False), tuple(sorted(int(c) for c in a.pair)))
        support = sorted(classes[a.support_idx])
    G = sc.build_group(a.channel, gauge, a.pair, sigma, a.cell, support)
    xv = G["xvars"]; order = sc.order_for(a.pair, gauge)
    names = tuple("x%d" % j for j in xv) + ("t", "j")
    C = fmpq_mpoly_ctx.get(names); gens = C.gens(); imgs = list(gens[:4])
    def lift(P): return P.compose(*imgs, ctx=C)
    AL = alphabet_gi.GIAlphabet(C, cache_path=os.path.join(WORK, "ckpt", "gaussian_split_cache_%s.json" % "".join(names[:3])))
    E = fibr_gi.Engine(AL)
    c = Fr(1); D = {}
    for nm, p in G["atoms"].items():
        ci, Ei = AL.factor(lift(p), "input-atom:" + nm); c /= ci ** G["M"][nm]
        for i, e in Ei.items(): D[i] = D.get(i, 0) + e * G["M"][nm]
    NL = lift(G["N"]) * fmpq(c.numerator, c.denominator)
    ts = {(E.dkey(D), ()): NL}
    log("%s: cells %d, N %d monomials; atoms %s; ext %s; order %s|x%d; %.1fs RSS %.0f MB" % (tag, G["n_cells"], len(NL), G["M"], G["external_d_powers"], order[:2], order[2], time.time() - t00, rss_mb()))
    rec = {"tag": tag, "channel": a.channel, "sigma": a.sigma, "pair": a.pair, "gauge": gauge, "cell": a.cell, "support_idx": a.support_idx, "support_atoms": support, "order": ["x%d" % q for q in order],
           "slice": "d12=t^2, d34=1, d13=d24=t^2/4+3t/10+1/4, d14=d23=t^2/4-3t/10+1/4 (sigma acts on labels)", "dsigma_polys": G["dsigma_polys"],
           "n_cells": G["n_cells"], "N_terms": len(NL), "atoms_M": G["M"], "external_d_powers": G["external_d_powers"], "ext_factor_poly_t": str(G["ext_factor_poly"]), "stages": []}
    for k, var in ((1, "x%d" % order[0]), (2, "x%d" % order[1])):
        t0 = time.time(); c0 = time.process_time(); E.stats = {"factor_calls": 0, "new_denominators": {}}
        out = {}; keys = sorted(ts.keys(), key=lambda kk: repr(kk)); binfo = []
        for b, key in enumerate(keys):
            t1 = time.time(); o = E.integrate({key: ts[key]}, var, log=None, consolidate=True, tag="s%d b%d" % (k, b)); E.ts_add_into(out, o)
            binfo.append({"block": b, "weight_in": sum(len(w) for w in key[1]), "in_monomials": len(ts[key]), "out_keys": len(o), "out_monomials": E.ts_nmon(o), "wall_s": round(time.time() - t1, 2)})
            log("  stage %d block %d/%d: in %d -> out %d keys %d monomials, %.1fs, RSS %.0f MB" % (k, b + 1, len(keys), len(ts[key]), len(o), E.ts_nmon(o), time.time() - t1, rss_mb()))
        out = E.consolidate(out); cen = E.census(out)
        sc.ts_dump(out, os.path.join(cdir, "stage%d.pkl" % k), names, AL, extra={"tag": tag, "order": order, "ext_factor_poly_t": str(G["ext_factor_poly"]), "external_d_powers": G["external_d_powers"], "sigma": a.sigma, "channel": a.channel, "pair": a.pair, "gauge": gauge})
        rec["stages"].append({"stage": k, "var": var, "wall_s": round(time.time() - t0, 2), "cpu_s": round(time.process_time() - c0, 2), "out_census": {kk: v for kk, v in cen.items() if kk != "denominator_letters"},
                              "blocks": binfo, "rss_MB_after": round(rss_mb(), 1), "timers": {kk: round(v, 2) for kk, v in E.stats.items() if kk.startswith("t_")}, "gi_stats": dict(AL.stats)})
        log("stage %d (%s): %s wall %.1fs RSS %.0f MB" % (k, var, rec["stages"][-1]["out_census"], time.time() - t0, rss_mb()))
        ts = out
    jx = names.index("x%d" % order[2]); jt = names.index("t")
    def tdeg_den(Dk): return sum(e * int(AL.letters[i].degrees()[jt]) for i, e in Dk)
    byw = {}
    for (Dk, Wk), N in ts.items():
        w = sum(len(wd) for wd in Wk); d = byw.setdefault(w, {"keys": 0, "monomials": 0, "max_deg_x_num": 0, "max_deg_t_num": 0, "max_deg_t_den": 0})
        d["keys"] += 1; d["monomials"] += len(N); d["max_deg_x_num"] = max(d["max_deg_x_num"], int(N.degrees()[jx])); d["max_deg_t_num"] = max(d["max_deg_t_num"], int(N.degrees()[jt])); d["max_deg_t_den"] = max(d["max_deg_t_den"], tdeg_den(Dk))
    rec["final_census_by_weight"] = {str(k): v for k, v in sorted(byw.items())}
    xl = {}
    for i, L in enumerate(AL.letters):
        dg = int(L.degrees()[jx])
        if dg >= 1: xl[str(L)] = {"idx": i, "deg_x": dg, "deg_t": int(L.degrees()[jt]), "has_j": int(L.degrees()[names.index("j")]) > 0}
    used = set(i for (Dk, Wk) in ts for i, e in Dk) | set(i for (Dk, Wk) in ts for wd in Wk for h in wd for i, e in h[1])
    rec["last_variable_letters_used(deg_x>=1)"] = {s: v for s, v in xl.items() if v["idx"] in used}
    rec["nonlinear_last_variable_letters_used"] = {s: v for s, v in xl.items() if v["idx"] in used and v["deg_x"] >= 2}
    rec["n_letters_alphabet"] = len(AL.letters); rec["alphabet"] = [str(L) for L in AL.letters]; rec["gi_stats"] = dict(AL.stats)
    rec["hlog_letters"] = sorted(set(E.hl_str(h) for (Dk, Wk) in ts for wd in Wk for h in wd))
    rec["wall_s_total"] = round(time.time() - t00, 1); rec["cpu_s_total"] = round(time.process_time(), 1); rec["peak_rss_MB"] = round(rss_mb(), 1); rec["checkpoint_dir"] = cdir
    slice_prov.write_receipt(rpath, rec, inputs=[os.path.join(HERE, "fibr_gi.py"), os.path.join(HERE, "alphabet_gi.py"), os.path.join(HERE, "slice_common.py"), os.path.join(sc.GATE, "cells", "%s_g%d.pkl" % (a.channel, gauge))])
    log("wrote %s ; nonlinear last-variable letters: %s ; wall %.1fs RSS %.0f MB" % (rpath, {s[:50]: (v["deg_x"], v["deg_t"], v["has_j"]) for s, v in rec["nonlinear_last_variable_letters_used"].items()}, time.time() - t00, rss_mb()))
if __name__ == "__main__":
    main()

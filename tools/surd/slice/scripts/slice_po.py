"""Per-group numeric reference ('piece oracle', GATE route): the 3-fold integral of a first-pair group's cells at slice d(t) by
gate/scripts/oracle_core (exact InstanceRF of the sub-cells, inner residues in Arb, 2-D exp-sinh), two charts/routes A=(g4;o1,i2|a3), B=(g1;o2,i3|a4).
Cache: ckpt/po/<channel>_s<sigma>_p<pair>_t<num>_<den>.json (write-once)."""
import os, sys, json, time
import slice_prov, slice_common as sc
from fractions import Fraction as Fr
from flint import arb, ctx
import oracle_core as CO, oracle_cells as OC, piece_select
WORK = slice_prov.WORK; PODIR = os.path.join(WORK, "ckpt", "po")
_CELLS = {}
def cells(ch, g):
    if (ch, g) not in _CELLS: _CELLS[(ch, g)] = OC.build_cells(ch, g, verbose=False)
    return _CELLS[(ch, g)]
def piece_oracle(channel, sigma, pair, tval, levels=("1/18", "1/60"), digits=34, force=False):
    tval = Fr(tval); fn = os.path.join(PODIR, "%s_s%s_p%s_t%s_%s.json" % (channel, sigma, pair, tval.numerator, tval.denominator))
    if os.path.exists(fn) and not force: return json.load(open(fn))
    t0 = time.time(); sig = tuple(int(c) for c in sigma); pr = tuple(sorted(int(c) for c in pair))
    dsig = sc.sigma_d(sc.d_of_t(tval), sig); ho, hi = Fr(levels[0]), Fr(levels[1]); rr = {}
    for r, (g, o, i, an) in (("A", (4, 1, 2, 3)), ("B", (1, 2, 3, 4))):
        ctx.prec = 384
        sub = piece_select.select_cells(cells(channel, g), pr, None); rf = CO.InstanceRF(sub, dsig); rev = CO.RouteEval(rf, o, i, an)
        t1 = time.time(); out = CO.integrate_instance(rev, ho, hi, 4.8, digits, prec0=384); w = time.time() - t1
        rr[r] = {"chart": g, "roles_o_i_a": [o, i, an], "n_cells": len(sub["cells"]), "value": out["value"].str(45, radius=True), "value_mid": out["value"].mid().str(50, radius=False), "wall_s": round(w, 1), "max_prec": out["stats"]["max_prec"]}
    ctx.prec = 512
    va, vb = arb(rr["A"]["value_mid"]), arb(rr["B"]["value_mid"]); rel = abs(va - vb) / abs(va) if va != 0 else abs(va - vb)
    rec = {"channel": channel, "sigma": sigma, "pair": pair, "t": str(tval), "d_sigma": {"%d%d" % k: str(v) for k, v in dsig.items()}, "routes": rr, "A_vs_B_rel": rel.mid().str(5, radius=False),
           "levels": list(levels), "wall_s": round(time.time() - t0, 1), "stamp_utc": slice_prov.stamp(), "producer": slice_prov.producer(inputs=[os.path.join(sc.GATE, "cells", "%s_g%d.pkl" % (channel, g)) for g in (4, 1)])}
    slice_prov.wdir("ckpt", "po"); tmp = fn + ".tmp%d" % os.getpid(); json.dump(rec, open(tmp, "w"), indent=1); os.replace(tmp, fn)
    return rec
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--channel", required=True); ap.add_argument("--sigma", required=True); ap.add_argument("--pair", required=True); ap.add_argument("--t", required=True)
    a = ap.parse_args(); slice_prov.log_pid("slice_po")
    for tv in a.t.split(","):
        r = piece_oracle(a.channel, a.sigma, a.pair, Fr(tv)); print(tv, r["routes"]["A"]["value"], "A_vs_B", r["A_vs_B_rel"], "wall", r["wall_s"])

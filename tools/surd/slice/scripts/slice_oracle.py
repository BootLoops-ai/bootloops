#!/usr/bin/env python3
"""S0: >=30-digit oracle for the LO collinear E4C on the dipole slice at rational t, per channel and per jet (n_f = 5), built on the
oracle machinery of gate/scripts/oracle.py in d_ab mode (task_run: exact cell assembly per (channel, chart, sigma), inner integral by residues in Arb,
2-D tensor exp-sinh, TWO ROUTES A=(g4;o1,i2|a3), B=(g1;o2,i3|a4); certified digits = -log10 |A-B|/|A|), with the labeled distance matrix taken
DIRECTLY from the slice: d12=t^2, d34=1, d13=d24=q+(t), d14=d23=q-(t) (slice_common.d_of_t).
Normalization of the closed forms of this package:   F_ch(t) = w_flavor(nf=5) * S_ch * sum_{sigma in S4} I_sigma(d(t))   (i.e. oracle.py's
G4_channel with xL^3 replaced by d34^3 = 1; the xL^3-normalized value = max(1,t^2)^3 * F(t); F(1/t) = t^6 F(t) exactly by the relabeling (13)(24)).
Checkpoints: ckpt/oracle/<tkey>.jsonl (resume = rerun).  Writes SLICE_S0_ORACLE_t<tkey>.json per t and SLICE_S0_ORACLE.json (table)."""
import os, sys, json, time, argparse, multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc
from fractions import Fraction as F
import oracle as ORA, oracle_cells as OC, oracle_core as CO
from flint import arb, ctx
WORK = slice_prov.WORK; CK = os.path.join(WORK, "ckpt", "oracle")
CHANNELS = OC.QUARK_JET + OC.GLUON_JET + ["n4"]
def tkey(t): t = F(t); return "t%s_%s" % (t.numerator, t.denominator)
def load_ckpt(path):
    done = {}
    if os.path.exists(path):
        for line in open(path):
            line = line.strip()
            if not line: continue
            r = json.loads(line); done[(r["channel"], tuple(r["sigma"]), r["route"], r["ho"], r["hi"])] = r
    return done
def _digits(rel): return ORA._digits(rel)
def _run_task(targs):
    t, args = targs; return t, ORA.task_run(args)
def assemble(t, done, classes, chans, nf=F(5)):
    ctx.prec = 384; out = {"channels": {}, "jets": {}}
    t = F(t); xL = max(sc.d_of_t(t).values())
    for ch in chans:
        w = OC.flavor_weight(ch, nf) * OC.SFAC[ch]; wa = arb(w.numerator) / arb(w.denominator)
        per = {}
        for r in ("A", "B"):
            RP = ORA.route_params(r, ch)
            if not all((ch, tuple(sig), r, str(RP["ho"]), str(RP["hi"])) in done for (sig, mult, dsig) in classes): continue
            tot = arb(0); tq = 0.0
            for (sig, mult, dsig) in classes:
                rec = done[(ch, tuple(sig), r, str(RP["ho"]), str(RP["hi"]))]; tot += mult * ORA._ball(rec); tq += rec["t_quad"] + rec["t_build"] + rec["t_load"]
            per[r] = (wa * tot, tq, [str(RP["ho"]), str(RP["hi"])])
        if not per: continue
        rec = {"per_route": {r: v[0].str(50, radius=True) for r, v in per.items()}, "wall_s": {r: round(v[1], 1) for r, v in per.items()}, "levels": {r: v[2] for r, v in per.items()},
               "flavor_weight_x_sfac": str(w), "per_sigma_class_routeA": {}}
        RP = ORA.route_params("A", ch)
        for (sig, mult, dsig) in classes:
            k = (ch, tuple(sig), "A", str(RP["ho"]), str(RP["hi"]))
            if k in done: rec["per_sigma_class_routeA"]["".join(map(str, sig))] = {"mult": mult, "I_sigma": ORA._ball(done[k]).str(40, radius=True)}
        prim = "A" if "A" in per else "B"; rec["primary_route"] = prim; rec["value"] = per[prim][0].mid().str(45, radius=False)
        if len(per) == 2:
            dif = abs(per["A"][0] - per["B"][0]); rel = dif / abs(per["A"][0]); rec["abs_diff"] = dif.str(5); rec["rel_diff"] = rel.str(5); rec["certified_digits"] = _digits(rel)
        rec["value_e4c_normalization(xL^3)"] = (per[prim][0] * (arb(xL.numerator) / arb(xL.denominator)) ** 3).mid().str(40, radius=False)
        out["channels"][ch] = rec
    for jn, lst in (("q_jet", OC.QUARK_JET), ("g_jet", OC.GLUON_JET)):
        if all(c in out["channels"] and "certified_digits" in out["channels"][c] for c in lst):
            A = arb(0); B = arb(0); err = arb(0)
            for c in lst:
                A += arb(out["channels"][c]["per_route"]["A"]); B += arb(out["channels"][c]["per_route"]["B"]); err += arb(out["channels"][c]["abs_diff"])
            out["jets"][jn] = {"value": A.mid().str(45, radius=False), "route_B": B.mid().str(45, radius=False), "certified_abs_err(sum |A-B|)": err.str(5), "certified_digits": _digits(err / abs(A)),
                               "value_e4c_normalization(xL^3)": (A * (arb(xL.numerator) / arb(xL.denominator)) ** 3).mid().str(40, radius=False)}
    out["xL=max d_ab"] = str(xL); out["d(t)"] = {"%d%d" % k: str(v) for k, v in sc.d_of_t(t).items()}
    return out
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--t", required=True, help="comma list of rational t"); ap.add_argument("--channels", default=",".join(CHANNELS))
    ap.add_argument("--nproc", type=int, default=2); ap.add_argument("--digits", type=int, default=30); ap.add_argument("--prec0", type=int, default=384); ap.add_argument("--assemble-only", dest="ao", action="store_true")
    ap.add_argument("--t-major", dest="tmajor", action="store_true", help="run tasks t-major in the order of --t (default: channel-major)")
    a = ap.parse_args(); slice_prov.log_pid("slice_oracle"); os.makedirs(CK, exist_ok=True)
    chans = a.channels.split(","); ts = [F(x) for x in a.t.split(",")]
    tasks = []; meta = {}
    for t in ts:
        d = sc.d_of_t(t); classes = CO.sigma_classes(d); ck = os.path.join(CK, tkey(t) + ".jsonl"); done = load_ckpt(ck); meta[t] = (classes, ck)
        assert len(classes) in (6, 12), len(classes)
        for ch in chans:
            for (sig, mult, dsig) in classes:
                for r in ("A", "B"):
                    RP = ORA.route_params(r, ch); key = (ch, tuple(sig), r, str(RP["ho"]), str(RP["hi"]))
                    if key in done: continue
                    tasks.append((t, (ch, tuple(sig), mult, [[list(k), str(vv)] for k, vv in dsig.items()], r, a.digits, a.prec0, (str(RP["ho"]), str(RP["hi"])))))
    order = {"g_gggg": 0, "g_qbggq": 1, "n4": 2, "q_qbqgq": 3, "g_qbqqbq": 4, "q_gggq": 5, "g_qbpqpqbq": 6, "q_qbpqpgq": 7}
    tasks.sort(key=(lambda x: (ts.index(x[0]), order.get(x[1][0], 9))) if a.tmajor else (lambda x: (order.get(x[1][0], 9), x[0])))
    print("slice oracle: t =", [str(t) for t in ts], "channels", chans, "tasks to run", len(tasks), "pid", os.getpid(), flush=True)
    T0 = time.time()
    if not a.ao and tasks:
        def sink(t, r):
            with open(meta[t][1], "a") as f: f.write(json.dumps(r) + "\n")
            print("%7.0fs done t=%s %s sigma=%s %s I=%s quad %.0fs maxprec %d rss %s" % (time.time() - T0, t, r["channel"], "".join(map(str, r["sigma"])), r["route"], r["value_ball"][:36], r["t_quad"], r["stats"]["max_prec"], r.get("maxrss_kb")), flush=True)
        if a.nproc > 1:
            with mp.Pool(a.nproc, maxtasksperchild=4) as pool:
                for t, r in pool.imap_unordered(_run_task, tasks): sink(t, r)     # checkpoint each result as soon as it lands
        else:
            for t, args in tasks: sink(t, ORA.task_run(args))
    table = {}
    for t in ts:
        classes, ck = meta[t]; done = load_ckpt(ck); res = assemble(t, done, classes, chans)
        res.update({"t": str(t), "n_sigma_classes": len(classes), "sigma_classes": [["".join(map(str, s)), m] for s, m, d_ in classes], "checkpoint": ck, "n_records": len(done),
                    "normalization": "F_ch(t) = w_flavor(nf=5) * S_ch * sum_sigma I_sigma(d(t)), d34 = 1 (NO xL^3); xL^3-normalized value = max(1,t^2)^3 * F"})
        out = os.path.join(WORK, "SLICE_S0_ORACLE_%s.json" % tkey(t))
        slice_prov.write_receipt(out, res, inputs=[ck] + [os.path.join(sc.GATE, "cells", "%s_g%d.pkl" % (c, g)) for c in chans for g in (4, 1)])
        table[str(t)] = {"channels": {c: {"value": v["value"], "certified_digits": v.get("certified_digits"), "e4c_norm": v["value_e4c_normalization(xL^3)"]} for c, v in res["channels"].items()}, "jets": res["jets"]}
        print("t=%s" % t, json.dumps(table[str(t)], indent=0)[:1500], flush=True)
    slice_prov.write_receipt(os.path.join(WORK, "SLICE_S0_ORACLE.json"), {"table": table, "per_t_receipts": [os.path.join(WORK, "SLICE_S0_ORACLE_%s.json" % tkey(t)) for t in ts]}, inputs=[meta[t][1] for t in ts])
if __name__ == "__main__":
    main()

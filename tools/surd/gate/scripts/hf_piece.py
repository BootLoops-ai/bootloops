"""Build (and optionally run) a HyperFLINT request for one PIECE = (channel, sigma, pair-group or 'all')
at a split-signature rational shape, chart x_g=1, order (k, j, l): expr = N/prod atoms^M with the last
variable rescaled x_l = lam*y (keeps rational last-variable letters symbolic inside the engine).
Group assignment rule: a cell belongs to the first vertex pair {k,j}
(lexicographic) contained in no 3-particle quadric and no partial-sum plane of its support.
The library part (first_pair, group_cells) needs nothing external; the command line (--run) drives the HyperFLINT engine of the
sibling package tools/subtropica: set SUBTROPICA_HF_BIN (engine wrapper) and SUBTROPICA_MZV_DATA (its MZV reduction table).
Request/response/meta files go under $SURD_WORK/hf_piece/."""
import os, sys, json, re, itertools, argparse, subprocess, time, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import oracle_cells as OC, oracle_core as CO, gate_prov
GATE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HF = os.environ.get("SUBTROPICA_HF_BIN", "")          # HyperFLINT engine wrapper (see tools/subtropica/GUIDE.md); needed only by main()
MZV = os.environ.get("SUBTROPICA_MZV_DATA", "")       # its MZV reduction table

def sset(nm): return frozenset(int(c) for c in nm[1:]) if nm[0] == "s" else frozenset(int(c) for c in nm[3:])
def first_pair(dk):
    atoms = [nm for nm, pw in dk if (nm[0] == "s" and len(nm) == 4) or nm.startswith("ell")]
    S = [sset(a) for a in atoms if a[0] == "s"]; T = [sset(a) for a in atoms if a.startswith("ell")]
    for p in itertools.combinations([1, 2, 3, 4], 2):
        if not any(frozenset(p) <= x for x in S) and not any(frozenset(p) <= x for x in T):
            return p
    return None

def group_cells(cellres, pair):
    if pair == "all":
        return cellres
    sub = dict(cellres); sub["cells"] = [(dk, conv) for dk, conv in cellres["cells"] if first_pair(dk) == pair]
    return sub

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--channel", required=True); ap.add_argument("--shape", required=True, help="u,ub,v,vb rationals (split-signature frame w=(u,1,0,v), wb=(ub,1,0,vb))")
    ap.add_argument("--sigma", default="1234"); ap.add_argument("--pair", default="all", help="kj or all")
    ap.add_argument("--gauge", type=int, required=True); ap.add_argument("--order", required=True, help="e.g. 341 = integrate x3,x4 then x1 (last, rescaled)")
    ap.add_argument("--tag", required=True); ap.add_argument("--run", action="store_true"); ap.add_argument("--check-div", dest="checkdiv", action="store_true")
    ap.add_argument("--alg", action="store_true", help="algebraic_letters tier (Wm/Wp for last-variable quadratics)")
    ap.add_argument("--nint", type=int, default=3, help="number of variables HyperFLINT integrates (3 = full; 2 = leave the last variable symbolic)")
    a = ap.parse_args()
    if not MZV or (a.run and not HF): sys.exit("hf_piece.py: set SUBTROPICA_MZV_DATA (and SUBTROPICA_HF_BIN for --run) to a built HyperFLINT engine of tools/subtropica")
    u, ub, v, vb = (F(x) for x in a.shape.split(","))
    d = CO.dmatrix_split(u, ub, v, vb)
    sig = tuple(int(c) for c in a.sigma)
    dsig = CO.sigma_dmatrix(d, sig)
    cellres = OC.build_cells(a.channel, a.gauge, verbose=False)
    pair = "all" if a.pair == "all" else tuple(sorted(int(c) for c in a.pair))
    sub = group_cells(cellres, pair)
    rf = CO.InstanceRF(sub, dsig)
    order = [int(c) for c in a.order]
    last = order[-1]
    Ns = str(rf.N)
    den = " * ".join("(%s)^%d" % (str(rf.atom_polys[nm]), rf.M[nm]) for nm in rf.atom_names)
    expr = "(%s)/(%s)" % (Ns, den)
    if a.nint == 3:
        expr = "lam*(" + re.sub(r"\bx%d\b" % last, "(lam*y%d)" % last, expr) + ")"
        vars_int = ["x%d" % k for k in order[:-1]] + ["y%d" % last]
        allv = vars_int + ["lam"]
    else:
        vars_int = ["x%d" % k for k in order[:a.nint]]
        allv = ["x%d" % k for k in order] 
    req = {"op": "hyperflint", "expr": expr, "vars_int": vars_int, "vars": allv, "mzv_data_path": MZV,
           "check_divergences": bool(a.checkdiv), "canonical_emission": True}
    if a.alg:
        req["algebraic_letters"] = True
    T3 = gate_prov.wdir("hf_piece")
    rq = os.path.join(T3, "hf_request_%s.json" % a.tag)
    json.dump(req, open(rq, "w"))
    meta = {"channel": a.channel, "sigma": list(sig), "pair_group": a.pair, "gauge": a.gauge, "order": order, "nint": a.nint,
            "shape_split": [str(u), str(ub), str(v), str(vb)], "d_sigma": {"%d%d" % k: str(vv) for k, vv in dsig.items()},
            "n_cells": len(sub["cells"]), "N_terms": len(rf.N), "N_degrees": [int(x) for x in rf.N.degrees()],
            "atoms": {nm: {"poly": str(rf.atom_polys[nm]), "M": int(rf.M[nm])} for nm in rf.atom_names}, "expr_chars": len(expr),
            "request": rq, "request_sha256": gate_prov.sha(rq), "check_divergences": bool(a.checkdiv), "algebraic_letters": bool(a.alg)}
    if a.run:
        rs = os.path.join(T3, "hf_response_%s.json" % a.tag); er = os.path.join(T3, "hf_stderr_%s.log" % a.tag)
        t0 = time.time()
        env = dict(os.environ, OMP_NUM_THREADS="1", HF_MAX_THREADS_PER_CALL="1", HF_LR_TIME_BUDGET_S=os.environ.get("HF_LR_TIME_BUDGET_S", "3600"))
        with open(rq) as fin, open(rs, "w") as fout, open(er, "w") as ferr:
            # ulimit -s unlimited: the CLI's regex JSON parser overflows the default 8 MB stack on requests >~28 kB (rc 139)
            p = subprocess.run(["bash", "-c", "ulimit -s unlimited; exec /usr/bin/time -v nice -n 10 %s eval-json" % HF], stdin=fin, stdout=fout, stderr=ferr, env=env)
        meta["hf_rc"] = p.returncode; meta["hf_wall_s"] = round(time.time() - t0, 1)
        tv = open(er).read()
        m1 = re.search(r"Maximum resident set size \(kbytes\): (\d+)", tv); m2 = re.search(r"Elapsed \(wall clock\).*: (.*)", tv)
        meta["hf_max_rss_kb"] = int(m1.group(1)) if m1 else None; meta["hf_elapsed"] = m2.group(1).strip() if m2 else None
        try:
            resp = json.loads(open(rs).read().strip().splitlines()[-1])
            meta["hf_keys"] = [k for k in resp if k != "result"]; meta["hf_n_terms"] = len(resp.get("result", [])) if "result" in resp else None
            meta["hf_flags"] = {k: resp[k] for k in ("divergent", "failed", "error", "reason") if k in resp}
        except Exception as e:
            meta["hf_parse_error"] = str(e)
        meta["response"] = rs
    meta["producer"] = gate_prov.producer(inputs=[rq, os.path.join(GATE, "cells", "%s_g%d.pkl" % (a.channel, a.gauge))], modules=["hf_piece.py", "oracle_cells.py", "oracle_core.py", "gate_prov.py"])
    json.dump(meta, open(os.path.join(T3, "hf_piece_%s.meta.json" % a.tag), "w"), indent=1)
    print(json.dumps({k: meta[k] for k in meta if k not in ("atoms", "producer", "d_sigma")}, indent=1))

if __name__ == "__main__":
    main()

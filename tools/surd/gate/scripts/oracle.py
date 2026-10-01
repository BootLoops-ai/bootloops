#!/usr/bin/env python3
""">=30-digit oracle for the LO QCD collinear four-point energy correlator (E4C).

  oracle.py --channel q_jet|g_jet|<channel> --shape 'u=3/7+2/5j;v=-1/3+4/9j' --digits 30 [--nproc 3]

Object (conventions: nf=5, CF=4/3, CA=3, TF=1/2, eps=0):
  G4_channel(shape) = xL^3 (1/S!) w_flavor sum_{sigma in S4} int_simplex d^4xi delta(1-sum xi)
                      (xi1 xi2 xi3 xi4)^2 <P_channel>(xi; s_ij = xi_i xi_j d_{sigma(i)sigma(j)}) / s_1234^3,
  shape = the four points {0, 1, u, v} of the plane (u, v complex with RATIONAL real/imaginary
  parts), d_ab = |w_a - w_b|^2, xL = max d_ab.  q_jet = q_qbpqpgq + q_qbqgq + q_gggq,
  g_jet = g_qbpqpqbq + g_qbqqbq + g_qbggq + g_gggg (flavor weights (nf-1),1,1 ; nf(nf-1)/2,nf,nf,1
  and identical-particle factors 1,1/2,1 ; 1,1/4,1/2,1/2 as in oracle_cells.flavor_weight / SFAC).

Method (oracle_core.py): Cheng-Wu chart x_g=1 (P is z-homogeneous, exponent k=6+h of
sum(xi) verified), exact assembly of each S4 instance as ONE rational function with multilinear
denominator atoms, innermost integral in closed form by residues (Arb ball arithmetic), remaining
2-D integral by the centered tensor exp-sinh rule.  TWO ROUTES = two different charts/analytic
variables/step sizes:  A = (gauge x4=1, analytic x3, outer x1, inner x2, steps 1/20,1/40),
B = (gauge x1=1, analytic x4, outer x2, inner x3, steps 1/18,1/36)  [--levelA/--levelB override].  Certified digits =
-log10 |A-B|/|A| (each route's own Arb radius is far below).  Checkpoints: every (shape, channel,
sigma-class, route) result is appended to $SURD_WORK/oracle_ckpt/<shapekey>.jsonl; rerunning the same
command resumes.
"""
import os, sys, json, time, argparse, hashlib, subprocess, multiprocessing as mp
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import oracle_cells as OC
import oracle_core as CO
import gate_prov
from flint import arb, ctx

GATE = os.path.dirname(HERE)
CKDIR = os.path.join(gate_prov.WORK, "oracle_ckpt")          # created in main() at run time, never at import
ROUTES = {"A": dict(g=4, o=1, i=2, a=3, ho=F(1, 20), hi=F(1, 40), S=4.8),
          "B": dict(g=1, o=2, i=3, a=4, ho=F(1, 18), hi=F(1, 36), S=4.8)}
# per-channel DE step policy (measured convergence ladders): channels whose tensor exp-sinh error decays
# slowly get finer grids so that BOTH routes are below ~1e-32.  (outer step, inner step): the inner
# (moving-singularity) direction converges slower (q_gggq at the square shape, route A: outer 1/16 already
# <1e-37, inner 1/32 -> 4e-26, 1/40 -> ~1e-31, 1/48 -> ~1e-33, 1/64 -> reference), hence the asymmetric
# 'slow' grids.
LEVELS = {"fast": {"A": (F(1, 20), F(1, 40)), "B": (F(1, 18), F(1, 36))},
          "gl": {"A": (F(1, 18), F(1, 36)), "B": (F(1, 16), F(1, 32))},
          "slow": {"A": (F(1, 18), F(1, 60)), "B": (F(1, 16), F(1, 54))},
          "slow2": {"A": (F(1, 24), F(1, 72)), "B": (F(1, 20), F(1, 64))}}
# gluon sigma=id probes at a second generic shape: (1/16,1/32) already 3.5e-35 (g_qbpqpqbq, g_qbqqbq),
# 1.6e-36 (g_qbggq); q_gggq is the slow channel (inner direction).  Slow INNER convergence is
# (instance, chart)-dependent (inner 1/32 -> 1e-26..1e-27, inner 1/54..1/64 -> <1e-37; outer 1/16 always
# <1e-37), so every channel except the two verified-fast quark channels uses the asymmetric fine-inner grids.
CHANNEL_SPEED = {"q_qbpqpgq": "fast", "q_qbqgq": "fast", "q_gggq": "slow",
                 "g_qbpqpqbq": "slow", "g_qbqqbq": "slow", "g_qbggq": "slow", "g_gggg": "slow"}
OVERRIDE = {}


def route_params(route, ch):
    R = dict(ROUTES[route])
    ho, hi = LEVELS[CHANNEL_SPEED.get(ch, "slow")][route]
    R["ho"], R["hi"] = ho, hi
    if route in OVERRIDE:
        R["ho"], R["hi"] = OVERRIDE[route]
    return R


def parse_cplx(s):
    s = s.strip().replace(" ", "").replace("i", "j")
    # forms: a+bj, a-bj, bj, a
    if "j" not in s:
        return (F(s), F(0))
    body = s[:-1]
    # split at the last + or - that is not leading and not after '/'
    k = max(body.rfind("+"), body.rfind("-"))
    while k > 0 and body[k - 1] == "/":
        k = max(body.rfind("+", 0, k), body.rfind("-", 0, k))
    if k <= 0:
        return (F(0), F(body) if body not in ("", "+", "-") else F(body + "1"))
    re, im = body[:k], body[k:]
    if im in ("+", "-"):
        im = im + "1"
    return (F(re), F(im))


def parse_shape(s):
    kv = dict(p.split("=") for p in s.replace(" ", "").split(";") if p)
    u = parse_cplx(kv["u"]); v = parse_cplx(kv["v"])
    return u, v


def shape_key(u, v):
    def c(z):
        return "%s%s%sj" % (z[0], "+" if z[1] >= 0 else "-", abs(z[1]))
    txt = "u=%s;v=%s" % (c(u), c(v))
    return txt, hashlib.sha256(txt.encode()).hexdigest()[:12]


def arb_to_str(x, n=60):
    return x.str(n, radius=True)


def task_run(args):
    """one (channel, sigma-class, route) -> dict (JSON-able).  Runs in a worker process."""
    ch, sig, mult, dsig_items, route, digits, prec0, lv = args
    R = dict(ROUTES[route]); R["ho"], R["hi"] = F(lv[0]), F(lv[1])
    dsig = {tuple(k): F(vv) for k, vv in dsig_items}
    t0 = time.time()
    cellres = OC.build_cells(ch, R["g"], verbose=False)
    t1 = time.time()
    ctx.prec = prec0
    rf = CO.InstanceRF(cellres, dsig)
    rev = CO.RouteEval(rf, R["o"], R["i"], R["a"])
    t2 = time.time()
    res = CO.integrate_instance(rev, R["ho"], R["hi"], R["S"], digits + 4, prec0=prec0)
    t3 = time.time()
    ctx.prec = prec0
    val = res["value"]
    try:
        import resource
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    except Exception:
        rss = None
    return {"channel": ch, "sigma": list(sig), "mult": mult, "dsig": dsig_items, "route": route,
            "roles": {k: R[k] for k in ("g", "o", "i", "a")}, "ho": str(R["ho"]), "hi": str(R["hi"]), "S": R["S"],
            "value_mid": val.mid().str(70), "value_rad": val.rad().str(5), "value_ball": val.str(60, radius=True),
            "N_terms": int(rf.nN), "N_degrees": [int(x) for x in rf.degN], "atoms_M": {k: int(vv) for k, vv in rf.M.items()},
            "t_load": round(t1 - t0, 2), "t_build": round(t2 - t1, 2), "t_quad": round(t3 - t2, 2),
            "stats": res["stats"], "maxrss_kb": rss, "prec0": prec0,
            "stamp_utc": subprocess.check_output(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip()}


def load_symmetries():
    """candidate relabeling symmetries per channel (oracle_cells.RELABELING_SYMMETRIES); used only for the symmetry-partner
    certificate, which compares independently computed instances and so also tests the candidates."""
    return {ch: [tuple(t) for t in v] for ch, v in OC.RELABELING_SYMMETRIES.items()}


def sym_certificate(ch, classes, done, route, RP, taus):
    """Symmetry-partner certificate: for an exact relabeling symmetry tau of |M|^2 the instances sigma and
    sigma o tau have equal integrals but are computed from different d-matrices (different charts/grids).
    Returns (sum_sigma mult*max_tau|I_sigma - I_sigma.tau|  as arb, n_pairs, max_rel) or None."""
    if not taus:
        return None
    keyof = {}
    vals = {}
    for (sig, mult, dsig) in classes:
        k = tuple(dsig[p] for p in CO.PAIRS)
        keyof[tuple(sig)] = k
        rec = done.get((ch, tuple(sig), route, str(RP["ho"]), str(RP["hi"])))
        if rec is None:
            return None
        vals[k] = _ball(rec)
    # map any sigma to its class key
    def key_sigma(sig):
        return tuple(CO.sigma_dmatrix_from_classes(sig) for _ in [0]) if False else None
    err = arb(0); npairs = 0; maxrel = arb(0)
    # need d for arbitrary sigma: rebuild from the class representative structure: classes carry dsig for reps only;
    # recover the point-distance matrix d from the identity representative
    drep = None
    for (sig, mult, dsig) in classes:
        inv = {sig[i]: i + 1 for i in range(4)}      # point -> particle
        drep = {}
        for a, b in CO.PAIRS:                        # d[(p,q)] for points p<q
            pa, pb = inv[a], inv[b]
            drep[(a, b)] = dsig[(pa, pb) if pa < pb else (pb, pa)]
        break
    for (sig, mult, dsig) in classes:
        k = keyof[tuple(sig)]
        best = None
        for tau in taus:
            st = tuple(sig[tau[i] - 1] for i in range(4))
            k2 = tuple(CO.sigma_dmatrix(drep, st)[p] for p in CO.PAIRS)
            if k2 == k or k2 not in vals:
                continue
            dlt = abs(vals[k] - vals[k2])
            best = dlt if best is None else (dlt if float(dlt.upper()) > float(best.upper()) else best)
            npairs += 1
            rel = dlt / abs(vals[k])
            if float(rel.upper()) > float(maxrel.upper()):
                maxrel = rel
        if best is None:
            return None          # some class has no non-trivial partner (symmetric shape): no certificate
        err += mult * best
    return err, npairs, maxrel


def load_ckpt(path):
    done = {}
    if os.path.exists(path):
        for line in open(path):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            done[(r["channel"], tuple(r["sigma"]), r["route"], r["ho"], r["hi"])] = r
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--channel", required=True, help="q_jet | g_jet | channel name | comma list")
    ap.add_argument("--shape", required=True, help="'u=3/7+2/5j;v=-1/3+4/9j' (points {0,1,u,v}); or 'split:u=..;ub=..;v=..;vb=..' (census-frame split-signature rational point, all d_ab>0)")
    ap.add_argument("--digits", type=int, default=30)
    ap.add_argument("--nproc", type=int, default=1)
    ap.add_argument("--routes", default="A,B")
    ap.add_argument("--prec0", type=int, default=384)
    ap.add_argument("--nf", default="5")
    ap.add_argument("--assemble-only", dest="assemble_only", action="store_true", help="do not compute, only assemble cached results")
    ap.add_argument("--speed", default=None, help="force level class fast|gl|slow|slow2 for all channels")
    ap.add_argument("--levelA", default=None, help="outer,inner DE steps for route A, e.g. '1/16,1/32'")
    ap.add_argument("--levelB", default=None, help="same for route B (default 1/15,1/30)")
    ap.add_argument("--out", default=None, help="write result JSON here (default oracle_ckpt/<key>_result.json)")
    a = ap.parse_args()
    for r, lv in (("A", a.levelA), ("B", a.levelB)):
        if lv:
            ho, hi = lv.split(","); OVERRIDE[r] = (F(ho), F(hi))
    if a.speed:
        for ch_ in CHANNEL_SPEED:
            CHANNEL_SPEED[ch_] = a.speed
    split = None
    if a.shape.replace(" ", "").startswith("split:"):
        # split-signature point of the frame w=(u,1,0,v), wb=(ub,1,0,vb): 'split:u=..;ub=..;v=..;vb=..' (u,ub,v,vb in Q)
        kv = dict(p.split("=") for p in a.shape.replace(" ", "")[6:].split(";") if p)
        split = tuple(F(kv[k]) for k in ("u", "ub", "v", "vb"))
        stxt = "split:u=%s;ub=%s;v=%s;vb=%s" % split
        skey = hashlib.sha256(stxt.encode()).hexdigest()[:12]
    else:
        u, v = parse_shape(a.shape)
        stxt, skey = shape_key(u, v)
    if a.channel == "q_jet":
        chans = list(OC.QUARK_JET)
    elif a.channel == "g_jet":
        chans = list(OC.GLUON_JET)
    else:
        chans = a.channel.split(",")
    if split is None:
        pts = CO.shape_points(u, v)
        d = CO.dmatrix(pts)
    else:
        pts = []
        d = CO.dmatrix_split(*split)
    if any(x <= 0 for x in d.values()):
        sys.exit("degenerate/non-Euclidean shape: some d_ab <= 0: %s" % d)
    xL = max(d.values())
    classes = CO.sigma_classes(d)
    gate_prov.wdir("oracle_ckpt")
    ck = os.path.join(CKDIR, "%s.jsonl" % skey)
    done = load_ckpt(ck)
    routes = a.routes.split(",")
    tasks = []
    for ch in chans:
        for (sig, mult, dsig) in classes:
            for r in routes:
                RP = route_params(r, ch)
                key = (ch, tuple(sig), r, str(RP["ho"]), str(RP["hi"]))
                if key in done:
                    continue
                tasks.append((ch, tuple(sig), mult, [[list(k), str(vv)] for k, vv in dsig.items()], r, a.digits, a.prec0, (str(RP["ho"]), str(RP["hi"]))))
    print("shape %s [key %s]: %d sigma-classes, channels %s, routes %s; %d tasks to run (%d cached)"
          % (stxt, skey, len(classes), chans, routes, len(tasks), len(done)), flush=True)
    T0 = time.time()
    if a.assemble_only:
        tasks = []
    if tasks:
        # big channels first for load balance
        order = {"g_gggg": 0, "g_qbggq": 1, "q_qbqgq": 2, "g_qbqqbq": 3, "q_gggq": 4, "g_qbpqpqbq": 5, "q_qbpqpgq": 6}
        tasks.sort(key=lambda t: order.get(t[0], 9))
        def sink(r):
            with open(ck, "a") as f:
                f.write(json.dumps(r) + "\n")
            print("  done %s sigma=%s route %s  I=%s  quad %.0fs build %.0fs maxprec %d" %
                  (r["channel"], r["sigma"], r["route"], r["value_ball"][:42], r["t_quad"], r["t_build"], r["stats"]["max_prec"]), flush=True)
        if a.nproc > 1:
            with mp.Pool(a.nproc) as pool:
                for r in pool.imap_unordered(task_run, tasks):
                    sink(r)
        else:
            for t in tasks:
                sink(task_run(t))
    done = load_ckpt(ck)
    # assemble
    ctx.prec = a.prec0
    nf = F(a.nf)
    out = {"shape": stxt, "shape_key": skey, "points": [[str(p[0]), str(p[1])] for p in pts],
           "d": {"%d%d" % k: str(vv) for k, vv in d.items()}, "xL": str(xL), "n_sigma_classes": len(classes),
           "nf": str(nf), "digits_target": a.digits, "routes": {r: {k: str(vv) for k, vv in ROUTES[r].items() if k not in ("ho", "hi")} for r in routes},
           "levels_per_channel": {ch: {r: [str(x) for x in (route_params(r, ch)["ho"], route_params(r, ch)["hi"])] for r in routes} for ch in chans},
           "channels": {}, "jets": {}, "script": os.path.abspath(__file__)}
    pref_x = arb(xL.numerator) / arb(xL.denominator)
    jet_tot = {}
    ch_err = {}
    SYMS = load_symmetries()
    for ch in chans:
        w = OC.flavor_weight(ch, nf) * OC.SFAC[ch] * xL ** 3
        wa = arb(w.numerator) / arb(w.denominator)
        per_route = {}
        tq = {}
        # assemble every route that is fully cached for this channel (in the order A, B); the first is primary
        avail = []
        for r in ("A", "B"):
            RP = route_params(r, ch)
            if all((ch, tuple(sig), r, str(RP["ho"]), str(RP["hi"])) in done for (sig, mult, dsig) in classes):
                avail.append(r)
        if not avail:
            print("channel %s incomplete" % ch); continue
        routes = avail
        for r in routes:
            tot = arb(0); tqq = 0.0; tb = 0.0
            RP = route_params(r, ch)
            for (sig, mult, dsig) in classes:
                rec = done[(ch, tuple(sig), r, str(RP["ho"]), str(RP["hi"]))]
                tot += mult * _ball(rec)
                tqq += rec["t_quad"]; tb += rec["t_build"] + rec["t_load"]
            per_route[r] = wa * tot
            tq[r] = {"quad_s": round(tqq, 1), "build_s": round(tb, 1), "levels": [str(RP["ho"]), str(RP["hi"])], "roles": {k: RP[k] for k in ("g", "o", "i", "a")}}
        rec = {"per_route": {r: per_route[r].str(50, radius=True) for r in per_route}, "timing": tq,
               "flavor_weight_x_sfac_x_xL3": str(w)}
        errs = []
        if len(per_route) >= 2:
            A, B = per_route[routes[0]], per_route[routes[1]]
            diff = abs(A - B); rel = diff / abs(A)
            rec["abs_diff"] = diff.str(5); rec["rel_diff"] = rel.str(5)
            rec["certified_digits_routes"] = _digits(rel)
            errs.append(diff)
        symc = sym_certificate(ch, classes, done, routes[0], route_params(routes[0], ch), SYMS.get(ch, []))
        if symc is not None:
            e_abs, npairs, maxrel = symc
            e_abs = wa * e_abs
            rec["sym_certificate"] = {"route": routes[0], "taus": [list(t) for t in SYMS.get(ch, [])], "n_pairs": npairs,
                                      "sum_mult_absdiff": e_abs.str(5), "max_pair_rel_diff": maxrel.str(5),
                                      "certified_digits_sym": _digits(e_abs / abs(per_route[routes[0]]))}
            errs.append(e_abs)
        if errs:
            # the certificate reported: route-vs-route when available, else the symmetry-partner bound
            e = errs[0]
            rec["certified_abs_err"] = e.str(5)
            rec["certified_digits"] = _digits(e / abs(per_route[routes[0]]))
            rec["certificate_kind"] = "routes A-vs-B" if len(per_route) >= 2 else "symmetry partners (route %s)" % routes[0]
            ch_err[ch] = e
        rec["value"] = per_route[routes[0]].mid().str(45, radius=False)
        rec["primary_route"] = routes[0]
        out["channels"][ch] = rec
        jet_tot.setdefault("primary", arb(0)); jet_tot["primary"] += per_route[routes[0]]
        for r in per_route:
            jet_tot.setdefault(r, arb(0))
            jet_tot[r] += per_route[r]
    jetname = {"q": "q_jet", "g": "g_jet"}
    if set(chans) == set(OC.QUARK_JET) or set(chans) == set(OC.GLUON_JET):
        jn = "q_jet" if set(chans) == set(OC.QUARK_JET) else "g_jet"
        if all(c in out["channels"] for c in chans):
            common = [r for r in ("A", "B") if all(r in out["channels"][c]["per_route"] for c in chans)]
            rec = {"per_route": {r: jet_tot[r].str(50, radius=True) for r in ["primary"] + common},
                   "primary_routes": {c: out["channels"][c]["primary_route"] for c in chans}}
            if all(len(out["channels"][c]["per_route"]) >= 2 for c in chans):
                A, B = jet_tot["A"], jet_tot["B"]
                rel = abs(A - B) / abs(A)
                rec["rel_diff_routes"] = rel.str(5); rec["certified_digits_routes"] = _digits(rel)
            if all(c in ch_err for c in chans):
                etot = arb(0)
                for c in chans:
                    etot += ch_err[c]
                rec["certified_abs_err"] = etot.str(5)
                rec["certified_digits"] = _digits(etot / abs(jet_tot["primary"]))
                rec["certificate_kinds"] = {c: out["channels"][c]["certificate_kind"] for c in chans}
            rec["value"] = jet_tot["primary"].mid().str(45, radius=False)
            rec["wall_total_quad_s"] = round(sum(t["quad_s"] for c in chans for t in out["channels"][c]["timing"].values()), 1)
            out["jets"][jn] = rec
    out["stamp_utc"] = subprocess.check_output(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip()
    out["producer"] = gate_prov.producer(inputs=[ck] + [os.path.join(GATE, "cells", "%s_g%d.pkl" % (c, g)) for c in chans for g in (4, 1)])
    out["wall_this_invocation_s"] = round(time.time() - T0, 1)
    op = a.out or os.path.join(CKDIR, "%s_%s_result.json" % (skey, a.channel.replace(",", "_")))
    json.dump(out, open(op, "w"), indent=1)
    # print summary
    print("\n=== E4C LO collinear oracle: shape %s (points 0, 1, u, v), nf=%s" % (stxt, nf))
    for ch, rec in out["channels"].items():
        print("%-11s value %s" % (ch, rec["value"]))
        for r, vv in rec["per_route"].items():
            print("            route %s %s" % (r, vv))
        if "rel_diff" in rec:
            print("            |A-B|/|A| = %s  -> digits (routes): %d" % (rec["rel_diff"], rec["certified_digits_routes"]))
        if "sym_certificate" in rec:
            sc = rec["sym_certificate"]
            print("            symmetry partners: %d pairs, max pair rel diff %s -> digits (sym): %d" % (sc["n_pairs"], sc["max_pair_rel_diff"], sc["certified_digits_sym"]))
        if "certified_digits" in rec:
            print("            CERTIFIED DIGITS: %d  [%s]" % (rec["certified_digits"], rec["certificate_kind"]))
    for jn, rec in out["jets"].items():
        print("%-11s value %s" % (jn, rec["value"]))
        for r, vv in rec["per_route"].items():
            print("            route %s %s" % (r, vv))
        if "rel_diff_routes" in rec:
            print("            |A-B|/|A| = %s  -> digits (routes): %d" % (rec["rel_diff_routes"], rec["certified_digits_routes"]))
        if "certified_digits" in rec:
            print("            CERTIFIED DIGITS: %d" % rec["certified_digits"])
    print("result JSON:", op)


def _ball(rec):
    """stored (mid string, radius string) -> arb ball"""
    return arb(rec["value_mid"]) + arb(0, arb(rec["value_rad"]))


def _digits(rel):
    import math
    try:
        up = float(rel.upper()) if hasattr(rel, "upper") else float(rel.mid()) + float(rel.rad())
    except Exception:
        up = float(str(rel.mid())) + float(str(rel.rad()))
    if up <= 0:
        return 99
    return int(math.floor(-math.log10(up)))


if __name__ == "__main__":
    main()

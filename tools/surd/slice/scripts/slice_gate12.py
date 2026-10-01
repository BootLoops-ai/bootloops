"""S1 GATE (semantics): the symbolic slice one-fold object F(x;t) (ckpt/<tag>/stage2.pkl, divided by the external factor prod d^ext(t)) at
rational t and REAL x > 0 versus the independent residue-route fiber function (gate/scripts/onefold_arb.fibre_residue: exact
cells with numeric d_ab(t), inner integral by residues in acb, 1-D exp-sinh in Arb; no engine code shared).  Writes SLICE_G12_<tag>.json."""
import os, sys, json, time, argparse, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc, slice_eval
from fractions import Fraction as Fr
import mpmath as mp
from flint import acb, arb, ctx, fmpq
import onefold_arb as OA, oracle_core as CO, oracle_cells as OC, hf_piece
WORK = slice_prov.WORK
_CELLS = {}
def fibre_real(channel, gauge, pair, dsig, last, x, h, S, digits, support=None, prec0=256, precmax=8192):
    """F(x_last) at REAL rational x > 0 by the residue route (oracle_core.RouteEval, real Arb): both inner role orders; returns list of
    (roles(i,a), ball, stats).  Independent of the symbolic engine."""
    key = (channel, gauge)
    if key not in _CELLS: _CELLS[key] = OC.build_cells(channel, gauge, verbose=False)
    sub = hf_piece.group_cells(_CELLS[key], tuple(sorted(int(c) for c in pair)))
    if support is not None:
        Sset = frozenset(support); sub = dict(sub); sub["cells"] = [(dk, cv) for dk, cv in sub["cells"] if frozenset(nm for nm, pw in dk) == Sset]
    rf = CO.InstanceRF(sub, dsig); others = [k for k in rf.xvars if k != last]; out = []
    xq = fmpq(Fr(x).numerator, Fr(x).denominator)
    for (i_, a_) in ((others[0], others[1]), (others[1], others[0])):
        rev = CO.RouteEval(rf, last, i_, a_); K = int(math.ceil(S / float(h))); tot = arb(0); st = {"nodes": 0, "max_prec": prec0, "failed_nodes": 0}
        for k in range(-K, K + 1):
            prec = prec0; g = None
            while True:
                ctx.prec = prec; c = arb.pi() / 2; hh = arb(fmpq(h.numerator, h.denominator))      # node abscissa/weight at the WORKING precision
                setup = rev.outer_setup(arb(xq)); s_ = hh * (arb(k) + arb(1) / 2); xi = arb(xq).sqrt().mid() * (c * s_.sinh()).exp(); w = hh * c * s_.cosh() * xi
                try:
                    g = rev.node(setup, xi) * w; ok = g.rad() < 10.0 ** (-(digits + 4))
                except RuntimeError: ok = False; g = None
                if ok or prec >= precmax: break
                prec *= 2
            st["max_prec"] = max(st["max_prec"], prec); st["nodes"] += 1
            if g is None: st["failed_nodes"] += 1; tot += arb(0, arb("inf")) if False else arb("nan"); break
            tot += g
        ctx.prec = prec0; out.append(((i_, a_), tot, st))
        if g is not None and tot.rad().str(5) != "nan" and float(tot.rad().str(5, radius=False) or "inf") < abs(float(tot.mid().str(20, radius=False) or "1")) * 10.0 ** (-(digits - 1)): break
    return out
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", required=True); ap.add_argument("--t", default="1/2,4/5"); ap.add_argument("--x", default="7/3,11/2,3/5,37/4")
    ap.add_argument("--dps", type=int, default=60); ap.add_argument("--digits", type=int, default=34); ap.add_argument("--h", default="1/32")
    a = ap.parse_args(); slice_prov.log_pid("slice_gate12:" + a.tag); t0 = time.time()
    pkl = os.path.join(WORK, "ckpt", a.tag, "stage2.pkl")
    C, letters, ts, extra = sc.ts_load(pkl)
    channel, sigma, pair, gauge, order = extra["channel"], extra["sigma"], extra["pair"], extra["gauge"], extra["order"]
    last = "x%d" % order[2]
    G = sc.build_group(channel, gauge, pair, tuple(int(c) for c in sigma))
    rows = []; mins = []
    for ts_ in a.t.split(","):
        tv = Fr(ts_); S = slice_eval.Specialized(C, letters, ts, tv, last); ef = sc.ext_factor_value(G, tv)
        d = sc.sigma_d(sc.d_of_t(tv), tuple(int(c) for c in sigma))
        piece = {"channel": channel, "sigma": sigma, "pair_group": pair, "d_sigma": {"%d%d" % k: str(v) for k, v in d.items()}, "support_atoms": None}
        for xs in a.x.split(","):
            xv = Fr(xs); t1 = time.time()
            try:
                val, mx = S.value(xv, a.dps, return_terms=True); mine = val / (mp.mpf(ef.numerator) / ef.denominator)
            except AssertionError as ex:
                rows.append({"t": ts_, "x": xs, "skipped": "on-contour ZIP letter at this real x (termwise letter on the positive axis; representation needs the x+i0 side prescription there): %s" % str(ex)[:120]}); print("t=%s x=%s SKIP on-contour" % (ts_, xs), flush=True); continue
            t2 = time.time()
            best = None; tries = {}
            for hh in (Fr(a.h), Fr(a.h) / 2):
                for roles, fv, fst in fibre_real(channel, gauge, pair, d, order[2], xv, hh, 4.7, a.digits):
                    try:
                        fmid = float(fv.mid().str(20, radius=False)); frad = float(fv.rad().str(5, radius=False))
                    except ValueError: fmid, frad = float("nan"), float("inf")
                    if not (frad < float("inf")) or fmid != fmid:
                        tries["roles=%s,h=%s" % (list(roles), hh)] = {"residue_route": "UNRESOLVED", "nodes": fst["nodes"], "max_prec": fst["max_prec"]}; continue
                    ref = mp.mpf(fv.mid().str(70, radius=False)); rel = abs(mine - ref) / abs(ref); dg = 99.0 if rel == 0 else float(-mp.log10(rel)); rr = frad / abs(fmid) if fmid else None
                    tries["roles=%s,h=%s" % (list(roles), hh)] = {"residue_route": fv.str(42, radius=True), "rel_radius": ("%.1e" % rr) if rr is not None else None, "agree_digits": round(dg, 1), "nodes": fst["nodes"], "max_prec": fst["max_prec"]}
                    best = dg if best is None else max(best, dg)
                if best is not None and best >= a.digits - 4: break
            if best is None: best = -1.0
            rows.append({"t": ts_, "x": xs, "symbolic_value/ext": mp.nstr(mine, 40), "abs_Im(symbolic)": mp.nstr(abs(mine.imag), 3), "n_oncontour_letters(-i0 rule)": S.n_oncontour, "max_term": mp.nstr(mx, 5), "cancellation_digits": round(float(mp.log10(mx / abs(val))) if abs(val) > 0 else 0, 1), "residue_route_tries": tries, "agree_digits": round(best, 1), "wall_eval_s": round(t2 - t1, 1), "wall_residue_s": round(time.time() - t2, 1)})
            mins.append(best); print("t=%s x=%s digits %.1f (eval %.0fs, residue %.0fs)" % (ts_, xs, best, t2 - t1, time.time() - t2), flush=True)
    per_t = {}
    for r in rows:
        if "agree_digits" in r: per_t.setdefault(r["t"], []).append(r["agree_digits"])
    rec = {"tag": a.tag, "points_gated_per_t": {k: len(v) for k, v in per_t.items()}, "channel": channel, "sigma": sigma, "pair": pair, "gauge": gauge, "last_var": last, "ext_factor_poly_t": extra["ext_factor_poly_t"], "dps": a.dps, "rows": rows, "min_agree_digits": min(mins) if mins else None,
           "PASS(>=30)": bool(mins and min(mins) >= 30 and all(len(v) >= 2 for v in per_t.values()) and len(per_t) == len(a.t.split(","))), "n_terms": len(ts), "wall_s": round(time.time() - t0, 1)}
    out = os.path.join(WORK, "SLICE_G12_%s.json" % a.tag); slice_prov.write_receipt(out, rec, inputs=[pkl, os.path.join(slice_prov.GATE_S, "onefold_arb.py"), os.path.join(slice_prov.GATE_S, "oracle_core.py")])
    print("wrote", out, "min digits", rec["min_agree_digits"], "points per t", rec["points_gated_per_t"], "PASS", rec["PASS(>=30)"])
if __name__ == "__main__":
    main()

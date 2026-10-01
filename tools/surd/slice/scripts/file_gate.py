"""FILE gate: the DELIVERABLE files out/E4C_LO_{QCD,N4}_dipole_slice_*.json.gz evaluated by assemble.evaluate_file (independent of the
group pickles: canonical keys, roots by polyroots, hyperlogs by hpath with the -i0 rule) at every t for which an S0-oracle record exists
(SLICE_S0_ORACLE_t*.json from slice_oracle.py: certified two-route Arb values) and, optionally (--e4c), versus an external float engine e4c.py
found under $SURD_E4C_ENGINE (Gauss n=16,24; its xL^3-normalized value = max(1,t^2)^3 F; not part of this package).
Oracle points T > 1 are compared through the exact relabeling symmetry F(T) = s^6 F(s), s = 1/T (the closed forms are constructed on 0<s<1).
Per-(file,t) values are cached in ckpt/file_eval/<label>__t<p>_<q>.json (resume = rerun; --shard k/n splits the (file,t) list over processes).
Writes SLICE_FILE_GATE.json with --receipt (reads all caches)."""
import os, sys, json, time, glob, argparse, resource
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc, assemble
from fractions import Fraction as Fr
import mpmath as mp
import oracle_cells as OC
WORK = slice_prov.WORK; DATA = slice_prov.DATA; OUT = os.path.join(DATA, "out"); CK = os.path.join(WORK, "ckpt", "file_eval")
FILES = {"quark": "E4C_LO_QCD_dipole_slice_quark.json.gz", "gluon": "E4C_LO_QCD_dipole_slice_gluon.json.gz", "n4": "E4C_LO_N4_dipole_slice_n4.json.gz"}
for ch_ in OC.QUARK_JET + OC.GLUON_JET: FILES[ch_] = "E4C_LO_QCD_dipole_slice_%s.json.gz" % ch_
def tkey(t): t = Fr(t); return "t%s_%s" % (t.numerator, t.denominator)
def cache_path(lab, s): return os.path.join(CK, "%s__%s.json" % (lab, tkey(s)))
def oracle_table():
    orc = {}
    for f in slice_prov.rglob("SLICE_S0_ORACLE_t*.json"):          # S0 oracle receipts: WORK if present there, else DATA
        d = json.load(open(f)); orc[Fr(d["t"])] = (d, f)
    return orc
def eval_points(orc):
    """closed-form evaluation points s in (0,1]: s = t for oracle t <= 1, s = 1/t for oracle t > 1"""
    pts = {}
    for t in orc:
        s = t if t <= 1 else 1 / t; pts.setdefault(s, []).append(t)
    return pts
def wpts(t):
    t = float(t); w1 = (t / 2) * complex(3, 4) / 5; return [w1, -w1, complex(-0.5, 0), complex(0.5, 0)]
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--labels", default="quark,gluon,n4"); ap.add_argument("--dps", type=int, default=60); ap.add_argument("--shard", default="0/1")
    ap.add_argument("--s", default="", help="comma list of evaluation points s (default: all oracle points)"); ap.add_argument("--receipt", action="store_true"); ap.add_argument("--e4c", action="store_true", help="also run the external float engine $SURD_E4C_ENGINE/e4c.py at oracle t (cached ckpt/file_eval/e4c.json)")
    ap.add_argument("--e4c-t", dest="e4ct", default="", help="restrict e4c to these oracle t"); ap.add_argument("--e4c-n", dest="e4cn", default="16,24"); ap.add_argument("--e4c-n-gluon", dest="e4cng", default="16"); ap.add_argument("--reverse", action="store_true")
    a = ap.parse_args(); slice_prov.log_pid("file_gate"); os.makedirs(CK, exist_ok=True); orc = oracle_table(); pts = eval_points(orc)
    if a.s: pts = {Fr(x): pts.get(Fr(x), [Fr(x)]) for x in a.s.split(",")}
    k, n = [int(x) for x in a.shard.split("/")]
    jobs = [(lab, s) for lab in a.labels.split(",") if lab in FILES for s in sorted(pts)]; jobs = [j for i, j in enumerate(jobs) if i % n == k]
    if a.reverse: jobs = jobs[::-1]
    if not a.receipt:
        docs = {}
        for lab, s in jobs:
            cp = cache_path(lab, s)
            p = os.path.join(OUT, FILES[lab])
            if not os.path.exists(p): print("missing file", p, flush=True); continue
            if os.path.exists(cp):
                cprev = json.load(open(cp))
                if cprev.get("dps", 0) >= a.dps and cprev.get("file_sha256") == slice_prov.sha(p): print("have", lab, s, flush=True); continue
            if lab not in docs:
                t0 = time.time(); import gzip; docs[lab] = json.load(gzip.open(p, "rt")); print("loaded %s (%d terms) %.0fs RSS %.0f MB" % (lab, len(docs[lab]["terms"]), time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), flush=True)
            doc = docs[lab]; t0 = time.time(); v, byw = assemble.evaluate_file(p, s, a.dps, by_weight=True, doc=doc); wall = time.time() - t0
            rec = {"label": lab, "file": p, "file_sha256": slice_prov.sha(p), "complete": doc.get("complete"), "s": str(s), "dps": a.dps, "value": mp.nstr(v.real, 50), "abs_Im": mp.nstr(abs(v.imag), 5),
                   "by_weight(Re,|Im|)": {str(w): [mp.nstr(x.real, 30), mp.nstr(abs(x.imag), 5)] for w, x in sorted(byw.items())}, "wall_s": round(wall, 1), "peak_rss_MB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1), "stamp_utc": slice_prov.stamp(),
                   "n_terms": len(doc["terms"])}
            tmp = cp + ".tmp%d" % os.getpid(); json.dump(rec, open(tmp, "w"), indent=1); os.replace(tmp, cp)
            print("%s s=%s value %s |Im| %s (%.0fs)" % (lab, s, mp.nstr(v.real, 30), mp.nstr(abs(v.imag), 3), wall), flush=True)
    if a.e4c:
        eng = os.environ.get("SURD_E4C_ENGINE") or sys.exit("file_gate.py --e4c needs SURD_E4C_ENGINE (directory of the external float engine e4c.py)"); sys.path.insert(0, eng); import e4c
        ep = os.path.join(CK, "e4c.json"); have = json.load(open(ep)) if os.path.exists(ep) else {}
        e4ts = [Fr(x) for x in a.e4ct.split(",")] if a.e4ct else sorted(orc)
        for t in e4ts:
            for ch in OC.QUARK_JET + OC.GLUON_JET:
                for nn in [int(x) for x in (a.e4cn if ch in OC.QUARK_JET else a.e4cng).split(",")]:
                    key = "%s|%s|n%d" % (t, ch, nn)
                    if key in have: continue
                    t0 = time.time(); val = e4c.G4_channel(wpts(t), ch, n=nn, nf=5.0)
                    have[key] = {"t": str(t), "channel": ch, "n": nn, "G4_channel(e4c, xL^3 convention)": float(val), "wall_s": round(time.time() - t0, 1), "stamp_utc": slice_prov.stamp()}
                    tmp = ep + ".tmp%d" % os.getpid(); json.dump(have, open(tmp, "w"), indent=1); os.replace(tmp, ep); print(key, float(val), "%.0fs" % (time.time() - t0), flush=True)
    if a.receipt:
        mp.mp.dps = 60; out = {"method": __doc__, "per_oracle_t": {}, "summary": {}}
        e4 = json.load(open(os.path.join(CK, "e4c.json"))) if os.path.exists(os.path.join(CK, "e4c.json")) else {}
        mins = {}
        for t in sorted(orc):
            d, rf = orc[t]; s = t if t <= 1 else 1 / t; fac = Fr(1) if t <= 1 else s ** 6; row = {"s(evaluation point)": str(s), "symmetry_factor s^6": str(fac), "oracle_receipt": rf, "functions": {}}
            xL = max(sc.d_of_t(t).values())
            for lab in FILES:
                cp = cache_path(lab, s)
                if not os.path.exists(cp): continue
                c = json.load(open(cp)); val = mp.mpf(c["value"]) * fac
                rec = {"closed_form(file) x factor": mp.nstr(val, 42), "abs_Im": c["abs_Im"], "complete": c["complete"], "file_sha256": c["file_sha256"][:16], "eval_wall_s": c["wall_s"], "dps": c["dps"], "n_terms": c["n_terms"]}
                if lab in ("quark", "gluon"):
                    oj = d["jets"].get("q_jet" if lab == "quark" else "g_jet")
                    if oj: ref = mp.mpf(oj["value"]); rec["oracle"] = oj["value"]; rec["oracle_certified_digits"] = oj.get("certified_digits")
                    else: ref = None
                    chs = OC.QUARK_JET if lab == "quark" else OC.GLUON_JET
                    for nn in (16, 24):
                        if all("%s|%s|n%d" % (t, ch, nn) in e4 for ch in chs):
                            ev = sum(e4["%s|%s|n%d" % (t, ch, nn)]["G4_channel(e4c, xL^3 convention)"] for ch in chs) / float(xL) ** 3
                            rec["e4c.py n=%d (/xL^3)" % nn] = "%.10g" % ev; rec["e4c n=%d rel.dev" % nn] = "%.3g" % float(abs(val - ev) / abs(val)) if val != 0 else None
                else:
                    oc = d["channels"].get(lab)
                    if oc: ref = mp.mpf(oc["value"]); rec["oracle"] = oc["value"]; rec["oracle_certified_digits"] = oc.get("certified_digits")
                    else: ref = None
                    for nn in (16, 24):
                        key = "%s|%s|n%d" % (t, lab, nn)
                        if key in e4: ev = e4[key]["G4_channel(e4c, xL^3 convention)"] / float(xL) ** 3; rec["e4c.py n=%d (/xL^3)" % nn] = "%.10g" % ev; rec["e4c n=%d rel.dev" % nn] = "%.3g" % float(abs(val - ev) / abs(val)) if val != 0 else None
                if ref is not None:
                    rel = abs(val - ref) / abs(ref); rec["agree_digits_vs_oracle"] = round(float(-mp.log10(rel)), 1) if rel > 0 else 99.0
                    if c["complete"]: mins[lab] = min(mins.get(lab, 999), rec["agree_digits_vs_oracle"])
                row["functions"][lab] = rec
            out["per_oracle_t"][str(t)] = row
            print("t=%s (s=%s):" % (t, s), {lab: (r.get("agree_digits_vs_oracle"), r["complete"]) for lab, r in row["functions"].items()}, flush=True)
        out["summary"] = {"min_agree_digits_vs_oracle(complete files only)": mins, "PASS(>=30 at every oracle t, complete)": {lab: (v >= 30) for lab, v in mins.items()}, "n_oracle_t": len(orc), "oracle_t": [str(t) for t in sorted(orc)]}
        print(out["summary"])
        slice_prov.write_receipt(os.path.join(WORK, "SLICE_FILE_GATE.json"), out, inputs=slice_prov.rglob("SLICE_S0_ORACLE_t*.json") + [os.path.join(OUT, f) for f in FILES.values() if os.path.exists(os.path.join(OUT, f))] + sorted(glob.glob(os.path.join(CK, "*.json"))))
if __name__ == "__main__":
    main()

"""The COMPLETE jet closed forms (the deliverable files) evaluated at t* -+ h, h = 1e-3 and 1e-5, on both sides of every Table-1 root t* in (0,1)
(values from file_gate.py's cache ckpt/file_eval; assemble.evaluate_file, dps 60).  For a function analytic at t*: odd part D(h) = F(t*+h) - F(t*-h) = 2F'h + O(h^3)
(so D(1e-3)/D(1e-5) = 100), even part E(h) = (F(t*+h) + F(t*-h))/2 = F(t*) + F''h^2/2 (so |E(1e-3) - E(1e-5)|/|F| ~ 1e-6 x O(1)); a (t-t*)^(-3/2), (t-t*)^(-1/2)
or (t-t*)^(1/2) component on this sheet would instead give a one-sided imaginary part, D ratio 10 or blow-up, and |E(1e-3)-E(1e-5)|/|F| ~ 3e-2.
Writes SLICE_S5_NEARROOT.json."""
import os, sys, json, glob, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, letters as LT
from fractions import Fraction as Fr
import mpmath as mp, sympy as sp
WORK = slice_prov.WORK; CK = os.path.join(WORK, "ckpt", "file_eval"); T = sp.Symbol("t"); mp.mp.dps = 60
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--labels", default="quark,gluon"); a = ap.parse_args(); slice_prov.log_pid("near_root")
    roots = []
    for nm, pol in LT.TABLE1.items():
        for r in sp.Poly(pol, T).nroots(n=30):
            if abs(sp.im(r)) < 1e-20 and 0 < sp.re(r) < 1: roots.append((nm, Fr(str(sp.re(r))[:20]) if "." in str(sp.re(r)) else Fr(sp.re(r))))
    out = {"method": __doc__, "functions": {}}
    for lab in a.labels.split(","):
        cache = {}
        for f in glob.glob(os.path.join(CK, "%s__t*.json" % lab)):
            d = json.load(open(f)); cache[Fr(d["s"])] = d
        rows = []
        for nm, r in roots:
            rec = {"table1": nm, "t*": "%.12f" % float(r), "points": {}}
            vals = {}
            for h in (Fr(1, 1000), Fr(1, 100000)):
                for sg in (1, -1):
                    x = Fr(round((r + sg * h) * 10 ** 12), 10 ** 12); c = cache.get(x)
                    if c: vals[(h, sg)] = mp.mpf(c["value"]); rec["points"]["%s%s" % ("+" if sg > 0 else "-", h)] = {"t": str(x), "F": c["value"][:36], "abs_Im": c["abs_Im"], "wall_s": c["wall_s"], "file_sha256": c["file_sha256"][:12]}
            if len(vals) == 4:
                h1, h2 = Fr(1, 1000), Fr(1, 100000)
                D1 = vals[(h1, 1)] - vals[(h1, -1)]; D2 = vals[(h2, 1)] - vals[(h2, -1)]; E1 = (vals[(h1, 1)] + vals[(h1, -1)]) / 2; E2 = (vals[(h2, 1)] + vals[(h2, -1)]) / 2
                Fs = abs(E2); rec.update({"D(1e-3)": mp.nstr(D1, 12), "D(1e-5)": mp.nstr(D2, 12), "D_ratio(expect 100 if analytic)": mp.nstr(D1 / D2, 8) if D2 != 0 else None, "F'(t*) estimate": mp.nstr(D2 / (2 * h2), 15),
                            "E(1e-3)": mp.nstr(E1, 20), "E(1e-5)=F(t*) to ~1e-10": mp.nstr(E2, 20), "|E(1e-3)-E(1e-5)|/|F| (expect ~1e-6 x F''/2F if analytic, ~3e-2 for a sqrt-type cusp)": mp.nstr(abs(E1 - E2) / Fs, 5),
                            "max abs_Im over the 4 points": max(float(mp.mpf(p["abs_Im"])) for p in rec["points"].values())})
                ratio = float(D1 / D2) if D2 != 0 else 0.0; ed = float(abs(E1 - E2) / Fs)
                rec["verdict"] = "ANALYTIC-consistent (C^1: odd part linear in h to %.2g%%, even part quadratic: rel. change %.1e; Im = 0)" % (abs(ratio / 100 - 1) * 100, ed) if (abs(ratio / 100 - 1) < 0.05 and ed < 1e-3) else "NOT analytic-consistent: D ratio %s, even-part change %.2e" % (mp.nstr(D1 / D2, 6) if D2 != 0 else None, ed)
            elif len(vals) == 0: rec["verdict"] = "not probed: this Table-1 letter does not occur in the %s representation (letters table: no X-letter class, no denominator)" % lab
            else: rec["verdict"] = "incomplete (%d of 4 points evaluated)" % len(vals)
            rows.append(rec); print(lab, nm, rec["t*"][:8], rec.get("verdict"), flush=True)
        out["functions"][lab] = rows
    slice_prov.write_receipt(os.path.join(WORK, "SLICE_S5_NEARROOT.json"), out, inputs=sorted(glob.glob(os.path.join(CK, "*.json"))))
if __name__ == "__main__":
    main()

"""Guarded command-line evaluator of the assembled slice closed forms: value at rational t with an honest precision report.
  python3 evaluate.py --file quark --t 3/5 [--dps 60] [--digits 30]
* calls assemble.evaluate_file(path, t, dps); reports F, |Im F| (the error monitor: the exact value is real) and the estimated digits
  min(dps - 5, -log10(|Im F|/|F|)); if fewer than --digits, re-evaluates at dps + 30 (up to twice);
* at a representation-singular rational t (assemble.evaluate_file raises ZeroDivisionError, e.g. t = 3/5 for the QCD files: a quadratic letter's
  leading X^2 coefficient and a rational-point denominator vanish while F is analytic) it evaluates at t -+ delta (delta = 10^-(digits+4))
  and returns the mean, whose error is O(delta^2 F'') + the two points' own precision, reported.
Files: labels quark, gluon, n4, q_qbpqpgq, q_qbqgq, q_gggq, g_qbpqpqbq, g_qbqqbq, g_qbggq, g_gggg or a path."""
import os, sys, json, gzip, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, assemble
from fractions import Fraction as Fr
import mpmath as mp
OUT = os.path.join(slice_prov.DATA, "out")
FILES = {"quark": "E4C_LO_QCD_dipole_slice_quark.json.gz", "gluon": "E4C_LO_QCD_dipole_slice_gluon.json.gz", "n4": "E4C_LO_N4_dipole_slice_n4.json.gz"}
for ch_ in ("q_qbpqpgq", "q_qbqgq", "q_gggq", "g_qbpqpqbq", "g_qbqqbq", "g_qbggq", "g_gggg"): FILES[ch_] = "E4C_LO_QCD_dipole_slice_%s.json.gz" % ch_
def digits_of(v, dps):
    if v == 0: return 0.0
    if v.imag == 0: return float(dps - 5)
    return float(min(dps - 5, -mp.log10(abs(v.imag) / abs(v))))
def evaluate(doc_or_path, t, dps=60, digits=30, max_escalations=2, doc=None):
    import hashlib
    def _sha(pth):
        h = hashlib.sha256()
        with open(pth, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""): h.update(chunk)
        return h.hexdigest()
    path = doc_or_path; t = Fr(t); rep = {"t": str(t), "file": path, "file_sha256": _sha(path), "script": os.path.abspath(__file__), "script_sha256": _sha(os.path.abspath(__file__)), "assemble_sha256": _sha(assemble.__file__), "attempts": []}
    for k in range(max_escalations + 1):
        d = dps + 30 * k
        try:
            t0 = time.time(); v = assemble.evaluate_file(path, t, d, doc=doc); dg = digits_of(v, d)
            rep["attempts"].append({"dps": d, "value": mp.nstr(v.real, 50), "abs_Im": mp.nstr(abs(v.imag), 5), "digits_est": round(dg, 1), "wall_s": round(time.time() - t0, 1)})
            if dg >= digits: rep.update({"value": mp.nstr(v.real, int(dg)), "digits_est": round(dg, 1), "mode": "direct"}); return rep
        except (ZeroDivisionError, AssertionError) as ex:
            rep["attempts"].append({"dps": d, "note": "direct evaluation raised %s: representation-singular t" % repr(ex)[:80]})
            # representation-singular rational t: symmetric limit (F(t+delta)+F(t-delta))/2, error ~ F'' delta^2 / 2.  delta must keep nearly-colliding
            # complex roots >= 1e-12 off the integration path (hpath's polygonal detour), so escalate delta from 10^-(digits+4) upwards until hpath accepts.
            last = None
            for dexp in (digits + 4, 24, 20, 16, 12, 8):
                delta = Fr(1, 10 ** dexp); d2 = max(d, 2 * dexp + 20)
                try:
                    vp = assemble.evaluate_file(path, t + delta, d2, doc=doc); vm = assemble.evaluate_file(path, t - delta, d2, doc=doc)
                except (ZeroDivisionError, AssertionError) as ex2: last = repr(ex2)[:120]; rep["attempts"].append({"dps": d2, "delta": str(delta), "note": "symmetric limit failed: %s" % last}); continue
                v = (vp + vm) / 2; scale = max(abs(vp), abs(vm), abs(vp - vm) / (2 * delta)); dmp = mp.mpf(delta.numerator) / delta.denominator; dg_delta = float(-mp.log10(dmp ** 2 * (1 + scale / max(abs(v), mp.mpf(10) ** -300))))
                dg = min(digits_of(vp, d2), digits_of(vm, d2), dg_delta)
                rep["attempts"].append({"dps": d2, "delta": str(delta), "note": "representation-singular t (%s): symmetric limit from t -+ delta; error <= O(F'' delta^2)" % ex.__class__.__name__, "value": mp.nstr(v.real, 50), "abs_Im": mp.nstr(max(abs(vp.imag), abs(vm.imag)), 5), "F(t+delta)-F(t-delta)": mp.nstr(vp.real - vm.real, 8), "digits_est": round(dg, 1)})
                rep.update({"value": mp.nstr(v.real, int(max(min(dg, 60), 5))), "digits_est": round(dg, 1), "mode": "symmetric limit (delta=%s)" % delta}); return rep
            rep.update({"value": None, "digits_est": 0, "mode": "failed: %s" % last}); return rep
    rep.update({"value": rep["attempts"][-1]["value"], "digits_est": rep["attempts"][-1]["digits_est"], "mode": "direct (target digits not reached; raise --dps)"}); return rep
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--file", required=True); ap.add_argument("--t", required=True); ap.add_argument("--dps", type=int, default=60); ap.add_argument("--digits", type=int, default=30)
    a = ap.parse_args(); p = a.file if os.path.exists(a.file) else os.path.join(OUT, FILES[a.file])
    doc = json.load(gzip.open(p, "rt")); r = evaluate(p, a.t, a.dps, a.digits, doc=doc); print(json.dumps(r, indent=1))

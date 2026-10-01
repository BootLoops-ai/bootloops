"""Are the RATIONAL Table-1 slice letters P(t) log-letters of an assembled slice function, per transcendental weight and per symbol slot?
For a prime P of Q[t] with a simple real root t* in (0,1) and the real place PP above it reached as t -> t*+0 (roots continued by proximity),
the functional  Psi_p = sum_T c_T(t0) * v_PP(e_{T,p}) * prod_{q != p} log|e_{T,q}(t0)|  (c_T = exact algebraic coefficient embedded at the
transcendental point t0 = e/7, e_{T,q} = the atoms of the Goncharov symbol of term T as in symtest.py, v_PP = order of vanishing at t*,
measured as the slope d log|e| / d log h over h = 1e-5, 1e-7, 1e-9, snapped to k/2 when consistent over both decades-pairs) is a linear functional on K (x) (K^x)^{(x) w} that
vanishes on every tensor whose slot p carries no PP-adic content.  Psi_p != 0  =>  log P(t) IS a letter in slot p of the weight-w symbol;
all Psi_p == 0 (to 25 digits of the scale) => absent (up to an accidental zero at t0; a second point t0' = pi/9 is also reported).
Atoms with v != 0 at a Table-1 root: rho_1 - rho_2 of the quadratic X-letters whose discriminant contains P (v = 1/2: the ramified atom that the
sqrt-odd test cannot see, its odd part being torsion), and any rational/algebraic atom that happens to vanish at t* (found numerically over ALL
atoms of ALL terms: the scan is global, not restricted to class terms).  Writes SLICE_S5_VALTEST_<label>.json."""
import os, sys, json, gzip, time, argparse, collections, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, symtest as SY, letters as LT
from fractions import Fraction as Fr
import mpmath as mp, sympy as sp
WORK = slice_prov.WORK; DATA = slice_prov.DATA; OUT = os.path.join(DATA, "out"); T = sp.Symbol("t")
FILES = {"quark": "E4C_LO_QCD_dipole_slice_quark.json.gz", "gluon": "E4C_LO_QCD_dipole_slice_gluon.json.gz", "n4": "E4C_LO_N4_dipole_slice_n4.json.gz"}
for ch_ in ("q_qbpqpgq", "q_qbqgq", "q_gggq", "g_qbpqpqbq", "g_qbqqbq", "g_qbggq", "g_gggg"): FILES[ch_] = "E4C_LO_QCD_dipole_slice_%s.json.gz" % ch_

class NumTracked(SY.Num):
    """Num at t* + h whose root labeling is matched by proximity to a reference Num (same letters at a nearby t)"""
    def __init__(self, doc, tval, dps, ref=None):
        super().__init__(doc, tval, dps); self.ref = ref
    def roots(self, L):
        if L not in self._roots:
            rr = SY.Num.roots(self, L)
            if self.ref is not None:
                ref = self.ref.roots(L); used = set(); perm = []
                def chord(a_, b_): return abs(a_ - b_) / mp.sqrt((1 + abs(a_) ** 2) * (1 + abs(b_) ** 2))     # metric on the Riemann sphere: roots escaping to oo stay matched
                for r0 in ref:
                    j = min((k for k in range(len(rr)) if k not in used), key=lambda k: chord(rr[k], r0)); used.add(j); perm.append(rr[j])
                self._roots[L] = perm
        return self._roots[L]
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--label", required=True); ap.add_argument("--dps", type=int, default=60); ap.add_argument("--weights", default="1,2,3")
    ap.add_argument("--h", default="1e-5,1e-7,1e-9"); ap.add_argument("--no-control", dest="control", action="store_false"); a = ap.parse_args(); slice_prov.log_pid("valtest:" + a.label); t00 = time.time()
    path = os.path.join(OUT, FILES.get(a.label, a.label)); doc = json.load(gzip.open(path, "rt")); mp.mp.dps = a.dps
    weights = [int(x) for x in a.weights.split(",")]; hs = [mp.mpf(x) for x in a.h.split(",")]; h1, h2, h3 = hs[0], hs[1], hs[2]
    # Table-1 primes with a real root in (0,1)
    primes = []
    for nm, pol in LT.TABLE1.items():
        P = sp.Poly(pol, T)
        for r in P.nroots(n=50):
            if abs(sp.im(r)) < 1e-40 and 0 < sp.re(r) < 1: primes.append((nm, str(P.as_expr()), mp.mpf(str(sp.re(r)))))
    if a.control: primes.append(("CONTROL(t-1): detector liveness, atoms DO degenerate at t=1", "t - 1", mp.mpf(1)))
    print("Table-1 primes with a root in (0,1):", [(nm, mp.nstr(r, 10)) for nm, P, r in primes], flush=True)
    # symbols of ALL terms per weight (cached builder); atoms collected
    SB = SY.SymBuilder(); per_w = {w: [] for w in weights}; t0 = time.time()
    for tm in doc["terms"]:
        wgt = sum(len(wz) for wz in tm["Z"]) + sum(len(c[1]) for c in tm["const"])
        if wgt not in per_w: continue
        S = {(): 1}
        for wz in tm["Z"]: S = SY.sym_mul(S, SB.sym_Z(wz))
        x0s = "%d/%d" % tuple(tm.get("x0", [11, 2]))
        for c in tm["const"]: S = SY.sym_mul(S, SB.sym_G0(c[1], x0s) if c[0] == "G0" else SB.sym_Z0(c[1], x0s))
        if S: per_w[wgt].append((tm, S))
    atoms = set()
    for w, lst in per_w.items():
        for tm, S in lst:
            for key in S:
                for e in key:
                    for at, n in e: atoms.add(at)
    print("symbols built: %s terms per weight, %d distinct atoms, %.0fs" % ({w: len(l) for w, l in per_w.items()}, len(atoms), time.time() - t0), flush=True)
    N0 = SY.Num(doc, mp.e / 7, a.dps); N0b = SY.Num(doc, mp.pi / 9, a.dps)
    def atom_val(at, N):
        if at[0] == "D": return N.point(json.loads(at[1])) - N.point(json.loads(at[2]))
        if at[0] == "P": return N.point(json.loads(at[1]))
        if at[0] == "1P": return 1 + N.point(json.loads(at[1]))
        if at[0] == "XP": return mp.mpf(Fr(at[2]).numerator) / Fr(at[2]).denominator - N.point(json.loads(at[1]))
        if at[0] == "X0": return mp.mpc(mp.mpf(Fr(at[1]).numerator) / Fr(at[1]).denominator)
        def hval(hs): h_, x0_ = json.loads(hs); return N.hl(h_, mp.mpf(Fr(x0_).numerator) / Fr(x0_).denominator)
        if at[0] == "H": return hval(at[1])
        if at[0] == "1H": return 1 + hval(at[1])
        if at[0] == "DH": return hval(at[1]) - hval(at[2])
        if at[0] == "MH": return hval(at[1])
        raise ValueError(at)
    out = {"label": a.label, "file": path, "file_sha256": slice_prov.sha(path), "complete": doc.get("complete"), "method": __doc__, "primes": {}}
    la0 = {at: mp.log(abs(atom_val(at, N0))) for at in atoms}; la0b = {at: mp.log(abs(atom_val(at, N0b))) for at in atoms}
    for nm, Pstr, tstar in primes:
        t1 = time.time()
        N1 = NumTracked(doc, tstar + h1, a.dps); N2 = NumTracked(doc, tstar + h2, a.dps, ref=N1); N3 = NumTracked(doc, tstar + h3, a.dps, ref=N2)
        val = {}; raw = {}
        for at in atoms:
            v1 = abs(atom_val(at, N1)); v2 = abs(atom_val(at, N2)); v3 = abs(atom_val(at, N3))
            if v1 == 0 or v2 == 0 or v3 == 0: raw[at] = None; continue
            s12 = (mp.log(v1) - mp.log(v2)) / (mp.log(h1) - mp.log(h2)); s23 = (mp.log(v2) - mp.log(v3)) / (mp.log(h2) - mp.log(h3)); k = mp.nint(2 * s23) / 2
            raw[at] = (float(s12), float(s23))
            # genuine zero/pole at t*: consistent power law over four decades; a near-coincidence (small nonzero limit) has s23 -> 0
            if k != 0 and abs(s23 - k) < mp.mpf("0.03") and abs(s12 - k) < mp.mpf("0.15"): val[at] = Fr(int(mp.nint(2 * s23)), 2)
        vat = sorted(((repr(at)[:160], str(v), "%.4f,%.4f" % raw[at]) for at, v in val.items()), key=lambda x: x[0])
        amb = sorted(("%.4f,%.4f" % r_, repr(at)[:120]) for at, r_ in raw.items() if r_ is not None and at not in val and abs(r_[1]) > 0.03)
        rec = {"prime": Pstr, "t*": mp.nstr(tstar, 20), "n_atoms_with_nonzero_valuation": len(val), "atoms_with_nonzero_valuation(atom, v, raw slope)": vat[:60], "atoms_with_ambiguous_slope(not within 0.02 of k/2)": amb[:20], "weights": {}}
        print("== %s (%s) t*=%s: %d atoms with v != 0, %d ambiguous (%.0fs)" % (nm, Pstr, mp.nstr(tstar, 8), len(val), len(amb), time.time() - t1), flush=True)
        for w, lst in per_w.items():
            res = {}
            for tag, N, la in (("t0=e/7", N0, la0), ("t0'=pi/9", N0b, la0b)):
                N._roots = N._roots  # keep cache
                Psi = [mp.mpc(0) for _ in range(w)]; scale = [mp.mpf(0) for _ in range(w)]; nterm = 0
                for tm, S in lst:
                    # quick skip: no atom of this symbol has v != 0
                    if not any(at in val for key in S for e in key for at, n in e): continue
                    c = mp.mpc(0)
                    for cs_ in tm["coef"]: c += N.gt(cs_)
                    for Lg, s_, ex in tm["gens"]: c *= N.roots(Lg)[s_] ** ex
                    nterm += 1
                    for key, mult in S.items():
                        logs = [sum(n * la[at] for at, n in e) for e in key]
                        for p in range(w):
                            vp = sum(n * val[at] for at, n in key[p] if at in val)
                            if vp == 0: continue
                            rest = mp.mpf(1)
                            for q in range(w):
                                if q != p: rest *= logs[q]
                            contrib = c * mult * rest; Psi[p] += contrib * (mp.mpf(vp.numerator) / vp.denominator); scale[p] = max(scale[p], abs(contrib))
                res[tag] = {"n_terms_contributing": nterm, "per_slot": {str(p + 1): {"Psi": mp.nstr(Psi[p], 10), "scale": mp.nstr(scale[p], 5), "ratio": (mp.nstr(abs(Psi[p]) / scale[p], 5) if scale[p] > 0 else None)} for p in range(w)}}
            slots_nz = sorted(set(int(p_) for tag in res for p_, x in res[tag]["per_slot"].items() if x["ratio"] is not None and float(x["ratio"]) > 1e-25))
            slots_tested = sorted(set(int(p_) for tag in res for p_, x in res[tag]["per_slot"].items() if x["ratio"] is not None))
            verdict = ("no atom with P-adic content in any weight-%d tensor: log P ABSENT trivially" % w) if not slots_tested else (("log P PRESENT in slot(s) %s" % slots_nz) if slots_nz else "log P CANCELS in every slot (tested slots %s, to >= 25 digits of scale at both t0)" % slots_tested)
            rec["weights"][str(w)] = {"functionals": res, "slots_with_content": slots_nz, "slots_tested": slots_tested, "verdict": verdict}
            print("   weight %d: %s   %s" % (w, verdict, {tag: {p_: x["ratio"] for p_, x in res[tag]["per_slot"].items()} for tag in res}), flush=True)
        out["primes"][nm] = rec
    out["wall_s"] = round(time.time() - t00, 1)
    slice_prov.write_receipt(os.path.join(WORK, "SLICE_S5_VALTEST_%s.json" % a.label), out, inputs=[path])
    print("wrote SLICE_S5_VALTEST_%s.json %.0fs" % (a.label, time.time() - t00))
if __name__ == "__main__":
    main()

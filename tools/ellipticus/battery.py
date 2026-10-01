"""ellipticus.battery — positive-control battery: reproduce reference
values through the tool's own code paths (fail-closed; receipts JSON).

Legs (reference values live in fixtures/battery/, committed with the
package; every leg checks against values that were NOT produced by the code
path under test):

L1  march engine end-to-end on the two synthetic families in
    fixtures/extended_systems.json (quartic, perfect-square fibre; cubic,
    square-x-linear fibre): polyform must EXACTLY match the pinned exact
    systems (fixtures/system_polyform.json, regenerable by polyform()),
    then moment kernels and dlog-letter words at three endpoints are
    checked against fixtures/battery/march_words.json — independent
    tanh-sinh quadrature of the defining integrals and of the documented
    letter conventions, never the march code.
L2  period/cycle-moment layer: the Legendre-K a-cycle period against the
    CLOSED FORM 2*K(1/4), and oval cycle moments against direct-quadrature
    oracles (fixtures/battery/cycle_controls.json).
L3  complex-pair contour verb: Om/X1/X2 cycle values around a conjugate
    root pair against an independent branch-cut segment-integral oracle
    (fixtures/battery/contour_controls.json).  Overall cycle orientation
    is convention: one global sign is recorded, never silently flipped.
L4  sunrise q-engine: rows 09/10/11 held-out engine oracles (bundled in
    vendor/) via qengine.sunrise_J at dps 100, + the equal-mass classical
    control.  Bar: > 95d at every gate point.
L5  ellred member: standard-elliptic reduction charts and atoms — quartic
    F/E/Pi closed atoms, cubic 3-real-root moments, analytic 1r2c reducer +
    Carlson-route agreement + derivative certification — against direct
    tanh-sinh quadrature of the defining integrals and the pinned modulus/
    chart strings in fixtures/battery/ellred_controls.json.

Reference-value routing: legs L1-L3 and L5 read their reference files under
$ELLIPTICUS_ARTIFACTS (former name $EMPL_EVAL_ARTIFACTS also read) when
set, else under fixtures/battery/ beside this
module; a leg whose files are absent prints a named SKIP and the battery
result rests on the remaining legs (L4 needs no reference files).

Run:  python3 -m ellipticus.battery      (tools/ on PYTHONPATH)
      python3 battery.py                (from this directory)
      (add --pilot for a low-dps timing pass)
Writes the receipts JSON next to the module unless ELLIPTICUS_BATTERY_OUT (or
its former name EMPL_BATTERY_OUT) is set.
"""
import json
import os
import sys
import time
from fractions import Fraction

if __name__ == "__main__" and not __package__:
    # run as a plain script (python3 battery.py): re-enter as the package
    # module so the relative imports below resolve, then stop here.
    import runpy
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    runpy.run_module("ellipticus.battery", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

import mpmath as mp
import sympy as sp

from . import QuarticCurve, MomentFamily, shared_digits, qengine, env

_PKG = os.path.dirname(os.path.abspath(__file__))
# Reference-value root for legs L1-L3 (fixtures/battery/ ships with the
# package; ELLIPTICUS_ARTIFACTS, or its former name EMPL_EVAL_ARTIFACTS,
# overrides to any root with the same layout).
_ART = env("ELLIPTICUS_ARTIFACTS", "EMPL_EVAL_ARTIFACTS",
           os.path.join(_PKG, "fixtures", "battery"))
# Family/system fixtures shared with the selftest (ELLIPTICUS_W3 / EMPL_EVAL_W3
# overrides).
_SYS = env("ELLIPTICUS_W3", "EMPL_EVAL_W3", os.path.join(_PKG, "fixtures"))

REC = {"legs": {}, "opened_utc": time.strftime("%FT%TZ", time.gmtime())}


def _have_reference(leg, *paths):
    """True when every reference file a leg needs is present; otherwise
    print the loud named SKIP (naming ELLIPTICUS_ARTIFACTS) and record it."""
    missing = [p for p in paths if not os.path.exists(p)]
    if not missing:
        return True
    REC["legs"][leg] = dict(
        status="SKIP",
        reason="reference data not present",
        missing=missing,
        cure=("set ELLIPTICUS_ARTIFACTS to a reference-data root, or place "
              "the files under " + os.path.join(_PKG, "fixtures", "battery") +
              " (layout: GUIDE.md ENV section)"))
    print(f"{leg}: SKIP — reference data not present "
          f"(missing: {', '.join(missing)}). Set ELLIPTICUS_ARTIFACTS to a "
          f"reference-data root, or place the files under "
          f"{os.path.join(_PKG, 'fixtures', 'battery')} (layout: GUIDE.md).")
    return False


def _frac(s):
    return Fraction(s)


# ------------------------------------------------------------------- L1
def leg1_march(dps=30):
    ref = os.path.join(_ART, "march_words.json")
    sysfix = os.path.join(_SYS, "extended_systems.json")
    pinfix = os.path.join(_SYS, "system_polyform.json")
    if not _have_reference("L1_march_words", ref, sysfix, pinfix):
        return True
    t0 = time.time()
    fx = json.load(open(ref))
    es = json.load(open(sysfix))
    pins = json.load(open(pinfix))
    window = (_frac(es["window"][0]), _frac(es["window"][1]))
    param = es["param"]
    endpoints = [_frac(e) for e in fx["endpoints"]]
    fams = {
        "quartic": (MomentFamily(es["P4_coeffs_x_desc"], window=window,
                                 kmax=2, param=param), 80),
        "cubic": (MomentFamily(es["cubic_coeffs_x_desc"], window=window,
                               kmax=1, param=param), 50),
    }
    bar_k = fx["bar_digits"]["kernels"]
    bar_w = fx["bar_digits"]["words"]
    rows = {}
    worst_k = worst_w = 999.0
    pin_ok = True
    for lab, (fam, extra) in fams.items():
        pf = fam.polyform()
        pin_ok &= (pf["Den"] == pins[lab]["Den"]
                   and pf["NumM"] == pins[lab]["NumM"]
                   and pf["n"] == pins[lab]["n"])
        words = [[w[0], (w[1][0], int(w[1][1]))] for w in fx["words"][lab]]
        res = fam.eval_words(words, endpoints, dps, wdps_extra=extra)
        with mp.workdps(dps + extra + 20):
            for ep in endpoints:
                o = fx["oracles"][lab][str(ep)]
                got = res[str(ep)]
                for k in range(fam.kmax + 1):
                    d = shared_digits(got["kernels"][k],
                                      mp.mpf(o["kernels"][k]))
                    rows[f"{lab}:z={ep}:I{k}"] = round(d, 1)
                    worst_k = min(worst_k, d)
                for j in range(len(words)):
                    d = shared_digits(got["words"][j], mp.mpf(o["words"][j]))
                    rows[f"{lab}:z={ep}:w{j}"] = round(d, 1)
                    worst_w = min(worst_w, d)
    ok = pin_ok and worst_k >= bar_k and worst_w >= bar_w
    REC["legs"]["L1_march_words"] = dict(
        dps=dps, rows=rows, worst_kernels=round(worst_k, 1),
        worst_words=round(worst_w, 1), bars={"kernels": bar_k, "words": bar_w},
        polyform_pins="EXACT match" if pin_ok else "MISMATCH",
        pass_=bool(ok), secs=round(time.time() - t0, 1))
    print(f"L1 march words: kernels worst {worst_k:.1f}d (bar {bar_k}), "
          f"words worst {worst_w:.1f}d (bar {bar_w}) vs independent "
          f"quadrature oracles, polyform pins "
          f"{'EXACT' if pin_ok else 'MISMATCH'} "
          f"[{'PASS' if ok else 'FAIL'}]  ({time.time()-t0:.0f}s)")
    return ok


# ------------------------------------------------------------------- L2
def leg2_cycles(dps=45):
    ref = os.path.join(_ART, "cycle_controls.json")
    if not _have_reference("L2_cycle_moments", ref):
        return True
    t0 = time.time()
    fx = json.load(open(ref))
    bar = fx["bar_digits"]
    rows = {}
    ok_all = True
    # closed-form period control
    pc = fx["period_control"]
    cur = QuarticCurve([_frac(c) for c in pc["quartic_coeffs_desc"]])
    with mp.workdps(dps + 30):
        per, meta = cur.period_oval(dps)
        d = shared_digits(per / 2, mp.mpf(pc["line"]))
    rows["period_vs_closed_form"] = round(d, 1)
    rows["period_route"] = meta.get("route")
    ok_all &= d >= bar
    # quadrature-oracle cycle moments
    mc = fx["moment_control"]
    curM = QuarticCurve([_frac(c) for c in mc["quartic_coeffs_desc"]])
    with mp.workdps(dps + 30):
        moms, metas = curM.cycle_moments((0, 1, 2), dps)
        for k in range(3):
            d = shared_digits(moms[k] / 2, mp.mpf(mc["halves"][f"I{k}"]))
            rows[f"I{k}"] = round(d, 1)
            ok_all &= d >= bar
    REC["legs"]["L2_cycle_moments"] = dict(dps=dps, rows=rows, bar=bar,
                                           pass_=bool(ok_all),
                                           secs=round(time.time() - t0, 1))
    print(f"L2 cycle layer: {rows} (bar {bar}) "
          f"[{'PASS' if ok_all else 'FAIL'}]  ({time.time()-t0:.0f}s)")
    return ok_all


# ------------------------------------------------------------------- L3
def leg3_contour(dps=50):
    ref = os.path.join(_ART, "contour_controls.json")
    if not _have_reference("L3_pair_contour", ref):
        return True
    t0 = time.time()
    fx = json.load(open(ref))
    bar = fx["bar_digits"]
    cur = QuarticCurve([_frac(c) for c in fx["quartic_coeffs_desc"]])
    with mp.workdps(dps + 30):
        rts = cur.roots(dps)
        pair = [i for i, rr in enumerate(rts) if mp.re(rr) < 0]
        if len(pair) != 2:
            raise AssertionError("conjugate pair not found")
        vals, meta = cur.cycle_pair_contour(tuple(pair), (0, 1, 2), dps)
        rows = {}
        ok_all = True
        sign = None
        for k, nm in enumerate(("Om", "X1", "X2")):
            o = fx["oracle"][nm]
            bv = mp.mpc(mp.mpf(o["re"]), mp.mpf(o["im"]))
            mine = vals[k]
            if sign is None:
                sign = 1 if abs(mine - bv) < abs(mine + bv) else -1
            d = shared_digits(sign * mine, bv)
            rows[nm] = round(d, 1)
            ok_all &= d >= bar
        rows["orientation_sign_vs_oracle"] = sign
        rows["contour_meta"] = {k: str(v) for k, v in meta.items()}
    REC["legs"]["L3_pair_contour"] = dict(dps=dps, rows=rows, bar=bar,
                                          pass_=bool(ok_all),
                                          secs=round(time.time() - t0, 1))
    print(f"L3 pair contour Om/X1/X2: "
          f"{ {k: rows[k] for k in ('Om', 'X1', 'X2')} } sign={sign} "
          f"(bar {bar}) [{'PASS' if ok_all else 'FAIL'}]  "
          f"({time.time()-t0:.0f}s)")
    return ok_all


# ------------------------------------------------------------------- L4
def leg4_sunrise(dps=100):
    t0 = time.time()
    rows = {}
    ok_all = True
    # rows data: bundled beside the package (ELLIPTICUS_ROWS_DIR overrides)
    from .qengine import _ROWS_DIR as RSD
    for row, fn in ((9, "row09_data.json"), (10, "row10_data.json"),
                    (11, "row11_data.json")):
        data = json.load(open(os.path.join(RSD, fn)))
        ms = [Fraction(m) for m in data["masses_sq"]]
        worst = 999
        for gp in data["gate_points"]:
            tq = Fraction(int(gp["t_num"]), int(gp["t_den"]))
            J, bound, diag = qengine.sunrise_J(tq, ms, dps)
            with mp.workdps(dps + 10):
                d = shared_digits(J, mp.mpf(gp["oracle_mid"]))
            worst = min(worst, d)
        rows[f"row{row:02d}"] = round(worst, 1)
        ok_all &= worst > min(95, dps - 5)
    em, eref, agree = qengine.sunrise_equal_mass_control(Fraction(-3),
                                                         dps=dps)
    rows["equal_mass_control"] = int(agree)
    ok_all &= agree >= dps - 8
    REC["legs"]["L4_sunrise_q"] = dict(dps=dps, rows=rows,
                                       pass_=bool(ok_all),
                                       secs=round(time.time() - t0, 1))
    print(f"L4 sunrise q-engine: {rows} [{'PASS' if ok_all else 'FAIL'}]  "
          f"({time.time()-t0:.0f}s)")
    return ok_all


# ------------------------------------------------------------------- L5
def leg5_ellred(dps=45):
    """ellred member: reproduce the reduction-gate class through the code
    paths at this home (ellipticus.ellred), against direct quadrature and the
    pinned reference strings in fixtures/battery/ellred_controls.json.
    a) quartic curve E (P4_at, cells A/B): chart modulus vs the pinned printed
       m string (bar = width-3) + closed F/E/Pi atoms vs direct quadrature
       (bar dps-12).
    b) cubic 3-real-root (E_alg + Kasin-image classes, cells A/B):
       Jmoments_cubic J0..J4 vs direct quadrature (bar min(30, dps-6): 30 at
       the default dps, scaled down for the low-dps --pilot pass) + pinned m.
    c) 1r2c E2 103849 cubic: analytic reduce_1r2c J0..J2 vs quadrature
       (bar 28) + Carlson R_F route agreement on J0 (bar 30) + pinned m
       string (bar 30) + J1 derivative certification."""
    ref = os.path.join(_ART, "ellred_controls.json")
    if not _have_reference("L5_ellred", ref):
        return True
    t0 = time.time()
    from . import ellred as EL
    ENG, CUB = EL.ellred_engine, EL.ellred_cubic
    R2C, CAR, EK = EL.ellred_1r2c, EL.ellred_carlson, EL.ellred_kernels
    ctrl = json.load(open(ref))
    rows, ok_all = [], True

    def _agree_pinned(val, s):
        width = len(s.replace("-", "").replace(".", "").lstrip("0"))
        return shared_digits(val, mp.mpf(s)), width

    with mp.workdps(dps + 15):
        # ---- a) quartic atoms, cells A/B ----
        for cell, tv, zv, x1f, x2f in [("A", "3/8", "3/10", "0.46", "0.54"),
                                       ("B", "1/4", "1/5", "0.66", "0.74")]:
            P4n = EK.P4_at(tv, zv)
            ch = ENG.chart(P4n)
            dm, wm = _agree_pinned(ch["m"], ctrl["quartic_E"][cell]["m"])
            x1, x2 = mp.mpf(x1f) ** 2, mp.mpf(x2f) ** 2
            z1, z2 = ENG.zof(ch, x1), ENG.zof(ch, x2)
            zlo, zhi = min(z1, z2), max(z1, z2)
            ch["_zwin"] = (zlo, zhi)
            ENG._calibrate_E(ch, zwin=(zlo, zhi))
            gz = lambda z: ENG.gz(ch, z)
            Fc = ch["jac"] / mp.sqrt(ch["A"]) * (ENG.atomF(ch, z2) - ENG.atomF(ch, z1))
            Fr = mp.quad(lambda t: 1 / mp.sqrt(EL.ell_normal.G_pval(P4n, t)), [x1, x2])
            if abs(Fc + Fr) < abs(Fc - Fr):
                Fc = -Fc
            Ec = ENG.atomE_closed(ch, zhi) - ENG.atomE_closed(ch, zlo)
            Er = mp.quad(lambda z: (ch["a2"] * z * z + ch["b2"]) / gz(z), [zlo, zhi])
            c = mp.mpf(2)
            Pc = ENG.atomPi_closed(ch, c, zhi) - ENG.atomPi_closed(ch, c, zlo)
            Pr = mp.quad(lambda z: 1 / ((z * z - c) * gz(z)), [zlo, zhi])
            dF, dE, dPi = (shared_digits(Fc, Fr), shared_digits(Ec, Er),
                           shared_digits(Pc, Pr))
            okq = dm >= wm - 3 and min(dF, dE, dPi) >= dps - 12
            ok_all &= okq
            rows.append(f"quartic {cell}: m {dm}/{wm}d F {dF:.0f} E {dE:.0f} "
                        f"Pi {dPi:.0f} [{'PASS' if okq else 'FAIL'}]")

        # ---- b) cubic 3-real-root moments ----
        for key, win in [("E_alg_A", ("0.46", "0.54")),
                         ("E_alg_B", ("0.66", "0.74")),
                         ("Kasin_img_A", ("0.46", "0.54")),
                         ("Kasin_img_B", ("0.66", "0.74"))]:
            b = ctrl["cubic"][key]
            cube = sp.Poly(sp.sympify(b["cube"]), sp.Symbol("x")).all_coeffs()
            C3n = [mp.mpf(str(c)) for c in reversed(cube)]
            while len(C3n) < 4:
                C3n.append(mp.mpf(0))
            x1, x2 = mp.mpf(win[0]) ** 2, mp.mpf(win[1]) ** 2
            J, meta = CUB.Jmoments_cubic(C3n, x1, x2, kmax=4)
            dm, wm = _agree_pinned(meta["m"], b["m"])
            # E_alg cubes are NEGATIVE on (e1,e2) in the pinned convention:
            # values are pure imaginary; compare complex-aware with ONE global
            # orientation sign fixed on J0 (orientation is convention, as L3).
            worst, sgn = 99, mp.mpf(1)
            for k, Jk in enumerate(J):
                ref_q = mp.quad(lambda t, k=k: t ** k /
                                mp.sqrt(mp.mpc(CUB.C3val(C3n, t))), [x1, x2])
                Jk = mp.mpc(Jk)
                if k == 0 and abs(Jk + ref_q) < abs(Jk - ref_q):
                    sgn = mp.mpf(-1)
                err = abs(sgn * Jk - ref_q) / abs(ref_q)
                worst = min(worst, float(-mp.log10(err)) if err > 0 else 99)
            cbar = min(30, dps - 6)   # 30 at the default dps; dps-relative
                                     # so the low-dps --pilot pass is judged
                                     # at its own working precision
            okc = dm >= wm - 3 and worst >= cbar
            ok_all &= okc
            rows.append(f"cubic {key} ({b['region']}): m {dm}/{wm}d "
                        f"J0..4 worst {worst:.0f} (bar {cbar}) "
                        f"[{'PASS' if okc else 'FAIL'}]")

        # ---- c) 1r2c analytic + Carlson route agreement (E2 103849) ----
        xx = sp.Symbol("x")
        cube = sp.expand(xx * (103849 * xx ** 2 + 25950 * xx + 5625))
        C3n = [mp.mpf(str(c)) for c in reversed(sp.Poly(cube, xx).all_coeffs())]
        while len(C3n) < 4:
            C3n.append(mp.mpf(0))
        x1, x2 = mp.mpf("0.46") ** 2, mp.mpf("0.54") ** 2
        res = R2C.reduce_1r2c(C3n, x1, x2, kmax=2)
        dm, wm = _agree_pinned(res["m"], ctrl["kasin_E2"]["modulus_m"])
        worst = 99
        for k, Jk in enumerate(res["J"]):
            ref_q = mp.quad(lambda t, k=k: t ** k /
                            mp.sqrt(CUB.C3val(C3n, t)), [x1, x2])
            worst = min(worst, shared_digits(mp.mpc(Jk).real, ref_q))
        resC = CAR.reduce_1r2c(C3n, x1, x2, kmax=1)
        dRF = shared_digits(mp.mpc(res["J"][0]).real,
                            mp.mpc(resC["J"][0]).real)
        cert1 = mp.mpf(res["J1_certify"])
        ok1 = (dm >= 30 and worst >= 28 and dRF >= 30 and
               cert1 <= mp.mpf(10) ** (-25))
        ok_all &= ok1
        rows.append(f"1r2c E2: m {dm}/{wm}d J0..2 worst {worst:.0f} "
                    f"Carlson-agree {dRF:.0f} J1_cert {mp.nstr(cert1, 3)} "
                    f"[{'PASS' if ok1 else 'FAIL'}]")

    REC["legs"]["L5_ellred"] = dict(dps=dps, rows=rows, pass_=bool(ok_all),
                                    secs=round(time.time() - t0, 1))
    print("L5 ellred:")
    for r in rows:
        print("    " + r)
    print(f"L5 {'PASS' if ok_all else 'FAIL'}  ({time.time()-t0:.0f}s)")
    return ok_all


def main():
    pilot = "--pilot" in sys.argv
    t0 = time.time()
    ok = True
    ok &= leg2_cycles(dps=45)
    ok &= leg3_contour(dps=50)
    ok &= leg4_sunrise(dps=100 if not pilot else 50)
    ok &= leg5_ellred(dps=45 if not pilot else 25)
    ok &= leg1_march(dps=30)
    skipped = [k for k, v in REC["legs"].items()
               if v.get("status") == "SKIP"]
    REC["all_pass"] = bool(ok)
    REC["skipped"] = skipped
    REC["secs_total"] = round(time.time() - t0, 1)
    out = env("ELLIPTICUS_BATTERY_OUT", "EMPL_BATTERY_OUT",
              os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "battery_results.json"))
    json.dump(REC, open(out, "w"), indent=1)
    print(f"BATTERY {'PASS' if ok else 'FAIL'}"
          + (f"  (skipped: {', '.join(skipped)})" if skipped else "")
          + f"  total {REC['secs_total']}s -> {out}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

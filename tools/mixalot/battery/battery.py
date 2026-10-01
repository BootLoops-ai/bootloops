"""MIXALOT battery — exercises the VENDORED engines on synthetic and
closed-form inputs, with mutation controls that must fail. Every check prints
PASS/FAIL; any FAIL => exit 1. Record -> battery/BATTERY.txt.

  M1  pins integrity (42 files).
  M2  MUTATION CONTROL: tampered vendor copy (1 byte) must refuse to compute.
  M3  cross-engine agreement grid (4 independent engines, exact equality).
  M4  LSX flagships: coin10 exact pair; swiss exact pair + 3,892,097-term
      support gate [--full]; LSX 5.5 3x3 N=132 held-out-prime residue vs the
      recorded 143/262-digit fraction [--full].
  M5  segmentation engine on a SYNTHETIC planted corpus (built in this
      file, no data dependence): the planted 2-segment structure is
      recovered (g*=2), and the exact-rational DP route cross-checks the
      float DP route to <= 5e-12 on every g.
  M6  blend evidence: both recorded GATE_FLIP_FROZEN_COMP cases, v1 == blind.
  M7  DP limit: Z_DPM((2,1),1)=5/48; g*(Z_g-Z_DPM)=-1/48 at g=64/128/256;
      three larger recorded rationals.
  M8  scaling row: bigk2 k=100 N=10000 lnZ=-41336.613834 (6dp string).
  M9  planted-seam power gate: mic_power arms vs recorded (gen 0/30 detected,
      half 0/30, language-swap 10/10) [--full; ~minutes].
  M10 comparator mutation controls: every exact comparator must FAIL on a
      perturbed value (guards against vacuous comparisons).
  M11-M19 member legs: MC-suite calibration (M11), input hardening (M12),
      worked examples (M13), and each member's own selftest/gate run from a
      clean cwd -- etienne_evaluate (M14), telegraph_evaluate (M15), the
      boxwalk subpackage via `python3 -m mixalot.boxwalk selftest` (M16,
      14 legs, its record routed to scratch), swap_route (M17), the
      Dirichlet emitter (M18), the eco/pilot gate_engines gate (M19).
The FAST battery is self-contained (synthetic + closed-form inputs built in
this file). The [--full] replay legs that compare against archived run
receipts SKIP by name unless their env-pointed read-only roots are set.
"""
import json
import os
import shutil
import subprocess
import sys
import time
from fractions import Fraction as F
from math import comb, log

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import mixalot  # noqa: E402
from mixalot.engines import load  # noqa: E402

FULL = "--full" in sys.argv
# Reference-data roots (env-pointed, read-only). Legs that replay recorded
# reference receipts SKIP by name when the roots are not configured.
PILOT = os.environ.get("MIXALOT_REFDATA_PILOT", "")   # pilot-engine result tree
SWEEP = os.environ.get("MIXALOT_REFDATA_SWEEP", "")   # census sweep tree (data/counts_*.json)

lines, fails = [], []


def check(name, ok, detail=""):
    s = "PASS" if ok else "FAIL"
    lines.append(f"{s}  {name}  {detail}")
    if not ok:
        fails.append(name)
    print(lines[-1], flush=True)


def info(msg):
    lines.append(f"INFO  {msg}")
    print(lines[-1], flush=True)


def skip(name, why):
    lines.append(f"SKIP  {name}  ({why})")
    print(lines[-1], flush=True)


def main():
    t_all = time.time()
    # ---------------- M1 pins
    rep = mixalot.verify(quiet=True)
    check("M1 pins integrity", len(rep["vendor_ok"]) == 42,
          f"{len(rep['vendor_ok'])} pinned; drift={len(rep['source_drift'])}")

    # ---------------- M2 vendor-tamper mutation control
    scratch = os.path.join(HERE, "_tamper_scratch")
    shutil.rmtree(scratch, ignore_errors=True)
    shutil.copytree(os.path.join(ROOT, "mixalot"),
                    os.path.join(scratch, "mixalot"))
    tgt = os.path.join(scratch, "mixalot", "vendor", "jeff", "bigg.py")
    b = open(tgt, "rb").read()
    open(tgt, "wb").write(b[:-1] + bytes([b[-1] ^ 1]))  # flip one bit, surely
    p = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); import mixalot; mixalot.verify()"
         % scratch], capture_output=True, text=True)
    check("M2 tampered vendor refuses (mutation control)",
          p.returncode != 0 and "VendorTamperError" in (p.stderr or ""),
          f"rc={p.returncode}")
    shutil.rmtree(scratch, ignore_errors=True)

    # ---------------- M2b pyc-shadow control: a poisoned same-size .pyc
    # planted in vendor __pycache__ must have NO effect — the loader
    # compiles from pinned source, never the cache.
    scratch = os.path.join(HERE, "_pyc_scratch")
    shutil.rmtree(scratch, ignore_errors=True)
    shutil.copytree(os.path.join(ROOT, "mixalot"),
                    os.path.join(scratch, "mixalot"))
    import py_compile
    real = os.path.join(scratch, "mixalot", "vendor", "pilot",
                        "closed_form_1var.py")
    poison_src = os.path.join(scratch, "_poison.py")
    body = open(real).read().replace("def Z_closed(", "def _Z_closed_real(")
    body += "\n\ndef Z_closed(u0, u1):\n    return 0\n"
    open(poison_src, "w").write(body)
    pyc_dir = os.path.join(scratch, "mixalot", "vendor", "pilot",
                           "__pycache__")
    os.makedirs(pyc_dir, exist_ok=True)
    tag = sys.implementation.cache_tag
    pyc_path = os.path.join(pyc_dir, f"closed_form_1var.{tag}.pyc")
    py_compile.compile(poison_src, cfile=pyc_path, doraise=True)
    # forge the header's (mtime, size) to match the REAL source so the
    # stdlib loader would accept the poisoned cache as fresh
    st = os.stat(real)
    import struct as _st
    pb = bytearray(open(pyc_path, "rb").read())
    pb[8:12] = _st.pack("<I", int(st.st_mtime) & 0xFFFFFFFF)
    pb[12:16] = _st.pack("<I", st.st_size & 0xFFFFFFFF)
    open(pyc_path, "wb").write(pb)
    p = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, %r); from mixalot.engines import "
         "load; v = load('closed_form_1var').Z_closed(10, 7); "
         "print(v); assert v != 0, 'POISONED PYC EXECUTED'"
         % scratch], capture_output=True, text=True)
    check("M2b pyc-shadow control: poisoned cache has NO effect",
          p.returncode == 0 and "8025569/2144534071680" in p.stdout,
          f"rc={p.returncode} out={p.stdout.strip()[:40]}")
    shutil.rmtree(scratch, ignore_errors=True)

    # ---------------- M3 cross-engine grid
    zs, cf = load("zseries"), load("closed_form_1var")
    w1, bg = load("w1_collapsed"), load("bigg")
    grid = [[10, 7], [5, 5], [12, 3], [7, 1]]
    ok = all(zs.z_1var_exact(U) == cf.Z_closed(*U) == w1.Z_g2(U)
             == bg.Z_bigg(U, 2) for U in grid)
    ok &= all(bg.Z_bigg(U, g) == w1.Z_general(U, g)
              for U in ([2, 1], [1, 2, 1]) for g in (2, 3))
    check("M3 cross-engine exact agreement (zseries/closed_form/w1/bigg)",
          ok, f"{len(grid)}+4 cells")

    # ---------------- M4 LSX flagships
    lsx = load("lsx_direct")
    t0 = time.time()
    m10 = lsx.Model.coin(lsx.COIN10_U)
    z10 = lsx.Z_phi(m10)
    check("M4a coin10 flagship == COIN10_REF (vendored, in-process)",
          z10 == lsx.COIN10_REF, f"{time.time()-t0:.2f}s")
    if FULL:
        t0 = time.time()
        msw = lsx.Model.from_table(lsx.SWISS_U)
        Zsw, support = lsx.Z_phi(msw, return_support=True)
        ok = (support == lsx.SWISS_SUPPORT
              and Zsw == F(lsx.SWISS_NUM, lsx.SWISS_DEN))
        check("M4b swiss flagship exact pair + 3,892,097-term support "
              "(vendored, in-process)", ok,
              f"support={support} {time.time()-t0:.0f}s")
        # LSX 5.5 held-out-prime residue vs recorded fraction
        if PILOT:
            t0 = time.time()
            rf = PILOT + "/LSX55_EXACT_RESULT.json"
            if not os.path.exists(rf):   # older reference trees use the previous name
                rf = PILOT + "/SCHIZO_EXACT_RESULT.json"
            sr = json.load(open(rf))
            num, den = int(sr["num"]), int(sr["den"])
            ok = (len(str(num)), len(str(den))) == (143, 262)
            check("M4c recorded LSX 5.5 fraction shape 143/262 digits", ok)
            pprime = 33554393
            want = num % pprime * pow(den % pprime, -1, pprime) % pprime
            zt = zs.Z_table_modp([[43, 16, 3], [6, 11, 10], [9, 18, 16]],
                                 pprime)
            check("M4d LSX 5.5 held-out-prime residue (engine vs recorded "
                  "fraction)", zt == want,
                  f"p={pprime} residue={zt} want={want} "
                  f"{time.time()-t0:.0f}s")
        else:
            skip("M4c/M4d recorded LSX 5.5 legs", "MIXALOT_REFDATA_PILOT unset")

    # ---------------- M5 segmentation engine on a SYNTHETIC planted corpus.
    # 12 units over a 6-symbol vocabulary with a planted seam after unit 6:
    # units 1-6 concentrate on symbols 1-3, units 7-12 on symbols 4-6 (with
    # small cross-leakage so nothing is degenerate). Deterministic — no data
    # files, nothing recorded; every expectation below is recomputed live.
    sv = load("seg_v1")
    a = F(1, 2)
    Xs = ([[7, 4, 2, 1, 0, 0], [6, 5, 1, 0, 1, 0], [8, 3, 2, 1, 0, 1],
           [7, 5, 1, 0, 1, 0], [6, 4, 3, 1, 0, 0], [8, 4, 1, 0, 1, 1]]
          + [[1, 0, 0, 7, 4, 2], [0, 1, 1, 6, 5, 1], [1, 0, 0, 8, 3, 3],
             [0, 1, 0, 7, 5, 2], [1, 0, 1, 6, 4, 3], [0, 0, 1, 8, 4, 2]])
    LWs = sv._logW_table(Xs, a)
    logZ = {}
    for g in range(1, 5):
        logZ[g] = sv.dp_float(Xs, g, a, LWs) - log(comb(len(Xs) - 1, g - 1))
    gstar = max(logZ, key=logZ.get)
    check("M5a planted 2-segment corpus: g* = 2 recovered", gstar == 2,
          " ".join(f"g{g}:{logZ[g]:.4f}" for g in sorted(logZ)))
    # exact-vs-float cross-route: the exact-rational DP recomputes every
    # cell's evidence and the float route must agree to <= 5e-12.
    mx_delta = 0.0
    for g in sorted(logZ):
        ex = sv.dp_exact(Xs, g, a)
        fl = sv.dp_float(Xs, g, a, LWs)
        import math as mm
        n_, d_ = ex.numerator, ex.denominator
        shn = max(0, n_.bit_length() - 500)
        shd = max(0, d_.bit_length() - 500)
        exf = (mm.log(n_ >> shn) + shn * mm.log(2)
               - mm.log(d_ >> shd) - shd * mm.log(2))
        mx_delta = max(mx_delta, abs(exf - fl))
    check("M5b exact-vs-float cross-route <= 5e-12 (synthetic, g=1..4)",
          mx_delta <= 5e-12, f"max delta {mx_delta:.3e}")

    # ---------------- M6 blend
    fc, fb = load("frozen_comp_v1"), load("frozen_comp_blind")
    p1 = [[F(7, 25), F(9, 25), F(9, 25)], [F(9, 17), F(5, 17), F(3, 17)]]
    b1, u1 = [F(9), F(3)], [1, 3, 0]
    zv = fc.evidence(p1, b1, u1)
    zb = fb.evidence_blind(p1, b1, u1)
    want1 = F(11721623169, 848260156250)
    check("M6a blend recorded case (v1==blind==recorded)",
          zv == zb == want1, str(zv))
    p2w = [[F(1, 3), F(7, 18), F(2, 9), F(1, 18)],
           [F(1, 9), F(1, 2), F(5, 18), F(1, 9)],
           [F(1, 2), F(1, 10), F(3, 10), F(1, 10)]]
    b2 = [F(5, 4), F(3), F(2)]
    u2 = [2, 0, 2, 2]
    bf_v = fc.presence_bayes_factor(p2w, b2, p2w[:2], b2[:2], u2)
    bf_b = fb.presence_bayes_factor_blind(p2w, b2, p2w[:2], b2[:2], u2)
    want2 = F(27542119402189, 8590418784375)
    check("M6b presence BF recorded case", bf_v == bf_b == want2, str(bf_v))

    # ---------------- M7 DP limit
    gf = load("w4_blind_gf")
    ok = gf.Zdpm_gf([2, 1], F(1)) == F(5, 48)
    check("M7a Z_DPM((2,1),1) = 5/48", ok)
    ok = all(g * (gf.Zg_gf([2, 1], F(1), g) - F(5, 48)) == F(-1, 48)
             for g in (64, 128, 256))
    check("M7b g*(Z_g - Z_DPM) = -1/48 exactly at g=64/128/256", ok)
    ok = (gf.Zdpm_gf([4, 3], F(1)) == F(29881, 5806080)
          and gf.Zdpm_gf([4, 3], F(22, 7)) == F(7238939, 1151479800)
          and gf.Zdpm_gf([3, 3, 2], F(1, 3))
          == F(12398482193, 223612316928000))
    check("M7c larger recorded DP rationals (3 instances)", ok)

    # ---------------- M8 scaling row
    bigk2 = load("bigk2")
    import random as _r
    _r.seed(42)
    k_, N_ = 100, 10000
    cuts = sorted(_r.sample(range(1, N_ + k_), k_ - 1))
    U = [b2_ - a2_ - 1 for a2_, b2_ in zip([0] + cuts, cuts + [N_ + k_])]
    t0 = time.time()
    Zv = bigk2.Z_fast2(U)
    dt = time.time() - t0

    def lnint(x):
        sh = max(0, x.bit_length() - 500)
        return log(x >> sh) + sh * 0.6931471805599453
    lnZ = lnint(Zv.numerator) - lnint(Zv.denominator)
    check("M8 scaling row k=100 N=10000 lnZ=-41336.613834",
          f"{lnZ:.6f}" == "-41336.613834", f"lnZ={lnZ:.6f} {dt:.2f}s")

    # ---------------- M9 planted-seam power (full)
    if FULL and not SWEEP:
        skip("M9 planted-seam power gate", "MIXALOT_REFDATA_SWEEP unset")
    if FULL and SWEEP:
        mp = load("mic_power")

        def chapter_matrix(book, channel):   # wrapper: the vendored
            # production_sweep.chapter_matrix body with the data path
            # parameterized to the env-pointed reference tree (read-only)
            X, vocab, chs = sv.load_matrix(
                SWEEP + f"/data/counts_{book}_A.json", channel)
            return X, vocab, [f"ch{c}" for c in chs]
        Xm, _, _ = chapter_matrix("Mic", "closed")
        sizes = [sum(u) for u in Xm]
        Xg, gv, _ = chapter_matrix("Gen", "closed")

        def pooled(Xx, rows):
            k = len(Xx[0])
            out = [0] * k
            for r in rows:
                for i, c in enumerate(Xx[r]):
                    out[i] += c
            return out
        pA = pooled(Xg, range(0, 11))
        pB = pooled(Xg, range(11, 50))
        wA = [x / sum(pA) for x in pA]
        wB = [x / sum(pB) for x in pB]
        t0 = time.time()
        gen = mp.arm("gen", sizes, wA, wB, len(gv), 30, 1000)
        wB2 = [(x + y) / 2 for x, y in zip(wA, wB)]
        half = mp.arm("half", sizes, wA, wB2, len(gv), 30, 2000)
        Xd, dv, chd = chapter_matrix("Dan", "closed")
        heb = pooled(Xd, [chd.index(f"ch{c}") for c in (1, 8, 9, 10, 11, 12)])
        ara = pooled(Xd, [chd.index(f"ch{c}") for c in (3, 4, 5, 6, 7)])
        wH = [x / sum(heb) for x in heb]
        wR = [x / sum(ara) for x in ara]
        swap = mp.arm("swap", sizes, wH, wR, len(dv), 10, 3000)
        check("M9 planted-seam power arms == recorded (0/30, 0/30, 10/10)",
              (gen["detected"], half["detected"], swap["detected"])
              == (0, 0, 10),
              f"gen={gen['detected']}/30 half={half['detected']}/30 "
              f"swap={swap['detected']}/10 {time.time()-t0:.0f}s")

    # ---------------- M11 exact-vs-MC gate cell (the comparison gate)
    estm = load("estimators")
    truth = lnfrac_ln = None
    zt = load("closed_form_1var").Z_closed(50, 50)
    import math as _mm
    n_, d_ = zt.numerator, zt.denominator
    shn = max(0, n_.bit_length() - 500)
    shd = max(0, d_.bit_length() - 500)
    truth = (_mm.log(n_ >> shn) + shn * _mm.log(2)
             - (_mm.log(d_ >> shd) + shd * _mm.log(2)))
    check("M11a exact truth == recorded GATE_ESTIMATORS truth_logZ",
          abs(truth - (-71.07939015729893)) < 1e-10, f"{truth:.11f}")
    rmc = estm.run("m1", 2, [50, 50], "nested", seed=1)
    bit = rmc["logZ_hat"] == -71.02035658502652
    zsc = abs((rmc["logZ_hat"] - truth) / rmc["err_est"])
    check("M11b vendored MC nested seed=1: calibrated vs exact (|z|<4)",
          zsc < 4, f"logZ_hat={rmc['logZ_hat']:.11f} z={zsc:.2f} "
          f"recorded-bit-repro={bit}")

    # ---------------- M12 input hardening (integer counts / positive priors /
    # measure name / size guards / g=1 shortcut / shape guard)
    import mixalot.core as mc
    hard = []
    try:
        mixalot.Z([1.5, 2.5], 1)
        hard.append("noninteger-counts-no-raise")
    except mc.NonIntegerCountError:
        pass
    try:
        mixalot.gstar([4, 1, 3, 2], gmax=2, prior={1: -1, 2: 2})
        hard.append("negative-prior-no-raise")
    except ValueError:
        pass
    try:
        mixalot.Z([10, 10], 2, measure="banana")
        hard.append("measure-no-raise")
    except ValueError:
        pass
    try:
        mixalot.Z([3_000_000, 5], 2)
        hard.append("size-guard-no-raise")
    except ValueError:
        pass
    t0 = time.time()
    ok1 = mixalot.Z([10 ** 12], 1) == 1 and time.time() - t0 < 1.0
    if not ok1:
        hard.append("g1-shortcut")
    try:
        mc.p_g1([[1, 2], [3, 4]])
        hard.append("shape-guard-no-raise")
    except ValueError:
        pass
    check("M12 input hardening (counts/priors/measure/size-guard/"
          "g1-shortcut/shape-guard)",
          not hard, f"failed={hard}")

    import tempfile
    # ---------------- M14 etienne-exact member:
    # vendored Etienne (2005) evaluator selftest (its own worked-example gate,
    # certified upstream vs the independent Hoppe/CRP urn DP + Etienne SA3)
    import subprocess as _sp
    _ep = os.path.join(os.path.dirname(os.path.abspath(mixalot.__file__)),
                       "vendor", "eco", "etienne_evaluate.py")
    _r = _sp.run([sys.executable, _ep, "--selftest"], capture_output=True,
                 text=True, cwd=tempfile.mkdtemp(prefix="m14_"), timeout=300)
    check("M14 etienne_evaluate --selftest rc=0 (vendored member)",
          _r.returncode == 0, (_r.stdout + _r.stderr).strip().splitlines()[-1]
          if (_r.stdout or _r.stderr) else "no output")

    # ---------------- M15 telegraph member: vendored
    # scrna evaluator selftest (its own certified reference gate)
    _tp = os.path.join(os.path.dirname(os.path.abspath(mixalot.__file__)),
                       "vendor", "scrna", "telegraph_evaluate.py")
    _tr = _sp.run([sys.executable, _tp, "--selftest"], capture_output=True,
                  text=True, cwd=tempfile.mkdtemp(prefix="m15_"), timeout=540)
    check("M15 telegraph_evaluate --selftest rc=0 (vendored member)",
          _tr.returncode == 0, (_tr.stdout + _tr.stderr).strip().splitlines()[-1]
          if (_tr.stdout or _tr.stderr) else "no output")

    # ---------------- M16 boxwalk member (mixalot/boxwalk/ subpackage): the
    # exact box-moment contiguity route's own 14-leg selftest, run through
    # the documented CLI `python3 -m mixalot.boxwalk selftest` from a clean
    # cwd with tools/mixalot on PYTHONPATH; its SELFTEST.json record is
    # routed into this leg's scratch dir, never the tracked tree.
    _m16_tmp = tempfile.mkdtemp(prefix="m16_")
    _m16_env = dict(os.environ,
                    # relative => lands in the scratch cwd; keeps the
                    # recorded TOTAL line free of machine paths
                    BOXWALK_SELFTEST_OUT="SELFTEST.json",
                    PYTHONPATH=ROOT + (os.pathsep + os.environ["PYTHONPATH"]
                                       if os.environ.get("PYTHONPATH")
                                       else ""))
    _br = _sp.run([sys.executable, "-m", "mixalot.boxwalk", "selftest"],
                  capture_output=True, text=True, cwd=_m16_tmp, timeout=900,
                  env=_m16_env)
    check("M16 boxwalk selftest rc=0 (member, python -m mixalot.boxwalk)",
          _br.returncode == 0 and "TOTAL: PASS" in _br.stdout,
          (_br.stdout.strip().splitlines()[-1] if _br.stdout else
           ((_br.stderr.strip().splitlines() or ["no output"])[-1])))

    # ---------------- M17 swap_route member: the k=2 float route's own
    # tiny-N float-vs-exact smoke test against the exact gf twins. The
    # verdict LINE is checked, not just rc (the script prints its verdict
    # either way and does not set an exit status).
    try:
        import numpy  # noqa: F401
        _np_ok = True
    except ImportError:
        _np_ok = False
    if _np_ok:
        _sw = os.path.join(os.path.dirname(os.path.abspath(mixalot.__file__)),
                           "vendor", "jeff", "swap_route.py")
        _swr = _sp.run([sys.executable, _sw, "selftest"], capture_output=True,
                       text=True, cwd=tempfile.mkdtemp(prefix="m17_"),
                       timeout=300)
        check("M17 swap_route selftest verdict PASS (k=2 float vs exact)",
              _swr.returncode == 0 and "SELFTEST PASS" in _swr.stdout,
              (_swr.stdout.strip().splitlines()[-1] if _swr.stdout
               else "no output"))
    else:
        skip("M17 swap_route selftest", "numpy absent (float route needs it)")

    # ---------------- M18 dirichlet emitter: the closed-form emitter's own
    # exact gate vs the brute-force Dirichlet engine, incl. the polynomial
    # branch (b0>A / b1>B numerator factors of G) and the uniform reduction
    _fe = os.path.join(os.path.dirname(os.path.abspath(mixalot.__file__)),
                       "vendor", "pilot", "formula_emitter_dirichlet.py")
    _fer = _sp.run([sys.executable, _fe], capture_output=True, text=True,
                   cwd=tempfile.mkdtemp(prefix="m18_"), timeout=300)
    check("M18 dirichlet emitter gate rc=0 incl. polynomial branch "
          "(vendored member)",
          _fer.returncode == 0 and "DIRICHLET EMITTER: PASS" in _fer.stdout,
          (_fer.stdout.strip().splitlines()[-1] if _fer.stdout
           else "no output"))

    # ---------------- M19 eco/pilot member: gate_engines EXECUTES its G1/G2
    # gate at load (exact product-tree == independent urn oracle; ball
    # containment/width; print-only). A gate failure exits via sys.exit(1)
    # inside the exec'd source, so SystemExit must be caught as a FAIL.
    t0 = time.time()
    try:
        load("gate_engines")
        check("M19 eco/pilot member gate (gate_engines G1/G2 at load)",
              True, f"{time.time()-t0:.2f}s")
    except SystemExit as e:
        check("M19 eco/pilot member gate (gate_engines G1/G2 at load)",
              False, f"SystemExit({e.code}) {time.time()-t0:.2f}s")

    # ---------------- M13 worked examples run as shipped, clean cwd,
    # print-only; the needles are live computations on the vendored engines
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = subprocess.run([sys.executable,
                            os.path.join(ROOT, "examples",
                                         "worked_examples.py")],
                           capture_output=True, text=True, cwd=td)
        wrote = os.listdir(td)
        need = ["-41336.613834", "-1/48", "PASS=True", "-71.07939015730"]
        missing = [n for n in need if n not in p.stdout]
        check("M13 worked examples: clean-cwd run, print-only, live "
              "computations present",
              p.returncode == 0 and not wrote and not missing,
              f"rc={p.returncode} wrote={wrote} missing={missing}")

    # ---------------- M10 comparator mutation controls
    muts = [
        ("blend", zv + F(1, 10**30) == want1),
        ("dp", gf.Zdpm_gf([2, 1], F(1)) == F(5, 48) + F(1, 10**30)),
        ("scaling", f"{lnZ + 1e-4:.6f}" == "-41336.613834"),
        # a perturbed float route must FAIL the M5b cross-route gate
        ("seg-crossroute", (mx_delta + 1e-3) <= 5e-12),
    ]
    bad = [n for n, hit in muts if hit]
    check("M10 comparator mutation controls all FAIL as required",
          not bad, f"vacuous={bad}")

    mode = "FULL" if FULL else "FAST"
    verdict = (f"OVERALL PASS [{mode}] "
               f"({sum(1 for L in lines if L.startswith('PASS'))} checks, "
               f"{time.time()-t_all:.0f}s)"
               if not fails else f"OVERALL FAIL [{mode}]: {', '.join(fails)}")
    lines.append(verdict)
    print(verdict, flush=True)
    # MIXALOT_BATTERY_OUT routes the record directory off the tracked tree
    # (the registered battery entry sets it to a mktemp dir); unset, the
    # record is written beside this script.
    out_dir = os.environ.get("MIXALOT_BATTERY_OUT") or HERE
    with open(os.path.join(out_dir, "BATTERY.txt"), "w") as f:
        f.write("MIXALOT battery record\n" + "\n".join(lines) + "\n")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())

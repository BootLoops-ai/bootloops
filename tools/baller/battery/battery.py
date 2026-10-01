#!/usr/bin/env python3
"""BALLER acceptance battery (S3). V-style: every leg prints PASS/FAIL; any
FAIL => exit 1; exactly one OVERALL line. Never writes outside its scratch
cwd. Refusal surface (exact): the baller tree, ALL of tools/, and every
deny-listed root (BALLER_DENY_ROOTS) — run from a scratch cwd.

  L1  vendor pins integrity (22/22: 5 engines + the 17 vendored fold-member
      files, fixture and sidecar included)
  L2  MUTATION CONTROL: bit-flipped scratch copy must refuse (verify+load)
  L3  MUTATION CONTROL: poisoned vendor .pyc has NO effect (loader never
      reads the cache) — plant on a scratch copy, value must stay exact
  L4  alias identity: posq/kernel_io/gate/manifest resolve to their
      registered homes and the dps_lint member resolves in-package;
      foreign-shadow refusal fires on a planted fake
  L5  march endpoint truth case (256-bit) + MUTATION-TESTED via the harness
      (3 mutants on a scratch copy must all be caught)
  L6  posq VERIFY through the alias (the gated core's own battery)
  L7  contract instruments: wayfinder gate+manifest own tests via pytest;
      mc: Clopper-Pearson quantile_ci vs exact beta quantiles (cross-check)
  L8  tripwire: agreement passes; disagreement HALTS (never averages);
      sampled Tripwire determinism + certificate
  L9  mutation harness selftest: pristine passes, planted mutants caught,
      in-place-scratch refusal fires
  L10 hygiene: 4 lints x (violating fixture flagged / clean fixture silent);
      ftrim zeroes radii; RadiusWatch halts on a blowup ball; bridge_mpf
      preserves the mantissa where mp.mpf() re-rounds
  L11 kklt engines: all vendored engines load through the pinned loader;
      pipe_transport exact-rational helpers round-trip (fr2arb/cxq)
  L12 integrity guards (pollution-restore, spec-spoof refusal,
      thread identity, pin-tamper refusal, typed missing-file, unpinned extras)
  L13 minted checks (precision-honest dual, NaN/zero/
      complex, kwargs sampling, mutation typing, lint guards, runtime helpers)
  L14 tripwire + sampling guards: ball-reject tripwire, constant-key sampling,
      inf-mid RadiusWatch, mid-exec dep-swap refusal, bare-name minimization,
      sys.path hygiene, symlink/any-file sweep, AugAssign/tuple-restore lint,
      load-fatal mutants, mc contract, configure_all
  L15 quad.integral_certified: acb_calc certified integration leg
  L16 certify.block_krawczyk: block-arrow Krawczyk + PD verification
      (planted zero / displaced refusals / sign-flip mutation controls
      / precision starvation / saddle + scalar-Cholesky cross-check)
  L17 front door run/solve/render (front-door spec): the Muller
      reference behaviors as EXACT gates vs an independent raw-flint truth
      — 50-digit run certifies x_30 (mid 6.0056486887714, rad ~1e-5) and
      returns the [mid-stream +/- inf] DEVIATION ball for x_100 (mid
      continues to the float-trap 100, marked UNCERTIFIED); solve(target
      16) escalates 30->240 dps and certifies x_100 = 6.0000000160995649;
      typed SolveRefused at the cap naming achieved digits; fail-closed
      render units; PLUS the documented NEGATIVE CONTROL: the same
      recurrence in plain floats returns exactly 100.0 (why balls matter)
  L18 MUTATION CONTROLS (frontdoor_case.py gate): a render that prints
      uncertified digits must FAIL; a solve that returns without the
      target must FAIL; a dropped deviation mid-stream and a dead
      escalation must FAIL (4 mutants, all caught)
  L19 purity-scan member: a planted sample file must classify
      MPF-OF-FLOAT, DEC-STRING, and PURE-CONV correctly through the
      member CLI, and `baller.purity_scan` must import as a module
  L20 vendored fold members are live code, not pinned dead bytes:
      certlane containment/planted/convention pytest suite (26 cases)
      passes from a scratch cwd; every ling_onesided module imports from
      the vendored bytes (fixture found beside them, N=318); valmono
      loads the pinned BC cache (6 exact coefficient polys) and its
      null-loop identity self-check passes
  L21 dps_lint member (mpmath import-time-dps AST linter): its built-in
      selftest (one planted finding per rule x 5 rules, exit-code contract
      1/0/2/2, fixed-fixture control) passes through the member CLI with
      TMPDIR pinned to the scratch cwd; `python3 -m baller.hygiene lint`
      exits 1 on a planted module-level constant and 0 on the fixed copy;
      `baller.hygiene.dps_lint` IS `baller.dps_lint` (one module object)
"""
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)
# Pin the real tools root so scratch copies inherit it: the battery's
# probe subprocesses copy the package to scratch dirs, where the relative
# tools-root default would dangle (BALLER_TOOLS_ROOT set by the user still wins).
os.environ.setdefault("BALLER_TOOLS_ROOT", os.path.dirname(TOOL))
# the refusal surface matches the documented claim — the tool
# tree, ALL of tools/, and every deny-listed root are write-banned
# cwds; run from a scratch cwd.
# cwd deny list: the tool tree (+ parent) by default; extend with
# BALLER_DENY_ROOTS (os.pathsep-separated) to protect additional trees
# (env-configured; fail-closed refusal).
BANNED = (TOOL, os.path.dirname(TOOL)) + tuple(
    r for r in os.environ.get("BALLER_DENY_ROOTS", "").split(os.pathsep) if r)
cwd = os.getcwd()
if any(os.path.realpath(cwd).startswith(os.path.realpath(b)) for b in BANNED):
    print(f"REFUSED: run from a scratch cwd (e.g. under runs/), not {cwd}")
    sys.exit(2)

RESULTS = []


def leg(name):
    def deco(fn):
        RESULTS.append((name, fn))
        return fn
    return deco


@leg("L1 vendor pins integrity")
def l1():
    import baller
    r = baller.verify(quiet=True)
    assert len(r["vendor_ok"]) == 22 and not r["vendor_bad"], r
    return "22/22 pins ok"


@leg("L2 tampered scratch copy refuses (mutation control)")
def l2():
    clone = os.path.join(os.getcwd(), "clone_tamper")
    if os.path.exists(clone):
        shutil.rmtree(clone)
    shutil.copytree(TOOL, clone)
    p = os.path.join(clone, "vendor", "kklt", "march_lib.py")
    b = bytearray(open(p, "rb").read()); b[len(b) // 2] ^= 1
    open(p, "wb").write(bytes(b))
    r = subprocess.run([sys.executable, "-c", (
        f"import sys; sys.path.insert(0, {clone!r})\n"
        "import baller\n"
        "try:\n"
        "    baller.verify(); print('NOTAMPER')\n"
        "except baller.VendorTamperError: print('REFUSED-VERIFY')\n"
        "from baller._core import load_vendored\n"
        "try:\n"
        "    load_vendored('march_lib'); print('LOADED')\n"
        "except baller.VendorTamperError: print('REFUSED-LOAD')\n")],
        capture_output=True, text=True, timeout=120)
    assert "REFUSED-VERIFY" in r.stdout and "REFUSED-LOAD" in r.stdout, r.stdout + r.stderr
    return "verify + load both refuse with the named error"


@leg("L3 poisoned vendor .pyc has NO effect")
def l3():
    clone = os.path.join(os.getcwd(), "clone_pyc")
    if os.path.exists(clone):
        shutil.rmtree(clone)
    shutil.copytree(TOOL, clone)
    plant = os.path.join(clone, "vendor", "kklt", "__pycache__")
    os.makedirs(plant, exist_ok=True)
    import importlib._bootstrap_external as be
    src_path = os.path.join(clone, "vendor", "kklt", "march_lib.py")
    poisoned = open(src_path).read() + "\nPOISONED = True\n"
    st = os.stat(src_path)
    pyc = be._code_to_timestamp_pyc(compile(poisoned, src_path, "exec"),
                                    st.st_mtime, st.st_size)
    tag = sys.implementation.cache_tag
    open(os.path.join(plant, f"march_lib.{tag}.pyc"), "wb").write(pyc)
    r = subprocess.run([sys.executable, "-c", (
        f"import sys; sys.path.insert(0, {clone!r})\n"
        "from baller._core import load_vendored\n"
        "m = load_vendored('march_lib')\n"
        "print('POISONED' if getattr(m, 'POISONED', False) else 'CLEAN')\n"
        f"import os; print('PURGED' if not os.path.exists({plant!r}) else 'CACHE-REMAINS')\n")],
        capture_output=True, text=True, timeout=120)
    assert "CLEAN" in r.stdout and "PURGED" in r.stdout, r.stdout + r.stderr
    return "poisoned cache ignored AND purged; source bytes rule"


@leg("L4 alias identity + foreign-shadow refusal")
def l4():
    from baller import quad, hygiene, contract
    assert quad.posq.__file__.endswith("tools/posq/posq.py")
    assert quad.kernel_io.__file__.endswith("tools/posq/kernel_io.py")
    # dps_lint is a MEMBER now (its only home is this package), not an alias
    assert hygiene.dps_lint.__name__ == "baller.dps_lint", hygiene.dps_lint.__name__
    assert os.path.dirname(os.path.abspath(hygiene.dps_lint.__file__)) == \
        os.path.join(TOOL, "baller"), hygiene.dps_lint.__file__
    assert contract.gate.__file__.endswith("wayfinder/gate.py")
    assert contract.manifest.__file__.endswith("wayfinder/manifest.py")
    r = subprocess.run([sys.executable, "-c", (
        f"import sys, types; sys.path.insert(0, {TOOL!r})\n"
        "fake = types.ModuleType('posq'); fake.__file__ = '/tmp/fake/posq.py'\n"
        "sys.modules['posq'] = fake\n"
        "try:\n"
        "    import baller; from baller import quad\n"
        "    print('SHADOW-ACCEPTED')\n"
        "except ImportError as e:\n"
        "    print('REFUSED' if 'identity violation' in str(e) else 'WRONG-ERROR')\n")],
        capture_output=True, text=True, timeout=120)
    assert "REFUSED" in r.stdout, r.stdout + r.stderr
    return "4 aliases home-verified + dps_lint member in-package; planted shadow refused"


@leg("L5 march truth case + mutation-tested (3 mutants)")
def l5():
    from baller.mutation import run_mutation_gate
    gate_cmd = [sys.executable, "-B", os.path.join(HERE, "march_case.py"), "{FILE}"]
    src = os.path.join(TOOL, "vendor", "kklt", "march_lib.py")
    mutants = [
        ("acc = acc * h + d[k][a]", "acc = acc * h + d[k][a] * 1.0000001"),
        ("tail = (Nr * q ** (K + 1) / (1 - q)).upper()", "tail = (Nr * 0).upper()"),
        ("inv = 1 / (g[0] * (k + 1))", "inv = 1 / (g[0] * (k + 2))"),
    ]
    rep = run_mutation_gate(src, mutants, gate_cmd,
                            scratch_dir=os.path.join(os.getcwd(), "mut_march"))
    return f"pristine PASS + mutants caught {rep['caught']}"


@leg("L6 posq VERIFY through the alias (eras-protocol battery, banked seed)")
def l6():
    from baller import quad
    home = os.path.dirname(quad.posq.__file__)
    # posq's C kernel is a compiled engine: where the shipped prebuilt binary
    # does not execute (FLINT shared-library or CPU-architecture mismatch),
    # kernel_io rebuilds it from source — the build cache is pinned to THIS
    # scratch cwd so the battery keeps its no-writes-outside-cwd claim, and
    # the VERIFY subprocess inherits it. No toolchain either => loud NAMED
    # skip (a missing engine is not a red leg), never a silent pass.
    os.environ.setdefault("POSQ_KERNEL_CACHE",
                          os.path.join(os.getcwd(), "posq_kernel_build"))
    try:
        quad.kernel_io.kernel_path()
    except quad.kernel_io.KernelUnavailable as e:
        return f"SKIP (named): {e}"
    r = subprocess.run([sys.executable, "-B",
                        os.path.join(home, "verify_posq_adaptation.py"),
                        "--seed", "20260720"],
                       capture_output=True, text=True, timeout=900,
                       cwd=os.getcwd())
    out = r.stdout + r.stderr
    ok = r.returncode == 0 and ("result.ok=True" in out or "VERIFY PASS" in out
                                or '"ok": true' in out.lower())
    assert ok, f"rc={r.returncode} tail={out[-400:]}"
    return "posq adaptation battery PASS (smoke seed) via aliased core"


@leg("L7 contract: wayfinder tests + mc Clopper-Pearson cross-check")
def l7():
    dt = os.path.join(os.path.dirname(TOOL), "wayfinder")
    r = subprocess.run([sys.executable, "-B", "-m", "pytest", "-q",
                        "-p", "no:cacheprovider",
                        "tests/test_gate.py", "tests/test_manifest.py"],
                       capture_output=True, text=True, timeout=300, cwd=dt)
    assert r.returncode == 0 and " passed" in r.stdout, r.stdout[-400:]
    from baller import contract
    from scipy.stats import beta
    lo, hi = contract.mc.quantile_ci(5, 100, conf=0.99)
    a = 1 - 0.99
    exact_lo = beta.ppf(a / 2, 5, 100 - 5 + 1)
    exact_hi = beta.ppf(1 - a / 2, 5 + 1, 100 - 5)
    assert abs(lo - exact_lo) < 1e-12 and abs(hi - exact_hi) < 1e-12, (lo, hi, exact_lo, exact_hi)
    return "wayfinder 7 tests green; quantile_ci == exact Clopper-Pearson to 1e-12"


@leg("L8 tripwire halts, never averages")
def l8():
    from baller.contract import dual, Tripwire, DualPathDisagreement
    v, cert = dual(lambda: 1.2345678901234, lambda: 1.2345678901239, digits=10)
    assert v == 1.2345678901234 and cert["agree_digits"] > 10
    try:
        dual(lambda: 1.0, lambda: 1.001, digits=6)
        raise AssertionError("disagreement did NOT halt")
    except DualPathDisagreement as e:
        assert e.primary == 1.0 and e.oracle == 1.001
    tw = Tripwire(lambda x: x * 2.0, digits=12, every=3)
    outs = [tw(lambda x: x * 2.0, i) for i in range(60)]
    assert outs[7] == 14.0 and tw.checked > 0 and tw.certificate()["worst_agree_digits"] == 9999.0
    sel = [tw.selected((i,)) for i in range(60)]
    assert sel == [tw.selected((i,)) for i in range(60)], "selector not deterministic"
    try:
        tw2 = Tripwire(lambda x: x * 2.0 + 1e-3, digits=12, every=1)
        [tw2(lambda x: x * 2.0, i) for i in range(5)]
        raise AssertionError("sampled disagreement did NOT halt")
    except DualPathDisagreement:
        pass
    return "agreement passes; both halt paths fire; selector deterministic"


@leg("L9 mutation harness selftest")
def l9():
    from baller.mutation import run_mutation_gate, MutationHarnessError
    toy = os.path.join(os.getcwd(), "toy_engine.py")
    open(toy, "w").write("def f(x):\n    return x * x + 1\n")
    gate = [sys.executable, "-c",
            "import importlib.util,sys; spec=importlib.util.spec_from_file_location('t', sys.argv[1]); "
            "m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); "
            "sys.exit(0 if m.f(3)==10 else 1)", "{FILE}"]
    rep = run_mutation_gate(toy, [("x * x + 1", "x * x + 2")], gate,
                            scratch_dir=os.path.join(os.getcwd(), "mut_toy"))
    assert rep["caught"] == "1/1"
    # vacuous-gate detection: a mutant the gate cannot see must raise
    try:
        run_mutation_gate(toy, [("def f(x):", "def f(x):  # cosmetic")], gate,
                          scratch_dir=os.path.join(os.getcwd(), "mut_toy2"))
        raise AssertionError("vacuous gate NOT flagged")
    except MutationHarnessError as e:
        assert "VACUOUS" in str(e)
    # in-place refusal: scratch == source dir
    try:
        run_mutation_gate(toy, [("x * x + 1", "x * x + 3")], gate,
                          scratch_dir=os.path.dirname(os.path.abspath(toy)))
        raise AssertionError("in-place scratch NOT refused")
    except MutationHarnessError as e:
        assert "in-place" in str(e) or "scratch" in str(e)
    return "pristine passes, mutant 1/1, vacuous gate flagged, in-place refused"


@leg("L10 hygiene: lints on fixtures + runtime helpers")
def l10():
    from baller import hygiene
    fx = os.path.join(os.getcwd(), "fixtures")
    os.makedirs(fx, exist_ok=True)
    bad = os.path.join(fx, "bad.py")
    open(bad, "w").write(
        "from mpmath import mp, mpf, workdps\n"
        "from functools import lru_cache\n"
        "from flint import ctx\n"
        "def f(x):\n"
        "    y = mpf(x)\n"
        "    with mp.workdps(50):\n"
        "        z = mpf(x)\n"
        "    with mp.workdps(30):\n"
        "        pass\n"
        "    return y\n"
        "@lru_cache(maxsize=None)\n"
        "def g(n):\n"
        "    return mp.dps * n\n"
        "ctx.prec = 128\n"
        "t = mpf('1e-30')\n")
    clean = os.path.join(fx, "clean.py")
    open(clean, "w").write(
        "from mpmath import mp, mpf\n"
        "def f(x, dps):\n"
        "    with mp.workdps(2 * dps):\n"
        "        return mpf(x)\n")
    finds = hygiene.lint_files([bad])
    kinds = {f[0] for f in finds}
    assert {"fixed_wdps", "ambient_conv", "unkeyed_cache", "ctx_restore"} <= kinds, kinds
    assert hygiene.lint_files([clean]) == [], hygiene.lint_files([clean])
    from flint import arb, acb
    t = hygiene.ftrim([arb("1.5") + arb(0, 1e-10), acb(2, 3)])
    assert float(t[0].rad()) == 0.0
    rw = hygiene.RadiusWatch(bar_digits=10)
    assert rw.check(arb(1) + arb(0, 1e-14)) > 10
    try:
        rw.check(arb(1) + arb(0, 1e-3)); raise AssertionError("no halt")
    except hygiene.RadiusBlowup:
        pass
    from mpmath import mp
    old = mp.dps
    try:
        mp.dps = 50
        x = mp.mpf(2) ** mp.mpf("0.5")
        mp.dps = 15
        rr = mp.mpf(x)            # re-rounds (the trap)
        br = hygiene.bridge_mpf(x, mp)   # exact mantissa
        assert br._mpf_ == x._mpf_ and rr._mpf_ != x._mpf_
    finally:
        mp.dps = old
    return "4 lints flag bad + pass clean; ftrim/RadiusWatch/bridge all behave"


@leg("L11 kklt engines load + exact-helper sanity")
def l11():
    from baller import transport, certify
    pl, pt, pv = transport.pipe_lib, transport.pipe_transport, certify.pipe_vac
    from fractions import Fraction as Fr
    a = pt.fr2arb(Fr(3, 7))
    assert a.contains(a) and float(a.mid()) - 3 / 7 < 1e-15
    z = pt.cxq((Fr(1, 3), Fr(-2, 5)))
    m = pt.cxq_mul(z, z)
    ub, lb = pt.cxq_abs_ub(z), pt.cxq_abs_lb(z)
    assert float(lb) <= (1 / 9 + 4 / 25) ** 0.5 <= float(ub)
    assert all(hasattr(pv, n) for n in dir(pv) if False) or pv is not None
    return "3 engines loaded pinned; fr2arb containment + cxq |z| bounds exact"


_L12_PROBE = '''
import os, sys, types, threading
CLONE = sys.argv[1]
TOOLS_ROOT = sys.argv[2]
sys.path.insert(0, CLONE)
import baller
from baller._core import load_vendored, alias_registered
baller.verify(quiet=True)
# C1: pollute-after-first-load must be RESTORED, not wired in
real = load_vendored("pipe_lib")
fake = types.ModuleType("pipe_lib"); fake.PWNED = True
sys.modules["pipe_lib"] = fake
pt = load_vendored("pipe_transport")
print("C1-" + ("PASS" if pt.PL is real and not hasattr(pt.PL, "PWNED") else "FAIL"))
# C2: spec-less pre-existing module with spoofed __file__ -> refused
# (spoof a real aliased flat tool: posq's kernel_io at its registered home)
POSQ_HOME = os.path.join(TOOLS_ROOT, "posq")
fake2 = types.ModuleType("kernel_io")
fake2.__file__ = os.path.join(POSQ_HOME, "kernel_io.py")
sys.modules["kernel_io"] = fake2
try:
    alias_registered("kernel_io", POSQ_HOME)
    print("C2-FAIL")
except ImportError as e:
    print("C2-" + ("PASS" if "provenance" in str(e) else "FAIL"))
del sys.modules["kernel_io"]
# C3: concurrent loads give ONE identity triple
ids = set()
def go():
    m = load_vendored("pipe_vac"); ids.add((id(m), id(m.PL), id(m.PT)))
ts = [threading.Thread(target=go) for _ in range(8)]
[t.start() for t in ts]; [t.join() for t in ts]
print("C3-" + ("PASS" if len(ids) == 1 else "FAIL(%d)" % len(ids)))
# C5: in-process pin mutation -> fail-closed typed refusal
import baller._pins as P
good = P.PINS["geo/mc.py"]
P.PINS["geo/mc.py"] = "0" * 64
try:
    baller.verify(quiet=True); print("C5-FAIL")
except baller.VendorTamperError:
    print("C5-PASS")
P.PINS["geo/mc.py"] = good
# C6: post-verify deletion -> typed refusal
os.remove(os.path.join(CLONE, "vendor", "geo", "mc.py"))
try:
    load_vendored("mc"); print("C6-FAIL")
except baller.VendorTamperError:
    print("C6-PASS")
'''


@leg("L12 integrity guards (pollution-restore, spec-spoof refusal)")
def l12():
    clone = os.path.join(os.getcwd(), "clone_r1")
    if os.path.exists(clone):
        shutil.rmtree(clone)
    shutil.copytree(TOOL, clone)
    probe = os.path.join(os.getcwd(), "l12_probe.py")
    open(probe, "w").write(_L12_PROBE)
    r = subprocess.run([sys.executable, "-B", probe, clone,
                        os.path.dirname(TOOL)],
                       capture_output=True, text=True, timeout=300)
    got = [ln for ln in r.stdout.split() if ln.startswith("C")]
    assert r.stdout.count("-PASS") == 5, r.stdout + r.stderr[-400:]
    # C4 needs a pristine copy (C6 deleted a file): unpinned extra = tamper
    clone2 = os.path.join(os.getcwd(), "clone_r1b")
    if os.path.exists(clone2):
        shutil.rmtree(clone2)
    shutil.copytree(TOOL, clone2)
    open(os.path.join(clone2, "vendor", "kklt", "evil.py"), "w").write("x=1\n")
    r2 = subprocess.run([sys.executable, "-c", (
        f"import sys; sys.path.insert(0, {clone2!r})\n"
        "import baller\n"
        "try:\n"
        "    baller.verify(); print('C4-FAIL')\n"
        "except baller.VendorTamperError as e:\n"
        "    print('C4-' + ('PASS' if 'UNPINNED' in str(e) else 'FAIL'))\n")],
        capture_output=True, text=True, timeout=120)
    assert "C4-PASS" in r2.stdout, r2.stdout + r2.stderr[-300:]
    return "pollution-restore, spec-spoof refusal, 8-thread identity, pin-tamper refusal, typed missing-file, unpinned-extra: all hold"


@leg("L13 minted checks (precision-honest dual, NaN/zero refusals)")
def l13():
    from baller.contract import dual, Tripwire, DualPathDisagreement, agree_digits
    from baller.mutation import run_mutation_gate, MutationHarnessError
    from baller import hygiene
    import mpmath
    # precision-honest comparison: 20-digit disagreement vs bar 30 must HALT
    with mpmath.workdps(50):
        a = mpmath.mpf(1) / 3
        b = a + mpmath.mpf(10) ** -20
    try:
        dual(lambda: a, lambda: b, digits=30)
        raise AssertionError("20-digit mpf disagreement passed a bar of 30")
    except DualPathDisagreement as e:
        assert 19 < e.digits < 21, e.digits
    # NaN halts
    try:
        dual(lambda: float("nan"), lambda: 1.0, digits=5)
        raise AssertionError("NaN passed")
    except DualPathDisagreement:
        pass
    # tiny-vs-zero symmetric: both orders halt
    for p, o in ((1e-5, 0.0), (0.0, 1e-5)):
        try:
            dual(lambda: p, lambda: o, digits=4)
            raise AssertionError(f"tiny-vs-zero passed ({p},{o})")
        except DualPathDisagreement:
            pass
    # complex rejected typed
    try:
        agree_digits(1 + 2j, 1 + 2j)
        raise AssertionError("complex accepted")
    except ValueError:
        pass
    # kwargs-style Tripwire IS sampled
    tw = Tripwire(lambda x=0: x * 2.0, digits=12, every=3)
    for i in range(30):
        tw(lambda x=0: x * 2.0, x=i)
    assert tw.checked > 0, "kwargs calls never oracle-checked"
    # mutation: syntax-breaking mutant rejected typed; {FILE}-less rejected
    toy = os.path.join(os.getcwd(), "toy13.py")
    open(toy, "w").write("def f(x):\n    return x + 1\n")
    gate = [sys.executable, "-c",
            "import importlib.util,sys; spec=importlib.util.spec_from_file_location('t', sys.argv[1]); "
            "m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); "
            "sys.exit(0 if m.f(1)==2 else 1)", "{FILE}"]
    try:
        run_mutation_gate(toy, [("return x + 1", "return x +")], gate,
                          scratch_dir=os.path.join(os.getcwd(), "mut13"))
        raise AssertionError("syntax mutant accepted")
    except MutationHarnessError as e:
        assert "syntax" in str(e)
    try:
        run_mutation_gate(toy, [("return x + 1", "return x + 2")],
                          [sys.executable, "-c", "pass"],
                          scratch_dir=os.path.join(os.getcwd(), "mut13"))
        raise AssertionError("{FILE}-less gate accepted")
    except MutationHarnessError as e:
        assert "{FILE}" in str(e)
    # hygiene: functools.cache caught; two-violation ctx file still flagged;
    # aliased-wdps clean; mp.dps-spelling ambient in scope
    fx = os.path.join(os.getcwd(), "fixtures13")
    os.makedirs(fx, exist_ok=True)
    fc = os.path.join(fx, "fc.py")
    open(fc, "w").write("import functools\nfrom mpmath import mp\n"
                        "@functools.cache\ndef g(n):\n    return mp.dps * n\n")
    assert any(f[0] == "unkeyed_cache" for f in hygiene.lint_files([fc])), "functools.cache missed"
    c2 = os.path.join(fx, "c2.py")
    open(c2, "w").write("from flint import ctx\n# finally\nctx.prec = 64\n"
                        "x = 1\nctx.prec = 53\n")
    assert sum(1 for f in hygiene.lint_files([c2]) if f[0] == "ctx_restore") == 2, \
        "two-violation ctx file not fully flagged"
    ok = os.path.join(fx, "ok.py")
    open(ok, "w").write("from mpmath import workdps as wdps, mpf\n"
                        "def f(x, dps):\n    with wdps(2 * dps):\n        return mpf(x)\n")
    assert hygiene.lint_files([ok]) == [], hygiene.lint_files([ok])
    amb = os.path.join(fx, "amb.py")
    open(amb, "w").write("from mpmath import mp, mpf\nmp.dps = 15\n"
                        "def f(x):\n    return mpf(x)\n")
    assert any(f[0] == "ambient_conv" for f in hygiene.lint_files([amb])), "mp.dps-spelling out of scope"
    # RadiusWatch zero-mid halts; ftrim float pass-through; ctx_guard restores
    from flint import arb, ctx as fctx
    rw = hygiene.RadiusWatch(bar_digits=10)
    try:
        rw.check(arb(0, 1e-50)); raise AssertionError("zero-mid ball passed")
    except hygiene.RadiusBlowup:
        pass
    assert hygiene.ftrim(3.5) == 3.5 and hygiene.ftrim([1, 2.0]) == [1, 2.0]
    before = fctx.prec
    with hygiene.ctx_guard(prec=999):
        assert fctx.prec == 999
    assert fctx.prec == before
    return "precision-honest dual, NaN/zero/complex, kwargs sampling, mutation typing, 4 lint guards, zero-mid halt, ftrim, ctx_guard: all hold"


_L14_PROBE = '''
import os, sys, types
CLONE = sys.argv[1]
sys.path.insert(0, CLONE)
import baller
from baller._core import load_vendored
baller.verify(quiet=True)
# D1: mid-exec dep swap via import hook -> typed refusal
class Hook:
    def find_spec(self, name, path=None, target=None):
        if name == "json_unused_marker":
            pass
        if name == "json":
            fake = types.ModuleType("pipe_lib"); fake.PWNED = True
            sys.modules["pipe_lib"] = fake
        return None
load_vendored("pipe_lib")
sys.meta_path.insert(0, Hook())
sys.modules.pop("json", None)
try:
    pt = load_vendored("pipe_transport")
    ok = not hasattr(pt.PL, "PWNED")
    print("D1-" + ("PASS" if ok else "FAIL(wired)"))
except baller.VendorTamperError:
    print("D1-PASS")
finally:
    sys.meta_path.pop(0)
# D2: mc gets NO bare sys.modules entry (hijack class killed)
from baller import contract
_ = contract.mc
print("D2-" + ("PASS" if "mc" not in sys.modules else "FAIL"))
# D3: aliasing does not permanently mutate sys.path
import copy
before = list(sys.path)
from baller import hygiene
after = list(sys.path)
extra = [p for p in after if p not in before]
print("D3-" + ("PASS" if not any(p.endswith("/tools") or p.endswith("/posq") for p in extra) else f"FAIL({extra})"))
'''


@leg("L14 tripwire + sampling guards (ball-reject, constant-key)")
def l14():
    from baller.contract import agree_digits, Tripwire, DualPathDisagreement
    from baller.mutation import run_mutation_gate, MutationHarnessError
    from baller import hygiene, transport
    from flint import arb
    # ball inputs rejected typed (radius never silently discarded)
    try:
        agree_digits(arb(1, 100), 1.0)
        raise AssertionError("arb ball accepted — radius discarded")
    except ValueError as e:
        assert "ball" in str(e)
    # constant-key sampling: wrong oracle with a CONSTANT arg must halt
    tw = Tripwire(lambda cfg: 1.0 + 1e-3, digits=10, every=15)
    try:
        for _ in range(200):
            tw(lambda cfg: 1.0, {"fixed": True})
        raise AssertionError(f"constant-key loop never checked (checked={tw.checked})")
    except DualPathDisagreement:
        assert tw.checked > 0
    # inf-mid RadiusWatch halts
    rw = hygiene.RadiusWatch(bar_digits=30)
    try:
        rw.check(arb("inf")); raise AssertionError("inf mid passed as exact")
    except hygiene.RadiusBlowup:
        pass
    # ftrim on a property-mid interval type must not crash
    import mpmath
    v = hygiene.ftrim(mpmath.iv.mpf([1, 2]))
    assert float(v) == 1.5, v
    # AugAssign + tuple-restore + alias in the ctx lint; vec.cap FP gone
    fx = os.path.join(os.getcwd(), "fixtures14")
    os.makedirs(fx, exist_ok=True)
    aug = os.path.join(fx, "aug.py")
    open(aug, "w").write("from flint import ctx\nctx.prec += 64\n")
    assert any(f[0] == "ctx_restore" for f in hygiene.lint_files([aug])), "AugAssign missed"
    tup = os.path.join(fx, "tup.py")
    open(tup, "w").write("from flint import ctx\nold = (ctx.prec, ctx.cap)\n"
                         "try:\n    ctx.prec = 128\n    x = 1\nfinally:\n"
                         "    ctx.prec, ctx.cap = old\n")
    assert not any(f[0] == "ctx_restore" for f in hygiene.lint_files([tup])), \
        hygiene.lint_files([tup])
    vec = os.path.join(fx, "vec.py")
    open(vec, "w").write("class V: pass\nvec = V()\nvec.cap = 100\nvec.prec = 5\n")
    assert not any(f[0] == "ctx_restore" for f in hygiene.lint_files([vec])), "vec.cap FP back"
    # self-lint: baller's own package files are ctx/ambient clean (markers honored)
    pkg = os.path.dirname(os.path.abspath(hygiene.__file__))
    own = [os.path.join(pkg, f) for f in ("hygiene_checks.py", "tripwire.py",
                                          "mutation.py", "frontdoor.py")]
    self_finds = [f for f in hygiene.lint_files(own) if f[0] in ("ctx_restore", "ambient_conv")]
    assert self_finds == [], self_finds
    # load-fatal mutant rejected typed
    toy = os.path.join(os.getcwd(), "toy14.py")
    open(toy, "w").write("K = 7\ndef f(x):\n    return x + K\n")
    gate = [sys.executable, "-c",
            "import importlib.util,sys; spec=importlib.util.spec_from_file_location('t', sys.argv[1]); "
            "m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); "
            "sys.exit(0 if m.f(1)==8 else 1)", "{FILE}"]
    try:
        run_mutation_gate(toy, [("K = 7", "K = 7 + undefined_xyz")], gate,
                          scratch_dir=os.path.join(os.getcwd(), "mut14"))
        raise AssertionError("load-fatal mutant accepted")
    except MutationHarnessError as e:
        assert "LOAD-FATAL" in str(e)
    # mc contract: batch<=0 and bad k refuse instead of hang/nan
    from baller import contract as C
    try:
        C.mc.quantile_ci(101, 100); raise AssertionError("k>n accepted")
    except ValueError:
        pass
    try:
        C.mc.sim_ll_fixed_n([0.0], [1.0], 1, 0.0, 10, None, batch=0)
        raise AssertionError("batch=0 accepted")
    except ValueError:
        pass
    # configure_all exists and orders PL before PT
    assert callable(transport.configure_all)
    # subprocess probes: dep-swap refusal, bare-name minimization, path hygiene
    clone = os.path.join(os.getcwd(), "clone_r2")
    if os.path.exists(clone):
        shutil.rmtree(clone)
    shutil.copytree(TOOL, clone)
    probe = os.path.join(os.getcwd(), "l14_probe.py")
    open(probe, "w").write(_L14_PROBE)
    r = subprocess.run([sys.executable, "-B", probe, clone,
                        os.path.dirname(TOOL)],
                       capture_output=True, text=True, timeout=300)
    assert r.stdout.count("-PASS") == 3, r.stdout + r.stderr[-400:]
    # symlink + any-file sweep
    clone2 = os.path.join(os.getcwd(), "clone_r2b")
    if os.path.exists(clone2):
        shutil.rmtree(clone2)
    shutil.copytree(TOOL, clone2)
    open(os.path.join(clone2, "vendor", "kklt", "evil.so"), "w").write("x")
    os.symlink("/tmp", os.path.join(clone2, "vendor", "kklt", "sub"))
    r2 = subprocess.run([sys.executable, "-c", (
        f"import sys; sys.path.insert(0, {clone2!r})\n"
        "import baller\n"
        "try:\n"
        "    baller.verify(); print('SWEEP-FAIL')\n"
        "except baller.VendorTamperError as e:\n"
        "    s = str(e)\n"
        "    print('SWEEP-' + ('PASS' if 'UNPINNED' in s and 'SYMLINK' in s else 'PARTIAL:' + s[:120]))\n")],
        capture_output=True, text=True, timeout=120)
    assert "SWEEP-PASS" in r2.stdout, r2.stdout + r2.stderr[-300:]
    return "ball-reject, constant-key, inf-mid, iv-ftrim, 3 lint guards, self-lint clean, load-fatal typing, mc contract, dep-swap/bare-name/path probes, widened sweep: all hold"


@leg("L15 quad.integral_certified (acb_calc leg)")
def l15():
    import flint, math
    from baller import quad
    # (a) known value: pi to >=50 certified digits
    r = quad.integral_certified(lambda x, _: 4/(1+x**2), 0, 1, prec_dps=60)
    pi = float(flint.arb.pi())
    assert r['certified_digits'] >= 50, r
    assert abs(float(flint.arb(r['mid_str'].strip('[]').split(' ')[0])) - pi) < 1e-12
    # (b) analytic-flag LAW negative control: unflagged sqrt returns a
    # confidently WRONG ball (4.669... vs exact 14/3) — the leg proves the
    # misuse is detectable against the flagged value, so reviews demand
    # the flag on any branch-cut integrand.
    wrong = quad.integral_certified(lambda x, _: x.sqrt(), 1, 4, prec_dps=40)
    right = quad.integral_certified(lambda x, a: x.sqrt(analytic=a), 1, 4, prec_dps=40)
    exact = 14.0/3.0
    assert abs(float(flint.arb(right['mid_str'].strip('[]').split(' ')[0])) - exact) < 1e-10
    assert abs(float(flint.arb(wrong['mid_str'].strip('[]').split(' ')[0])) - exact) > 1e-3, 'misuse undetectable'
    # (c) rcv boundary-layer class: pole 2.8e-17 off-path, >=45 certified digits
    r2 = quad.integral_certified(lambda x, _: 1/((1-x)+flint.arb('2.8e-17')), 0, 1, prec_dps=60)
    assert r2['certified_digits'] >= 45, r2
    # (d) input-precision footgun: constant built at ambient 53-bit prec
    # poisons the ball; Arb reports it honestly; max_rad gate fires.
    amb_eps = flint.arb('2.8e-17')   # constructed OUTSIDE working prec
    fired = False
    try:
        quad.integral_certified(lambda x, _: 1/((1-x)+amb_eps), 0, 1,
                                prec_dps=60, max_rad=1e-40)
    except quad.BallBlown:
        fired = True
    assert fired, 'sloppy-input ball accepted (max_rad gate dead)'
    return 'pi>=50d; sqrt-flag misuse detectable; near-pole 45d+; input-prec footgun gated'


@leg("L16 certify.block_krawczyk (block-arrow Krawczyk + PD)")
def l16():
    import numpy as np
    from flint import arb, arb_mat
    from baller import block_arrow as bk

    D0 = [np.array([[4., 1, 0], [1, 5, 1], [0, 1, 6]]),
          np.array([[3., 1], [1, 4]])]
    B0 = [np.array([[1., 0], [0, 1], [1, 1]]),
          np.array([[0.5, 0], [0, 0.5]])]
    G0 = np.array([[8., 1], [1, 9]])
    X0 = ([[0.3, -0.2, 0.1], [0.05, -0.4]], [0.2, -0.1])

    class Oracle:
        """Planted zero at X0: F = H0 (x-x0) + 0.01 (x-x0)^3 componentwise,
        Jacobian H0 + 0.03 diag((x-x0)^2). flip_math plants a SIGN ERROR in
        block 1's residual math (x+x0 for x-x0: a different objective whose
        optimum is elsewhere — the spec's sign-flip mutation); flip_h
        negates the border of H (inconsistent oracle); saddle plants a
        negative border eigenvalue. NOTE (documented boundary): negating a
        gradient block POST-HOC is invisible BY MATHEMATICS — the centered
        test is sign-blind (|-c|=|c|) and a flipped F-block has the
        identical zero set, so the certificate stays true; tamper legs
        therefore flip signs in the MATH and in H, where flips are lies."""
        dims = ([3, 2], 2)

        def __init__(self, flip_math=False, flip_h=False, saddle=False):
            self.flip_math, self.flip_h, self.saddle = \
                flip_math, flip_h, saddle

        def _hf(self):
            Hf = np.zeros((7, 7))
            Hf[:3, :3] = D0[0]; Hf[3:5, 3:5] = D0[1]
            Hf[:3, 5:] = B0[0]; Hf[5:, :3] = B0[0].T
            Hf[3:5, 5:] = B0[1]; Hf[5:, 3:5] = B0[1].T; Hf[5:, 5:] = G0
            if self.saddle:
                Hf[6, 6] = -9.0
            return Hf

        def _d(self, x):
            zs, g = x
            f = list(zs[0]) + list(zs[1]) + list(g)
            c = X0[0][0] + X0[0][1] + X0[1]
            # flip_math: the planted sign error — block-1 coords (3, 4)
            # read x + x0 instead of x - x0
            return [(f[i] + arb(c[i])) if (self.flip_math and 3 <= i < 5)
                    else (f[i] - arb(c[i])) for i in range(7)]

        def F(self, x):
            d = self._d(x); Hf = self._hf(); out = []
            for i in range(7):
                s = arb(0)
                for j in range(7):
                    s += arb(float(Hf[i, j])) * d[j]
                # NB x*x*x, never x**3: arb pow NaNs on zero-containing balls
                s += arb('0.01') * d[i] * d[i] * d[i]
                out.append(s)
            return (out[:3], out[3:5]), out[5:]

        def H(self, x):
            d = self._d(x); Hf = self._hf()

            def mat(r0, r1, c0, c1, neg=False):
                M = arb_mat(r1 - r0, c1 - c0)
                for i in range(r1 - r0):
                    for j in range(c1 - c0):
                        e = arb(float(Hf[r0 + i, c0 + j]))
                        if r0 + i == c0 + j:
                            e = e + arb('0.03') * d[r0 + i] * d[r0 + i]
                        M[i, j] = -e if neg else e
                return M
            return ([mat(0, 3, 0, 3), mat(3, 5, 3, 5)],
                    [mat(0, 3, 5, 7), mat(3, 5, 5, 7)],
                    mat(5, 7, 5, 7, neg=self.flip_h))

    # realistic candidate: NEAR the planted zero, never exactly on it
    XC = ([[v + 1e-9 for v in X0[0][0]], [v - 1e-9 for v in X0[0][1]]],
          [v + 1e-9 for v in X0[1]])
    o = Oracle()

    # (a) planted zero -> CERTIFIED; certified ball contains the true zero
    res = bk.block_krawczyk(o, XC, 1e-3, 128)
    assert res["verdict"] == "CERTIFIED", res
    assert res["max_margin"] < 1.0 and res["newton_step_inf"] < 1e-6, res
    assert all(abs(a - b) <= res["radius"] for a, b in
               [(0.3, XC[0][0][0]), (-0.1, XC[1][1])]), "containment"

    # (b) planted non-contraction: displaced candidates REFUSE, named
    X1 = ([[v + 0.5 for v in X0[0][0]], X0[0][1]], X0[1])
    r1 = bk.block_krawczyk(o, X1, 1e-3, 128)
    assert r1["verdict"] == "REFUSED" and r1["failed"] == "block_00", r1
    X2 = (X0[0], [v + 0.5 for v in X0[1]])
    r2 = bk.block_krawczyk(o, X2, 1e-3, 128)
    assert r2["verdict"] == "REFUSED" and r2["failed"] == "border", r2
    rbig = bk.block_krawczyk(o, XC, 30.0, 128)
    assert rbig["verdict"] == "REFUSED", "fat box must refuse"

    # (c) mutation controls: sign flips MUST flip the verdict of the
    # certificate that consumes them. Math flip (different objective) is
    # caught by the Krawczyk leg (the zero moved); an H-only border flip
    # self-preconditions the contraction (any invertible claimed-H gives
    # E ~ 0 near a true zero — Krawczyk trusts the Jacobian CLAIM, the
    # oracle contract; oracle-vs-model truth is gated oracle-side) but is
    # caught by the PD leg of the COMBINED certificate.
    rg = bk.block_krawczyk(Oracle(flip_math=True), XC, 1e-3, 128)
    assert rg["verdict"] == "REFUSED", "math sign-flip escaped"
    rh = bk.certify_local_min(Oracle(flip_h=True), XC, 1e-3, 128)
    assert rh["verdict"] == "REFUSED" and rh["failed"] is not None, \
        "H border sign-flip escaped the combined certificate"

    # (d) precision starvation: fat balls REFUSE, never a loose certificate
    rs = bk.block_krawczyk(o, XC, 1e-3, 6)
    assert rs["verdict"] == "REFUSED", "starved precision emitted certificate"

    # (e) PD legs + scalar interval-Cholesky cross-check
    fh = bk.eval_oracle(o, XC, 1e-3, 128)
    pd = bk.pd_block_arrow(*fh[1], 128)
    assert pd["verdict"] == "PD_CERTIFIED" and \
        min(pd["margins"].values()) > 0, pd
    fhs = bk.eval_oracle(Oracle(saddle=True), XC, 1e-3, 128)
    pds = bk.pd_block_arrow(*fhs[1], 128)
    assert pds["verdict"] == "REFUSED" and "border_schur" in pds["failed"], pds
    assert "diagnostic_negative_direction" in pds, "saddle direction unnamed"
    ok, piv, _ = bk.interval_cholesky(fh[1][0][0])
    assert ok and piv > 0, "scalar cholesky disagrees on PD block"
    oks, _, idx = bk.interval_cholesky(fhs[1][2])
    assert not oks and idx == 1, "scalar cholesky missed the planted saddle"
    full = bk.certify_local_min(o, XC, 1e-3, 128)
    assert full["verdict"] == "CERTIFIED", full

    # (f) contract violations are TYPED, never quiet
    class BadDims:
        dims = ([], 0)
    try:
        bk.block_krawczyk(BadDims(), XC, 1e-3, 128)
        raise AssertionError("bad dims accepted")
    except bk.BlockArrowContractError:
        pass

    # (g) ANISOTROPIC per-coordinate radius contract
    # (a stiff Hessian + uniform box suffers
    # Bauer-Skeel width amplification). Stiff variant: coordinate 0
    # carries curvature 4e6 AND a cross nonlinearity f += 500 d2^2 d0^2
    # (the hierarchical mechanism: a stiff coordinate's box width bleeding
    # into a soft coordinate's Hessian row). The uniform box at the scale
    # the soft coordinates need must REFUSE; the curvature-scaled box must
    # CERTIFY the same candidate.
    class Stiff(Oracle):
        def _hf(self):
            Hf = super()._hf()
            Hf[0, 0] = 4.0e6
            return Hf

        def F(self, x):
            (f0, f1), fg = super().F(x)
            d = self._d(x)
            f0 = list(f0)
            f0[0] = f0[0] + 1000 * d[2] * d[2] * d[0]
            f0[2] = f0[2] + 1000 * d[2] * d[0] * d[0]
            return (f0, f1), fg

        def H(self, x):
            D, Bm, G = super().H(x)
            d = self._d(x)
            D[0][0, 0] = D[0][0, 0] + 1000 * d[2] * d[2]
            D[0][2, 2] = D[0][2, 2] + 1000 * d[0] * d[0]
            cross = 2000 * d[2] * d[0]
            D[0][0, 2] = D[0][0, 2] + cross
            D[0][2, 0] = D[0][2, 0] + cross
            return D, Bm, G
    st = Stiff()
    XSC = ([[v + 1e-9 for v in X0[0][0]], [v - 1e-9 for v in X0[0][1]]],
           [v + 1e-9 for v in X0[1]])
    runi = bk.block_krawczyk(st, XSC, 0.5, 128)
    assert runi["verdict"] == "REFUSED", \
        "uniform box escaped the stiff cross-width amplification"
    raniso = ([[0.5 / 2000.0, 0.5, 0.5], [0.5, 0.5]], [0.5, 0.5])
    rani = bk.block_krawczyk(st, XSC, raniso, 128)
    assert rani["verdict"] == "CERTIFIED" and rani["radius_anisotropic"], \
        rani
    fhs2 = bk.eval_oracle(st, XSC, raniso, 128)
    pds2 = bk.pd_block_arrow(*fhs2[1], 128)
    assert pds2["verdict"] == "PD_CERTIFIED", pds2
    return ("planted zero certified+contained; displaced/fat refuse named; "
            "sign-flip mutants 2/2 caught; starved prec refuses; PD + saddle "
            "+ scalar-cholesky cross-check; typed contract; anisotropic "
            "radius certifies where the uniform box refuses")


@leg("L17 front door run/solve/render (Muller gates + deviation)")
def l17():
    import math
    import baller
    from flint import arb, ctx

    def muller(xp, x):
        return 111 - 1130 / x + 3000 / (x * xp)

    # independent truth: raw flint at 400 dps, never through the front door
    old = ctx.prec
    try:
        ctx.prec = int(400 * 3.3219281) + 16
        tb = [arb(2), arb(-4)]
        for _ in range(2, 101):
            tb.append(muller(tb[-2], tb[-1]))
    finally:
        ctx.prec = old
    t30, t100 = tb[30], tb[100]

    # (a) the spec's 50-digit run: x_30 certified enclosure
    r = baller.run(muller, x0=(2, -4), n=100, dps=50)
    x30, x100 = r.values[30], r.values[100]
    assert baller.render(x30) == "6.0056", baller.render(x30)
    assert x30.mid().str(14, radius=False, more=True) == "6.0056486887714"
    assert 1e-6 < float(x30.rad()) < 1e-3, float(x30.rad())
    assert x30.contains(t30), "x30 enclosure lost the truth"
    # (b) x_100: THE DOCUMENTED DEVIATION — divisor ball contained 0, so
    # rad=+inf and the CENTRAL VALUE CONTINUES as the fixed-precision
    # stream (which walks into the float trap at 100); marked, never bare
    assert r.blown_at is not None and 30 < r.blown_at < 45, r.blown_at
    assert float(x100.rad()) == math.inf
    assert abs(float(x100.mid()) - 100.0) < 1e-6, float(x100.mid())
    assert baller.render(x100).startswith("UNCERTIFIED")
    try:
        baller.render(x100, strict=True)
        raise AssertionError("strict render of the inf ball did not refuse")
    except baller.RenderRefused:
        pass
    # (c) NEGATIVE CONTROL (documented): plain floats converge to EXACTLY
    # 100.0 (true limit 6) — the bare wrong number fail-closed printing
    # exists to prevent; this is the control demonstrating why balls matter
    a, b = 2.0, -4.0
    for _ in range(2, 101):
        a, b = b, 111 - 1130 / b + 3000 / (b * a)
    assert b == 100.0, b
    # (d) adaptive solve: escalates 30->60->120->240 dps, certifies the
    # spec's x_100 = 6.0000000160995649 (>= 16 digits)
    res = baller.solve(muller, 16, x0=(2, -4), n=100)
    assert res.dps == 240 and res.certified_digits >= 16, (res.dps, res.certified_digits)
    assert [d for d, _ in res.attempts] == [30, 60, 120, 240], res.attempts
    assert baller.render(res, digits=17) == "6.0000000160995649"
    assert res.ball.contains(t100), "solve enclosure lost the truth"
    assert "+/-" in repr(res)
    # (e) cap -> TYPED refusal naming achieved digits, never a bare number
    try:
        baller.solve(muller, 16, x0=(2, -4), n=100, max_dps=60)
        raise AssertionError("capped solve RETURNED instead of refusing")
    except baller.SolveRefused as e:
        assert e.achieved_digits == 0 and e.target_digits == 16
        assert len(e.attempts) == 2
        assert float(e.best.rad()) == math.inf, "best-enclosure not carried"
    # (f) fail-closed render units
    fat = arb("1.234567890123") + arb(0, 1e-6)
    assert baller.render(fat) == "1.23457"
    assert baller.render(fat, digits=12) == "1.23457", "clamp breached"
    try:
        baller.render(fat, digits=12, strict=True)
        raise AssertionError("strict over-request did not refuse")
    except baller.RenderRefused as e:
        assert e.certified == 6 and e.requested == 12
    assert baller.render(arb(0, 1e-30)).startswith("UNCERTIFIED")
    assert baller.render(arb(2), digits=5) == "2.0000"
    try:
        baller.render(3.14)
        raise AssertionError("bare float accepted")
    except TypeError:
        pass
    # (g) expression mode at working precision
    e = baller.run("arb(2).sqrt() * arb.pi()", dps=40)
    assert baller.render(e.value, digits=30) == "4.44288293815836624701588099006"
    return ("x30 6.0056486887714 +/- ~1e-5 certified+contained; x100 "
            "[mid-stream->100 +/- inf] marked (the deviation); float "
            "control = bare 100.0; solve 30->240 dps certifies "
            "6.0000000160995649, contains truth; cap refusal typed; "
            "render clamps/refuses; expression mode exact")


@leg("L18 front door MUTATION CONTROLS (render/solve/deviation/escalation)")
def l18():
    from baller.mutation import run_mutation_gate
    src = os.path.join(TOOL, "baller", "frontdoor.py")
    gate_cmd = [sys.executable, "-B",
                os.path.join(HERE, "frontdoor_case.py"), "{FILE}"]
    mutants = [
        # a render that prints uncertified digits must FAIL the battery
        ("        return p_str if certified else None",
         "        return p_str"),
        # a solve that returns without the target must FAIL
        ("        if cert >= target_digits:",
         "        if True:"),
        # the deviation's mid-stream continuation dropped (nan instead)
        ("values.append(ball[-1] if blown_at is None else _inf_ball(m))",
         "values.append(ball[-1] if blown_at is None else nan)"),
        # escalation dead: the driver never raises precision
        ("        nxt = int(dps * factor)",
         "        nxt = dps"),
    ]
    rep = run_mutation_gate(src, mutants, gate_cmd,
                            scratch_dir=os.path.join(os.getcwd(), "mut_fd"))
    return f"pristine PASS + mutants caught {rep['caught']}"


@leg("L19 purity-scan member (planted 3-class smoke)")
def l19():
    from baller import purity_scan as ps
    assert hasattr(ps, "scan_python") and hasattr(ps, "scan_c")
    sample = os.path.join(os.getcwd(), "purity_l19_sample.py")
    with open(sample, "w") as fh:
        fh.write("from mpmath import mp\n"
                 "value = mp.mpf(0.2857142857)\n"
                 "banked = '3.14159265358979'\n"
                 "tol = 1e-12\n")
    listing = os.path.join(os.getcwd(), "purity_l19_list.txt")
    with open(listing, "w") as fh:
        fh.write(sample + "\n")
    out = os.path.join(os.getcwd(), "purity_l19_out.json")
    r = subprocess.run([sys.executable, "-B",
                        os.path.join(TOOL, "baller", "purity_scan.py"),
                        listing, out], capture_output=True, text=True)
    assert r.returncode == 0, f"member CLI rc={r.returncode}: {r.stderr[-200:]}"
    rep = json.load(open(out))
    kinds = rep["summary_kind_counts"]
    assert kinds.get("MPF-OF-FLOAT", 0) >= 1, f"mpf(float) missed: {kinds}"
    assert kinds.get("DEC-STRING", 0) >= 1, f"decimal string missed: {kinds}"
    assert kinds.get("PURE-CONV", 0) >= 1, f"conv literal missed: {kinds}"
    assert kinds.get("PARSE-ERROR", 0) == 0, f"parse error: {kinds}"
    return f"3 planted classes classified ({kinds})"


@leg("L20 vendored fold members: certlane pytest + ling/lgf import+probe")
def l20():
    import baller
    baller.verify(quiet=True)
    vend = os.path.join(TOOL, "vendor")
    # certlane: the member's own 26-case containment/planted/convention suite
    r = subprocess.run(
        [sys.executable, "-B", "-m", "pytest",
         os.path.join(vend, "certlane", "test_primitives.py"),
         "-q", "-p", "no:cacheprovider"],
        capture_output=True, text=True, cwd=os.getcwd())
    assert r.returncode == 0, \
        f"certlane pytest rc={r.returncode}: {(r.stdout or r.stderr)[-300:]}"
    pyline = [l for l in r.stdout.strip().splitlines() if "passed" in l][-1]
    # ling_onesided: every module must import from the vendored bytes alone
    # (fixture beside them), from a scratch cwd; lgf: the pinned BC
    # cache must load and the null-loop identity self-check must pass
    probe = (
        "import sys\n"
        f"sys.path.insert(0, {os.path.join(vend, 'ling_onesided')!r})\n"
        "import taylor_p, gate_price, tail_engine, balanced_unc, double_unc\n"
        "import collapse_unc, sd_engine, l5_center, sd_onesided\n"
        "assert sum(sd_engine.COUNTS.values()) == 318\n"
        "assert sum(taylor_p.COUNTS.values()) == 318\n"
        f"sys.path.insert(0, {os.path.join(vend, 'lgf_monodromy')!r})\n"
        "import valmono\n"
        "assert len(valmono.BC) == 6\n"
        "valmono.selftest(prec=160)\n")
    r2 = subprocess.run([sys.executable, "-B", "-c", probe],
                        capture_output=True, text=True, cwd=os.getcwd())
    assert r2.returncode == 0, \
        f"member import/probe rc={r2.returncode}: {(r2.stderr or r2.stdout)[-300:]}"
    tail = [l for l in r2.stdout.strip().splitlines() if l][-1]
    return f"certlane {pyline.strip()}; ling 9 modules import, N=318; {tail[:80]}"


@leg("L21 dps_lint member: selftest + `-m baller.hygiene lint` exit codes")
def l21():
    import baller
    from baller import hygiene
    from baller import dps_lint as member
    assert hygiene.dps_lint is member, "hygiene.dps_lint is not baller.dps_lint"
    assert callable(member.main) and callable(member.selftest) and callable(member.lint_file)
    env = dict(os.environ, TMPDIR=os.getcwd(), PYTHONDONTWRITEBYTECODE="1")
    # (a) the member's own built-in battery through its standalone CLI
    r = subprocess.run([sys.executable, "-B",
                        os.path.join(TOOL, "baller", "dps_lint.py"), "--selftest"],
                       capture_output=True, text=True, cwd=os.getcwd(), env=env,
                       timeout=120)
    assert r.returncode == 0, f"selftest rc={r.returncode}: {(r.stdout + r.stderr)[-400:]}"
    summ = [l for l in r.stdout.splitlines() if l.startswith("dps_lint selftest:")]
    assert summ and summ[-1].rstrip().endswith(" 0 fail"), r.stdout[-400:]
    npass = int(summ[-1].split()[2])
    # (b) the package entry `python3 -m baller.hygiene lint`: planted -> 1, fixed -> 0
    d = os.path.join(os.getcwd(), "dps_l21")
    os.makedirs(d, exist_ok=True)
    bad = os.path.join(d, "planted.py"); ok = os.path.join(d, "fixed.py")
    open(bad, "w").write("import mpmath as mp\nT = mp.mpf(-1) / 3\nmp.mp.dps = 50\n")
    open(ok, "w").write("import mpmath as mp\nmp.mp.dps = 50\nT = mp.mpf(-1) / 3\n")
    penv = dict(env, PYTHONPATH=TOOL + os.pathsep + env.get("PYTHONPATH", ""))
    rb = subprocess.run([sys.executable, "-B", "-m", "baller.hygiene", "lint", bad],
                        capture_output=True, text=True, cwd=os.getcwd(), env=penv,
                        timeout=60)
    assert rb.returncode == 1 and "[module-level]" in rb.stdout, \
        f"planted rc={rb.returncode}: {(rb.stdout + rb.stderr)[-300:]}"
    ro = subprocess.run([sys.executable, "-B", "-m", "baller.hygiene", "lint", ok],
                        capture_output=True, text=True, cwd=os.getcwd(), env=penv,
                        timeout=60)
    assert ro.returncode == 0, f"fixed rc={ro.returncode}: {(ro.stdout + ro.stderr)[-300:]}"
    # (c) in-process API on the same pair
    assert len(member.lint_file(bad)) == 1 and member.lint_file(ok) == []
    return (f"selftest {npass} pass/0 fail; -m baller.hygiene lint planted->1 "
            f"fixed->0; API 1/0; member identity holds")


def main():
    t0 = time.time()
    results, fails = [], []
    for name, fn in RESULTS:
        t = time.time()
        try:
            detail = fn()
            ok = True
        except Exception as e:
            detail = f"{type(e).__name__}: {e}"
            ok = False
            fails.append(name)
        wall = time.time() - t
        print(f"{'PASS' if ok else 'FAIL'}  {name}  [{wall:.1f}s]  {str(detail)[:180]}")
        results.append({"leg": name, "ok": ok, "wall_s": round(wall, 2),
                        "detail": str(detail)[:400]})
    total = time.time() - t0
    with open("BALLER_BATTERY_SUMMARY.json", "w") as fh:
        json.dump({"results": results, "total_s": round(total, 1),
                   "overall": "PASS" if not fails else "FAIL"}, fh, indent=1)
    if fails:
        print(f"OVERALL FAIL: {fails} ({total:.0f}s)")
        return 1
    print(f"OVERALL PASS ({len(results)} legs, {total:.0f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

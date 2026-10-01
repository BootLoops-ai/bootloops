#!/usr/bin/env python3
"""pmflow smoke battery — what runs without a built amflow_cli:

  1. CLI surface: pmflow.py --help exits 0 and lists all six subcommands;
     the solve subparser exposes --anchors/--n-eps.
  2. ratio salvage on a synthetic AMFLOW_DUMP_EPS_GRID log: every master at a
     node carries a common node-erratic wild factor (the documented failure
     mode of raw eps-fits); the planted rational eps^0 ratios
     m1/m0 = -2/5 and m2/m0 = 22/7 must come back EXACTLY through
     gf_eps0_ratios.py — the factor cancels in the ratio and the Lagrange
     eps->0 limit is exact on a polynomial ratio.
  3. fail-closed: gf_eps0_ratios.py --selftest without $GF_SELFTEST_LOG must
     refuse loudly (the reference log is not shipped).
  4. anchored closure solve on a synthetic 2-master fixed-point fixture:
     gf_solve.py --anchors "1:J63" must reproduce the planted solution
     (anchored master = J63(eps), null-space master closed to 2*J63(eps) by
     the second-flow rows) with an all-consistent row census.
  5. --anchors refusals: unknown FORM, sector index outside the slot count,
     sector absent from the fpA basis, malformed entry, and the built-in
     reference-family default on a foreign family all refuse loudly with
     named errors (never a bare traceback from basis.index).
  6. GRAVITYFLOW_FIXEDPOINT surface: the detect-side regex parses the
     engine guard's structured abort line (family + depth) and rejects the
     CUTREGION line; when the engine source ships in-tree, the emit
     literals in amfsystem.cpp must match what the regex consumes.
  7. engine-source syntax gate: the patched amfsystem.cpp passes
     `g++ -fsyntax-only` — skipped BY NAME when the engine tree, g++, or
     the flint headers are absent.

The engine legs (detect/discover/respond/solve/inject/map against a built
amflow_cli) are not run here.  Exit 0 only if every non-skipped leg behaves
as documented.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
fails = []


def grade(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        fails.append(name)


# 1. CLI surface
r = subprocess.run([sys.executable, os.path.join(HERE, "pmflow.py"), "--help"],
                   capture_output=True, text=True)
grade("pmflow CLI surface (--help, six subcommands)",
      r.returncode == 0 and all(s in r.stdout for s in
      ("detect", "discover", "respond", "solve", "inject", "map")))
r = subprocess.run([sys.executable, os.path.join(HERE, "pmflow.py"),
                    "solve", "--help"], capture_output=True, text=True)
grade("solve subparser exposes --anchors/--n-eps",
      r.returncode == 0 and "--anchors" in r.stdout and "--n-eps" in r.stdout)

# 2. synthetic ratio salvage with planted rational eps^0 ratios
mp.mp.dps = 130
eps = [mp.mpf(1) / d for d in (100, 128, 160, 200, 256, 320, 400, 512)]
ratios = {1: lambda e: mp.mpf(-2) / 5 + 3 * e + 7 * e ** 2,   # -> -2/5
          2: lambda e: mp.mpf(22) / 7 - e + 5 * e ** 3}       # -> 22/7
tmp = tempfile.mkdtemp(prefix="pmflow_selftest_")
log = os.path.join(tmp, "synthetic_eps_grid.log")
with open(log, "w") as fh:
    for j, e in enumerate(eps):
        fh.write(f"[EPS_GRID] eps[{j}] = {mp.nstr(e, 120)}\n")
    for j, e in enumerate(eps):
        wild = mp.exp(mp.sin(17 * (j + 1)) * 20)  # node-erratic common factor
        fh.write(f"[EPS_GRID] M[0][{j}] = {mp.nstr(wild, 120)}\n")
        for row, f in ratios.items():
            fh.write(f"[EPS_GRID] M[{row}][{j}] = {mp.nstr(f(e) * wild, 120)}\n")
out = os.path.join(tmp, "table.json")
r = subprocess.run([sys.executable, os.path.join(HERE, "gf_eps0_ratios.py"),
                    "--log", log, "--rows", "0-2", "--n-eps", "8",
                    "--dps", "120", "--out", out],
                   capture_output=True, text=True)
tab = json.load(open(out)) if r.returncode == 0 and os.path.isfile(out) else {}
got = {k: (tab.get(k, {}).get("p"), tab.get(k, {}).get("q")) for k in ("m1", "m2")}
ok = r.returncode == 0 and got["m1"] == ("-2", "5") and got["m2"] == ("22", "7")
grade("ratio salvage recovers the planted rationals exactly", ok,
      f"(m1/m0={got['m1'][0]}/{got['m1'][1]}, m2/m0={got['m2'][0]}/{got['m2'][1]})")
if not ok:
    print((r.stdout + r.stderr)[-800:])

# 3. fail-closed selftest without the (unshipped) reference log
env = {k: v for k, v in os.environ.items() if k != "GF_SELFTEST_LOG"}
r = subprocess.run([sys.executable, os.path.join(HERE, "gf_eps0_ratios.py"),
                    "--selftest"], capture_output=True, text=True,
                   cwd=tmp, env=env)
grade("--selftest fails closed without $GF_SELFTEST_LOG",
      r.returncode != 0 and ("FileNotFoundError" in r.stderr
                             or "No such file" in r.stderr))

# 4. anchored closure solve on a synthetic 2-master fixed-point fixture.
#    Family "toy", 2 slots, basis [(1,0),(1,1)].  One J63-kind probe key
#    (unit column on master 0) and one master-kind key on index 1 whose
#    unit column makes master 1 a free direction; the second flow's
#    A_alpha = [[1,0],[1,1/2]] pins M1 = 2*M0 through the row census.
#    With --anchors "1:J63" (sector 1 -> corner (1,0) -> index 0) the
#    solution is M0 = J63(eps), M1 = 2*J63(eps), all rows consistent.
mp.mp.dps = 60
toy = tempfile.mkdtemp(prefix="pmflow_selftest_anchor_")
toy_probes = os.path.join(toy, "probes")
os.makedirs(toy_probes)
toy_eps = [mp.mpf(1) / 100, mp.mpf(1) / 50, mp.mpf(3) / 100]
json.dump({"family": "toy",
           "preferred": ["toy|1|0", "toy|1|1"],
           "sub_masters": ["toy|1|0", "toy|1|1"],
           "A": [[["1", "0"], ["1", "0.5"]] for _ in toy_eps]},
          open(os.path.join(toy, "fp_A.json"), "w"))
json.dump(["toyend|0|1", "toyend|1|1"],
          open(os.path.join(toy_probes, "keys.json"), "w"))
json.dump({"toyend|0|1": ["J63", None], "toyend|1|1": ["master", 1]},
          open(os.path.join(toy, "circular_map.json"), "w"))
for k, col in enumerate(([1, 0], [0, 1])):
    with open(os.path.join(toy_probes, f"log_{k:02d}.log"), "w") as fh:
        for j, e in enumerate(toy_eps):
            fh.write(f"[EPS_GRID] eps[{j}] = {mp.nstr(e, 50)}\n")
        for m_ in (0, 1):
            for j in range(len(toy_eps)):
                fh.write(f"[EPS_GRID] M[{m_}][{j}] = "
                         f"{mp.nstr(mp.mpf(col[m_]), 30)}\n")
toy_out = os.path.join(toy, "fp_boundary.json")


def run_solve(*extra):
    return subprocess.run(
        [sys.executable, os.path.join(HERE, "gf_solve.py"),
         "--fpA", os.path.join(toy, "fp_A.json"), "--probes", toy_probes,
         "--cmap", os.path.join(toy, "circular_map.json"), "--out", toy_out,
         "--n-eps", "3", "--dps", "60"] + list(extra),
        capture_output=True, text=True)


r = run_solve("--anchors", "1:J63")
anchored_ok = r.returncode == 0 and os.path.isfile(toy_out)
if anchored_ok:
    tab = json.load(open(toy_out))

    def j63(e):
        return -(mp.mpc(0, 1) / 8) * mp.pi ** mp.mpf('-1.5') \
            * mp.gamma(e - mp.mpf('0.5')) ** 3

    for j, e in enumerate(toy_eps):
        for key, mult in (("toy|1|0", 1), ("toy|1|1", 2)):
            v = mp.mpc(mp.mpf(tab[key][j]["re"]), mp.mpf(tab[key][j]["im"]))
            want = mult * j63(e)
            if abs(v - want) / abs(want) > mp.mpf("1e-40"):
                anchored_ok = False
grade("--anchors closure solve reproduces the planted solution "
      "(J63, 2*J63; census all-consistent)", anchored_ok)
if not anchored_ok:
    print((r.stdout + r.stderr)[-800:])

# 5. --anchors refusals: named errors, never a bare traceback.
refusals = (
    (("--anchors", "1:FOO"), "unknown anchor form"),
    (("--anchors", "7:J63"), "does not fit a 2-slot family"),
    (("--anchors", "2:J63"), "not among the 2 fpA preferred masters"),
    (("--anchors", "garbage"), "malformed entry"),
    ((), "does not fit a 2-slot family"),   # reference default, foreign family
)
ref_ok = True
for extra, needle in refusals:
    r = run_solve(*extra)
    if not (r.returncode != 0 and needle in r.stderr
            and "Traceback" not in r.stderr):
        ref_ok = False
        print(f"  refusal miss for {extra!r}: rc={r.returncode} "
              f"stderr={r.stderr[-200:]!r}")
grade("--anchors refuses loudly (bad form/sector/basis/spec + foreign "
      "default)", ref_ok)

# 6. GRAVITYFLOW_FIXEDPOINT surface: detect-side regex vs the engine emit.
sys.path.insert(0, HERE)
import pmflow  # noqa: E402

fp_line = ("AMFSystem::setup: GRAVITYFLOW_FIXEDPOINT: boundary family "
           "'pmfam_b0_r0_f0' at depth 1 reproduces its parent's "
           "propagator set (15 propagators identical after conservation)")
m_fp = pmflow.FIXEDPOINT_RE.search(fp_line)
cut_line = ("AMFSystem::build_boundary: GRAVITYFLOW_CUTREGION: family "
            "'pmfam' carries a cut quadratic propagator")
parse_ok = (m_fp is not None and m_fp.group(1) == "pmfam_b0_r0_f0"
            and m_fp.group(2) == "1"
            and pmflow.FIXEDPOINT_RE.search(cut_line) is None
            and pmflow.CUTREGION_RE.search(fp_line) is None)
grade("FIXEDPOINT_RE parses the guard line (family+depth), rejects "
      "CUTREGION", parse_ok)

# Engine source tree: the AMFlow.cpp fork lives in the sibling repository
# amflow-cpp. AMFLOW_CPP_SRC names a checkout of it; when
# unset, the default is a checkout beside this repository
# (<repo root>/../amflow-cpp). Legs 6 and 7 skip by name when absent.
ENGINE_DIR = os.environ.get("AMFLOW_CPP_SRC") or os.path.normpath(
    os.path.join(HERE, "..", "..", "..", "amflow-cpp"))
ENGINE_CPP = os.path.join(ENGINE_DIR, "src", "pipeline", "amfsystem.cpp")
if os.path.isfile(ENGINE_CPP):
    src = open(ENGINE_CPP, errors="replace").read()
    # The emit is a stream expression; pin the literals the regex consumes
    # around the family name and the depth, plus the debug bypass knob.
    lit_ok = ("GRAVITYFLOW_FIXEDPOINT: boundary " in src
              and "family '" in src and "' at depth " in src
              and "AMFLOW_ALLOW_FIXEDPOINT_RECURSION" in src)
    grade("engine emit literals match the parse surface", lit_ok)
else:
    print("[SKIP] engine emit literal pin: engine source not present "
          f"({ENGINE_CPP} not found; set AMFLOW_CPP_SRC to a checkout of "
          "the sibling repository amflow-cpp)")

# 7. engine-source syntax gate (skips BY NAME when toolchain/tree absent).
gxx = shutil.which("g++")
flint_ok = any(os.path.isfile(os.path.join(d, "flint", "fmpz.h"))
               for d in ("/usr/include", "/usr/local/include"))
if not os.path.isdir(ENGINE_DIR):
    print("[SKIP] engine syntax gate: engine source tree absent "
          f"({ENGINE_DIR}; set AMFLOW_CPP_SRC to a checkout of "
          "the sibling repository amflow-cpp)")
elif gxx is None:
    print("[SKIP] engine syntax gate: g++ not on PATH")
elif not flint_ok:
    print("[SKIP] engine syntax gate: flint headers not found")
else:
    r = subprocess.run(
        [gxx, "-std=c++17", "-fsyntax-only", "-Iinclude",
         "-I/usr/include/flint", "-I/usr/local/include/flint",
         os.path.join("src", "pipeline", "amfsystem.cpp")],
        capture_output=True, text=True, cwd=ENGINE_DIR)
    grade("engine amfsystem.cpp passes g++ -fsyntax-only", r.returncode == 0)
    if r.returncode != 0:
        print(r.stderr[-1200:])

if fails:
    print(f"OVERALL: FAIL {fails}")
    sys.exit(1)
print("OVERALL: PASS (smoke; engine legs need a built amflow_cli)")

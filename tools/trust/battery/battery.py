#!/usr/bin/env python3
"""TRUST acceptance battery (assembly-package pattern; every wall MEASURED
before wiring).

V-style: refuses tool/protected/parked cwds; runs from run scratch; OVERALL
states COVERAGE. Every leg mutation-controlled (break the gated thing in a
scratch copy -> the leg's own grade FAILS; receipts in
MUTATION_RECEIPTS.json).

ENV PERIMETER: a PYTHONPATH stdlib-
shadow (pins-aware hashlib) greened a tampered vendor — battery.py imports
hashlib before any check, so PYTHONPATH/PYTHONSTARTUP/PYTHONHOME worlds are
REFUSED at the preamble (rc=2, nothing runs; membership test —
set-but-empty vars refuse too); all PYTHON* env is scrubbed for children and
every python child launches -E -P -s -B (-s means the
user-writable user site is never site-processed, so its .pth files never
execute; children needing sympy/flint get the user-site dir appended at the
TAIL of sys.path by an inert bootstrap — CHILD_BOOT; -I/-S call sites keep
the plain maximum-isolation launch). The SAME refusal runs at
trust/__init__ IMPORT time (typed EnvPoisonError; -E worlds exempt — the
interpreter ignored the vars), so the package's in-process code paths are no
longer the battery's asymmetric blind side (T8 probes both directions).
trust/__init__ ALSO authenticates the security-critical
stdlib modules post-import (typed StdlibShadowError) — the env-FREE
sys.path[0]/cwd shadow class (T8 probes it with a script-dir pins-aware
hashlib + a gate-off mutation control). That authentication
is against STDLIB roots ONLY (BASE-installation sysconfig stdlib/platstdlib
+ the site-segment belt): site-packages is NOT a stdlib home, so the
site-.pth reorder class — a 2-line .pth reorders sys.path at interpreter
startup, before any user code — refuses too (T8 probes it with a scratch
venv + .pth plant + a roots-widened mutation control).

Legs:
  T0 PACKAGE pins — trust/*.py sha-pinned BATTERY-SIDE, verified BEFORE any
     trust import (the battery-side pin pattern); pycache purged; the
     trust/ DIRECTORY must contain the pinned files and NOTHING else
     (dirs/symlinks/strays refuse); TOOL root + battery/ swept for importable
     strays; vendor tree walked against the pin table SINGLE-SOURCED from
     the sha-verified trust/_pins.py (ast-parsed pre-import, nothing
     executed; an unauthenticated pin source refuses by name).
  T1 vendor pins 4/4 + LINKED homes + LINKED_SHAS
     advisory baselines (verify() compares and reports
     linked_drift — loud, never fatal); the imported trust.PINS cross-
     checked == the T0-parsed table (an import that resolves a different pin
     table than the verified file refuses here).
  T2 tamper control (T1's mutation control): 6 scratch COPIES of the package,
     child-mode verify — byteflip / unpinned-extra / symlink-swap (byte-
     identical content!) / pycache-content each raise VendorTamperError;
     pristine copy PASSES (positive control); byteflipped copy also refuses
     at run_vendored LAUNCH; linkdrift copy (flipped LINKED_SHAS baseline)
     verifies WITHOUT raising but REPORTS the drift (the advisory promise).
  T3 lp-syz restore-verify (shipped scope): the engines' BUILT-IN control
     stages run child-mode and are graded battery-side with NO external
     receipts — tadpole closed-form control (11 cols/60 rels, value-level
     vs the env-pointed truth table; NAMED SKIP when TRUST_LPSYZ_KTAB is
     unset), planted-violation positive control ctrl45 (EXACTLY ONE
     violation, named row supported on master columns only, 0 free
     non-masters) and clean-port negative control ctrl452 (0 violations by
     design). Structural counts (relations/columns/free) follow the
     Singular syzygy basis, so they are RECORDED, not pinned; mutants are
     graded against the SAME-RUN pristine receipts. The source control
     receipts and the sec63 value-level reference are not included in the
     release, so sec63 is a NAMED SKIP.
     MUTATIONS (the planted lp-syz traps, wired): Aut(G)-quotient-skip in
     scratch lp_syz_431.py -> ctrl45 design grade FAILS (violations 1->0,
     measured) AND ctrl452 differs from pristine (cols 203->378 measured);
     descent-sign flip in the port descent block = the MEASURED BLINDSPOT
     demo (receipts field-identical to pristine on BOTH shipped controls —
     the class only a value-level grade catches; recorded informational;
     its retired catcher is the T3m1 named skip); grader-teeth dict
     mutations.
  T4 joint three-check fixture (pilot B + fixture receipts): all pilot/
     fixture receipts + the kira table sha-pinned battery-side and the
     agreement bars parsed FIELD-EXACT (84/84 CERTIFIED, 0 CLAIM_MISMATCH,
     2 primes; 42/42 EXACT at 2 rational points, sigma=-1, 0 violations);
     LIVE bounded re-executions: P1 {2,5,6,7} tower 7/7 vs kira (4.0 s
     measured) + fresh adapter re-parse and 84-witness re-verify through
     BOTH receipt cores (0.6 s measured; emit-dir witness rollup sha-pinned)
     + INDEPENDENT stdlib re-parse spot-check: a
     sha-seeded FOREIGN subsample of SYSTEM_*.gz rows re-derived battery-
     side (gzip text -> ast-based modular ints, NO loader/CoeffEvaluator
     import) must equal the adapter rows the dual-core verify used — the
     named STRATA-LOADER SUBSTRATE caveat's wired mitigation.
     MUTATIONS: kira-table coefficient flip in a scratch table copy ->
     single-target compare refuses/fails; witness lam tamper -> verify
     FAILS; bar-grader teeth; adapter-row value flip -> re-parse comparator
     FAILS.
  T5 witness v1.0 round-trip (pilot C as a leg): 14 lp_syz witnesses
     (7 targets x 2 primes) emitted fresh (41.5 s measured), verified by
     BOTH cores; SCHEMA-BYTE assertions (R5: v1.0 UNCHANGED) — key sets vs
     the v1.0 spec AND vs a pinned pilot-B eliminator witness (lp_syz keys
     == eliminator keys minus optional `labels`); receipt_version/kind/
     identity byte-exact; syzygy provenance ONLY in system.source.
     MUTATIONS: emit-refusal (perturbed claim -> ValueError, nothing
     written) + planted faults lam/c/fingerprint 3/3 MUST-FAIL.
     Stub-detector floor is CPU-time, not wall (wall floors
     calibrated near the loaded wall threw a false red on a quiet box —
     6.4-41.5 s swing measured for identical work).
  T6 FRACTION-ORACLE leg (R3): pure-stdlib child (python3 -I -S), loads
     fraction_oracle.py BY PATH with its own sha assert, then: leg1 anchor
     reduce_column == flint-route reduction (15.5 s measured); leg2 FOREIGN
     rational point (103/31,7/13): pure-Fraction RE-ELIMINATION of the
     tower's exact rows -> I-space vs kira (37.2 s measured in-process
     prior; SCOPE-HONEST wording: the tower build + master/
     free classification + kira-side eval are shared with the standard
     flint route — the head-to-head FLINT cross-gate is leg1); leg3
     sha-seeded FOREIGN row subsample (k=25) exact->mod-p == the
     witness-system rows; leg4 sys.modules sweep — NOTHING beyond stdlib +
     the path-loaded module.
     MUTATION: corrupted flint reduction in a scratch export copy -> child
     leg1 equality FAILS (~15 s).
  T7 identity-front integrity: strata/receipt fronts resolve to the LIVE
     registered homes; cores byte-identical; the ONE-LINEAGE caveat +
     arbitration table asserted present in trust.strata.__doc__; the
     STRATA-LOADER SUBSTRATE caveat (check 2 as executed) asserted present
     in trust.receipt.__doc__ AND MANUAL.md.
     MUTATIONS: pre-imported shadow strata -> alias refuses (identity
     violation); sys.modules spoof with in-tree __file__ but no matching
     spec -> spec-authentication refuses (child-mode, both).
  T7C lineage census (the SCOPE-promised independence leg): import-graph census (ast, battery-side) — check-2 CORE
     stdlib-only; check-3 engine files (vendored lpsyz + oracle_k2disp +
     fraction_oracle + witness_bridge) import NO strata/loader/
     fp_eliminate/ibplapper code; the adapters/strata.py -> loader +
     fp_eliminate edge (the check-2 as-executed substrate) is FOUND and
     DOCUMENTED as the named caveat, never presented as disjoint; sha
     census: lineage file-sets pairwise byte-disjoint (sole designed
     exception: the two byte-identical receipt cores, ONE lineage).
     MUTATION: injected fp_eliminate import in a scratch oracle copy ->
     census flags it.
  T8 refusal battery: env-less TRUST_OUT_ROOT refuses TYPED (leaf
     OutputRootError, names the var); forbidden VALUES refuse (relative
     paths, parked-archive/home/protected trees, the tools tree, the data
     root, /, undeclared-neutral —
     ALLOWLIST posture, symlink-resolved escape); positive
     control: scratch value ACCEPTED; RELOCATED package copy refuses its
     OWN tree but accepts run scratch (the guard travels);
     run_vendored cwd escape refuses; unknown engine -> KeyError; junk
     fixture -> reduce_and_compare raises TYPED (ValueError; the old
     0-entries success-shaped path is CLOSED), absent target -> KeyError;
     missing file -> FileNotFoundError. Import-perimeter probes: in-process trust
     import in a PYTHONPATH world refuses typed EnvPoisonError (child, no
     -E) with the -E exemption as positive control (F1); a crafted evil
     kira table refuses ValueError WITHOUT executing its payload (F2,
     sympify retired for a locked rational-function grammar); env-less and
     out-of-root witness emission refuse OutputRootError, nothing written
     (F5); reserved provenance keys (producer/schema_note) refuse ValueError
     (F4) with an in-root stamped-producer positive control. Stdlib-identity
     probes: env-FREE script-dir/cwd stdlib shadow refuses typed
     StdlibShadowError naming the offending path, clean-dir positive
     control, gate-off copy mutation (F1); parse-coverage refusals — term
     lines/entries the grammar cannot consume refuse typed NAMING line
     numbers instead of silently certifying the parsed subset (dropterm/
     dropent/spaced-payload legs + the pinned 42-entry table as positive
     control + coverage-off copy mutation) (F2); in-grammar resource bombs
     refuse fast + typed (2^10^18 in ms, over-cap exponent, div0-at-build,
     eval-point pole; in-cap positive control + cap-off copy mutation)
     (F3). Site/.pth probes: a .pth-reordered site-packages hashlib plant
     (scratch venv) refuses typed StdlibShadowError naming the site dir,
     clean-venv positive control, roots-widened copy mutation (F1); the
     rebuild bombs INSIDE the R3 caps — 18KB Mult chain, 4KB div chain,
     22KB symbolic-exponent merge, 44KB Add flood, eager Pow-over-Mul-head
     — all refuse typed on CPU floors, eval-side merged-exponent refusal,
     budget-off copy mutation (F2); byte-alphabet refusals — U+2192 arrow/
     fullwidth-digit/lookalike-famname tables refuse typed naming byte
     offset + line BEFORE parsing, real non-ASCII file refuses, CRLF
     positive control, alphabet-off copy mutation showing the silent
     2-entry subset certification (F3).

DECLARED-UNWIRED (counted, with the measured number — legs over the 10-min
budget are declared, not wired):
  U1 pilot-B emit table-mode re-run: 997.5 s measured LOADED
     (eliminator solve 576.9+415.1 s of it; unloaded prior 2.2-2.3
     s/(row,prime)) — receipts sha-pinned + bars parsed (T4), 84-witness
     dual-core re-verify wired live instead (0.6 s).
  U2 FULL {1..7} 42-target lp_syz re-run at 2 points: 269.2 + 369.5 s
     measured (638.7 s combined) — receipts sha-pinned + bars parsed
     (T4), P1 sub-tower wired live instead (4.0 s).
  U3 strata full 4-slice 2-prime dictionary certification: reference run
     1493-3400 s/slice class — identity-linked; strata's own battery stays
     the certification for its bytes.
  U4 wrong-table MUST-FAIL: reference run 6/6 wrong + 6/6 controls in the
     receipt battery — linked by identity, not duplicated here.
  U5 lp_syz_prod execution: 8-var top-node syz RIGHT-CENSORED 1800 s hard
     — files pinned (T1); execution not wired into this battery.
  U6 adversarial hardening: separately gated before final acceptance.
  U7 m2-symbolic-route witnesses (sec503 class): winnow reference run
     (solve 176.7 s, 8/8) — fixture here is numeric-mass/(d,eta)-rational;
     scope pin restated.

Run: python3 battery/battery.py <scratch-dir>
Wall: ~30 s measured on a fast quiet workstation (T3 ~4 s, T4 ~9 s, T8 ~8 s
incl. the venv leg); plan minutes loaded.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

# --- env perimeter  --------------------------------
# A PYTHONPATH/PYTHONHOME/PYTHONSTARTUP world can shadow the stdlib (pins-
# aware hashlib PoC greened a tampered vendor). battery.py has ALREADY
# imported hashlib by the time any check could run, so the only sound
# response is refusal — scrubbing after the fact cannot un-import a shadow.
# MEMBERSHIP, not truthiness — set-but-empty PYTHONPATH refuses
# too (its path-inertness is undocumented interpreter behavior).
_BAD_ENV = [k for k in ("PYTHONPATH", "PYTHONSTARTUP", "PYTHONHOME")
            if k in os.environ]
if _BAD_ENV:
    print(f"REFUSED: {_BAD_ENV} set in the environment — the stdlib-shadow "
          "class defeats sha pinning (TRUST BLOCKING). Unset and "
          "rerun from a clean environment.")
    sys.exit(2)
# children must not inherit interpreter-config surprises either
for _k in [k for k in os.environ if k.startswith("PYTHON")
           and k != "PYTHONDONTWRITEBYTECODE"]:
    os.environ.pop(_k, None)

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
PKG = os.path.join(TOOL, "trust")
# scratch copies resolve the linked receipt root through this env (default: the real
# receipt MEMBER tree beside the package, <tool>/receipt) so the linked-drift
# advisory stays definable there (a scratch copy holds trust/ + vendor/ only).
os.environ.setdefault("TRUST_RECEIPT_ROOT", os.path.join(TOOL, "receipt"))
VENDOR = os.path.join(TOOL, "vendor")

# ---------------------------------------------------------------- guards
# cwd deny list: the tool tree (+ parent) and home dirs by default; extend
# with TRUST_BATTERY_DENY_ROOTS (colon-separated) to protect additional
# trees (env-configured; fail-closed refusal behavior unchanged).
BANNED = (TOOL, os.path.dirname(TOOL), "/home") + tuple(
    r for r in os.environ.get("TRUST_BATTERY_DENY_ROOTS", "").split(":") if r)
_cwd = os.path.realpath(os.getcwd())
if any(_cwd == os.path.realpath(b) or _cwd.startswith(os.path.realpath(b) + os.sep)
       for b in BANNED):
    print(f"REFUSED: run from a scratch dir, not {_cwd}")
    sys.exit(2)

OUT = sys.argv[1] if len(sys.argv) > 1 else None
if not OUT:
    print("usage: battery.py <scratch-dir>")
    sys.exit(2)
OUT = os.path.realpath(os.path.abspath(OUT))
if any(OUT == os.path.realpath(b) or OUT.startswith(os.path.realpath(b) + os.sep)
       for b in BANNED):
    print(f"REFUSED: scratch dir {OUT} is inside a banned tree")
    sys.exit(2)
os.makedirs(OUT, exist_ok=True)

# --- battery-side pins (the battery pins the package it is about to execute;
# battery.py cannot pin itself — it is the operator-run trust root. Values are
# re-emitted from the shipped bytes by script whenever a pinned file changes,
# never typed.)
PACKAGE_PINS = {
    "__init__.py": "e32c30b5ab44ea79dc3c92b0058f6c3079960dd965e762b0c6daa3d4fc8ed8bf",
    "_core.py": "68759af57d2b5b67ee6b2a3341c2db358b859052d0e22635beae549e2ab7c1dc",
    # The vendor engine set (4 files) is pinned inside _pins.py; the battery
    # consumes that table from the sha-verified _pins.py itself (see
    # _vendor_pins_from_package), so there is one vendor-pin surface only.
    "_pins.py": "6049d82d607d4f3bd7ad41fb304799299893404b39f09f099af78957540c18ca",
    "fraction_oracle.py": "ce0cd79078c74f1841e4ebf3416fe27297b067ace6132e0d5c1aadf525c6b3e0",
    "oracle_k2disp.py": "860f737a41cd2151d6efd440ff54dbc182044f40a0e4bd87e0d970d616dad766",
    "receipt.py": "f9ead65b5b934607f25abd1ca74a2431c7807d04e38989ef09308032f6e251e9",
    "strata.py": "164187bbaafb1ac9aefebe7246b1b30239a3a43e8ac969020b39924a77ba102a",
    "witness_bridge.py": "fc2554756764fb4f5da5256ed601c331dd5befb1d4441f72afdbf67917418f37",
}
# vendored bytes: SINGLE-SOURCED from trust/_pins.py, the package's own pin
# table. The battery formerly carried a second, hand-copied vendor pin dict
# here; two surfaces can drift (a vendor re-pin moves _pins.py while the
# battery copy lags), and T0 then fails tamper-shaped on pristine shipped
# bytes. One pin surface: T0 verifies _pins.py's own
# bytes against the battery-side PACKAGE_PINS sha above, then ast-parses its
# PINS literal (nothing executed — the pre-import property holds) and walks
# the vendor tree against THAT table. Fail-closed: an unverified or
# unparseable pin source returns None and T0_vendor_pins_battery_side FAILS
# by name — a pin table the battery cannot authenticate is never consumed.


def _vendor_pins_from_package():
    """PINS from the sha-verified trust/_pins.py, or None (refuse)."""
    import ast
    src_path = os.path.join(PKG, "_pins.py")
    try:
        if sha(src_path) != PACKAGE_PINS["_pins.py"]:
            return None
        tree = ast.parse(open(src_path, "rb").read())
    except (OSError, SyntaxError, ValueError):
        return None
    for node in tree.body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == "PINS"):
            try:
                pins = ast.literal_eval(node.value)
            except ValueError:
                return None
            if (isinstance(pins, dict) and pins
                    and all(isinstance(k, str) and isinstance(v, str)
                            and len(v) == 64 and set(v) <= set("0123456789abcdef")
                            for k, v in pins.items())):
                return pins
            return None
    return None
# The source control receipts (CONTROL*.json, sec63_final.json) are reference
# stores that do not ship, so they are not pinned here. T3 grades the engines'
# built-in control stages self-contained
# (design facts + same-run pristine differentials) instead.
# joint-fixture reference receipts (T4 bars) + the kira truth table.
# TRUST_RUNS_ROOT points at a reference run tree (pilotA/pilotB/
# pilotC/fixture layout); TRUST_KIRA_TABLE and TRUST_ART_DIR point at the
# kira truth table and the SYSTEM-dump dir. These reference receipts are not
# included in the release; absent files fail loudly (env-configured).
_TRUN = os.environ.get("TRUST_RUNS_ROOT", "")
RECEIPT_PINS = {
    "PILOT_A": (os.path.join(_TRUN, "pilotA/PILOT_A.json"),
                "bbac3a818b0f4a43525ea1c9d726b163d0fb4263ae0a13047cb9b27a900f4e21"),
    "PILOT_B": (os.path.join(_TRUN, "pilotB/PILOT_B_legs12.json"),
                "7252d471613fdbb8e4669ee59e79db9d2e0f5846497e19db9a4af29295b151a1"),
    "PILOT_C": (os.path.join(_TRUN, "pilotC/PILOT_C.json"),
                "e70309e32b96fa521692d5113a0c6639b61c56a8d21097cca8b52c2fd206510e"),
    "FULL_pt1": (os.path.join(_TRUN, "fixture/FULL_1to7.json"),
                 "28479e3f5c0bb80df6eb8b9b5408bd73684bdd071bcbfa27a2cda9089e1a37e9"),
    "FULL_pt2": (os.path.join(_TRUN, "fixture/FULL_1to7_pt2.json"),
                 "d005a586f6803e52e1f2c80515d10f98d75c691f93e2b14f3a55a969ffbf7b9b"),
    "P1": (os.path.join(_TRUN, "fixture/P1_2567.json"),
           "956239535e52da72699aa7a408665e7489159f80f70c488f52994915a1499286"),
    "FIXTURE_DECISION": (os.path.join(_TRUN, "FIXTURE_DECISION.md"),
                         # (23-min composition stated; 30/30 cite rows)
                         "2373a3f5c87879ef272d6d83fc7e25929764645178aaffcb5ffbb5153c88f118"),
    "ADAPTER_PROBE": (os.path.join(_TRUN, "ADAPTER_PROBE.md"),
                      "88c57c96bbdc07c5023427ae6ec899e073fd55772f613165db19110c3c598e14"),
    "EMIT_REPORT": (os.path.join(_TRUN, "pilotB/emit/EMIT_REPORT.json"),
                    "1eeb40aafa87618b04ee9c657887446e657c5de75d586d0103a9a694d200a1b7"),
    "KIRA_TABLE": (os.environ.get("TRUST_KIRA_TABLE", ""),
                   "5f75ee981569e780d81f7e60d31fed800fc1a59a33dfba0b5b3d519194eef8bc"),
}
EMIT_DIR = os.path.join(_TRUN, "pilotB/emit")
EMIT_WITNESS_ROLLUP = ("3f2f8f02f36b22a007924219275256b1f5d02faf7f1ce18b7f126a430f26a83e", 84)

FAM = "lbl3m2L_k2disp"
PRIMES = [2147483647, 2147483629]
SLICE_D0, SLICE_E0 = 1234577, 87654321
ART = os.environ.get("TRUST_ART_DIR", "")

V10_TOP_KEYS = {"receipt_version", "kind", "identity", "p", "point", "family",
                "target", "c", "labels", "lam", "system"}
V10_REQUIRED = {"receipt_version", "kind", "identity", "p", "target", "c", "lam"}
V10_IDENTITY = "sum_i lam[i]*R_i == e_target - sum_m c[m]*e_m (mod p)"

FAILED = []
RESULTS = {"battery": "TRUST", "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
           "scratch": OUT, "legs": {}, "walls_s": {}}
MUTATIONS = []
T0_START = time.time()


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check(name, cond, extra=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name} {extra}", flush=True)
    RESULTS["legs"][name] = {"pass": bool(cond), "extra": str(extra)[:300]}
    if not cond:
        FAILED.append(name)


def leg_wall(leg, t0):
    RESULTS["walls_s"][leg] = round(time.time() - t0, 1)


def save():
    json.dump(RESULTS, open(os.path.join(OUT, "BATTERY_REPORT.json"), "w"),
              indent=1, default=str)
    json.dump(MUTATIONS, open(os.path.join(OUT, "MUTATION_RECEIPTS.json"), "w"),
              indent=1, default=str)


def mutation(mid, leg, target, mut, fixture, wall, caught, note=""):
    MUTATIONS.append({"id": mid, "leg": leg, "target": target, "mutation": mut,
                      "fixture": fixture, "wall_s": round(wall, 1),
                      "caught": bool(caught), "note": note})
    print(f"[{'PASS' if caught else 'FAIL'}] {mid} (mutation control) "
          f"caught={caught} {note[:120]}", flush=True)
    if not caught:
        FAILED.append(mid)


# child bootstrap (mirrors trust._core._CHILD_BOOT): under -s the
# user-writable user site is never site-processed — its .pth files (arbitrary
# code at interpreter startup, the site-reorder vector) never run.
# Children that need the third-party substrate (sympy/python-flint, when
# they live in the user site) get that ONE directory appended at the TAIL of
# sys.path — stdlib and system site keep precedence, the append executes
# nothing.
CHILD_BOOT = (
    "import sys, site, runpy\n"
    "p = sys.argv.pop(1)\n"
    "sys.argv[0] = p\n"
    "u = site.getusersitepackages()\n"
    "if u and u not in sys.path:\n"
    "    sys.path.append(u)\n"
    "runpy.run_path(p, run_name='__main__')\n"
)


def child(args, timeout=120, env_extra=None):
    """Child launcher: PYTHON* env scrubbed  and,
    for python children, -E -P -s -B forced in code so a shadow dir on an
    inherited PYTHONPATH can never reach a child (-E), the script-dir/cwd
    sys.path prepend is off (-P — forced in code, no longer
    call-site discipline; run scratch holds battery-written files), the
    user site is never site-processed (-s — user-site .pth
    files never execute; sympy/flint ride CHILD_BOOT's inert tail append)
    and no bytecode is written (-B; with -E the env var form is ignored).
    Call sites passing -I or -S keep the plain direct launch — those are
    the pure-stdlib children (T6) whose whole point is NO site dirs at all."""
    env = {k: v for k, v in os.environ.items()
           if not (k.startswith("PYTHON") and k != "PYTHONDONTWRITEBYTECODE")}
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if env_extra:
        env.update(env_extra)
    if args and args[0] == sys.executable:
        if "-I" in args or "-S" in args:
            flags = [f for f in ("-E", "-P", "-B") if f not in args]
            args = [args[0]] + flags + list(args[1:])
        else:
            rest = [a for a in args[1:] if a not in ("-E", "-P", "-B", "-s")]
            args = [args[0], "-E", "-P", "-s", "-B", "-c", CHILD_BOOT] + rest
    return subprocess.run(args, capture_output=True, text=True,
                          timeout=timeout, env=env, cwd=OUT)


# =========================== T0: package pins (PRE-IMPORT) ==================
t0 = time.time()
print("== T0 package pins (pre-import, battery-side) ==", flush=True)
# purge regenerable bytecode caches in the package + battery dirs
for d in (PKG, HERE):
    pc = os.path.join(d, "__pycache__")
    if os.path.isdir(pc):
        shutil.rmtree(pc, ignore_errors=True)

# trust/ directory: file-set-exact vs PACKAGE_PINS — nothing else, no dirs,
# no symlinks, shas match
ok_set = True
details = []
entries = sorted(os.listdir(PKG))
if set(entries) != set(PACKAGE_PINS):
    ok_set = False
    details.append(f"file-set {sorted(set(entries) ^ set(PACKAGE_PINS))}")
for e in entries:
    p = os.path.join(PKG, e)
    if os.path.islink(p) or not os.path.isfile(p):
        ok_set = False
        details.append(f"{e}: not a regular file")
    elif e in PACKAGE_PINS and sha(p) != PACKAGE_PINS[e]:
        ok_set = False
        details.append(f"{e}: sha mismatch")
check("T0_package_fileset_exact_and_pinned", ok_set, "; ".join(details)[:200]
      or f"{len(PACKAGE_PINS)}/8 pinned, nothing else")

# TOOL root: exactly GUIDE.md + MANUAL.md + selftest.sh + the member/package
# dirs (battery, receipt, strata, trust, vendor); nothing importable strays
root_files = sorted(os.listdir(TOOL))
# GUIDE.md is part of the pinned fileset
ok_root = (set(root_files) == {"GUIDE.md", "MANUAL.md", "battery", "receipt",
                               "selftest.sh", "strata", "trust", "vendor"}
           and not any(os.path.islink(os.path.join(TOOL, e)) for e in root_files))
check("T0_tool_root_fileset", ok_root, str(root_files))
bat_files = sorted(e for e in os.listdir(HERE) if e != "__pycache__")
check("T0_battery_dir_fileset", bat_files == ["battery.py"], str(bat_files))

# vendor tree: walked battery-side against the single-sourced pin table
# (PINS ast-parsed from the sha-verified trust/_pins.py; None = refuse)
VENDOR_PINS = _vendor_pins_from_package()
if VENDOR_PINS is None:
    check("T0_vendor_pins_battery_side", False,
          "REFUSED: vendor pin source trust/_pins.py failed its battery-side "
          "sha gate or its PINS table did not parse — an unauthenticated pin "
          "table is never consumed")
else:
    ok_v = True
    vdetails = []
    seen = set()
    for root, dirs, files in os.walk(VENDOR):
        for d in list(dirs):
            if os.path.islink(os.path.join(root, d)):
                ok_v = False
                vdetails.append(f"symlinked dir {d}")
        for f in files:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, VENDOR)
            if "__pycache__" in rel:
                ok_v = False
                vdetails.append(f"pycache content {rel}")
                continue
            seen.add(rel)
            if os.path.islink(p):
                ok_v = False
                vdetails.append(f"symlink {rel}")
            elif rel not in VENDOR_PINS:
                ok_v = False
                vdetails.append(f"unpinned {rel}")
            elif sha(p) != VENDOR_PINS[rel]:
                ok_v = False
                vdetails.append(f"sha mismatch {rel}")
    if seen != set(VENDOR_PINS):
        ok_v = False
        vdetails.append(f"missing {sorted(set(VENDOR_PINS) - seen)}")
    check("T0_vendor_pins_battery_side", ok_v, "; ".join(vdetails)[:200]
          or f"{len(VENDOR_PINS)}/{len(VENDOR_PINS)} vs the _pins.py table")
leg_wall("T0", t0)
save()
if FAILED:
    print(f"\nOVERALL: FAIL at T0 — refusing to import the package: {FAILED}")
    sys.exit(1)

sys.path.insert(0, TOOL)
import trust  # noqa: E402
from trust import oracle_k2disp as okd  # noqa: E402
from trust import witness_bridge as wb  # noqa: E402
from trust import fraction_oracle as fo  # noqa: E402
from trust import receipt as trr  # noqa: E402
from fractions import Fraction as F  # noqa: E402

os.environ["TRUST_OUT_ROOT"] = os.path.join(OUT, "out")

# =========================== T1: vendor pins ================================
t0 = time.time()
print("== T1 vendor pins + LINKED ==", flush=True)
check("T1_pins_crosscheck_battery_vs_package",
      dict(trust.PINS) == VENDOR_PINS,
      "imported trust.PINS == the T0 ast-parsed table from the sha-verified "
      "_pins.py (single-sourced; an import resolving a different pin table "
      "than the verified file refuses here)")
rep = trust.verify(quiet=True)
check("T1_verify_vendor_ok", len(rep["vendor_ok"]) == 4, f"{len(rep['vendor_ok'])}/4")
check("T1_linked_homes_present", not rep["linked_missing"],
      str(rep["linked_missing"])[:150])
# LINKED_SHAS advisory baselines exist for every LINKED home and
# verify() actually compares (linked_drift key). Drift itself is ADVISORY by
# design (live tools advance independently) — reported here, never a FAIL.
check("T1_linked_shas_baselined",
      set(trust.LINKED_SHAS) == set(trust.LINKED)
      and all(trust.LINKED_SHAS[k] for k in trust.LINKED_SHAS)
      and "linked_drift" in rep,
      f"{sum(len(v) for v in trust.LINKED_SHAS.values())} baselines across "
      f"{sorted(trust.LINKED_SHAS)}; current drift: "
      f"{rep['linked_drift'] or 'none'}")
if rep["linked_drift"]:
    print(f"[INFO] LINKED drift (advisory, not a FAIL): {rep['linked_drift']}",
          flush=True)
leg_wall("T1", t0)
save()

# =========================== T2: tamper control (scratch copies) ============
t0 = time.time()
print("== T2 tamper control (scratch copies, child-mode) ==", flush=True)
CLONES = os.path.join(OUT, "t2_clone")
shutil.rmtree(CLONES, ignore_errors=True)
VERIFY_CHILD = os.path.join(OUT, "t2_verify_child.py")
open(VERIFY_CHILD, "w").write("""\
import json, sys
sys.path.insert(0, sys.argv[1])
import trust
try:
    rep = trust.verify(quiet=True)
    print(json.dumps({"raised": None,
                      "n_drift": len(rep.get("linked_drift", []))}))
except Exception as e:
    print(json.dumps({"raised": type(e).__name__, "msg": str(e)[:200]}))
""")


def make_clone(tag):
    dst = os.path.join(CLONES, tag)
    os.makedirs(dst)
    shutil.copytree(PKG, os.path.join(dst, "trust"),
                    ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(VENDOR, os.path.join(dst, "vendor"),
                    ignore=shutil.ignore_patterns("__pycache__"))
    return dst


def clone_verify(tag):
    r = child([sys.executable, "-P", VERIFY_CHILD, os.path.join(CLONES, tag)])
    try:
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}


c = make_clone("pristine")
v = clone_verify("pristine")
check("T2_pristine_clone_verifies", v.get("raised") is None, str(v)[:120])

c = make_clone("byteflip")
with open(os.path.join(c, "vendor/lpsyz/lp_syz.py"), "ab") as f:
    f.write(b"\n# tampered\n")
v = clone_verify("byteflip")
check("T2_byteflip_refused", v.get("raised") == "VendorTamperError", str(v)[:120])

c = make_clone("unpinned")
open(os.path.join(c, "vendor/lpsyz/evil.py"), "w").write("x = 1\n")
v = clone_verify("unpinned")
check("T2_unpinned_extra_refused", v.get("raised") == "VendorTamperError", str(v)[:120])

c = make_clone("symlink")
tgt = os.path.join(c, "vendor/lpsyz/lp_syz_431.py")
side = os.path.join(c, "identical_bytes.py")
shutil.move(tgt, side)
os.symlink(side, tgt)  # byte-identical content, still must refuse
v = clone_verify("symlink")
check("T2_symlink_swap_refused_despite_identical_bytes",
      v.get("raised") == "VendorTamperError", str(v)[:120])

c = make_clone("pycache")
pc = os.path.join(c, "vendor/lpsyz/__pycache__")
os.makedirs(pc)
open(os.path.join(pc, "lp_syz.cpython-312.pyc"), "wb").write(b"\x00fake")
v = clone_verify("pycache")
check("T2_pycache_content_refused", v.get("raised") == "VendorTamperError", str(v)[:120])

# mutation control: a flipped LINKED_SHAS baseline in a scratch copy must
# be REPORTED by verify() (linked_drift non-empty) WITHOUT raising — advisory
# loud, never a brick (the promise the three docstrings make, now implemented)
c = make_clone("linkdrift")
_pp = os.path.join(c, "trust", "_pins.py")
_ptxt = open(_pp).read()
_base = trust.LINKED_SHAS["receipt"]["core.py"]
assert _base in _ptxt, "LINKED_SHAS baseline anchor missing in _pins.py"
open(_pp, "w").write(_ptxt.replace(_base, "0" * 64))
v = clone_verify("linkdrift")
check("T2_linked_drift_advisory_loud_not_fatal",
      v.get("raised") is None and (v.get("n_drift") or 0) >= 1,
      f"{v.get('n_drift')} drift entr(ies) reported, verify() did NOT raise")

# launch-path refusal on the byteflipped copy (valid env, tampered bytes)
LAUNCH_CHILD = os.path.join(OUT, "t2_launch_child.py")
open(LAUNCH_CHILD, "w").write("""\
import json, os, sys
sys.path.insert(0, sys.argv[1])
os.environ["TRUST_OUT_ROOT"] = sys.argv[2]
import trust
try:
    trust.run_vendored("lp_syz_431", ["--stage", "control452", "--D0", "6",
                                      "--out", "x.json"], timeout=60)
    print(json.dumps({"raised": None}))
except Exception as e:
    print(json.dumps({"raised": type(e).__name__, "msg": str(e)[:200]}))
""")
r = child([sys.executable, "-P", LAUNCH_CHILD, os.path.join(CLONES, "byteflip"),
           os.path.join(OUT, "t2_launch_out")])
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}
check("T2_launch_refused_on_tampered_clone",
      v.get("raised") == "VendorTamperError", str(v)[:120])
leg_wall("T2", t0)
save()

# =========================== T3: lp-syz restore-verify ======================
t0 = time.time()
print("== T3 restore-verify (shipped control stages, self-graded) ==", flush=True)
# The engines' BUILT-IN control stages, graded battery-side with NO external
# receipts: ctrl45 is the planted-violation positive control (a known rank
# deficiency among the 5 forced masters — the run must find EXACTLY ONE
# violation, named over master columns only) and ctrl452 the clean-port
# negative control (0 violations by design). Structural counts follow the
# Singular syzygy basis, so they are recorded, never pinned; the mutation
# controls grade mutants against the SAME-RUN pristine receipts. The former
# value-level legs (sec63 details, source control receipts) consumed
# reference stores that do not ship and are NAMED SKIPS.


def skip(name, why):
    print(f"[SKIP] {name} — {why}", flush=True)
    RESULTS["legs"][name] = {"pass": None, "skip": why}


def grade_tadpole(g):
    return (g.get("sigma") == -1 and g.get("match_sigma_-1") is True
            and g.get("n_violations") == 0 and g.get("n_free_nonmaster") == 0
            and g.get("n_columns") == 11 and g.get("n_relations") == 60)


def grade_452(g):
    # clean-port negative control: 0 violations BY DESIGN (engine's own
    # stage comment + control_pass verdict), and the solve really happened
    return (g.get("n_violations") == 0 and not g.get("viol_rows_named")
            and g.get("control_pass") is True
            and (g.get("n_relations") or 0) > 0 and (g.get("n_columns") or 0) > 0)


def grade_45(g):
    # planted-violation positive control: EXACTLY ONE violation, its named
    # row supported on master columns only, no free non-masters
    vr = g.get("viol_rows_named")
    return (g.get("n_violations") == 1 and g.get("n_free_nonmaster") == 0
            and g.get("control_pass") is False
            and bool(vr) and len(vr) == 1 and len(vr[0]) >= 2
            and all(k.startswith("M(") for k in vr[0]))


_T3_STRUCT = ("n_relations", "n_columns", "n_violations", "n_free_nonmaster")


def _t3_struct(g):
    return ", ".join(f"{k}={g.get(k)}" for k in _T3_STRUCT)


T3_RUNS = [
    ("T3_tadpole_control", "lp_syz",
     ["--stage", "control", "--D0", "4", "--out", "T3_tadpole.json"], grade_tadpole),
    ("T3_ctrl452_negative", "lp_syz_431",
     ["--stage", "control452", "--D0", "6", "--out", "T3_ctrl452.json"], grade_452),
    ("T3_ctrl45_positive", "lp_syz_431",
     ["--stage", "control", "--D0", "6", "--out", "T3_ctrl45.json"], grade_45),
]
skip("T3_sec63_tower", "value-level reference receipt is not included in "
     "the release; leg skipped (the ctrl45/"
     "ctrl452 planted controls above are the shipped-scope legs)")
t3_good = {}
for name, eng, args, grade in T3_RUNS:
    if eng == "lp_syz" and not os.environ.get("TRUST_LPSYZ_KTAB"):
        skip(name, "TRUST_LPSYZ_KTAB unset — the lp_syz truth table is "
             "a reference file not included; set the var to run the "
             "value-level tadpole control")
        continue
    # stale-readback guard: unlink the expected output BEFORE the
    # child spawns — a reused scratch must never grade last run's bytes
    outp = os.path.join(os.environ["TRUST_OUT_ROOT"], dict(zip(args, args[1:]))["--out"])
    try:
        os.unlink(outp)
    except FileNotFoundError:
        pass
    t1 = time.time()
    r = trust.run_vendored(eng, args, timeout=900)
    wall = round(time.time() - t1, 1)
    if r.returncode != 0 or not os.path.isfile(outp):
        check(name, False, f"rc={r.returncode} {r.stderr[-150:]}")
        continue
    g = json.load(open(outp))
    t3_good[name] = g
    check(name, grade(g), f"{wall} s, self-graded; {_t3_struct(g)} "
          f"(counts recorded, basis-dependent — not pinned)")
# wall floor (stub detector, re-tuned when the sec63 leg retired): the two
# ported control solves cannot be sub-1 s with a real Singular; the
# mutation controls below are the primary anti-stub evidence
check("T3_wall_floor", (time.time() - t0) >= 1.0,
      f"{round(time.time() - t0, 1)} s >= 1 s")

# --- T3 mutation controls (the planted lp-syz traps, walls measured) --------
MUT = os.path.join(OUT, "t3_mut")
os.makedirs(MUT, exist_ok=True)
src_q = open(trust.vendored_path("lp_syz_431")).read()
# the port carries the identical 3-line descent block, so the descent-sign
# trap mutates the SAME engine the shipped controls exercise
OLD_DESC = ("                            key = self.canon(Sc, e)\n"
            "                            if key not in self.colid: ok = False; break\n"
            "                            row[self.colid[key]] = row.get(self.colid[key], Fraction(0)) + c")
assert src_q.count(OLD_DESC) == 1, "descent-block anchor not unique in pinned lp_syz_431.py"
open(os.path.join(MUT, "mut_descent.py"), "w").write(
    src_q.replace(OLD_DESC, OLD_DESC.replace("Fraction(0)) + c", "Fraction(0)) - c")))
OLD_AUT = "    def _aut_min(self, R, m):\n        best = tuple(m)"
assert src_q.count(OLD_AUT) == 1, "_aut_min anchor not unique in pinned lp_syz_431.py"
open(os.path.join(MUT, "mut_autskip.py"), "w").write(
    src_q.replace(OLD_AUT, "    def _aut_min(self, R, m):\n"
                  "        return tuple(m)  # MUTATION: Aut(G) quotient skipped\n"
                  "        best = tuple(m)"))


def run_mutant(script, args, timeout=600):
    outp = os.path.join(MUT, script + ".out.json")
    # stale-readback guard: pre-delete so a crashed child can never
    # be graded against the PREVIOUS run's file
    try:
        os.unlink(outp)
    except FileNotFoundError:
        pass
    t1 = time.time()
    r = child([sys.executable, "-P", os.path.join(MUT, script)] + args +
              ["--out", outp], timeout=timeout)
    wall = time.time() - t1
    try:
        g = json.load(open(outp))
    except Exception:
        g = {"_err": r.stderr[-150:]}
    # child rc recorded in every mutation receipt: disambiguates caught-by-
    # grade from caught-by-crash 
    g["_child_rc"] = r.returncode
    return g, wall


_p45 = t3_good.get("T3_ctrl45_positive")
_p452 = t3_good.get("T3_ctrl452_negative")


def _differs(g, p):
    """Mutant vs SAME-RUN pristine receipt: any graded structural field or
    the named violation rows differ. Basis-controlled by construction (same
    box, same Singular, minutes apart) — no pinned counts needed."""
    if not p:
        return False
    return (any(g.get(k) != p.get(k) for k in _T3_STRUCT)
            or g.get("viol_rows_named") != p.get("viol_rows_named"))


# m1: the descent-sign trap's wired catcher was the sec63 value-level
# details grade against the internal reference receipt — retired with it.
skip("T3m1_descent_sign_valuelevel_catcher",
     "catcher was the sec63 value-level details grade vs a "
     "reference receipt (not included); the shipped controls are "
     "measured-blind to this mutation — the T3m1b demo records it live")
# m1b: measured blindspot demo (informational, NOT a battery check): the
# descent-sign flip leaves the shipped control receipts FIELD-IDENTICAL to
# the same-run pristine ones — a support/structure grade cannot see it;
# only a value-level grade can (the documented footgun)
g, w = run_mutant("mut_descent.py", ["--stage", "control452", "--D0", "6"])
_same452 = _p452 is not None and not _differs(g, _p452)
MUTATIONS.append({"id": "T3m1b_descent_sign_ctrl452_blindspot", "leg": "T3",
                  "target": "lp_syz_431.py (port descent block)",
                  "mutation": "descent sign flip",
                  "fixture": "ctrl452 vs same-run pristine",
                  "wall_s": round(w, 1), "caught": not _same452,
                  "expected_blindspot": True,
                  "note": "MEASURED: mutant ctrl452 receipt field-identical "
                          "to pristine (struct + viol rows) — documents the "
                          "footgun; a value-level grade is the only catcher "
                          f"[field_identical={_same452}]"})
print("[INFO] T3m1b descent-sign ctrl452 blindspot demo recorded "
      f"(field-identical={_same452}; caught=False EXPECTED, {round(w, 1)} s)",
      flush=True)
# m2: Aut(G)-quotient skip -> the ctrl45 design grade must FAIL (the planted
# violation disappears) and ctrl452 must DIFFER from same-run pristine
g, w = run_mutant("mut_autskip.py", ["--stage", "control", "--D0", "6"])
mutation("T3m2_autskip_caught_by_ctrl45", "T3", "lp_syz_431.py",
         "_aut_min returns m unquotiented",
         "ctrl45 design grade + pristine diff", w,
         (not grade_45(g)) or _differs(g, _p45),
         f"violations {g.get('n_violations')} vs 1 by design; cols "
         f"{g.get('n_columns')} vs pristine {(_p45 or {}).get('n_columns')} "
         f"[child_rc={g.get('_child_rc')}]")
g, w = run_mutant("mut_autskip.py", ["--stage", "control452", "--D0", "6"])
mutation("T3m2b_autskip_caught_by_ctrl452", "T3", "lp_syz_431.py",
         "_aut_min returns m unquotiented", "ctrl452 vs same-run pristine", w,
         (not grade_452(g)) or _differs(g, _p452),
         f"cols {g.get('n_columns')} vs pristine "
         f"{(_p452 or {}).get('n_columns')} [child_rc={g.get('_child_rc')}]")
# m3: grader teeth — corrupted result dicts must not grade PASS
if _p45:
    bad = dict(_p45)
    bad["n_violations"] = 0
    mutation("T3m3_grader_teeth_ctrl45", "T3", "result dict (scratch copy)",
             "n_violations 1->0", "grade_45", 0.0, not grade_45(bad))
if "T3_tadpole_control" in t3_good:
    bad = dict(t3_good["T3_tadpole_control"])
    bad["n_violations"] = 1
    mutation("T3m3b_grader_teeth_tadpole", "T3", "result dict (scratch copy)",
             "n_violations 0->1", "grade_tadpole", 0.0, not grade_tadpole(bad))
leg_wall("T3", t0)
save()

# =========================== T4: joint three-check fixture ==================
t0 = time.time()
print("== T4 joint three-check fixture (bars pinned + live bounded) ==", flush=True)
recs = {}
ok_pins = True
for k, (path, want) in RECEIPT_PINS.items():
    if not os.path.isfile(path) or sha(path) != want:
        ok_pins = False
        print(f"    receipt pin MISMATCH: {k}", flush=True)
    elif path.endswith(".json"):
        recs[k] = json.load(open(path))
check("T4_receipts_sha_pinned", ok_pins, f"{len(RECEIPT_PINS)} receipts/inputs")
import glob  # noqa: E402

wfiles = sorted(glob.glob(os.path.join(EMIT_DIR, "w_*.json")))
roll = hashlib.sha256("\n".join(
    os.path.basename(w) + " " + sha(w) for w in wfiles).encode()).hexdigest()
check("T4_emit_witness_rollup_pinned",
      (roll, len(wfiles)) == EMIT_WITNESS_ROLLUP, f"{len(wfiles)} witnesses")


def bars_pilot_b(g):
    return (g.get("B1_verdict") == "PASS" and g.get("B1_certified") == "84/84"
            and g.get("B1_row_statuses") == {"CERTIFIED": 84}
            and g.get("B1_n_anchored") == 42
            and g.get("B2_cores_byte_identical") is True
            and g.get("B2", {}).get("n") == 84
            and g.get("B2", {}).get("pass_a") == 84
            and g.get("B2", {}).get("pass_b") == 84
            and g.get("B2", {}).get("tamper_must_fail") is True)


def bars_full(g, d0, e0):
    det = g.get("details") or []
    return (g.get("n_targets") == 42 and g.get("n_match") == 42
            and g.get("sigma") == -1 and g.get("n_violations") == 0
            and g.get("d0") == d0 and g.get("eta0") == e0
            and len(det) == 42 and all(d.get("match") is True for d in det))


if ok_pins:
    check("T4_bars_pilotB_emit_84of84", bars_pilot_b(recs["PILOT_B"]),
          "84/84 CERTIFIED, 0 CLAIM_MISMATCH, 2 primes; dual-core 84+84")
    check("T4_bars_full42_pt1", bars_full(recs["FULL_pt1"], "97/23", "5/7"),
          "42/42 EXACT, sigma=-1, 0 viol (269.2 s reference run)")
    check("T4_bars_full42_pt2", bars_full(recs["FULL_pt2"], "61/17", "3/11"),
          "42/42 EXACT at 2nd point (369.5 s reference run)")
    check("T4_bars_pilotA", recs["PILOT_A"].get("verdict") == "4/4 legs REPRODUCED")
    pc = recs["PILOT_C"]
    check("T4_bars_pilotC", pc.get("verdict") == "PASS"
          and pc.get("verify", {}).get("n") == 14
          and pc.get("verify", {}).get("pass_tools_core") == 14
          and pc.get("verify", {}).get("pass_winnow_core") == 14
          and all(pc.get("planted_faults", {}).values()))

# --- live L1: check-3 bounded — P1 {2,5,6,7} tower vs kira (4.0 s measured)
TABLE = RECEIPT_PINS["KIRA_TABLE"][0]
ents = okd.parse_kira_table(TABLE)
p1_targets = []
for tt in sorted(ents):
    u = set(okd.support(tt))
    for k in ents[tt]:
        u.update(okd.support(k))
    if u <= {2, 5, 6, 7}:
        p1_targets.append(tt)
check("T4_live_P1_target_census", len(p1_targets) == 7, f"{len(p1_targets)} targets")
t1 = time.time()
r1 = okd.reduce_and_compare(TABLE, F(97, 23), F(5, 7), targets=p1_targets,
                            verbose=False, engine_tag="batteryA")
p1 = recs.get("P1", {})
check("T4_live_P1_7of7_exact_vs_kira",
      r1["sigma"] == -1 and r1["n_match"] == 7 and r1["n_violations"] == 0
      and r1["n_relations"] == p1.get("n_relations")
      and r1["n_columns"] == p1.get("n_columns"),
      f"{round(time.time() - t1, 1)} s; 437 rels/182 cols == pinned P1 receipt")
tw = r1["_tower"]

# --- live L2: check-2 bounded — fresh adapter re-parse, 84 witnesses through
# BOTH receipt cores (0.6 s measured)
core_a = trr.core()
core_b = trr.winnow_core()
okbi, _shas = trr.byte_identity()
check("T4_cores_byte_identical", okbi)
RECEIPT_DIR = trust.LINKED["receipt"]  # the receipt member (<tool>/receipt)
sys.path.insert(0, RECEIPT_DIR)
from adapters import strata as ad  # noqa: E402

sys.path.remove(RECEIPT_DIR)
t1 = time.time()
n = pa = pb = 0
rows_by_p = {}
for p in PRIMES:
    sysd = ad.load_rows(ART, FAM, p, SLICE_D0, SLICE_E0)
    rows = sysd["rows"]
    rows_by_p[p] = rows
    fp = core_a.system_fingerprint(rows, p)
    for q in sorted(glob.glob(os.path.join(EMIT_DIR, f"w_{FAM}_{p}_*.json"))):
        wa = core_a.load_witness(q)
        oka, _ = core_a.verify_row(rows, wa, check_fingerprint=fp)
        wbv = core_b.load_witness(q)
        okb, _ = core_b.verify_row(rows, wbv, check_fingerprint=fp)
        n += 1
        pa += bool(oka)
        pb += bool(okb)
check("T4_live_reverify_84_dual_core", n == 84 and pa == 84 and pb == 84,
      f"{pa}+{pb}/84+84 in {round(time.time() - t1, 1)} s (own re-parse, fp-checked)")

# --- live L3 : INDEPENDENT stdlib re-parse spot-check ----
# The adapter rows above come from strata's OWN loader/CoeffEvaluator (the
# named STRATA-LOADER SUBSTRATE caveat, trust/receipt.py docstring).
# Re-derive a sha-seeded FOREIGN subsample of rows DIRECTLY from the
# SYSTEM_*.gz dumps with battery-side code only: gzip text parse + ast-walk
# modular arithmetic. NO loader import, NO eval(), NEVER floats (a
# CoeffEvaluator float-division bug class is exactly what this must not share).
import ast  # noqa: E402
import gzip  # noqa: E402
import re  # noqa: E402

_COEFF_OK = re.compile(r"^[0-9deta+\-*/^() ]+$")


def _indep_lit_int(nd):
    if isinstance(nd, ast.Constant) and isinstance(nd.value, int):
        return nd.value
    if isinstance(nd, ast.UnaryOp) and isinstance(nd.op, ast.USub):
        return -_indep_lit_int(nd.operand)
    raise ValueError("non-literal exponent")


def _indep_eval_modp(expr, p, d0, e0):
    if not _COEFF_OK.match(expr):
        raise ValueError(f"unexpected coeff chars: {expr!r}")
    env = {"d": d0 % p, "eta": e0 % p}

    def ev(nd):
        if isinstance(nd, ast.Expression):
            return ev(nd.body)
        if isinstance(nd, ast.Constant) and isinstance(nd.value, int):
            return nd.value % p
        if isinstance(nd, ast.Name):
            return env[nd.id]
        if isinstance(nd, ast.UnaryOp):
            v = ev(nd.operand)
            return (-v) % p if isinstance(nd.op, ast.USub) else v
        if isinstance(nd, ast.BinOp):
            if isinstance(nd.op, ast.Pow):
                return pow(ev(nd.left), _indep_lit_int(nd.right), p)
            a, b = ev(nd.left), ev(nd.right)
            if isinstance(nd.op, ast.Add):
                return (a + b) % p
            if isinstance(nd.op, ast.Sub):
                return (a - b) % p
            if isinstance(nd.op, ast.Mult):
                return (a * b) % p
            if isinstance(nd.op, ast.Div):
                if b == 0:
                    raise ZeroDivisionError(f"denominator hit mod {p}: {expr!r}")
                return (a * pow(b, p - 2, p)) % p
        raise ValueError(f"unsupported node in {expr!r}")

    return ev(ast.parse(expr.replace("^", "**"), mode="eval"))


def _indep_seed_indices(n, k, material):
    """battery-side sha-chain FOREIGN sampler (deliberately NOT package code)."""
    if k >= n:
        return list(range(n))
    out, seen = [], set()
    hh = hashlib.sha256(material.encode()).digest()
    while len(out) < k:
        for i in range(0, len(hh) - 3, 4):
            v = int.from_bytes(hh[i:i + 4], "big") % n
            if v not in seen:
                seen.add(v)
                out.append(v)
                if len(out) == k:
                    break
        hh = hashlib.sha256(hh).digest()
    return sorted(out)


def _indep_structural(art_dir, family):
    """SYSTEM_*.gz -> rows as [(coeff_str, col)], trivial-sector terms
    dropped, NO evaluation (bounded design point: only the sampled rows get
    evaluated)."""
    tpath = os.path.join(art_dir, "sectormappings", family, "trivialsector")
    triv = set()
    if os.path.exists(tpath):
        triv = {int(t) for t in re.findall(r"\d+", open(tpath).read())}
    sfiles = sorted(glob.glob(os.path.join(art_dir, "tmp", family,
                                           "SYSTEM_*.gz")))
    rows = []
    for f in sfiles:
        lines = gzip.open(f, "rt").read().split("\n")
        i, n = 0, len(lines)
        while i < n:
            if lines[i] != "Eq":
                i += 1
                continue
            nt = int(lines[i + 2])
            row = []
            for j in range(nt):
                parts = lines[i + 3 + j].split()
                if len(parts) != 6:
                    raise ValueError(f"bad term line {lines[i + 3 + j]!r}")
                if int(parts[3]) in triv:
                    continue
                row.append((parts[1], int(parts[2])))
            i += 3 + nt
            if row:
                rows.append(row)
    return rows, sfiles


def _indep_eval_row(struct_row, p):
    row = {}
    for coeff, wcol in struct_row:
        cv = _indep_eval_modp(coeff, p, SLICE_D0, SLICE_E0)
        if cv:
            row[wcol] = (row.get(wcol, 0) + cv) % p
            if not row[wcol]:
                del row[wcol]
    return row


t1 = time.time()
struct_rows, _sfiles = _indep_structural(ART, FAM)
n_ok = all(len(struct_rows) == len(rows_by_p[p]) for p in PRIMES)
check("T4_indep_reparse_rowcount", n_ok,
      f"{len(struct_rows)} structural rows == loader row count, both primes "
      f"(cancelled_to_empty==0 asserted by the adapter)")
K_SUB = 25
sub = _indep_seed_indices(len(struct_rows), K_SUB,
                          "trustT4-indep|" + "|".join(sha(f) for f in _sfiles))
ok_sub = n_ok
bad = []
for p in PRIMES:
    for i in sub:
        if _indep_eval_row(struct_rows[i], p) != rows_by_p[p][i]:
            ok_sub = False
            bad.append((p, i))
check("T4_indep_reparse_subsample_foreign", ok_sub and len(sub) == K_SUB,
      f"k={K_SUB} sha-seeded FOREIGN rows x {len(PRIMES)} primes == the "
      f"loader rows the dual-core verify used (NO loader import), "
      f"{round(time.time() - t1, 1)} s" + (f"; bad {bad[:3]}" if bad else ""))

# --- T4 mutation controls -----------------------------------------------
# m1: kira-table coefficient flip in a SCRATCH copy -> compare must FAIL
txt = open(TABLE).read()
OLD_ENT = ("lbl3m2L_k2disp[0,0,0,0,0,1,2,0,0] -> \n"
           " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*((d-2)/(2*eta+10))")
assert txt.count(OLD_ENT) == 1, "kira-table mutation anchor not unique"
MUT_TABLE = os.path.join(OUT, "t4_kira_target_MUT.m")
open(MUT_TABLE, "w").write(txt.replace(OLD_ENT, OLD_ENT.replace("(d-2)", "(d-3)")))
t1 = time.time()
rm = okd.reduce_and_compare(MUT_TABLE, F(97, 23), F(5, 7),
                            targets=[(0, 0, 0, 0, 0, 1, 2, 0, 0)],
                            verbose=False, engine_tag="batteryM")
mutation("T4m1_kira_table_coeff_flip", "T4", "kira_target.m (scratch copy)",
         "(d-2) -> (d-3) in the smallest entry", "single-target compare",
         time.time() - t1, rm["sigma"] is None and rm["n_match"] == 0,
         "match=False both sigmas on the mutated claim")
# m2: witness lam tamper -> dual-core verify must FAIL
p = PRIMES[0]
q = sorted(glob.glob(os.path.join(EMIT_DIR, f"w_{FAM}_{p}_*.json")))[0]
wt = core_a.load_witness(q)
wt["lam_val"][0] = (wt["lam_val"][0] + 1) % p
fp = core_a.system_fingerprint(rows_by_p[p], p)
tok, _ = core_a.verify_row(rows_by_p[p], wt, check_fingerprint=fp)
mutation("T4m2_witness_lam_tamper", "T4", "pilot-B witness (in-memory copy)",
         "lam_val[0] += 1", "core verify_row", 0.0, not tok)
# m3: bar-grader teeth
bad = json.loads(json.dumps(recs["PILOT_B"]))
bad["B1_certified"] = "83/84"
mutation("T4m3_bar_grader_teeth", "T4", "PILOT_B receipt (scratch copy)",
         "B1_certified 84/84 -> 83/84", "bars_pilot_b", 0.0, not bars_pilot_b(bad))
# m4: adapter-row value flip -> the independent re-parse comparator must FAIL
_p0, _i0 = PRIMES[0], sub[0]
_tampered = dict(rows_by_p[_p0][_i0])
_kk = next(iter(_tampered))
_tampered[_kk] = (_tampered[_kk] + 1) % _p0
mutation("T4m4_indep_reparse_comparator_teeth", "T4",
         "adapter row (in-memory copy)", "row[k] += 1 mod p",
         "battery-side stdlib re-parse", 0.0,
         _indep_eval_row(struct_rows[_i0], _p0) != _tampered)
leg_wall("T4", t0)
save()

# =========================== T5: witness v1.0 round-trip ====================
t0 = time.time()
_t5_cpu0 = time.process_time()
print("== T5 witness v1.0 round-trip (pilot C as a leg; schema-byte) ==", flush=True)
# the emit path honors the output-root law — T5 emits under the
# validated TRUST_OUT_ROOT (out-of-root emission now refuses typed, see T8)
T5_DIR = os.path.join(os.environ["TRUST_OUT_ROOT"], "t5_witnesses")
os.makedirs(T5_DIR, exist_ok=True)
rows_modp_by_p = {}
nw = 0
emit_fail = None
for p in PRIMES:
    rows_modp = okd.tower_relation_rows_modp(tw, p)
    if len(rows_modp) != r1["n_relations"]:
        emit_fail = f"mod-{p} row drop: {len(rows_modp)} != {r1['n_relations']}"
        break
    rows_modp_by_p[p] = rows_modp
    for nu in p1_targets:
        K = okd.support(nu)
        mK = tuple(nu[i - 1] - 1 for i in K)
        cid = tw.colid[tw.canon(K, mK)]
        red = tw.expand(cid)
        c_modp = {col: wb.frac_mod(v, p) for col, v in red.items() if v}
        wpath = os.path.join(T5_DIR, f"w_lpsyz_{p}_{'_'.join(map(str, nu))}.json")
        prov = {"family": FAM, "tower_top": r1["top"],
                "engine": "lp_syz (pinned vendored, Singular syz + tower descent)",
                "space": "Lee-Pomeransky J-space column ids (tower colid)",
                "d0": "97/23", "eta0": "5/7",
                "relation_rows": "lp_syz syzygy-descent rows mod p, node order"}
        try:
            wb.emit_lpsyz_witness(wpath, rows_modp, cid, c_modp, p, family=FAM,
                                  point={"d": wb.frac_mod(F(97, 23), p),
                                         "eta": wb.frac_mod(F(5, 7), p)},
                                  target_label=list(nu), provenance=prov)
            nw += 1
        except ValueError as e:
            emit_fail = str(e)[:150]
            break
check("T5_emit_14_witnesses", emit_fail is None and nw == 14,
      emit_fail or f"14 witnesses, {round(time.time() - t0, 1)} s")
# stub-detector floor on CPU time, NOT wall (re-tuned: 2.7 s CPU measured
# on a fast box — the old 3.0 s floor false-redded real work there; a stub
# is near-zero CPU, so 1.0 s keeps the detector's teeth)
_t5_cpu = time.process_time() - _t5_cpu0
check("T5_emit_cpu_floor", _t5_cpu >= 1.0,
      f"{round(_t5_cpu, 1)} s CPU >= 1 s (14 solve_lam eliminations cannot "
      f"be sub-1 s CPU; wall {round(time.time() - t0, 1)} s)")

n = pa = pb = 0
for p in PRIMES:
    if p not in rows_modp_by_p:
        continue
    fp = core_a.system_fingerprint(rows_modp_by_p[p], p)
    for q in sorted(glob.glob(os.path.join(T5_DIR, f"w_lpsyz_{p}_*.json"))):
        wa = core_a.load_witness(q)
        oka, _ = core_a.verify_row(rows_modp_by_p[p], wa, check_fingerprint=fp)
        wbv = core_b.load_witness(q)
        okb, _ = core_b.verify_row(rows_modp_by_p[p], wbv, check_fingerprint=fp)
        n += 1
        pa += bool(oka)
        pb += bool(okb)
check("T5_verify_dual_core_14of14", n == 14 and pa == 14 and pb == 14,
      f"{pa}+{pb}/14+14 through BOTH cores unchanged")

# schema-byte assertions (R5: v1.0 UNCHANGED, no fork)
elim_keys = set(json.load(open(sorted(glob.glob(
    os.path.join(EMIT_DIR, "w_*.json")))[0])).keys())
ok_schema = True
sdetails = []
for q in sorted(glob.glob(os.path.join(T5_DIR, "w_lpsyz_*.json"))):
    w = json.load(open(q))
    keys = set(w.keys())
    if not (V10_REQUIRED <= keys <= V10_TOP_KEYS):
        ok_schema = False
        sdetails.append(f"{os.path.basename(q)}: keys {sorted(keys)}")
    if keys != elim_keys - {"labels"}:
        ok_schema = False
        sdetails.append(f"{os.path.basename(q)}: != eliminator keys minus labels")
    if (str(w.get("receipt_version")) != "1.0" or w.get("kind") != "lambda-witness"
            or w.get("identity") != V10_IDENTITY):
        ok_schema = False
        sdetails.append(f"{os.path.basename(q)}: version/kind/identity bytes")
    if (not set(w["target"]) <= {"col", "label"}
            or set(w["lam"]) != {"n_rows", "idx", "val"}
            or not set(w.get("system", {})) <= {"n_rows", "fingerprint", "source"}):
        ok_schema = False
        sdetails.append(f"{os.path.basename(q)}: nested key sets")
    try:
        prov = json.loads(w["system"]["source"])
        if prov.get("producer") != "trust.lp_syz":
            raise ValueError("producer")
    except Exception as e:
        ok_schema = False
        sdetails.append(f"{os.path.basename(q)}: system.source provenance ({e})")
check("T5_schema_v10_byte_assertions", ok_schema, "; ".join(sdetails)[:250]
      or "keys == v1.0 spec == eliminator-producer minus optional labels; "
         "provenance rides system.source ONLY")

# mutation controls: emit-refusal + planted faults
p = PRIMES[0]
if p in rows_modp_by_p:
    nu = p1_targets[0]
    K = okd.support(nu)
    cid = tw.colid[tw.canon(K, tuple(nu[i - 1] - 1 for i in K))]
    red = tw.expand(cid)
    c_bad = {col: wb.frac_mod(v, p) for col, v in red.items() if v}
    k0 = next(iter(c_bad))
    c_bad[k0] = (c_bad[k0] + 1) % p
    refuse_path = os.path.join(T5_DIR, "REFUSED_never_written.json")
    try:
        wb.emit_lpsyz_witness(refuse_path, rows_modp_by_p[p], cid, c_bad, p,
                              family=FAM, point={"d": 1, "eta": 1})
        caught = False
    except ValueError:
        caught = not os.path.exists(refuse_path)
    mutation("T5m1_emit_refuses_false_claim", "T5", "claimed reduction (copy)",
             "c[k0] += 1 mod p", "emit_lpsyz_witness", 0.0, caught,
             "ValueError raised, nothing written")
    fp = core_a.system_fingerprint(rows_modp_by_p[p], p)
    q = sorted(glob.glob(os.path.join(T5_DIR, f"w_lpsyz_{p}_*.json")))[0]
    w1 = core_a.load_witness(q)
    w1["lam_val"][0] = (w1["lam_val"][0] + 1) % p
    f1, _ = core_a.verify_row(rows_modp_by_p[p], w1, check_fingerprint=fp)
    w2 = core_a.load_witness(q)
    k0 = next(iter(w2["c"]))
    w2["c"][k0] = (w2["c"][k0] + 1) % p
    f2, _ = core_a.verify_row(rows_modp_by_p[p], w2, check_fingerprint=fp)
    w3 = core_a.load_witness(q)
    f3, _ = core_a.verify_row(rows_modp_by_p[p], w3,
                              check_fingerprint="deadbeef" * 8)
    mutation("T5m2_planted_faults_3of3", "T5", "emitted witness (in-memory)",
             "lam+1 / c+1 / fingerprint", "core verify_row", 0.0,
             (not f1) and (not f2) and (not f3),
             f"lam:{not f1} c:{not f2} fp:{not f3}")
leg_wall("T5", t0)
save()

# =========================== T6: fraction-oracle leg (R3) ===================
t0 = time.time()
print("== T6 fraction-oracle leg (pure-stdlib child, FOREIGN subsample) ==",
      flush=True)
T6 = os.path.join(OUT, "t6")
os.makedirs(T6, exist_ok=True)


def fr_s(v):
    return f"{v.numerator}/{v.denominator}"


rows_exact = [row for R in tw.node_order for row in tw.rel_by_node[R]]
rows_ser = [{str(c): fr_s(v) for c, v in row.items()} for row in rows_exact]
# FOREIGN target: sha-seeded over the serialized target list (never "first")
tmat = "trustT6-target|" + hashlib.sha256(
    json.dumps([list(t) for t in p1_targets]).encode()).hexdigest()
nu_f = p1_targets[fo._seed_indices(len(p1_targets), 1, tmat)[0]]
K = okd.support(nu_f)
cid_f = tw.colid[tw.canon(K, tuple(nu_f[i - 1] - 1 for i in K))]
masters_free = ([tw.colid[c] for c in tw.cols[tw.first_master_col:]]
                + [tw.colid[c] for c in tw.free_nonmaster])
red_flint = {c: v for c, v in tw.expand(cid_f).items() if v}
p6 = PRIMES[0]
if p6 not in rows_modp_by_p:
    rows_modp_by_p[p6] = okd.tower_relation_rows_modp(tw, p6)
rmat = "trustT6-rows|" + hashlib.sha256(
    json.dumps(rows_ser, sort_keys=True).encode()).hexdigest()
sub_idx = fo._seed_indices(len(rows_ser), 25, rmat)
# FOREIGN rational point (pilot-C R3 point; never used by flint fits)
t1 = time.time()
rf = okd.reduce_and_compare(TABLE, F(103, 31), F(7, 13), targets=[nu_f],
                            verbose=False, engine_tag="batteryF")
check("T6_foreign_point_flint_route_1of1",
      rf["sigma"] == -1 and rf["n_match"] == 1,
      f"foreign (103/31,7/13) tower {round(time.time() - t1, 1)} s")
twf = rf["_tower"]
rows_f = [row for R in twf.node_order for row in twf.rel_by_node[R]]
cid_ff = twf.colid[twf.canon(K, tuple(nu_f[i - 1] - 1 for i in K))]
mf = ([twf.colid[c] for c in twf.cols[twf.first_master_col:]]
      + [twf.colid[c] for c in twf.free_nonmaster])
ct = twf.conv(twf.cols[cid_ff], -1)
conv_ratio = {str(c): fr_s(twf.conv(twf.cols[c], -1) / ct) for c in mf}
kmap = {}
for knu, kex in ents[nu_f].items():
    Kk = okd.support(knu)
    kc = twf.canon(Kk, tuple(knu[i - 1] - 1 for i in Kk))
    kmap[twf.colid[kc]] = kmap.get(twf.colid[kc], F(0)) + okd.eval_coeff(
        kex, F(103, 31), F(7, 13))
export = {
    "anchor": {"rows": rows_ser, "cid": cid_f, "masters": masters_free,
               "red_flint": {str(c): fr_s(v) for c, v in red_flint.items()},
               "p": p6,
               "rows_modp": [{str(c): v for c, v in row.items()}
                             for row in rows_modp_by_p[p6]],
               "sub_k": 25, "sub_idx": sub_idx},
    "foreign": {"rows": [{str(c): fr_s(v) for c, v in row.items()}
                         for row in rows_f],
                "cid": cid_ff, "masters": mf, "conv_ratio": conv_ratio,
                "kmap": {str(c): fr_s(v) for c, v in kmap.items() if v},
                "foreign_target": list(nu_f)},
}
EXP = os.path.join(T6, "t6_export.json")
json.dump(export, open(EXP, "w"))
T6_CHILD = os.path.join(T6, "t6_child.py")
open(T6_CHILD, "w").write('''\
"""T6 child: pure-stdlib Fraction-oracle leg. Runs under python3 -I -S; loads
fraction_oracle.py BY PATH with a sha assert; sweeps sys.modules at the end —
anything beyond stdlib + the path-loaded module is a FAIL."""
import hashlib, importlib.util, json, sys, time
from fractions import Fraction

sys.dont_write_bytecode = True   # -I implies -E: env var form is inert; never
                                 # write pycache into the PACKAGE tree
fo_path, fo_sha, exp_path, out_path = sys.argv[1:5]
legs = sys.argv[5] if len(sys.argv) > 5 else "all"
src = open(fo_path, "rb").read()
assert hashlib.sha256(src).hexdigest() == fo_sha, "fraction_oracle sha mismatch"
spec = importlib.util.spec_from_file_location("fraction_oracle_t6", fo_path)
fo = importlib.util.module_from_spec(spec)
sys.modules["fraction_oracle_t6"] = fo
spec.loader.exec_module(fo)


def F(s):
    n, _, d = s.partition("/")
    return Fraction(int(n), int(d or 1))


exp = json.load(open(exp_path))
res = {}
A = exp["anchor"]
rows = [{int(c): F(v) for c, v in r.items()} for r in A["rows"]]
t0 = time.time()
c0 = time.process_time()   # (T5-floor) class, T6 instance: the
                           # stub-floor must be load-insensitive (CPU, not wall)
red = fo.reduce_column(rows, A["cid"], A["masters"], max_rows=600, max_cols=600)
res["leg1_cpu_s"] = round(time.process_time() - c0, 2)
res["leg1_wall_s"] = round(time.time() - t0, 1)
red_flint = {int(c): F(v) for c, v in A["red_flint"].items()}
res["leg1_equal"] = (red == red_flint)
if legs == "all":
    B = exp["foreign"]
    rowsf = [{int(c): F(v) for c, v in r.items()} for r in B["rows"]]
    t0 = time.time()
    redf = fo.reduce_column(rowsf, B["cid"], B["masters"],
                            max_rows=600, max_cols=600)
    mine_I = {}
    for c, v in (redf or {}).items():
        w = v * F(B["conv_ratio"][str(c)])
        if w:
            mine_I[c] = w
    kmap = {int(c): F(v) for c, v in B["kmap"].items()}
    res["leg2_wall_s"] = round(time.time() - t0, 1)
    res["leg2_equal"] = (mine_I == kmap)
    # leg3: sha-seeded FOREIGN row subsample, exact -> mod p == witness system
    p = A["p"]
    mat = "trustT6-rows|" + hashlib.sha256(
        json.dumps(A["rows"], sort_keys=True).encode()).hexdigest()
    idx = fo._seed_indices(len(rows), A["sub_k"], mat)
    res["leg3_idx_match"] = (idx == A["sub_idx"])
    ok3 = True
    for i in idx:
        got = {}
        for c, v in rows[i].items():
            den = v.denominator % p
            if den == 0:
                ok3 = False
                break
            val = (v.numerator % p) * pow(den, p - 2, p) % p
            if val:
                got[c] = val
        want = {int(c): w for c, w in A["rows_modp"][i].items()}
        if got != want:
            ok3 = False
    res["leg3_modp_consistent"] = ok3
offend = []
for name in list(sys.modules):
    top = name.split(".")[0]
    if top in sys.stdlib_module_names:
        continue
    if top in ("fraction_oracle_t6", "__main__"):
        continue
    offend.append(name)
res["stdlib_offenders"] = offend
res["stdlib_only"] = not offend
json.dump(res, open(out_path, "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "stdlib_offenders"}))
''')
T6_OUTJ = os.path.join(T6, "t6_result.json")
try:
    os.unlink(T6_OUTJ)          # stale-readback guard
except FileNotFoundError:
    pass
t1 = time.time()
rc = child([sys.executable, "-I", "-S", T6_CHILD, os.path.join(PKG, "fraction_oracle.py"),
            PACKAGE_PINS["fraction_oracle.py"], EXP, T6_OUTJ, "all"], timeout=600)
t6w = round(time.time() - t1, 1)
try:
    t6r = json.load(open(T6_OUTJ))
except Exception:
    t6r = {}
check("T6_child_ran", rc.returncode == 0 and bool(t6r),
      f"{t6w} s; rc={rc.returncode} {rc.stderr[-120:] if rc.returncode else ''}")
check("T6_leg1_flint_vs_fraction_equal", t6r.get("leg1_equal") is True,
      f"anchor point, {t6r.get('leg1_wall_s')} s pure-Fraction RREF")
check("T6_leg2_fraction_reelim_vs_kira_equal", t6r.get("leg2_equal") is True,
      f"(103/31,7/13) Fraction RE-ELIMINATION of the tower's exact rows vs "
      f"kira, {t6r.get('leg2_wall_s')} s, target {list(nu_f)} — tower build/"
      f"classification + kira eval SHARED with the flint route ; the head-to-head FLINT cross-gate is leg1")
check("T6_leg3_foreign_row_subsample", t6r.get("leg3_idx_match") is True
      and t6r.get("leg3_modp_consistent") is True,
      "k=25 sha-seeded rows: exact -> mod p == the witness-system rows")
check("T6_stdlib_only_sweep", t6r.get("stdlib_only") is True,
      str(t6r.get("stdlib_offenders", "?"))[:150])
check("T6_cpu_floor", (t6r.get("leg1_cpu_s") or 0) >= 0.3,
      f"leg1 {t6r.get('leg1_cpu_s')} s CPU >= 0.3 s (wall {t6r.get('leg1_wall_s')} s "
      f"reported, not floored — the T5 CPU-floor rule extended to T6; "
      f"0.96 s CPU measured on a fast box, so the floor sits well below "
      f"real work and far above a stub)")
# mutation: corrupted flint reduction in a SCRATCH export copy -> leg1 FAILS
mexp = json.loads(json.dumps(export))
k0 = next(iter(mexp["anchor"]["red_flint"]))
nn, _, dd = mexp["anchor"]["red_flint"][k0].partition("/")
mexp["anchor"]["red_flint"][k0] = f"{int(nn) + 1}/{dd or 1}"
MEXP = os.path.join(T6, "t6_export_MUT.json")
json.dump(mexp, open(MEXP, "w"))
MOUTJ = os.path.join(T6, "t6_result_MUT.json")
try:
    os.unlink(MOUTJ)            # stale-readback guard
except FileNotFoundError:
    pass
t1 = time.time()
rc = child([sys.executable, "-I", "-S", T6_CHILD, os.path.join(PKG, "fraction_oracle.py"),
            PACKAGE_PINS["fraction_oracle.py"], MEXP, MOUTJ, "1"], timeout=300)
try:
    mres = json.load(open(MOUTJ))
except Exception:
    mres = {}
mutation("T6m1_corrupted_flint_reduction", "T6", "t6 export (scratch copy)",
         "red_flint[k0] numerator += 1", "child leg1 equality",
         time.time() - t1, mres.get("leg1_equal") is False,
         f"pure-Fraction leg refuses the corrupted flint claim "
         f"[child_rc={rc.returncode}]")
leg_wall("T6", t0)
save()

# =========================== T7: identity-front integrity ===================
t0 = time.time()
print("== T7 identity fronts (live homes; one-lineage caveat) ==", flush=True)
smod = trust.strata.load()
check("T7_strata_front_resolves_live",
      os.path.abspath(smod.__file__).startswith(trust.LINKED["strata"] + os.sep)
      or os.path.dirname(os.path.abspath(smod.__file__)) == trust.LINKED["strata"],
      smod.__file__)
check("T7_receipt_core_identity",
      os.path.abspath(core_a.__file__) ==
      os.path.join(trust.LINKED["receipt"], "core.py"), core_a.__file__)
check("T7_winnow_core_identity",
      os.path.abspath(core_b.__file__) == os.path.join(
          trust.LINKED["ibplapper"], "ibplapper", "receipt", "core.py"),
      core_b.__file__ + " (exact-path, never a substring check)")
wmod = trr.winnow()
check("T7_winnow_pkg_identity",
      os.path.abspath(wmod.__file__).startswith(trust.LINKED["ibplapper"] + os.sep),
      wmod.__file__)
okbi2, _ = trr.byte_identity()
check("T7_cores_byte_identical_recheck", okbi2)
doc = trust.strata.__doc__ or ""
check("T7_one_lineage_caveat_in_docstring",
      "ONE-LINEAGE CAVEAT" in doc and "never present as independence" in doc
      and "VERBATIM LIFTS" in doc, "strata/winnow one-lineage, import-level")
check("T7_arbitration_table_in_docstring",
      "ARBITRATION TABLE (fixed" in doc and "END-TO-END certification" in doc
      and "INDEPENDENT exact oracle" in doc, "fixed design decision packaged, not re-litigated")
doc_r = trust.receipt.__doc__ or ""
_manual_txt = open(os.path.join(TOOL, "MANUAL.md")).read()
check("T7_check2_substrate_caveat_present",
      "STRATA-LOADER SUBSTRATE CAVEAT" in doc_r
      and "Never present check 2 as executed" in doc_r
      and "STRATA-LOADER SUBSTRATE CAVEAT" in _manual_txt,
      "check-2 as-executed caveat named in trust.receipt "
      "docstring AND MANUAL.md, parallel to the FLINT caveat")

# mutation controls: both shadow classes must REFUSE (child-mode)
T7D = os.path.join(OUT, "t7_fake")
os.makedirs(os.path.join(T7D, "strata"), exist_ok=True)
open(os.path.join(T7D, "strata", "__init__.py"), "w").write("SHADOW = True\n")
SH1 = os.path.join(OUT, "t7_shadow1.py")
open(SH1, "w").write("""\
import json, sys
fake_root, tool = sys.argv[1], sys.argv[2]
sys.path.insert(0, fake_root)
import strata  # the SHADOW, pre-imported
sys.path.insert(0, tool)
import trust
try:
    trust.strata.load()
    print(json.dumps({"raised": None}))
except ImportError as e:
    print(json.dumps({"raised": "ImportError",
                      "violation": "identity violation" in str(e)}))
""")
r = child([sys.executable, "-P", SH1, T7D, TOOL])
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}
mutation("T7m1_preimported_shadow_refused", "T7", "sys.path (child scratch)",
         "fake strata pkg pre-imported from scratch dir", "alias_identity",
         0.0, v.get("raised") == "ImportError" and v.get("violation") is True,
         str(v)[:120])
SH2 = os.path.join(OUT, "t7_shadow2.py")
open(SH2, "w").write("""\
import json, sys, types
tool = sys.argv[1]
sys.path.insert(0, tool)
import trust
m = types.ModuleType("strata")
m.__file__ = trust.LINKED["strata"] + "/__init__.py"  # in-tree path, no spec
sys.modules["strata"] = m
try:
    trust.strata.load()
    print(json.dumps({"raised": None}))
except ImportError as e:
    print(json.dumps({"raised": "ImportError",
                      "violation": "identity violation" in str(e)}))
""")
r = child([sys.executable, "-P", SH2, TOOL])
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}
mutation("T7m2_spec_spoof_refused", "T7", "sys.modules (child)",
         "ModuleType with in-tree __file__ but no import spec",
         "alias_identity spec-authentication", 0.0,
         v.get("raised") == "ImportError" and v.get("violation") is True,
         str(v)[:120])
leg_wall("T7", t0)
save()

# ============== T7C: lineage census (the SCOPE independence leg) ===============
# SCOPE.md: "Battery ships an import-graph + sha census proving the three
# verifier lineages disjoint." HONESTLY scoped: the CORES are
# disjoint (censused below); the checks AS EXECUTED share the strata loader
# on check-2's row path — that edge is FOUND and DOCUMENTED, never denied.
t0 = time.time()
print("== T7C lineage census (import-graph + sha, as-executed) ==", flush=True)
_STDLIB = set(sys.stdlib_module_names)
_STRATA_LINEAGE = {"strata", "loader", "fp_eliminate", "corpus_bench",
                   "coverage", "weights", "bank", "engine", "ibplapper"}


def _imports_of(path):
    tops = set()
    for nd in ast.walk(ast.parse(open(path).read())):
        if isinstance(nd, ast.Import):
            for a in nd.names:
                tops.add(a.name.split(".")[0])
        elif isinstance(nd, ast.ImportFrom):
            tops.add("." if nd.level else (nd.module or "").split(".")[0])
    return tops


_CORE2 = os.path.join(trust.LINKED["receipt"], "core.py")
_CORE2W = os.path.join(trust.LINKED["ibplapper"], "ibplapper", "receipt", "core.py")
_AD2 = os.path.join(trust.LINKED["receipt"], "adapters", "strata.py")
_c2_tops = _imports_of(_CORE2)
check("T7C_check2_core_stdlib_only", _c2_tops <= _STDLIB,
      f"receipt core imports {sorted(_c2_tops)} — all stdlib")

_C3_FILES = ([trust.vendored_path(nm) for nm in
              ("lp_syz", "lp_syz_431", "lp_syz_prod", "bessel_oracle2")]
             + [os.path.join(PKG, f) for f in
                ("oracle_k2disp.py", "fraction_oracle.py", "witness_bridge.py")])
_bad3 = {}
for _f in _C3_FILES:
    _hit = _imports_of(_f) & _STRATA_LINEAGE
    if _hit:
        _bad3[os.path.basename(_f)] = sorted(_hit)
check("T7C_check3_no_eliminator_imports", not _bad3,
      str(_bad3) if _bad3 else
      "7 files: zero strata/loader/fp_eliminate/ibplapper imports "
      "(sympy/flint = the declared FLINT-substrate caveat; '.' = trust-internal)")

_ad_tops = _imports_of(_AD2)
_edge = {"loader", "fp_eliminate"} <= _ad_tops
check("T7C_adapter_edge_found_and_documented",
      _edge and "STRATA-LOADER SUBSTRATE CAVEAT" in (trust.receipt.__doc__ or "")
      and "STRATA-LOADER SUBSTRATE CAVEAT" in _manual_txt,
      "adapters/strata.py DOES import loader+fp_eliminate (check-2 as "
      "executed) — the census DOCUMENTS the edge as the named caveat; "
      "independent stdlib re-parse spot-check wired at T4")

_sh3 = {sha(f) for f in _C3_FILES}
_sh2 = {sha(f) for f in (_CORE2, _CORE2W, _AD2)}
_shS = {sha(os.path.join(trust.LINKED["strata"], f))
        for f in ("loader.py", "fp_eliminate.py")
        if os.path.isfile(os.path.join(trust.LINKED["strata"], f))}
check("T7C_sha_census_pairwise_disjoint",
      not (_sh3 & _sh2) and not (_sh3 & _shS) and not (_sh2 & _shS)
      and len(_shS) == 2,
      "check-3 / check-2 / strata file bytes pairwise distinct (designed "
      "exception counted once: the two receipt cores are byte-identical — "
      "ONE lineage by construction, T7_cores_byte_identical)")

# mutation: an injected eliminator import in a scratch oracle copy must flag
_T7C_MUT = os.path.join(OUT, "t7c_mut_oracle.py")
open(_T7C_MUT, "w").write(
    open(os.path.join(PKG, "oracle_k2disp.py")).read()
    + "\nfrom fp_eliminate import eliminate_fast  # MUTATION\n")
_hitm = _imports_of(_T7C_MUT) & _STRATA_LINEAGE
mutation("T7Cm1_injected_eliminator_import_flagged", "T7C",
         "oracle_k2disp.py (scratch copy)",
         "append 'from fp_eliminate import eliminate_fast'",
         "_imports_of census", 0.0, bool(_hitm), f"flags {sorted(_hitm)}")
leg_wall("T7C", t0)
save()

# =========================== T8: refusal battery ============================
t0 = time.time()
print("== T8 refusal battery (typed, value-validated) ==", flush=True)
ENV_BAK = os.environ.pop("TRUST_OUT_ROOT")


def typed_refusal(value):
    """Returns (leaf_type_name, message) from out_root() with env=value."""
    if value is not None:
        os.environ["TRUST_OUT_ROOT"] = value
    try:
        trust.out_root()
        return None, "accepted"
    except Exception as e:
        return type(e).__name__, str(e)
    finally:
        os.environ.pop("TRUST_OUT_ROOT", None)


tname, msg = typed_refusal(None)
check("T8_envless_refusal_typed", tname == "OutputRootError"
      and "TRUST_OUT_ROOT" in msg, f"{tname}: {msg[:80]}")

# --- witness-emit path honors the output-root law -------------
# (TRUST_OUT_ROOT is popped right now: env-less emission must refuse typed
# with NOTHING written — the poc_b class dead)
_rows_pb = [{5: 1, 7: 100}]                     # R0 = e5 + 100*e7
_c_pb = {7: PRIMES[0] - 100}                    # claim e5 == (p-100)*e7 == R0
_wp = os.path.join(OUT, "t8_envless_emit.json")
try:
    wb.emit_lpsyz_witness(_wp, _rows_pb, 5, _c_pb, PRIMES[0], family=FAM,
                          point={"d": 1, "eta": 1})
    _ok_we = False
except trust.OutputRootError:
    _ok_we = not os.path.exists(_wp)
except Exception:
    _ok_we = False
check("T8_emit_envless_outroot_refused", _ok_we,
      "env-less emit_lpsyz_witness refuses OutputRootError, nothing written "
      "(the emit path was the package's only out-root-blind writing path)")

# --- in-process trust import refuses a PYTHONPATH world -------
T8_ENV = os.path.join(OUT, "t8_envpoison_child.py")
open(T8_ENV, "w").write("""\
import json, sys
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
try:
    import trust
    print(json.dumps({"raised": None}))
except Exception as e:
    print(json.dumps({"raised": type(e).__name__, "msg": str(e)[:150]}))
""")
_shadow_dir = os.path.join(OUT, "t8_empty_shadow")
os.makedirs(_shadow_dir, exist_ok=True)
_penv = {k: v for k, v in os.environ.items() if not k.startswith("PYTHON")}
_penv["PYTHONPATH"] = _shadow_dir
_penv["PYTHONDONTWRITEBYTECODE"] = "1"
# deliberately RAW subprocess, NO -E: the poisoned-interpreter world itself
r = subprocess.run([sys.executable, T8_ENV, TOOL], capture_output=True,
                   text=True, timeout=60, env=_penv, cwd=OUT)
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}
check("T8_import_env_refusal_typed",
      v.get("raised") == "EnvPoisonError" and "PYTHONPATH" in v.get("msg", ""),
      f"{v.get('raised')}: in-process trust import refuses the PYTHONPATH "
      f"world ")
# positive control: same world under -E — the interpreter IGNORED the vars,
# the environment did not poison it, import must succeed (battery children)
r = subprocess.run([sys.executable, "-E", T8_ENV, TOOL], capture_output=True,
                   text=True, timeout=60, env=_penv, cwd=OUT)
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}
check("T8_import_env_E_exemption_positive_control", v.get("raised") is None,
      "same PYTHONPATH world under -E imports fine (sys.flags."
      "ignore_environment exemption — children stay unaffected)")

# --- crafted evil table refuses WITHOUT executing -------------
_EVIL = os.path.join(OUT, "t8_evil_table.m")
_PROOF = os.path.join(OUT, "t8_SYMPIFY_EXEC_PROOF")
try:
    os.unlink(_PROOF)
except FileNotFoundError:
    pass
open(_EVIL, "w").write(
    "{\n"
    "lbl3m2L_k2disp[0,0,0,0,0,1,2,0,0] -> \n"
    " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*((d-2)/(2*eta+10) + "
    "0*__import__('pathlib').Path(" + repr(_PROOF) + ").write_text('code ran'))\n"
    "}\n")
try:
    okd.parse_kira_table(_EVIL)
    _ok_ev = False
except ValueError:
    _ok_ev = not os.path.exists(_PROOF)
except Exception:
    _ok_ev = False
check("T8_evil_table_refused_without_execution", _ok_ev,
      "breaker-class payload in a coefficient: ValueError from the locked "
      "rational-function grammar, NO proof file written (sympify retired "
      "from the untrusted-table path)")
# Synthetic probe values (generic stand-ins for the protected
# classes — parked archive, home, protected trees, the tools tree, a bare
# data root, /, and an undeclared-neutral dir; allowlist posture means
# every one must refuse regardless of the specific prefix).
FORBIDDEN_PROBES = {
    "relative": "rel/path", "parked_archive": "/srv/parked_archive/zz",
    "tools_tree": os.path.join(TOOL, "x"), "home": "/home/zz",
    "protected_a": "/srv/protected_a/x",
    "protected_b": "/srv/protected_b/x",
    # values the old denylist ACCEPTED — now refused (allowlist
    # posture; this list is no longer a mirror of the implementation's list)
    "protected_c": "/srv/protected_c/x_probe",
    "protected_d": "/srv/protected_d/x_probe",
    "protected_e": "/srv/protected_e/x_probe",
    "denied_tree": "/srv/denied_tree/x_probe",
    "data_root": "/data/x_probe",
    "fs_root": "/",
    "undeclared_neutral": "/opt/trust_x_probe",
}
ok_f = True
fdet = []
for k, val in FORBIDDEN_PROBES.items():
    tname, msg = typed_refusal(val)
    if tname != "OutputRootError":
        ok_f = False
        fdet.append(f"{k}: {tname}")
check("T8_forbidden_value_refusals_13of13", ok_f, "; ".join(fdet) or
      "relative/parked/tools/home/protected-a..e/"
      "denied-tree/data-root/fs-root/undeclared all refused typed")
link = os.path.join(OUT, "t8_link")
if not os.path.islink(link):
    os.symlink("/srv/parked_archive", link)
tname, msg = typed_refusal(os.path.join(link, "foo"))
check("T8_symlink_resolved_escape_refused", tname == "OutputRootError",
      f"{tname} (realpath resolves under the parked-archive probe root)")
os.environ["TRUST_OUT_ROOT"] = ENV_BAK
check("T8_valid_root_accepted", trust.out_root() == os.path.realpath(ENV_BAK),
      "positive control: scratch root accepted")
# --- (root set): out-of-root emit path refuses, nothing written
_wp2 = os.path.join(OUT, "t8_escape_emit.json")   # OUT is the PARENT of root
try:
    wb.emit_lpsyz_witness(_wp2, _rows_pb, 5, _c_pb, PRIMES[0], family=FAM,
                          point={"d": 1, "eta": 1})
    _ok_esc = False
except trust.OutputRootError:
    _ok_esc = not os.path.exists(_wp2)
except Exception:
    _ok_esc = False
check("T8_emit_path_escape_refused", _ok_esc,
      "witness path outside TRUST_OUT_ROOT refuses OutputRootError")
# --- reserved provenance keys refuse; producer stamp holds ----
_wp3 = os.path.join(trust.out_root(), "t8_forged_prov.json")
try:
    wb.emit_lpsyz_witness(_wp3, _rows_pb, 5, _c_pb, PRIMES[0], family=FAM,
                          point={"d": 1, "eta": 1},
                          provenance={"producer": "strata.fp_eliminate"})
    _ok_fp = False
except ValueError:
    _ok_fp = not os.path.exists(_wp3)
except Exception:
    _ok_fp = False
check("T8_emit_reserved_provenance_refused", _ok_fp,
      "provenance {'producer': ...} refuses ValueError (forged-lineage "
      "guard), nothing written")
_wp4 = os.path.join(trust.out_root(), "t8_emit_ok.json")
try:
    os.unlink(_wp4)
except FileNotFoundError:
    pass
_w4, _ = wb.emit_lpsyz_witness(_wp4, _rows_pb, 5, _c_pb, PRIMES[0], family=FAM,
                               point={"d": 1, "eta": 1},
                               provenance={"note": "T8 positive control"})
_src4 = json.loads(_w4["system"]["source"])
check("T8_emit_positive_control_producer_stamped",
      os.path.exists(_wp4) and _src4.get("producer") == "trust.lp_syz"
      and _src4.get("note") == "T8 positive control",
      "in-root emit with clean provenance succeeds; bridge stamp intact, "
      "caller keys carried")
# the own-tree guard must TRAVEL with a relocated copy — and a
# relocated copy must NOT go spuriously red on run scratch (the baseline-red
# class that muddies attribution)
T8_RELOC = os.path.join(OUT, "t8_reloc_child.py")
open(T8_RELOC, "w").write("""\
import json, os, sys
clone, out = sys.argv[1:3]
sys.path.insert(0, clone)
os.environ["TRUST_OUT_ROOT"] = out
import trust
try:
    r = trust.out_root()
    print(json.dumps({"raised": None, "root": r}))
except Exception as e:
    print(json.dumps({"raised": type(e).__name__, "msg": str(e)[:150]}))
""")
_reloc = os.path.join(CLONES, "pristine")


def _reloc_probe(out_val):
    r = child([sys.executable, "-P", T8_RELOC, _reloc, out_val])
    try:
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}


v = _reloc_probe(os.path.join(_reloc, "vendor", "x_probe"))
check("T8_relocated_own_tree_refused", v.get("raised") == "OutputRootError",
      f"{v.get('raised')}: relocated copy refuses ITS OWN tree even under an "
      f"allowlisted prefix (the guard travels)")
v = _reloc_probe(os.path.join(OUT, "t8_reloc_ok"))
check("T8_relocated_scratch_accepted", v.get("raised") is None,
      "relocated copy accepts run scratch (no spurious relocation red)")
try:
    trust.run_vendored("lp_syz_431", ["--stage", "control452"],
                       cwd=os.path.join(OUT, "t8_cwd_escape"), timeout=30)
    ok_cwd = False
except trust.OutputRootError:
    ok_cwd = True
except Exception:
    ok_cwd = False
check("T8_run_cwd_escape_refused", ok_cwd, "cwd outside TRUST_OUT_ROOT")
try:
    trust.vendored_path("no_such_engine")
    ok_u = False
except KeyError:
    ok_u = True
try:
    trust.load_vendored("no_such_engine")
    ok_u = False
except KeyError:
    pass
except Exception:
    ok_u = False
check("T8_unknown_engine_typed", ok_u, "KeyError from vendored_path AND load_vendored")
# fixture note: ADAPTER_PROBE.md (the old junk fixture) contains
# em-dashes, so it now refuses at the BYTE-ALPHABET gate — asserted below as
# a real-file specimen; the 0-entry/junk-driver legs ride an ASCII junk file.
_JUNK = os.path.join(OUT, "t8_ascii_junk.m")
open(_JUNK, "w").write("this is not a kira table\nno records here at all\n")
check("T8_unknown_fixture_no_entries",
      len(okd.parse_kira_table(_JUNK)) == 0,
      "ASCII non-table file parses to 0 entries at the PARSE level (nothing "
      "silently invented; the driver-level hole is closed below)")
# 0: the vacuous success-shaped dict path is CLOSED — junk table
# and requested-but-absent targets refuse TYPED at the driver
try:
    okd.reduce_and_compare(_JUNK, F(97, 23), F(5, 7))
    ok_j = False
except ValueError:
    ok_j = True
except Exception:
    ok_j = False
check("T8_junk_fixture_typed_refusal", ok_j,
      "reduce_and_compare(non-table) raises ValueError — no sigma=1/"
      "n_violations=0 success shape from a 0-entry parse")
try:
    okd.reduce_and_compare(RECEIPT_PINS["KIRA_TABLE"][0], F(97, 23), F(5, 7),
                           targets=[(9,) * 9])
    ok_a = False
except KeyError:
    ok_a = True
except Exception:
    ok_a = False
check("T8_absent_target_typed_refusal", ok_a,
      "requested-but-absent target raises KeyError (silent filter removed)")
try:
    okd.parse_kira_table(os.path.join(OUT, "no_such_table.m"))
    ok_m = False
except FileNotFoundError:
    ok_m = True
check("T8_missing_fixture_typed", ok_m, "FileNotFoundError")

# --- cwd/script-dir stdlib shadow (env-FREE) refuses typed ----
# sys.path[0] is the script's dir in script mode and '' (cwd) under -c/REPL/
# stdin — NO poisoned env var, so an env gate has nothing to fire on; a
# pins-aware passthrough hashlib.py planted there could green a tampered
# vendor. Under test: post-import stdlib IDENTITY vs the
# interpreter-owned roots — identity, not vector enumeration.
T8_SHADOW = os.path.join(OUT, "t8_cwd_shadow")
os.makedirs(T8_SHADOW, exist_ok=True)
open(os.path.join(T8_SHADOW, "hashlib.py"), "w").write(f"""\
# passthrough stdlib shadow (battery probe): everything WORKS, provenance lies
import importlib.util as _iu
_spec = _iu.spec_from_file_location("_real_hashlib", {hashlib.__file__!r})
_real = _iu.module_from_spec(_spec)
_spec.loader.exec_module(_real)
sha256 = _real.sha256
def __getattr__(n):
    return getattr(_real, n)
""")
T8_SHADOW_BODY = """\
import json, sys
sys.dont_write_bytecode = True
sys.path.insert(1, sys.argv[1])
import hashlib
out = {"shadowed": "t8_cwd_shadow" in (getattr(hashlib, "__file__", "") or "")}
try:
    import trust
    out["raised"] = None
except Exception as e:
    out["raised"] = type(e).__name__
    out["msg"] = str(e)[:220]
print(json.dumps(out))
"""
T8_SHADOW_PROBE = os.path.join(T8_SHADOW, "t8_shadow_probe.py")
open(T8_SHADOW_PROBE, "w").write(T8_SHADOW_BODY)


def _shadow_probe(script, pkg_root):
    # RAW subprocess ON PURPOSE: child() forces -P, which closes the very
    # script-dir sys.path[0] channel this leg probes.
    r = subprocess.run([sys.executable, script, pkg_root],
                       capture_output=True, text=True, timeout=120,
                       env=dict(os.environ), cwd=OUT)
    try:
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}


v = _shadow_probe(T8_SHADOW_PROBE, TOOL)
check("T8_cwd_shadow_import_refused_typed",
      v.get("shadowed") is True and v.get("raised") == "StdlibShadowError"
      and "hashlib" in v.get("msg", "") and "t8_cwd_shadow" in v.get("msg", ""),
      f"{v.get('raised')}: script-dir hashlib shadow, ZERO env vars — trust "
      f"import refuses typed naming the offending path ")
T8_CLEAN = os.path.join(OUT, "t8_cwd_clean")
os.makedirs(T8_CLEAN, exist_ok=True)
T8_CLEAN_PROBE = os.path.join(T8_CLEAN, "t8_shadow_probe.py")
open(T8_CLEAN_PROBE, "w").write(T8_SHADOW_BODY)
v = _shadow_probe(T8_CLEAN_PROBE, TOOL)
check("T8_cwd_clean_import_positive_control",
      v.get("shadowed") is False and v.get("raised") is None,
      "same probe from a clean script dir imports fine (no spurious red)")
# mutation control: gate-call removed in a scratch copy -> the shadow world
# imports CLEAN (greened) — the refusal above is the gate's doing
_t0m = time.time()
_c1 = make_clone("r3_gateoff")
_ip = os.path.join(_c1, "trust", "__init__.py")
_itxt = open(_ip).read()
_anchor = ("_stdlib_identity_gate()   "
           "# GATE-CALL (battery mutation anchor)")
assert _anchor in _itxt, "gate-call anchor missing in trust/__init__.py"
open(_ip, "w").write(_itxt.replace(_anchor,
                                   "pass   # MUTANT: F1 gate call removed"))
v = _shadow_probe(T8_SHADOW_PROBE, _c1)
mutation("T8m_f1_gate_removed", "T8", "copy r3_gateoff trust/__init__.py",
         "stdlib-identity gate call removed",
         "script-dir pins-aware hashlib shadow world", time.time() - _t0m,
         v.get("shadowed") is True and v.get("raised") is None,
         "gate-off copy imports CLEAN through the shadow world — the "
         "refusal is the gate's doing, not an accident of the fixture")

# --- parse coverage — out-of-grammar table content REFUSES ----
# instead of being silently DROPPED (a table is never certified against a
# parsed SUBSET of its bytes; refusals name the offending line numbers)
_T8HDR = "lbl3m2L_k2disp[0,0,0,0,0,1,2,0,0] -> \n"
_T8T1 = " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*((d-2)/(2*eta+10))\n"
_dt = os.path.join(OUT, "t8_dropterm.m")
open(_dt, "w").write("{\n" + _T8HDR + _T8T1
                     + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]* (999)\n}\n")


def _parse_refusal(path):
    try:
        ents = okd.parse_kira_table(path)
        return None, f"PARSED {len(ents)} entries"
    except ValueError as e:
        return "ValueError", str(e)
    except Exception as e:
        return type(e).__name__, str(e)


tname, msg = _parse_refusal(_dt)
check("T8_dropterm_coverage_refused_named_line",
      tname == "ValueError" and "line(s) [4]" in msg,
      f"{tname}: {msg[:110]}")
_de = os.path.join(OUT, "t8_dropent.m")
open(_de, "w").write("{\n" + _T8HDR + _T8T1 + ",\n"
                     "lbl3m2L_k2disp[0,0,0,0,0,2,1,0,0] ->\n"
                     " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*(7)\n}\n")
tname, msg = _parse_refusal(_de)
check("T8_dropent_coverage_refused_named_line",
      tname == "ValueError" and "2 raw '->' records" in msg
      and "line(s) [5]" in msg, f"{tname}: {msg[:110]}")
_PR3 = os.path.join(OUT, "t8_R3_EXEC_PROOF")
try:
    os.unlink(_PR3)
except FileNotFoundError:
    pass
_ev = os.path.join(OUT, "t8_evil_spaced.m")
open(_ev, "w").write(
    "{\n" + _T8HDR + _T8T1
    + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]* ((d-2)/(2*eta+10) + "
    "0*__import__('pathlib').Path(" + repr(_PR3) + ").write_text('ran'))\n}\n")
tname, msg = _parse_refusal(_ev)
check("T8_evil_spaced_payload_refused_not_dropped",
      tname == "ValueError" and not os.path.exists(_PR3),
      "the R2 payload rewritten with '* (' spacing REFUSES typed (a permissive "
      "parser silently drops it and certifies the subset); nothing executed")
_t0p = time.time()
check("T8_coverage_positive_real_table_42",
      len(okd.parse_kira_table(RECEIPT_PINS["KIRA_TABLE"][0])) == 42,
      f"pinned 42-entry table full-coverage parse OK "
      f"({time.time() - _t0p:.1f} s — the assertion rejects nothing real)")
# mutation control: both coverage refusals disabled in a scratch copy ->
# dropterm parses to 1 term and the spaced payload is silently dropped
_t0m = time.time()
_c2 = make_clone("r3_coverageoff")
_op = os.path.join(_c2, "trust", "oracle_k2disp.py")
_otxt = open(_op).read()
_aA = ("if len(ents) != len(raw_recs):   "
       "# GUARD COVERAGE (battery anchor A)")
_aB = "if bad_lines:                # GUARD COVERAGE (battery anchor B)"
assert _aA in _otxt and _aB in _otxt, "F2 coverage anchors missing"
open(_op, "w").write(_otxt.replace(_aA, "if False:   # MUTANT anchor A")
                     .replace(_aB, "if False:   # MUTANT anchor B"))
T8_COVER_CHILD = os.path.join(OUT, "t8_cover_child.py")
open(T8_COVER_CHILD, "w").write("""\
import json, sys
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
from trust import oracle_k2disp as okd
out = {}
for tag, path in (("dropterm", sys.argv[2]), ("evil", sys.argv[3])):
    try:
        ents = okd.parse_kira_table(path)
        out[tag] = {"refused": False,
                    "n_terms": max(len(t) for t in ents.values())}
    except ValueError:
        out[tag] = {"refused": True}
print(json.dumps(out))
""")
r = child([sys.executable, "-P", T8_COVER_CHILD, _c2, _dt, _ev], timeout=180)
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"err": (r.stderr or r.stdout)[-200:]}
mutation("T8m_f2_coverage_removed", "T8",
         "copy r3_coverageoff oracle_k2disp.py",
         "both parse-coverage refusals disabled",
         "dropterm + spaced-payload tables", time.time() - _t0m,
         v.get("dropterm", {}).get("refused") is False
         and v.get("dropterm", {}).get("n_terms") == 1
         and v.get("evil", {}).get("refused") is False
         and v.get("evil", {}).get("n_terms") == 1
         and not os.path.exists(_PR3),
         "coverage-off copy silently certifies the 1-term SUBSET of both "
         "tables — the refusal is the coverage gate's doing")

# --- in-grammar resource bombs + div0 refuse FAST + typed -----
_bm = os.path.join(OUT, "t8_bomb.m")
open(_bm, "w").write("{\n" + _T8HDR
                     + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]"
                       "*(2^1000000000000000000)\n}\n")
_t0b = time.time()
tname, msg = _parse_refusal(_bm)
_wallb = time.time() - _t0b
check("T8_pow_bomb_refused_fast_typed",
      tname == "ValueError" and "cap" in msg and _wallb < 5.0,
      f"60-byte 2^10^18 bomb: {tname} in {_wallb * 1000:.0f} ms (an uncapped "
      f"eval is timeout-killed at 30 s under RLIMIT_AS 2 GiB)")
_oc = os.path.join(OUT, "t8_overcap.m")
open(_oc, "w").write("{\n" + _T8HDR
                     + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*(2^10001)\n}\n")
tname, msg = _parse_refusal(_oc)
check("T8_pow_overcap_refused_typed", tname == "ValueError" and "cap" in msg,
      f"{tname}: exponent 10001 > 10^4 cap")
_dv = os.path.join(OUT, "t8_div0.m")
open(_dv, "w").write("{\n" + _T8HDR
                     + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*(1/0)\n}\n")
tname, msg = _parse_refusal(_dv)
check("T8_div0_refused_typed_at_build",
      tname == "ValueError" and "zero denominator" in msg,
      f"{tname}: {msg[:90]} (zoo never constructed)")
_pc = os.path.join(OUT, "t8_pow_ok.m")
open(_pc, "w").write("{\n" + _T8HDR
                     + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]"
                       "*((2^100*d^3+1)/(3*eta+1))\n}\n")
try:
    _e = okd.parse_kira_table(_pc)
    _vv = okd.eval_coeff(_e[(0, 0, 0, 0, 0, 1, 2, 0, 0)]
                         [(0, 0, 0, 0, 0, 1, 1, 0, 0)], F(3, 2), F(1, 2))
    _okp = _vv == (F(2) ** 100 * F(27, 8) + 1) / F(5, 2)
except Exception:
    _okp = False
check("T8_pow_within_cap_positive_control", _okp,
      "in-cap numeric+symbolic powers parse AND evaluate exactly")
try:
    okd.eval_coeff(okd._safe_rational_expr("1/(eta-5)"), F(97, 23), F(5))
    _okpole = False
except ValueError:
    _okpole = True
except Exception:
    _okpole = False
check("T8_eval_pole_refused_typed", _okpole,
      "denominator vanishing AT the evaluation point refuses ValueError "
      "(never an untyped TypeError via zoo)")
# mutation control: caps + div0 check removed in a scratch copy -> overcap
# parses fine and 1/0 parses to zoo without raising
_t0m = time.time()
_c3 = make_clone("r3_capoff")
_op3 = os.path.join(_c3, "trust", "oracle_k2disp.py")
_otxt3 = open(_op3).read()
_aC = "_MAX_POW_EXP = 10 ** 4"
_aD = "if b.is_zero:   # GUARD DIV0 (battery mutation anchor)"
assert _aC in _otxt3 and _aD in _otxt3, "F3 anchors missing"
open(_op3, "w").write(_otxt3.replace(_aC, "_MAX_POW_EXP = 10 ** 30")
                      .replace(_aD, "if False:   # MUTANT div0 removed"))
T8_CAP_CHILD = os.path.join(OUT, "t8_cap_child.py")
open(T8_CAP_CHILD, "w").write("""\
import json, sys
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
from trust import oracle_k2disp as okd
out = {}
try:
    okd.parse_kira_table(sys.argv[2])
    out["overcap"] = {"refused": False}
except ValueError:
    out["overcap"] = {"refused": True}
try:
    ents = okd.parse_kira_table(sys.argv[3])
    cx = list(list(ents.values())[0].values())[0]
    out["div0"] = {"refused": False, "zoo": "zoo" in str(cx)}
except ValueError:
    out["div0"] = {"refused": True}
print(json.dumps(out))
""")
r = child([sys.executable, "-P", T8_CAP_CHILD, _c3, _oc, _dv], timeout=180)
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"err": (r.stderr or r.stdout)[-200:]}
mutation("T8m_f3_caps_removed", "T8", "copy r3_capoff oracle_k2disp.py",
         "exponent cap lifted to 10^30 + div0 build check removed",
         "overcap (2^10001) + div0 (1/0) tables", time.time() - _t0m,
         v.get("overcap", {}).get("refused") is False
         and v.get("div0", {}).get("refused") is False
         and v.get("div0", {}).get("zoo") is True,
         "cap-off copy evaluates the over-cap power and silently builds "
         "zoo — both refusals are the gate's doing")

# --- site-packages .pth reorder (trusted-root shadow) ---------
# site.py executes .pth 'import' lines at interpreter startup, BEFORE any
# user code — a 2-line .pth in site-packages reorders sys.path so a planted
# hashlib WINS import from a dir the old gate whitelisted (R4 W2 greened a
# tampered vendor). Cure under test: STDLIB-roots-only authentication
# (BASE-installation sysconfig vars) + the site-segment belt.
T8_VENV = os.path.join(OUT, "t8_venv")
if not os.path.isdir(T8_VENV):
    subprocess.run([sys.executable, "-m", "venv", "--without-pip", T8_VENV],
                   capture_output=True, timeout=120)
_VPY = os.path.join(T8_VENV, "bin", "python3")
import glob as _glob
_VSP = _glob.glob(os.path.join(T8_VENV, "lib", "python3.*",
                               "site-packages"))[0]
open(os.path.join(_VSP, "hashlib.py"), "w").write(f"""\
# passthrough stdlib shadow (battery probe): everything WORKS, provenance lies
import importlib.util as _iu
_spec = _iu.spec_from_file_location("_real_hashlib", {hashlib.__file__!r})
_real = _iu.module_from_spec(_spec)
_spec.loader.exec_module(_real)
sha256 = _real.sha256
def __getattr__(n):
    return getattr(_real, n)
""")
open(os.path.join(_VSP, "r4_reorder.py"), "w").write(
    f"import sys\nsys.path.insert(0, {_VSP!r})\n")
open(os.path.join(_VSP, "r4reorder.pth"), "w").write("import r4_reorder\n")
T8_VPROBE = os.path.join(T8_VENV, "t8_vprobe.py")
open(T8_VPROBE, "w").write("""\
import json, sys
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
import hashlib
out = {"shadowed": "site-packages" in (getattr(hashlib, "__file__", "") or "")}
try:
    import trust
    out["raised"] = None
except Exception as e:
    out["raised"] = type(e).__name__
    out["msg"] = str(e)[:220]
print(json.dumps(out))
""")


def _venv_probe(pkg_root):
    # RAW venv python ON PURPOSE: default site processing must run the .pth
    # (that IS the vector); no -P/-s/-E.
    r = subprocess.run([_VPY, T8_VPROBE, pkg_root], capture_output=True,
                       text=True, timeout=120, env=dict(os.environ), cwd=OUT)
    try:
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"raised": "CHILD-ERROR", "msg": (r.stderr or r.stdout)[-200:]}


v = _venv_probe(TOOL)
check("T8_sitepth_shadow_refused_typed",
      v.get("shadowed") is True and v.get("raised") == "StdlibShadowError"
      and "site-packages" in v.get("msg", ""),
      f"{v.get('raised')}: .pth-reordered site-packages hashlib plant, ZERO "
      f"env vars — trust import refuses typed naming the site dir "
      f"")
# mutation control: roots widened back to the default (venv) scheme + the
# site-segment belt removed in a scratch copy -> the same world imports
# CLEAN — the refusal above is the gate's doing
_t0m = time.time()
_c4 = make_clone("r4_rootsoff")
_ip4 = os.path.join(_c4, "trust", "__init__.py")
_itxt4 = open(_ip4).read()
_aE = ("            p = _sc.get_path(_key, _scheme, vars=_bv)   "
       "# GUARD STDLIB-ROOTS (battery mutation anchor)")
_aF = ('        if any(_seg in _parts for _seg in ("site-packages", '
       '"dist-packages")):   # GUARD SITE-SEG (battery mutation anchor)')
assert _aE in _itxt4 and _aF in _itxt4, "STDLIB-ROOTS/SITE-SEG anchors missing"
open(_ip4, "w").write(
    _itxt4.replace(_aE, "            p = _sc.get_path(_key)   # MUTANT: "
                        "default scheme, in-venv platstdlib whitelisted")
    .replace(_aF, "        if False:   # MUTANT: site-segment belt removed"))
v = _venv_probe(_c4)
mutation("T8m_f1r4_siteroots_widened", "T8", "copy r4_rootsoff trust/__init__.py",
         "stdlib roots -> default (venv) scheme + site-segment belt removed",
         ".pth-reordered venv site-packages hashlib shadow world",
         time.time() - _t0m,
         v.get("shadowed") is True and v.get("raised") is None,
         "roots-widened copy imports CLEAN through the .pth-reordered "
         "shadow world — the refusal is the gate's doing, not the fixture")
# clean control: plants removed -> the venv imports trust fine (no spurious
# red in venv worlds; venvs share the base stdlib)
for _f in ("hashlib.py", "r4_reorder.py", "r4reorder.pth"):
    try:
        os.unlink(os.path.join(_VSP, _f))
    except FileNotFoundError:
        pass
shutil.rmtree(os.path.join(_VSP, "__pycache__"), ignore_errors=True)
v = _venv_probe(TOOL)
check("T8_sitepth_clean_venv_positive_control",
      v.get("shadowed") is False and v.get("raised") is None,
      "same venv with plants removed imports trust fine")

# --- aggregate budgets — rebuild bombs INSIDE the R3 caps -----
# The R3 budgets bounded single Pow results; Mult/Div chains, symbolic-
# exponent merges, eager Pow-over-Mul-heads and Add floods rebuilt the
# refused values from in-cap pieces (measured 6-60+ s, multi-GiB). All must
# now refuse typed on CPU-time floors (wall floors false-red on a loaded
# box — the lesson).


def _timed_refusal(path):
    c0 = time.process_time()
    tname, msg = _parse_refusal(path)
    return tname, msg, time.process_time() - c0


_r4mc = os.path.join(OUT, "t8_r4_multchain.m")
open(_r4mc, "w").write("{\n" + _T8HDR
                       + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*("
                       + "*".join(["2^9999"] * 2000) + ")\n}\n")
tname, msg, _cpu = _timed_refusal(_r4mc)
check("T8_multchain_bomb_refused_fast_typed",
      tname == "ValueError" and "CUMULATIVE" in msg and _cpu < 1.0,
      f"18KB 2^9999-x2000 Mult chain: {tname} in {_cpu * 1000:.0f} ms CPU "
      f"(a permissive build COMPLETES this in 6.3 s / 2.0 GiB RSS — the value the caps refuse "
      f"as (2^9999)^2000)")
_r4dc = os.path.join(OUT, "t8_r4_divchain.m")
open(_r4dc, "w").write("{\n" + _T8HDR
                       + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*(1/"
                       + "/".join(["3^9999"] * 400) + ")\n}\n")
tname, msg, _cpu = _timed_refusal(_r4dc)
check("T8_divchain_bomb_refused_fast_typed",
      tname == "ValueError" and "CUMULATIVE" in msg and _cpu < 1.0,
      f"4KB div chain: {tname} in {_cpu * 1000:.0f} ms CPU (an uncapped eval: 12.5 s)")
_r4sp = os.path.join(OUT, "t8_r4_sympow.m")
open(_r4sp, "w").write("{\n" + _T8HDR
                       + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*("
                       + "*".join(["(d+1)^9999"] * 2000) + ")\n}\n")
tname, msg, _cpu = _timed_refusal(_r4sp)
check("T8_sympow_merge_bomb_refused_fast_typed",
      tname == "ValueError" and "ACCUMULATED" in msg and _cpu < 1.0,
      f"22KB (d+1)^9999-x2000 merge (exponent 19,998,000 >> 1e4 cap): "
      f"{tname} in {_cpu * 1000:.0f} ms CPU (uncapped: parse cheap, eval "
      f"TIMEOUT-KILLED at 60 s — reachable via reduce_and_compare)")
_r4af = os.path.join(OUT, "t8_r4_addflood.m")
open(_r4af, "w").write("{\n" + _T8HDR
                       + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*("
                       + "+".join(f"d^{i}" for i in range(1, 5001)) + ")\n}\n")
tname, msg, _cpu = _timed_refusal(_r4af)
check("T8_addflood_refused_fast_typed",
      tname == "ValueError" and "op-units" in msg and _cpu < 2.0,
      f"44KB 5000-term Add flood: {tname} in {_cpu * 1000:.0f} ms CPU "
      f"(uncapped: 48.6 s quadratic build)")
_r4ep = os.path.join(OUT, "t8_r4_eagerpow.m")
open(_r4ep, "w").write("{\n" + _T8HDR
                       + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]"
                         "*((2^9999*d)^9999)\n}\n")
tname, msg, _cpu = _timed_refusal(_r4ep)
check("T8_eagerpow_mulhead_refused_typed",
      tname == "ValueError" and "estimate" in msg and _cpu < 1.0,
      f"20-byte (2^9999*d)^9999: {tname} in {_cpu * 1000:.0f} ms CPU (sympy "
      f"eagerly distributes the power over the Mul's numeric head — 1e8-bit "
      f"integer at construction; found while curing, same class)")
# eval-side soft budget: a merged-exponent expr arriving AT eval refuses
# typed instead of hanging (eval_coeff is a budgeted Fraction walk now)
import sympy as _sp
try:
    okd.eval_coeff((_sp.Symbol("d") + 1) ** _sp.Integer(19998000),
                   F(97, 23), F(5, 7))
    _ok_ev = False
except ValueError:
    _ok_ev = True
except Exception:
    _ok_ev = False
check("T8_eval_softbudget_refused_typed", _ok_ev,
      "merged-exponent expr at eval_coeff refuses ValueError (internal "
      "soft budget — never a 60 s hang; sympy .subs retired from this code path)")
# mutation control: all four budget caps lifted + the exp walk removed in a
# scratch copy -> the mult chain COMPLETES and the merged exponent parses
# clean — the measured DoS band reopens (the refusals are the caps' doing)
_t0m = time.time()
_c5 = make_clone("r4_budgetoff")
_op5 = os.path.join(_c5, "trust", "oracle_k2disp.py")
_otxt5 = open(_op5).read()
_r4anchors = [
    ("_MAX_NUM_BITS = 2 * 10 ** 6      # GUARD NUM-BITS (battery mutation anchor)",
     "_MAX_NUM_BITS = 2 * 10 ** 60      # MUTANT"),
    ("_MAX_TOTAL_BITS = 2 * 10 ** 6    # GUARD TOTAL-BITS (battery mutation anchor)",
     "_MAX_TOTAL_BITS = 2 * 10 ** 60    # MUTANT"),
    ("_MAX_TOTAL_OPS = 10 ** 5         # GUARD TOTAL-OPS (battery mutation anchor)",
     "_MAX_TOTAL_OPS = 10 ** 15         # MUTANT"),
    ("_MAX_AST_NODES = 5 * 10 ** 4     # GUARD AST-NODES (battery mutation anchor)",
     "_MAX_AST_NODES = 10 ** 12     # MUTANT"),
    ("for _pw in ex.atoms(sp.Pow):   # GUARD EXP-WALK (battery mutation anchor)",
     "for _pw in ():   # MUTANT exp walk off"),
]
for _a, _m in _r4anchors:
    assert _a in _otxt5, f"anchor missing: {_a[:40]}"
    _otxt5 = _otxt5.replace(_a, _m)
open(_op5, "w").write(_otxt5)
_r4mx = os.path.join(OUT, "t8_r4_mut_multchain.txt")
open(_r4mx, "w").write("*".join(["2^9999"] * 400))
_r4sx = os.path.join(OUT, "t8_r4_mut_sympow.txt")
open(_r4sx, "w").write("*".join(["(d+1)^9999"] * 50))
T8_BUD_CHILD = os.path.join(OUT, "t8_bud_child.py")
open(T8_BUD_CHILD, "w").write("""\
import json, sys
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
from trust import oracle_k2disp as okd
import sympy as sp
out = {}
try:
    e = okd._safe_rational_expr(open(sys.argv[2]).read().replace("^", "**"))
    out["mult"] = {"refused": False, "is_number": bool(e.is_Number)}
except ValueError:
    out["mult"] = {"refused": True}
try:
    e2 = okd._safe_rational_expr(open(sys.argv[3]).read().replace("^", "**"))
    mx = max([abs(int(p.exp)) for p in e2.atoms(sp.Pow)
              if p.exp.is_Integer] or [0])
    out["sympow"] = {"refused": False, "max_exp": mx}
except ValueError:
    out["sympow"] = {"refused": True}
print(json.dumps(out))
""")
r = child([sys.executable, "-P", T8_BUD_CHILD, _c5, _r4mx, _r4sx], timeout=180)
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"err": (r.stderr or r.stdout)[-200:]}
mutation("T8m_f2r4_budgets_removed", "T8", "copy r4_budgetoff oracle_k2disp.py",
         "all four aggregate caps lifted + merged-exponent walk removed",
         "2^9999-x400 Mult chain + (d+1)^9999-x50 merge", time.time() - _t0m,
         v.get("mult", {}).get("refused") is False
         and v.get("mult", {}).get("is_number") is True
         and v.get("sympow", {}).get("refused") is False
         and v.get("sympow", {}).get("max_exp") == 499950,
         "budget-off copy COMPLETES the 4e6-bit rebuild and parses the "
         "499,950 merged exponent clean — the refusals are the caps' doing")

# --- byte-alphabet gate — lookalike bytes refuse BEFORE parse -
# A record written with U+2192, a fullwidth digit or a lookalike family name
# matched NEITHER coverage counter and was silently dropped (the R3 'whole
# byte stream' contract violated). The raw table must be ASCII-printable +
# newline; anything else refuses typed naming offset/line.
_T8T1B = " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*((d))\n"
_GOOD2 = ("{\n" + _T8HDR + _T8T1 + ",\n"
          + "lbl3m2L_k2disp[0,0,0,0,0,2,1,0,0] -> \n"
          + " + lbl3m2L_k2disp[0,0,0,0,0,1,1,0,0]*(7)\n" + ",\n")
_lookalikes = {
    "arrow": "lbl3m2L_k2disp[3,1,0,0,0,0,1,0,0] → \n" + _T8T1B + ",\n",
    "fullwidth": "lbl3m2L_k2disp[３,1,0,0,0,0,1,0,0] -> \n"
                 + _T8T1B + ",\n",
    "famname": "lbl3m2L_k2disp".replace("l", "ⅼ", 1)
               + "[4,1,0,0,0,0,1,0,0] -> \n" + _T8T1B + ",\n",
}
ok_la = True
_ladet = []
for _tag, _evil in _lookalikes.items():
    _p = os.path.join(OUT, f"t8_r4_lookalike_{_tag}.m")
    with open(_p, "w", newline="") as _fh:
        _fh.write(_GOOD2 + _evil + "}\n")
    tname, msg = _parse_refusal(_p)
    if not (tname == "ValueError" and "offset" in msg and "line" in msg
            and "alphabet" in msg):
        ok_la = False
        _ladet.append(f"{_tag}: {tname} {msg[:60]}")
check("T8_alphabet_lookalikes_refused_named_offset", ok_la,
      "; ".join(_ladet) or
      "U+2192 arrow / fullwidth digit / lookalike famname all refuse typed "
      "naming byte offset + line BEFORE parsing (a permissive parser silently "
      "drops the evil entry from the certified set)")
tname, msg = _parse_refusal(RECEIPT_PINS["ADAPTER_PROBE"][0])
check("T8_alphabet_nonascii_real_file_refused",
      tname == "ValueError" and "alphabet" in msg,
      "a real non-ASCII file (ADAPTER_PROBE.md, em-dashes) refuses at the "
      "byte gate — the old 0-entry junk parse now refuses even earlier")
_p_lf = os.path.join(OUT, "t8_r4_crlf_base.m")
with open(_p_lf, "w", newline="") as _fh:
    _fh.write(_GOOD2 + "}\n")
_p_cr = os.path.join(OUT, "t8_r4_crlf.m")
with open(_p_cr, "w", newline="") as _fh:
    _fh.write((_GOOD2 + "}\n").replace("\n", "\r\n"))
try:
    _e_lf = okd.parse_kira_table(_p_lf)
    _e_cr = okd.parse_kira_table(_p_cr)
    _ok_cr = (len(_e_lf) == 2 and sorted(_e_lf) == sorted(_e_cr))
except Exception:
    _ok_cr = False
check("T8_alphabet_crlf_positive_control", _ok_cr,
      "CRLF table parses identically to its LF twin (CR allowed + "
      "normalized; the R4 probe found no CRLF differential and the gate "
      "must not create one)")
# mutation control: alphabet gate removed in a scratch copy -> the arrow
# table PARSES to 2 entries, the evil third entry silently dropped (the
# exact R4-F3 differential: both coverage counters blind to it)
_t0m = time.time()
_c6 = make_clone("r4_alphaoff")
_op6 = os.path.join(_c6, "trust", "oracle_k2disp.py")
_otxt6 = open(_op6).read()
_aG = ("    _mba = _TABLE_BAD_BYTE.search(raw)   "
       "# GUARD ALPHABET (battery mutation anchor)")
assert _aG in _otxt6, "anchor missing"
open(_op6, "w").write(_otxt6.replace(
    _aG, "    _mba = None   # MUTANT: alphabet gate removed"))
T8_ALPHA_CHILD = os.path.join(OUT, "t8_alpha_child.py")
open(T8_ALPHA_CHILD, "w").write("""\
import json, sys
sys.dont_write_bytecode = True
sys.path.insert(0, sys.argv[1])
from trust import oracle_k2disp as okd
try:
    ents = okd.parse_kira_table(sys.argv[2])
    print(json.dumps({"refused": False, "n_entries": len(ents)}))
except ValueError:
    print(json.dumps({"refused": True}))
""")
r = child([sys.executable, "-P", T8_ALPHA_CHILD, _c6,
           os.path.join(OUT, "t8_r4_lookalike_arrow.m")], timeout=180)
try:
    v = json.loads(r.stdout.strip().splitlines()[-1])
except Exception:
    v = {"err": (r.stderr or r.stdout)[-200:]}
mutation("T8m_f3r4_alphabet_removed", "T8", "copy r4_alphaoff oracle_k2disp.py",
         "byte-alphabet gate removed",
         "U+2192 lookalike-arrow table", time.time() - _t0m,
         v.get("refused") is False and v.get("n_entries") == 2,
         "alphabet-off copy silently certifies the 2-entry SUBSET, the "
         "arrow entry dropped invisibly to both coverage counters — the "
         "refusal is the gate's doing")
leg_wall("T8", t0)
save()

# =========================== OVERALL + COVERAGE =============================
UNWIRED = [
    ("U1 pilot-B emit table-mode re-run",
     "997.5 s measured LOADED (eliminator 576.9+415.1 s of it; unloaded "
     "prior 2.2-2.3 s/(row,prime))",
     "receipts sha-pinned + bars parsed (T4); 84-witness dual-core re-verify "
     "wired live (0.6 s)"),
    ("U2 FULL {1..7} 42-target lp_syz re-run, 2 points",
     "269.2 + 369.5 s measured (638.7 s combined, over the 10-min leg budget)",
     "receipts sha-pinned + bars parsed (T4); P1 sub-tower wired live (4 s)"),
    ("U3 strata full 4-slice 2-prime dictionary certification",
     "reference run 1493-3400 s/slice class (measured band)",
     "identity-linked (T7); strata's own battery certifies its bytes"),
    ("U4 wrong-table MUST-FAIL",
     "reference run 6/6 wrong + 6/6 controls (receipt battery)",
     "linked by identity, not duplicated in TRUST"),
    ("U5 lp_syz_prod execution",
     "8-var top-node syz RIGHT-CENSORED 1800 s hard kill",
     "files pinned (T0/T1); execution not wired into this battery"),
    ("U6 adversarial hardening rounds >=2 + confirmatory",
     "gated before final acceptance",
     "separately gated, before final acceptance"),
    ("U7 m2-symbolic-route witnesses (sec503 class)",
     "winnow reference run: solve 176.7 s, 8/8 witnesses",
     "fixture here is numeric-mass/(d,eta)-rational; scope pin restated"),
]
RESULTS["unwired"] = UNWIRED
n_skips = sum(1 for v in RESULTS["legs"].values() if v.get("pass") is None)
n_checks = len(RESULTS["legs"]) - n_skips
n_mut = sum(1 for m in MUTATIONS if not m.get("expected_blindspot"))
n_mut_ok = sum(1 for m in MUTATIONS
               if m["caught"] and not m.get("expected_blindspot"))
RESULTS["overall"] = "PASS" if not FAILED else f"FAIL {FAILED}"
RESULTS["walls_s"]["TOTAL"] = round(time.time() - T0_START, 1)
save()
print("\n== COVERAGE ==", flush=True)
print(f"wired: {n_checks} checks across T0-T8 + {n_mut} mutation controls "
      f"({n_mut_ok} caught) + 1 recorded blindspot demo (T3m1b, expected) "
      f"+ {n_skips} named skips (reference sources not included).")
print("declared-unwired (measured numbers):")
for name, num, standin in UNWIRED:
    print(f"  {name}: {num} -> {standin}")
print(f"\nper-leg walls (s): " + ", ".join(
    f"{k} {v}" for k, v in RESULTS["walls_s"].items()))
print(f"\nOVERALL: {'ALL PASS' if not FAILED else 'FAIL'} — "
      f"{n_checks} checks, {n_mut_ok}/{n_mut} mutations caught, "
      f"{n_skips} named skips, wall {RESULTS['walls_s']['TOTAL']} s"
      + ("" if not FAILED else f"; failed: {FAILED}"))
sys.exit(0 if not FAILED else 1)

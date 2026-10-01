"""trust._core — integrity (fail-closed vendor pins)
+ the two loading mechanisms (the house assembly pattern for packages that
front vendored engines and live tools):

  load_vendored(name)   — restored lp-syz engines executed via
                          exec(compile(source)) with a per-load source-hash
                          check; the bytecode cache is never consulted and
                          never written (pycache-clobbers-sealed-bytes footgun).
  run_vendored(name,..) — CHILD-MODE runner (the supported compute mode for
                          engine RUNS; children launch -E -P -s -B
                          env-scrubbed — -s means the
                          user-writable user site is never site-processed,
                          so its .pth files never execute and it can never
                          be reordered ahead of the stdlib; the engines'
                          third-party substrate (sympy/python-flint, when
                          installed in the user site) is reached by an
                          explicit TAIL append that executes nothing).
                          In-process load_vendored/driver use is supported
                          ONLY behind the trust/__init__ import-time gates:
                          the env refusal (a raw-PYTHONPATH
                          world refuses typed before any pin is trusted) AND
                          the post-import stdlib-identity check (a stdlib module reached through
                          sys.path[0]/cwd OR a .pth-reordered site dir,
                          neither needing an env var, refuses typed
                          StdlibShadowError; roots are sysconfig stdlib/
                          platstdlib ONLY, never site-packages). The
                          identity check authenticates __file__ — spoofable
                          by an already-hostile interpreter, so it is the
                          honest-import perimeter, not a sandbox. NAMED
                          RESIDUALS: a .pth in the ROOT-owned
                          system site-packages executes arbitrary code at
                          every interpreter startup (children included —
                          root write access required); the host process's
                          own user-site .pth runs before the import gates
                          (same-user write = outside the pin perimeter).
                          Fail-closed output root: refuses to run without
                          TRUST_OUT_ROOT and VALUE-validates it
                          (the house fail-closed output-root rule).
  alias_identity(...)   — live tools (strata, ibplapper, the receipt member)
                          imported IN PLACE by identity from their registered
                          homes — never copied, never shadowed (the alias
                          contract shared with tools/popcorn).
"""
import hashlib
import importlib
import os
import shutil
import subprocess
import sys
import threading
import types

from ._pins import PINS, LINKED, LINKED_SHAS

_PKG = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(_PKG)
_VENDOR = os.path.join(TOOL, "vendor")
_verified = {"ok": False}
_cache = {}
_lock = threading.RLock()

# vendored-engine short names -> vendor-relative paths
VENDORED = {
    "lp_syz": "lpsyz/lp_syz.py",
    "lp_syz_431": "lpsyz/lp_syz_431.py",
    "lp_syz_prod": "lpsyz/lp_syz_prod.py",
    "bessel_oracle2": "lpsyz/bessel_oracle2.py",
}

# Output-root VALUE validation: ALLOWLIST posture — the root must resolve
# under a declared scratch tree; the denylist stays as belt-and-braces and
# carries the package's OWN tree derived from its actual location, so a
# relocated copy still refuses to write into itself (a hardcoded
# canonical prefix did not travel with the package).
# Both lists are env-configurable (colon-separated): TRUST_ALLOWED_OUT_ROOTS
# replaces the allowlist; TRUST_DENY_ROOTS extends the denylist. The generic
# defaults keep the same fail-closed posture on any machine (public-port
# machine-specific prefixes come from env, never the source).
_ALLOWED_OUT_PREFIXES = tuple(
    p for p in os.environ.get("TRUST_ALLOWED_OUT_ROOTS", "").split(":") if p
) or (
    "/tmp",                              # scratch
    "/var/tmp",
)
_FORBIDDEN_OUT_PREFIXES = tuple(
    p for p in os.environ.get("TRUST_DENY_ROOTS", "").split(":") if p
) + (
    "/home",
    os.path.realpath(TOOL),                # this package's own tree, wherever it lives
)

# Env keys scrubbed from every child launch (a
# PYTHONPATH stdlib-shadow — pins-aware hashlib PoC — greened a tampered
# vendor; children are launched -E AND with these keys removed).
_ENV_SCRUB = ("PYTHONPATH", "PYTHONSTARTUP", "PYTHONHOME", "PYTHONEXECUTABLE",
              "PYTHONUSERBASE", "PYTHONCASEOK", "PYTHONNOUSERSITE")


class VendorTamperError(RuntimeError):
    pass


class OutputRootError(RuntimeError):
    pass


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(quiet=False):
    """FAIL-CLOSED on any vendored mismatch/missing/unpinned/symlink/pycache
    content (VendorTamperError); identity-linked live tools
    checked advisory-loud (drift reported, never raises — live tools advance
    independently)."""
    report = {"vendor_ok": [], "vendor_bad": [],
              "linked_missing": [], "linked_drift": []}
    for rel, vsha in PINS.items():
        vp = os.path.join(_VENDOR, rel)
        if not os.path.isfile(vp) or os.path.islink(vp):
            report["vendor_bad"].append((rel, "MISSING-OR-SYMLINK", vsha))
            continue
        try:
            got = _sha(vp)
        except OSError as e:
            report["vendor_bad"].append(
                (rel, f"UNREADABLE({e.__class__.__name__})", vsha))
            continue
        if got == vsha:
            report["vendor_ok"].append(rel)
        else:
            report["vendor_bad"].append((rel, got, vsha))
    # tamper sweep: unpinned files / symlinks / pycache content under vendor/
    for root, dirs, files in os.walk(_VENDOR):
        for d in list(dirs):
            p = os.path.join(root, d)
            if os.path.islink(p):
                report["vendor_bad"].append(
                    (os.path.relpath(p, _VENDOR), "SYMLINKED-DIR", "-"))
                dirs.remove(d)
            elif d == "__pycache__":
                for f2 in os.listdir(p):
                    report["vendor_bad"].append(
                        (os.path.join(os.path.relpath(root, _VENDOR), d, f2),
                         "PYCACHE-CONTENT", "-"))
                dirs.remove(d)
        for f in files:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, _VENDOR)
            if os.path.islink(p):
                report["vendor_bad"].append((rel, "SYMLINK", "-"))
            elif rel not in PINS:
                report["vendor_bad"].append((rel, "UNPINNED", "-"))
    # advisory: identity-linked live homes present + LINKED_SHAS baseline
    # compare — loud, never fatal: live tools advance
    # independently; re-baseline LINKED_SHAS with a recorded reason
    for name, home in LINKED.items():
        if not os.path.isdir(home):
            report["linked_missing"].append((name, home))
            continue
        for rel, base in LINKED_SHAS.get(name, {}).items():
            lp = os.path.join(home, rel)
            try:
                if not os.path.isfile(lp):
                    report["linked_missing"].append((name, lp))
                elif _sha(lp) != base:
                    report["linked_drift"].append((name, rel))
            except OSError:
                report["linked_missing"].append((name, lp + " (unreadable)"))
    if report["vendor_bad"]:
        _verified["ok"] = False
        raise VendorTamperError(
            "TRUST vendor integrity FAILED (tool refuses to compute): "
            + "; ".join(f"{r}: {str(g)[:16]} != pin {str(pn)[:16]}"
                        for r, g, pn in report["vendor_bad"]))
    _verified["ok"] = True
    if not quiet and (report["linked_missing"] or report["linked_drift"]):
        print("TRUST ADVISORY: linked homes differ or missing "
              "(vendored bytes remain the certified copy):")
        for n, h in report["linked_missing"]:
            print(f"  missing linked home: {n} at {h}")
        for n, rel in report["linked_drift"]:
            print(f"  LINKED DRIFT (advisory): {n}/{rel} no longer matches the "
                  f"recorded baseline — the identity front tracks a moved "
                  f"target; re-baseline _pins.LINKED_SHAS with a recorded "
                  f"reason if the move is legitimate")
    return report


def require_verified():
    if not _verified["ok"]:
        verify(quiet=True)


def _purge_pycache():
    for root, dirs, _files in os.walk(_VENDOR):
        for d in list(dirs):
            if d == "__pycache__":
                shutil.rmtree(os.path.join(root, d), ignore_errors=True)
                dirs.remove(d)


def vendored_path(name):
    return os.path.join(_VENDOR, VENDORED[name])


def load_vendored(name):
    """Load a vendored engine by short name: source-hash-checked at load time,
    exec(compile(source)); bytecode cache never an input. The engines are
    self-contained scripts (no intra-vendor imports) — no bare-name seeding.
    NOTE: engines guard their CLI under __main__, so exec is import-safe."""
    with _lock:
        if name in _cache:
            mod = _cache[name]
            sys.modules[f"trust.vendor.{name}"] = mod
            return mod
        rel = VENDORED[name]
        require_verified()
        _purge_pycache()
        path = os.path.join(_VENDOR, rel)
        if not os.path.isfile(path) or os.path.islink(path):
            raise VendorTamperError(
                f"engine {name}: vendored path is not a regular file — refusing")
        src = open(path, "rb").read()
        got = hashlib.sha256(src).hexdigest()
        if got != PINS[rel]:
            raise VendorTamperError(
                f"engine {name}: source bytes sha {got[:16]} != pin "
                f"{PINS[rel][:16]} at load time")
        old_dwb = sys.dont_write_bytecode
        sys.dont_write_bytecode = True
        try:
            mod = types.ModuleType(f"trust.vendor.{name}")
            mod.__file__ = path
            sys.modules[mod.__name__] = mod
            pre_path = list(sys.path)
            try:
                exec(compile(src, path, "exec"), mod.__dict__)
            finally:
                sys.path[:] = pre_path
        finally:
            sys.dont_write_bytecode = old_dwb
        _cache[name] = mod
        return mod


def load_vendored_instance(name, tag):
    """Load a PRIVATE instance of a vendored engine (sha-checked, exec'd
    fresh, NOT the shared cache) so a front may configure family globals
    (build_G/GRAM/PROPS) without mutating the module every other consumer
    sees. Registered as trust.vendor.<name>__<tag>."""
    with _lock:
        rel = VENDORED[name]
        require_verified()
        _purge_pycache()
        path = os.path.join(_VENDOR, rel)
        src = open(path, "rb").read()
        got = hashlib.sha256(src).hexdigest()
        if got != PINS[rel]:
            raise VendorTamperError(
                f"engine {name}: source bytes sha {got[:16]} != pin at load time")
        old_dwb = sys.dont_write_bytecode
        sys.dont_write_bytecode = True
        try:
            mod = types.ModuleType(f"trust.vendor.{name}__{tag}")
            mod.__file__ = path
            sys.modules[mod.__name__] = mod
            pre_path = list(sys.path)
            try:
                exec(compile(src, path, "exec"), mod.__dict__)
            finally:
                sys.path[:] = pre_path
        finally:
            sys.dont_write_bytecode = old_dwb
        return mod


def out_root():
    """Fail-closed output root: TRUST_OUT_ROOT env, VALUE-validated."""
    root = os.environ.get("TRUST_OUT_ROOT")
    if not root:
        raise OutputRootError(
            "TRUST_OUT_ROOT is not set — trust refuses env-less compute "
            "(fail-closed output-root law). Point it at a scratch run directory.")
    if not os.path.isabs(root):
        raise OutputRootError(f"TRUST_OUT_ROOT={root!r} is not absolute — refusing")
    rp = os.path.realpath(root)
    for pref in _FORBIDDEN_OUT_PREFIXES:
        if rp == pref or rp.startswith(pref + os.sep):
            raise OutputRootError(
                f"TRUST_OUT_ROOT={root!r} resolves under forbidden prefix "
                f"{pref!r} (protected tree: source/parked/tool/published) — refusing")
    if not any(rp == a or rp.startswith(a + os.sep)
               for a in _ALLOWED_OUT_PREFIXES):
        raise OutputRootError(
            f"TRUST_OUT_ROOT={root!r} (realpath {rp!r}) is outside every "
            f"allowed scratch prefix {_ALLOWED_OUT_PREFIXES} — allowlist "
            f"posture: point it at a permitted scratch root "
            f"(TRUST_ALLOWED_OUT_ROOTS extends the allowlist)")
    os.makedirs(rp, exist_ok=True)
    return rp


# child bootstrap: launched as `python3 -E -P -s -B -c BOOT
# <engine path> <args...>`. Under -s the user site is NEVER site-processed —
# its .pth files (arbitrary code at interpreter startup, executed BEFORE any
# script line; the site-reorder vector) never run. The engines'
# third-party substrate (sympy/python-flint) is user-site-installed on this
# box, so the bootstrap appends that ONE directory at the TAIL of sys.path:
# the stdlib and the system site keep precedence (a planted hashlib.py there
# stays inert) and the append executes nothing. runpy runs the engine with
# __name__ == "__main__" and __file__ = the engine path — script semantics,
# minus the script-dir sys.path prepend (-P holds under -c too).
_CHILD_BOOT = (
    "import sys, site, runpy\n"
    "p = sys.argv.pop(1)\n"
    "sys.argv[0] = p\n"
    "u = site.getusersitepackages()\n"
    "if u and u not in sys.path:\n"
    "    sys.path.append(u)\n"
    "runpy.run_path(p, run_name='__main__')\n"
)


def run_vendored(name, args, timeout=None, cwd=None, env_extra=None):
    """CHILD-MODE runner: verify pins, then run the vendored engine as
    `python3 -E -P -s -B -c <bootstrap> <vendored path> ...` with cwd under
    the validated TRUST output root. -P: no script-dir/cwd sys.path prepend
    -E + _ENV_SCRUB: env-derived interpreter config
    (PYTHONPATH / PYTHONHOME / PYTHONSTARTUP) neutralized — the stdlib-shadow
    class (TRUST BLOCKING). -s (F1): the user-writable user
    site is never site-processed — its .pth files never execute; sympy/flint
    access rides the bootstrap's inert TAIL append (see _CHILD_BOOT). NOTE
    the earlier rationale for omitting -E was WRONG: the user site comes from
    site.py, not PYTHONPATH, so sympy/python-flint import fine under -E
    (verified live). -B: no bytecode writes (with -E the
    PYTHONDONTWRITEBYTECODE env var is ignored). Residual (named): .pth files
    in the ROOT-owned system site-packages still execute at child startup.
    Returns subprocess.CompletedProcess."""
    require_verified()
    _purge_pycache()
    path = vendored_path(name)
    got = _sha(path)
    if got != PINS[VENDORED[name]]:
        raise VendorTamperError(f"engine {name}: sha changed between verify and launch")
    root = out_root()
    rundir = cwd or root
    rr, rt = os.path.realpath(rundir), os.path.realpath(root)
    if rr != rt and not rr.startswith(rt + os.sep):
        raise OutputRootError(f"run cwd {rundir!r} escapes TRUST_OUT_ROOT — refusing")
    os.makedirs(rundir, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if k not in _ENV_SCRUB}
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if env_extra:
        env.update(env_extra)
    return subprocess.run([sys.executable, "-E", "-P", "-s", "-B", "-c",
                          _CHILD_BOOT, path] + list(args),
                          cwd=rundir, env=env, capture_output=True,
                          text=True, timeout=timeout)


# sys.path ROOTS that contain the importable top-level packages (identity
# imports only — inserted transiently, identity preserved via sys.modules).
# strata ships as the member tree beside this package (TOOL/strata) and the
# front resolves there by default; TRUST_STRATA_ROOT overrides. Other
# linked-tool roots come from the environment (alias_identity refuses
# unregistered/unset names fail-closed).
IMPORT_ROOTS = {
    k: v for k, v in {
        "strata": (os.environ.get("TRUST_STRATA_ROOT", "")
                   or (TOOL if os.path.isdir(os.path.join(TOOL, "strata"))
                       else "")),
        "ibplapper": os.environ.get("TRUST_IBPLAPPER_ROOT", ""),
    }.items() if v
}


def alias_identity(name, root=None):
    """Import module/package `name` from sys.path root `root` IN PLACE by
    identity — shared module object with every legacy importer; refuses a
    shadowing same-named module from elsewhere (the house alias_registered
    contract, incl. spec authentication)."""
    with _lock:
        top = name.split(".")[0]
        # membership BEFORE abspath — os.path.abspath("")
        # returns cwd, which would make the typed refusal unreachable and
        # fall open to cwd imports for unregistered names.
        reg = root if root is not None else IMPORT_ROOTS.get(top)
        if not reg:
            raise ImportError(
                f"no registered import root for {name!r} — refusing "
                f"(fail-closed: unregistered names never import from cwd)")
        root = os.path.abspath(reg)
        mod = sys.modules.get(name)
        pre_existing = mod is not None
        if mod is None:
            added = root not in sys.path
            if added:
                sys.path.insert(0, root)
            try:
                mod = importlib.import_module(name)
            finally:
                if added and root in sys.path:
                    sys.path.remove(root)
        f = os.path.abspath(getattr(mod, "__file__", "") or "")
        pkg_dir = os.path.join(root, top)
        ok_loc = (f.startswith(pkg_dir + os.sep)
                  or f == pkg_dir + ".py"
                  or os.path.dirname(f) == pkg_dir)
        if not ok_loc or not os.path.isfile(f):
            raise ImportError(
                f"trust identity violation: module {name!r} resolved to "
                f"{f or None!r}, not inside {pkg_dir} — a same-named module "
                f"is shadowing the registered tool (or a namespace-package "
                f"fallthrough); refusing to alias it.")
        if pre_existing:
            spec = getattr(mod, "__spec__", None)
            origin = os.path.abspath(getattr(spec, "origin", "") or "")
            if spec is None or origin != f:
                raise ImportError(
                    f"trust identity violation: pre-existing module {name!r} "
                    f"has no import spec matching its __file__ — provenance "
                    f"unverifiable; refusing to alias it.")
        return mod

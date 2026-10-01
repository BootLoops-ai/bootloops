"""baller._core — integrity (fail-closed vendor pins)
+ the two loading mechanisms:

  load_vendored(name)  — vendored engines executed via exec(compile(source))
                         with a per-load source-hash check; the bytecode cache
                         is never consulted and never written.
  alias_registered(..) — registered tools (posq core, wayfinder
                         gate/manifest) imported IN PLACE by identity from
                         their own homes (their receipts pin their paths;
                         BALLER fronts, never moves).
"""
import hashlib
import importlib
import os
import shutil
import sys
import threading
import types

from ._pins import PINS

_PKG = os.path.dirname(os.path.abspath(__file__))
_VENDOR = os.path.join(os.path.dirname(_PKG), "vendor")
# Registered-tools root: the aliased tools (posq, wayfinder) ship BESIDE
# this package in the repository's tools/ directory, so
# the root defaults to this package's grandparent. A checkout that keeps the
# aliased tools elsewhere sets BALLER_TOOLS_ROOT.
TOOLS = os.environ.get("BALLER_TOOLS_ROOT") or os.path.dirname(os.path.dirname(_PKG))
_verified = {"ok": False}
_cache = {}
_lock = threading.RLock()   # concurrent first-loads would shatter engine
# identity (distinct module objects per thread, incoherent cache) —
# every load/alias runs under this lock

VENDORED = {
    "pipe_transport": "kklt/pipe_transport.py",
    "pipe_vac": "kklt/pipe_vac.py",
    "pipe_lib": "kklt/pipe_lib.py",
    "march_lib": "kklt/march_lib.py",
    "mc": "geo/mc.py",
}

# intra-vendor bare-name imports (grep-verified at vendor time): deps preload
# through THIS loader and register under the bare name, so `import pipe_lib`
# inside a vendored engine resolves to the pinned bytes — never to sys.path.
_PRELOAD = {
    "pipe_transport": ["pipe_lib"],
    "pipe_vac": ["pipe_lib", "pipe_transport"],
}
# ONLY names that vendored engines actually import bare get a sys.modules
# bare entry — never process-wide aliases (registering every engine would
# hijack hyper-generic names: a user's own ./mc.py could silently bind the
# vendored geo mc or vice versa)
_NEEDS_BARE = {"pipe_lib", "pipe_transport"}


class VendorTamperError(RuntimeError):
    pass


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(quiet=False):
    """FAIL-CLOSED on any vendored mismatch/missing file (VendorTamperError)."""
    report = {"vendor_ok": [], "vendor_bad": [], "vendor_unpinned": []}
    for rel, vsha in PINS.items():
        vp = os.path.join(_VENDOR, rel)
        if not os.path.isfile(vp):
            report["vendor_bad"].append((rel, "MISSING", vsha))
            continue
        try:
            got = _sha(vp)
        except OSError as e:               # typed on the vendor side too
            report["vendor_bad"].append((rel, f"UNREADABLE({e.__class__.__name__})", vsha))
            continue
        if got != vsha:
            report["vendor_bad"].append((rel, got, vsha))
        else:
            report["vendor_ok"].append(rel)
    # ANY unpinned file under vendor/ is a tamper class (never .py-only —
    # planted .pyc/.so/.pth/case-variants must trip the sweep), and SYMLINKS
    # anywhere under vendor/ break self-containment (a dir symlink would hide
    # payloads from the walk entirely).
    for root, dirs, files in os.walk(_VENDOR):
        for d in list(dirs):
            p = os.path.join(root, d)
            if os.path.islink(p):
                report["vendor_bad"].append(
                    (os.path.relpath(p, _VENDOR), "SYMLINKED-DIR", "-"))
                dirs.remove(d)
            elif d == "__pycache__":
                dirs.remove(d)            # purged at every load, never read
        for f in files:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, _VENDOR)
            if os.path.islink(p):
                report["vendor_bad"].append((rel, "SYMLINK", "-"))
            elif rel not in PINS:
                report["vendor_bad"].append((rel, "UNPINNED", "-"))
    if report["vendor_bad"]:
        _verified["ok"] = False
        raise VendorTamperError(
            "BALLER vendor integrity FAILED (tool refuses to compute): "
            + "; ".join(f"{r}: {g[:16] if g != 'MISSING' else g} != pin {p[:16]}"
                        for r, g, p in report["vendor_bad"]))
    _verified["ok"] = True
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


def load_vendored(name):
    """Load a vendored engine by short name: source-hash-checked at load time,
    exec(compile(source)) — the bytecode cache is never an input. Thread-safe
    (one lock for all loads); NOTE (documented semantics): an engine already
    loaded in-process keeps being served from the in-memory cache even if a
    LATER on-disk verify() fails — its bytes were verified at its own load."""
    with _lock:
        if name in _cache:
            mod = _cache[name]
            # cache-bypass shadow guard: the fast path must not return
            # without restoring sys.modules — post-first-load pollution of
            # sys.modules[name] got wired into later engines' bare imports.
            # The verified cached module is authoritative: restore the dotted
            # name always, the bare name only where engines need it, and
            # say so when a foreign occupant is displaced (the miss path
            # refuses typed; a silent clobber here would invert that contract).
            sys.modules[f"baller.vendor.{name}"] = mod
            if name in _NEEDS_BARE:
                prior = sys.modules.get(name)
                if prior is not None and prior is not mod:
                    print(f"BALLER ADVISORY: displacing foreign sys.modules"
                          f"[{name!r}] ({getattr(prior, '__file__', '?')}) "
                          f"with the pinned vendored engine", file=sys.stderr)
                sys.modules[name] = mod
            return mod
        rel = VENDORED[name]
        require_verified()
        _purge_pycache()
        path = os.path.join(_VENDOR, rel)
        if not os.path.isfile(path) or os.path.islink(path):
            # the path must be a regular file: a FIFO planted there would
            # make bare open() HANG
            # (blocking is not an OSError); isfile is False for FIFOs and
            # symlinked stand-ins — refuse typed before touching it
            raise VendorTamperError(
                f"engine {name}: vendored path is not a regular file "
                f"(missing/replaced/special) — refusing")
        try:
            src = open(path, "rb").read()
        except OSError as e:                       # typed refusal
            raise VendorTamperError(
                f"engine {name}: vendored file MISSING/unreadable at load "
                f"time ({e}) — incomplete or tampered tree") from e
        got = hashlib.sha256(src).hexdigest()
        if got != PINS[rel]:
            raise VendorTamperError(
                f"engine {name}: source bytes sha {got[:16]} != pin "
                f"{PINS[rel][:16]} at load time")
        old_dwb = sys.dont_write_bytecode
        sys.dont_write_bytecode = True
        try:
            for dep in _PRELOAD.get(name, []):
                load_vendored(dep)
            prior = sys.modules.get(name)
            if prior is not None and getattr(prior, "__name__", "") != f"baller.vendor.{name}":
                raise VendorTamperError(
                    f"engine {name}: a foreign module already occupies sys.modules[{name!r}] "
                    f"({getattr(prior, '__file__', '?')}) — refusing to shadow or be shadowed.")
            mod = types.ModuleType(f"baller.vendor.{name}")
            mod.__file__ = path
            sys.modules[mod.__name__] = mod
            if name in _NEEDS_BARE:       # bare-name capture ONLY where
                sys.modules[name] = mod   # engines import it bare
            code = compile(src, path, "exec")   # resolve here, never via sys.path
            exec(code, mod.__dict__)
            # an import hook firing DURING exec (any uncached import
            # in the engine body) could swap a dep out of sys.modules between
            # preload and the engine's own bare import — verify dep identity
            # AFTER exec, both in sys.modules and as actually bound in the
            # engine's namespace.
            for dep in _PRELOAD.get(name, []):
                want = _cache.get(dep)
                if sys.modules.get(dep) is not want:
                    raise VendorTamperError(
                        f"engine {name}: sys.modules[{dep!r}] was swapped "
                        f"during exec (import-hook interference) — refusing")
                # the engine's actually-BOUND alias must BE the pinned object
                for v in mod.__dict__.values():
                    if (isinstance(v, types.ModuleType)
                            and getattr(v, "__name__", "")
                            in (dep, f"baller.vendor.{dep}")
                            and v is not want):
                        raise VendorTamperError(
                            f"engine {name}: bound dep {dep!r} is not the "
                            f"pinned module (swapped mid-exec) — refusing")
        finally:
            sys.dont_write_bytecode = old_dwb
        _cache[name] = mod
        return mod


def alias_registered(name, home):
    """Import module `name` from directory `home` IN PLACE by identity —
    shared module object with every legacy importer; refuses a shadowing
    same-named module from elsewhere (same loader contract as tools/popcorn).

    Authentication (__file__ alone is spoofable by one attribute):
    a pre-existing sys.modules entry must carry a real import __spec__ whose
    origin agrees with its __file__, and that file must exist on disk in
    `home`. THREAT BOUNDARY (documented in MANUAL): an in-process actor who
    controls sys.modules can forge a full ModuleSpec too — no in-process
    check can beat an attacker who can equally monkeypatch baller itself;
    this check catches every accidental shadow and naive plant."""
    with _lock:
        home = os.path.abspath(home)
        mod = sys.modules.get(name)
        pre_existing = mod is not None
        if mod is None:
            # sys.path is restored after the import — a permanent insert
            # at [0] would shadow user modules process-wide
            # (identity is preserved via sys.modules, not via sys.path)
            added = home not in sys.path
            if added:
                sys.path.insert(0, home)
            try:
                mod = importlib.import_module(name)
            finally:
                if added and home in sys.path:
                    sys.path.remove(home)
        f = os.path.abspath(getattr(mod, "__file__", "") or "")
        if os.path.dirname(f) != home or not os.path.isfile(f):
            raise ImportError(
                f"baller identity violation: module {name!r} resolved to {f}, "
                f"not a file in {home} — a same-named module is shadowing the "
                f"registered tool; refusing to alias it.")
        if pre_existing:
            spec = getattr(mod, "__spec__", None)
            origin = os.path.abspath(getattr(spec, "origin", "") or "")
            if spec is None or origin != f:
                raise ImportError(
                    f"baller identity violation: pre-existing module {name!r} "
                    f"has no import spec matching its __file__ (spec origin "
                    f"{origin or None!r} vs {f!r}) — provenance unverifiable; "
                    f"refusing to alias it.")
        return mod

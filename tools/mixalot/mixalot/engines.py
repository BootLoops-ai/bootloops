"""mixalot.engines — the vendored evidence engines, importable by name.

SECURITY MODEL: the pins guard SOURCE
bytes, so the loader must never trust Python's bytecode cache — a poisoned
same-size .pyc in vendor/__pycache__ would otherwise execute while verify()
reports clean. Therefore load():
  1. deletes any __pycache__ under vendor/ (idempotent, every call),
  2. re-hashes the engine's SOURCE bytes against _pins.PINS (fail-closed,
     per-load — not just at verify() time),
  3. executes exec(compile(source_bytes)) — the cache is never consulted,
  4. holds sys.dont_write_bytecode=True during execution so sibling imports
     cannot mint new cache files,
  5. pre-loads vendored siblings under their PLAIN module names so an
     engine's own `import sibling` (after its hardcoded sys.path insert)
     resolves to the vendored, source-verified module via sys.modules —
     never via the path system.
Vendor bytes are never edited; wrappers parameterize cwd/output paths.
"""
import hashlib
import os
import shutil
import sys
import types

from ._pins import PINS
from ._verify import VendorTamperError, require_verified

_VENDOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor")
_DIRS = [os.path.join(_VENDOR, d)
         for d in ("pilot", "jeff", "census", "blend", "tools", "eco/pilot")]
_cache = {}

_ENGINES = {
    "telegraph_evaluate": "scrna/telegraph_evaluate.py",
    "etienne_evaluate": "eco/etienne_evaluate.py",
    "lsx_direct": "pilot/lsx_direct.py",
    "zseries": "pilot/zseries.py",
    "closed_form_1var": "pilot/closed_form_1var.py",
    "formula_emitter": "pilot/formula_emitter.py",
    "formula_emitter_dirichlet": "pilot/formula_emitter_dirichlet.py",
    "lsx55_exact": "pilot/lsx55_exact.py",
    "bigg": "jeff/bigg.py",
    "bigk2": "jeff/bigk2.py",
    "biggen": "jeff/biggen.py",
    "f_nu": "jeff/f_nu.py",
    "w1_collapsed": "jeff/w1_collapsed.py",
    "annihilator": "tools/annihilator.py",
    "estimators": "jeff/estimators.py",
    "w1_brute": "jeff/w1_brute.py",
    "bigk": "jeff/bigk.py",
    "w4_blind_gf": "jeff/w4_blind_gf.py",
    "w4_dpm_limit": "jeff/w4_dpm_limit.py",
    "swap_route": "jeff/swap_route.py",
    "seg_v1": "census/seg_v1.py",
    "seg_blind": "census/seg_blind.py",
    "mic_power": "census/mic_power.py",
    "production_sweep": "census/production_sweep.py",
    "frozen_comp_v1": "blend/frozen_comp_v1.py",
    "frozen_comp_blind": "blend/frozen_comp_blind.py",
    # Etienne pilot engine cluster member (eco/pilot/). LOADABLE set only —
    # gate_rf / gate_worked_example / verify_* are pinned-but-not-loadable
    # (top-level scripts; gate_rf reads data by relative path at import).
    # gate_engines RUNS its G1/G2 gate at load (~0.2 s, print-only).
    "ball_engine": "eco/pilot/ball_engine.py",
    "etienne_oracle": "eco/pilot/etienne_oracle.py",
    "phase2_certify": "eco/pilot/phase2_certify.py",
    "phase2_multisample_eqI": "eco/pilot/phase2_multisample_eqI.py",
    "multisample": "eco/pilot/multisample.py",
    "multisample_ball": "eco/pilot/multisample_ball.py",
    "hier_kron": "eco/pilot/hier_kron.py",
    "hierarchical3": "eco/pilot/hierarchical3.py",
    "hier_ball": "eco/pilot/hier_ball.py",
    "hier_eval": "eco/pilot/hier_eval.py",
    "rf_engine": "eco/pilot/rf_engine.py",
    "gate_engines": "eco/pilot/gate_engines.py",
}

_PRELOAD = {
    "mic_power": ["seg_v1", "production_sweep"],
    "production_sweep": ["seg_v1"],
    "bigk2": ["bigk"],
    "w1_collapsed": ["w1_brute"],
    "zseries": ["lsx_direct"],
    "ball_engine": ["etienne_oracle"],
    "phase2_certify": ["ball_engine"],
    "phase2_multisample_eqI": ["multisample_ball", "phase2_certify"],
    "rf_engine": ["ball_engine"],
    "hier_kron": ["hierarchical3"],
    "hier_ball": ["hierarchical3"],
    "hier_eval": ["hierarchical3"],
    "gate_engines": ["etienne_oracle", "ball_engine"],
}


def purge_pycache():
    """Delete every __pycache__ under vendor/ (the cache is never
    trusted and never allowed to accumulate)."""
    n = 0
    for root, dirs, _files in os.walk(_VENDOR):
        for d in list(dirs):
            if d == "__pycache__":
                shutil.rmtree(os.path.join(root, d), ignore_errors=True)
                dirs.remove(d)
                n += 1
    return n


def load(name):
    """Load a vendored engine: source-hash-checked against the pins at load
    time, executed via exec(compile(source)) — the bytecode cache is never
    consulted and never written."""
    if name in _cache:
        return _cache[name]
    if name not in _ENGINES:
        raise KeyError(f"unknown engine {name!r}; have {sorted(_ENGINES)}")
    require_verified()
    purge_pycache()
    rel = _ENGINES[name]
    path = os.path.join(_VENDOR, rel)
    src = open(path, "rb").read()
    pinned = PINS[rel][1]
    got = hashlib.sha256(src).hexdigest()
    if pinned is None or got != pinned:
        raise VendorTamperError(
            f"engine {name}: source bytes sha {got[:16]} != pin "
            f"{(pinned or 'UNSTAMPED')[:16]} at load time")
    old_path = list(sys.path)
    old_dwb = sys.dont_write_bytecode
    sys.path[:0] = _DIRS
    sys.dont_write_bytecode = True
    try:
        for dep in _PRELOAD.get(name, []):
            load(dep)
        mod = types.ModuleType(f"mixalot.vendor.{name}")
        mod.__file__ = path
        sys.modules[mod.__name__] = mod
        sys.modules[name] = mod          # plain-name cache beats any sys.path
        code = compile(src, path, "exec")
        exec(code, mod.__dict__)
    finally:
        sys.path[:] = old_path
        sys.dont_write_bytecode = old_dwb
    _cache[name] = mod
    return mod


def __getattr__(name):
    if name in _ENGINES:
        return load(name)
    raise AttributeError(name)

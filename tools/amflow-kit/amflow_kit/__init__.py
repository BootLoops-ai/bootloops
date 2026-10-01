"""amflow_kit — the importable half of amflow-kit, the house kit around the
AMFlow.cpp fork (the sibling repository amflow-cpp).

The kit directory tools/amflow-kit/ holds three command-line output gates
(amflow_output_lint.py, amflow_smoke.sh, compare_amflow_json.py) and this
package, whose two modules are libraries with their own command lines:

  amflow_kit.keypred   offline predictor of the fork's IBP-cache key, hashing
                       the canonical jobs.yaml form the fork itself hashes
                       (python3 -m amflow_kit.keypred key DIR, --selftest)
  amflow_kit.memfence  the memory-fence policy and post-spawn environment
                       readback shared by every launcher of amflow_cli
                       (pmflow, numkin's shift_opt, your own)

Import from a sibling tools/<pkg>/ module with the sys.path hop
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "amflow-kit")); from amflow_kit import memfence
which resolves to tools/amflow-kit/amflow_kit/memfence.py in both the repo and
a mirror tree laid out as tools/<pkg>/. Stdlib only. See ../GUIDE.md,
../KEYPRED.md and ../MEMFENCE.md.
"""
import importlib

__all__ = ["keypred", "memfence"]
__version__ = "1.0.0"


def __getattr__(name):
    # members load on first touch (PEP 562), so `python3 -m amflow_kit.keypred`
    # runs the module once, as __main__, without a prior import of it here
    if name in __all__:
        mod = importlib.import_module("." + name, __name__)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

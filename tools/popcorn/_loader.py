"""popcorn._loader — import-in-place by identity.

The engine modules are plain top-level modules (bare-name imports of each
other), shipped flat beside this file. Each popcorn.* submodule aliases the
module objects obtained by importing them in place (sys.path insert +
bare-name import). Because Python caches by name in sys.modules, POPCORN and
any direct importer share the SAME module object — one source of truth, zero
duplicated math.

The identity assert below is the guard that makes aliasing honest: if a
same-named module from some OTHER tree is already loaded (name collisions on
generic names like certquad do happen), POPCORN refuses to alias it rather
than silently re-exporting foreign code.
"""
import importlib
import sys
from pathlib import Path

from ._pins import PKG_DIR


def import_in_place(name, home=PKG_DIR):
    """Import module `name` from directory `home` (Path), asserting identity.

    Returns the module object shared with every other importer of `name` in
    this process. Raises ImportError if the resolved file does not live in
    `home` (a shadowing module got there first — refuse, never alias)."""
    home = Path(home).resolve()
    mod = sys.modules.get(name)
    if mod is None:
        p = str(home)
        if p not in sys.path:
            sys.path.insert(0, p)
        mod = importlib.import_module(name)
    f = Path(getattr(mod, "__file__", "") or "").resolve()
    if f.parent != home:
        raise ImportError(
            f"popcorn identity violation: module {name!r} resolved to {f}, not "
            f"{home} — a same-named module is shadowing the shipped original; "
            f"refusing to alias it (fold-by-identity would be dishonest).")
    return mod

"""lockpick.ring15 — the ring-basis member: conductor-15 / Q(sqrt(-15))
Chowla-Selberg CM constant ring (27 members, dps 210+120 cross-checked).

THIS member is lockpick's live reach (pslq_gate --selftest asserts member
integrity vs the PINS.json sha table).

Submodules:

  ring15_core   — characters chi_{-3,-4,5,-15}, Dirichlet L (mpmath),
                  Lerch/Chowla-Selberg periods Omega_D, Dedekind eta,
                  CM periods, disc -4 lemniscatic controls.
  ring15_hecke  — exact integer a_n of the two weight-3 level-15 CM
                  newforms (Grossencharacter; 3-way eta-product/PARI
                  validation).
  ring15_lfun   — L(f,s), s=1,2,3, both forms: split Hecke integrals,
                  Fricke matrix fit, t-invariance check.
  ring15        — orchestrator (builds RING15.json at two dps;
                  rebuild runs write under RING15_OUT, default cwd —
                  launch from a scratch dir).

Dependency note: Eichler.jl (shipped in this repo at upgrades/Eichler.jl)
exports its own dirichlet_L / cusp_dictionary; ring15_core.dirichlet_L here
is an independent mpmath implementation, not a copy of Eichler.jl's.
"""

__all__ = ["ring15", "ring15_core", "ring15_hecke", "ring15_lfun"]


def __getattr__(name):
    if name in __all__:
        import importlib
        return importlib.import_module("." + name, __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(__all__))

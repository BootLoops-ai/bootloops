"""popcorn.dominance — certified selection-with-dominance evaluators (h != 1/2).

E(S,h) is an entire q-deformation of the Kummer line (heat relation).
Aliases, in place and by identity (see popcorn._loader):
  dominance_oracle   — certified 2-fold oracle (R_certified, E_dominant);
                       rc-coded selftest (battery --full).
  dominance_qseries  — stable q-series evaluator (E_dominant_qseries);
                       certify by a two-dps agreement gate (the battery
                       checks it against the pinned reference values).

Self-consistency radii are self-consistency, not proven balls — any document
quoting dominance/mixture enclosures must carry that distinction.
"""
from ._loader import import_in_place

dominance_oracle = import_in_place("dominance_oracle")
dominance_qseries = import_in_place("dominance_qseries")

E_dominant = dominance_oracle.E_dominant
E_dominant_qseries = dominance_qseries.E_dominant_qseries

__all__ = ["dominance_oracle", "dominance_qseries",
           "E_dominant", "E_dominant_qseries"]

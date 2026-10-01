"""popcorn.enclosure — certified enclosures of the TRANSIENT selected SFS.

Aliases, in place and by identity (see popcorn._loader), loaded lazily so
`import popcorn.enclosure` needs nothing; the first attribute access loads
the engine, which needs python-flint (Arb ball arithmetic):
  transient_enclosure_engine — the expected unfolded sample SFS E_n(i; t),
        i = 1..n-1 (theta = 1), under genic selection S <= 0 and a
        piecewise-constant size history [(rho, T), ...] after an ancestral
        Wright equilibrium, as arb balls: the moment system of the forward
        diffusion truncated at M >= n with a proved two-sided tail closure
        (LOWER w_{M+1} := 0 and UPPER w_{M+1} := w_M bracket the true
        moments because the truncated matrix is Metzler), exact fmpq
        matrices, rigorous arb_mat exponential + solve per epoch, exact
        integer projection to sample size n. certified_cell(M, S, n, prec,
        epochs) -> report dict with 'entry_balls', 'enclosure', 'gates';
        equilibrium_entry / weq_ball (certified stationary objects);
        float_cell (FLOAT register, scipy) and xcheck (FLOAT-HP register,
        mpmath) as implementation checks; moments_cell (optional `moments`
        cross-check, imported only there); shipped reference cells under
        reference/enclosure/ (load_reference_cell, entry_ball).

REGISTER: CERTIFIED-ENCLOSURE, MODEL-CONDITIONAL. Every entry is an interval
guaranteed to contain the value of the moment system of the stated
diffusion (genic selection, S <= 0, piecewise-constant N(t), theta = 1):
the two truncation brackets enclose the untruncated moments by a proved
lemma and the reported ball is their hull, so truncation is inside the
width (which self-reports it) and is independently controlled by comparing
two truncations (the shipped M = 200 and M = 500 cells agree to all 30
recorded digits). It is not a bound on model error, it does not cover
S > 0, dominance or linkage, and the float/mpmath/moments routes in the
same engine are implementation checks that must never be quoted as
certified. Use it to certify a float transient solver
(popcorn.transient) at spot cells, exactly as the certified stationary
engine (popcorn.sfs) certifies float stationary paths. Conventions:
S = 4 N_ref s (dadi/moments gamma = S/2), time in 2 N_ref generations,
rho = N/N_ref, theta = 1; S integer, epoch entries flint.fmpq (epochs_fmpq
below converts ints, Fractions and 'p/q' strings). certified_cell sets the
process-wide flint precision ctx.prec and leaves it set.
"""
from ._loader import import_in_place

_ENGINE = "transient_enclosure_engine"
_NAMES = ("certified_cell", "equilibrium_entry", "weq_ball", "project_entry",
          "rat_rows", "to_arb_mat", "bvec", "qarb", "float_cell", "xcheck",
          "moments_cell", "run_cell", "EPOCHS", "THETA", "ANC_RHO",
          "REFERENCE_DIR", "REFERENCE_CELLS", "load_reference_cell",
          "load_reference_xcheck", "entry_ball")

__all__ = [_ENGINE, "epochs_fmpq"] + list(_NAMES)


def epochs_fmpq(pairs):
    """[(rho, T), ...] with int / Fraction / 'p/q' string / fmpq entries ->
    the [(fmpq, fmpq), ...] history certified_cell expects (exact rationals;
    a binary float is refused rather than silently rationalized)."""
    from fractions import Fraction
    from flint import fmpq

    def q(v):
        if isinstance(v, fmpq):
            return v
        if isinstance(v, float):
            raise TypeError("epochs_fmpq: pass exact rationals (int, Fraction, "
                            "'p/q' string or fmpq), not binary floats")
        f = Fraction(v)
        return fmpq(f.numerator, f.denominator)
    return [(q(r), q(t)) for r, t in pairs]


def __getattr__(name):
    if name == _ENGINE or name in _NAMES:      # lazy: needs python-flint
        mod = import_in_place(_ENGINE)
        globals()[_ENGINE] = mod
        for k in _NAMES:
            globals()[k] = getattr(mod, k)
        return globals()[name]
    raise AttributeError(f"module 'popcorn.enclosure' has no attribute {name!r}")

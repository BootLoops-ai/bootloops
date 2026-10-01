"""baller.certify — rigorous existence/uniqueness certification.

  pipe_vac       the generalized Krawczyk / interval-Newton engine with
                 forward-AD over acb and the Groenwall-Cauchy leg (VENDORED;
                 the model-specific sympy PF-ideal half is not part of this
                 package — see vendor header): certified zero-trapping +
                 unique-optimum proofs.
  block_krawczyk — generic
                 block-arrow structured Krawczyk — blockwise approximate
                 inverse + interval Schur complement on the border, strict-
                 interior contraction semantics — plus interval block-arrow
                 PD verification (pd_block_arrow / interval_cholesky) and
                 the combined certify_local_min front door. Module:
                 baller.block_arrow (battery leg L16).
"""
from ._core import load_vendored


def __getattr__(name):
    if name == "pipe_vac":
        mod = load_vendored("pipe_vac")
        globals()[name] = mod
        return mod
    if name == "block_krawczyk":
        from baller import block_arrow as mod
        globals()[name] = mod
        return mod
    raise AttributeError(f"module 'baller.certify' has no attribute {name!r}")


__all__ = ["pipe_vac", "block_krawczyk"]

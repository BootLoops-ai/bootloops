"""frobenius_boundary — reusable Frobenius-branch boundary machinery.

Convention: M_k(var) ~ var^(L*d/2 + ai_k + lambda_j), lambda_j = eig(D1),
D1 = R_Delta - diag(L*d/2 + ai).  See README.md (section 'Convention').
"""
import mpmath as mp
if mp.mp.dps < 50:
    mp.mp.dps = 50          # set FIRST (import-dps footgun)

from .core import (build_poly_DE, spectrum, classify, exclude_strata,
                   branch_series, phi_matrix)

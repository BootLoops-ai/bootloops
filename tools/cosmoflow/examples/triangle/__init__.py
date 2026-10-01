"""cosmoflow.examples.triangle -- the WORKED EXAMPLE of the cosmoflow package:
the FRW 1-loop 3-site (triangle) elliptic sub-sector of arXiv:2408.16386 and
its tree positive control, the two-site chain F(X1,X2,Y;eps) of 2312.05303.

  oracle.py -- K(m)-split nested-quadrature oracle for the 9 triangle masters
               e1..e9 (pinned reference e6(2,1/2,c=1,eps=0)), the two-site
               F/F1/F_twosite_certified evaluators and the GATE-ONLY print
               form F_twosite_print342; eval_master(n_s in {2,3}, ...).
  maxcut.py -- residue-period sampler at q_{G12}=0, the paper's L_2, its
               analytic periods varpi0/varpi1 and the symbolic L_2 proof.

General (any site graph) entry points live one level up: cosmoflow.alphabet
(alphabet_graph), cosmoflow.polytope (baikov_B, q_subgraph).
"""
from .oracle import (eval_master, eval_all_triangle, eval_master_lspace,
                     MASTERS, FINITE_EPS0,
                     F_twosite, F1_twosite, F_twosite_certified,
                     F_twosite_print342)
from .maxcut import (period_residue, L2_paper_coeffs, K2_modulus,
                     varpi0, varpi1, prove_L2_varpi0_symbolic)

__all__ = [
    'eval_master', 'eval_all_triangle', 'eval_master_lspace',
    'MASTERS', 'FINITE_EPS0',
    'F_twosite', 'F1_twosite', 'F_twosite_certified', 'F_twosite_print342',
    'period_residue', 'L2_paper_coeffs', 'K2_modulus',
    'varpi0', 'varpi1', 'prove_L2_varpi0_symbolic',
]

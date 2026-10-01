"""numkin kit — the "numkin" technique family: escape kits for
stuck/oversized kira/AMFlow reductions at the same pipeline stage.
  numkin_sweep.sh   — generation half: sed-substitute ONE numeric kinematic
                      value into already-generated SYSTEM_*.gz files, run
                      cheap (N-1)-var FireFly per point (shell entry point;
                      the flat path tools/numkin_sweep.sh is an exec-wrapper
                      here).
  numkin_harvest    — reconstruction half: exact bivariate (d, var) Pade
                      reassembly across the point grid, >=2 held-out verify.
  shift_opt         — per-sector loop-momentum-shift optimizer for AMFlow
                      eta-DE walls (unimodular affine shifts, ISP recompletion,
                      probe/harvest verbs).
  basisland         — basis-constrained exact function landing (verbs
                      rank|solve): #exact-samples ~ dimension, independent of
                      #variables, mixed point+functional rows, surplus-row
                      certificate.

Flat entry points at tools/ (numkin_sweep.sh, numkin_harvest.py,
shift_opt.py, basisland forms) forward here.
NOTE: numkin_harvest performs a load-bearing sys.path.insert of
tools/frobenius_boundary (kira_parse) on import — see its header comment.
This __init__ deliberately imports nothing (submodules pull flint/sympy;
keep `import numkin` cheap and side-effect-free).
"""

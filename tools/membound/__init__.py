"""membound — soft-region (omega->0) boundary constants as frequency-space
Bessel-kernel integrals for worldline / post-Minkowskian-type expansions.

eps-layer moments of user-supplied two- and three-frequency omega-cores
(`python3 -m membound --spec core.json`), assembled into boundary vectors in
the |omega|^{-k eps} layer basis, with retarded or advanced i0+ per
frequency factor ('ret'/'adv') and a phase-only diagnostic flag 'fey' that
is NOT the Feynman propagator prescription. Worked example and validation
gate: the post-Minkowskian memory family (gate_cM.py reproduces c_M = 1).
See GUIDE.md; battery: selftest.py.
"""

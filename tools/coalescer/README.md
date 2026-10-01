# coalescer — spectral-projector coalescence-coefficient extractor

Extracts the connection coefficient c_α of a holomorphic solution of a Fuchsian
operator onto a fractional-power (finite-monodromy) Frobenius branch Φ_α at another
singular point, generic over order / exponent set / eigenvalue set.  Two doors onto
one kernel: `coalescer.py --op FILE.json` takes YOUR operator L = Σ_j P_j(z) θ_z^j
(exact coefficients) and the exact Taylor coefficients of the solution at z = 0
(schema: GUIDE.md INPUTS; `dump_pf_data.py` writes the layout); `coalescer.py
--masses 1,1,1,9` is the worked example, the L-loop banana at threshold, with the
operator and the BFKNS seed built for you.  Projector (kills rank-n unipotent
log-tower + all other semisimple eigenspaces exactly):

    P_λ = (M₀−I)^n · ∏_{μ≠λ}(M₀−μI) / [(λ−1)^n · ∏_{μ≠λ}(λ−μ)]

**Validation (the banana worked example):** K3(1,1,1,9) c_{3/2} = −√3/(36π) at 195d (Arb ball-cert), 217d
same code path; CY₃(1,1,1,1,16) c_{5/4} = −1/(√(2π)Γ(¼)²), c_{7/4} =
−5Γ(¼)²/(384√2π^{5/2}) at ≥180d, product = 5/(768π³). Two-method (Route A
mpmath/s-coord ↔ Route B Arb/z-coord) gate 72–80d cross, 260d Route-B vs closed.

**Footgun:** anchor s_b must satisfy s_b > 1/min|z_sing| (MUM convergence
radius), NOT 1.5·max|z_sing| — Arb ball stays tight on the *wrong* sum if
seeded outside the disk.  The CLI checks this from the operator it builds
before any transport and refuses a violating anchor by name (rc 3); `--sb`,
`--sb2`, `--lift` override the defaults, `--check-only` prints the check, and
`--anchor-law seed|roots` picks the reading of min|z_sing| (seed: the leading
threshold, floor thr — exact for the multinomial seed; roots: every root of the
leading polynomial — conservative).  `coalescer.py --help` states the law with
the CY₄ (1,1,1,1,1,25) case as the example.


## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

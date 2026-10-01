# sec53_qorbit -- a 2x2 block with an irreducible quadratic letter, as a Counterweight control set

A 2x2 block of a two-loop differential equation (a gg -> Z gamma family, one sector, two masters) with an
irreducible quadratic letter q = 25x^2 - 43x + 25 (disc -651) in the chart 25s + 43 = 25(x + 1/x)
(s = (25x^2 - 43x + 25)/(25x); the block's own square root sqrt((25s-7)(25s+93)) = 25(x-1)(x+1)/x is rationalized).
The input blocks and the two pinned eps-forms are listed with their sha256 in MANIFEST.json.

What the fixture IS: a POSITIVE control for the rational-point Lee balances (the engine reaches a rational eps-form:
ok = true; the four letters x-1, x, x+1, q are read off out_engine_b53_chart1_x.json) and a NEGATIVE control for the
orbit-balance gate (the q-orbit residue is scalar, -eps I, odd part 0 -- pinned data in sec53_qorbit_block.json /
sec53_qorbit_epsform.json, not recomputed by the test; the gate must stay closed and no orbit move is needed). Two exact
eps-forms are pinned: a manual Lee-balance construction (out_lee_b53_chart1_x.json) and the engine's own output
(out_engine_b53_chart1_x.json); both verify by the exact point verifier vendored beside them (verify_epsform_point.py:
T A T^-1 + T' T^-1 == eps * Atilde entrywise in exact Fractions at two rational (eps, y) points, Atilde eps-free).
What it is NOT: a specimen of the orbit MOVE (no non-scalar orbit residue exists in this block in any rational chart of
this family).

Expected verdicts (each ASSERTED by test/test_sec53_qorbit.jl; the parenthesized reasons are explanatory, not asserted):
- b53_chart1_x.json      -> engine ok = true; verifier PASS on the fresh output and on out_lee_b53_chart1_x.json
- b53_raw_s.json         -> engine STOP by name 'no eps^0-reducing balance found and constant eps-decoupling failed ...' (half-integer exponents at 25s-7, 25s+93)
- b53_landau_x.json      -> engine STOP by name (the same; this chart adds 25x^2-143x+25 and keeps the half-integer exponents)
- planted_*.json         -> verifier FAIL rc 3 (a dropped balance; a mutated Atilde; eps inside Atilde)
- the row-major reading (--no-transpose) of the fresh output and of both pinned forms -> verifier FAIL rc 3 (a control of the reading convention)

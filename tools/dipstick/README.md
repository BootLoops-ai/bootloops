# DIPSTICK — pre-compute triage, unified

Check the integral BEFORE you commit the farm.

Verbs (CLI `./dipstick`, thin dispatcher — execs the members, no rewrites):

| verb      | member              | regime / job                                          |
|-----------|---------------------|-------------------------------------------------------|
| order     | pf_rank.jl          | Griffiths–Dwork mod-p PF order + holonomic rank, n<=4; 2-prime mandatory; saturation-refusal gate built in |
| count     | critcount.py        | n>=5 (closes order's NOT-FOR wall): Lagrange/ML-degree critical system -> msolve count = rank. RANK probe, NOT a solver |
| cycletype | cycletype.py        | Frobenius cycle-type statistics of the critical points; assumes squarefree reduction |
| nonres    | nonres.py           | GKZ facet + beta-pairing nonresonance certificate; certifies only at the tested beta |
| audit     | resaudit2.py        | msolve-INDEPENDENT exact resultant audit of a 2-roamer-variable .ms critical system (`--fast` = evaluate-interpolate variant for n>=20 lines); 2-variable `count` results should always be cross-checked through it |
| flatchi   | flatchi.py / flatchi_aff.py | exact flat-lattice Mobius chi for hyperplane arrangements: char poly, chi_top of the projective complement, Zaslavsky bounded count (with a .ms file: affine route) |
| regions   | regions/regions.jl  | expansion-by-regions completeness certificate: lower facets of Newt(U+F) under the small-parameter weight = the regions; volume-tiling + chi-additivity sums certify a user region list and NAME the missing scaling vector. julia + Oscar (engine-gated; skips by name without them) |

> **CHAR-P COEFFICIENT FOOTGUN:** msolve (0.6.5) silently garbles char-p
> input coefficients larger than 2^31 — wrong parse, rc 0, no warning; the
> solved system is then NOT the emitted system. EVERY .ms emitter MUST
> reduce every coefficient mod p to [0, p) at emission and refuse to emit
> otherwise. Integer-cleared exact-arithmetic systems are exactly the
> exposed class. After wiring the fix anywhere, re-run one known-answer
> wrong-then-right probe pair and one emitting-scale control before
> trusting any new verdict.

`tools/pf_rank.jl` forwards here (compat shim), providing
`pf_probe`/`pf_probe_gen`/`pf_probe_esc`.

Validation: `regress_order.jl` (7/7, run against member AND shim routes) +
`battery.py` (reproductions vs the reference receipts in `fixtures/` beside it;
its `regions` leg runs `regions/selftest.jl`, the three-input control battery
of the regions member).

Not here: series-only inputs (series_pf.py), LCLM factorization.
Manual: `GUIDE.md` in this directory.

## Credits

The master-integral counts of `order` and `regions` are computed by the Leviathan
engine (`upgrades/leviathan`), a Julia port of the Euler-characteristic-drop method of
Vsevolod Chestnov, Giulio Crisanti and Mathieu Giroux (arXiv:2606.29612), in part
translated from their Mathematica code and from DiscKosky by Giulio Crisanti, Luke
Lippstreu, Andrew J. McLeod and Maria Polackova. Dipstick calls that engine at run time
and contains none of its lines; see `upgrades/leviathan/NOTICE`. The region enumeration
of `regions` implements the Newton-polytope construction of A. Pak and A. Smirnov
(arXiv:1011.4863); the caveat on Glauber-type and potential regions is that of Jantzen,
Smirnov and Smirnov (arXiv:1206.0546). The credits of the other verbs are in the CREDIT
paragraph of `GUIDE.md`.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# Surd — format of the assembled closed-form files

`slice/scripts/assemble.py` writes one gzipped JSON document per channel or jet, `out/E4C_LO_{QCD,N4}_dipole_slice_<label>.json.gz` under the output
root. This page is the format specification, the evaluation and precision notes, and the census of the reference files for the LO collinear
four-point energy correlator (E4C) on the dipole slice (those files are large, 4.8-86 MB each, and are not included in the package; the census
lets a rebuilt file be checked against them).

## Object and normalization

Slice: detector positions w1 = (t/2)(3+4i)/5, w2 = -w1, w3 = -1/2, w4 = 1/2 (cos phi = 3/5); labeled distances d12 = t^2, d34 = 1,
d13 = d24 = t^2/4 + 3t/10 + 1/4, d14 = d23 = t^2/4 - 3t/10 + 1/4 (all polynomial in t).

Normalization: F_ch(t) = w_flavor(n_f = 5) * S_ch * sum_{sigma in S4} I_sigma(d(t)), with I_sigma the Cheng-Wu simplex integral of
`gate/scripts/oracle_cells.py` (flavor weights and identical-particle factors S_ch in `oracle_cells.flavor_weight` / `SFAC`); the quark and gluon
jet files are the sums over their channels. The xL^3-normalized value of `gate/scripts/oracle.py` (xL = max d_ab = max(1, t^2)) is
max(1,t^2)^3 F(t). The relabeling (13)(24) gives F(1/t) = t^6 F(t) exactly; files are constructed and checked on 0 < t < 1.

## The `format` block (identical in every file)

| key | meaning |
|---|---|
| `value` | F(t) = sum over `terms` of coef(t) * prod_{(L,s,e) in gens} rho_{L,s}(t)^e * prod const * prod_{w in Z} Z(w) |
| `coef` | LIST of summands, each `[A, B, D]` = coefficient lists (increasing powers of t) of (A(t) + i B(t))/D(t) in Q(i)(t); the coefficient is their sum |
| `gens` | `[L, s, e]`: rho_{L,s} = the s-th root in X of the letter polynomial L(X; t, j = i), raised to the integer power e (any fixed bijection slots <-> roots, used consistently everywhere) |
| points | `['Q', [A,B,D]]` = the exact point (A + iB)/D in Q(i)(t); `['A', L, s]` = rho_{L,s} |
| `Z` | Z(w) = ZIP_reg[w] = shuffle-regularized G(w_1, ..., w_n; oo) from 0 (ZeroInfPeriod semantics: log-divergences at 0 and oo dropped) = sum over letters a of {1 (coef -1)} + {a/(1+a) (coef +1, absent if a = -1)} of G(...; 1) |
| `const` | `['Z0', word]` = ZIP_reg of the word of hyperlog letters c * prod L^e (rational functions of X, t, j) evaluated at X = x0 = 11/2; `['G0', v]` = G_0(v; x0) = iterated integral from 0 to x0 with letters = points (G(0; x0) = log x0) |
| `branches` | every letter that lies exactly on the positive real integration path is `letter - i0` (the path passes above it); all other hyperlogs on their principal branch |
| `letters` | letter polynomial keys are JSON lists of `[[deg_X, deg_t, deg_j], [p, q]]` monomials; the top-level `letters` object gives their X-coefficients as Q(i)(t) elements (`GT` serialization of `slice/scripts/k3field.py`) |

Besides `format`, `letters` and `terms`, a file carries `object`, `channel`/`jet`, `complete` (all groups present), `missing_groups`, `x0`,
`census` and `groups` (one entry per contributing support group with its scale factor).

## Evaluation and precision

`slice/scripts/assemble.py: evaluate_file(path, Fraction t, dps)` returns an mpmath complex number; the exact value is REAL, so |Im F| of the
returned value is the error monitor. Digits delivered at working precision dps = D are about D minus the termwise cancellation depth of the
representation at that t: ~10 generically, ~25 near t = 4/5, up to ~30 within ~1e-2 of t = 0.4110 and t = 0.8110 (the roots of 15t^2 -+ 6t - 5:
poles of rational hyperlog points; F itself is analytic there) and near the Table-1 roots 0.8178, 0.9511 (example: q_qbqgq at t = 163/200 gives
34.2 digits at dps 60 and 37.4 at dps 80). Raise dps there, or use `slice/scripts/evaluate.py --file <label or path> --t p/q --digits N`, which
escalates dps until |Im F|/|F| meets the target and reports the digits actually delivered.

Representation-singular rational t: `evaluate_file` raises `ZeroDivisionError` where a letter's leading X-coefficient, a rational-point
denominator or a coefficient denominator vanishes exactly (t = 3/5 for 7 of the 9 QCD files: quark, gluon, q_qbqgq, q_gggq, g_qbqqbq, g_qbggq,
g_gggg; n4, q_qbpqpgq, g_qbpqpqbq evaluate directly; t = 5/8 for the files carrying the 8t -+ 5 class; F is analytic there). `evaluate.py`
returns the symmetric limit (F(t + delta) + F(t - delta))/2 there, escalating delta until the nearly colliding roots clear the 1e-12 path-distance
guard of `hpath.py`, and reports the O(F'' delta^2) error.

Slot invariance: the value is invariant under any permutation of the root slots of any letter (the term store lives in the symmetric
splitting algebra); only `symtest.py`'s two-point bookkeeping needs a definite slot bijection, which it fixes by canonical order plus continuation.

Per file, the representation-singular rational t and the irrational t near which precision is lost (raise dps):

| file | singular rational t | precision loss near t ~ |
|---|---|---|
| quark | 3/5, 5/8 | 0.06449216, 0.41101009, 0.81101009, 0.81776661, 0.95109994 |
| gluon | 3/5, 5/8 | 0.06449216, 0.41101009, 0.49394157, 0.60854676, 0.70710678, 0.81101009, 0.81776661, 0.95109994 |
| n4 | 3/5 | 0.06449216, 0.41101009, 0.81101009 |
| q_qbpqpgq | 3/5 | 0.06449216, 0.41101009, 0.81101009 |
| q_qbqgq | 3/5, 5/8 | 0.06449216, 0.41101009, 0.81101009, 0.81776661, 0.95109994 |
| q_gggq | 3/5 | 0.06449216, 0.41101009, 0.81101009 |
| g_qbpqpqbq | 3/5 | 0.06449216, 0.41101009, 0.81101009 |
| g_qbqqbq | 3/5 | 0.06449216, 0.41101009, 0.81101009 |
| g_qbggq | 3/5, 5/8 | 0.06449216, 0.41101009, 0.81101009, 0.81776661, 0.95109994 |
| g_gggg | 3/5, 5/8 | 0.06449216, 0.41101009, 0.49394157, 0.60854676, 0.70710678, 0.81101009, 0.81776661, 0.95109994 |

## Census of the reference files

| file | sha256[:16] | MB | groups | terms | distinct Z | constants | points | algebraic X-letters | root generators | max deg_t (num, den) | weight histogram |
|---|---|---|---|---|---|---|---|---|---|---|---|
| n4 | `ec709955cf60e52a` | 4.8 | 72 | 21348 | 2959 | 682 | 107 | 6 | 10 | [168, 166] | {0: 1, 1: 49, 2: 1961, 3: 19337} |
| g_gggg | `5942411214f5c410` | 85.7 | 36 | 40042 | 5167 | 706 | 171 | 36 | 40 | [424, 416] | {0: 1, 1: 49, 2: 2072, 3: 37920} |
| g_qbggq | `fa3a2ab01e3dda93` | 59.6 | 24 | 28284 | 3987 | 566 | 135 | 22 | 26 | [426, 417] | {0: 1, 1: 47, 2: 1716, 3: 26520} |
| g_qbpqpqbq | `db0c99b3a8310f0a` | 22.5 | 36 | 13086 | 1725 | 423 | 87 | 6 | 10 | [356, 348] | {0: 1, 1: 47, 2: 1330, 3: 11708} |
| g_qbqqbq | `0516d448467d420f` | 20.6 | 18 | 11422 | 1572 | 370 | 76 | 5 | 8 | [366, 358] | {0: 1, 1: 39, 2: 1140, 3: 10242} |
| gluon | `09197bf72a040454` | 81.5 | 114 | 41034 | 5167 | 706 | 171 | 36 | 40 | [438, 429] | {0: 1, 1: 49, 2: 2072, 3: 38912} |
| q_gggq | `d9b13df8008969ff` | 60.2 | 72 | 23160 | 3427 | 681 | 111 | 6 | 10 | [352, 346] | {0: 1, 1: 49, 2: 2013, 3: 21097} |
| q_qbpqpgq | `6e7799319b3bd80b` | 40.0 | 48 | 18181 | 2614 | 620 | 107 | 6 | 10 | [264, 260] | {0: 1, 1: 49, 2: 1807, 3: 16324} |
| q_qbqgq | `b5302d792e01aff5` | 51.9 | 36 | 27873 | 3903 | 612 | 135 | 18 | 22 | [354, 348] | {0: 1, 1: 48, 2: 1856, 3: 25968} |
| quark | `574e885b27223e08` | 58.2 | 156 | 31117 | 4243 | 698 | 135 | 18 | 22 | [362, 356] | {0: 1, 1: 49, 2: 2047, 3: 29020} |

terms = distinct hyperlog monomials after exact collection; every file is complete (all support groups present).

## The non-linear X-letters (their roots are hyperlog points)

36 distinct polynomials L(X; t) of degree 2 or 3 in X occur across the files (j = i has been substituted; these are real on the slice):

| deg_X | L(X; t) | files |
|---|---|---|
| 3 | 25 - 30*t + 25*t^2 - 75*X*t + 118*X*t^2 - 15*X*t^3 - 15*X^2*t + 118*X^2*t^2 - 75*X^2*t^3 + 25*X^3*t^2 - 30*X^3*t^3 + 25*X^3*t^4 | n4, g_gggg, g_qbggq, g_qbpqpqbq, g_qbqqbq, gluon, q_gggq, q_qbpqpgq, q_qbqgq, quark |
| 3 | 25 + 30*t + 25*t^2 + 75*X*t + 118*X*t^2 + 15*X*t^3 + 15*X^2*t + 118*X^2*t^2 + 75*X^2*t^3 + 25*X^3*t^2 + 30*X^3*t^3 + 25*X^3*t^4 | n4, g_gggg, g_qbggq, g_qbpqpqbq, g_qbqqbq, gluon, q_gggq, q_qbpqpgq, q_qbqgq, quark |
| 2 | 25 + 14*t^2 + 25*t^4 - 150*X - 28*X*t^2 + 50*X*t^4 + 25*X^2 + 14*X^2*t^2 + 25*X^2*t^4 | n4, g_gggg, g_qbggq, g_qbpqpqbq, g_qbqqbq, gluon, q_gggq, q_qbpqpgq, q_qbqgq, quark |
| 2 | 25 + 14*t^2 + 25*t^4 + 50*X - 28*X*t^2 - 150*X*t^4 + 25*X^2 + 14*X^2*t^2 + 25*X^2*t^4 | n4, g_gggg, g_qbggq, g_qbpqpqbq, g_qbqqbq, gluon, q_gggq, q_qbpqpgq, q_qbqgq, quark |
| 3 | 25*t^2 - 30*t^3 + 25*t^4 - 15*X*t + 118*X*t^2 - 75*X*t^3 - 75*X^2*t + 118*X^2*t^2 - 15*X^2*t^3 + 25*X^3 - 30*X^3*t + 25*X^3*t^2 | n4, g_gggg, g_qbggq, g_qbpqpqbq, gluon, q_gggq, q_qbpqpgq, q_qbqgq, quark |
| 3 | 25*t^2 + 30*t^3 + 25*t^4 + 15*X*t + 118*X*t^2 + 75*X*t^3 + 75*X^2*t + 118*X^2*t^2 + 15*X^2*t^3 + 25*X^3 + 30*X^3*t + 25*X^3*t^2 | n4, g_gggg, g_qbggq, g_qbpqpqbq, g_qbqqbq, gluon, q_gggq, q_qbpqpgq, q_qbqgq, quark |
| 2 | -1 + 1*t^2 - 2*X - 1*X^2 + 1*X^2*t^2 | g_gggg, gluon |
| 2 | -1 + 1*t^2 + 2*X*t^2 - 1*X^2 + 1*X^2*t^2 | g_gggg, gluon |
| 2 | -10 + 10*t^2 - 5*X - 6*X*t - 5*X*t^2 + 5*X^2 + 6*X^2*t + 5*X^2*t^2 | g_gggg, g_qbggq, gluon |
| 2 | -10 + 10*t^2 - 5*X + 6*X*t - 5*X*t^2 + 5*X^2 - 6*X^2*t + 5*X^2*t^2 | g_gggg, g_qbggq, gluon |
| 2 | -5 - 6*t - 5*t^2 + 5*X + 6*X*t + 5*X*t^2 - 10*X^2 + 10*X^2*t^2 | g_gggg, g_qbggq, gluon, q_qbqgq, quark |
| 2 | -5 + 6*t - 5*t^2 - 15*X + 18*X*t - 15*X*t^2 - 20*X^2 + 12*X^2*t | g_gggg, gluon, q_qbqgq, quark |
| 2 | -5 + 6*t - 5*t^2 + 5*X - 6*X*t + 5*X*t^2 - 10*X^2 + 10*X^2*t^2 | g_gggg, g_qbggq, gluon, q_qbqgq, quark |
| 2 | 10 - 10*t^2 - 5*X - 6*X*t - 5*X*t^2 + 5*X^2 + 6*X^2*t + 5*X^2*t^2 | g_gggg, gluon |
| 2 | 10 - 10*t^2 - 5*X + 6*X*t - 5*X*t^2 + 5*X^2 - 6*X^2*t + 5*X^2*t^2 | g_gggg, gluon |
| 2 | 20 - 12*t + 15*X - 18*X*t + 15*X*t^2 + 5*X^2 - 6*X^2*t + 5*X^2*t^2 | g_gggg, g_qbggq, gluon |
| 2 | 20 - 12*t + 15*X - 30*X*t + 15*X*t^2 - 12*X^2*t + 20*X^2*t^2 | g_gggg, gluon |
| 2 | 20 + 12*t + 15*X + 18*X*t + 15*X*t^2 + 5*X^2 + 6*X^2*t + 5*X^2*t^2 | g_gggg, g_qbggq, gluon |
| 2 | 20 + 12*t + 15*X + 30*X*t + 15*X*t^2 + 12*X^2*t + 20*X^2*t^2 | g_gggg, gluon |
| 2 | 5 - 3*t - 6*X*t + 5*X*t^2 + 5*X^2*t^2 | g_gggg, g_qbggq, gluon |
| 2 | 5 - 6*t + 5*t^2 - 5*X + 6*X*t - 5*X*t^2 - 10*X^2 + 10*X^2*t^2 | g_gggg, gluon, q_qbqgq, quark |
| 2 | 5 - 6*t + 5*t^2 + 15*X - 18*X*t + 15*X*t^2 - 12*X^2*t + 20*X^2*t^2 | g_gggg, g_qbggq, gluon, q_qbqgq, quark |
| 2 | 5 + 3*t + 6*X*t + 5*X*t^2 + 5*X^2*t^2 | g_gggg, g_qbggq, gluon |
| 2 | 5 + 6*t + 5*t^2 - 5*X - 6*X*t - 5*X*t^2 - 10*X^2 + 10*X^2*t^2 | g_gggg, gluon, q_qbqgq, quark |
| 2 | 5 + 6*t + 5*t^2 + 15*X + 18*X*t + 15*X*t^2 + 20*X^2 + 12*X^2*t | g_gggg, gluon, q_qbqgq, quark |
| 2 | 5 + 6*t + 5*t^2 + 15*X + 18*X*t + 15*X*t^2 + 12*X^2*t + 20*X^2*t^2 | g_gggg, g_qbggq, gluon, q_qbqgq, quark |
| 2 | 5 + 5*X - 6*X*t - 3*X^2*t + 5*X^2*t^2 | g_gggg, g_qbggq, gluon, q_qbqgq, quark |
| 2 | 5 + 5*X + 6*X*t + 3*X^2*t + 5*X^2*t^2 | g_gggg, g_qbggq, gluon, q_qbqgq, quark |
| 2 | -12*t + 20*t^2 + 15*X - 18*X*t + 15*X*t^2 + 5*X^2 - 6*X^2*t + 5*X^2*t^2 | g_gggg, gluon |
| 2 | -3*t + 5*t^2 + 5*X - 6*X*t + 5*X^2 | g_gggg, g_qbggq, gluon |
| 2 | 12*t - 20*t^2 - 15*X + 30*X*t - 15*X*t^2 - 20*X^2 + 12*X^2*t | g_gggg, gluon |
| 2 | 12*t + 20*t^2 + 15*X + 18*X*t + 15*X*t^2 + 5*X^2 + 6*X^2*t + 5*X^2*t^2 | g_gggg, gluon |
| 2 | 12*t + 20*t^2 + 15*X + 30*X*t + 15*X*t^2 + 20*X^2 + 12*X^2*t | g_gggg, gluon |
| 2 | 3*t + 5*t^2 + 5*X + 6*X*t + 5*X^2 | g_gggg, g_qbggq, gluon |
| 2 | -5*t^2 + 6*X*t - 5*X*t^2 - 5*X^2 + 3*X^2*t | g_gggg, g_qbggq, gluon, q_qbqgq, quark |
| 2 | 5*t^2 + 6*X*t + 5*X*t^2 + 5*X^2 + 3*X^2*t | g_gggg, g_qbggq, gluon, q_qbqgq, quark |


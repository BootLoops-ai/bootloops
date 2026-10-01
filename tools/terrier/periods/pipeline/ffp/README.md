# periods/pipeline/ffp — finite-field fingerprint engines

Engine code behind the facade `periods/ffp.py` (see its docstring for the
front door, the int64 prime law and the list of importable gates; battery
`periods/selftest_ffp.py`).

| file | content |
|---|---|
| `field.py`, `ext_field.py` | F_q index arithmetic (q = p, p^2) and batched F_{p^m} arithmetic |
| `hv_count.py` | AESZ34 / Hulek-Verrill torus fibre counts over F_q for all psi in one pass |
| `brute.py` | brute-force validators of the counts at small q |
| `weil.py`, `quartic.py`, `modsqrt.py` | exact Weil-circle test, rank-2 split of the Frobenius quartic, Tonelli-Shanks |
| `groundtruth_p19.py` | the p = 19 Frobenius table of arXiv:1912.06146 for AESZ34, replayed as a gate |
| `anchors.py` | independent arithmetic anchors: a_p of the elliptic curve 14a1 by direct count, loaders for the newform tables below |
| `pilot2.py` | AESZ34 pass-set driver at one prime (E-law, k-classes, F_{p^2} side, Weil windows, quartic split) |
| `scale/bcm*.py` | BCM finite-hypergeometric a_p(z) engine (Gauss sums), its calibration, validation battery, production pass-sets and control gate |
| `scale/check_candidate.py`, `scale/checker_gate.py` | candidate checker (per-prime IN/OUT + AND verdict) and its regression gate |
| `scale/passsets/*.json` | stored pass-sets: AESZ34 at 14 primes 53..127, AESZ4 and AESZ11 at 12 primes 61..337 |
| `banks/*.jsonl` | reference gate outputs (`BCM_GATE.jsonl`, `CHECKER_GATE.jsonl`) that the battery regenerates and compares |

## Data sources

The files `nf_14_4_a_a.json`, `nf_34_2_b_a.json` and `scale/nf_180_4_a_e.json`
hold Hecke eigenvalues a_p of the classical newforms with LMFDB labels
14.4.a.a, 34.2.b.a and 180.4.a.e, taken from the L-functions and Modular
Forms Database (The LMFDB Collaboration, https://www.lmfdb.org). The Hecke
eigenvalue tables in this directory are taken from the LMFDB (The L-functions
and Modular Forms Database, https://www.lmfdb.org), whose data is licensed
CC BY-SA 4.0 (https://www.lmfdb.org/license).

* `nf_14_4_a_a.json` — {"p": a_p} for all 168 primes p < 1000 (p = 2 ... 997);
  weight 4, level 14. Used with the point count of the elliptic curve 14a1
  (`anchors.ap_14a`) as the rank-2 anchor T1 of the AESZ34 pass-sets.
* `nf_34_2_b_a.json` — {"p": a_p} at the nine odd primes p < 100 that split
  in Q(sqrt 17): 13, 19, 43, 47, 53, 59, 67, 83, 89; weight 2, level 34.
* `scale/nf_180_4_a_e.json` — {"label", "source", "note", "ap": {"p": a_p}}
  for the 23 primes 5 <= p <= 97 (a_2 = a_3 = 0 since 4, 9 | 180); weight 4,
  level 180. Used by `scale/bcm_gate.py` to anchor the beta eigenvalue of the
  AESZ11 control point z = -1/432.

Literature anchors used by the gates: the p = 19 Frobenius-polynomial table
for AESZ34 of Candelas, de la Ossa, Elmi and van Straten, arXiv:1912.06146
(`groundtruth_p19.py`); the rank-2 attractor points of AESZ4 (z = -1/5832)
and AESZ11 (z = -1/432) with their modular-form partners as given in
arXiv:2203.09426 (`scale/bcm_passsets.py`, `scale/bcm_gate.py`).

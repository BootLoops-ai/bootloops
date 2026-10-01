# Longhand — independent numerical ground truth for Feynman integrals, done the slow honest way

When a fast method hands you a number for a Feynman integral — a differential-equation
transport, an auxiliary-mass-flow run, a fitted ansatz — the question that decides whether
you may build on it is whether a computation sharing none of that machinery agrees, and to
how many digits, with an error statement that was measured rather than assumed. Longhand
is the package that answers it, by the two slow routes that need the least trust.

**Route A — the arbitrary-precision parametric evaluator** (`hiprec.py`, `feynman.py`,
`parametric.py`; manual [`PARAMETRIC.md`](PARAMETRIC.md)). Evaluates the Feynman-parametric
representation of a finite (eps^0) integral directly in arbitrary precision in the
deep-Euclidean region, after integrating out every parameter that enters F linearly in
closed form. Spectral to 40+ digits when the effective dimension is ≤3 (the one-loop
massive box: two independent reductions agree to 46-47 digits at dps=45/60; pinned 40-digit
reference `0.1245570572327092697104024441545542480922` at s=−1, t=−1/3, m²=1, M5=2), a few
digits with a statistical bar by parallel CBC-lattice QMC at irreducible dimension 4-12, and
any finite projective integrand written in the UF-spec JSON schema. No pySecDec needed except
to turn a propagator list into U and F.

```python
from fractions import Fraction
from longhand import box, integrate_spec            # tools/ on sys.path
box(-1, Fraction(-1, 3), 1, 2, dps=35)               # exact kinematics in, 35 digits out
integrate_spec('toy_qmc_spec.json', N=8009, n_shifts=8, dps=25, nproc=8)   # {'value', 'error', 'digits', ...}
```

**Route B — the pySecDec disteval driver** (`longhand.disteval`, command line
`disteval/pysecdec_point.py`; manual [`DISTEVAL.md`](DISTEVAL.md)). Runs a compiled pySecDec
`loop_package` at one point as three disteval stages — two lattice sizes and the
median-lattice second integrator — split into chunks of sectors that are resumable
checkpoints, then compares the stages per regulator order and states the digits by object:
the minimum over the B-C pair and their bars, with a planted-fail control that proves the
compare can see a wrong digit and a foreign control against coarse independent runs of the
same library. All eps orders, double precision (pySecDec's ceiling), pySecDec's numerics
(GPL-3.0, `pip install pySecDec`, never bundled), our staging and checks.

```
cd tools/longhand/disteval
python3 pysecdec_point.py --disteval-dir PKG/disteval --stage A --workers 8 --work WORK     # then --stage B, --stage C
python3 pysecdec_point.py --stage compare --results A=WORK/stage_A/RESULT_A.json B=... C=... --out WORK/COMPARE.json
python3 pysecdec_point.py --stage check --work WORK --compare WORK/COMPARE.json --out WORK/GATE.json
python3 pysecdec_point.py --stage check --example        # replay the worked example's records: 21 checks by name
```

**Worked example** ([`examples/`](examples/README.md)): the generic-mass five-point
double-pentagon top-sector scalar (332 sectors, 968 kernels) at a Euclidean point through
route B — eps^0 = 8.246341248953e-01 ± 7.3e-08, digit class 7, leading pole to 9 digits —
whose stage results, compare, acceptance check and memory readbacks ship as sha256-pinned
records the battery replays.

**Battery:** `python3 selftest.py` (both routes; about a minute with pySecDec installed,
~13 s without — the pySecDec-dependent legs skip by name); `--full` for route A's
multi-process tier. **Requirements:** `mpmath`, `sympy`, `numpy`, `gmpy2`, `pytest`;
pySecDec optional (GPL-3.0, user-installed). Everything else — when to use which route,
the UF-spec schema, every flag, the acceptance criteria and the known hazards — is in
[`GUIDE.md`](GUIDE.md).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# nestor.dispersion — exponentially convergent singularity-subtracted dispersion quadrature

Member of the `nestor` package (`tools/nestor/dispersion/`); package guide: [`../GUIDE.md`](../GUIDE.md).

Reusable tool for any **2-particle-reducible self-energy-insertion** dispersion integral

```
        I = (1/pi) ∫ rho(w') K(w') dw'
```

where `rho = Im Sigma(w')` is an **expensive** spectral density carrying threshold
non-analyticities (a turn-on at the lowest threshold, a cusp/jump at each higher normal threshold),
and `K(w')` is a **smooth analytic** kernel (e.g. the one-loop box `Box1(s,t; M²=w')`).

## The problem it solves

A plain **Gauss–Legendre** dispersion quadrature converges only **polynomially** when `rho` has a
threshold singularity — e.g. `(w'−1)·log(w'−1)` (a 2-particle turn-on with a massless line) or a
finite **jump** at a higher cut. On the LBL3SE benchmark this floored the assembly at ~6 digits with
180 nodes.

Two fixes, in order of impact:

1. **Use tanh-sinh (double-exponential) with panel breaks exactly at the thresholds.** tanh-sinh
   natively resolves *integrable* endpoint singularities (`(w'−w*)^α`, `log`, …) with exponential
   convergence; a panel boundary placed *at* a finite **jump** makes each side analytic. This alone
   turns 6 d → tens of digits. (Surrogate: GL 3,5,6,7 d at n=20..160 → tanh-sinh 15,29,99 d at
   level 3,4,5.)

2. **Singularity subtraction** (this tool) for the cases tanh-sinh alone does not nail and to cut the
   number of expensive `rho` evaluations: on a short sub-panel adjacent to the threshold `w*`,
   subtract the known leading singular model
   ```
        S(w') = Σ_j c_j (w'−w*)^{α_j} [log(w'−w*)]^{m_j}
   ```
   so `rho − S` is analytic (tanh-sinh exponential), and **add the subtracted piece back in closed
   form**: expand the analytic kernel `K(w') = Σ_n K_n (w'−w*)^n` and use the exact power-log moments
   ```
        ∫₀^L u^{α} (log u)^m du = M_{α,m}(L)   (closed form, see moment_powerlog)
   ```
   so `∫ S·K = Σ c_j K_n M_{α_j+n, m_j}(L)`. The `K_n` series truncates exponentially (analytic
   kernel); the sub-panel length `L` must be **smaller than the kernel's analyticity radius**.

## API (`disp_sub.py`; `from nestor.dispersion import ...` with `tools/` on the path)

- `moment_powerlog(p, m, L)` — closed form of `∫₀^L u^p (log u)^m du`, `m ∈ {0,1,2}`, `p > −1`.
- `kernel_taylor(K, wstar, nmax, radius, dps, method='cheb')` — local Taylor coefficients of the
  analytic kernel. Default `'cheb'` (Chebyshev collocation) uses **real** evaluations only (the box
  is real-only); `'quad'` is Cauchy-integral (needs a complex-capable kernel).
- `addback_endpoint(coeffs, Kc, L, side)` — closed-form `∫ S·K` over a one-threshold sub-panel.
- `disp_subtracted(rho, K, panels, dps, kernel_nmax, maxdegree)` — the full assembly. Each panel:
  `{'a','b','thresh': None | ('left'|'right', w*, [(c,α,m),…]), 'sub_width', 'radius', 'map':'tail'}`.
  A pure jump at an interior threshold needs **no model** — just give it a panel boundary
  (`thresh=None` each side).

## Tests

All four legs run inside `python3 -m nestor.selftest` (D1–D4); the files are also pytest-collectable
(`python3 -m pytest tools/nestor/tests -q -p no:cacheprovider` from the repository root).

- `../tests/test_dispersion_moments.py` — closed-form moments vs adaptive quad (≥30 d) + add-back on
  a model kernel (33 d).
- `../tests/test_dispersion_surrogate.py` — a fast analytic surrogate with the **exact** LBL3SE
  singularity structure (`(w'−1)log(w'−1)` turn-on + `√(w'−9)` onset). Asserted legs: the level-5
  subtracted assembly must match a high-precision tanh-sinh reference to ≥30 digits, and an
  add-back-sign mutation control must collapse the agreement (the closed-form add-back is
  load-bearing). Run directly (`python3 tools/nestor/tests/test_dispersion_surrogate.py`) it also
  prints the GL-polynomial vs tanh-sinh/subtracted-exponential convergence ladder.

## Worked kernel

`../examples/box1_dilog.py` — the one-loop box (three equal masses + one large corner mass) in
1D dilogarithmic form with Vieta-stabilized small roots; `make_K(s, t, m2, dps)` returns a memoized
kernel `K(w')` ready for `disp_subtracted` panels. See the package guide for its region and gates.

## LBL3SE application

Measured singular structure: w'=1 turn-on
`Im Σ ~ (w'−1)(−2π·log(w'−1) + B)` (subtract the `−2π(w'−1)log(w'−1)`); w'=9 elliptic 3-massive cut
is a **finite jump** (equal-mass sunrise disc is analytic-with-jump at threshold) → panel break, no
model.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

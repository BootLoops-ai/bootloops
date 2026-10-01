# nullspace_fast

**What:** int64 numpy GF(p) RREF (`nullspace_mod_fast`) + polynomial-coefficient
recurrence finder (`find_recurrence_fast`) for integer sequences mod p (p < 2^31).
Vectorized row elimination gives 10-100x over `sympy.Matrix.nullspace` for mod-p
nullspaces; the recurrence finder sweeps (order r, coeff-degree d) candidates in
increasing cost order, solves the (r+1)(d+1)-column ansatz with `overdet` extra
rows, and verifies the found relation on every available row before returning.

**Reference scale:** a phylogenetics GKZ application — found the order-24, degree-19 contiguity recurrence of the full
4-taxon Jukes-Cantor marginal likelihood from a 550-term mod-p Z-sequence in
204 s (2-prime stable, p=2147483629 / 2147483587), where the sympy route was
ETA-prohibitive. Same role as `sparse_cascade.py`'s GF(p) incremental nullspace,
but for dense recurrence ansatz matrices rather than sparse symbol cascades.

**API:**
- `nullspace_mod_fast(A, p) -> list | None` — one nullvector of A mod p (or None if full column rank).
- `find_recurrence_fast(Zs, p, rmax=40, dmax=40, overdet=5, verbose=False) -> (r, d, v) | None` —
  smallest-cost recurrence `sum_{j=0}^{r} c_j(n) Z_{n+j} = 0 mod p` with
  `c_j(n) = sum_i v[j*(d+1)+i] n^i`; `v` is the flattened coefficient vector.

**Caveats:** p must satisfy p^2 * max-matrix-dim < 2^63 (int64 products before
the `% p`; the shipped PRIMES ~2^31 are fine for nrows up to ~few thousand).
Single nullvector only (last free column), not the full nullspace basis. Always
confirm a found (r,d) at a second prime before claiming an order (2-prime rule).

Run `python3 test_nullspace_fast.py` (a few seconds; includes the
`annihilator.py --selftest` exit-status gate and its mutation control).

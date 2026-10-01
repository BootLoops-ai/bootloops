"""ellipticus.gmtel — gmtel member: certified Gauss-Manin q-telescoper.

Exact Picard-Fuchs operators of K3 PENCIL periods derived directly from a
Weierstrass model (exact 2-parameter fiber Gauss-Manin + creative telescoping
over the base q with an explicit certificate row-vector), with exact
certificates.  The Torelli-grade PF fingerprint instrument; the GEOMETRIC half
of a two-route pattern whose data half is a series + annihilator fit
(tools/annihilator, pf_from_series).

Modules:
  famgen_pf_routeA — the telescoper (route A, geometric): exact fiber GM
                     matrices M_q/M_s -> mod-p telescoper shape search (flint
                     nmod_poly, 3 primes) -> canonical rational-function
                     reconstruction + CRT/Farey lift -> EXACT certificate
                     re-solve over Q[q] at rational s points -> x-chart gates
                     vs route-B operator and the reference L5 operator (gauge
                     twist exhibited).
  famgen_pf_routeB — route B (data): exact moments -> annihilator
                     pf_from_series (two-prime (r,s) pin, exact Q-nullspace,
                     held-out terms) + reference-L5 comparator.  Rides along
                     so the member's battery reproduces the full two-route
                     control; the general-purpose data route is
                     tools/annihilator itself.
  famgen_common    — shared exact operator algebra (canonical D-forms,
                     primitive normalization, theta-form/chart translation,
                     1/x chart inversion, moments, reference-L5 loader).
  gmtel_selftest   — the member's own battery: two routes on the Watson/HLY
                     (2,2,2) family at (A,B,C)=(1,4,16) must land EXACTLY on
                     the reference L5 operator up to the exhibited x^(-1/2)
                     twist.

Run modules as scripts (they sys.path themselves):
  python3 tools/ellipticus/gmtel/gmtel_selftest.py            # full battery
  python3 tools/ellipticus/gmtel/gmtel_selftest.py --pilot    # fast smoke
Env: FAMGEN_L5 (reference L5 operator; default the packaged copy in
fixtures/gmtel/L5_theta.txt), FAMGEN_TOOLS (dir holding the annihilator
engine module; default the sibling tools/annihilator/ in this tree).
"""

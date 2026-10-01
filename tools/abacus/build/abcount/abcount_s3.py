"""abcount stage S3 — the general (2,2) layer.

O_K = Z[sqrt(-5)] action verification (exact rho receipts: rho^2 = -5,
lattice preservation, Rosati compatibility, analytic C-linearity + the
holomorphic-eigenvalue sign pin s), polarization-type extraction (Smith /
Frobenius normal form of the alternating form — the N2 predicate), the
B4/B4s counting routes (Path 1 = PARI ellap / torsion over the residue
fields; Path 2 = hecke_desk-style exact Z[sqrt(-5)] Hecke evaluation, NO
PARI), the B4s rho/-rho detector, and the B5 Weil-restriction routes
(PARI hyperellcharpoly over F_p / F_p^2 vs pure-python tower counts).

Seeds: battery/SEEDS.json S3 revision (sha quoted by run_s3.py).
Sign law: sign_convention s1 = "P" is the pin of record (Q6); the SINGLE
conversion site remains abcount_s0.convert_sign_convention — nothing in
this module converts; the sign-pin receipt here only VERIFIES and REFUSES
(FAIL-SIGN-PIN), it never relabels.
"""

import json
import math
import subprocess

from flint import acb, arb, ctx

# ----------------------------------------------------------------------------
# exact integer matrices (small, dense, python ints)
# ----------------------------------------------------------------------------

def mat_mul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)]
            for i in range(n)]


def mat_scal(c, A):
    return [[c * x for x in row] for row in A]


def mat_add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def mat_T(A):
    return [list(r) for r in zip(*A)]


def mat_eye(n, c=1):
    return [[c if i == j else 0 for j in range(n)] for i in range(n)]


def block_diag(blocks):
    n = sum(len(b) for b in blocks)
    M = [[0] * n for _ in range(n)]
    off = 0
    for b in blocks:
        for i, row in enumerate(b):
            for j, x in enumerate(row):
                M[off + i][off + j] = x
        off += len(b)
    return M


RHO_BLOCK = [[0, -5], [1, 0]]     # SEEDS B4.rho_8x8 per-factor block, basis {1, rho}
E2_BLOCK = [[0, 1], [-1, 0]]      # principal symplectic form per factor


def rho_product_matrix(tamper=None):
    """8x8 diagonal O_K-action; tamper = (i, j, delta) 0-based adds delta (N3)."""
    M = block_diag([RHO_BLOCK] * 4)
    if tamper is not None:
        i, j, d = tamper
        M[i][j] += d
    return M


def product_polarization_E(divisors=(1, 1, 1, 1)):
    """Alternating form; per factor [[0, d], [-d, 0]] (principal iff all d = 1)."""
    return block_diag([[[0, d], [-d, 0]] for d in divisors])


# ----------------------------------------------------------------------------
# O_K-action receipts (exact integer predicates; named, machine-readable)
# ----------------------------------------------------------------------------

def ok_action_receipts(rho, E):
    """Exact integer predicates on (rho, E).  Returns (ok, receipt) with the
    FIRST failing predicate named in receipt['failing_predicate']."""
    n = len(rho)
    preds = []
    r2 = mat_mul(rho, rho)
    preds.append(("OK-action-rho-square-eq-minus5",
                  r2 == mat_eye(n, -5),
                  {"claim": "rho^2 = -5*Id (exact integer)",
                   "residual_nonzero_entries":
                       [[i, j, r2[i][j] - (-5 if i == j else 0)]
                        for i in range(n) for j in range(n)
                        if r2[i][j] != (-5 if i == j else 0)][:8]}))
    preds.append(("OK-action-lattice-preservation",
                  all(isinstance(x, int) for row in rho for x in row),
                  {"claim": "rho is an integer matrix on the pinned lattice basis"}))
    m1 = mat_mul(mat_T(rho), mat_mul(E, rho))
    preds.append(("OK-action-rosati-scaling",
                  m1 == mat_scal(5, E),
                  {"claim": "E(rho x, rho y) = 5 E(x, y): rho^T E rho = 5E"}))
    m2 = mat_add(mat_mul(mat_T(rho), E), mat_mul(E, rho))
    preds.append(("OK-action-rosati-antiadjoint",
                  all(x == 0 for row in m2 for x in row),
                  {"claim": "E(rho x, y) + E(x, rho y) = 0 (Rosati(rho) = -rho)"}))
    failing = next((name for name, ok, _ in preds if not ok), None)
    return failing is None, {
        "predicates": [{"name": nm, "ok": bool(ok), **det} for nm, ok, det in preds],
        "failing_predicate": failing}


# ----------------------------------------------------------------------------
# polarization type: Smith normal form of the alternating integer form (N2)
# ----------------------------------------------------------------------------

def smith_divisors(M):
    """Elementary divisors (SNF diagonal, each dividing the next) of an integer
    matrix — full exact integer algorithm, small matrices."""
    A = [list(r) for r in M]
    n, m = len(A), len(A[0])
    divs = []
    r = 0
    while r < min(n, m):
        piv = None
        for i in range(r, n):
            for j in range(r, m):
                if A[i][j] != 0 and (piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            break
        i0, j0 = piv
        A[r], A[i0] = A[i0], A[r]
        for row in A:
            row[r], row[j0] = row[j0], row[r]
        while True:
            done = True
            for i in range(r + 1, n):
                if A[i][r] % A[r][r] != 0:
                    q = A[i][r] // A[r][r]
                    for j in range(m):
                        A[i][j] -= q * A[r][j]
                    A[r], A[i] = A[i], A[r]
                    done = False
            for j in range(r + 1, m):
                if A[r][j] % A[r][r] != 0:
                    q = A[r][j] // A[r][r]
                    for i in range(n):
                        A[i][j] -= q * A[i][r]
                    for i in range(n):
                        A[i][r], A[i][j] = A[i][j], A[i][r]
                    done = False
            if done:
                break
        for i in range(r + 1, n):
            q = A[i][r] // A[r][r]
            for j in range(m):
                A[i][j] -= q * A[r][j]
        for j in range(r + 1, m):
            q = A[r][j] // A[r][r]
            for i in range(n):
                A[i][j] -= q * A[i][r]
        divs.append(abs(A[r][r]))
        r += 1
    # enforce divisibility chain
    for i in range(len(divs)):
        for j in range(i + 1, len(divs)):
            g = math.gcd(divs[i], divs[j])
            l = divs[i] * divs[j] // g
            divs[i], divs[j] = g, l
    return sorted(divs)


def polarization_type_check(E, declared):
    """Named predicate: alternating integrality + SNF pair structure + declared
    elementary divisors.  Returns (ok, receipt)."""
    n = len(E)
    alt_ok = all(E[i][j] == -E[j][i] for i in range(n) for j in range(n)) \
        and all(E[i][i] == 0 for i in range(n))
    divs = smith_divisors(E)
    pairs_ok = (len(divs) == n and
                all(divs[2 * k] == divs[2 * k + 1] for k in range(n // 2)))
    symp_type = [divs[2 * k] for k in range(n // 2)] if pairs_ok else None
    ok = alt_ok and pairs_ok and symp_type == list(declared)
    return ok, {
        "predicate": "polarization-type-frobenius-normal-form",
        "alternating_integral": alt_ok,
        "snf_divisors": divs,
        "snf_pair_structure_ok": pairs_ok,
        "symplectic_type": symp_type,
        "declared_elementary_divisors": list(declared),
        "match": symp_type == list(declared),
        "failing_predicate": None if ok else "polarization-type-frobenius-normal-form"}


# ----------------------------------------------------------------------------
# analytic layer: period matrix, C-linearity of rho, sign pin s (balls)
# ----------------------------------------------------------------------------

def b4_period_matrix(conjugated=False):
    """Pi = (I_4 | tau), tau = i*sqrt(5)*I_4 under the PINNED embedding
    rho -> +i*sqrt(5) (basis per factor {1, rho}, lattice O_K).  The N5
    control conjugates the complex structure: Pi -> conj(Pi)."""
    s5 = arb(5).sqrt()
    z = acb(0, s5) if not conjugated else acb(0, -s5)
    Pi = [[acb(0)] * 8 for _ in range(4)]
    for k in range(4):
        Pi[k][2 * k] = acb(1)
        Pi[k][2 * k + 1] = z
    return Pi


def analytic_rep_receipts(Pi, M):
    """From the raw input (Pi 4x8 balls, M 8x8 integers): compute the analytic
    representation A with A*Pi = Pi*M, verify C-linearity on ALL columns,
    A^2 = -5, and extract the holomorphic eigenvalue sign s (rho ~ s*i*sqrt5).
    Returns (ok, s, receipt); s = 0 when indeterminate."""
    PiM = [[sum(Pi[i][t] * M[t][j] for t in range(8)) for j in range(8)]
           for i in range(4)]
    # first 4 columns of Pi form P1; here P1 = permutation-selected identity-like
    # block (B4 layout: columns 2k are the '1' generators).  Solve A from those.
    idx1 = [2 * k for k in range(4)]
    idx2 = [2 * k + 1 for k in range(4)]
    # A * P1 = PiM[:, idx1]; with P1 = I (B4 layout) A is direct:
    A = [[PiM[i][j] for j in idx1] for i in range(4)]
    # verify on the OTHER 4 columns: A * Pi[:, idx2] must overlap PiM[:, idx2]
    lin_ok = True
    for j in idx2:
        for i in range(4):
            lhs = sum(A[i][t] * Pi[t][j] for t in range(4))
            if not bool((lhs - PiM[i][j]).contains(acb(0))):
                lin_ok = False
    A2 = [[sum(A[i][t] * A[t][j] for t in range(4)) for j in range(4)]
          for i in range(4)]
    sq_ok = all(bool((A2[i][j] - acb(-5 if i == j else 0)).contains(acb(0)))
                for i in range(4) for j in range(4))
    # diagonal-scalar receipt + sign
    s5 = arb(5).sqrt()
    diag_ok = all(bool(A[i][j].contains(acb(0))) for i in range(4)
                  for j in range(4) if i != j)
    re0_ok = all(bool(A[k][k].real.contains(arb(0))) for k in range(4))
    s = 0
    if diag_ok and re0_ok:
        if all(bool(A[k][k].imag > 0) and bool((A[k][k].imag - s5).contains(arb(0)))
               for k in range(4)):
            s = 1
        elif all(bool(A[k][k].imag < 0) and bool((A[k][k].imag + s5).contains(arb(0)))
                 for k in range(4)):
            s = -1
    ok = lin_ok and sq_ok and s != 0
    return ok, s, {
        "predicates": [
            {"name": "rho-C-linearity(A*Pi=Pi*M)", "ok": bool(lin_ok)},
            {"name": "rho-analytic-square-eq-minus5", "ok": bool(sq_ok)},
            {"name": "rho-holomorphic-eigenvalue-scalar", "ok": bool(diag_ok and re0_ok)}],
        "holomorphic_eigenvalue": "s*i*sqrt(5), s = %+d" % s if s else "INDETERMINATE",
        "failing_predicate": None if ok else (
            "rho-C-linearity(A*Pi=Pi*M)" if not lin_ok else
            "rho-analytic-square-eq-minus5" if not sq_ok else
            "rho-holomorphic-eigenvalue-scalar")}


def sign_pin_receipt(s, declared_s1):
    """Q6 sign-pin verification (VERIFY/REFUSE only; the single conversion site
    stays abcount_s0.convert_sign_convention).  P lane requires s = +1."""
    want = 1 if declared_s1 == "P" else -1
    ok = (s == want)
    return ok, {
        "predicate": "sign-convention-pin(s1)",
        "declared_s1": declared_s1,
        "required_holomorphic_sign": want,
        "input_holomorphic_sign": s,
        "verdict": "OK" if ok else "FAIL-SIGN-PIN",
        "law": "verification only; NO conversion here (single site = abcount_s0.convert_sign_convention)"}


# ----------------------------------------------------------------------------
# gp helper (Path-1 side only; NEVER used by Path-2 desk routes)
# ----------------------------------------------------------------------------

def gp_run(script, timeout=300):
    try:
        r = subprocess.run(["gp", "-q", "-f"], input=script + "\nquit()\n",
                           capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError:
        raise RuntimeError(
            "PARI/GP binary `gp` not found on PATH — the S3 battery's "
            "B4/B4s/B5 truth-side legs drive it directly (named refusal; "
            "see GUIDE.md Dependencies)")
    if r.returncode != 0:
        raise RuntimeError("gp failed rc=%d: %s" % (r.returncode, r.stderr[:500]))
    return [ln for ln in r.stdout.strip().splitlines() if ln]


# ----------------------------------------------------------------------------
# B4 / B4s desk data (exact residue bookkeeping)
# ----------------------------------------------------------------------------

J_A, J_B = 632000, 282880          # j(E20) = J_A + J_B*sqrt5


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return 1 if r == 1 else -1


def sqrt5_residues(p):
    return sorted(x for x in range(p) if (x * x - 5) % p == 0)


def e20_j_residues(p):
    """{sqrt5_residue: j mod P} for the two degree-1 primes of F above p."""
    return {r: (J_A + J_B * r) % p for r in sqrt5_residues(p)}


def ellap_e20_at(p, jres):
    """Path 1 per-prime a_P.  Standard pinned model where nonsingular
    (j(1728-j) != 0 mod p); else the ellfromj model with the SEEDS
    twist-invariance license (inert-in-K primes only).  Returns (a_P, receipt)."""
    degenerate = (jres * (1728 - jres)) % p == 0
    if not degenerate:
        out = gp_run(
            "j=Mod(%d,%d); A=3*j*(1728-j); B=2*j*(1728-j)^2;"
            "E=ellinit([A,B]); print(ellap(E),\" \",lift(E.j));" % (jres, p))
        ap_s, jchk = out[0].split()
        assert int(jchk) == jres, "model j mismatch"
        return int(ap_s), {"model": "pinned standard model reduced at P",
                           "j_mod_P": jres, "ellap": int(ap_s)}
    # degenerate: only legal at inert-in-K primes (a_P twist-invariant there)
    if legendre(-5, p) != -1:
        raise RuntimeError("degenerate-model at split prime: no licensed route")
    out = gp_run("E=ellinit(ellfromj(Mod(%d,%d))); print(ellap(E),\" \",lift(E.j));"
                 % (jres, p))
    ap_s, jchk = out[0].split()
    assert int(jchk) == jres
    return int(ap_s), {
        "model": "ellfromj (standard model singular: j(1728-j) = 0 mod p)",
        "twist_invariance_license":
            "p inert in K => supersingular => a_P = 0 for EVERY twist "
            "(SEEDS B4.model_degeneracy_note + derivation_p11_twist_invariance)",
        "j_mod_P": jres, "ellap": int(ap_s)}


def cornacchia_a2_5b2(p):
    """Exact: the (a, b), a, b > 0, with a^2 + 5 b^2 = p (p split in K,
    principal); brute force — desk arithmetic, no PARI."""
    for b in range(1, math.isqrt(p // 5) + 1):
        r = p - 5 * b * b
        if r > 0:
            a = math.isqrt(r)
            if a * a == r:
                return a, b
    raise RuntimeError("no representation a^2+5b^2 = %d" % p)


def hecke_ok_eval(p, u):
    """Path 2 — hecke_desk O_K extension: exact Z[sqrt(-5)] evaluation of the
    Frobenius data of E20 at the primes of F above p.  NO PARI.
    inert-in-K: a_P = 0 closed form.  split (29-class): Frob_P = u*pi at the
    PINNED K-prime, u*pibar at the conjugate; output (a_P, b_P) per prime in
    the pinned basis {1, rho}."""
    chi = legendre(-5, p)
    if chi == -1:
        return {"case": "inert-in-K", "a_P": 0,
                "receipt": {
                    "euler_criterion": "(-5|%d) = -1" % p,
                    "closed_form": "inert => supersingular reduction => a_P = 0 mod p; "
                                   "|a_P| <= 2*sqrt(p) < p => a_P = 0 exactly",
                    "pari_free": True}}
    a, b = cornacchia_a2_5b2(p)
    # pinned K-prime contains pi = a + b*rho (SEEDS B4s.pinned_K_prime receipt)
    # Frob_P = u*(a + b*rho) = (a_P/2) + b_P*rho with a_P = 2*u*a, b_P = u*b
    return {"case": "split-in-K", "pi": [a, b],
            "at_pinned": {"a_P": 2 * u * a, "b_P": u * b,
                          "Frob": "u*(%d + %d*rho)" % (a, b)},
            "at_conjugate": {"a_P": 2 * u * a, "b_P": -u * b,
                             "Frob": "u*(%d - %d*rho)" % (a, b)},
            "receipt": {"cornacchia": "%d^2 + 5*%d^2 = %d (exact)" % (a, b, p),
                        "unit_pin_u": u,
                        "norm_checks": {"N(pi-1)": (a - 1) ** 2 + 5 * b * b,
                                        "N(pi+1)": (a + 1) ** 2 + 5 * b * b},
                        "pari_free": True}}


# ----------------------------------------------------------------------------
# L-factor assembly for B4/B4s (exact integer polynomials)
# ----------------------------------------------------------------------------

def poly_mul(f, g):
    h = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            h[i + j] += a * b
    return h


def poly_pow(f, n):
    r = [1]
    for _ in range(n):
        r = poly_mul(r, f)
    return r


def induced_lp_b4(a_P_by_prime, p):
    """L_p(B4/F, T) = prod over P|p of (1 - a_P T + p T^2)^4 (all f(P) = 1);
    ascending coeffs; plus the count of record prod (p + 1 - a_P)^4."""
    L = [1]
    count = 1
    for a in a_P_by_prime:
        L = poly_mul(L, poly_pow([1, -a, p], 4))
        count *= (p + 1 - a) ** 4
    return L, count


# ----------------------------------------------------------------------------
# B4s Path 1: Frobenius matrix on E20[3] over F_29 + torsion congruence
# ----------------------------------------------------------------------------

def frobenius_matrix_on_E3(p, jres):
    """Explicit 2x2 Frobenius matrix on E[3] (basis from the two rational
    x-roots of the 3-division polynomial; points over F_p^2), via gp.
    Returns (F (2x2 int mod 3), receipt)."""
    script = """
p=%d; j=Mod(%d,p); A=3*j*(1728-j); B=2*j*(1728-j)^2;
E=ellinit([A,B]); d3=elldivpol(E,3);
fa=factor(d3); rs=[]; for(i=1,matsize(fa)[1], if(poldegree(fa[i,1])==1, rs=concat(rs,[-polcoef(fa[i,1],0)/polcoef(fa[i,1],1)])));
if(#rs!=2, print("ERR nroots ",#rs));
t=ffgen([p,2],'t);
EE=ellinit([A,B]*t^0);
x1=rs[1]*t^0; y1=sqrt(x1^3+A*x1+B);
x2=rs[2]*t^0; y2=sqrt(x2^3+A*x2+B);
P1=[x1,y1]; P2=[x2,y2];
ok1=(ellmul(EE,P1,3)==[0]); ok2=(ellmul(EE,P2,3)==[0]);
fr(P)=[P[1]^p,P[2]^p];
cf(Q)=my(res=[-1,-1]); for(a=0,2, for(b=0,2, if(elladd(EE,ellmul(EE,P1,a),ellmul(EE,P2,b))==Q, res=[a,b]))); res;
c1=cf(fr(P1)); c2=cf(fr(P2));
r1=(y1^p==y1); \\\\ P1 rational over F_p iff y1 fixed
print(ok1," ",ok2," ",c1[1]," ",c1[2]," ",c2[1]," ",c2[2]," ",lift(rs[1])," ",lift(rs[2])," ",r1);
""" % (p, jres)
    out = gp_run(script)
    tok = out[-1].split()
    ok1, ok2 = tok[0] == "1", tok[1] == "1"
    a, b, c, d = (int(tok[2]), int(tok[3]), int(tok[4]), int(tok[5]))
    F = [[a, c], [b, d]]
    return F, {"order3_checks": [ok1, ok2],
               "psi3_rational_xroots": [tok[6], tok[7]],
               "basis": "P1 = (x1, y1), P2 = (x2, y2) over F_%d^2" % p,
               "P1_rational_over_Fp": tok[8] == "1",
               "frobenius_matrix_columns_FP1_FP2": F}


def torsion_congruence_receipt(F, a_P, p, u, b_pin, b_conj, m=3):
    """SEEDS B4s.torsion_check: R := 2^{-1} (u*F - 3I) mod m must satisfy the
    exact rho relations; charpoly of F ties to the ellap-side a_P; the claimed
    per-prime b_P signs must be coherent (b_pin = -b_conj, |b| = 2)."""
    inv2 = pow(2, -1, m)
    R = [[(inv2 * (u * F[i][j] - (3 if i == j else 0))) % m for j in range(2)]
         for i in range(2)]
    trF = (F[0][0] + F[1][1]) % m
    detF = (F[0][0] * F[1][1] - F[0][1] * F[1][0]) % m
    cp_ok = (trF == a_P % m) and (detF == p % m)
    R2 = [[(R[i][0] * R[0][j] + R[i][1] * R[1][j]) % m for j in range(2)]
          for i in range(2)]
    rho_sq_ok = R2 == [[(-5) % m, 0], [0, (-5) % m]]
    trR = (R[0][0] + R[1][1]) % m
    detR = (R[0][0] * R[1][1] - R[0][1] * R[1][0]) % m
    tr_ok, det_ok = trR == 0, detR == 5 % m
    nonscalar = not (R[0][1] == 0 and R[1][0] == 0 and R[0][0] == R[1][1])
    signs_ok = (b_pin == -b_conj and abs(b_pin) == 2)
    ok = cp_ok and rho_sq_ok and tr_ok and det_ok and nonscalar and signs_ok
    return ok, {
        "m": m,
        "R_eq_halfinv_uF_minus_3": R,
        "charpoly_ties_ellap": {"tr_F_mod_m": trF, "a_P_mod_m": a_P % m,
                                "det_F_mod_m": detF, "p_mod_m": p % m, "ok": cp_ok},
        "rho_relations": {"R^2=-5I": rho_sq_ok, "tr=0": tr_ok,
                          "det=5": det_ok, "nonscalar": nonscalar},
        "claimed_bP_coherence": {"b_pinned": b_pin, "b_conjugate": b_conj,
                                 "flip_ok": signs_ok},
        "scope": "congruence-structure receipt: Frob_P = u*(3 + b_P*rho) verified on "
                 "E[3] up to the rho vs -rho orientation, which is supplied by the "
                 "ANALYTIC sign pin s (input rho/complex structure) + the Path-2 "
                 "conjugation-aware CM desk truth (committed desk derivation)",
        "failing": None if ok else "torsion-congruence"}


# ----------------------------------------------------------------------------
# B5: tower counts (pure python, NO PARI) + PARI hyperellcharpoly (Path 1)
# ----------------------------------------------------------------------------

def b5_counts_inert(p):
    """C_K: y^2 = x^5 + s x + 1 at inert p: counts over F_q (q = p^2) and
    F_{q^2}; F_q = F_p[s]/(s^2 - c), c = -5 mod p (nonresidue).  Returns
    (n1, n2, charpoly monic desc over F_q, receipt)."""
    c = (-5) % p
    assert legendre(c, p) == -1
    q = p * p
    MUL = [[0] * q for _ in range(q)]
    for a in range(q):
        a0, a1 = a % p, a // p
        row = MUL[a]
        for b in range(a, q):
            b0, b1 = b % p, b // p
            v = (a0 * b0 + c * a1 * b1) % p + p * ((a0 * b1 + a1 * b0) % p)
            row[b] = v
            MUL[b][a] = v

    def addq(a, b):
        return (a % p + b % p) % p + p * ((a // p + b // p) % p)

    S = p                      # the element s = (0, 1)
    sqq = set(MUL[a][a] for a in range(1, q))
    n1 = 1
    for x in range(q):
        x2 = MUL[x][x]; x4 = MUL[x2][x2]; x5 = MUL[x4][x]
        v = addq(addq(x5, MUL[S][x]), 1)
        n1 += 1 if v == 0 else (2 if v in sqq else 0)
    d = next(i for i in range(2, q) if i not in sqq)
    sq2 = set()
    for A in range(q):
        rowA = MUL[A]
        for B in range(q):
            sq2.add(addq(rowA[A], MUL[d][MUL[B][B]]) * q + addq(rowA[B], rowA[B]))
    n2 = 1
    for A in range(q):
        for B in range(q):
            a2 = addq(MUL[A][A], MUL[d][MUL[B][B]]); b2 = addq(MUL[A][B], MUL[A][B])
            a4 = addq(MUL[a2][a2], MUL[d][MUL[b2][b2]]); b4 = addq(MUL[a2][b2], MUL[a2][b2])
            a5 = addq(MUL[a4][A], MUL[d][MUL[b4][B]]); b5 = addq(MUL[a4][B], MUL[b4][A])
            va = addq(addq(a5, MUL[S][A]), 1); vb = addq(b5, MUL[S][B])
            n2 += 1 if va == 0 and vb == 0 else (2 if (va * q + vb) in sq2 else 0)
    t1 = q + 1 - n1
    t2 = q * q + 1 - n2
    e1 = t1
    assert (t1 * t1 - t2) % 2 == 0
    e2 = (t1 * t1 - t2) // 2
    cp = [1, -e1, e2, -q * e1, q * q]
    return n1, n2, cp, {
        "field": "F_%d[s]/(s^2 - %d), then quadratic ext by nonsquare %d" % (p, c, d),
        "derivation": "pure-python squares-set character counts + power sums (NO PARI)"}


def b5_counts_split(p, sval):
    """Split-prime reduction y^2 = x^5 + sval*x + 1 over F_p; counts over F_p,
    F_p^2 by the same pure-python method (B2 pattern)."""
    sqp = set(x * x % p for x in range(1, p))
    n1 = 1
    for x in range(p):
        v = (pow(x, 5, p) + sval * x + 1) % p
        n1 += 1 if v == 0 else (2 if v in sqp else 0)
    nr = next(i for i in range(2, p) if i not in sqp)
    mul = lambda a, b: ((a[0] * b[0] + nr * a[1] * b[1]) % p,
                        (a[0] * b[1] + a[1] * b[0]) % p)
    els = [(i % p, i // p) for i in range(p * p)]
    sq = set(mul(e, e) for e in els if e != (0, 0))
    n2 = 1
    for x in els:
        x2 = mul(x, x); x4 = mul(x2, x2); x5 = mul(x4, x)
        v = ((x5[0] + sval * x[0] + 1) % p, (x5[1] + sval * x[1]) % p)
        n2 += 1 if v == (0, 0) else (2 if v in sq else 0)
    t1 = p + 1 - n1
    t2 = p * p + 1 - n2
    e1 = t1
    e2 = (t1 * t1 - t2) // 2
    return n1, n2, [1, -e1, e2, -p * e1, p * p]


def hyperellcharpoly_ff(p, deg, sval_or_c):
    """Path 1: PARI hyperellcharpoly of y^2 = x^5 + s x + 1 over F_p (deg=1,
    s = sval) or F_p^2 (deg=2, s = sqrt(c), c = -5 mod p).  Returns monic-desc
    integer coeff list [1, c3, c2, c1, c0]."""
    if deg == 1:
        out = gp_run("print(Vec(hyperellcharpoly(Mod(1,%d)*(x^5+%d*x+1))));"
                     % (p, sval_or_c))
    else:
        out = gp_run("t=ffgen([%d,2],'t); s=sqrt(%d*t^0);"
                     "print(Vec(hyperellcharpoly(x^5+s*x+1)));" % (p, sval_or_c))
    v = [int(x) for x in out[-1].strip("[]").split(",")]
    assert len(v) == 5 and v[0] == 1
    return v


def lp_from_charpoly_T2(cp):
    """L_P(J, T^2) ascending deg-8 coeffs from monic-desc quartic charpoly
    (L ascending coeffs = the monic-desc list, then spread into T^2)."""
    L4 = cp                      # ascending L coeffs == monic-desc charpoly list
    L8 = [0] * 9
    for i, v in enumerate(L4):
        L8[2 * i] = v
    return L8


def lp_split_product(cp_a, cp_b):
    """L_(pi)(T) * L_(pibar)(T) ascending deg-8 coeffs."""
    return poly_mul(cp_a, cp_b)

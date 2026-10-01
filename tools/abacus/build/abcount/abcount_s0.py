"""abcount — S0 dimension-1 elliptic core (stage S0 of the sha-pinned spec
M3T_REGISTRATION.md sec 4).

SIGN CONVENTION (Q6, S0PRE_RECEIPT.md — present from the first lines of this module):
  GLOBAL pin of record: "P" (the r40 object of record, sig (4,2)).  "P(-1)" is
  PERMITTED as engine-native INTERNAL representation ONLY, and every artifact
  carries a mandatory `sign_convention` field.  Boundary conversion happens at a
  SINGLE site: the ONE function `convert_sign_convention` below — no per-stage
  ad-hoc convention flips anywhere else in this codebase.
  s1 = complex-structure convention ("P" | "P(-1)").  In dimension 1 the flip is
  complex conjugation of the complex structure: tau -> -conj(tau) (again in H_1),
  period balls -> conjugate balls.  s2 = K-embedding orientation: S0 has no
  O_K-action (dimension 1); the field is CARRIED per contract and declared
  "N/A-dim1", never dropped.

Q3 bindings (S0PRE_RECEIPT.md, riding every output):
  (c) every output receipt carries the ACTUAL input lattice's own invariants
      (polarization elementary divisors, dimension, exact model invariants
      c4/c6/disc/j when model provenance exists, lattice covolume ball) plus the
      declared isogeny data (here: none).
  (d) the tool NEVER infers equivalence to the pinned lattice P: field
      `pinned_lattice_claim` is emitted "NONE" unless the pinned lattice itself
      is the input or a separately-receipted stage-level isometry is supplied
      (neither exists in S0).

Verdict classes (closed list, TOTAL — registration sec 1, review item C9):
  OK | FAIL-<named-check> | UNDECIDED-PRECISION | UNDECIDED-PRECONDITION |
  STOP-WALL | FAIL-INTERNAL.

Certified-arithmetic substrate law: ball arithmetic is python-flint (Arb) —
existing certified substrate, NOT rebuilt here.  The desk tail bound below is an
independent PROOF-BACKED cross-enclosure, not a rebuild of Arb.

Desk theta tail bound (dimension-1 instance of the C7 Gaussian-comparison shape;
the full genus-4 inequality is an S1 obligation — manuals/abcount.md at S1):
  Let y = Im(tau) > 0 (the smallest — only — eigenvalue of Im tau in H_1) and
  r = |q| = exp(-pi*y) < 1, q = exp(i*pi*tau).
  theta3(0,tau) = 1 + 2*sum_{n>=1} q^{n^2};  theta4 likewise with sign (-1)^n;
  theta2(0,tau) = 2*exp(i*pi*tau/4)*sum_{n>=0} q^{n(n+1)}.
  For n >= N+1 the term ratio r^{(n+1)^2-n^2} = r^{2n+1} <= r^{2N+3} < 1, so by
  geometric comparison (each Gaussian term dominated by a geometric series):
    |theta3/theta4 tail past n=N|          <= 2 r^{(N+1)^2} / (1 - r^{2N+3}),
    |theta2 tail: 2|q^{1/4}| sum_{n>N}|    <= 2 r^{1/4} r^{(N+1)(N+2)} / (1 - r^{2N+4}).
  Bounds evaluated in Arb upper-bound arithmetic; the desk enclosure is
  [partial-sum ball] +/- [tail bound] and must OVERLAP the Arb enclosure —
  disjoint enclosures = STOP (theta mismatch).
"""

import math
import sys
import time
import resource

from flint import acb, arb, ctx, fmpq

# ----------------------------------------------------------------------------
# Q6: sign-convention constants + THE single boundary-conversion function
# ----------------------------------------------------------------------------

SIGN_S1_VALUES = ("P", "P(-1)")
GLOBAL_PIN_S1 = "P"          # pin of record; P(-1) internal-only
S2_DIM1 = "N/A-dim1"         # carried, never dropped (no K-action in S0)


def convert_sign_convention(artifact):
    """THE single receipted boundary-conversion site (Q6).

    artifact: dict with `sign_convention` {"s1":..., "s2":...} and optional
    `tau` (acb), `periods` (pair of acb).  Flips s1 P <-> P(-1) by conjugating
    the complex structure: tau -> -conj(tau); periods -> entrywise conj.
    Returns (converted_artifact, conversion_receipt).  No other function in
    this codebase may flip a convention.
    """
    s1 = artifact["sign_convention"]["s1"]
    if s1 not in SIGN_S1_VALUES:
        raise ValueError("FAIL-SIGN-PIN: unknown s1 value %r" % (s1,))
    new_s1 = "P(-1)" if s1 == "P" else "P"
    out = dict(artifact)
    out["sign_convention"] = dict(artifact["sign_convention"], s1=new_s1)
    receipt = {"site": "convert_sign_convention (single site, Q6)",
               "from_s1": s1, "to_s1": new_s1,
               "s2": artifact["sign_convention"]["s2"]}
    if artifact.get("tau") is not None:
        out["tau"] = -artifact["tau"].conjugate()
        receipt["tau_map"] = "tau -> -conj(tau)"
    if artifact.get("periods") is not None:
        out["periods"] = tuple(w.conjugate() for w in artifact["periods"])
        receipt["period_map"] = "w -> conj(w)"
    return out, receipt


# ----------------------------------------------------------------------------
# exact model invariants (Q3c: the ACTUAL input's own invariants, exact integers)
# ----------------------------------------------------------------------------

def weierstrass_invariants(a1, a2, a3, a4, a6):
    """Exact b2,b4,b6,b8,c4,c6,disc and j = c4^3/disc from the input model."""
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    num = b2 * b6 - b4 * b4
    if num % 4 != 0:
        raise ValueError("FAIL-INVARIANTS: b8 not integral")
    b8 = num // 4
    c4 = b2 * b2 - 24 * b4
    c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
    disc = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    assert 1728 * disc == c4 ** 3 - c6 ** 2, "c4^3 - c6^2 != 1728*disc"
    if disc == 0:
        raise ValueError("FAIL-INVARIANTS: singular model (disc = 0)")
    return {"b2": b2, "b4": b4, "b6": b6, "b8": b8,
            "c4": c4, "c6": c6, "disc": disc,
            "j_num": c4 ** 3, "j_den": disc}


# ----------------------------------------------------------------------------
# tau construction: period balls -> reduced tau in H_1 (SL2(Z) receipt)
# ----------------------------------------------------------------------------

def tau_from_periods(w1, w2):
    """tau = w1/w2 (swap if needed for certified Im > 0), then SL2(Z)-reduce
    toward the standard fundamental domain.  j is SL2(Z)-invariant, so theta/j
    at the reduced tau computes the lattice's j.  Returns (tau, receipt)."""
    tau = w1 / w2
    swapped = False
    if not bool(tau.imag > 0):
        tau = w2 / w1
        swapped = True
        if not bool(tau.imag > 0):
            raise ValueError("UNDECIDED-PRECISION: cannot certify Im(tau) > 0")
    a, b, c, d = 1, 0, 0, 1
    moves = 0
    while moves < 200:
        nn = int(round(float(tau.real.mid())))
        if nn != 0:
            tau = tau - nn
            a, b = a - nn * c, b - nn * d
        if bool(tau.abs_upper() < 1):
            tau = -1 / tau
            a, b, c, d = -c, -d, a, b
            moves += 1
            continue
        break
    if not bool(tau.imag > 0):
        raise ValueError("UNDECIDED-PRECISION: reduction lost Im(tau) > 0")
    return tau, {"swapped_basis": swapped, "sl2_inversions": moves,
                 "sl2_matrix_abcd": [a, b, c, d],
                 "tau_mid": str(complex(tau)),
                 "im_tau_positive_certified": True}


# ----------------------------------------------------------------------------
# theta nulls: route ARB (acb_modular, certified) + route DESK (partial sum +
# proved tail bound, C7 shape in dim 1)
# ----------------------------------------------------------------------------

def theta_nulls_arb(tau):
    """Certified theta nulls via Arb acb_modular_theta (z = 0)."""
    _t1, t2, t3, t4 = acb(0).modular_theta(tau)
    return t2, t3, t4


def theta_nulls_desk(tau, N):
    """Desk enclosure: truncation at N + PROVED tail bound (module docstring).
    Independent of Arb's acb_modular_theta code path (shares only the ball
    add/mul/exp substrate, per the substrate law)."""
    ipi = acb(0, arb.pi())
    q = (ipi * tau).exp()
    r = arb(q.abs_upper())          # certified upper bound for |q|
    if not bool(r < 1):
        raise ValueError("UNDECIDED-PRECISION: |q| < 1 not certified")
    t3 = acb(1)
    t4 = acb(1)
    for n in range(1, N + 1):
        qn2 = q ** (n * n)
        t3 = t3 + 2 * qn2
        t4 = t4 + 2 * ((-1) ** n) * qn2
    q14 = (ipi * tau / 4).exp()     # e^{i pi tau / 4}, principal by construction
    s = acb(0)
    for n in range(0, N + 1):
        s = s + q ** (n * (n + 1))
    t2 = 2 * q14 * s
    tail34 = 2 * r ** ((N + 1) ** 2) / (1 - r ** (2 * N + 3))
    r14 = arb(q14.abs_upper())
    tail2 = 2 * r14 * r ** ((N + 1) * (N + 2)) / (1 - r ** (2 * N + 4))
    def widen(ball, t):
        pad = arb("0 +/- 1") * t          # rigorous [-t_up, t_up] enclosure
        return ball + acb(pad, pad)
    return (widen(t2, tail2), widen(t3, tail34), widen(t4, tail34),
            {"cutoff_N": N, "abs_q_upper": r.str(10),
             "tail34_bound": tail34.str(6), "tail2_bound": tail2.str(6),
             "inequality": "geometric/Gaussian-comparison (module docstring; C7 dim-1 shape)"})


def balls_overlap(x, y):
    return bool((x - y).contains(acb(0)))


# ----------------------------------------------------------------------------
# j from theta nulls + recognition against the actual input's exact j
# ----------------------------------------------------------------------------

def j_from_theta(t2, t3):
    lam = (t2 / t3) ** 4
    num = (lam * lam - lam + 1) ** 3
    den = (lam * (1 - lam)) ** 2
    return 256 * num / den


def ball_contains_rational(ball, p, q):
    """Rigorous: ball contains p/q  <=>  (ball*q - p) contains 0 (q > 0 exact)."""
    assert q > 0
    return bool((ball * q - p).contains(acb(0)))


def j_recognition(j_ball, j_num, j_den):
    """Recognition receipt (held-out-gate discipline): (i) the exact rational j
    of the ACTUAL input model lies inside the certified ball; (ii) positive
    controls: nearby wrong targets (j+1, j-1, -j) must be EXCLUDED."""
    if j_den < 0:
        j_num, j_den = -j_num, -j_den
    inside = ball_contains_rational(j_ball, j_num, j_den)
    controls = {}
    for name, wnum in (("j_plus_1", j_num + j_den),
                       ("j_minus_1", j_num - j_den),
                       ("j_negated", -j_num)):
        controls[name] = ("EXCLUDED"
                         if not ball_contains_rational(j_ball, wnum, j_den)
                         else "NOT-EXCLUDED")
    ok = inside and all(v == "EXCLUDED" for v in controls.values())
    return ok, {"exact_j_in_ball": inside, "negative_controls": controls,
                "j_ball_mid": str(complex(j_ball)),
                "j_ball_rad_re": j_ball.real.rad().str(6)}


# ----------------------------------------------------------------------------
# point counts at p: route C (PARI ellap) vs route DESK (naive count, no PARI)
# ----------------------------------------------------------------------------

def ap_naive(a1, a2, a3, a4, a6, p):
    """Independent desk count over F_p (pure integer arithmetic, NO PARI):
    complete y in y^2 + (a1 x + a3) y = rhs(x); for odd p the y-count at x is
    1 + chi((a1 x + a3)^2 + 4 rhs) (chi = Euler-criterion Legendre symbol,
    chi(0) := 0 giving exactly 1 solution).  Plus the point at infinity.
    a_p = p + 1 - #E(F_p)."""
    if p % 2 == 0:
        raise ValueError("UNDECIDED-PRECONDITION: desk route coded for odd p")
    npts = 1
    for x in range(p):
        rhs = (x * x * x + a2 * x * x + a4 * x + a6) % p
        dy = ((a1 * x + a3) ** 2 + 4 * rhs) % p
        if dy == 0:
            npts += 1
        else:
            npts += 2 if pow(dy, (p - 1) // 2, p) == 1 else 0
    return p + 1 - npts, npts


def ap_pari(model, p):
    """Route C: PARI ellap via cypari2 (padic_counters row pattern)."""
    import cypari2
    pari = cypari2.Pari()
    e = pari.ellinit(list(model))
    return int(pari.ellap(e, p))


def pari_periods(model, dps):
    """Input generation: period lattice of the actual input model at declared
    dps, imported as balls with declared per-entry radius 1e-(dps-5)."""
    import cypari2
    pari = cypari2.Pari()
    pari.set_real_precision(dps + 10)
    bits = int((dps + 12) * 3.33)
    e = pari.ellinit(list(model), precision=bits)
    om = pari.ellperiods(e, precision=bits)
    rad = arb("0 +/- 1e-%d" % (dps - 5))
    def to_ball(z):
        re_s = str(pari.real(z)).replace(" E", "e").replace("E", "e")
        im_s = str(pari.imag(z)).replace(" E", "e").replace("E", "e")
        return acb(arb(re_s) + rad, arb(im_s) + rad)
    return to_ball(om[0]), to_ball(om[1]), {"per_entry_radius": "1e-%d" % (dps - 5),
                                            "source": "PARI ellperiods, dps+10 guard digits"}


# ----------------------------------------------------------------------------
# Weil-box integer isolation layer (Builds item 2, dimension-1 instance)
# ----------------------------------------------------------------------------

def weil_box_isolate_g1(enclosure, p):
    """FULL ENUMERATION of integers a with |a| <= 2 sqrt(p) intersected with the
    enclosure ball (never per-coefficient rounding).  Returns (verdict, cands):
    OK (unique) | FAIL-WEIL-BOX-EMPTY | UNDECIDED-PRECISION."""
    A = math.isqrt(4 * p)
    cands = [a for a in range(-A, A + 1) if enclosure.contains(acb(a))]
    if len(cands) == 1:
        return "OK", cands
    if len(cands) == 0:
        return "FAIL-WEIL-BOX-EMPTY", cands
    return "UNDECIDED-PRECISION", cands


# ----------------------------------------------------------------------------
# wall (S0 row: minutes wallclock, <= 2 GB, dps <= 60)
# ----------------------------------------------------------------------------

def maxrss_kib():
    """Peak RSS of this process in KiB.  getrusage ru_maxrss is KiB on Linux
    but BYTES on macOS — normalize here so the wall compares real units."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1024 if sys.platform == "darwin" else rss


class Wall:
    def __init__(self, seconds=300, mem_gb=2.0, dps=60):
        if ctx.dps > dps:
            raise RuntimeError("STOP-WALL: ctx.dps %d exceeds S0 wall %d" % (ctx.dps, dps))
        self.t0 = time.time()
        self.seconds = seconds
        self.mem_gb = mem_gb

    def check(self, where):
        el = time.time() - self.t0
        rss_gb = maxrss_kib() / (1024 ** 2)
        if el > self.seconds:
            raise RuntimeError("STOP-WALL: wallclock %.1fs > %ds at %s"
                               % (el, self.seconds, where))
        if rss_gb > self.mem_gb:
            raise RuntimeError("STOP-WALL: rss %.2f GB > %.1f GB at %s"
                               % (rss_gb, self.mem_gb, where))
        return {"elapsed_s": round(el, 3), "rss_gb": round(rss_gb, 4)}

"""abcount — S1 products extension (stage S1 of the sha-pinned spec
M3T_REGISTRATION.md sec 4).

EXTENDS abcount_s0 (S0 module reused verbatim, banked sha 5e85503f...); the Q6
single conversion site REMAINS abcount_s0.convert_sign_convention — this module
introduces NO second conversion site (S0PRE_RECEIPT.md binding).

S1 scope (pinned spec row): battery member B1 (block-diagonal genus-4 tau);
ONE acb_theta genus-4 call cost measured FIRST (timed pilot, one member, one
prime); independent tail-bound cross-check on that one call (Q1 / review item C7 —
inequality stated and proved in manuals/abcount.md, quoted in the S1 receipt;
pinned criterion: acb_theta claimed radius <= desk tail bound at matched
truncation, ANY violation = STOP).

C7 TAIL INEQUALITY (manuals/abcount.md sec 1, quoted here):
  For tau in H_g, lambda = lambda_min(Im tau), any half-integer characteristic
  [alpha;beta] (alpha, beta in {0,1/2}^g), and the box partial sum S_N over
  ||n||_inf <= N:
    |theta[alpha;beta](0,tau) - S_N| <= TailBound(g, N, lambda)
        = g * S_tail(N, lambda) * M(lambda)^(g-1),
    S_tail(N, lambda) = 2 exp(-pi lambda (N+1/2)^2) / (1 - exp(-pi lambda (2N+1))),
    M(lambda)         = 1 + 2 exp(-pi lambda / 4) / (1 - exp(-pi lambda)).
  Proof: Rayleigh (term modulus <= exp(-pi lambda ||n+alpha||^2)) + union bound
  over the coordinate exceeding N + separable 1-D geometric-comparison lemmas
  L1/L2 — manuals/abcount.md. Block-diagonal tau: lambda_min = min_k Im tau_k.

Genus-4 theta engine: the EXISTING acb_theta wrapper
Software/Eichler/src/siegel.jl (`siegel_theta_all` = FLINT acb_theta_all,
Kieffer 2023) — registration Calls list; values cross the bridge as EXACT
dyadic mantissa/exponent integer pairs in both directions (no decimal
rounding anywhere on the bridge).

Route D (B1): per-factor a_p -> exact integer product assembly of the four
degree-2 L-factors into the degree-8 L_p (abcount-internal, NO PARI in the
assembly layer) -> Weil-box FULL-ENUMERATION isolation (registration C3).
"""

import itertools
import math
import os
import shutil
import subprocess
import time

from flint import acb, arb, ctx, fmpz

import abcount_s0 as ab

# ABACUS_JULIA overrides the julia executable (an absolute path or a name
# looked up on PATH); the default is whatever `julia` is on PATH.
JULIA = os.environ.get("ABACUS_JULIA", "julia")
# Public port: the Eichler julia project (siegel.jl wrapper) is located via
# the environment; the julia-call sites refuse loudly when it is unset.
EICHLER_PROJECT = os.environ.get("ABACUS_EICHLER_PROJECT", "")
SIEGEL_JL = (os.path.join(EICHLER_PROJECT, "src", "siegel.jl")
             if EICHLER_PROJECT else "")


def require_eichler():
    """Named refusal: the Eichler siegel.jl wrapper (used as-is, never edited)
    must be pointed to."""
    if not EICHLER_PROJECT:
        raise RuntimeError(
            "ABACUS_EICHLER_PROJECT unset — path to the Eichler julia "
            "project (src/siegel.jl wrapper) is required")


def require_julia():
    """Named refusal: a julia executable must resolve — ABACUS_JULIA if set,
    else `julia` on PATH."""
    if shutil.which(JULIA) is None:
        raise RuntimeError(
            "julia executable %r does not resolve — set ABACUS_JULIA to a "
            "julia executable (an absolute path or a name on PATH) or put "
            "julia on PATH" % JULIA)
THETA_BRIDGE_JL = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "theta_g4_bridge.jl")

# ----------------------------------------------------------------------------
# FLINT characteristic bookkeeping, genus 4
# ----------------------------------------------------------------------------

G = 4


def char_decode(idx):
    """FLINT 2g-bit encoding a_1..a_g b_1..b_g, MSB first (siegel.jl docs)."""
    bits = [(idx >> (2 * G - 1 - j)) & 1 for j in range(2 * G)]
    return tuple(bits[:G]), tuple(bits[G:])


def char_is_even(idx):
    a, b = char_decode(idx)
    return sum(x * y for x, y in zip(a, b)) % 2 == 0


def genus1_theta_char(a, b, t2, t3, t4):
    """theta[a;b](0,tau) for a,b in {0,1} from the S0 genus-1 nulls:
    (0,0)=theta3, (0,1)=theta4, (1,0)=theta2, (1,1)=theta1(0)=0 EXACT."""
    if (a, b) == (0, 0):
        return t3
    if (a, b) == (0, 1):
        return t4
    if (a, b) == (1, 0):
        return t2
    return acb(0)


# ----------------------------------------------------------------------------
# C7 desk side: 1-D characteristic partial sums + the PROVED TailBound
# ----------------------------------------------------------------------------

def theta1d_char_partial(tau_k, a, b, N):
    """S^(k)_N(a,b) = sum_{n=-N..N} exp(pi*i*[(n+a/2)^2 tau_k + (n+a/2) b]).
    Ball arithmetic; independent of acb_theta's code path (substrate law)."""
    ipi = acb(0, arb.pi())
    s = acb(0)
    for n in range(-N, N + 1):
        half = acb(2 * n + a) / 2
        s = s + (ipi * (half * half * tau_k + half * b)).exp()
    return s


def desk_theta_g4_all(taus, N):
    """Box partial sums for ALL 256 genus-4 characteristics of the
    block-diagonal tau = diag(taus): the box sum over ||n||_inf <= N
    factorizes EXACTLY into the product of 1-D characteristic sums."""
    per_factor = []
    for tau_k in taus:
        per_factor.append({(a, b): theta1d_char_partial(tau_k, a, b, N)
                           for a in (0, 1) for b in (0, 1)})
    out = []
    for idx in range(1 << (2 * G)):
        a, b = char_decode(idx)
        v = acb(1)
        for k in range(G):
            v = v * per_factor[k][(a[k], b[k])]
        out.append(v)
    return out


def tail_bound_c7(g, N, lam_lb):
    """TailBound(g, N, lambda) of the manual, evaluated in ball arithmetic on
    a certified LOWER bound lam_lb of lambda_min(Im tau); the returned arb is
    a certified upper bound (TailBound is decreasing in lambda).
    Returns (bound_upper: arb exact, receipt)."""
    lam = arb(lam_lb)
    if not bool(lam > 0):
        raise ValueError("UNDECIDED-PRECONDITION: lambda_min lower bound not positive")
    pi = arb.pi()
    s_tail = 2 * (-pi * lam * arb(2 * N + 1) ** 2 / 4).exp() \
        / (1 - (-pi * lam * (2 * N + 1)).exp())
    m_full = 1 + 2 * (-pi * lam / 4).exp() / (1 - (-pi * lam).exp())
    bound = g * s_tail * m_full ** (g - 1)
    ub = bound.upper()
    receipt = {
        "inequality_quoted": ("|theta[a;b](0,tau) - S_N| <= g * S_tail(N,lambda) * M(lambda)^(g-1); "
                              "S_tail = 2 exp(-pi lambda (N+1/2)^2)/(1 - exp(-pi lambda (2N+1))); "
                              "M = 1 + 2 exp(-pi lambda/4)/(1 - exp(-pi lambda)); "
                              "lambda = lambda_min(Im tau) [manuals/abcount.md sec 1, PROVED: "
                              "Rayleigh + union bound + geometric-comparison lemmas L1/L2]"),
        "g": g, "N": N, "lambda_lower_bound": lam_lb.str(20),
        "tail_bound_upper": ub.str(12),
    }
    return ub, receipt


# ----------------------------------------------------------------------------
# exact dyadic bridge helpers (python <-> julia, no decimal rounding)
# ----------------------------------------------------------------------------

def _arb_exact_from_man_exp(man, e):
    """Exact arb = man * 2^e (integers). Temporarily boosts precision so the
    fmpz -> arb conversion is exact."""
    man = int(man)
    e = int(e)
    if man == 0:
        return arb(0)
    old = ctx.prec
    try:
        ctx.prec = max(old, man.bit_length() + 8)
        x = arb(fmpz(man))
        y = x * (arb(2) ** e if e >= 0 else 1 / (arb(2) ** (-e)))
    finally:
        ctx.prec = old
    # exactness receipt: dyadic scaling of an exact integer stays exact
    if not y.is_exact():
        raise RuntimeError("FAIL-INTERNAL: dyadic import not exact")
    return y


def _man_exp_of_mid(x):
    """(man, exp) of the EXACT midpoint of arb x."""
    m = x.mid()
    if m.is_zero():
        return 0, 0
    man, e = m.man_exp()
    return int(man), int(e)


def dyadic_truncate_64(x):
    """Deterministic 64-bit dyadic truncation of the midpoint of arb x
    (SEEDS.json S1_C7_pilot pilot_input_rule)."""
    man, e = _man_exp_of_mid(x)
    if man == 0:
        return 0, 0
    bits = abs(man).bit_length()
    if bits > 64:
        sh = bits - 64
        man = man >> sh if man > 0 else -((-man) >> sh)
        e += sh
    return man, e


def julia_theta_all_g4(tau_entries, prec, workdir, tag):
    """ONE genus-4 acb_theta call through the Eichler wrapper siegel.jl
    (used as-is, never edited).

    tau_entries: list of 4 dicts {re_man, re_exp, im_man, im_exp,
    rad_man, rad_exp} (exact dyadics; rad added to both parts via add_error).
    Returns dict: values (256 acb balls, enclosing the wrapper output
    exactly), claimed_radii (256 arb exact upper radii = max of re/im radius),
    tau_echo entries, call_seconds, load_seconds, maxrss_gb, prec.
    """
    inp = os.path.join(workdir, "theta_in_%s.txt" % tag)
    outp = os.path.join(workdir, "theta_out_%s.txt" % tag)
    with open(inp, "w") as f:
        f.write("prec %d\n" % prec)
        f.write("g %d\n" % G)
        for e in tau_entries:
            f.write("block %d %d %d %d %d %d\n" % (
                e["re_man"], e["re_exp"], e["im_man"], e["im_exp"],
                e["rad_man"], e["rad_exp"]))
    t0 = time.time()
    require_eichler()
    require_julia()
    r = subprocess.run(
        [JULIA, "--project=%s" % EICHLER_PROJECT, THETA_BRIDGE_JL, inp, outp],
        capture_output=True, text=True, timeout=3300,
        env=dict(os.environ, ABACUS_SIEGEL_JL=SIEGEL_JL))
    wall_subprocess = time.time() - t0
    if r.returncode != 0:
        raise RuntimeError("FAIL-INTERNAL: julia bridge rc=%d stderr=%s"
                           % (r.returncode, r.stderr[-2000:]))
    values, claimed, echo = [None] * 256, [None] * 256, []
    call_ns = load_ns = maxrss = None
    with open(outp) as f:
        for line in f:
            w = line.split()
            if not w:
                continue
            if w[0] == "theta":
                idx = int(w[1])
                re_m = _arb_exact_from_man_exp(w[2], w[3])
                im_m = _arb_exact_from_man_exp(w[4], w[5])
                re_r = _arb_exact_from_man_exp(w[6], w[7])
                im_r = _arb_exact_from_man_exp(w[8], w[9])
                pad_r = arb("0 +/- 1") * re_r
                pad_i = arb("0 +/- 1") * im_r
                values[idx] = acb(re_m + pad_r, im_m + pad_i)
                claimed[idx] = re_r.max(im_r)
            elif w[0] == "tau_echo":
                echo.append([int(v) for v in w[1:]])
            elif w[0] == "call_ns":
                call_ns = int(w[1])
            elif w[0] == "load_ns":
                load_ns = int(w[1])
            elif w[0] == "maxrss_bytes":
                maxrss = int(w[1])
    if any(v is None for v in values) or call_ns is None or maxrss is None:
        raise RuntimeError("FAIL-INTERNAL: julia bridge output incomplete")
    return {"values": values, "claimed_radii": claimed, "tau_echo": echo,
            "call_seconds": call_ns / 1e9, "load_seconds": load_ns / 1e9,
            "subprocess_wall_seconds": wall_subprocess,
            "maxrss_gb": maxrss / (1024 ** 3), "prec": prec}


# ----------------------------------------------------------------------------
# route-D assembly (abcount-internal; exact integers; NO PARI in this layer)
# ----------------------------------------------------------------------------

def assemble_lp_deg8(a_list, p):
    """Exact integer product of the four degree-2 factors
    L(T) = prod_i (1 - a_i T + p T^2) -> coefficients c_0..c_8 (c_0 = 1).
    Charpoly convention: P(T) = T^8 L(1/T), so the Weil-box tuple is
    (c_1, c_2, c_3, c_4) and #B(F_p) = P(1) = sum c_k."""
    c = [1]
    for a in a_list:
        f = [1, -int(a), int(p)]
        c = [sum(c[j] * f[k - j] for j in range(len(c)) if 0 <= k - j < 3)
             for k in range(len(c) + 2)]
    assert len(c) == 9 and c[0] == 1
    return c


def functional_equation_receipt(c, p):
    """C3 inequality set part (ii): a_(8-k) = p^(4-k) a_k for k = 0..3."""
    rows = {"k=%d" % k: (c[8 - k] == p ** (4 - k) * c[k]) for k in range(4)}
    return all(rows.values()), rows


def weil_box_isolate_deg8(coeff_enclosures, p):
    """Registration C3 isolation semantics: FULL ENUMERATION of integer
    tuples (a_1..a_4) in the Weil box |a_k| <= C(8,k) p^(k/2) intersected
    with the enclosure, NEVER per-coefficient rounding.
    coeff_enclosures: 4 acb balls (zero-radius for integer-exact route D).
    Returns (verdict, tuples, receipt)."""
    per_k = []
    for k in range(1, 5):
        A = math.isqrt(math.comb(8, k) ** 2 * p ** k)   # floor(C(8,k) p^(k/2))
        ball = coeff_enclosures[k - 1]
        cands = [a for a in range(-A, A + 1) if ball.contains(acb(a))]
        per_k.append({"k": k, "box_halfwidth": A, "n_candidates": len(cands),
                      "cands": cands})
        if len(cands) > 4000:
            return ("UNDECIDED-PRECISION", [],
                    {"reason": "coefficient a_%d enclosure admits %d integers"
                     % (k, len(cands)),
                     "needed_dps_estimate": _needed_dps(ball),
                     "per_k": [{kk: vv for kk, vv in d.items() if kk != "cands"}
                               for d in per_k]})
    tuples = list(itertools.product(*[d["cands"] for d in per_k]))
    receipt = {"per_k": [{kk: vv for kk, vv in d.items() if kk != "cands"}
                         for d in per_k],
               "n_tuples_enumerated": len(tuples),
               "semantics": "FULL enumeration box-intersect-enclosure (C3); never rounding"}
    if len(tuples) == 1:
        return "OK", tuples, receipt
    if len(tuples) == 0:
        return "FAIL-WEIL-BOX-EMPTY", tuples, receipt
    receipt["needed_dps_estimate"] = max(_needed_dps(b) for b in coeff_enclosures)
    return "UNDECIDED-PRECISION", tuples, receipt


def _needed_dps(ball):
    rad = ball.real.rad().max(ball.imag.rad())
    if bool(rad == 0):
        return 0
    return int(-float(rad.log() / arb(10).log()).real) + 5 if bool(rad > 0) else 0


def count_from_charpoly(c):
    """#B(F_p) = P(1) = sum_k c_k (monic charpoly P(T) = T^8 L(1/T))."""
    return sum(c)

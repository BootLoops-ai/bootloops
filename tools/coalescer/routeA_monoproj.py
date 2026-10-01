#!/usr/bin/env python3
r"""
Route A — MONODROMY-PROJECTION extraction of c_alpha for the L-loop banana
non-unipotent threshold, GENERIC ORDER.  Same code path runs the K3 positive
control (must reproduce c_{3/2} = -sqrt(3)/(36 pi)) and the CY3 target.

Spectral projector (the K3 trick lifted): with threshold monodromy M_0,
distinct eigenvalues E = {1} U {e^{2 pi i alpha} : fractional alpha}, the
projector onto eigenvalue lambda is

    P_lambda = (M_0 - I)^nilp * prod_{mu in E\{1,lambda}} (M_0 - mu I)
               -------------------------------------------------------
               (lambda - 1)^nilp * prod_{mu in E\{1,lambda}} (lambda - mu)

where nilp = rank of the unipotent block (>= its nilpotent index).  (M-I)^nilp
ANNIHILATES the unipotent log-tower exactly; the remaining linear factors
separate the semisimple eigenspaces.

Two doors onto one kernel:
  run_operator(tag, Pj, order, seed_a, ...)  -- GENERIC: any Fuchsian operator
      L = sum_j P_j(z) theta_z^j given as exact polynomial coefficients, plus the
      exact (rational) Taylor coefficients a_n of the holomorphic solution
      f(z) = sum_n a_n z^n at z = 0 whose connection coefficients onto the
      fractional-exponent Frobenius branches at z = infinity are wanted.
      load_operator_json() reads such an operator + seed from a JSON file (the
      layout dump_pf_data.py writes; schema in GUIDE.md INPUTS).
  run(tag, msq, ...)  -- the banana worked example: builds the L-loop banana
      Picard-Fuchs operator for the mass tuple msq (find_minimal_op) and the
      BFKNS series seed, then calls run_operator.  Behaviour unchanged.
"""
import sys, os, json, time, argparse
from fractions import Fraction
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pf_engine import (find_minimal_op, indicial_threshold, theta_t_terms,
                       frob_power_branch, phi_theta_state, wms_series_mp,
                       build_ode_s, transport_adaptive, digits, z)


def loop_waypoints(r, npts=32, ccw=True):
    sgn = 1 if ccw else -1
    return [r * mp.expj(sgn * 2 * mp.pi * k / npts) for k in range(npts + 1)]


def matpow(M, k):
    R = mp.eye(M.rows)
    for _ in range(k):
        R = R * M
    return R


# ---------------------------------------------------------------------------
# Generic front door: exact seed coefficients -> theta_z-state; operator from
# coefficient lists; JSON loader.
# ---------------------------------------------------------------------------
def exact_coeff(c):
    """One user-given seed coefficient -> Fraction.  Accepts int, Fraction,
    sympy Rational, 'p/q' or integer strings, and [p, q] pairs."""
    if isinstance(c, bool):
        raise TypeError("boolean is not a coefficient")
    if isinstance(c, int):
        return Fraction(c)
    if isinstance(c, Fraction):
        return c
    if isinstance(c, (list, tuple)) and len(c) == 2:
        return Fraction(int(c[0]), int(c[1]))
    if isinstance(c, str):
        return Fraction(c.strip())
    if isinstance(c, sp.Rational):
        return Fraction(int(c.p), int(c.q))
    if isinstance(c, float):
        if c != int(c):
            raise TypeError(f"inexact float seed coefficient {c!r}: give an int, 'p/q' string or [p, q] pair")
        return Fraction(int(c))
    raise TypeError(f"unsupported seed coefficient type {type(c).__name__}: {c!r}")


def _coeff_mp(c):
    """Seed coefficient -> mp number at the CURRENT mp.mp.dps.  mp numbers pass
    through unchanged (the banana door feeds wms_series_mp floats, which keeps
    its arithmetic bit-identical); everything else goes through exact_coeff."""
    if isinstance(c, (mp.mpf, mp.mpc)):
        return c
    f = exact_coeff(c)
    if f.denominator == 1:
        return mp.mpf(f.numerator)
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)


def make_series_state(seed_a, order, Nser=None):
    """theta_z-state of the seed solution f(z) = sum_n a_n z^n at z = -1/s.

    seed_a : the coefficients a_0..a_N (exact: int / Fraction / 'p/q' / [p, q];
             or mp numbers).  Nser truncates to a_0..a_Nser (default: all).
    Returns state(s) -> column [theta_z^k f]_{k=0..order-1}, theta_z^k z^n = n^k z^n
    (the same loop as pf_engine.make_period_state, which is this function
    applied to the BFKNS banana series)."""
    if Nser is None:
        Nser = len(seed_a) - 1
    if Nser < 1 or Nser > len(seed_a) - 1:
        raise ValueError(f"Nser={Nser} out of range for a seed of {len(seed_a)} coefficients")
    a_mp = [_coeff_mp(c) for c in seed_a[:Nser + 1]]

    def state(sv):
        sv = mp.mpc(sv)
        zz = mp.mpf(-1) / sv
        Y = mp.matrix(order, 1)
        zn = mp.mpf(1)
        for n in range(len(a_mp)):
            term = a_mp[n] * zn
            nk = mp.mpf(1)
            for k in range(order):
                Y[k, 0] += nk * term
                nk *= n
            zn *= zz
        return Y
    return state


def op_from_coeffs(theta_coeffs):
    """{j: [c_j0, c_j1, ...]} or [[...], ...] (P_j(z) integer/rational coefficients,
    LOW -> HIGH power of z; operator L = sum_j P_j(z) theta_z^j) -> (Pj dict of sympy
    polynomials in pf_engine.z with a common integer normalization, order)."""
    if isinstance(theta_coeffs, dict):
        rows = {int(k): v for k, v in theta_coeffs.items()}
    else:
        rows = {j: v for j, v in enumerate(theta_coeffs)}
    order = max(rows)
    fr = {j: [exact_coeff(c) for c in rows.get(j, [])] for j in range(order + 1)}
    L = 1
    for j in fr:
        for c in fr[j]:
            L = L * c.denominator // sp.igcd(L, c.denominator)
    Pj = {}
    for j in range(order + 1):
        Pj[j] = sp.expand(sum((sp.Integer(int(c * L)) * z ** k for k, c in enumerate(fr[j])), sp.Integer(0)))
    if Pj[order] == 0:
        raise ValueError(f"leading coefficient P_{order}(z) is zero")
    return Pj, order


def op_degz(Pj):
    return max(sp.Poly(p, z).degree() for p in Pj.values() if p != 0)


OP_JSON_KEYS = ("theta_form_Pj", "theta_coeffs", "Pj")
SEED_JSON_KEYS = ("bfkns_a", "seed")


def load_operator_json(path):
    """Read an operator + seed JSON (schema: GUIDE.md INPUTS; dump_pf_data.py writes it).
    Returns dict(Pj, order, degz, seed (list of Fraction), alpha (str or None), label,
    seed_radius (float or None), knobs (dict of optional sb/sb2/lift/sdec/nthr/nser))."""
    with open(path, encoding="utf-8") as fh:
        J = json.load(fh)
    key = next((k for k in OP_JSON_KEYS if k in J), None)
    if key is None:
        raise ValueError(f"{path}: no operator: give one of {OP_JSON_KEYS} "
                         "(P_j(z) coefficient lists, low->high power of z, j = 0..order)")
    Pj, order = op_from_coeffs(J[key])
    if "order" in J and int(J["order"]) != order:
        raise ValueError(f"{path}: 'order'={J['order']} but the operator table has j = 0..{order}")
    skey = next((k for k in SEED_JSON_KEYS if k in J), None)
    if skey is None:
        raise ValueError(f"{path}: no seed: give one of {SEED_JSON_KEYS} "
                         "(exact Taylor coefficients a_0, a_1, ... of the holomorphic solution at z = 0)")
    seed = [exact_coeff(c) for c in J[skey]]
    if len(seed) < 2:
        raise ValueError(f"{path}: the seed needs at least a_0, a_1")
    alpha = J.get("alpha")
    if alpha is not None:
        alpha = str(sp.Rational(str(alpha)))
    label = J.get("label") or J.get("tag") or (f"msq={tuple(J['msq'])}" if "msq" in J else os.path.basename(path))
    R = J.get("seed_radius")
    knobs = {k: J[k] for k in ("sb", "sb2", "lift", "sdec", "nthr", "nser") if k in J}
    return dict(Pj=Pj, order=order, degz=op_degz(Pj), seed=seed, alpha=alpha, label=str(label),
                seed_radius=(None if R is None else float(sp.Rational(str(R)))), knobs=knobs, path=path)


def run(tag, msq, dps, sb, sb2, lift, sdec_list, Nser, Nthr, frac, nloop):
    """The banana worked example (thin wrapper, behaviour unchanged): exact PF
    operator of the L-loop banana with squared masses msq via find_minimal_op,
    BFKNS multinomial series as the seed, then the generic kernel."""
    mp.mp.dps = dps
    print(f"\n[{tag}] Route A monodromy-projection, msq={msq}, dps={dps}")
    t0 = time.time()

    # ---- exact PF operator ----
    order, degz, Pj, resid, _ = find_minimal_op(list(msq), ord_max=10, degz_max=18)
    print(f"  PF order={order} degz={degz} residual={resid}")
    # ---- BFKNS seed (mp floats, as before: make_series_state passes them through) ----
    seed = wms_series_mp(list(msq), Nser)
    out = run_operator(tag, Pj, order, seed, dps=dps, sb=sb, sb2=sb2, lift=lift,
                       sdec_list=sdec_list, Nthr=Nthr, frac=frac, nloop=nloop, Nser=Nser,
                       seed_name="BFKNS seed", t0=t0, header=False)
    res = dict(tag=tag, msq=list(msq), dps=dps, order=order, degz=degz)
    res.update({k: v for k, v in out.items() if k not in res})
    return res


def run_operator(tag, Pj, order, seed_a, dps, sb, sb2, lift, sdec_list, Nthr, frac, nloop,
                 Nser=None, label=None, seed_name="seed series", t0=None, header=True):
    """GENERIC kernel.  Pj: {j: sympy polynomial in z} for L = sum_j P_j(z) theta_z^j
    (op_from_coeffs builds it from coefficient lists); seed_a: the exact Taylor
    coefficients a_n of the holomorphic solution f = sum a_n z^n at z = 0 (the point
    the anchors s_b = -1/z_b sit near); the fractional-exponent Frobenius branches
    Phi_alpha at z = infinity (local coordinate t = 1/z = -s) are read from the
    operator's own indicial equation, and c_alpha = the coefficient of f on Phi_alpha
    (branch arg t = -pi at Euclidean s_dec > 0) is returned for every fractional alpha.
    Anchors sb, sb2 must lie inside the seed's convergence region (s_b > 1/R): the
    caller checks that (coalescer.py does, by name, before calling)."""
    mp.mp.dps = dps
    if t0 is None:
        t0 = time.time()
    if Nser is None:
        Nser = len(seed_a) - 1
    if header:
        print(f"\n[{tag}] Route A monodromy-projection, {label or 'operator'}, dps={dps}")
        print(f"  operator order={order} degz={op_degz(Pj)} (given), seed terms={Nser}")

    # ---- indicial data at the threshold z = infinity ----
    ind_t, exps_t, _ = indicial_threshold(Pj)
    frac_exps = sorted([k for k in exps_t if not k.is_integer])
    int_exps = [k for k in exps_t for _ in range(exps_t[k]) if k.is_integer]
    nilp = len(int_exps)  # >= nilpotent index of unipotent block
    if not frac_exps:
        raise ValueError("no fractional exponent at z = infinity: nothing to project onto "
                         f"(indicial roots {dict(exps_t)})")
    print(f"  threshold exponents: {dict(exps_t)}  fractional={frac_exps}  nilp_rank={nilp}")
    eigs = {sp.Rational(a): sp.simplify(sp.exp(2*sp.pi*sp.I*a)) for a in frac_exps}
    print(f"  semisimple eigenvalues: {eigs}")

    # ---- companion ODE in s, singularities ----
    num, den, svar, sings = build_ode_s(Pj, order)
    sings_re = sorted(set(round(float(s_.real), 4) for s_ in sings if abs(s_.imag) < 1e-6))
    print(f"  s-singularities (real): {sings_re}")

    # ---- seed state at the deep anchors (inside the seed's disk), truncation self-check ----
    state_fn = make_series_state(seed_a, order, Nser)
    Yb = state_fn(mp.mpf(sb)); Yb2 = state_fn(mp.mpf(sb2))
    Nchk = Nser - 200 if Nser > 300 else max(1, (2 * Nser) // 3)
    Yb_chk = make_series_state(seed_a, order, Nchk)(mp.mpf(sb))
    seed_d = digits(Yb[0, 0], Yb_chk[0, 0])
    print(f"  {seed_name} Nser={Nser} convergence at s_b={sb}: {seed_d:.1f} d  [{time.time()-t0:.1f}s]")

    ordr = int(dps / (-mp.log10(frac))) + 12
    print(f"  adaptive Taylor: ordr={ordr}, frac={frac}")

    # ---- threshold Frobenius: clean power branches Phi_alpha for each fractional alpha ----
    terms = theta_t_terms(Pj, order)
    a_alpha = {}
    for alpha in frac_exps:
        a, obstr = frob_power_branch(terms, alpha, N=Nthr)
        assert obstr is None, f"obstruction at n={obstr} for rho={alpha}"
        a_alpha[alpha] = a
    print(f"  Frobenius branches built (Nthr={Nthr}): {[str(a) for a in frac_exps]}  [{time.time()-t0:.1f}s]")

    results = {}
    cvals = {str(alpha): {} for alpha in frac_exps}
    for s_dec in sdec_list:
        print(f"\n  >>> s_dec = {float(s_dec)}")
        wp = [mp.mpc(sb, 0), mp.mpc(sb, lift), mp.mpc(s_dec, lift), mp.mpc(s_dec, 0)]
        Yd = transport_adaptive(num, den, Yb, wp, ordr, sings, frac=frac)
        wp2 = [mp.mpc(sb2, 0), mp.mpc(sb2, lift), mp.mpc(s_dec, lift), mp.mpc(s_dec, 0)]
        Yd2 = transport_adaptive(num, den, Yb2, wp2, ordr, sings, frac=frac)
        gate = min(digits(Yd[k, 0], Yd2[k, 0]) for k in range(order))
        print(f"      two-anchor transport gate at s_dec: {gate:.1f} d  [{time.time()-t0:.1f}s]")

        # ---- monodromy matrix M_0: transport identity CCW around |s|=s_dec ----
        loop = loop_waypoints(s_dec, npts=nloop, ccw=True)
        I_ = mp.eye(order)
        M = transport_adaptive(num, den, I_, loop, ordr, sings, frac=frac, ncols=order)
        try:
            E, _ = mp.eig(M)
            eigdist = {str(a): float(min(abs(e - mp.mpc(complex(eigs[a]))) for e in E))
                       for a in frac_exps}
        except Exception:
            eigdist = None
        detM = mp.det(M)
        # nilpotent-index probe: ||(M-I)^k|| Frobenius norm
        MmI = M - I_
        nilprobe = []
        Mk = mp.eye(order)
        for k in range(1, nilp + 2):
            Mk = Mk * MmI
            nilprobe.append(float(mp.mnorm(Mk, 'F')))
        print(f"      det(M)={mp.nstr(detM,20)}  eig-dist={eigdist}")
        print(f"      ||(M-I)^k||_F, k=1..{nilp+1}: {['%.3g'%x for x in nilprobe]}")

        # ---- spectral projectors P_lambda for each fractional alpha ----
        eigs_mp = {a: mp.mpc(complex(eigs[a])) for a in frac_exps}
        sdec_out = {}
        for alpha in frac_exps:
            lam = eigs_mp[alpha]
            # P = (M-I)^nilp * prod_{mu!=lam} (M-mu I) / [(lam-1)^nilp * prod (lam-mu)]
            P = matpow(MmI, nilp)
            denom = (lam - 1) ** nilp
            for beta in frac_exps:
                if beta == alpha:
                    continue
                mu = eigs_mp[beta]
                P = P * (M - mu * I_)
                denom *= (lam - mu)
            P = P / denom
            PY = P * Yd; PY2 = P * Yd2
            phi = phi_theta_state(a_alpha[alpha], alpha, s_dec, order, branch='lower', N=Nthr)
            c_rows = [PY[k, 0] / phi[k, 0] for k in range(order)]
            c_rows2 = [PY2[k, 0] / phi[k, 0] for k in range(order)]
            row_cons = min(digits(c_rows[0], c_rows[k]) for k in range(1, order))
            anc_cons = min(digits(c_rows[k], c_rows2[k]) for k in range(order))
            c = c_rows[0]
            print(f"      c_{{{alpha}}} = {mp.nstr(c, min(50, dps-10))}")
            print(f"        row-consistency ({order} components agree): {row_cons:.1f} d   "
                  f"two-anchor: {anc_cons:.1f} d")
            cvals[str(alpha)][str(float(s_dec))] = c
            sdec_out[str(alpha)] = dict(
                c=mp.nstr(c, dps - 5), c_re=mp.nstr(mp.re(c), dps - 5), c_im=mp.nstr(mp.im(c), dps - 5),
                row_consistency_d=round(row_cons, 1), two_anchor_d=round(anc_cons, 1),
            )
        results[str(float(s_dec))] = dict(
            transport_gate_d=round(gate, 1), det_M=mp.nstr(detM, 30),
            eig_dist=eigdist, nilprobe=nilprobe, c_alpha=sdec_out,
        )

    # ---- s_dec-independence ----
    sdec_ind = {}
    for alpha in frac_exps:
        keys = list(cvals[str(alpha)].keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                sdec_ind[f"{alpha}:{keys[i]}_vs_{keys[j]}"] = round(
                    digits(cvals[str(alpha)][keys[i]], cvals[str(alpha)][keys[j]]), 1)
    print(f"\n  s_dec-independence: {sdec_ind}  [{time.time()-t0:.1f}s total]")

    out = dict(
        tag=tag, label=label, dps=dps, order=order, degz=op_degz(Pj),
        ind_threshold=str(ind_t), frac_exps=[str(a) for a in frac_exps], nilp=nilp,
        eigs={str(a): str(eigs[a]) for a in frac_exps},
        sb=sb, sb2=sb2, lift=lift, Nser=Nser, Nthr=Nthr, ordr=ordr, frac=frac,
        seed_convergence_d=round(seed_d, 1),
        per_sdec=results,
        sdec_independence=sdec_ind,
        elapsed_s=round(time.time()-t0, 1),
        c_alpha_best={str(a): mp.nstr(cvals[str(a)][list(cvals[str(a)].keys())[0]], dps - 5)
                      for a in frac_exps},
    )
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", choices=["K3", "CY3", "both"], default="both")
    ap.add_argument("--dps", type=int, default=60)
    ap.add_argument("--out", default="ROUTEA.json")
    ap.add_argument("--probe", action="store_true", help="small/fast settings for ETA")
    args = ap.parse_args()

    cfgs = []
    if args.target in ("K3", "both"):
        cfgs.append(("K3_control", (1,1,1,9), dict(sb=80, sb2=150, lift=10,
                     sdec_list=[mp.mpf('0.15'), mp.mpf('0.25')], Nser=1200, Nthr=120)))
    if args.target in ("CY3", "both"):
        cfgs.append(("CY3_target", (1,1,1,1,16), dict(sb=100, sb2=200, lift=20,
                     sdec_list=[mp.mpf('0.15'), mp.mpf('0.25'), mp.mpf('0.35')], Nser=1400, Nthr=150)))

    if args.probe:
        for c in cfgs:
            c[2].update(sdec_list=c[2]['sdec_list'][:1], Nser=500, Nthr=40)

    OUT = {}
    for tag, msq, kw in cfgs:
        OUT[tag] = run(tag, msq, args.dps, frac=0.3, nloop=32, **kw)
    json.dump(OUT, open(os.path.join(HERE, args.out), "w"), indent=1)
    print(f"\nWROTE {args.out}")

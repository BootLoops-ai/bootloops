#!/usr/bin/env python3
"""ffp.py (periods wing) — finite-field fingerprint engines for the
all-z Frobenius scan: HV convolution oracle + BCM Gauss-sum engine +
candidate checker / pass-set I/O + the validation batteries as importable
gates. Engine code lives in pipeline/ffp/; this facade only
routes, guards and packages — no fingerprint logic is duplicated.

FRONT DOOR — int64 prime law (known pitfall: mod-p block-16 int64 overflow):
every engine here accumulates in numpy int64 (HV convolution steps,
mu-row bincounts) or rounds complex Gauss sums against integer targets;
the validated envelope is p < 2^31 AND the production battery range
(HV p <= 127, BCM p <= 337). Primes >= 2^31 are REFUSED
(BigPrimeRefusal) unless the caller passes allow_big_prime=True and
owns the overflow analysis. The refusal is the tool remembering the
31-bit block-16 overflow: garbage modulus mimics huge-degree functions.

Engines (delegation only):
  hv_torus_counts / hv_twisted12_counts   AESZ34/HV torus fiber counts
     over F_q for ALL psi in one pass (pipeline/ffp/hv_count.py)
  hv_passset_AESZ34   production pass-set at one prime (pilot2 law:
     E-law blind a, k-class, F_p^2 side, class Weil-windows, quartic
     split; the in-law asserts are the integrity gate — do NOT weaken)
  bcm_gauss_field / bcm_h_coeffs / bcm_h_all_fp / bcm_passset   BCM
     finite-hypergeometric a_p(z) engine, complex-exact Gauss sums with
     float residual certified < 0.05 before rounding (scale/bcm*.py);
     unit-root congruence asserted at EVERY z in production pass-sets
  load_passsets / check / parse_spec / uvd_minpoly   candidate checker:
     per-prime IN/OUT + AND verdict, split primes only for quadratic z,
     skips are named (inert/ramified/bad-den/singular), survivor law
     >= MIN_CHECKABLE checkable primes (scale/check_candidate.py)
Validation gates (importable; the batteries that earned the oracle):
  gate_brute_small_q          brute-force count equality at tiny q
  gate_unit_root_congruence   trace4(z) == truncated series mod p, all z
  gate_weil_integrality       b-integrality + exact Weil circle + split
  gate_groundtruth_p19        CDEvS 1912.06146 p=19 table replay
  checker_gate_battery / bcm_gate_battery   full regression batteries
     (T1/T2 ALL-IN + 9 controls AND-fail; K4/K11 12/12 + negatives) —
     pure replay, no file writes
Battery: selftest_ffp.py; reference receipts in pipeline/ffp/banks/."""
import glob, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
FFP_DIR = os.path.join(HERE, "pipeline", "ffp")
SCALE_DIR = os.path.join(FFP_DIR, "scale")
PASSSET_DIR = os.path.join(SCALE_DIR, "passsets")
BANK_DIR = os.path.join(FFP_DIR, "banks")
for _d in (FFP_DIR, SCALE_DIR):
    if _d not in sys.path:
        sys.path.insert(0, _d)

INT64_PRIME_CEILING = 1 << 31
MIN_CHECKABLE = 5


class BigPrimeRefusal(RuntimeError):
    """Raised by the front door for p >= 2^31 without allow_big_prime."""


def require_int64_prime(p, allow_big_prime=False):
    """Front door: int64 law. Returns p; refuses p >= 2^31 unless the
    caller explicitly owns the overflow analysis (allow_big_prime)."""
    if p >= INT64_PRIME_CEILING and not allow_big_prime:
        raise BigPrimeRefusal(
            f"REFUSED: p = {p} >= 2^31. The HV/BCM engines accumulate in "
            "int64 and are validated only for p < 2^31 (production: HV "
            "p<=127, BCM p<=337). A 31-bit prime overflowing an int64 "
            "convolution yields a garbage modulus that mimics huge-degree "
            "functions (a measured pitfall, mod-p block-16). Pass "
            "allow_big_prime=True ONLY with your own overflow analysis.")
    return p


# ------------------------------------------------ HV convolution oracle
def hv_torus_counts(p, deg, nvar=5, allow_big_prime=False):
    """N_torus(psi) for ALL psi in F_q, q = p^deg (hv_count law)."""
    require_int64_prime(p, allow_big_prime)
    from hv_count import torus_counts
    return torus_counts(p, deg, nvar)


def hv_twisted12_counts(p, nvar_fixed=3, allow_big_prime=False):
    """(12)-twisted count over F_p (one conjugate pair + fixed vars)."""
    require_int64_prime(p, allow_big_prime)
    from hv_count import twisted12_counts
    return twisted12_counts(p, nvar_fixed)


def hv_passset_AESZ34(p, allow_big_prime=False):
    """Production AESZ34 pass-set at prime p (pilot2 E-law + k-class +
    N2 windows + quartic split; the in-law asserts ARE the integrity
    gate). Returns the pass-set dict."""
    require_int64_prime(p, allow_big_prime)
    from pilot2 import compute, test_point
    t0 = time.time()
    try:
        ctx = compute(p)
    except AssertionError as e:
        return dict(op="AESZ34", p=p, status="BAD-PRIME-LAW-ASSERT",
                    detail=repr(e))
    fires, alphas = [], {}
    for f in ctx["smooth"]:
        r = test_point(ctx, f)
        if r["fired"]:
            fires.append(int(f))
            alphas[str(f)] = sorted({h["alpha"] for h in r["hits"]})
    n = len(ctx["smooth"])
    return dict(op="AESZ34", p=p, status="OK",
                oracle="HV-torus pilot2 v1 (E-law+windows+quartic)",
                n_smooth=n, sing=sorted(int(x) for x in ctx["sing"]),
                fire_count=len(fires), rate=round(len(fires) / n, 4),
                fires=fires, alphas=alphas,
                windows={str(k): list(v) for k, v in ctx["wins"].items()},
                includes_fp2_pass=True,
                compute_seconds=round(time.time() - t0, 1))


# ------------------------------------------------ BCM Gauss-sum engine
def bcm_gauss_field(p, deg, allow_big_prime=False):
    """GaussField(p, deg) with the front door applied (complex-exact
    Gauss-sum table; g(triv) == -1 and |g| == sqrt(q) asserted)."""
    require_int64_prime(p, allow_big_prime)
    from bcm import GaussField
    return GaussField(p, deg)


def bcm_h_coeffs(gf, alpha):
    from bcm import h_coeffs
    return h_coeffs(gf, alpha)


def bcm_h_all_fp(gf, C, conj=False):
    from bcm import h_all_fp
    return h_all_fp(gf, C, conj)


def bcm_passset(op, hname, ctrl, p, alpha=None, M=None,
                allow_big_prime=False):
    """One BCM production pass-set (bcm_passsets.run VERBATIM: FIRE :=
    b integral & weil_circle & split_12; unit-root congruence asserted
    at every z; float residual certified < 0.05). ctrl = (num, den) of
    the positive-control point. hname keys bcm_validate.OPS; pass
    (alpha, M) to register a NEW hypergeometric datum under hname
    ((q-1)*a_i integrality asserted downstream). No file writes."""
    require_int64_prime(p, allow_big_prime)
    from bcm_validate import OPS
    if alpha is not None:
        OPS[hname] = (list(alpha), int(M))
    import bcm_passsets
    return bcm_passsets.run(op, hname, ctrl, p)


# --------------------------------------- candidate checker / pass-set I/O
def load_passsets(op, dirpath=None):
    """{p: {fires, sing, alphas}} from <dirpath>/<op>.p*.json (status OK
    only). Default dirpath = packaged pipeline/ffp/scale/passsets."""
    dirpath = dirpath or PASSSET_DIR
    ps = {}
    for fn in sorted(glob.glob(os.path.join(dirpath, f"{op}.p*.json"))):
        d = json.load(open(fn))
        if d.get("status") == "OK":
            ps[d["p"]] = dict(fires=set(d["fires"]), sing=set(d["sing"]),
                              alphas=d.get("alphas", {}))
    return ps


def check(cand, passsets):
    """Per-prime IN/OUT + AND verdict (check_candidate.check law:
    quadratic z checked at split primes only, both Galois reductions
    must fire; named skips; survivor needs >= MIN_CHECKABLE=5)."""
    import check_candidate as _cc
    return _cc.check(cand, passsets)


def parse_spec(kind, z=None, minpoly=None):
    import check_candidate as _cc
    return _cc.parse_spec(kind, z=z, minpoly=minpoly)


def uvd_minpoly(u, v, d):
    import check_candidate as _cc
    return _cc.uvd_minpoly(u, v, d)


def residues_mod_p(cand, p):
    import check_candidate as _cc
    return _cc.residues_mod_p(cand, p)


# ------------------------------------ validation batteries (importable)
def gate_brute_small_q(torus_cases=((7, 1), (11, 1), (13, 1), (5, 2)),
                       twist_ps=(5, 7, 11)):
    """Brute-force count equality at tiny q (brute.py battery). Exact
    numpy array equality per case; receipt dict with 'pass'."""
    import numpy as np
    from hv_count import torus_counts, twisted12_counts
    from brute import brute_torus, brute_twisted12
    out = {}
    for (p, deg) in torus_cases:
        _K, F = torus_counts(p, deg)
        out[f"torus_p{p}_deg{deg}"] = bool(np.array_equal(
            F, brute_torus(p, deg)))
    for p in twist_ps:
        out[f"twist12_p{p}"] = bool(np.array_equal(
            twisted12_counts(p), brute_twisted12(p)))
    out["pass"] = all(out.values())
    return out


def gate_unit_root_congruence(hname, p, allow_big_prime=False):
    """d=4 unit-root congruence: -a(z) == sum_{n<p} A_n z^n mod p at ALL
    smooth z, a = -H_p(Mz) (bcm_validate V2 production convention)."""
    require_int64_prime(p, allow_big_prime)
    import numpy as np
    from bcm import GaussField, h_coeffs, h_all_fp
    from bcm_validate import OPS, series_trunc_mod_p, horner_all
    alpha, M = OPS[hname]
    gf = GaussField(p, 1)
    H1r = h_all_fp(gf, h_coeffs(gf, alpha))
    resid = float(np.max(np.abs(H1r[1:] - np.round(H1r.real[1:]))))
    H1 = np.round(H1r.real).astype(np.int64)
    S = horner_all(series_trunc_mod_p(alpha, M, p), p)
    bad, nchk = [], 0
    for z in range(1, p):
        t = M * z % p
        if t in (0, 1):
            continue
        nchk += 1
        if (int(H1[t]) - S[z]) % p:          # -a == S[z], a = -H1[t]
            bad.append(z)
    return dict(gate="unit-root-congruence", hname=hname, p=p,
                float_resid=resid, n_checked=nchk, bad_z=bad,
                **{"pass": (not bad) and resid < 0.05})


def gate_weil_integrality(a, s2, p):
    """Weight-3 quartic acceptance chain on (a, s2): b = (a^2-s2)/(2p)
    integral, exact Weil circle, rank-2 split_12. Receipt dict."""
    from weil import weil_circle
    from quartic import split_12
    if (a * a - s2) % (2 * p):
        return dict(gate="weil-b-integrality", p=p, a=a, s2=s2,
                    b_integral=False, weil=False, splits=[], fire=False)
    b = (a * a - s2) // (2 * p)
    w = bool(weil_circle(a, b, p))
    sp = split_12(a, b, p) if w else []
    return dict(gate="weil-b-integrality", p=p, a=a, s2=s2, b=b,
                b_integral=True, weil=w, splits=sp, fire=bool(sp))


def gate_groundtruth_p19():
    """Replay CDEvS 1912.06146 p=19 AESZ34 table: all smooth rows pass
    weil_circle; b+40 perturbation must fail most rows; every 'fac' row's
    (alpha, beta) is recovered by split_12."""
    import groundtruth_p19 as gt
    from weil import weil_circle
    from quartic import split_12
    smooth = [(f, gt.ab(f)) for f in sorted(gt.ROWS) if gt.ab(f)]
    n = len(smooth)
    n_weil = sum(bool(weil_circle(a, b, gt.P)) for _f, (a, b) in smooth)
    n_pert = sum(bool(weil_circle(a, b + 40, gt.P))
                 for _f, (a, b) in smooth)
    fac_bad = [f for f in sorted(gt.ROWS) if gt.ROWS[f][0] == "fac"
               and (gt.ROWS[f][1], gt.ROWS[f][2])
               not in split_12(*gt.ab(f), gt.P)]
    return dict(gate="groundtruth-p19", n_smooth=n, n_weil_pass=n_weil,
                n_perturbed_pass=n_pert, fac_rows_missed=fac_bad,
                **{"pass": n_weil == n and n_pert < n and not fac_bad})


def checker_gate_battery(dirpath=None):
    """AESZ34 checker regression (checker_gate.py CASES, replay only —
    NO file writes): T1/T2 ALL-IN at >= 5 checkable primes, 9 pilot
    controls fail the AND. Returns dict(verdict, rows)."""
    from checker_gate import CASES
    ps = load_passsets("AESZ34", dirpath)
    rows, ok = [], True
    for name, cand, expect in CASES:
        r = check(cand, ps)
        good = (r["and_pass"] and r["n_checkable"] >= MIN_CHECKABLE) \
            if expect == "ALL-IN" else (not r["and_pass"])
        ok &= good
        rows.append(dict(name=name, expected=expect,
                         gate="PASS" if good else "FAIL",
                         n_in=r["n_in"], n_checkable=r["n_checkable"]))
    return dict(gate="CHECKER-REGRESSION", n_primes=len(ps),
                verdict="PASS" if ok else "FAIL", rows=rows)


def bcm_gate_battery(dirpath=None):
    """BCM positive/negative control battery (bcm_gate.py law, replay
    only — NO file writes): AESZ4/AESZ11 control points fire at ALL
    pass-set primes with alpha matched to the identified wt-2 curve
    (and AESZ11 beta to 180.4.a.e where archived); pilot rational
    controls must fail the AND. Returns dict(verdict, rows)."""
    from fractions import Fraction
    from bcm_gate import ap_curve, CURVES, NF180, NEG
    dirpath = dirpath or PASSSET_DIR
    rows, ok = [], True
    for op, ctrl in (("AESZ4", Fraction(-1, 5832)),
                     ("AESZ11", Fraction(-1, 432))):
        ps = load_passsets(op, dirpath)
        r = check(dict(kind="rational", z=ctrl), ps)
        amatch = bmatch = True
        for p in sorted(ps):
            d = json.load(open(os.path.join(dirpath, f"{op}.p{p}.json")))
            sp = d["control"].get("splits", [])
            amatch &= any(al == ap_curve(*CURVES[op], p) for al, _be in sp)
            if op == "AESZ11" and p in NF180:
                bmatch &= any(be == NF180[p] for _al, be in sp)
        g = r["and_pass"] and r["n_checkable"] == len(ps) and amatch \
            and bmatch
        ok &= g
        rows.append(dict(gate=f"{op}-POSITIVE-CONTROL", z=str(ctrl),
                         verdict="PASS" if g else "FAIL",
                         alpha_match=amatch, beta_match=bmatch,
                         n_in=r["n_in"], n_checkable=r["n_checkable"]))
        for z in NEG:
            rn = check(dict(kind="rational", z=Fraction(z)), ps)
            ok &= not rn["and_pass"]
            rows.append(dict(gate=f"{op}-NEG-{z}", verdict="PASS"
                             if not rn["and_pass"] else "FAIL"))
    return dict(gate="BCM-GATE", verdict="PASS" if ok else "FAIL",
                rows=rows)

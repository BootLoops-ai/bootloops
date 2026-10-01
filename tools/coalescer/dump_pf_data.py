#!/usr/bin/env python3
"""Dump exact PF operator (theta + d/dz forms), Frobenius branch coeffs, and
BFKNS series for the independent Julia/Arb route.  Generic in (msq, order, degz).

Writes two files: <out>.jl (Route B's data) and <out>.json.  The JSON is also the
operator + seed layout that `coalescer.py --op FILE.json` reads (GUIDE.md INPUTS):
theta_form_Pj = {j: P_j(z) integer coefficients, low -> high power of z} for
L = sum_j P_j(z) theta_z^j, order, bfkns_a = the exact seed coefficients a_n of
the holomorphic solution at z = 0, frac_exps, label.

    python3 dump_pf_data.py --msq 1,1,1,9 --nthr 250 --nser 300 --out pf_k3_data.jl
(a bare --out name lands beside this file, where routeB_monoproj.jl looks; a path
with a directory is used as given).  Importable: dump(msq, nthr, nser, out)."""
import sys, os, json, argparse
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pf_engine import (find_minimal_op, theta_t_terms, frob_power_branch,
                       wms_series_int, indicial_threshold, indicial_mum, z)

HERE = os.path.dirname(os.path.abspath(__file__))


def dump(msq, nthr=250, nser=300, out="pf_cy3_data.jl", verbose=True):
    """Build the banana PF operator for squared masses msq and write <out> (.jl) and
    the matching .json; returns (jl_path, json_path, summary dict)."""
    msq = [int(x) for x in msq]
    log = print if verbose else (lambda *a, **k: None)
    jl_path = out if (os.path.isabs(out) or os.path.dirname(out)) else os.path.join(HERE, out)
    json_path = os.path.splitext(jl_path)[0] + ".json"

    order, degz, Pj, resid, _ = find_minimal_op(msq, ord_max=10, degz_max=18)
    ind_t, exps_t, _ = indicial_threshold(Pj)
    ind_m, exps_m, _ = indicial_mum(Pj)
    frac_exps = sorted([k for k in exps_t if not k.is_integer])
    log(f"msq={msq} order={order} degz={degz} residual={resid}")
    log(f"threshold exps={dict(exps_t)}  fractional={frac_exps}")

    # theta_z-form: P_j(z) coeffs
    theta_coeffs = {}
    for j in range(order + 1):
        cs = sp.Poly(Pj.get(j, 0), z).all_coeffs()[::-1]
        theta_coeffs[j] = [int(c) for c in cs]

    # d/dz-form via Stirling-2nd: theta^j = sum_k S2(j,k) z^k D^k
    S2 = [[0] * (order + 1) for _ in range(order + 1)]
    S2[0][0] = 1
    for n in range(1, order + 1):
        for k in range(1, n + 1):
            S2[n][k] = k * S2[n - 1][k] + S2[n - 1][k - 1]
    ddz = {k: sp.Integer(0) for k in range(order + 1)}
    for j in range(order + 1):
        for k in range(j + 1):
            ddz[k] += Pj.get(j, 0) * S2[j][k] * z ** k
    ddz_coeffs = {}
    for k in range(order + 1):
        cs = sp.Poly(sp.expand(ddz[k]), z).all_coeffs()[::-1]
        ddz_coeffs[k] = [int(c) for c in cs]

    # exact-Q Frobenius power branches at threshold
    terms = theta_t_terms(Pj, order)
    abranch = {}
    for alpha in frac_exps:
        a, obstr = frob_power_branch(terms, alpha, N=nthr)
        assert obstr is None, f"obstruction at n={obstr} for rho={alpha}"
        abranch[str(alpha)] = [[int(sp.Rational(a[n]).p), int(sp.Rational(a[n]).q)] for n in sorted(a)]
        log(f"  Phi_{alpha}: a_0..a_3 = {[a[n] for n in range(4)]}")

    # BFKNS series (exact ints)
    a_int = wms_series_int(msq, nser)

    # emit Julia data file
    with open(jl_path, "w") as f:
        f.write(f"# auto-generated: msq={msq} order={order} degz={degz}\n")
        f.write(f"const ORDER = {order}\n")
        f.write(f"const MSQ = {msq}\n")
        f.write(f"const SUM_MSQ = {sum(msq)}\n")
        f.write("const DDZ_COEFFS = [\n")
        for k in range(order + 1):
            f.write(f"  big.({ddz_coeffs[k]}),\n")
        f.write("]\n")
        f.write("const THETA_COEFFS = [\n")
        for j in range(order + 1):
            f.write(f"  big.({theta_coeffs[j]}),\n")
        f.write("]\n")
        f.write(f"const FRAC_EXPS = [{', '.join(f'{sp.Rational(a).p}//{sp.Rational(a).q}' for a in frac_exps)}]\n")
        for i, alpha in enumerate(frac_exps):
            nums = [pq[0] for pq in abranch[str(alpha)]]
            dens = [pq[1] for pq in abranch[str(alpha)]]
            f.write(f"const A{i}_NUM = big.({nums})\n")
            f.write(f"const A{i}_DEN = big.({dens})\n")
        f.write(f"const BFKNS_A = big.({[int(x) for x in a_int]})\n")

    # also JSON (the --op layout: theta_form_Pj + order + bfkns_a (+ label, frac_exps, seed_radius))
    thr = sp.nsimplify(sum(sp.sqrt(m) for m in msq) ** 2)
    summary = dict(
        label=f"banana msq={tuple(msq)}",
        msq=msq, order=order, degz=degz, residual=str(resid),
        ind_threshold=str(ind_t), exps_threshold={str(k): int(v) for k, v in exps_t.items()},
        ind_mum=str(ind_m), exps_mum={str(k): int(v) for k, v in exps_m.items()},
        frac_exps=[str(a) for a in frac_exps],
        seed_radius=str(sp.Rational(1) / thr) if thr.is_rational else str(sp.N(1 / thr, 17)),
        theta_form_Pj=theta_coeffs, ddz_form_pk=ddz_coeffs,
        leading_factored=str(sp.factor(Pj[order])),
        abranch=abranch, bfkns_a=[int(x) for x in a_int],
    )
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=1)
    log(f"wrote {jl_path} and {json_path}")
    return jl_path, json_path, summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--msq", default="1,1,1,1,16")
    ap.add_argument("--nthr", type=int, default=250)
    ap.add_argument("--nser", type=int, default=300)
    ap.add_argument("--out", default="pf_cy3_data.jl")
    args = ap.parse_args()
    dump([int(x) for x in args.msq.split(",")], nthr=args.nthr, nser=args.nser, out=args.out)

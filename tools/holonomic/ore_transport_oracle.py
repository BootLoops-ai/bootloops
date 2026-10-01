#!/usr/bin/env python3
r"""ore_transport_oracle.py -- a code-independent SECOND IMPLEMENTATION of the house period transports:
SageMath ore_algebra `numerical_transition_matrix` (rigorous Arb balls) on OUR operators along OUR paths,
emitting a receipt with the ball radii and the digits of agreement against the house values.

extends: tools/holonomic (the RELOCATED-SAGEENV pattern: the pinned SageMath 10.7 + ore_algebra 0.5 venv
         invoked by its bin/python3.11 DIRECTLY, sage.all imported first, algorithm='naive' mandatory,
         the IC convention DETERMINED per session by holonomic_transport.determine_ic_convention) and
         tools/amflow-kit/compare_amflow_json.py (the pair form: every value a ball "[mid +/- rad]",
         matched_digits on the midpoints, the radius kept beside, OVERALL PASS/FAIL at --min-digits).
         INDEPENDENCE RULE: this wrapper is the external oracle the house transports are gated AGAINST;
         it imports nothing from detransport / famhar / baller / the served evaluators.

RUN CONTEXT (mandatory; the manual's law):
  SE=<the pinned sageenv root (SAGEENV_PIN.md; HOLONOMIC_SAGEENV)>
  nice -n 10 ionice -c3 prlimit --as=<bytes> $SE/bin/python3.11 ore_transport_oracle.py --spec SPEC.json --out RECEIPT.json
  (HOLONOMIC_DIR = the dir holding holonomic_transport.py, default: this file's own dir or ../holonomic;
   COMPARE_AMFLOW_JSON = the amflow-kit comparator, default: ../amflow-kit/compare_amflow_json.py.)

THE SPEC (JSON):
  name            : a label
  operator        : {"form": "theta_z_Pj", "chart": "z=-1/s", "order": r, "degz": d, "Pj": {"j": [p_j0, p_j1, ...]},
                     "source": {"path","sha256"}}      -- the theta_z-form of record (P_j(z) low -> high),
                    the operator L = sum_j P_j(z) theta_z^j rewritten in the s-chart: theta_z = -s d/ds,
                    s^d P_j(-1/s) polynomial;  or
                    {"form": "ddz_coefD_s", "coefD": [[p_k(s) low->high] for k], ...}   -- a d/ds form in s;  or
                    {"form": "theta_z_Pj", "chart": "z", ...}  -- the theta form in the variable z itself
                    (theta_z = z d/dz), optionally "compose_left_D": true (D o L, the homogeneous operator of an
                    inhomogeneous L y = const).
  variable        : the name of the transport variable ("s" or "z")
  seed            : {"kind": "bfkns_multinomial", "msq": [...], "at": "100", "N": 1200}  -- the exact
                    multinomial-squared period series sum_n c_n x^{-n} (x = s) summed EXACTLY over Q with a
                    RIGOROUS tail bound |c_n| <= thr^n, thr = (sum sqrt M_i)^2, added as the ball radius;  or
                    {"kind": "explicit", "at": "1/128", "values": ["..."], "radii": ["..."]}  -- the jet
                    [f, f', ..., f^(r-1)] at the base point as decimal strings (+ radii).
  path            : ["100", "50+50*I", "5+12*I", "1+3*I", "1", "0"]  -- exact points (rationals / Gaussian
                    rationals), the base point first; the last point may be a regular singular point of L
                    (a Frobenius landing: the coordinates come back in ore_algebra's canonical local basis).
  landing         : {"kind": "local_basis", "fractional": {"5/4": {"phase": "exp(i*pi*5/4)"}, ...}}  -- pick the
                    local-basis coordinates whose leading monomial is x^alpha and multiply by the phase (the
                    t = e^{-i pi} s convention of record: s^alpha = e^{i pi alpha} t^alpha);  or
                    {"kind": "regular_point"}  -- the transported jet itself is the result.
  targets         : {"5/4": {"label": ..., "value": "<decimal re>", "im": "<decimal im>|0", "digits": N,
                             "source": {...}}, ...}  -- the house values (by object) for the agreement;
                    for a regular landing the keys are "0", "1", ... (the jet components).
  eps, prec       : the engine tolerance (e.g. "1e-80") and the ball precision in bits (default 430).
  min_digits      : the pair-form bar (default 30).

WHAT THE RECEIPT CARRIES: engine versions (must be the pin 10.7 / 0.5), the IC convention determined this
session, the operator's sha + order + leading coefficient's real roots vs the path, the seed's exact tail bound,
the path, eps, prec, wall and maxrss, every transported entry as "[mid +/- rad]" (re, im), per target:
our midpoint, the ball radius as digits (-log10(rad/|mid|)), matched_digits vs the house value on the
midpoints, the certified agreement = min(matched, ball digits), PASS/FAIL at min_digits; OVERALL.
--smoke-toy runs two engine-free-of-our-data toys in seconds and exits (cos(1) from Ds^2 + 1 on [0 -> 1];
the fractional landing x^{1/2} from x Dx - 1/2 on [1 -> 0], coordinate 1 on the local basis).
EXIT: 0 OVERALL PASS; 1 a compared component below the bar (named); 2 usage / spec error; 3 the engine
versions are not the pin (refused before any transport); 4 a source sha pin mismatch (refused by name).
"""
import argparse
import hashlib
import importlib.util
import json
import os
import re
import resource
import subprocess
import sys
import time
from fractions import Fraction as Fr

EXIT_FAIL, EXIT_USAGE, EXIT_ENGINE, EXIT_PIN = 1, 2, 3, 4
SAGE_PIN, ORE_PIN = "10.7", "0.5"
_HERE = os.path.dirname(os.path.abspath(__file__))
# defaults: this file is meant to live in tools/holonomic/ beside holonomic_transport.py; amflow-kit is its sibling package
HOLONOMIC_DIR_DEFAULT = os.environ.get("HOLONOMIC_DIR", _HERE if os.path.exists(os.path.join(_HERE, "holonomic_transport.py")) else os.path.join(_HERE, "..", "holonomic"))
COMPARE_DEFAULT = os.environ.get("COMPARE_AMFLOW_JSON", os.path.join(os.path.dirname(HOLONOMIC_DIR_DEFAULT), "amflow-kit", "compare_amflow_json.py"))


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def load_by_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- exact period series
def multinomial_squared_series(msq, N):
    """c_n = sum_{|k|=n} (n!/prod k_i!)^2 prod M_i^{k_i} as exact ints, n = 0..N (via the product of
    sum_k M^k x^k / k!^2 and n!^2)."""
    from math import factorial
    prod = [Fr(1)]
    for M in msq:
        ser = [Fr(M) ** k / Fr(factorial(k)) ** 2 for k in range(N + 1)]
        new = [Fr(0)] * (N + 1)
        for i, a in enumerate(prod):
            if a == 0:
                continue
            for j in range(N + 1 - i):
                new[i + j] += a * ser[j]
        prod = new
    out = []
    for n in range(N + 1):
        c = prod[n] * Fr(factorial(n)) ** 2
        assert c.denominator == 1
        out.append(int(c))
    return out


def seed_jet_exact(msq, sb, N, order):
    """[f^(k)(sb)]_{k<order} of f(s) = sum_n c_n s^{-n} as exact rationals, with the rigorous tail bound
    T_k <= sb^{-k} (N+1+k)^k (thr/sb)^{N+1} / (1 - rho_k), rho_k = e^{k/(N+1)} thr/sb < 1 (asserted),
    from |c_n| <= thr^n, thr = (sum sqrt M_i)^2 (Cauchy-Schwarz on the multinomial sum)."""
    import math
    c = multinomial_squared_series(msq, N)
    m = [math.isqrt(M) for M in msq]
    assert all(mi * mi == M for mi, M in zip(m, msq))
    thr = Fr(sum(m) ** 2)
    sb = Fr(sb)
    assert sb > thr, f"seed point {sb} inside the threshold disc (thr = {thr})"
    jets, tails = [], []
    for k in range(order):
        tot = Fr(0)
        for n in range(N + 1):
            ff = 1
            for i in range(k):
                ff *= (-n - i)
            tot += Fr(c[n] * ff) / sb ** (n + k)
        rho = math.exp(k / (N + 1)) * float(thr / sb)
        assert rho < 1
        T = float(sb) ** (-k) * (N + 1 + k) ** k * float(thr / sb) ** (N + 1) / (1 - rho)
        jets.append(tot)
        tails.append(T)
    return jets, tails, thr, c[:6]


# --------------------------------------------------------------------------- operator construction
def build_operator(opspec, var, Pol, D):
    """Return the OreAlgebra operator per opspec (see the module docstring)."""
    x = Pol.gen()
    form = opspec["form"]
    if form == "ddz_coefD_s":
        L = 0
        for k, co in enumerate(opspec["coefD"]):
            L += Pol([Fr(str(cf)) for cf in co]) * D ** k
        return L
    if form == "theta_z_Pj":
        order = int(opspec["order"])
        Pj = {int(j): [Fr(str(cf)) for cf in co] for j, co in opspec["Pj"].items()}
        chart = opspec.get("chart", "z=-1/s")
        if chart == "z=-1/s":
            degz = int(opspec.get("degz", max(len(co) - 1 for co in Pj.values())))
            th = -x * D                         # theta_z = -s d/ds
            L = 0
            for j, co in Pj.items():
                # s^degz P_j(-1/s) = sum_k p_jk (-1)^k s^{degz-k}
                Q = Pol({degz - k: cf * (-1) ** k for k, cf in enumerate(co)})
                L += Q * th ** j
        elif chart == "z":
            th = x * D                          # theta_z = z d/dz
            L = 0
            for j, co in Pj.items():
                L += Pol(co) * th ** j
        else:
            raise ValueError(f"unknown chart {chart!r}")
        if opspec.get("compose_left_D"):
            L = D * L
        return L
    raise ValueError(f"unknown operator form {form!r}")


def parse_point(s, QQ, QQbar, I):
    """'100' | '1/4' | '50+50*I' | '-3+3*I' -> an exact Sage number."""
    s = str(s).replace(" ", "")
    if "I" not in s:
        return QQ(Fr(s))
    m = re.fullmatch(r"([+-]?[0-9/.]+)?([+-][0-9/.]*)\*?I", s)
    if not m:
        raise ValueError(f"cannot parse path point {s!r}")
    re_part = Fr(m.group(1)) if m.group(1) else Fr(0)
    im_txt = m.group(2)
    im_part = Fr(im_txt) if im_txt not in ("+", "-") else Fr(im_txt + "1")
    return QQbar(QQ(re_part)) + QQbar(QQ(im_part)) * I


def flip_digit(s, k):
    out = list(s)
    seen = 0
    for i, ch in enumerate(out):
        if ch.isdigit() and (seen or ch != "0"):
            seen += 1
            if seen == k:
                out[i] = str((int(ch) + 1) % 10)
                return "".join(out)
    raise ValueError("string too short")


def ball_str(b, digits=60):
    """'[mid +/- rad]' for a real ball (str(b) at the ball precision is exact-enough; keep the engine's form)."""
    return str(b)


def leading_exponent(expansion_str):
    """The exponent alpha of the leading monomial x^alpha (or x^alpha*log(x)^k) of a local basis expansion."""
    s = expansion_str.replace(" ", "")
    m = re.match(r"^([-+]?)(?:\(?([0-9]+(?:/[0-9]+)?)\)?\*)?([a-zA-Z]\w*)?(?:\^\(?([-0-9/]+)\)?)?", s)
    # ore_algebra prints e.g. '1 + ...', 's + ...', 's^(5/4) + ...', 'log(s)*s + ...'
    if s.startswith("1+") or s == "1" or s.startswith("1-"):
        return Fr(0)
    mm = re.search(r"([a-zA-Z]\w*)\^\(?(-?[0-9]+(?:/[0-9]+)?)\)?", s.split("+")[0].split("-")[0] if not s.startswith("-") else s[1:].split("+")[0])
    if mm:
        return Fr(mm.group(2))
    head = s.split("+")[0]
    if re.fullmatch(r"(log\(\w+\)(\^\d+)?\*)?\w+", head) or re.fullmatch(r"\w+(\*log\(\w+\)(\^\d+)?)?", head):
        return Fr(1)
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description="ore_algebra transport oracle (see the module docstring)")
    ap.add_argument("--spec")
    ap.add_argument("--out", help="receipt JSON (O_EXCL)")
    ap.add_argument("--smoke-toy", action="store_true")
    ap.add_argument("--eps", default=None)
    ap.add_argument("--prec", type=int, default=None)
    ap.add_argument("--min-digits", type=float, default=None)
    ap.add_argument("--holonomic-dir", default=os.environ.get("HOLONOMIC_DIR", HOLONOMIC_DIR_DEFAULT))
    ap.add_argument("--compare-tool", default=COMPARE_DEFAULT)
    ap.add_argument("--planted-target", default=None, metavar="KEY", help="the control: flip the 25th significant digit of the target KEY's value (must FAIL by name)")
    a = ap.parse_args(argv)
    if not (a.spec or a.smoke_toy):
        ap.error("--spec SPEC.json or --smoke-toy")
    t_import = time.time()
    from sage.all import QQ, ZZ, PolynomialRing, RealBallField, ComplexBallField, QQbar, I  # noqa: E402 (sage.all FIRST)
    from ore_algebra import OreAlgebra  # noqa: E402
    H = load_by_path("holonomic_transport", os.path.join(a.holonomic_dir, "holonomic_transport.py"))
    C = load_by_path("compare_amflow_json", a.compare_tool)
    sv, ov = H.engine_versions()
    t_import = time.time() - t_import
    rec = {"receipt": "ore_transport_oracle", "stamp_utc": subprocess.run(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"], capture_output=True, text=True, check=True).stdout.strip(),
           "engine": {"sage": sv, "ore_algebra": ov, "pin": [SAGE_PIN, ORE_PIN], "python": sys.version.split()[0], "sys_prefix": sys.prefix},
           "house": {"holonomic_transport.py": {"path": os.path.join(a.holonomic_dir, "holonomic_transport.py"), "sha256": sha256_of(os.path.join(a.holonomic_dir, "holonomic_transport.py"))},
                     "compare_amflow_json.py": {"path": a.compare_tool, "sha256": sha256_of(a.compare_tool)}},
           "import_wall_s": round(t_import, 2)}
    print(f"engine SageMath {sv} / ore_algebra {ov} (pin {SAGE_PIN} / {ORE_PIN}); sys.prefix {sys.prefix}; imports {t_import:.2f} s")
    if not (sv.startswith(SAGE_PIN) and ov == ORE_PIN):
        sys.stderr.write(f"REFUSED (exit {EXIT_ENGINE}): engine {sv}/{ov} is not the pin {SAGE_PIN}/{ORE_PIN}\n")
        return EXIT_ENGINE
    t0 = time.time()
    conv = H.determine_ic_convention(prec=430)
    rec["ic_convention"] = conv
    print(f"IC convention determined this session: {conv} ({time.time()-t0:.2f} s)")

    if a.smoke_toy:
        prec = a.prec or 430
        RBF, CBF = RealBallField(prec), ComplexBallField(prec)
        Pol = PolynomialRing(QQ, "x"); x = Pol.gen()
        A = OreAlgebra(Pol, "Dx"); Dx = A.gen()
        toys = {}
        t1 = time.time()
        L1 = Dx ** 2 + 1
        M = L1.numerical_transition_matrix([QQ(0), QQ(1)], 1e-60, algorithm="naive")
        v = H.mat_vec(M, H.ic_vector([RBF(1), RBF(0)], conv))
        cos1 = RBF(1).cos()
        ok1 = bool(v[0].contains_exact(cos1) if hasattr(v[0], "contains_exact") else abs((v[0] - cos1).mid()) < 1e-50)
        toys["cos1"] = {"value": str(v[0]), "reference": str(cos1), "contains": ok1, "wall_s": round(time.time() - t1, 3)}
        print(f"toy 1: cos(1) by Ds^2+1 on [0 -> 1]: {v[0]} contains cos(1): {ok1}")
        t1 = time.time()
        L2 = x * Dx - QQ(1) / 2
        M2 = L2.numerical_transition_matrix([QQ(1), QQ(0)], 1e-60, algorithm="naive")
        mons = [str(m) for m in L2.local_basis_monomials(QQ(0))] if hasattr(L2, "local_basis_monomials") else [str(e) for e in L2.local_basis_expansions(QQ(0), order=1)]
        v2 = H.mat_vec(M2, [RBF(1)])
        ok2 = abs((v2[0] - 1).mid()) < 1e-50 and float(v2[0].rad()) < 1e-50
        toys["frac_landing"] = {"local_basis": mons, "coordinate": str(v2[0]), "expected": "1 on x^(1/2)", "PASS": bool(ok2), "wall_s": round(time.time() - t1, 3)}
        print(f"toy 2: x Dx - 1/2 on [1 -> 0]: local basis {mons}; coordinate {v2[0]} (expect 1): {'PASS' if ok2 else 'FAIL'}")
        rec["smoke_toy"] = toys
        rec["maxrss_MB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
        rec["wall_s"] = round(time.time() - t0, 3)
        ok = ok1 and ok2
        rec["OVERALL"] = "PASS" if ok else "FAIL"
        if a.out:
            with open(a.out, "x") as f:
                json.dump(rec, f, indent=1)
            print("receipt:", a.out)
        print("OVERALL:", rec["OVERALL"])
        return 0 if ok else EXIT_FAIL

    spec = json.load(open(a.spec))
    rec["spec"] = {"path": a.spec, "sha256": sha256_of(a.spec), "name": spec.get("name")}
    eps = float(a.eps or spec.get("eps", "1e-40"))
    prec = int(a.prec or spec.get("prec", 430))
    min_digits = float(a.min_digits if a.min_digits is not None else spec.get("min_digits", 30))
    var = spec.get("variable", "s")
    RBF, CBF = RealBallField(prec), ComplexBallField(prec)
    Pol = PolynomialRing(QQ, var); x = Pol.gen()
    A = OreAlgebra(Pol, "D" + var); D = A.gen()
    opspec = spec["operator"]
    if opspec.get("source") and os.path.exists(opspec["source"].get("path", "")):
        got = sha256_of(opspec["source"]["path"])
        if got != opspec["source"].get("sha256"):
            sys.stderr.write(f"REFUSED (exit {EXIT_PIN}): operator source {opspec['source']['path']} sha256 {got[:16]} != pin {str(opspec['source'].get('sha256'))[:16]}\n")
            return EXIT_PIN
        rec["operator_source_verified"] = {"path": opspec["source"]["path"], "sha256": got}
    L = build_operator(opspec, var, Pol, D)
    order = L.order()
    lead = L.leading_coefficient()
    rec["operator"] = {"form": opspec["form"], "chart": opspec.get("chart"), "order": int(order), "degree": int(L.degree()),
                       "leading_coefficient_real_roots": [float(r) for r, _m in lead.roots(RealBallField(53).base() if False else __import__("sage.all", fromlist=["RR"]).RR)]}
    print(f"operator: order {order}, degree {L.degree()}; leading coefficient real roots {rec['operator']['leading_coefficient_real_roots']}")
    # seed
    sd = spec["seed"]
    at = QQ(Fr(str(sd["at"])))
    if sd["kind"] == "bfkns_multinomial":
        import math as _m
        thr_guess = (sum(_m.isqrt(M) for M in sd["msq"])) ** 2
        N_eps = int(-_m.log(eps) / _m.log(float(Fr(str(sd["at"])) / thr_guess))) + 30   # tail bound below eps
        N = min(int(sd.get("N", 1200)), N_eps) if sd.get("N_cap_by_eps", True) else int(sd.get("N", 1200))
        jets, tails, thr, head = seed_jet_exact(sd["msq"], Fr(str(sd["at"])), N, order)
        ics_der = [RBF(QQ(j)).add_error(RBF(t)) for j, t in zip(jets, tails)]
        rec["seed"] = {"kind": "bfkns_multinomial", "msq": sd["msq"], "at": str(sd["at"]), "N": N, "threshold": str(thr), "c_head": head,
                       "tail_bounds": [f"{t:.3e}" for t in tails], "jet_radii_after_wrap": [str(b.rad()) for b in ics_der]}
        print(f"seed: exact multinomial series at {sd['at']} with N = {N} terms; c_n head {head}; tail bounds {[f'{t:.2e}' for t in tails]}")
    elif sd["kind"] == "explicit":
        vals = sd["values"]
        rads = sd.get("radii", ["0"] * len(vals))
        ics_der = [RBF(v).add_error(RBF(r)) for v, r in zip(vals, rads)]
        rec["seed"] = {"kind": "explicit", "at": str(sd["at"]), "values": vals, "radii": rads}
    else:
        ap.error("seed.kind must be bfkns_multinomial | explicit")
    ics = H.ic_vector(ics_der, conv)
    path = [parse_point(p, QQ, QQbar, I) for p in spec["path"]]
    assert path[0] == at, f"path must start at the seed point {at}, got {path[0]}"
    rec["path"] = [str(p) for p in spec["path"]]
    rec["eps"], rec["prec_bits"] = eps, prec
    print(f"path {spec['path']} eps {eps} prec {prec} bits")
    t1 = time.time()
    M = L.numerical_transition_matrix(path, eps, algorithm="naive")
    t_tm = time.time() - t1
    v = H.mat_vec(M, ics)
    rec["transition_wall_s"] = round(t_tm, 3)
    rec["transported"] = [{"re": str(CBF(e).real()), "im": str(CBF(e).imag())} for e in v]
    print(f"transition matrix: {t_tm:.2f} s; max entry radius {max(float(CBF(e).rad()) for e in v):.3e}")
    landing = spec.get("landing", {"kind": "regular_point"})
    results = {}
    if landing["kind"] == "local_basis":
        end = path[-1]
        try:
            mons = [str(m) for m in L.local_basis_monomials(end)]
        except Exception:  # noqa: BLE001
            mons = [str(e) for e in L.local_basis_expansions(end, order=1)]
        rec["local_basis_monomials"] = mons
        print("local basis monomials at the landing point:", mons)
        for alpha, fspec in landing["fractional"].items():
            al = Fr(alpha)
            idx = [i for i, m in enumerate(mons) if leading_exponent(m) == al]
            if len(idx) != 1:
                sys.stderr.write(f"cannot identify the local basis element x^{alpha} among {mons} (matches {idx})\n")
                return EXIT_USAGE
            coord = CBF(v[idx[0]])
            phase = CBF(0)
            if fspec.get("phase", "1") == "1":
                phase = CBF(1)
            else:
                m = re.fullmatch(r"exp\(i\*pi\*([0-9/]+)\)", fspec["phase"].replace(" ", ""))
                if not m:
                    ap.error("phase must be 1 or exp(i*pi*p/q)")
                q = Fr(m.group(1))
                phase = (CBF.pi() * CBF(0, 1) * CBF(QQ(q))).exp()
            val = coord * phase
            results[alpha] = {"basis_index": idx[0], "monomial": mons[idx[0]], "coordinate": {"re": str(coord.real()), "im": str(coord.imag())},
                              "phase": fspec.get("phase", "1"), "value": {"re": str(val.real()), "im": str(val.imag())}, "_ball": val}
    else:
        for k, e in enumerate(v):
            results[str(k)] = {"value": {"re": str(CBF(e).real()), "im": str(CBF(e).imag())}, "_ball": CBF(e)}
    # the pair form (compare_amflow_json's digit rule on the midpoints, the radius kept beside)
    import mpmath as mp
    mp.mp.dps = int(prec * 0.30103) + 20
    fails, compared = [], 0
    for key, tgt in spec.get("targets", {}).items():
        if a.planted_target == key:
            tgt = dict(tgt)
            tgt["value"] = flip_digit(tgt["value"], 25)
            tgt["label"] = tgt.get("label", "") + " [PLANTED: digit 25 flipped]"
            rec["planted_target"] = key
        if key not in results:
            results[key] = {"note": "no transported value for this key"}
            continue
        b = results[key]["_ball"]
        mid_re, mid_im = mp.mpf(str(b.real().mid())), mp.mpf(str(b.imag().mid()))
        rad = float(b.rad())
        ref_re = mp.mpf(str(tgt["value"]))
        ref_im = mp.mpf(str(tgt.get("im", "0")))
        d_re = C.matched_digits(ref_re, mid_re)
        d_re = float(d_re) if d_re is not None else float("inf")
        scale = max(float(abs(mid_re)), float(abs(mid_im)), 1e-300)
        ball_d = -mp.log10(mp.mpf(rad) / scale) if rad > 0 else mp.inf
        ball_d = float(ball_d)
        T = C.zero_threshold(mp.mpf(rad), mp.mpf(0), None)
        im_status = "NUMERICAL ZERO" if (abs(mid_im) <= max(T, mp.mpf(10) ** (-(int(tgt.get("digits", 30)) - 2))) and abs(ref_im) <= mp.mpf(10) ** (-(int(tgt.get("digits", 30)) - 2))) else None
        d_im = None if im_status else C.matched_digits(ref_im, mid_im)
        cap = float(tgt.get("digits", 10 ** 6))
        certified = min(d_re, ball_d, cap)
        ok = certified >= min_digits and (im_status is not None or (d_im is not None and float(d_im) >= min_digits))
        compared += 1
        results[key].update({"target": {k2: v2 for k2, v2 in tgt.items()}, "matched_digits_re": round(d_re, 2),
                             "im": im_status or f"matched_digits_im {float(d_im):.2f}", "ball_digits": round(ball_d, 2),
                             "certified_agreement_d": round(certified, 2), "bar": min_digits, "PASS": bool(ok)})
        print(f"  {key}: ours {results[key]['value']['re'][:70]}... rad {rad:.3e} ({ball_d:.1f} d); house {tgt.get('label','')} -> re {d_re:.2f} d, im {im_status or ('%.2f d' % float(d_im))}; certified {certified:.2f} d vs bar {min_digits}: {'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"{key} vs {tgt.get('label','target')}: certified {certified:.2f} d < {min_digits}")
    for r in results.values():
        r.pop("_ball", None)
    rec["results"] = results
    rec["OVERALL"] = "PASS" if (compared and not fails) else ("VACUOUS (0 compared)" if not compared else "FAIL")
    rec["fails"] = fails
    rec["wall_s"] = round(time.time() - t0, 3)
    rec["maxrss_MB"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    print(f"OVERALL: {rec['OVERALL']} ({compared} compared); wall {rec['wall_s']} s; maxrss {rec['maxrss_MB']} MB")
    if a.out:
        with open(a.out, "x") as f:
            json.dump(rec, f, indent=1)
        print("receipt:", a.out)
    return 0 if rec["OVERALL"] == "PASS" else EXIT_FAIL


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
r"""membound acceptance battery.

Default tier: exact and closed-form checks of every layer of the
machinery at reduced precision — Laurent algebra, the automatic sector-weight
derivation vs the hand-derived special case, the per-flag order-eps phase
table (ret/adv/fey),
sign-domain classing, the single-frequency Mellin closed form vs direct
quadrature (including the kappa2 finite-difference layer), the engine's
order-0 moments vs exact/regression targets, the Gamma-prefactor pole
assembly vs the closed-form pole, the boundary-vector output contract, and
the 3-frequency octant machinery (domain partition, exact coordinate and
Jacobian identities, coupler phase flags, order-0 engine vs the
exact 3/8 moment), and the user-core front door (`python3 -m membound
--spec examples/core_example.json` through the process farm vs the exact
order-0 value -1/15 and vs the serial library call).

--full: additionally runs the complete reduced-precision c_M validation gate
(gate_cM.py at dps 11: all six omega-cores through order eps^2, Laurent
assembly, IBP, spurious-pole cancellation, c_M = 1, PSLQ). Minutes, not
seconds; farmed over min(12, cpu) processes.

Exit 0 = all legs pass; exit 1 = any leg fails. If mpmath is missing the
battery SKIPs by name (exit 0) stating what to install.
"""
import argparse
import json
import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import mpmath as mp
except ImportError:
    print("SKIP membound battery: python package 'mpmath' not installed "
          "(pip install mpmath)")
    sys.exit(0)

from membound.core import (KCache, OmegaCoreSpec, make_integrand,
                           domain_classes, domains3, _coords3, _lone_pair,
                           compute_core)
from membound.spec import I1M_CORE, I3M_CORE, REGISTRY
from membound.laurent import L
from membound import assemble as asm
from membound.cli import load_spec

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE_SPEC = os.path.join(HERE, "examples", "core_example.json")

DPS = 9
PANELS = ("0", "0.5", "2", "8", "inf")
RESULTS = []


def leg(name):
    def wrap(fn):
        def run(*a, **kw):
            t0 = time.time()
            try:
                ok, detail = fn(*a, **kw)
            except Exception as ex:  # a crashed leg is a failed leg
                ok, detail = False, f"EXCEPTION {type(ex).__name__}: {ex}"
            dt = time.time() - t0
            print(f"  [{'PASS' if ok else 'FAIL'}] {name:<24} {detail}  ({dt:.1f}s)")
            RESULTS.append(ok)
            return ok
        return run
    return wrap


def digits(x, ref):
    """Agreement digits of x against a nonzero reference."""
    err = abs(x - ref) / abs(ref)
    return mp.inf if err == 0 else -mp.log10(err)


@leg("laurent-algebra")
def leg_laurent():
    a = L(-1, [1, 2, 3])            # eps^-1 + 2 + 3 eps
    b = L(0, [4, 5])                # 4 + 5 eps
    p = a * b                       # 4 eps^-1 + 13 + 22 eps + 15 eps^2
    s = a + b                       # eps^-1 + 6 + 8 eps
    checks = [
        p.coeff(-1) == 4, p.coeff(0) == 13, p.coeff(1) == 22, p.coeff(2) == 15,
        s.coeff(-1) == 1, s.coeff(0) == 6, s.coeff(1) == 8,
        (2 * a).coeff(0) == 4,
    ]
    t = L.taylor(mp.exp, 4)         # 1 + eps + eps^2/2 + ...
    checks.append(digits(t.coeff(2), mp.mpf(1) / 2) > DPS - 3)
    return all(checks), f"{sum(checks)}/{len(checks)} identities"


@leg("weight-engine-B1")
def leg_b1(kc):
    """Automatic order-1 weight vs the hand-derived -(3 l1 + 3 l2 + l3)."""
    x, y = mp.mpf("0.7"), mp.mpf("1.3")
    f = make_integrand(I1M_CORE, ("same", 1, 1), 1, kc)
    got = f(x, y)
    l1, l2, l3 = mp.log(x), mp.log(y), mp.log(x + y)
    base = x * y * mp.besselk(0, x) * mp.besselk(0, y) * mp.besselk(0, x + y)
    want = (3 * l1 + 3 * l2 + l3) * base   # re[(-1)*(A+iB)] = -A, A = -(3l1+3l2+l3)
    d = digits(got, want)
    return d > DPS - 3, f"{mp.nstr(d, 3)} digits vs hand formula"


@leg("prescription-phase")
def leg_phase(kc):
    """Per-flag order-eps phase table: order-1 imaginary weight -B in a mixed-sign sector.

    ret: B = 0 (opposite signs cancel); adv: B = 0; fey: B = 2 pi (constant
    phase -pi/2 per factor; the phase-only substitution, not a Feynman
    propagator, see core.py).
    """
    x, y = mp.mpf("0.31"), mp.mpf("2.1")
    outs = {}
    for presc in ("ret", "adv", "fey"):
        spec = OmegaCoreSpec(
            name=f"t-{presc}",
            ret=[(1, 1, 2, presc), (2, 1, 2, presc)],
            klines=[1, 2, 3], proj="im")
        f = make_integrand(spec, ("oppgt", 1, -1), 1, kc)
        outs[presc] = f(x, y)
    # coords in 'oppgt': |w2|=x, |w3|=y, |w1|=x+y; base carries sign s2 = -1
    base = -(x + y) * x * (mp.besselk(0, x + y) * mp.besselk(0, x)
                           * mp.besselk(0, y))
    want_fey = -2 * mp.pi * base           # im[(-1)*(A+iB)] = -B, B_fey = 2 pi
    tol = mp.mpf(10) ** (-(DPS - 3))
    ok = (abs(outs["ret"]) < tol * abs(base)
          and abs(outs["adv"]) < tol * abs(base)
          and digits(outs["fey"], want_fey) > DPS - 3)
    return ok, "ret/adv null, fey = 2pi weight"


@leg("domain-classes")
def leg_domains(kc):
    for order in (0, 1):
        cls = domain_classes(I1M_CORE, order, kc)
        if sum(m for _, m, _ in cls) != 6 or not (1 <= len(cls) <= 6):
            return False, f"order {order}: bad class partition {cls}"
    return True, "6 sign-domains partitioned, multiplicities conserved"


@leg("mellin-closed-form")
def leg_mellin(kc):
    """Single-frequency Mellin law vs direct quadrature, orders 0-2.

    Int w^{s0-1} |w|^{-e} K_e(w) dw  =  2^{s-2} G((s+e)/2) G((s-e)/2),
    s = s0 - e (s1 = 0).  Order 2 exercises the kappa2 FD layer.
    """
    s0 = mp.mpf(3)
    closed = asm.mellin_K(s0, 0, N=3)
    brk = [mp.mpf(0), mp.mpf(1), mp.mpf(5), mp.inf]
    n0 = mp.quad(lambda w: w ** (s0 - 1) * kc.K0(w), brk, maxdegree=6)
    n1 = mp.quad(lambda w: -w ** (s0 - 1) * mp.log(w) * kc.K0(w),
                 brk, maxdegree=6)
    n2 = mp.quad(lambda w: w ** (s0 - 1) * (mp.log(w) ** 2 / 2 * kc.K0(w)
                                            + kc.kap2(w)), brk, maxdegree=6)
    d0 = digits(n0, closed.coeff(0))
    d1 = digits(n1, closed.coeff(1))
    d2 = digits(n2, closed.coeff(2))
    ok = d0 > DPS - 3 and d1 > DPS - 3 and d2 > 6
    return ok, (f"orders 0/1/2 agree {mp.nstr(d0, 3)}/{mp.nstr(d1, 3)}/"
                f"{mp.nstr(d2, 3)} digits")


def order0_moment(spec, kc, brk):
    tot = mp.mpf(0)
    for rep, mult, _ in domain_classes(spec, 0, kc):
        f = make_integrand(spec, rep, 0, kc)
        tot += mult * mp.quad(f, brk, brk, maxdegree=3)
    return tot / (2 * mp.pi) ** 2


@leg("engine-order0-exact")
def leg_order0(kc, brk, store):
    """2-D engine end to end at order eps^0.

    j1^(0) has the exact value 1/30 (arXiv:2601.16256 memory sector).
    j3^(0) is checked against -8/105, a regression reference: this engine's
    own dps-16 run agrees with that rational to 15 digits; it is a drift
    detector, not an independently derived value.
    """
    j10 = order0_moment(I1M_CORE, kc, brk)
    j30 = order0_moment(I3M_CORE, kc, brk)
    store["j10"], store["j30"] = j10, j30
    d1 = digits(j10, mp.mpf(1) / 30)
    d3 = digits(j30, mp.mpf(-8) / 105)
    return (d1 > 6 and d3 > 6,
            f"j1(0) vs 1/30: {mp.nstr(d1, 3)}d; j3(0) vs -8/105: "
            f"{mp.nstr(d3, 3)}d (regression)")


@leg("pole-assembly")
def leg_pole(store):
    """Gamma prefactor x order-0 moment vs the closed-form leading pole."""
    zero = {0: store["j10"], 1: mp.mpf(0), 2: mp.mpf(0)}
    I1 = asm.prefactor_I1() * asm.core_to_L(zero)
    target = 1 / (15 * (8 * mp.pi) ** 4)
    d = digits(I1.coeff(-1), target)
    # structural smoke: IBP combination lands at eps^-4 with finite coeffs
    I3 = asm.prefactor_I3() * asm.core_to_L(
        {0: store["j30"], 1: mp.mpf(0), 2: mp.mpf(0)})
    I2 = asm.ibp_I2(I1, I3)
    ok = d > 6 and I2.nmin == -4 and mp.isfinite(I2.coeff(-4))
    return ok, f"I1[eps^-1] vs 1/(15(8pi)^4): {mp.nstr(d, 3)} digits; IBP smoke"


@leg("boundary-vector-contract")
def leg_contract():
    ent = REGISTRY[("P", "memory", "t40-static")]
    bv = asm.boundary_vector(
        "P", "memory", "t40-static", ent["prescription"],
        layers=[{"k": ent["layer_k"], "slots": [
            {"name": "I1M", "closure": "numeric(core)+gamma(prefactor)",
             "eps_series": {"-1": "1.0"}, "digits": 6.0, "notes": "battery"}]}],
        meta={"dps": DPS})
    fd, path = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    try:
        asm.emit(path, bv)
        with open(path) as fh:
            back = json.load(fh)
    finally:
        os.unlink(path)
    ok = (back["contract"] == "membound.boundary_vector.v1"
          and back["layer_basis"] == "|omega|^{-k*eps}"
          and back["layers"][0]["k"] == 7)
    return ok, "emit/load round-trip, contract fields, k=7 layer"


@leg("octant-domains")
def leg_domains3(kc):
    """3-frequency sign-octant decomposition: partition structure + classing."""
    doms = domains3()
    tags = [d[0] for d in doms]
    octants = {tuple(d[1:]) for d in doms}
    ok = (len(doms) == 14 and tags.count('same') == 2
          and tags.count('gt') == 6 and tags.count('lt') == 6
          and len(octants) == 8)
    for oc in octants:
        got = sorted(d[0] for d in doms if tuple(d[1:]) == oc)
        want = ['same'] if len(set(oc)) == 1 else ['gt', 'lt']
        ok = ok and got == want
    # classing conserves the 14 domains and never merges across panel kinds
    spec = OmegaCoreSpec("t3cls", ret=[(4, 1, 2, 'ret')], klines=[1, 2, 3],
                         proj='re', n_freq=3)
    for order in (0, 1):
        cls = domain_classes(spec, order, kc)
        if sum(m for _, m, _ in cls) != 14 or not (1 <= len(cls) <= 14):
            return False, f"order {order}: bad class partition"
        for _, _, members in cls:
            if len({'u' if d[0] == 'gt' else 'w' for d in members}) != 1:
                return False, "panel kinds merged into one class"
    return ok, "14 octant-domains (2 same + 6 gt + 6 lt), classes conserved"


def _invert3(dom, a):
    """Inverse of the octant coordinate map; None when the point falls on
    the other piece of a mixed octant."""
    tag = dom[0]
    if tag == 'same':
        return a[0], a[1], a[2]
    lone, i, j = _lone_pair(dom[1], dom[2], dom[3])
    al, ai, aj = a[lone - 1], a[i - 1], a[j - 1]
    if tag == 'gt':
        y = ai + aj - al
        return (al, y, ai / (ai + aj)) if y > 0 else None
    z = al - ai - aj
    return (ai, aj, z) if z > 0 else None


@leg("octant-coords")
def leg_coords3():
    """Exact geometry of the octant decomposition.

    Fixed signed points (w1, w2, w3): exactly one piece of the point's octant
    contains it, the forward map reproduces the magnitudes, and the coupler
    identity s1|w1| + s2|w2| + s3|w3| = s4|w4| holds to working precision.
    The 'gt' Jacobian is checked against a central-difference determinant
    (the map is quadratic in the coordinates, so central differences are
    exact up to rounding).
    """
    pts = [("0.9", "-0.4", "1.7"), ("0.2", "-1.5", "0.3"),
           ("-0.6", "-0.8", "-0.1"), ("1.1", "0.5", "-1.55")]
    tol = mp.mpf(10) ** (-(DPS - 2))
    checks = []
    for pw in pts:
        w = [mp.mpf(t) for t in pw]
        s = tuple(1 if t > 0 else -1 for t in w)
        a = [abs(t) for t in w]
        w4 = w[0] + w[1] + w[2]
        hits = []
        for d in domains3():
            if tuple(d[1:]) != s:
                continue
            xyz = _invert3(d, a)
            if xyz is not None:
                hits.append((d, xyz))
        checks.append(len(hits) == 1)
        if not hits:
            continue
        d, (x, y, z) = hits[0]
        a1, a2, a3, a4, s4, jac = _coords3(d[0], d[1], d[2], d[3], x, y, z)
        checks.append(max(abs(a1 - a[0]), abs(a2 - a[1]), abs(a3 - a[2])) < tol)
        checks.append(abs(s4 * a4 - w4) < tol)     # coupler identity
    # Jacobian of the simplex-parametrized piece vs finite differences
    x0, y0, u0, h = mp.mpf("0.7"), mp.mpf("1.3"), mp.mpf("0.45"), mp.mpf("0.001")
    M = mp.matrix(3, 3)
    for c, dxyz in enumerate([(h, 0, 0), (0, h, 0), (0, 0, h)]):
        up = _coords3('gt', 1, 1, -1, x0 + dxyz[0], y0 + dxyz[1], u0 + dxyz[2])
        dn = _coords3('gt', 1, 1, -1, x0 - dxyz[0], y0 - dxyz[1], u0 - dxyz[2])
        for r in range(3):
            M[r, c] = (up[r] - dn[r]) / (2 * h)
    checks.append(digits(abs(mp.det(M)), x0 + y0) > DPS - 3)
    return all(checks), f"{sum(checks)}/{len(checks)} geometry identities"


@leg("weight-engine-coupler")
def leg_weight3(kc):
    """Order-1 weight with a retarded factor ON THE COUPLER, mixed octants.

    Spec (0+ - i w4)^{1-2e} x K-lines 1,2,3, proj 're': the order-1 value is
    B * base with B = -2 phi(s4), and s4 = -1 on the (-,-,+) 'gt' piece,
    +1 on (+,+,-). The s4-odd base makes ret s4-even (both octants equal the
    hand formula), adv the global sign flip, and fey (the phase-only flag)
    the s4-follower: three distinguishable phase rules through the coupler sign.
    """
    x, y, u = mp.mpf("0.7"), mp.mpf("1.3"), mp.mpf("0.45")
    vals = {}
    for presc in ("ret", "adv", "fey"):
        spec = OmegaCoreSpec(f"t4-{presc}", ret=[(4, 1, 2, presc)],
                             klines=[1, 2, 3], proj='re', n_freq=3)
        for octn, dom in (("pp", ('gt', 1, 1, -1)), ("mm", ('gt', -1, -1, 1))):
            vals[presc, octn] = make_integrand(spec, dom, 1, kc)(x, y, u)
    a1, a2, a3, a4, s4, jac = _coords3('gt', 1, 1, -1, x, y, u)
    hand = (mp.pi * a4 * jac * mp.besselk(0, a1) * mp.besselk(0, a2)
            * mp.besselk(0, a3))
    checks = [
        digits(vals["ret", "pp"], hand) > DPS - 3,
        digits(vals["ret", "mm"], hand) > DPS - 3,
        digits(vals["adv", "pp"], -hand) > DPS - 3,
        digits(vals["adv", "mm"], -hand) > DPS - 3,
        digits(vals["fey", "pp"], hand) > DPS - 3,
        digits(vals["fey", "mm"], -hand) > DPS - 3,
    ]
    return all(checks), f"{sum(checks)}/6 coupler phase routings vs hand formula"


@leg("engine3-order0-exact")
def leg_order0_3():
    """3-frequency engine end to end at order eps^0, reduced precision.

    Core K0(|w1|) K0(|w2|) K0(|w3|) w4^2 with w4 = w1 + w2 + w3: the exact
    moment is J^(0) = 3/8 (expand w4^2; odd moments vanish; the Mellin law
    closes the rest). The uniform-octant class alone is held to a sharper bar
    against its separable closed form 3 A2 A0^2 + 6 A1^2 A0 (A_n from
    mellin_K). The 14-domain sum is quadrature-truncation-limited at
    maxdegree 2, but any coordinate, Jacobian, multiplicity, or measure
    defect shifts it at O(1), so the 2-digit bar is decisive.
    """
    spec = OmegaCoreSpec("T3", ret=[], klines=[1, 2, 3], poly={4: 2},
                         proj='re', n_freq=3)
    with mp.workdps(6):
        kc6 = KCache(6)
        brk = [mp.mpf(0), mp.inf]
        ubrk = [mp.mpf(0), mp.mpf(1)]
        f = make_integrand(spec, ('same', 1, 1, 1), 0, kc6)
        same = mp.quad(f, brk, brk, brk, maxdegree=2)
        A0 = asm.mellin_K(1, 0, N=1).coeff(0)
        A1 = asm.mellin_K(2, 0, N=1).coeff(0)
        A2 = asm.mellin_K(3, 0, N=1).coeff(0)
        d_same = digits(same, 3 * A2 * A0 ** 2 + 6 * A1 ** 2 * A0)
        tot = mp.mpf(0)
        for rep, mult, _ in domain_classes(spec, 0, kc6):
            f = make_integrand(spec, rep, 0, kc6)
            third = ubrk if rep[0] == 'gt' else brk
            tot += mult * mp.quad(f, brk, brk, third, maxdegree=2)
        d_full = digits(tot / (2 * mp.pi) ** 3, mp.mpf(3) / 8)
    return (d_same > mp.mpf("3.5") and d_full > 2,
            f"same-class {mp.nstr(d_same, 3)}d vs closed form; "
            f"J(0) vs 3/8: {mp.nstr(d_full, 3)}d (14-domain sum)")


def _run_cli(argv):
    """`python3 -m membound ...` as a user would run it (cwd = tools/)."""
    import subprocess
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, "-B", "-m", "membound"] + argv,
                          cwd=os.path.dirname(HERE), env=env,
                          capture_output=True, text=True)


@leg("cli-user-spec")
def leg_cli():
    """The --spec front door on the shipped non-registry example core.

    (a) known value: `python3 -m membound --spec examples/core_example.json`
        through the process farm (--jobs 2) at maxdegree 3 gives J^(0)
        against the exact -1/15 (Fourier route, see the example's comment;
        the same route yields the registry 1/30); the JSON record round-trips
        the spec and carries layer_k = 4.
    (b) plumbing: at a cheap setting (maxdegree 2, panels 0,1,inf) the farm's
        J^(0), J^(1) agree with the serial library call compute_core on the
        spec loaded by load_spec (same nodes, same arithmetic: spec loading,
        task farming, class multiplicities and the measure are what is tested).
    """
    checks, notes = [], []
    with open(EXAMPLE_SPEC) as fh:
        raw = json.load(fh)
    # (a) known value through the farm
    fd, out_a = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    fd, out_b = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    try:
        r = _run_cli(["--spec", EXAMPLE_SPEC, "--dps", str(DPS),
                      "--maxdegree", "3", "--panels", "0,2,inf",
                      "--order-max", "0", "--jobs", "2", "--quiet",
                      "--json-out", out_a])
        if r.returncode != 0:
            return False, f"CLI rc={r.returncode}: {r.stderr.strip()[-200:]}"
        with open(out_a) as fh:
            rec = json.load(fh)
        d0 = digits(mp.mpf(rec["J"]["0"]), mp.mpf(-1) / 15)
        checks += [rec["contract"] == "membound.core_moments.v1",
                   rec["layer_k"] == 4,
                   rec["spec"]["ret"] == raw["ret"],
                   rec["spec"]["klines"] == raw["klines"],
                   "J^(0) = " in r.stdout,
                   d0 > 6]
        notes.append(f"J(0) vs -1/15: {mp.nstr(d0, 3)}d")
        # (b) farm vs serial library at a cheap setting, orders 0 and 1
        r = _run_cli(["--spec", EXAMPLE_SPEC, "--dps", str(DPS),
                      "--maxdegree", "2", "--panels", "0,1,inf",
                      "--order-max", "1", "--jobs", "4", "--quiet",
                      "--json-out", out_b])
        if r.returncode != 0:
            return False, f"CLI rc={r.returncode}: {r.stderr.strip()[-200:]}"
        with open(out_b) as fh:
            rec = json.load(fh)
        spec = load_spec(EXAMPLE_SPEC)
        brk = [mp.mpf(0), mp.mpf(1), mp.inf]
        lib = compute_core(spec, DPS, brk, 2, log=lambda *a, **k: None,
                           orders=(0, 1))
        dd = [digits(mp.mpf(rec["J"][str(k)]), lib[k]) for k in (0, 1)]
        checks += [d > DPS - 3 for d in dd]
        notes.append(f"farm vs library J(0)/J(1): {mp.nstr(dd[0], 3)}/"
                     f"{mp.nstr(dd[1], 3)}d")
    finally:
        for pth in (out_a, out_b):
            if os.path.exists(pth):
                os.unlink(pth)
    return all(checks), f"{sum(checks)}/{len(checks)} " + "; ".join(notes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="also run the complete c_M validation gate "
                         "(gate_cM.py --dps 11; minutes)")
    args = ap.parse_args()

    t0 = time.time()
    mp.mp.dps = DPS
    print(f"# membound battery  dps={DPS}")
    kc = KCache(DPS)
    brk = [mp.inf if p == "inf" else mp.mpf(p) for p in PANELS]
    store = {}

    leg_laurent()
    leg_b1(kc)
    leg_phase(kc)
    leg_domains(kc)
    leg_mellin(kc)
    leg_order0(kc, brk, store)
    if "j10" in store:
        leg_pole(store)
    leg_contract()
    leg_domains3(kc)
    leg_coords3()
    leg_weight3(kc)
    leg_order0_3()
    leg_cli()

    ok = all(RESULTS)
    print(f"# default tier: {sum(RESULTS)}/{len(RESULTS)} legs pass  "
          f"wall={time.time() - t0:.1f}s")

    if args.full and ok:
        import subprocess
        gate = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "gate_cM.py")
        print("# --full: c_M validation gate (dps 11)")
        rc = subprocess.call([sys.executable, "-B", gate,
                              "--dps", "11", "--maxdegree", "4"])
        ok = ok and rc == 0

    print(f"# BATTERY: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

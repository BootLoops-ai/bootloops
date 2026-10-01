"""dist_cert test battery (module spec in each file's docstring).

LAW: NO census data. Inputs are (a) the pinned loci, (b) SYNTHETIC
candidates generated from the fixed plant-rule seeds (one PLANT_RULE JSON
per geometry under TERRIER_PLANTS_DIR; reference data not included in the
package), (c) synthetic metric models. Test set:
  1. pinned-locus points return distance 0 (EXACT);
  2. plants (fixed rules) come back CERTIFIED off-locus;
  3. two-chart consistency on synthetic points where both charts are defined.
Run: python3 tests/test_dist_cert.py   (single core, minutes-class).
"""
import json
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ball import CBall, RIv, ln_enclosure, sqrt_enclosure
import lattice_exact as lx
import zchart as zc
import wp_chart as wp

F = Fraction
# plant-rule directory from the environment; unset => loud refusal (fail-closed).
PLANTS_DIR = os.environ.get("TERRIER_PLANTS_DIR", "")
if not PLANTS_DIR:
    raise SystemExit("REFUSE: env TERRIER_PLANTS_DIR unset — must point at "
                     "the plant-rule directory (reference data not "
                     "included in the package)")
PASS = []


def check(name, cond):
    PASS.append((name, bool(cond)))
    print(("PASS " if cond else "FAIL ") + name)
    return cond


def plant_seed(geometry):
    with open(f"{PLANTS_DIR}/PLANT_RULE_{geometry}.json") as fh:
        return int(json.load(fh)["seed"], 16)


def test_ball_core():
    lo, hi = sqrt_enclosure(F(4))
    check("ball.sqrt exact at 4", lo == hi == 2)
    lo, hi = sqrt_enclosure(F(0))
    check("ball.sqrt exact at 0", lo == hi == 0)
    lo, hi = sqrt_enclosure(F(2))
    check("ball.sqrt(2) certified", lo * lo < 2 < hi * hi and hi - lo <= F(2, 2**63))
    llo, lhi = ln_enclosure(F(2))
    check("ball.ln(2) certified",
          llo < F(6931472, 10**7) and lhi > F(6931471, 10**7)
          and F(69, 100) < llo and lhi < F(70, 100))
    check("zchart sqrt2 pin", zc._SQRT2_LO**2 < 2 < zc._SQRT2_HI**2)
    d = CBall(F(1, 3), 0, 0)
    check("ball.dist exact 0 on identical exact points",
          (d - CBall(F(1, 3))).abs_iv() .lo == 0
          and (d - CBall(F(1, 3))).abs_iv().hi == 0)


def test_root_enumeration_named_frozen():
    # Complete-enumeration validation on the NAMED G3 frozen sublattices:
    # ADE root counts A17=306, A7+A7=112, D4+D4=48.
    check("G3-1 frozen A17 root count 306",
          len(lx.roots(lx.frozen_gram_G3(1))) == 306)
    check("G3-2 frozen A7+A7 root count 112",
          len(lx.roots(lx.frozen_gram_G3(2))) == 112)
    check("G3-3 frozen D4+D4 root count 48",
          len(lx.roots(lx.frozen_gram_G3(3))) == 48)


def test_exact_metric_pinned_locus_zero():
    # ON-locus: parent A2+A1; flux = the A1 root; flux-perp = A2 has roots.
    G = lx.direct_sum(lx.gram_An(2), lx.gram_An(1))
    d, wit = lx.dist_cert_exact(G, [[0, 0, 1]])
    ok = d == 0 and wit is not None and lx.inner(G, wit, wit) == -2 \
        and lx.inner(G, wit, [0, 0, 1]) == 0
    check("G2/G3 pinned-locus point -> distance 0 EXACT + root witness", ok)
    # Frozen-quotient rule: same candidate, but frozen basis = that A2
    # => every perp root lies in the frozen span => OFF (unfrozen-sector law).
    d2, cert2 = lx.dist_cert_exact(G, [[0, 0, 1]],
                                   frozen_basis=[[1, 0, 0], [0, 1, 0]])
    check("frozen-sublattice roots do NOT fire the predicate",
          d2 == "OFF" and cert2["roots_examined"] == 6)
    # CM stratum: exact finite-list membership.
    d3, w3 = lx.dist_cert_exact(G, [[1, 1, 1]], cm_key="tau=i",
                                frozen_cm_list=("tau=i",))
    check("CM finite-list membership -> distance 0 EXACT",
          d3 == 0 and w3 == ("CM", "tau=i"))


def test_exact_metric_plants_certified_off():
    # G3 plant rule (fixed seed): PCG64(seed) draws from a SYNTHETIC
    # candidate list, ACCEPT iff the exact predicate fails, stop at 2.
    from numpy.random import Generator, PCG64
    rng = Generator(PCG64(plant_seed("G3")))
    G = lx.direct_sum(lx.gram_An(2), lx.gram_An(1))  # synthetic parent
    cands = [[a, b, c] for a in range(-2, 3) for b in range(-2, 3)
             for c in range(-2, 3) if (a, b, c) != (0, 0, 0)]
    order = rng.permutation(len(cands))
    plants = []
    for idx in order:
        f = cands[int(idx)]
        on, _, _ = lx.on_locus_root(G, [f])
        if not on:  # rule: accept iff OFF by the exact predicate
            plants.append(f)
        if len(plants) == 2:
            break
    ok = len(plants) == 2
    for f in plants:
        d, cert = lx.dist_cert_exact(G, [f])
        ok = ok and d == "OFF" and "roots_examined" in cert
    check("G3 plants (fixed seed) -> CERTIFIED off-locus x2", ok)


def test_zchart_pinned_locus_zero():
    for geometry, loci in (("G4", zc.PINNED_LOCI_G4),
                           ("G5", zc.PINNED_LOCI_G5)):
        for s in loci:
            if s.minpoly is None:  # exact rational pinned point
                q = CBall(s.ball.re, s.ball.im, 0)
                d = zc.dist_to_locus(q, s)
                check(f"{geometry} {s.name} query=locus -> RIv[0,0] EXACT",
                      d.lo == 0 and d.hi == 0)
            else:  # algebraic pinned point: same-point identity path
                d = zc.dist_to_locus(s.ball, s, same_point=True)
                check(f"{geometry} {s.name} same algebraic point -> 0 EXACT",
                      d.lo == 0 and d.hi == 0)
    # minpoly sanity for (11 +- 5 sqrt2)/2: 4z^2 - 44z + 71 = 0 on enclosure
    for s in zc.PINNED_LOCI_G5:
        if s.minpoly:
            c, b, a = s.minpoly
            z = RIv(s.ball.re - s.ball.rad, s.ball.re + s.ball.rad)
            val = z * z * RIv(a) + z * RIv(b) + RIv(c)
            check(f"G5 {s.name} minpoly enclosure straddles 0",
                  val.lo <= 0 <= val.hi)
    check("G5 apparent singularity 7/4 EXCLUDED by convention pin",
          all(s.ball.re != F(7, 4) for s in zc.PINNED_LOCI_G5))


def test_zchart_verdicts_and_gates():
    eps0 = zc.EPS0["G4"]
    # certified ON (non-exact): center 1/25 + eps0/2, radius eps0/10
    v = CBall(F(1, 25) + eps0 / 2, 0, eps0 / 10)
    verdict, d, s = zc.classify([v], "G4")
    check("G4 ON within eps0 (certified upper bound)",
          verdict == "ON" and d.hi <= eps0 and s.name == "conifold_1/25")
    # validity gate FIRST: radius > eps0/10 => N_failed(ball-blowup)
    v = CBall(F(1, 25), 0, eps0 / 9)
    check("eps0-validity gate: rad > eps0/10 -> FAILED ball-blowup",
          zc.classify([v], "G4")[:2] == ("FAILED", "ball-blowup"))
    # straddle: enclosure crosses eps0 => indeterminate, counts N_failed
    v = CBall(F(1, 25) + eps0, 0, eps0 / 10)
    check("straddling enclosure -> FAILED indeterminate",
          zc.classify([v], "G4")[:2] == ("FAILED", "indeterminate"))
    # duality-quotient law: min over orbit representatives
    far, near = CBall(F(1, 2), 0, 0), CBall(F(1, 9) + eps0 / 2, 0, eps0 / 20)
    verdict, d, s = zc.classify([far, near], "G4")
    check("quotient dist = min over orbit reps",
          verdict == "ON" and s.name == "conifold_1/9")


def test_ball_metric_plants_certified_off():
    # G4/G5 plant rule (fixed seed): PCG64(seed); accept iff
    # float-proxy distance to every s in S > 10*eps0; plant success =
    # CERTIFIED off-locus (must-fail definition).
    from numpy.random import Generator, PCG64
    for geometry, loci in (("G4", zc.PINNED_LOCI_G4),
                           ("G5", zc.PINNED_LOCI_G5)):
        eps0 = zc.EPS0[geometry]
        rng = Generator(PCG64(plant_seed(geometry)))
        plants, tried = [], 0
        while len(plants) < 2 and tried < 10000:
            tried += 1
            z = F(int(rng.integers(-3000, 3000)), 1000)  # synthetic grid draw
            zf = float(z)
            prox = min(abs(zf - float(s.ball.re)) for s in loci)
            if prox > 10 * float(eps0):  # fixed acceptance rule
                plants.append(CBall(z, 0, eps0 / 10))
        ok = len(plants) == 2
        for p in plants:
            verdict, d, _ = zc.classify([p], geometry)
            ok = ok and verdict == "OFF" and d.lo > eps0
        check(f"{geometry} plants (fixed seed) -> CERTIFIED off-locus x2", ok)


def test_two_chart_consistency():
    # Test 3: synthetic points where BOTH charts are defined.
    # g == 1: WP == |dz| exactly; g == 4: WP == 2|dz|.
    a, b = CBall(F(1, 25) + F(1, 500)), CBall(F(1, 9) - F(1, 700))
    zd = (b - a).abs_iv()
    for c in (1, 4):
        ok, w, ref = wp.two_chart_consistency(a, b, wp.ConstantG(c), zd)
        check(f"two-chart consistency, constant g={c}: WP overlaps ref", ok)
        check(f"  WP enclosure tight, g={c}",
              w.hi - w.lo <= ref.hi - ref.lo + F(1, 2**40))
    # complex segment too (im component nonzero)
    a2, b2 = CBall(0, F(1, 10)), CBall(F(1, 5), F(-1, 10))
    zd2 = (b2 - a2).abs_iv()
    ok, w, ref = wp.two_chart_consistency(a2, b2, wp.ConstantG(1), zd2)
    check("two-chart consistency on complex synthetic segment", ok)


def test_wp_log_singular_no_cutoff():
    # Pinned regularization: log-divergent g integrable, NO cutoff.
    # g(t) = ln(1/(1-t)) on [0,1]: true WP length = Gamma(3/2) = sqrt(pi)/2.
    w = wp.wp_distance(CBall(0), CBall(1), wp.LogG(1, 0, F(1, 256)),
                       n_panels=255)
    ref_lo, ref_hi = F(8862269, 10**7), F(8862270, 10**7)  # sqrt(pi)/2
    check("WP log-singular head: enclosure contains sqrt(pi)/2 (finite)",
          w.lo < ref_lo and w.hi > ref_hi and w.hi < 2)
    check("WP log-singular head: certified width < 0.2", w.hi - w.lo < F(1, 5))
    # refinement monotonicity: finer head+panels -> tighter upper bound
    w2 = wp.wp_distance(CBall(0), CBall(1), wp.LogG(1, 0, F(1, 1024)),
                        n_panels=1023)
    check("WP refinement tightens the enclosure", w2.hi - w2.lo < w.hi - w.lo)


def main():
    test_ball_core()
    test_root_enumeration_named_frozen()
    test_exact_metric_pinned_locus_zero()
    test_exact_metric_plants_certified_off()
    test_zchart_pinned_locus_zero()
    test_zchart_verdicts_and_gates()
    test_ball_metric_plants_certified_off()
    test_two_chart_consistency()
    test_wp_log_singular_no_cutoff()
    n_fail = sum(1 for _, ok in PASS if not ok)
    print(f"\n{len(PASS) - n_fail}/{len(PASS)} checks passed")
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()

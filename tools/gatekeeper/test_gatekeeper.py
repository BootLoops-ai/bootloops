#!/usr/bin/env python3
"""
Self-tests for gatekeeper.

  T1  integrity scan catches a cross-tag duplicate (the duplicate-oracle
      bug, miniature)
  T2  held-out CV recovers exact rationals on a known function to ≥40 digits
  T3  Lyndon basis: Duval gives the Witt-formula count; transition matrix
      is square and invertible
  T4  ill-conditioned synthetic: plain LSQ+PSLQ FAILS, auto_condition (LLL)
      RECOVERS the exact integer relation  (the cond~1e36 wall, mini)
  T5  pslq_with_lll finds a 6-term integer relation directly
  T5b amf_laurent: coefficients[].order is ABSOLUTE (no lo+order double count)
  T6  ball-form values (AMFlow "[mid +/- rad]", verbatim fixture) with a
      nonzero midpoint are NOT flagged suspicious_zero / NOT quarantined;
      amf_arb reads them identically to the original split idiom
  T7  true zeros MUST flag: "0", "0.0", "0e0", "[0 +/- 1e-60]",
      "[0.000000 +/- 1e-60]", Arb's zero-midpoint "[+/- 1e-60]"
  T8  a ball with radius > |mid| is consistent_with_zero: its own field,
      NOT suspicious_zero, CLEAN by default, quarantines only under
      strict_zero_balls / --strict-zero-balls
  T9  CLI end-to-end (python -m gatekeeper.oracle_integrity on temp dirs):
      exit 0 CLEAN for T6's fixture, exit 2 for T7's zeros, exit 0 then 2
      for T8's ball without/with --strict-zero-balls
  T10 AMFlow out.json form (result[i].coefficients balls, integral.indices):
      one tagged dump + one out.json dump with two finite-ball entries ->
      n_files 2, n_records 3, no suspicious_zero, CLEAN (API and CLI)
  T11 an all-zero out.json entry is flagged by its ENTRY tag
      '<file-tag>#<indices>' (integral.family fallback / top-level tag);
      an out.json whose result list is empty stays 'empty coeffs'
  T12 the same numerics under two different indices across two out.json
      files -> ONE duplicate group (tags = the two entry tags, the second
      file quarantined); the same entry re-run in two files -> no duplicate
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

from gatekeeper import (scan_integrity, heldout_certify,      # noqa: E402
                            load_oracle_values, lyndon_basis,
                            pslq_with_lll)
from gatekeeper.heldout_cv import _duval_lyndon_words         # noqa: E402


# --------------------------------------------------------------------- T1
def test_integrity_catches_duplicate():
    tmp = tempfile.mkdtemp(prefix="oh_integ_")
    try:
        def w(name, tag, idx, s, t, re):
            json.dump({"tag": tag, "indices": idx, "s": s, "t": t,
                       "leading": -3,
                       "coeffs": {"-1": {"re": re, "im": "0"}}},
                      open(os.path.join(tmp, name), "w"))

        # 3 genuinely distinct
        w("m00_p1.json", "m00", [1, 0, 1], "-3", "-5", "1.111111111111111")
        w("m01_p1.json", "m01", [0, 1, 1], "-3", "-5", "2.222222222222222")
        w("m03_p1.json", "m03", [1, 1, 1], "-3", "-5", "3.333333333333333")
        # the bug: m02 was mis-collected as m00 — same numerics, DIFFERENT tag
        w("m02_p1.json", "m02", [1, 0, 0], "-3", "-5", "1.111111111111111")
        # an all-zero (eps_order footgun)
        w("m04_p1.json", "m04", [0, 0, 1], "-3", "-5", "0")

        rep = scan_integrity([tmp])
        assert rep.n_files == 5, rep.n_files
        assert len(rep.duplicates) == 1, rep.duplicates
        dup = rep.duplicates[0]
        assert set(dup["tags"]) == {"m00", "m02"}, dup["tags"]
        assert len(rep.suspicious_zero) == 1
        assert rep.suspicious_zero[0]["tag"] == "m04"
        assert rep.verdict == "QUARANTINE_NEEDED"
        # quarantine_set must drop exactly the *second* duplicate + the zero
        q = rep.quarantine_set()
        assert len(q) == 2
        print("T1 integrity-dedup: PASS  "
              f"(dup tags={dup['tags']}, zero={rep.suspicious_zero[0]['tag']})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------- T2
def test_heldout_cv_recovers_known():
    """Target:  f(s,t) = zeta(2)*log(-s) + 3*log(-t)^2  at weight 2.
    Basis = [log(-s), log(-t)^2].  6 Euclidean points."""
    mp.mp.dps = 100
    z2 = mp.zeta(2)
    pts = [(-3, -5), (-3, -7), (-3, -11), (-2, -13), (-5, -17), (-7, -19)]

    tmp = tempfile.mkdtemp(prefix="oh_cv_")
    try:
        for i, (s, t) in enumerate(pts):
            val = z2 * mp.log(-mp.mpf(s)) + 3 * mp.log(-mp.mpf(t)) ** 2
            json.dump({"tag": "demo", "indices": [1, 1],
                       "s": str(s), "t": str(t), "leading": 0,
                       "coeffs": {"2": {"re": mp.nstr(val, 90), "im": "0"}}},
                      open(os.path.join(tmp, f"demo_p{i}.json"), "w"))

        values, kin, _ = load_oracle_values([tmp])
        assert "demo" in values and len(values["demo"]) == 6

        names = ["zeta2*log(-s)", "log(-t)^2"]

        def basis_fn(pk, s, t):
            return [z2 * mp.log(-s), mp.log(-t) ** 2]

        res = heldout_certify(values, kin, weight=2,
                              basis_fn=basis_fn, basis_names=names,
                              gate=40, dps=100)
        r = res["demo"]
        assert r.certified, f"min held-out digits = {r.min}"
        assert r.min >= 40, r.min
        assert r.pslq_rationals == ["1", "3"], r.pslq_rationals
        print(f"T2 heldout-CV: PASS  (min={r.min:.1f}d  mean={r.mean:.1f}d  "
              f"rats={r.pslq_rationals})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------- T3
def _witt(k, n):
    """Necklace/Witt count: (1/n) Σ_{d|n} μ(n/d) k^d."""
    from sympy import mobius, divisors
    return sum(int(mobius(n // d)) * k ** d for d in divisors(n)) // n


def test_lyndon_basis_counts():
    # weight-3 over {a,b}: Witt(2,3) = (8-2)/3 = 2  →  words {aab, abb}
    w3 = list(_duval_lyndon_words(["a", "b"], 3))
    assert len(w3) == 2 == _witt(2, 3), w3
    assert set("".join(w) for w in w3) == {"aab", "abb"}, w3
    # spot-check a few more against Witt
    for k, n, expect in [(2, 1, 2), (2, 2, 1), (2, 4, 3), (2, 5, 6),
                         (3, 1, 3), (3, 2, 3), (3, 3, 8), (3, 4, 18)]:
        got = len(list(_duval_lyndon_words(list("abcde"[:k]), n)))
        assert got == expect == _witt(k, n), (k, n, got, expect)
    # full Lyndon shuffle-basis at weight 3 over {a,b}: dimension must be
    # 2^3 = 8 and the transition matrix must be square + invertible (Radford)
    names, words, T, _ = lyndon_basis(["a", "b"], 3)
    assert len(words) == 8 and len(names) == 8, (len(names), len(words))
    Tm = mp.matrix([[float(x) for x in row] for row in T])
    assert abs(mp.det(Tm)) > 0.5, "transition matrix singular?"
    print(f"T3 Lyndon/Witt: PASS  (w3 over 2 letters: {len(w3)} Lyndon words; "
          f"full shuffle basis dim=8, |det T|={int(abs(mp.det(Tm)))})")


# --------------------------------------------------------------------- T4
def test_illcond_lll_recovers():
    """The conditioning wall in miniature, through the public heldout_certify API.

    Six O(1) basis functions over 12 points, built from TWO independent
    primitives g,h plus a 1e-20 perturbation so the 6×12 evaluation matrix
    is full-rank but with cond ~1e36 (the classic pathological
    construction {1, 1+ε, 1+2ε, ...}).  Target = 3B0 - 2B1 + 5B3, an
    EXACT small-integer combo.

      auto_cond=False : normal-equations LSQ at dps=60 floors (held-out
                        digits ≈ 0, PSLQ returns junk) — conditioning, not
                        data, is the wall.
      auto_cond=True  : LLL on the column lattice finds the small-integer
                        near-dependencies, cond collapses, every held-out
                        fold certifies, and pslq_with_lll rationalises.
    """
    DPS = 60
    mp.mp.dps = DPS
    nB, nP = 6, 12
    g = lambda x: mp.log(1 + x)                          # noqa: E731
    h = lambda x: mp.polylog(2, x / 2)                   # noqa: E731
    # B_j = g + j·h/5   →  six O(1) columns, EXACT rank-2; cond → ∞ at dps.
    # (this is the "{1, 1+ε, 1+2ε,...}" pathological basis with ε→0;
    #  the real-world analogue is shuffle-redundant polylog products.)
    def Bj(j, x):
        return g(x) + (mp.mpf(j) / 5) * h(x)

    pts = [mp.mpf(k) / 13 for k in range(2, 2 + nP)]
    c_true = [3, -2, 0, 5, 0, 0]
    Btab, values, kin = {}, {"ill": {}}, {}
    for x in pts:
        row = [Bj(j, x) for j in range(nB)]
        tgt = sum(c_true[j] * row[j] for j in range(nB))
        pk = f"{mp.nstr(x, 12)}|0"
        Btab[pk] = row
        values["ill"][pk] = {2: mp.mpc(tgt, 0)}
        kin[pk] = (x, mp.mpf(0))
    names = [f"B{j}" for j in range(nB)]
    basis_fn = lambda pk, s=None, t=None: Btab[pk]       # noqa: E731

    # ---- BEFORE: auto_cond OFF (reproduce old behaviour) -------------
    r_off = heldout_certify(values, kin, weight=2, basis_fn=basis_fn,
                            basis_names=names, gate=30, dps=DPS,
                            auto_cond=False)["ill"]
    # ---- AFTER: auto_cond ON -----------------------------------------
    r_on = heldout_certify(values, kin, weight=2, basis_fn=basis_fn,
                           basis_names=names, gate=30, dps=DPS,
                           auto_cond=True, cond_target=1e6)["ill"]

    cb = mp.mpf(r_on.cond_before)
    ca = mp.mpf(r_on.cond_after)
    assert cb > mp.mpf("1e30"), f"setup not ill-conditioned: {r_on.cond_before}"
    # BEFORE: either LSQ is singular (min=0) or it limps through but PSLQ
    # cannot rationalise the garbage coefficients — that is the conditioning floor.
    assert r_off.pslq_rationals is None, \
        f"plain path unexpectedly rationalised: {r_off.pslq_rationals}"
    assert r_on.lll_applied and ca < mp.mpf("1e6"), \
        f"LLL didn't help: {r_on.cond_before} -> {r_on.cond_after}"
    assert r_on.certified and r_on.min >= 30, \
        f"LLL path failed to certify (min={r_on.min}d)"
    assert r_on.pslq_rationals is not None, "rationalisation failed"
    assert r_on.pslq_method in ("lll", "exact/lll", "exact"), r_on.pslq_method
    # verify the returned rational representative reproduces the target
    pk0 = next(iter(Btab))
    rec = mp.mpf(0)
    for j, s in enumerate(r_on.pslq_rationals):
        p, _, q = s.partition("/")
        rec += (mp.mpf(p) / mp.mpf(q or 1)) * Btab[pk0][j]
    err = abs(rec - values["ill"][pk0][2].real)
    assert err < mp.mpf(10) ** -25, \
        f"recovered rats don't reproduce target (err={mp.nstr(err,3)})"

    print("T4 ill-cond LLL (before/after via heldout_certify): PASS")
    print(f"    BEFORE  cond={r_on.cond_before:>10s}  "
          f"pslq_rationals={r_off.pslq_rationals}  "
          f"<-- PSLQ FLOORS (conditioning, not data)")
    print(f"    AFTER   cond={r_on.cond_after:>10s}  "
          f"held-out min={r_on.min:5.1f}d  pslq_method={r_on.pslq_method}  "
          f"rats={r_on.pslq_rationals}  <-- RECOVERS")
    if r_on.note:
        print(f"            note: {r_on.note}")


# --------------------------------------------------------------------- T5
def test_pslq_with_lll_multibasis():
    """Direct LLL relation:  7·ζ(3) - 2·π²log2 + 0·log³2 - 3·Li₃(1/2)
       has NO relation, but  8·Li₃(1/2) - 7·ζ(3) + π²log2·... — use a
       known one:  35·ζ(3) - 2·π²·log2·? ... simpler: synthetic exact."""
    mp.mp.dps = 80
    b = [mp.zeta(3), mp.pi ** 2 * mp.log(2), mp.log(2) ** 3,
         mp.mpf(1), mp.catalan * mp.log(2)]
    target = 7 * b[0] - 2 * b[1] + 11 * b[2] - 5 * b[3] + 3 * b[4]
    rel, method = pslq_with_lll([target] + b, prec=60, max_coeff_bits=20)
    assert rel is not None and method == "lll", (rel, method)
    # normalise sign so rel[0] > 0
    if rel[0] < 0:
        rel = [-r for r in rel]
    # expect a0=1 (or divides), and -a_i/a0 = true coeffs
    a0 = rel[0]
    got = [-r / a0 for r in rel[1:]]
    assert got == [7, -2, 11, -5, 3], (rel, got)
    print(f"T5 pslq_with_lll: PASS  (method={method}, rel={rel})")


# -------------------------------------------------------------------- T5b
def test_amf_laurent_absolute_order():
    """coefficients[].order is ABSOLUTE; leading_order is redundant min().
    Guards the lo+c['order'] double-count: the offset must be applied once."""
    from gatekeeper import amf_laurent
    entry = {"leading_order": -2,
             "coefficients": [
                 {"order": -2, "value": {"re": "[1.5 +/- 1e-50]", "im": "0"}},
                 {"order": -1, "value": {"re": "[2.5 +/- 1e-50]", "im": "0"}},
                 {"order":  0, "value": {"re": "[3.5 +/- 1e-50]", "im": "0"}},
             ]}
    L = amf_laurent(entry)
    assert set(L) == {-2, -1, 0}, f"expected {{-2,-1,0}}, got {set(L)}"
    assert abs(L[-2].real - 1.5) < 1e-10
    print(f"T5b amf_laurent absolute order: PASS  (orders={sorted(L)})")


# ------------------------------------------------------- T6-T9 (ball forms)
FIXTURE_BALL = os.path.join(os.path.dirname(HERE), "fixtures", "gatekeeper",
                            "amflow_ball_nonzero.json")
_OI_MOD = ["-m", "gatekeeper.oracle_integrity"]


def _write_dump(tmp, name, tag, coeffs):
    json.dump({"tag": tag, "indices": [1, 1, 1], "s": "-3", "t": "-5",
               "leading": min(int(k) for k in coeffs), "coeffs": coeffs},
              open(os.path.join(tmp, name), "w"))


def _old_split_idiom(s):
    """The pre-fix amf_arb body, verbatim, as the identity oracle for T6."""
    s = str(s)
    if s.startswith("["):
        s = s.split(" +/- ")[0].lstrip("[")
    return mp.mpf(s)


def test_ballform_nonzero_not_flagged():
    from gatekeeper import amf_arb, amf_ball
    fx = json.load(open(FIXTURE_BALL))
    balls = [v["re"] for v in fx["coeffs"].values() if v["re"].startswith("[")]
    assert len(balls) >= 10, f"fixture lost its ball strings ({len(balls)})"
    assert all(" +/- " in b for b in balls)
    # (i) the scanner: a real nonzero ball-form dump is NOT suspicious_zero
    tmp = tempfile.mkdtemp(prefix="oh_ball_")
    try:
        shutil.copy(FIXTURE_BALL, os.path.join(tmp, "g70_p0.json"))
        rep = scan_integrity([tmp])
        assert rep.n_files == 1
        assert rep.suspicious_zero == [], rep.suspicious_zero
        assert rep.consistent_with_zero == [], rep.consistent_with_zero
        assert rep.verdict == "CLEAN", rep.verdict
        assert rep.quarantine_set() == set()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    # (ii) the shared parser keeps the full midpoint (no float truncation)
    mp.mp.dps = 130
    mid_s, rad_s = amf_ball(balls[0])
    assert mid_s == balls[0].split(" +/- ")[0].lstrip("["), mid_s
    assert rad_s == balls[0].split(" +/- ")[1].rstrip("]"), rad_s
    assert len(mid_s.replace("-", "").replace(".", "")) >= 100
    # (iii) amf_arb is IDENTICAL to the original split idiom on every form
    for s in balls + ["0", "0.25", "-1.0e-5", "[1.5 +/- 1e-50]"]:
        assert amf_arb(s) == _old_split_idiom(s), s
    # (iv) forms the old idiom could not read: Arb zero-midpoint + interval
    assert amf_ball("[+/- 1.01e-60]") == ("0", "1.01e-60")
    assert amf_ball("1.5 +/- 1e-9") == ("1.5", "1e-9")
    m, r = amf_ball("[0.24, 0.26]")
    assert abs(mp.mpf(m) - mp.mpf("0.25")) < mp.mpf("1e-15"), m
    assert abs(mp.mpf(r) - mp.mpf("0.01")) < mp.mpf("1e-15"), r
    print(f"T6 ball-form nonzero NOT flagged: PASS  ({len(balls)} verbatim "
          f"AMFlow balls, midpoint kept to {len(mid_s)-2} digits, "
          f"amf_arb == old idiom on {len(balls)+4} strings, verdict=CLEAN)")


ZERO_FORMS = ["0", "0.0", "0e0", "[0 +/- 1e-60]", "[0.000000 +/- 1e-60]",
              "[+/- 1e-60]"]


def test_true_zeros_flag():
    tmp = tempfile.mkdtemp(prefix="oh_zero_")
    try:
        for i, z in enumerate(ZERO_FORMS):
            _write_dump(tmp, f"z{i}.json", f"z{i}",
                        {"-2": {"re": z, "im": "0"}, "-1": {"re": "0", "im": z}})
        # control: one nonzero plain-float file in the same dir must survive
        _write_dump(tmp, "ok.json", "ok", {"-1": {"re": "1.25", "im": "0"}})
        rep = scan_integrity([tmp])
        flagged = sorted(s["tag"] for s in rep.suspicious_zero)
        assert flagged == [f"z{i}" for i in range(len(ZERO_FORMS))], flagged
        assert rep.verdict == "QUARANTINE_NEEDED"
        q = rep.quarantine_set()
        assert len(q) == len(ZERO_FORMS) and not any(
            f.endswith("ok.json") for f in q), q
        print(f"T7 true zeros MUST flag: PASS  ({len(ZERO_FORMS)} forms "
              f"{ZERO_FORMS} all suspicious_zero; plain-float control kept)")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


CWZ_BALL = "[1e-40 +/- 1e-30]"


def test_consistent_with_zero_status():
    tmp = tempfile.mkdtemp(prefix="oh_cwz_")
    try:
        _write_dump(tmp, "c0.json", "c0",
                    {"-2": {"re": "1.5", "im": "0"},
                     "-1": {"re": CWZ_BALL, "im": "0"}})
        rep = scan_integrity([tmp])
        assert rep.suspicious_zero == [], rep.suspicious_zero
        assert len(rep.consistent_with_zero) == 1, rep.consistent_with_zero
        c = rep.consistent_with_zero[0]
        assert (c["tag"], c["order"], c["part"]) == ("c0", "-1", "re"), c
        assert mp.mpf(c["mid"]) == mp.mpf("1e-40") and \
            mp.mpf(c["rad"]) == mp.mpf("1e-30"), c
        assert rep.verdict == "CLEAN" and rep.quarantine_set() == set()
        # strict mode promotes it
        rep_s = scan_integrity([tmp], strict_zero_balls=True)
        assert rep_s.suspicious_zero == []
        assert len(rep_s.consistent_with_zero) == 1
        assert rep_s.verdict == "QUARANTINE_NEEDED"
        assert rep_s.quarantine_set() == {os.path.join(tmp, "c0.json")}
        # a tight ball (rad << |mid|) is neither
        _write_dump(tmp, "t0.json", "t0",
                    {"-1": {"re": "[1e-40 +/- 1e-60]", "im": "0"}})
        rep_t = scan_integrity([os.path.join(tmp, "t0.json")])
        assert rep_t.suspicious_zero == [] and rep_t.consistent_with_zero == []
        print(f"T8 consistent_with_zero: PASS  ({CWZ_BALL} -> own field "
              f"(mid={c['mid']}, rad={c['rad']}), not suspicious_zero; "
              f"default verdict CLEAN, strict -> QUARANTINE_NEEDED)")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _cli(dirs, *flags):
    tmp_rep = tempfile.mktemp(prefix="oh_rep_", suffix=".json")
    p = subprocess.run([sys.executable, *_OI_MOD, *dirs, "--report", tmp_rep,
                        *flags], cwd=os.path.dirname(HERE),
                       capture_output=True, text=True, timeout=120)
    rep = json.load(open(tmp_rep)) if os.path.exists(tmp_rep) else None
    if os.path.exists(tmp_rep):
        os.remove(tmp_rep)
    return p.returncode, p.stdout + p.stderr, rep


def test_cli_end_to_end():
    tmpA = tempfile.mkdtemp(prefix="oh_cliA_")
    tmpB = tempfile.mkdtemp(prefix="oh_cliB_")
    tmpC = tempfile.mkdtemp(prefix="oh_cliC_")
    try:
        shutil.copy(FIXTURE_BALL, os.path.join(tmpA, "g70_p0.json"))
        for i, z in enumerate(ZERO_FORMS):
            _write_dump(tmpB, f"z{i}.json", f"z{i}", {"-1": {"re": z, "im": "0"}})
        _write_dump(tmpC, "c0.json", "c0", {"-1": {"re": CWZ_BALL, "im": "0"}})

        rcA, outA, repA = _cli([tmpA])
        assert rcA == 0 and "VERDICT      : CLEAN" in outA, (rcA, outA)
        assert repA["verdict"] == "CLEAN" and repA["suspicious_zero"] == []
        rcB, outB, repB = _cli([tmpB])
        assert rcB == 2 and "QUARANTINE_NEEDED" in outB, (rcB, outB)
        assert len(repB["suspicious_zero"]) == len(ZERO_FORMS)
        rcC, outC, repC = _cli([tmpC])
        assert rcC == 0 and repC["verdict"] == "CLEAN" and \
            len(repC["consistent_with_zero"]) == 1 and \
            repC["strict_zero_balls"] is False, (rcC, outC)
        assert "consistent-0 : 1 ball(s)" in outC, outC
        rcS, outS, repS = _cli([tmpC], "--strict-zero-balls")
        assert rcS == 2 and repS["verdict"] == "QUARANTINE_NEEDED" and \
            repS["strict_zero_balls"] is True, (rcS, outS)
        # no file moved (no --quarantine given)
        assert os.path.exists(os.path.join(tmpC, "c0.json"))
        print(f"T9 CLI end-to-end: PASS  (ball fixture rc={rcA} CLEAN; "
              f"{len(ZERO_FORMS)} zero forms rc={rcB}; cwz rc={rcC} default, "
              f"rc={rcS} --strict-zero-balls)")
    finally:
        for d in (tmpA, tmpB, tmpC):
            shutil.rmtree(d, ignore_errors=True)


# ------------------------------------------- T10-T12 (AMFlow out.json form)
_OUT_OPTIONS = {"d0": "2", "working_pre": 200, "x_order": 100}
OUT_BALL_A = "[1.25000000000000000000000000000000000000000000 +/- 3e-45]"
OUT_BALL_B = "[-3.86248153090527578720902875031351820625741058 +/- 8.06e-45]"
OUT_BALL_C = "[9.10431585514040789453326108942443432117447190 +/- 3.07e-45]"
OUT_BALL_D = "[0.5 +/- 1e-50]"


def _out_entry(indices, re_by_order, family="fam"):
    """One result[] entry in amflow-cpp's solve_integrals shape:
    {coefficients: [{order, value: {im, re}}], integral: {family, indices},
    leading_order}."""
    return {"coefficients": [{"order": o, "value": {"im": "0", "re": s}}
                             for o, s in sorted(re_by_order.items())],
            "integral": {"family": family, "indices": list(indices)},
            "leading_order": min(re_by_order)}


def _write_out(tmp, name, entries, **top):
    """An AMFlow out.json dump {mode, options, result: [...]} (+ top-level
    keys, e.g. tag=...)."""
    d = {"mode": "solve_integrals", "options": dict(_OUT_OPTIONS),
         "result": entries}
    d.update(top)
    json.dump(d, open(os.path.join(tmp, name), "w"), indent=1)


def test_outjson_two_form_fixture():
    tmp = tempfile.mkdtemp(prefix="oh_out2f_")
    try:
        _write_dump(tmp, "tagged_p0.json", "m00",
                    {"-1": {"re": "1.5", "im": "0"}})
        _write_out(tmp, "out_p0.json",
                   [_out_entry([1, 0, 1], {-2: OUT_BALL_A, -1: OUT_BALL_B,
                                           0: OUT_BALL_C}),
                    _out_entry([0, 1, 1], {-1: OUT_BALL_D, 0: OUT_BALL_A})])
        rep = scan_integrity([tmp])
        assert rep.n_files == 2, rep.n_files
        assert rep.n_records == 3, rep.n_records
        assert rep.suspicious_zero == [], rep.suspicious_zero
        assert rep.consistent_with_zero == [], rep.consistent_with_zero
        assert rep.duplicates == [] and rep.ok_count == 3, (rep.duplicates,
                                                            rep.ok_count)
        assert rep.verdict == "CLEAN" and rep.quarantine_set() == set()
        rc, out, repj = _cli([tmp])
        assert rc == 0 and "VERDICT      : CLEAN" in out, (rc, out)
        assert repj["n_files"] == 2 and repj["n_records"] == 3, repj
        assert repj["suspicious_zero"] == [] and repj["verdict"] == "CLEAN"
        print("T10 out.json two-form fixture: PASS  (1 tagged + 1 out.json "
              f"of 2 entries: n_files={rep.n_files} n_records={rep.n_records} "
              f"suspicious_zero={len(rep.suspicious_zero)} verdict={rep.verdict}; "
              f"CLI rc={rc})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_outjson_zero_entry_and_empty_result():
    tmp = tempfile.mkdtemp(prefix="oh_outz_")
    try:
        _write_out(tmp, "out_z.json",
                   [_out_entry([1, 0, 1], {-1: OUT_BALL_A}),
                    _out_entry([0, 1, 1], {-2: "[0 +/- 1e-60]",
                                           -1: "[+/- 1e-60]", 0: "0"})])
        _write_out(tmp, "out_empty.json", [])
        _write_out(tmp, "out_tagged.json",
                   [_out_entry([1, 1, 0], {0: "0"})], tag="g7")
        rep = scan_integrity([tmp])
        assert rep.n_files == 3 and rep.n_records == 4, (rep.n_files,
                                                         rep.n_records)
        flagged = {s["tag"]: s["reason"] for s in rep.suspicious_zero}
        assert len(rep.suspicious_zero) == 3 and len(flagged) == 3, flagged
        # the all-zero entry, by its entry tag (integral.family fallback)
        assert flagged["fam#0,1,1"].startswith("all coeffs |.|<"), flagged
        # a top-level tag wins over the family
        assert flagged["g7#1,1,0"].startswith("all coeffs |.|<"), flagged
        # the empty result list: one file-level record, "empty coeffs"
        assert flagged[None] == "empty coeffs", flagged
        files = {os.path.basename(s["file"]) for s in rep.suspicious_zero}
        assert files == {"out_z.json", "out_empty.json", "out_tagged.json"}
        assert rep.verdict == "QUARANTINE_NEEDED"
        # quarantine is file-granular: the flagged entry moves its file
        assert {os.path.basename(f) for f in rep.quarantine_set()} == files
        print("T11 out.json zero entry + empty result: PASS  "
              f"(flagged by entry tag: {sorted(str(k) for k in flagged)})")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_outjson_duplicate_across_files():
    tmpA = tempfile.mkdtemp(prefix="oh_outdupA_")
    tmpB = tempfile.mkdtemp(prefix="oh_outdupB_")
    try:
        same = {-2: OUT_BALL_A, -1: OUT_BALL_B}
        _write_out(tmpA, "f1.json",
                   [_out_entry([1, 0, 1], same),
                    _out_entry([0, 1, 1], {-1: OUT_BALL_D, 0: OUT_BALL_C})])
        # f2's entry [1,1,0] carries f1's [1,0,1] numerics: the mis-collection
        _write_out(tmpA, "f2.json",
                   [_out_entry([1, 1, 1], {0: OUT_BALL_C}),
                    _out_entry([1, 1, 0], same)])
        rep = scan_integrity([tmpA])
        assert rep.n_files == 2 and rep.n_records == 4
        assert len(rep.duplicates) == 1, rep.duplicates
        g = rep.duplicates[0]
        assert g["n"] == 2 and set(g["tags"]) == {"fam#1,0,1", "fam#1,1,0"}, g
        assert {os.path.basename(f) for f in g["files"]} == {"f1.json", "f2.json"}
        assert rep.ok_count == 2 and rep.suspicious_zero == []
        assert rep.verdict == "QUARANTINE_NEEDED"
        assert {os.path.basename(f) for f in rep.quarantine_set()} == {"f2.json"}
        # control: the SAME entry (same indices, same numerics) in two files
        # is the re-run class — no duplicate, as for tagged files
        _write_out(tmpB, "f1.json", [_out_entry([1, 0, 1], same)])
        _write_out(tmpB, "f1_again.json", [_out_entry([1, 0, 1], same)])
        rep_b = scan_integrity([tmpB])
        assert rep_b.n_records == 2 and rep_b.duplicates == [], rep_b.duplicates
        assert rep_b.ok_count == 2 and rep_b.verdict == "CLEAN"
        print("T12 out.json duplicate across files: PASS  "
              f"(1 group, tags={sorted(g['tags'])}, n={g['n']}; "
              f"re-run control: 0 groups, ok_count={rep_b.ok_count})")
    finally:
        shutil.rmtree(tmpA, ignore_errors=True)
        shutil.rmtree(tmpB, ignore_errors=True)


if __name__ == "__main__":
    test_integrity_catches_duplicate()
    test_heldout_cv_recovers_known()
    test_lyndon_basis_counts()
    test_illcond_lll_recovers()
    test_pslq_with_lll_multibasis()
    test_amf_laurent_absolute_order()
    test_ballform_nonzero_not_flagged()
    test_true_zeros_flag()
    test_consistent_with_zero_status()
    test_cli_end_to_end()
    test_outjson_two_form_fixture()
    test_outjson_zero_entry_and_empty_result()
    test_outjson_duplicate_across_files()
    print("\nALL TESTS PASS")

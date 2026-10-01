#!/usr/bin/env python3
"""CLINCH acceptance battery (legs L1-L8). V-style: every leg prints PASS/FAIL; any FAIL =>
exit 1; exactly one OVERALL line. Run from a scratch cwd (refuses the
tool tree and the protected roots; extend via CLINCH_DENY_ROOTS).

  L1  planted-optimum: synthetic hierarchical problem, EXACT rational
      optimum (Fraction elimination) -> CERTIFIED-INTERIOR and the
      certified box contains the exact optimum. Plus: a corner
      candidate with the bound declared -> CERTIFIED-CORNER with strict
      KKT margins; the same candidate WITHOUT the bound -> fat-ball
      REFUSAL (the true-but-weak-statement control); an unported space ->
      OUT-OF-SCOPE.
  L2  planted-saddle: flipped latent-prior curvature -> the Krawczyk leg
      still certifies the stationary point, the PD leg REFUSES naming the
      block and the (float-diagnostic) negative direction.
  L3  THE FLAT-VALLEY CONTROL (v3.6 unpinned vs v3.7 pinned): receipt-
      consuming leg over the centered-space (oracle_v36, auxiliary-border
      KKT extension) results — the stall optimum REFUSED pointing at the
      advance-level direction, the pinned KKT system Krawczyk-certified,
      identity/FD gates in bar. (The reduced-Hessian PD statement for
      centered spaces is proven self-contained by L8.)
  L4  oracle-tamper/mutation controls: (a) synthetic math-sign-flip at the
      pristine certificate's radius -> REFUSED; (b) bit-flipped scratch
      copy -> clinch.verify() refuses typed; (c) poisoned __pycache__ in
      a copy -> NO effect (purged at import; source bytes rule); (d) a
      sign flip planted in the REAL fit13 oracle inputs moves the arb nll
      far outside the receipted gate value -> the identity gate catches
      model drift.
  L5  precision honesty: synthetic certificates at dps 30/50/80 agree
      (same verdict, margins within 10x); a starved run (prec=8 bits)
      REFUSES — never a loose certificate.
  L6  flagship receipts: clinch_fit13_null / clinch_fit17_null CERT.json
      verified — direct certificate an HONEST NAMED REFUSAL (the
      under-polish finding), polished-center certificate
      CERTIFIED-INTERIOR with contraction margin < 1 and PD margins > 0,
      oracle identity gates adjudicated ok, ball-center shas + statements
      present. (Producer: battery/flagship.py; certificates are minutes —
      the battery verifies the stored receipts rather than recomputing.)
  L8  centered-space reduced-Hessian PD (the rank-4 congruence
      extension, clinch.kkt_pd): the toy CENTERED fixture (mean-zero
      centering, rank-2 auxiliary-border KKT extension) gets the FULL
      certificate — Krawczyk KKT point + PD_REDUCED_CERTIFIED with the
      certified inertia, exact rational KKT point contained in the box.
      Controls: a flipped border prior -> typed kkt-signature refusal
      (one negative direction too many, counts named); a flipped latent
      prior -> midpoint-Cholesky refusal naming the block.
"""
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)
sys.path.insert(0, HERE)
sys.path.insert(0, os.environ.get("CLINCH_BALLER_DIR")
                or os.path.normpath(os.path.join(TOOL, "..", "baller")))

# protected roots are env-extensible (colon-separated CLINCH_DENY_ROOTS);
# the tool tree itself is always refused; machine-specific prefixes come
# from the environment, never the source
BANNED = (TOOL, os.path.dirname(TOOL)) + tuple(
    p for p in os.environ.get("CLINCH_DENY_ROOTS", "").split(":") if p)
cwd = os.getcwd()
if any(os.path.realpath(cwd).startswith(os.path.realpath(b)) for b in BANNED):
    print(f"REFUSED: run from a scratch cwd, not {cwd}")
    sys.exit(2)

# reference results tree (CLINCH_REFERENCE_DIR) env-required for the
# receipt-reading legs (named skip/refusal in the legs when unset; no
# machine-local default)
RCPT = (os.environ.get("CLINCH_REFERENCE_DIR")
        or os.environ.get("CLINCH_RECEIPTS_ROOT")
        or os.environ.get("CLINCH_LANE_RECEIPTS") or "")
RESULTS = []


def leg(name):
    def deco(fn):
        RESULTS.append((name, fn))
        return fn
    return deco


def _setup():
    import numpy as np
    import clinch
    clinch.verify(quiet=True)
    from synthetic import SynthOracle, exact_optimum, NZ, NG
    from clinch import engine
    zs_ex, g_ex = exact_optimum()
    theta = [float(v) for v in zs_ex[0] + zs_ex[1] + g_ex]
    return np, clinch, SynthOracle, (zs_ex, g_ex), theta, engine, NZ, NG


@leg("L1 planted optimum + corner + out-of-scope")
def l1():
    from fractions import Fraction
    np, clinch, SynthOracle, (zs_ex, g_ex), theta, engine, NZ, NG = _setup()
    res = engine.certify(SynthOracle(), theta, prec=128)
    assert res["verdict"] == "CERTIFIED-INTERIOR", res["verdict"]
    ex = zs_ex[0] + zs_ex[1] + g_ex
    rmin = Fraction(res["radius_min"]).limit_denominator(10 ** 18)
    assert all(abs(Fraction(theta[i]) - ex[i]) < rmin
               for i in range(len(ex))), "exact optimum not contained"
    # corner: bound g_0 above its unconstrained optimum; solve the
    # reduced system by masked Newton; declared bound -> CERTIFIED-CORNER
    from clinch.mask import MaskedOracle
    from flint import arb
    b0 = float(g_ex[0]) + 0.25
    base = SynthOracle()
    mo = MaskedOracle(base, {5: b0})
    th = np.array(theta); th[5] = b0
    for _ in range(60):
        ctr = ([[arb(float(th[t])) for t in idx] for idx in mo.part.blocks],
               [arb(float(th[t])) for t in mo.part.border_idx])
        Fz, Fg = mo.F(ctr)
        D, B, G = mo.H(ctr)
        idxs = ([t for idx in mo.part.blocks for t in idx]
                + list(mo.part.border_idx))
        n = len(idxs)

        def midnp(M):
            Mm = M.mid()
            return np.array([[float(Mm[i, j]) for j in range(M.ncols())]
                             for i in range(M.nrows())])
        o0, o1 = len(mo.part.blocks[0]), len(mo.part.blocks[1])
        H = np.zeros((n, n))
        Fv = np.array([float(v.mid()) for blk in Fz for v in blk]
                      + [float(v.mid()) for v in Fg])
        H[:o0, :o0] = midnp(D[0]); H[o0:o0 + o1, o0:o0 + o1] = midnp(D[1])
        H[:o0, o0 + o1:] = midnp(B[0]); H[o0 + o1:, :o0] = midnp(B[0]).T
        H[o0:o0 + o1, o0 + o1:] = midnp(B[1])
        H[o0 + o1:, o0:o0 + o1] = midnp(B[1]).T
        H[o0 + o1:, o0 + o1:] = midnp(G)
        step = np.linalg.solve(H, Fv)
        for i, t in enumerate(idxs):
            th[t] -= step[i]
        if np.max(np.abs(step)) < 1e-15:
            break
    lo = np.full(len(th), -np.inf); lo[5] = b0
    rc = engine.certify(base, list(th), bounds=(lo, None), prec=128)
    assert rc["verdict"] == "CERTIFIED-CORNER", rc["verdict"]
    assert rc["kkt"]["ok"] and rc["kkt"]["margins"][5]["grad_lower"] > 0
    # wrong-reason control: no bound declared -> the fat-ball law refuses
    rw = engine.certify(base, list(th), prec=128)
    assert rw["verdict"] == "REFUSED" and rw["failed"] == "fat-ball", rw
    # out-of-scope: unported space refuses TYPED via certify_guarded
    from clinch.oracle_v31 import OracleRefusal

    def make_bad():
        raise OracleRefusal("check_space", "not the v3.1 space (fixture)")
    ro = engine.certify_guarded(make_bad, theta)
    assert ro["verdict"] == "OUT-OF-SCOPE", ro
    return ("interior certified + exact-optimum containment; corner "
            "certified with strict KKT; undeclared-bound fat-ball refusal; "
            "out-of-scope typed")


@leg("L2 planted saddle -> PD refusal naming the direction")
def l2():
    np, clinch, SynthOracle, _, theta, engine, NZ, NG = _setup()
    from flint import arb
    osad = SynthOracle(saddle_block=1)
    th = np.array(theta)
    part = osad.part
    for _ in range(80):
        ctr = ([[arb(float(th[t])) for t in idx] for idx in part.blocks],
               [arb(float(th[t])) for t in part.border_idx])
        Fz, Fg = osad.F(ctr)
        D, B, G = osad.H(ctr)

        def midnp(M):
            Mm = M.mid()
            return np.array([[float(Mm[i, j]) for j in range(M.ncols())]
                             for i in range(M.nrows())])
        n = sum(NZ) + NG
        H = np.zeros((n, n))
        Fv = np.array([float(v.mid()) for blk in Fz for v in blk]
                      + [float(v.mid()) for v in Fg])
        H[:3, :3] = midnp(D[0]); H[3:5, 3:5] = midnp(D[1])
        H[:3, 5:] = midnp(B[0]); H[5:, :3] = midnp(B[0]).T
        H[3:5, 5:] = midnp(B[1]); H[5:, 3:5] = midnp(B[1]).T
        H[5:, 5:] = midnp(G)
        step = np.linalg.solve(H, Fv)
        for i in range(n):
            th[i] -= step[i]
        if np.max(np.abs(step)) < 1e-14:
            break
    res = engine.certify(osad, list(th), prec=128)
    assert res["verdict"] == "REFUSED", res["verdict"]
    assert res["krawczyk"]["verdict"] == "CERTIFIED", \
        "stationary point should still certify"
    assert "block_01" in str(res["failed"]) or \
        "midpoint_cholesky" in str(res["failed"]), res["failed"]
    assert "diagnostic_negative_direction" in (res.get("pd") or {}), \
        "negative direction unnamed"
    return ("saddle: existence+uniqueness certified, PD REFUSED naming "
            "the flipped block + float-diagnostic direction")


@leg("L3 flat-valley control (v3.6 unpinned vs v3.7 pinned)")
def l3():
    """Receipt-consuming leg (producers: battery/l3_run.py -> L3_RESULTS.json,
    battery/l3_gates.py -> L3_GATES.json, expected under
    $CLINCH_REFERENCE_DIR/clinch_l3/). The centered-space oracle is
    clinch.oracle_v36 — the auxiliary-border KKT extension (m = ubar +
    multiplier per centered channel; the v3.7 pin as a KKT border row).
    The stall optimum is the reference fit's stage-1 checkpoint
    (CKPT_STAGE1.npz, gmax 0.026999); the reference fit's final theta
    (THETA_STAGE1.npz) is a later tail-polished point and is receipted as
    an addendum."""
    bdir = os.path.join(RCPT, "clinch_l3")
    d = json.load(open(os.path.join(bdir, "L3_RESULTS.json")))
    g = json.load(open(os.path.join(bdir, "L3_GATES.json")))
    # v36 oracle identity gates: float64 register at reference-noise
    # level (the v3.1-established class) + FD gate incl. aux columns
    assert g["fit16"]["max_rel"] <= 5e-9, g["fit16"]["max_rel"]
    assert g["fit18a"]["max_rel"] <= 5e-9, g["fit18a"]["max_rel"]
    assert max(g["fd"].values()) <= 1e-6, g["fd"]
    # (a) the v3.6 UNPINNED space at its fitted (stall) optimum:
    #     REFUSED, pointing at the advance-level direction. MEASURED
    #     register: the Newton correction at the stall is ~1 curvature-
    #     sigma and concentrates on the ADVANCE-LEVEL machinery — the
    #     a_adv intercepts and the sp_adv pooling scale log_sig[3] (the
    #     level-redistribution pair V37_MECHANISM names) + advance-
    #     channel slopes/latents.
    st = d["fit16_stall_analysis"]
    assert d["fit16_stall_at_coverage"]["verdict"] == "REFUSED"

    def _adv_family(nm):
        return (nm.split("[")[0] in ("a_adv", "b_gro", "g_gro", "w_gro",
                                     "z_sp_adv", "z_gen_adv", "z_fam_adv")
                or nm in ("log_sig[3]", "log_sig[4]", "log_sig[5]",
                          "aux[1]"))
    assert all(_adv_family(nm) for nm, _ in st["top_scaled"]), \
        st["top_scaled"]
    assert st["scaled_inf"] >= 0.5, st   # the walked valley: ~1 sigma out
    # (b) the SAME model under the v3.7 pin: the contraction defect
    #     shrinks by ORDERS OF MAGNITUDE (common-radius register)
    comp = d["comparison"]
    assert comp["defect_ratio_common_radius"] >= 1e2, comp
    # (c) the pinned KKT system (deep-polished center) carries a full
    #     Krawczyk certificate (existence+uniqueness of the KKT point)
    assert d["fit18a_cert"]["verdict"] == "CERTIFIED" \
        and d["fit18a_cert"]["max_margin"] < 1.0, d["fit18a_cert"]
    return (f"stall optimum REFUSED at {st['scaled_inf']:.2f} sigma, top "
            f"step coords all advance-channel "
            f"({', '.join(nm for nm, _ in st['top_scaled'])}); defect "
            f"ratio {comp['defect_ratio_common_radius']:.3g} at common "
            f"radius; fit18a pinned-KKT Krawczyk CERTIFIED margin "
            f"{d['fit18a_cert']['max_margin']:.3g} at ru="
            f"{d['fit18a_cert']['r_unit']:.2g}; findings receipted: "
            f"final reference theta already tail-polished ("
            f"its own refusal 1.2e5 = under-polish class), fit18a active "
            f"q-rail clip census; PD-for-centered-spaces = the "
            f"reduced-Hessian leg (clinch.kkt_pd, proven by L8)")


class HonestRemainder(Exception):
    pass


@leg("L4 oracle-tamper / mutation controls")
def l4():
    np, clinch, SynthOracle, _, theta, engine, NZ, NG = _setup()
    from baller.certify import block_krawczyk as bk
    # (a) math sign flip at the pristine certificate's own radius
    res = engine.certify(SynthOracle(), theta, prec=128)
    part = SynthOracle.part
    center = ([[theta[t] for t in part.blocks[0]],
               [theta[t] for t in part.blocks[1]]],
              [theta[t] for t in part.border_idx])
    mut = bk.block_krawczyk(SynthOracle(flip_math=True), center,
                            res["radius"], 128)
    assert mut["verdict"] == "REFUSED", "math sign-flip escaped"
    # (b) bit-flipped scratch copy refuses typed
    clone = os.path.join(os.getcwd(), "clinch_clone_tamper")
    if os.path.exists(clone):
        shutil.rmtree(clone)
    shutil.copytree(TOOL, clone)
    p = os.path.join(clone, "clinch", "oracle_v31.py")
    b = bytearray(open(p, "rb").read()); b[len(b) // 2] ^= 1
    open(p, "wb").write(bytes(b))
    r = subprocess.run([sys.executable, "-c", (
        f"import sys; sys.path.insert(0, {clone!r})\n"
        "import clinch\n"
        "try:\n"
        "    clinch.verify(); print('NOTAMPER')\n"
        "except clinch.ClinchTamperError: print('REFUSED-TYPED')\n")],
        capture_output=True, text=True, timeout=120)
    assert "REFUSED-TYPED" in r.stdout, r.stdout + r.stderr
    # (c) poisoned pyc has NO effect (purged at import)
    clone2 = os.path.join(os.getcwd(), "clinch_clone_pyc")
    if os.path.exists(clone2):
        shutil.rmtree(clone2)
    shutil.copytree(TOOL, clone2)
    plant = os.path.join(clone2, "clinch", "__pycache__")
    os.makedirs(plant, exist_ok=True)
    import importlib._bootstrap_external as be
    srcp = os.path.join(clone2, "clinch", "jets.py")
    poisoned = open(srcp).read() + "\nPOISONED = True\n"
    st = os.stat(srcp)
    pyc = be._code_to_timestamp_pyc(compile(poisoned, srcp, "exec"),
                                    st.st_mtime, st.st_size)
    tag = sys.implementation.cache_tag
    open(os.path.join(plant, f"jets.{tag}.pyc"), "wb").write(pyc)
    r2 = subprocess.run([sys.executable, "-c", (
        f"import sys; sys.path.insert(0, {clone2!r})\n"
        "import clinch\n"
        "from clinch import jets\n"
        "print('POISONED' if getattr(jets, 'POISONED', False) else 'CLEAN')\n")],
        capture_output=True, text=True, timeout=120)
    assert "CLEAN" in r2.stdout, r2.stdout + r2.stderr
    # (d) REAL-oracle input tamper vs the receipted identity-gate value
    gate_p = os.path.join(RCPT, "clinch_fit13_null", "CERT.json")
    cert = json.load(open(gate_p))
    nll_pinned = cert["oracle_identity_gates"]["gradient"]["registered"][
        "nll_arb_mid"]
    from clinch import adapter_v31 as A
    from clinch.oracle_v31 import ModelV31Oracle
    from baller.hygiene import ctx_guard
    from flint import arb
    arrs, consts, th13, prov = A.load("fit13")
    arrs_bad = dict(arrs)
    arrs_bad["X"] = -np.asarray(arrs["X"])          # planted sign flip
    o_bad = ModelV31Oracle(arrs_bad, consts)
    with ctx_guard(prec=192):
        thv = [arb(float(t)) for t in th13]
        nll_bad, _, _ = o_bad.evaluate(thv, order2=False)
        dev = abs(float(nll_bad.mid()) - nll_pinned)
    assert dev > 1.0, f"sign-flipped oracle within gate bar (dev={dev})"
    return ("math-flip refused at fixed radius; scratch-copy tamper refused "
            "typed; poisoned pyc inert; real-oracle sign flip moves nll "
            f"{dev:.3g} vs receipted gate value — identity gate catches it")


@leg("L5 precision honesty: dps ladder agrees; starved prec refuses")
def l5():
    np, clinch, SynthOracle, _, theta, engine, NZ, NG = _setup()
    outs = []
    for dps in (30, 50, 80):
        r = engine.certify(SynthOracle(), theta,
                           prec=int(dps / 0.30103) + 1)
        assert r["verdict"] == "CERTIFIED-INTERIOR", (dps, r["verdict"])
        outs.append(r["krawczyk"]["max_margin"])
    assert max(outs) <= 10 * max(min(outs), 1e-30) or max(outs) < 1e-3, \
        f"dps ladder margins disagree: {outs}"
    rs = engine.certify(SynthOracle(), theta, prec=8)
    assert rs["verdict"] == "REFUSED", "starved precision emitted certificate"
    return f"dps 30/50/80 certified (margins {['%.2g' % m for m in outs]}); prec=8 refused"


@leg("L6 flagship receipts: fit13_null + fit17_null")
def l6():
    out = []
    for tag in ("fit13", "fit17"):
        p = os.path.join(RCPT, f"clinch_{tag}_null", "CERT.json")
        assert os.path.exists(p), f"{p} missing — run battery/flagship.py"
        c = json.load(open(p))
        # direct certificate: an honest NAMED refusal (the polish finding)
        assert c["verdict"] == "REFUSED" and c["failed"], \
            f"{tag}: direct verdict {c['verdict']}"
        assert c["ball_center_sha16"] and c["statement"]
        g = c["oracle_identity_gates"]["gradient"]
        assert g["adjudication_ok"] and g["adjudication_bar"] == 1e-12, g
        assert c["oracle_identity_gates"]["fdhess"]["ok"]
        # polished-center certificate: CERTIFIED-INTERIOR, margins in bar
        pc = c["polished_center_certificate"]
        assert pc["verdict"] == "CERTIFIED-INTERIOR", pc["verdict"]
        assert pc["max_margin"] < 1.0 and pc["radius"] > 0
        assert min(pc["pd_margins"].values()) > 0
        assert pc["distance_inf_from_receipted"] < 1e-3
        out.append(f"{tag}: direct REFUSED({c['failed']}), polished "
                   f"CERTIFIED-INTERIOR margin {pc['max_margin']:.3g} "
                   f"dist {pc['distance_inf_from_receipted']:.3g}")
    return "; ".join(out)


@leg("L7 engine exhaustion path: all rungs raise -> typed refusal "
     "(etienne F6)")
def l7():
    """The etienne-deployment F6 finding: when EVERY radius rung raises
    an oracle refusal, the engine must return REFUSED(all-rungs-failed)
    with the per-rung failure classes listed — never crash on a
    best=None dereference."""
    np, clinch, SynthOracle, _, theta, engine, NZ, NG = _setup()
    from clinch.oracle_v31 import OracleRefusal

    class RailOracle(SynthOracle):
        """Works at the tight center; raises on every widened box."""

        def H(self, x):
            zs, g = x
            if any(float(v.rad()) > 0.0 for blk in zs for v in blk) or \
                    any(float(v.rad()) > 0.0 for v in g):
                raise OracleRefusal("test.rail",
                                    "planted per-rung refusal")
            return super().H(x)

    res = engine.certify(RailOracle(), theta, prec=128)
    assert res["verdict"] == "REFUSED", res["verdict"]
    assert res["failed"] == "all-rungs-failed", res["failed"]
    classes = [r["krawczyk"].get("failed") for r in res["rungs"]]
    assert classes and all(c == "oracle:test.rail" for c in classes), \
        classes
    assert "oracle:test.rail" in res["reason"], res["reason"]
    return (f"{len(classes)} rungs all raised; typed refusal "
            f"'all-rungs-failed' with per-rung classes listed — no crash")


@leg("L8 centered space: reduced-Hessian PD (rank-4 congruence "
     "extension)")
def l8():
    """The extended-KKT PD statement (oracle_v36's SEMANTICS NOTE): the
    Krawczyk leg certifies the KKT point; clinch.kkt_pd certifies the
    REDUCED Hessian PD by the signed-congruence inertia certificate. The
    toy centered fixture runs the same code path at rank 2 that the
    centered model spaces ride at rank 4/5."""
    from fractions import Fraction
    np, clinch, SynthOracle, _, _, engine, NZ, NG = _setup()
    from synthetic import CenteredSynthOracle, exact_centered_optimum
    # (a) planted centered optimum: FULL certificate, both legs
    zs, g, m, mu = exact_centered_optimum()
    ex = zs[0] + zs[1] + g + [m, mu]
    n_ext = len(ex)
    th = [float(v) for v in ex]
    res = engine.certify(CenteredSynthOracle(), th, prec=128)
    assert res["verdict"] == "CERTIFIED-INTERIOR", res["verdict"]
    pd = res["pd"]
    assert pd["verdict"] == "PD_REDUCED_CERTIFIED", pd
    assert pd["inertia"] == [n_ext - 1, 1, 0] and pd["n_mult"] == 1 \
        and pd["reduced_dim"] == n_ext - 2, pd
    assert min(pd["margins"].values()) > 0, pd["margins"]
    assert "constrained" in res["reason"], res["reason"]
    rmin = Fraction(res["radius_min"]).limit_denominator(10 ** 18)
    assert all(abs(Fraction(th[i]) - ex[i]) < rmin
               for i in range(n_ext)), "exact KKT point not contained"
    # (b) control: indefiniteness in a BORDER direction of the reduced
    #     space -> the typed kkt-signature refusal (the inertia COUNT,
    #     the code path L2's Cholesky control never reaches); the KKT
    #     point itself stays unique, so Krawczyk still certifies
    zs2, g2, m2, mu2 = exact_centered_optimum(gflip=-6)
    th2 = [float(v) for v in zs2[0] + zs2[1] + g2 + [m2, mu2]]
    r2 = engine.certify(CenteredSynthOracle(gflip=-6.0), th2, prec=128)
    assert r2["verdict"] == "REFUSED" and r2["failed"] == \
        "kkt-signature", (r2["verdict"], r2["failed"])
    assert r2["krawczyk"]["verdict"] == "CERTIFIED", \
        "KKT point should still certify"
    assert r2["pd"]["negative_count"] == 2 \
        and r2["pd"]["expected_negative_count"] == 1, r2["pd"]
    # (c) control: latent-block indefiniteness -> midpoint-Cholesky
    #     refusal naming the block (same register as L2, new leg)
    zs3, g3, m3, mu3 = exact_centered_optimum(saddle_block=1)
    th3 = [float(v) for v in zs3[0] + zs3[1] + g3 + [m3, mu3]]
    r3 = engine.certify(CenteredSynthOracle(saddle_block=1), th3,
                        prec=128)
    assert r3["verdict"] == "REFUSED" and str(r3["failed"]).startswith(
        "midpoint_cholesky:block_01"), (r3["verdict"], r3["failed"])
    assert "diagnostic_negative_direction" in (r3.get("pd") or {}), \
        "negative direction unnamed"
    return (f"centered KKT fixture CERTIFIED-INTERIOR with "
            f"PD_REDUCED_CERTIFIED: inertia {tuple(pd['inertia'])}, "
            f"reduced dim {pd['reduced_dim']}, min margin "
            f"{min(pd['margins'].values()):.3g}, exact KKT point "
            f"contained; flipped-border control -> kkt-signature "
            f"(2 negatives vs 1 allowed); flipped-latent control -> "
            f"midpoint-Cholesky refusal naming block_01")


def main():
    t0 = time.time()
    results, fails, remainders = [], [], []
    for name, fn in RESULTS:
        t = time.time()
        try:
            detail = fn()
            ok, status = True, "PASS"
        except HonestRemainder as e:
            detail = str(e)
            ok, status = True, "REMAINDER"
            remainders.append(name)
        except Exception as e:
            detail = f"{type(e).__name__}: {e}"
            ok, status = False, "FAIL"
            fails.append(name)
        wall = time.time() - t
        print(f"{status}  {name}  [{wall:.1f}s]  {str(detail)[:220]}")
        results.append({"leg": name, "status": status,
                        "wall_s": round(wall, 2),
                        "detail": str(detail)[:500]})
    total = time.time() - t0
    with open("CLINCH_BATTERY_SUMMARY.json", "w") as fh:
        json.dump({"results": results, "total_s": round(total, 1),
                   "remainders": remainders,
                   "overall": "FAIL" if fails else
                   ("PASS-WITH-REMAINDER" if remainders else "PASS")}, fh,
                  indent=1)
    if fails:
        print(f"OVERALL FAIL: {fails} ({total:.0f}s)")
        return 1
    if remainders:
        print(f"OVERALL PASS-WITH-REMAINDER ({len(results)} legs, "
              f"remainder: {remainders}) ({total:.0f}s) — STAGED until "
              f"the remainder lands")
        return 0
    print(f"OVERALL PASS ({len(results)} legs, {total:.0f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

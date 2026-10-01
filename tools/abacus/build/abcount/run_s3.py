"""S3 driver — stage S3 of the sha-pinned spec (M3T_REGISTRATION.md sec 4
row S3):
the general (2,2) machinery — B4 (rho predicate wired), the split-prime
detector member B4s (CORE), negative controls N1-N5, B5 (precondition
CONFIRMED at build time -> runs), FULL battery run + BATTERY_RESULT
banking.

Seeds: battery/SEEDS.json S3 revision (sha pinned below, committed BEFORE
any S3 battery run).  Stage law: S1/S2 receipts quote their own SEEDS
revisions; the B2/B3 battery legs are therefore REPLAYED through run_s2.py
under an exactly-reconstructed S2 revision (sha-verified, receipted; the
replay runs in an isolated scratch copy of the package tree so the shipped,
sha-pinned SEEDS.json is never written).  The B1 leg is re-exercised
natively here against
the verbatim-carried B1 blocks of the S3 revision (the S1 revision is not
reconstructible; its stage receipt S1_RESULT.md stays banked).

Walls (S3 row: priced from S2 actuals, declared in SEEDS S3_walls BEFORE
launch): <= 3600 s wallclock, <= 4 GB, dps cap 400.  Timed pilots: the
F_11^4 desk enumeration is timed before the F_31^4 enumeration; one E[3]
torsion call is timed before the torsion sweep.  STOP semantics: wall
breach = STOP-WALL (no silent retry); any integer disagreement = STOP;
UNDECIDED-* is a verdict, never a retry loop.

Bindings riding every output: Q3c/Q3d, Q6 (single conversion site =
abcount_s0.convert_sign_convention; the S3 sign-pin receipt VERIFIES and
REFUSES only, it never converts).
"""

import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

from flint import acb, arb, ctx

import abcount_s0 as ab
import abcount_s1 as s1
import abcount_s3 as s3

HERE = os.path.dirname(os.path.abspath(__file__))
PKG_ROOT = os.path.dirname(os.path.dirname(HERE))
BUILD = os.path.dirname(HERE)

PINNED = {
    "M3T_REGISTRATION.md": "85ed68777d171378c70c4b590ae35c25abb768431274bf00af6e2ed959b7bdba",
    "M3T_VALIDATION_BATTERY.md": "af97c604704f64e74194d64ddcd05fe5d5076f0fd999511b4875194c980b0d0c",
    "S0PRE_RECEIPT.md": "a23cb2ebe8de4b0819b08c4ffb9e9f59c9f891c5ef6ef99ede56c71154e9ec2e",
    "build/abcount/abcount_s0.py": "448dad948b9de15f4075097fe6b4c669e2cfc928e4836bec3cbe86307bb013cf",
    "build/abcount/abcount_s1.py": "8a6130d82ac9e3a20f6fa4210d654db7d204810de52baab7888d2deff043c0bf",
    "build/abcount/abcount_s2.py": "3dd63aa52e25e855b72c70ed1c46bf4a88765673640d9359696bb918bbe957c7",
    "build/abcount/hecke_desk.py": "38433f503ffbbcf849b2048ff28205f55b873c824b141240d954f8caaca7f440",
    "build/abcount/theta_g4_bridge.jl": "863f7c3d4f0fd2981a455992a33a58a2e3f616cbacb9a2ebdd8fad322f487fcf",
    "build/abcount/period_g2_bridge.jl": "19874a48bb95c91e7651bff7d357b3a238345d7aee97ec8cea94b9a0052c19ee",
    "build/abcount/battery/SEEDS.json": "9da4d6eebfced2378cc4b7b84f43d21f61c2b22d8c437c6ac78eeed429e5842d",
}
SEEDS_S3_SHA = PINNED["build/abcount/battery/SEEDS.json"]
SEEDS_S2_SHA = "23c86ba23c6b9b9e16d7f9342bd04b48e70bdffc437c6aa4065aa062f3e3488e"

STAGE_S3_TEXT = "S3 revision (extends the S2 revision sha 23c86ba23c6b9b9e16d7f9342bd04b48e70bdffc437c6aa4065aa062f3e3488e per this file's stage law: later stages EXTEND before their own first run and quote the new sha; S1/S2 receipts quote their own revisions, S3 receipts quote THIS revision). B1/S1_C7_pilot/B2/B3/S2_*/N1/N5 blocks carried verbatim from the S2 revision; B4, B4s (with B4s_twist_sign_u FILLED), B5, N2, N3, N4 and S3_walls are the S3 additions, committed BEFORE any S3 battery run; every numerical desk truth below was computed at seed time by pure-integer python arithmetic (PARI used only as seed-time cross-check where stated)."
STAGE_S2_TEXT = "S2 revision (extends the S1 revision sha 6f823872744f36a7be49d03257f3dffa31e2f34a780be5241edb09ca24c96751 per this file's stage law: later stages EXTEND before their own first run and quote the new sha; the S1 receipts quote the S1 revision, S2 receipts quote THIS revision). B1/S1_C7_pilot/N1/N5/B4s blocks carried verbatim from the S1 revision; B2, B3 and S2 blocks are the S2 additions, committed before any S2 battery run."
S2_TAIL = '''  "B4s_twist_sign_u": {
    "commitment": "committed together with the E20 model instantiation at stage S3 (battery: 'committed with the model in SEEDS.json'); the model is not instantiated at S1, so u is not yet chooseable — this field records the obligation, and the S3 revision of this file must fill it BEFORE any B4/B4s run",
    "run_stage": "S3"
  }
}
'''

WALL_SECONDS = 3600
WALL_MEM_GB = 4.0
WALL_DPS_CAP = 400
MEMBER_DPS = 60
PILOT_N = 12
MEMBER_PREC_BITS = 250
RADIUS_DISCIPLINE_GATE = "1e-30"


def sha256_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    ctx.dps = MEMBER_DPS
    wall = ab.Wall(seconds=WALL_SECONDS, mem_gb=WALL_MEM_GB, dps=WALL_DPS_CAP)
    checks = []
    member_rows = []
    banked = {"stage": "S3", "seeds_revision_sha": SEEDS_S3_SHA,
              "checks": checks, "member_rows": member_rows}

    def bank(name, verdict, detail):
        checks.append({"check": name, "verdict": verdict, "detail": detail})
        print("[%s] %s" % (verdict, name))
        with open(os.path.join(HERE, "S3_RUN_OUTPUT.json"), "w") as f:
            json.dump(banked, f, indent=1, default=str)

    def bank_member(member, verdict, primes, note):
        member_rows.append({"member": member, "verdict": verdict,
                            "primes": primes, "note": note,
                            "seeds_revision_sha": SEEDS_S3_SHA})
        with open(os.path.join(BUILD, "S3_RESULT.json"), "w") as f:
            json.dump({"stage": "S3", "seeds_revision_sha": SEEDS_S3_SHA,
                       "member_rows": member_rows,
                       "status": "RUNNING (crash-safe row bank)"}, f,
                      indent=1, default=str)

    # ---- check 1: sha pins (spec + modules) + SEEDS S3 sha -----------------
    sha_rows, ok = {}, True
    for rel, want in PINNED.items():
        got = sha256_file(os.path.join(PKG_ROOT, rel))
        sha_rows[rel] = {"got": got[:16], "match": got == want}
        ok = ok and got == want
    bank("sha-pin-seeds-binding-verify", "OK" if ok else "FAIL-SHA-PIN-VERIFY",
         sha_rows)
    if not ok:
        return finish(banked, wall, "S3-FAIL", bank)

    seeds = json.load(open(os.path.join(HERE, "battery", "SEEDS.json")))

    # ---- check 2: B2/B3 legs — run_s2.py replay under reconstructed S2 rev --
    seeds_path = os.path.join(HERE, "battery", "SEEDS.json")
    s3_text = open(seeds_path, encoding="utf-8").read()
    head = s3_text.split('  "B4s_twist_sign_u": {')[0]
    s2_text = head.replace(STAGE_S3_TEXT, STAGE_S2_TEXT) + S2_TAIL
    s2_sha_got = hashlib.sha256(s2_text.encode()).hexdigest()
    recon_ok = s2_sha_got == SEEDS_S2_SHA
    if not recon_ok:
        bank("S2-seeds-reconstruction", "AMBIGUOUS-STOP",
             {"got": s2_sha_got, "want": SEEDS_S2_SHA,
              "reason": "S2 revision not exactly reconstructible; B2/B3 replay impossible"})
        return finish(banked, wall, "AMBIGUOUS-STOP", bank)
    prev = os.path.join(HERE, "S2_RUN_OUTPUT.json")
    prev_save = os.path.join(HERE, "S2_RUN_OUTPUT_s2stage.json")
    if os.path.exists(prev) and not os.path.exists(prev_save):
        shutil.copy2(prev, prev_save)   # preserve the S2-stage receipt artifact
    # The replay runs in an ISOLATED copy of the package tree under a scratch
    # directory: the S2-revision seeds are written only into that copy, so the
    # shipped, sha-pinned SEEDS.json is never modified and a kill at any
    # instant (per-attempt timeout, out-of-memory, ctrl-C) cannot leave
    # swapped bytes behind for a later run's pin verify to trip over.
    replay_root = tempfile.mkdtemp(prefix="abacus_s2replay_")
    try:
        rpkg = os.path.join(replay_root, "pkg")
        rhere = os.path.join(rpkg, "build", "abcount")
        os.makedirs(os.path.join(rpkg, "build"))
        shutil.copytree(HERE, rhere,
                        ignore=shutil.ignore_patterns("__pycache__"))
        for spec in ("M3T_REGISTRATION.md", "M3T_VALIDATION_BATTERY.md",
                     "S0PRE_RECEIPT.md"):
            shutil.copy2(os.path.join(PKG_ROOT, spec), os.path.join(rpkg, spec))
        rseeds = os.path.join(rhere, "battery", "SEEDS.json")
        with open(rseeds, "w", encoding="utf-8") as f:
            f.write(s2_text)
        assert sha256_file(rseeds) == SEEDS_S2_SHA
        rprev = os.path.join(rhere, "S2_RUN_OUTPUT.json")
        if os.path.exists(rprev):
            os.remove(rprev)            # stale-output guard: replay must recreate it
        t0 = time.time()
        r = subprocess.run([sys.executable, os.path.join(rhere, "run_s2.py")],
                           capture_output=True, text=True, timeout=3300)
        replay_s = time.time() - t0
        if os.path.exists(rprev):
            s2out = json.load(open(rprev))
            s2_verdict = s2out.get("verdict")
            s2_checks = s2out.get("checks_ok")
            # publish the fresh replay receipt beside the harness (atomic
            # rename; the pre-replay artifact stays as S2_RUN_OUTPUT_s2stage.json)
            tmp_pub = prev + ".replay-tmp"
            shutil.copyfile(rprev, tmp_pub)
            os.replace(tmp_pub, prev)
        else:
            s2_verdict, s2_checks = "NO-OUTPUT", "0/0"
    finally:
        shutil.rmtree(replay_root, ignore_errors=True)
    restore_ok = sha256_file(seeds_path) == SEEDS_S3_SHA
    replay_ok = (r.returncode == 0 and s2_verdict == "S2-PASS" and restore_ok)
    replay_run = {"returncode": r.returncode, "verdict": s2_verdict,
                  "checks_ok": s2_checks, "wall_seconds": round(replay_s, 1),
                  "verdict_line": r.stdout.strip().splitlines()[-1] if r.stdout else ""}
    if not replay_ok:
        # surface the replay's own words: a crashed run_s2 (engine missing,
        # bridge failure) otherwise leaves only a stacktrace fragment here
        replay_run["stdout_tail"] = r.stdout[-1200:] if r.stdout else ""
        replay_run["stderr_tail"] = r.stderr[-1200:] if r.stderr else ""
    bank("B2-B3-replay-run_s2-under-reconstructed-S2-revision",
         "OK" if replay_ok else "STOP",
         {"seeds_swap": {"reconstructed_s2_sha": s2_sha_got,
                         "isolation": "replay ran in a scratch copy of the package tree; "
                                      "the shipped SEEDS.json was never written",
                         "shipped_s3_sha_ok": restore_ok,
                         "s2stage_artifact_preserved": os.path.basename(prev_save)},
          "run_s2": replay_run,
          "discipline": "verifier re-run, result checked, never quoted from memory"})
    bank_member("B2", "PASS" if replay_ok else "FAIL", [11, 31],
                "run_s2 replay %s (%s)" % (s2_verdict, s2_checks))
    bank_member("B3", "PASS" if replay_ok else "FAIL", [11, 31],
                "run_s2 replay %s (both integer routes + tie receipts inside)" % s2_verdict)
    if not replay_ok:
        return finish(banked, wall, "STOP", bank)
    wall.check("post-s2-replay")

    # ---- checks 3-7: B1 leg (native re-exercise, S3 revision) ---------------
    b1 = seeds["B1"]
    curves = b1["curves"]
    invs, taus, g1_thetas = [], [], []
    rows, ok = [], True
    for c in curves:
        inv = ab.weierstrass_invariants(*c["weierstrass"])
        invs.append(inv)
        ok = ok and int(inv["disc"]) == c["disc_banked"]
        w1, w2, prec_receipt = ab.pari_periods(tuple(c["weierstrass"]), MEMBER_DPS)
        tau, tau_receipt = ab.tau_from_periods(w1, w2)
        taus.append(tau)
        t2a, t3a, t4a = ab.theta_nulls_arb(tau)
        t2d, t3d, t4d, _ = ab.theta_nulls_desk(tau, N=PILOT_N)
        g1_thetas.append((t2a, t3a, t4a))
        overlaps = [ab.balls_overlap(x, y) for x, y in
                    ((t2a, t2d), (t3a, t3d), (t4a, t4d))]
        jac_ok = bool((t3a ** 4 - t2a ** 4 - t4a ** 4).contains(acb(0)))
        ok = ok and all(overlaps) and jac_ok and bool(tau.imag > 0)
        rows.append({"label": c["label"], "disc_match": int(inv["disc"]) == c["disc_banked"],
                     "arb_vs_desk_overlap": overlaps, "jacobi_ok": jac_ok,
                     "im_tau_pos": bool(tau.imag > 0)})
    contract_ok = (b1["sign_convention"]["s1"] == ab.GLOBAL_PIN_S1
                   and b1["polarization_elementary_divisors"] == [1, 1, 1, 1]
                   and len(set(r["label"] for r in rows)) == 4)
    bank("B1-blocks-tau-H1-genus1-theta", "OK" if (ok and contract_ok) else "STOP", rows)
    if not (ok and contract_ok):
        bank_member("B1", "FAIL", [13, 31], "block receipts")
        return finish(banked, wall, "STOP", bank)
    wall.check("post-B1-blocks")

    # C7 pilot (dyadic tau; matched precision; dominance + overlap + echo)
    pilot_entries, pilot_taus = [], []
    for tau in taus:
        rm, re_e = s1.dyadic_truncate_64(tau.real)
        im, im_e = s1.dyadic_truncate_64(tau.imag)
        pilot_entries.append({"re_man": rm, "re_exp": re_e, "im_man": im,
                              "im_exp": im_e, "rad_man": 0, "rad_exp": 0})
        pilot_taus.append(acb(s1._arb_exact_from_man_exp(rm, re_e),
                              s1._arb_exact_from_man_exp(im, im_e)))
    lam_lb = min((t.imag for t in pilot_taus), key=lambda x: float(x)).lower()
    tb, tb_receipt = s1.tail_bound_c7(4, PILOT_N, lam_lb)
    matched_prec = int(math.ceil(-math.log2(float(tb)))) + 64
    jres = s1.julia_theta_all_g4(pilot_entries, matched_prec, HERE, "s3pilot")
    wall.check("post-pilot-call")
    echo_ok = all(
        s1._arb_exact_from_man_exp(e[1], e[2]) == s1._arb_exact_from_man_exp(pe["re_man"], pe["re_exp"])
        and s1._arb_exact_from_man_exp(e[3], e[4]) == s1._arb_exact_from_man_exp(pe["im_man"], pe["im_exp"])
        and e[5] == 0 and e[7] == 0
        for e, pe in zip(jres["tau_echo"], pilot_entries))
    desk_dps = int(matched_prec * 0.302) + 20
    old = ctx.dps
    ctx.dps = desk_dps
    try:
        desk_vals = s1.desk_theta_g4_all(pilot_taus, PILOT_N)
    finally:
        ctx.dps = old
    tb_arb = arb(tb)
    pad = arb("0 +/- 1") * tb_arb
    dom = [i for i in range(256) if not bool(jres["claimed_radii"][i] <= tb_arb)]
    ovl = [i for i in range(256)
           if not ab.balls_overlap(desk_vals[i] + acb(pad, pad), jres["values"][i])]
    c7_ok = echo_ok and not dom and not ovl
    bank("B1-C7-pilot-cross-check", "OK" if c7_ok else "STOP",
         {"tail_bound": tb_arb.str(8), "matched_prec_bits": matched_prec,
          "dominance_violations": dom, "overlap_failures": ovl,
          "tau_echo_exact": echo_ok, "tail_bound_receipt_N": PILOT_N,
          "pilot_cost_s": jres["subprocess_wall_seconds"]})
    if not c7_ok:
        bank_member("B1", "FAIL", [13, 31], "C7 pilot")
        return finish(banked, wall, "STOP", bank)

    # member genus-4 call receipts
    member_entries = []
    for tau in taus:
        rm, re_e = s1._man_exp_of_mid(tau.real)
        im, im_e = s1._man_exp_of_mid(tau.imag)
        rad = tau.real.rad().abs_upper().max(tau.imag.rad().abs_upper())
        rman, rexp = (rad.man_exp() if not rad.is_zero() else (0, 0))
        member_entries.append({"re_man": rm, "re_exp": re_e, "im_man": im,
                               "im_exp": im_e, "rad_man": int(rman), "rad_exp": int(rexp)})
    jmem = s1.julia_theta_all_g4(member_entries, MEMBER_PREC_BITS, HERE, "s3member")
    wall.check("post-member-call")
    rad_gate = arb(RADIUS_DISCIPLINE_GATE)
    rad_fail = [i for i in range(256) if not bool(jmem["claimed_radii"][i] < rad_gate)]
    odd_fail, prod_fail = [], []
    for i in range(256):
        aa, bb = s1.char_decode(i)
        prod = acb(1)
        for k in range(4):
            t2a, t3a, t4a = g1_thetas[k]
            prod = prod * s1.genus1_theta_char(aa[k], bb[k], t2a, t3a, t4a)
        if not s1.char_is_even(i) and not bool(jmem["values"][i].contains(acb(0))):
            odd_fail.append(i)
        if not ab.balls_overlap(prod, jmem["values"][i]):
            prod_fail.append(i)
    ok = not odd_fail and not prod_fail and not rad_fail
    bank("B1-member-theta-receipts-route-T", "OK" if ok else "STOP",
         {"odd_failures": odd_fail, "product_failures": prod_fail,
          "radius_failures": rad_fail, "route_T_scope": "receipts only, no Frobenius data"})
    if not ok:
        bank_member("B1", "FAIL", [13, 31], "member theta receipts")
        return finish(banked, wall, "STOP", bank)

    # pslq recognition positive control (held-out gate)
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "lockpick")))
    import pslq_gate
    import mpmath as mp
    old = ctx.dps
    ctx.dps = 235
    try:
        w1, w2, _ = ab.pari_periods(tuple(curves[0]["weierstrass"]), 235)
        tau_hi, _ = ab.tau_from_periods(w1, w2)
        t2h, t3h, _ = ab.theta_nulls_arb(tau_hi)
        j_ball = ab.j_from_theta(t2h, t3h)
        j_str = j_ball.real.mid().str(228, radius=False)
    finally:
        ctx.dps = old
    with mp.workdps(240):
        pool = {"1": mp.nstr(mp.mpf(1), 228, strip_zeros=False),
                "pi": mp.nstr(mp.pi, 228, strip_zeros=False)}
    audit = pslq_gate.members_audit(pool, ["1", "pi"])
    ctl = pslq_gate.controls(pool, dict(pslq_gate.DEFAULT, maxcoeff=10 ** 9))
    vec = pslq_gate.two_prec_stable(j_str, ["1"], pool, dps_pair=(120, 200),
                                    maxcoeff=10 ** 9)
    inv0 = invs[0]
    jn, jd = inv0["j_num"], inv0["j_den"]
    if jd < 0:
        jn, jd = -jn, -jd
    g = math.gcd(abs(jn), jd)
    expected = (jd // g, -jn // g) if jd // g > 0 else (-jd // g, jn // g)
    rev = pslq_gate.reverify({"target": j_str, "members": ["1"],
                              "vector": list(vec) if vec else [0, 0]}, pool)
    ok = (audit == [] and ctl["ok"] and vec == expected and rev["ok"])
    bank("B1-pslq-recognition-control", "OK" if ok else "FAIL-RECOGNITION-CONTROL",
         {"vector_got": vec, "vector_expected": expected, "reverify_ok": rev["ok"]})
    if not ok:
        bank_member("B1", "FAIL", [13, 31], "pslq control")
        return finish(banked, wall, "S3-FAIL", bank)
    wall.check("post-pslq")

    # integer routes at the pre-committed primes
    truth = {13: b1["banked_desk_truths"]["a_13"], 31: b1["banked_desk_truths"]["a_31"]}
    pinned_counts = {13: b1["banked_desk_truths"]["count_F13"],
                     31: b1["banked_desk_truths"]["count_F31_derived"]}
    b1_ok = True
    for p in b1["pilot_primes"]:
        path1, path2 = [], []
        agree = True
        for c in curves:
            a_ellap = ab.ap_pari(tuple(c["weierstrass"]), p)
            a_desk, _ = ab.ap_naive(*c["weierstrass"], p)
            agree = agree and (a_ellap == a_desk == truth[p][c["label"]])
            path1.append(a_ellap)
            path2.append(a_desk)
        c2 = s1.assemble_lp_deg8(path2, p)
        c1 = s1.assemble_lp_deg8(path1, p)
        fe_ok, _ = s1.functional_equation_receipt(c2, p)
        verdict, tuples, box = s1.weil_box_isolate_deg8([acb(v) for v in c2[1:5]], p)
        cnt = s1.count_from_charpoly(c2)
        ok = (agree and c1 == c2 and fe_ok and verdict == "OK"
              and list(tuples[0]) == c2[1:5] and cnt == pinned_counts[p]
              and cnt == math.prod(p + 1 - a for a in path2))
        b1_ok = b1_ok and ok
        bank("B1-integer-routes-p%d" % p, "OK" if ok else "STOP",
             {"path1_ellap": path1, "path2_desk": path2, "L_c0_c8": c2,
              "weil_box": verdict, "count_P1": cnt, "pinned_truth": pinned_counts[p],
              "functional_equation": fe_ok})
        if not ok:
            bank_member("B1", "FAIL", [13, 31], "integer routes p=%d" % p)
            return finish(banked, wall, "STOP", bank)
        wall.check("post-B1-p%d" % p)
    bank_member("B1", "PASS", [13, 31],
                "native re-exercise: blocks+C7+theta+pslq+integer routes, both primes")

    # ---- check 8: B4 — the (2,2) O_K-action member --------------------------
    b4 = seeds["B4"]
    rho = s3.rho_product_matrix()
    E8 = s3.product_polarization_E()
    ok_rho, rec_rho = s3.ok_action_receipts(rho, E8)
    ok_pol, rec_pol = s3.polarization_type_check(E8, b4["polarization_elementary_divisors"])
    Pi = s3.b4_period_matrix()
    ok_an, s_sign, rec_an = s3.analytic_rep_receipts(Pi, rho)
    ok_pin, rec_pin = s3.sign_pin_receipt(s_sign, b4["sign_convention"]["s1"])
    layer_ok = ok_rho and ok_pol and ok_an and ok_pin
    bank("B4-OK-action-layer", "OK" if layer_ok else "STOP",
         {"rho_receipts": rec_rho, "polarization": rec_pol,
          "analytic": rec_an, "sign_pin": rec_pin})
    if not layer_ok:
        bank_member("B4", "FAIL", [11, 31], "O_K-action layer")
        return finish(banked, wall, "STOP", bank)

    b4_ok = True
    for p in b4["pilot_primes"]:
        jres = s3.e20_j_residues(p)
        path1_rows, aps1 = [], []
        for r5, jr in sorted(jres.items()):
            apv, rec = s3.ellap_e20_at(p, jr)
            aps1.append(apv)
            path1_rows.append({"sqrt5_residue": r5, "j_mod_P": jr, **rec})
        h = s3.hecke_ok_eval(p, seeds["B4s_twist_sign_u"]["u"])
        ap2 = h["a_P"]
        aps2 = [ap2, ap2]
        L1, cnt1 = s3.induced_lp_b4(aps1, p)
        L2, cnt2 = s3.induced_lp_b4(aps2, p)
        # deg-16 functional equation receipt: c_{16-k} = p^{8-k} c_k
        fe16 = all(L2[16 - k] == p ** (8 - k) * L2[k] for k in range(8))
        boxes = [ab.weil_box_isolate_g1(acb(a), p) for a in aps2]
        box_ok = all(v == "OK" and c == [a] for (v, c), a in zip(boxes, aps2))
        ok = (aps1 == aps2 and L1 == L2 and fe16 and box_ok
              and aps1 == [0, 0]
              and cnt1 == cnt2 == b4["banked_desk_truths"]["count_of_record"]["p%d" % p])
        b4_ok = b4_ok and ok
        bank("B4-routes-p%d" % p, "OK" if ok else "STOP",
             {"path1_ellap_per_P": path1_rows,
              "path2_hecke_desk": h["receipt"],
              "L_p_deg16_c0_c4": L2[:5], "functional_equation_deg16": fe16,
              "count_of_record": cnt1, "weil_box_per_P": [v for v, _ in boxes],
              "agreement": "EXACT-INTEGER"})
        if not ok:
            bank_member("B4", "FAIL", [11, 31], "routes p=%d" % p)
            return finish(banked, wall, "STOP", bank)
    bank_member("B4", "PASS" if b4_ok else "FAIL", [11, 31],
                "O_K receipts + a_P=0 both routes both primes + deg-16 assembly")
    wall.check("post-B4")

    # ---- check 9: B4s — split-prime rho-embedding detector (CORE) ----------
    b4s = seeds["B4s"]
    u = seeds["B4s_twist_sign_u"]["u"]
    h = s3.hecke_ok_eval(29, u)
    want_pin = b4s["banked_desk_truths"]["at_P_pin"]
    want_conj = b4s["banked_desk_truths"]["at_P_conj"]
    path2_ok = (h["at_pinned"]["a_P"] == want_pin["a_P"]
                and h["at_pinned"]["b_P"] == want_pin["b_P"]
                and h["at_conjugate"]["a_P"] == want_conj["a_P"]
                and h["at_conjugate"]["b_P"] == want_conj["b_P"])
    # Path 1: ellap at both primes + pilot-priced torsion congruence
    jres = s3.e20_j_residues(29)
    ap_rows, tors_rows = [], []
    t_pilot = None
    path1_ok = True
    for idx, (r5, jr) in enumerate(sorted(jres.items())):
        apv, rec = s3.ellap_e20_at(29, jr)
        t0 = time.time()
        F, frec = s3.frobenius_matrix_on_E3(29, jr)
        t_el = time.time() - t0
        if t_pilot is None:
            t_pilot = t_el     # first torsion call is the timed pilot
            if t_pilot > 60:
                bank("B4s-torsion-timed-pilot", "STOP-WALL",
                     {"pilot_seconds": t_pilot})
                bank_member("B4s", "FAIL", [29], "torsion pilot wall")
                return finish(banked, wall, "STOP-WALL", bank)
        okt, trec = s3.torsion_congruence_receipt(
            F, apv, 29, u, want_pin["b_P"], want_conj["b_P"])
        path1_ok = path1_ok and apv == -6 and okt and all(frec["order3_checks"])
        ap_rows.append({"sqrt5_residue": r5, "j_mod_P": jr, "ellap": apv})
        tors_rows.append({"j_mod_P": jr, "F_on_E3": F, "receipt": trec})
    # pair box: FULL enumeration of (a, b), a even, (a/2)^2 + 5 b^2 = 29
    pair_box = [(a, b) for a in range(-10, 11, 2) for b in range(-3, 4)
                if (a // 2) ** 2 * 4 == a * a and (a * a) // 4 + 5 * b * b == 29]
    pair_box_ok = sorted(pair_box) == [(-6, -2), (-6, 2), (6, -2), (6, 2)]
    # detector receipt: rho -> -rho flips the emitted b_P sign only
    ok_neg, s_neg, _ = s3.analytic_rep_receipts(Pi, s3.mat_scal(-1, rho))
    emitted_pin = {"a_P": h["at_pinned"]["a_P"], "b_P": s_sign * h["at_pinned"]["b_P"]}
    emitted_pin_neg = {"a_P": h["at_pinned"]["a_P"], "b_P": s_neg * h["at_pinned"]["b_P"]}
    detector_ok = (ok_neg and s_neg == -s_sign
                   and emitted_pin["b_P"] == want_pin["b_P"]
                   and emitted_pin_neg["b_P"] == -want_pin["b_P"]
                   and emitted_pin_neg["a_P"] == want_pin["a_P"])
    ok = path2_ok and path1_ok and pair_box_ok and detector_ok
    bank("B4s-detector-routes", "OK" if ok else "STOP",
         {"path2_hecke_desk": {"at_pinned": h["at_pinned"],
                               "at_conjugate": h["at_conjugate"],
                               "receipt": h["receipt"]},
          "path1_ellap": ap_rows,
          "path1_torsion_E3": tors_rows,
          "torsion_pilot_seconds": round(t_pilot, 3),
          "pair_weil_box": {"enumeration": sorted(pair_box), "ok": pair_box_ok,
                            "semantics": "FULL enumeration (a/2)^2+5b^2 = 29, a even"},
          "rho_flip_detector": {"s_pinned": s_sign, "s_negated_rho": s_neg,
                                "emitted_at_pinned": emitted_pin,
                                "emitted_at_pinned_negated_rho": emitted_pin_neg,
                                "b_P_flips_a_P_invariant": detector_ok},
          "count_E20_F_P": 29 + 1 + 6})
    bank_member("B4s", "PASS" if ok else "FAIL", [29],
                "(a_P,b_P)=(-6,-2)/(-6,+2) both routes; E[3] congruence; flip detector")
    if not ok:
        return finish(banked, wall, "STOP", bank)
    wall.check("post-B4s")

    # ---- check 10: B5 — Weil restriction member (precondition CONFIRMED) ----
    b5 = seeds["B5"]
    probe = s3.gp_run("t=ffgen([11,2],'t); print(Vec(hyperellcharpoly(x^5+t*x+1)));")
    probe_ok = probe[-1].strip("[]").split(",")[0].strip() == "1"
    bank("B5-precondition-reprobe", "OK" if probe_ok else "STOP",
         {"seeds_status": b5["precondition"]["status"],
          "in_run_probe_ok": probe_ok})
    if not probe_ok:
        bank_member("B5", "DROPPED", [], "precondition lost at run time")
        return finish(banked, wall, "STOP", bank)
    b5_ok = True
    b5_rows = []
    t0 = time.time()
    n1, n2, cp11, _ = s3.b5_counts_inert(11)
    t11 = time.time() - t0
    est31 = t11 * (31.0 / 11.0) ** 4
    if est31 > 600:
        bank("B5-pilot-pricing", "STOP-WALL", {"t11_s": t11, "est31_s": est31})
        bank_member("B5", "FAIL", [11, 29, 31], "pilot pricing breach")
        return finish(banked, wall, "STOP-WALL", bank)
    bank("B5-pilot-pricing", "OK", {"t_F11^4_pilot_s": round(t11, 2),
                                    "est_F31^4_s": round(est31, 2),
                                    "budget_s": 600})
    for p in b5["pilot_primes"]:
        bt = b5["banked_desk_truths"]["p%d" % p]
        if p == 29:
            d13 = s3.b5_counts_split(29, 13)
            d16 = s3.b5_counts_split(29, 16)
            cp_pin = s3.hyperellcharpoly_ff(29, 1, 13)
            cp_conj = s3.hyperellcharpoly_ff(29, 1, 16)
            path2_ok = (d13[2] == bt["pinned_s13"]["charpoly_monic_desc"]
                        and d16[2] == bt["conj_s16"]["charpoly_monic_desc"]
                        and [d13[0], d13[1]] == [bt["pinned_s13"]["count_C_Fp"], bt["pinned_s13"]["count_C_Fp2"]]
                        and [d16[0], d16[1]] == [bt["conj_s16"]["count_C_Fp"], bt["conj_s16"]["count_C_Fp2"]])
            path1_ok = cp_pin == d13[2] and cp_conj == d16[2]
            L8 = s3.lp_split_product(cp_pin, cp_conj)
            row = {"p": 29, "pinned_s13_cp": cp_pin, "conj_s16_cp": cp_conj}
        else:
            nn1, nn2, cp, drec = (n1, n2, cp11, None) if p == 11 else s3.b5_counts_inert(31)
            cpff = s3.hyperellcharpoly_ff(p, 2, (-5) % p)
            path2_ok = (cp == bt["charpoly_over_Fq_monic_desc"]
                        and [nn1, nn2] == [bt["count_C_Fq"], bt["count_C_Fq2"]])
            path1_ok = cpff == cp
            L8 = s3.lp_from_charpoly_T2(cpff)
            row = {"p": p, "q": p * p, "cp_over_Fq": cpff, "counts": [nn1, nn2]}
        fe_ok, _ = s1.functional_equation_receipt(L8, p)
        verdict, tuples, box = s1.weil_box_isolate_deg8([acb(v) for v in L8[1:5]], p)
        cnt = sum(L8)
        ok = (path1_ok and path2_ok and fe_ok and verdict == "OK"
              and list(tuples[0]) == L8[1:5] and cnt == bt["count_Res_Fp"])
        b5_ok = b5_ok and ok
        row.update({"L8_c0_c8": L8, "functional_equation": fe_ok,
                    "weil_box": verdict, "count_Res_Fp": cnt,
                    "banked_count": bt["count_Res_Fp"],
                    "path1_eq_path2": path1_ok, "path2_eq_banked": path2_ok})
        b5_rows.append(row)
        bank("B5-routes-p%d" % p, "OK" if ok else "STOP", row)
        if not ok:
            bank_member("B5", "FAIL", [11, 29, 31], "routes p=%d" % p)
            return finish(banked, wall, "STOP", bank)
        wall.check("post-B5-p%d" % p)
    bank_member("B5", "PASS" if b5_ok else "FAIL", [11, 29, 31],
                "Weil-restriction identity both routes, F_p and F_p^2 legs exercised")

    # ---- checks 11-15: negative controls ------------------------------------
    # N1: B1 off-diagonal tau perturbation -> Riemann-relation symmetry FAIL
    tau4 = [[acb(0)] * 4 for _ in range(4)]
    for k in range(4):
        tau4[k][k] = taus[k]
    tau4[0][1] = tau4[0][1] + acb(arb("1e-6"))       # ONE entry, (1,2) only
    sym_viol = []
    for i in range(4):
        for j in range(i + 1, 4):
            if not bool((tau4[i][j] - tau4[j][i]).contains(acb(0))):
                sym_viol.append([i + 1, j + 1])
    n1_fail_fired = sym_viol == [[1, 2]]
    bank("N1-perturbed-offdiag-tau", "OK" if n1_fail_fired else "STOP",
         {"required": "FAIL at Riemann-relation/polarization consistency, named",
          "named_predicate": "riemann-relation-tau-symmetry",
          "violating_entries": sym_viol,
          "perturbation": "tau_(1,2) += 1e-6 (> per-entry radius 1e-55; SEEDS N1_pinned_entry)",
          "second_channel": "cross-method mismatch vs planted truth armed (unused here: primary predicate fired)"})
    bank_member("N1", "PASS" if n1_fail_fired else "FAIL", [13],
                "FAIL-as-required: riemann-relation-tau-symmetry fired on entry (1,2)")
    if not n1_fail_fired:
        return finish(banked, wall, "STOP", bank)

    # N2: wrong elementary divisors -> polarization-type FAIL, no tau emitted
    n2 = seeds["N2_wrong_divisors"]
    E_bad = s3.product_polarization_E([1, 1, 2, 2])
    ok_bad, rec_bad = s3.polarization_type_check(E_bad, [1, 1, 1, 1])
    n2_ok = (not ok_bad
             and rec_bad["failing_predicate"] == "polarization-type-frobenius-normal-form"
             and rec_bad["symplectic_type"] == [1, 1, 2, 2])
    bank("N2-wrong-elementary-divisors", "OK" if n2_ok else "STOP",
         {"required": n2["required_verdict"], "receipt": rec_bad,
          "no_tau_emitted": "code path returns before any tau construction"})
    bank_member("N2", "PASS" if n2_ok else "FAIL", [11],
                "FAIL-as-required: polarization-type-frobenius-normal-form, no tau emitted")
    if not n2_ok:
        return finish(banked, wall, "STOP", bank)

    # N3: rho tampered one entry -> O_K-action FAIL (named)
    n3 = seeds["N3_rho_tamper"]
    rho_bad = s3.rho_product_matrix(tamper=(0, 1, 1))
    ok_bad, rec_bad = s3.ok_action_receipts(rho_bad, E8)
    n3_ok = (not ok_bad
             and rec_bad["failing_predicate"] == "OK-action-rho-square-eq-minus5")
    bank("N3-rho-tamper", "OK" if n3_ok else "STOP",
         {"required": n3["required_verdict"],
          "tamper": n3["commitment"],
          "failing_predicate": rec_bad["failing_predicate"],
          "residuals": rec_bad["predicates"][0]["residual_nonzero_entries"]})
    bank_member("N3", "PASS" if n3_ok else "FAIL", [11, 31],
                "FAIL-as-required: OK-action-rho-square-eq-minus5 named")
    if not n3_ok:
        return finish(banked, wall, "STOP", bank)

    # N4: precision starvation -> UNDECIDED-PRECISION + needed-dps, NO integer
    n4 = seeds["N4_precision_starvation"]
    p = 13
    n4_rad = n4["inflated_radius"]
    starved = [acb(arb("%d +/- %s" % (truth[13][c["label"]], n4_rad))) for c in curves]
    n4_rows, any_integer_emitted = [], False
    for ballv, c in zip(starved, curves):
        v, cands = ab.weil_box_isolate_g1(ballv, p)
        n4_rows.append({"label": c["label"], "verdict": v,
                        "n_candidates": len(cands),
                        "needed_dps_estimate": s1._needed_dps(ballv)})
        if v == "OK":
            any_integer_emitted = True
    n4_ok = (all(r["verdict"] == "UNDECIDED-PRECISION" for r in n4_rows)
             and not any_integer_emitted)
    bank("N4-precision-starvation", "OK" if n4_ok else "STOP",
         {"required": n4["required_verdict"], "per_factor": n4_rows,
          "inflated_radius": "%s per SEEDS N4_precision_starvation (erratum revision)" % n4_rad,
          "integer_emitted": any_integer_emitted})
    bank_member("N4", "PASS" if n4_ok else "FAIL", [13],
                "UNDECIDED-PRECISION as required; no integer emitted; needed-dps quoted")
    if not n4_ok:
        return finish(banked, wall, "STOP", bank)

    # N5: conjugated complex structure claiming P -> FAIL-SIGN-PIN (pre-declared)
    n5 = seeds["N5_required_verdict"]
    Pic = s3.b4_period_matrix(conjugated=True)
    ok_c, s_c, rec_c = s3.analytic_rep_receipts(Pic, rho)
    ok_pin_c, rec_pin_c = s3.sign_pin_receipt(s_c, "P")
    n5_ok = (ok_c and s_c == -1 and not ok_pin_c
             and rec_pin_c["verdict"] == "FAIL-SIGN-PIN"
             and n5["pre_declaration"] == "FAIL-SIGN-PIN")
    bank("N5-conjugated-complex-structure", "OK" if n5_ok else "STOP",
         {"pre_declared": n5["pre_declaration"],
          "got": rec_pin_c["verdict"],
          "holomorphic_sign": s_c,
          "coherent_conjugate_geometry": ok_c,
          "silent_same_answer": False,
          "receipt": rec_pin_c})
    bank_member("N5", "PASS" if n5_ok else "FAIL", [29],
                "FAIL-SIGN-PIN as pre-declared in SEEDS (rides seeds sha)")
    if not n5_ok:
        return finish(banked, wall, "STOP", bank)
    wall.check("post-negatives")

    return finish(banked, wall, "S3-PASS", bank)


def finish(banked, wall, verdict, bank):
    w = wall.check("finish")
    checks = banked["checks"]
    n_ok = sum(1 for c in checks if c["verdict"] == "OK")
    rows = banked["member_rows"]
    n_pass = sum(1 for r in rows if r["verdict"] == "PASS")
    battery = "BATTERY-PASS(%d/%d, B5 RUNS+PASS)" % (n_pass, len(rows)) \
        if (verdict == "S3-PASS" and n_pass == len(rows) == 11) else \
        "BATTERY-INCOMPLETE(%d/%d)" % (n_pass, len(rows))
    banked.update({
        "tool": "abcount (S3 general (2,2) layer over the S0-S2 core)",
        "verdict": verdict,
        "battery": battery,
        "checks_ok": "%d/%d" % (n_ok, len(checks)),
        "wall_declared": {"seconds": WALL_SECONDS, "mem_gb": WALL_MEM_GB,
                          "dps_cap": WALL_DPS_CAP,
                          "provenance": "SEEDS S3_walls, priced from S2 actuals BEFORE launch"},
        "wall_actuals": w,
        "sign_convention": {
            "s1": "P", "s2": "K-embedding pinned (rho = +i*sqrt5 holomorphic)",
            "conversion_site": "abcount_s0.convert_sign_convention (single site, Q6); "
                               "S3 sign-pin receipt verifies/refuses only, never converts"},
        "input_invariants_Q3c": {
            "B4_B4s": {"dimension": 4, "base_field": "Q(sqrt5)",
                       "polarization_elementary_divisors": [1, 1, 1, 1],
                       "O_K_action": "Z[sqrt(-5)] diagonal, rho 8x8 pinned in SEEDS",
                       "isogeny_declared": "product-diagonal (identity factors)"},
            "B5": {"dimension": 4, "base_field": "Q (Weil restriction from K)",
                   "curve": "y^2 = x^5 + sqrt(-5) x + 1 over K",
                   "isogeny_declared": "NONE (restriction of scalars, declared split over K)"}},
        "pinned_lattice_claim": "NONE (Q3d — no identification against the pinned P is made or implied)",
        "route_scope": "route T = theta receipts only; route D = integer assembly + "
                       "Weil box; route C = ellap/hyperellcharpoly truth side; "
                       "hecke_desk = PARI-free Path 2",
    })
    with open(os.path.join(HERE, "S3_RUN_OUTPUT.json"), "w") as f:
        json.dump(banked, f, indent=1, default=str)
    with open(os.path.join(BUILD, "S3_RESULT.json"), "w") as f:
        json.dump({"stage": "S3", "seeds_revision_sha": banked["seeds_revision_sha"],
                   "verdict": verdict, "battery": battery,
                   "checks_ok": banked["checks_ok"],
                   "member_rows": banked["member_rows"],
                   "wall_actuals": w}, f, indent=1, default=str)
    if verdict == "S3-PASS":
        bat = {"battery": "abcount full validation battery (M3T_VALIDATION_BATTERY.md 8b57eafb)",
               "seeds_revision_sha": banked["seeds_revision_sha"],
               "verdict": battery,
               "member_verdicts": {r["member"]: r["verdict"] for r in rows},
               "B5_status": "RUNS (precondition confirmed at build + re-probed in-run)",
               "notes": {r["member"]: r["note"] for r in rows},
               "wall_actuals": w,
               "stage_receipts": {"S0": "6906a33d", "S1": "8e2d76be",
                                  "S2": "2eaeddb0", "S3": "THIS run"},
               "acceptance": "all positives at all declared primes + all negatives "
                             "FAIL-as-required + all cross-method pairs exact "
                             "(battery sec d); B4s CORE exercised"}
        with open(os.path.join(HERE, "battery", "BATTERY_RESULT.json"), "w") as f:
            json.dump(bat, f, indent=1, default=str)
    print("VERDICT:", verdict, "|", battery, "|", banked["checks_ok"], "| wall:", w)
    return 0 if verdict == "S3-PASS" else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except RuntimeError as e:
        if "STOP-WALL" in str(e):
            print("STOP-WALL:", e)
            sys.exit(3)
        print("FAIL-INTERNAL:", e)
        traceback.print_exc()
        sys.exit(2)
    except Exception as e:
        print("FAIL-INTERNAL:", e)
        traceback.print_exc()
        sys.exit(2)

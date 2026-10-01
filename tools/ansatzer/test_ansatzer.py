#!/usr/bin/env python3
"""
Tests for ansatzer.py.

  T1. trivial 2-letter {x,y}: hand-checkable integrability (w2: 4→3, w3: 8→4).
  T2. P4 J23 (3-letter rational): validated integrable tower (3,7,15,31),
      first-entry residual (2,4,8,16).
  T3. P4 J25 (4-letter, parity-odd): validated w2 16→12→4→3 chain.
  T4. synthetic 7-letter alphabet (invented dlog-independent letters with a
      formal parity grading and channel tags): w1 cuts are hand-checkable
      (first=|FE|=5, last=|FE∩LE|=3); w2/w3 exact_Q dims pinned as regression
      values, with and without the Steinmann channel cut; strong collapse
      (verdict not FULL).
  T5. landau_alphabet.py-output adapter: load box1l alphabet.json directly.
  T6. Steinmann cut: hand-checkable 2-letter S/T-channel case (w2: 3→2,
      w3: 4→2), exact_Q vs mod-p agreement on the 7-letter alphabet, honesty
      flags (applied / vacuous / disabled).
  T7. per-orbit input adapter: recursive_Kw accepts a fixture path, a
      landau-alphabet-style dict, and a dict without entry lists (defaults to
      ALL letters), all reproducing the hand-checked trivial-2-letter K_w.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ansatzer import (load_alphabet, predict_collapse, format_table,
                                symbol_space_dims, integrability_cut,
                                recursive_Kw)

ALPH = os.path.join(HERE, "alphabets")
FAIL = []


def check(name, got, want, tol="=="):
    ok = (got == want) if tol == "==" else (got <= want)
    tag = "OK " if ok else "FAIL"
    print(f"  [{tag}] {name}: got {got}  want {tol}{want}")
    if not ok:
        FAIL.append(name)


def T1_trivial():
    print("\n[T1] trivial 2-letter {x,y}")
    a = load_alphabet(os.path.join(ALPH, "trivial_2var.json"))
    raw = symbol_space_dims(a, 3)
    check("raw w2", raw[2], 4); check("raw w3", raw[3], 8)
    _, _, d2, _ = integrability_cut(a, 2, method="exact_Q")
    _, _, d3, _ = integrability_cut(a, 3, method="exact_Q")
    # dlog x ∧ dlog y ≠ 0 ⇒ w2: kill (xy−yx) → 3;  w3: integrable = 4 (sym tensors)
    check("integrable w2", d2, 3)
    check("integrable w3", d3, 4)


def T2_P4_J23():
    print("\n[T2] P4 J23 (3-letter rational, FE={s,s+t})")
    rep = predict_collapse(os.path.join(ALPH, "p4_J23.json"),
                           weight_max=4, method="exact_Q")
    print(format_table(rep))
    # validated counts: integrable tower (1,3,7,15,31); first-entry residual (1,2,4,8,16)
    for w, (ig, fe) in {1: (3, 2), 2: (7, 4), 3: (15, 8), 4: (31, 16)}.items():
        check(f"w{w} integrable", rep["per_weight"][w]["integrability"], ig)
        check(f"w{w} first-entry", rep["per_weight"][w]["first"], fe)


def T3_P4_J25():
    print("\n[T3] P4 J25 (4-letter, parity-odd, FE={s,t,O_A3})")
    rep = predict_collapse(os.path.join(ALPH, "p4_J25.json"),
                           weight_max=4, method="exact_Q")
    print(format_table(rep))
    # J25 validated counts: w2 raw=16 integ=12 parity=4 first=3; w3 64→32→12→8; w4 256→80→32→20.
    # In this cut order (integ→first→last→parity) the FINAL residual should match.
    check("w2 integrable", rep["per_weight"][2]["integrability"], 12)
    check("w2 RESIDUAL",   rep["per_weight"][2]["RESIDUAL"], 3)
    check("w3 integrable", rep["per_weight"][3]["integrability"], 32)
    check("w3 RESIDUAL",   rep["per_weight"][3]["RESIDUAL"], 8)
    check("w4 integrable", rep["per_weight"][4]["integrability"], 80)
    check("w4 RESIDUAL",   rep["per_weight"][4]["RESIDUAL"], 20)


def T4_synth7():
    print("\n[T4] synthetic 7-letter alphabet (collapse expected)")
    rep = predict_collapse(os.path.join(ALPH, "synth_7letter.json"),
                           weight_max=3, method="exact_Q")
    print(format_table(rep))
    # w1 is hand-checkable: the 5 first-entry letters are multiplicatively
    # independent, so first = 5; FE∩LE = {u, v, omuv}, so last = 3; all
    # three are parity-even and w1 has no adjacent pair, so w1 RESIDUAL = 3.
    check("w1 first-entry = |FE| = 5", rep["per_weight"][1]["first"], 5)
    check("w1 RESIDUAL = |FE∩LE| = 3", rep["per_weight"][1]["RESIDUAL"], 3)
    # The fixture carries channel tags (u/omu/opu→A, v/omv/opv→B, omuv→C), so
    # the Steinmann cut is active.  w2 steinmann is hand-checkable at word
    # level: same-channel FE×LE words = {u,omu}×{u,opu} + {v,omv}×{v,opv} +
    # {omuv}×{omuv} = 4+4+1 = 9.  Other dims: exact-Q regression pins
    # (engine correctness is checked on the hand-derivable cases in T1/T6).
    check("w2 integrable", rep["per_weight"][2]["integrability"], 37)
    check("w2 steinmann", rep["per_weight"][2]["steinmann"], 9)
    check("w2 RESIDUAL", rep["per_weight"][2]["RESIDUAL"], 5)
    check("w3 integrable", rep["per_weight"][3]["integrability"], 175)
    check("w3 steinmann", rep["per_weight"][3]["steinmann"], 25)
    check("w3 RESIDUAL", rep["per_weight"][3]["RESIDUAL"], 13)
    check("verdict not FULL", "FULL" not in rep["verdict"], True)
    # Without the channel block the same alphabet reproduces the untagged
    # pins, and Steinmann residuals can only be ≤ the untagged ones.
    with open(os.path.join(ALPH, "synth_7letter.json")) as fh:
        d = json.load(fh)
    d.pop("channel")
    rep0 = predict_collapse(d, weight_max=3, method="exact_Q")
    check("w2 RESIDUAL (no channel)", rep0["per_weight"][2]["RESIDUAL"], 9)
    check("w3 RESIDUAL (no channel)", rep0["per_weight"][3]["RESIDUAL"], 31)
    for w in (2, 3):
        check(f"w{w} steinmann ≤ untagged residual",
              rep["per_weight"][w]["RESIDUAL"],
              rep0["per_weight"][w]["RESIDUAL"], tol="<=")


def T5_adapter():
    print("\n[T5] landau_alphabet.py output adapter (box1l)")
    # box1l = the one-loop box family (physics tag, not a machine). The fixture
    # is a landau_alphabet.py output; leg SKIPs when none is supplied.
    p = os.environ.get("COLLAPSE_T5_ALPHABET", "")
    if not p or not os.path.exists(p):
        print("  [SKIP] no alphabet.json supplied (set COLLAPSE_T5_ALPHABET)")
        return
    rep = predict_collapse(p, weight_max=3, method="exact_Q")
    check("box1l n_letters", rep["n_letters"], 2)
    check("box1l w2 integrable", rep["per_weight"][2]["integrability"], 3)


def T6_steinmann():
    print("\n[T6] Steinmann adjacent-channel cut")
    with open(os.path.join(ALPH, "trivial_2var.json")) as fh:
        d = json.load(fh)
    # Hand-checkable: tag x→S, y→T (overlapping channels).  The w2 integrable
    # space is span{xx, yy, xy+yx}; Steinmann kills the words xy and yx, so
    # only span{xx, yy} survives → 2.  At w3 every mixed word has an adjacent
    # x,y pair, so the 4-dim integrable space drops to span{xxx, yyy} → 2.
    d_st = dict(d); d_st["channel"] = {"x": ["S"], "y": ["T"]}
    rep = predict_collapse(d_st, weight_max=3, method="exact_Q")
    check("w2 integrable", rep["per_weight"][2]["integrability"], 3)
    check("w2 steinmann (kill xy,yx)", rep["per_weight"][2]["steinmann"], 2)
    check("w3 integrable", rep["per_weight"][3]["integrability"], 4)
    check("w3 steinmann (xxx,yyy only)", rep["per_weight"][3]["steinmann"], 2)
    check("honesty: applied",
          rep["honesty"]["steinmann"].startswith("applied"), True)
    # bare-string tags are accepted (one tag, not splintered per character)
    a = load_alphabet({**d_st, "channel": {"x": "CH1", "y": "CH1"}})
    check("bare-string tag parses whole", a.channel[0], {"CH1"})
    # untagged alphabet → cut vacuous; disabled → SKIPPED, untagged dims back
    rep0 = predict_collapse(d, weight_max=2, method="exact_Q")
    check("honesty: vacuous",
          rep0["honesty"]["steinmann"].startswith("vacuous"), True)
    check("vacuous w2 RESIDUAL", rep0["per_weight"][2]["RESIDUAL"], 3)
    repoff = predict_collapse(d_st, weight_max=3, method="exact_Q",
                              steinmann=False)
    check("disabled honesty: SKIPPED",
          repoff["honesty"]["steinmann"].startswith("SKIPPED"), True)
    check("disabled w3 RESIDUAL", repoff["per_weight"][3]["RESIDUAL"], 4)
    # exact_Q vs mod-p agreement on the tagged 7-letter alphabet: the mod-p
    # word pre-filter and the exact projection must land on the same residual.
    repm = predict_collapse(os.path.join(ALPH, "synth_7letter.json"),
                            weight_max=3, method="modp")
    check("modp w2 steinmann words", repm["per_weight"][2]["steinmann"], 9)
    check("modp w2 RESIDUAL", repm["per_weight"][2]["RESIDUAL"], 5)
    check("modp w3 steinmann words", repm["per_weight"][3]["steinmann"], 25)
    check("modp w3 RESIDUAL", repm["per_weight"][3]["RESIDUAL"], 13)


def T7_per_orbit_adapter():
    print("\n[T7] per-orbit input adapter (recursive_Kw via load_alphabet)")
    # Reference: trivial 2-letter, FE=LE=all → K_w equals the integrable dims
    # T1 checks by hand (w1: 2, w2: 3, w3: 4).
    want = {0: 1, 1: 2, 2: 3, 3: 4}
    r = recursive_Kw(os.path.join(ALPH, "trivial_2var.json"),
                     wmax=3, npts_max=24)
    check("fixture path Kw_LE", r["Kw_LE"], want)
    check("converged", all(c[0] for c in r["converged"].values()), True)
    # landau-alphabet-style output (alphabet/invariants, no letters key)
    r2 = recursive_Kw({"alphabet": ["x", "y"], "invariants": ["x", "y"]},
                      wmax=3, npts_max=24)
    check("landau-style dict Kw_LE", r2["Kw_LE"], want)
    # collapse schema without entry lists → first/last default to ALL letters
    r3 = recursive_Kw({"name": "t", "vars": ["x", "y"],
                       "letters": {"x": "x", "y": "y"}}, wmax=3, npts_max=24)
    check("entry-defaults Kw_LE", r3["Kw_LE"], want)


if __name__ == "__main__":
    T1_trivial()
    T2_P4_J23()
    T3_P4_J25()
    T4_synth7()
    T5_adapter()
    T6_steinmann()
    T7_per_orbit_adapter()
    print("\n" + ("ALL PASS" if not FAIL else f"FAILED: {FAIL}"))
    sys.exit(1 if FAIL else 0)

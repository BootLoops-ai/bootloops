#!/usr/bin/env python3
# derive_posq_sentences.py -- sentence layer for QUARTET_POSQ_TABLE.json,
# BY CODE, mirroring the heirloom derive_sentences.py law:
#   gap(winner, rivals) = winner lnZ LOWER - max over rivals lnZ UPPER
#   displays: lower endpoints / lnBF floors ROUND_FLOOR, uppers ROUND_CEILING
#   TIMESTAMP-OUT: no wall-clock key in pinned bytes (provenance -> stdout);
#   regeneration with unchanged inputs is byte-identical; sha256 printed.
import json, math, sys, time, hashlib, os
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topos

# production input roots from the environment (data not shipped);
# unset => loud SKIP at use (fail-closed, no machine-local defaults)
RESULTS = os.environ.get("POSQ_RESULTS_DIR", "")
SD = os.environ.get("POSQ_SD_DIR", "")
def _require_root(v, name):
    if not v:
        raise SystemExit(f"SKIP: env {name} unset — points at the production "
                         "input root (not shipped)")
    return v

def rdn(x, d):
    return str(Decimal(x).quantize(Decimal(1).scaleb(-d), rounding=ROUND_FLOOR))

def rup(x, d):
    return str(Decimal(x).quantize(Decimal(1).scaleb(-d), rounding=ROUND_CEILING))

# rounding-direction unit regression (heirloom law), runs on every generation
assert rdn(-46.0437, 2) == "-46.05" and rup(-43.2095, 2) == "-43.20"
assert rdn(14.14365, 4) == "14.1436" and rup(-45.898, 2) == "-45.89"
for _v in (0.15, -33.333, 14.1437, -1e-7):
    assert float(rdn(_v, 4)) <= _v <= float(rup(_v, 4)), "outward law violated"

def ival(lo, hi, d=2):
    return f"[{rdn(lo, d)}, {rup(hi, d)}]"

NEG_INF = float("-inf")

def derive(table_path, out_path=None):
    tab = json.load(open(table_path))
    rows = {r["topo_index"]: r for r in tab["per_topology"]}
    assert set(rows) == set(range(15)), "need all 15 topology rows"
    LO = {i: (rows[i]["lnZ"][0] if rows[i]["lnZ"][0] is not None else NEG_INF)
          for i in rows}
    HI = {i: rows[i]["lnZ"][1] for i in rows}
    win = max(rows, key=lambda i: LO[i])
    assert LO[win] > NEG_INF, "winner must have a certified lower endpoint"

    def gap_to(rivals):
        rivals = [i for i in rivals if i != win]
        if not rivals:
            return None
        b = max(rivals, key=lambda i: HI[i])
        v = LO[win] - HI[b]
        return {"certified_lnBF_at_least": v,
                "certified_lnBF_display_floor": rdn(v, 4),
                "certified_log10BF_display_floor": rdn(v / math.log(10), 4),
                "binding_rival": {"topo_index": b,
                                  "topology": rows[b]["topology"],
                                  "lnZ_upper": HI[b]},
                "rival_set_size": len(rivals)}

    out = {
        "derived_from": table_path,
        "rule": "every gap = winner lnZ LOWER bound minus max-over-rivals lnZ "
                "UPPER bound (certified direction); binding rival = largest "
                "upper bound; display: lower endpoints and bound floors "
                "ROUND_FLOOR, upper endpoints ROUND_CEILING (heirloom "
                "derive_sentences.py law, mirrored by code)",
        "model": "posq CTMC: binary CTMC, strict clock, Exp(lambda=10) "
                 "increments, pi1~U(0,1), uncorrected likelihood "
                 "(PINNED_SPEC.md); certified two-sided lnZ_unc enclosures",
        "winner": {"topo_index": win, "topology": rows[win]["topology"],
                   "lnZ": [LO[win], HI[win]],
                   "lnZ_display": ival(LO[win], HI[win])},
        "runner_up_margin": gap_to([i for i in rows if i != win]),
        "split_classes": {},
    }
    wcls = next(c for c, m in topos.SPLIT_CLASSES.items() if win in m)
    out["winner_split_class"] = wcls
    for cname, members in topos.SPLIT_CLASSES.items():
        entry = {"members": sorted(members),
                 "class_best_lnZ_lower": max(LO[i] for i in members),
                 "class_max_lnZ_upper": max(HI[i] for i in members)}
        if cname != wcls:
            b = max(members, key=lambda i: HI[i])
            v = max(LO[i] for i in topos.SPLIT_CLASSES[wcls]) - HI[b]
            entry["certified_gap_from_winner_class"] = {
                "certified_lnBF_at_least": v,
                "certified_lnBF_display_floor": rdn(v, 4),
                "binding_rival": {"topo_index": b,
                                  "topology": rows[b]["topology"]}}
        out["split_classes"][cname] = entry
    # evidence-weighted total (uniform 1/15); missing lowers contribute 0
    mx = max(HI.values())
    lb = math.log(sum(math.exp(LO[i] - mx) for i in rows
                      if LO[i] > NEG_INF) / 15) + mx - 1e-9
    ub = math.log(sum(math.exp(HI[i] - mx) for i in rows) / 15) + mx + 1e-9
    out["evidence_weighted_total_lnZbar"] = [lb, ub]
    # winner posterior weight under uniform 1/15 (certified directions)
    wl = math.exp(LO[win] - mx) / (math.exp(LO[win] - mx)
                                   + sum(math.exp(HI[i] - mx)
                                         for i in rows if i != win))
    wu = math.exp(HI[win] - mx) / (math.exp(HI[win] - mx)
                                   + sum(math.exp(LO[i] - mx)
                                         for i in rows
                                         if i != win and LO[i] > NEG_INF))
    out["winner_posterior_weight"] = [wl - 1e-12, min(wu + 1e-12, 1.0)]
    # cross-model check vs heirloom SD MAP (the robustness sentence)
    sd = json.load(open(f"{_require_root(SD, 'POSQ_SD_DIR')}/DERIVED_SENTENCES.json"))
    sd_win = sd["winner"]
    agree = (sd_win["topo_index"] == win)
    strict = LO[win] > max(HI[i] for i in rows if i != win)
    out["cross_model_check"] = {
        "posq_ctmc_map": {"topo_index": win,
                          "topology": rows[win]["topology"],
                          "strictly_certified_map": bool(strict)},
        "heirloom_sd_map": {"topo_index": sd_win["topo_index"],
                            "topology": sd_win["topology"]},
        "agree": bool(agree)}
    g_ru = out["runner_up_margin"]
    sents = [
        f"MAP topology (posq CTMC): {rows[win]['topology']} with lnZ_unc in "
        f"{ival(LO[win], HI[win])} (certified enclosure).",
        f"Certified lnBF >= {g_ru['certified_lnBF_display_floor']} over the "
        f"runner-up (binding rival {g_ru['binding_rival']['topology']}, by "
        f"UPPER bound).",
    ]
    for cname in sorted(out["split_classes"]):
        e = out["split_classes"][cname]
        if "certified_gap_from_winner_class" in e:
            g = e["certified_gap_from_winner_class"]
            sents.append(
                f"Certified lnBF >= {g['certified_lnBF_display_floor']} of "
                f"the winner class ({wcls}) over split class {cname} "
                f"(binding rival {g['binding_rival']['topology']}).")
    sents.append(
        f"Evidence-weighted total lnZbar in {ival(lb, ub)} under the "
        f"registered uniform 1/15.")
    sents.append(
        f"Cross-model: the posq CTMC MAP topology "
        f"{rows[win]['topology']} "
        f"{'AGREES with' if agree else 'DIFFERS from'} the heirloom SD MAP "
        f"{sd_win['topology']} -- "
        + ("the two independently-registered models select the same rooted "
           "topology." if agree else
           "the two independently-registered models select different rooted "
           "topologies."))
    out["sentences"] = sents
    print("generated:", time.strftime("%Y-%m-%d %H:%M:%S %Z"),
          "(stdout provenance only; NOT in the pinned bytes)")
    if out_path:
        js = json.dumps(out, indent=1)
        open(out_path, "w").write(js)
        # byte-stability: regenerate in-memory and compare
        assert json.dumps(out, indent=1) == js
        sha = hashlib.sha256(js.encode()).hexdigest()
        print(f"sha256({os.path.basename(out_path)}) = {sha}")
    print(json.dumps(out["sentences"], indent=1))
    return out

if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else f"{_require_root(RESULTS, 'POSQ_RESULTS_DIR')}/QUARTET_POSQ_TABLE.json"
    dst = sys.argv[2] if len(sys.argv) > 2 else f"{_require_root(RESULTS, 'POSQ_RESULTS_DIR')}/QUARTET_POSQ_SENTENCES.json"
    o1 = derive(src, dst)
    o2 = derive(src, None)
    assert json.dumps(o1, indent=1) == json.dumps(o2, indent=1), \
        "double-derivation byte-identity FAILED"
    print("double derivation byte-identical: PASS")

"""Assemble the verdict table + gate table from the CERT.json receipts
(mechanical; the report prose is written by the operator around it)."""
import json
import os
import sys

sys.dont_write_bytecode = True

from . import RECEIPTS_DIR


def row(tag, c, product):
    kw = c.get("contraction_margins") or {}
    pd = c.get("pd_margins") or {}
    mx = max(kw.values()) if kw else None
    pdm = min(pd.values()) if pd else None
    rmin = c.get("dims") or {}
    rad = c.get("radius")
    return dict(tag=tag, product=product, verdict=c["verdict"],
                radius=rad, max_margin=mx, min_pd=pdm,
                failed=c.get("failed"),
                n=(c.get("dims") or {}).get("n_total"),
                gates_ok=_gates_ok(c))


def _gates_ok(c):
    g = c.get("oracle_identity_gates") or {}
    oks = []
    for k, v in g.items():
        if isinstance(v, dict) and "ok" in v and v["ok"] is not None:
            oks.append(bool(v["ok"]))
    return all(oks) if oks else None


def collect():
    rows = []
    for d in sorted(os.listdir(RECEIPTS_DIR)):
        p = os.path.join(RECEIPTS_DIR, d, "CERT.json")
        if os.path.isfile(p):
            c = json.load(open(p))
            rows.append(row(d, c, c.get("product", "direct")))
        pp = os.path.join(RECEIPTS_DIR, d, "polished", "CERT.json")
        if os.path.isfile(pp):
            c = json.load(open(pp))
            rows.append(row(d + "/polished", c, "polished_center"))
    return rows


def fmt(v, spec="{:.3g}"):
    return "—" if v is None else spec.format(v)


def main():
    rows = collect()
    print("| fit | product | verdict | radius (u) | worst contraction "
          "margin | min PD margin | gates |")
    print("|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['tag']} | {r['product']} | **{r['verdict']}** | "
              f"{fmt(r['radius'])} | {fmt(r['max_margin'])} | "
              f"{fmt(r['min_pd'])} | "
              f"{'ok' if r['gates_ok'] else r['gates_ok']} |")


if __name__ == "__main__":
    main()

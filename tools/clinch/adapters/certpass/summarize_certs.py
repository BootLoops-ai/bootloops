#!/usr/bin/env python3
"""summarize_certs — collect every banked CERT.json into one table
(markdown rows + a JSON master) for the addendum."""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "receipts")


def fmt(v, g="{:.3g}"):
    return "-" if v is None else g.format(v)


def main():
    rows = []
    for p in sorted(glob.glob(os.path.join(R, "*", "CERT.json"))):
        c = json.load(open(p))
        tag = c["model"]
        kw = c.get("contraction_margins") or {}
        pd = c.get("pd_margins") or {}
        shift = c.get("shift_vs_banked") or {}
        rows.append(dict(
            tag=tag, verdict=c["verdict"],
            failed=c.get("failed"),
            radius=c.get("radius"), radius_min=c.get("radius_min"),
            worst_kw=(max(kw.values()) if kw else None),
            worst_pd=(min(pd.values()) if pd else None),
            dims=(c.get("dims") or {}).get("n_total"),
            dist=shift.get("distance_inf"),
            f_drop=shift.get("f_drop", shift.get("nll_drop")),
            stamp=c.get("stamp"), path=os.path.relpath(p, HERE)))
    sp = os.path.join(R, "r4c_selection", "SCOPE_ADJUDICATION.json")
    if os.path.exists(sp):
        c = json.load(open(sp))
        rows.append(dict(tag="r4c_selection (T3)", verdict=c["verdict"],
                         failed=None, radius=None, radius_min=None,
                         worst_kw=None, worst_pd=None, dims=None,
                         dist=None, f_drop=None, stamp=c["stamp"],
                         path=os.path.relpath(sp, HERE)))
    out = os.path.join(R, "CERT_TABLE.json")
    json.dump(rows, open(out, "w"), indent=1)
    print(f"| target | p | verdict | worst Krawczyk margin | worst PD "
          f"margin | radii | shift vs banked (inf) | nll drop |")
    print("|---|---|---|---|---|---|---|---|")
    for r in rows:
        rad = (f"[{fmt(r['radius_min'])}, {fmt(r['radius'])}]"
               if r["radius"] is not None else "-")
        print(f"| {r['tag']} | {r['dims']} | {r['verdict']}"
              f"{('('+r['failed']+')') if r['failed'] else ''} | "
              f"{fmt(r['worst_kw'])} | {fmt(r['worst_pd'])} | {rad} | "
              f"{fmt(r['dist'])} | {fmt(r['f_drop'])} |")
    print("->", out)


if __name__ == "__main__":
    main()

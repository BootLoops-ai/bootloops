#!/usr/bin/env python3
"""RECEIPT — reduction-certificate inspection: every table ships with a receipt.

Per-(kinematic point, prime) lambda-multiplier certificates for IBP/Laporta
reduction tables: verify (independent, ~1 ms/row), emit (retrofit any
existing table), detect (false-master/structural audit).

    receipt.py verify --witness ...  [--system-jsonl S | --adapter strata ...]
    receipt.py emit   --adapter strata --system-dir D --staging STG \
                      --family F --table kira_target.m --slices p:d:eta[,...]
    receipt.py detect [--adapter strata ... | --system-jsonl S --masters M]

rc: 0 = all rows PASS / no alarms; 1 = any row FAILED / alarm; 2 = malformed
input. See README.md + WITNESS_FORMAT.md.

The tool name is a constant below; the directory is rename-trivial (all
internal paths are __file__-relative, all imports are same-directory).
"""
import argparse
import glob
import json
import os
import sys
import time

TOOL_NAME = "RECEIPT"

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import core                                    # noqa: E402
from core import (MalformedInput, RC_FAIL, RC_MALFORMED, RC_OK,  # noqa: E402
                  load_system_jsonl, load_witness, verify_many,
                  write_witness)


def _die_malformed(msg):
    print(json.dumps({"tool": TOOL_NAME, "rc": RC_MALFORMED, "error": msg}))
    sys.exit(RC_MALFORMED)


def _expand_witness_args(items):
    paths = []
    for it in items:
        if os.path.isdir(it):
            paths.extend(sorted(glob.glob(os.path.join(it, "*.json"))))
        else:
            hits = sorted(glob.glob(it))
            if not hits:
                _die_malformed(f"no witness files match {it!r}")
            paths.extend(hits)
    return paths


def _parse_slices(s):
    out = []
    for part in s.split(","):
        p, d0, e0 = (int(x) for x in part.split(":"))
        out.append((p, d0, e0))
    return out


def _strata_rows_for(args, p, d0, e0):
    from adapters import strata
    sysd = strata.load_rows(args.system_dir, args.family, p, d0, e0)
    return sysd


def _load_table_json(path, p):
    try:
        raw = json.load(open(path))
        return {int(t): {int(m): int(v) % p for m, v in row.items()}
                for t, row in raw.items()}
    except (OSError, json.JSONDecodeError, TypeError, ValueError) as e:
        raise MalformedInput(f"{path}: unreadable table ({e})")


# ------------------------------------------------------------------- verify
def cmd_verify(args):
    try:
        paths = _expand_witness_args(args.witness)
        if not paths:
            _die_malformed("no witness files given")
        wits = []
        legacy = [q for q in paths if os.path.basename(q).startswith("c_")]
        v1 = [q for q in paths if q not in set(legacy)]
        for q in v1:
            wits.append(load_witness(q))
    except MalformedInput as e:
        _die_malformed(str(e))

    # group by (p, point): each group verifies against its own system eval
    def point_of(w):
        pt = (w["meta"].get("point") or {})
        return (w["p"], pt.get("d"), pt.get("eta"))

    groups = {}
    for w in wits:
        groups.setdefault(point_of(w), []).append(w)

    all_results, n_pass, n_rows = [], 0, 0
    t_verify = 0.0
    for (p, d0, e0), group in sorted(groups.items()):
        if args.system_jsonl:
            try:
                rows = load_system_jsonl(args.system_jsonl)
            except MalformedInput as e:
                _die_malformed(str(e))
            table = (_load_table_json(args.table, p) if args.table else None)
        elif args.adapter == "strata":
            if d0 is None or e0 is None:
                _die_malformed("witness lacks point{d,eta} for strata adapter")
            sysd = _strata_rows_for(args, p, d0, e0)
            rows = sysd["rows"]
            table = (_load_table_json(args.table, p) if args.table else None)
        else:
            _die_malformed("need --system-jsonl or --adapter strata")
        # legacy pairs need the adapter's masters map
        if legacy:
            if args.adapter != "strata":
                _die_malformed("legacy c_*.json witnesses need --adapter strata")
            from adapters import strata as _st
            info = getattr(cmd_verify, "_info_cache", None)
            if info is None:
                sector_of = {w: s for w, s in sysd["sector_of"].items()}
                _sysd, info, _closed = _st.stage_partition(
                    args.system_dir, args.staging, args.family, p, d0, e0,
                    sector_of)
                cmd_verify._info_cache = info
            for q in list(legacy):
                lw = _st.read_legacy_pair(q, info["masters_map"], len(rows))
                if lw["p"] == p:
                    group.append(lw)
                    legacy.remove(q)
        t0 = time.time()
        rc, rep = verify_many(rows, group, table=table)
        t_verify += time.time() - t0
        all_results.extend(rep["results"])
        n_pass += rep["n_pass"]
        n_rows += rep["n_rows"]

    if legacy:
        # legacy witnesses at a prime no v1 witness pinned: verify per file
        from adapters import strata as _st
        for q in legacy:
            meta = json.load(open(q))
            p, d0, e0 = int(meta["p"]), int(meta["d0"]), int(meta["eta0"])
            sysd = _strata_rows_for(args, p, d0, e0)
            info = getattr(cmd_verify, "_legacy_info", {}).get((p, d0, e0))
            if info is None:
                _s, info, _c = _st.stage_partition(
                    args.system_dir, args.staging, args.family, p, d0, e0,
                    sysd["sector_of"])
                cmd_verify._legacy_info = getattr(
                    cmd_verify, "_legacy_info", {})
                cmd_verify._legacy_info[(p, d0, e0)] = info
            lw = _st.read_legacy_pair(q, info["masters_map"],
                                      len(sysd["rows"]))
            t0 = time.time()
            rc, rep = verify_many(sysd["rows"], [lw])
            t_verify += time.time() - t0
            all_results.extend(rep["results"])
            n_pass += rep["n_pass"]
            n_rows += rep["n_rows"]

    ok = (n_rows > 0 and n_pass == n_rows)
    out = {"tool": TOOL_NAME, "cmd": "verify", "n_rows": n_rows,
           "n_pass": n_pass, "n_fail": n_rows - n_pass, "all_pass": ok,
           "verify_wall_s": round(t_verify, 4),
           "verify_ms_per_row": (round(1000 * t_verify / n_rows, 3)
                                 if n_rows else None),
           "results": all_results}
    if args.report:
        with open(args.report, "w") as fh:
            json.dump(out, fh, indent=1)
    print(json.dumps(out if not args.quiet else
                     {k: v for k, v in out.items() if k != "results"},
                     indent=1))
    sys.exit(RC_OK if ok else (RC_FAIL if n_rows else RC_MALFORMED))


# --------------------------------------------------------------------- emit
def cmd_emit(args):
    import numpy as np
    from emitter import FpSystem, emit_for_target

    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    if args.adapter != "strata":
        _emit_generic(args, rng)
        return
    from adapters import strata
    from adapters.strata import CoeffEvaluator, parse_oracle_table

    try:
        slices = _parse_slices(args.slices)
    except ValueError:
        _die_malformed(f"bad --slices {args.slices!r}")
    struct_rows, exprs, sector_of = strata.parse_structural(
        args.system_dir, args.family)
    p1, d1, e1 = slices[0]
    sysd_a, info, closed_a = strata.stage_partition(
        args.system_dir, args.staging, args.family, p1, d1, e1, sector_of)
    assert len(info["pivots"]) == len(struct_rows) and \
        info["n_unclassified"] == 0, "column partition not square/complete"
    ev_a = CoeffEvaluator(p1, d1, e1)
    orc = parse_oracle_table(args.table)
    targets, skip = strata.anchor_table(orc, info, closed_a, sector_of,
                                        ev_a, p1)
    if args.max_rows:
        targets = targets[:args.max_rows]
    del closed_a, sysd_a
    masters_map = info["masters_map"]
    midx2w = {idx: w for w, idx in masters_map.items()}
    master_cols = set(masters_map)

    summary = {"tool": TOOL_NAME, "cmd": "emit", "family": args.family,
               "n_oracle_rows": len(orc), "n_anchored": len(targets),
               "skips": skip, "slices": []}
    any_fail = False
    for (p, d0, e0) in slices:
        ev = CoeffEvaluator(p, d0, e0)
        vals = np.array([ev(e) for e in exprs], np.int64)
        rows_terms = [[(w, int(vals[k])) for k, w in row if vals[k]]
                      for row in struct_rows]
        S = FpSystem(rows_terms, info["pivots"], sorted(sector_of), p)
        fp_rows = [dict(rt) for rt in rows_terms]
        fpr = core.system_fingerprint(fp_rows, p)
        minpoly = None
        recs = []
        for tgt, terms, w_t in targets:
            claimed = strata.oracle_c(terms, ev, midx2w, p)
            if isinstance(claimed, str):
                recs.append({"target": list(tgt), "status": f"ORACLE_{claimed}"})
                any_fail = True
                continue
            rec, lam, mp_new = emit_for_target(
                S, w_t, master_cols, rng, minpoly=minpoly, claimed=claimed)
            if mp_new is not None:
                minpoly = mp_new
            rec["target"] = list(tgt)
            if lam is not None and rec["status"] == "CERTIFIED":
                nz = np.flatnonzero(lam)
                wp = os.path.join(
                    args.out, f"w_{args.family}_{p}_{w_t}.json")
                write_witness(
                    wp, p=p, target_col=w_t, c=rec["c"],
                    lam_idx=[int(i) for i in nz],
                    lam_val=[int(lam[i]) for i in nz],
                    n_rows=S.n_eqs, family=args.family,
                    point={"d": d0, "eta": e0}, target_label=list(tgt),
                    labels={w: list(masters_map[w]) for w in rec["c"]},
                    system={"n_rows": S.n_eqs, "fingerprint": fpr,
                            "source": os.path.abspath(args.system_dir)})
                rec["witness"] = os.path.basename(wp)
            else:
                any_fail = True
            rec.pop("c", None)
            recs.append(rec)
        n_cert = sum(1 for r in recs if r["status"] == "CERTIFIED")
        summary["slices"].append(
            {"slice": [p, d0, e0], "n_rows": len(recs),
             "n_certified": n_cert,
             "n_claim_mismatch": sum(1 for r in recs
                                     if r["status"] == "CLAIM_MISMATCH"),
             "solve_wall_total_s": round(sum(r.get("solve_wall_s", 0)
                                             for r in recs), 2),
             "rows": recs})
    with open(os.path.join(args.out, "EMIT_REPORT.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps({**{k: v for k, v in summary.items()
                         if k != "slices"},
                      "slices": [{k: v for k, v in sl.items() if k != "rows"}
                                 for sl in summary["slices"]]}, indent=1))
    sys.exit(RC_FAIL if any_fail else RC_OK)


def _emit_generic(args, rng):
    """Generic emit: JSONL system + pivot list + one target col."""
    import numpy as np
    from emitter import FpSystem, emit_for_target
    try:
        rows = load_system_jsonl(args.system_jsonl)
        pivots = json.load(open(args.pivots))
        p = int(args.p)
        masters = set(json.load(open(args.masters))) if args.masters else \
            {c for r in rows for c in r} - set(pivots)
        claimed = None
        if args.table:
            tb = _load_table_json(args.table, p)
            claimed = tb.get(int(args.target_col))
    except (MalformedInput, OSError, json.JSONDecodeError, ValueError) as e:
        _die_malformed(str(e))
    rows_terms = [list(r.items()) for r in rows]
    S = FpSystem(rows_terms, pivots, sorted({c for r in rows for c in r}), p)
    rec, lam, _ = emit_for_target(S, int(args.target_col), masters, rng,
                                  claimed=claimed)
    if lam is not None and rec["status"] == "CERTIFIED":
        nz = np.flatnonzero(lam)
        wp = os.path.join(args.out, f"w_generic_{p}_{args.target_col}.json")
        write_witness(wp, p=p, target_col=int(args.target_col), c=rec["c"],
                      lam_idx=[int(i) for i in nz],
                      lam_val=[int(lam[i]) for i in nz], n_rows=S.n_eqs,
                      system={"n_rows": S.n_eqs,
                              "fingerprint": core.system_fingerprint(rows, p)})
        rec["witness"] = wp
    print(json.dumps({"tool": TOOL_NAME, "cmd": "emit", **{
        k: v for k, v in rec.items() if k != "c"}}, indent=1))
    sys.exit(RC_OK if rec["status"] == "CERTIFIED" else RC_FAIL)


# ------------------------------------------------------------------- detect
def cmd_detect(args):
    from detector import detect
    if args.adapter == "strata":
        from adapters import strata
        try:
            slices = _parse_slices(args.slices)
        except ValueError:
            _die_malformed(f"bad --slices {args.slices!r}")
        p, d0, e0 = slices[0]
        sysd = strata.load_rows(args.system_dir, args.family, p, d0, e0)
        sector_of = sysd["sector_of"]
        # declared masters: staged census (preferred ordinals + survivor
        # buckets need the partition; the DECLARED set here = the loader's
        # master-class weights + optional extra columns from --masters)
        masters = set(sysd["masters"])
        shell = set()
        if args.staging:
            _s, info, _c = strata.stage_partition(
                args.system_dir, args.staging, args.family, p, d0, e0,
                sector_of)
            masters = set(info["masters_map"])
            shell = set(info["shell"])
        if args.masters:
            masters = set(json.load(open(args.masters)))
        order = {w: w for w in sector_of if w not in masters}
        rows = sysd["rows"]
        table = _load_table_json(args.table, p) if args.table else None
        targets = json.load(open(args.targets)) if args.targets else ()
    else:
        try:
            rows = load_system_jsonl(args.system_jsonl)
            p = int(args.p)
            masters = set(json.load(open(args.masters)))
            shell = set(json.load(open(args.shell))) if args.shell else ()
            table = _load_table_json(args.table, p) if args.table else None
            targets = json.load(open(args.targets)) if args.targets else ()
        except (MalformedInput, OSError, json.JSONDecodeError,
                TypeError, ValueError) as e:
            _die_malformed(str(e))
        order = None
    t0 = time.time()
    alarm, report = detect(rows, p, masters, order=order, shell=shell,
                           targets=targets, table=table)
    report.update({"tool": TOOL_NAME, "cmd": "detect", "p": p,
                   "wall_s": round(time.time() - t0, 3)})
    if args.report:
        with open(args.report, "w") as fh:
            json.dump(report, fh, indent=1)
    print(json.dumps(report, indent=1))
    sys.exit(RC_FAIL if alarm else RC_OK)


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(
        prog=TOOL_NAME.lower(), description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(g):
        g.add_argument("--adapter", choices=["strata"], default=None)
        g.add_argument("--system-dir", help="kira artifact dir (adapter)")
        g.add_argument("--staging", help="staged config dir (adapter)")
        g.add_argument("--family", help="family name (adapter)")
        g.add_argument("--system-jsonl", help="generic system snapshot")
        g.add_argument("--table", help="claimed table (JSON col-keyed)")
        g.add_argument("--report", help="write full JSON report here")

    v = sub.add_parser("verify", help="verify witnesses (solver-agnostic)")
    common(v)
    v.add_argument("--witness", nargs="+", required=True,
                   help="witness file(s)/glob(s)/dir(s); c_*.json = legacy")
    v.add_argument("--quiet", "-q", action="store_true")

    e = sub.add_parser("emit", help="emit lambda-witnesses (retrofit driver)")
    common(e)
    e.add_argument("--slices", help="p:d:eta[,p:d:eta...] (adapter)")
    e.add_argument("--out", required=True)
    e.add_argument("--max-rows", type=int, default=0)
    e.add_argument("--seed", type=int, default=20260706)
    e.add_argument("--pivots", help="JSON pivot col list (generic)")
    e.add_argument("--masters", help="JSON master col list (generic)")
    e.add_argument("--target-col", help="target column id (generic)")
    e.add_argument("--p", help="prime (generic)")

    d = sub.add_parser("detect", help="false-master/structural audit")
    common(d)
    d.add_argument("--slices", help="p:d:eta (adapter; first slice used)")
    d.add_argument("--masters", help="JSON declared master col list")
    d.add_argument("--shell", help="JSON shell col list (generic)")
    d.add_argument("--targets", help="JSON demanded target col list")
    d.add_argument("--p", help="prime (generic)")

    args = ap.parse_args()
    {"verify": cmd_verify, "emit": cmd_emit, "detect": cmd_detect}[args.cmd](args)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""gatekeeper.registry — tower/values oracle index for per-point oracle farms.

Home of record for the flat tool `tools/oracle_registry.py`.
The flat path stays as an attribute-complete script-mode shim that delegates
here; the engine, CLI verbs (--rebuild/--point/--dupes/--selftest) and the
ORACLE_REG_* env overrides are UNCHANGED from the original registration.

Scans values_*.json / *towers*.json under the configured roots, indexes
(s,t) point, goal, eps orders, integral
lists + sha1 (dedup: byte-identical towers filed under two master labels
inflate the apparent master coverage — the duplicate-oracle gotcha),
and master coverage via combos.json (masters 0-12) / combos_full41.json
(masters 13-40 identity combos).  Index: oracle/REGISTRY.json under the registry root.
(append-safe rebuild: vanished files kept, marked stale).

Usage:
  python -m gatekeeper.registry --rebuild
  python -m gatekeeper.registry --point -2,-3 --masters 13..40 --min-goal 50
  python -m gatekeeper.registry --point -4,-12/7 --masters 0..12
  python -m gatekeeper.registry --dupes
  python -m gatekeeper.registry --selftest
  # legacy call form (shim, identical behaviour):
  python3 tools/oracle_registry.py --dupes

Env configuration — point the engine at your own oracle tree:
  ORACLE_REG_ROOTS     colon-separated root globs to scan
                       (default: oracle_dumps/*)
  ORACLE_REG_REGISTRY  registry json path (default: oracle/REGISTRY.json)
  ORACLE_REG_COMBOS    combos dir (default: dirname of the registry;
                       combos.json rides with the registry)
Without ORACLE_REG_ROOTS pointing at real dumps, --selftest coverage
checks fail environmentally — that is expected on a bare clone.
"""
import argparse, glob, hashlib, json, os, sys, time
from fractions import Fraction

ROOT_GLOBS = (os.environ["ORACLE_REG_ROOTS"].split(":")
              if os.environ.get("ORACLE_REG_ROOTS") else
              ["oracle_dumps/*"])  # override via ORACLE_REG_ROOTS
REGISTRY = os.environ.get("ORACLE_REG_REGISTRY", "oracle/REGISTRY.json")
COMBOS_DIR = os.environ.get("ORACLE_REG_COMBOS", os.path.dirname(REGISTRY))
MAX_BYTES = 200 * 1024 * 1024
MAX_DEPTH = 4

def frac(x):
    try: return str(Fraction(str(x)))
    except Exception: return str(x)

def sha1_of(obj):
    return hashlib.sha1(json.dumps(obj, sort_keys=True).encode()).hexdigest()[:16]

def load_combos():
    maps = {}
    p = os.path.join(COMBOS_DIR, "combos.json")
    if os.path.exists(p):
        for m, combo in json.load(open(p))["combos"].items():
            maps[int(m)] = [tuple(t[1]) for t in combo]
    p = os.path.join(COMBOS_DIR, "combos_full41.json")
    if os.path.exists(p):
        for i, idx in enumerate(json.load(open(p))["unique"]):
            maps[13 + i] = [tuple(idx)]
    return maps

def parse_file(path, combo_maps):
    rec = {"path": path, "mtime": os.path.getmtime(path), "size": os.path.getsize(path)}
    try: d = json.load(open(path))
    except Exception as e:
        rec.update(format="unreadable", error=str(e)[:120], points=[]); return rec
    pts = []
    def values_point(d):
        ints = [list(i) for i in d.get("integrals", [])]
        iset = set(tuple(i) for i in ints)
        cov = sorted(m for m, need in combo_maps.items() if need and all(t in iset for t in need))
        eo = sorted({int(o) for r in d.get("result", []) for o in r.get("coeffs", {})})
        return {"s": frac(d["s"]), "t": frac(d["t"]), "goal": int(d["goal"]),
                "goal_src": "field", "eps_order": d.get("eps_order"),
                "eps_orders": eo, "n_integrals": len(ints), "integrals": ints,
                "int_sha1": sha1_of(sorted(ints)), "masters_covered": cov}
    def towers_point(d):
        mk = sorted(int(k) for k in d["masters"])
        lens = [len(v) for tw in d["masters"].values() for v in
                (tw.values() if isinstance(tw, dict) else [])
                if isinstance(v, str)]
        eo = sorted({int(o) for tw in d["masters"].values() if isinstance(tw, dict) for o in tw})
        return {"s": frac(d["s"]), "t": frac(d["t"]),
                "goal": (min(lens) - 2 if lens else None), "goal_src": "strlen",
                "eps_orders": eo, "n_integrals": len(mk), "integrals": None,
                "int_sha1": sha1_of(mk), "masters_covered": mk}
    if isinstance(d, dict) and {"s", "t", "goal", "integrals", "result"} <= set(d):
        rec["format"] = "values"; pts.append(values_point(d))
    elif isinstance(d, dict) and {"s", "t", "masters"} <= set(d):
        rec["format"] = "towers"; pts.append(towers_point(d))
    elif isinstance(d, dict) and d and all(isinstance(v, dict) and {"s", "t", "masters"} <= set(v) for v in d.values()):
        rec["format"] = "towers_multi"
        for tag, blk in d.items():
            p = towers_point(blk); p["block"] = tag; pts.append(p)
    else:
        rec["format"] = "unknown"
    rec["points"] = pts
    return rec

def rebuild():
    t0 = time.time(); combo_maps = load_combos()
    found = []
    for g in ROOT_GLOBS:
        for root in sorted(glob.glob(g)):
            base_depth = root.rstrip("/").count("/")
            if os.path.isfile(root): continue
            for dp, dns, fns in os.walk(root):
                if dp.count("/") - base_depth >= MAX_DEPTH: dns[:] = []
                for fn in fns:
                    if fn.endswith(".json") and (fn.startswith("values_") or "towers" in fn or "master_tower" in fn):
                        p = os.path.join(dp, fn)
                        try:
                            if os.path.getsize(p) > MAX_BYTES: continue
                        except OSError: continue
                        found.append(p)
    old = {}
    if os.path.exists(REGISTRY):
        try: old = {e["path"]: e for e in json.load(open(REGISTRY))["files"]}
        except Exception: old = {}
    for e in old.values(): e["stale"] = not os.path.exists(e["path"])
    for p in sorted(set(found)):
        e = old.get(p)
        if e and not e.get("stale") and e.get("mtime") == os.path.getmtime(p) and e.get("points"):
            e["stale"] = False; continue
        old[p] = parse_file(p, combo_maps); old[p]["stale"] = False
    out = {"built": time.time(), "scan_seconds": round(time.time() - t0, 2),
           "n_files": len(found), "files": sorted(old.values(), key=lambda e: e["path"])}
    os.makedirs(os.path.dirname(REGISTRY), exist_ok=True)
    tmp = REGISTRY + ".tmp"; json.dump(out, open(tmp, "w"), indent=1); os.replace(tmp, REGISTRY)
    print("rebuilt %s: %d files, %.2fs scan" % (REGISTRY, len(found), out["scan_seconds"]))
    return out

def load_reg():
    if not os.path.exists(REGISTRY): sys.exit("no registry; run --rebuild first")
    return json.load(open(REGISTRY))

def query(point, masters, min_goal, quiet=False):
    s, t = (frac(x) for x in point.split(","))
    reg = load_reg(); cover = {}
    for e in reg["files"]:
        if e.get("stale"): continue
        for p in e["points"]:
            if p["s"] == s and p["t"] == t and p.get("goal") is not None and p["goal"] >= min_goal:
                for m in p["masters_covered"]:
                    cover.setdefault(m, []).append((p["goal"], e["path"]))
    missing = [m for m in masters if m not in cover]
    if not quiet:
        print("point (%s,%s)  min_goal=%d  masters %d..%d" % (s, t, min_goal, masters[0], masters[-1]))
        for m in masters:
            if m in cover:
                g, f = max(cover[m])
                print("  m%-3d best_goal=%-4d %s (+%d more)" % (m, g, f, len(cover[m]) - 1))
        print("MISSING (%d): %s" % (len(missing), missing if missing else "none"))
    return cover, missing

def dupes(quiet=False):
    reg = load_reg(); groups = {}
    for e in reg["files"]:
        if e.get("stale"): continue
        for p in e["points"]:
            tag = "#" + p["block"] if p.get("block") else ""
            groups.setdefault((p["s"], p["t"], p["int_sha1"]), []).append(e["path"] + tag)
    bad = {k: sorted(set(v)) for k, v in groups.items() if len(set(v)) > 1}
    if not quiet:
        print("duplicate (point, integral-hash) groups: %d" % len(bad))
        for (s, t, h), fs in sorted(bad.items()):
            print("  (%s,%s) sha1=%s:" % (s, t, h))
            for f in fs: print("    %s" % f)
    return bad

def parse_masters(spec):
    a, b = spec.split(".."); return list(range(int(a), int(b) + 1))

def selftest():
    ok = True
    def chk(name, cond):
        nonlocal ok; ok &= bool(cond); print("%s %s" % ("PASS" if cond else "FAIL", name))
    t0 = time.time(); rebuild(); dt = time.time() - t0
    chk("rebuild under 60s (%.1fs)" % dt, dt < 60)
    cov, miss = query("-2,-3", list(range(13)), 110, quiet=True)
    chk("(-2,-3) 13-block covered at goal>=110", not miss)
    chk("(-2,-3) 13-block includes values_P0hi",
        any("values_P0hi" in f for m in cov for _, f in cov[m]))
    cov2, miss2 = query("-4,-12/7", list(range(13)), 30, quiet=True)
    chk("(-4,-12/7) 13-block covered (values_PCRV)", not miss2 and
        any("values_PCRV" in f for m in cov2 for _, f in cov2[m]))
    _, miss3 = query("-2,-3", list(range(13, 41)), 1, quiet=True)
    chk("full-41 at (-2,-3) MISSING (P0f41b still running): %d/28 missing" % len(miss3), len(miss3) == 28)
    print("--- dupes output ---"); dupes()
    print("SELFTEST", "PASS" if ok else "FAIL"); return ok

def main():
    # allow `--point -2,-3` (leading-dash value) by folding it into --point=...
    av = sys.argv[1:]; sys.argv = sys.argv[:1]
    while av:
        x = av.pop(0)
        if x in ("--point", "--masters", "--min-goal") and av:
            x += "=" + av.pop(0)
        sys.argv.append(x)
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true")
    ap.add_argument("--point"); ap.add_argument("--masters", default="0..40")
    ap.add_argument("--min-goal", type=int, default=1)
    ap.add_argument("--dupes", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest: sys.exit(0 if selftest() else 1)
    if a.rebuild: rebuild()
    if a.point: query(a.point, parse_masters(a.masters), a.min_goal)
    if a.dupes: dupes()
    if not (a.rebuild or a.point or a.dupes): ap.print_help()

if __name__ == "__main__":
    main()

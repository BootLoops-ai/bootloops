#!/usr/bin/env python3
"""corpus_bench.py — STRATA end-to-end certification driver.

Per family:
  A. kira Generate of the solve box (or --reuse-artifacts), then at every
     pre-registered slice (p, d, eta): Pass-2 stratified Route-B F_p solve
     (B2) with per-stratum banking into T, assembly strictly from T read-back.
  B. dictionary run: kira full-box select_mandatory_list + triangular +
     back-substitution + kira2math; weight<->integral dictionary by
     reduction-vector matching at the ANCHOR slice, held out on the others;
     doubles as a full-table cross-check (kira algebra vs Route-B).
  C. fresh identities from the INDEPENDENT vacuum IBP generator (ibp_gen),
     translated to weight space via the dictionary; coverage-counted
     residual gate (c_min) + planted-fault mutation trials (scratch copies).
  D. oracle row-match vs every cached corpus instance (anchor at slice 1,
     held out on slices 2-4; zero-normalization applied; 100 pct or VOID).

NOT-VALID-OFF-GRID: all banked values are slice-numeric.
"""
import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coverage import (assemble_closed, coverage_gate,          # noqa: E402
                      planted_fault_trials)
from fp_eliminate import b2b3_compare, stratified_solve        # noqa: E402
from ibp_gen import identities, parse_family_spec, seed_box    # noqa: E402
from interface_table import InterfaceTable                     # noqa: E402
from loader import (CoeffEvaluator, load_system,               # noqa: E402
                    load_trivial_sectors, parse_masters_file)
from weightdict import (MultiEvalSome, config_has_eta,         # noqa: E402
                        enumerate_box, kira_table_path, match_dictionary,
                        run_dict_kira, run_dict_kira_numeric,
                        run_dict_kira_numeric_persector, stream_eval_table,
                        stream_rows)
from weights import popcount                                   # noqa: E402


def maxrss_kb(who=resource.RUSAGE_SELF):
    """Peak RSS in KiB: getrusage ru_maxrss is KiB on Linux, BYTES on macOS —
    normalize so the report's *_kb keys mean what they say on both."""
    rss = resource.getrusage(who).ru_maxrss
    return rss // 1024 if sys.platform == "darwin" else rss


JOBS_TMPL = """jobs:
 - reduce_sectors:
    reduce:
     - {{topologies: [{fam}], sectors: [{secs}], r: {r}, s: {s}}}
    select_integrals:
      select_mandatory_recursively:
        - {{topologies: [{fam}], sectors: [{secs}], r: {r}, s: {s}, d: {d}}}
    preferred_masters: preferred
    integral_ordering: 5
    run_initiate: true
    run_firefly: false
"""

# --gen-select list: provenance-faithful reproduction of production
# pins that selected via an explicit mandatory target list (sec1015 pin
# jobs.yaml shape) instead of the recursive box select.
JOBS_TMPL_LIST = """jobs:
 - reduce_sectors:
    reduce:
     - {{topologies: [{fam}], sectors: [{secs}], r: {r}, s: {s}, d: {d}}}
    select_integrals:
      select_mandatory_list:
        - [{fam}, target]
    preferred_masters: preferred
    integral_ordering: 5
    run_initiate: true
    run_firefly: false
"""


def run_kira(workdir, fam, staging, box, tag, gen_select="recursive"):
    os.makedirs(workdir, exist_ok=True)
    shutil.copytree(os.path.join(staging, "config"),
                    os.path.join(workdir, "config"), dirs_exist_ok=True)
    shutil.copy(os.path.join(staging, "preferred"),
                os.path.join(workdir, "preferred"))
    if gen_select == "list":
        shutil.copy(os.path.join(staging, "gen_target"),
                    os.path.join(workdir, "target"))
    tmpl = JOBS_TMPL_LIST if gen_select == "list" else JOBS_TMPL
    with open(os.path.join(workdir, "jobs.yaml"), "w") as fh:
        fh.write(tmpl.format(
            fam=fam, secs=",".join(map(str, box["sectors"])),
            r=box["r"], s=box["s"], d=box["d"]))
    t0 = time.time()
    with open(os.path.join(workdir, f"kira_{tag}.log"), "w") as lg:
        rc = subprocess.call(["nice", "-n", "10", "kira", "--parallel=4",
                              "jobs.yaml"], cwd=workdir, stdout=lg,
                             stderr=subprocess.STDOUT, timeout=900)
    wall = round(time.time() - t0, 2)
    assert rc == 0, f"kira rc={rc} in {workdir} (see kira_{tag}.log)"
    return wall


# ------------------------------------------------------ oracle table parsing
import re                                                       # noqa: E402

_ROW_HEAD = re.compile(r"^\s*(\w+)\[([-\d, ]+)\]\s*->\s*(.*)$")
_TERM = re.compile(r"([+-])?\s*(\w+)\[([-\d, ]+)\](?:\*\((.*)\))?\s*,?\s*$")


def parse_oracle_table(path):
    """kira_target.m -> list of (target_idx, [(sign, master_idx, expr), ...])."""
    rows, cur = [], None
    for raw in open(path):
        line = raw.strip()
        if line in ("{", "}", ""):
            continue
        if line == ",":
            cur = None
            continue
        m = _ROW_HEAD.match(line)
        if m:
            tgt = tuple(int(x) for x in m.group(2).split(","))
            cur = (tgt, [])
            rows.append(cur)
            rest = m.group(3).strip().rstrip(",")
            if rest in ("", "+", "0"):
                continue
            line = rest
        if cur is None:
            raise ValueError(f"orphan line in {path}: {raw!r}")
        body = line.rstrip(",").strip()
        if body in ("0", "+ 0", "+0", "- 0", "-0") or not body:
            continue
        t = _TERM.match(body)
        if not t:
            raise ValueError(f"unparsed term in {path}: {raw!r}")
        sign = -1 if t.group(1) == "-" else 1
        midx = tuple(int(x) for x in t.group(3).split(","))
        cur[1].append((sign, midx, t.group(4) or "1"))
    return rows


def planted_fault_trials_xcheck(banked, order, p, wdict, kira_rows, masters,
                                n_trials=3, seed=20260706):
    """Planted-fault trials for s=0 boxes (no in-box IBP seed
    exists, so the identity gate is structurally N/A). Detector = the
    INDEPENDENT kira full-table: a fault planted in a banked row must move
    some matched closed row off its kira reduction vector. Scratch copies
    only — never the originals."""
    import random
    rng = random.Random(seed)
    matched = sorted((t, w) for t, w in wdict.items() if w in banked)
    if not matched:
        return {"trials": [], "detection_rate": None, "all_detected": False,
                "note": "no matched banked rows to mutate"}
    trials = []
    for tno in range(n_trials):
        _, u = matched[rng.randrange(len(matched))]
        mutated = {k: dict(v) for k, v in banked.items()}   # scratch copy
        row = mutated[u]
        if row:
            c = rng.choice(sorted(row))
            row[c] = (row[c] + 1 + rng.randrange(p - 2)) % p or 1
        else:
            mc = rng.choice(sorted(masters))
            row[mc] = 1 + rng.randrange(p - 1)
        closed_m = assemble_closed(mutated, order, p)
        detected = any(closed_m.get(w, {}) != kira_rows.get(t, {})
                       for t, w in matched)
        trials.append({"trial": tno, "mutated_pivot": u, "detected": detected})
    return {"detector": "fulltable_kira_crosscheck", "trials": trials,
            "detection_rate": (sum(tr["detected"] for tr in trials)
                               / max(1, len(trials))),
            "all_detected": all(tr["detected"] for tr in trials)}


def idx_sector(idx):
    return sum(1 << i for i, a in enumerate(idx) if a > 0)


def idx_class(idx):
    return sum(a - 1 for a in idx if a > 0), -sum(a for a in idx if a < 0)


def eval_vec(terms, ev, midx2w, p):
    """(sign, midx, expr) list -> NEGATED F_p vector over master weights
    (matching the subs convention pivot + rest = 0), or ('denhit'|'nomaster').
    Masters keyed by full index TUPLE (multi-master sectors exist)."""
    vec = {}
    for sign, midx, expr in terms:
        try:
            v = ev(expr)
        except ZeroDivisionError:
            return "denhit"
        mw = midx2w.get(midx)
        if mw is None:
            return "nomaster"
        vec[mw] = (vec.get(mw, 0) - sign * v) % p
    return {k: v for k, v in vec.items() if v}


def _eval_vec_multi(terms, mev, midx2w, ps):
    """Streamed-mode analogue of eval_vec at ALL slices in one term pass,
    preserving the batch reference's PER-SLICE first-failure semantics
    (eval before master lookup within each term; a denominator hit censors
    only its own slice; a missing master censors every slice still alive).
    -> list per slice: vec dict | 'denhit' | 'nomaster'."""
    n = len(ps)
    vecs = [{} for _ in range(n)]
    state = [None] * n
    active = list(range(n))
    for sign, midx, expr in terms:
        if not active:
            break
        vals = mev.eval_some(expr, active)
        nxt = []
        for si in active:
            if vals[si] is None:
                state[si] = "denhit"
            else:
                nxt.append(si)
        active = nxt
        mw = midx2w.get(midx)
        if mw is None:
            for si in active:
                state[si] = "nomaster"
            active = []
            break
        for si in active:
            vecs[si][mw] = (vecs[si].get(mw, 0) - sign * vals[si]) % ps[si]
    return [state[si] if state[si] is not None
            else {k: v for k, v in vecs[si].items() if v}
            for si in range(n)]


def oracle_compare(mode, cdir, slices, closed_by_slice, kira_rows_by_slice,
                   midx2w, trivial, fam_masters_idx, class_of, sector_of,
                   shell_cols, ledger_path, fam, wdict):
    """Stage-D corpus oracle comparison.

    mode='batch': the legacy reference — materializes each
    kira_target.m per slice (parse_oracle_table + per-slice CoeffEvaluator +
    linear closed-row scan). Measured pathology: killed at 133 min /
    +112 GB RSS on the 343 MB sec1015 stage-1 table. Retained for
    regression only.

    mode='stream': single pass per table (stream_rows), one compile per
    distinct expr across all slices (MultiEvalSome, bounded cache),
    signature-INDEXED anchor matching, shell-diagnosis class set
    precomputed. Ledger lines are buffered and emitted slice-major in the
    batch order, so verdicts AND the ledger are byte-identical to batch
    (regression-gated). Duplicate targets within
    one instance table would make batch heldout statuses order-dependent —
    stream asserts they are absent (verified file-only on the banked
    reference tables).

    Returns (inst_verdicts, dict_stats)."""
    instances = sorted(os.listdir(cdir)) if os.path.isdir(cdir) else []
    anchor_map = {}
    dict_stats = {"anchored": 0, "ambiguous": 0, "class_mismatch": 0,
                  "unmatched": 0, "master_identity": 0, "zero_normalized": 0,
                  "anchor_vs_dictionary_disagree": 0, "shell_blocked": 0,
                  "dictrow_verified": 0, "dictrow_mismatch": 0}
    SHELL_SENTINEL = object()
    DICTROW_SENTINEL = object()
    inst_verdicts = {h: {"rows": {}} for h in instances}

    if mode == "batch":
        ledger = open(ledger_path, "a")
        for si, (p, d0, e0) in enumerate(slices):
            ev = CoeffEvaluator(p, d0, e0)
            closed = closed_by_slice[si]
            for h in instances:
                idir = os.path.join(cdir, h)
                iv = inst_verdicts[h]
                if si == 0:
                    imast = [i for i, _ in parse_masters_file(
                        os.path.join(idir, "masters"))]
                    iv["masters_subset_ok"] = set(imast) <= fam_masters_idx
                    iv["n_instance_masters"] = len(imast)
                ktm = os.path.join(idir, "kira_target.m")
                orows = parse_oracle_table(ktm) if os.path.exists(ktm) else []
                iv.setdefault("n_oracle_rows", len(orows))
                for tgt, terms in orows:
                    key = str(tgt)
                    vec = eval_vec(terms, ev, midx2w, p)
                    if vec == "denhit":
                        status = "denominator_hit_reported"
                    elif vec == "nomaster":
                        status = "MASTER_NOT_IN_RUN"
                    elif idx_sector(tgt) in trivial and not vec:
                        status = "zero_normalized"
                        dict_stats["zero_normalized"] += si == 0
                    elif tgt in fam_masters_idx:
                        status = ("match_master_identity"
                                  if midx2w.get(tgt) is not None
                                  and vec == {midx2w[tgt]: p - 1}
                                  else "MISMATCH")
                        dict_stats["master_identity"] += si == 0
                    elif si == 0:
                        cands = [w for w, r in closed.items() if r == vec]
                        if not cands:
                            kvec0 = kira_rows_by_slice[si].get(tgt)
                            cls = idx_class(tgt)
                            t_tgt = sum(1 for a in tgt if a > 0)
                            if kvec0 is not None:
                                if kvec0 == vec:
                                    status = "dictrow_verified"
                                    dict_stats["dictrow_verified"] += 1
                                    anchor_map[tgt] = DICTROW_SENTINEL
                                else:
                                    status = "MISMATCH"
                                    dict_stats["dictrow_mismatch"] += 1
                            else:
                                shellish = [
                                    w for w, r in closed.items()
                                    if class_of.get(w) == cls
                                    and popcount(sector_of[w]) == t_tgt
                                    and any(c in shell_cols for c in r)]
                                if shellish:
                                    status = "shell_blocked_reported"
                                    dict_stats["shell_blocked"] += 1
                                    anchor_map[tgt] = SHELL_SENTINEL
                                else:
                                    status = "UNMATCHED_AT_ANCHOR"
                                    dict_stats["unmatched"] += 1
                        else:
                            cls = idx_class(tgt)
                            good = [w for w in cands
                                    if class_of.get(w) == cls]
                            wsel = min(good) if good else min(cands)
                            anchor_map[tgt] = wsel
                            dict_stats["anchored"] += 1
                            dict_stats["ambiguous"] += len(cands) > 1
                            dict_stats["class_mismatch"] += not good
                            dict_stats["anchor_vs_dictionary_disagree"] += (
                                wdict.get(tgt) is not None
                                and closed.get(wdict[tgt]) != vec)
                            status = "anchored"
                    else:
                        wsel = anchor_map.get(tgt)
                        if wsel is DICTROW_SENTINEL:
                            kv = kira_rows_by_slice[si].get(tgt)
                            status = ("match_dictrow_heldout" if kv == vec
                                      and kv is not None else "MISMATCH")
                        else:
                            status = ("shell_blocked_reported"
                                      if wsel is SHELL_SENTINEL else
                                      "NO_ANCHOR" if wsel is None else
                                      "match_heldout"
                                      if closed.get(wsel) == vec
                                      else "MISMATCH")
                    iv["rows"].setdefault(key, []).append(status)
                    ledger.write(json.dumps(
                        {"family": fam, "instance": h, "target": key,
                         "slice": [p, d0, e0], "status": status}) + "\n")
            ledger.flush()
        ledger.close()
        return inst_verdicts, dict_stats

    # ---- stream mode: one pass per table, slice-major ledger emission
    assert mode == "stream", mode
    ps = [s[0] for s in slices]
    mev = MultiEvalSome(slices)
    closed0 = closed_by_slice[0]
    sig_index = {}
    for w, r in closed0.items():
        sig_index.setdefault(tuple(sorted(r.items())), []).append(w)
    shell_touch_cls = {(popcount(sector_of[w]), class_of.get(w))
                       for w, r in closed0.items()
                       if any(c in shell_cols for c in r)}
    buf = [[] for _ in slices]
    for h in instances:
        idir = os.path.join(cdir, h)
        iv = inst_verdicts[h]
        imast = [i for i, _ in parse_masters_file(
            os.path.join(idir, "masters"))]
        iv["masters_subset_ok"] = set(imast) <= fam_masters_idx
        iv["n_instance_masters"] = len(imast)
        ktm = os.path.join(idir, "kira_target.m")
        n_rows_h = 0
        seen_h = set()
        for tgt, terms in (stream_rows(ktm)
                           if os.path.exists(ktm) else ()):
            n_rows_h += 1
            assert tgt not in seen_h, \
                f"duplicate target {tgt} in {ktm} (order-dependent in batch)"
            seen_h.add(tgt)
            key = str(tgt)
            vres = _eval_vec_multi(terms, mev, midx2w, ps)
            for si, (p, d0, e0) in enumerate(slices):
                vec = vres[si]
                if vec == "denhit":
                    status = "denominator_hit_reported"
                elif vec == "nomaster":
                    status = "MASTER_NOT_IN_RUN"
                elif idx_sector(tgt) in trivial and not vec:
                    status = "zero_normalized"
                    dict_stats["zero_normalized"] += si == 0
                elif tgt in fam_masters_idx:
                    status = ("match_master_identity"
                              if midx2w.get(tgt) is not None
                              and vec == {midx2w[tgt]: p - 1}
                              else "MISMATCH")
                    dict_stats["master_identity"] += si == 0
                elif si == 0:
                    cands = sig_index.get(tuple(sorted(vec.items())), [])
                    if not cands:
                        kvec0 = kira_rows_by_slice[0].get(tgt)
                        cls = idx_class(tgt)
                        t_tgt = sum(1 for a in tgt if a > 0)
                        if kvec0 is not None:
                            if kvec0 == vec:
                                status = "dictrow_verified"
                                dict_stats["dictrow_verified"] += 1
                                anchor_map[tgt] = DICTROW_SENTINEL
                            else:
                                status = "MISMATCH"
                                dict_stats["dictrow_mismatch"] += 1
                        elif (t_tgt, cls) in shell_touch_cls:
                            status = "shell_blocked_reported"
                            dict_stats["shell_blocked"] += 1
                            anchor_map[tgt] = SHELL_SENTINEL
                        else:
                            status = "UNMATCHED_AT_ANCHOR"
                            dict_stats["unmatched"] += 1
                    else:
                        cls = idx_class(tgt)
                        good = [w for w in cands
                                if class_of.get(w) == cls]
                        wsel = min(good) if good else min(cands)
                        anchor_map[tgt] = wsel
                        dict_stats["anchored"] += 1
                        dict_stats["ambiguous"] += len(cands) > 1
                        dict_stats["class_mismatch"] += not good
                        dict_stats["anchor_vs_dictionary_disagree"] += (
                            wdict.get(tgt) is not None
                            and closed0.get(wdict[tgt]) != vec)
                        status = "anchored"
                else:
                    wsel = anchor_map.get(tgt)
                    if wsel is DICTROW_SENTINEL:
                        kv = kira_rows_by_slice[si].get(tgt)
                        status = ("match_dictrow_heldout" if kv == vec
                                  and kv is not None else "MISMATCH")
                    else:
                        status = ("shell_blocked_reported"
                                  if wsel is SHELL_SENTINEL else
                                  "NO_ANCHOR" if wsel is None else
                                  "match_heldout"
                                  if closed_by_slice[si].get(wsel) == vec
                                  else "MISMATCH")
                iv["rows"].setdefault(key, []).append(status)
                buf[si].append(json.dumps(
                    {"family": fam, "instance": h, "target": key,
                     "slice": [p, d0, e0], "status": status}) + "\n")
        iv.setdefault("n_oracle_rows", n_rows_h)
    ledger = open(ledger_path, "a")
    for si in range(len(slices)):
        for line in buf[si]:
            ledger.write(line)
        ledger.flush()
    ledger.close()
    return inst_verdicts, dict_stats


# ----------------------------------------------------------------- benchmark
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--staging", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--slices", required=True,
                    help="comma list p:d:eta; FIRST is the anchor slice")
    ap.add_argument("--reuse-artifacts", default=None)
    ap.add_argument("--b2b3", action="store_true")
    ap.add_argument("--c-min", type=int, default=2)
    ap.add_argument("--dict-timeout", type=int, default=1800)
    ap.add_argument("--skip-dict", action="store_true",
                    help="pre-declared degradation 3a: skip the kira "
                    "dictionary/cross-check stage (oracle-only "
                    "certification, ships PARTIAL)")
    ap.add_argument("--reuse-dict", default=None,
                    help="existing SYMBOLIC dict workdir (skip the kira dict "
                    "run, reuse its kira_target.m; consumed streaming unless "
                    "--dict-mode batch)")
    ap.add_argument("--dict-mode", default="numeric-slice",
                    choices=["numeric-slice", "stream", "batch"],
                    help="dict-table cost mode: numeric-slice = one eta-pinned "
                    "kira run per distinct slice eta, streamed (default; "
                    "measured 1493 s vs 2107 s stream vs ~2800 s/17.3 GB "
                    "batch on lp1disp); stream = one symbolic (d,eta) run, "
                    "single-pass streaming eval (use for >~3 eta points); "
                    "batch = legacy materializing path (regression only)")
    ap.add_argument("--dict-scope", default="span",
                    choices=["span", "core", "targets"],
                    help="dict target-list scope: span = +1-shell "
                    "box (r+1,s+1, full dots window; legacy default — "
                    "measured 1.16M tuples on a production family, prohibitive); "
                    "core = the solve box itself (widened dots window; "
                    "covers every in-box lead incl. all corpus targets; "
                    "shell-touched rows are excluded from the crosscheck "
                    "scope and counted, per the shell accounting); "
                    "targets = ONLY the corpus instances' oracle tuples "
                    "(mini-dict: independent kira route on exactly the "
                    "oracle rows; crosscheck ships PARTIAL)")
    ap.add_argument("--gen-select", default="recursive",
                    choices=["recursive", "list"],
                    help="Generate selection mode: recursive = box-recursive "
                    "(default); list = select_mandatory_list of the staging "
                    "'gen_target' file (provenance-faithful reproduction of "
                    "production pins, e.g. sec1015)")
    ap.add_argument("--oracle-mode", default="stream",
                    choices=["stream", "batch"],
                    help="corpus-oracle comparison mode: "
                    "stream = single-pass per-table streaming eval, bounded "
                    "RSS, one compile per distinct expr, signature-indexed "
                    "anchor matching (default; regression-gated "
                    "byte-identical vs batch on the banked sh + "
                    "sec1015-prefix audits); batch = the legacy materializing "
                    "reference (measured: killed at 133 min / +112 GB "
                    "RSS on the 343 MB sec1015 stage-1 table)")
    ap.add_argument("--dict-staging", default="cell",
                    choices=["cell", "persector"],
                    help="numeric-slice dict staging: "
                    "cell = one kira run per eta on the whole cell "
                    "(reference; measured cell-reduce-bound on sh) ; "
                    "persector = one kira leg per target sector, r-budget "
                    "per leg (see --dict-rbudget), smallest sectors first")
    ap.add_argument("--dict-rbudget", default="trim",
                    choices=["trim", "box"],
                    help="persector leg r-budget: trim = max target r of "
                    "the leg (regression-gated vs banked micro-probe rows); "
                    "box = the full dict-box r (fallback)")
    ap.add_argument("--prereg-shape", default=None)
    ap.add_argument("--prereg-front", type=int, default=None)
    ap.add_argument("--coverage-mode", default="auto",
                    choices=["auto", "xcheck"],
                    help="auto = identity-coverage gate when in-box IBP "
                    "seeds exist (small/medium rungs); xcheck = FORCE the "
                    "composite full-table-crosscheck path (at "
                    "flagship scale the banked-dependency closure of the "
                    "identity gate is unpriced — forcing xcheck is a "
                    "pre-declared, reported scope statement, not a pass)")
    args = ap.parse_args()

    fam = args.family
    os.makedirs(args.out, exist_ok=True)
    box = json.load(open(os.path.join(args.staging, "box.json")))
    slices = [tuple(int(x) for x in tok.split(":"))
              for tok in args.slices.split(",")]
    report = {"family": fam, "box": box, "slices": slices,
              "grid_stamp": "NOT-VALID-OFF-GRID: values defined only at the "
                            "listed (p,d,eta) slices", "stages": {}}
    t_all = time.time()

    # ---------------- A. solve system + per-slice stratified solve
    if args.reuse_artifacts:
        art = args.reuse_artifacts
        report["stages"]["kira_solve"] = {"reused": art}
    else:
        art = os.path.join(args.workdir, "solve")
        report["stages"]["kira_solve"] = {
            "wall_s": run_kira(art, fam, args.staging, box["solve"], "solve",
                               gen_select=args.gen_select),
            "gen_select": args.gen_select}

    closed_by_slice, banked_by_slice, per_slice = [], [], []
    sysd0 = rows0 = None
    for si, (p, d0, e0) in enumerate(slices):
        t0 = time.time()
        sysd = load_system(art, fam, p, d0, e0)
        rows, order = sysd["rows"], sysd["order"]
        stratum_of, masters = sysd["stratum_of"], sysd["masters"]
        if si == 0:
            sysd0, rows0 = sysd, rows
            shape = {}
            for r in rows:
                lead = max((c for c in r if c in order), default=None,
                           key=lambda c: order.get(c, -1))
                if lead is not None:
                    shape[stratum_of[lead]] = shape.get(stratum_of[lead], 0) + 1
            report["measured_selected_by_t"] = shape
            report["loader_counters"] = sysd["counters"]
        tbl = InterfaceTable(os.path.join(args.out, f"T_{fam}"),
                             meta=None if si else {"family": fam})
        node = f"p{p}_d{d0}_e{e0}"
        subs_all, leftovers, stats_per = stratified_solve(
            rows, order, p, stratum_of,
            policy=os.environ.get("STRATA_POLICY", "B2"),
            forbid=set(masters), table=tbl, cell_key="tiny", node_key=node)
        banked = {}
        for t in sorted(set(stratum_of.values()), reverse=True):
            kn = f"{node}/t{t}"
            if tbl.has("tiny", kn, family=fam):
                for lead, terms in tbl.get("tiny", kn, family=fam):
                    banked[lead] = dict(terms)
        assert set(banked) == set(subs_all), "T round-trip lost pivots"
        closed = assemble_closed(banked, order, p)
        closed_by_slice.append(closed)
        banked_by_slice.append(banked)
        if si == 0:
            leftovers0 = [dict(r) for r in leftovers]
        per_slice.append({
            "slice": [p, d0, e0], "n_pivots": len(closed),
            "n_leftover_master_relations": len(leftovers),
            "stats_per_stratum": {str(k): v for k, v in stats_per.items()},
            "solve_wall_s": round(time.time() - t0, 2)})
        if si == 0 and args.b2b3:
            shadow = {w: sum(c) for w, c in sysd["class_of"].items() if c}
            per_slice[0]["b2b3"] = b2b3_compare(
                rows, order, p, shadow, set(masters), tag=f"{fam}_tiny_real",
                out_jsonl=os.path.join(args.out, "b2b3_compare.jsonl"))

    # -------- survivor-master labeling (non-preferred masters have ordinary
    # weights and can never be pivoted on, so they SURVIVE elimination)
    masters = dict(sysd0["masters"])            # preferred/ordinal masters
    order = sysd0["order"]
    surv = set()
    for closed in closed_by_slice:
        surv |= {c for r in closed.values() for c in r if c not in masters}
    surv |= {c for r in leftovers0 for c in r if c not in masters}
    # SHELL columns (kira dump-closure gap, measured): mandatory-
    # recursive SYSTEM dumps carry terms one shell beyond the select box
    # whose closing equations kira only regenerates internally (never
    # dumped; escalating the box just moves the shell). They are NOT
    # masters; rows referencing them are reductions modulo shell and are
    # certified only against the corpus oracle, not the full-table check.
    sb0 = box["solve"]
    shell_cols = set()
    for w in list(surv):
        cls = sysd0["class_of"].get(w)
        t_w = popcount(sysd0["sector_of"][w])
        if cls and (cls[0] > sb0["d"] or cls[1] > sb0["s"]
                    or t_w + cls[0] > sb0["r"]):
            shell_cols.add(w)
            surv.discard(w)
    used_idx = {m["indices"] for m in masters.values() if m["indices"]}
    buckets = {}
    for idx, sec in sysd0["file_masters"]:
        if idx in used_idx:
            continue
        buckets.setdefault((sec, idx_class(idx)), []).append(idx)
    label_ambig = label_missing = 0
    for w in sorted(surv):
        key = (sysd0["sector_of"][w], sysd0["class_of"].get(w))
        b = buckets.get(key, [])
        if not b:
            label_missing += 1
            masters[w] = {"sector": key[0], "indices": None}
        else:
            if len(b) > 1:
                label_ambig += 1
            masters[w] = {"sector": key[0], "indices": b.pop(0)}
    file_unmatched = sum(len(b) for b in buckets.values())
    shell_touched0 = {w for w, r in closed_by_slice[0].items()
                      if any(c in shell_cols for c in r)}
    report["shell"] = {
        "n_shell_cols": len(shell_cols),
        "n_rows_shell_touched": len(shell_touched0),
        "shell_col_classes": sorted(
            {(sysd0["sector_of"][w],) + tuple(sysd0["class_of"][w] or ())
             for w in shell_cols})}
    report["master_census"] = {
        "masters_file_n": len(sysd0["file_masters"]),
        "ordinal_masters": len(sysd0["masters"]),
        "survivor_masters": len(surv),
        "label_ambiguous": label_ambig, "label_missing": label_missing,
        "file_masters_unmatched": file_unmatched,
        "match": (label_missing == 0 and file_unmatched == 0
                  and len(sysd0["masters"]) + len(surv)
                  == len(sysd0["file_masters"]))}
    report["fam_masters_idx"] = sorted(
        [list(i) for i, _ in sysd0["file_masters"]])
    midx2w = {m["indices"]: w for w, m in masters.items()
              if m["indices"] is not None}
    for si in range(len(per_slice)):
        per_slice[si]["n_unsolved_nonmaster_cols"] = len(
            {c for r in closed_by_slice[si].values() for c in r}
            - set(masters))
    try:
        spec = parse_family_spec(
            open(os.path.join(args.staging, "config",
                              "integralfamilies.yaml")).read(),
            open(os.path.join(args.staging, "config",
                              "kinematics.yaml")).read())
        n_idx = len(spec["u"])
    except Exception as exc:  # noqa: BLE001 — spec is only needed for the
        # identity generator, which is structurally N/A on s=0 boxes;
        # non-canonical propagator forms (lp1disp composite
        # "(3*l^2+...)/(3)") raise Assertion/Syntax/Attribute errors here.
        spec = None
        n_idx = len(sysd0["file_masters"][0][0])
        report["spec_parse_error"] = f"{type(exc).__name__}: {exc}"[:200]
    sectors_nontriv = sorted({sysd0["sector_of"][w] for w in sysd0["order"]
                              if w not in masters}
                             | {m["sector"] for m in masters.values()})

    # ---------------- B. dictionary run + full-table cross-check
    # +1-shell box: selected equations may carry terms one IBP step outside
    # the seed box (r+1 / s+1 / d+1), and those columns are solved too —
    # the dictionary must cover the whole SPAN, not just the box
    # (vac3_b0_r0_f0 finding).
    # d=None (dots-window fix, measured falsification x3):
    # dots enumerate to the full r-budget rmax - t per sector. The old cap
    # solve.d+2 missed deep-dot pivots of the recursive selection (dots 6-7
    # at t=2 on k2disp: 232 clean rows class-absent from the kira table).
    dict_dir = os.path.join(args.workdir, "dict")
    # --dict-scope: span (legacy +1-shell box), core (solve box,
    # widened dots), targets (corpus oracle tuples only, mini-dict).
    if args.dict_scope == "core":
        dbox = {"sectors": box["solve"]["sectors"], "r": box["solve"]["r"],
                "s": box["solve"]["s"], "d": None}
    else:
        dbox = {"sectors": box["solve"]["sectors"],
                "r": box["solve"]["r"] + 1,
                "s": box["solve"]["s"] + 1, "d": None}
    dict_target_tuples = None
    if args.dict_scope == "targets":
        cdir0 = os.path.join(args.staging, "corpus")
        tset = set()
        for h in (sorted(os.listdir(cdir0))
                  if os.path.isdir(cdir0) else []):
            ktm0 = os.path.join(cdir0, h, "kira_target.m")
            if os.path.exists(ktm0):
                for tgt0, _terms0 in parse_oracle_table(ktm0):
                    tset.add(tgt0)
        dict_target_tuples = sorted(tset)
    try:
        if args.skip_dict:
            raise RuntimeError("dict stage skipped by --skip-dict "
                               "(pre-declared degradation 3a)")
        n_dict_eval_skipped = 0
        dict_mode = args.dict_mode
        if dict_mode == "numeric-slice" and not config_has_eta(args.staging):
            # eta-free family: numeric-slice degenerates to the single
            # symbolic run consumed streaming.
            dict_mode = "stream"
        if args.dict_scope == "targets" and dict_mode != "numeric-slice":
            raise RuntimeError("--dict-scope targets requires numeric-slice")
        # numeric-slice REUSE — a previous leg's per-eta pinned
        # tables (dict_eta*/DICT_ETA.json) are d-symbolic and reusable at
        # any d-slices with the same eta set.
        reuse_numeric = None
        if args.reuse_dict and dict_mode == "numeric-slice":
            import glob as _glob
            stamps = sorted(
                _glob.glob(os.path.join(args.reuse_dict, "dict_eta*",
                                        "DICT_ETA.json"))
                + _glob.glob(os.path.join(args.reuse_dict, "dict_eta*",
                                          "leg_sec*", "DICT_ETA.json")))
            if stamps:
                # A stamp dir may be a single cell table OR a set of
                # persector leg dirs (leg_sec*/DICT_ETA.json) — one eta maps
                # to a LIST of table dirs, merged (targets disjoint by
                # construction: legs partition tuples by sector).
                reuse_numeric = {}
                for sp in stamps:
                    st = json.load(open(sp))
                    reuse_numeric.setdefault(st["eta"], []).append(
                        os.path.dirname(sp))
        if dict_mode == "numeric-slice" and reuse_numeric is not None:
            kira_rows_by_slice = [None] * len(slices)
            n_targets = None
            for si, (_p, _d0, e0) in enumerate(slices):
                if e0 not in reuse_numeric:
                    raise RuntimeError(
                        f"reuse-dict numeric-slice: no table for eta {e0}")
            groups = {}
            for si, (_p, _d0, e0) in enumerate(slices):
                groups.setdefault(e0, []).append(si)
            for e0, sidx in groups.items():
                krs_m, n_targets = None, 0
                for wd in reuse_numeric[e0]:
                    tab = kira_table_path(wd, fam)
                    st = json.load(open(os.path.join(wd, "DICT_ETA.json")))
                    n_targets += st["n_targets"]
                    krs, _n_rows, n_skip = stream_eval_table(
                        tab, [slices[i] for i in sidx], midx2w)
                    n_dict_eval_skipped += n_skip
                    if krs_m is None:
                        krs_m = krs
                    else:
                        for k in range(len(krs)):
                            dup = set(krs_m[k]) & set(krs[k])
                            assert not dup, f"persector legs overlap: {dup}"
                            krs_m[k].update(krs[k])
                for k, si in enumerate(sidx):
                    kira_rows_by_slice[si] = krs_m[k]
            common = set(kira_rows_by_slice[0])
            for kr in kira_rows_by_slice[1:]:
                common &= set(kr)
            for kr in kira_rows_by_slice:
                for tgt in set(kr) - common:
                    kr.pop(tgt)
            report["stages"]["kira_dict"] = {
                "mode": "numeric-slice", "reused": args.reuse_dict,
                "n_targets": n_targets,
                "n_eval_skipped": n_dict_eval_skipped}
        elif dict_mode == "numeric-slice" and not args.reuse_dict:
            # numeric-slice default: one dict kira
            # per DISTINCT slice eta (numkin pattern, Fermat backend),
            # each table consumed streaming at its own slices.
            groups = []
            for si, (_p, _d0, e0) in enumerate(slices):
                for g in groups:
                    if g["eta"] == e0:
                        g["sidx"].append(si)
                        break
                else:
                    groups.append({"eta": e0, "sidx": [si]})
            kira_rows_by_slice = [None] * len(slices)
            walls, sizes, n_targets = [], [], None
            for gi, g in enumerate(groups):
                wd = os.path.join(args.workdir, f"dict_eta{gi}")
                if args.dict_staging == "persector":
                    tuples_g = (dict_target_tuples
                                if dict_target_tuples is not None else
                                enumerate_box(sectors_nontriv, n_idx,
                                              dbox["r"], dbox["s"],
                                              dbox.get("d")))
                    legs = run_dict_kira_numeric_persector(
                        wd, fam, args.staging, n_idx, g["eta"], tuples_g,
                        box_r=dbox["r"], box_s=dbox["s"],
                        timeout_per_leg=args.dict_timeout,
                        r_budget=args.dict_rbudget)
                    walls.append(round(sum(l["wall_s"] for l in legs), 2))
                    n_targets = sum(l["n_targets"] for l in legs)
                    sizes.append(sum(os.path.getsize(l["table"])
                                     for l in legs))
                    report["stages"].setdefault("kira_dict_legs", {})[
                        f"eta{gi}"] = [
                        {k: l[k] for k in ("sector", "t", "r", "s",
                                           "n_targets", "wall_s")}
                        for l in legs]
                    krs_m = None
                    for l in legs:
                        krs, _n_rows, n_skip = stream_eval_table(
                            l["table"], [slices[i] for i in g["sidx"]],
                            midx2w)
                        n_dict_eval_skipped += n_skip
                        if krs_m is None:
                            krs_m = krs
                        else:
                            for k in range(len(krs)):
                                dup = set(krs_m[k]) & set(krs[k])
                                assert not dup, \
                                    f"persector legs overlap: {dup}"
                                krs_m[k].update(krs[k])
                    for k, si in enumerate(g["sidx"]):
                        kira_rows_by_slice[si] = krs_m[k]
                    continue
                wall, n_targets = run_dict_kira_numeric(
                    wd, fam, args.staging, dbox, sectors_nontriv, n_idx,
                    g["eta"], timeout=args.dict_timeout,
                    target_tuples=dict_target_tuples)
                walls.append(wall)
                tab = kira_table_path(wd, fam)
                sizes.append(os.path.getsize(tab))
                krs, _n_rows, n_skip = stream_eval_table(
                    tab, [slices[i] for i in g["sidx"]], midx2w)
                n_dict_eval_skipped += n_skip
                for k, si in enumerate(g["sidx"]):
                    kira_rows_by_slice[si] = krs[k]
            common = set(kira_rows_by_slice[0])
            for kr in kira_rows_by_slice[1:]:
                common &= set(kr)
            for kr in kira_rows_by_slice:  # cross-table skip consistency
                for tgt in set(kr) - common:
                    kr.pop(tgt)
            report["stages"]["kira_dict"] = {
                "mode": "numeric-slice", "wall_s": round(sum(walls), 2),
                "wall_s_per_eta": walls, "n_targets": n_targets,
                "table_bytes_per_eta": sizes,
                "scope": args.dict_scope,
                "n_eval_skipped": n_dict_eval_skipped}
        else:
            if args.reuse_dict:
                dict_dir = args.reuse_dict
                report["stages"]["kira_dict"] = {"reused": dict_dir,
                                                 "mode": dict_mode}
            else:
                wall, n_targets = run_dict_kira(dict_dir, fam, args.staging,
                                                dbox, sectors_nontriv, n_idx,
                                                timeout=args.dict_timeout)
                report["stages"]["kira_dict"] = {"wall_s": wall,
                                                 "n_targets": n_targets,
                                                 "mode": dict_mode}
            tab = kira_table_path(dict_dir, fam)
            if dict_mode == "batch":        # legacy materializing reference
                ktable = parse_oracle_table(tab)
                kira_rows_by_slice = []
                skipped_tuples = set()
                for (p, d0, e0) in slices:
                    ev = CoeffEvaluator(p, d0, e0)
                    kr = {}
                    for tgt, terms in ktable:
                        v = eval_vec(terms, ev, midx2w, p)
                        if v in ("denhit", "nomaster"):
                            n_dict_eval_skipped += 1
                            skipped_tuples.add(tgt)
                            continue
                        kr[tgt] = v
                    kira_rows_by_slice.append(kr)
                if skipped_tuples:  # drop at ALL slices or none
                    for kr in kira_rows_by_slice:
                        for tgt in skipped_tuples:
                            kr.pop(tgt, None)
            else:                           # stream: single-pass eval
                kira_rows_by_slice, _n_rows, n_dict_eval_skipped = \
                    stream_eval_table(tab, slices, midx2w)
            report["stages"]["kira_dict"]["n_eval_skipped"] = \
                n_dict_eval_skipped
        dict_ok = True
    except Exception as exc:                 # pre-declared degradation 3a:
        # dictionary/cross-check stage failed or timed out — continue with
        # oracle-only certification (anchor matching needs only midx2w);
        # crosscheck + planted-fault detector unavailable => certification
        # ships PARTIAL, never silently passed.
        report["stages"]["kira_dict"] = {"failed": str(exc)[:200]}
        dict_ok = False
    if not dict_ok:
        kira_rows_by_slice = [{} for _ in slices]
        wdict, dstats = {}, {"skipped_dict_stage": True}
        report["dictionary_run"] = dstats
        report["fulltable_crosscheck_pass"] = False
        report["certification_level"] = "oracle_only_dict_skipped"
        zero_closed = set()
    else:
        wdict, dstats = match_dictionary(kira_rows_by_slice, closed_by_slice,
                                         midx2w, [s[0] for s in slices])
        dstats["n_closed_rows"] = len(closed_by_slice[0])
        dstats["closed_rows_unmatched"] = (len(closed_by_slice[0])
                                           - dstats["matched"])
        # ZERO closed rows (b18_r1_f1 finding): a pivot whose row
        # is empty at ALL slices is a genuine zero reduction; zero vectors
        # carry no signature, so match_dictionary structurally cannot pair
        # them (its zero-tuple skip is correct). Certified as a THIRD
        # category: all-slice-zero (4 slices, 2 primes) + per-class count
        # vs kira's own in-solve-span zero tuples (must not exceed them).
        zero_closed = {w for w in closed_by_slice[0]
                       if all(not cs.get(w) for cs in closed_by_slice)}
        sb_solve = box["solve"]

        def in_solve_span(sec, cls):
            t = popcount(sec)
            return (cls is not None and cls[0] <= sb_solve["d"]
                    and cls[1] <= sb_solve["s"]
                    and t + cls[0] <= sb_solve["r"])

        zk_class = {}
        for tgt in kira_rows_by_slice[0]:
            if any(kr.get(tgt) for kr in kira_rows_by_slice):
                continue
            sec, cls = idx_sector(tgt), idx_class(tgt)
            if in_solve_span(sec, cls):
                zk_class[(sec, cls)] = zk_class.get((sec, cls), 0) + 1
        zc_class = {}
        for w in zero_closed:
            key = (sysd0["sector_of"][w], sysd0["class_of"].get(w))
            zc_class[key] = zc_class.get(key, 0) + 1
        zero_class_gate = all(zc_class[k] <= zk_class.get(k, 0)
                              for k in zc_class)
        dstats["zero_closed_rows"] = len(zero_closed)
        dstats["zero_class_gate_pass"] = zero_class_gate
        dstats["zero_class_exact_equality"] = (zc_class == {
            k: v for k, v in zk_class.items() if k in zc_class}
            and zero_class_gate)
        report["dictionary_run"] = dstats
        # DUPLICATE-CLASS accounting (k2disp finding): symmetric/
        # degenerate integrals give IDENTICAL reduction vectors (measured:
        # 7,135 rows -> 3,203 distinct vectors, classes up to 79). The 1:1
        # matcher can only consume as many rows per class as kira presents
        # tuples; the excess duplicates are nevertheless verified by the
        # SAME independent kira vector — certified here by explicit
        # all-slice equality with a matched representative.
        matched_ws = set(wdict.values())
        rep_of_sig = {}
        for w in matched_ws:
            row0 = closed_by_slice[0].get(w)
            if row0:
                rep_of_sig.setdefault(tuple(sorted(row0.items())), w)
        n_dup_verified = 0
        # SURVIVOR-blocked rows: support includes a box-edge survivor
        # master that the LARGER dict box reduces further — structurally
        # unmatchable (basis mismatch at the box edge), reported like shell.
        survivor_touched0 = {
            w for w, r in closed_by_slice[0].items()
            if w not in shell_touched0 and any(c in surv for c in r)}
        for w, row0 in closed_by_slice[0].items():
            if (w in matched_ws or w in shell_touched0 or w in zero_closed
                    or w in survivor_touched0 or not row0):
                continue
            w0 = rep_of_sig.get(tuple(sorted(row0.items())))
            if w0 is not None and all(
                    cs.get(w) == cs.get(w0) for cs in closed_by_slice):
                n_dup_verified += 1
        dstats["dup_of_matched_verified_all_slices"] = n_dup_verified
        dstats["survivor_blocked_rows"] = len(survivor_touched0)
        # pass = every CLEAN NONZERO Route-B row matched a kira row at
        # anchor AND survived held-out (directly or via its duplicate
        # class). Shell-touched rows cannot match by construction (kira
        # rows never contain shell columns) and zero rows cannot match by
        # construction (no signature) — together with verified duplicates
        # they must account exactly for the unmatched count; zero rows must
        # additionally pass the class-count gate vs kira's independent zeros.
        report["fulltable_crosscheck_pass"] = (
            dstats["heldout_mismatch"] == 0
            and dstats["closed_rows_unmatched"]
            == (len(shell_touched0 - zero_closed) + len(zero_closed)
                + n_dup_verified + len(survivor_touched0))
            and zero_class_gate)

    # ---------------- C. independent fresh identities + coverage gate
    trivial = load_trivial_sectors(art, fam)
    sb = box["solve"]
    seeds = (seed_box(sectors_nontriv, n_idx, sb["r"], sb["s"], sb["d"])
             if spec is not None and args.coverage_mode == "auto" else [])
    report["coverage_mode_flag"] = args.coverage_mode
    report["fresh_identities"] = {
        "n_seeds": len(seeds),
        "n_ops": (len(spec["loops"]) * len(spec["mom"]))
        if spec is not None else None,
        "spec_available": spec is not None}
    for w, m in masters.items():        # masters are legitimate identity terms
        if m["indices"]:
            wdict.setdefault(tuple(m["indices"]), w)
    inv_dict = {w: t for t, w in wdict.items()}
    smax_box = sb["s"]
    if not seeds:
        # s=0 solve box (measured finding): NO in-box IBP seed
        # exists — every identity image raises s to 1 and leaves the box, so
        # the whole box is the boundary shell and the identity-coverage gate
        # is STRUCTURALLY N/A. Certification for these families is the
        # composite: full-table kira cross-check (100% row coverage, all
        # slices, 2 primes) + corpus oracle + planted faults detected
        # against the independent kira table. Reported, never hidden.
        report["coverage_mode"] = "fulltable_xcheck_s0_box"
        for si in range(len(slices)):
            per_slice[si]["coverage_gate"] = {
                "structurally_na_s0_box": True,
                "n_identities_raw": 0,
                "coverage_source": "fulltable_kira_crosscheck",
                "pass": True,       # no identity claim is made
                "composite_pass": report["fulltable_crosscheck_pass"]}
        per_slice[0]["planted_faults"] = planted_fault_trials_xcheck(
            banked_by_slice[0], order, slices[0][0], wdict,
            kira_rows_by_slice[0], masters)
    else:
        report["coverage_mode"] = "identity_coverage"
    for si, (p, d0, e0) in enumerate(slices) if seeds else ():
        fresh, n_raw, n_unusable = [], 0, 0
        for _a, _i, _v, irow in identities(spec, seeds, p, d0, e0):
            n_raw += 1
            wrow, ok = {}, True
            for tup, cv in irow.items():
                if idx_sector(tup) in trivial:
                    continue
                w = wdict.get(tup)
                if w is None:
                    ok = False
                    break
                wrow[w] = (wrow.get(w, 0) + cv) % p
            if not ok:
                n_unusable += 1
                continue
            wrow = {k: v for k, v in wrow.items() if v}
            if wrow:
                fresh.append(wrow)
        gate, touch = coverage_gate(fresh, closed_by_slice[si],
                                    banked_by_slice[si], masters, p,
                                    c_min=args.c_min,
                                    excluded_cols=frozenset(shell_cols))
        gate["n_identities_raw"] = n_raw
        gate["n_identities_unusable"] = n_unusable
        # Boundary-shell split: in-box IBP identities cannot
        # reach the outermost s-shell with full multiplicity (the (0,smax)
        # corners are structurally unreachable: covering them requires
        # s=smax seeds whose identities leave the box). The s-shell is
        # therefore certified by the independent full-table cross-check,
        # the interior by the coverage gate. Composite, reported as such.
        uncov = [u for u in banked_by_slice[si]
                 if touch.get(u, 0) < args.c_min]
        bshell = sum(1 for u in uncov
                     if inv_dict.get(u) is not None
                     and -sum(a for a in inv_dict[u] if a < 0) == smax_box)
        n_sht = sum(1 for u in uncov
                    if u in shell_touched0 and u not in zero_closed)
        n_zc = sum(1 for u in uncov
                   if u in zero_closed and inv_dict.get(u) is None)
        gate[f"uncovered_boundary_shell_s{smax_box}"] = bshell
        gate["uncovered_shell_touched"] = n_sht
        gate["uncovered_zero_closed"] = n_zc      # certified by zero-class gate
        gate["uncovered_interior"] = len(uncov) - bshell - n_sht - n_zc
        gate["composite_pass"] = (gate["n_nonzero_residual"] == 0
                                  and gate["uncovered_interior"] == 0)
        per_slice[si]["coverage_gate"] = gate
        if si == 0:
            per_slice[0]["planted_faults"] = planted_fault_trials(
                fresh, banked_by_slice[0], order, masters, p,
                n_trials=3, c_min=args.c_min,
                excluded_cols=frozenset(shell_cols))

    # ---------------- D. corpus oracle comparison (oracle_compare: stream =
    # default, batch = the legacy materializing reference)
    cdir = os.path.join(args.staging, "corpus")
    instances = sorted(os.listdir(cdir)) if os.path.isdir(cdir) else []
    fam_masters_idx = {tuple(x) for x in map(tuple, report["fam_masters_idx"])}
    _orss0 = maxrss_kb()
    _ot0 = time.time()
    inst_verdicts, dict_stats = oracle_compare(
        args.oracle_mode, cdir, slices, closed_by_slice, kira_rows_by_slice,
        midx2w, trivial, fam_masters_idx, sysd0["class_of"],
        sysd0["sector_of"], shell_cols, args.ledger, fam, wdict)
    report["oracle_mode"] = args.oracle_mode
    # Record the oracle-audit stage wall + RSS increment (the
    # streamed-comparator record datum; zero behavior change).
    report["oracle_stage"] = {
        "wall_s": round(time.time() - _ot0, 2),
        "maxrss_incr_kb": maxrss_kb() - _orss0}

    report["per_slice"] = per_slice
    report["anchor_dictionary"] = dict_stats
    bad = ("MISMATCH", "MASTER_NOT_IN_RUN", "UNMATCHED_AT_ANCHOR", "NO_ANCHOR")
    n_rows = n_ok = n_shell_rows = n_dictroute = 0
    for h, iv in inst_verdicts.items():
        ok = iv.get("masters_subset_ok", False)
        n_dr_inst = 0
        for sts in iv["rows"].values():
            n_rows += 1
            if any(s == "shell_blocked_reported" for s in sts):
                n_shell_rows += 1       # unverifiable from the dumped
                continue                # system — reported, NOT a pass
            row_ok = not any(s in bad for s in sts)
            if row_ok and any(s in ("dictrow_verified",
                                    "match_dictrow_heldout") for s in sts):
                # Row verified by the INDEPENDENT dict route (kira on
                # the fresh-generated staging), not by a Route-B closed row.
                # Counted separately — a value verification, not a solve
                # coverage claim.
                n_dictroute += 1
                n_dr_inst += 1
            else:
                n_ok += row_ok
            ok = ok and row_ok
        iv["n_rows_dictroute_verified"] = n_dr_inst
        iv["instance_pass"] = ok
    report["instances"] = inst_verdicts
    gates_ok = all(sl["coverage_gate"]["pass"] for sl in per_slice)
    composite_ok = (all(sl["coverage_gate"]["composite_pass"]
                        for sl in per_slice)
                    and report["fulltable_crosscheck_pass"])
    faults_ok = per_slice[0]["planted_faults"]["all_detected"]
    report["summary"] = {
        "n_instances": len(instances),
        "n_instances_pass": sum(iv["instance_pass"]
                                for iv in inst_verdicts.values()),
        "n_oracle_rows": n_rows, "n_oracle_rows_ok": n_ok,
        "n_oracle_rows_dictroute_verified": n_dictroute,
        "n_oracle_rows_shell_blocked": n_shell_rows,
        "n_oracle_rows_verifiable": n_rows - n_shell_rows,
        "row_match_pct_verifiable": (
            100.0 * (n_ok + n_dictroute) / (n_rows - n_shell_rows))
        if n_rows - n_shell_rows else None,
        "VOID": ((n_rows - n_shell_rows) > 0
                 and n_ok + n_dictroute != n_rows - n_shell_rows),
        "fulltable_crosscheck_pass": report["fulltable_crosscheck_pass"],
        "coverage_gate_all_slices_pass": gates_ok,
        "certification_composite_pass": composite_ok and faults_ok,
        "planted_faults_all_detected": faults_ok,
        "wall_s_total": round(time.time() - t_all, 2),
        "max_rss_kb_self": maxrss_kb(),
        "max_rss_kb_children_kira": maxrss_kb(resource.RUSAGE_CHILDREN),
    }
    if args.prereg_shape:
        pred = {int(k): v for k, v in json.loads(args.prereg_shape).items()}
        meas = report.get("measured_selected_by_t", {})
        report["P_SHAPE"] = {
            str(t): {"pred": pred.get(t), "meas": meas.get(t, 0),
                     "dev_pct": (100.0 * (meas.get(t, 0) - pred[t]) / pred[t])
                     if pred.get(t) else None,
                     "within_30pct": pred.get(t) is not None and
                     abs(meas.get(t, 0) - pred[t]) <= 0.3 * pred[t]}
            for t in sorted(set(pred) | set(meas))}
    if args.prereg_front is not None:
        peak = max(st["peak_live_rows"]
                   for sl in per_slice
                   for st in sl["stats_per_stratum"].values())
        report["P_FRONT"] = {"bound": args.prereg_front,
                             "measured_peak_live_rows": peak,
                             "pass": peak <= args.prereg_front}
    with open(os.path.join(args.out, f"FAMILY_REPORT_{fam}.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    print(json.dumps(report["summary"]))
    return 1 if report["summary"]["VOID"] else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""blade_family -- TURNKEY one-command CLI for the Wolfram-free Blade pipeline.

A user (or driver script) points at THIS file (or `python3 -m blade.blade_family`)
and never needs to touch the stage modules.  Three subcommands:

  prepare <family_spec.json> --workdir W [--max-round-seconds 1800]
          [--bank-dir DIR] [--bisect] [--nthreads N]
      Runs the full chain on a FRESH workdir W:
        extension -> database -> scheme -> ansatz -> search (escalation with a
        measured STOP-line) -> grow-prime rational export -> held-out oracle
        gate -> banks a self-describing family dir (btform.json manifest).
      PARTIAL closures bank honestly: closed blocks are usable via `probe`,
      open blocks are listed in btform.json["open_works"] with their measured
      stop reasons.  Exit codes: 0 = all blocks closed; 2 = partial closure
      banked; 1/3 = typed failure (ERROR.json names the artifact; nothing
      half-banked).

  probe <family_dir> --points FILE [--prime P] [--batch] [--out FILE]
      BTOracle round-trip on a banked family dir: reduces the closed targets
      to master-coefficient VALUES mod P at the given points.  FILE is JSON:
      either  [[x1,...,eps], ...]  (coords in btform.json["params"] order,
      eps LAST) or {"prime": P, "points": [[...], ...]}.  Values go to stdout
      as JSON (and --out FILE).  Works at ANY prime (exact backend).

  status <family_dir>
      Manifest + gate summary: loads the BT form (sha256 provenance gates
      fire here), prints family/params/closed/open blocks/relation counts.
      Exit 0 iff the bank loads clean.

family_spec JSON (prepare)
--------------------------
{
  "family": "db",                              # REQUIRED integral-family head
  "search_params": ["t", "eps"],               # REQUIRED, eps LAST
  "targets": [[1,1,1,1,1,1,1,-3,0], ...],      # REQUIRED (database row order)
  "masters": [[0,1,0,0,1,0,1,0,0], ...],       # database order contract;
                                               # omitted => derived (sorted)
                                               # from the kira table RHS
  # ---- database lane (auto-picked from what you provide) -------------------
  "kira_table": "/path/kira_targets_ext.m",    # lane A: kira2math symbolic
                                               # table, evaluated mod p
  "kira_config_dir": "/path/kira",             # lane A': dir with
                                               # results/*/kira_*.m (the table
                                               # is DISCOVERED; if none exists
                                               # you get a typed refusal --
                                               # run Kira with a kira2math
                                               # export first)
  "database_dir": "/path/database",            # lane B: banked data_p*/red_*
                                               # pools (fixed points; growth
                                               # capped at the pool)
  "fflow_system": "/path/sys.json",            # lane C: raw IBP(-like) sparse
                                               # system in the fflow JSON
                                               # grammar, reduced per point by
                                               # the fflow sparse eliminator
                                               # (fflowcli learn/evalmany; no
                                               # Kira anywhere).  Masters must
                                               # occupy the HIGHEST column ids
                                               # (left-to-right RREF pivots).
  "fflow_columns": "/path/columns.json",       # lane C REQUIRED companion:
                                               # {"columns": [[...], ...]} (or
                                               # a bare list; a path or the
                                               # inline list) mapping system
                                               # column id -> integral indices
  "fflow_homog": true,                         # lane C: homogeneous system
                                               # (default true)
  # Providing none of the four lane keys above is a typed refusal.
  # ---- kira-table coefficient environment ----------------------------------
  "coeff_env": {"t": "t", "d": "4-2*eps"},     # every symbol appearing in the
                                               # table's coefficients ->
                                               # expression in search_params;
                                               # default: identity on
                                               # search_params + d=4-2*eps
  # ---- extension (optional) -------------------------------------------------
  "zero_sectors": [[1,0,...], ...],            # scaleless sectors (upstream
                                               # maximal-cut output); without
                                               # it the extension record is a
                                               # SUPERSET (no elimination)
  # ---- knobs (all optional, defaults shown) ---------------------------------
  "n_primes": 3,                               # fit primes (>= 3: 2 CRT + 1
                                               # held-out verify minimum)
  "max_primes": 8,                             # grow-prime export ceiling
  "held_out_prime": 9223372036854775507,       # independent oracle prime;
                                               # default: last BIG_UINT prime
                                               # (reserved, never fitted)
  "seed": 20260709,                            # point-sampling seed
  "scheme_mode": "single_block",               # or "full" (BLGenerateJob
                                               # multi-block scheme)
  "scheme_order": "given",                     # single_block: "given"|"weight"
  "recipe": {"ansatz_size": 7, "data_size": null,   # null => 10/30/100/300/
             "level_step": 1, "nmax": 20,           # 1000 by #params
             "break_rounds": 10,
             "var_groups": [["t","eps"]],           # default: one group =
             "var_weights": [[1,1]],                # search_params, weight 1
             "intmode": ["dimension"], "cut": []}
}

Gate discipline (inherited from the stage modules, nothing new trusted):
artifact gates only; rc recorded never gating; CRT export verifies EVERY
entry at >= 1 held-out prime and grows the prime set on overflow (fail-closed
ReconstructionOverflow at the ceiling); the banked relations additionally
pass a held-out relation oracle at a prime never used in search/fit/verify
(kira lane: the independent Kira-expression oracle; fflow-system lane: the
system's own sparse solve), with a synthetic-truth mutation control.  Every stage
wall is measured into PREPARE_LEDGER.json.  The escalation STOP-line refuses
to LAUNCH a round whose projected wall (measured ratio extrapolation) exceeds
--max-round-seconds, and kills (exact PID) any round that overruns it.

Banked family dir layout (self-describing; `probe`/`status` need nothing else):
  btform.json        manifest: params/labels/target+master gids/closed works/
                     relations sha256s/open_works + envelope notes/prepare walls
  relations.json     exact-rational relations (CRT+Wang, held-out verified)
  run_tree/search/   full search tree (database pools, schemes, fit tables,
                     kinematics) => enables the ssolve C backend of the probe

Typical lane usage:
  python3 tools/blade/blade_family.py prepare spec.json --workdir runs/<p>/<tag>/fam
  python3 tools/blade/blade_family.py status  runs/<p>/<tag>/fam/family/<name>
  python3 tools/blade/blade_family.py probe   runs/<p>/<tag>/fam/family/<name> \
          --points pts.json --prime 9223372036854775783
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import shutil
import sys
import time
from typing import Callable, Dict, List, Optional, Sequence, Tuple

if __package__ in (None, ""):                      # run as a plain script
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from blade import formats, search as bs                    # noqa: E402
    from blade.ansatz import FamilyAnsatz, integral_dimension  # noqa: E402
    from blade.extension import Extension                      # noqa: E402
    from blade.pipeline import BladeGateError                  # noqa: E402
    from blade.probe import BTForm, BTOracle                   # noqa: E402
    from blade.scheme import (RedTable, build_scheme,          # noqa: E402
                                  single_block_scheme, write_scheme)
    from blade.formats import Scheme, big_uint_primes          # noqa: E402
else:
    from . import formats, search as bs
    from .ansatz import FamilyAnsatz, integral_dimension
    from .extension import Extension
    from .pipeline import BladeGateError
    from .probe import BTForm, BTOracle
    from .scheme import (RedTable, build_scheme, single_block_scheme,
                         write_scheme)
    from .formats import Scheme, big_uint_primes


# --------------------------------------------------------------------------
# typed errors
# --------------------------------------------------------------------------

class SpecError(BladeGateError):
    """The family_spec JSON is missing/inconsistent -- fix the spec."""


class LaneUnavailable(BladeGateError):
    """The spec does not provide enough to build a database -- the error
    message names what to produce first (and with which banked recipe)."""


class NothingClosed(BladeGateError):
    """Escalation closed ZERO blocks -- nothing bankable.  The per-work stop
    reasons are in PREPARE_LEDGER.json."""


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# --------------------------------------------------------------------------
# spec loading / lane resolution
# --------------------------------------------------------------------------

_REQUIRED = ("family", "search_params", "targets")
_DATASIZE_BY_NVARS = {1: 10, 2: 30, 3: 100, 4: 300}


def load_spec(path: str) -> dict:
    if not os.path.exists(path):
        raise SpecError(f"family_spec not found: {path}")
    try:
        spec = json.load(open(path))
    except json.JSONDecodeError as e:
        raise SpecError(f"{path}: not valid JSON ({e})")
    for k in _REQUIRED:
        if k not in spec:
            raise SpecError(f"{path}: required key {k!r} missing "
                            f"(required: {_REQUIRED})")
    if not spec["search_params"] or spec["search_params"][-1] != "eps":
        raise SpecError(f"{path}: search_params must end in 'eps' "
                        f"(BLSearchParameter order contract); got "
                        f"{spec['search_params']}")
    nix = {len(t) for t in spec["targets"]}
    if len(nix) != 1:
        raise SpecError(f"{path}: targets have mixed index lengths {nix}")
    if "masters" in spec:
        nix |= {len(m) for m in spec["masters"]}
        if len(nix) != 1:
            raise SpecError(f"{path}: master/target index lengths differ")
    return spec


def resolve_lane(spec: dict) -> Tuple[str, str]:
    """-> ("kira-table", table_path) | ("banked-db", database_dir) |
    ("fflow-system", system_path)."""
    if spec.get("kira_table"):
        p = spec["kira_table"]
        if not os.path.exists(p):
            raise SpecError(f"kira_table not found: {p}")
        return "kira-table", p
    if spec.get("kira_config_dir"):
        d = spec["kira_config_dir"]
        hits = sorted(glob.glob(os.path.join(d, "results", "*", "kira_*.m"))
                      + glob.glob(os.path.join(d, "results", "*", "*.m")))
        hits = [h for h in hits if os.path.getsize(h) > 2]
        if not hits:
            raise LaneUnavailable(
                f"kira_config_dir {d}: no kira2math results table found under "
                f"results/*/;  run Kira with a kira2math export first and "
                f"re-run prepare, or pass kira_table directly")
        return "kira-table", hits[0]
    if spec.get("database_dir"):
        d = spec["database_dir"]
        pools = sorted(glob.glob(os.path.join(d, "data_p*")))
        if not pools:
            raise LaneUnavailable(f"database_dir {d}: no data_p*/ pools")
        return "banked-db", d
    if spec.get("fflow_system"):
        p = spec["fflow_system"]
        if not os.path.exists(p):
            raise SpecError(f"fflow_system not found: {p}")
        if not spec.get("fflow_columns"):
            raise SpecError(
                "fflow_system needs the fflow_columns companion (system "
                "column id -> integral index map; {'columns': [[...], ...]} "
                "or a bare list, as a path or inline)")
        return "fflow-system", p
    raise LaneUnavailable(
        "spec provides none of kira_table / kira_config_dir / database_dir "
        "/ fflow_system -- provide one (kira2math table, Kira config dir "
        "with results, banked red_* pools, or a raw fflow JSON sparse "
        "system plus fflow_columns)")


def load_fflow_columns(doc) -> List[Tuple[int, ...]]:
    """fflow_columns: a path to JSON ({'columns': [...]} or a bare list) or
    the inline list itself -> list of integral index tuples, one per system
    column id."""
    if isinstance(doc, str):
        if not os.path.exists(doc):
            raise SpecError(f"fflow_columns file not found: {doc}")
        doc = json.load(open(doc))
    if isinstance(doc, dict):
        doc = doc.get("columns")
    if not isinstance(doc, list) or not doc:
        raise SpecError("fflow_columns: expected {'columns': [[...], ...]} "
                        "or a nonempty list of integral index tuples")
    return [tuple(int(x) for x in c) for c in doc]


def make_env_fn(search_params: Sequence[str],
                coeff_env: Dict[str, str]) -> Callable:
    """Compile the coefficient-symbol environment: every symbol used by the
    kira table's coefficients -> expression in search_params, evaluated in
    GF(p) (integer literals are exact; '/' is modular division)."""
    compiled = {}
    for name, expr in coeff_env.items():
        compiled[name] = compile(str(expr).replace("^", "**"),
                                 f"<coeff_env:{name}>", "eval")

    def env_fn(coords: Sequence[int], p: int) -> dict:
        base = {v: bs._GF(c, p) for v, c in zip(search_params, coords)}
        ns = dict(base)
        ns["__builtins__"] = {}
        out = {}
        for name, code in compiled.items():
            v = eval(code, ns)
            out[name] = v if isinstance(v, bs._GF) else bs._GF(v, p)
        return out

    return env_fn


def default_coeff_env(search_params: Sequence[str]) -> Dict[str, str]:
    env = {p: p for p in search_params}
    if "eps" in search_params:
        env["d"] = "4-2*eps"
    return env


def derive_masters(table_path: str, family: str,
                   targets: Sequence[Tuple[int, ...]]) -> List[Tuple[int, ...]]:
    """Masters = every family[...] tuple on the RHS of a target rule, sorted
    ascending (deterministic database-order convention when the spec does not
    pin one)."""
    import re
    jb = re.compile(re.escape(family) + r"\[([0-9,\- ]+)\]")
    txt = open(table_path).read()
    tset = {tuple(t) for t in targets}
    masters = set()
    for chunk in re.split(r",\s*(?=" + re.escape(family) + r"\[)",
                          txt.strip()[1:-1]):
        if "->" not in chunk:
            continue
        lhs_txt, rhs_txt = chunk.split("->", 1)
        m = jb.match(lhs_txt.strip())
        if not m:
            continue
        lhs = tuple(int(x) for x in m.group(1).split(","))
        if lhs not in tset:
            continue
        for mm in jb.findall(rhs_txt):
            masters.add(tuple(int(x) for x in mm.split(",")))
    if not masters:
        raise SpecError(f"{table_path}: no target rule matched -- wrong "
                        f"family head {family!r} or targets not in table")
    return sorted(masters)


class _DBView:
    """Duck-typed subset of SearchFamily that make_kira_grower needs, so the
    database can be grown BEFORE the scheme/ansatz exist."""

    def __init__(self, searchdir: str, datanames: List[str]):
        self.searchdir = searchdir
        self.datanames = datanames
        self.log = print

    def dbdir(self, pid: int) -> str:
        return os.path.join(self.searchdir, "database", self.datanames[pid])


# --------------------------------------------------------------------------
# prepare
# --------------------------------------------------------------------------

def cmd_prepare(args) -> int:
    t00 = time.time()
    workdir = os.path.abspath(args.workdir)
    # FRESH-workdir contract, checked BEFORE any write (a reused tree carries
    # stale flag/config/pool state that silently corrupts a rerun -- measured:
    # stale flag==1 skips escalation over a reset 2-point pool; and the
    # refusal must never clobber the previous run's ledger).
    stale = [s for s in ("search", "family", "relations.json",
                         "PREPARE_LEDGER.json") if
             os.path.exists(os.path.join(workdir, s))]
    if stale:
        os.makedirs(workdir, exist_ok=True)
        msg = (f"workdir {workdir} is not fresh ({', '.join(stale)} exist) "
               f"-- prepare refuses stale state; use a new --workdir tag "
               f"(expert resume lives in the stage modules, not the turnkey)")
        with open(os.path.join(workdir, "ERROR.json"), "w") as f:
            json.dump({"type": "BladeGateError", "message": msg,
                       "workdir": workdir}, f, indent=1)
        print(f"PREPARE FAILED [BladeGateError]: {msg}", file=sys.stderr)
        return 1
    os.makedirs(workdir, exist_ok=True)
    ledger = {"cli": "blade_family.py", "cmd": "prepare",
              "spec": os.path.abspath(args.spec),
              "spec_sha256": None, "lane": None, "gates": [], "walls": {},
              "works": {}, "stop_line_seconds": args.max_round_seconds}
    lpath = os.path.join(workdir, "PREPARE_LEDGER.json")

    def bank_ledger():
        with open(lpath, "w") as f:
            json.dump(ledger, f, indent=1, default=str)

    def gate(name: str, ok: bool, detail: str = ""):
        ledger["gates"].append({"gate": name, "pass": bool(ok),
                                "detail": detail})
        print(f"GATE {'PASS' if ok else 'FAIL'}: {name} {detail}")
        if not ok:
            bank_ledger()
            raise BladeGateError(f"gate failed: {name} -- {detail} "
                                 f"(ledger: {lpath})")

    try:
        return _prepare_inner(args, workdir, ledger, lpath, bank_ledger,
                              gate, t00)
    except BladeGateError as e:
        err = {"type": type(e).__name__, "message": str(e),
               "ledger": lpath, "workdir": workdir}
        with open(os.path.join(workdir, "ERROR.json"), "w") as f:
            json.dump(err, f, indent=1)
        bank_ledger()
        print(f"PREPARE FAILED [{type(e).__name__}]: {e}", file=sys.stderr)
        print(f"artifacts: {os.path.join(workdir, 'ERROR.json')}",
              file=sys.stderr)
        return 3 if isinstance(e, NothingClosed) else 1


def _prepare_inner(args, workdir, ledger, lpath, bank_ledger, gate,
                   t00) -> int:
    spec = load_spec(args.spec)
    ledger["spec_sha256"] = _sha256(args.spec)
    lane, lane_path = resolve_lane(spec)
    ledger["lane"] = {"kind": lane, "path": lane_path}
    print(f"[prepare] lane: {lane} <- {lane_path}")

    family = spec["family"]
    params = tuple(spec["search_params"])
    targets = [tuple(t) for t in spec["targets"]]
    seed = int(spec.get("seed", 20260709))
    n_primes = int(spec.get("n_primes", 3))
    if n_primes < 3:
        raise SpecError("n_primes must be >= 3 (2 CRT + 1 held-out verify)")

    # ---- stage: extension ---------------------------------------------------
    t0 = time.time()
    ext = Extension(zero_sectors=[tuple(z) for z in
                                  spec.get("zero_sectors", [])])
    extended = ext.generate_scheme(targets)
    if ext.zero_sectors:
        # BLEliminateZeroSectors semantics (the maximal-cut recipe):
        # without this the extension record keeps every cut-invalid
        # staircase integral (measured on a 13-cut maximal-cut family:
        # 528,384 -> 70).
        # The single_block scheme path never consumed the superset, so this
        # only fixes the RECORD (and its wall/JSON size), not the relations.
        extended = ext.eliminate_zero_sectors(extended)
    ext_doc = {"family": family, "targets": [list(t) for t in targets],
               "n_extended": len(extended),
               "zero_sectors_supplied": "zero_sectors" in spec,
               "note": ("zero_sectors not supplied: extended set is a "
                        "SUPERSET (no scaleless elimination)"
                        if "zero_sectors" not in spec else ""),
               "extended": [list(v) for v in extended]}
    ext_path = os.path.join(workdir, "extension", "extended_integrals.json")
    os.makedirs(os.path.dirname(ext_path), exist_ok=True)
    with open(ext_path, "w") as f:
        json.dump(ext_doc, f, indent=1)
    ledger["walls"]["extension"] = time.time() - t0
    gate("extension covers targets", set(targets) <= set(extended),
         f"{len(extended)} extended integrals -> {ext_path}")

    # ---- stage: database ----------------------------------------------------
    t0 = time.time()
    searchdir = os.path.join(workdir, "search")
    dbroot = os.path.join(searchdir, "database")
    os.makedirs(dbroot, exist_ok=True)
    P = big_uint_primes()
    kc = None
    fsys = None
    if lane == "kira-table":
        masters = ([tuple(m) for m in spec["masters"]] if "masters" in spec
                   else derive_masters(lane_path, family, targets))
        if "masters" not in spec:
            ledger["masters_derived"] = ("sorted ascending from kira-table "
                                         "RHS (spec pinned no order)")
        coeff_env = spec.get("coeff_env", default_coeff_env(params))
        ledger["coeff_env"] = coeff_env
        kc = bs.KiraCoeffs(lane_path, family, targets, masters,
                           env_fn=make_env_fn(params, coeff_env))
        NT, NM = len(targets), len(masters)
        covered = {i for (i, _j) in kc.cexpr}
        gate("kira table covers every target", covered == set(range(NT)),
             f"missing target rows: {sorted(set(range(NT)) - covered)}")
        groups = [[i for i in range(NT) if (i, m) in kc.cexpr] + [NT + m]
                  for m in range(NM)]
        posi = [i for g in groups for i in g]
        part = [0]
        for g in groups:
            part.append(part[-1] + len(g))
        datanames = [f"data_p{pid}" for pid in range(n_primes)]
        for pid in range(n_primes):
            db = formats.Database(prime=P[pid], npara=len(params),
                                  nint=NT + NM, nmaster=NM, nentry=len(posi),
                                  nps=0, part=part, posi=posi, ps=[], data=[],
                                  pid_raw=f"{pid}\n", cons_raw="{}\n")
            db.write(os.path.join(dbroot, datanames[pid]))
        grower0 = bs.make_kira_grower(_DBView(searchdir, datanames), kc,
                                      seed=seed)
        grower0(2)                       # scheme needs 2 sample points
    elif lane == "fflow-system":
        fsys = bs.FflowSystem(
            lane_path, load_fflow_columns(spec["fflow_columns"]), targets,
            masters=([tuple(m) for m in spec["masters"]]
                     if "masters" in spec else None),
            homog=bool(spec.get("fflow_homog", True)),
            nthreads=args.nthreads,
            timeout=float(args.max_round_seconds),
            scratch=os.path.join(workdir, "fflow"))
        masters = fsys.masters
        if "masters" not in spec:
            ledger["masters_derived"] = ("sorted ascending from the learned "
                                         "independent set (spec pinned no "
                                         "order)")
        NT, NM = len(targets), len(masters)
        ledger["fflow_learn"] = {"nparsout": fsys.nparsout,
                                 "npars": fsys.npara, "neqs": fsys.neqs}
        # target coverage + independents==masters already gated typed inside
        # FflowSystem._learn; this gate pins the parameter-count contract
        gate("fflow system npars == len(search_params)",
             fsys.npara == len(params),
             f"{NT} targets, {NM} masters, nparsout={fsys.nparsout}, "
             f"system npars={fsys.npara} vs search_params {len(params)}")
        stored = fsys.measure_pattern(P[0], npoints=2, seed=seed)
        gate("fflow pattern probe: every target row stores >= 1 master",
             all(any((i, m) in stored for m in range(NM))
                 for i in range(NT)),
             f"{len(stored)} stored (target,master) pairs at 2 probe points")
        groups = [[i for i in range(NT) if (i, m) in stored] + [NT + m]
                  for m in range(NM)]
        posi = [i for g in groups for i in g]
        part = [0]
        for g in groups:
            part.append(part[-1] + len(g))
        datanames = [f"data_p{pid}" for pid in range(n_primes)]
        for pid in range(n_primes):
            db = formats.Database(prime=P[pid], npara=len(params),
                                  nint=NT + NM, nmaster=NM, nentry=len(posi),
                                  nps=0, part=part, posi=posi, ps=[], data=[],
                                  pid_raw=f"{pid}\n", cons_raw="{}\n")
            db.write(os.path.join(dbroot, datanames[pid]))
        grower0 = bs.make_fflow_grower(_DBView(searchdir, datanames), fsys,
                                       seed=seed)
        grower0(2)                       # scheme needs 2 sample points
    else:                                # banked-db
        pools = sorted(glob.glob(os.path.join(lane_path, "data_p*")),
                       key=lambda p: int(p.rsplit("data_p", 1)[1]))
        datanames = [os.path.basename(p) for p in pools]
        for p, dn in zip(pools, datanames):
            shutil.copytree(p, os.path.join(dbroot, dn))
        if "masters" not in spec:
            raise SpecError("banked-db lane needs explicit 'masters' "
                            "(database row-order contract)")
        masters = [tuple(m) for m in spec["masters"]]
        NT, NM = len(targets), len(masters)
        hdr = formats.read_ints(os.path.join(dbroot, datanames[0],
                                             "red_common"))
        gate("banked database header matches spec",
             hdr[1] == len(params) and hdr[2] == NT + NM and hdr[3] == NM,
             f"red_common={hdr} vs npara={len(params)} nint={NT+NM} nm={NM}")
        if spec.get("kira_table"):
            pass                          # kc already impossible here (lane)
    all_ints = targets + masters
    labels = [family + "[" + ",".join(map(str, t)) + "]" for t in all_ints]
    ledger["walls"]["database"] = time.time() - t0
    ledger["database"] = {"lane": lane, "datanames": datanames,
                          "primes": [formats.read_ints(os.path.join(
                                         dbroot, dn, "red_common"))[0]
                                     for dn in datanames]}

    # ---- stage: scheme --------------------------------------------------------
    t0 = time.time()
    rt, rt2 = RedTable.pair(os.path.join(dbroot, datanames[0]))
    mode = spec.get("scheme_mode", "single_block")
    if mode == "single_block":
        res = single_block_scheme(all_ints, NM, rt, rt2,
                                  order=spec.get("scheme_order", "given"))
    elif mode == "full":
        res = build_scheme(all_ints, NM, rt, rt2)
    else:
        raise SpecError(f"scheme_mode {mode!r} not in "
                        f"{{'single_block','full'}}")
    jobdir = os.path.join(searchdir, "job0")
    write_scheme(res, jobdir)
    nworks = len(res.blocks)
    ledger["walls"]["scheme"] = time.time() - t0
    ledger["scheme"] = {"mode": mode, "n_blocks": nworks,
                        "sizes": res.sizes(), "notes": res.notes}
    gate("scheme produced >= 1 block", nworks >= 1, res.summary())

    # ---- stage: ansatz recipe + search ----------------------------------------
    rc = dict(spec.get("recipe", {}))
    data_size = rc.get("data_size") or _DATASIZE_BY_NVARS.get(len(params),
                                                              1000)
    fam = FamilyAnsatz(
        search_params=params,
        var_groups=rc.get("var_groups", [list(params)]),
        var_weights=rc.get("var_weights", [[1] * len(params)]),
        intmode=rc.get("intmode", ["dimension"]),
        cut=rc.get("cut", []),
        nmax=int(rc.get("nmax", 20)))
    dims = [integral_dimension(x) for x in all_ints]
    work_weights = [[dims[gid] for gid in b.ids] for b in res.blocks]
    cfg = bs.SearchFamily(
        searchdir=searchdir, job="job0", fam=fam, work_weights=work_weights,
        datanames=list(datanames), nthreads=args.nthreads,
        ansatz_size=int(rc.get("ansatz_size", 7)),
        level_base=int(rc.get("level_base", 0)),
        level_step=int(rc.get("level_step", 1)),
        data_size=int(data_size),
        break_rounds=int(rc.get("break_rounds", 10)),
        datasize_policy=rc.get("datasize_policy", "reactive"),
        round_timeout=float(args.max_round_seconds),
        stop_line_seconds=float(args.max_round_seconds))
    if kc is not None:
        cfg.grower = bs.make_kira_grower(cfg, kc, seed=seed)
    elif fsys is not None:
        cfg.grower = bs.make_fflow_grower(cfg, fsys, seed=seed)

    t0 = time.time()
    outcomes = [bs.escalate_work(cfg, w) for w in range(nworks)]
    ledger["walls"]["escalation"] = time.time() - t0
    for o in outcomes:
        ledger["works"][o.work] = {
            "closed": o.closed, "stop_reason": o.stop_reason,
            "rounds": [vars(r) for r in o.rounds]}
    closed = [o.work for o in outcomes if o.closed]
    open_works = [o for o in outcomes if not o.closed]
    print(f"[prepare] escalation: {len(closed)}/{nworks} blocks closed "
          f"({ledger['walls']['escalation']:.1f}s)")
    if not closed:
        raise NothingClosed(
            f"0/{nworks} blocks closed; stop reasons: "
            + "; ".join(f"w{o.work}: {o.stop_reason}" for o in open_works))

    # ---- stage: fit at every prime + grow-prime rational export ---------------
    # Multi-block footgun (measured here, failmode gate F): each work's
    # escalation set_data_number()s pool 0 to ITS OWN want, so after driving
    # several works the logical nps can be far below what an earlier, bigger
    # block's fit needs (fitrel flag==0 self_check).  Fit at the FULL physical
    # pool (upstream fits at searchnumber; bisection shrinks afterwards).
    bs.set_data_number(cfg.dbdir(0), bs._file_rows(cfg.dbdir(0)))
    t0 = time.time()
    for pid in range(len(cfg.datanames)):
        for w in closed:
            shutil.rmtree(os.path.join(cfg.workdir(w), "fit",
                                       cfg.datanames[pid]),
                          ignore_errors=True)
            bs.run_fitrel_works(cfg, cfg.datanames[pid], w, w + 1)
    ledger["walls"]["fitrel"] = time.time() - t0

    held_prime = int(spec.get("held_out_prime", P[-1]))
    fit_primes = set(ledger["database"]["primes"])
    gate("held-out prime not a fit prime", held_prime not in fit_primes,
         str(held_prime))
    prime_adder = None
    if kc is not None:
        prime_adder = lambda: bs.add_fit_prime(          # noqa: E731
            cfg, kc, closed, seed, reserved=(held_prime,),
            timeout=float(args.max_round_seconds))
    elif fsys is not None:
        prime_adder = lambda: bs.add_fit_prime(          # noqa: E731
            cfg, None, closed, seed, reserved=(held_prime,),
            timeout=float(args.max_round_seconds),
            grower_factory=lambda pids: bs.make_fflow_grower(
                cfg, fsys, seed=seed, pids=pids))
    t0 = time.time()
    rel_path = os.path.join(workdir, "relations.json")
    doc = bs.export_relations(
        cfg, closed, labels, rel_path,
        extra_meta={"family": family,
                    "provenance": f"blade_family prepare, lane {lane}, "
                                  f"spec sha256 {ledger['spec_sha256'][:16]}"},
        prime_adder=prime_adder,
        max_primes=int(spec.get("max_primes", 8)))
    ledger["walls"]["export"] = time.time() - t0
    g1_total = sum(Scheme.read(cfg.workdir(w)).g1 for w in closed)
    gate("export: n_relations == sum g1 over closed works",
         doc["n_relations"] == g1_total,
         f"{doc['n_relations']} vs {g1_total}")

    # ---- stage: held-out oracle (kira / fflow-system lanes) -------------------
    if kc is not None:
        t0 = time.time()
        checks, fails = bs.heldout_relation_check(
            doc, kc, labels, NM, held_prime=held_prime, npoints=8,
            seed=seed + 1)
        ledger["walls"]["heldout"] = time.time() - t0
        gate("held-out relation oracle (fresh prime + mutation control)",
             fails == 0, f"{checks} checks, {fails} fails")
    elif fsys is not None:
        # honesty boundary: the rows replay the input system's OWN sparse
        # solve at a fresh prime -- this gates the search/fit/export chain,
        # not the input system (validate the system upstream)
        t0 = time.time()
        checks, fails = bs.heldout_relation_check(
            doc, fsys, labels, NM, held_prime=held_prime, npoints=8,
            seed=seed + 1)
        ledger["walls"]["heldout"] = time.time() - t0
        gate("held-out relation gate vs the system's own solve "
             "(fresh prime + mutation control)",
             fails == 0, f"{checks} checks, {fails} fails")
    else:
        ledger["heldout"] = ("SKIPPED: banked-db lane has no independent "
                             "coefficient oracle (provide kira_table too "
                             "to enable it)")
        print(f"[prepare] WARNING: {ledger['heldout']}")

    if args.bisect:
        t0 = time.time()
        fn = bs.fitnumber_bisection(cfg, closed, pid=0)
        ledger["walls"]["bisection"] = time.time() - t0
        doc["fitnumber"] = fn
        with open(rel_path, "w") as f:
            json.dump(doc, f, indent=1)

    # ---- stage: bank the self-describing family dir ---------------------------
    t0 = time.time()
    bank_dir = os.path.abspath(args.bank_dir) if args.bank_dir else \
        os.path.join(workdir, "family", family)
    if os.path.exists(bank_dir):
        raise BladeGateError(f"bank dir {bank_dir} already exists -- "
                             f"refusing to overwrite a banked family")
    os.makedirs(bank_dir)
    shutil.copy(rel_path, os.path.join(bank_dir, "relations.json"))
    shutil.copytree(searchdir, os.path.join(bank_dir, "run_tree", "search"))
    works_meta = []
    for w in closed:
        schw = Scheme.read(cfg.workdir(w))
        works_meta.append({"work": w, "g1_gids": schw.intid[:schw.g1]})
    closed_gids = [g for wm in works_meta for g in wm["g1_gids"]]
    manifest = {
        "family": family,
        "params": list(params),
        "labels": labels,
        "target_gids": list(range(NT)),
        "master_gids": list(range(NT, NT + NM)),
        "works": works_meta,
        "relations_files": ["relations.json"],
        "provenance": {"relations.json":
                       _sha256(os.path.join(bank_dir, "relations.json"))},
        "ssolve": {"search_tree": "run_tree/search", "job": "job0",
                   "datanames": list(cfg.datanames), "nworks": nworks},
        "status": "closed" if not open_works else "partial",
        "open_works": [{"work": o.work, "stop_reason": o.stop_reason}
                       for o in open_works],
        "targets": [list(t) for t in targets],
        "masters": [list(m) for m in masters],
        "prepare": {
            "spec_sha256": ledger["spec_sha256"],
            "lane": lane,
            "seed": seed,
            "held_out_prime": held_prime,
            "fit_primes": ledger["database"]["primes"],
            "stop_line_seconds": args.max_round_seconds,
            "walls": dict(ledger["walls"], TOTAL=time.time() - t00),
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
        "envelope": {
            "closed_target_gids": [g for g in range(NT) if g in closed_gids],
            "open_target_gids": [g for g in range(NT)
                                 if g not in closed_gids],
            "note": ("all blocks closed" if not open_works else
                     "PARTIAL closure: open blocks listed in open_works; "
                     "probe serves the closed blocks only"),
        },
    }
    with open(os.path.join(bank_dir, "btform.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    ledger["walls"]["bank"] = time.time() - t0
    ledger["bank_dir"] = bank_dir

    # load-gate the bank through the probe itself
    oracle = BTOracle(bank_dir)
    gate("banked family loads through BTOracle",
         set(closed_gids) <= set(oracle.form.closed_gids),
         f"closed gids {oracle.form.closed_gids}, load "
         f"{oracle.load_wall:.3f}s")

    ledger["walls"]["TOTAL"] = time.time() - t00
    ledger["status"] = manifest["status"]
    bank_ledger()
    print(f"[prepare] BANKED {manifest['status'].upper()}: {bank_dir}")
    print(f"[prepare] total wall {ledger['walls']['TOTAL']:.1f}s "
          f"(ledger: {lpath})")
    return 0 if not open_works else 2


# --------------------------------------------------------------------------
# probe
# --------------------------------------------------------------------------

def cmd_probe(args) -> int:
    try:
        oracle = BTOracle(args.family_dir)
        form = oracle.form
        raw = json.load(open(args.points))
        prime = args.prime
        pts_in = raw
        if isinstance(raw, dict):
            pts_in = raw["points"]
            prime = prime or raw.get("prime")
        prime = int(prime) if prime else big_uint_primes()[0]
        npar = len(form.params)
        pts = []
        for i, row in enumerate(pts_in):
            if len(row) != npar:
                raise BladeGateError(
                    f"{args.points}: point {i} has {len(row)} coords; family "
                    f"params are {form.params} (need {npar}, eps LAST)")
            pts.append(([int(c) for c in row[:-1]], int(row[-1])))
        t0 = time.time()
        results = oracle.reduce_many(pts, prime)
        wall = time.time() - t0
        out = {"family": form.family, "prime": prime,
               "params": form.params,
               "basis": [form.label(g) for g in form.independent_gids],
               "closed_targets": [form.label(g) for g in form.target_gids
                                  if g in form.closed_gids],
               "points": [list(k) + [e] for k, e in pts],
               "results": results,
               "wall_seconds": wall}
        txt = json.dumps(out, indent=1)
        print(txt)
        if args.out:
            with open(args.out, "w") as f:
                f.write(txt + "\n")
        return 0
    except (BladeGateError, OSError, KeyError, ValueError,
            json.JSONDecodeError) as e:
        json.dump({"type": type(e).__name__, "message": str(e),
                   "family_dir": args.family_dir}, sys.stderr, indent=1)
        sys.stderr.write("\n")
        return 1


# --------------------------------------------------------------------------
# status
# --------------------------------------------------------------------------

def cmd_status(args) -> int:
    fdir = os.path.abspath(args.family_dir)
    try:
        meta = json.load(open(os.path.join(fdir, "btform.json")))
        t0 = time.time()
        form = BTForm(fdir)               # sha256 provenance gates fire here
        wall = time.time() - t0
        nrel = sum(len(w.relations) for w in form.works.values())
        closed_t = [g for g in form.target_gids if g in form.closed_gids]
        open_t = [g for g in form.target_gids if g not in form.closed_gids]
        print(f"family     : {form.family}")
        print(f"dir        : {fdir}")
        print(f"params     : {form.params}")
        status = meta.get("status", "(pre-CLI bank: no status key)")
        print(f"status     : {status}")
        print(f"targets    : {len(form.target_gids)} total, "
              f"{len(closed_t)} closed, {len(open_t)} open")
        for g in closed_t:
            print(f"  CLOSED   {form.label(g)}")
        for g in open_t:
            print(f"  PENDING  {form.label(g)}")
        print(f"masters    : {len(form.master_gids)}"
              + (f" (+{len(form.extra_independents)} extra independents)"
                 if form.extra_independents else ""))
        print(f"works      : {sorted(form.works)} closed; relations {nrel}")
        for ow in meta.get("open_works", []):
            print(f"  OPEN w{ow['work']}: {ow['stop_reason']}")
        print(f"provenance : {len(meta['provenance'])} file(s) sha256-"
              f"verified in {wall:.3f}s")
        ss = ("available (" + str(meta["ssolve"]["datanames"]) + ")"
              if "ssolve" in meta else "not banked (exact backend only)")
        print(f"ssolve     : {ss}")
        if "prepare" in meta:
            pw = meta["prepare"].get("walls", {})
            print(f"prepare    : lane {meta['prepare'].get('lane')}, "
                  f"total {pw.get('TOTAL', float('nan')):.1f}s, "
                  f"created {meta['prepare'].get('created')}")
        print("GATES      : btform.json schema + sha256 provenance + block "
              "squareness + triangularizability all PASS (BTForm load)")
        return 0
    except (BladeGateError, OSError, KeyError, json.JSONDecodeError) as e:
        print(f"STATUS FAILED [{type(e).__name__}]: {e}", file=sys.stderr)
        return 1


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="blade_family",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepare", help="full chain: spec -> banked family dir")
    p.add_argument("spec", help="family_spec JSON (see module docstring)")
    p.add_argument("--workdir", required=True,
                   help="FRESH working dir (tree + ledger + bank land here)")
    p.add_argument("--max-round-seconds", type=float, default=1800.0,
                   help="STOP line: refuse a search round whose projected "
                        "wall exceeds this; also the hard per-round kill "
                        "timeout (default 1800)")
    p.add_argument("--bank-dir", default=None,
                   help="bank target (default <workdir>/family/<family>)")
    p.add_argument("--bisect", action="store_true",
                   help="run fitnumber bisection after export")
    p.add_argument("--nthreads", type=int, default=6)
    p.set_defaults(fn=cmd_prepare)

    p = sub.add_parser("probe", help="BTOracle round-trip on a banked family")
    p.add_argument("family_dir")
    p.add_argument("--points", required=True,
                   help="JSON: [[x1,..,eps],...] or {'prime':P,'points':[..]}")
    p.add_argument("--prime", type=int, default=None,
                   help="reduction prime (default: points-file prime, else "
                        "BIG_UINT_PRIMES[0])")
    p.add_argument("--batch", action="store_true",
                   help="accepted for compatibility; probe always batches")
    p.add_argument("--out", default=None, help="also write the JSON here")
    p.set_defaults(fn=cmd_probe)

    p = sub.add_parser("status", help="manifest + gate summary of a bank")
    p.add_argument("family_dir")
    p.set_defaults(fn=cmd_status)

    args = ap.parse_args(argv)
    try:
        return args.fn(args)
    except SpecError as e:
        print(f"SPEC ERROR: {e}", file=sys.stderr)
        return 1
    except LaneUnavailable as e:
        print(f"LANE UNAVAILABLE: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

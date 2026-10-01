"""blade.adaptive -- port of Blade's ADAPTIVE-SEARCH decision logic.

WHAT IS PORTED (the load-bearing decision logic, with .wl line cites):

  * powerlimit2weight (BLSearch/SearchRelation.wl:275): turn per-variable
    maximal single-variable-analytic relation powers into integer ansatz
    weights (weight_i ~ 1/limit_i, normalized so the HARDEST variable gets
    weight 1).  This is the "adaptive" VariableWeight semantics documented in
    SearchOptions.wl:39 ("Complexity of variable is determined by
    maximal-power of single-variable-analytic block-triangular relations").

  * getLowerLimit (SearchRelation.wl:143-150): the measured "complexity" of a
    probe job = (number of configs carrying tmp_nsol) - 1, maxed over works.
    The upstream off-by-one comparison (`If[part>max, max=part-1]`, i.e. a
    config COUNT compared against a config INDEX) is replicated verbatim.

  * single-variable probe searches (BLSearchRelations "adaptive" branch,
    SearchRelation.wl:288-296): for each variable v run
    BLSearchRelationsBasic[name, name, {{v}}, {{1}}, {intmode}, {}, ...,
    "SearchOnly"->True] -- a cheap search whose ONLY output is the closing
    level.  Probe artifacts are quarantined (their relations are relative to
    FROZEN other parameters and are never fitted or exported).

  * conservative-parameter freezing (BLIBPSystem/Database.wl:131-141): probe
    databases sample ONLY the probed variable; all other search parameters
    ($Conservative, SearchRelation.wl:174) are frozen to one shared value
    across the pool (`pslist[[All,cons]] = pslist[[1,cons]]`, Database.wl:140).

  * single-group main search (production default at nana = #params:
    SearchOptions.wl:227-235 sets VariableGroup -> ONE group of all params
    with "adaptive" weights and IntegralWeight "uniform"; the reference dbox
    recipes use "dimension", which we keep for ladder comparability).

WHAT IS **SKIPPED**, AND WHY (documented here):

  * BLAdaptiveSearch's nana-escalation loop + timeEstimate/BLTimeIBP
    (BLSearch/AdaptiveSearch.wl:39-63 and 69-74): that logic chooses HOW MANY
    parameters stay analytic (nana) by weighing BT-search cost against
    per-point sampling cost for the SemiBL production pipeline.  For w0
    closure the answer is forced: full symbolic parity requires ALL params
    analytic (nana = len(BLSearchParameter)), so the wall-clock model is not
    load-bearing.  Our timed-pilot ladder projections replace it as the cost gate.
  * multi-group sorting + cut determination (SearchRelation.wl:311-326):
    only active when VariableGroup has >1 group; the nana-max default is one
    group.  Held as a fallback lever, not ported.
  * blSolveLearning / average_time (SearchRelation.wl:249-250): SemiBL
    sampling-time estimate, not needed for closure.
  * probe round-limit: upstream Break=10 would $Failed on a slow direction;
    probes here allow more rounds (each is sub-second at probe sizes) --
    disclosed deviation, gated by the same artifact checks.

Everything here drives the existing gated machinery in search.py; no C code
is touched and no gate is bypassed.
"""
from __future__ import annotations

import math
import os
import shutil
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import formats
from .formats import Scheme, read_ints, write_fab_rows
from .ansatz import FamilyAnsatz
from .pipeline import BladeGateError
from .search import (SearchFamily, WorkOutcome, KiraCoeffs, escalate_work,
                     n_in_configs, _file_rows)

__all__ = ["powerlimit2weight", "get_lower_limit", "make_frozen_grower",
           "init_probe_database", "make_probe_job", "single_var_probe",
           "adaptive_weights"]


# --------------------------------------------------------------------------
# powerlimit2weight (SearchRelation.wl:275)
# --------------------------------------------------------------------------

def powerlimit2weight(limits: Sequence[int]) -> List[int]:
    """temp = limits*1000 /. 0->1; temp = LCM@@temp / temp;
    Floor[temp/Min[temp]].  (*1000 so a 0 limit maps to weight >> others*)"""
    temp = [l * 1000 if l * 1000 != 0 else 1 for l in limits]
    l = temp[0]
    for t in temp[1:]:
        l = l * t // math.gcd(l, t)
    temp = [l // t for t in temp]
    m = min(temp)
    return [t // m for t in temp]


# --------------------------------------------------------------------------
# getLowerLimit (SearchRelation.wl:143-150)
# --------------------------------------------------------------------------

def get_lower_limit(jobdir: str) -> int:
    """max over works of (#configs with tmp_nsol) - 1, with the verbatim
    upstream comparison `If[part > max, max = part-1]`."""
    mx = 0
    w = 0
    while os.path.isdir(os.path.join(jobdir, str(w))):
        part = 0
        while os.path.exists(os.path.join(jobdir, str(w), "config",
                                          str(part), "tmp_nsol")):
            part += 1
        if part > mx:
            mx = part - 1
        w += 1
    return mx


# --------------------------------------------------------------------------
# frozen-conservative database (Database.wl:131-141 + SearchRelation.wl:174)
# --------------------------------------------------------------------------

def init_probe_database(searchdir: str, dataname: str, src_dataname: str
                        ) -> str:
    """Create database/<dataname> with header/part/posi copied from
    src_dataname (same prime -- upstream probes run at PrimeNo 0 too,
    SearchRelation.wl:169) and an EMPTY pool.  Refuses to overwrite."""
    dbdir = os.path.join(searchdir, "database", dataname)
    src = os.path.join(searchdir, "database", src_dataname)
    if os.path.exists(dbdir):
        raise BladeGateError(f"init_probe_database: {dbdir} exists -- refusing")
    common = read_ints(os.path.join(src, "red_common"))
    os.makedirs(dbdir)
    shutil.copy(os.path.join(src, "red_part"), os.path.join(dbdir, "red_part"))
    shutil.copy(os.path.join(src, "red_posi"), os.path.join(dbdir, "red_posi"))
    write_fab_rows(os.path.join(dbdir, "red_common"),
                   [list(common[:-1]) + [0]])
    write_fab_rows(os.path.join(dbdir, "red_ps"), [])
    write_fab_rows(os.path.join(dbdir, "red_data"), [])
    return dbdir


def make_frozen_grower(cfg: SearchFamily, kc: KiraCoeffs, seed: int,
                       free_idx: int, pids: Sequence[int] = (0,)
                       ) -> Callable[[int], None]:
    """Grower with $Conservative semantics: every parameter EXCEPT free_idx is
    frozen to one shared value per pool (Database.wl:140
    `pslist[[All,cons]] = pslist[[1,cons]]`); only the probed variable
    varies.  Frozen values are drawn once per (seed, pid); the free
    coordinate is resampled on value-vanish (stored entries must be nonzero,
    FORMATS contract) while frozen values stay put, exactly as upstream
    resampling would keep the shared conservative row."""
    import random as _random

    def grow(n_total: int) -> None:
        for pid in pids:
            dbdir = cfg.dbdir(pid)
            common = read_ints(os.path.join(dbdir, "red_common"))
            prime, npara, nint, nmaster, nentry, _nps_log = common
            part = read_ints(os.path.join(dbdir, "red_part"))
            posi = read_ints(os.path.join(dbdir, "red_posi"))
            ps = formats._rows(read_ints(os.path.join(dbdir, "red_ps")),
                               npara, "red_ps")
            data = formats._rows(read_ints(os.path.join(dbdir, "red_data")),
                                 nentry, "red_data")
            if len(data) != len(ps):
                raise BladeGateError(f"{dbdir}: red_ps {len(ps)} rows != "
                                     f"red_data {len(data)} rows")
            frz_rng = _random.Random(f"{seed}:frozen:{pid}")
            frozen = [frz_rng.randrange(1, prime) for _ in range(npara)]
            if ps:                       # pool already started: freeze row 0
                frozen = list(ps[0])     # (Database.wl:136 oldps[[1,cons]])
            phys = len(ps)
            if phys >= n_total:
                continue
            nt = len(kc.targets)
            keys = sorted(kc.cexpr)
            for idx in range(phys, n_total):
                for attempt in range(200):
                    rng = _random.Random(f"{seed}:{pid}:{idx}:{attempt}")
                    coord = list(frozen)
                    coord[free_idx] = rng.randrange(1, prime)
                    vals = {}
                    ok = True
                    for (i, j) in keys:
                        v = kc.value(i, j, coord, prime)
                        if v == 0:
                            ok = False
                            break
                        vals[(i, j)] = v
                    if ok:
                        break
                else:
                    raise BladeGateError(
                        "frozen grower: 200 resamples of the free coordinate "
                        "all value-vanish -- frozen point is degenerate, "
                        "rerun with a new seed")
                row = []
                for m in range(nmaster):
                    for gi in posi[part[m]:part[m + 1]]:
                        row.append(vals[(gi, m)] if gi < nt else 1)
                ps.append(coord)
                data.append(row)
            write_fab_rows(os.path.join(dbdir, "red_ps"), ps)
            write_fab_rows(os.path.join(dbdir, "red_data"), data)
            common[-1] = len(ps)
            write_fab_rows(os.path.join(dbdir, "red_common"), [common])

    return grow


# --------------------------------------------------------------------------
# probe job scaffolding + the probe itself (SearchRelation.wl:288-296)
# --------------------------------------------------------------------------

def make_probe_job(searchdir: str, probe_job: str, src_job: str,
                   works: Sequence[int]) -> None:
    """Fresh job dir with the SAME schemes as src_job (BLGenerateJob[name]
    regenerates the family's scheme; the scheme is ansatz-independent, so a
    byte-copy of sch_* is the faithful equivalent).  Refuses to overwrite."""
    jobdir = os.path.join(searchdir, probe_job)
    if os.path.exists(jobdir):
        raise BladeGateError(f"make_probe_job: {jobdir} exists -- refusing")
    for w in works:
        Scheme.read(os.path.join(searchdir, src_job, str(w))).write(
            os.path.join(jobdir, str(w)))


def single_var_probe(searchdir: str, src_job: str, works: Sequence[int],
                     work_weights: List[List[int]],
                     search_params: Sequence[str], var_index: int,
                     intmode: str, kc: KiraCoeffs, seed: int,
                     nmax: int = 60, break_rounds: int = 60,
                     round_timeout: float = 900.0, nthreads: int = 8,
                     src_dataname: str = "data_p0",
                     log: Callable[[str], None] = print
                     ) -> Tuple[int, List[WorkOutcome], str]:
    """One adaptive-weight probe (SearchRelation.wl:294):
    BLSearchRelationsBasic[name, name, {{v}}, {{1}}, {intmode}, {},
    "SearchOnly"->True] with a frozen-conservative database.

    Returns (powerlimit = getLowerLimit(name), outcomes, probe_job).  The
    probe job/database are left on disk for audit but must NEVER be fitted
    or exported (their relations hold only at the frozen point)."""
    v = search_params[var_index]
    probe_job = f"probe_{v}"
    make_probe_job(searchdir, probe_job, src_job, works)
    init_probe_database(searchdir, probe_job, src_dataname)
    fam = FamilyAnsatz(search_params=tuple(search_params),
                       var_groups=[[v]], var_weights=[[1]],
                       intmode=[intmode], cut=[], nmax=nmax)
    cfg = SearchFamily(
        searchdir=searchdir, job=probe_job, fam=fam,
        work_weights=work_weights, datanames=[probe_job], nthreads=nthreads,
        ansatz_size=7,            # ansatzSize[{}] = 7  (SearchRelation.wl:140)
        level_base=0, level_step=1,
        data_size=10,             # dataSize for 1 var  (SearchRelation.wl:139)
        break_rounds=break_rounds,
        datasize_policy="structural", structural_margin=2.5,
        round_timeout=round_timeout, log=log)
    cfg.grower = make_frozen_grower(cfg, kc, seed=seed, free_idx=var_index,
                                    pids=(0,))
    outcomes = [escalate_work(cfg, w) for w in works]
    pl = get_lower_limit(os.path.join(searchdir, probe_job))
    return pl, outcomes, probe_job


def adaptive_weights(searchdir: str, src_job: str, works: Sequence[int],
                     work_weights: List[List[int]],
                     search_params: Sequence[str], intmode: str,
                     kc: KiraCoeffs, seed: int,
                     log: Callable[[str], None] = print, **probe_kw
                     ) -> Tuple[List[int], List[int], Dict[str, list]]:
    """The 'adaptive' branch of BLSearchRelations (SearchRelation.wl:288-305)
    for ONE group of all search params: probe every variable, then
    powerlimit2weight.  Returns (weights, powerlimits, outcome-dict).
    A probe that fails to close is a HARD error (upstream: $Failed)."""
    limits: List[int] = []
    all_out: Dict[str, list] = {}
    for i, v in enumerate(search_params):
        pl, outs, pjob = single_var_probe(
            searchdir, src_job, works, work_weights, search_params, i,
            intmode, kc, seed, log=log, **probe_kw)
        closed = all(o.closed for o in outs)
        log(f"[probe {v}] powerlimit={pl} closed={closed} "
            f"({', '.join(f'w{o.work}:{o.stop_reason or 'closed'}' for o in outs)})")
        if not closed:
            raise BladeGateError(
                f"probe {pjob} did not close all works -- upstream would "
                f"$Failed (SearchRelation.wl:294 via wl:201); "
                f"outcomes: {[(o.work, o.stop_reason) for o in outs]}")
        limits.append(pl)
        all_out[v] = [(o.work, o.closed,
                       [(r.level, round(r.wall, 2)) for r in o.rounds])
                      for o in outs]
    weights = powerlimit2weight(limits)
    log(f"[adaptive] Variables: {list(search_params)}, powerlimits: {limits},"
        f" Weights: {weights}  (SearchRelation.wl:305)")
    return weights, limits, all_out

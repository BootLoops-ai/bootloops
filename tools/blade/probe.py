"""blade.probe -- per-point numeric reduction oracle.

BTOracle(family_dir) loads a BANKED block-triangular (BT) form -- the
exact-rational relations exported by search.export_relations plus the block
scheme metadata -- and reduces target integrals to master-coefficient VALUES
mod p at numeric kinematic points, without Kira, without fits, without
reconstruction.  This is the "numeric-probe mode": the BT form is the one-time
asset; per point we only (1) evaluate the ansatz monomials, (2) solve the
tiny g1 x g1 block systems in block-triangular sequence.

family_dir layout (built by a bank step such as `blade_family prepare`):

    btform.json          -- the manifest (see BTForm.load for the schema)
    relations*.json      -- exact-rational relation exports (sha256-pinned
                            in btform.json["provenance"])
    run_tree/...         -- OPTIONAL full search tree (schemes, fit tables,
                            kinematics, database headers): enables the
                            ssolve-driven backend

Backends
--------
* exact  (default): pure-python GF(p) solve of the banked EXACT-RATIONAL
  relations.  Works at ANY prime (fresh primes included) and on partial
  families (closed blocks only; unknowns of unclosed blocks become declared
  independents).  This is the fresh-prime verification lane.
* ssolve: drives the shipped C stage `ssolve` on the banked run tree's
  per-prime FIT TABLES (the pipeline path fitrel->ssolve pays for during the
  search).  Restricted to fflow primes with banked fit tables; artifact gates
  follow blade.pipeline (rc never gates; bytes/dimensions do).

Gate discipline: every load/solve failure raises BladeGateError naming the
artifact; outputs are cross-gated (ssolve backend self-gates its element
ordering against the exact backend on every batch).
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

from .formats import (Database, EvalFile, PointsFile, Scheme, SystemFile,
                      big_uint_primes, flags_size)
from .pipeline import BladeGateError, DEFAULT_BIN_DIR


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------
# BT form loading
# --------------------------------------------------------------------------

@dataclass
class _Work:
    work: int
    g1_gids: List[int]                       # unknowns of this block
    relations: List[dict] = field(default_factory=list)
    # compiled: per relation, list of (gid, Fraction, exps tuple over params)
    terms: List[List[Tuple[int, Fraction, Tuple[int, ...]]]] = field(
        default_factory=list)


class BTForm:
    """btform.json schema (all gids are family-global integral ids):

    {
      "family": "db",
      "params": ["t", "eps"],              # BLSearchParameter order, eps LAST
      "labels": ["db[...]", ...],          # index == gid
      "target_gids": [0,1,2,3],            # the family's declared targets
      "master_gids": [4,...,11],
      "works": [{"work": 0, "g1_gids": [0,1,2,3]}, ...],   # CLOSED works only
      "relations_files": ["relations_w0.json", ...],
      "provenance": {"<file>": "<sha256>", ...},
      "ssolve": {                          # optional, enables ssolve backend
         "search_tree": "run_tree/search", "job": "job0",
         "datanames": ["data_p0", ...],    # index == fflow prime id
         "nworks": 1
      }
    }
    """

    def __init__(self, family_dir: str):
        self.dir = os.path.abspath(family_dir)
        mpath = os.path.join(self.dir, "btform.json")
        if not os.path.exists(mpath):
            raise BladeGateError(f"btform.json missing in {self.dir}")
        self.meta = json.load(open(mpath))
        for key in ("family", "params", "labels", "target_gids",
                    "master_gids", "works", "relations_files", "provenance"):
            if key not in self.meta:
                raise BladeGateError(f"btform.json: key {key!r} missing")
        self.family: str = self.meta["family"]
        self.params: List[str] = list(self.meta["params"])
        self.labels: List[str] = list(self.meta["labels"])
        self.target_gids: List[int] = list(self.meta["target_gids"])
        self.master_gids: List[int] = list(self.meta["master_gids"])

        # ---- provenance gate: every relations file byte-pinned -------------
        self.works: Dict[int, _Work] = {
            int(w["work"]): _Work(int(w["work"]), list(w["g1_gids"]))
            for w in self.meta["works"]}
        unknown_gids = set()
        for w in self.works.values():
            if unknown_gids & set(w.g1_gids):
                raise BladeGateError("btform: works share unknown gids")
            unknown_gids |= set(w.g1_gids)
        for rf in self.meta["relations_files"]:
            path = os.path.join(self.dir, rf)
            want = self.meta["provenance"].get(rf)
            if want is None:
                raise BladeGateError(f"btform: no provenance sha for {rf}")
            got = _sha256(path)
            if got != want:
                raise BladeGateError(
                    f"{path}: sha256 {got[:16]}... != banked {want[:16]}... "
                    f"(BT form corrupted or swapped)")
            doc = json.load(open(path))
            ra = doc.get("reduced_analytic")
            if ra is not None:
                raise BladeGateError(
                    f"{path}: reduced-analytic export (nana={ra.get('nana')}, "
                    f"pinned {ra.get('pinned_params')}) -- these relations "
                    f"hold ONLY on the pinned slice; a full-symbolic BT form "
                    f"must not be built from them (FAIL-CLOSED).  Per-point "
                    f"fitting in the pinned parameters is the consumer's job "
                    f"(search.plan_pin_instances / search.fit_at_pin).")
            if list(doc["search_params"]) != self.params:
                raise BladeGateError(f"{path}: search_params {doc['search_params']}"
                                     f" != btform params {self.params}")
            for rel in doc["relations"]:
                w = int(rel["work"])
                if w not in self.works:
                    raise BladeGateError(f"{path}: relation for undeclared work {w}")
                self.works[w].relations.append(rel)

        # ---- structural gates + compile ------------------------------------
        pindex = {p: i for i, p in enumerate(self.params)}
        for w in self.works.values():
            if len(w.relations) != len(w.g1_gids):
                raise BladeGateError(
                    f"work {w.work}: {len(w.relations)} relations != g1 = "
                    f"{len(w.g1_gids)} (block system not square)")
            for rel in w.relations:
                if rel["n_terms"] != len(rel["terms"]):
                    raise BladeGateError(f"work {w.work}: n_terms mismatch")
                compiled = []
                for t in rel["terms"]:
                    gid = int(t["global_id"])
                    if not (0 <= gid < len(self.labels)) \
                            or t["integral"] != self.labels[gid]:
                        raise BladeGateError(
                            f"work {w.work}: term gid {gid} / label "
                            f"{t['integral']!r} inconsistent with btform labels")
                    exps = [0] * len(self.params)
                    for pn, e in t["monomial"].items():
                        exps[pindex[pn]] = int(e)
                    compiled.append((gid, Fraction(t["coeff"]), tuple(exps)))
                w.terms.append(compiled)

        # ---- independents & solve order -------------------------------------
        solved_unknowns = unknown_gids
        referenced = {gid for w in self.works.values()
                      for rel in w.terms for gid, _, _ in rel}
        self.independent_gids: List[int] = list(self.master_gids) + sorted(
            g for g in referenced
            if g not in solved_unknowns and g not in self.master_gids)
        self.extra_independents = self.independent_gids[len(self.master_gids):]
        # order works so every referenced solved-elsewhere gid is ready
        order, ready = [], set(self.independent_gids)
        pending = dict(self.works)
        while pending:
            for wk in sorted(pending):
                w = pending[wk]
                need = {g for rel in w.terms for g, _, _ in rel
                        if g not in w.g1_gids}
                if need <= ready:
                    order.append(w)
                    ready |= set(w.g1_gids)
                    del pending[wk]
                    break
            else:
                raise BladeGateError(
                    f"BT works not triangularizable; pending {sorted(pending)}")
        self.solve_order: List[_Work] = order
        self.closed_gids: List[int] = [g for w in order for g in w.g1_gids]

    def label(self, gid: int) -> str:
        return self.labels[gid]


# --------------------------------------------------------------------------
# the oracle
# --------------------------------------------------------------------------

class BTOracle:
    """Per-point reduction oracle over a banked BT form.

    reduce_point(kin_point, eps, prime) -> {target_label: {basis_label: value}}
    where basis = masters (+ declared extra independents for partial families).
    Zero coefficients are omitted; values are ints in [0, prime).
    """

    def __init__(self, family_dir: str, bin_dir: str = DEFAULT_BIN_DIR):
        t0 = time.time()
        self.form = BTForm(family_dir)
        self.bin_dir = bin_dir
        self._mod_cache: Dict[int, list] = {}   # prime -> per work/rel coeffs
        self.load_wall = time.time() - t0

    # ---- exact backend ------------------------------------------------------

    def _coeffs_mod(self, p: int):
        """Reduce every rational coefficient mod p once per prime."""
        cached = self._mod_cache.get(p)
        if cached is not None:
            return cached
        out = {}
        for w in self.form.solve_order:
            rels = []
            for rel in w.terms:
                row = []
                for gid, fr, exps in rel:
                    den = fr.denominator % p
                    if den == 0:
                        raise BladeGateError(
                            f"prime {p} divides a coefficient denominator in "
                            f"work {w.work} -- pick another prime")
                    row.append((gid, fr.numerator % p * pow(den, -1, p) % p,
                                exps))
                rels.append(row)
            out[w.work] = rels
        self._mod_cache[p] = out
        return out

    def _coords(self, kin_point: Sequence[int], eps: int) -> List[int]:
        if len(kin_point) != len(self.form.params) - 1:
            raise BladeGateError(
                f"kin_point has {len(kin_point)} coords; family params are "
                f"{self.form.params} (eps passed separately)")
        if self.form.params[-1] != "eps":
            raise BladeGateError("btform params do not end in 'eps'")
        return list(kin_point) + [eps]

    def reduce_point(self, kin_point: Sequence[int], eps: int, prime: int,
                     want_gids: Optional[Sequence[int]] = None) -> dict:
        return self.reduce_many([(kin_point, eps)], prime, want_gids)[0]

    def reduce_many(self, points: Sequence[Tuple[Sequence[int], int]],
                    prime: int,
                    want_gids: Optional[Sequence[int]] = None) -> List[dict]:
        """Batch reduction (the eta-grid use case).  points = [(kin, eps)]."""
        p = prime
        form = self.form
        mods = self._coeffs_mod(p)
        basis = form.independent_gids
        bidx = {g: i for i, g in enumerate(basis)}
        nb = len(basis)
        if want_gids is None:
            want_gids = [g for g in form.target_gids if g in form.closed_gids]
            if not want_gids:
                want_gids = form.closed_gids
        for g in want_gids:
            if g not in form.closed_gids:
                raise BladeGateError(
                    f"gid {g} ({form.label(g)}) is not closed by the banked "
                    f"BT form (closed: {form.closed_gids})")
        results = []
        for kin, eps in points:
            coords = [c % p for c in self._coords(kin, eps)]
            vec: Dict[int, List[int]] = {}
            for w in form.solve_order:
                g1 = w.g1_gids
                uidx = {g: i for i, g in enumerate(g1)}
                n = len(g1)
                A = [[0] * n for _ in range(n)]
                B = [[0] * nb for _ in range(n)]
                for r, rel in enumerate(mods[w.work]):
                    for gid, c, exps in rel:
                        m = c
                        for x, e in zip(coords, exps):
                            if e:
                                m = m * pow(x, e, p) % p
                        ui = uidx.get(gid)
                        if ui is not None:
                            A[r][ui] = (A[r][ui] + m) % p
                        elif gid in bidx:
                            B[r][bidx[gid]] = (B[r][bidx[gid]] - m) % p
                        else:
                            vg = vec.get(gid)
                            if vg is None:
                                raise BladeGateError(
                                    f"work {w.work}: term gid {gid} neither "
                                    f"unknown, independent, nor solved")
                            Br = B[r]
                            for j in range(nb):
                                if vg[j]:
                                    Br[j] = (Br[j] - m * vg[j]) % p
                # Gauss-Jordan on [A | B] mod p
                for col in range(n):
                    piv = next((r for r in range(col, n) if A[r][col]), None)
                    if piv is None:
                        raise BladeGateError(
                            f"work {w.work}: block matrix singular mod {p} at "
                            f"point {coords} (bad point/prime -- resample)")
                    A[col], A[piv] = A[piv], A[col]
                    B[col], B[piv] = B[piv], B[col]
                    inv = pow(A[col][col], -1, p)
                    A[col] = [a * inv % p for a in A[col]]
                    B[col] = [b * inv % p for b in B[col]]
                    for r in range(n):
                        if r != col and A[r][col]:
                            f = A[r][col]
                            A[r] = [(a - f * ac) % p
                                    for a, ac in zip(A[r], A[col])]
                            B[r] = [(b - f * bc) % p
                                    for b, bc in zip(B[r], B[col])]
                for g in g1:
                    vec[g] = B[uidx[g]]
            out = {}
            for g in want_gids:
                out[form.label(g)] = {
                    form.label(basis[j]): vec[g][j]
                    for j in range(nb) if vec[g][j]}
            results.append(out)
        return results

    # ---- ssolve backend (drives the shipped C stage on banked fit tables) ---

    def _ssolve_meta(self):
        ss = self.form.meta.get("ssolve")
        if not ss:
            raise BladeGateError(
                f"{self.form.dir}: btform has no 'ssolve' section (fit-table "
                f"run tree not banked for this family)")
        tree = os.path.join(self.form.dir, ss["search_tree"])
        return ss, tree

    def reduce_many_ssolve(self, points: Sequence[Tuple[Sequence[int], int]],
                           pid: int, scratch: str, nthreads: int = 1,
                           self_gate: bool = True) -> List[dict]:
        """Reduce a batch by running `ssolve` on the banked per-prime fit
        tables (block sequence).  prime = BIG_UINT_PRIMES[pid]; the banked
        tree must hold fit tables for datanames[pid].  Artifact gates: eval
        byte size / header / per-row coord+prime echo; element ordering is
        self-gated against the exact backend on the first point of every
        batch (self_gate=False only for pure timing runs)."""
        ss, tree = self._ssolve_meta()
        form = self.form
        if not (0 <= pid < len(ss["datanames"])):
            raise BladeGateError(f"pid {pid} has no banked dataname")
        dataname = ss["datanames"][pid]
        jobdir = os.path.join(tree, ss["job"])
        # header/posi/part read raw: red_common nps may be LOGICALLY shrunk by
        # fitnumber bisection (set_data_number) -- Database.read refuses the
        # physical-pool mismatch (documented footgun, search.make_kira_grower)
        from .formats import read_ints
        dbdir = os.path.join(tree, "database", dataname)
        common = read_ints(os.path.join(dbdir, "red_common"))
        db_prime, db_npara, db_nint, db_nmaster = common[0], common[1], \
            common[2], common[3]
        db_part = read_ints(os.path.join(dbdir, "red_part"))
        db_posi = read_ints(os.path.join(dbdir, "red_posi"))
        P = big_uint_primes()
        if db_prime != P[pid]:
            raise BladeGateError(f"{dataname}: red_common prime {db_prime} != "
                                 f"BIG_UINT_PRIMES[{pid}]={P[pid]}")
        nt = db_nint - db_nmaster
        # stored (target, master) elements, ssolve output order (target-major)
        elements = sorted(
            (gi, m) for m in range(db_nmaster)
            for gi in db_posi[db_part[m]:db_part[m + 1]] if gi < nt)
        npo = len(elements)
        fs = flags_size(npo)
        nparsin = len(form.params)
        os.makedirs(scratch, exist_ok=True)
        # points file: coords, prime, all-outputs-needed flags
        mask = (1 << npo) - 1
        fwords = [(mask >> (64 * i)) & ((1 << 64) - 1) for i in range(fs)]
        rows = []
        for kin, eps in points:
            coords = [c % db_prime for c in self._coords(kin, eps)]
            rows.append(coords + [db_prime] + fwords)
        pf = PointsFile(n=nparsin + fs, rows=rows)
        ppath = os.path.join(scratch, f"probe_points_p{pid}.fflow")
        pf.write(ppath)
        # system file (absolute blockdir; kinematics resolved relative to it)
        sysf = SystemFile(db_prime, db_npara, db_nint, db_nmaster,
                          targets=[g for w in form.solve_order
                                   for g in w.g1_gids],
                          blockdirs=[os.path.abspath(jobdir)],
                          dataname=dataname)
        spath = os.path.join(scratch, f"probe_system_p{pid}.txt")
        sysf.write(spath)
        epath = os.path.join(scratch, f"probe_eval_p{pid}.fflow")
        if os.path.exists(epath):
            os.unlink(epath)
        ptsord = "".join(str(i) for i in range(1, nparsin + 1))
        argv = [os.path.join(self.bin_dir, "ssolve"), str(nparsin), str(pid),
                "0", str(len(rows)), str(nthreads), spath, ptsord, ppath,
                epath]
        cp = subprocess.run(argv, cwd=scratch, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=3600)
        log = cp.stdout.decode(errors="replace")
        if not os.path.exists(epath):
            raise BladeGateError(f"ssolve: {epath} missing (rc={cp.returncode}"
                                 f", rc lies anyway). log tail: {log[-400:]}")
        want = EvalFile.expected_bytes(nparsin, npo, len(rows))
        got = os.path.getsize(epath)
        if got != want:
            raise BladeGateError(f"ssolve: {epath}: {got} bytes != {want} "
                                 f"expected ({len(rows)} pts x {npo} outputs)")
        ev = EvalFile.read(epath)
        if ev.nparsin != nparsin or ev.nout_plus_flags != npo + fs \
                or len(ev.rows) != len(rows):
            raise BladeGateError(f"ssolve: {epath}: header inconsistent")
        results = []
        for i, r in enumerate(ev.rows):
            if r[:nparsin] != rows[i][:nparsin] or r[nparsin] != db_prime:
                raise BladeGateError(
                    f"ssolve: {epath}: row {i} coords/prime echo mismatch")
            outs = r[nparsin + 1 + fs:]
            out: Dict[str, Dict[str, int]] = {}
            for (gi, m), v in zip(elements, outs):
                if v:
                    out.setdefault(form.label(gi), {})[
                        form.label(db_nint - db_nmaster + m)] = v
            results.append(out)
        if self_gate and points:
            exact = self.reduce_many([points[0]], db_prime)[0]
            got0 = {k: results[0].get(k, {}) for k in exact}
            if got0 != exact:
                raise BladeGateError(
                    "ssolve backend disagrees with the exact-relation backend "
                    "on the batch's first point -- element-order or fit-table "
                    "corruption; refusing to return values")
        return results


__all__ = ["BTForm", "BTOracle", "BladeGateError"]

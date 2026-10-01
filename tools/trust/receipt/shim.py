#!/usr/bin/env python3
"""shim.py — minimal stable API for live runs (coordinated READ-ONLY;
this tool never writes into a live run directory).

Import surface (kept deliberately tiny and stable):

    sys.path.insert(0, "<path to>/tools/trust/receipt")
    from shim import ReceiptSession

    rs = ReceiptSession.from_strata(art_dir, staging_dir, family,
                                    p, d0, eta0)          # ~seconds
    rec, path = rs.emit(target_col, claimed=c_dict,       # ~2.3 s/row at
                        out_dir=my_out)                   #  reference scale
    ok, detail = rs.verify(path)                          # ~1 ms/row
    ok, detail = rs.verify(path, table_row=claimed_c)     # wrong-table mode
    alarm, report = rs.audit(targets=[...])               # structural audit

Or construct generically (any solver, no kira anywhere):

    rs = ReceiptSession(rows, pivots, masters, p, point={"d":..,"eta":..})

COST HONESTY (measured): emission is Wiedemann at ~n matvecs/target —
2.2-2.3 s/(row,prime) at the 7,135-equation reference scale, and a
MEASURED-THROUGHPUT LOWER BOUND of ~24 min/(target,prime) at large scale
(n >= 179,447). On a 397,813-equation system, per-row emission is for
PRE-REGISTERED SAMPLES of rows, not whole tables; whole-table certification
at that scale is priced out. Verification stays ~1 ms/row at any emitted
witness.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import core                                    # noqa: E402
from core import load_witness, verify_row, write_witness   # noqa: E402
from detector import detect                    # noqa: E402
from emitter import FpSystem, emit_for_target  # noqa: E402


class ReceiptSession:
    """One (system, point, prime) context: emit / verify / audit."""

    def __init__(self, rows, pivots, masters, p, point=None, family=None,
                 labels=None, shell=(), seed=20260706):
        """rows: list of dict col->val (mod p) — THE row order of record
        (lambda indexes it). pivots: square pivot column list (emit only;
        pass [] for verify/audit-only sessions). masters: master column set.
        labels: optional col->indices dict for human-readable witnesses."""
        import numpy as np
        self.rows = rows
        self.p = int(p)
        self.point = dict(point or {})
        self.family = family
        self.masters = set(masters)
        self.shell = set(shell)
        self.labels = dict(labels or {})
        self._rng = np.random.default_rng(seed)
        self._minpoly = None
        self._S = None
        if pivots:
            allcols = sorted({c for r in rows for c in r})
            self._S = FpSystem([list(r.items()) for r in rows],
                               list(pivots), allcols, self.p)
        self._fpr = None

    # ------------------------------------------------------------ builders
    @classmethod
    def from_strata(cls, art_dir, staging_dir, family, p, d0, eta0,
                    with_emit=True, seed=20260706):
        """Build from a kira/STRATA staged artifact dir (adapter path).
        with_emit=False skips the pivot-partition staging pass (verify-only
        sessions load in ~0.2 s at production registry scale)."""
        from adapters import strata
        sysd = strata.load_rows(art_dir, family, p, d0, eta0)
        if with_emit:
            _s, info, _c = strata.stage_partition(
                art_dir, staging_dir, family, p, d0, eta0,
                sysd["sector_of"])
            masters = set(info["masters_map"])
            labels = {w: list(idx) for w, idx in info["masters_map"].items()}
            shell = set(info["shell"])
            pivots = info["pivots"]
        else:
            masters = set(sysd["masters"])
            labels = {w: (m["indices"] or [])
                      for w, m in sysd["masters"].items()}
            shell, pivots = set(), []
        obj = cls(sysd["rows"], pivots, masters, p,
                  point={"d": d0, "eta": eta0}, family=family,
                  labels=labels, shell=shell, seed=seed)
        obj._sysd = sysd
        return obj

    # ---------------------------------------------------------------- emit
    def emit(self, target_col, claimed=None, out_dir=None, target_label=None):
        """Solve + write one witness. Returns (record, path or None).
        record["status"]: CERTIFIED | CLAIM_MISMATCH | SHELL_RESIDUAL |
        SOLVE_FAILED. CLAIM_MISMATCH = the claimed row is NOT implied by the
        system (wrong-table caught at emit time); no witness is written."""
        import numpy as np
        assert self._S is not None, "session built without pivots (emit off)"
        rec, lam, mp_new = emit_for_target(
            self._S, target_col, self.masters, self._rng,
            minpoly=self._minpoly, claimed=claimed)
        if mp_new is not None:
            self._minpoly = mp_new
        path = None
        if lam is not None and rec["status"] == "CERTIFIED" and out_dir:
            if self._fpr is None:
                self._fpr = core.system_fingerprint(self.rows, self.p)
            nz = np.flatnonzero(lam)
            fam = self.family or "generic"
            path = os.path.join(out_dir, f"w_{fam}_{self.p}_{target_col}.json")
            write_witness(
                path, p=self.p, target_col=target_col, c=rec["c"],
                lam_idx=[int(i) for i in nz],
                lam_val=[int(lam[i]) for i in nz],
                n_rows=len(self.rows), family=self.family, point=self.point,
                target_label=target_label,
                labels={w: self.labels[w] for w in rec["c"]
                        if w in self.labels},
                system={"n_rows": len(self.rows), "fingerprint": self._fpr})
            rec["witness"] = path
        return rec, path

    # -------------------------------------------------------------- verify
    def verify(self, witness, table_row=None, check_fingerprint=False):
        """witness: path or dict from core.load_witness. ~1 ms/row after the
        session holds the rows. Solver-agnostic (core.verify_row)."""
        wit = load_witness(witness) if isinstance(witness, str) else witness
        fpr = None
        if check_fingerprint:
            if self._fpr is None:
                self._fpr = core.system_fingerprint(self.rows, self.p)
            fpr = self._fpr
        return verify_row(self.rows, wit, table_row=table_row,
                          check_fingerprint=fpr)

    # --------------------------------------------------------------- audit
    def audit(self, targets=(), table=None):
        """Structural mis-reduction detector (see detector.py)."""
        return detect(self.rows, self.p, self.masters, shell=self.shell,
                      targets=targets, table=table)

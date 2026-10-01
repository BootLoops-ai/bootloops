#!/usr/bin/env python3
"""Certify-after-float hull membership + certificate reuse.

Decision chain (float64 never decides, every accepted answer is re-verified
in exact arithmetic):

  0. REUSE (warm starts across a refinement sequence / related LPs on the
     same hull):
       - last verified Farkas functional: ONE exact re-substitution over the
         full point set (shift trick) re-certifies it for the new query;
       - last verified basis (plus t-neighbor columns): ONE exact dual-simplex
         call restricted to those columns — for a nearby query this is a
         single exact basis solve; FEASIBLE on a subset is globally valid.
  1. PROPOSE: float64 LP (scipy HiGHS, phase-1 form — PROPOSER ONLY)
     proposes a support basis (feasible candidate) or a dual Farkas
     direction (infeasible candidate).
  2. VERIFY: infeasible candidate — one exact re-substitution of the proposed
     functional over the FULL point set (shift trick: with mn = min_i F(P_i)
     computed exactly, F - mn is a valid exact certificate independent of
     float error iff F(q) < mn exactly). Feasible candidate — exact solve on
     the proposed columns via stage 3 (float64 rank-collapses near-duplicate
     curve columns, so the raw support is often exactly rank-deficient; the
     subset simplex is the sound completion of the one-shot solve).
  3. SUBSET SIMPLEX: the exact dual-simplex (cone_lp_b) restricted to
     float-priced candidate columns; FEASIBLE on a subset is globally valid;
     INFEASIBLE is accepted only after an exact full-point-set Farkas scan
     (violating columns join the subset and the LP re-runs).
  4. FALLBACK: the unrestricted exact dual-simplex path, warm-seeded.

Exactness statement: every FEASIBLE answer returns lambda with lambda >= 0,
sum lambda = 1, sum lambda_i P_i = q re-verified in gmpy2.mpq; every
INFEASIBLE answer returns (w, w0) with <w,P_i> + w0 >= 0 for ALL points i and
<w,q> + w0 < 0 re-verified in gmpy2.mpq against the full point set.
Unverified proposals are discarded, never reported.
"""
import os
import sys
import numpy as np
from gmpy2 import mpq

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import cone_lp_b as cone_lp                   # the exact dual-simplex engine
                                              # (rref basis completion — see
                                              # cone_lp_b.py header note)
from cone_lp_b import LPStall, col_from_point # noqa: F401 (LPStall re-exported)
from fast_lp import _rationalize              # float -> small-denominator mpq

from scipy.optimize import linprog            # PROPOSER ONLY


def _F(w, w0, p):
    return sum((w[i] * p[i] for i in range(len(p))), mpq(0)) + w0


class CertHull:
    """Growing hull of exact simplex points, certify-after-float membership.

    decide(q) -> (status, payload, basis_or_None, route)
      'FEASIBLE':   payload = lam (full-length list of mpq), exact-verified
      'INFEASIBLE': payload = (w, w0), exact-verified affine Farkas functional
    route in {'reuse-basis','reuse-farkas','float-farkas','subset-simplex',
              'exact-fallback'}
    """

    SUBSET_ROUND_CAP = 60
    VIOL_BATCH = 8

    def __init__(self, points):
        self.points = []
        self._A_flt = None            # cached float64 matrix (d+1, N)
        self._icols = {}              # j -> (integer cone column, scale)
        self.last_basis = None        # verified feasible support (point indices)
        self.last_farkas = None       # {'w','w0','mn','k'} verified Farkas
        self.last_kind = None
        self.stats = {'reuse_basis': 0, 'reuse_farkas': 0, 'float_farkas': 0,
                      'subset_simplex': 0, 'fallback': 0, 'float_calls': 0,
                      'subset_rounds': 0}
        self.add_points(points)

    # ------------------------------------------------------------------ #
    def add_points(self, pts):
        self.points.extend(pts)
        self._A_flt = None

    def _icol(self, j):
        if j not in self._icols:
            self._icols[j] = col_from_point(self.points[j])
        return self._icols[j]

    def _sub_membership(self, cand, q, warm=None, log=None):
        """Exact phase-0 membership on the candidate subset (memoized integer
        columns). warm carries GLOBAL point indices and is mapped to subset
        positions. Returns cone_lp.membership_exact output on the subset."""
        sub = [self.points[j] for j in cand]
        cols = []
        scales = []
        pos = {}
        for k, j in enumerate(cand):
            g, D = self._icol(j)
            cols.append(g)
            scales.append(D)
            pos[j] = k
        warm_local = [pos[j] for j in (warm or []) if j in pos]
        return cone_lp.membership_exact(sub, q, warm=warm_local, log=log,
                                        _cols=cols, _scales=scales)

    # ------------------------- exact verification ---------------------- #
    def _accept_feasible(self, support, lam_s, q):
        """Independent exact re-substitution in gmpy2.mpq (re-substitution law); records the
        reuse basis. support/lam_s indexed over global point indices."""
        d = len(q)
        assert all(l >= 0 for l in lam_s)
        assert sum(lam_s) == 1
        for i in range(d):
            assert sum((lam_s[k] * self.points[j][i]
                        for k, j in enumerate(support)), mpq(0)) == mpq(q[i])
        lam = [mpq(0)] * len(self.points)
        for k, j in enumerate(support):
            lam[j] = lam_s[k]
        self.last_basis = list(support)
        self.last_kind = 'FEASIBLE'
        return lam

    def _verify_basis(self, idxs, q, log=None):
        """Reuse path: ONE exact dual-simplex call restricted to the stored
        basis (single exact basis solve when the old basis still carries q),
        widened to t-neighbor columns on failure. Returns full-length lam
        list or None."""
        for cand in (sorted(set(idxs)), sorted(self._neighbors(idxs))):
            if not cand:
                continue
            try:
                st, res, _basis = self._sub_membership(cand, q,
                                                       warm=list(idxs), log=log)
            except LPStall:
                continue
            if st != 'FEASIBLE':
                continue
            support = [cand[k] for k in range(len(cand)) if res[k] > 0]
            lam_s = [res[k] for k in range(len(cand)) if res[k] > 0]
            return self._accept_feasible(support, lam_s, q)
        return None

    def _verify_farkas(self, cert, q):
        """Shift-trick exact verification of a proposed affine functional
        against the FULL point set. cert = {'w','w0','mn','k'}; mn/k cache the
        exact minimum over points[:k]. Returns (w, w0_shifted) or None."""
        w, w0 = cert['w'], cert['w0']
        mn, k = cert['mn'], cert['k']
        for j in range(k, len(self.points)):
            v = _F(w, w0, self.points[j])
            if mn is None or v < mn:
                mn = v
        cert['mn'], cert['k'] = mn, len(self.points)
        if not (_F(w, w0, q) < mn):
            return None
        w0x = w0 - mn
        # ---- independent exact re-check, full point set ----
        assert _F(w, w0x, q) < 0
        for p in self.points:
            assert _F(w, w0x, p) >= 0
        self.last_farkas = cert
        self.last_kind = 'INFEASIBLE'
        return (w, w0x)

    # --------------------------- float proposer ------------------------ #
    def _float_propose(self, q):
        """scipy HiGHS phase-1 on the float64 shadow problem. Returns
        ('feas', support_idx_list) | ('infeas', (w, w0)) | None."""
        d = len(q)
        m = d + 1
        N = len(self.points)
        self.stats['float_calls'] += 1
        if self._A_flt is None or self._A_flt.shape[1] != N:
            self._A_flt = np.array(
                [[float(p[i]) for p in self.points] for i in range(d)]
                + [[1.0] * N])
        A = self._A_flt.copy()
        b = np.array([float(v) for v in q] + [1.0])
        sgn = [1] * m
        for i in range(m):
            if b[i] < 0:
                sgn[i] = -1
                b[i] = -b[i]
                A[i] = -A[i]
        c = np.concatenate([np.zeros(N), np.ones(m)])
        A_eq = np.hstack([A, np.eye(m)])
        try:
            res = linprog(c, A_eq=A_eq, b_eq=b, bounds=(0, None), method='highs')
        except Exception:
            return None
        if res.status != 0 or res.x is None:
            return None
        if res.fun < 1e-9:
            x = res.x[:N]
            supp = [int(j) for j in np.argsort(-x) if x[j] > 1e-13][:m]
            return ('feas', supp) if supp else None
        eq = getattr(res, 'eqlin', None)
        y = None if eq is None else eq.marginals
        if y is None:
            return None
        yu = [_rationalize(float(sgn[i] * y[i])) for i in range(m)]
        scale = max(abs(v) for v in yu)
        if scale == 0:
            return None
        yu = [v / scale for v in yu]
        w = [-yu[i] for i in range(d)]
        w0 = -yu[d]
        return 'infeas', (w, w0)

    # ------------------- subset dual simplex (stage 3) ------------------ #
    def _neighbors(self, idxs):
        out = set()
        for j in idxs:
            out |= {j - 1, j, j + 1}
        return {j for j in out if 0 <= j < len(self.points)}

    def _spread(self):
        N = len(self.points)
        d1 = len(self.points[0]) if self.points else 1
        step = max(1, N // (d1 + 1))
        return set(range(0, N, step))

    def _subset_simplex(self, q, seed, log=None):
        """Exact cone dual simplex (phase-0 engine, read-only) on a candidate
        column subset; FEASIBLE on a subset is globally valid; INFEASIBLE is
        accepted only after an exact full-set Farkas scan (violators join the
        subset). Returns ('FEASIBLE', lam) | ('INFEASIBLE', (w, w0)) | None."""
        d = len(q)
        N = len(self.points)
        cand = sorted(self._neighbors(seed) | set(self.last_basis or []))
        spread_added = False
        if len(cand) < d + 1:     # cannot possibly carry a full basis
            cand = sorted(set(cand) | self._spread())
            spread_added = True
        warm = list(seed)
        for _round in range(self.SUBSET_ROUND_CAP):
            self.stats['subset_rounds'] += 1
            try:
                st, res, _basis = self._sub_membership(cand, q, warm=warm,
                                                       log=log)
            except LPStall:
                if not spread_added:
                    cand = sorted(set(cand) | self._spread())
                    spread_added = True
                    continue
                return None
            if st == 'FEASIBLE':
                support = [cand[k] for k in range(len(cand)) if res[k] > 0]
                lam_s = [res[k] for k in range(len(cand)) if res[k] > 0]
                return ('FEASIBLE', self._accept_feasible(support, lam_s, q))
            w, w0 = res
            viol = []
            for j in range(N):
                if j in cand:
                    continue
                if _F(w, w0, self.points[j]) < 0:
                    viol.append(j)
                    if len(viol) >= self.VIOL_BATCH:
                        break
            if not viol:
                # exact full-set Farkas re-check through the common gate
                out = self._verify_farkas(
                    {'w': w, 'w0': w0, 'mn': None, 'k': 0}, q)
                assert out is not None, 'subset Farkas failed global re-check'
                return ('INFEASIBLE', out)
            cand = sorted(set(cand) | set(viol))
        return None

    # ------------------------------ decide ----------------------------- #
    def decide(self, q, log=None):
        # 0. certificate reuse (verified answers from earlier related LPs)
        if self.last_kind == 'INFEASIBLE' and self.last_farkas is not None:
            out = self._verify_farkas(self.last_farkas, q)
            if out is not None:
                self.stats['reuse_farkas'] += 1
                return 'INFEASIBLE', out, None, 'reuse-farkas'
        if self.last_basis:
            out = self._verify_basis(self.last_basis, q, log=log)
            if out is not None:
                self.stats['reuse_basis'] += 1
                return 'FEASIBLE', out, list(self.last_basis), 'reuse-basis'
        if self.last_kind != 'INFEASIBLE' and self.last_farkas is not None:
            out = self._verify_farkas(self.last_farkas, q)
            if out is not None:
                self.stats['reuse_farkas'] += 1
                return 'INFEASIBLE', out, None, 'reuse-farkas'
        # 1.+2. float64 proposal; exact one-shot verification for Farkas
        prop = self._float_propose(q)
        seed = []
        if prop is not None:
            pkind, data = prop
            if pkind == 'feas':
                seed = list(data)
            else:
                w, w0 = data
                for wc, w0c in ((w, w0), ([-v for v in w], -w0)):
                    out = self._verify_farkas(
                        {'w': wc, 'w0': w0c, 'mn': None, 'k': 0}, q)
                    if out is not None:
                        self.stats['float_farkas'] += 1
                        return 'INFEASIBLE', out, None, 'float-farkas'
        # 3. exact dual simplex on float-priced candidate subset
        sub = self._subset_simplex(q, seed, log=log)
        if sub is not None:
            self.stats['subset_simplex'] += 1
            st, payload = sub
            basis = list(self.last_basis) if st == 'FEASIBLE' else None
            return st, payload, basis, 'subset-simplex'
        # 4. unrestricted phase-0 exact dual simplex, warm-seeded
        self.stats['fallback'] += 1
        if log:
            log('    cert_lp: subset path stalled -> full exact fallback')
        warm = list(seed) + list(self.last_basis or [])
        for j in range(len(self.points)):   # memoize all integer columns
            self._icol(j)
        cols = [self._icols[j][0] for j in range(len(self.points))]
        scales = [self._icols[j][1] for j in range(len(self.points))]
        st, res, basis = cone_lp.membership_exact(self.points, q, warm=warm,
                                                  log=log, _cols=cols,
                                                  _scales=scales)
        if st == 'FEASIBLE':
            self.last_basis = [j for j in range(len(res)) if res[j] > 0]
            self.last_kind = 'FEASIBLE'
            return st, res, list(self.last_basis), 'exact-fallback'
        w, w0 = res
        self.last_farkas = {'w': w, 'w0': w0, 'mn': None, 'k': 0}
        self.last_kind = 'INFEASIBLE'
        return st, res, None, 'exact-fallback'

"""boxwalk.relations — exact relation calculus + peel + prime-batched refill.

Math (pilot p3_staircase, verified vs dense oracle; dimension-generalized):
for a polynomial field a = (a_1..a_n), F_e := a_e*X_e with X_e = x_e - x_e^2
(tangent to the box boundary), and for each class j the exact division
    G_j := sum_e a_e X_e d_e P_j = q_j*P_j + r_j        (flint divmod)
gives, for every base monomial gamma and count vector u:

  0 = sum_delta [ D(delta) + sum_e gamma_e F_e(delta-e_e)
                  + sum_j u_j q_j(delta) ] * M(gamma+delta; u)
      + sum_{j: r_j != 0} u_j sum_delta r_j(delta) * M(gamma+delta; u-e_j)

with D = sum_e d_e F_e.  Tangency to j <=> r_j = 0 <=> no down-shift into
u-e_j.  Coefficients are affine in (gamma, u): record the pivot program ONCE,
replay at every (u, prime).

FamTables/peel/refill follow pilot p3d (vectorized) + p3b (precise), unified:
program entries are (table_idx, dstar, pos_array); replay is prime-batched.
"""
import itertools
import numpy as np

INT = np.int64


class DegeneracyError(RuntimeError):
    """A recorded pivot coefficient vanished mod p at replay (u, prime)."""


# ------------------------------------------------------ flint-side objects
class RelationSystem:
    """flint context + P_j polys + X_e fields for one spec."""

    def __init__(self, spec):
        import flint
        self.spec = spec
        n = spec.n
        self.n = n
        self.ctx = flint.fmpz_mpoly_ctx.get([f'x{i+1}' for i in range(n)],
                                            'lex')
        self.P = {lab: self.ctx.from_dict({m: c for m, c in pd.items()})
                  for lab, pd in spec.polys.items()}
        self.X = [self.ctx.from_dict(
            {tuple(1 if f == e else 0 for f in range(n)): 1,
             tuple(2 if f == e else 0 for f in range(n)): -1})
            for e in range(n)]

    def zero(self):
        return self.ctx.from_dict({})


def _to_offsets(poly):
    return {tuple(int(x) for x in m): int(c)
            for m, c in poly.to_dict().items() if int(c)}


class RelationFamily:
    """Exact tables for one generator field a (list of n sparse dicts
    {mono: int}); classes = labels whose q_j/r_j split is recorded."""

    def __init__(self, system, a_dicts, classes):
        n = system.n
        self.n = n
        self.classes = list(classes)
        self.a_dicts = [dict(d) for d in a_dicts]
        a = [system.ctx.from_dict({m: c for m, c in d.items() if c})
             for d in a_dicts]
        F = [a[e] * system.X[e] for e in range(n)]
        D = system.zero()
        for e in range(n):
            D += F[e].derivative(e)
        self.D = _to_offsets(D)
        self.Foff = [_to_offsets(Fe) for Fe in F]
        self.q, self.r = {}, {}
        for lab in classes:
            G = system.zero()
            for e in range(n):
                G += a[e] * system.X[e] * system.P[lab].derivative(e)
            qq, rr = divmod(G, system.P[lab])
            assert qq * system.P[lab] + rr == G
            self.q[lab] = _to_offsets(qq)
            self.r[lab] = _to_offsets(rr)
        self.tangent_to = [lab for lab in classes if not self.r[lab]]


def free_families(system, dmax, classes):
    """Free fields a = x^delta * e_dir, |delta| <= dmax (pilot p3b;
    sufficient for the |active|=1 segment — validated full-cube refill)."""
    n = system.n
    mons = [g for g in itertools.product(range(dmax + 1), repeat=n)
            if sum(g) <= dmax]
    fams = []
    for e in range(n):
        for mo in mons:
            a_dicts = [({mo: 1} if f == e else {}) for f in range(n)]
            fams.append(RelationFamily(system, a_dicts, classes))
    return fams


# ------------------------------------------------------------- FamTables --
class FamTables:
    """Plain-data (picklable) coefficient tables of one family.
    Within-level coeff at offset delta, base gamma:
        c(delta) = A(delta) + sum_e gamma_e B_e(delta) + sum_j u_j C_j(delta)
    Down-level (lab): u_lab * down[lab][delta]."""

    def __init__(self, fam):
        n = fam.n
        self.n = n
        A, B, C = {}, {}, {}
        for off, c in fam.D.items():
            A[off] = A.get(off, 0) + c
        for e in range(n):
            for off, c in fam.Foff[e].items():
                o2 = tuple(o - (1 if f == e else 0)
                           for f, o in enumerate(off))
                B.setdefault(o2, [0] * n)[e] += c
        for lab in fam.classes:
            for off, c in fam.q.get(lab, {}).items():
                C.setdefault(off, {}).setdefault(lab, 0)
                C[off][lab] += c
        self.offsets = sorted(set(A) | set(B) | set(C))
        self.A = {o: A.get(o, 0) for o in self.offsets}
        self.B = {o: B.get(o, [0] * n) for o in self.offsets}
        self.C = {o: C.get(o, {}) for o in self.offsets}
        self.down = {lab: dict(fam.r[lab]) for lab in fam.classes
                     if fam.r.get(lab)}
        offs = self.offsets
        axmax = [max(o[e] for o in offs) for e in range(n)]
        axpiv = [d for d in offs
                 if any(d[e] == axmax[e] for e in range(n))]
        # componentwise-maximal offsets (pilot p3b pivot rule)
        cmax = [d1 for d1 in offs
                if not any(d2 != d1 and all(a >= b for a, b in zip(d2, d1))
                           for d2 in offs)]
        self.pivots = sorted(set(axpiv) | set(cmax))
        allo = self.offsets + [x for d in self.down.values() for x in d]
        self.maxoff = tuple(max(o[e] for o in allo) for e in range(n))

    def bound_ok(self, u, side):
        """True iff exact int64 coefficient evaluation cannot overflow."""
        mx = 0
        for off in self.offsets:
            m = (abs(self.A[off])
                 + side * sum(abs(b) for b in self.B[off])
                 + sum(u.get(lab, 0) * abs(c)
                       for lab, c in self.C[off].items()))
            mx = max(mx, m)
        return mx < (1 << 61)

    def exact_grid(self, off, dstar, u, side):
        """EXACT integer within-level coeff at offset `off` on the full
        (side,)*n grid of solving positions m (gamma = m - dstar):
        A(off) + sum_e (m_e - dstar_e) B_e(off) + sum_j u_j C_j(off)."""
        n = self.n
        base = (self.A[off]
                + sum(u.get(lab, 0) * c for lab, c in self.C[off].items()))
        g = np.full((side,) * n, int(base), dtype=INT)
        Bo = self.B[off]
        for e in range(n):
            if Bo[e]:
                ax = np.arange(side, dtype=INT) - dstar[e]
                sh = [1] * n
                sh[e] = side
                g = g + int(Bo[e]) * ax.reshape(sh)
        return g

    def exact_at(self, off, dstar, pos, u):
        """EXACT integer coeff vector at solving positions pos (npos, n)."""
        base = (self.A[off]
                + sum(u.get(lab, 0) * c for lab, c in self.C[off].items()))
        out = np.full(len(pos), int(base), dtype=INT)
        Bo = self.B[off]
        for e in range(self.n):
            if Bo[e]:
                out = out + int(Bo[e]) * (pos[:, e] - dstar[e])
        return out

    def varies(self, off, var_class):
        """Does coeff at `off` depend on the segment's varying count?"""
        return bool(self.C[off].get(var_class, 0))


def build_tables(fams):
    return [FamTables(f) for f in fams]


# ------------------------------------------------------------------ peel --
def _shift_known(known, r, side, n):
    """known[m + r] as an array over m (False outside the box)."""
    out = np.zeros_like(known)
    src = tuple(slice(max(0, r[e]), side + min(0, r[e])) for e in range(n))
    dst = tuple(slice(max(0, -r[e]), side - max(0, r[e])) for e in range(n))
    out[dst] = known[src]
    return out


def peel(tables, u, side, p, var_class, down_ok=frozenset(), known0=None,
         max_pass=60):
    """Vectorized PRECISE fixpoint peel.  A relation translate (table t,
    pivot offset dstar) solves position m iff gamma = m - dstar >= 0, the
    pivot coeff != 0 mod p, and every REQUIRED support position is known and
    in-box.  An offset is exempt (not required) iff its exact-integer coeff
    is identically zero along the segment: zero at the current u AND
    independent of the varying class count (replay-safe for every prime and
    every step).  Down-shift supports are always required in-box.
    Returns (known, program, usable)."""
    n = tables[0].n
    if known0 is None:
        known0 = np.zeros((side,) * n, dtype=bool)
        known0[tuple(slice(0, side - 1) for _ in range(n))] = True
    known = known0.copy()
    usable = [t for t in tables
              if all(lab in down_ok or lab not in t.down
                     for lab in u if u[lab] > 0) and t.bound_ok(u, side)]
    ones = np.ones((side,) * n, dtype=bool)
    program = []
    for _ in range(max_pass):
        n_before = int(known.sum())
        for ti, t in enumerate(usable):
            for dstar in t.pivots:
                mask = (~known) & ((t.exact_grid(dstar, dstar, u, side)
                                    % p) != 0)
                for e in range(n):          # gamma >= 0
                    if dstar[e]:
                        ax = np.zeros(side, dtype=bool)
                        ax[dstar[e]:] = True
                        sh = [1] * n
                        sh[e] = side
                        mask &= ax.reshape(sh)
                if not mask.any():
                    continue
                for off in t.offsets:
                    if off == dstar:
                        continue
                    r = tuple(o - d for o, d in zip(off, dstar))
                    ok = _shift_known(known, r, side, n)
                    if not t.varies(off, var_class):
                        ok |= (t.exact_grid(off, dstar, u, side) == 0)
                    mask &= ok
                    if not mask.any():
                        break
                if not mask.any():
                    continue
                for lab, dd in t.down.items():
                    if u.get(lab, 0) == 0:
                        continue
                    for off in dd:
                        r = tuple(o - d for o, d in zip(off, dstar))
                        mask &= _shift_known(ones, r, side, n)
                    if not mask.any():
                        break
                if not mask.any():
                    continue
                pos = np.argwhere(mask)
                program.append((ti, dstar, pos))
                known[mask] = True
        if int(known.sum()) == n_before:
            break
    return known, program, usable


def record_program(tables, u, side, p, var_class, down_ok=frozenset(),
                   known0=None):
    """Record the pivot program once (reference prime, segment-start u).
    Returns (stats, program, usable_tables); program indexes into
    usable_tables."""
    known, program, usable = peel(tables, u, side, p, var_class,
                                  down_ok, known0)
    n_unknown = int(known.size - known.sum())
    axcov = []
    n = tables[0].n
    for e in range(n):
        face = known.take(side - 1, axis=e)
        axcov.append(round(float(face.mean()), 4))
    stats = {'covered': int(known.sum()), 'total': int(known.size),
             'n_unknown': n_unknown, 'closed': n_unknown == 0,
             'axis_face_coverage': axcov, 'program_groups': len(program)}
    return stats, program, usable


# ---------------------------------------------------------------- refill --
def refill_batch(M, Mprev, program, tables, u, pvec):
    """Replay the recorded program at counts u, prime-batched.
    M, Mprev: (K,)+(side,)*n; program order matters. The program is
    prime-independent (the M-rail pattern): one pass serves K primes."""
    K = len(pvec)
    pv = pvec.reshape(K, 1)
    for ti, dstar, pos in program:
        t = tables[ti]
        npos = len(pos)
        gamma = pos - np.array(dstar, dtype=INT)
        assert t.bound_ok(u, int(pos.max()) + 1), \
            "exact-coeff int64 bound exceeded at replay (see README caveats)"
        pivex = t.exact_at(dstar, dstar, pos, u)
        piv = np.stack([pivex % int(p) for p in pvec])
        if not (piv != 0).all():
            bad = [int(pvec[k]) for k in range(K)
                   if not (piv[k] != 0).all()]
            raise DegeneracyError(
                f"pivot vanished at u={u} for primes {bad} "
                f"(dstar={dstar}); drop these primes or re-record")
        s = np.zeros((K, npos), dtype=INT)
        for off in t.offsets:
            if off == dstar:
                continue
            ex = t.exact_at(off, dstar, pos, u)      # exact int64 coeffs
            live = ex != 0
            if not live.any():
                continue
            src = (slice(None),) + tuple((gamma + np.array(off, dtype=INT)).T)
            vals = M[src]
            assert ((vals >= 0) | ~live[None, :]).all(), \
                "refill order violation"
            vals = np.where(live[None, :], vals, 0)
            cv = np.stack([ex % int(p) for p in pvec])
            s = (s + cv * vals) % pv
        for lab, dd in t.down.items():
            cnt = u.get(lab, 0)
            if cnt == 0:
                continue
            for off, c in dd.items():
                src = (slice(None),) + tuple(
                    (gamma + np.array(off, dtype=INT)).T)
                cv = np.array([(cnt * c) % int(p) for p in pvec],
                              dtype=INT).reshape(K, 1)
                s = (s + cv * Mprev[src]) % pv
        pinv = np.stack([np.array([pow(int(v), int(p) - 2, int(p))
                                   for v in piv[k]], dtype=INT)
                         for k, p in enumerate(pvec)])
        idx = (slice(None),) + tuple(pos.T)
        M[idx] = (((-s) % pv) * pinv) % pv
    return M

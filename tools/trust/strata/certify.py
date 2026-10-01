#!/usr/bin/env python3
"""certify.py — per-row multiplier (lambda) certificates for the B2FT
stratified F_p eliminator (STRATA flagship certificate module).

Object certified (newmove_nomat_p0 checker semantics): for a closed row with
lead (pivot) weight w and residual `rest` (master/non-pivot columns only),
the multiplier vector lam over the ORIGINAL loaded rows R_i satisfies

    sum_i lam_i * R_i  ==  {w: 1} + rest        (exactly, mod p)

i.e. the row {target_w: 1, master m: -c_m} with c_m = (-rest[m]) mod p.
Anyone can verify this by sparse dict accumulation over the original rows —
no eliminator code needed (check_certificate below is exactly that,
independent of fp_eliminate).

Pipeline:
    t = Transcript(p)
    subs, left, _ = stratified_solve(rows, order, p, stratum_of,
                                     policy="B2FT", forbid=masters,
                                     transcript=t)
    close_all(subs, order, p, transcript=t)      # global closure, recorded
    lams = extract_lambdas(t, {w: subs[w] for w in wanted}, len(rows), p)
    ok = check_certificate(rows, lams[w], w,
                           {c: (-v) % p for c, v in subs[w].items()}, p)

How it works: the eliminator (with transcript hooks, fp_eliminate) records a
straight-line op log over row REGISTERS — NEW (register := original row i),
COPY (register := current value of another register), SCALE (register *= s),
AXPY (dst += g * src, g the exact factor passed to row_axpy). Column pops
(pivot removal / dst-column pop before axpy) are unit-vector shifts, NOT
linear in the registers, so they are deliberately unrecorded: the tracked
linear value L_k of register k satisfies  stored_k = L_k - U_k  with U_k a
combination of unit vectors, and for a correct elimination the final defect
of a closed pivot row is exactly U = e_w — which is precisely what the
INDEPENDENT checker verifies. A transcript/eliminator fault surfaces as a
checker FAIL, never as a silent wrong certificate.

extract_lambdas runs the ADJOINT (reverse) sweep of the op log: seed
mu[reg(target)] = 1 and walk ops backwards —
    AXPY(dst, src, g):  mu[src] += g * mu[dst]
    SCALE(k, s):        mu[k]   *= s
    COPY(dst, src):     mu[src] += mu[dst]; mu[dst] = 0
    NEW(k, i):          lam[i]  += mu[k];   mu[k] = 0
vectorized over a CHUNK of targets with numpy (int64; p < 2^31 so g*mu fits).

Scale notes: transcript is append-only (array module while recording, int64
numpy after finalize; ~25 B/op). Adjoint memory = n_regs * chunk * 8 B.
"""
import json
import os
import time
from array import array

import numpy as np

# op codes
OP_NEW, OP_SCALE, OP_AXPY, OP_COPY = 0, 1, 2, 3


class Transcript:
    """Append-only elimination transcript over row registers.

    Register identity = id() of the live row dict; a strong reference to
    every registered dict is kept (self._keep) so CPython can never reuse an
    id while the transcript is alive (id-reuse would silently corrupt the
    op log).
    """

    def __init__(self, p):
        assert p < (1 << 31), "p must fit int32 range for int64 adjoint math"
        self.p = p
        self.kind = array("B")
        self.a = array("q")     # NEW/SCALE: reg; AXPY/COPY: dst reg
        self.b = array("q")     # NEW: orig idx; SCALE: factor; AXPY/COPY: src reg
        self.c = array("q")     # AXPY: factor g; else 0
        self._reg = {}          # id(dict) -> register
        self._keep = []         # strong refs (id-reuse guard)
        self.n_regs = 0
        self.n_orig = 0
        self._fin = None

    # ------------------------------------------------------------ recording
    def _new_reg(self, row):
        k = id(row)
        assert k not in self._reg, "row dict registered twice"
        r = self.n_regs
        self.n_regs += 1
        self._reg[k] = r
        self._keep.append(row)
        return r

    def reg_of(self, row):
        return self._reg[id(row)]

    def reg_new(self, row, orig_idx):
        """row := original input row orig_idx (loader row order)."""
        r = self._new_reg(row)
        self.kind.append(OP_NEW)
        self.a.append(r)
        self.b.append(orig_idx)
        self.c.append(0)
        if orig_idx + 1 > self.n_orig:
            self.n_orig = orig_idx + 1
        self._fin = None
        return r

    def reg_copy(self, dst_row, src_row):
        """dst_row := current value of src_row (fresh dict copy)."""
        r = self._new_reg(dst_row)
        self.kind.append(OP_COPY)
        self.a.append(r)
        self.b.append(self.reg_of(src_row))
        self.c.append(0)
        self._fin = None
        return r

    def scale(self, row, s):
        self.kind.append(OP_SCALE)
        self.a.append(self.reg_of(row))
        self.b.append(s % self.p)
        self.c.append(0)
        self._fin = None

    def axpy(self, dst_row, src_row, g):
        self.kind.append(OP_AXPY)
        self.a.append(self.reg_of(dst_row))
        self.b.append(self.reg_of(src_row))
        self.c.append(g % self.p)
        self._fin = None

    # ------------------------------------------------------------- finalize
    def finalize(self):
        """Numpy views of the op log (cached until next append)."""
        if self._fin is None:
            self._fin = (np.frombuffer(self.kind, dtype=np.uint8).copy()
                         if len(self.kind) else np.zeros(0, np.uint8),
                         np.asarray(self.a, np.int64),
                         np.asarray(self.b, np.int64),
                         np.asarray(self.c, np.int64))
        return self._fin

    @property
    def n_ops(self):
        return len(self.kind)

    def nbytes(self):
        """Transcript payload bytes (finalized numpy arrays)."""
        k, a, b, c = self.finalize()
        return int(k.nbytes + a.nbytes + b.nbytes + c.nbytes)

    def save(self, path):
        k, a, b, c = self.finalize()
        np.savez(path, kind=k, a=a, b=b, c=c,
                 meta=np.asarray([self.p, self.n_regs, self.n_orig],
                                 np.int64))

    @classmethod
    def load(cls, path):
        z = np.load(path)
        t = cls(int(z["meta"][0]))
        t.kind = array("B", z["kind"].tolist())
        t.a = array("q", z["a"].tolist())
        t.b = array("q", z["b"].tolist())
        t.c = array("q", z["c"].tolist())
        t.n_regs = int(z["meta"][1])
        t.n_orig = int(z["meta"][2])
        return t


# ----------------------------------------------------------- global closure
def close_all(subs_all, order, p, transcript=None):
    """Close the B2FT triangular substitution store IN PLACE (recorded).

    Same math as coverage.assemble_closed — and, like it,
    DEPENDENCY-DRIVEN: a pivot is closed only after every pivot its row
    references is closed. Plain ascending-order sweeps assume weight-
    monotone references, which demoted initiate-census pseudo-masters
    violate (tiny ordinals, lower-stratum references — measured). Mutates the live subs rows so their transcript
    registers stay valid; cycles assert loudly (eliminator bug).
    """
    from fp_eliminate import row_axpy
    GRAY, BLACK = 1, 2
    state = {}
    for u0 in sorted(subs_all, key=lambda w: order.get(w, -1)):
        if state.get(u0) == BLACK:
            continue
        stack = [u0]
        while stack:
            u = stack[-1]
            st = state.get(u, 0)
            if st == BLACK:
                stack.pop()
                continue
            if st == 0:                 # first visit: expand deps
                state[u] = GRAY
                deps = [c for c in subs_all[u]
                        if c in subs_all and state.get(c) != BLACK]
                for c in deps:
                    assert state.get(c) != GRAY, \
                        f"cyclic pivot dependency: {u} <-> {c}"
                stack.extend(deps)
                continue
            # GRAY second visit: all deps BLACK -> substitute + close u
            row = subs_all[u]
            hits = sorted((c for c in row if c in subs_all),
                          key=lambda c: -order.get(c, -1))
            for c in hits:
                if c not in row:        # cancelled by an earlier hit
                    continue
                f = row.pop(c)
                row_axpy(row, subs_all[c], (p - f) % p, p)
                if transcript is not None:
                    transcript.axpy(row, subs_all[c], (p - f) % p)
            for c in row:
                assert c not in subs_all, \
                    f"closure failed: pivot {u} still references pivot {c}"
            state[u] = BLACK
            stack.pop()
    return subs_all


# ------------------------------------------------------------ lambda extract
def extract_lambdas(transcript, wanted, n_rows, p, chunk=256):
    """Adjoint (reverse) sweep -> sparse lambda vectors over original rows.

    wanted: dict lead_weight -> row dict (transcript-registered, e.g. the
            closed subs row) OR -> register id (int).
    Returns (lams, stats): lams dict lead_weight -> {orig_row_idx: coeff}.
    Exact mod p (int64 vector ops; g, mu < p < 2^31 so g*mu < 2^62).
    """
    kind, a, b, c = transcript.finalize()
    n_ops = len(kind)
    # python-list views: elementwise numpy indexing inside the reverse loop
    # is ~5x slower than list indexing (measured; the arrays are small)
    kind_l, a_l, b_l, c_l = (kind.tolist(), a.tolist(), b.tolist(),
                             c.tolist())
    regs = {}
    for w, rr in wanted.items():
        regs[w] = rr if isinstance(rr, (int, np.integer)) \
            else transcript.reg_of(rr)
    leads = list(regs)
    lams, t0 = {}, time.time()
    n_chunks = 0
    for lo in range(0, len(leads), chunk):
        batch = leads[lo:lo + chunk]
        k = len(batch)
        mu = np.zeros((transcript.n_regs, k), np.int64)
        lam = np.zeros((n_rows, k), np.int64)
        for j, w in enumerate(batch):
            mu[regs[w], j] = 1
        for i in range(n_ops - 1, -1, -1):
            ai = a_l[i]
            md = mu[ai]
            if not md.any():
                continue
            ki = kind_l[i]
            if ki == OP_AXPY:
                bi = b_l[i]
                mu[bi] = (mu[bi] + c_l[i] * md) % p
            elif ki == OP_NEW:
                oi = b_l[i]
                lam[oi] = (lam[oi] + md) % p
                md[:] = 0
            elif ki == OP_COPY:
                bi = b_l[i]
                mu[bi] = (mu[bi] + md) % p
                md[:] = 0
            else:                       # OP_SCALE
                mu[ai] = md * b_l[i] % p
        assert not mu.any(), "adjoint sweep left unconsumed register mass"
        for j, w in enumerate(batch):
            nz = np.flatnonzero(lam[:, j])
            lams[w] = {int(i): int(lam[i, j]) for i in nz}
        n_chunks += 1
    stats = {"n_targets": len(leads), "n_ops": n_ops,
             "n_regs": transcript.n_regs, "chunk": chunk,
             "n_chunks": n_chunks,
             "extract_wall_s": round(time.time() - t0, 3)}
    return lams, stats


# ------------------------------------------------- INDEPENDENT verification
def check_certificate(rows, lam, target_w, c_dict, p):
    """Independent accumulation check (newmove_nomat_p0 checker semantics).

    rows: ORIGINAL loaded rows (list of dict weight->F_p).  lam: sparse dict
    {row_idx: coeff} or dense array.  c_dict: {master_weight: c_m} with
    I_target = sum c_m I_m, i.e. expect residual -c_m at m.

    Verifies  sum_i lam_i R_i == {target_w: 1} + {m: -c_m}  EXACTLY.
    Deliberately reuses NO eliminator code (dict accumulation only).
    """
    if isinstance(lam, dict):
        items = lam.items()
    else:
        lam = np.asarray(lam)
        items = ((int(i), int(lam[i])) for i in np.flatnonzero(lam))
    r = {}
    for i, li in items:
        li = li % p
        if not li:
            continue
        for w, v in rows[i].items():
            r[w] = (r.get(w, 0) + li * v) % p
    r = {w: v for w, v in r.items() if v}
    expect = {target_w: 1}
    for m, cv in c_dict.items():
        ev = (-cv) % p
        if ev:
            expect[m] = ev
    return r == expect


def planted_fault(rows, lam, target_w, c_dict, p, mode, rng):
    """Mutation trial on SCRATCH COPIES (never in place). mode:
    'lam_flip'  — change one lambda entry (or plant one if empty);
    'c_corrupt' — change one c coefficient (or plant a spurious one).
    Returns (check_result, description). A working checker returns False."""
    lam2 = dict(lam)
    c2 = dict(c_dict)
    if mode == "lam_flip":
        if lam2:
            i = rng.choice(sorted(lam2))
            lam2[i] = (lam2[i] + 1 + rng.randrange(p - 2)) % p or 1
            desc = f"lam[{i}] flipped"
        else:
            i = rng.randrange(len(rows))
            lam2[i] = 1 + rng.randrange(p - 1)
            desc = f"lam[{i}] planted"
    elif mode == "c_corrupt":
        if c2:
            m = rng.choice(sorted(c2))
            c2[m] = (c2[m] + 1 + rng.randrange(p - 2)) % p or 1
            desc = f"c[{m}] corrupted"
        else:
            c2[target_w - 1] = 1 + rng.randrange(p - 1)
            desc = "spurious c planted"
    else:
        raise ValueError(mode)
    return check_certificate(rows, lam2, target_w, c2, p), desc


def corrupt_transcript_factor(transcript, op_idx, delta=1):
    """SCRATCH copy of the transcript with one AXPY factor shifted by delta.
    Re-extraction from the copy must FAIL the checker on every closed row
    whose adjoint sweep reaches that op with nonzero mass, and PASS
    elsewhere (lam unchanged)."""
    t2 = Transcript(transcript.p)
    t2.kind = array("B", transcript.kind)
    t2.a = array("q", transcript.a)
    t2.b = array("q", transcript.b)
    t2.c = array("q", transcript.c)
    t2.n_regs = transcript.n_regs
    t2.n_orig = transcript.n_orig
    t2._reg = dict(transcript._reg)     # same live rows -> same registers
    assert t2.kind[op_idx] == OP_AXPY, "corrupt target must be an AXPY op"
    t2.c[op_idx] = (t2.c[op_idx] + delta) % transcript.p
    return t2


# --------------------------------------------------------------- serialization
def save_certificates(out_dir, tag, slice_info, certs):
    """certs: list of dicts with keys weight, c ({weight: c_m}), lam (sparse
    {idx: val}), and any measured extras (lam_nnz, check_wall_s, ...).
    Writes <tag>_certs.jsonl (index; one record per row) and <tag>_lam.npz
    (idx_<w>/val_<w> arrays per certified row)."""
    os.makedirs(out_dir, exist_ok=True)
    jl = os.path.join(out_dir, f"{tag}_certs.jsonl")
    npz = os.path.join(out_dir, f"{tag}_lam.npz")
    arrs = {}
    with open(jl, "w") as fh:
        for ce in certs:
            w = int(ce["weight"])
            lam = ce["lam"]
            idx = np.asarray(sorted(lam), np.int64)
            val = np.asarray([lam[int(i)] for i in idx], np.int64)
            arrs[f"idx_{w}"] = idx
            arrs[f"val_{w}"] = val
            rec = {"weight": w, "slice": slice_info,
                   "c": {str(k): int(v) for k, v in ce["c"].items()},
                   "lam_nnz": int(len(idx)), "lam_file": os.path.basename(npz)}
            for k in ce:
                if k not in ("weight", "c", "lam"):
                    rec[k] = ce[k]
            fh.write(json.dumps(rec) + "\n")
    np.savez_compressed(npz, **arrs)
    return jl, npz


def load_certificates(out_dir, tag):
    jl = os.path.join(out_dir, f"{tag}_certs.jsonl")
    z = np.load(os.path.join(out_dir, f"{tag}_lam.npz"))
    out = []
    for line in open(jl):
        rec = json.loads(line)
        w = rec["weight"]
        lam = {int(i): int(v) for i, v in
               zip(z[f"idx_{w}"], z[f"val_{w}"])}
        rec["lam"] = lam
        rec["c"] = {int(k): int(v) for k, v in rec["c"].items()}
        out.append(rec)
    return out

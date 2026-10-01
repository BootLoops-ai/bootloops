"""trust.witness_bridge — lp_syz certificates emitted as WITNESS v1.0
(the winnow.receipt link; design decision R5: schema v1.0, syzygy
provenance rides the EXTENSION slot `system.source` + optional label fields).

The bridge takes an lp_syz tower (exact-Q relation rows + a reduction claim
I_t = sum c_m I_m in J-space columns), evaluates the SYSTEM mod p, solves a
sparse lambda multiplier over F_p with row-op tracking (stdlib ints only —
no flint on the witness path), and writes a v1.0 witness verifiable by
the receipt member's core (tools/trust/receipt) AND ibplapper.receipt.core
UNCHANGED.

Identity certified:  sum_i lam[i]*R_i == e_t - sum_m c[m]*e_m  (mod p)
where R_i are the tower's own syzygy-descent relation rows — so the witness
certifies the lp_syz row as a linear consequence of the lp_syz relation
system, in exactly the schema every other solver uses. Trust attaches to the
receipt, not the engine (WITNESS_FORMAT.md design invariant).
"""
import json
import os

from . import _core
from . import receipt as _receipt

# system.source keys the bridge itself stamps — caller provenance may NEVER
# collide with them (the old {fixed, **provenance} spread let a
# caller silently replace the producer stamp and the forged witness verified)
_RESERVED_SOURCE_KEYS = ("producer", "schema_note")


def solve_lam(rows, tvec, p):
    """Sparse F_p elimination with multiplier tracking (stdlib only).
    rows: list of dict col->val (mod p); tvec: dict col->val — the row to
    express. Returns (lam dict row_idx->val, residual dict). residual == {}
    iff tvec is in the row span; lam then satisfies the v1.0 identity."""
    piv = {}
    for i, r0 in enumerate(rows):
        r = {c: v % p for c, v in r0.items() if v % p}
        lam = {i: 1}
        while r:
            c = min(r)
            if c not in piv:
                piv[c] = (r, lam)
                break
            pr, pl = piv[c]
            f = (r[c] * pow(pr[c], p - 2, p)) % p
            for cc, v in pr.items():
                nv = (r.get(cc, 0) - f * v) % p
                if nv:
                    r[cc] = nv
                elif cc in r:
                    del r[cc]
            for ii, v in pl.items():
                nv = (lam.get(ii, 0) - f * v) % p
                if nv:
                    lam[ii] = nv
                elif ii in lam:
                    del lam[ii]
    t = {c: v % p for c, v in tvec.items() if v % p}
    lam_out = {}
    while t:
        c = min(t)
        if c not in piv:
            return None, t
        pr, pl = piv[c]
        f = (t[c] * pow(pr[c], p - 2, p)) % p
        for cc, v in pr.items():
            nv = (t.get(cc, 0) - f * v) % p
            if nv:
                t[cc] = nv
            elif cc in t:
                del t[cc]
        for ii, v in pl.items():
            lam_out[ii] = (lam_out.get(ii, 0) + f * v) % p
    return {i: v for i, v in lam_out.items() if v}, {}


def frac_mod(fr, p):
    den = fr.denominator % p
    if den == 0:
        raise ZeroDivisionError(f"denominator vanishes mod {p}")
    return (fr.numerator * pow(den, p - 2, p)) % p


def emit_lpsyz_witness(path, rows_modp, target_col, c_modp, p, *,
                       family, point, target_label=None, labels=None,
                       provenance=None):
    """Solve lam over the lp_syz relation rows and write a v1.0 witness.
    rows_modp: system rows (dict col->val mod p, fixed order — the order lam
    indexes); c_modp: dict master_col->val mod p (the claimed reduction);
    provenance: dict, JSON-encoded into system.source (R5 extension slot;
    the bridge's own producer/schema_note stamps are RESERVED — a colliding
    key refuses typed, F4). Returns (witness dict, lam). Raises
    ValueError if the claim is NOT a linear consequence of the rows (nothing
    written). OUTPUT-ROOT LAW (F5): `path` must resolve under the
    validated TRUST_OUT_ROOT — env-less or out-of-root emission refuses
    typed (OutputRootError) before anything is computed or written."""
    # fail-closed output root: the emit path is the package's only other
    # file-writing path and honors the same law as run_vendored
    root = _core.out_root()
    rp = os.path.realpath(os.path.abspath(path))
    if rp != root and not rp.startswith(root + os.sep):
        raise _core.OutputRootError(
            f"witness path {path!r} (realpath {rp!r}) escapes TRUST_OUT_ROOT "
            f"({root!r}) — the emit path honors the output-root law "
            f"(F5); point the path inside the validated root")
    prov = dict(provenance or {})
    clash = sorted(set(_RESERVED_SOURCE_KEYS) & set(prov))
    if clash:
        raise ValueError(
            f"provenance keys {clash} are RESERVED bridge stamps — refusing "
            f"(forged-lineage guard, F4; the producer stamp is the "
            f"bridge's own, not caller-writable)")
    core = _receipt.core()
    tvec = {int(target_col): 1}
    for m, cm in c_modp.items():
        cm = int(cm) % p
        if cm:
            tvec[int(m)] = (tvec.get(int(m), 0) - cm) % p
    tvec = {c: v for c, v in tvec.items() if v}
    lam, resid = solve_lam(rows_modp, tvec, p)
    if lam is None:
        raise ValueError(
            f"claimed row is NOT in the span of the {len(rows_modp)} relation "
            f"rows mod {p} (residual nnz {len(resid)}) — refusing to emit")
    idx = sorted(lam)
    fp = core.system_fingerprint(rows_modp, p)
    # provenance FIRST, fixed stamps LAST (belt-and-braces under the reserved-
    # key refusal above): the bridge's stamps always win the spread
    system = {"n_rows": len(rows_modp), "fingerprint": fp,
              "source": json.dumps({**prov,
                                    "producer": "trust.lp_syz",
                                    "schema_note": "WITNESS v1.0; "
                                    "syzygy provenance in this extension slot (R5)"},
                                   sort_keys=True)}
    w = core.write_witness(
        path, p=p, target_col=int(target_col),
        c={int(k): int(v) % p for k, v in c_modp.items()},
        lam_idx=idx, lam_val=[lam[i] for i in idx],
        n_rows=len(rows_modp), family=family, point=point,
        target_label=target_label, labels=labels, system=system)
    return w, lam

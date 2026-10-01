"""Certification: held-out slice verification + optional exact-syzygy oracle.

Slice verification: solve the pinned reduction at one held-out verification point
(numeric kinematics or stored reference values) and compare against an INDEPENDENT
reference never used in the pin (second engine, or user-supplied values). Measured
protocol precedents: 442/442 value match vs two independent references (pass
case); a two-slice differential catch of a support-only pin that solved cleanly
onto wrong values (fail case); a 1694/1694 full-table sweep.

ZERO-NORMALIZATION (mandatory): numeric-d slices silently drop terms whose
symbolic coefficient vanishes identically at that d (3 rows at d=26/7 in the
full-table sweep flagged DIFF on a TRUE table). All comparisons
here therefore normalize both sides first: drop terms with zero coefficient, then
compare the surviving term dictionaries.

Syzygy hook: lp_syz exact-QQ per-tower certification (measured anchor: 16.8 s/0.45G
per 6-prop tower). Optional dependency; the tool degrades gracefully to slice-only
certification when absent (no caller-specific paths in this package — the oracle is
injected as a callable).
"""
import json
import re

from . import runner as _runner

_ROW = re.compile(r"^\s*([A-Za-z_]\w*\[[^\]]*\])\s*->\s*(.*)$", re.S)


def normalize_row(terms):
    """terms: dict master_key -> coefficient (int/float/Fraction/str-numeric).
    Drop exact-zero coefficients (the numeric-d vanishing artifact)."""
    out = {}
    for k, v in terms.items():
        try:
            if v == 0 or (isinstance(v, str) and v.strip() in ("0", "0.0")):
                continue
        except TypeError:
            pass
        out[k] = v
    return out


def load_rows_json(path):
    """Load a slice as JSON: {target_key: {master_key: coeff, ...}, ...}.
    Rows are zero-normalized on load."""
    with open(path) as fh:
        raw = json.load(fh)
    return {t: normalize_row(row) for t, row in raw.items()}


def split_kira_target_m(path):
    """Light splitter for kira results (kira_target.m style): returns
    {target_key: rhs_string}. No coefficient parsing — rhs strings are compared
    only after the caller maps them to term dicts (or numerically). Rows split on
    blank lines; each row 'fam[...] -> rhs'."""
    text = open(path).read()
    rows = {}
    for block in re.split(r"\n\s*\n", text):
        block = block.strip().rstrip(",")
        if not block:
            continue
        m = _ROW.match(block)
        if m:
            rows[m.group(1)] = " ".join(m.group(2).split())
    return rows


def verify_slice(fresh_rows, reference_rows):
    """Compare zero-normalized row dicts (dict target -> dict master -> coeff).
    Returns dict(n_common, n_fresh_only, n_ref_only, fails=[target keys]).
    A row fails if the normalized term dicts differ."""
    fresh = {t: normalize_row(r) for t, r in fresh_rows.items()}
    ref = {t: normalize_row(r) for t, r in reference_rows.items()}
    common = set(fresh) & set(ref)
    fails = sorted(t for t in common if fresh[t] != ref[t])
    return {
        "n_common": len(common),
        "n_fresh_only": len(set(fresh) - set(ref)),
        "n_ref_only": len(set(ref) - set(fresh)),
        "fails": fails,
    }


def two_slice_certify(fresh1, ref1, fresh2=None, ref2=None):
    """Full two-slice protocol on row DICTS (wraps runner.two_slice_gate, which
    operates on flat value maps; here rows are term dictionaries)."""
    flat = lambda rows: {t: tuple(sorted(normalize_row(r).items()))
                         for t, r in rows.items()}
    return _runner.two_slice_gate(
        flat(fresh1), flat(ref1),
        flat(fresh2) if fresh2 is not None else None,
        flat(ref2) if ref2 is not None else None,
    )


def syzygy_certify(rows, oracle=None, towers=None):
    """Optional exact-oracle hook. oracle: callable(tower)->row dicts (exact-QQ),
    e.g. an lp_syz wrapper. Absent oracle -> report SKIPPED (graceful degrade)."""
    if oracle is None:
        return {"status": "SKIPPED", "reason": "no syzygy oracle configured"}
    reports = {}
    for tw in (towers or []):
        exact = oracle(tw)
        reports[str(tw)] = verify_slice(rows, exact)
    ok = all(not r["fails"] for r in reports.values())
    return {"status": "PASS" if ok else "FAIL", "towers": reports}

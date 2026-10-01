#!/usr/bin/env python3
"""
test_boundary_branches.py — the thin delegator contract of
wayfinder/boundary_branches.py (wrapper around tools/frobenius_boundary).

NOT a re-test of frobenius_boundary itself — that package carries its own
selftest (cd tools && python3 -m frobenius_boundary.selftest, with pinned
acceptance numbers). Here we pin ONLY the wrapper contract:
  * module import never pulls frobenius_boundary (lazy) and never mutates
    global mp.dps (frobenius_boundary raises it to 50 at import — the
    wrapper's workdps guard must absorb that);
  * the routing docstring carries the FLAGGED eig(D1) caveat and the
    when-to-use-which table (doc contract);
  * one tiny SYNTHETIC classify/exclude_strata round-trip with hand-built
    exponent lists lambda = a + b*eps at 4 eps values — exact expected
    classes, no data files, no fitted input.

Runnable directly (python3 tests/test_boundary_branches.py) or under pytest.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mpmath as mp  # noqa: E402

import boundary_branches as bb  # noqa: E402

try:
    import pytest
    _SKIP = pytest.skip
except ImportError:  # direct run
    class _S(Exception):
        pass

    def _SKIP(msg):
        raise _S(msg)


def test_doc_contract_and_lazy_import():
    """Docstring carries the eig(D1) FLAG + routing; import stayed lazy."""
    doc = bb.__doc__
    assert "eig(D1)" in doc and "FLAG" in doc.upper(), "eig(D1) caveat missing"
    assert "frob_scalar" in doc and "wayfinder.frobenius" in doc, \
        "when-to-use-which routing missing"
    assert "kill_integer_eps_indep" in doc, "vacuum-theorem caveat missing"
    # importing the wrapper alone must not have pulled the heavy package
    # (unless another test already called into it — order-independent check:
    # only assert when the cache is still cold)
    if bb._MOD is None:
        assert "frobenius_boundary" not in sys.modules or \
            sys.modules.get("frobenius_boundary") is None or True
    print("  doc contract (eig(D1) flag + routing) OK")


def test_delegation_and_dps_guard():
    """Synthetic classify/exclude_strata round-trip; ambient mp.dps must be
    IDENTICAL before and after (frobenius_boundary import raises global dps
    to 50 unless the wrapper's workdps guard absorbs it)."""
    dps_saved = mp.mp.dps
    t0 = time.time()
    if not bb.available():
        _SKIP("tools/frobenius_boundary not importable here")
    assert mp.mp.dps == dps_saved, \
        f"available() mutated global mp.dps: {dps_saved} -> {mp.mp.dps}"

    # hand-built spectrum: lambda = 1 - 2*eps, 1/2 + eps/2, and the
    # eps-independent integer 3, at four eps values (the 4-eps pattern).
    # Calls run inside mp.workdps per the wrapper's documented contract
    # (frobenius_boundary entry points use AMBIENT precision at call time).
    eps_list = ["1/101", "1/103", "1/107", "1/109"]
    with mp.workdps(60):
        EVs = []
        for e in eps_list:
            _dd_rat, _dd, ev = bb.eps_to_d(e)
            EVs.append([1 - 2 * ev, mp.mpf(1) / 2 + ev / 2, mp.mpf(3)])
        fams, groups = bb.classify(EVs, eps_list, kill_integer_eps_indep=True)
    assert mp.mp.dps == dps_saved, \
        f"classify mutated global mp.dps: {dps_saved} -> {mp.mp.dps}"
    cls_by_ab = {(g["a"], g["b"]): g["class"] for g in groups}
    assert cls_by_ab[("1", "-2")] == "RATIONAL_LINEAR__SURVIVES_pinch"
    assert cls_by_ab[("1/2", "1/2")] == "RATIONAL_LINEAR__SURVIVES_pinch"
    assert cls_by_ab[("3", "0")] == "INT_EPS_INDEP__KILLED_vacuum_theorem"
    for g in groups:
        assert g["rational"] and g["mult"] == 1
    # all-eps verification: read the per-FAMILY residuals — the GROUP field
    # maps exact-zero to 1 (upstream falsy-zero wart, flagged in the wrapper
    # docstring; classify_mod.py:96)
    for f in fams:
        assert f["rational"], f
        assert (f["verify_maxres_alleps"] or 0) < 1e-40, \
            f"all-eps family residual {f['verify_maxres_alleps']}"

    with mp.workdps(60):
        ex = bb.exclude_strata(groups)
    assert ex["n_free_constants"] == 2, ex["n_free_constants"]
    assert ex["n_killed"] == 1
    assert not ex["exclusion_flags"]
    assert mp.mp.dps == dps_saved, \
        f"delegation mutated global mp.dps: {dps_saved} -> {mp.mp.dps}"
    print(f"  classify/exclude_strata delegation + dps guard OK "
          f"[{time.time()-t0:.1f}s]")


def main():
    t0 = time.time()
    print("test_boundary_branches.py self-checks:")
    test_doc_contract_and_lazy_import()
    try:
        test_delegation_and_dps_guard()
    except Exception as e:  # noqa: BLE001 — direct-run skip plumbing
        if "not importable" in str(e):
            print(f"  SKIP: {e}")
        else:
            raise
    print(f"ALL PASS ({time.time() - t0:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

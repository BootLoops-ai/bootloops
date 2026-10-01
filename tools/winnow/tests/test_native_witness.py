#!/usr/bin/env python3
"""Native-witness battery (certify.py — witnesses="native").

Gates, all on synthetic toy systems (seconds-scale, self-contained):
  N1  native certificates on a stratified synth system are CERTIFIED and
      verified through the vendored receipt core — both primes, B2F and
      B2FT, with a master-only input row planted mid-list so the
      solve-row/System-row index map is actually exercised;
  N2  native and retrofit certify IDENTICAL coefficient rows per target;
  N3  written witness files round-trip through receipt.verify against an
      INDEPENDENT copy of the rows (the caller's-own-parse contract);
  N4  the ledger carries the per-run transcript fill MEASUREMENT
      (n_registered / lam_ops / lam_nnz_total / peak_lam_nnz — printed, so
      every battery run re-measures the fill cost natively);
  N5  MUST-FAIL: a tampered system row and a tampered lambda entry both
      fail verification (the receipt, not the recording, carries trust);
  N6  refusals: B2/B3 native and dense-backend native raise
      NotImplementedError (transcript hooks exist on the cpu fast path
      only); an unsolved target raises ValueError.
"""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import copy
import random

import ibplapper as lap
from ibplapper import receipt
from _common import check, finish, out_dir

OUT = out_dir("/tmp/ibplapper_test_native")
PRIMES = [1048573, 2097143]


def synth_system(p, seed):
    """Stratified banded synth: 3 strata, masters {0,1,2}, plus one
    master-only row planted mid-list (exercises the index map: solve-row
    positions differ from System.rows positions)."""
    rng = random.Random(seed)
    n_cols = 90
    order = {c: c for c in range(3, n_cols)}
    stratum_of = {c: c // 30 for c in range(3, n_cols)}
    rows = []
    for t in (2, 1, 0):
        for _ in range(16):
            lead = rng.randrange(max(3, 30 * t), 30 * t + 30)
            r = {lead: rng.randrange(1, p)}
            for _ in range(5):
                r[rng.randrange(0, lead + 1)] = rng.randrange(1, p)
            rows.append(r)
    rows.insert(len(rows) // 2, {1: 5, 2: 7})   # master-only: pre-leftover
    return lap.System(rows, p, order, forbid={0, 1, 2},
                      stratum_of=stratum_of)


# ------------------------------------------------- N1/N2/N4: certify+verify
fill_lines = []
for p in PRIMES:
    for pol in ("B2F", "B2FT"):
        sys_ = synth_system(p, seed=41 + p % 97)
        base = lap.eliminate(sys_, lap.Schedule(policy=pol))
        targets = sorted(base.subs)
        check(f"N1_solved_pivots_nonempty::{pol}_p{p}", len(targets) >= 20,
              f"{len(targets)} pivots")
        wd = f"{OUT}/wit_{pol}_{p}"
        res = lap.eliminate(sys_, lap.Schedule(policy=pol),
                            witnesses="native", witness_targets=targets,
                            witness_dir=wd)
        ok_cert = all(res.witnesses[t]["status"] == "CERTIFIED"
                      for t in targets)
        ok_ver = all(res.witnesses[t].get("verify_pass") for t in targets)
        check(f"N1_native_all_certified::{pol}_p{p}", ok_cert)
        check(f"N1_native_all_verified::{pol}_p{p}", ok_ver)
        led = res.ledger["witnesses"]
        fillrec = led.get("transcript_fill") or {}
        check(f"N4_ledger_fill_measured::{pol}_p{p}",
              led["mode"] == "native"
              and led["n_certified"] == len(targets)
              and fillrec.get("n_registered", 0) > 0
              and fillrec.get("peak_lam_nnz", 0) > 0, fillrec)
        fill_lines.append((pol, p, len(sys_.rows), len(targets), fillrec))

        # N2: retrofit certifies the SAME coefficient rows
        res_r = lap.eliminate(sys_, lap.Schedule(policy=pol),
                              witnesses="retrofit", witness_targets=targets)
        agree = all(
            res_r.witnesses[t]["status"] == "CERTIFIED"
            and res_r.witnesses[t]["c"] == res.witnesses[t]["c"]
            for t in targets)
        check(f"N2_native_matches_retrofit::{pol}_p{p}", agree)

        # N3: file round-trip against an INDEPENDENT copy of the rows
        t0 = targets[len(targets) // 2]
        wp = res.witnesses[t0].get("witness_path")
        rows_indep = [dict(r) for r in sys_.rows]
        ok3, _ = receipt.verify(rows_indep, wp)
        check(f"N3_witness_file_roundtrip::{pol}_p{p}",
              bool(wp) and bool(ok3))

        # N5: MUST-FAIL — tamper a row the witness's lambda actually
        # references (a lambda certificate binds exactly the rows in its
        # support; rows outside it are certified by THEIR witnesses)
        wit0 = receipt.load_witness(wp)
        i0 = wit0["lam_idx"][0]
        rows_bad = [dict(r) for r in sys_.rows]
        c0 = sorted(rows_bad[i0])[0]
        rows_bad[i0][c0] = (rows_bad[i0][c0] + 1) % p or 1
        bad_ok, _ = receipt.verify(rows_bad, wp)
        check(f"N5_tampered_rows_fail::{pol}_p{p}", not bad_ok)
        # N5: MUST-FAIL — tampered lambda entry
        wit_bad = receipt.load_witness(wp)
        wit_bad = copy.deepcopy(wit_bad)
        wit_bad["lam_val"][0] = (wit_bad["lam_val"][0] + 1) % p or 1
        bad_ok2, _ = receipt.verify(rows_indep, wit_bad)
        check(f"N5_tampered_lambda_fails::{pol}_p{p}", not bad_ok2)

print("\ntranscript fill, measured this run "
      "(policy, p, n_rows, n_targets, fill):", flush=True)
for pol, p, nr, nt, f in fill_lines:
    print(f"  {pol} p={p} rows={nr} targets={nt} "
          f"registered={f.get('n_registered')} lam_ops={f.get('lam_ops')} "
          f"lam_nnz_total={f.get('lam_nnz_total')} "
          f"peak_lam_nnz={f.get('peak_lam_nnz')}", flush=True)

# ------------------------------------------------------------- N6: refusals
p = PRIMES[0]
sys_ = synth_system(p, seed=7)
try:
    lap.eliminate(sys_, lap.Schedule(policy="B2"), witnesses="native",
                  witness_targets=[89])
    check("N6_native_B2_refused", False, "no exception")
except NotImplementedError:
    check("N6_native_B2_refused", True)
try:
    lap.eliminate(sys_, lap.Schedule(policy="B2FT"), witnesses="native",
                  witness_targets=[89], backend="dense")
    check("N6_native_dense_refused", False, "no exception")
except NotImplementedError:
    check("N6_native_dense_refused", True)
try:
    lap.eliminate(sys_, lap.Schedule(policy="B2FT"), witnesses="native",
                  witness_targets=[1])          # a master, never a pivot
    check("N6_unsolved_target_refused", False, "no exception")
except ValueError:
    check("N6_unsolved_target_refused", True)

finish("test_native_witness")

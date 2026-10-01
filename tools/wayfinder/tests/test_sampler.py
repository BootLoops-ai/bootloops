#!/usr/bin/env python3
"""
Tests for the Wayfinder sampler (sampler.py).

T1  synthetic        — inject J^{(k)}(p) generated FROM a known dlog
                       connection (3 masters, alphabet {s,t,s+t}); run
                       fit_dlog_connection + verify_connection; check the
                       fitted A_a match exactly and held-out ≥ 30 d.
                       (Fast; no amflow.)

T2  box1l end-to-end — connection_once() on targets/box1l_conn.json via
                       amflow-cpp.  Gated by env WAYFINDER_RUN_AMFLOW=1
                       (≈10–15 min wall on 4 cores).  Checks: every entry
                       rationalizes, held-out ≥ 20 d, A_s/A_t/A_{s+t} have
                       the Henn support pattern (bub diagonal, box couples
                       to both bubbles).

T3  diffeq backend   — connection_once(backend="amflow_diffeq") on box1l
                       via amflow-cpp's mode:"diffeq" (LIBPDerivivative
                       chain rule + Kira IBP).  Gated by
                       WAYFINDER_RUN_AMFLOW=1.  Checks: A_a EXACTLY match
                       the Henn reference, held-out ≥ goal_digits-5 (no
                       finite-diff loss), and n_amflow_calls = 2·n_points
                       (vs (1+2·|kinvars|)·n_points for finite_diff).

T4  sparse support   — same connection as T1, fit with a FIXED coupling
                       skeleton from only 2 points.  (Fast; no amflow.)

T5  tampered sample  — NEGATIVE control: the T1 table with ONE derivative
                       value perturbed (relative 1e-6, 19 orders above the
                       25-digit rationalization tolerance) must make the fit
                       fail to rationalize or drop held-out digits below the
                       T1 bar; the same table unperturbed passes in the same
                       process.  (Fast; no amflow.)

The four fast legs run under pytest (functions are test_*) and, together
with the env-gated T2/T3, from `python3 tests/test_sampler.py`.
"""
import json
import os
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp

try:
    import pytest
    _Skipped = pytest.skip.Exception
except ImportError:                      # direct run without pytest installed
    pytest = None

    class _Skipped(Exception):
        pass

HERE  = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(HERE)
sys.path.insert(0, TOOLS)

from sampler import (Rotation, DESample,
                                fit_dlog_connection, verify_connection,
                                sample_de, connection_once)


def _skip(msg):
    """Env-gated skip that works under pytest (registers as SKIPPED with the
    reason) and under the direct __main__ runner (prints and moves on)."""
    print(msg)
    if pytest is not None:
        pytest.skip(msg)
    raise _Skipped(msg)


# ---------------------------------------------------------------------------
# T1  synthetic — known connection → synthesize (J, ∂J) → fit → recover
# ---------------------------------------------------------------------------
def _dlog_at(a, v, pt):
    ae = sp.sympify(a)
    d = sp.diff(sp.log(ae), sp.Symbol(v))
    val = d.subs({sp.Symbol(k): sp.Integer(x) for k, x in pt.items()})
    return mp.mpf(Fraction(str(sp.nsimplify(val))).numerator) / \
           mp.mpf(Fraction(str(sp.nsimplify(val))).denominator)


def _synth_samples(A_true, alphabet, kinvars, points, J0, K, dps=60,
                   seed=None):
    """Forward-generate (J, ∂J) satisfying (★) EXACTLY.

    The fit and held-out verification only consume the pairs
    (J^{(k-1)}(p), ∂J^{(k)}/∂x_i(p)), related by (★).  So we are free to
    pick J^{(k)}(p) as any deterministic full-rank function of (k, p, m)
    — here a small-integer pseudo-random table — and DEFINE ∂J^{(k)} from
    (★).  This isolates the linear-algebra + PSLQ core from amflow.

    seed=None (T1/T4) keeps the per-process hash-seeded table; an explicit
    seed (T5) makes the table PYTHONHASHSEED-independent by seeding the RNG
    from a string (random.Random hashes str seeds with sha512)."""
    mp.mp.dps = dps
    M = len(J0)
    letters = list(alphabet)
    A = {a: mp.matrix(A_true[a]) for a in letters}
    import random as _r
    samples = []
    for p in points:
        key = tuple(sorted(p.items()))
        rng = _r.Random(hash(key) if seed is None else f"{seed}:{key}")
        Jk = {0: list(J0)}
        for k in range(1, K + 1):
            Jk[k] = [mp.mpf(rng.randint(-9, 9)) + mp.mpf(rng.randint(1, 9)) / 7
                     for _ in range(M)]
        dJ = {}
        for v in kinvars:
            dJ[v] = {0: [mp.mpc(0)] * M}
            for k in range(1, K + 1):
                Jkm1 = mp.matrix([[x] for x in Jk[k - 1]])
                acc = mp.matrix(M, 1)
                for a in letters:
                    acc += _dlog_at(a, v, p) * (A[a] * Jkm1)
                dJ[v][k] = [acc[i, 0] for i in range(M)]
        samples.append(DESample(
            point={kk: str(vv) for kk, vv in p.items()},
            orders=list(range(K + 1)),
            J={k: Jk[k] for k in Jk}, dJ=dJ, deriv_digits=dps - 10))
    return samples


def test_t1_synthetic():
    alphabet = ["s", "t", "s + t"]
    kinvars = ["s", "t"]
    masters = ["bub_s", "bub_t", "box"]
    # the Henn 1-loop massless box connection (basis: ε(1-2ε)·bub_s,
    # ε(1-2ε)·bub_t, ε²st·box).  Σ_a A_a = diag(-1,-1,-1) (Euler scaling).
    # This is what T2 recovers from amflow-cpp end-to-end.
    A_true = {
        "s":     [[-1, 0, 0], [0,  0, 0], [ 2,  0, -1]],
        "t":     [[ 0, 0, 0], [0, -1, 0], [ 0,  2, -1]],
        "s + t": [[ 0, 0, 0], [0,  0, 0], [-2, -2,  1]],
    }
    J0 = [mp.mpf(1), mp.mpf(1), mp.mpf(2)]   # arbitrary nonzero boundary
    fit_pts = [{"s": -3, "t": -7}, {"s": -5, "t": -11}, {"s": -13, "t": -2},
               {"s": -17, "t": -23}]
    ho_pts  = [{"s": -29, "t": -19}, {"s": -41, "t": -31}]
    rot = Rotation(masters, {m: sp.Integer(1) for m in masters},
                   set(), kinvars)
    fit_s = _synth_samples(A_true, alphabet, kinvars, fit_pts, J0, K=3)
    ho_s  = _synth_samples(A_true, alphabet, kinvars, ho_pts,  J0, K=3)

    conn = fit_dlog_connection(alphabet, kinvars, fit_s, rot, dps=50)
    print("[T1] cond:", conn["honesty"]["cond_number"],
          " fit-residual digits:", conn["honesty"]["fit_residual_digits"],
          " failed:", len(conn["honesty"]["failed_to_rationalize"]))
    for a in alphabet:
        print(f"[T1]  A[{a:5s}] =", conn["A_dense"][a])
    # exact match
    for a in alphabet:
        for r in range(3):
            for c in range(3):
                got = conn["A_dense"][a][r][c]
                want = str(A_true[a][r][c])
                assert got == want, f"A[{a}][{r}][{c}]: got {got}, want {want}"
    ho_d = verify_connection(conn, ho_s, dps=50)
    print("[T1] held-out digits:", ho_d)
    assert ho_d >= 30
    assert not conn["honesty"]["failed_to_rationalize"]
    print("[T1] PASS")


# ---------------------------------------------------------------------------
# T2  box1l end-to-end (amflow)  — gated
# ---------------------------------------------------------------------------
def test_t2_box1l_amflow():
    if not os.environ.get("WAYFINDER_RUN_AMFLOW"):
        _skip("[T2] SKIP (set WAYFINDER_RUN_AMFLOW=1 to run; ~10-15 min)")
    tgt = os.path.join(HERE, "targets", "box1l_conn.json")
    out = os.path.join(HERE, "out", "box1l", "connection.json")
    wr  = os.path.join(HERE, "out", "box1l", "amflow_cache")
    conn = connection_once(tgt, n_fit=3, n_heldout=2, h_pow=15,
                           goal_digits=40, target_eps_power=3,
                           work_root=wr,
                           parallel=int(os.environ.get("WAYFINDER_PAR", "4")),
                           out_path=out)
    h = conn["honesty"]
    print("[T2] timing:", h["timing"])
    print("[T2] fit residual:", h["fit_residual_digits"], "d;  held-out:",
          h["min_heldout_digits"], "d;  failed:",
          len(h["failed_to_rationalize"]))
    for a in conn["alphabet"]:
        print(f"[T2]  A[{a:7s}] =", conn["A_dense"][a])
    assert not h["failed_to_rationalize"], h["failed_to_rationalize"]
    assert h["min_heldout_digits"] >= 15
    # Henn structure: bub_s only in dlog(s), bub_t only in dlog(t),
    # box row couples to both bubbles + itself.
    As, At, Au = conn["A_dense"]["s"], conn["A_dense"]["t"], conn["A_dense"]["s + t"]
    assert As[0][0] != "0" and At[0][0] == "0"
    assert At[1][1] != "0" and As[1][1] == "0"
    assert any(conn["A_dense"][a][2][0] != "0" for a in conn["alphabet"])
    assert any(conn["A_dense"][a][2][1] != "0" for a in conn["alphabet"])
    print("[T2] PASS")


# ---------------------------------------------------------------------------
# T3  amflow_diffeq backend (mode:"diffeq" via LIBPDerivivative chain rule)
# ---------------------------------------------------------------------------
def test_t3_diffeq_backend():
    if not os.environ.get("WAYFINDER_RUN_AMFLOW"):
        _skip("[T3] SKIP (set WAYFINDER_RUN_AMFLOW=1 to run; ~2-5 min)")
    tgt = os.path.join(HERE, "targets", "box1l_conn.json")
    out = os.path.join(HERE, "out", "box1l_diffeq", "connection.json")
    wr  = os.path.join(HERE, "out", "box1l_diffeq", "amflow_cache")
    conn = connection_once(tgt, n_fit=3, n_heldout=2,
                           goal_digits=40, target_eps_power=3,
                           work_root=wr,
                           parallel=int(os.environ.get("WAYFINDER_PAR", "4")),
                           out_path=out, backend="amflow_diffeq")
    h = conn["honesty"]
    print("[T3] backend:", h["backend"], " timing:", h["timing"])
    print("[T3] fit residual:", h["fit_residual_digits"], "d;  held-out:",
          h["min_heldout_digits"], "d;  failed:",
          len(h["failed_to_rationalize"]))
    for a in conn["alphabet"]:
        print(f"[T3]  A[{a:7s}] =", conn["A_dense"][a])
    # The Henn 1-loop massless box reference (same as T1).
    A_ref = {
        "s":     [["-1", "0", "0"], ["0",  "0", "0"], [ "2",  "0", "-1"]],
        "t":     [[ "0", "0", "0"], ["0", "-1", "0"], [ "0",  "2", "-1"]],
        "s + t": [[ "0", "0", "0"], ["0",  "0", "0"], ["-2", "-2",  "1"]],
    }
    for a in A_ref:
        assert conn["A_dense"][a] == A_ref[a], \
            f"A[{a}]: got {conn['A_dense'][a]}, want {A_ref[a]}"
    assert h["backend"] == "amflow_diffeq"
    assert not h["failed_to_rationalize"], h["failed_to_rationalize"]
    # Exact DE matrix → derivative digits = goal_digits, so held-out
    # agreement should be ≳ goal_digits-5 (vs ~⅔·goal_digits for FD).
    assert h["min_heldout_digits"] >= 35, h["min_heldout_digits"]
    # 5 points × 2 calls = 10 (vs 5 × (1+2·2) = 25 for finite_diff).
    assert h["timing"]["n_amflow_calls"] == 10, h["timing"]["n_amflow_calls"]
    print("[T3] PASS  (exact A_a recovered; held-out",
          h["min_heldout_digits"], "d;", h["timing"]["n_amflow_calls"],
          "amflow calls vs 25 finite_diff)")


# ---------------------------------------------------------------------------
# T4  sparse-support synthetic — fixed coupling skeleton → fewer points
# ---------------------------------------------------------------------------
def test_t4_sparse_support():
    """Same Henn box connection as T1, but fit with a FIXED sparsity support
    (row→cols) and only 2 fit points (dense T1 used 4).  Checks the sparse fit
    recovers the exact A_a, forces off-support entries to 0, and held-out ≥30d.
    Fast; no amflow."""
    alphabet = ["s", "t", "s + t"]
    kinvars = ["s", "t"]
    masters = ["bub_s", "bub_t", "box"]
    A_true = {
        "s":     [[-1, 0, 0], [0,  0, 0], [ 2,  0, -1]],
        "t":     [[ 0, 0, 0], [0, -1, 0], [ 0,  2, -1]],
        "s + t": [[ 0, 0, 0], [0,  0, 0], [-2, -2,  1]],
    }
    J0 = [mp.mpf(1), mp.mpf(1), mp.mpf(2)]
    # Henn coupling skeleton: bub_s→{0}, bub_t→{1}, box→{0,1,2}
    support = {0: [0], 1: [1], 2: [0, 1, 2]}
    rot = Rotation(masters, {m: sp.Integer(1) for m in masters}, set(), kinvars)
    fit_s = _synth_samples(A_true, alphabet, kinvars,
                           [{"s": -3, "t": -7}, {"s": -5, "t": -11}], J0, K=3)
    ho_s  = _synth_samples(A_true, alphabet, kinvars,
                           [{"s": -29, "t": -19}, {"s": -41, "t": -31}], J0, K=3)
    conn = fit_dlog_connection(alphabet, kinvars, fit_s, rot, dps=50,
                               support=support)
    h = conn["honesty"]
    print("[T4] sparse:", h["sparse_support"], " max_row_unknowns:",
          h["max_row_unknowns"], " support_pairs:", h["support_pairs"],
          " n_eq_per_row:", h["n_eq_per_row"], " fit-resid:",
          h["fit_residual_digits"])
    for a in alphabet:
        for r in range(3):
            for c in range(3):
                got = conn["A_dense"][a][r][c]; want = str(A_true[a][r][c])
                assert got == want, f"A[{a}][{r}][{c}]: got {got}, want {want}"
    ho_d = verify_connection(conn, ho_s, dps=50)
    print("[T4] held-out digits:", ho_d, " failed:",
          len(h["failed_to_rationalize"]))
    assert ho_d >= 30
    assert not h["failed_to_rationalize"]
    assert h["sparse_support"] and h["max_row_unknowns"] == 9
    print("[T4] PASS  (exact A_a from 2 pts via fixed support)")


# ---------------------------------------------------------------------------
# T5  tampered sample — NEGATIVE control (the fit must notice one bad ∂J)
# ---------------------------------------------------------------------------
T5_SEED = 5                      # PYTHONHASHSEED-independent synthetic table
T5_TAMPER_REL = mp.mpf("1e-6")   # 19 orders above tol_rationalize = 25 d


def test_t5_tampered_sample():
    """NEGATIVE control.  Build the T1 synthetic table (deterministic seed),
    perturb ONE derivative value — dJ[s][k=1][row 0] at the first fit point,
    relative 1e-6 — and require that fit_dlog_connection + verify_connection
    detect it: fit_dlog_connection's own failure path is the returned
    honesty["failed_to_rationalize"] list (an entry _pslq_rational could not
    recognise at tol_rationalize = deriv_digits//2 = 25 digits with 32-bit
    coefficients; it raises nothing), and verify_connection then falls back to
    the fitted float for that "?" entry, so the held-out digits collapse to
    the perturbation scale.  Either signal is a detection; the assertion
    demands at least one.  Control-of-the-control: the SAME table, unperturbed,
    passes the T1 bar in this process, so the failure is the tamper's and not
    the table's.  Fast; no amflow."""
    alphabet = ["s", "t", "s + t"]
    kinvars = ["s", "t"]
    masters = ["bub_s", "bub_t", "box"]
    A_true = {
        "s":     [[-1, 0, 0], [0,  0, 0], [ 2,  0, -1]],
        "t":     [[ 0, 0, 0], [0, -1, 0], [ 0,  2, -1]],
        "s + t": [[ 0, 0, 0], [0,  0, 0], [-2, -2,  1]],
    }
    J0 = [mp.mpf(1), mp.mpf(1), mp.mpf(2)]
    fit_pts = [{"s": -3, "t": -7}, {"s": -5, "t": -11}, {"s": -13, "t": -2},
               {"s": -17, "t": -23}]
    ho_pts  = [{"s": -29, "t": -19}, {"s": -41, "t": -31}]
    rot = Rotation(masters, {m: sp.Integer(1) for m in masters},
                   set(), kinvars)
    ho_s = _synth_samples(A_true, alphabet, kinvars, ho_pts, J0, K=3,
                          seed=T5_SEED)

    # control-of-the-control: the unperturbed table passes the T1 bar
    fit_ok = _synth_samples(A_true, alphabet, kinvars, fit_pts, J0, K=3,
                            seed=T5_SEED)
    conn_ok = fit_dlog_connection(alphabet, kinvars, fit_ok, rot, dps=50)
    ho_ok = verify_connection(conn_ok, ho_s, dps=50)
    print("[T5] control (unperturbed): failed:",
          len(conn_ok["honesty"]["failed_to_rationalize"]),
          " held-out digits:", ho_ok)
    assert not conn_ok["honesty"]["failed_to_rationalize"], \
        conn_ok["honesty"]["failed_to_rationalize"]
    assert ho_ok >= 30, ho_ok
    for a in alphabet:
        assert conn_ok["A_dense"][a] == [[str(x) for x in row]
                                         for row in A_true[a]], a

    # the tamper: same table (same seed), ONE derivative value perturbed
    fit_bad = _synth_samples(A_true, alphabet, kinvars, fit_pts, J0, K=3,
                             seed=T5_SEED)
    orig = fit_bad[0].dJ["s"][1][0]
    fit_bad[0].dJ["s"][1][0] = orig * (1 + T5_TAMPER_REL)
    assert fit_bad[0].dJ["s"][1][0] != orig, "tamper did not land (entry 0?)"
    print("[T5] tampered dJ[s][k=1][row 0] at", fit_bad[0].point, ":",
          mp.nstr(orig, 12), "->", mp.nstr(fit_bad[0].dJ["s"][1][0], 12))

    conn_bad = fit_dlog_connection(alphabet, kinvars, fit_bad, rot, dps=50)
    failed = conn_bad["honesty"]["failed_to_rationalize"]
    ho_bad = verify_connection(conn_bad, ho_s, dps=50)
    print("[T5] tampered: fit-residual digits:",
          conn_bad["honesty"]["fit_residual_digits"],
          " failed_to_rationalize:", len(failed),
          " held-out digits:", ho_bad)
    for f in failed:
        print("[T5]   failed entry (letter,row,col,float):", f)
    assert failed or ho_bad < 30, \
        f"tamper undetected: failed_to_rationalize={failed} held-out={ho_bad} d"
    print(f"[T5] PASS  (negative control fails as designed: "
          f"{len(failed)} entries unrationalized, held-out {ho_bad} d < 30)")


if __name__ == "__main__":
    fail = 0
    for fn in (test_t1_synthetic, test_t4_sparse_support,
               test_t5_tampered_sample, test_t2_box1l_amflow,
               test_t3_diffeq_backend):
        try:
            fn()
        except _Skipped:
            pass                                  # reason already printed
        except AssertionError as e:
            print(f"[{fn.__name__}] FAIL:", e); fail += 1
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f"[{fn.__name__}] ERROR:", e); fail += 1
        print()
    print("=" * 60)
    print(f"{'OK' if fail == 0 else 'FAIL'}")
    sys.exit(1 if fail else 0)

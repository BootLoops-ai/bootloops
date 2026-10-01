"""Receipt emission with budget-aware backend routing.

Measured failure this module exists to prevent: retrofit-witness receipts
routed through the dense-torch backend OOM-crashed three production cells — a
6.98 GB stratum buffer was requested under a standing 31G `ulimit -v` while
the nominal dense_budget (8 GB) said "fine". Root cause: a nominal budget
check that is not VA-headroom-aware — and a receipt pass has no business
dying on an allocator.

Contract of this module:
  * choose_backend() routes dense vs cpu by the LARGEST-STRATUM dense-buffer
    ESTIMATE against min(nominal budget, safety * measured VA headroom).
  * emit_receipts() NEVER dense-OOMs a receipt pass: any MemoryError /
    torch-allocator RuntimeError on a dense call is caught once and the call
    is re-run on the engine ("cpu") path, loudly recorded. A cpu-path failure
    still raises (fail-closed) — receipts are certification, not best-effort.
  * The returned record always NAMES its evidence (witness_dir) and carries
    the measured pass walls (pass1_wall_s, pass2_wall_s) so downstream
    projections are rate-fit from the tool's own record, never guessed
    (pass2 is NOT pass1-sized; it was measured at ~2-3x pass1 on the
    production systems — never project pass2 walls from pass1).

ibplapper is an injected dependency (same policy as the cli hooks); the
eliminate/load callables are parameters so tests exercise the routing without
a real solve.
"""
import os
import resource
import time

_DENSE_CELL_BYTES = 8          # torch int64/f64 cell
_DEFAULT_SAFETY = 0.5          # fraction of measured VA headroom dense may use


def va_headroom_bytes():
    """RLIMIT_AS minus current VmSize; None when the limit is unlimited."""
    soft, _hard = resource.getrlimit(resource.RLIMIT_AS)
    if soft in (resource.RLIM_INFINITY, -1):
        return None
    vmsize = None
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmSize:"):
                    vmsize = int(line.split()[1]) * 1024
                    break
    except OSError:
        return None
    if vmsize is None:
        return None
    return max(0, soft - vmsize)


def stratum_dense_estimates(system):
    """Per-stratum dense-buffer ESTIMATE: rows_in_stratum * distinct_cols * 8.

    system must expose .rows (iterable of col->coeff mappings) and
    .stratum_of (mapping col -> stratum id). Rows are grouped by the min
    stratum id among their eliminable columns (routing estimate only — the
    exact backend block layout is the backend's business; this only has to
    be the right order of magnitude, and errs conservative via the caller's
    safety factor)."""
    st_rows, st_cols = {}, {}
    stratum_of = system.stratum_of
    for row in system.rows:
        sts = [stratum_of[c] for c in row if c in stratum_of]
        if not sts:
            continue
        s = min(sts)
        st_rows[s] = st_rows.get(s, 0) + 1
        st_cols.setdefault(s, set()).update(row.keys())
    return {s: st_rows[s] * len(st_cols[s]) * _DENSE_CELL_BYTES
            for s in st_rows}


def choose_backend(system, budget_bytes, safety=_DEFAULT_SAFETY):
    """Return (backend, info). backend 'dense' only when the largest-stratum
    estimate fits within min(budget, safety * VA headroom)."""
    est = stratum_dense_estimates(system)
    max_est = max(est.values()) if est else 0
    headroom = va_headroom_bytes()
    cap = budget_bytes
    if headroom is not None:
        cap = min(cap, int(safety * headroom))
    backend = "dense" if max_est <= cap else "cpu"
    return backend, {
        "max_stratum_bytes_est": max_est,
        "budget_bytes": budget_bytes,
        "va_headroom_bytes": headroom,
        "effective_cap_bytes": cap,
        "safety": safety,
        "n_strata_est": len(est),
    }


def _is_alloc_failure(exc):
    if isinstance(exc, MemoryError):
        return True
    return isinstance(exc, RuntimeError) and (
        "allocate" in str(exc).lower() or "out of memory" in str(exc).lower())


def emit_receipts(system, witness_dir, n_targets=2, policy="B2FT",
                  budget_bytes=8 * 2**30, safety=_DEFAULT_SAFETY,
                  eliminate=None, schedule=None):
    """Two-pass retrofit receipts with routed backend + cpu fallback.

    system: ibplapper System (or duck-typed equivalent); witness_dir: where
    lambda-witness files go (ALWAYS named in the returned record).
    eliminate/schedule: injectable for tests; default = ibplapper's.

    Returns record dict: {n, ok, witness_dir, backend_requested,
    backend_used, backend_choice, pass1_wall_s, pass2_wall_s, fallbacks}.
    Raises on cpu-path failure (fail-closed)."""
    if eliminate is None or schedule is None:
        import ibplapper as lap
        eliminate = eliminate or lap.eliminate
        schedule = schedule or lap.Schedule
    backend, info = choose_backend(system, budget_bytes, safety)
    fallbacks = []

    def _run(**kw):
        nonlocal backend
        try:
            return eliminate(system, schedule(policy=policy),
                             backend=backend, device="cpu",
                             dense_budget=budget_bytes, **kw)
        except Exception as e:                       # noqa: BLE001
            if backend == "dense" and _is_alloc_failure(e):
                fallbacks.append({"from": "dense", "to": "cpu",
                                  "error": repr(e)[:200]})
                backend = "cpu"
                return eliminate(system, schedule(policy=policy),
                                 backend=backend, device="cpu",
                                 dense_budget=budget_bytes, **kw)
            raise                                     # fail-closed

    t0 = time.time()
    res0 = _run()
    tg = list(res0.subs)[:n_targets]
    res0 = None                                       # free pass1 (VA headroom)
    t1 = time.time()
    os.makedirs(witness_dir, exist_ok=True)
    res = _run(witnesses="retrofit", witness_targets=tg,
               witness_dir=witness_dir)
    t2 = time.time()
    ok = all(res.witnesses[t]["status"] == "CERTIFIED"
             and res.witnesses[t].get("verify_pass") for t in tg)
    return {"n": len(tg), "ok": bool(ok),
            "witness_dir": os.path.abspath(witness_dir),
            "backend_requested": ("dense" if info["max_stratum_bytes_est"]
                                  <= info["effective_cap_bytes"] else "cpu"),
            "backend_used": backend + ("-fallback" if fallbacks else ""),
            "backend_choice": info, "fallbacks": fallbacks,
            "pass1_wall_s": round(t1 - t0, 1),
            "pass2_wall_s": round(t2 - t1, 1)}

#!/usr/bin/env python3
r"""
PMflow — self-consistent auxiliary-mass-flow CLI for cut-eikonal (PM) families.

Design: tools/pmflow/DESIGN.md
Solver:  tools/pmflow/gf_solve.py  (do not modify)

Subcommands
-----------
detect   --json FAMILY.json --etac CSV
         Run amflow_cli briefly (timeout 900, AMFLOW_DEBUG_SCHEME=1, AMFLOW_ETAC_OVERRIDE=CSV).
         Parse "[scheme]" lines from the log; report boundary family names per level and
         whether successive levels have identical names-modulo-suffix (fixed point) or descend.

discover --json CHAIN.json --keys-out keys.json [--etac CSV] [--force-depth N]
         Zero-probe iteration (from eb_discover.sh): run with all-zero explicit_boundary,
         parse MISS keys, append, repeat (max 15) until solved; write final key list.
         DUAL-SURFACE: both error shapes are MISS surfaces —
           * "no Vacuum entry" master keys (classic surface, kind "vacuum";
             dual-shape regex: current binary prints the key BEFORE the phrase);
           * "GRAVITYFLOW_CUTREGION" guard aborts (cut-quadratic families, e.g. a cut
             graviton line): recorded as a kind "cutregion" family marker, then discover
             switches to the depth-0 leaf-Trivial injection surface
             (AMFLOW_FORCE_ENDING_DEPTH=0) and the masters found there carry kind "cut".
         The guard only fires on the normal path (forced ending depth >= 1 bypasses
         build_boundary), so cut-quadratic input JSONs get one normal-path probe first.
         Emitted key list: flat ["key", ...] when every key is kind "vacuum"
         (byte-compatible with the legacy flat format); otherwise tagged
         [{"key": ..., "kind": "vacuum"|"cut"|"cutregion"}, ...] so respond can route
         cut keys to the injection path and vacuum keys to the classic path.

respond  --json CHAIN.json --keys keys.json --dir PROBEDIR --par N
         Generate per-key unit-injection JSONs and run them N-way parallel with
         AMFLOW_DUMP_EPS_GRID=1. Each run timeout 1800.

solve    ...
         Thin wrapper delegating to gf_solve.py (imported, same args).

inject   --parent PARENT.json --table fp_boundary.json --depth 1 --out OUT.json
         Build parent JSON with explicit_boundary and run amflow_cli with
         AMFLOW_FORCE_ENDING_DEPTH.

map      --keys keys.json --families ending_families.json --fpA fp_A.json --out circular_map.json
         Build key -> b0-master-index mapping by matching used-prop content (normalized) to the
         basis family's prop list.  Parent prop content comes from --parent-json CHAIN.json
         (family.propagators; same cfg the other subcommands take) or an ending_families entry
         named like the parent.  Unmatched keys are printed; they default to ['J63', None]
         only with --foreign-as-j63.

All subcommands accept:
  --amflow-cli PATH   path to amflow_cli binary (default: $AMFLOW_CLI, else amflow_cli on PATH)
  --jemalloc-wrap PATH path to jemalloc wrapper (default: $AMFLOW_JEMALLOC_WRAP; none if unset)
  --cache PATH        AMFLOW_IBP_CACHE directory
"""

import argparse
import concurrent.futures
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "amflow-kit")); from amflow_kit import memfence  # noqa: E402  (tools/amflow-kit/amflow_kit/memfence.py)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
# Defaults come from the environment:
#   AMFLOW_CLI            amflow-cpp CLI binary (build from the sibling repository ../amflow-cpp)
#   AMFLOW_JEMALLOC_WRAP  optional jemalloc wrapper script; unset/empty = run the CLI bare
#   FERMATPATH            fer64 binary;  AMFLOW_IBP_CACHE  IBP cache dir
DEFAULT_CLI = os.environ.get("AMFLOW_CLI", "amflow_cli")
DEFAULT_JW = os.environ.get("AMFLOW_JEMALLOC_WRAP", "")
DEFAULT_FERMAT = os.environ.get("FERMATPATH", "fer64")
DEFAULT_IBP_CACHE = os.environ.get("AMFLOW_IBP_CACHE", "")
DEFAULT_ETAC = "0,0,0,0,0,0,-1,-1,0,0,0,0,0,0,0"

# --- discover MISS surfaces (dual-surface) ---------------------------------
# Current binary (amfsystem.cpp:963,996) prints the key BEFORE the phrase:
#   "ending system master 'FAM|i|...' has no Vacuum entry, no explicit_boundary
#    value, and reduction did not lower it"
# The older binary shape (key AFTER the phrase) is kept as a fallback — a
# single-shape regex is blind on the other shape's logs.
MISS_RES = (
    re.compile(r"ending system master '([^']+)' has no Vacuum entry"),
    re.compile(r"no Vacuum entry[^']*'([^']+)'"),
)
# CUTREGION guard abort (amfsystem.cpp build_boundary guard):
#   "... GRAVITYFLOW_CUTREGION: family 'FAM' carries a cut quadratic ..."
CUTREGION_RE = re.compile(r"GRAVITYFLOW_CUTREGION: family '([^']+)'")
# FIXEDPOINT guard abort (amfsystem.cpp setup() self-similarity guard —
# fires on plain runs when a boundary family reproduces its parent's
# propagator set):
#   "... GRAVITYFLOW_FIXEDPOINT: boundary family 'FAM' at depth D
#    reproduces its parent's propagator set ..."
# `detect` treats it as the engine-certified fixed-point verdict.
FIXEDPOINT_RE = re.compile(
    r"GRAVITYFLOW_FIXEDPOINT: boundary family '([^']+)' at depth (\d+)")
# Engine abort when AMFLOW_ETAC_OVERRIDE length != #propagators.
ETAC_MISMATCH_MARK = "AMFLOW_ETAC_OVERRIDE: length"


def _parse_miss(log_text):
    """Dual-shape 'no Vacuum entry' MISS key, or None."""
    for rex in MISS_RES:
        m = rex.search(log_text)
        if m:
            return m.group(1)
    return None


def _resolve_etac(args, tag):
    """--etac semantics:
       CSV  -> AMFLOW_ETAC_OVERRIDE=CSV;
       ''   -> engine-chosen etac (no override; what cut-quadratic families need);
       omitted -> DEPRECATED legacy 15-slot DEFAULT_ETAC, with a loud warning
                  (kept for backward compatibility with earlier runs)."""
    etac = getattr(args, "etac", None)
    if etac is None:
        print(f"[{tag}] DEPRECATION: --etac not given; falling back to the legacy "
              f"15-slot DEFAULT_ETAC ({DEFAULT_ETAC}). On any family whose "
              f"propagator count is not 15, amflow aborts with an "
              f"AMFLOW_ETAC_OVERRIDE length mismatch. "
              f"Pass --etac CSV explicitly, or --etac '' for the "
              f"engine-chosen etac.")
        return DEFAULT_ETAC
    return etac


def _etac_env(etac):
    """Env fragment for an etac choice; '' means no override (engine treats
    an empty AMFLOW_ETAC_OVERRIDE as unset, amfsystem.cpp:4253)."""
    return {} if etac == "" else {"AMFLOW_ETAC_OVERRIDE": etac}


def _load_key_entries(path):
    """Read a keys.json in either format:
       legacy flat  ["FAM|i|...", ...]                (all kind 'vacuum')
       tagged       [{"key":..., "kind":...}, ...]    (dual-surface discover)
    Returns dict key -> kind."""
    raw = json.loads(Path(path).read_text())
    entries = {}
    for e in raw:
        if isinstance(e, str):
            entries.setdefault(e, "vacuum")
        else:
            entries.setdefault(e["key"], e.get("kind", "vacuum"))
    return entries


def _emit_key_entries(path, entries):
    """Write keys.json. All-vacuum -> legacy flat sorted list (byte-compatible
    with earlier key lists and gf_solve); any cut/cutregion kind -> tagged list
    sorted by key, so respond can route by kind."""
    if all(k == "vacuum" for k in entries.values()):
        Path(path).write_text(json.dumps(sorted(entries), indent=0))
    else:
        Path(path).write_text(json.dumps(
            [{"key": k, "kind": entries[k]} for k in sorted(entries)],
            indent=0))


def _family_cut_quadratic(cfg):
    """Mirror of the binary's GRAVITYFLOW_CUTREGION criterion on the input
    JSON: any CUT propagator whose loop-total-degree is >= 2 (a cut graviton
    line (k1+k2-q)^2; eikonal 2u.k cut props are degree 1 and fine).
    Returns True/False, or None when the JSON cannot be classified (missing
    fields, sympy unavailable, unparsable propagator) — callers must treat
    None as 'unknown', keeping classic behavior. The binary's guard remains
    the ground truth; this only decides whether to SPEND a normal-path probe."""
    try:
        fam = cfg["family"]
        loops = list(fam["loops"])
        props = list(fam["propagators"])
        cut = list(fam.get("cut", []))
        if not any(cut):
            return False
        import sympy
        idents = set()
        for p in props:
            idents.update(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", p))
        syms = {name: sympy.Symbol(name) for name in idents}
        loop_syms = [syms[l] for l in loops if l in syms]
        if not loop_syms:
            return None
        for i, p in enumerate(props):
            if i >= len(cut) or not cut[i]:
                continue
            expr = sympy.expand(
                sympy.sympify(p.replace("^", "**"), locals=syms))
            if not (set(loop_syms) & expr.free_symbols):
                continue
            deg = sympy.Poly(expr, *loop_syms).total_degree()
            if deg >= 2:
                return True
        return False
    except Exception:
        return None

# Standard env block — DEFAULTS only:
# the caller's environment wins over every entry
# here (passthrough-with-defaults). A {**os.environ,
# **BASE_ENV} merge would silently discard exported staging vars (measured:
# AMFLOW_REDUCE_DCAP=1 in the shell arrived as 0 at the binary).
# Subcommand semantics (--cache, forced depth, etac, grid dump)
# still override the caller via the `extra` dict.
BASE_ENV = {
    "FERMATPATH": DEFAULT_FERMAT,
    "AMFLOW_KIRA_MEM_CAP_GB": "60",
    "AMFLOW_IBP_CACHE": DEFAULT_IBP_CACHE,
    "AMFLOW_KIRA_P_AUTO": "0",
    "AMFLOW_NEPS_OVERRIDE": "12",
    "AMFLOW_EPS_PARALLEL": "1",
    "AMFLOW_REDUCE_DCAP": "0",
    "AMFLOW_NO_AUTO_FIREFLY": "0",
    "OPENBLAS_NUM_THREADS": "2",
}


def _env(args, extra=None):
    """BASE_ENV defaults under the caller's environment (caller wins), then
    --cache override and the optional extra dict (subcommand semantics win)."""
    e = {**{k: v for k, v in BASE_ENV.items() if v}, **os.environ}
    if getattr(args, "cache", None):
        e["AMFLOW_IBP_CACHE"] = args.cache
    if extra:
        e.update(extra)
    # memfence: the cap in e is the memory BUDGET (BASE_ENV default 60 GiB,
    # caller/extra may override). Inside a finite-memory.max cgroup leaf the
    # helper exports '0' (the binary's no-cap value; the leaf is the fence);
    # otherwise it exports the 2.5 x budget ADDRESS-SPACE fuse, which the
    # binary applies as RLIMIT_AS to its Kira child — never the 1x budget.
    memfence.apply_policy_to_env(e, memfence.fuse_policy(e.get(memfence.CAP_ENV)))
    return e


def _cli(args):
    return getattr(args, "amflow_cli", DEFAULT_CLI)


def _jw(args):
    return getattr(args, "jemalloc_wrap", DEFAULT_JW)


def _run_fenced(argv, **kw):
    """subprocess.run with the memfence env readback. Popen the child, read
    its /proc/<pid>/environ — and, after a short settle, every descendant's
    (timeout(1) -> jemalloc wrapper -> the exec'd amflow_cli; a wrapper that
    drops an export shows only there) — against the cap fields of the env we
    passed (memfence.intended_from_env), write the readback beside the child's log
    as <log>.memfence.json (pmflow has no other per-run receipt surface), and
    on a mismatch abort BY PID — the pid Popen returned, cmdline-checked
    against argv, never a pattern — raising RuntimeError with the diff.
    Otherwise wait and return a CompletedProcess(argv, returncode). An
    optional timeout= kwarg is honoured as subprocess.run's would be, mapped
    to returncode 124 after the child is killed by pid."""
    timeout = kw.pop("timeout", None)
    env = kw.get("env") or os.environ
    p = subprocess.Popen(list(argv), **kw)
    rb = memfence.env_readback_tree(p.pid, memfence.intended_from_env(env))
    out = getattr(kw.get("stdout"), "name", None)
    side = (out + ".memfence.json") if isinstance(out, str) else None
    if side:
        Path(side).write_text(json.dumps(
            {"argv": list(argv), "pid": p.pid, "fence_mode": env.get(memfence.MODE_ENV),
             "cap_export": env.get(memfence.CAP_ENV), "readback": rb},
            indent=1, default=str))
    if rb["diff"]:
        abs_ = memfence.abort_tree_by_pid(p.pid, list(argv), reason="env readback mismatch")
        try:
            p.wait(timeout=5)
        except Exception:
            pass
        raise RuntimeError(
            f"memfence: env readback mismatch for pid {p.pid}: {rb['diff']} "
            f"(killed by pid, tree of {len(abs_)}: gone={[a['gone'] for a in abs_]}); see {side}")
    try:
        rc = p.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        memfence.abort_by_name(p.pid, list(argv), reason=f"wall clamp {timeout}s")
        p.wait()
        rc = 124
    return subprocess.CompletedProcess(list(argv), rc)


def _wrap_cmd(args, *tail):
    """[jemalloc-wrap] + tail, dropping the wrapper when none is configured."""
    return [c for c in ([_jw(args)] + list(tail)) if c]


def _run_clamped(argv, secs, **kw):
    """subprocess.run under a wall clamp. Prefers timeout(1) (absent on stock
    macOS; coreutils installs it as gtimeout) — identical to the historical
    behavior, rc=124 on expiry. Without either binary, falls back to
    subprocess's own timeout, mapped to the same rc=124 contract, so callers
    see one interface on every platform."""
    exe = shutil.which("timeout") or shutil.which("gtimeout")
    if exe:
        return _run_fenced([exe, str(secs)] + list(argv), **kw)
    return _run_fenced(list(argv), timeout=float(secs), **kw)


# ---------------------------------------------------------------------------
# Subcommand: detect
# ---------------------------------------------------------------------------

def cmd_detect(args):
    """
    Run amflow_cli briefly (timeout 900) with AMFLOW_DEBUG_SCHEME=1 and the given
    AMFLOW_ETAC_OVERRIDE; parse "[scheme]" lines from the log; report boundary family
    names per level and whether successive levels have identical names-modulo-suffix
    (fixed point) or descend.  Print verdict.
    """
    t0 = time.monotonic()
    family_json = Path(args.json).resolve()
    etac = args.etac
    work_dir = Path(args.work_dir).resolve() if args.work_dir else \
        family_json.parent / "wd_detect"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Inject work_dir into the probe JSON.
    probe = json.loads(family_json.read_text())
    probe["work_dir"] = str(work_dir)
    tmp = work_dir / "probe_detect_tmp.json"
    tmp.write_text(json.dumps(probe, indent=1))
    out_json = work_dir / "out_detect.json"
    log_path = work_dir / "log_detect.log"

    env = _env(args, {
        "AMFLOW_DEBUG_SCHEME": "1",
        "AMFLOW_ETAC_OVERRIDE": etac,
    })

    print(f"[detect] family={family_json.name}  etac={etac}")
    print(f"[detect] log → {log_path}")
    r = _run_clamped(
        _wrap_cmd(args, _cli(args), str(tmp), str(out_json)), 900,
        env=env,
        stdout=open(log_path, "w"),
        stderr=subprocess.STDOUT,
    )
    elapsed = time.monotonic() - t0
    print(f"[detect] amflow_cli exit={r.returncode}  wall={elapsed:.1f}s")

    # Parse [scheme] lines: extract family name + depth indicator.
    # Lines look like: [scheme] building boundary for pmfam_b0_r0_f0_b0_r0_f0_b1_r0_f0
    # or:              [scheme] depth=2: family pmfam_..._b2_r0_f0, scheme=Tradition
    log_text = log_path.read_text(errors="replace")
    scheme_lines = [ln for ln in log_text.splitlines() if "[scheme]" in ln.lower() or
                    "scheme" in ln.lower()]

    # Extract family names per boundary level.
    level_families = {}
    for ln in log_text.splitlines():
        # Match lines referencing pm-prefixed family names with depth info.
        m = re.search(r'depth[=: ]+(\d+)[^a-zA-Z].*?(pm\w+)', ln)
        if m:
            level_families.setdefault(int(m.group(1)), []).append(m.group(2))
        else:
            # Also match any [scheme] line that has a family name.
            m2 = re.search(r'\[scheme\].*?(pm[a-zA-Z0-9_]+)', ln, re.IGNORECASE)
            if m2:
                level_families.setdefault(-1, []).append(m2.group(1))

    print(f"[detect] boundary families per level:")
    levels = sorted(level_families)
    for lv in levels:
        fams = level_families[lv]
        print(f"  level {lv}: {fams}")

    # Check for fixed-point: are successive levels identical modulo b-suffix?
    def strip_b_suffix(name):
        """Strip the trailing _bN_rN_fN component(s) to get the base family name."""
        return re.sub(r'(_b\d+_r\d+_f\d+)+$', '', name)

    verdict = "UNKNOWN"
    fpm = FIXEDPOINT_RE.search(log_text)
    if fpm:
        # Engine-certified: the setup() self-similarity guard compared the
        # boundary family's propagator set against its parent and aborted
        # with the structured code — stronger than the name heuristics.
        print(f"[detect] engine guard GRAVITYFLOW_FIXEDPOINT: boundary "
              f"family {fpm.group(1)} reproduces its parent's propagator "
              f"set at depth {fpm.group(2)} (engine-certified)")
        verdict = "FIXED_POINT"
    elif len(levels) >= 2:
        # Compare consecutive levels
        fixed_point_count = 0
        for i in range(len(levels) - 1):
            fams_a = [strip_b_suffix(f) for f in level_families[levels[i]]]
            fams_b = [strip_b_suffix(f) for f in level_families[levels[i + 1]]]
            if set(fams_a) == set(fams_b):
                fixed_point_count += 1
        if fixed_point_count >= 1:
            verdict = "FIXED_POINT"
        else:
            verdict = "DESCENDING"
    elif len(levels) == 0:
        # Fall back: look for "identical" propagator structure in raw log.
        if "identical" in log_text.lower() or "fixed" in log_text.lower():
            verdict = "FIXED_POINT"
        else:
            verdict = "UNKNOWN (no [scheme] lines found; check log)"

    print(f"[detect] VERDICT: {verdict}")
    if verdict == "FIXED_POINT":
        print("[detect] → AMFlow recursion fixed point detected: vanilla AMFlow cannot produce")
        print("[detect]   boundary data. Proceed: discover → respond → solve → inject.")
    else:
        print("[detect] → Recursion descends; standard AMFlow may suffice.")
    print(f"[detect] wall time: {elapsed:.1f}s")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: discover
# ---------------------------------------------------------------------------

def cmd_discover(args):
    """
    Zero-probe iteration from eb_discover.sh, generalized; dual-surface.

    Run with all-zero explicit_boundary, parse the MISS surface, append,
    repeat (max 15), until solved.  Write final key list to --keys-out.

    Surfaces:
      * "no Vacuum entry" (dual-shape regex) -> master key, kind "vacuum" on
        the classic forced-depth surface, kind "cut" on the depth-0
        leaf-Trivial surface entered after a GRAVITYFLOW_CUTREGION abort.
      * "GRAVITYFLOW_CUTREGION" guard abort -> family marker, kind
        "cutregion"; discover switches to AMFLOW_FORCE_ENDING_DEPTH=0 and
        keeps going (the injection surface; measured behavior pattern).
        The guard only fires on the NORMAL path (forced ending bypasses
        build_boundary), so a cut-quadratic input JSON gets one normal-path
        probe before the classic loop.
    """
    t0 = time.monotonic()
    chain_json = Path(args.json).resolve()
    keys_out = Path(args.keys_out).resolve()

    # Start with empty key list (or existing if file present).
    if keys_out.exists():
        entries = _load_key_entries(keys_out)
        print(f"[discover] resuming with {len(entries)} existing keys from {keys_out}")
    else:
        entries = {}
        print(f"[discover] starting discovery from zero keys")

    etac = _resolve_etac(args, "discover")
    force_depth = args.force_depth
    cut_mode = any(k in ("cut", "cutregion") for k in entries.values())
    if cut_mode and force_depth != "0":
        print("[discover] resumed keys carry cut/cutregion kinds — staying on "
              "the leaf-Trivial injection surface (AMFLOW_FORCE_ENDING_DEPTH=0)")
        force_depth = "0"

    base = chain_json.parent
    wd_disc = base / "wd_disc"
    probe_disc = base / "probe_disc.json"
    out_disc = base / "out_disc.json"
    log_disc = base / "log_disc.log"

    import shutil

    def run_probe(extra_env, log_path):
        # Build the probe JSON with all-zero explicit_boundary over the
        # master keys known so far (family-level cutregion markers excluded).
        probe = json.loads(chain_json.read_text())
        if "amf_options" not in probe:
            probe["amf_options"] = {}
        probe["amf_options"]["explicit_boundary"] = {
            k: {"re": "0", "im": "0"}
            for k in sorted(entries) if entries[k] != "cutregion"
        }
        probe["work_dir"] = str(wd_disc)
        probe_disc.write_text(json.dumps(probe, indent=1))
        if wd_disc.exists():
            shutil.rmtree(str(wd_disc))
        step_t0 = time.monotonic()
        r = _run_clamped(
            _wrap_cmd(args, _cli(args), str(probe_disc), str(out_disc)),
            args.timeout,
            env=_env(args, extra_env),
            stdout=open(log_path, "w"),
            stderr=subprocess.STDOUT,
        )
        return r, time.monotonic() - step_t0

    # ---- GRAVITYFLOW_CUTREGION preflight + normal-path probe --------------
    if not cut_mode and getattr(args, "cut_probe", "auto") != "off":
        cq = _family_cut_quadratic(json.loads(chain_json.read_text()))
        if cq and force_depth == "0":
            cut_mode = True  # user pinned the injection surface directly
            print("[discover] cut-quadratic family with --force-depth 0: "
                  "keys will be tagged kind=cut (injection path)")
        elif cq:
            print("[discover] cut-quadratic CUT propagator in the input JSON — "
                  "running one normal-path probe (no forced ending) for the "
                  "GRAVITYFLOW_CUTREGION guard")
            log_cutprobe = base / "log_disc_cutprobe.log"
            r, wall = run_probe(_etac_env(etac), log_cutprobe)
            print(f"[discover] cut-probe: amflow exit={r.returncode}  "
                  f"wall={wall:.1f}s  log={log_cutprobe}")
            text = log_cutprobe.read_text(errors="replace")
            if ETAC_MISMATCH_MARK in text:
                print("[discover] ERROR: AMFLOW_ETAC_OVERRIDE length mismatch "
                      "(legacy 15-slot default on a non-15-prop family). "
                      "Pass --etac CSV for this family, or --etac '' for the "
                      "engine-chosen etac.")
                return 1
            cm = CUTREGION_RE.search(text)
            if cm:
                fam_key = cm.group(1)
                print(f"[discover] MISS (kind=cutregion): {fam_key} — no "
                      f"eta->inf region family exists; switching to the "
                      f"leaf-Trivial injection surface "
                      f"(AMFLOW_FORCE_ENDING_DEPTH=0)")
                entries.setdefault(fam_key, "cutregion")
                _emit_key_entries(keys_out, entries)
                cut_mode = True
                force_depth = "0"
            else:
                print("[discover] cut-probe did not fire the guard — "
                      "falling back to the classic surface")

    for it in range(1, args.max_iter + 1):
        r, step_elapsed = run_probe(
            {**_etac_env(etac), "AMFLOW_FORCE_ENDING_DEPTH": force_depth},
            log_disc)
        print(f"[discover] iter {it}: amflow exit={r.returncode}  wall={step_elapsed:.1f}s")

        log_text = log_disc.read_text(errors="replace")
        if ETAC_MISMATCH_MARK in log_text:
            print("[discover] ERROR: AMFLOW_ETAC_OVERRIDE length mismatch "
                  "(legacy 15-slot default on a non-15-prop family). "
                  "Pass --etac CSV for this family, or --etac '' for the "
                  "engine-chosen etac.")
            return 1
        cm = CUTREGION_RE.search(log_text)
        if cm:
            fam_key = cm.group(1)
            if fam_key not in entries:
                entries[fam_key] = "cutregion"
                _emit_key_entries(keys_out, entries)
            if not cut_mode:
                print(f"[discover] iter {it}: MISS (kind=cutregion): {fam_key} "
                      f"— switching to AMFLOW_FORCE_ENDING_DEPTH=0 and retrying")
                cut_mode = True
                force_depth = "0"
                continue
            print(f"[discover] ERROR: GRAVITYFLOW_CUTREGION fired on the "
                  f"depth-0 surface for {fam_key} — no further surface to "
                  f"probe; aborting")
            break
        miss = _parse_miss(log_text)
        if not miss:
            print(f"[discover] iter {it}: NO missing key — discovery complete")
            break
        kind = "cut" if cut_mode else "vacuum"
        print(f"[discover] iter {it}: adding key {miss} (kind={kind})")
        entries.setdefault(miss, kind)
        _emit_key_entries(keys_out, entries)
    else:
        print(f"[discover] WARNING: reached max_iter={args.max_iter} without converging")

    _emit_key_entries(keys_out, entries)
    n_vac = sum(1 for v in entries.values() if v == "vacuum")
    n_cut = sum(1 for v in entries.values() if v == "cut")
    n_mark = sum(1 for v in entries.values() if v == "cutregion")
    total = time.monotonic() - t0
    kinds = (f" ({n_vac} vacuum, {n_cut} cut, {n_mark} cutregion markers)"
             if (n_cut or n_mark) else "")
    print(f"[discover] wrote {len(entries)} keys{kinds} → {keys_out}  "
          f"total wall={total:.1f}s")
    # Return contract unchanged from the shipped tool: 0 even when max_iter
    # was reached (callers read the WARNING line / key list, not the rc).
    return 0


# ---------------------------------------------------------------------------
# Subcommand: respond
# ---------------------------------------------------------------------------

def cmd_respond(args):
    """
    Generate per-key unit-injection JSONs and run them N-way parallel with
    AMFLOW_DUMP_EPS_GRID=1 (env block = BASE_ENV). Each run timeout 1800.
    """
    t0 = time.monotonic()
    chain_json = Path(args.json).resolve()
    entries = _load_key_entries(args.keys)
    # Kind routing (dual-surface discover): vacuum keys take the
    # classic unit-injection probe path below; cut keys are injection-path
    # keys (their values come from closed-form providers via `pmflow
    # inject`, see the manual) and are NOT unit-probed; family-level
    # cutregion markers are provenance only. Legacy flat key lists are all
    # kind "vacuum" — behavior unchanged.
    master_keys = [k for k in sorted(entries) if entries[k] != "cutregion"]
    keys = [k for k in master_keys if entries[k] == "vacuum"]
    n_cutk = sum(1 for k in master_keys if entries[k] == "cut")
    n_mark = sum(1 for v in entries.values() if v == "cutregion")
    if n_cutk or n_mark:
        print(f"[respond] key list carries {n_cutk} cut keys and {n_mark} "
              f"cutregion markers: cut keys route to the injection path "
              f"(pmflow inject + closed-form boundary tables) and are "
              f"not unit-probed here; probing the {len(keys)} vacuum keys")
    probe_dir = Path(args.dir).resolve()
    probe_dir.mkdir(parents=True, exist_ok=True)
    n_par = args.par
    timeout = args.timeout
    etac = _resolve_etac(args, "respond")

    base_probe = json.loads(chain_json.read_text())
    if "amf_options" not in base_probe:
        base_probe["amf_options"] = {}

    # Generate one probe JSON per key (unit-injection: that key = 1, all others = 0).
    probe_paths = []
    for k, key in enumerate(keys):
        p = dict(base_probe)
        p["amf_options"] = dict(base_probe.get("amf_options", {}))
        p["amf_options"]["explicit_boundary"] = {
            kk: ({"re": "1", "im": "0"} if kk == key else {"re": "0", "im": "0"})
            for kk in master_keys
        }
        wd = probe_dir / f"wd_{k:02d}"
        wd.mkdir(exist_ok=True)
        p["work_dir"] = str(wd)
        pjson = probe_dir / f"probe_{k:02d}.json"
        pjson.write_text(json.dumps(p, indent=1))
        probe_paths.append((k, pjson))

    print(f"[respond] {len(keys)} probes → {probe_dir}  par={n_par}  timeout={timeout}s")

    env = _env(args, {
        **_etac_env(etac),
        "AMFLOW_FORCE_ENDING_DEPTH": "2",
        "AMFLOW_DUMP_EPS_GRID": "1",
    })

    def run_one(k_pjson):
        k, pjson = k_pjson
        out = probe_dir / f"out_{k:02d}.json"
        log = probe_dir / f"log_{k:02d}.log"
        step_t0 = time.monotonic()
        r = _run_clamped(
            _wrap_cmd(args, _cli(args), str(pjson), str(out)), timeout,
            env=env,
            stdout=open(log, "w"),
            stderr=subprocess.STDOUT,
        )
        elapsed = time.monotonic() - step_t0
        ok = log.exists() and "EPS_GRID" in log.read_text(errors="replace")
        print(f"[respond] probe {k:02d}: exit={r.returncode}  "
              f"EPS_GRID={'YES' if ok else 'NO'}  wall={elapsed:.1f}s")
        return ok

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_par) as ex:
        results = list(ex.map(run_one, probe_paths))

    n_ok = sum(results)
    total = time.monotonic() - t0
    print(f"[respond] {n_ok}/{len(keys)} probes wrote EPS_GRID output  "
          f"total wall={total:.1f}s")
    # Write keys.json next to the logs so gf_solve can find them (kind tags
    # preserved when present; legacy flat lists stay flat).
    kj = probe_dir / "keys.json"
    _emit_key_entries(kj, entries)
    print(f"[respond] keys.json written → {kj}")
    return 0 if n_ok == len(keys) else 1


# ---------------------------------------------------------------------------
# Subcommand: solve
# ---------------------------------------------------------------------------

def cmd_solve(args):
    """
    Thin wrapper delegating to gf_solve.py (imported, same args).
    """
    t0 = time.monotonic()
    # Import gf_solve from the same directory.
    here = Path(__file__).parent
    sys.path.insert(0, str(here))
    import gf_solve
    rc = gf_solve.solve(
        fpA_path=args.fpA,
        probes_dir=args.probes,
        cmap_path=args.cmap,
        out_path=args.out,
        extra_zero=getattr(args, "extra_zero", []),
        dps=getattr(args, "dps", 120),
        n_eps=getattr(args, "n_eps", 12),
        # None = gf_solve's built-in reference-family anchors; a spec
        # string ("SECTOR:FORM,...") is decoded against the fpA basis.
        anchors=getattr(args, "anchors", None),
    )
    elapsed = time.monotonic() - t0
    print(f"[solve] wall time: {elapsed:.1f}s")
    return rc


# ---------------------------------------------------------------------------
# Subcommand: inject
# ---------------------------------------------------------------------------

def cmd_inject(args):
    """
    Build the parent JSON with explicit_boundary and run amflow_cli with
    AMFLOW_FORCE_ENDING_DEPTH.
    """
    t0 = time.monotonic()
    parent_json = Path(args.parent).resolve()
    table_json = Path(args.table).resolve()
    out_json = Path(args.out).resolve()
    depth = args.depth

    parent = json.loads(parent_json.read_text())
    boundary = json.loads(table_json.read_text())

    if "amf_options" not in parent:
        parent["amf_options"] = {}
    parent["amf_options"]["explicit_boundary"] = boundary

    tmp = out_json.parent / (out_json.stem + "_injected.json")
    tmp.write_text(json.dumps(parent, indent=1))

    log_path = out_json.parent / (out_json.stem + "_inject.log")
    env = _env(args, {"AMFLOW_FORCE_ENDING_DEPTH": str(depth)})

    print(f"[inject] parent={parent_json.name}  depth={depth}  → {log_path}")
    r = _run_fenced(
        _wrap_cmd(args, _cli(args), str(tmp), str(out_json)),
        env=env,
        stdout=open(log_path, "w"),
        stderr=subprocess.STDOUT,
    )
    elapsed = time.monotonic() - t0
    print(f"[inject] amflow_cli exit={r.returncode}  wall={elapsed:.1f}s")
    return r.returncode


# ---------------------------------------------------------------------------
# Subcommand: map
# ---------------------------------------------------------------------------

def cmd_map(args):
    """
    Build key -> b0-master-index (circular_map.json) by matching used-prop content
    (normalized strings) to the basis family's prop list.

    Logic (used-prop content matching — the map compares prop CONTENT,
    never just computes norm_props; a content-blind match is 0/68 on the
    reference inputs):
      - Load keys from keys.json; load ending_families.json (family -> {props, masters}).
      - Load fp_A.json for the preferred master basis (parent family); parent
        prop CONTENT from --parent-json (family.propagators) or an
        ending_families entry named like the parent.
      - For each key, extract the family name prefix (everything before the first |).
      - Name matches first: key in the parent preferred basis directly, or in
        sub_masters (mapped through perm_al), keeps working as before.
      - Otherwise CONTENT match: the used-prop signature — multiset of
        (normalized prop, power) over slots with nonzero key index, family-name
        blind — equals a preferred master's signature ⇒ ['master', idx].
        Signature collisions on the preferred side are fail-closed (those keys
        stay unmatched, loudly).
      - If the family is not in ending_families, or no signature matches:
        flag as unmatched.  With --foreign-as-j63, emit ['J63', None].
    """
    t0 = time.monotonic()
    _map_entries = _load_key_entries(args.keys)
    keys = [k for k in sorted(_map_entries) if _map_entries[k] != "cutregion"]
    _n_mark = len(_map_entries) - len(keys)
    if _n_mark:
        print(f"[map] {_n_mark} family-level cutregion markers in the key "
              f"list — provenance only, skipped")
    families = json.loads(Path(args.families).read_text())
    fpa = json.loads(Path(args.fpA).read_text())
    out_path = Path(args.out).resolve()

    # Build preferred basis lookup: key-string -> index.
    preferred = fpa["preferred"]
    parent_fam = fpa["family"]
    # Each preferred entry: "FAMILY|i1|i2|..." — strip the family prefix for matching.
    pref_vecs = {}
    for idx, p in enumerate(preferred):
        parts = p.split("|")
        vec = tuple(int(x) for x in parts[1:])
        fam_prefix = parts[0]
        pref_vecs[(fam_prefix, vec)] = idx

    # sub_masters name-lookup (alpha-flow basis -> preferred positions via
    # perm_al, built the same way as gf_solve).
    sub_masters = fpa.get("sub_masters", [])
    sub_vecs = {}
    for idx2, sm in enumerate(sub_masters):
        sm_parts = sm.split("|")
        sub_vecs[(sm_parts[0], tuple(int(x) for x in sm_parts[1:]))] = idx2
    basis = [tuple(int(x) for x in p.split("|")[1:]) for p in preferred]
    pos_al = {v: i for i, v in enumerate(basis)}
    perm_al = [pos_al[tuple(int(x) for x in sm.split("|")[1:])]
               for sm in sub_masters]

    # Parent prop CONTENT for used-prop matching (fp_A carries names only).
    parent_props = None
    if getattr(args, "parent_json", None):
        parent_props = json.loads(Path(args.parent_json).read_text())[
            "family"]["propagators"]
    elif parent_fam in families:
        parent_props = families[parent_fam].get("props")
    if parent_props is None:
        print(f"[map] ERROR: parent family '{parent_fam}' prop content is in "
              f"neither --parent-json nor ending_families — used-prop content "
              f"matching cannot run. Pass --parent-json CHAIN.json "
              f"(family.propagators).")
        return 2

    def norm_prop(s):
        return re.sub(r'\s+', '', s.lower())

    def used_sig(props, vec):
        """Used-prop content signature: multiset of (normalized prop, power)
        over the slots with nonzero key index — family-name blind."""
        return tuple(sorted((norm_prop(props[s]), vec[s])
                            for s in range(len(vec)) if vec[s] != 0))

    # Preferred-master signatures; collisions are fail-closed (a colliding
    # signature can never be matched — the key stays unmatched, loudly).
    pref_sigs, dup_sigs = {}, set()
    for idx, vec in enumerate(basis):
        sig = used_sig(parent_props, vec)
        if sig in pref_sigs:
            dup_sigs.add(sig)
        else:
            pref_sigs[sig] = idx
    if dup_sigs:
        print(f"[map] WARNING: {len(dup_sigs)} content-ambiguous preferred "
              f"signatures — keys matching one stay unmatched (fail-closed)")

    unmatched = []
    cmap = {}
    for key in keys:
        parts = key.split("|")
        fam_name = parts[0]
        vec = tuple(int(x) for x in parts[1:]) if len(parts) > 1 else ()

        # Try direct match in the parent preferred basis (fam_name matches parent_fam).
        if (fam_name, vec) in pref_vecs:
            cmap[key] = ["master", pref_vecs[(fam_name, vec)]]
            continue

        # Try matching sub_masters (alpha-flow basis, mapped to preferred positions).
        if (fam_name, vec) in sub_vecs:
            cmap[key] = ["master", perm_al[sub_vecs[(fam_name, vec)]]]
            continue

        # Used-prop CONTENT match into the parent b0 basis (ending families).
        fam_info = families.get(fam_name)
        if fam_info is None:
            # Family not in ending_families at all.
            unmatched.append(key)
            continue
        props = fam_info.get("props", [])
        if len(vec) != len(props):
            unmatched.append(key)
            continue
        sig = used_sig(props, vec)
        if sig in pref_sigs and sig not in dup_sigs:
            cmap[key] = ["master", pref_sigs[sig]]
        else:
            # Content-foreign (genuine cut-vacuum leaf) or ambiguous.
            unmatched.append(key)

    if unmatched:
        print(f"[map] {len(unmatched)} unmatched keys (not in preferred basis or ending_families):")
        for uk in unmatched:
            print(f"  {uk}")
        if args.foreign_as_j63:
            print("[map] --foreign-as-j63: defaulting unmatched to ['J63', None]")
            for uk in unmatched:
                cmap[uk] = ["J63", None]
        else:
            print("[map] Re-run with --foreign-as-j63 to default unmatched keys to J63.")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(cmap, indent=0))
    elapsed = time.monotonic() - t0
    print(f"[map] wrote {len(cmap)} entries → {out_path}  wall={elapsed:.1f}s")
    masters_count = sum(1 for v in cmap.values() if v[0] == "master")
    j63_count = sum(1 for v in cmap.values() if v[0] == "J63")
    print(f"[map] master: {masters_count}  J63: {j63_count}  unresolved: {len(unmatched)}")
    return 0 if not (unmatched and not args.foreign_as_j63) else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_common(p):
    p.add_argument("--amflow-cli", default=DEFAULT_CLI, metavar="PATH",
                   help="path to amflow_cli binary")
    p.add_argument("--jemalloc-wrap", default=DEFAULT_JW, metavar="PATH",
                   help="path to jemalloc wrapper script")
    p.add_argument("--cache", default=None, metavar="PATH",
                   help="AMFLOW_IBP_CACHE directory (overrides default)")


def build_parser():
    p = argparse.ArgumentParser(
        prog="pmflow",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # detect
    d = sub.add_parser("detect",
                       help="Run AMFLOW_DEBUG_SCHEME=1 + etac; report fixed-point verdict")
    _add_common(d)
    d.add_argument("--json", required=True, metavar="FAMILY.json",
                   help="AMFlow input JSON for the b0 family")
    d.add_argument("--etac", default=DEFAULT_ETAC, metavar="CSV",
                   help="AMFLOW_ETAC_OVERRIDE comma-separated vector")
    d.add_argument("--work-dir", default=None, metavar="DIR")
    d.set_defaults(func=cmd_detect)

    # discover
    dsc = sub.add_parser("discover",
                         help="Zero-probe key discovery loop (from eb_discover.sh)")
    _add_common(dsc)
    dsc.add_argument("--json", required=True, metavar="CHAIN.json",
                     help="AMFlow input JSON for the chain family")
    dsc.add_argument("--keys-out", required=True, metavar="keys.json",
                     help="output file for discovered key list")
    dsc.add_argument("--max-iter", type=int, default=15, metavar="N",
                     help="maximum discovery iterations (default: 15)")
    dsc.add_argument("--etac", default=None, metavar="CSV",
                     help="AMFLOW_ETAC_OVERRIDE (length must equal #props). "
                          "'' = engine-chosen etac (cut-quadratic families). "
                          "Omitted = DEPRECATED legacy 15-slot default, with "
                          "a warning (backward compatibility only)")
    dsc.add_argument("--force-depth", default="2", metavar="N",
                     help="AMFLOW_FORCE_ENDING_DEPTH for probes (default 2; "
                          "GRAVITYFLOW_CUTREGION families auto-switch to 0)")
    dsc.add_argument("--cut-probe", default="auto", choices=("auto", "off"),
                     help="'auto' (default): cut-quadratic input JSONs get one "
                          "normal-path probe so the GRAVITYFLOW_CUTREGION "
                          "guard can fire; 'off' restores the classic loop "
                          "exactly")
    dsc.add_argument("--timeout", type=int, default=1200, metavar="SEC",
                     help="per-iteration amflow_cli timeout (default: 1200)")
    dsc.set_defaults(func=cmd_discover)

    # respond
    rsp = sub.add_parser("respond",
                         help="Generate + run unit-injection probes N-way parallel")
    _add_common(rsp)
    rsp.add_argument("--json", required=True, metavar="CHAIN.json",
                     help="AMFlow input JSON for the chain family")
    rsp.add_argument("--keys", required=True, metavar="keys.json",
                     help="key list (from discover)")
    rsp.add_argument("--dir", required=True, metavar="PROBEDIR",
                     help="directory for probe JSONs and log files")
    rsp.add_argument("--par", type=int, required=True, metavar="N",
                     help="number of parallel amflow_cli processes")
    rsp.add_argument("--etac", default=None, metavar="CSV",
                     help="AMFLOW_ETAC_OVERRIDE (length must equal #props). "
                          "'' = engine-chosen etac. Omitted = DEPRECATED legacy "
                          "15-slot default, with a warning")
    rsp.add_argument("--timeout", type=int, default=1800, metavar="SEC",
                     help="per-probe timeout seconds (default: 1800)")
    rsp.set_defaults(func=cmd_respond)

    # solve
    sv = sub.add_parser("solve",
                        help="Assemble fixed-point system → explicit_boundary JSON (via gf_solve.py)")
    _add_common(sv)
    sv.add_argument("--fpA", required=True, metavar="fp_A.json",
                    help="second-flow response (AMFLOW_FIXEDPOINT_PROBE)")
    sv.add_argument("--probes", required=True, metavar="DIR",
                    help="directory with log_XX.log files + keys.json")
    sv.add_argument("--cmap", required=True, metavar="circular_map.json",
                    help="key -> ['master', idx] | ['J63', None]")
    sv.add_argument("--out", required=True, metavar="fp_boundary.json",
                    help="output: per-eps explicit_boundary table")
    sv.add_argument("--extra-zero", nargs="*", default=[], metavar="KEY",
                    help="additional keys emitted as exact zeros")
    sv.add_argument("--anchors", default=None, metavar="SPEC",
                    help="comma-separated SECTOR:FORM analytic anchors, e.g. "
                         "'63:J63,111:J63' (the default). SECTOR is the "
                         "sector index of a corner master in the fpA "
                         "preferred basis (bit i of SECTOR = index slot i); "
                         "FORM names the analytic value (built-in: J63, the "
                         "cut-tadpole closed form). REQUIRED in practice on "
                         "any family that is not the 15-slot reference "
                         "family — the default anchors only fit that one")
    sv.add_argument("--n-eps", type=int, default=12, metavar="N",
                    help="eps-grid length of the probe logs (default: 12)")
    sv.add_argument("--dps", type=int, default=120, metavar="N",
                    help="mpmath precision (default: 120)")
    sv.set_defaults(func=cmd_solve)

    # inject
    inj = sub.add_parser("inject",
                         help="Build parent JSON + run amflow_cli with FORCE_ENDING_DEPTH")
    _add_common(inj)
    inj.add_argument("--parent", required=True, metavar="PARENT.json",
                     help="AMFlow input JSON for the parent (cut) family")
    inj.add_argument("--table", required=True, metavar="fp_boundary.json",
                     help="explicit_boundary table (output of solve)")
    inj.add_argument("--depth", type=int, default=1, metavar="N",
                     help="AMFLOW_FORCE_ENDING_DEPTH (default: 1)")
    inj.add_argument("--out", required=True, metavar="OUT.json",
                     help="destination for the AMFlow solve_integrals output")
    inj.set_defaults(func=cmd_inject)

    # map
    mp_ = sub.add_parser("map",
                         help="Build key->b0-index circular_map.json")
    _add_common(mp_)
    mp_.add_argument("--keys", required=True, metavar="keys.json",
                     help="key list (from discover or respond)")
    mp_.add_argument("--families", required=True, metavar="ending_families.json",
                     help="ending family prop/master metadata")
    mp_.add_argument("--fpA", required=True, metavar="fp_A.json",
                     help="second-flow response (preferred basis)")
    mp_.add_argument("--parent-json", default=None, metavar="CHAIN.json",
                     help="AMFlow input JSON carrying the parent family's "
                          "propagators (family.propagators) for used-prop "
                          "content matching; omitted = look the parent family "
                          "up in ending_families.json")
    mp_.add_argument("--out", required=True, metavar="circular_map.json",
                     help="output: key -> ['master', idx] or ['J63', None]")
    mp_.add_argument("--foreign-as-j63", action="store_true",
                     help="default unmatched keys to ['J63', None]")
    mp_.set_defaults(func=cmd_map)

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()

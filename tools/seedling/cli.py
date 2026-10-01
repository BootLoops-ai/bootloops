"""seedling CLI.

Subcommands:
  pin       parse targets -> support + chain floor -> emit pinned jobs.yaml
            (default DRY: writes jobs.yaml + prints the cell; --execute stages
            in a fresh cgroup via runner)
  run       the full pipeline:
            config+targets -> support (with optional sector-tree closure
            under family symmetry maps via --sector-maps) -> census ->
            MARGINS (registered M2
            model via margins.py, or --margins=const baseline arm) ->
            per-sector schedule (pinner.build_schedule) -> jobs.yaml +
            run-dir assembly -> preflight (bulkmodel) -> [--execute: staged
            kira in a fresh cgroup -> SEL/masters extraction -> two-slice
            gate -> de-census gate -> receipt emission (injected hook,
            retrofit mode)] -> ledger JSONL. Without --execute: DRY-RUN
            prints the full plan + preflight memory prediction.
            Slice row files and the receipt emitter are INJECTED (JSON
            paths / hook module) — no caller-specific paths inside this
            package; gate steps without inputs are ledgered as SKIPPED with
            the reason, never silently omitted.
  preflight predict staged bulk + RSS from the seed-box model (census file)
  certify   two-slice differential gate on JSON row dicts (zero-normalized)
  census    DE-census gate: masters + eta-props + gated rows -> uncovered subset
            + mini-kira spec emission
  ledger    print a run's JSONL telemetry ledger
  audit     family preflight: full autopermutation group of an integral
            family (external-leg crossings ∘ C ∘ affine unimodular loop maps
            that permute the denominator set; anchor-solved search, every
            generator re-verified by explicit sympy substitution). Emits
            <fam>_AUT.json + .kira seeds + magic_relations drop-in.
  identity  same-family check for paper-transcribed families via the affine
            scalar-product-dictionary / second-Symanzik-F route (audit.py
            member owns the parser; no-arg form = the built-in wpair
            regression fixture). rc 0 PASS / 1 FAIL-loud, mismatch named.

Safety defaults: nothing executes kira unless --execute is passed AND a cgroup
spec is given; margins obey the hard floors (meff >= 1; s >= chain floor —
chain-floor rule) unless --unsafe-s0 is passed explicitly (logged).
"""
import argparse
import json
import os
import re
import shutil
import sys
import time

from . import support, pinner, bulkmodel, runner, certify


def _family_from_config(config_dir):
    """Read family name + top sectors from a kira config dir."""
    import yaml
    with open(os.path.join(config_dir, "integralfamilies.yaml")) as fh:
        doc = yaml.safe_load(fh)
    fam = doc["integralfamilies"][0]
    return fam["name"], list(fam.get("top_level_sectors", [])), len(fam["propagators"])


def _parse_margins(text):
    vals = {"r": 0, "s": 0, "d": 0}
    for part in text.split(","):
        k, v = part.split(":")
        vals[k.strip()] = int(v)
    return (vals["r"], vals["s"], vals["d"])


_MEM_MAX_RE = re.compile(r"max|[0-9]+[KkMmGgTt]?")
_CPU_MAX_RE = re.compile(r"(max|[0-9]+)( +[0-9]+)?")


def _cgroup_spec(text):
    """Parse `--cgroup NAME,MEM_MAX,CPU_MAX` into runner's cgroup_cfg.

    Returns (cfg, None) or (None, reason). Refused — with NOTHING built —
    are: an empty or illegal NAME (runner.validate_cgroup_name: non-empty,
    ^[A-Za-z0-9_.-]+$, never a dot-only component such as `..`), an empty
    MEM_MAX or CPU_MAX, and a value that is not a cgroup v2 memory.max /
    cpu.max literal (the values are interpolated into a shell script that
    runs under sudo). The cgroup.procs write discipline is documented in
    runner.build_cgroup_script. An empty NAME used to resolve to the cgroup
    ROOT and was refused only by accident (the root has no memory.max)."""
    parts = text.split(",", 2)
    if len(parts) != 3:
        return None, f"--cgroup {text!r}: need NAME,MEM_MAX,CPU_MAX"
    name, mem, cpu = (p.strip() for p in parts)
    try:
        runner.validate_cgroup_name(name)
    except ValueError as e:
        return None, f"--cgroup {text!r}: {e}"
    if not mem:
        return None, f"--cgroup {text!r}: MEM_MAX is empty"
    if not cpu:
        return None, f"--cgroup {text!r}: CPU_MAX is empty"
    if not _MEM_MAX_RE.fullmatch(mem):
        return None, (f"--cgroup {text!r}: MEM_MAX {mem!r} is not a memory.max "
                      f"literal (max or N[KMGT])")
    if not _CPU_MAX_RE.fullmatch(cpu):
        return None, (f"--cgroup {text!r}: CPU_MAX {cpu!r} is not a cpu.max "
                      f"literal (max | QUOTA [PERIOD])")
    return {"name": name, "mem_max": mem, "cpu_max": cpu}, None


def _refuse_bad_cgroup(a):
    """rc 2 + one stderr line when a given --cgroup spec is illegal; None
    when it is absent or valid. Runs BEFORE anything is built."""
    if not a.cgroup:
        return None
    _, why = _cgroup_spec(a.cgroup)
    if why:
        print(f"refusing: {why} (nothing built)", file=sys.stderr)
        return 2
    return None


def cmd_pin(a):
    rc = _refuse_bad_cgroup(a)
    if rc:
        return rc
    name, tops, _ = _family_from_config(a.config)
    targets = support.parse_targets(a.targets)
    g = support.global_support(targets)
    s_floor = 0 if a.unsafe_s0 else None
    r, s, d = pinner.pinned_cell(g, targets, tops, _parse_margins(a.margin), s_floor)
    os.makedirs(a.out, exist_ok=True)
    text = pinner.render_jobs_yaml(name, tops, r, s, d,
                                   preferred_name=a.preferred)
    with open(os.path.join(a.out, "jobs.yaml"), "w") as fh:
        fh.write(text)
    print(f"pin: family={name} tops={tops} support=({g['r']},{g['s']},{g['d']}) "
          f"cell=({r},{s},{d}){' UNSAFE-S0' if a.unsafe_s0 else ''} -> {a.out}/jobs.yaml")
    if a.execute:
        if not a.cgroup:
            print("refusing to execute without --cgroup NAME,MEM,CPUMAX", file=sys.stderr)
            return 2
        cg_cfg, why = _cgroup_spec(a.cgroup)
        if why:                      # unreachable after the early refusal; fail closed
            raise ValueError(why)
        rec, _ = runner.run_phase("stage", a.out, a.kira_cmd, cg_cfg,
                                  os.path.join(a.out, "ledger.jsonl"))
        print(json.dumps(rec))
        return 0 if rec["exit_status"] == 0 else 1
    return 0


_TGT_KEY = re.compile(r"^\s*[A-Za-z_]\w*\[([-\d,\s]*)\]")


_GATE_PHASES = ("two_slice_gate", "de_census", "receipts")
_EVIDENCE_KEYS = ("evidence", "evidence_path", "witness_dir", "receipt_dir",
                  "paths", "minikira_dir")


def _check_gate_bookkeeping(rec):
    """Gate-bookkeeping assertion: every GATE emission must name its
    receipt/evidence file(s); a SKIPPED gate must name its reason.
    Fail-closed: a gate record we cannot point back to evidence is a bug."""
    if rec.get("phase") not in _GATE_PHASES:
        return
    if rec.get("verdict") == "SKIPPED":
        if not rec.get("reason"):
            raise ValueError(f"gate bookkeeping: SKIPPED {rec['phase']} "
                             f"record must carry a 'reason'")
        return
    if not any(rec.get(k) for k in _EVIDENCE_KEYS):
        raise ValueError(
            f"gate bookkeeping: {rec['phase']} record names no evidence "
            f"file (need one of {_EVIDENCE_KEYS}); record keys: "
            f"{sorted(rec)}")


def _ledger_write(path, rec):
    rec = dict(rec)
    rec.setdefault("ts", time.time())
    _check_gate_bookkeeping(rec)
    with open(path, "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return rec


def _n_ops_from_config(config_dir, n_loops):
    """n_ops = L*(L+E); E = #independent external momenta (incoming - 1)."""
    import yaml
    try:
        with open(os.path.join(config_dir, "kinematics.yaml")) as fh:
            kin = yaml.safe_load(fh)
        n_ext = len(kin["kinematics"]["incoming_momenta"]) - 1
    except (OSError, KeyError, TypeError):
        n_ext = 3          # recorded in the ledger either way
    return n_loops * (n_loops + n_ext), n_ext


def _load_census(nts_path, tops):
    cap = max(tops) if tops else None
    census = {}
    for line in open(nts_path):
        parts = line.split()
        if len(parts) == 2 and (cap is None or int(parts[0]) <= cap):
            census[int(parts[0])] = int(parts[1])
    return census


def _parse_row_keys(keys):
    """'fam[i1,...]' strings -> index tuples (for de-census / escalation)."""
    out = []
    for k in keys:
        m = _TGT_KEY.match(k)
        if m and m.group(1).strip():
            out.append(tuple(int(x) for x in m.group(1).split(",")))
    return out


def cmd_run(a):
    """FAIL-CLOSED wrapper (known pitfall: a loud print with rc 0):
    ANY exception in ANY stage -> ledger phase=failed (stage name +
    traceback) + rc 1. A gate that cannot run is a FAILURE, never a skip."""
    rc = _refuse_bad_cgroup(a)       # illegal --cgroup: rc 2 before anything is built
    if rc:
        return rc
    os.makedirs(a.out, exist_ok=True)
    ledger = os.path.join(a.out, "ledger.jsonl")
    stage = {"name": "init"}
    try:
        return _run_pipeline(a, ledger, stage)
    except Exception as e:
        import traceback
        _ledger_write(ledger, {"phase": "failed", "stage": stage["name"],
                               "error": repr(e),
                               "traceback": traceback.format_exc()[-2000:]})
        print(f"seedling run FAILED-CLOSED at stage {stage['name']!r}: "
              f"{e!r} (ledger phase=failed)", file=sys.stderr)
        return 1


_LAPPER_HOME = os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))


def _ensure_dep(mod_name, subdir):
    """Locate a Lapper-family dependency robustly; clear error naming it."""
    try:
        __import__(mod_name)
        return
    except ImportError:
        pass
    cand = os.path.join(_LAPPER_HOME, subdir)
    if os.path.isdir(cand) and cand not in sys.path:
        sys.path.insert(0, cand)
    try:
        __import__(mod_name)
    except ImportError as e:
        raise RuntimeError(
            f"seedling run: missing Lapper-family dependency {mod_name!r} "
            f"(not importable, and no package under {cand})") from e


def _load_hook(path):
    """Load an injected hook module; hooks may use ibplapper (dep-guarded
    HERE so a broken env fails loudly with the dep named, not with a bare
    ModuleNotFoundError inside the hook)."""
    _ensure_dep("ibplapper", "ibplapper")
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_seedling_hook_" + os.path.basename(path).replace(".", "_"), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_pipeline(a, ledger, stage):
    from . import margins as margins_mod, rankcheck, escalate
    stage["name"] = "config"
    name, tops, n_props = _family_from_config(a.config)
    import yaml
    with open(os.path.join(a.config, "integralfamilies.yaml")) as fh:
        n_loops = len(yaml.safe_load(fh)["integralfamilies"][0]["loop_momenta"])
    targets = support.parse_targets(a.targets)
    if getattr(a, "sector_maps", None):
        # sector-tree closure: fold every symmetry image of the targets into
        # the support before pinning, so mapped-to sectors get target-grade
        # cells (undecidable images are ledgered, never silently dropped)
        stage["name"] = "support_closure"
        smaps = rankcheck.parse_sector_maps(a.sector_maps, n_props)
        targets, srep = support.symmetry_closure(targets, smaps)
        _ledger_write(ledger, {"phase": "support_closure",
                               "maps_path": a.sector_maps, **srep})
        print(f"support closure: {srep['n_input']} targets + "
              f"{srep['n_added']} symmetry images "
              f"({len(srep['undecidable'])} undecidable)")
        stage["name"] = "config"
    tab = support.support_table(targets)
    g = support.global_support(targets)
    nts = a.nts or os.path.join(os.path.dirname(os.path.abspath(a.config)),
                                "sectormappings", name, "nonTrivialSector")
    census = _load_census(nts, tops)

    # --- margins (registered model or const arm; hard floors) --------------
    stage["name"] = "margins"
    mrep = margins_mod.predict_margins(
        census, tab, targets, tops, mode=a.margins, n_sp=n_props,
        labels_path=(a.labels or margins_mod.DEFAULT_LABELS))
    # schedule s-floor: every cell s >= 1 (chain-floor rule; target cells
    # additionally carry their own target-s via build_schedule's
    # max(e["s"], s_floor) — exactly the validated staging recipe).
    # chain_s_floor (global-pin semantics) is recorded for provenance;
    # s=0 only via explicit unsafe flag.
    s_used = 0 if a.unsafe_s0 else 1
    _ledger_write(ledger, {"phase": "margins", "mode": mrep["mode"],
                           "per_class": mrep["per_class"],
                           "s_floor": mrep["s_floor"], "s_used": s_used,
                           "unsafe_s0": bool(a.unsafe_s0),
                           "provenance": mrep["provenance"]})
    if a.unsafe_s0:
        print("WARNING: --unsafe-s0 (UNSAFE: drops the chain floor) — logged", file=sys.stderr)

    # --- per-sector schedule (pinner) ---------------------------------------
    stage["name"] = "schedule"
    mm = {t: v["margin"] for t, v in mrep["per_class"].items()}

    def closure_cell(t):
        m = mm.get(t, mrep["b_const"])   # unseen class: B-const (conservative)
        return (t + m, s_used, m)

    groups = pinner.build_schedule(tab, census, closure_cell, s_floor=s_used)
    ytxt = pinner.render_schedule_jobs_yaml(
        name, groups, preferred_name=("preferred" if a.preferred else None))
    with open(os.path.join(a.out, "jobs.yaml"), "w") as fh:
        fh.write(ytxt)
    shutil.copy(a.targets, os.path.join(a.out, "target"))
    if a.preferred:
        shutil.copy(a.preferred, os.path.join(a.out, "preferred"))
    dst_cfg = os.path.join(a.out, "config")
    if not os.path.isdir(dst_cfg):
        shutil.copytree(a.config, dst_cfg)
    _ledger_write(ledger, {"phase": "schedule",
                           "groups": {f"{c}": secs
                                      for c, secs in sorted(groups.items())},
                           "n_sectors": sum(len(v) for v in groups.values())})

    # --- preflight (bulkmodel; estimate, ~2x band) ---------------------------
    stage["name"] = "preflight"
    n_ops, n_ext = _n_ops_from_config(a.config, n_loops)
    B = 0
    for (r, s, d), secs in groups.items():
        for sec in secs:
            B += n_ops * bulkmodel.seeds_per_sector(
                bin(sec).count("1"), n_props, r, s, d)
    pf = {"bulk_eqns": B, "n_ops": n_ops, "L": n_loops, "E": n_ext,
          "select_floor_GB": round(bulkmodel.C_SEL * B, 3),
          "generate_peak_GB": round(bulkmodel.A_GEN + bulkmodel.C_GEN * B, 3)}
    _ledger_write(ledger, {"phase": "preflight", **pf})
    print(f"run: family={name} tops={tops} support=({g['r']},{g['s']},{g['d']}) "
          f"margins={a.margins} cells={sorted(groups)} s_used={s_used}")
    print(f"preflight: bulk={B:,} eqns, gen-peak ~{pf['generate_peak_GB']} GB "
          f"(2x band) -> {a.out}/jobs.yaml")
    if not a.execute:
        print("DRY-RUN (no --execute): plan + preflight emitted; "
              "nothing launched")
        return 0

    # --- staged kira (fresh cgroup; runner) ----------------------------------
    stage["name"] = "stage"
    if not a.cgroup:
        print("refusing to execute without --cgroup NAME,MEM,CPUMAX",
              file=sys.stderr)
        return 2
    cg_cfg, why = _cgroup_spec(a.cgroup)
    if why:                          # unreachable after cmd_run's early refusal; fail closed
        raise ValueError(why)
    rec, _ = runner.run_phase("stage", a.out, a.kira_cmd, cg_cfg, ledger)
    esc = escalate.EscalationState()
    rc = 0
    if rec["exit_status"] != 0:
        action, newm = escalate.next_margins(esc, escalate.LEFTOVER_TARGETS)
        _ledger_write(ledger, {"phase": "escalation",
                               "signal": escalate.LEFTOVER_TARGETS,
                               "action": action, "new_margins": list(newm)})
        print(f"stage FAILED (rc={rec['exit_status']}); escalation: "
              f"{action} -> margins {newm}", file=sys.stderr)
        return 1

    # --- SEL/masters extraction ----------------------------------------------
    stage["name"] = "extract"
    klog = os.path.join(a.out, "kira.log")
    masters = rankcheck.parse_kira_masters(klog) if os.path.exists(klog) else []
    sel = None
    if os.path.exists(klog):
        sels = [int(m.group(1)) for m in
                re.finditer(r"Number of selected equations to reduce:\s*(\d+)",
                            open(klog).read())]
        sel = max(sels) if sels else None
    with open(os.path.join(a.out, "masters.json"), "w") as fh:
        json.dump([list(m[0]) for m in masters], fh)
    _ledger_write(ledger, {"phase": "extract", "n_masters": len(masters),
                           "selected_eqns": sel})
    print(f"extract: masters={len(masters)} selected_eqns={sel}")

    # --- two-slice gate (mandatory; inputs injected) --------------------------
    stage["name"] = "slice_hook"
    if a.slice_hook and not (a.fresh and a.ref):
        hook = _load_hook(a.slice_hook)
        paths = hook.make(a.out)      # (fresh, ref, fresh2, ref2) JSON paths
        a.fresh, a.ref, a.fresh2, a.ref2 = paths
        _ledger_write(ledger, {"phase": "slice_hook", "paths": list(paths)})
    stage["name"] = "two_slice_gate"
    if a.fresh and a.ref:
        f1, r1 = certify.load_rows_json(a.fresh), certify.load_rows_json(a.ref)
        f2 = certify.load_rows_json(a.fresh2) if a.fresh2 else None
        r2 = certify.load_rows_json(a.ref2) if a.ref2 else None
        v = certify.two_slice_certify(f1, r1, f2, r2)
        _ledger_write(ledger, {"phase": "two_slice_gate", **v,
                               "evidence": {"fresh": a.fresh, "ref": a.ref,
                                            "fresh2": a.fresh2,
                                            "ref2": a.ref2}})
        print(f"two-slice gate: {v['verdict']}")
        if v["verdict"] != runner.GATE_PASS:
            rc = 1                 # ANY non-PASS = uncertified (fail-closed)
        if v["verdict"] == runner.GATE_FAIL_PIN:
            sig = escalate.classify_gate_failure(
                _parse_row_keys(v["fails1"]))
            action, newm = escalate.next_margins(esc, sig)
            _ledger_write(ledger, {"phase": "escalation", "signal": sig,
                                   "action": action,
                                   "new_margins": list(newm)})
            print(f"GATE FAIL_PIN; escalation: {sig} -> {action} margins "
                  f"{newm}", file=sys.stderr)
            rc = 1
        gated_keys = [] if v["verdict"] != runner.GATE_PASS else list(f1)
    else:
        _ledger_write(ledger, {"phase": "two_slice_gate",
                               "verdict": "SKIPPED",
                               "reason": "no --fresh/--ref slice JSONs "
                                         "supplied (gate remains MANDATORY "
                                         "before consuming the table)"})
        print("two-slice gate: SKIPPED (no slice JSONs) — table NOT certified")
        gated_keys = []
        rc = max(rc, 3)

    # --- de-census gate --------------------------------------------------------
    stage["name"] = "de_census"
    if a.eta_props:
        eta = [int(x) for x in a.eta_props.split(",")]
        drep = runner.de_census_gate([m[0] for m in masters], eta,
                                     _parse_row_keys(gated_keys))
        _ledger_write(ledger, {"phase": "de_census",
                               "census_n": drep["census_n"],
                               "gated_n": drep["gated_n"],
                               "master_self_n": drep["master_self_n"],
                               "n_uncovered": len(drep["uncovered"]),
                               "evidence": {
                                   "masters": os.path.join(a.out,
                                                           "masters.json"),
                                   "ledger": ledger}})
        print(f"de-census: {drep['census_n']} rows, "
              f"{len(drep['uncovered'])} uncovered")
        if drep["uncovered"]:
            spec = runner.minikira_spec(name, tops, drep["uncovered"])
            mk = os.path.join(a.out, "minikira")
            os.makedirs(mk, exist_ok=True)
            for fn, key in (("target", "target"), ("jobs.yaml", "jobs_yaml")):
                with open(os.path.join(mk, fn), "w") as fh:
                    fh.write(spec[key])
            action, _ = escalate.next_margins(esc,
                                              escalate.DE_CENSUS_UNCOVERED)
            _ledger_write(ledger, {"phase": "escalation",
                                   "signal": escalate.DE_CENSUS_UNCOVERED,
                                   "action": action,
                                   "minikira_cell": list(spec["cell"])})
            print(f"de-census uncovered -> minikira spec (cell "
                  f"{spec['cell']}) -> {mk}/ [{action}]")
    else:
        _ledger_write(ledger, {"phase": "de_census", "verdict": "SKIPPED",
                               "reason": "no --eta-props supplied"})

    # --- receipt emission (retrofit mode; injected hook) -----------------------
    stage["name"] = "receipts"
    if a.receipt_hook:
        hook = _load_hook(a.receipt_hook)
        rrec = hook.emit(a.out)
        _ledger_write(ledger, {"phase": "receipts", **rrec})
        print(f"receipts: {rrec}")
        if not rrec.get("ok", False):
            rc = max(rc, 1)
    else:
        _ledger_write(ledger, {"phase": "receipts", "verdict": "SKIPPED",
                               "reason": "no --receipt-hook supplied"})
    _ledger_write(ledger, {"phase": "done", "rc": rc})
    return rc


def cmd_preflight(a):
    census = []
    for line in open(a.census):
        parts = line.split()
        if len(parts) == 2 and int(parts[0]) <= a.top:
            census.append(int(parts[1]))
    r, s, d = (int(x) for x in a.cell.split(","))
    pf = bulkmodel.preflight(census, a.nsp, a.nops, r, s, d)
    print(json.dumps({"sectors": len(census), "cell": [r, s, d], **pf}, indent=1))
    return 0


def cmd_certify(a):
    fresh1 = certify.load_rows_json(a.fresh)
    ref1 = certify.load_rows_json(a.ref)
    fresh2 = certify.load_rows_json(a.fresh2) if a.fresh2 else None
    ref2 = certify.load_rows_json(a.ref2) if a.ref2 else None
    v = certify.two_slice_certify(fresh1, ref1, fresh2, ref2)
    print(json.dumps(v))
    return 0 if v["verdict"] == runner.GATE_PASS else 1


def cmd_census(a):
    masters = [tuple(m) for m in json.load(open(a.masters))]
    gated = [tuple(g) for g in json.load(open(a.gated))] if a.gated else []
    eta_props = [int(x) for x in a.eta_props.split(",")]
    rep = runner.de_census_gate(masters, eta_props, gated)
    print(json.dumps({k: v if k != "uncovered" else [list(u) for u in v]
                      for k, v in rep.items()}))
    if rep["uncovered"] and a.emit_minikira:
        name, tops, _ = _family_from_config(a.config)
        spec = runner.minikira_spec(name, tops, rep["uncovered"])
        os.makedirs(a.emit_minikira, exist_ok=True)
        for fn, key in (("target", "target"), ("jobs.yaml", "jobs_yaml")):
            with open(os.path.join(a.emit_minikira, fn), "w") as fh:
                fh.write(spec[key])
        print(f"minikira spec (cell {spec['cell']}) -> {a.emit_minikira}/")
    return 0 if not rep["uncovered"] or a.emit_minikira else 1


def cmd_ledger(a):
    path = os.path.join(a.rundir, "ledger.jsonl")
    if not os.path.exists(path):
        print(f"no ledger at {path}", file=sys.stderr)
        return 1
    for line in open(path):
        print(line.rstrip())
    return 0


def cmd_audit(a):
    """Aut-mode family preflight (the audit.py member)."""
    from . import audit as _aud
    _aud.run(a.config, out_dir=a.out, max_shift=a.max_shift,
             use_cache=not a.no_cache, report_crossings=a.report_crossings,
             n_denom=a.n_denom)
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # `identity` passthrough: the audit member owns this parser — no-arg =
    # the built-in wpair regression fixture, flags = the GENERAL
    # affine-dictionary mode (its argparse handles -h).
    if argv[:1] == ["identity"]:
        from . import audit as _aud
        if len(argv) == 1:
            _aud.identity_check_wpair()
            return 0
        return _aud.identity_check_general(argv[1:])
    p = argparse.ArgumentParser(prog="seedling",
        description="Adaptive certified seeding wrapper for Kira IBP reduction")
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("pin", help="emit pinned jobs.yaml from target support")
    sp.add_argument("--config", required=True, help="kira config dir")
    sp.add_argument("--targets", required=True)
    sp.add_argument("--out", required=True)
    sp.add_argument("--margin", default="r:0,s:0,d:0")
    sp.add_argument("--preferred", default=None)
    sp.add_argument("--unsafe-s0", action="store_true")
    sp.add_argument("--execute", action="store_true")
    sp.add_argument("--cgroup", default=None, help="NAME,MEM_MAX,CPU_MAX")
    sp.add_argument("--kira-cmd", default="kira --parallel=4 jobs.yaml")
    sp.set_defaults(fn=cmd_pin)

    sp = sub.add_parser("run", help="full pipeline: margins -> schedule -> "
                                    "stage -> gates -> receipts -> ledger")
    sp.add_argument("--config", required=True, help="kira config dir")
    sp.add_argument("--targets", required=True)
    sp.add_argument("--out", required=True)
    sp.add_argument("--margins", default="m2", choices=("m2", "const"),
                    help="m2 = registered model; const = B-const baseline arm")
    sp.add_argument("--labels", default=None)
    sp.add_argument("--nts", default=None,
                    help="nonTrivialSector file (default: sibling "
                         "sectormappings/<fam>/nonTrivialSector)")
    sp.add_argument("--sector-maps", default=None,
                    help="kira sectorSymmetries/sectorRelations file: close "
                         "the target support under the family symmetry maps "
                         "before pinning (support.symmetry_closure; images "
                         "land in the mapped-to sectors)")
    sp.add_argument("--preferred", default=None,
                    help="preferred-masters file to copy into the run dir")
    sp.add_argument("--unsafe-s0", action="store_true")
    sp.add_argument("--execute", action="store_true")
    sp.add_argument("--cgroup", default=None, help="NAME,MEM_MAX,CPU_MAX")
    sp.add_argument("--kira-cmd", default="kira --parallel=4 jobs.yaml")
    sp.add_argument("--fresh", default=None, help="slice-1 rows JSON")
    sp.add_argument("--ref", default=None, help="slice-1 reference rows JSON")
    sp.add_argument("--fresh2", default=None)
    sp.add_argument("--ref2", default=None)
    sp.add_argument("--eta-props", default=None,
                    help="comma 0-based eta-carrying propagator positions")
    sp.add_argument("--receipt-hook", default=None,
                    help="python file with emit(rundir)->dict (retrofit "
                         "receipt emission; injected, caller-side)")
    sp.add_argument("--slice-hook", default=None,
                    help="python file with make(rundir)->(fresh,ref,fresh2,"
                         "ref2) JSON paths (two-slice inputs from the "
                         "staged output; injected, caller-side)")
    sp.set_defaults(fn=cmd_run)

    sp = sub.add_parser("preflight", help="predict bulk + RSS for a cell")
    sp.add_argument("--census", required=True, help="kira nonTrivialSector file")
    sp.add_argument("--top", type=int, required=True, help="top sector id filter")
    sp.add_argument("--nsp", type=int, required=True)
    sp.add_argument("--nops", type=int, required=True)
    sp.add_argument("--cell", required=True, help="r,s,d")
    sp.set_defaults(fn=cmd_preflight)

    sp = sub.add_parser("certify", help="two-slice gate on JSON row dicts")
    sp.add_argument("--fresh", required=True)
    sp.add_argument("--ref", required=True)
    sp.add_argument("--fresh2", default=None)
    sp.add_argument("--ref2", default=None)
    sp.set_defaults(fn=cmd_certify)

    sp = sub.add_parser("census", help="DE-census gate + minikira emission")
    sp.add_argument("--masters", required=True, help="JSON list of index lists")
    sp.add_argument("--eta-props", required=True, help="comma 0-based positions")
    sp.add_argument("--gated", default=None, help="JSON list of gated rows")
    sp.add_argument("--config", default=None)
    sp.add_argument("--emit-minikira", default=None, help="output dir")
    sp.set_defaults(fn=cmd_census)

    sp = sub.add_parser("ledger", help="dump run telemetry")
    sp.add_argument("--rundir", required=True)
    sp.set_defaults(fn=cmd_ledger)

    sp = sub.add_parser("audit", help="family autopermutation preflight "
                                      "(crossings ∘ C ∘ loop maps; Family "
                                      "Audit fold)")
    sp.add_argument("config", help="kira config dir OR amflow input JSON")
    sp.add_argument("--out", default=None,
                    help="output dir (default: alongside config)")
    sp.add_argument("--max-shift", type=int, default=1,
                    help="coefficient bound for |M_ij|,|A_ij| (default 1): a "
                         "valid hit above it is STILL recorded, the receipt "
                         "sets search_cap_hit (flags, does not filter)")
    sp.add_argument("--n-denom", type=int, default=None,
                    help="override denominator/ISP split")
    sp.add_argument("--no-cache", action="store_true")
    sp.add_argument("--report-crossings", action="store_true",
                    help="also report leg perms that PERMUTE invariants")
    sp.set_defaults(fn=cmd_audit)

    # `identity` is dispatched BEFORE this parser (the audit member owns its
    # argparse); this stub exists so it shows in --help.
    sub.add_parser("identity", add_help=False,
                   help="same-family check (affine SP-dictionary / 2nd-"
                        "Symanzik F); no-arg = built-in wpair fixture")

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())

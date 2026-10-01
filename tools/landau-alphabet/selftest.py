#!/usr/bin/env python3
r"""selftest.py -- validation battery, two tiers.

Default (no flag): FAST code-check battery.  Real checks with known answers on
the package's own code paths, a few seconds of wall on a laptop:

  * full face-engine run on the 1-loop massless box (9 faces): alphabet must be
    exactly {s, t}, every face resolved, complete flag true, honesty tallies
    consistent, first entry {s, t};
  * full face-engine run on the 2-loop equal-mass sunrise: the normal and
    pseudo thresholds s-m2 and s-9*m2 must be in the alphabet, no timeouts;
  * incremental flush: the on-disk output JSON must equal the in-memory result;
  * the silent-drop auditor (audit_faces.py): a clean box run must pass
    (exit 0); the same log doctored to claim complete over a timed-out face
    must exit 2 with SILENT DROP; a legacy hand-written COMPLETE_* that
    records 2 strata must exit 2 against the 54-face census that the auditor
    computes from the K4 graph spec alone (face-lattice enumeration only --
    no engine run on the K4 graph in this tier);
  * the Groebner escalation backend (pld_groebner_backend via msolve or
    Singular) forced on the sunrise leading face: must resolve and return the
    known thresholds s-m2 and s-9*m2 (named SKIP if neither engine is
    installed);
  * the mixed quadratic+linear engine (linear_landau.py): Symanzik U/F of the
    1-loop PM box against the hand computation, its Landau alphabet exactly
    {q2, y-1, y+1} and {x, x-1, x+1} in the x chart; the saturated-Groebner
    verify_letter filter must refute the spurious y on the 3PM H face and
    confirm y-1, y+1 on the worldline sub-face; the letter-loss regression
    face must yield y-1, y+1 (resultant self-annihilation guard);
  * the MAIN engine's letter-loss guard (discriminant_letters with the
    Groebner escalation masked off): the box1l leading face must keep BOTH
    s and t, and the self-annihilation reproducer face must yield y-1, y+1
    with status resolved;
  * the topbox spec's topology (graph-theoretic, no engine run): the graph
    with its external legs joined to a common apex must be NON-planar
    (exact Wagner check -- K5/K3,3-minor search over connected branch
    sets), the massive edges must close into a single 4-ring, the massless
    edges must form an open 3-chain between two opposite ring vertices
    with the externals at the chain interiors + the remaining ring
    vertices, and the spec must load with loop number 2.  Guards against
    a planar/broken-ring mis-specification, which fails every prong of
    this leg.

--full: the original engine-driving suites plus the bridge fallback:

  * test_landau_alphabet.py -- the K4-graph legs: 54-face short-timeout
    honesty run, msolve/Singular escalation of the leading K4 face, auditor
    on real engine output (minutes of wall; this is the tier the K4 graph
    lives in);
  * PYTHONHASHSEED=1 test_linear_landau.py -- mixed engine, 8/8 (2PM box,
    3PM IY/H, Apery gamma=3 anomalous threshold);
  * the pld_bridge.jl HC-only numerical fallback on box1l: positive control
    (all torus witnesses match the sympy backend's own letters, exit 0) and
    negative control (t withheld from the candidates must be CAUGHT, exit 2;
    named SKIP without julia or HomotopyContinuation.jl);
  * the Julia engine suite (LandauAlphabet.jl): Pkg.instantiate() then
    test/runtests.jl and test/test_cuspvals.jl (named SKIP if julia is not
    on PATH).

Exits nonzero on failure.  Missing msolve/Singular and missing julia are
named skips -- they are external engines, not requirements of the Python
face engine.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

FAIL = 0


def _leg(tag, ok, detail):
    global FAIL
    print(f"[{tag}] {detail}")
    if ok:
        print("  PASS")
    else:
        print("  FAIL")
        FAIL += 1


def _skip(tag, why):
    print(f"[{tag}] SKIP ({why})")


def _has_letters(got, expected_strs):
    """True iff every expected expression appears (up to a unit) in `got`.

    Letters are re-sympified through their string form: engine letters carry
    spec-local Symbol assumptions and would not cancel against plain symbols."""
    import sympy as sp
    got_e = [sp.sympify(str(g)) for g in got]
    for e in expected_strs:
        e = sp.sympify(e)
        if not any(sp.cancel(e / g).is_number for g in got_e if g != 0):
            return False
    return True


# --------------------------------------------------------------------------- #
#  fast tier
# --------------------------------------------------------------------------- #
def fast_box_engine():
    """Face engine on the 1-loop box: exact alphabet + honesty contract."""
    from landau_alphabet import load_spec, run
    gs = load_spec(os.path.join(HERE, "graph_specs", "box1l.json"))
    log = run(gs, out_path=None, face_timeout=20, global_minutes=5)
    h = log["honesty"]
    ok = (set(log["alphabet"]) == {"s", "t"}
          and log["complete"] is True
          and h["bounded_timeout"] == 0 and h["bounded_partial"] == 0
          and h["missing_faces"] == 0
          and h["attempted_faces"] == h["expected_faces"]
          and "s" in log["first_entry"] and "t" in log["first_entry"])
    _leg("box engine", ok,
         f"box1l alphabet={log['alphabet']} complete={log['complete']} "
         f"resolved={h['resolved']}/{h['n_faces']} first_entry={log['first_entry']}")
    return log


def fast_sunrise_engine():
    """Face engine on the 2-loop sunrise: normal + pseudo thresholds."""
    from landau_alphabet import load_spec, run
    gs = load_spec(os.path.join(HERE, "graph_specs", "sunrise2l.json"))
    log = run(gs, out_path=None, face_timeout=20, global_minutes=5)
    ok = (_has_letters(log["alphabet"], ["s - m2", "s - 9*m2"])
          and log["honesty"]["bounded_timeout"] == 0)
    _leg("sunrise engine", ok,
         f"sunrise2l alphabet={log['alphabet']} (need s-m2 and s-9*m2, no timeouts)")


def fast_flush():
    """Incrementally flushed output JSON equals the in-memory result."""
    from landau_alphabet import load_spec, run
    tmp = os.path.join(tempfile.mkdtemp(prefix="lb_selftest_"), "box1l.json")
    gs = load_spec(os.path.join(HERE, "graph_specs", "box1l.json"))
    log = run(gs, out_path=tmp, face_timeout=20, global_minutes=5)
    with open(tmp) as fh:
        ondisk = json.load(fh)
    ok = (ondisk["alphabet"] == log["alphabet"]
          and ondisk["honesty"]["n_faces"] == log["honesty"]["n_faces"])
    _leg("flush", ok, f"on-disk output == in-memory result ({tmp})")


def fast_auditor():
    """audit_faces.py: exit 0 on a clean log, exit 2 on both silent-drop lies."""
    from landau_alphabet import load_spec, run
    audit = os.path.join(HERE, "audit_faces.py")
    box_spec = os.path.join(HERE, "graph_specs", "box1l.json")
    k4_spec = os.path.join(HERE, "graph_specs", "k4box3l.json")
    gs = load_spec(box_spec)
    log = run(gs, out_path=None, face_timeout=20, global_minutes=5)
    # (0) positive control: the honest log must pass
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(log, fh, default=str)
        clean = fh.name
    r0 = subprocess.run([sys.executable, audit, clean, "--spec", box_spec],
                        capture_output=True, text=True)
    # (1) tally-lie: summary claims complete but a face record is timed_out
    log["faces"][-1]["face_status"] = "timed_out"
    log["complete"] = True
    log["honesty"]["complete"] = True
    log["honesty"]["bounded_timeout"] = 0
    log["honesty"]["unresolved_faces"] = []
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(log, fh, default=str)
        lie = fh.name
    r1 = subprocess.run([sys.executable, audit, lie, "--spec", box_spec],
                        capture_output=True, text=True)
    # (2) legacy COMPLETE_* recording 2 strata vs the true 54-face K4 lattice
    #     (census computed from the graph spec alone -- no engine run)
    legacy = {"task": "x", "graph": "K4 prose", "FACES": [1, 2, 3],
              "COMPLETE_LANDAU_VARIETY": {"strata": {"a": 1, "b": 2}},
              "honesty": ["all covered"]}
    with tempfile.NamedTemporaryFile("w", suffix="COMPLETE_k4box3l.json",
                                     delete=False) as fh:
        json.dump(legacy, fh)
        leg = fh.name
    r2 = subprocess.run([sys.executable, audit, leg, "--spec", k4_spec],
                        capture_output=True, text=True)
    ok = (r0.returncode == 0
          and r1.returncode == 2 and "SILENT DROP" in r1.stdout
          and r2.returncode == 2 and "SILENT DROP" in r2.stdout
          and "of 54" in r2.stdout)
    _leg("auditor", ok,
         f"clean rc={r0.returncode} (need 0); tally-lie rc={r1.returncode} "
         f"(need 2+SILENT DROP); legacy-vs-54-face-census rc={r2.returncode} "
         f"(need 2+SILENT DROP)")


def fast_groebner():
    """Forced Groebner escalation on the sunrise leading face vs known thresholds."""
    import landau_alphabet as L
    if L._GB is None:
        why = ("disabled by PLD_GROEBNER=off"
               if os.environ.get("PLD_GROEBNER") == "off"
               else "neither msolve nor Singular installed -- install either "
                    "to enable")
        _skip("groebner", f"{why}; escalation backend not exercised")
        return
    gs = L.load_spec(os.path.join(HERE, "graph_specs", "sunrise2l.json"))
    sg = L.symanzik_subgraph(gs, set(gs.ALLE))          # leading face
    fv = [gs.ev[e] for e in gs.ALLE]
    g = L._groebner_letters(gs, sg["F"], fv, 60)
    ok = (g is not None and _has_letters(g[0], ["s - m2", "s - 9*m2"]))
    _leg("groebner", ok,
         "sunrise leading face via "
         f"{g[1] if g else 'FAILED'}: letters={sorted(map(str, g[0])) if g else None} "
         "(need s-m2 and s-9*m2)")


def fast_mixed_engine():
    """Mixed quadratic+linear engine: exact 2PM-box alphabet + spurious filter."""
    import sympy as sp
    import linear_landau as LL
    y, q2, x = sp.symbols("y q2 x")

    def aset(ls):
        return {sp.factor(l) for l in ls}

    # Symanzik U/F of the 1-loop PM box vs the hand computation
    U, F, a, _ = LL.symanzik_mixed(*LL.family_2PM_box())
    a1, a2, a3, a4 = a
    ok_uf = (sp.expand(U - (a1 + a2)) == 0
             and sp.expand(F - (a1*a2*q2 - a3**2 - 2*y*a3*a4 - a4**2)) == 0)
    # its Landau alphabet, exactly, in both charts
    res = LL.landau_locus_mixed(*LL.family_2PM_box())
    ok_alpha = (aset(res["alphabet"]) == {q2, y - 1, y + 1}
                and aset(LL.pm_alphabet_x(res["alphabet"])) == {x, x - 1, x + 1})
    # verify_letter must refute the spurious y on the 3PM H face (0,3,4,5,6)
    # and confirm y-1, y+1 on the rung+worldline sub-face
    U3, F3, a3v, _ = LL.symanzik_mixed(*LL.family_3PM_H())
    Ff = sp.expand(F3.subs({a3v[1]: 0, a3v[2]: 0}))
    Fwl = sp.expand(F3.subs({a3v[0]: 0, a3v[1]: 0, a3v[2]: 0, a3v[3]: 0}))
    wl = [a3v[4], a3v[5], a3v[6]]
    ok_ver = (LL.verify_letter(Ff, [a3v[0], a3v[3], a3v[4], a3v[5], a3v[6]],
                               y, [y, q2]) is False
              and LL.verify_letter(Fwl, wl, y - 1, [y, q2]) is True
              and LL.verify_letter(Fwl, wl, y + 1, [y, q2]) is True)
    # letter-loss regression: resultant self-annihilation must not lose y-+1
    b2, b5, b7, b8 = sp.symbols("a2 a5 a7 a8", positive=True)
    P = -b2*b5*(b7**2 + 2*y*b7*b8 + b8**2)
    letters, status = LL._discriminant_letters(
        P, [b2, b5, b7, b8], {y, q2}, [b2, b5, b7, b8], face_timeout=30)
    got = aset(letters)
    ok_loss = (y - 1) in got and (y + 1) in got
    _leg("mixed engine", ok_uf and ok_alpha and ok_ver and ok_loss,
         f"2PM box U/F exact={ok_uf} alphabet exact={ok_alpha} "
         f"spurious-y refuted + y-+1 confirmed={ok_ver} "
         f"letter-loss guard={ok_loss} (status {status})")


def fast_letterloss_main():
    """Letter-loss guard in the MAIN engine's bounded resultant chain.

    discriminant_letters must strip pure-x monomial content each round and
    harvest kinematic factors from every intermediate generator set, so
    resultant self-annihilation cannot silently empty a face's letter set
    while the face reports 'resolved'.  Runs with the Groebner escalation
    masked off so the pure-sympy chain itself is exercised."""
    import landau_alphabet as L
    saved = L._GB
    L._GB = None            # exercise the bounded chain, no escalation
    try:
        # box1l leading face: the chain must keep BOTH s and t on its own
        gs = L.load_spec(os.path.join(HERE, "graph_specs", "box1l.json"))
        sg = L.symanzik_subgraph(gs, set(gs.ALLE))
        fv = [gs.ev[e] for e in gs.ALLE]
        lets, st = L.discriminant_letters(gs, sg["F"], fv, 20)
        ok_box = _has_letters(lets, ["s", "t"]) and st == "resolved"
        # self-annihilation reproducer (mirror of the mixed-engine regression):
        # shared pure-x monomial content across generators
        spec = {"name": "loss8", "edges": [[1, 2, 0]] * 8,
                "nodes": {"1": [1], "2": [2]},
                "kinematics": {"invariants": ["y"], "channels": {"1": "y"}}}
        gs2 = L.GraphSpec(spec)
        y = gs2.invariants[0]
        x2, x5, x7, x8 = (gs2.ev[e] for e in ("x2", "x5", "x7", "x8"))
        P = -x2 * x5 * (x7**2 + 2 * y * x7 * x8 + x8**2)
        lets2, st2 = L.discriminant_letters(gs2, P, [x2, x5, x7, x8], 30)
        ok_rep = _has_letters(lets2, ["y - 1", "y + 1"]) and st2 == "resolved"
    finally:
        L._GB = saved
    _leg("letter-loss (main engine)", ok_box and ok_rep,
         f"box1l leading face without escalation: {sorted(map(str, lets))} ({st}, "
         f"need s AND t); self-annihilation face: {sorted(map(str, lets2))} "
         f"({st2}, need y-1 and y+1)")


def _graph_is_planar(verts, pairs):
    """Exact planarity via Wagner's theorem: planar iff no K5 and no K3,3 minor.

    Brute-force minor search over connected branch sets (bitmask BFS).  With
    six branch sets pairwise disjoint each has at most n-5 vertices (n-4 for
    the five K5 sets), so the candidate pool stays tiny for the shipped
    graph specs (<= 8 vertices including an external apex)."""
    idx = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    adj = [0] * n
    for a, b in pairs:
        if a != b:
            adj[idx[a]] |= 1 << idx[b]
            adj[idx[b]] |= 1 << idx[a]

    def connected(mask):
        seen = front = mask & -mask
        while front:
            nxt = 0
            m = front
            while m:
                v = (m & -m).bit_length() - 1
                m &= m - 1
                nxt |= adj[v] & mask
            front = nxt & ~seen
            seen |= front
        return seen == mask

    def reach(mask):
        r = 0
        m = mask
        while m:
            v = (m & -m).bit_length() - 1
            m &= m - 1
            r |= adj[v]
        return r & ~mask

    def minor(k, bipartite):
        cap = n - k + 1                     # each branch set leaves >=1 vertex
        cand = [m for m in range(1, 1 << n)
                if bin(m).count("1") <= cap and connected(m)]
        rs = {m: reach(m) for m in cand}

        def grow(chosen, used, start):
            if len(chosen) == k:
                return True
            need = (chosen[:3] if bipartite and len(chosen) >= 3
                    else (chosen if not bipartite else []))
            for i in range(start if not (bipartite and len(chosen) == 3) else 0,
                           len(cand)):
                m = cand[i]
                if m & used:
                    continue
                if all(rs[c] & m for c in need):
                    if grow(chosen + [m], used | m, i + 1):
                        return True
            return False

        return grow([], 0, 0)

    return not (minor(5, False) or minor(6, True))


def _topbox_topology_ok(spec_path):
    """Graph-theoretic contract of the topbox spec (no engine run).

    (a) NON-planar once the external legs are joined to a common apex
        (planarity of a diagram is planarity of graph + apex vertex);
    (b) the massive edges close into a single ring (here of length 4);
    (c) the massless edges form an open chain between two NON-adjacent
        (opposite) ring vertices, the chain interior carries external legs,
        and the two ring vertices off the chain carry the other externals;
    (d) the spec loads and its full-edge face has loop number 2.
    Returns (ok, detail)."""
    from landau_alphabet import load_spec, symanzik_subgraph
    with open(spec_path) as fh:
        spec = json.load(fh)
    edges = []
    for e in spec["edges"]:
        v1, v2 = int(e[0]), int(e[1])
        m2 = e[2] if len(e) == 3 else 0
        edges.append((v1, v2, m2 not in (0, "0")))
    verts = sorted({v for a, b, _ in edges for v in (a, b)})
    ext_verts = sorted(int(v) for v, legs in spec.get("nodes", {}).items()
                       if legs)

    # (a) planarity with the external apex
    apex = max(verts) + 1
    pairs = [(a, b) for a, b, _ in edges] + [(apex, v) for v in ext_verts]
    nonplanar = not _graph_is_planar(verts + [apex], pairs)

    # (b) massive ring: every touched vertex has massive-degree 2, connected
    mass = [(a, b) for a, b, heavy in edges if heavy]
    mdeg = {}
    for a, b in mass:
        mdeg[a] = mdeg.get(a, 0) + 1
        mdeg[b] = mdeg.get(b, 0) + 1
    ring_verts = sorted(mdeg)
    ring_ok = (len(mass) == 4 and len(ring_verts) == 4
               and all(d == 2 for d in mdeg.values()))
    if ring_ok:                             # single cycle, not two 2-cycles
        comp = {ring_verts[0]}
        grew = True
        while grew:
            grew = False
            for a, b in mass:
                if (a in comp) != (b in comp):
                    comp |= {a, b}
                    grew = True
        ring_ok = comp == set(ring_verts)

    # (c) massless open chain between two opposite ring vertices
    light = [(a, b) for a, b, heavy in edges if not heavy]
    ldeg = {}
    for a, b in light:
        ldeg[a] = ldeg.get(a, 0) + 1
        ldeg[b] = ldeg.get(b, 0) + 1
    chain_ends = sorted(v for v, d in ldeg.items() if d == 1)
    chain_mid = sorted(v for v, d in ldeg.items() if d == 2)
    chain_ok = (ring_ok and len(light) == 3
                and len(chain_ends) == 2 and len(chain_mid) == 2
                and all(v in ring_verts for v in chain_ends)
                and not any(set(e) == set(chain_ends) for e in mass)
                and all(v not in ring_verts for v in chain_mid)
                and sorted(chain_mid
                           + [v for v in ring_verts if v not in chain_ends])
                    == ext_verts)

    # (d) the engine can load it; loop number 2 on the full-edge face
    gs = load_spec(spec_path)
    loops = symanzik_subgraph(gs, set(gs.ALLE))["L"]

    ok = nonplanar and ring_ok and chain_ok and loops == 2
    detail = (f"nonplanar(graph+apex)={nonplanar} massive-4-ring={ring_ok} "
              f"massless-chain+externals={chain_ok} loops={loops}")
    return ok, detail


def fast_topbox_topology():
    """topbox spec: non-planar topology + closed massive ring (graph check)."""
    ok, detail = _topbox_topology_ok(
        os.path.join(HERE, "graph_specs", "topbox.json"))
    _leg("topbox topology", ok, detail)


def run_fast():
    t0 = time.time()
    try:
        import sympy                                    # noqa: F401
    except ImportError:
        print("FAIL: requirement sympy is not installed (pip install sympy) -- "
              "the face engine cannot run")
        return 1
    fast_box_engine()
    fast_sunrise_engine()
    fast_flush()
    fast_auditor()
    fast_groebner()
    fast_mixed_engine()
    fast_letterloss_main()
    fast_topbox_topology()
    verdict = "FAIL" if FAIL else "PASS"
    print(f"fast battery {verdict}  ({time.time() - t0:.1f}s wall; K4 engine "
          f"runs, the full pytest suites and the Julia engine suite are "
          f"behind --full)")
    return 1 if FAIL else 0


# --------------------------------------------------------------------------- #
#  --full tier: the original engine-driving suites, unchanged
# --------------------------------------------------------------------------- #
def _full_cmd(tag, cmd, cwd=HERE, env_extra=None):
    global FAIL
    print(f"\n=== [{tag}] {' '.join(cmd)}")
    env = dict(os.environ, **(env_extra or {}))
    r = subprocess.run(cmd, cwd=cwd, env=env)
    if r.returncode == 0:
        print(f"=== [{tag}] PASS")
    else:
        print(f"=== [{tag}] FAIL (exit {r.returncode})")
        FAIL += 1


def full_pld_bridge_fallback():
    """pld_bridge.jl's HomotopyContinuation-only numerical fallback, box1l.

    The fallback witnesses INTERIOR torus pinches of each face's G = U + F.
    On box1l the one such pinch sits on the second-type u-channel surface
    s+t (x = (1/s, 1/t, 1/s, 1/t) forces s+t = 0) -- a surface the DE-level
    alphabet {s, t} legitimately prunes; the endpoint-class letters s and t
    have no interior pinch and are the sympy chain's territory.  Positive
    control: candidates = the sympy letters plus s+t -> every witness
    matches, s+t is witnessed, exit 0.  Negative control: candidates = the
    sympy output alone -> the second-type witness must be CAUGHT (exit 2).
    Named SKIP without julia or without HomotopyContinuation.jl."""
    import shutil
    if shutil.which("julia") is None:
        _skip("pld bridge HC fallback", "julia not on PATH -- HC fallback "
              "not exercised")
        return
    probe = subprocess.run(["julia", "--startup-file=no", "-e",
                            "import HomotopyContinuation"],
                           capture_output=True, text=True)
    if probe.returncode != 0:
        _skip("pld bridge HC fallback", "HomotopyContinuation.jl not installed "
              "in the default julia env -- HC fallback not exercised")
        return
    print("\n=== [pld bridge HC fallback] box1l witnesses vs the sympy backend "
          "(HC probe + two bridge runs; each pays julia+HC startup)", flush=True)
    from landau_alphabet import load_spec, run
    spec = os.path.join(HERE, "graph_specs", "box1l.json")
    bridge = os.path.join(HERE, "pld_bridge.jl")
    tmpd = tempfile.mkdtemp(prefix="lb_bridge_")
    sympy_out = os.path.join(tmpd, "box1l_sympy.json")
    log_s = run(load_spec(spec), out_path=sympy_out, face_timeout=20,
                global_minutes=5)
    # positive control: sympy letters + the known second-type letter s+t
    cand = os.path.join(tmpd, "box1l_candidates.json")
    with open(cand, "w") as fh:
        json.dump({"alphabet": list(log_s["alphabet"]) + ["s + t"]}, fh)
    outj = os.path.join(tmpd, "box1l_hc.json")
    r = subprocess.run(["julia", "--startup-file=no", bridge, spec, outj,
                        "--fallback", "--check", cand, "--seed", "0"],
                       capture_output=True, text=True)
    ok = False
    detail = f"positive control rc={r.returncode} (need 0)"
    if r.returncode == 0 and os.path.exists(outj):
        with open(outj) as fh:
            log = json.load(fh)
        wit = sum(len(f.get("t_witnesses", [])) for f in log["faces"])
        ok = (log["complete"] is False
              and not log["honesty"]["unmatched_witnesses"]
              and wit > 0
              and "s + t" in log["alphabet"])
        detail += (f"; witnesses={wit} matched-alphabet={log['alphabet']} "
                   f"(need s + t witnessed) complete={log['complete']} (need false)")
    else:
        detail += f"; stderr tail: {r.stderr.strip()[-300:]}"
    # negative control: the sympy alphabet ALONE must be caught as missing
    # the second-type surface (exit 2)
    r2 = subprocess.run(["julia", "--startup-file=no", bridge, spec,
                         os.path.join(tmpd, "box1l_hc2.json"),
                         "--fallback", "--check", sympy_out, "--seed", "0"],
                        capture_output=True, text=True)
    ok2 = (r2.returncode == 2 and "CROSS-VALIDATION CATCH" in r2.stdout)
    _leg("pld bridge HC fallback", ok and ok2,
         detail + f"; negative control (sympy letters alone) rc={r2.returncode} "
                  f"(need 2 + CROSS-VALIDATION CATCH)")


def run_full():
    t0 = time.time()
    _full_cmd("face engine suite", [sys.executable,
                                    os.path.join(HERE, "test_landau_alphabet.py")])
    _full_cmd("mixed engine suite", [sys.executable,
                                     os.path.join(HERE, "test_linear_landau.py")],
              env_extra={"PYTHONHASHSEED": "1"})
    full_pld_bridge_fallback()
    import shutil
    if shutil.which("julia") is None:
        _skip("julia engine suite", "julia not on PATH -- LandauAlphabet.jl "
              "suite not exercised; install Julia >= 1.11 to enable")
    else:
        jl = os.path.join(HERE, "LandauAlphabet.jl")
        _full_cmd("julia instantiate",
                  ["julia", f"--project={jl}", "-e",
                   "using Pkg; Pkg.instantiate()"])
        _full_cmd("julia engine suite",
                  ["julia", f"--project={jl}",
                   os.path.join(jl, "test", "runtests.jl")])
        _full_cmd("julia cusp values",
                  ["julia", f"--project={jl}",
                   os.path.join(jl, "test", "test_cuspvals.jl")])
    verdict = "FAIL" if FAIL else "PASS"
    print(f"\nfull battery {verdict}  ({time.time() - t0:.1f}s wall)")
    return 1 if FAIL else 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python3 selftest.py",
        description="landau-alphabet validation battery")
    ap.add_argument("--full", action="store_true",
                    help="run the original engine-driving suites (K4 face "
                         "engine + auditor, mixed engine 8/8, Julia engine "
                         "suite) instead of the fast code-check battery")
    args = ap.parse_args(argv)
    sys.exit(run_full() if args.full else run_fast())


if __name__ == "__main__":
    main()

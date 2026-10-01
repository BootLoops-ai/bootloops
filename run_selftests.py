#!/usr/bin/env python3
"""run_selftests.py — verify a BootLoops clone in one command.

Walks every package under tools/, reads its declared verification class from
tools/README.md (selftest / partial / smoke / data-gated), looks up its battery in
tools/BATTERIES.json, runs it, and prints a class-aware table. Designed for a cold clone: legs that
need an external engine or gated data are expected to skip or fail closed with
a named error — that is not breakage.

    python3 run_selftests.py                # everything
    python3 run_selftests.py mixalot eras   # named packages only
    python3 run_selftests.py --par 8        # parallel (batteries must not
                                            #   write outside their package)
    python3 run_selftests.py --timeout 600  # per-attempt timeout (default 300)

Exit code: 1 if any package FAILs (a data-gated REFUSED is not a failure), else 0. The full per-package
record (commands tried, exit codes, output tails) lands in selftest_results.json.

Each package's canonical battery command and working directory are pinned in
tools/BATTERIES.json — the runner executes exactly that entry, once. There is
no discovery and no guessing: a package absent from the manifest is reported
no-manifest-entry with the fix named. Adding a package to the repo means
adding its battery line to the manifest, same as its row in tools/README.md.
Prerequisites and how to read the results table: INSTALL.md,
"Verifying an installation".
"""
import argparse, contextlib, glob, json, os, re, signal, subprocess, sys, threading, time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, "tools")
STANDARD_S = 20   # a package's default battery completes in ~20 s on a laptop
                  # (tools/README.md, "Time standard"); heavier tiers go behind flags
# Julia batteries share one depot: concurrent precompiles/loads race its locks
# and fail nondeterministically at high --par, so they run one at a time.
JULIA_LOCK = threading.Lock()
def classes():
    out = {}
    for line in open(os.path.join(TOOLS, "README.md"), encoding="utf-8"):
        m = re.match(r'\| \`([A-Za-z0-9_.-]+)/\` \|.*\| ([a-z-]+)', line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


MANIFEST = {}
_mp = os.path.join(TOOLS, "BATTERIES.json")
if os.path.isfile(_mp):
    MANIFEST = json.load(open(_mp))


def run_pkg(pkg, cls, timeout):
    d = os.path.join(TOOLS, pkg)
    if not os.path.isdir(d):
        return {"class": cls, "status": "missing-dir"}
    if pkg not in MANIFEST:
        # no guessing: every package declares its battery in tools/BATTERIES.json
        return {"class": cls, "status": "no-manifest-entry",
                "fix": f"add a line for {pkg!r} to tools/BATTERIES.json (cmd + cwd)"}
    m = MANIFEST[pkg]
    cands = [("manifest", m["cmd"])]
    forced_cwd = d if m.get("cwd") == "package" else HERE
    uses_julia = ("julia" in m["cmd"] or m.get("serial")
                  or os.path.isfile(os.path.join(d, "Project.toml")))
    gate = JULIA_LOCK if uses_julia else contextlib.nullcontext()
    # Julia-component packages: make the committed environment concrete first
    if os.path.isfile(os.path.join(d, "Project.toml")):
        with JULIA_LOCK:
            subprocess.run("julia --project=. -e 'using Pkg; Pkg.instantiate()'",
                           shell=True, cwd=d, capture_output=True, timeout=1800)
    pypath = [TOOLS, d, HERE]
    if os.path.isdir(os.path.join(d, "src")):        # src-layout packages
        pypath.insert(0, os.path.join(d, "src"))
    env = dict(os.environ, PYTHONPATH=":".join(pypath + [os.environ.get("PYTHONPATH", "")]))
    entry, usage_only = {"class": cls, "tried": []}, True
    attempts = [(k, c, forced_cwd) for k, c in cands]
    for kind, cmd, cwd in attempts:
      with gate:
        t0 = time.monotonic()
        try:
            # own process group per battery: some batteries fork daemons that
            # inherit the pipes — on timeout the WHOLE group must die, or the
            # pipe never reaches EOF and the runner hangs (ops/turnstile, macOS)
            p = subprocess.Popen(cmd, shell=True, cwd=cwd, env=env, text=True,
                                 stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, start_new_session=True)
            try:
                out, _ = p.communicate(timeout=timeout)
                rc = p.returncode
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    p.kill()
                try:
                    out, _ = p.communicate(timeout=15)
                except Exception:
                    out = ""
                rc = -9
                out = (out or "") + f"\nTIMEOUT after {timeout}s (process group killed)"
            finally:
                # a stray double-forked survivor must not hold the NEXT battery hostage
                try:
                    os.killpg(p.pid, signal.SIGKILL)
                except Exception:
                    pass
            tail = (out or "")[-1200:]
        except OSError as ex:
            rc, tail = -8, f"could not run: {ex}"
        entry["secs"] = round(time.monotonic() - t0, 1)
        entry["tried"].append({"kind": kind, "cmd": cmd, "cwd": cwd, "rc": rc,
                               "wall_s": round(time.monotonic() - t0, 1), "tail": tail})
        entry["wall_s"] = round(time.monotonic() - t0, 1)
        if not (rc == 2 and re.search(r'usage:|required:|unrecognized arguments', tail)):
            usage_only = False
        if rc == 0:
            entry.update(status="PASS", via=cmd)
            return entry
    if usage_only:
        entry["status"] = "entry-not-found"
    elif cls == "data-gated":
        entry["status"] = "REFUSED (by design)"   # data-gated: loud refusal is the contract
    else:
        entry["status"] = "FAIL"
    return entry


def _overtime(e):
    # soft warning only — never changes status or exit code
    secs = e.get("secs")
    if secs and secs > 3 * STANDARD_S:
        return f"  over the ~{STANDARD_S} s battery standard (took {secs:.0f}s)"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("packages", nargs="*", help="restrict to these packages")
    ap.add_argument("--par", type=int, default=1, help="parallel workers (default 1)")
    ap.add_argument("--timeout", type=int, default=300)
    args = ap.parse_args()
    cls = classes()
    todo = args.packages or sorted(cls, key=lambda p: (cls[p] != "selftest", p))
    unknown = [p for p in todo if p not in cls]
    if unknown:
        sys.exit(f"not in tools/README.md index: {unknown}")
    results = {}

    def one(p):
        return p, run_pkg(p, cls[p], args.timeout)

    if args.par > 1:
        import concurrent.futures as cf
        with cf.ThreadPoolExecutor(max_workers=args.par) as ex:
            it = ex.map(one, todo)
            for p, e in it:
                results[p] = e
                print(f"{e['status']:>18}  [{cls[p]:>10}]  {p}{_overtime(e)}", flush=True)
    else:
        for p in todo:
            p, e = one(p)
            results[p] = e
            print(f"{e['status']:>18}  [{cls[p]:>10}]  {p}{_overtime(e)}", flush=True)

    json.dump(results, open(os.path.join(HERE, "selftest_results.json"), "w"), indent=1)
    counts = {}
    for e in results.values():
        counts[f"{e['class']}/{e['status']}"] = counts.get(f"{e['class']}/{e['status']}", 0) + 1
    print(json.dumps(counts, indent=1, sort_keys=True))
    hard = [p for p, e in results.items() if e["status"] in ("FAIL", "no-manifest-entry")]
    if hard:
        print(f"failures (every class must run green or refuse by design from a cold clone): {sorted(hard)}")
        sys.exit(1)


if __name__ == "__main__":
    main()

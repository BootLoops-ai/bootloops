#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Part of BootLoops 1.0 (repository-root
# LICENSE and NOTICE).
"""ffcapital selftest — smoke checks, the REAL synthetic control battery,
then the designed data-gate refusal.

The ff_save-decoding gates run against a user-supplied FireFly `ff_save/`
state; no public fixtures ship (GUIDE.md BATTERY). But the recon_symbolic_d
synthetic-truth control (gate d) needs only a slice-name grid and the FireFly
prime table, both synthesizable here — so it runs FOR REAL:
  1. `--help` exits 0 on both harvest members (the CLIs import and parse),
  2. `recon_symbolic_d.py control` with a scratch --outdir aborts with the
     named unset CONFIG (the fail-closed refusal fires as designed),
  3. control legs (skipped BY NAME if python-flint is absent):
     a. reference-style grid (contains d=118/29): full battery must PASS,
     b. generic grid (no d=118/29): the grid-parametrized plant must be
        picked from the grid and the full battery must PASS,
     c. too-small grid: the named CONTROL refusal must fire (nonzero rc),
        never a crash.
Without user ff_save data the exit is always nonzero: the final line names
the data gate for the harvest members. Point the member CLIs at your save
(GUIDE.md INVOKE) to run the full gate battery.
"""
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def _run(args):
    return subprocess.run([sys.executable] + args, cwd=HERE,
                          capture_output=True, text=True, timeout=120)


def _is_prime_u64(n):
    """Deterministic Miller-Rabin for n < 2^64 (fixed witness set)."""
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _write_prime_fixture(path, count=40):
    """The largest 63-bit primes, generated (never hand-typed); the member's
    own firefly_primes() assert gates the table head."""
    ps = []
    n = 2 ** 63 - 1
    while len(ps) < count:
        if _is_prime_u64(n):
            ps.append(n)
        n -= 2
    with open(path, "w") as f:
        f.write("// selftest fixture: the 300 largest 63-bit primes table,\n"
                "// truncated to the first %d entries.\n" % count)
        f.write("static const uint64_t primes[] = {"
                + ", ".join(f"{q}uLL" for q in ps) + "};\n")


def _write_cov(path, names):
    with gzip.open(path, "wt") as f:
        f.write(json.dumps({"fn": 0, "tag": "SYNTH", "nd": 4, "dd": 3,
                            "slices": names}) + "\n")


def _control(cov, helper, outdir):
    return _run([os.path.join(HERE, "recon_symbolic_d.py"), "control",
                 "--coverage", cov, "--ff-helper", helper,
                 "--outdir", outdir])


def control_legs(scratch):
    """Legs 3a-3c. Returns 0 on pass, 1 on any FAIL."""
    helper = os.path.join(scratch, "helper_fixture.cpp")
    _write_prime_fixture(helper)
    # d = 4 - 2a/b per slice name; the reference-style grid contains
    # node_eps_-1_29 (d=118/29) and node_eps_-3_26 (d=55/13)
    grids = {
        "reference": ["node_eps_-1_29", "node_eps_-3_26"]
                     + [f"node_eps_{k}_17" for k in range(-9, 10)],
        "generic": [f"node_eps_{k}_17" for k in range(-10, 11)],
        "small": [f"node_eps_{k}_17" for k in range(5)],
    }
    for name, names in grids.items():
        cov = os.path.join(scratch, f"cov_{name}.jsonl.gz")
        _write_cov(cov, names)
        out = os.path.join(scratch, f"out_{name}")
        os.makedirs(out, exist_ok=True)
        r = _control(cov, helper, out)
        msg = r.stdout + r.stderr
        if name == "small":
            if r.returncode == 0 or "CONTROL: grid has" not in msg:
                print("CONTROL FAIL: too-small grid did not hit the named "
                      f"refusal (rc={r.returncode})")
                print(msg.strip()[-400:])
                return 1
            print(f"control ok: {name} grid refused by name "
                  f"(rc={r.returncode})")
            continue
        if r.returncode != 0 or "CONTROL VERDICT: PASS" not in msg:
            print(f"CONTROL FAIL: {name} grid rc={r.returncode}")
            print(msg.strip()[-600:])
            return 1
        with open(os.path.join(out, "control_result.json")) as f:
            res = json.load(f)
        plant = res["degrees"]["planted_cancellation"]
        if name == "generic" and "118/29" in plant:
            print("CONTROL FAIL: generic grid still uses the hard-wired "
                  f"plant ({plant})")
            return 1
        print(f"control ok: {name} grid PASS "
              f"({len(res['cases'])} cases; {plant})")
    return 0


def main():
    for member in ("harvest_ffsave_eta.py", "harvest_ffsave_2var.py"):
        r = _run([os.path.join(HERE, member), "--help"])
        if r.returncode != 0:
            print(f"SMOKE FAIL: {member} --help rc={r.returncode}")
            print((r.stdout + r.stderr)[-400:])
            return 1
        print(f"smoke ok: {member} --help rc=0")
    scratch = tempfile.mkdtemp(prefix="ffcapital_selftest_")
    r = _run([os.path.join(HERE, "recon_symbolic_d.py"), "control",
              "--outdir", scratch])
    msg = (r.stdout + r.stderr).strip()
    if r.returncode == 0 or not re.search(r"CONFIG: [A-Z0-9_]+ unset", msg):
        print("SMOKE FAIL: control did not refuse with a named unset CONFIG")
        print(msg[-400:])
        return 1
    print(f"smoke ok: control refuses fail-closed ({msg.splitlines()[-1]})")
    try:
        import flint  # noqa: F401
    except ImportError:
        print("SKIP (by name): synthetic control legs — python-flint absent "
              "(recon_symbolic_d hard-requires it)")
    else:
        if control_legs(scratch):
            return 1
    print("DATA GATE: no ff_save/ state supplied — the harvest-member gates "
          "run against a user-supplied FireFly ff_save/ reconstruction state "
          "and no public fixtures ship (see GUIDE.md INVOKE). The "
          "recon_symbolic_d synthetic control battery above is the "
          "data-free part.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

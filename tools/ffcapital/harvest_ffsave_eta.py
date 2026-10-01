#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Part of BootLoops 1.0 (repository-root
# LICENSE and NOTICE). Reads the plain-text state files FireFly writes; contains
# no FireFly code (THIRD_PARTY.md B2).
# ffcapital member — 1-var salvage (numeric-d eta-slices).
# Siblings: harvest_ffsave_2var.py (2-var), recon_symbolic_d.py (cross-slice
# symbolic-d). All roots via CLI/--config; fail-closed on unset CONFIG.
r"""harvest_ffsave_eta.py — farm harvest shape (A): extract per-fn
reconstructed 1-var eta rational functions DIRECTLY from a stopped node's
ff_save (no job completion / kira_target.m needed).

Ground truth format: RatReconst::save_state(),
FireFly source RatReconst.cpp:2537-2716 (Klappert-Klein-Lange)
(g_ni/g_di lines = "e1 e2 <num> <den>" — the ALREADY-RECONSTRUCTED rational
coefficients; n=2 vars [d,eta] with d numerically substituted in the SYSTEM
files => e1 == 0 for every monomial, asserted). A fn with is_done=1 has its
COMPLETE function in g_ni/g_di: f(eta) = sum(g_ni)/sum(g_di).

Independent validation (no fabrication): ff_save/validation.gz is written by
Reconstructor (Reconstructor.hpp:2157ff) as
  line 1:   the n black-box input values (mod the prime current at write time)
  line i+2: the BLACK-BOX result of fn i at that input (item order = tag order)
These are black-box probe values, independent of the reconstruction state.
For each extracted done fn we require  N(x2)/D(x2) == stored  (mod p), with the
prime index auto-detected once (consistent across fns) from the FireFly table.

Modes:
  census  : all states -> jsonl {fn, tag, done, np(exact CRT depth), degN, degD}
  extract : given --fns (or --all-done) -> jsonl with exact rational coeffs
            {fn, tag, done, degN, degD, N: {e: "p/q"}, D: {e: "p/q"}, validated}

usage:
  harvest_ffsave_eta.py census  --ffsave <dir> --out census.jsonl [--workers 16]
  harvest_ffsave_eta.py extract --ffsave <dir> --out funcs.jsonl
                                (--fns 0,17,26683 | --all-done) [--workers 16]
                                [--no-validate]
"""
import argparse, glob, gzip, json, os, sys
from fractions import Fraction
from multiprocessing import Pool

# FireFly 63-bit prime table head (ReconstHelper.cpp:24ff), enough for 40 primes
PRIMES = [
    9223372036854775783, 9223372036854775643, 9223372036854775549,
    9223372036854775507, 9223372036854775433, 9223372036854775421,
    9223372036854775417, 9223372036854775399, 9223372036854775351,
    9223372036854775337, 9223372036854775291, 9223372036854775279,
    9223372036854775259, 9223372036854775181, 9223372036854775159,
    9223372036854775139, 9223372036854775097, 9223372036854775073,
    9223372036854775057, 9223372036854774959, 9223372036854774937,
    9223372036854774917, 9223372036854774893, 9223372036854774797,
    9223372036854774739, 9223372036854774713, 9223372036854774679,
    9223372036854774629, 9223372036854774587, 9223372036854774571,
    9223372036854774559, 9223372036854774511, 9223372036854774509,
    9223372036854774499, 9223372036854774451, 9223372036854774413,
    9223372036854774341, 9223372036854774319, 9223372036854774307,
    9223372036854774277]

KEYWORDS = {"combined_prime", "tag_name", "is_done", "max_deg_num", "max_deg_den",
            "individual_degrees_num", "individual_degrees_den", "need_prime_shift",
            "normalizer_deg", "normalize_to_den", "normalizer_den_num",
            "shifted_max_num_eqn", "shift", "sub_num", "sub_den", "zero_degs_num",
            "zero_degs_den", "g_ni", "g_di", "combined_ni", "combined_di",
            "combined_primes_ni", "combined_primes_di", "interpolations"}

FFDIR = None  # set in main; module-global so Pool workers inherit it


def n_primes_exact(cp):
    k = 0
    while cp > 1:
        if k >= len(PRIMES) or cp % PRIMES[k]:
            return -1
        cp //= PRIMES[k]
        k += 1
    return k


def state_path(fn):
    hits = glob.glob(os.path.join(FFDIR, "states", f"{fn}_*.gz"))
    if len(hits) != 1:
        raise RuntimeError(f"fn {fn}: {len(hits)} state files (want exactly 1): {hits}")
    return hits[0]


def parse_state(path):
    sections, cur = {}, None
    with gzip.open(path, "rt") as f:
        first = f.readline().rstrip("\n")
        if first == "ZERO":
            return {"zero": True}
        cur = first
        sections[cur] = []
        for line in f:
            line = line.rstrip("\n")
            if line in KEYWORDS:
                cur = line
                sections[cur] = []
            elif line.strip():
                sections[cur].append(line)
    return sections


def poly_from_g(lines, fn, side):
    """g_ni/g_di lines 'e1 e2 num den' -> {e2: Fraction}. Asserts e1==0."""
    P = {}
    for ln in lines:
        p = ln.split()
        assert len(p) == 4, f"fn {fn} {side}: monomial line has {len(p)} fields (want 4): {ln!r}"
        e1, e2 = int(p[0]), int(p[1])
        assert e1 == 0, f"fn {fn} {side}: e1={e1} != 0 — d-dependence present, NOT a numeric-d slice"
        assert e2 not in P, f"fn {fn} {side}: duplicate exponent {e2}"
        P[e2] = Fraction(int(p[2]), int(p[3]))
    return P


def extract_one(fn):
    sec = parse_state(state_path(fn))
    if sec.get("zero"):
        return {"fn": fn, "zero": True, "done": 1, "N": {}, "D": {"0": "1"}}
    done = int(sec["is_done"][0])
    N = poly_from_g(sec.get("g_ni", []), fn, "g_ni")
    D = poly_from_g(sec.get("g_di", []), fn, "g_di")
    return {"fn": fn, "tag": sec["tag_name"][0], "done": done,
            "np": n_primes_exact(int(sec["combined_prime"][0])),
            "degN": max(N) if N else 0, "degD": max(D) if D else 0,
            "N": {str(e): f"{c.numerator}/{c.denominator}" for e, c in sorted(N.items())},
            "D": {str(e): f"{c.numerator}/{c.denominator}" for e, c in sorted(D.items())}}


def census_one(path):
    fn = int(os.path.basename(path).split("_")[0])
    sec = parse_state(path)
    if sec.get("zero"):
        return {"fn": fn, "zero": True}
    gn = sec.get("g_ni", [])
    gd = sec.get("g_di", [])
    em = lambda ls: max((int(l.split()[1]) for l in ls), default=0)
    return {"fn": fn, "tag": sec["tag_name"][0], "done": int(sec["is_done"][0]),
            "np": n_primes_exact(int(sec["combined_prime"][0])),
            "mdn": int(sec["max_deg_num"][0]), "mdd": int(sec["max_deg_den"][0]),
            "degN": em(gn), "degD": em(gd), "ngn": len(gn), "ngd": len(gd),
            "ncn": len(sec.get("combined_ni", [])), "ncd": len(sec.get("combined_di", []))}


# ---------------- validation against ff_save/validation.gz (black-box values)
def load_validation():
    vals = []
    with gzip.open(os.path.join(FFDIR, "validation.gz"), "rt") as f:
        pt = [int(x) for x in f.readline().split()]
        for ln in f:
            ln = ln.strip()
            if ln:
                vals.append(int(ln))
    return pt, vals


def eval_mod(P, x, p):
    """{e2: Fraction} at eta=x mod p."""
    r = 0
    for e, c in P.items():
        t = (c.numerator % p) * pow(c.denominator % p, p - 2, p) % p
        r = (r + t * pow(x, e, p)) % p
    return r


def detect_prime(recs, pt, vals):
    """Find the unique prime index under which extracted done fns reproduce
    validation.gz. Uses up to 8 nonzero-den done fns; requires unanimity."""
    sample = [r for r in recs if r.get("done") and not r.get("zero")][:8]
    if not sample:
        return None
    cands = []
    for pi, p in enumerate(PRIMES):
        x2 = pt[1] % p
        ok = True
        for r in sample:
            N = {int(e): Fraction(v) for e, v in r["N"].items()}
            D = {int(e): Fraction(v) for e, v in r["D"].items()}
            dv = eval_mod(D, x2, p)
            if dv == 0 or eval_mod(N, x2, p) * pow(dv, p - 2, p) % p != vals[r["fn"]] % p:
                ok = False
                break
        if ok:
            cands.append(pi)
    return cands


def main():
    global FFDIR
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["census", "extract"])
    ap.add_argument("--ffsave", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fns", default=None, help="comma list of fn indices")
    ap.add_argument("--all-done", action="store_true")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--no-validate", action="store_true")
    a = ap.parse_args()
    FFDIR = a.ffsave

    if a.mode == "census":
        files = sorted(glob.glob(os.path.join(FFDIR, "states", "*.gz")))
        print(f"[census] {len(files)} state files", file=sys.stderr)
        with Pool(a.workers) as pool, open(a.out, "w") as out:
            for rec in pool.imap_unordered(census_one, files, chunksize=128):
                out.write(json.dumps(rec) + "\n")
        print(f"[census] wrote {a.out}", file=sys.stderr)
        return

    if a.fns:
        fns = [int(x) for x in a.fns.split(",")]
    elif a.all_done:
        files = glob.glob(os.path.join(FFDIR, "states", "*.gz"))
        with Pool(a.workers) as pool:
            cen = pool.map(census_one, files, chunksize=128)
        fns = sorted(r["fn"] for r in cen if r.get("zero") or r.get("done"))
        print(f"[extract] {len(fns)} done fns of {len(files)}", file=sys.stderr)
    else:
        ap.error("extract needs --fns or --all-done")

    with Pool(a.workers) as pool:
        recs = pool.map(extract_one, fns, chunksize=64)

    if not a.no_validate:
        pt, vals = load_validation()
        assert len(vals) == len(glob.glob(os.path.join(FFDIR, "states", "*.gz"))), \
            "validation.gz entry count != number of state files"
        cands = detect_prime(recs, pt, vals)
        assert cands and len(cands) == 1, f"prime autodetect ambiguous/failed: {cands}"
        pi = cands[0]
        p = PRIMES[pi]
        x2 = pt[1] % p
        nv = nz = 0
        failures = []
        for r in recs:
            if r.get("zero") or not r.get("done"):
                continue
            N = {int(e): Fraction(v) for e, v in r["N"].items()}
            D = {int(e): Fraction(v) for e, v in r["D"].items()}
            dv = eval_mod(D, x2, p)
            if dv == 0:
                r["validated"] = False
                r["validation_error"] = "denominator vanishes at validation point"
                failures.append(r["fn"])
                print(f"[extract] VALIDATION FAILURE fn {r['fn']}: denominator "
                      f"vanishes at validation point (prime idx {pi})", file=sys.stderr)
                continue
            got = eval_mod(N, x2, p) * pow(dv, p - 2, p) % p
            if got != vals[r["fn"]] % p:
                r["validated"] = False
                r["validation_error"] = f"blackbox mismatch @prime{pi}"
                failures.append(r["fn"])
                print(f"[extract] VALIDATION FAILURE fn {r['fn']}: blackbox "
                      f"MISMATCH (prime idx {pi})", file=sys.stderr)
                continue
            r["validated"] = f"blackbox@prime{pi}"
            nv += 1
        print(f"[extract] validated {nv} done fns against validation.gz "
              f"(prime index {pi} = {p}); {nz} zero; {len(failures)} FAILED", file=sys.stderr)

    with open(a.out, "w") as out:
        for r in recs:
            out.write(json.dumps(r) + "\n")
    print(f"[extract] wrote {a.out} ({len(recs)} fns)", file=sys.stderr)
    if not a.no_validate and failures:
        print(f"[extract] EXIT NONZERO: {len(failures)} validation failures: "
              f"{failures[:50]}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()

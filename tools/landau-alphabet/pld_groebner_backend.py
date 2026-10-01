#!/usr/bin/env python3
r"""
pld_groebner_backend.py -- fast Groebner backend for the PLD / Landau Alphabet
per-face discriminant (Landau-letter) computation.

WHY THIS EXISTS
---------------
The per-face Landau letters are the discriminant of the face Symanzik
polynomial P(x; params):  the locus in kinematic-parameter space where the
Landau equations

      P = 0,    dP/dx_i = 0   (all Feynman params x_i),     x_i =/= 0  (torus)

have a common solution.  Eliminating the x_i gives an ideal in the kinematic
ring whose generators (factored) are the letters.

The stock `landau_alphabet.py` does this with a bounded sympy
resultant-chain.  That chain EXPLODES on heavy faces -- on the 3-loop K4
graph the leading codim-1 face is degree-267, and the four 5-edge "3m+2gamma"
faces time out (>900 s) -- because pairwise resultants blow up the
intermediate coefficients.  The honest fix is a real Groebner elimination
engine.

WHAT THIS DOES
--------------
Elimination ideal via **msolve** (Berthomieu-Eder-Safey El Din), the fastest
open F4/FGLM implementation over Q.  We avoid msolve's `-S` f4sat path (which
does not emit a printable elimination GB) and instead saturate the ideal
against the product of Feynman parameters by the classical RABINOWITSCH trick
(adjoin t with t*prod(x_i) - 1) and read off the elimination ideal with
`-e ELIM -g 2` (modular grevlex GB; eliminated vars dropped automatically).

The letters are the irreducible factors (in the kinematic ring) of the
elimination-ideal generators.  We factor with sympy AFTER msolve has done the
heavy elimination, so factoring sees only the small, fully-reduced result.

FALLBACK ORDER (each guarded, honest status reported):
  1. msolve  (system /usr/bin/msolve or msolve_jll artifact)   [PREFERRED]
  2. Singular  (system /usr/bin/Singular: eliminate() on the saturated ideal)
  3. sympy groebner (sp.groebner(..., order='lex'))            [last resort]

Both saturate the ideal against the Feynman-parameter torus, so x_i=0
spurious branches are removed -- matching the bounded resultant engine's
torus restriction.

PUBLIC API
----------
    face_letters(F, fvars, kinvars, *, char0=True, threads=8, timeout=900,
                 backend="auto") -> dict
        F        : sympy expr, the face F-polynomial (or U-polynomial)
        fvars    : list of sympy symbols, the Feynman parameters to eliminate
        kinvars  : list of sympy symbols, the kinematic ring to keep
        returns  : {"letters": [...sympy...], "letters_str": [...],
                    "backend": ..., "status": "resolved"|"timeout"|"error",
                    "secs": float, "elim_gens": [...], "note": ...}

CLI
---
    python pld_groebner_backend.py FACE_5EDGE_UF.json --out OUT.json \
           [--field s,t,m2] [--threads 8] [--timeout 900] [--backend auto]

reads the FACE_*_UF.json schema {face: {U, F, F_m2_1, edges, masses, ...}}
and computes letters for every face, dehomogenizing the chart automatically.
"""
from __future__ import annotations
import sympy as sp
import json, os, sys, time, shutil, subprocess, tempfile, argparse, glob

# ---------------------------------------------------------------------------
# locate msolve
# ---------------------------------------------------------------------------
def _find_msolve() -> str | None:
    m = shutil.which("msolve")
    if m:
        return m
    # msolve_jll artifact
    for p in glob.glob(os.path.expanduser(
            "~/.julia/artifacts/*/bin/msolve")):
        if os.access(p, os.X_OK):
            return p
    return None

_MSOLVE = _find_msolve()
_SINGULAR = shutil.which("Singular")


# ---------------------------------------------------------------------------
# sympy -> msolve text
# ---------------------------------------------------------------------------
def _ms_poly(expr) -> str:
    # msolve SILENTLY misparses the sympy sstr rational form "c*mono/den"
    # (probe: "36*a1^2/77-1" -> GB {36*a1^2-1}, denominator dropped, rc=0),
    # so coefficients are cleared to integers before writing: clear
    # denominators (unit scaling of each generator; ideal and variety
    # unchanged) so the written string is always integer-coefficient.
    e = sp.expand(expr)
    if e.free_symbols:
        try:
            e = sp.Poly(e, *sorted(e.free_symbols, key=str)).clear_denoms()[1].as_expr()
        except Exception:
            pass
    return sp.printing.str.sstr(sp.expand(e)).replace("**", "^").replace(" ", "")


def _write_msolve_input(path, allvars, polys):
    with open(path, "w") as f:
        f.write(",".join(str(v) for v in allvars) + "\n")
        f.write("0\n")                       # characteristic 0 (rational)
        f.write(",\n".join(_ms_poly(p) for p in polys) + "\n")


# ---------------------------------------------------------------------------
# parse msolve -g 2 GB output  ->  list of sympy polys in the kept vars
# ---------------------------------------------------------------------------
def _parse_msolve_gb(text, keptvars):
    # output: a header block of #-lines, then "[g1, \n g2, \n ... ]:"
    loc = {str(v): v for v in keptvars}
    # grab the bracketed body
    lb = text.find("[")
    rb = text.rfind("]")
    if lb < 0 or rb < 0 or rb <= lb:
        return []
    body = text[lb + 1: rb]
    gens = []
    for chunk in body.split(","):
        c = chunk.strip().strip("\n")
        if not c:
            continue
        c = c.replace("^", "**")
        try:
            gens.append(sp.sympify(c, locals=loc))
        except Exception:
            pass
    return gens


# ---------------------------------------------------------------------------
# backend: msolve  (Rabinowitsch saturation + elimination, modular grevlex GB)
# ---------------------------------------------------------------------------
def _msolve_elim(gens, fvars, kinvars, threads, timeout):
    if _MSOLVE is None:
        raise RuntimeError("msolve not found")
    t = sp.Symbol("rab_t")
    sat = t * sp.prod(fvars) - 1
    polys = list(gens) + [sat]
    # eliminated block FIRST: [rab_t, fvars...] then kinvars
    elimvars = [t] + list(fvars)
    allvars = elimvars + list(kinvars)
    nelim = len(elimvars)
    with tempfile.TemporaryDirectory(prefix="pld_ms_") as td:
        inp = os.path.join(td, "in.ms")
        out = os.path.join(td, "out.ms")
        _write_msolve_input(inp, allvars, polys)
        # Cap law: every msolve exec goes through the capped wrapper (the
        # dipstick member) — one cap door, never a raw exec.  MSOLVE_CAP_GB
        # defaults to 30 here; an explicit environment value wins, and the
        # wrapper's refusal (rc=2) surfaces as backend-unavailable, falling
        # through to the Singular leg.
        wrapper = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               os.pardir, "dipstick", "run_msolve_capped.sh")
        env = dict(os.environ)
        env.setdefault("MSOLVE_CAP_GB", "30")
        env["MSOLVE_BIN"] = _MSOLVE
        cmd = ["bash", wrapper, "-f", inp, "-e", str(nelim), "-g", "2",
               "-t", str(threads), "-o", out]
        subprocess.run(cmd, timeout=timeout, check=True, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        text = open(out).read() if os.path.exists(out) else ""
    return _parse_msolve_gb(text, kinvars)


# ---------------------------------------------------------------------------
# backend: Singular  (eliminate() on the saturated ideal)
# ---------------------------------------------------------------------------
def _singular_elim(gens, fvars, kinvars, timeout):
    if _SINGULAR is None:
        raise RuntimeError("Singular not found")
    t = sp.Symbol("rab_t")
    allvars = [t] + list(fvars) + list(kinvars)
    polys = list(gens) + [t * sp.prod(fvars) - 1]
    varstr = ",".join(str(v) for v in allvars)
    def sg(e):  # singular poly text
        return sp.printing.str.sstr(sp.expand(e)).replace("**", "^").replace(" ", "")
    elim_prod = "*".join(str(v) for v in ([t] + list(fvars)))
    idgens = ",".join(sg(p) for p in polys)
    script = f"""
ring R = 0,({varstr}),dp;
ideal I = {idgens};
ideal E = eliminate(I, {elim_prod});
print(E);
quit;
"""
    with tempfile.TemporaryDirectory(prefix="pld_sg_") as td:
        sp_path = os.path.join(td, "s.sing")
        open(sp_path, "w").write(script)
        res = subprocess.run([_SINGULAR, "-q", sp_path], timeout=timeout,
                             capture_output=True, text=True)
    loc = {str(v): v for v in kinvars}
    gens_out = []
    for line in res.stdout.splitlines():
        line = line.strip().rstrip(",")
        if not line or line.startswith("//"):
            continue
        try:
            e = sp.sympify(line.replace("^", "**"), locals=loc)
            if e.free_symbols and e.free_symbols <= set(kinvars):
                gens_out.append(e)
        except Exception:
            pass
    return gens_out


# ---------------------------------------------------------------------------
# backend: sympy groebner  (last resort)
# ---------------------------------------------------------------------------
def _sympy_elim(gens, fvars, kinvars, timeout):
    t = sp.Symbol("rab_t")
    polys = list(gens) + [t * sp.prod(fvars) - 1]
    order_vars = [t] + list(fvars) + list(kinvars)
    G = sp.groebner(polys, *order_vars, order="lex")
    out = []
    for g in G.exprs:
        fs = g.free_symbols
        if fs and fs <= set(kinvars):
            out.append(g)
    return out


# ---------------------------------------------------------------------------
# factor elimination-ideal generators into letters (kinematic ring only)
# ---------------------------------------------------------------------------
def _letters_from_gens(gens, kinvars):
    letters = []
    seen = []
    for g in gens:
        g = sp.expand(g)
        if g == 0 or g.is_number:
            continue
        for fpow in sp.Mul.make_args(sp.factor(g)):
            base = fpow.as_base_exp()[0]
            base = sp.factor(base)
            if base.is_number:
                continue
            if not (base.free_symbols & set(kinvars)):
                continue
            # strip numeric content
            c, _ = base.as_coeff_Mul()
            if c.is_number and c not in (1, 0):
                base = sp.factor(sp.expand(base / c))
            # dedupe up to scalar
            dup = False
            for s in seen:
                r = sp.cancel(base / s)
                if r.is_number:
                    dup = True
                    break
            if not dup:
                seen.append(base)
                letters.append(base)
    return letters


# ---------------------------------------------------------------------------
# PUBLIC: face_letters
# ---------------------------------------------------------------------------
def face_letters(F, fvars, kinvars, *, threads=8, timeout=900,
                 backend="auto", dehom=True):
    """Discriminant letters of face polynomial F by Groebner elimination."""
    F = sp.expand(F)
    fvars = [v for v in fvars if v in F.free_symbols]
    if not fvars:
        return {"letters": [], "letters_str": [], "backend": "none",
                "status": "no-free-params", "secs": 0.0, "elim_gens": []}
    # dehomogenize: set last Feynman param chart -> 1 (projective scaling
    # invariance of the Landau equations; matches the bounded engine).
    elim = fvars
    if dehom and len(fvars) >= 2:
        chart = fvars[-1]
        F = sp.expand(F.subs(chart, 1))
        elim = fvars[:-1]
    gens = [F] + [sp.expand(sp.diff(F, v)) for v in elim]
    order = (["msolve", "singular", "sympy"] if backend == "auto"
             else [backend])
    t0 = time.time()
    last_err = None
    for be in order:
        try:
            if be == "msolve":
                eg = _msolve_elim(gens, elim, kinvars, threads, timeout)
            elif be == "singular":
                eg = _singular_elim(gens, elim, kinvars, timeout)
            elif be == "sympy":
                eg = _sympy_elim(gens, elim, kinvars, timeout)
            else:
                continue
            letters = _letters_from_gens(eg, kinvars)
            return {"letters": letters,
                    "letters_str": [sp.printing.str.sstr(l) for l in letters],
                    "backend": be, "status": "resolved",
                    "secs": round(time.time() - t0, 2),
                    "elim_gens": [sp.printing.str.sstr(g) for g in eg]}
        except subprocess.TimeoutExpired:
            last_err = f"{be}:timeout({timeout}s)"
            return {"letters": [], "letters_str": [], "backend": be,
                    "status": "timeout", "secs": round(time.time() - t0, 2),
                    "elim_gens": [], "note": last_err}
        except Exception as e:
            last_err = f"{be}:{type(e).__name__}:{e}"
            continue
    return {"letters": [], "letters_str": [], "backend": "none",
            "status": "error", "secs": round(time.time() - t0, 2),
            "elim_gens": [], "note": last_err}


# ---------------------------------------------------------------------------
# CLI: run over a FACE_*_UF.json
# ---------------------------------------------------------------------------
def _run_facefile(path, fieldvars, threads, timeout, backend, out_path):
    data = json.load(open(path))
    kinvars = [sp.Symbol(s, real=True) for s in fieldvars]
    loc = {s: v for s, v in zip(fieldvars, kinvars)}
    results = {}
    for face, fd in data.items():
        # collect Feynman params x1.. from the polynomial
        Fexpr = sp.sympify(str(fd["F_m2_1"] if "F_m2_1" in fd else fd["F"]),
                           locals=loc)
        fvars = sorted([s for s in Fexpr.free_symbols
                        if str(s).startswith("x")], key=lambda z: str(z))
        sys.stderr.write(f"[{face}] F-elim over {fvars} keep {kinvars} ...\n")
        rF = face_letters(Fexpr, fvars, kinvars, threads=threads,
                          timeout=timeout, backend=backend)
        sys.stderr.write(f"   F: {rF['status']} via {rF['backend']} "
                         f"({rF['secs']}s) letters={rF['letters_str']}\n")
        rec = {"F": rF}
        # also the U-polynomial (second-type / at-infinity) if present
        if "U" in fd:
            Uexpr = sp.sympify(str(fd["U"]), locals=loc)
            uv = sorted([s for s in Uexpr.free_symbols
                         if str(s).startswith("x")], key=lambda z: str(z))
            rU = face_letters(Uexpr, uv, kinvars, threads=threads,
                              timeout=timeout, backend=backend)
            rec["U"] = rU
            sys.stderr.write(f"   U: {rU['status']} via {rU['backend']} "
                             f"({rU['secs']}s) letters={rU['letters_str']}\n")
        results[face] = rec
        if out_path:
            json.dump(results, open(out_path, "w"), indent=2, default=str)
    # union alphabet
    alpha = []
    for face, rec in results.items():
        for k in ("F", "U"):
            for l in rec.get(k, {}).get("letters_str", []):
                if l not in alpha:
                    alpha.append(l)
    summary = {"file": path, "field": fieldvars, "backend_requested": backend,
               "msolve": _MSOLVE, "singular": _SINGULAR,
               "alphabet": alpha, "faces": results}
    if out_path:
        json.dump(summary, open(out_path, "w"), indent=2, default=str)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description="Fast Groebner backend for PLD/Landau face letters.")
    ap.add_argument("facefile", help="FACE_*_UF.json")
    ap.add_argument("--out", default=None)
    ap.add_argument("--field", default="s,t,m2", help="kinematic ring, comma-separated")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--backend", default="auto", choices=["auto", "msolve", "singular", "sympy"])
    args = ap.parse_args(argv)
    field = [s.strip() for s in args.field.split(",") if s.strip()]
    out = args.out or os.path.splitext(args.facefile)[0] + "_GROEBNER_LETTERS.json"
    s = _run_facefile(args.facefile, field, args.threads, args.timeout,
                      args.backend, out)
    print(f"\nmsolve={_MSOLVE}  singular={_SINGULAR}")
    print(f"alphabet ({len(s['alphabet'])} letters): {s['alphabet']}")
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

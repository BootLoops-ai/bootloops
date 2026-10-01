#!/usr/bin/env python3
r"""
sampler.py — the Wayfinder "dlog-connection-once" sampler.

THE INSIGHT
===========
In a canonical (UT) basis the differential-equation matrix is

        d J  =  ε · ( Σ_a  A_a · d log a )  · J ,                     (★)

with the {A_a} CONSTANT ℚ-matrices and {a} the symbol alphabet.  Symbolic IBP
recovers (★) by reconstructing every entry of A(x,ε) as a multivariate rational
function (FireFly over ≳9 variables — measured ~30 min/node on a
two-loop five-point family).  But (★) has
only  n_letters × M²  RATIONAL unknowns.  So instead:

    1.  alphabet {a₁,…,aₙ}                    ← Landau Alphabet (already built)
    2.  K cheap NUMERIC oracle samples        ← amflow-cpp solve_integrals
    3.  per-point, per-x_i, per-ε^k slice of (★) is a LINEAR equation in the
        rational entries of {A_a}             ← stack → least squares
    4.  PSLQ/LLL-rationalize each entry        ← gatekeeper.pslq_with_lll
    5.  verify on held-out points

never running symbolic FireFly.

NUMERIC-DE BACKEND
==================
Two backends:

  finite_diff   — central difference of `mode:"solve_integrals"` Laurent
                  values at  p ± h·e_i.  (1+2·|kinvars|) amflow calls per
                  point, ~⅔·goal_digits effective derivative precision.

  amflow_diffeq — `mode:"diffeq"` (src/api/run_json.cpp), which now wraps
                  `ibp::diffeq()` with the upstream-MMA momentum-derivative
                  chain rule (`LIBPDerivivative`) ported into
                  `src/ibp/libp_deriv.cpp` so replacement-defined kinematic
                  invariants (s, t, sᵢⱼ) are supported.  2 amflow calls per
                  point (1× solve_integrals + 1× diffeq), full goal_digits
                  derivative precision (exact-in-ℚ DE matrix · numeric I).

PUBLIC API
==========
    Rotation(spec)                            — UT rotation (1/LS · ε^p, per master)
    sample_de(target, points, ...)            → list[DESample]   (finite-diff backend)
    fit_dlog_connection(alphabet, ...)        → dict (A_a + honesty fields)
    verify_connection(conn, ...)              → min held-out digits
    connection_once(target_json_path, ...)    — one command, writes connection.json

OUTPUT  connection.json
=======================
    {
      "alphabet": [...], "kinvars": [...], "masters": [...],
      "A": { "<letter>": [[i,j,"p/q"], ...], ... },          # sparse
      "A_dense": { "<letter>": [["p/q",...], ...], ... },
      "rotation": {...},
      "honesty": { n_points_used, cond_number, min_heldout_digits,
                   failed_to_rationalize, ut_check, backend, ... }
    }

This is the schema `tools/ansatzer/` consumes via `--connection`
(see `ansatzer.py:load_connection`) and the schema `transport`
will read.

PROVENANCE
==========
    alphabet                       : tools/landau-alphabet/landau_alphabet.py
    UT rotation (LS, default)      : tools/canonical_form/canonical_form.py
    PSLQ/LLL rationalize, lstsq    : tools/gatekeeper/heldout_cv.py
    amflow-cpp solve_integrals     : Software/AMFlow-port/.../run_json.cpp:718
    libp_deriv scope note          : Software/AMFlow-port/.../libp_deriv.hpp:26
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from fractions import Fraction

import mpmath as mp
import sympy as sp

HERE  = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(HERE)
sys.path.insert(0, TOOLS)
sys.path.insert(0, os.path.join(TOOLS, "gatekeeper"))
sys.path.insert(0, os.path.join(TOOLS, "canonical_form"))
# public repo layout: canonical_form ships as a member of Counterweight
sys.path.insert(0, os.path.join(TOOLS, "counterweight", "canonical_form"))
sys.path.insert(0, os.path.join(TOOLS, "landau-alphabet"))

try:
    from heldout_cv import (_pslq_rational, _lstsq, _cond_svd,   # type: ignore
                            auto_condition, pslq_with_lll)
    import canonical_form as cf                                   # type: ignore
except ImportError as _e:
    raise ImportError(
        "wayfinder.sampler needs `heldout_cv` (gatekeeper) and "
        "`canonical_form` (ships with the Counterweight tool) "
        "importable — add their directories to PYTHONPATH or place them as "
        "siblings under tools/") from _e

AMFLOW_CLI = os.environ.get(
    "AMFLOW_CLI",
    "amflow_cli",  # set AMFLOW_CLI to your amflow-cpp build (the sibling repository amflow-cpp)
)
DEFAULT_IBP_CACHE = os.environ.get("AMFLOW_IBP_CACHE", "")

_ARB = re.compile(r"^\s*\[\s*([^\s\]]*)\s*\+/-\s*([^\s\]]+)\s*\]\s*$")


def _sp_to_mpc(e, dps: int) -> mp.mpc:
    """sympy number → mp.mpc at `dps` digits, EXACTLY for Rationals (never via
    nsimplify, which rounds near-integers and would destroy finite-diff
    shifts)."""
    if e == 0:
        return mp.mpc(0)
    if getattr(e, "is_Rational", False):
        return mp.mpc(mp.mpf(int(e.p)) / mp.mpf(int(e.q)))
    # algebraic / float: sympy high-precision N, then split re/im
    with mp.workdps(dps + 10):
        v = sp.N(e, dps + 10)
        return mp.mpc(mp.mpmathify(sp.re(v)), mp.mpmathify(sp.im(v)))


def _frac(v) -> sp.Rational:
    f = Fraction(str(v))
    return sp.Rational(f.numerator, f.denominator)


# ===========================================================================
#  Rotation — diagonal  J_m = c_m(x, ε) · I_m^{raw},   c_m polynomial in ε.
#
#  The general canonical rotation is  J = T(x, ε) · I  with T rational in ε.
#  In practice the per-master factor is  c_m(x) · (polynomial in ε): the
#  ε^{p_m} weight shift, times any (1-2ε)-type factors that absorb sub-sector
#  Γ-ratios (the standard Henn move for massless bubbles).  We store c_m as a
#  sympy expression in {kinvars ∪ {eps}}; `eval_eps_poly` evaluates the
#  kinematics and returns the ε-Laurent coefficients {j: mpc}.
# ===========================================================================
EPS = sp.Symbol("eps")


@dataclass
class Rotation:
    masters:   list[str]
    coeff:     dict[str, sp.Expr]      # tag → symbolic c_m(x, ε)
    elliptic:  set                      # tags whose LS is elliptic / unknown
    kinvars:   list[str]

    @classmethod
    def from_spec(cls, spec: dict, masters: list[dict], kinvars: list[str]):
        """Build from a target['rotation'] block:
            {tag: {"coeff": "<sympy str in kinvars and eps>"}}
        Legacy keys `eps_power` (multiply by ε^p) and `eps_poly` (list of
        ε-coeffs from ε⁰ upward) are also accepted."""
        L = spec.get("_default_eps_power", 1)
        c, ell = {}, set()
        for m in masters:
            t = m["tag"]
            ent = spec.get(t, {})
            ce = sp.sympify(str(ent.get("coeff", 1)), locals={"eps": EPS})
            if "eps_poly" in ent:
                ce = ce * sum(sp.Rational(str(a)) * EPS**j
                              for j, a in enumerate(ent["eps_poly"]))
            elif EPS not in ce.free_symbols:
                ce = ce * EPS ** int(ent.get("eps_power", L))
            c[t] = sp.expand(ce)
            if ent.get("elliptic"):
                ell.add(t)
        return cls([m["tag"] for m in masters], c, ell, list(kinvars))

    @classmethod
    def from_canonical_form(cls, family_spec: dict, masters: list[dict],
                            kinvars: list[str]):
        """Default: diagonal LS rotation from tools/canonical_form/.  Elliptic
        sectors are FLAGGED (the dlog ansatz only covers the polylog block)."""
        fam = cf.Family(family_spec)
        LS = cf.leading_singularities(fam, masters)
        c, ell = {}, set()
        for m in masters:
            t = m["tag"]; r = LS.get(t)
            if r is None or r.kind not in ("rational", "sqrt") or r.expr is None:
                c[t] = EPS ** fam.L
                ell.add(t)
            else:
                c[t] = sp.simplify(1 / r.expr) * EPS ** fam.L
        return cls([m["tag"] for m in masters], c, ell, list(kinvars))

    def eval_eps_poly(self, tag: str, point: dict, dps: int) -> dict[int, mp.mpc]:
        """Return {j: c_m^{(j)}(point)} so that c_m(point, ε) = Σ_j c^{(j)} ε^j."""
        with mp.workdps(dps):
            subs = {sp.Symbol(k): _frac(v) for k, v in point.items()}
            ce = sp.expand(self.coeff[tag].subs(subs))
            poly = sp.Poly(ce, EPS)
            return {int(j): _sp_to_mpc(a, dps) for (j,), a in poly.terms()}

    def to_json(self):
        return {t: {"coeff": str(self.coeff[t]),
                    "elliptic": (t in self.elliptic)}
                for t in self.masters}


# ===========================================================================
#  DESample — one kinematic point's rotated Laurent data + finite-diff ∂J
# ===========================================================================
@dataclass
class DESample:
    point: dict[str, str]                       # exact-rational strings
    orders: list[int]                            # ε-orders present (rotated)
    J:  dict[int, list]                          # order → length-M list[mpc]
    dJ: dict[str, dict[int, list]]               # var → order → length-M list[mpc]
    deriv_digits: int                            # effective digits in dJ
    basis_diag: dict = field(default_factory=dict)  # Kira↔user embedding diag

    def J_mat(self, k):  return mp.matrix([[v] for v in self.J[k]])
    def dJ_mat(self, v, k): return mp.matrix([[x] for x in self.dJ[v][k]])


# ===========================================================================
#  amflow-cpp wrapper  (mode: solve_integrals)
# ===========================================================================
def _build_amflow_input(target: dict, masters: list[dict], pt: dict[str, str],
                        goal: int, K: int, work_dir: str) -> dict:
    amf_opts = copy.deepcopy(target.get("amf_options") or {})
    bb = amf_opts.setdefault("blackbox", {})
    nv = bb.setdefault("numeric_values", {})
    for k, v in pt.items():
        nv[str(k)] = str(v)
    if DEFAULT_IBP_CACHE and "ibp_cache_dir" not in bb:
        bb["ibp_cache_dir"] = DEFAULT_IBP_CACHE
    return {
        "mode": "solve_integrals",
        "options": target.get("amflow_engine_options",
                              {"silent_mode": True, "chop_pre": 20}),
        "family": target["family"],
        "integrals": [{"indices": m["indices"]} for m in masters],
        "goal_digits": goal,
        "target_eps_power": K,
        "work_dir": work_dir,
        "amf_options": amf_opts,
    }


def _run_amflow(cfg: dict, work_root: str, tag: str, timeout_s: int) -> dict:
    in_p  = os.path.join(work_root, "in",  tag + ".json")
    out_p = os.path.join(work_root, "out", tag + ".json")
    log_p = os.path.join(work_root, "log", tag + ".log")
    for d in (os.path.dirname(in_p), os.path.dirname(out_p),
              os.path.dirname(log_p)):
        os.makedirs(d, exist_ok=True)
    # cache key: physics inputs only (not work_dir, which is per-run)
    key_cfg = {k: v for k, v in cfg.items() if k != "work_dir"}
    sha = hashlib.sha1(json.dumps(key_cfg, sort_keys=True).encode()).hexdigest()
    if os.path.exists(out_p):
        try:
            prev = json.load(open(out_p))
            if prev.get("_input_sha") == sha:
                return prev
        except Exception:
            pass
    with open(in_p, "w") as f:
        json.dump(cfg, f, indent=2)
    with open(log_p, "a") as lf:
        rc = subprocess.run([AMFLOW_CLI, in_p, out_p],
                            stdout=lf, stderr=subprocess.STDOUT,
                            timeout=timeout_s).returncode
    if rc != 0:
        raise RuntimeError(f"amflow_cli exit {rc} (see {log_p})")
    res = json.load(open(out_p))
    res["_input_sha"] = sha
    json.dump(res, open(out_p, "w"))
    return res


def _parse_ball(s: str) -> mp.mpf:
    s = str(s).strip()
    m = _ARB.match(s)
    return mp.mpf(m.group(1) or "0") if m else mp.mpf(s)


def _laurent(amf_out: dict, masters: list[dict]) -> dict[str, dict[int, mp.mpc]]:
    out = {}
    for ent, m in zip(amf_out["result"], masters):
        co = {}
        for c in ent.get("coefficients", []):
            o = int(c["order"])
            co[o] = mp.mpc(_parse_ball(c["value"]["re"]),
                           _parse_ball(c["value"]["im"]))
        out[m["tag"]] = co
    return out


def _rotate(laurent: dict, rotation: Rotation, point: dict, dps: int,
            order_max: int) -> dict[int, list]:
    """J_m^{(k)} = Σ_j c_m^{(j)}(point) · I_m^{(k-j)}   (ε-polynomial conv).
    Returns {k: [J_m^{(k)}]_m for k=0..order_max}."""
    out: dict[int, list] = {}
    polys = {t: rotation.eval_eps_poly(t, point, dps) for t in rotation.masters}
    for k in range(0, order_max + 1):
        row = []
        for t in rotation.masters:
            v = mp.mpc(0)
            for j, cj in polys[t].items():
                v += cj * laurent[t].get(k - j, mp.mpc(0))
            row.append(v)
        out[k] = row
    return out


# ===========================================================================
#  sample_de — finite-difference backend
# ===========================================================================
def _shift_point(pt: dict, var: str, h: Fraction) -> dict[str, str]:
    out = dict(pt)
    out[var] = str(Fraction(str(pt[var])) + h)
    return out


def sample_de(target: dict, points: list[dict], *, rotation: Rotation,
              kinvars: list[str], h: Fraction = Fraction(1, 10**18),
              goal_digits: int = 50, target_eps_power: int = 4,
              work_root: str | None = None, parallel: int = 4,
              timeout_s: int = 1800,
              backend: str = "finite_diff") -> list[DESample]:
    """Sample the rotated DE at each point.

    backend = "finite_diff":
        For each point p run amflow-cpp solve_integrals at p and at p ± h·e_i
        for every kinematic variable x_i, central-difference the rotated
        Laurent coefficients.  Effective derivative digits ≈
        min(2·log10(1/h), goal_digits - log10(1/h)).

    backend = "amflow_diffeq":
        For each point p run amflow-cpp `mode:"diffeq"` (exact-in-ℚ DE matrix
        in `d`, via the LIBPDerivivative chain rule + Kira IBP) and one
        `mode:"solve_integrals"` for the raw Laurent values; assemble
        ∂I/∂x_i = M_{x_i}(d=4-2ε)·I as an exact ε-Laurent product, then
        rotate.  2 amflow calls/point, derivative digits = goal_digits.
    """
    if backend == "amflow_diffeq":
        return _sample_de_amflow_diffeq(
            target, points, rotation=rotation, kinvars=kinvars,
            goal_digits=goal_digits, target_eps_power=target_eps_power,
            work_root=work_root, parallel=parallel, timeout_s=timeout_s)

    masters = target["masters"]
    work_root = work_root or tempfile.mkdtemp(prefix="wayfinder_samp_")
    os.makedirs(work_root, exist_ok=True)
    dps = goal_digits + 20
    mp.mp.dps = dps

    # Build the call list:  for each point p, calls = {center, ±h per var}
    calls: list[tuple[str, dict]] = []
    for ip, p in enumerate(points):
        p = {k: str(v) for k, v in p.items()}
        calls.append((f"p{ip}_c", p))
        for v in kinvars:
            calls.append((f"p{ip}_{v}_p", _shift_point(p, v, +h)))
            calls.append((f"p{ip}_{v}_m", _shift_point(p, v, -h)))

    def _one(tag, pt):
        wd = os.path.join(work_root, "work", tag)
        cfg = _build_amflow_input(target, masters, pt, goal_digits,
                                  target_eps_power, wd)
        out = _run_amflow(cfg, work_root, tag, timeout_s)
        return tag, pt, _laurent(out, masters)

    raw: dict[str, tuple[dict, dict]] = {}
    t0 = time.time()
    if parallel <= 1:
        for tag, pt in calls:
            raw[tag] = _one(tag, pt)[1:]
    else:
        with ThreadPoolExecutor(max_workers=parallel) as ex:
            futs = {ex.submit(_one, tag, pt): tag for tag, pt in calls}
            for f in as_completed(futs):
                tag, pt, L = f.result()
                raw[tag] = (pt, L)
    t_sample = time.time() - t0

    # Assemble DESamples (rotate, then central-diff the ROTATED Laurent)
    hf = mp.mpf(h.numerator) / mp.mpf(h.denominator)
    deriv_digits = int(min(2 * (-mp.log10(hf)),
                           goal_digits - (-mp.log10(hf))) - 2)
    samples: list[DESample] = []
    K = target_eps_power
    for ip, p0 in enumerate(points):
        pc, Lc = raw[f"p{ip}_c"]
        Jc = _rotate(Lc, rotation, pc, dps, K)
        dJ: dict[str, dict[int, list]] = {}
        for v in kinvars:
            pp, Lp = raw[f"p{ip}_{v}_p"]
            pm, Lm = raw[f"p{ip}_{v}_m"]
            Jp = _rotate(Lp, rotation, pp, dps, K)
            Jm = _rotate(Lm, rotation, pm, dps, K)
            dJ[v] = {k: [(Jp[k][m] - Jm[k][m]) / (2 * hf)
                         for m in range(len(rotation.masters))]
                     for k in range(0, K + 1)}
        samples.append(DESample(point={k: str(v) for k, v in p0.items()},
                                orders=list(range(K + 1)),
                                J=Jc, dJ=dJ, deriv_digits=deriv_digits))
    sys.stderr.write(f"[sample_de] {len(calls)} amflow calls in "
                     f"{t_sample:.1f}s  (~{t_sample/max(len(calls),1):.1f}s/call), "
                     f"∂-digits≈{deriv_digits}\n")
    return samples


# ---------------------------------------------------------------------------
#  amflow_diffeq backend
# ---------------------------------------------------------------------------
_D = sp.Symbol("d")


def _build_diffeq_input(target: dict, masters: list[dict], pt: dict[str, str],
                        kinvars: list[str], work_dir: str) -> dict:
    amf_opts = copy.deepcopy(target.get("amf_options") or {})
    bb = amf_opts.setdefault("blackbox", {})
    nv = bb.setdefault("numeric_values", {})
    for k, v in pt.items():
        nv[str(k)] = str(v)
    bb["work_dir"] = work_dir
    if DEFAULT_IBP_CACHE and "ibp_cache_dir" not in bb:
        bb["ibp_cache_dir"] = DEFAULT_IBP_CACHE
    return {
        "mode": "diffeq",
        "family": target["family"],
        "integrals": [{"indices": m["indices"]} for m in masters],
        "variables": list(kinvars),
        "amf_options": amf_opts,
    }


def _de_entry_laurent(expr_str: str, dps: int, K: int) -> dict[int, mp.mpc]:
    """Parse an `ibp::diffeq` Mfrac entry (rational in `d`, possibly `eta`)
    and return its ε-Laurent coefficients {j: c_j} for j ∈ [j₀..K] under
    d → 4-2ε.  `eta` is mapped to 0 (the family is not η-injected; it only
    appears because the reduction context names it)."""
    if expr_str in ("0", "(0)"):
        return {}
    e = sp.sympify(expr_str.replace("^", "**"),
                   locals={"d": _D, "eta": sp.Integer(0)})
    e = sp.together(e.subs(_D, 4 - 2 * EPS))
    # Laurent series about ε=0.  K+2 terms is plenty (entries are O(ε⁰)
    # or at worst a single 1/ε from a degenerate Gram — never seen here).
    ser = sp.series(e, EPS, 0, K + 2).removeO()
    poly = sp.Poly(sp.expand(ser), EPS)
    out: dict[int, mp.mpc] = {}
    with mp.workdps(dps):
        for (j,), c in poly.terms():
            out[int(j)] = _sp_to_mpc(c, dps)
    return out


def _sample_de_amflow_diffeq(target: dict, points: list[dict], *,
                             rotation: Rotation, kinvars: list[str],
                             goal_digits: int, target_eps_power: int,
                             work_root: str | None, parallel: int,
                             timeout_s: int) -> list[DESample]:
    masters = target["masters"]
    work_root = work_root or tempfile.mkdtemp(prefix="wayfinder_samp_de_")
    os.makedirs(work_root, exist_ok=True)
    dps = goal_digits + 20
    mp.mp.dps = dps
    K = target_eps_power

    # Build the call list:  for each point p, calls = {solve_integrals, diffeq}
    calls: list[tuple[str, str, dict]] = []
    for ip, p in enumerate(points):
        p = {k: str(v) for k, v in p.items()}
        wd_si = os.path.join(work_root, "work", f"p{ip}_si")
        wd_de = os.path.join(work_root, "work", f"p{ip}_de")
        calls.append((f"p{ip}_si", "si",
                      _build_amflow_input(target, masters, p, goal_digits, K, wd_si)))
        calls.append((f"p{ip}_de", "de",
                      _build_diffeq_input(target, masters, p, kinvars, wd_de)))

    def _one(tag, kind, cfg):
        return tag, kind, _run_amflow(cfg, work_root, tag, timeout_s)

    raw: dict[str, dict] = {}
    t0 = time.time()
    if parallel <= 1:
        for tag, kind, cfg in calls:
            raw[tag] = _one(tag, kind, cfg)[2]
    else:
        with ThreadPoolExecutor(max_workers=parallel) as ex:
            futs = {ex.submit(_one, *c): c[0] for c in calls}
            for f in as_completed(futs):
                tag, _kind, out = f.result()
                raw[tag] = out
    t_sample = time.time() - t0

    # Assemble DESamples.
    idx_key = lambda lst: tuple(int(x) for x in lst)
    user_idx = {idx_key(m["indices"]): m["tag"] for m in masters}
    samples: list[DESample] = []
    for ip, p0 in enumerate(points):
        pc = {k: str(v) for k, v in p0.items()}
        L_raw = _laurent(raw[f"p{ip}_si"], masters)          # tag → {ord: mpc}
        de = raw[f"p{ip}_de"]
        kira_masters = de["masters"]
        Mk = len(kira_masters)
        # ---- Projection / embedding of the Kira DE basis into the user basis.
        #
        # The diffeq matrix M_v is (Mk × Mk) over Kira's OWN master basis,
        # which is the genuine, IBP-irreducible basis of the family.  The
        # user's target['masters'] is an INDEPENDENT list of integrals we asked
        # solve_integrals to evaluate.  These two sets are NOT size-matched in
        # general:
        #
        #   * Kira ⊋ user  (measured case: Kira=218, a stale target supplied 27):
        #       Kira couples each row to columns we never solved for.  Those
        #       columns' values are unavailable, so their contribution to
        #       ∂I/∂v cannot be assembled — we DROP them and flag the affected
        #       rows as having an incomplete (sub-basis) derivative rather than
        #       crashing (the earlier implementation crashed here with a
        #       RuntimeError).
        #   * user ⊋ Kira  (some user masters are reducible / not in the Kira
        #       basis): those user masters have no DE row; handled below by
        #       `kira_row_of.get(tg) is None`.
        #
        # The correspondence is by Baikov index vector (the sector mask + dot
        # powers), which is the canonical Kira↔user master identity.  We embed
        # by index match and synthesise placeholder tags (`__kira<c>`) for Kira
        # masters absent from the user basis, so all Mk columns stay addressable
        # while `solved_cols` records which actually have a solved value.
        kira_tags: list[str] = []          # length Mk: user tag or "__kira<c>"
        solved_cols: list[bool] = []       # length Mk: do we have L_raw value?
        unsolved_kira: list[list[int]] = []
        for c, km in enumerate(kira_masters):
            t = user_idx.get(idx_key(km["indices"]))
            if t is None:
                kira_tags.append(f"__kira{c}")
                solved_cols.append(False)
                unsolved_kira.append([int(x) for x in km["indices"]])
            else:
                kira_tags.append(t)
                solved_cols.append(True)
        kira_tag_set = set(t for t, ok in zip(kira_tags, solved_cols) if ok)
        missing = [m["tag"] for m in masters if m["tag"] not in kira_tag_set]
        if missing:
            sys.stderr.write(f"[sample_de/amflow_diffeq] WARN: {len(missing)} "
                             f"user master(s) are reducible (not in the Kira "
                             f"basis); their rows are filled by direct rotation "
                             f"of solve_integrals + zero coupling. "
                             f"First few: {missing[:8]}\n")
        if unsolved_kira:
            sys.stderr.write(
                f"[sample_de/amflow_diffeq] NOTE: Kira DE basis has {Mk} "
                f"masters but target['masters'] supplied {len(masters)}; "
                f"{len(unsolved_kira)} Kira column(s) have no solved value and "
                f"are PROJECTED OUT of the derivative assembly.  Rows that "
                f"couple to them are flagged incomplete (see "
                f"DESample.basis_diag).  For a FULL, sound connection the "
                f"target must supply the complete Kira basis (here {Mk} "
                f"masters) in the target's own conn_target json.\n")

        # ∂I_r/∂v |_{kira order} = Σ_c M_v[r][c](ε) ⊗_ε I_c  (Laurent conv.)
        # then rotate per user master:  J_m = c_m(x,ε)·I_m,
        #   ∂J_m/∂v = (∂_v c_m)·I_m + c_m·(∂_v I_m).
        # We need ∂_v I_m (the m-th user master's raw derivative), which is
        # row r=kira_tags.index(m) of M_v · I_kira.
        #
        # Laurent bookkeeping: solve_integrals returns I^{(k)} for k from
        # leading_order through K.  M_v entries are O(ε⁰) (rational in d),
        # so the convolution preserves leading orders.  We need rotated
        # orders 0..K, and rotation.coeff[m] is a polynomial in ε starting
        # at ε^{≥1} typically, so I^{(k)} for k down to (−max rotation deg)
        # are consumed — exactly the orders solve_integrals produces.
        polys = {t: rotation.eval_eps_poly(t, pc, dps) for t in rotation.masters}
        # ∂_v c_m(x,ε): differentiate the symbolic coeff, evaluate kinematics.
        dpolys = {v: {} for v in kinvars}
        subs = {sp.Symbol(k): _frac(vv) for k, vv in pc.items()}
        for tg in rotation.masters:
            ce = rotation.coeff[tg]
            for v in kinvars:
                dc = sp.expand(sp.diff(ce, sp.Symbol(v)).subs(subs))
                if dc == 0:
                    dpolys[v][tg] = {}
                else:
                    pl = sp.Poly(dc, EPS)
                    dpolys[v][tg] = {int(j): _sp_to_mpc(a, dps)
                                     for (j,), a in pl.terms()}

        # Pre-expand DE entries once per (v, r, c).
        de_laur = {v: [[_de_entry_laurent(de["diffeq"][v][r][c], dps, K + 4)
                        for c in range(Mk)] for r in range(Mk)]
                   for v in kinvars}

        # Raw ∂I_m/∂v (per user tag) as Laurent dict {ord: mpc}.  Only the
        # FIRST occurrence of a user tag in the Kira basis defines its row
        # (indices are unique, so this is a 1-1 embedding for solved columns).
        kira_row_of: dict[str, int] = {}
        for r, t in enumerate(kira_tags):
            kira_row_of.setdefault(t, r)
        ord_hi = K
        # Per-user-master row diagnostic: did we drop any live coupling because
        # the daughter column was an unsolved (projected-out) Kira master?
        incomplete_rows: dict[str, int] = {}
        dI: dict[str, dict[str, dict[int, mp.mpc]]] = {v: {} for v in kinvars}
        for tg in rotation.masters:
            r = kira_row_of.get(tg)
            for v in kinvars:
                acc: dict[int, mp.mpc] = {}
                if r is not None:
                    for c in range(Mk):
                        Mvc = de_laur[v][r][c]
                        if not Mvc:
                            continue
                        if not solved_cols[c]:
                            # Live coupling to an unsolved Kira master: its
                            # value is unavailable, so this term is dropped.
                            # Record it so the row is not silently wrong.
                            incomplete_rows[tg] = incomplete_rows.get(tg, 0) + 1
                            continue
                        Ic = L_raw[kira_tags[c]]
                        for jM, aM in Mvc.items():
                            for jI, aI in Ic.items():
                                k = jM + jI
                                if k > ord_hi:
                                    continue
                                acc[k] = acc.get(k, mp.mpc(0)) + aM * aI
                dI[v][tg] = acc
        if incomplete_rows:
            top = sorted(incomplete_rows.items(), key=lambda kv: -kv[1])[:8]
            sys.stderr.write(
                f"[sample_de/amflow_diffeq] WARN: {len(incomplete_rows)} "
                f"master row(s) couple to projected-out (unsolved) Kira "
                f"daughters; their ∂I/∂v is a SUB-BASIS derivative.  "
                f"Worst rows (dropped-couplings): {top}\n")

        # Rotated J and ∂J (orders 0..K).
        Jc = _rotate(L_raw, rotation, pc, dps, K)
        dJ: dict[str, dict[int, list]] = {}
        for v in kinvars:
            dJ[v] = {}
            for k in range(0, K + 1):
                row = []
                for tg in rotation.masters:
                    val = mp.mpc(0)
                    # c_m ⊗ ∂I_m
                    for j, cj in polys[tg].items():
                        val += cj * dI[v][tg].get(k - j, mp.mpc(0))
                    # (∂_v c_m) ⊗ I_m
                    for j, cj in dpolys[v][tg].items():
                        val += cj * L_raw[tg].get(k - j, mp.mpc(0))
                    row.append(val)
                dJ[v][k] = row
        samples.append(DESample(
            point=pc, orders=list(range(K + 1)),
            J=Jc, dJ=dJ, deriv_digits=goal_digits,
            basis_diag={
                "kira_basis_size": Mk,
                "user_basis_size": len(masters),
                "projected_out_kira_cols": len(unsolved_kira),
                "reducible_user_masters": missing,
                "incomplete_rows": dict(incomplete_rows),
                "full_basis": (not unsolved_kira and not incomplete_rows),
            }))

    sys.stderr.write(f"[sample_de/amflow_diffeq] {len(calls)} amflow calls in "
                     f"{t_sample:.1f}s  (~{t_sample/max(len(calls),1):.1f}s/call), "
                     f"∂-digits={goal_digits} (exact DE)\n")
    return samples


# ===========================================================================
#  fit_dlog_connection — row-by-row linear system + PSLQ-rationalize
# ===========================================================================
def _dlog_grad(alphabet: list, kinvars: list[str], point: dict,
               dps: int) -> dict:
    """Return G[a][v] = ∂log(a)/∂x_v |_{point}  as mpc."""
    subs = {sp.Symbol(k): _frac(v) for k, v in point.items()}
    G = {}
    with mp.workdps(dps):
        for a in alphabet:
            ae = sp.sympify(a)
            G[str(a)] = {}
            for v in kinvars:
                d = sp.together(sp.diff(sp.log(ae), sp.Symbol(v))).subs(subs)
                G[str(a)][v] = _sp_to_mpc(d, dps)
    return G


def _ut_check(samples: list[DESample], rotation: Rotation) -> dict:
    """Canonical-form sanity:  J^{(0)} should be POINT-INDEPENDENT and
    ∂J^{(0)}/∂x_i ≈ 0.  Returns per-master spread digits and a verdict."""
    M = len(rotation.masters)
    spread = {}
    for m, t in enumerate(rotation.masters):
        vals = [s.J[0][m] for s in samples]
        ref = vals[0]
        sp_ = max((abs(v - ref) for v in vals[1:]), default=mp.mpf(0))
        d = int(-mp.log10(sp_ / max(abs(ref), mp.mpf(1)))) if sp_ > 0 else 999
        # ∂J^{(0)} magnitude
        dmax = max(abs(s.dJ[v][0][m]) for s in samples for v in s.dJ)
        spread[t] = {"J0": mp.nstr(ref, 12), "spread_digits": d,
                     "max_dJ0": mp.nstr(dmax, 4)}
    ok = all(v["spread_digits"] >= 6 for v in spread.values())
    return {"per_master": spread, "ok": ok,
            "note": ("rotated leading order is constant ⇒ rotation is UT-ready"
                     if ok else
                     "rotated leading order VARIES across points — diagonal LS "
                     "rotation is NOT canonical for the flagged masters; the "
                     "fit will proceed but those rows may not rationalize.")}


def fit_dlog_connection(alphabet: list, kinvars: list[str],
                        samples: list[DESample], rotation: Rotation, *,
                        dps: int = 60, tol_rationalize: int | None = None,
                        max_coeff_bits: int = 32,
                        use_auto_condition: bool = True,
                        support: dict | None = None) -> dict:
    r"""Fit the constant ℚ-matrices {A_a} from sampled rotated Laurent data.

    The order-k slice of (★) reads, row-by-row,

        ∂J_r^{(k)}/∂x_i |_p  =  Σ_a Σ_c  (A_a)_{rc} · ∂log(a)/∂x_i|_p · J_c^{(k-1)}(p)

    so for each row r the unknowns are the n_letters·M numbers
    {(A_a)_{rc}}_{a,c}, and each (sample p, variable x_i, order k≥1) gives
    one equation.  Solve by least squares (mpmath, full precision), then
    PSLQ/LLL-recognise each entry as p/q.

    SPARSE-SUPPORT MODE  (`support` kwarg)
    --------------------------------------
    The canonical (UT) DE matrix is extremely sparse: master r couples to only
    a handful of daughters c (block-/sector-triangular), and each coupling
    carries only the letters of the (r,c) sub-topology.  Passing

        support = {r: [c0, c1, ...], ...}            # row → live block-cols
        support = {r: {c0: ["letter", ...], ...}}    # row → col → live letters

    restricts each row's fit to ONLY the live (letter, col) unknowns.  This
    collapses the binding point-count from ~2·nL·M/(nV·K) (every entry free) to
    ~2·max_row_unknowns/(nV·K) — a handful of points instead of tens.  Entries
    outside the support are forced to 0 (and verified ≈0 by the held-out check).
    The support is a SAFE SUPERSET when taken from the raw-Laporta coupling
    skeleton: a diagonal LS rotation cannot create new couplings, only kill
    them, so the UT support ⊆ raw support.

    Returns
    -------
    dict with keys
        A_dense   : {letter: [["p/q",..],..]}      (M×M, "?" where PSLQ failed)
        A         : {letter: [[i,j,"p/q"],...]}    (sparse, only nonzero)
        A_float   : {letter: [[float,..],..]}      (the fitted floats)
        honesty   : {n_points_used, n_eq_per_row, cond_number,
                     residual_digits, failed_to_rationalize, ut_check}
    """
    mp.mp.dps = dps
    masters = rotation.masters
    M = len(masters)
    letters = [str(a) for a in alphabet]
    nL = len(letters)
    n_unk = nL * M
    li = {a: i for i, a in enumerate(letters)}

    # Normalise the support into  row → list[(letter_idx, col)]  global-column
    # indices (into the dense nL·M unknown vector, col = letter_idx*M + c).
    row_cols: dict[int, list[int]] | None = None
    if support is not None:
        row_cols = {}
        for r in range(M):
            ent = support.get(r, support.get(str(r)))
            cols: list[int] = []
            if ent is None:
                cols = []                                   # decoupled row
            elif isinstance(ent, dict):                     # col → letters
                for c, las in ent.items():
                    c = int(c)
                    for a in las:
                        cols.append(li[str(a)] * M + c)
            else:                                           # list of cols, all letters
                for c in ent:
                    c = int(c)
                    for ai in range(nL):
                        cols.append(ai * M + c)
            row_cols[r] = sorted(set(cols))

    if tol_rationalize is None:
        tol_rationalize = max(8, min(s.deriv_digits for s in samples) // 2)

    # Precompute dlog gradients per sample point
    G = {id(s): _dlog_grad(alphabet, kinvars, s.point, dps) for s in samples}

    ut = _ut_check(samples, rotation)

    # Row-by-row design matrix  (shared for every r:  columns = (a,c),
    # row = (p,i,k) → entry = ∂log(a)/∂x_i|_p · J_c^{(k-1)}(p))
    design_rows: list[list] = []
    rhs_tags: list[tuple] = []         # (sample_idx, var, k)
    for si, s in enumerate(samples):
        Gp = G[id(s)]
        for v in kinvars:
            for k in s.orders:
                if k < 1:
                    continue
                row = [mp.mpc(0)] * n_unk
                for ai, a in enumerate(letters):
                    g = Gp[a][v]
                    if g == 0:
                        continue
                    for c in range(M):
                        row[ai * M + c] = g * s.J[k - 1][c]
                design_rows.append(row)
                rhs_tags.append((si, v, k))
    n_eq = len(design_rows)
    # In dense mode every row carries n_unk unknowns; in sparse mode the binding
    # figure is the widest row support.  Guard on the relevant one.
    max_row_unk = (n_unk if row_cols is None
                   else max((len(v) for v in row_cols.values()), default=0))
    if n_eq < max_row_unk:
        raise RuntimeError(
            f"under-determined: {n_eq} equations < {max_row_unk} unknowns "
            f"(densest row); need more points / ε-orders / kinvars")

    # Build mp design matrix once (real part — A_a are real ℚ, samples Euclidean)
    Amat = mp.matrix(n_eq, n_unk)
    for i, row in enumerate(design_rows):
        for j, v in enumerate(row):
            Amat[i, j] = mp.re(v)

    # Solve per row r
    A_float = {a: [[mp.mpf(0)] * M for _ in range(M)] for a in letters}
    A_dense = {a: [["?"] * M for _ in range(M)] for a in letters}
    A_sparse: dict[str, list] = {a: [] for a in letters}
    failed: list = []
    resid_digits_min = 999
    cond_reported = None

    # In sparse mode each row may project onto a DIFFERENT column subset, so the
    # conditioning / LLL reduction must be redone per row.  In dense mode it is
    # shared (computed once, lazily below).
    shared = None       # (cond0, U, Aeff, cond_after) for the dense full design

    def _prep(sub: "mp.matrix"):
        c0 = _cond_svd(sub)
        if use_auto_condition and c0 > 1e8:
            U_, _a, c1_, Aeff_ = auto_condition(sub, target_cond=1e6, _log=False)
            return c0, U_, Aeff_, c1_
        return c0, None, sub, c0

    for r in range(M):
        if row_cols is not None:
            cols = row_cols[r]                       # restricted unknowns
            if not cols:                              # decoupled row → all zero
                resid_digits_min = min(resid_digits_min, 999)
                continue
            sub = mp.matrix(n_eq, len(cols))
            for jj, gc in enumerate(cols):
                for i in range(n_eq):
                    sub[i, jj] = Amat[i, gc]
            cond0_r, U, Aeff, cond_after_r = _prep(sub)
            full_sub = sub
        else:
            if shared is None:
                shared = _prep(Amat)
            cond0_r, U, Aeff, cond_after_r = shared
            cols = list(range(n_unk))
            full_sub = Amat
        if cond_reported is None:
            cond_reported = (cond0_r, cond_after_r)

        b = mp.matrix(n_eq, 1)
        for i, (si, v, k) in enumerate(rhs_tags):
            b[i, 0] = mp.re(samples[si].dJ[v][k][r])
        # solve on (possibly LLL-reduced) basis
        if U is not None:
            cU, _ = _lstsq(Aeff, b)
            csol = U * cU      # back to original (restricted) column ordering
        else:
            csol, _ = _lstsq(Aeff, b)
        # residual (against the restricted design)
        res = full_sub * csol - b
        rn = mp.sqrt(sum(abs(res[i, 0])**2 for i in range(n_eq)))
        bn = mp.sqrt(sum(abs(b[i, 0])**2 for i in range(n_eq)))
        rd = int(-mp.log10(rn / max(bn, mp.mpf(1)))) if rn > 0 else 999
        resid_digits_min = min(resid_digits_min, rd)
        # scatter the restricted solution back to (letter, col) and rationalize
        for jj, gc in enumerate(cols):
            a = letters[gc // M]
            col = gc % M
            x = csol[jj, 0]
            A_float[a][r][col] = x
            q, _ = _pslq_rational(mp.mpf(x), tol_rationalize,
                                  max_coeff_bits=max_coeff_bits)
            if q is None:
                failed.append([a, r, col, mp.nstr(x, 16)])
                A_dense[a][r][col] = "?"
            else:
                A_dense[a][r][col] = q
                if q not in ("0", "0/1"):
                    A_sparse[a].append([r, col, q])
        # entries outside the support stay "0" in A_dense (sparse mode)
        if row_cols is not None:
            live = set(cols)
            for ai in range(nL):
                for col in range(M):
                    if ai * M + col not in live:
                        A_dense[letters[ai]][r][col] = "0"

    cond0 = cond_reported[0] if cond_reported else mp.mpf(1)
    cond_after = cond_reported[1] if cond_reported else mp.mpf(1)

    return {
        "alphabet": letters,
        "kinvars": list(kinvars),
        "masters": list(masters),
        "rotation": rotation.to_json(),
        "A": A_sparse,
        "A_dense": A_dense,
        "A_float": {a: [[mp.nstr(x, 20) for x in row] for row in M_]
                    for a, M_ in A_float.items()},
        "honesty": {
            "backend": "finite_diff",
            "n_points_used": len(samples),
            "n_eq_per_row": n_eq,
            "n_unknowns_per_row": n_unk,
            "sparse_support": row_cols is not None,
            "max_row_unknowns": (max_row_unk if row_cols is not None else n_unk),
            "support_pairs": (sum(len(v) for v in row_cols.values())
                              if row_cols is not None else None),
            "cond_number": mp.nstr(cond0, 4),
            "cond_after_LLL": mp.nstr(cond_after, 4),
            "tol_rationalize_digits": tol_rationalize,
            "fit_residual_digits": resid_digits_min,
            "failed_to_rationalize": failed,
            "ut_check": ut,
            "elliptic_masters_flagged": sorted(rotation.elliptic),
        },
    }


# ===========================================================================
#  verify_connection — held-out residual
# ===========================================================================
def verify_connection(conn: dict, samples_heldout: list[DESample], *,
                      dps: int = 60) -> int:
    """Return the minimum agreement digits of (★) on the held-out samples,
    using the RATIONALIZED A_a (so any '?' entry uses the float fallback)."""
    mp.mp.dps = dps
    letters = conn["alphabet"]
    kinvars = conn["kinvars"]
    masters = conn["masters"]
    M = len(masters)
    # build numeric A_a (rational where available, float fallback)
    A = {}
    for a in letters:
        Ma = mp.matrix(M, M)
        for r in range(M):
            for c in range(M):
                q = conn["A_dense"][a][r][c]
                if q != "?":
                    if "/" in q:
                        p_, q_ = q.split("/"); Ma[r, c] = mp.mpf(p_) / mp.mpf(q_)
                    else:
                        Ma[r, c] = mp.mpf(q)
                else:
                    Ma[r, c] = mp.mpf(conn["A_float"][a][r][c])
        A[a] = Ma
    min_d = 999
    for s in samples_heldout:
        Gp = _dlog_grad(letters, kinvars, s.point, dps)
        for v in kinvars:
            for k in s.orders:
                if k < 1:
                    continue
                lhs = mp.matrix([[mp.re(x)] for x in s.dJ[v][k]])
                Jkm1 = mp.matrix([[mp.re(x)] for x in s.J[k - 1]])
                rhs = mp.matrix(M, 1)
                for a in letters:
                    rhs += mp.re(Gp[a][v]) * (A[a] * Jkm1)
                diff = lhs - rhs
                dn = max(abs(diff[i, 0]) for i in range(M))
                ln = max(abs(lhs[i, 0]) for i in range(M)) or mp.mpf(1)
                d = int(-mp.log10(dn / ln)) if dn > 0 else 999
                min_d = min(min_d, d)
    return min_d


# ===========================================================================
#  connection_once — one command
# ===========================================================================
def _gen_points(kinvars: list[str], n: int, seed: int = 0,
                lo: int = -29, hi: int = -3) -> list[dict]:
    """Generic Euclidean rational points (small negative integers, distinct)."""
    rnd = random.Random(seed)
    pts, seen = [], set()
    while len(pts) < n:
        p = {v: rnd.randint(lo, hi) for v in kinvars}
        # avoid degenerate alphabet zeros (any letter = 0 at p)
        key = tuple(sorted(p.items()))
        if key in seen:
            continue
        seen.add(key)
        pts.append(p)
    return pts


def _alphabet_from_target(target: dict) -> list[str]:
    if "alphabet" in target:
        a = target["alphabet"]
        return list(a.values()) if isinstance(a, dict) else list(a)
    # try Landau Alphabet on the graph block
    try:
        import landau_alphabet as LB                       # type: ignore
        gs = LB.load_spec(target["graph"])
        log = LB.run(gs, out_path=None, face_timeout=30)
        return list(log["alphabet"])
    except Exception as e:
        raise RuntimeError(
            f"no 'alphabet' in target and Landau Alphabet failed: {e}")


def connection_once(target_json: str, *, n_fit: int | None = None,
                    n_heldout: int = 2, h_pow: int = 18,
                    goal_digits: int = 50, target_eps_power: int = 4,
                    work_root: str | None = None, parallel: int = 4,
                    out_path: str | None = None, seed: int = 0,
                    timeout_s: int = 1800,
                    backend: str = "finite_diff",
                    support: dict | None = None) -> dict:
    """alphabet → masters → sample K points → fit → verify → connection.json.

    `support` (or target['support']): sparse coupling skeleton row→[cols] or
    row→{col:[letters]} (see fit_dlog_connection).  When given, n_fit is sized
    from the DENSEST row, not nL·M — the connection-route point-count collapse.
    """
    target = json.load(open(target_json)) if isinstance(target_json, str) \
        else dict(target_json)
    name = target.get("name", "family")
    masters = target["masters"]
    kinvars = list(target.get("kinvars")
                   or target["graph"]["kinematics"]["invariants"])
    alphabet = _alphabet_from_target(target)
    M, nL, nV = len(masters), len(alphabet), len(kinvars)
    if support is None:
        support = target.get("support")

    # Rotation: explicit > canonical_form default
    if "rotation" in target:
        rot = Rotation.from_spec(target["rotation"], masters, kinvars)
    else:
        rot = Rotation.from_canonical_form(target["family"], masters, kinvars)

    # K ≳ n·|kinvars| / (|kinvars|·|orders|) + slack.
    # Dense: row-unknowns = nL·M, equations/row = n_fit·nV·K_ord.  Choose n_fit
    # so eqs ≥ 2·unknowns.  Sparse: unknowns/row = widest row support.
    K_ord = target_eps_power
    if n_fit is None:
        if support is not None:
            def _rowdeg(ent):
                if ent is None:
                    return 0
                if isinstance(ent, dict):
                    return sum(len(v) for v in ent.values())
                return len(ent) * nL
            max_unk = max((_rowdeg(support.get(r, support.get(str(r))))
                           for r in range(M)), default=0)
            n_fit = max(3, -(-2 * max_unk // (nV * K_ord)))
        else:
            n_fit = max(3, -(-2 * nL * M // (nV * K_ord)))
    pts = (target.get("points") or [])[:]
    need = n_fit + n_heldout
    if len(pts) < need:
        pts += _gen_points(kinvars, need - len(pts), seed=seed)
    fit_pts, ho_pts = pts[:n_fit], pts[n_fit:n_fit + n_heldout]

    work_root = work_root or os.path.join(
        HERE, "out", name, time.strftime("%Y%m%d-%H%M%S"))
    sys.stderr.write(
        f"[connection_once] {name}: M={M} masters, nL={nL} letters, "
        f"kinvars={kinvars}; fit on {n_fit} pts (+{n_heldout} held-out)\n"
        f"  work_root = {work_root}\n")

    h = Fraction(1, 10 ** h_pow)
    t0 = time.time()
    samples = sample_de(target, fit_pts + ho_pts, rotation=rot,
                        kinvars=kinvars, h=h, goal_digits=goal_digits,
                        target_eps_power=target_eps_power,
                        work_root=work_root, parallel=parallel,
                        timeout_s=timeout_s, backend=backend)
    t_sample = time.time() - t0
    fit_s, ho_s = samples[:n_fit], samples[n_fit:]

    conn = fit_dlog_connection(alphabet, kinvars, fit_s, rot,
                               dps=goal_digits + 10, support=support)
    conn["honesty"]["backend"] = backend
    ho_d = verify_connection(conn, ho_s, dps=goal_digits + 10)
    n_calls = (len(samples) * (1 + 2 * nV) if backend == "finite_diff"
               else len(samples) * 2)
    conn["honesty"]["min_heldout_digits"] = ho_d
    conn["honesty"]["n_heldout"] = len(ho_s)
    conn["honesty"]["heldout_points"] = [s.point for s in ho_s]
    conn["honesty"]["timing"] = {
        "sample_s": round(t_sample, 1),
        "n_amflow_calls": n_calls,
        "per_call_s": round(t_sample / max(1, n_calls), 2),
        "vs_symbolic_firefly": (
            f"symbolic IBP over {nV} vars would reconstruct {M}×{M} entries × "
            f"{nV} matrices as multivariate rationals (FireFly); for many-scale "
            f"families that is ~30 min/node.  This route used "
            f"{n_calls} cheap numeric IBPs."),
    }
    conn["name"] = name

    out_path = out_path or os.path.join(HERE, "out", name, "connection.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(conn, f, indent=2)
    sys.stderr.write(f"[connection_once] wrote → {out_path}\n"
                     f"  fit residual: {conn['honesty']['fit_residual_digits']} d, "
                     f"held-out: {ho_d} d, "
                     f"failed-to-rationalize: "
                     f"{len(conn['honesty']['failed_to_rationalize'])}\n")
    return conn


# ===========================================================================
#  CLI
# ===========================================================================
def _main(argv=None):
    ap = argparse.ArgumentParser(
        description="dlog-connection-once: fit constant A_a from numeric DE samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", help="target JSON (family + masters + alphabet)")
    ap.add_argument("--n-fit", type=int, default=None)
    ap.add_argument("--n-heldout", type=int, default=2)
    ap.add_argument("--h-pow", type=int, default=18,
                    help="finite-diff step h = 10^{-h_pow}")
    ap.add_argument("--goal-digits", type=int, default=50)
    ap.add_argument("--eps-power", type=int, default=4)
    ap.add_argument("--parallel", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--work-root", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--backend", default="finite_diff",
                    choices=("finite_diff", "amflow_diffeq"))
    ap.add_argument("--support", default=None,
                    help="JSON file: sparse coupling skeleton row→[cols] or "
                         "row→{col:[letters]} (collapses the fit point-count)")
    a = ap.parse_args(argv)
    sup = None
    if a.support:
        raw = json.load(open(a.support))
        sup = raw.get("support", raw)
    connection_once(a.target, n_fit=a.n_fit, n_heldout=a.n_heldout,
                    h_pow=a.h_pow, goal_digits=a.goal_digits,
                    target_eps_power=a.eps_power, work_root=a.work_root,
                    parallel=a.parallel, out_path=a.out, seed=a.seed,
                    timeout_s=a.timeout, backend=a.backend, support=sup)


if __name__ == "__main__":
    _main()

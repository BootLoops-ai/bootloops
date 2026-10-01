"""qinvert.solver — profile-enumeration engine for quantile-statistic inversion.

Given a dataset (sorted Fractions) and a list of published statistic targets,
find minimal-size delta families (k_add additions with value windows +
k_rem removals of released values) that make every target recompute EXACTLY.

Method:

1. Index-shift lemma: a delta of k = k_add + k_rem entries moves the rank of
   any fixed value by at most k, so the order statistic at position j of the
   new array is either an added value or an original value with index within
   [j - k_add, j + k_rem].  All candidate enumeration below is windowed by
   this lemma (window pad k + 4).
2. For each (k_add, k_rem) and each removal base, the targets pin "blocks":
   point constraints x'_(j) = v and weighted-average constraints
   w_j x'_(j) + w_{j+1} x'_(j+1) = v derived from each statistic's
   order-stat support at n' = n - k_rem + k_add (targets.support).  The
   averaging/interpolation case follows an n mod 4 parity
   argument: for QNTLDEF=5 at p = 1/4, 3/4 the two-index (averaging) branch
   fires iff 4 | n', and the general rule per definition comes straight from
   support(n').
   Blocks pinned to the SAME value are merged into one span constraint
   (every covered order stat equals v) — scenarios with two blocks pinned
   to the same value are enumerated, never silently pruned.
3. Per block, enumerate realizations: counts of adds pinned AT critical
   values plus an interval of allowed counts strictly below the block
   (linear inequalities in exact counts against the base's cumulative
   counts).  Weighted constraints admit two branches: both order stats equal
   v, or a straddle pair (s, s') with w_j s + w_{j+1} s' = v and no merged
   value strictly between.
4. Assemble realizations left-to-right with filler adds in the open regions
   between blocks; filler values carry WINDOWS (any register point strictly
   inside the region works — the constraints depend only on region counts).
5. Verification gate (authoritative): every witness is materialized and every
   target re-evaluated exactly via targets.evaluate_targets; failures are
   dropped.  The gate makes the enumeration sound; completeness is bounded by
   the loud caps below (see README "Honest limits").

Exactness: fractions.Fraction everywhere; no float participates in any
accept/reject decision (there are no float-semantic sites in this tool).
Multiprocessing: spawn context only — fork+BLAS deadlocks.
"""

import math
import multiprocessing as mp
from bisect import bisect_left, bisect_right
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement

from .targets import (Interval, LinearCombo, Mean, OrderStat, Quantile,
                      Statistic, TukeyFence, as_fraction, evaluate_targets,
                      invert_fence_pair, targets_exact)

__all__ = ["solve", "solve_many", "DEFAULT_CAPS"]

# ---------------------------------------------------------------------------
# loud caps: enumeration is complete only up to these; every hit is recorded
# in result["caps_hit"] so a truncated search can never masquerade as a proof
# of emptiness.
# ---------------------------------------------------------------------------
DEFAULT_CAPS = {
    "grid_cap": 800,          # register multiples per candidate window
    "scenario_cap": 20000,    # candidate pinned-value scenarios per (k_add,k_rem)
    "assembly_budget": 300000,  # filler/realization iterations per add-solve
    "rem_cands_cap": 48,      # distinct removal candidate values
    "rem_bases_cap": 20000,   # removal multisets per (k, k_rem)
    "probe_budget": 20000,    # verify calls in the probe fallback
    "probe_k_max": 3,         # probe fallback only for k_add <= this
}

WINDOW_PAD = 4  # index-shift lemma window padding beyond k


# ---------------------------------------------------------------------------
# base counts (sorted distinct values + cumulative counts)
# ---------------------------------------------------------------------------

class BaseCounts:
    def __init__(self, sorted_vals):
        self.vals = list(sorted_vals)
        self.n = len(self.vals)
        dv, cum = [], []
        for v in self.vals:
            if not dv or v != dv[-1]:
                dv.append(v)
                cum.append(1 if not cum else cum[-1] + 1)
            else:
                cum[-1] += 1
        self.dv, self.cum = dv, cum

    def cnt_le(self, v):
        i = bisect_right(self.dv, v)
        return self.cum[i - 1] if i else 0

    def cnt_lt(self, v):
        i = bisect_left(self.dv, v)
        return self.cum[i - 1] if i else 0

    def cnt_at(self, v):
        return self.cnt_le(v) - self.cnt_lt(v)


# ---------------------------------------------------------------------------
# small exact helpers
# ---------------------------------------------------------------------------

def _floor(x):
    return x.numerator // x.denominator


def _ceil(x):
    return -((-x).numerator // (-x).denominator)


def _addable_fn(reg, vlo, vhi):
    def addable(v):
        if vlo is not None and v < vlo:
            return False
        if vhi is not None and v > vhi:
            return False
        if reg is not None and (v / reg).denominator != 1:
            return False
        return True
    return addable


def _component_quantiles(stat):
    """Flatten a statistic into its order-stat components (Quantile/OrderStat
    objects) for windowing purposes.  Dense components yield nothing."""
    if isinstance(stat, Quantile) or isinstance(stat, OrderStat):
        return [stat]
    if isinstance(stat, TukeyFence):
        return [stat.q_lo, stat.q_hi]
    if isinstance(stat, LinearCombo):
        out = []
        for _, s in stat.terms:
            out.extend(_component_quantiles(s))
        return out
    return []


# ---------------------------------------------------------------------------
# filler windows and canonical values
# ---------------------------------------------------------------------------

def _canonical_in(lo_b, hi_b, lo_strict, hi_strict, reg):
    """A canonical value in the window (None if empty).  Bounds may be None
    (unbounded).  With a register, the canonical value is a register multiple;
    without, a rational interior point."""
    if reg is not None:
        if hi_b is not None:
            m = _floor(hi_b / reg)
            v = m * reg
            if hi_strict and v == hi_b:
                v -= reg
            if lo_b is not None and (v < lo_b or (lo_strict and v == lo_b)):
                return None
            return v
        if lo_b is not None:
            m = _ceil(lo_b / reg)
            v = m * reg
            if lo_strict and v == lo_b:
                v += reg
            return v
        return F(0)
    # no register: any rational
    if lo_b is not None and hi_b is not None:
        if lo_b > hi_b or (lo_b == hi_b and (lo_strict or hi_strict)):
            return None
        return (lo_b + hi_b) / 2
    if hi_b is not None:
        return hi_b - 1 if hi_strict else hi_b
    if lo_b is not None:
        return lo_b + 1 if lo_strict else lo_b
    return F(0)


def _window_int_range(win, reg):
    """Integer register-unit range [lo_u, hi_u] (None = unbounded) of a
    filler window (lo_b, hi_b, lo_strict, hi_strict)."""
    lo_b, hi_b, lo_s, hi_s = win
    lo_u = hi_u = None
    if lo_b is not None:
        lo_u = _ceil(lo_b / reg)
        if lo_s and lo_u * reg == lo_b:
            lo_u += 1
    if hi_b is not None:
        hi_u = _floor(hi_b / reg)
        if hi_s and hi_u * reg == hi_b:
            hi_u -= 1
    return lo_u, hi_u


# ---------------------------------------------------------------------------
# block realizations
# (generalized from membership_stacker._quartile_realizations)
# ---------------------------------------------------------------------------

def _real_point(bc, j, v, k, addable):
    """Realizations of x'_(j) = v with up to k adds (special case of
    _real_span with j_lo == j_hi)."""
    return _real_span(bc, j, j, v, k, addable)


def _real_span(bc, j_lo, j_hi, v, k, addable):
    """Realizations of x'_(i) = v for EVERY i in [j_lo, j_hi] with up to k
    adds.  Each realization:
    {"crits": [(value, count>0), ...], "lo_rng": (min, max) allowed adds
     strictly below `below`, "below": block bottom, "top": block top}.

    Exact counting (with e added copies AT v and c adds strictly below v):
      x'_(j_lo) >= v  <=>  c <= j_lo - 1 - cnt_lt(v)
      x'_(j_hi) <= v  <=>  c >= j_hi - cnt_le(v) - e
    A span of several indices may need up to k fresh copies of v when v is
    scarce in the base, so e ranges over 0..k — scarce-in-base entries are
    legal (planted-truth receipts tests/test_all.py::TestSameValueBlocks)."""
    outs = []
    for e in range(0, k + 1):
        if e and not addable(v):
            break
        lo_min = j_hi - bc.cnt_le(v) - e
        lo_max = j_lo - 1 - bc.cnt_lt(v)
        if lo_max < 0 or lo_min > k:
            continue
        outs.append({"crits": [(v, e)] if e else [],
                     "lo_rng": (max(0, lo_min), lo_max),
                     "below": v, "top": v})
    return outs


def _real_wavg(bc, j, wj, v, k, addable, reg):
    """Realizations of wj*x'_(j) + (1-wj)*x'_(j+1) = v.

    Branch (a): both order stats equal v (the g-independent branch).
    Branch (b): straddle pair s < v < s' with wj*s + (1-wj)*s' = v and no
    merged value strictly between s and s'; then exactly j merged values
    are <= s.  Candidate s: register steps below v, existing values below v
    in the index-shift window, and values derived from existing t > v via
    s = (v - (1-wj)*t)/wj."""
    wj1 = 1 - wj
    outs = []
    # (a) both equal v; e ranges 0..k for the same reason as _real_span
    # (v may be scarce in the base and the region below blocked)
    for e in range(0, k + 1):
        if e and not addable(v):
            break
        lo_min = j + 1 - bc.cnt_le(v) - e
        lo_max = j - 1 - bc.cnt_lt(v)
        if lo_max < 0 or lo_min > k:
            continue
        outs.append({"crits": [(v, e)] if e else [],
                     "lo_rng": (max(0, lo_min), lo_max),
                     "below": v, "top": v})
    # (b) straddles
    s_cands = set()
    if reg is not None:
        base = _floor(v / reg) * reg
        s = base if base < v else base - reg
        for _ in range(4):
            s_cands.add(s)
            s -= reg
    i_hi = bisect_left(bc.dv, v)
    for i in range(max(0, i_hi - (k + WINDOW_PAD)), i_hi):
        s_cands.add(bc.dv[i])
    i_lo2 = bisect_right(bc.dv, v)
    for i in range(i_lo2, min(len(bc.dv), i_lo2 + k + WINDOW_PAD)):
        t = bc.dv[i]
        s = (v - wj1 * t) / wj
        if s < v:
            s_cands.add(s)
    for s in sorted(s_cands, reverse=True):
        sp = (v - wj * s) / wj1
        if not (s < v < sp):
            continue
        if bc.cnt_lt(sp) != bc.cnt_le(s):
            continue  # merged data strictly between s and s'
        s_ex = bc.cnt_at(s) > 0
        sp_ex = bc.cnt_at(sp) > 0
        for es in range(0, 3):
            if es and not addable(s):
                break
            if not (s_ex or es):
                continue
            for esp in range(0, 3):
                if esp and not addable(sp):
                    break
                if not (sp_ex or esp):
                    continue
                c_lo = j - bc.cnt_le(s) - es
                if c_lo < 0 or c_lo > k:
                    continue
                if bc.cnt_lt(s) + c_lo > j - 1:
                    continue
                if bc.cnt_le(sp) + c_lo + es + esp < j + 1:
                    continue
                crits = [(s, es)] if es else []
                if esp:
                    crits.append((sp, esp))
                outs.append({"crits": crits, "lo_rng": (c_lo, c_lo),
                             "below": s, "top": sp})
    return outs


# ---------------------------------------------------------------------------
# candidate values for a free (gridded) quantile
# (generalized from membership_stacker._q_pair_candidates.window_vals)
# ---------------------------------------------------------------------------

def _window_values(bc, quantile, n_new, k, reg, caps, caps_hit):
    """Candidate exact VALUES the quantile could take on the new array,
    windowed by the index-shift lemma."""
    sup = quantile.support(n_new)
    j = sup[0][0]
    i_lo = max(0, j - k - WINDOW_PAD - 1)
    i_hi = min(bc.n - 1, j + k + WINDOW_PAD - 1)
    if bc.n == 0:
        return set()
    vlo_w, vhi_w = bc.vals[i_lo], bc.vals[i_hi]
    out = set()
    d0 = bisect_left(bc.dv, vlo_w)
    d1 = bisect_right(bc.dv, vhi_w)
    exist = bc.dv[d0:d1]
    out.update(exist)
    if len(sup) == 2:
        wj = sup[0][1]
        wj1 = sup[1][1]
        for a, b in zip(exist, exist[1:]):
            out.add(wj * a + wj1 * b)
    if reg is not None:
        # register (sub)lattice: interpolated quantiles of register-valued
        # adds live on reg / lcm(weight denominators)
        den = 1
        if len(sup) == 2:
            den = sup[0][1].denominator
            d2 = sup[1][1].denominator
            den = den * d2 // math.gcd(den, d2)
        if den > 8:
            caps_hit.add(f"grid-sublattice-den-{den}-capped-at-reg")
            den = 1
        step = reg / den
        start = _ceil(vlo_w / step)
        stop = _floor(vhi_w / step)
        if stop - start + 1 <= caps["grid_cap"]:
            out.update(m * step for m in range(start, stop + 1))
        else:
            caps_hit.add("grid_cap")
    else:
        # register=None: the free-quantile grid holds data values and their
        # interpolations ONLY — a structural sub-lattice limit (README
        # "Honest limits").  Recorded loudly so emptiness in register-None
        # one-sided fence modes is distinguishable from exhaustiveness by
        # caps_hit alone (receipt: tests/test_all.py::
        # TestLoudness.test_register_none_fence_grid_is_loud).
        caps_hit.add("register-none-grid-data-values-only")
    return out


def _tail_values(bc, side, k, reg, vlo, vhi, published):
    """Candidate values for min (side='lo') or max (side='hi') on the new
    array: nearby existing values, register steps beyond the extreme, the
    bound itself, and the published value."""
    out = set()
    if bc.n:
        if side == "lo":
            for i in range(0, min(len(bc.dv), k + WINDOW_PAD)):
                out.add(bc.dv[i])
            if reg is not None:
                for t in range(1, k + WINDOW_PAD):
                    out.add(bc.dv[0] - t * reg)
        else:
            for i in range(max(0, len(bc.dv) - k - WINDOW_PAD), len(bc.dv)):
                out.add(bc.dv[i])
            if reg is not None:
                for t in range(1, k + WINDOW_PAD):
                    out.add(bc.dv[-1] + t * reg)
    if vlo is not None:
        out.add(vlo)
    if vhi is not None:
        out.add(vhi)
    if published is not None:
        out.add(published)
    return out


# ---------------------------------------------------------------------------
# scenario generation: pinned (quantile -> required value) maps
# ---------------------------------------------------------------------------

def _scenarios(tg, bc, k_add, n_new, reg, vlo, vhi, caps, caps_hit):
    """Return (scenarios, mean_req).  Each scenario is a list of block dicts
    {"kind": "point"|"wavg", "j": idx, "wj": weight, "v": value} sorted by j,
    or None overall if the targets are inconsistent at this n'."""
    qmap = {}      # (p, definition) -> required value
    omap = {}      # 1-based index at n_new -> required value
    mean_req = None
    fence_groups = {}
    gridded_choices = []   # list of lists of {qkey/okey: value} alternatives

    def put_q(key, v):
        if key in qmap and qmap[key] != v:
            return False
        qmap[key] = v
        return True

    for stat, pub in tg:
        if isinstance(pub, Interval):
            continue  # verification-only
        if isinstance(stat, Quantile):
            if not put_q((stat.p, stat.definition), pub):
                return [], None
        elif isinstance(stat, OrderStat):
            i = stat.index(n_new)
            if i in omap and omap[i] != pub:
                return [], None
            omap[i] = pub
        elif isinstance(stat, Mean):
            mean_req = pub
        elif isinstance(stat, TukeyFence):
            g = fence_groups.setdefault(stat.pair_key(), {
                "lo": None, "hi": None, "cap_lo": None, "cap_hi": None,
                "stat": stat})
            g[stat.side] = pub
            g["cap_lo" if stat.side == "lo" else "cap_hi"] = stat.cap
        elif isinstance(stat, LinearCombo):
            inv = _invertible_combo(stat)
            if inv is None:
                caps_hit.add("linearcombo-verify-only")
                continue
            (c1, os1), (c2, os2) = inv
            i1, i2 = os1.index(n_new), os2.index(n_new)
            alts = []
            side1 = "lo" if i1 <= n_new // 2 else "hi"
            side2 = "lo" if i2 <= n_new // 2 else "hi"
            for s in sorted(_tail_values(bc, side1, k_add, reg, vlo, vhi,
                                         None)):
                t = (pub - c1 * s) / c2
                alts.append({("idx", i1): s, ("idx", i2): t})
            for t in sorted(_tail_values(bc, side2, k_add, reg, vlo, vhi,
                                         None)):
                s = (pub - c2 * t) / c1
                alts.append({("idx", i1): s, ("idx", i2): t})
            gridded_choices.append(alts)
        else:
            caps_hit.add(f"undrivable-target-{type(stat).__name__}")

    for key, g in fence_groups.items():
        mult, definition, p_lo, p_hi = key
        q1r, q3r, mode = invert_fence_pair(g["lo"], g["hi"], mult,
                                           g["cap_lo"], g["cap_hi"])
        if mode == "both":
            if not put_q((p_lo, definition), q1r):
                return [], None
            if not put_q((p_hi, definition), q3r):
                return [], None
            continue
        m = mult
        q_lo_stat = Quantile(p_lo, definition)
        q_hi_stat = Quantile(p_hi, definition)
        alts = []
        seen_pairs = set()

        def add_pair(q1, q3):
            if q1 <= q3 and (q1, q3) not in seen_pairs:
                seen_pairs.add((q1, q3))
                alts.append({("q", p_lo, definition): q1,
                             ("q", p_hi, definition): q3})

        if mode == "hi_only":
            U = as_fraction(g["hi"])
            for q3 in sorted(_window_values(bc, q_hi_stat, n_new, k_add, reg,
                                            caps, caps_hit)):
                add_pair(((1 + m) * q3 - U) / m, q3)
            for q1 in sorted(_window_values(bc, q_lo_stat, n_new, k_add, reg,
                                            caps, caps_hit)):
                add_pair(q1, (U + m * q1) / (1 + m))
        elif mode == "lo_only":
            L = as_fraction(g["lo"])
            for q1 in sorted(_window_values(bc, q_lo_stat, n_new, k_add, reg,
                                            caps, caps_hit)):
                add_pair(q1, ((1 + m) * q1 - L) / m)
            for q3 in sorted(_window_values(bc, q_hi_stat, n_new, k_add, reg,
                                            caps, caps_hit)):
                add_pair((L + m * q3) / (1 + m), q3)
        else:  # 'none': both sides capped -> inequality only; verification
            continue
        gridded_choices.append(alts)

    # cartesian product of gridded choices with the fixed qmap/omap
    scenarios = [dict()]
    for alts in gridded_choices:
        nxt = []
        for sc in scenarios:
            for alt in alts:
                d = dict(sc)
                ok = True
                for kk, vv in alt.items():
                    if kk in d and d[kk] != vv:
                        ok = False
                        break
                    d[kk] = vv
                if ok:
                    nxt.append(d)
                if len(nxt) > caps["scenario_cap"]:
                    caps_hit.add("scenario_cap")
                    break
            if len(nxt) > caps["scenario_cap"]:
                break
        scenarios = nxt
        if not scenarios:
            return [], mean_req

    # materialize block lists
    out = []
    for sc in scenarios:
        blocks = []
        ok = True
        local_q = dict(qmap)
        local_o = dict(omap)
        for kk, vv in sc.items():
            if kk[0] == "q":
                key = (kk[1], kk[2])
                if key in local_q and local_q[key] != vv:
                    ok = False
                    break
                local_q[key] = vv
            else:
                i = kk[1]
                if i in local_o and local_o[i] != vv:
                    ok = False
                    break
                local_o[i] = vv
        if not ok:
            continue
        for (p, definition), v in local_q.items():
            sup = Quantile(p, definition).support(n_new)
            if len(sup) == 1:
                blocks.append({"kind": "point", "j": sup[0][0], "wj": F(1),
                               "v": v})
            else:
                blocks.append({"kind": "wavg", "j": sup[0][0],
                               "wj": sup[0][1], "v": v})
        for i, v in local_o.items():
            blocks.append({"kind": "point", "j": i, "wj": F(1), "v": v})
        if not blocks:
            continue  # nothing drivable in this scenario
        blocks.sort(key=lambda b: (b["j"], b["v"]))
        # disjointness of supports (overlaps are skipped loudly)
        used = set()
        for b in blocks:
            idxs = {b["j"]} if b["kind"] == "point" else {b["j"], b["j"] + 1}
            if used & idxs:
                same = [bb for bb in blocks
                        if bb["j"] == b["j"] and bb["kind"] == b["kind"]]
                if len(same) > 1 and all(bb["v"] == b["v"] for bb in same):
                    continue
                caps_hit.add("overlapping-supports-skipped")
                ok = False
                break
            used |= idxs
        if not ok:
            continue
        # merge exact duplicates
        dedup = []
        for b in blocks:
            if dedup and dedup[-1]["j"] == b["j"] and \
                    dedup[-1]["kind"] == b["kind"] and dedup[-1]["v"] == b["v"] \
                    and dedup[-1]["wj"] == b["wj"]:
                continue
            dedup.append(b)
        # merge runs of blocks pinned to the SAME value into one "span"
        # constraint x'_(i) = v over every covered index: order-statistic
        # monotonicity forces every index between two equal-pinned blocks to
        # the same value, and a wavg block in such a run can only realize its
        # both-equal branch (a straddle end sits strictly off v, contradicting
        # the equal-pinned neighbor).  Equal-quantile escape: when
        # neighboring blocks pin the SAME value the span merges
        # (order-statistic monotonicity) — without the merge the strict
        # left-to-right chaining prune in _assemble.choose would drop every
        # such scenario SILENTLY (planted-truth receipts
        # tests/test_all.py::TestSameValueBlocks).
        merged = []
        for b in dedup:
            if merged and merged[-1]["v"] == b["v"]:
                p = merged[-1]
                p_hi = p["j_hi"] if p["kind"] == "span" else (
                    p["j"] + 1 if p["kind"] == "wavg" else p["j"])
                b_hi = b["j"] + 1 if b["kind"] == "wavg" else b["j"]
                merged[-1] = {"kind": "span", "j": p["j"],
                              "j_hi": max(p_hi, b_hi), "wj": F(1),
                              "v": b["v"]}
            else:
                merged.append(b)
        out.append(merged)
    return out, mean_req


def _invertible_combo(stat):
    """A LinearCombo is drivable iff it has exactly two OrderStat terms with
    nonzero coefficients (e.g. Range = max - min)."""
    if not isinstance(stat, LinearCombo) or stat.dense:
        return None
    terms = [(c, s) for c, s in stat.terms if c != 0]
    if len(terms) == 2 and all(isinstance(s, OrderStat) for _, s in terms):
        return terms
    return None


# ---------------------------------------------------------------------------
# assembly: realizations x fillers -> witness add-multisets
# ---------------------------------------------------------------------------

def _assemble(bc, k_add, blocks, reg, vlo, vhi, mean_req, tg, seen, found,
              max_w, caps, caps_hit, budget):
    addable = _addable_fn(reg, vlo, vhi)
    reals_lists = []
    for b in blocks:
        if b["kind"] == "point":
            rl = _real_point(bc, b["j"], b["v"], k_add, addable)
        elif b["kind"] == "span":
            rl = _real_span(bc, b["j"], b["j_hi"], b["v"], k_add, addable)
        else:
            rl = _real_wavg(bc, b["j"], b["wj"], b["v"], k_add, addable, reg)
        if not rl:
            return budget
        reals_lists.append(rl)

    m = len(reals_lists)

    def emit(chosen, fcounts, fill_vals, gaps):
        specs = []
        for i, r in enumerate(chosen):
            if fcounts[i]:
                specs += [{"value": fill_vals[i], "window": gaps[i],
                           "role": "filler"}] * fcounts[i]
            for v, e in r["crits"]:
                specs += [{"value": v, "window": None, "role": "crit"}] * e
        if fcounts[m]:
            specs += [{"value": fill_vals[m], "window": gaps[m],
                       "role": "filler"}] * fcounts[m]
        if mean_req is not None:
            base_sum = sum(bc.vals)
            need = mean_req * (bc.n + k_add) - base_sum
            specs = _mean_adjust(specs, need, reg)
            if specs is None:
                caps_hit.add("mean-unadjustable-witness-dropped")
                return
        values = sorted(s["value"] for s in specs)
        key = tuple(values)
        if key in seen:
            return
        seen.add(key)
        after = sorted(bc.vals + values)
        if targets_exact(after, tg):
            found.append({"adds": specs, "n_new": len(after)})

    def choose(bi, chosen, b):
        if len(found) >= max_w or b[0] <= 0:
            return
        if bi == m:
            fill(chosen, b)
            return
        for r in reals_lists[bi]:
            b[0] -= 1
            if b[0] <= 0:
                caps_hit.add("assembly_budget")
                return
            if chosen and not (chosen[-1]["top"] < r["below"]):
                # Equality after the same-value merge means a straddle end
                # touching a NEIGHBORING pinned value (different targets,
                # different required values, shared boundary value).  The
                # chained count accounting does not model shared-boundary
                # adds, so these realizations are skipped — LOUDLY, as a
                # structural limit (receipt: tests/
                # test_all.py::TestLoudness.test_touching_straddle_is_loud).
                # Strict inequality (values out of order) is a genuine
                # infeasibility and stays silent.
                if chosen[-1]["top"] == r["below"]:
                    caps_hit.add("touching-realizations-skipped")
                continue
            choose(bi + 1, chosen + [r], b)

    def fill(chosen, b):
        crit_counts = [sum(e for _, e in r["crits"]) for r in chosen]
        if sum(crit_counts) > k_add:
            return
        gaps = []
        gaps.append((vlo, chosen[0]["below"], False, True))
        for i in range(m - 1):
            gaps.append((chosen[i]["top"], chosen[i + 1]["below"], True, True))
        gaps.append((chosen[-1]["top"], vhi, True, False))
        fill_vals = [_canonical_in(*g, reg) for g in gaps]
        crit_rest = [0] * (m + 1)
        for i in range(m - 1, -1, -1):
            crit_rest[i] = crit_rest[i + 1] + crit_counts[i]

        def frec(i, below_count, used, fcounts):
            if len(found) >= max_w or b[0] <= 0:
                return
            b[0] -= 1
            if b[0] <= 0:
                caps_hit.add("assembly_budget")
                return
            if i == m:
                f_top = k_add - used
                if f_top < 0:
                    return
                if f_top > 0 and fill_vals[m] is None:
                    return
                emit(chosen, fcounts + [f_top], fill_vals, gaps)
                return
            a, bb = chosen[i]["lo_rng"]
            fmin = max(0, a - below_count)
            fmax = min(bb - below_count, k_add - used - crit_rest[i])
            for f in range(fmin, fmax + 1):
                if f > 0 and fill_vals[i] is None:
                    break
                frec(i + 1, below_count + f + crit_counts[i],
                     used + f + crit_counts[i], fcounts + [f])

        frec(0, 0, 0, [])

    b = [budget]
    choose(0, [], b)
    return b[0]


def _mean_adjust(specs, need_sum, reg):
    """Shift filler values inside their windows so sum(adds) == need_sum.
    Crits are pinned.  With a register the fillers move on integer register
    units (exact greedy on interval slack); without, a single filler absorbs
    the whole deficit if its window permits.  Returns new specs or None."""
    cur = sum(s["value"] for s in specs)
    deficit = need_sum - cur
    if deficit == 0:
        return specs
    movable = [i for i, s in enumerate(specs) if s["role"] == "filler"]
    if not movable:
        return None
    if reg is not None:
        du = deficit / reg
        if du.denominator != 1:
            return None
        du = du.numerator
        out = [dict(s) for s in specs]
        for i in movable:
            if du == 0:
                break
            lo_u, hi_u = _window_int_range(out[i]["window"], reg)
            cu = out[i]["value"] / reg
            assert cu.denominator == 1
            cu = cu.numerator
            if du > 0:
                room = du if hi_u is None else min(du, hi_u - cu)
            else:
                room = du if lo_u is None else max(du, lo_u - cu)
            if (du > 0 and room > 0) or (du < 0 and room < 0):
                out[i]["value"] = (cu + room) * reg
                du -= room
        return out if du == 0 else None
    # no register: put the whole deficit on one filler if its window allows
    for i in movable:
        lo_b, hi_b, lo_s, hi_s = specs[i]["window"]
        v = specs[i]["value"] + deficit
        if lo_b is not None and (v < lo_b or (lo_s and v == lo_b)):
            continue
        if hi_b is not None and (v > hi_b or (hi_s and v == hi_b)):
            continue
        out = [dict(s) for s in specs]
        out[i]["value"] = v
        return out
    return None


def _mean_only_adds(bc, k_add, mean_req, reg, vlo, vhi):
    """Direct witnesses when Mean is the only driving target: k_add values,
    register-valid within [vlo, vhi], summing to mean_req*n' - base_sum."""
    need = mean_req * (bc.n + k_add) - sum(bc.vals)
    if k_add == 0:
        return []
    if reg is not None:
        t = need / reg
        if t.denominator != 1:
            return []
        t = t.numerator
        lo_u = _ceil(vlo / reg) if vlo is not None else None
        hi_u = _floor(vhi / reg) if vhi is not None else None
        if lo_u is not None and t < k_add * lo_u:
            return []
        if hi_u is not None and t > k_add * hi_u:
            return []
        q, r = divmod(t, k_add)
        units = [q + 1] * r + [q] * (k_add - r)
        # equal split lies within [lo_u, hi_u] whenever the sum does
        vals = [u * reg for u in units]
    else:
        v = need / k_add
        if vlo is not None and v < vlo:
            return []
        if vhi is not None and v > vhi:
            return []
        vals = [v] * k_add
    return [[{"value": v, "window": None, "role": "mean"} for v in vals]]


# ---------------------------------------------------------------------------
# probe fallback (interval-only target sets, small k)
# ---------------------------------------------------------------------------

def _probe_adds(bc, k_add, tg, reg, vlo, vhi, caps, caps_hit, seen, found,
                max_w):
    if k_add > caps["probe_k_max"]:
        caps_hit.add("probe-k-cap")
        return
    caps_hit.add("probe-fallback-used")
    n_new = bc.n + k_add
    addable = _addable_fn(reg, vlo, vhi)
    probes = set()
    for stat, _ in tg:
        for comp in _component_quantiles(stat):
            sup = comp.support(n_new)
            j = sup[0][0]
            for i in range(max(0, j - k_add - WINDOW_PAD - 1),
                           min(bc.n, j + k_add + WINDOW_PAD)):
                probes.add(bc.vals[i])
    if bc.n:
        probes |= {bc.vals[0], bc.vals[-1]}
        if reg is not None:
            probes |= {bc.vals[0] - reg, bc.vals[-1] + reg}
    if vlo is not None:
        probes.add(vlo)
    if vhi is not None:
        probes.add(vhi)
    probes = sorted(v for v in probes if addable(v))
    budget = caps["probe_budget"]
    for combo in combinations_with_replacement(probes, k_add):
        budget -= 1
        if budget <= 0:
            caps_hit.add("probe_budget")
            return
        values = sorted(combo)
        key = tuple(values)
        if key in seen:
            continue
        seen.add(key)
        after = sorted(bc.vals + values)
        if targets_exact(after, tg):
            found.append({"adds": [{"value": v, "window": None,
                                    "role": "probe"} for v in values],
                          "n_new": len(after)})
            if len(found) >= max_w:
                return


# ---------------------------------------------------------------------------
# removal bases
# ---------------------------------------------------------------------------

def _removal_candidates(xs, tg, k, reg, caps, caps_hit):
    """Distinct candidate removal values: index-shift windows around every
    target's support positions plus one representative per outside region.
    Complete up to region-equivalence for pure order-stat targets (removing
    any value in a region shifts all pinned indices identically); for dense
    targets (Mean) every distinct value matters, so all are included up to
    the loud cap."""
    n = len(xs)
    bc = BaseCounts(xs)
    dense = any(isinstance(s, Mean) or (isinstance(s, LinearCombo) and s.dense)
                for s, _ in tg)
    if dense:
        if len(bc.dv) <= caps["rem_cands_cap"]:
            return list(bc.dv)
        caps_hit.add("rem-cands-dense-truncated")
    idxs = set()
    for stat, _ in tg:
        for comp in _component_quantiles(stat):
            try:
                sup = comp.support(n)
            except ValueError:
                # OrderStat index invalid at the CURRENT n (e.g. j > n): it
                # contributes no removal window; removal only shrinks n, so
                # the target is handled by the _solve_adds infeasibility
                # guard at each n'.
                continue
            for j, _w in sup:
                for i in range(max(1, j - k - WINDOW_PAD),
                               min(n, j + k + WINDOW_PAD) + 1):
                    idxs.add(i)
    cands = {xs[i - 1] for i in idxs}
    # region representatives: extremes + a midpoint index per gap
    cands.add(xs[0])
    cands.add(xs[-1])
    sidx = sorted(idxs)
    prev = 1
    for j in sidx + [n]:
        if j - prev > 1:
            cands.add(xs[(prev + j) // 2 - 1])
        prev = j
    cands = sorted(cands)
    if len(cands) > caps["rem_cands_cap"]:
        caps_hit.add("rem_cands_cap")
        cands = cands[:caps["rem_cands_cap"]]
    return cands


def _removal_bases(xs, tg, k_rem, k_total, reg, caps, caps_hit):
    """Yield (removed_values_tuple, reduced_sorted_list)."""
    if k_rem == 0:
        yield (), xs
        return
    avail = Counter(xs)
    cands = _removal_candidates(xs, tg, k_total, reg, caps, caps_hit)
    budget = caps["rem_bases_cap"]
    for rem in combinations_with_replacement(cands, k_rem):
        budget -= 1
        if budget <= 0:
            caps_hit.add("rem_bases_cap")
            return
        c = Counter(rem)
        if any(c[v] > avail[v] for v in c):
            continue
        reduced = list(xs)
        for v in rem:
            reduced.remove(v)
        yield rem, reduced


# ---------------------------------------------------------------------------
# add-only solve on one base
# ---------------------------------------------------------------------------

def _solve_adds(base_sorted, k_add, tg, reg, vlo, vhi, caps, caps_hit,
                max_w):
    """Witness add-multisets of size exactly k_add on the given base.
    Returns list of {"adds": [spec...], "n_new": int} (exact-verified)."""
    # An OrderStat whose index falls outside [1, n'] is unsatisfiable at this
    # array size: no dataset of size n' has that order statistic, so this
    # (k_add, k_rem) combination is INFEASIBLE — an exact skip, not a
    # truncation (no cap note); the k loop proceeds to sizes where the index
    # exists — never a ValueError out of solve() (receipt
    # tests/test_all.py::TestOrderStatOutOfRange).
    n_probe = len(base_sorted) + k_add
    for stat, _pub in tg:
        for comp in _component_quantiles(stat):
            if isinstance(comp, OrderStat):
                try:
                    comp.index(n_probe)
                except ValueError:
                    return []
    if k_add == 0:
        if targets_exact(base_sorted, tg):
            return [{"adds": [], "n_new": len(base_sorted)}]
        return []
    bc = BaseCounts(base_sorted)
    n_new = bc.n + k_add
    found = []
    seen = set()
    scenarios, mean_req = _scenarios(tg, bc, k_add, n_new, reg, vlo, vhi,
                                     caps, caps_hit)
    if not scenarios and mean_req is not None:
        # Mean is the only driving target
        for specs in _mean_only_adds(bc, k_add, mean_req, reg, vlo, vhi):
            values = sorted(s["value"] for s in specs)
            after = sorted(bc.vals + values)
            if targets_exact(after, tg):
                found.append({"adds": specs, "n_new": len(after)})
        return found[:max_w]
    if not scenarios:
        # no driving equalities at all (interval-only targets): probe
        _probe_adds(bc, k_add, tg, reg, vlo, vhi, caps, caps_hit, seen,
                    found, max_w)
        return found[:max_w]
    budget = caps["assembly_budget"]
    for blocks in scenarios:
        if len(found) >= max_w or budget <= 0:
            break
        budget = _assemble(bc, k_add, blocks, reg, vlo, vhi, mean_req, tg,
                           seen, found, max_w, caps, caps_hit, budget)
    return found[:max_w]


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------

def solve(values, targets, register=None, value_lo=None, value_hi=None,
          k_cap=8, rem_cap=0, max_witnesses=6, caps=None):
    """Minimal-k delta families making every target exact.

    values      : released data (Fractions / ints / decimal strings)
    targets     : [(Statistic, published)] — published is an exact value
                  (Fraction/int/decimal string) or an Interval
    register    : value lattice for ADDED entries (e.g. F(1,100) for 2dp
                  data); None = any rational
    value_lo/hi : allowed range for added values (caps on scores)
    k_cap       : max k_add + k_rem searched (loud, in result)
    rem_cap     : max removals per family
    max_witnesses: families returned per minimal k

    Returns {"k": minimal k or None, "families": [...], "caps_hit": [...],
             "exact_already": bool, "n": len(values)}.
    Each family: {"k_add", "k_rem", "removed": [Fraction...],
                  "adds": [{"value", "window", "role"}...], "n_new",
                  "verified": True}.
    Windows are (lo, hi, lo_strict, hi_strict) with None = unbounded: any
    register point inside the window keeps every ORDER-STAT target exact
    (region-count argument; the canonical value is always exact-verified,
    and the tests re-verify an alternate in-window point).
    """
    caps = dict(DEFAULT_CAPS, **(caps or {}))
    reg = as_fraction(register) if register is not None else None
    vlo = as_fraction(value_lo) if value_lo is not None else None
    vhi = as_fraction(value_hi) if value_hi is not None else None
    xs = sorted(as_fraction(v) for v in values)
    tg = [(s, pub if isinstance(pub, Interval) else as_fraction(pub))
          for s, pub in targets]
    caps_hit = set()
    result = {"k": None, "families": [], "caps_hit": caps_hit,
              "exact_already": False, "n": len(xs),
              "k_cap": k_cap, "rem_cap": rem_cap}
    if targets_exact(xs, tg):
        result.update(k=0, exact_already=True,
                      families=[{"k_add": 0, "k_rem": 0, "removed": [],
                                 "adds": [], "n_new": len(xs),
                                 "verified": True}])
        result["caps_hit"] = sorted(caps_hit)
        return result
    for k in range(1, k_cap + 1):
        fams = []
        for k_rem in range(0, min(k, rem_cap) + 1):
            k_add = k - k_rem
            for rem_vals, base in _removal_bases(xs, tg, k_rem, k, reg,
                                                 caps, caps_hit):
                ws = _solve_adds(base, k_add, tg, reg, vlo, vhi, caps,
                                 caps_hit, max_witnesses - len(fams))
                for w in ws:
                    fams.append({"k_add": k_add, "k_rem": k_rem,
                                 "removed": list(rem_vals),
                                 "adds": w["adds"], "n_new": w["n_new"],
                                 "verified": True})
                if len(fams) >= max_witnesses:
                    break
            if len(fams) >= max_witnesses:
                break
        if fams:
            result["k"] = k
            result["families"] = fams
            break
    if result["k"] is None:
        # The k loop exhausted its depth without a witness.  Recorded LOUDLY:
        # an empty result can never certify emptiness beyond the searched
        # depth, so every witness-free search carries this note and the
        # README's caps_hit contract holds for k_cap/rem_cap too (receipt:
        # tests/test_all.py::TestLoudness).
        caps_hit.add("k_cap-exhausted")
    result["caps_hit"] = sorted(caps_hit)
    return result


def _solve_job(job):
    """Worker for solve_many (must be module-level for spawn pickling)."""
    res = solve(job["values"], job["targets"],
                register=job.get("register"),
                value_lo=job.get("value_lo"), value_hi=job.get("value_hi"),
                k_cap=job.get("k_cap", 8), rem_cap=job.get("rem_cap", 0),
                max_witnesses=job.get("max_witnesses", 6),
                caps=job.get("caps"))
    res["job_id"] = job.get("job_id")
    return res


def solve_many(jobs, workers=2):
    """Parallel solve over independent jobs.  spawn context ONLY —
    fork+BLAS deadlocks."""
    if workers <= 1 or len(jobs) <= 1:
        return [_solve_job(j) for j in jobs]
    ctx = mp.get_context("spawn")
    with ctx.Pool(min(workers, len(jobs))) as pool:
        return pool.map(_solve_job, jobs)

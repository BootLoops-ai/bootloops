"""ansatz.py -- faithful port of Blade's BLSearch/PolynomialAnsatz.wl.

Generates the STRUCTURAL search inputs consumed by the C pipeline:
  kinematics/<job>/{kin_common,kin_table}      (reader: src/search/kinematics.c:14)
  <job>/<workid>/config/<k>/{in_nvar,in_var}   (reader: src/search/template.c:105)

Authoritative semantics: BLSearch/PolynomialAnsatz.wl (line numbers cited below).
Authoritative FORMATS/orderings: src/search/{kinematics.c,template.c}.
Emission is byte-identical to the reference trees (fabricator row style via
formats.write_fab_rows).

Ordering contracts ported (each pinned by the search-replay gate + mutation tests):
  O1  kin_table = degree blocks for composite index ii = 0..n ASCENDING
      (toInternalAnsatz, PolynomialAnsatz.wl:65,80: KinTable = Join@@monos).
  O2  within a degree block: FrobeniusSolve order = ascending lexicographic in
      the flattened variable-group order (PolynomialAnsatz.wl:54).
  O3  multi-group: Outer[Join, sols_1, ..., sols_g] -- group 1 varies SLOWEST,
      row = concat of per-group exponent tuples (PolynomialAnsatz.wl:54).
  O4  kin_table columns = BLSearchParameter order (NOT the vars flatten order);
      vars are scattered into their BLSearchParameter positions, parameters not
      in vars get all-zero columns ("conservative", PolynomialAnsatz.wl:62-72).
  O5  composite degree index = mixed radix, group 1 = LEAST significant digit,
      radix cut[i]+1 for group i+1, last group unbounded
      (fromDigits/integerDigits, PolynomialAnsatz.wl:44,48).
  O6  config k is CUMULATIVE: blocks of composite level s for all s digit-wise
      dominated by level n_k, in DESCENDING numeric s (distributeAnsatz
      PolynomialAnsatz.wl:153-155: sub tuples <= ansatz, Sort[#,Greater]&);
      all-"uniform" intmode instead takes ALL s = n_k..0 (PolynomialAnsatz.wl:115-118).
  O7  within a level: localIntegralID ascending (Range order after Select,
      PolynomialAnsatz.wl:141-147); within an integral: monomials in kin_table
      order of its exact-degree block.
  O8  in_var rows = {localIntegralID, kinTableRowID} both 0-based
      (validateAnsatz subtracts {0,1} from WL's 1-based KinTableID,
      PolynomialAnsatz.wl:214); C absorb heuristic detects the leading-degree
      block boundary as the first DECREASE of localIntegralID (template.c:149).
  O9  "dimension" intmode: integral weight = integralDimension, normalized by
      subtracting the min OVER THE BLOCK's integrals (PolynomialAnsatz.wl:122-124
      -- per work/part, NOT global); the weight is subtracted from EVERY
      "dimension" group's degree budget (PolynomialAnsatz.wl:126-128).
  O10 integralDimension of a bare BL integral = -(sum of propagator powers)
      (PolynomialAnsatz.wl:92).  Dressed/extended integrals need the caller to
      supply weights (USInts/EXInts substitution is upstream-side).

CAVEAT: the multi-group (cut != []) orderings O3/O5 and the all-uniform branch
of O6 have NO banked upstream artifact; they are pinned only by hand-derived
unit expectations in selftest_multigroup() (transcribed from the .wl formulas).
Single-group orderings are gated byte-identical against phase0 (db) and phase2
(dbox) banked trees.
"""
from __future__ import annotations

import itertools
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .formats import AnsatzConfig, Kinematics

__all__ = [
    "frobenius_solve", "from_digits", "integer_digits", "polynomial_ansatz",
    "build_kinematics", "distribute_ansatz", "integral_dimension",
    "FamilyAnsatz", "selftest_multigroup",
]


# --------------------------------------------------------------------------
# Basic combinatorics
# --------------------------------------------------------------------------

def frobenius_solve(weights: Sequence[int], total: int) -> List[Tuple[int, ...]]:
    """All nonnegative integer solutions of sum_i weights[i]*e_i == total,
    in Mathematica FrobeniusSolve canonical order = ascending lexicographic.
    (PolynomialAnsatz.wl:54 relies on this order for the kin_table rows.)"""
    if total < 0:
        return []
    k = len(weights)
    out: List[Tuple[int, ...]] = []

    def rec(i: int, rem: int, prefix: Tuple[int, ...]) -> None:
        if i == k - 1:
            w = weights[i]
            if rem % w == 0:
                out.append(prefix + (rem // w,))
            return
        w = weights[i]
        for e in range(rem // w + 1):          # ascending (O2)
            rec(i + 1, rem - e * w, prefix + (e,))

    if k:
        rec(0, total, ())
    elif total == 0:
        out.append(())
    return out


def from_digits(powers: Sequence[int], cut: Sequence[int]) -> int:
    """PolynomialAnsatz.wl:44-45 fromDigits.  powers = per-group degrees
    (p_1..p_g), cut = (c_1..c_{g-1}) bounds all but the LAST group.  Composite
    index = p_1 + (c_1+1)*(p_2 + (c_2+1)*(... p_g)); group 1 least significant
    (O5).  Any negative component -> -1."""
    if any(p < 0 for p in powers):
        return -1
    val = powers[-1]
    for p, c in zip(reversed(powers[:-1]), reversed(list(cut))):
        val = val * (c + 1) + p
    return val


def integer_digits(n: int, cut: Sequence[int]) -> List[int]:
    """Inverse of from_digits (PolynomialAnsatz.wl:48): composite index ->
    per-group degrees, length len(cut)+1, group 1 first."""
    digits: List[int] = []
    for c in cut:
        digits.append(n % (c + 1))
        n //= (c + 1)
    digits.append(n)
    return digits


def polynomial_ansatz(var_weights: Sequence[Sequence[int]],
                      cut: Sequence[int], n: int) -> List[Tuple[int, ...]]:
    """PolynomialAnsatz.wl:52-55.  Exponent rows (flattened variable order) of
    the exact-composite-degree-n block.  Group 1 solutions vary slowest (O3)."""
    powers = integer_digits(n, cut)
    per = [frobenius_solve(w, p) for w, p in zip(var_weights, powers)]
    return [tuple(e for sol in combo for e in sol)
            for combo in itertools.product(*per)]


def integral_dimension(indices: Sequence[int]) -> int:
    """Bare BL integral branch of integralDimension (PolynomialAnsatz.wl:92):
    -(sum of propagator powers).  Dressed integrals: caller supplies weights."""
    return -sum(indices)


# --------------------------------------------------------------------------
# KinTable (toInternalAnsatz, PolynomialAnsatz.wl:64-84)
# --------------------------------------------------------------------------

def build_kinematics(search_params: Sequence[str],
                     var_groups: Sequence[Sequence[str]],
                     var_weights: Sequence[Sequence[int]],
                     cut: Sequence[int], nmax: int
                     ) -> Tuple[List[Tuple[int, ...]], Dict[int, List[int]]]:
    """Return (kin_table rows, leading_term).  kin_table row j = exponent
    vector of monomial j in BLSearchParameter column order (O4); leading_term[s]
    = 0-based kin_table row ids of the exact-composite-degree-s block."""
    flat = [v for g in var_groups for v in g]
    posi = []
    for v in flat:
        if v not in search_params:
            raise ValueError(f"toInternalAnsatz unknown pattern: {v!r} "
                             f"not in BLSearchParameter {tuple(search_params)}")
        posi.append(search_params.index(v))     # PolynomialAnsatz.wl:67
    npar = len(search_params)

    table: List[Tuple[int, ...]] = []
    leading: Dict[int, List[int]] = {}
    for ii in range(nmax + 1):                  # ascending blocks (O1)
        rows = []
        for flatrow in polynomial_ansatz(var_weights, cut, ii):
            r = [0] * npar
            for p, e in zip(posi, flatrow):     # scatter (O4), wl:72
                r[p] = e
            rows.append(tuple(r))
        leading[ii] = list(range(len(table), len(table) + len(rows)))
        table.extend(rows)
    if len(set(table)) != len(table):
        raise ValueError("KinTable rows not unique; KinTableID would clash")
    return table, leading


# --------------------------------------------------------------------------
# distributeAnsatz (PolynomialAnsatz.wl:112-159)
# --------------------------------------------------------------------------

def distribute_ansatz(ints_weight: Sequence[int], intmode: Sequence[str],
                      cut: Sequence[int], n: int,
                      leading_term: Dict[int, List[int]]
                      ) -> List[List[Tuple[int, int]]]:
    """Cumulative configs for levels 0..n; each config = list of
    (localIntegralID, kinTableRowID), both 0-based (O8)."""
    nints = len(ints_weight)
    ngrp = len(cut) + 1
    if len(intmode) != ngrp:
        raise ValueError("intmode length must equal len(cut)+1 group count")

    if all(m == "uniform" for m in intmode):    # wl:114-119
        out: List[List[Tuple[int, int]]] = []
        prev: List[Tuple[int, int]] = []
        for i in range(n + 1):
            blk = [(k, r) for k in range(nints) for r in leading_term[i]]
            cur = blk + prev                    # new level FIRST (O6)
            out.append(cur)
            prev = cur
        return out

    # "dimension" path: normalize weights over THIS block (wl:122-124, O9)
    wmin = min(ints_weight)
    weight = [w - wmin for w in ints_weight]
    dim_groups = [g for g, m in enumerate(intmode) if m == "dimension"]  # wl:126
    power = [[weight[k] if g in dim_groups else 0 for g in range(ngrp)]
             for k in range(nints)]             # wl:127-128

    lt: Dict[int, List[Tuple[int, int]]] = {}
    out = []
    for i in range(n + 1):
        ansatz = integer_digits(i, cut)         # wl:134
        entries: List[Tuple[int, int]] = []
        for k in range(nints):                  # localid ascending (O7)
            comp = from_digits([a - p for a, p in zip(ansatz, power[k])], cut)
            if comp < 0:                        # wl:144
                continue
            entries.extend((k, r) for r in leading_term[comp])  # wl:147
        lt[i] = entries
        # dominated sub-levels, descending composite (wl:153-155, O6)
        sub = sorted((from_digits(tpl, cut)
                      for tpl in itertools.product(*(range(a + 1) for a in ansatz))),
                     reverse=True)
        out.append([e for s in sub for e in lt[s]])
    return out


# --------------------------------------------------------------------------
# Family-level driver
# --------------------------------------------------------------------------

@dataclass
class FamilyAnsatz:
    """One family's ansatz recipe = the arguments of BLGenerateAnsatzForJob
    (PolynomialAnsatz.wl:238) plus the BLSearchParameter order."""
    search_params: Tuple[str, ...]              # BLSearchParameter order
    var_groups: List[List[str]]                 # 'vars'
    var_weights: List[List[int]]                # 'varweight'
    intmode: List[str]                          # per group
    cut: List[int]                              # len = groups-1
    nmax: int                                   # 'n'
    _kin: Optional[Tuple[List[Tuple[int, ...]], Dict[int, List[int]]]] = field(
        default=None, repr=False)

    def _built(self):
        if self._kin is None:
            self._kin = build_kinematics(self.search_params, self.var_groups,
                                         self.var_weights, self.cut, self.nmax)
        return self._kin

    def kinematics(self) -> Kinematics:
        table, _ = self._built()
        return Kinematics(len(table), len(self.search_params),
                          [list(r) for r in table])

    def kin_prefix(self, nmax: int) -> Kinematics:
        """Smaller-degree variant: blocks 0..nmax are a PREFIX of the full
        table (O1), so row ids in any config with level <= nmax stay valid."""
        if nmax > self.nmax:
            raise ValueError("prefix beyond built nmax")
        table, leading = self._built()
        upto = leading[nmax][-1] + 1 if leading[nmax] else 0
        return Kinematics(upto, len(self.search_params),
                          [list(r) for r in table[:upto]])

    def block_configs(self, ints_weight: Sequence[int],
                      levels: Optional[Sequence[int]] = None
                      ) -> List[AnsatzConfig]:
        """Configs for ONE block (work) given its integrals' RAW dimensions
        (integral_dimension values; normalization happens inside, per block).
        levels: config k = cumulative ansatz at levels[k]; default = 0..nmax
        (upstream writes every level, writeAnsatzForJob wl:206-208; the
        escalation driver may skip levels)."""
        if levels is None:
            levels = list(range(self.nmax + 1))
        if max(levels) > self.nmax:
            raise ValueError("level beyond kin_table nmax")
        _, leading = self._built()
        allcfg = distribute_ansatz(ints_weight, self.intmode, self.cut,
                                   max(levels), leading)
        return [AnsatzConfig(list(allcfg[lv])) for lv in levels]

    def write_job(self, search_dir: str, job: str,
                  blocks: Sequence[Tuple[Sequence[int], Optional[Sequence[int]]]]
                  ) -> None:
        """Emit kinematics/<job>/ + <job>/<workid>/config/<k>/ for each block
        (writeKinematics wl:226-230 + writeAnsatzForJob wl:191-209).
        blocks[workid] = (ints_weight, levels_or_None)."""
        self.kinematics().write(os.path.join(search_dir, "kinematics", job))
        for workid, (w, levels) in enumerate(blocks):
            for k, cfg in enumerate(self.block_configs(w, levels)):
                cfg.write(os.path.join(search_dir, job, str(workid),
                                       "config", str(k)))


# --------------------------------------------------------------------------
# Hand-derived multi-group pins (no banked upstream artifact exists for
# cut != []; expectations transcribed by hand from the .wl formulas)
# --------------------------------------------------------------------------

def selftest_multigroup() -> None:
    # -- O5 mixed radix: vars {{x},{y}}, cut={2}: composite = px + 3*py
    assert integer_digits(7, [2]) == [1, 2]
    assert from_digits([1, 2], [2]) == 7
    assert from_digits([2, -1], [2]) == -1
    tab, lead = build_kinematics(("x", "y"), [["x"], ["y"]], [[1], [1]], [2], 7)
    assert tab == [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2)]

    # -- O6/O9 dimension path, mixed intmode, 2 ints weights (0,1):
    cfg = distribute_ansatz([0, 1], ["uniform", "dimension"], [2], 7, lead)
    assert cfg[7] == [(0, 7), (1, 4), (0, 6), (1, 3), (0, 4), (1, 1),
                      (0, 3), (1, 0), (0, 1), (0, 0)]
    # composite 2 = (2,0) is NOT digit-dominated by 4 = (1,1): absent at lvl 4
    assert all(r != 2 for _, r in cfg[4])

    # -- O6 all-uniform branch DOES include non-dominated level 2 at lvl 4
    cfgu = distribute_ansatz([0], ["uniform", "uniform"], [2], 4, lead)
    assert cfgu[4] == [(0, 4), (0, 3), (0, 2), (0, 1), (0, 0)]

    # -- O3 Outer order + O4 scatter with permuted BLSearchParameter
    tab2, _ = build_kinematics(("y1", "x1", "y2", "x2"),
                               [["x1", "x2"], ["y1", "y2"]],
                               [[1, 1], [1, 1]], [1], 3)
    assert tab2 == [(0, 0, 0, 0),
                    (0, 0, 0, 1), (0, 1, 0, 0),
                    (0, 0, 1, 0), (1, 0, 0, 0),
                    (0, 0, 1, 1), (1, 0, 0, 1), (0, 1, 1, 0), (1, 1, 0, 0)]

    # -- O2 general weights: FrobeniusSolve[{1,2},4] ascending lex
    assert frobenius_solve([1, 2], 4) == [(0, 2), (2, 1), (4, 0)]


if __name__ == "__main__":
    selftest_multigroup()
    print("ansatz.py multigroup selftest: OK")

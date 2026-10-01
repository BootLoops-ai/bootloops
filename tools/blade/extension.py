"""blade.extension — Blade's integral-extension combinatorics, ported from .wl.

Every method cites the .wl lines it ports (Blade, MIT; notice retained):

  blade/BLSearch/IntegralExtension.wl
      sectorDown            :34-39   sectorDownCollect      :41-43
      globalExtendG         :50-56   OperatorExtendG        :63-76
      BLGenerateScheme      :84-90   getTopSector           :97
      divideSubFamilyBasic  :108-130
  blade/BLSearch/Integrals.wl
      BLIntPropagators :31  BLIntDots :32  BLIntRank :33
      BLSectorOrSubsectorQ :36  BLSector :40  BLDenominator :41
      BLLowerDenom :44  IntegralWeight :57-62  BLSortIntegrals :66
      BLEliminateZeroSectors :70-73 (reads zero sectors from
      $MaximalCutDirectory/results/zerosectors — the .wl does NOT compute
      them either; callers supply them here, exactly as upstream supplies a
      file produced by the maximal-cut stage)
  blade/BLGeneral/AuxiliaryFunctions.wl
      frobeniusSolve2 :92-95 (WL FrobeniusSolve on unit weights = all
      nonnegative compositions in ascending lexicographic order; this is
      the order contract of the KinTable generation)

Defaults transcribed from blade/BLSearch/SearchOptions.wl:149
  Options[BLSetSchemeOptions] = {"OperatorExtendGQ"->True,
  "MinimalSchemeRank"->1, ...}  and :85 "IntegralOrdering"->1.

Representation: an integral is a plain tuple of ints (the BL[fam,{...}]
power list; the family head is implicit — one family per Extension instance,
like the per-family upstream globals).  A "sector"/"denominator" is the same
kind of tuple.  BLExtraIntDerivDen general-integrand branches
(IntegralExtension.wl:42-43) are NOT ported: standard Feynman families have
Length[BLExtraIntDerivDen]==1 (DefineFamily.wl:91) so the branch is dead for
them; a guard raises if a caller ever needs it.

Ordering fidelity note (Integrals.wl:57-62): IntegralWeight keys 4-7 are
(!ExtMappedQ[int], !MappedQ[int], Min[powers], int).  In the search context
ExtMappedQ/MappedQ are undefined symbols (comment :65 "useless"), so WL's
SortBy tie-breaks by the canonical order of the unevaluated expressions,
which for same-family integrals reduces to canonical (elementwise) order of
the power lists.  We therefore tie-break by the power tuple directly; keys
6-7 can then never fire and are omitted.
"""

from __future__ import annotations

from itertools import product
from typing import Dict, FrozenSet, Iterable, List, Sequence, Tuple

Ints = Tuple[int, ...]


def props(v: Ints) -> int:
    """BLIntPropagators (Integrals.wl:31): number of positive powers."""
    return sum(1 for a in v if a > 0)


def dots(v: Ints) -> int:
    """BLIntDots (Integrals.wl:32): Total@Select[v-1, #>0&]."""
    return sum(a - 1 for a in v if a > 1)


def rank(v: Ints) -> int:
    """BLIntRank (Integrals.wl:33): -Total@Select[v, #<0&]."""
    return -sum(a for a in v if a < 0)


def sector(v: Ints) -> Ints:
    """BLSector (Integrals.wl:40)."""
    return tuple(1 if a > 0 else 0 for a in v)


def denominator(v: Ints) -> Ints:
    """BLDenominator (Integrals.wl:41): keep positive powers, zero the rest."""
    return tuple(a if a > 0 else 0 for a in v)


def lower_denoms(v: Ints) -> List[Ints]:
    """BLLowerDenom (Integrals.wl:44): Tuples[Range[0,Max[0,a]]&/@v].
    WL Tuples iterates the LAST slot fastest == itertools.product."""
    return [t for t in product(*[range(0, max(0, a) + 1) for a in v])]


def sector_or_subsector_q(s1: Ints, s2: Ints) -> bool:
    """BLSectorOrSubsectorQ (Integrals.wl:36): no component of s1-s2 positive."""
    return not any(a - b > 0 for a, b in zip(s1, s2))


def get_top_sector(ints: Sequence[Ints]) -> Ints:
    """getTopSector (IntegralExtension.wl:97): componentwise Max>0 flag."""
    return tuple(1 if max(col) > 0 else 0 for col in zip(*ints))


def integrals_in(ints: Iterable[Ints]) -> List[Ints]:
    """BLIntegralsIn (Integrals.wl:67): DeleteDuplicates, keep first order."""
    seen, out = set(), []
    for v in ints:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def sort_key(v: Ints, ordering: int = 1):
    """IntegralWeight[ordering] (Integrals.wl:57-62); tie-break by the power
    tuple (see module docstring on the ExtMappedQ/MappedQ keys)."""
    p, d, r = props(v), dots(v), rank(v)
    if ordering == 1:
        head = (-p, -d, -r)
    elif ordering == 2:
        head = (-p, -r, -d)
    elif ordering == 3:
        head = (-p, -d - r, -d, -r)
    elif ordering == 4:
        head = (-p, -d - r, -r, -d)
    else:
        raise ValueError(f"unknown integral ordering {ordering}")
    return head + (v,)


def frobenius_solve_unit(nslots: int, n: int) -> List[Ints]:
    """frobeniusSolve2[ConstantArray[1,nslots], n] (AuxiliaryFunctions.wl:92-95):
    {} weights: {{}} if n==0 else {}; else WL FrobeniusSolve = all nonnegative
    solutions of Total==n in ascending lexicographic order."""
    if nslots == 0:
        return [()] if n == 0 else []
    if nslots == 1:
        return [(n,)]
    out: List[Ints] = []
    for first in range(n + 1):
        out += [(first,) + rest for rest in frobenius_solve_unit(nslots - 1, n - first)]
    return out


class Extension:
    """Port of the IntegralExtension.wl staircase for ONE family.

    zero_sectors plays the role of results/zerosectors read by
    BLEliminateZeroSectors (Integrals.wl:71); upstream computes it in the
    maximal-cut stage, we accept it as input."""

    def __init__(self, minimal_scheme_rank: int = 1, operator_extend_gq: bool = True,
                 integral_ordering: int = 1,
                 zero_sectors: Iterable[Ints] = ()) -> None:
        # defaults: SearchOptions.wl:149 and :85
        self.minimal_scheme_rank = minimal_scheme_rank
        self.operator_extend_gq = operator_extend_gq
        self.integral_ordering = integral_ordering
        self.zero_sectors: FrozenSet[Ints] = frozenset(tuple(s) for s in zero_sectors)

    # -- Integrals.wl -------------------------------------------------------
    def sort_integrals(self, ints: Iterable[Ints]) -> List[Ints]:
        """BLSortIntegrals (Integrals.wl:66)."""
        return sorted(ints, key=lambda v: sort_key(v, self.integral_ordering))

    def eliminate_zero_sectors(self, ints: Iterable[Ints]) -> List[Ints]:
        """BLEliminateZeroSectors (Integrals.wl:70-73): zero-sector integrals
        become 0 and drop out of any subsequent BLIntegralsIn."""
        return [v for v in ints if sector(v) not in self.zero_sectors]

    # -- IntegralExtension.wl ------------------------------------------------
    def sector_down(self, sect: Ints, n: int) -> List[Ints]:
        """sectorDown (IntegralExtension.wl:34-39): rank-n integrals of the
        generalized sector `sect` (dots allowed, no ISP): distribute -n over
        the zero slots in frobeniusSolve order."""
        posi = [i for i, a in enumerate(sect) if a == 0]
        out = []
        for sol in frobenius_solve_unit(len(posi), n):
            tmp = list(sect)
            for i, e in zip(posi, sol):
                tmp[i] = -e
            out.append(tuple(tmp))
        return out

    def sector_down_collect(self, sect: Ints, n: int) -> List[Ints]:
        """sectorDownCollect (IntegralExtension.wl:41-43): ranks n,n-1,...,0
        (Reverse[Range[0,n]]).  General-integrand branch (:42-43) not ported —
        guard lives in the caller contract (single integrand)."""
        out: List[Ints] = []
        for m in range(n, -1, -1):
            out += self.sector_down(sect, m)
        return out

    def global_extend_g(self, target: Sequence[Ints]) -> List[Tuple[Ints, int]]:
        """globalExtendG (IntegralExtension.wl:50-56): every lower denominator
        of every target denominator, with rank Max[$MinimalSchemeRank, ranks
        of targets in the SAME SECTOR]; keep rank>=0.  Union sorts denominators
        canonically (ascending lex for same-length tuples)."""
        gpsec: Dict[Ints, List[Ints]] = {}
        for t in target:
            gpsec.setdefault(sector(t), []).append(t)

        def getrank(bas: Ints) -> int:
            return max([self.minimal_scheme_rank]
                       + [rank(t) for t in gpsec.get(sector(bas), [])])

        basis = integrals_in(
            ld for d in sorted(set(denominator(t) for t in target))
            for ld in lower_denoms(d))
        return [(b, getrank(b)) for b in basis if getrank(b) >= 0]

    def operator_extend_g(self, target: Sequence[Ints]) -> List[Tuple[Ints, int]]:
        """OperatorExtendG (IntegralExtension.wl:63-76): the r/p_{+d} staircase.
        For each target denominator with rank n, every lower denominator ld
        gets rank props(ld) - (props(den) - n); unionscheme (:67-68) keeps the
        max rank per denominator; keep rank>=0."""
        gpden: Dict[Ints, List[Ints]] = {}
        for t in target:
            gpden.setdefault(denominator(t), []).append(t)

        def getrank(bas: Ints) -> int:
            return max([-1] + [rank(t) for t in gpden.get(denominator(bas), [])])

        def ncircleddash(sect: Ints, n: int) -> List[Tuple[Ints, int]]:
            return [(ld, props(ld) - (props(sect) - n)) for ld in lower_denoms(sect)]

        def unionscheme(s1, s2):
            if not s1:                      # unionscheme[{},a_]:=a  (:67)
                return list(s2)
            merged, best = [], {}
            for pair in s1 + s2:            # DeleteDuplicates + GatherBy (:68)
                if pair in best.setdefault(pair[0], set()):
                    continue
                best[pair[0]].add(pair)
            order = integrals_in(p[0] for p in s1 + s2)
            for den in order:   # SortBy[{-rank}][[1]] = max rank; ties are
                # impossible after DeleteDuplicates (same den+rank = same pair)
                merged.append(max(best[den], key=lambda q: q[1]))
            return merged

        basis = sorted(set(denominator(t) for t in target))   # Union (:70)
        scheme: List[Tuple[Ints, int]] = []
        for b in basis:                                        # Do (:72-74)
            scheme = unionscheme(scheme, ncircleddash(b, getrank(b)))
        return [p for p in scheme if p[1] >= 0]                # Select (:75)

    def generate_scheme(self, target: Sequence[Ints]) -> List[Ints]:
        """BLGenerateScheme (IntegralExtension.wl:84-90): Join the two schemes
        (OperatorExtendG only if $OperatorExtendGQ), expand each (denominator,
        rank) with sectorDownCollect, DeleteDuplicates, assert the targets are
        covered (:88)."""
        target = list(target)
        scheme = (self.operator_extend_g(target) if self.operator_extend_gq else [])
        scheme = scheme + self.global_extend_g(target)
        full = integrals_in(v for (sect, n) in scheme
                            for v in self.sector_down_collect(sect, n))
        if not set(target) <= set(full):
            raise RuntimeError(
                "error: current scheme is not sufficient to reduce target "
                "(IntegralExtension.wl:88)")
        return full

    def divide_sub_family_basic(self, target: Sequence[Ints],
                                dividelevel=float("inf")
                                ) -> List[Tuple[Ints, List[Ints]]]:
        """divideSubFamilyBasic (IntegralExtension.wl:108-130) -> list of
        (top sector, extended integral list) per sub-family.  Upstream stores
        them in tarinfo[sysid, "top"/"target"]; Counter == len(result).
        GenerateScheme input filter and Complement-resort per :113-118."""
        ints = self.sort_integrals(integrals_in(self.eliminate_zero_sectors(target)))
        out: List[Tuple[Ints, List[Ints]]] = []
        if not ints:
            return out
        toplevel = props(ints[0])
        while ints and toplevel - props(ints[0]) < dividelevel:
            top = sector(ints[0])
            r0, d0 = rank(ints[0]), dots(ints[0])
            sel = [v for v in ints
                   if sector(v) == top and rank(v) <= r0 and dots(v) <= d0]
            ext = integrals_in(self.eliminate_zero_sectors(self.generate_scheme(sel)))
            out.append((top, ext))
            covered = set(ext)
            ints = self.sort_integrals(v for v in ints if v not in covered)
        if ints:                                     # remainder block (:121-124)
            ext = integrals_in(self.eliminate_zero_sectors(self.generate_scheme(ints)))
            out.append((get_top_sector(ints), ext))
        return out

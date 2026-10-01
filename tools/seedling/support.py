"""Target-support measurement.

Kira conventions: for an index vector nu over n_sp scalar products,
  t = #(nu_i > 0)            (propagators on; sector id = sum 2^i over nu_i>0)
  r = sum of positive nu_i   (>= t)
  s = sum |negative nu_i|    (ISP degree)
  d = r - t                  (dots)

Target file format (kira select_mandatory_list): lines "family[i1, i2, ...]",
blank lines ignored.

Scope: global and per-sector support maxima, plus sector-tree closure of the
target list under family symmetry maps (`symmetry_closure`; consumes a kira
sectorSymmetries/sectorRelations file via rankcheck.parse_sector_maps).
The global-maxima path matches the deployed manual pin recipe (golden-case
validated); closure feeds the same support_table/global_support consumers.
"""
from dataclasses import dataclass
import re

from . import rankcheck

_LINE = re.compile(r"^\s*([A-Za-z_]\w*)\s*\[([^\]]*)\]\s*$")


@dataclass(frozen=True)
class Target:
    family: str
    indices: tuple

    @property
    def sector(self):
        return sum(1 << i for i, n in enumerate(self.indices) if n > 0)

    @property
    def t(self):
        return sum(1 for n in self.indices if n > 0)

    @property
    def r(self):
        return sum(n for n in self.indices if n > 0)

    @property
    def s(self):
        return sum(-n for n in self.indices if n < 0)

    @property
    def d(self):
        return self.r - self.t


def parse_targets(path):
    """Return list[Target]. Raises ValueError on malformed non-blank lines."""
    out = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip():
                continue
            m = _LINE.match(line)
            if not m:
                raise ValueError(f"{path}:{ln}: unparseable target line: {line!r}")
            fam, body = m.group(1), m.group(2)
            idx = tuple(int(x) for x in body.split(",")) if body.strip() else ()
            out.append(Target(fam, idx))
    return out


def symmetry_closure(targets, maps):
    """Close a target list under family symmetry maps (sector-tree closure).

    maps: rankcheck.parse_sector_maps output ({src, tgt, perm} dicts from a
    kira sectorSymmetries/sectorRelations file). Every map whose src sector
    matches a target's sector is applied FORWARD (the direction kira applies
    it — the reduction demand lands in the image sector), iterated to a
    fixpoint: the returned list holds the input targets plus every reachable
    image, so support_table/global_support on it give the CLOSED support and
    mapped-to sectors get target-grade pin cells instead of closure cells.

    Images the map file leaves undecidable (rankcheck.apply_perm -> None:
    a nonzero index at an unmapped position) are REPORTED, never silently
    dropped — over-pinning is safe, a silently missing image sector is not.
    An image whose computed sector differs from the map's recorded tgt is
    added under its computed sector (the only reading consistent with its
    indices) and reported as a sector_mismatch (inconsistent map line).

    Returns (closed_targets, report). report keys: n_input, n_maps, n_added,
    added ([family, indices, sector] per new target), undecidable
    ({indices, src, tgt} per unresolvable image), sector_mismatch
    ({indices, image, map_tgt, image_sector}).
    """
    closed = list(targets)
    seen = {(t.family, t.indices) for t in closed}
    report = {"n_input": len(closed), "n_maps": len(maps), "n_added": 0,
              "added": [], "undecidable": [], "sector_mismatch": []}
    frontier = list(closed)
    while frontier:
        nxt = []
        for tg in frontier:
            for mp in maps:
                if mp["src"] != tg.sector:
                    continue
                if len(mp["perm"]) != len(tg.indices):
                    raise ValueError(
                        f"symmetry_closure: map perm length "
                        f"{len(mp['perm'])} != index length "
                        f"{len(tg.indices)} (map {mp['src']} -> {mp['tgt']})")
                img = rankcheck.apply_perm(tg.indices, mp["perm"])
                if img is None:
                    rec = {"indices": list(tg.indices), "src": mp["src"],
                           "tgt": mp["tgt"]}
                    if rec not in report["undecidable"]:
                        report["undecidable"].append(rec)
                    continue
                new = Target(tg.family, img)
                if new.sector != mp["tgt"]:
                    rec = {"indices": list(tg.indices), "image": list(img),
                           "map_tgt": mp["tgt"], "image_sector": new.sector}
                    if rec not in report["sector_mismatch"]:
                        report["sector_mismatch"].append(rec)
                key = (new.family, new.indices)
                if key in seen:
                    continue
                seen.add(key)
                closed.append(new)
                report["added"].append([new.family, list(img), new.sector])
                nxt.append(new)
        frontier = nxt
    report["n_added"] = len(report["added"])
    return closed, report


def support_table(targets):
    """Per-sector support: dict sector_id -> dict(r=, s=, d=, t=, n=count)."""
    tab = {}
    for tg in targets:
        e = tab.setdefault(tg.sector, {"r": 0, "s": 0, "d": 0, "t": tg.t, "n": 0})
        e["r"] = max(e["r"], tg.r)
        e["s"] = max(e["s"], tg.s)
        e["d"] = max(e["d"], tg.d)
        e["n"] += 1
    return tab


def global_support(targets):
    """Global maxima: dict(r=, s=, d=, t=, n=, families=list)."""
    if not targets:
        raise ValueError("empty target list")
    return {
        "r": max(t.r for t in targets),
        "s": max(t.s for t in targets),
        "d": max(t.d for t in targets),
        "t": max(t.t for t in targets),
        "n": len(targets),
        "families": sorted({t.family for t in targets}),
    }

"""Independent master-basis rank check.

Motivating exhibits (both measured, opposite signs, same silence):
  - a basis SHORT by 2 — an S4 symmetry kira's detection misses; certified
    vs a Bessel-moment oracle to 1e-70. Deficits from MISSED symmetries cannot be
    found from kira's own maps (kira missed them) — they need an EXTERNAL
    symmetry oracle (`seedling audit`/Feynson-class), injected here as a callable.
  - a basis LONG by 1 — reducible master retained, selection-dependent.
    Redundancy under kira's OWN maps IS detectable file-side: two masters
    in one sector related by a self-symmetry permutation.
No production engine performs either check; the solver emits a consistent-looking
table in both failure modes.

Checks (any-of):
  1. expected-count comparison (user-supplied: Azurite-class count, literature,
     or syzygy dimension via the certify hook)
  2. self-symmetry redundancy scan over kira's sectormappings (detects
     retained-redundant-master LONG bases)
  3. external-symmetry-oracle hook (detects missed-symmetry SHORT bases)

kira file formats parsed (observed, kira 3.x):
  masters block in kira.log:  'fam[i1,...,in]  # <sector>'
  sectorSymmetries/sectorRelations line: '<src> placeholder*6 1 <tgt> <t> 0
     {mom-rules} {aux-rule} <n_props image positions, -1=off> 1 0'
"""
import re

_MASTER = re.compile(r"^\s*([A-Za-z_]\w*)\[([-\d,\s]+)\]\s*#\s*(\d+)\s*$")


def parse_kira_masters(log_path):
    """Extract the master list from kira.log. Returns list of (indices, sector)."""
    masters, inblock = [], False
    for line in open(log_path):
        if "Number of master integrals:" in line:
            inblock = True
            continue
        if inblock:
            m = _MASTER.match(line)
            if m:
                idx = tuple(int(x) for x in m.group(2).split(","))
                masters.append((idx, int(m.group(3))))
            elif masters:
                break
    return masters


def parse_sector_maps(path, n_props):
    """Parse kira sectorSymmetries/sectorRelations. Returns list of dicts
    {src, tgt, perm} where perm[i] = image position of propagator i (-1 = off)."""
    maps = []
    for line in open(path):
        if "{" not in line:
            continue
        head, rest = line.split("{", 1)
        toks = head.split()
        try:
            src, tgt = int(toks[0]), int(toks[8])
        except (ValueError, IndexError):
            continue
        tail = rest.rsplit("}", 1)[-1].split()
        ints = [int(x) for x in tail if re.fullmatch(r"-?\d+", x)]
        if len(ints) < n_props + 2:
            continue
        perm = ints[:n_props]           # last two tokens are flags
        maps.append({"src": src, "tgt": tgt, "perm": perm})
    return maps


def apply_perm(indices, perm):
    """Apply a propagator-position permutation to an index vector.

    Returns None when the image is not determined by the map: kira's sector maps
    record image positions only for the on-propagators of the source sector
    (perm[i] = -1 elsewhere), so any nonzero index at an unmapped position —
    an on-propagator mapped off OR an ISP power whose image scalar product is
    not tracked — makes the image undecidable from this file alone. (Treating
    those as zero fabricates false redundancy pairs, e.g. corner vs ISP master.)
    """
    out = [0] * len(indices)
    for i, p in enumerate(perm):
        if p >= 0:
            out[p] = indices[i]
        elif indices[i] != 0:
            return None
    return tuple(out)


def redundancy_scan(masters, maps):
    """#11-class check: masters related by a SELF-symmetry of their own sector.
    Returns list of (master_a, master_b, src_sector) suspicious pairs."""
    by_sector = {}
    for idx, sec in masters:
        by_sector.setdefault(sec, set()).add(idx)
    pairs = []
    for mp in maps:
        if mp["src"] != mp["tgt"]:
            continue
        pool = by_sector.get(mp["src"], ())
        for idx in pool:
            img = apply_perm(idx, mp["perm"])
            if img and img != idx and img in pool:
                a, b = sorted((idx, img))
                if (a, b, mp["src"]) not in pairs:
                    pairs.append((a, b, mp["src"]))
    return pairs


def rank_report(masters, expected_count=None, maps=None, external_orbits=None):
    """Assemble the rank verdict.

    expected_count: independent master count (Azurite/syzygy/literature).
    maps: parse_sector_maps output (enables redundancy scan).
    external_orbits: iterable of sets of index tuples from an EXTERNAL symmetry
      oracle; two masters in one orbit = redundancy kira cannot see; an orbit
      whose members are all absent is fine; missed-symmetry deficits show up as
      orbit-inconsistent value structure and need the oracle by construction.
    Verdict: OK / LONG / SHORT / SUSPECT with evidence lists.
    """
    n = len(masters)
    report = {"found": n, "expected": expected_count, "verdict": "OK",
              "redundant_pairs": [], "external_redundant": []}
    if maps:
        report["redundant_pairs"] = redundancy_scan(masters, maps)
    if external_orbits:
        have = {idx for idx, _ in masters}
        for orbit in external_orbits:
            hit = sorted(set(orbit) & have)
            if len(hit) > 1:
                report["external_redundant"].append(hit)
    if expected_count is not None and n < expected_count:
        report["verdict"] = "SHORT"          # missed-symmetry class
    elif (expected_count is not None and n > expected_count) or \
            report["redundant_pairs"] or report["external_redundant"]:
        report["verdict"] = "LONG"           # master-#11 class
    elif expected_count is None and not maps and not external_orbits:
        report["verdict"] = "SUSPECT"        # nothing independent supplied
    return report

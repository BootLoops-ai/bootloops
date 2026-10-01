# topos.py -- pinned 15-topology enumeration (EXACTLY sd_engine.py order) +
# per-topology exponent vectors over the 16 shared pattern VALUES.
# Value index y = int("".join(bits), 2) in TREE-SLOT order of the evaluator:
#   caterpillar base (((s1,s2),s3),s4): ab_fiber pats (a,b,c,d) = slots 1..4
#   balanced case Q(old1, old2, young1, young2): bal_fiber pats = slots 1..4
# Two INDEPENDENT constructions of the exponent vectors (index arithmetic vs
# string permutation) are asserted equal on import.
import itertools, json, os

BASE = os.environ.get("POSQ_STAGE0", os.path.dirname(os.path.abspath(__file__)))
# The production pattern-count data is NOT shipped; the topology enumeration
# below needs no data file. Point POSQ_STAGE0 at a directory containing
# QUARTET_COLLAPSE.json to enable the exponent-vector functions (else they
# refuse by name).
_QC_PATH = f"{BASE}/QUARTET_COLLAPSE.json"
QC = json.load(open(_QC_PATH)) if os.path.exists(_QC_PATH) else None
LEAVES = "ABCD"


def _require_counts():
    if QC is None:
        raise SystemExit(
            "REFUSE: POSQ-PRODUCTION-DATA-ABSENT -- QUARTET_COLLAPSE.json "
            "is not shipped; point POSQ_STAGE0 at a directory containing it "
            "to use the per-topology exponent vectors")
    return PAT_COUNTS

TOPOLOGIES = []                     # pinned order: 12 caterpillars, 3 balanced
for c1, c2 in itertools.combinations(LEAVES, 2):
    rest = [x for x in LEAVES if x not in (c1, c2)]
    for k in (0, 1):
        TOPOLOGIES.append(("cat", (c1, c2, rest[k], rest[1 - k])))
TOPOLOGIES += [("bal", ("A", "B", "C", "D")), ("bal", ("A", "C", "B", "D")),
               ("bal", ("A", "D", "B", "C"))]
assert len(TOPOLOGIES) == 15
CAT_IDX = list(range(12))
BAL_IDX = [12, 13, 14]

def topo_name(t):
    kind, (a, b, c, d) = t
    return f"((({a},{b}),{c}),{d})" if kind == "cat" else f"(({a},{b}),({c},{d}))"

# split classes (derive_sentences.py law, same index sets)
SPLIT_CLASSES = {
    "AB|CD": {0, 1, 10, 11, 12},
    "AC|BD": {2, 3, 8, 9, 13},
    "AD|BC": {4, 5, 6, 7, 14},
}

# observed pattern counts keyed by (bA,bB,bC,bD) bit tuple (None when the
# production data file is absent — see _require_counts above)
PAT_COUNTS = (None if QC is None else
              {tuple(int(ch) for ch in p): int(n)
               for p, n in QC["pattern_counts"].items()})
if PAT_COUNTS is not None:
    assert sum(PAT_COUNTS.values()) == 318 and len(PAT_COUNTS) == 14

def _slots(perm):
    """leaf letters -> tuple index positions (A=0..D=3)."""
    return tuple(LEAVES.index(ch) for ch in perm)

def expvec_cat(ti):
    """16-vector: exponent of value y for caterpillar topology ti."""
    _require_counts()
    kind, perm = TOPOLOGIES[ti]
    assert kind == "cat"
    sl = _slots(perm)
    v = [0] * 16
    for x, n in PAT_COUNTS.items():
        y = (x[sl[0]] << 3) | (x[sl[1]] << 2) | (x[sl[2]] << 1) | x[sl[3]]
        v[y] += n
    return v

def expvec_bal_cases(ti):
    """two 16-vectors (case1: first cherry older; case2: second older)."""
    _require_counts()
    kind, (a, b, c, d) = TOPOLOGIES[ti]
    assert kind == "bal"
    out = []
    for order in ((a, b, c, d), (c, d, a, b)):
        sl = _slots(order)
        v = [0] * 16
        for x, n in PAT_COUNTS.items():
            y = (x[sl[0]] << 3) | (x[sl[1]] << 2) | (x[sl[2]] << 1) | x[sl[3]]
            v[y] += n
        out.append(v)
    return out

# ---- independent construction #2: string permutation on pattern keys ----
def _expvec_str(perm):
    v = [0] * 16
    pos = {ch: i for i, ch in enumerate(LEAVES)}
    for pstr, n in QC["pattern_counts"].items():
        s = "".join(pstr[pos[ch]] for ch in perm)
        v[int(s, 2)] += int(n)
    return v

# the two independent constructions are asserted equal on import whenever
# the production data is present (skipped when it is absent — see above)
if PAT_COUNTS is not None:
    for _ti in CAT_IDX:
        assert expvec_cat(_ti) == _expvec_str(TOPOLOGIES[_ti][1]), _ti
        assert sum(expvec_cat(_ti)) == 318
    for _ti in BAL_IDX:
        _kind, (_a, _b, _c, _d) = TOPOLOGIES[_ti]
        _v1, _v2 = expvec_bal_cases(_ti)
        assert _v1 == _expvec_str((_a, _b, _c, _d)) and _v2 == _expvec_str((_c, _d, _a, _b))
        assert sum(_v1) == 318 and sum(_v2) == 318
    assert expvec_cat(0) == _expvec_str("ABCD")  # base topology = identity

def used_values(vecs):
    """sorted list of value indices with nonzero exponent in any vector."""
    u = set()
    for v in vecs:
        u |= {m for m in range(16) if v[m]}
    return sorted(u)

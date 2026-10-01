#!/usr/bin/env python3
"""Fixed-seed SYNTHETIC generator + orbit-family builder (tests only).

No census data anywhere; this file builds synthetic Sp(4,Z) monodromy
generators (symplectic transvections) and synthetic orbit families with
PROVABLE ground truth, consumed only by test_dedup.py with the fixed seeds
in TEST_SEEDS.json.

Ground-truth guarantees:
- distinct-orbit: families are seeded at bases with pairwise DISTINCT
  budgets |f^T Sigma h|; the budget is an exact group invariant, so
  cross-family merges are impossible (proof-level negative control).
- same-orbit: every member is base * (random word), witnessed by
  construction (proof-level positive control).
Path norms are capped at generation so the engines' ball caps
provably contain the connecting paths (recoverability by construction).
"""
import random
import witness

def make_gens(seed, n=4, count=2):
    """Symplectic transvections x -> x + c*(x^T Sigma v)*v, exact-checked."""
    rng = random.Random(seed)
    S = witness.sig(n)
    gens = []
    while len(gens) < count:
        v = tuple(rng.randrange(-1, 2) for _ in range(n))
        if not any(v):
            continue
        c = rng.choice([1, -1])
        Sv = witness.mv(S, v)
        M = [[(1 if i == j else 0) + c * v[i] * Sv[j] for j in range(n)]
             for i in range(n)]
        if all(M[i][j] == (1 if i == j else 0) for i in range(n)
               for j in range(n)):
            continue
        if max(abs(x) for r in M for x in r) > 2 or M in gens:
            continue
        if not witness.sympl_ok(M, S):
            raise AssertionError("transvection not symplectic (bug)")
        gens.append(M)
    return gens

def alphabet(n_gens):
    return ["S", "s", "T", "t"] + ["%s%d" % (p, i) for i in range(n_gens)
                                   for p in ("M", "m")]

def _n2(st):
    return sum(x * x for x in st[0]) + sum(x * x for x in st[1])

def _l1(st):
    return sum(abs(x) for x in st[0]) + sum(abs(x) for x in st[1])

def make_family(rng, gens, tab, b, size, cfg):
    """One orbit family at budget b: base + (size-1) random-word images.

    Every accepted member's generating path stays inside path_n2_cap and
    path_l1_cap, so both engines' search balls contain the way back."""
    n = cfg["n"]
    base = ((b,) + (0,) * (n - 1), (0,) * (n // 2) + (1,) + (0,) * (n // 2 - 1))
    toks = alphabet(len(gens))
    mem, seen, tries = [base], {base}, 0
    while len(mem) < size and tries < 40000:
        tries += 1
        st, p2, p1 = base, _n2(base), _l1(base)
        for _ in range(rng.randrange(1, cfg["wordlen_max"] + 1)):
            st = witness.step(st, rng.choice(toks), tab)
            p2, p1 = max(p2, _n2(st)), max(p1, _l1(st))
        if st in seen or _n2(st) > cfg["member_n2_cap"]:
            continue
        if p2 > cfg["path_n2_cap"] or p1 > cfg["path_l1_cap"]:
            continue
        seen.add(st)
        mem.append(st)
    if len(mem) < size:
        raise RuntimeError("family generation exhausted at budget %d" % b)
    return mem

def build_battery(seed_fam, seed_shuf, gens, cfg):
    """All families + plants, shuffled. Returns (states, truth, tags, labels).

    truth = ground-truth partition (list of index lists); plants are their
    own single-member families at unique budgets (own orbit class by the
    budget invariant), tagged 'plant' and inserted PRE-dedup."""
    rng = random.Random(seed_fam)
    tab = witness.make_table(gens, cfg["n"])
    rows = []
    for b in cfg["honest_budgets"]:
        size = rng.randrange(cfg["family_size_range"][0],
                             cfg["family_size_range"][1] + 1)
        lab = "on" if b <= cfg["on_locus_label_budget_max"] else "off"
        for st in make_family(rng, gens, tab, b, size, cfg):
            rows.append((st, b, "honest", lab))
    for b in cfg["plant_budgets"]:
        n = cfg["n"]
        base = ((b,) + (0,) * (n - 1),
                (0,) * (n // 2) + (1,) + (0,) * (n // 2 - 1))
        rows.append((base, b, "plant", "off"))
    random.Random(seed_shuf).shuffle(rows)
    states = [r[0] for r in rows]
    tags = [r[2] for r in rows]
    labels = [r[3] for r in rows]
    fams = {}
    for i, r in enumerate(rows):
        fams.setdefault(r[1], []).append(i)
    truth = sorted(sorted(c) for c in fams.values())
    return states, truth, tags, labels

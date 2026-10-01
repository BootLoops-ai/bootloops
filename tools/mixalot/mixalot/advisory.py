"""Exact-frontier advisories for users sizing a run.

MIXALOT's own exact evaluator needs no estimator; these advisories exist for
users who want to know how far the exact routes reach before they consider a
sampler. Every line is quoted from the relevant engine docstring (measured
practical frontiers), not a theorem.
"""

# exact-evaluator practical frontier, quoted from the engine docstring
# (bigg.py, MEMORY/PRACTICAL FRONTIER section; 31 GB ulimit basis):
EXACT_FRONTIER = {
    2: "g=2: 1-D prefix DP — N in the 1e5 class is easy (time-bound)",
    3: "g=3: N=3000 ~1 GB, N=6000 ~5 GB, N=1e4 ~25-30 GB (edge; do not "
       "attempt on a 31 GB limit)",
    4: "g>=4 generic path: N ~ 500-700 memory edge; time binds first",
}


def advise(N=None, g=None):
    """Return the advisory block: the exact-frontier line for the requested
    g (all of them when g is None), plus a sizing note for the requested N.
    Purely a quotation layer — nothing computed."""
    out = {"exact_frontier": (EXACT_FRONTIER.get(g)
                              if g else list(EXACT_FRONTIER.values()))}
    if N is not None:
        out["note_for_N"] = (
            f"requested N={N}: if the exact evaluator's frontier (above) "
            "covers it, use the exact path and ignore estimators entirely; "
            "otherwise validate any sampler against the exact routes on the "
            "largest N they still reach (the battery's exact-vs-MC leg is "
            "the pattern).")
    return out

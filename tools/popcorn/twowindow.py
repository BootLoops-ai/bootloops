"""popcorn.twowindow — two-window coalescent simulation harness for the
two-site frequency spectrum (msprime).

Aliases, in place and by identity (see popcorn._loader):
  twowindow_sim — run_two_window(d_bp, n_reps, seed, ...): two 1-bp
                  branch-AFS windows a distance d apart on one sequence
                  (crossover, gene conversion with geometric tracts, any
                  msprime model: sweeps via sweep_model, Beta / Dirac
                  coalescents, demography), accumulating the ordered
                  M[i-1][j-1] += afs_i(A) afs_j(B) over replicates;
                  run_two_locus(rho_scaled, n_reps, seed, ...): two sites at
                  population-scaled recombination rho = 4 Ne m;
                  run_probe_track(...): many same-window probes along one
                  long recombining sequence; pool_pairs / pair_keys (the
                  diagonal-doubled pair vector, 28 entries at n = 8),
                  block_jackknife (delete-one-block SEs), layout /
                  margin_for / window_edges / sim_kwargs_for (geometry),
                  seed_from_name / root_seed (seed law), n_afs_distinct
                  (measured decorrelation fraction and vacuity check).
                  numpy; msprime + tskit imported on first simulation.

Register: SIMULATION (Monte Carlo estimates with jackknife standard errors;
nothing exact, nothing certified). The exact comparand for a completely
linked arm (rho = 0, or no process between the windows) is popcorn.twosfs
(exact_linked_pooled below; the n = 8 Kingman values also ship under
reference/twowindow/); popcorn.twolocus is the float-validated expectation
at rho > 0. Label every number from this module as simulated when it
appears beside exact or certified ones.
"""
import json

from ._loader import import_in_place
from ._pins import PKG_DIR

twowindow_sim = import_in_place("twowindow_sim")

# Convenience re-exports.
run_two_window = twowindow_sim.run_two_window
run_two_locus = twowindow_sim.run_two_locus
run_probe_track = twowindow_sim.run_probe_track
sweep_model = twowindow_sim.sweep_model
layout = twowindow_sim.layout
margin_for = twowindow_sim.margin_for
window_edges = twowindow_sim.window_edges
sim_kwargs_for = twowindow_sim.sim_kwargs_for
pair_keys = twowindow_sim.pair_keys
pair_counts = twowindow_sim.pair_counts
pool_pairs = twowindow_sim.pool_pairs
jackknife_cov = twowindow_sim.jackknife_cov
block_jackknife = twowindow_sim.block_jackknife
probe_positions = twowindow_sim.probe_positions
probe_weights = twowindow_sim.probe_weights
seed_from_name = twowindow_sim.seed_from_name
root_seed = twowindow_sim.root_seed
keys28 = twowindow_sim.keys28
v28_from_M = twowindow_sim.v28_from_M
block_jackknife_v28 = twowindow_sim.block_jackknife_v28
N_HAP = twowindow_sim.N_HAP
NE = twowindow_sim.NE
MARGIN_BP = twowindow_sim.MARGIN_BP

REFERENCE_DIR = PKG_DIR / "reference" / "twowindow"


def load_reference(n_hap=8):
    """The shipped exact Kingman linked-pair moments at n_hap (only 8 is
    shipped): dict with keys [(i,j)], E_Li, E_LiLj and pooled as Fractions,
    plus the raw JSON under 'doc'."""
    from fractions import Fraction
    p = REFERENCE_DIR / f"linked_kingman_n{int(n_hap)}.json"
    if not p.is_file():
        raise FileNotFoundError(f"no shipped reference for n_hap = {n_hap} ({p.name})")
    with open(p) as fh:
        doc = json.load(fh)
    keys = [tuple(k) for k in doc["keys"]]
    return {"keys": keys,
            "E_Li": [Fraction(x) for x in doc["E_Li"]],
            "E_LiLj": {k: Fraction(doc["E_LiLj"][f"{k[0]},{k[1]}"]) for k in keys},
            "pooled": [Fraction(doc["pooled"][f"{k[0]},{k[1]}"]) for k in keys],
            "doc": doc}


def _exact_engine():
    """popcorn.twosfs if importable, else the flat twosfs_engine module;
    ImportError naming the missing piece otherwise."""
    try:
        from . import twosfs as T                       # package alias
        T.lambda_moments, T.kingman_rate                # noqa: B018 (touch)
        return T
    except (ImportError, AttributeError):
        pass
    try:
        return import_in_place("twosfs_engine")         # flat engine beside this file
    except ImportError as e:
        raise ImportError("the exact linked-pair comparand needs popcorn.twosfs "
                          "(twosfs_engine.py + lambda_exact.py)") from e


def exact_linked_pooled(n_hap=N_HAP, lam=None):
    """Exact comparand for a completely linked arm: (keys, pooled) with
    pooled[k] = E[L_i L_j] / sum_{i<=j} E[L_i L_j] as Fractions in pair_keys
    order — the expectation of pool_pairs(M) when both windows carry one
    genealogy — under the Kingman coalescent (lam None) or any
    Lambda-coalescent rate function lam(b, k) from popcorn.twosfs /
    popcorn.lambda_coalescent. Computed live by popcorn.twosfs."""
    T = _exact_engine()
    n = int(n_hap)
    lam = T.kingman_rate() if lam is None else lam
    _u, v = T.lambda_moments(n, lam)
    keys = pair_keys(n)
    tot = sum(v[k] for k in keys)
    return keys, [v[k] / tot for k in keys]


__all__ = ["twowindow_sim", "run_two_window", "run_two_locus",
           "run_probe_track", "sweep_model", "layout", "margin_for",
           "window_edges", "sim_kwargs_for", "pair_keys", "pair_counts",
           "pool_pairs", "jackknife_cov", "block_jackknife",
           "probe_positions", "probe_weights", "seed_from_name", "root_seed",
           "keys28", "v28_from_M", "block_jackknife_v28", "N_HAP", "NE",
           "MARGIN_BP", "REFERENCE_DIR", "load_reference",
           "exact_linked_pooled"]

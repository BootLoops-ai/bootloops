#!/usr/bin/env python3
"""twowindow_sim — two-window coalescent simulation harness for the two-site
frequency spectrum (msprime, branch mode).

Register: SIMULATION. Every spectrum returned by this module is a Monte
Carlo estimate with Monte Carlo error (a delete-one-block jackknife standard
error is supplied on request). Nothing here is exact or certified; the exact
comparand at zero recombination is the linked-pair engine (twosfs_engine /
popcorn.twosfs), and popcorn.twolocus gives the float-validated expectation
at any recombination rate. numpy at import; msprime (with tskit) is imported
only when a simulation is actually requested, so the pooling, jackknife and
geometry helpers work without it.

Object
------
For a sample of n haploid sequences let afs_i(W), i = 1..n-1, be the
branch-mode allele frequency spectrum of the local genealogy in a 1-bp
window W: the total branch length (in generations) subtending exactly i of
the n samples (tskit allele_frequency_spectrum with mode="branch",
polarised=True, span_normalise=False). For two windows A and B on the same
simulated sequence the harness accumulates, over independent replicates,
the ORDERED matrix

    M[i-1][j-1] += afs_i(A) * afs_j(B),        i, j = 1..n-1,

so M / replicates estimates E[T_i^A T_j^B], the second moments whose
normalized form is the expected two-site frequency spectrum of one
infinite-sites mutation in each window (small-theta limit). With no process
acting between the windows (no crossover, no gene conversion) both windows
carry the SAME genealogy, afs(A) == afs(B) bit for bit, and M estimates the
completely linked object E[L_i L_j] that twosfs_engine computes exactly;
crossover or gene conversion between the windows interpolates M toward the
product of independent genealogies, outer(E[L_i], E[L_j]). Each replicate
whose two windows return bitwise-different spectra is counted in
n_afs_distinct: a directly measured decorrelation fraction, and the
harness's own vacuity check — an arm with no between-window process MUST
return n_afs_distinct == 0, otherwise the geometry or the simulator call is
wrong and every other number from the run is void.

Two geometries
--------------
run_two_window(d_bp, ...): one contiguous sequence, window A at position
    pa = margin, window B at pb = pa + d_bp, sequence_length = pb + 1 + margin
    (layout()), with `margin` bp of flanking sequence on each side so that
    gene-conversion tracts covering a window are not edge-truncated; margin
    law margin_for(tract) = max(MARGIN_BP, 10 * mean tract length) unless an
    explicit margin is passed. Between-window processes: crossover at a
    per-bp rate (recomb_rate, default 0), gene conversion at a per-bp
    initiation rate with geometric tract lengths of a given mean (msprime's
    gene_conversion_rate / gene_conversion_tract_length), and any msprime
    ancestry model or model list (sweeps via sweep_model(), Beta- and
    Dirac-coalescents via msprime.BetaCoalescent / DiracCoalescent,
    piecewise demography via extra_sim_kwargs={"demography": ...}).
run_two_locus(rho_scaled, ...): a pure two-locus system, sites 0 and 1 on a
    sequence of length 2 with crossover probability m = rho_scaled / (4 Ne)
    per generation between them, placed as msprime.RateMap(position=[0, 1,
    2], rate=[0, m]). In msprime's discrete-genome coordinates a breakpoint
    at integer position x separates [0, x) from [x, L) and is driven by the
    rate-map mass on [x, x+1): mass on [1, 2) controls the link between
    sites 0 and 1, mass on [0, 1) is inert (a breakpoint at 0 separates
    nothing). rho_scaled is therefore the population-scaled recombination
    rate 4 Ne m between the two sites, the rho of popcorn.twolocus.
    rho_scaled == 0 is the completely linked arm (vacuity check applies).
run_probe_track(...): one long recombining sequence with many 1-bp probe
    windows; per probe p the SAME-window product M_p[i][j] += afs_i(p) *
    afs_j(p) (the completely linked object of the local tree at p), with
    trapezoid bp-weights for pooling over probes. The single-window object
    is blind to anything acting BETWEEN sites (the marginal genealogy at one
    site is the same coalescent whatever the recombination elsewhere); it is
    the right geometry for processes that distort the local tree itself,
    e.g. a selective sweep at a given distance from the probe.

Pooling convention and errors
-----------------------------
pair_keys(n) lists the index pairs (i, j), 1 <= i <= j <= n-1, in the order
(1,1),(1,2),...,(1,n-1),(2,2),...,(n-1,n-1) — n(n-1)/2 of them, 28 at the
default n = 8 — the same order as twosfs_engine.flatten and
twolocus_engine.pool_pairs. pool_pairs(M) reads the symmetrized matrix on
those pairs, 2*M_ii on the diagonal and M_ij + M_ji off it, normalized to
unit sum; for the completely linked arm its expectation is E[L_i L_j] /
sum_{i<=j} E[L_i L_j] (the exact values normalized the same way). It is a
pooling convention, not the class law of an unordered pair of sites (which
carries M_ii, not 2 M_ii, on the diagonal): compare only against vectors
pooled the same way. block_jackknife(M_blocks) splits the replicates into
equal consecutive blocks (n_blocks in the run_* calls), pools each block's
un-normalized pair counts, and returns the pooled vector of the total with
the delete-one-block jackknife standard error of the normalized composition,
    Sigma = (B-1)/B * sum_k (q_(-k) - qbar)(q_(-k) - qbar)^T,
q_(-k) the composition with block k deleted; se = sqrt(diag Sigma). The SE is
itself an estimate whose relative noise is about 1/sqrt(2(B-1)) (B = 10:
~24 %, B = 40: ~11 %): quote the block count with every z-score.

Seeds and batches: `seed` is an int root seed or a str hashed to one by
seed_from_name (blake2b, 4 bytes). All msprime seeds are drawn from
numpy.random.default_rng(root) in code order, one per batch of `batch`
replicates (msprime num_replicates), so a run is reproducible from (seed,
n_reps, batch, n_blocks) at a fixed msprime version. A batch in which
msprime raises (e.g. a sweep-trajectory failure) is counted in `failed` and
skipped, never silently dropped (check failed == 0; msprime's last message
is returned as last_error); a run in which EVERY replicate failed has no
estimate and raises RuntimeError with that message.

Units and conventions: msprime's — population_size = Ne is the diploid
effective size when ploidy = 2, times and branch lengths in generations
(pooled vectors are scale free); samples = number of `ploidy`-ploid
individuals, n = samples * ploidy haploid sequences (default 4 x 2 = 8);
crossover and gene-conversion initiation rates per bp per generation.
sweep_model(position, s, Ne) builds msprime.SweepGenicSelection with
msprime's fitness parametrization W_bb = 1, W_Bb = 1 + s/2, W_BB = 1 + s;
in this package's convention (S = 4 Ne s' with genic fitnesses
1 : 1+s' : 1+2s') that is S = 2 Ne s.

Measured cost (one core, msprime with its C library): two-locus arm at n = 8
about 0.2 ms per replicate at rho = 0 and 0.35 ms at rho = 3; two-window arm
with d = 100 bp, margin 300 bp about 0.15 ms per replicate with no process
and 0.7 ms with gene conversion covering the windows; peak memory under
100 MB. Relative standard errors of the pooled 28-vector at n = 8: about
2 % for the largest cells, 5 % median and 11 % for the smallest at 4,000
replicates; half that at 16,000 (1/sqrt(replicates)).

Scope: neutral or swept single population (anything msprime.sim_ancestry
accepts through `model` / extra_sim_kwargs); expected second moments of the
two-window branch spectra and the decorrelation fraction; no mutation
simulation (branch mode is the infinite-sites expectation), no scoring or
fitting of the resulting vectors. Multiple-merger arms meant for comparison
with the exact Lambda-coalescent engine must run at ploidy = 1 (samples =
n): msprime's BetaCoalescent / DiracCoalescent at ploidy p > 1 merge up to
2p groups simultaneously (a Xi-coalescent), which is a different process
(measured at n = 8, Beta(3/2), 16,000 replicates: ploidy 1 agrees with the
exact Beta spectrum at max|z| 2.9 and rejects Kingman at 41 SE; ploidy 2
sits 38 SE from the exact Lambda-Beta spectrum).

CLI:  python3 twowindow_sim.py [rho] [n_reps]
  runs the two-locus arm at rho = 0 (vacuity check printed) and at the given
  rho (default 1.0) with n_reps replicates (default 2000) in 10 blocks, and
  prints the pooled vectors with jackknife SEs. Writes nothing.

References
  Baumdicker, F. et al. (2022). Efficient ancestry and mutation simulation
    with msprime 1.0. Genetics 220, iyab229.  (The simulator.)
  Ralph, P., Thornton, K. & Kelleher, J. (2020). Efficiently summarizing
    relationships in large samples: a general duality between statistics of
    genealogies and genomes. Genetics 215, 779-797.  (Branch-mode statistics,
    the allele frequency spectrum of a genealogy.)
  Hudson, R. R. (1983). Properties of a neutral allele model with intragenic
    recombination. Theor. Popul. Biol. 23, 183-201; Griffiths, R. C. (1981).
    Neutral two-locus multiple allele models with recombination. Theor.
    Popul. Biol. 19, 169-186.  (The two-locus ancestral process.)
  Wiuf, C. & Hein, J. (2000). The coalescent with gene conversion. Genetics
    155, 451-462.  (Gene conversion in the coalescent.)
  Braverman, J. M., Hudson, R. R., Kaplan, N. L., Langley, C. H. & Stephan,
    W. (1995). The hitchhiking effect on the site frequency spectrum of DNA
    polymorphisms. Genetics 140, 783-796; Kern, A. D. & Schrider, D. R.
    (2016). Discoal: flexible coalescent simulations with selection.
    Bioinformatics 32, 3839-3841.  (The structured-coalescent sweep model
    msprime.SweepGenicSelection implements.)
  Fu, Y.-X. (1995). Statistical properties of segregating sites. Theor.
    Popul. Biol. 48, 172-197.  (Exact Kingman E[L_i L_j], the rho = 0
    comparand.)
"""
import math
import time
from hashlib import blake2b

import numpy as np

# ---------------------------------------------------------------- defaults
SAMPLES = 4                    # individuals sampled
PLOIDY = 2
N_HAP = SAMPLES * PLOIDY       # haploid sample size n (classes i = 1..n-1)
NE = 1.0e4                     # msprime population_size (diploid Ne at ploidy 2)
MARGIN_BP = 3_000              # flanking sequence beyond each window, bp
GC_TRACT_BP = 300.0            # default mean gene-conversion tract length, bp
BATCH = 5_000                  # msprime num_replicates per seed batch (two-window, two-locus)
BATCH_TRACK = 500              # per seed batch for the long-sequence probe track
R_BP = 1.0e-8                  # per-bp per-generation crossover rate (probe-track default)
TRACK_LENGTH = 2_000_000       # probe-track sequence length, bp
SWEEP_DT = 1e-6                # sweep trajectory time step (msprime dt)

# probe-track offsets from the sequence midpoint, bp (23 probes)
PROBE_OFFSETS = tuple(sorted(
    {0, 1_000, 3_000, 10_000, 30_000, 100_000, 200_000, 300_000,
     450_000, 600_000, 750_000, 900_000}
    | {-o for o in (1_000, 3_000, 10_000, 30_000, 100_000, 200_000,
                    300_000, 450_000, 600_000, 750_000, 900_000)}))


def _require_msprime():
    try:
        import msprime
    except ImportError as e:                       # named refusal, not a crash at import
        raise ImportError("twowindow_sim needs msprime (and tskit) to simulate: "
                          "pip install msprime") from e
    return msprime


# ------------------------------------------------------------------- seeds
def seed_from_name(name):
    """str -> 32-bit root seed (blake2b, 4-byte digest, big-endian)."""
    h = blake2b(name.encode(), digest_size=4)
    return int.from_bytes(h.digest(), "big")


def root_seed(seed):
    """int -> itself; str -> seed_from_name(str)."""
    if isinstance(seed, str):
        return int(seed_from_name(seed))
    if isinstance(seed, (int, np.integer)) and not isinstance(seed, bool):
        return int(seed)
    raise TypeError("seed must be an int root seed or a str seed name")


# ------------------------------------------------------ pooling and errors
def pair_keys(n_hap=N_HAP):
    """[(i, j) for 1 <= i <= j <= n_hap-1] in the order (1,1),(1,2),...,
    (n-1,n-1): n(n-1)/2 pairs (28 at n = 8)."""
    n = int(n_hap)
    return [(i, j) for i in range(1, n) for j in range(i, n)]


def pair_counts(M):
    """Un-normalized pooled pair vector of an ordered (n-1)x(n-1) matrix:
    2*M_ii on the diagonal, M_ij + M_ji off it, in pair_keys order."""
    Q = np.asarray(M, dtype=float)
    if Q.ndim == 3 and Q.shape[0] == 1:
        Q = Q[0]
    if Q.ndim != 2 or Q.shape[0] != Q.shape[1] or Q.shape[0] < 1:
        raise ValueError(f"expected a square (n-1)x(n-1) matrix, got {Q.shape}")
    n = Q.shape[0] + 1
    return np.array([2.0 * Q[i - 1, j - 1] if i == j
                     else Q[i - 1, j - 1] + Q[j - 1, i - 1]
                     for (i, j) in pair_keys(n)])


def pool_pairs(M):
    """Ordered (n-1)x(n-1) sum matrix -> normalized pair vector (pair_counts
    divided by its sum). A pooling convention (diagonal doubled): compare
    only against vectors pooled the same way."""
    v = pair_counts(M)
    return v / v.sum()


def jackknife_cov(block_vectors):
    """Delete-one-block jackknife of a pooled composition.

    block_vectors: sequence (or dict) of un-normalized count vectors, one
    per block. Returns (qhat, Sigma, B): qhat = total / total.sum(),
    Sigma = (B-1)/B * sum_k (q_(-k) - qbar)(q_(-k) - qbar)^T with q_(-k) the
    normalized composition with block k deleted, B the number of blocks."""
    if isinstance(block_vectors, dict):
        block_vectors = [block_vectors[k] for k in sorted(block_vectors)]
    Vs = np.array([np.asarray(v, dtype=float) for v in block_vectors])
    if Vs.ndim != 2 or Vs.shape[0] < 2:
        raise ValueError("need >= 2 blocks of equal length for a jackknife")
    Vtot = Vs.sum(axis=0)
    qhat = Vtot / Vtot.sum()
    reps = []
    for k in range(Vs.shape[0]):
        Vd = Vtot - Vs[k]
        reps.append(Vd / Vd.sum())
    reps = np.array(reps)
    Gb = len(reps)
    qbar = reps.mean(axis=0)
    Dv = reps - qbar
    Sigma = (Gb - 1.0) / Gb * (Dv.T @ Dv)
    return qhat, Sigma, Gb


def block_jackknife(M_blocks):
    """Delete-one-block jackknife of the pooled pair vector: returns
    (pooled, se) with pooled == pool_pairs(sum of blocks) and se =
    sqrt(diag Sigma) from jackknife_cov over the blocks' pair_counts.
    Needs >= 2 blocks. The SE is an estimate; its own relative noise is
    about 1/sqrt(2(B-1)) — quote the block count with every z."""
    B = [np.asarray(b, dtype=float) for b in M_blocks]
    nb = len(B)
    if nb < 2:
        raise ValueError("need >= 2 blocks for a jackknife")
    shp = B[0].shape
    for b in B:
        if b.ndim != 2 or b.shape != shp or shp[0] != shp[1]:
            raise ValueError(f"every block must be the same square matrix shape, got {b.shape}")
    qhat, Sigma, Gb = jackknife_cov([pair_counts(b) for b in B])
    v_ref = pool_pairs(sum(B))
    if not np.allclose(qhat, v_ref, rtol=0.0, atol=1e-12):
        raise RuntimeError("jackknife total and pool_pairs(sum of blocks) disagree")
    se = np.sqrt(np.clip(np.diag(Sigma), 0.0, None))
    return v_ref, se


# ----------------------------------------------------------------- geometry
def margin_for(gc_tract=None, margin_bp=MARGIN_BP):
    """Margin law: max(margin_bp, 10 * mean tract length) bp of sequence
    beyond each window (residual edge truncation of a geometric tract is
    then below e^-10). gc_tract None or 0 -> margin_bp."""
    if not gc_tract:
        return int(margin_bp)
    return int(max(int(margin_bp), math.ceil(10.0 * float(gc_tract))))


def layout(d_bp, margin=None):
    """(pa, pb, L): 1-bp windows at pa = margin and pb = pa + d_bp on a
    sequence of length L = pb + 1 + margin. margin None -> MARGIN_BP."""
    d_bp = int(d_bp)
    if d_bp < 1:
        raise ValueError("d_bp must be >= 1")
    margin = MARGIN_BP if margin is None else int(margin)
    if margin < 1:
        raise ValueError("margin must be >= 1 bp")
    pa = margin
    pb = pa + d_bp
    L = pb + 1 + margin
    return pa, pb, L


def window_edges(pa, pb, L):
    """tskit window edges for 1-bp windows at pa and pb on [0, L), and the
    two row indices of those windows in the windowed AFS: (edges, ra, rb)."""
    w = np.array([0.0, pa, pa + 1, pb, pb + 1, float(L)])
    return w, 1, 3


def sim_kwargs_for(L, Ne=NE, samples=SAMPLES, ploidy=PLOIDY, gc_rate=0.0,
                   gc_tract=None, recomb_rate=0.0, model=None,
                   extra_sim_kwargs=None):
    """msprime.sim_ancestry keyword arguments for a two-window arm: crossover
    rate as given (default 0), gene conversion only when gc_rate > 0 (then
    gc_tract, the mean tract length in bp, is required), `model` for sweeps /
    Beta / Dirac / model lists, extra_sim_kwargs for e.g. demography=... .
    `samples` is an integer count of `ploidy`-ploid individuals (the runners
    read n = samples * ploidy from it; msprime sample-set forms are not
    supported here). num_replicates and random_seed are owned by the
    batch/seed law and are refused here."""
    kw = dict(samples=int(samples), ploidy=int(ploidy),
              population_size=float(Ne),
              recombination_rate=recomb_rate, sequence_length=L)
    if gc_rate and float(gc_rate) > 0.0:
        if not gc_tract or float(gc_tract) <= 0.0:
            raise ValueError("gc_rate > 0 requires gc_tract (mean bp) > 0")
        kw["gene_conversion_rate"] = float(gc_rate)
        kw["gene_conversion_tract_length"] = float(gc_tract)
    if model is not None:
        kw["model"] = model
    if extra_sim_kwargs:
        for k in extra_sim_kwargs:
            if k in ("num_replicates", "random_seed"):
                raise ValueError(f"{k} is owned by the batch/seed law")
        kw.update(extra_sim_kwargs)
    return kw


def sweep_model(position, s, Ne=NE, ploidy=PLOIDY, dt=SWEEP_DT,
                start_frequency=None, end_frequency=None):
    """msprime.SweepGenicSelection at `position` with selection coefficient s
    in msprime's parametrization (W_bb = 1, W_Bb = 1 + s/2, W_BB = 1 + s;
    S = 2 Ne s in this package's S = 4 Ne s' convention). Defaults: a hard
    sweep from one copy, start_frequency = 1/(ploidy Ne), to fixation,
    end_frequency = 1 - 1/(ploidy Ne); trajectory step dt. Compose as
    model=[sweep_model(...), msprime.StandardCoalescent()] (sweep ending at
    sampling time) or [msprime.StandardCoalescent(duration=tau), sweep,
    msprime.StandardCoalescent()] (ended tau generations ago). See msprime's
    documentation of SweepGenicSelection for the model's own limits
    (single population, no size change during the sweep, start_frequency)."""
    msprime = _require_msprime()
    lo = 1.0 / (ploidy * Ne)
    return msprime.SweepGenicSelection(
        position=float(position),
        start_frequency=float(lo if start_frequency is None else start_frequency),
        end_frequency=float((1.0 - lo) if end_frequency is None else end_frequency),
        s=float(s), dt=dt)


# ------------------------------------------------------- two-window runner
def _n_hap(sim_kwargs):
    return int(sim_kwargs["samples"]) * int(sim_kwargs.get("ploidy", PLOIDY))


def _run_2w(sim_kwargs, n_reps, rng, windows, ra, rb, M, batch):
    """Two-window batched runner. M has shape (1, n-1, n-1); accumulates the
    ordered outer product afs(A) x afs(B). Returns (done, failed,
    n_afs_distinct, last_error)."""
    msprime = _require_msprime()
    n_hap = _n_hap(sim_kwargs)
    done = failed = distinct = 0
    last_error = None
    while done + failed < n_reps:
        nb = min(batch, n_reps - done - failed)
        seed = int(rng.integers(1, 2**31 - 1))
        try:
            reps = msprime.sim_ancestry(
                num_replicates=nb, random_seed=seed, **sim_kwargs)
            for ts in reps:
                afs = ts.allele_frequency_spectrum(
                    mode="branch", windows=windows, polarised=True,
                    span_normalise=False)
                A = afs[ra, 1:n_hap]
                B = afs[rb, 1:n_hap]
                M[0] += np.outer(A, B)
                if not np.array_equal(A, B):
                    distinct += 1
            done += nb
        except Exception as e:                  # counted, never dropped
            failed += nb
            last_error = f"{type(e).__name__}: {e}"
    return done, failed, distinct, last_error


def _run_blocks(sim_kwargs, windows, ra, rb, n_reps, n_blocks, rng, batch):
    if n_blocks < 1 or n_reps < n_blocks or n_reps % n_blocks:
        raise ValueError("n_reps must be a positive multiple of n_blocks")
    n_hap = _n_hap(sim_kwargs)
    if n_hap < 2:
        raise ValueError("need at least 2 haploid samples (samples * ploidy)")
    per = n_reps // n_blocks
    blocks, done, failed, distinct = [], 0, 0, 0
    dist_blocks = []
    last_error = None
    for _ in range(n_blocks):
        M = np.zeros((1, n_hap - 1, n_hap - 1))
        d, f, x, err = _run_2w(sim_kwargs, per, rng, windows, ra, rb, M,
                               int(batch))
        blocks.append(M[0])
        dist_blocks.append(int(x))
        done += d
        failed += f
        distinct += x
        last_error = err or last_error
    _all_failed_guard(done, failed, last_error)
    return blocks, done, failed, distinct, dist_blocks, last_error


def _all_failed_guard(done, failed, last_error):
    """Failed batches are counted, not raised (a sweep trajectory can fail
    now and then); but a run in which EVERY replicate failed carries no
    estimate at all and almost always means a bad simulator argument, so
    that case raises with msprime's last message."""
    if done == 0 and failed > 0:
        raise RuntimeError(f"all {failed} replicates failed inside msprime "
                           f"(last error: {last_error})")


def _finish(out, blocks, done, failed, distinct, dist_blocks, last_error,
            pool, t0):
    M = sum(blocks)
    n_hap = M.shape[0] + 1
    out.update({
        "M": M,
        "M_blocks": blocks if len(blocks) > 1 else None,
        "n_afs_distinct_blocks": dist_blocks if len(blocks) > 1 else None,
        "done": int(done), "failed": int(failed), "last_error": last_error,
        "n_afs_distinct": int(distinct),
        "afs_distinct_fraction": distinct / max(done, 1),
        "keys": pair_keys(n_hap),
        "pooled": None, "pooled_se": None,
        "register": "simulation (msprime branch-mode Monte Carlo; Monte "
                    "Carlo error applies)",
        "object": "M[i-1][j-1] = sum over replicates of afs_i(window A) * "
                  "afs_j(window B), branch-mode polarized 1-bp windows, "
                  "classes i, j = 1..n-1, ORDERED (pool_pairs symmetrizes)",
    })
    if pool and done > 0:
        if len(blocks) > 1:
            v, se = block_jackknife(blocks)
            out["pooled"], out["pooled_se"] = v, se
        else:
            out["pooled"] = pool_pairs(M)
    out["elapsed_s"] = time.time() - t0
    return out


def run_two_window(d_bp, n_reps, seed, Ne=NE, samples=SAMPLES, gc_rate=0.0,
                   gc_tract=None, recomb_rate=0.0, model=None, batch=BATCH,
                   n_blocks=1, margin=None, pool=True, ploidy=PLOIDY,
                   extra_sim_kwargs=None):
    """Contiguous-sequence two-window arm. Two 1-bp branch-AFS windows at
    separation d_bp on a sequence with `margin` bp beyond each window
    (default: margin_for(gc_tract) when gc_rate > 0, else MARGIN_BP);
    accumulates the ORDERED sum M[i-1][j-1] += afs_i(A) * afs_j(B) over
    n_reps msprime replicates (batch/seed law in the module docstring;
    batches that raise are counted in `failed`, never dropped).

    seed: int root seed, or str -> seed_from_name(str).
    Returns dict(M, M_blocks, pooled, pooled_se (jackknife over n_blocks if
    > 1), keys, done, failed, n_afs_distinct, afs_distinct_fraction,
    params, ...). An arm with no between-window process (gc_rate 0,
    recomb_rate 0) MUST return n_afs_distinct == 0 (vacuity check)."""
    t0 = time.time()
    if margin is None:
        margin = margin_for(gc_tract if (gc_rate and gc_rate > 0) else None)
    pa, pb, L = layout(d_bp, margin)
    windows, ra, rb = window_edges(pa, pb, L)
    sim_kwargs = sim_kwargs_for(L, Ne=Ne, samples=samples, ploidy=ploidy,
                                gc_rate=gc_rate, gc_tract=gc_tract,
                                recomb_rate=recomb_rate, model=model,
                                extra_sim_kwargs=extra_sim_kwargs)
    root = root_seed(seed)
    rng = np.random.default_rng(root)
    blocks, done, failed, distinct, dist_blocks, last_error = _run_blocks(
        sim_kwargs, windows, ra, rb, int(n_reps), int(n_blocks), rng, batch)
    out = {
        "geometry": "contiguous two-window",
        "seed_root": root,
        "seed_law": "all msprime seeds drawn from default_rng(seed_root) in "
                    "code order, one per batch",
        "params": {
            "n_haploid": int(samples) * int(ploidy), "samples": int(samples),
            "ploidy": int(ploidy), "Ne": float(Ne),
            "gene_conversion_rate_per_bp": float(gc_rate or 0.0),
            "gene_conversion_tract_mean_bp":
                float(gc_tract) if (gc_rate and gc_rate > 0) else None,
            "crossover_rate": (float(recomb_rate)
                               if isinstance(recomb_rate, (int, float))
                               else repr(recomb_rate)),
            "model": None if model is None else repr(model),
            "separation_bp": int(d_bp), "window_a_bp": pa, "window_b_bp": pb,
            "sequence_length": L, "margin_bp": int(margin),
            "margin_law": f"max({MARGIN_BP}, 10 * mean tract) bp unless an "
                          "explicit margin is passed",
            "probe_weights_bp": [1.0],
            "batch": int(batch), "n_blocks": int(n_blocks),
            "n_replicates_requested": int(n_reps),
            "extra_sim_kwargs": sorted(extra_sim_kwargs or {}),
        },
    }
    return _finish(out, blocks, done, failed, distinct, dist_blocks,
                   last_error, pool, t0)


def run_two_locus(rho_scaled, n_reps, seed, Ne=NE, samples=SAMPLES,
                  model=None, batch=BATCH, n_blocks=1, pool=True,
                  ploidy=PLOIDY, extra_sim_kwargs=None):
    """Pure two-locus arm: sites 0 and 1 on a sequence of length 2 with
    crossover probability m = rho_scaled / (4 Ne) per generation between
    them, as msprime.RateMap(position=[0, 1, 2], rate=[0, m]) (the mass on
    [1, 2) drives the breakpoint at position 1, which separates site 0 from
    site 1; mass on [0, 1) would be inert). rho_scaled = 4 Ne m is the
    coalescent-scaled rate at ploidy 2 (pair coalescence rate 1/(2 Ne) per
    generation); m is always placed as rho_scaled / (4 Ne), so at ploidy p
    the coalescent-scaled rate is rho_scaled * p / 2. rho_scaled == 0 is
    the completely linked arm: n_afs_distinct must be 0. Same return dict
    as run_two_window."""
    msprime = _require_msprime()
    t0 = time.time()
    rho = float(rho_scaled)
    if rho < 0:
        raise ValueError("rho_scaled must be >= 0")
    windows = np.array([0.0, 1.0, 2.0])
    kw = dict(samples=int(samples), ploidy=int(ploidy),
              population_size=float(Ne))
    if rho > 0:
        m = rho / (4.0 * float(Ne))
        kw["recombination_rate"] = msprime.RateMap(
            position=[0.0, 1.0, 2.0], rate=[0.0, m])
    else:
        m = 0.0
        kw["sequence_length"] = 2
    if model is not None:
        kw["model"] = model
    if extra_sim_kwargs:
        for k in extra_sim_kwargs:
            if k in ("num_replicates", "random_seed"):
                raise ValueError(f"{k} is owned by the batch/seed law")
        kw.update(extra_sim_kwargs)
    root = root_seed(seed)
    rng = np.random.default_rng(root)
    blocks, done, failed, distinct, dist_blocks, last_error = _run_blocks(
        kw, windows, 0, 1, int(n_reps), int(n_blocks), rng, batch)
    out = {
        "geometry": "pure two-locus",
        "seed_root": root,
        "seed_law": "all msprime seeds drawn from default_rng(seed_root) in "
                    "code order, one per batch",
        "params": {
            "n_haploid": int(samples) * int(ploidy), "samples": int(samples),
            "ploidy": int(ploidy), "Ne": float(Ne),
            "rho_scaled": rho, "inter_site_mass_per_gen": m,
            "rate_map": "RateMap(position=[0,1,2], rate=[0, m]): mass on "
                        "[1,2) controls the 0-1 link" if rho > 0 else None,
            "sequence_length": 2, "model": None if model is None else repr(model),
            "probe_weights_bp": [1.0],
            "batch": int(batch), "n_blocks": int(n_blocks),
            "n_replicates_requested": int(n_reps),
            "extra_sim_kwargs": sorted(extra_sim_kwargs or {}),
        },
    }
    return _finish(out, blocks, done, failed, distinct, dist_blocks,
                   last_error, pool, t0)


# --------------------------------------------------------- probe track
def probe_positions(sequence_length=TRACK_LENGTH, offsets=PROBE_OFFSETS):
    """Integer probe positions midpoint + offsets on [1, L-1)."""
    L = int(sequence_length)
    pos = (L // 2 + np.asarray(sorted(offsets), dtype=int)).astype(int)
    if pos.min() < 1 or pos.max() + 1 >= L or len(np.unique(pos)) != len(pos):
        raise ValueError("probe offsets must be distinct and land inside (0, L-1)")
    return pos


def probe_weights(offsets=PROBE_OFFSETS):
    """Trapezoid bp-widths of the sorted probe offsets (each probe owns half
    the gap to each neighbor; the two end probes are extended by half the
    terminal gap). A uniform-in-bp pooling weight; overall normalization is
    irrelevant for normalized spectra."""
    o = np.asarray(sorted(offsets), dtype=float)
    if len(o) == 1:
        return np.array([1.0])
    mid = (o[1:] + o[:-1]) / 2.0
    lo = np.concatenate([[o[0] - (o[1] - o[0]) / 2.0], mid])
    hi = np.concatenate([mid, [o[-1] + (o[-1] - o[-2]) / 2.0]])
    return hi - lo


def _track_windows(sequence_length, positions):
    L = int(sequence_length)
    P = np.asarray(positions, dtype=int)
    w = np.unique(np.concatenate([[0], P, P + 1, [L]])).astype(float)
    rows = np.searchsorted(w, P)
    return w, rows


def _accumulate(ts, rows, windows, M, n_hap):
    afs = ts.allele_frequency_spectrum(
        mode="branch", windows=windows, polarised=True, span_normalise=False)
    A = afs[rows, 1:n_hap]                      # (P, n-1) classes 1..n-1
    M += np.einsum("pi,pj->pij", A, A)


def _run_batches(sim_kwargs, n_reps, rng, rows, windows, M, batch):
    """Fixed-model arms: batched replicates. Returns (done, failed, last_error)."""
    msprime = _require_msprime()
    n_hap = _n_hap(sim_kwargs)
    done = failed = 0
    last_error = None
    while done + failed < n_reps:
        nb = min(batch, n_reps - done - failed)
        seed = int(rng.integers(1, 2**31 - 1))
        try:
            reps = msprime.sim_ancestry(
                num_replicates=nb, random_seed=seed, **sim_kwargs)
            for ts in reps:
                _accumulate(ts, rows, windows, M, n_hap)
            done += nb
        except Exception as e:                  # counted, never dropped
            failed += nb
            last_error = f"{type(e).__name__}: {e}"
    _all_failed_guard(done, failed, last_error)
    return done, failed, last_error


def _run_singles(model_fn, sim_kwargs, n_reps, rng, rows, windows, M):
    """Per-replicate randomized arms: one simulation per draw. model_fn(rng)
    -> (model, extra) must consume its parameter draws from rng BEFORE the
    seed draw, every call, so failures do not desynchronize the stream."""
    msprime = _require_msprime()
    n_hap = _n_hap(sim_kwargs)
    done = failed = 0
    draws = []
    last_error = None
    for _ in range(n_reps):
        model, extra = model_fn(rng)
        draws.append(extra)
        seed = int(rng.integers(1, 2**31 - 1))
        try:
            ts = msprime.sim_ancestry(random_seed=seed, model=model,
                                      **sim_kwargs)
            _accumulate(ts, rows, windows, M, n_hap)
            done += 1
        except Exception as e:                  # counted, never dropped
            failed += 1
            last_error = f"{type(e).__name__}: {e}"
    _all_failed_guard(done, failed, last_error)
    return done, failed, draws, last_error


def run_probe_track(n_reps, seed, sequence_length=TRACK_LENGTH,
                    probe_offsets=PROBE_OFFSETS, recomb_rate=R_BP, Ne=NE,
                    samples=SAMPLES, ploidy=PLOIDY, model=None, model_fn=None,
                    batch=BATCH_TRACK, extra_sim_kwargs=None):
    """Single-sequence probe track: one recombining sequence of length L
    with 1-bp probes at midpoint + probe_offsets; per probe p accumulates
    the SAME-window product M[p][i-1][j-1] += afs_i(p) * afs_j(p) (the
    completely linked object of the local tree at p) over n_reps replicates.
    Either a fixed `model` (batched) or a per-replicate `model_fn(rng) ->
    (model, extra)` (e.g. a sweep at a random position; its draws are
    returned under "draws"). Returns dict(M_per_probe (P, n-1, n-1),
    probe_positions, probe_weights_bp, pooled_per_probe, pooled (weights
    applied), done, failed, params, ...)."""
    t0 = time.time()
    if model is not None and model_fn is not None:
        raise ValueError("pass either model or model_fn, not both")
    L = int(sequence_length)
    pos = probe_positions(L, probe_offsets)
    windows, rows = _track_windows(L, pos)
    kw = sim_kwargs_for(L, Ne=Ne, samples=samples, ploidy=ploidy,
                        recomb_rate=recomb_rate, model=model,
                        extra_sim_kwargs=extra_sim_kwargs)
    n_hap = _n_hap(kw)
    root = root_seed(seed)
    rng = np.random.default_rng(root)
    M = np.zeros((len(pos), n_hap - 1, n_hap - 1))
    draws = None
    if model_fn is None:
        done, failed, last_error = _run_batches(kw, int(n_reps), rng, rows,
                                                windows, M, int(batch))
    else:
        kw.pop("model", None)
        done, failed, draws, last_error = _run_singles(
            model_fn, kw, int(n_reps), rng, rows, windows, M)
    wts = probe_weights(probe_offsets)
    out = {
        "geometry": "single-sequence probe track",
        "seed_root": root,
        "seed_law": "all msprime seeds drawn from default_rng(seed_root) in "
                    "code order (one per batch, or one per replicate after "
                    "model_fn's draws)",
        "params": {
            "n_haploid": n_hap, "samples": int(samples), "ploidy": int(ploidy),
            "Ne": float(Ne), "crossover_rate": (float(recomb_rate)
                                                if isinstance(recomb_rate, (int, float))
                                                else repr(recomb_rate)),
            "sequence_length": L, "model": None if model is None else repr(model),
            "probe_positions": pos.tolist(),
            "probe_offsets": [int(o) for o in sorted(probe_offsets)],
            "probe_weights_bp": wts.tolist(),
            "batch": int(batch) if model_fn is None else 1,
            "n_replicates_requested": int(n_reps),
            "extra_sim_kwargs": sorted(extra_sim_kwargs or {}),
        },
        "M_per_probe": M,
        "probe_positions": pos,
        "probe_weights_bp": wts,
        "done": int(done), "failed": int(failed), "last_error": last_error,
        "draws": draws,
        "keys": pair_keys(n_hap),
        "register": "simulation (msprime branch-mode Monte Carlo; Monte "
                    "Carlo error applies)",
        "object": "M[p][i-1][j-1] = sum over replicates of afs_i(p) * "
                  "afs_j(p), branch-mode polarized 1-bp probe windows, "
                  "classes i, j = 1..n-1 (E[L_i L_j] of the local genealogy)",
    }
    if done > 0:
        out["pooled_per_probe"] = np.array([pool_pairs(M[p]) for p in range(len(pos))])
        out["pooled"] = pool_pairs(np.einsum("p,pij->ij", wts, M))
    else:
        out["pooled_per_probe"] = None
        out["pooled"] = None
    out["elapsed_s"] = time.time() - t0
    return out


# aliases at the default n = 8, where the pooled vector has 28 entries
def keys28():
    """pair_keys(8)."""
    return pair_keys(8)


v28_from_M = pool_pairs
block_jackknife_v28 = block_jackknife

__all__ = [
    "SAMPLES", "PLOIDY", "N_HAP", "NE", "MARGIN_BP", "GC_TRACT_BP", "BATCH",
    "BATCH_TRACK", "R_BP", "TRACK_LENGTH", "SWEEP_DT", "PROBE_OFFSETS",
    "seed_from_name", "root_seed", "pair_keys", "pair_counts", "pool_pairs",
    "jackknife_cov", "block_jackknife", "margin_for", "layout",
    "window_edges", "sim_kwargs_for", "sweep_model", "run_two_window",
    "run_two_locus", "probe_positions", "probe_weights", "run_probe_track",
    "keys28", "v28_from_M", "block_jackknife_v28",
]


if __name__ == "__main__":
    import sys
    rho = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    n_reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    np.set_printoptions(precision=4, suppress=True, linewidth=120)
    r0 = run_two_locus(0.0, n_reps, seed="twowindow_sim.demo.rho0", n_blocks=10, batch=200)
    print(f"two-locus rho=0, n={r0['params']['n_haploid']}, {r0['done']} replicates "
          f"({r0['failed']} failed) in {r0['elapsed_s']:.2f}s; n_afs_distinct = "
          f"{r0['n_afs_distinct']} (must be 0)")
    print("pooled  :", r0["pooled"])
    print("se (B=10):", r0["pooled_se"])
    r1 = run_two_locus(rho, n_reps, seed="twowindow_sim.demo.rho", n_blocks=10, batch=200)
    print(f"two-locus rho={rho:g}: {r1['done']} replicates in {r1['elapsed_s']:.2f}s; "
          f"afs_distinct_fraction = {r1['afs_distinct_fraction']:.4f}")
    print("pooled  :", r1["pooled"])
    print("se (B=10):", r1["pooled_se"])

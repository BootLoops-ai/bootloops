#!/usr/bin/env python3
"""synth_truths — planted truths for site-frequency-spectrum estimators.

Two things live in this file:

1. An exact (closed-form, float) expected branch-length spectrum under a
   piecewise-constant population size,

       E[L_i],  i = 1..n-1,   epochs = [(t_start_gen, N_chrom), ...],

   computed from Tavare's (1984) ancestral-process transition function
   P(A_n(tau) = k) = sum_j c_{kj} exp(-j(j-1) tau / 2) integrated epoch by
   epoch, and E[L_i] = sum_k k P(a level-k line subtends i leaves) E[T_k]
   (the combinatorial weight is C(n-i-1, k-2)/C(n-1, k-1); Fu 1995). Time is
   in generations, ploidy-1 convention: with N chromosomes each pair
   coalesces at rate 1/N per generation, so tau = integral dt / N(t). Under a
   constant size N this is E[L_i] = 2N/i exactly, which the built-in
   analytic_selfchecks() asserts. The alternating sum is numerically safe for
   the small n it is used at here (n <= 10 or so); it is FLOAT, not certified.
   match_epoch_sizes() uses it to fit the epoch sizes of one history so that
   its normalized expected spectrum matches another history's at several n
   (an "SFS-matched pair": two different histories that a single-site
   spectrum cannot tell apart at those n; cf. Myers, Fefferman & Patterson
   2008 on histories with identical spectra).

2. Seeded msprime generators (Kelleher, Etheridge & McVean 2016; Baumdicker
   et al. 2022) that PLANT a truth and write synthetic data in the formats a
   two-site / SFS pipeline consumes, for named truth classes:

     TK1   Kingman, 3 epochs (bottleneck): sizes 20000 / 3000 / 16000
           chromosomes changing at 500 and 1500 generations
     TK2   Kingman, 4 epochs with change points 0/150/700/2500 generations
           and sizes FITTED (match_epoch_sizes) so that its normalized
           expected spectrum at n = 4 and 6 equals TK1's to <= 7.3e-5
           relative: the SFS-matched partner of TK1
     B13   Beta(2-alpha, alpha) coalescent, alpha = 1.3  (Schweinsberg 2003)
     B17   Beta coalescent, alpha = 1.7
     D005  Dirac (psi) coalescent, psi = 0.05, c = 10  (Eldon & Wakeley 2006,
           in msprime's DiracCoalescent(psi, c) parametrization)
     D020  Dirac coalescent, psi = 0.2, c = 10
   all at ploidy 1 (a population of NPOP = 20000 chromosomes for the
   multiple-merger classes). Further classes and arms can be declared in the
   pins file (see make_pins / class_params / arm_spec): kind
   'kingman_epochs' with 'epochs', kind 'beta' with 'alpha' and 'N', kind
   'dirac' with 'psi', 'c' and 'N'.

   Planted CLASS truth (run_truth): for each n in TRUTH_NS = (4, 6), from
   `reps` independent single-locus msprime trees (branch mode, no mutation),
       q_ij = E[L_i L_j] / E[L_tot^2]   and   xi_i = E[L_i] / E[L_tot]
   as ratio-of-means estimates with delete-one-jackknife Monte Carlo standard
   errors per cell, plus E[L_tot]. q_ij is the normalized second-moment
   matrix of the branch-length spectrum, i.e. the two-site frequency
   spectrum of a pair of completely linked sites (distance -> 0) up to
   normalization; because the single-locus law does not depend on the
   recombination rate, one truth per class serves every arm of that class.

   ARMS (run_arm): class x {recombination ON, OFF} over a REGION_L = 5 Mb
   region, plus one gene-conversion arm TK1_rON_GC (1 Mb, mean tract 300 bp,
   non-crossover initiation rate = 2 x the crossover rate). Each arm draws
   N_DIP = 500 "diploid" individuals = 1000 exchangeable chromosomes paired
   (2k, 2k+1) (a ploidy-1 coalescent, NOT a diploid multiple-merger model;
   unphased genotype statistics are blind to the pairing), calibrates the
   per-bp mutation rate mu from the class's own mean total branch length at
   that sample size so that the expected density of segregating sites is
   TARGET_DENSITY = 3.7e-3 per bp (Watterson 1975: the density a sample of
   1000 chromosomes has under mu = 1.25e-8 in a constant population of
   2 x 10^4 chromosomes, 2 N mu H_999), sets the crossover
   rate r = mu x (r_map / MU_HUMAN) so that r/mu keeps the ratio of a given
   map-average rate to the human mutation rate, simulates ancestry with
   recombination (and gene conversion), drops infinite-sites binary
   mutations, and writes
     <ARM>.gt.tsv.gz   VCF-like genotype table: '##' meta lines, then
                       #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S0001..
                       with unphased GT 'a/b' (REF = ancestral = 0), contig
                       'chr21S', 1-based positions, random REF/ALT letters;
                       multi-hit and non-segregating sites are dropped and
                       counted
     <ARM>.anc.tsv.gz  the TRUE polarization table 'pos ref alt anc conf'
                       (anc == ref by construction, conf = high), the schema
                       popcorn.ancestral writes for real data
     arm_<ARM>.json    parameters (mu, r, gene-conversion rate, calibration),
                       the four seeds, site counts, file sizes and hashes
                       (the .gz container hash AND the mtime-free payload
                       hash, see gz_payload_sha256).
   For a constant Kingman population of 2 x 10^4 chromosomes at the default
   500 individuals the calibration returns mu ~= 1.25e-8 and r ~= r_map;
   other histories rescale mu (TK1: ~1.9e-8), and for the multiple-merger
   classes the absolute per-generation rates are synthetic (B13: ~5.6e-7):
   only the population-scaled intensities (site density, r/mu) are
   meaningful there.

   Seeds: every random stream is derived from one base seed by
   derive_seed(base, tag) = sha256('<base>:<tag>')[:4] -> 1..2^31-2 with
   tags 'truth:<cls>:n<n>' and 'arm:<arm>:{calib,ancestry,mutation,refalt}'
   (the reduced-scale selftest prefixes 'selftest:'), so any single object
   can be regenerated alone. The base seed, the map-average rate, the fitted
   TK2 sizes and the class/arm tables are written ONCE to a pins file
   (make_pins) before anything is generated; every generator takes the pins
   as input.

REGISTER: SIMULATION / DATA-UTILITY. The planted q_ij and xi_i are Monte
Carlo expectations known to their quoted jackknife SE (20000 replicates for
the multiple-merger classes and 2 x 10^6 for TK1/TK2 in the shipped
reference set), not exact values; the only closed-form object is
expected_branch_sfs and it is double precision. Never quote a planted q_ij
as exact, and label it when it appears beside popcorn's exact or certified
numbers. Bit-for-bit regeneration of the shipped reference set is a property
of the recorded simulator versions (pins.json "versions"); under other
msprime/numpy versions only agreement within the Monte Carlo errors is
expected, and the battery says which mode it checked.

Reference data (reference/planted/, read-only; regenerable from pins.json):
  pins.json              base seed, seed law, pinned constants, map-average
                         rate used for r/mu, fitted TK2 sizes with the
                         achieved match, class and arm tables, simulator
                         versions the set was generated with
  truth_<CLS>.json       planted q_ij / xi_i / SEs / E[L_tot] at n = 4, 6 for
                         the six classes (TK2 also carries its two-site
                         distance from TK1: at these n it is within the
                         Monte Carlo error, max 1.8 SE at 2 x 10^6 replicates)
  arm_TK1_rON_GC.json    the gene-conversion arm at production scale:
                         parameters, seeds, counts, payload hashes
  selftest_expected.json the reduced-scale selftest's expected numbers
                         (exact E[L_i] for TK1 at n = 6, the seeded z-check,
                         the subscale arm's counts and payload hashes)

Cost (single core): expected_branch_sfs is instantaneous; a class truth at
20000 replicates ~2 s per n; TK1/TK2 at 2 x 10^6 replicates ~4 min per n; an
arm at 500 individuals is dominated by writing the table (~2-3 s per Mb, up
to ~16 s for a 5 Mb Dirac-class arm); the reduced selftest < 1 s. Memory:
the jackknife holds reps x (n-1)^2 doubles (400 MB at 2 x 10^6, n = 6).

Scope: neutral, single population, piecewise-constant sizes or the two
multiple-merger families above; infinite-sites binary mutations; flat
recombination and gene-conversion rates per arm; truths are CLEAN (no
genotype error, no mis-polarization: perturbations belong downstream).
Nothing here reads real data except an optional genetic-map bedGraph for
its length-weighted mean rate (measure_map_mean).

CLI (outputs only under --out; nothing is written next to this file):
  python3 synth_truths.py list
  python3 synth_truths.py check                      # analytic checks vs reference, no msprime
  python3 synth_truths.py pin   --out DIR [--seed S] [--map BEDGRAPH --chrom chr21 | --r-per-bp R] [--extra JSON]
  python3 synth_truths.py truth --cls B13 --out DIR [--pins FILE] [--reps 20000]
  python3 synth_truths.py arm   --arm TK1_rON_GC --out DIR [--pins FILE] [--n-dip 500] [--length L]
  python3 synth_truths.py all   --out DIR [--pins FILE] [--reps 20000]
  python3 synth_truths.py selftest [--out DIR] [--pins FILE]   # reduced-scale end to end, needs msprime
--pins defaults to the shipped reference pins, so `truth --cls B13 --out D`
regenerates reference/planted/truth_B13.json from its recorded seeds.

Requirements: numpy (always); scipy (match_epoch_sizes, hence `pin`);
msprime + tskit (truths, arms, selftest). The analytic functions and
`check` need numpy only.

References
  Tavare, S. (1984). Line-of-descent and genealogical processes, and their
    applications in population genetics models. Theor. Popul. Biol. 26,
    119-164.
  Fu, Y.-X. (1995). Statistical properties of segregating sites. Theor.
    Popul. Biol. 48, 172-197.
  Watterson, G. A. (1975). On the number of segregating sites in genetical
    models without recombination. Theor. Popul. Biol. 7, 256-276.
  Polanski, A. & Kimmel, M. (2003). New explicit expressions for relative
    frequencies of single-nucleotide polymorphisms with application to
    statistical inference on population growth. Genetics 165, 427-436.
  Myers, S., Fefferman, C. & Patterson, N. (2008). Can one learn history
    from the allelic spectrum? Theor. Popul. Biol. 73, 342-348.
  Schweinsberg, J. (2003). Coalescent processes obtained from supercritical
    Galton-Watson processes. Stochastic Process. Appl. 106, 107-139.
  Eldon, B. & Wakeley, J. (2006). Coalescent processes when the distribution
    of offspring number among individuals is highly skewed. Genetics 172,
    2621-2633.
  Kelleher, J., Etheridge, A. M. & McVean, G. (2016). Efficient coalescent
    simulation and genealogical analysis for large sample sizes. PLoS
    Comput. Biol. 12, e1004842.
  Baumdicker, F. et al. (2022). Efficient ancestry and mutation simulation
    with msprime 1.0. Genetics 220, iyab229.
"""

import argparse
import gzip
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
import time

import numpy as np

# ----------------------------------------------------------------------------
# paths + pinned constants
# ----------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.join(HERE, "reference", "planted")
PINS_REF = os.path.join(REFERENCE_DIR, "pins.json")

REGISTER = ("simulation: ploidy-1 Kingman/Lambda coalescents (msprime); "
            "diploid individuals = exchangeable-paired chromosomes")

MU_HUMAN = 1.25e-8          # per bp per gen; the r/mu-ratio denominator
TARGET_DENSITY = 3.7e-3     # kept-SNV target per bp (Watterson n=1000,
                            # Ne_dip=1e4, mu=1.25e-8 -> 2*2e4*H_999*mu)
N_DIP = 500                 # diploid individuals per arm
REGION_L = 5_000_000        # bp, region length of the standard arms
GC_REGION_L = 1_000_000     # bp, region length of the gene-conversion arm
GC_TRACT = 300.0            # bp, mean gene-conversion tract length
GC_NCO_FACTOR = 2.0         # non-crossover initiation rate = 2 x crossover r
DIRAC_C = 10.0              # msprime DiracCoalescent c for the Dirac classes
NPOP = 20000.0              # ploidy-1 chromosome-population size (= Ne_dip 1e4)
TRUTH_REPS_DEFAULT = 20000  # branch-mode replicates per class per n
TRUTH_NS = (4, 6)           # sample sizes of the planted class truths
CALIB_REPS = 8              # single-tree reps for the mu calibration
FALLBACK_R_MAP = 1.2e-8     # per bp per gen, used when no map rate is given

# T-K histories (sizes are ploidy-1 chromosome counts; diploid-equivalent /2)
TK1_EPOCHS = [(0.0, 20000.0), (500.0, 3000.0), (1500.0, 16000.0)]
TK2_TIMES = [0.0, 150.0, 700.0, 2500.0]   # change-points; sizes fitted at pin
TK2_LOGN_BOUNDS = (math.log(500.0), math.log(200000.0))
TK2_X0 = math.log(12000.0)

CLASSES = ("TK1", "TK2", "B13", "B17", "D005", "D020")

ARMS = {}
for _cls in CLASSES:
    ARMS[_cls + "_rON"] = {"cls": _cls, "recomb": True, "gc": False,
                           "L": REGION_L}
    ARMS[_cls + "_rOFF"] = {"cls": _cls, "recomb": False, "gc": False,
                            "L": REGION_L}
ARMS["TK1_rON_GC"] = {"cls": "TK1", "recomb": True, "gc": True,
                      "L": GC_REGION_L}

BASES = ("A", "C", "G", "T")


# ----------------------------------------------------------------------------
# small utilities
# ----------------------------------------------------------------------------

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gz_payload_sha256(path):
    """sha256 of the decompressed table of a .gz written by run_arm: every
    line except '##' meta lines (i.e. the '#CHROM'/'pos' header and the data
    rows). Unlike the file hash it does not depend on the gzip container's
    mtime, so it identifies a regenerated table."""
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b"##"):
                continue
            h.update(line)
    return h.hexdigest()


def derive_seed(base_seed, tag):
    """Deterministic per-purpose seed: sha256(base:tag) -> 1..2^31-2."""
    d = hashlib.sha256(f"{base_seed}:{tag}".encode()).digest()
    return int.from_bytes(d[:4], "big") % (2**31 - 2) + 1


def write_json(path, obj):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1)
        fh.write("\n")
    return path


def load_pins(path=None):
    """The pins a generator runs under (dict). Default: the shipped reference
    pins (reference/planted/pins.json); pass the pins.json that make_pins
    wrote to generate under your own base seed / map rate / classes."""
    path = path or PINS_REF
    if not os.path.exists(path):
        raise FileNotFoundError(
            "%s missing: run 'pin --out DIR' first (the base seed and the "
            "fitted sizes must be written BEFORE anything is generated)"
            % path)
    with open(path) as fh:
        return json.load(fh)


def load_reference_truth(cls):
    """The shipped planted truth of class `cls` (reference/planted/truth_<cls>.json)."""
    with open(os.path.join(REFERENCE_DIR, f"truth_{cls}.json")) as fh:
        return json.load(fh)


def load_reference_arm(arm="TK1_rON_GC"):
    """The shipped arm record (reference/planted/arm_<arm>.json)."""
    with open(os.path.join(REFERENCE_DIR, f"arm_{arm}.json")) as fh:
        return json.load(fh)


def versions():
    """Interpreter and library versions (None for a library that is absent);
    recorded in pins.json because bit-for-bit regeneration depends on them."""
    out = {"python": sys.version.split()[0], "numpy": np.__version__}
    for name in ("scipy", "msprime", "tskit"):
        try:
            out[name] = __import__(name).__version__
        except ImportError:
            out[name] = None
    return out


# ----------------------------------------------------------------------------
# exact branch-expected 1-SFS under piecewise-constant N (ploidy-1 units)
# Tavare (1984) ancestral-process alternating sum; n <= 6 here, so stable.
# ----------------------------------------------------------------------------

def _rising(a, b):
    out = 1.0
    for m in range(b):
        out *= a + m
    return out


def _falling(a, b):
    out = 1.0
    for m in range(b):
        out *= a - m
    return out


def tavare_coefs(n):
    """P(A_n(tau)=k) = sum_j coef * exp(-lam*tau); returns {k: [(coef, lam)]}.
    tau = standard coalescent time (pair rate 1)."""
    coefs = {}
    for k in range(2, n + 1):
        terms = []
        for j in range(k, n + 1):
            lam = j * (j - 1) / 2.0
            c = (2 * j - 1) * ((-1.0) ** (j - k))
            c *= _rising(k, j - 1) * _falling(n, j)
            c /= math.factorial(k) * math.factorial(j - k) * _rising(n, j)
            terms.append((c, lam))
        coefs[k] = terms
    return coefs


def expected_branch_sfs(epochs, n):
    """E[L_i] (generations) for i=1..n-1 under piecewise-constant N(t).
    epochs: [(t_start_gen, N_chrom)] with t_start[0]=0; last epoch infinite.
    ploidy-1 convention: pair coalescence rate 1/N per generation."""
    coefs = tavare_coefs(n)
    # tau boundaries of the epochs (standard-time rescale tau = int dt/N)
    taus = [0.0]
    for m in range(1, len(epochs)):
        dt = epochs[m][0] - epochs[m - 1][0]
        taus.append(taus[-1] + dt / epochs[m - 1][1])
    # E[T_k] in generations = sum_epochs N_m * int_{tau_m}^{tau_{m+1}} P_k
    ETk = {}
    for k in range(2, n + 1):
        tot = 0.0
        for m, (_, Nm) in enumerate(epochs):
            ta = taus[m]
            tb = taus[m + 1] if m + 1 < len(taus) else None
            for c, lam in coefs[k]:
                ea = math.exp(-lam * ta)
                eb = math.exp(-lam * tb) if tb is not None else 0.0
                tot += Nm * c / lam * (ea - eb)
        ETk[k] = tot
    EL = np.zeros(n - 1)
    for i in range(1, n):
        acc = 0.0
        for k in range(2, n + 1):
            p = (math.comb(n - i - 1, k - 2) / math.comb(n - 1, k - 1)
                 if 0 <= k - 2 <= n - i - 1 else 0.0)
            acc += k * p * ETk[k]
        EL[i - 1] = acc
    return EL


def analytic_selfchecks():
    """Exact-machinery guards: constant-N 1/i law and n=2 survival."""
    n = 6
    N = 12345.0
    EL = expected_branch_sfs([(0.0, N)], n)
    for i in range(1, n):
        assert abs(EL[i - 1] - 2.0 * N / i) < 1e-6 * N, (i, EL[i - 1])
    # n=2: P(A=2)(tau) must be exp(-tau)
    (c, lam), = tavare_coefs(2)[2]
    assert abs(c - 1.0) < 1e-12 and abs(lam - 1.0) < 1e-12
    return True


def match_epoch_sizes(target_epochs, times, ns=TRUTH_NS):
    """Coarse moment-match: fit the epoch sizes of a history with change
    points `times` so that its normalized branch-expected 1-SFS at n in ns
    matches that of `target_epochs`. Returns (sizes, report). Needs scipy."""
    from scipy.optimize import least_squares

    targets = {}
    for n in ns:
        xi = expected_branch_sfs(target_epochs, n)
        targets[n] = xi / xi.sum()

    def resid(logNs):
        epochs = [(times[m], math.exp(logNs[m]))
                  for m in range(len(times))]
        out = []
        for n in ns:
            xi = expected_branch_sfs(epochs, n)
            xi = xi / xi.sum()
            out.extend(((xi - targets[n]) / targets[n]).tolist())
        return np.array(out)

    x0 = np.full(len(times), TK2_X0)
    fit = least_squares(resid, x0, bounds=(TK2_LOGN_BOUNDS[0],
                                           TK2_LOGN_BOUNDS[1]),
                        xtol=1e-14, ftol=1e-14, gtol=1e-14)
    sizes = [float(math.exp(v)) for v in fit.x]
    epochs2 = [(times[m], sizes[m]) for m in range(len(times))]
    report = {"converged": bool(fit.success), "cost": float(fit.cost),
              "achieved": {}}
    for n in ns:
        xi1 = expected_branch_sfs(target_epochs, n)
        xi1 = xi1 / xi1.sum()
        xi2 = expected_branch_sfs(epochs2, n)
        xi2 = xi2 / xi2.sum()
        report["achieved"][f"n{n}"] = {
            "xi_target": xi1.tolist(), "xi_fitted": xi2.tolist(),
            "max_abs_rel_diff": float(np.max(np.abs(xi2 - xi1) / xi1))}
    return sizes, report


# ----------------------------------------------------------------------------
# msprime model plumbing
# ----------------------------------------------------------------------------

def class_params(cls, pins):
    """Model of truth class `cls`: the six built-in classes, else a class
    declared under pins['classes'] ({'kind': 'kingman_epochs', 'epochs':
    [[t_gen, N_chrom], ...]} | {'kind': 'beta', 'alpha': a, 'N': N} |
    {'kind': 'dirac', 'psi': p, 'c': c, 'N': N})."""
    if cls == "TK1":
        return {"kind": "kingman_epochs", "epochs": TK1_EPOCHS}
    if cls == "TK2":
        sizes = pins["tk2_fitted_sizes"]
        return {"kind": "kingman_epochs",
                "epochs": [(TK2_TIMES[m], sizes[m])
                           for m in range(len(TK2_TIMES))]}
    if cls == "B13":
        return {"kind": "beta", "alpha": 1.3, "N": NPOP}
    if cls == "B17":
        return {"kind": "beta", "alpha": 1.7, "N": NPOP}
    if cls == "D005":
        return {"kind": "dirac", "psi": 0.05, "c": DIRAC_C, "N": NPOP}
    if cls == "D020":
        return {"kind": "dirac", "psi": 0.2, "c": DIRAC_C, "N": NPOP}
    spec = (pins or {}).get("classes", {}).get(cls)
    if spec is not None and spec.get("kind") in ("kingman_epochs", "beta",
                                                 "dirac"):
        return dict(spec)
    raise ValueError("unknown truth class %r (built-in: %s; or declare it "
                     "under pins['classes'])" % (cls, ", ".join(CLASSES)))


def arm_spec(arm, pins=None):
    """{'cls', 'recomb', 'gc', 'L'} of arm `arm`: the built-in arms, else an
    arm declared under pins['arms'] with the same keys."""
    if arm in ARMS:
        return ARMS[arm]
    spec = (pins or {}).get("arms", {}).get(arm)
    if spec is not None and {"cls", "recomb", "gc", "L"} <= set(spec):
        return spec
    raise ValueError("unknown arm %r (built-in: %s; or declare it under "
                     "pins['arms'] with keys cls, recomb, gc, L)"
                     % (arm, ", ".join(sorted(ARMS))))


def sim_kwargs(cp):
    """msprime.sim_ancestry kwargs for a class-params dict (ploidy=1)."""
    import msprime
    if cp["kind"] == "kingman_epochs":
        dem = msprime.Demography()
        dem.add_population(name="p0", initial_size=cp["epochs"][0][1])
        for t, Nv in cp["epochs"][1:]:
            dem.add_population_parameters_change(time=t, initial_size=Nv,
                                                 population="p0")
        return {"demography": dem}
    if cp["kind"] == "beta":
        return {"population_size": cp["N"],
                "model": msprime.BetaCoalescent(alpha=cp["alpha"])}
    if cp["kind"] == "dirac":
        return {"population_size": cp["N"],
                "model": msprime.DiracCoalescent(psi=cp["psi"], c=cp["c"])}
    raise ValueError(cp)


def branch_class_lengths(tree, n):
    """L_i (i=1..n-1): total branch length subtending i of the n samples."""
    import tskit
    L = np.zeros(n - 1)
    for u in tree.nodes():
        pa = tree.parent(u)
        if pa == tskit.NULL:
            continue
        k = tree.num_samples(u)
        if 1 <= k <= n - 1:
            L[k - 1] += tree.time(pa) - tree.time(u)
    return L


# ----------------------------------------------------------------------------
# pin step (seeds + params written BEFORE any generation)
# ----------------------------------------------------------------------------

def measure_map_mean(map_path, chrom="chr21"):
    """Length-weighted mean recombination rate (cM/Mb -> per bp per
    generation) of chromosome `chrom` in a bedGraph genetic map
    ('chrom start end rate_cM_per_Mb' per line, plain or .gz)."""
    if not map_path:
        raise ValueError("no genetic map given: pass map_path / --map PATH "
                         "(bedGraph 'chrom start end cM/Mb'), or give the "
                         "rate directly (r_per_bp / --r-per-bp)")
    s = w = 0.0
    rows = 0
    opener = gzip.open if str(map_path).endswith(".gz") else open
    with opener(map_path, "rt") as fh:
        for line in fh:
            f = line.split()
            if f[0] != chrom:
                continue
            span = float(f[2]) - float(f[1])
            s += span * float(f[3])
            w += span
            rows += 1
    if w <= 0:
        raise ValueError("no rows for chromosome %r in %s" % (chrom, map_path))
    mean_cm_mb = s / w
    return {"source": os.path.basename(str(map_path)),
            "source_sha256": sha256_file(map_path),
            "chrom": chrom, "rows": rows, "bp_covered": int(w),
            "mean_cM_per_Mb": mean_cm_mb,
            "r_map_per_bp": mean_cm_mb * 1e-8,
            "usage": "mean rate only: synthetic arms use a FLAT rate "
                     "r = mu * r_map_per_bp / mu_human"}


def make_pins(out_dir, seed=None, map_path=None, chrom="chr21",
              r_per_bp=None, classes=None, arms=None):
    """Write out_dir/pins.json (never overwritten: an existing file is
    returned as is) and return the pins dict. seed: base seed (default: 8
    random bytes). The map-average rate comes from `map_path` (bedGraph,
    see measure_map_mean), else from `r_per_bp`, else FALLBACK_R_MAP.
    classes / arms: extra {name: spec} tables merged into the pins (see
    class_params / arm_spec for the spec keys). Needs scipy (TK2 fit)."""
    path = os.path.join(out_dir, "pins.json")
    os.makedirs(out_dir, exist_ok=True)
    if os.path.exists(path):
        print("pins already written (never regenerated):", path)
        return load_pins(path)
    analytic_selfchecks()
    base_seed = int(seed) if seed is not None else \
        int.from_bytes(os.urandom(8), "big")
    if map_path:
        map_info = measure_map_mean(map_path, chrom)
    else:
        rate = float(r_per_bp) if r_per_bp is not None else FALLBACK_R_MAP
        map_info = {"source": "given" if r_per_bp is not None else "default",
                    "r_map_per_bp": rate,
                    "usage": "synthetic arms use a FLAT rate "
                             "r = mu * r_map_per_bp / mu_human"}
    tk2_sizes, match_report = match_epoch_sizes(TK1_EPOCHS, TK2_TIMES)
    all_classes = {c: class_params(c, {"tk2_fitted_sizes": tk2_sizes})
                   for c in CLASSES}
    for name, spec in (classes or {}).items():
        if name in all_classes:
            raise ValueError("class %r is built in" % name)
        all_classes[name] = class_params(name, {"classes": {name: spec}})
    all_arms = dict(ARMS)
    for name, spec in (arms or {}).items():
        if name in all_arms:
            raise ValueError("arm %r is built in" % name)
        spec = arm_spec(name, {"arms": {name: spec}})
        if spec["cls"] not in all_classes:
            raise ValueError("arm %r: unknown class %r" % (name, spec["cls"]))
        all_arms[name] = spec
    pins = {
        "generator": "popcorn synth_truths",
        "versions": versions(),
        "register": REGISTER,
        "base_seed": base_seed,
        "seed_law": "sha256('<base_seed>:<tag>')[:4] -> 1..2^31-2; tags: "
                    "truth:<cls>:n<n>, arm:<arm>:{calib,ancestry,mutation,"
                    "refalt}; selftest stream prefixes tags with 'selftest:'",
        "pinned": {
            "mu_human": MU_HUMAN, "target_density": TARGET_DENSITY,
            "n_dip": N_DIP, "region_L": REGION_L, "gc_region_L": GC_REGION_L,
            "gc_tract_bp": GC_TRACT, "gc_nco_factor": GC_NCO_FACTOR,
            "dirac_c": DIRAC_C, "npop_chrom": NPOP,
            "truth_reps_default": TRUTH_REPS_DEFAULT, "truth_ns": list(TRUTH_NS),
            "tk1_epochs": TK1_EPOCHS, "tk2_times": TK2_TIMES},
        "recomb_law": "r_arm = mu_arm * (r_map_per_bp / mu_human); "
                      "mu_arm = target_density / mean single-tree total "
                      "branch length at the arm's n (CALIB_REPS=%d reps); "
                      "for constant-N=2e4 Kingman at n=1000 this recovers "
                      "mu ~= 1.25e-8 and r ~= r_map; multi-epoch and MMC "
                      "classes rescale mu (and r with it) to hold the site "
                      "density and the r/mu ratio — the population-scaled "
                      "pins" % CALIB_REPS,
        "map": map_info,
        "tk2_fitted_sizes": tk2_sizes,
        "tk2_match": match_report,
        "classes": all_classes,
        "arms": all_arms,
        "note": "written BEFORE any truth/arm generation; never regenerated",
    }
    write_json(path, pins)
    print("pinned:", path)
    for n in TRUTH_NS:
        print("  tk2 match n=%d max|rel diff| = %.3e" %
              (n, match_report["achieved"][f"n{n}"]["max_abs_rel_diff"]))
    return pins


# ----------------------------------------------------------------------------
# planted class truths: q_ij at n=4,6 from branch-mode replicates
# ----------------------------------------------------------------------------

def run_truth(cls, pins, reps, out_dir, tag_prefix=""):
    """Planted truth of class `cls` from `reps` msprime replicates per n in
    TRUTH_NS under `pins`; writes out_dir/truth_<cls>.json unless out_dir is
    None; returns the record."""
    import msprime
    cp = class_params(cls, pins)
    kwargs = sim_kwargs(cp)
    blocks = {}
    t_all = time.time()
    for n in TRUTH_NS:
        seed = derive_seed(pins["base_seed"], f"{tag_prefix}truth:{cls}:n{n}")
        t0 = time.time()
        X = np.empty((reps, n - 1))
        it = msprime.sim_ancestry(samples=n, ploidy=1, num_replicates=reps,
                                  random_seed=seed, **kwargs)
        for idx, ts in enumerate(it):
            X[idx] = branch_class_lengths(ts.first(), n)
        tot = X.sum(axis=1)
        S2 = X.T @ X
        T2 = float(np.sum(tot ** 2))
        S1 = X.sum(axis=0)
        T1 = float(tot.sum())
        q = S2 / T2
        xi = S1 / T1
        # delete-one jackknife for q and xi
        qj = np.empty((reps, n - 1, n - 1))
        xj = np.empty((reps, n - 1))
        for r in range(reps):
            o = np.outer(X[r], X[r])
            qj[r] = (S2 - o) / (T2 - tot[r] ** 2)
            xj[r] = (S1 - X[r]) / (T1 - tot[r])
        fac = (reps - 1) / reps
        q_se = np.sqrt(fac * np.sum((qj - qj.mean(axis=0)) ** 2, axis=0))
        xi_se = np.sqrt(fac * np.sum((xj - xj.mean(axis=0)) ** 2, axis=0))
        blocks[f"n{n}"] = {
            "n": n, "reps": reps, "seed": seed,
            "q_ij": q.tolist(), "q_ij_mc_se": q_se.tolist(),
            "xi_1sfs": xi.tolist(), "xi_1sfs_mc_se": xi_se.tolist(),
            "E_Ltot_gen": T1 / reps,
            "wall_s": round(time.time() - t0, 3),
            "index_law": "row/col index m is derived-count i=m+1; "
                         "q_ij = E[L_i L_j]/E[L_tot^2], single-locus "
                         "(no-recomb) marginal = the d->0 planted truth"}
    rec = {"generator": "popcorn synth_truths",
           "kind": "class_truth", "cls": cls,
           "model": {k: v for k, v in class_params(cls, pins).items()},
           "register": pins.get("register", REGISTER),
           "mc": blocks, "wall_s_total": round(time.time() - t_all, 3),
           "subscale": bool(reps < TRUTH_REPS_DEFAULT)}
    if reps < TRUTH_REPS_DEFAULT:
        rec["subscale_note"] = ("reps below production default %d — selftest "
                                "register only" % TRUTH_REPS_DEFAULT)
    # SFS-matched-pair evidence: quantify the 2-SFS difference vs TK1
    if cls == "TK2" and out_dir is not None:
        tk1_path = os.path.join(out_dir, "truth_TK1.json")
        if os.path.exists(tk1_path):
            with open(tk1_path) as fh:
                tk1 = json.load(fh)
            vs = {}
            for n in TRUTH_NS:
                b1 = tk1["mc"].get(f"n{n}")
                b2 = blocks[f"n{n}"]
                if b1 is None:
                    continue
                d = np.abs(np.array(b2["q_ij"]) - np.array(b1["q_ij"]))
                se = np.sqrt(np.array(b2["q_ij_mc_se"]) ** 2 +
                             np.array(b1["q_ij_mc_se"]) ** 2)
                vs[f"n{n}"] = {"max_abs_dq": float(d.max()),
                               "max_z": float((d / se).max()),
                               "z_matrix": (d / se).tolist()}
            rec["vs_TK1_2sfs"] = vs
    if out_dir is not None:
        path = write_json(os.path.join(out_dir, f"truth_{cls}.json"), rec)
        print("truth record:", path)
    return rec


# ----------------------------------------------------------------------------
# arm generation: genotype TSV + TRUE polarization TSV
# ----------------------------------------------------------------------------

def run_arm(arm, pins, out_truths, out_records, n_dip=None, L=None,
            tag_prefix=""):
    """Generate arm `arm` under `pins`: writes <arm>.gt.tsv.gz and
    <arm>.anc.tsv.gz into out_truths and arm_<arm>.json into out_records
    (the two may be the same directory); n_dip / L override the arm's
    individuals / region length (reduced-scale runs). Returns the record."""
    import msprime
    spec = arm_spec(arm, pins)
    cls = spec["cls"]
    n_dip = n_dip or N_DIP
    L = L or spec["L"]
    n_hap = 2 * n_dip
    cp = class_params(cls, pins)
    kwargs = sim_kwargs(cp)
    base = pins["base_seed"]
    seeds = {k: derive_seed(base, f"{tag_prefix}arm:{arm}:{k}")
             for k in ("calib", "ancestry", "mutation", "refalt")}
    timings = {}

    # 1) mu calibration: mean single-tree total branch length at this n
    t0 = time.time()
    it = msprime.sim_ancestry(samples=n_hap, ploidy=1,
                              num_replicates=CALIB_REPS,
                              random_seed=seeds["calib"], **kwargs)
    btots = np.array([ts.first().total_branch_length for ts in it])
    bhat = float(btots.mean())
    mu = TARGET_DENSITY / bhat
    rratio = pins["map"]["r_map_per_bp"] / MU_HUMAN
    r = mu * rratio if spec["recomb"] else 0.0
    timings["calib_s"] = round(time.time() - t0, 3)

    # 2) ancestry
    t0 = time.time()
    gc_kwargs = {}
    if spec["gc"]:
        gc_kwargs = {"gene_conversion_rate": GC_NCO_FACTOR * r,
                     "gene_conversion_tract_length": GC_TRACT}
    ts = msprime.sim_ancestry(samples=n_hap, ploidy=1, sequence_length=L,
                              recombination_rate=r,
                              random_seed=seeds["ancestry"],
                              **gc_kwargs, **kwargs)
    timings["ancestry_s"] = round(time.time() - t0, 3)

    # 3) mutations (binary model; ancestral '0', derived '1')
    t0 = time.time()
    mts = msprime.sim_mutations(ts, rate=mu,
                                model=msprime.BinaryMutationModel(),
                                random_seed=seeds["mutation"])
    timings["mutations_s"] = round(time.time() - t0, 3)

    # 4) write the tables
    t0 = time.time()
    rng = np.random.default_rng(seeds["refalt"])
    os.makedirs(out_truths, exist_ok=True)
    gt_path = os.path.join(out_truths, f"{arm}.gt.tsv.gz")
    anc_path = os.path.join(out_truths, f"{arm}.anc.tsv.gz")
    n_multi = n_nonseg = n_kept = n_single = n_double = 0
    sample_ids = ["S%04d" % (k + 1) for k in range(n_dip)]
    with gzip.open(gt_path, "wt") as gfh, gzip.open(anc_path, "wt") as afh:
        gfh.write("##fileformat=VCF-like_TSV_synthetic\n")
        gfh.write("##source=synth_truths.py arm=%s cls=%s (SYNTHETIC — no "
                  "real data)\n" % (arm, cls))
        gfh.write("##contig=<ID=chr21S,length=%d>\n" % int(L))
        gfh.write("##arm_record=arm_%s.json\n" % arm)
        gfh.write("##polarization_truth=%s\n" % os.path.basename(anc_path))
        gfh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(sample_ids) + "\n")
        afh.write("pos\tref\talt\tanc\tconf\n")
        for var in mts.variants():
            site = var.site
            if len(site.mutations) != 1:
                n_multi += 1
                continue
            didx = var.alleles.index(site.mutations[0].derived_state)
            g = (var.genotypes == didx).astype(np.int8)
            dc = int(g.sum())
            if dc < 1 or dc > n_hap - 1:
                n_nonseg += 1
                continue
            n_kept += 1
            n_single += dc == 1
            n_double += dc == 2
            pos = int(site.position) + 1
            ref = BASES[rng.integers(4)]
            alt = rng.choice([b for b in BASES if b != ref])
            row = "\t".join("%d/%d" % (min(g[2 * k], g[2 * k + 1]),
                                       max(g[2 * k], g[2 * k + 1]))
                            for k in range(n_dip))
            gfh.write("chr21S\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n"
                      % (pos, ref, alt, row))
            afh.write("%d\t%s\t%s\t%s\thigh\n" % (pos, ref, alt, ref))
    timings["write_s"] = round(time.time() - t0, 3)

    rec = {"generator": "popcorn synth_truths", "kind": "arm",
           "arm": arm, "cls": cls,
           "class_truth_record": f"truth_{cls}.json",
           "model": cp, "register": pins.get("register", REGISTER),
           "params": {"n_dip": n_dip, "n_hap": n_hap, "region_L": int(L),
                      "mu_per_bp_gen": mu, "r_per_bp_gen": r,
                      "recomb_on": spec["recomb"], "gc_on": spec["gc"],
                      "gc_rate": GC_NCO_FACTOR * r if spec["gc"] else 0.0,
                      "gc_tract_bp": GC_TRACT if spec["gc"] else None,
                      "target_density": TARGET_DENSITY,
                      "calib_bhat_gen": bhat,
                      "calib_btot_cv": float(btots.std() / bhat),
                      "r_over_mu": rratio if spec["recomb"] else 0.0},
           "seeds": seeds,
           "counts": {"num_trees": int(ts.num_trees),
                      "sites_mutated_total": int(mts.num_sites),
                      "dropped_multi_mutation": n_multi,
                      "dropped_non_segregating": n_nonseg,
                      "kept_biallelic_snv": n_kept,
                      "singletons": int(n_single),
                      "doubletons": int(n_double),
                      "achieved_density_per_bp": n_kept / L},
           "files": {"gt": {"path": gt_path, "sha256": sha256_file(gt_path),
                            "payload_sha256": gz_payload_sha256(gt_path),
                            "bytes": os.path.getsize(gt_path)},
                     "anc": {"path": anc_path,
                             "sha256": sha256_file(anc_path),
                             "payload_sha256": gz_payload_sha256(anc_path),
                             "bytes": os.path.getsize(anc_path)}},
           "timings_s": timings,
           "polarization_note": "anc == REF by construction (TRUE "
                                "polarization); mis-polarization, if "
                                "wanted, is applied downstream",
           "subscale": bool(n_dip != N_DIP or L != spec["L"])}
    path = write_json(os.path.join(out_records, f"arm_{arm}.json"), rec)
    print("arm record:", path)
    return rec


# ----------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------

def reference_checks(pins=None):
    """Analytic checks against the shipped reference pins (no msprime):
    analytic_selfchecks(); expected_branch_sfs for TK1 and for TK2 (fitted
    sizes) reproduces the recorded normalized spectra of pins['tk2_match']
    to 1e-12; TK1 at n = 6 reproduces selftest_expected.json 'exact_EL_gen'.
    Returns a dict of the worst relative deviations; raises AssertionError
    on failure."""
    pins = pins or load_pins()
    analytic_selfchecks()
    out = {}
    tk2 = [(TK2_TIMES[m], pins["tk2_fitted_sizes"][m])
           for m in range(len(TK2_TIMES))]
    worst = 0.0
    for n in TRUTH_NS:
        ach = pins["tk2_match"]["achieved"][f"n{n}"]
        for epochs, key in ((TK1_EPOCHS, "xi_target"), (tk2, "xi_fitted")):
            xi = expected_branch_sfs(epochs, n)
            xi = xi / xi.sum()
            rel = max(abs(a - b) / abs(b) for a, b in zip(xi.tolist(), ach[key]))
            worst = max(worst, rel)
    assert worst <= 1e-12, "tk2_match spectra not reproduced: %.3g" % worst
    out["tk2_match_rel"] = worst
    exp_path = os.path.join(REFERENCE_DIR, "selftest_expected.json")
    if os.path.exists(exp_path):
        with open(exp_path) as fh:
            zb = json.load(fh)["zcheck_exact_vs_mc"]
        exact = expected_branch_sfs(TK1_EPOCHS, zb["n"])
        rel = max(abs(a - b) / abs(b)
                  for a, b in zip(exact.tolist(), zb["exact_EL_gen"]))
        assert rel <= 1e-12, "exact E[L_i](TK1, n=6) not reproduced: %.3g" % rel
        out["exact_EL_rel"] = rel
    return out


def selftest(out_dir=None, pins=None):
    """Reduced-scale end-to-end (needs msprime): analytic checks -> seeded
    exact-vs-MC z-check (TK1, n = 6, 2000 replicates, |z| < 4 required) ->
    subscale class truth (TK1, 300 replicates) -> the gene-conversion arm at
    40 individuals over 300 kb (recombination + gene conversion + mutation +
    both writers), all on the 'selftest:' seed stream. Files go under out_dir
    (a temporary directory, removed afterwards, when None); returns the
    record, whose 'pass' is the verdict."""
    import msprime
    t_all = time.time()
    pins = pins or load_pins()
    analytic_selfchecks()
    tmp = None
    if out_dir is None:
        tmp = out_dir = tempfile.mkdtemp(prefix="planted_selftest_")
    os.makedirs(out_dir, exist_ok=True)
    try:
        # A) exact matcher vs msprime branch-mode MC (validates BOTH the
        # Tavare machinery and the ploidy/time conventions of the sim path)
        t0 = time.time()
        n, reps = 6, 2000
        seed = derive_seed(pins["base_seed"], "selftest:zcheck:TK1:n6")
        kwargs = sim_kwargs(class_params("TK1", pins))
        X = np.empty((reps, n - 1))
        it = msprime.sim_ancestry(samples=n, ploidy=1, num_replicates=reps,
                                  random_seed=seed, **kwargs)
        for idx, ts in enumerate(it):
            X[idx] = branch_class_lengths(ts.first(), n)
        exact = expected_branch_sfs(TK1_EPOCHS, n)
        mc_mean = X.mean(axis=0)
        mc_se = X.std(axis=0, ddof=1) / math.sqrt(reps)
        z = (mc_mean - exact) / mc_se
        zcheck = {"n": n, "reps": reps, "seed": seed,
                  "exact_EL_gen": exact.tolist(), "mc_EL_gen": mc_mean.tolist(),
                  "mc_se": mc_se.tolist(), "z": z.tolist(),
                  "max_abs_z": float(np.max(np.abs(z))),
                  "pass": bool(np.max(np.abs(z)) < 4.0),
                  "wall_s": round(time.time() - t0, 3)}

        # B) subscale class truth (TK1, reps=300 — marked subscale)
        truth_rec = run_truth("TK1", pins, 300, out_dir,
                              tag_prefix="selftest:")

        # C) smallest arm end-to-end (GC arm exercises recomb + gene
        # conversion + mutation + both output writers), reduced n and L
        arm_rec = run_arm("TK1_rON_GC", pins, out_dir, out_dir,
                          n_dip=40, L=300_000, tag_prefix="selftest:")

        rec = {"generator": "popcorn synth_truths", "kind": "selftest",
               "versions": versions(),
               "analytic_selfchecks": "PASS (constant-N 2N/i exact; n=2 Tavare)",
               "zcheck_exact_vs_mc": zcheck,
               "subscale_truth_TK1": {
                   "reps": 300, "seeds": {k: b["seed"] for k, b in
                                         truth_rec["mc"].items()},
                   "q_ij": {k: b["q_ij"] for k, b in truth_rec["mc"].items()},
                   "wall_s": truth_rec["wall_s_total"]},
               "subscale_arm_TK1_rON_GC": {
                   "n_dip": 40, "L": 300000,
                   "seeds": arm_rec["seeds"],
                   "counts": arm_rec["counts"],
                   "payload_sha256": {k: arm_rec["files"][k]["payload_sha256"]
                                      for k in ("gt", "anc")},
                   "timings_s": arm_rec["timings_s"]},
               "wall_s_total": round(time.time() - t_all, 3),
               "pass": bool(zcheck["pass"] and
                            arm_rec["counts"]["kept_biallelic_snv"] > 0),
               "out_dir": None if tmp else out_dir,
               "note": "selftest register: reduced reps/n/L, separate seed "
                       "stream ('selftest:' tags); production truths/arms "
                       "use the defaults (reps>=%d, n_dip=%d)"
                       % (TRUTH_REPS_DEFAULT, N_DIP)}
        if tmp is None:
            write_json(os.path.join(out_dir, "selftest.json"), rec)
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)
    return rec


# ----------------------------------------------------------------------------
# commands
# ----------------------------------------------------------------------------

def cmd_list(args):
    pins = load_pins(args.pins) if getattr(args, "pins", None) else None
    print("classes:", ", ".join(CLASSES))
    arms = dict(ARMS)
    if pins:
        extra = [c for c in pins.get("classes", {}) if c not in CLASSES]
        if extra:
            print("classes (pins):", ", ".join(extra))
        arms.update({a: s for a, s in pins.get("arms", {}).items()
                     if a not in ARMS})
    for a, s in arms.items():
        print("arm %-12s cls=%-4s recomb=%-5s gc=%-5s L=%d"
              % (a, s["cls"], s["recomb"], s["gc"], s["L"]))


def cmd_check(args):
    out = reference_checks(load_pins(args.pins) if args.pins else None)
    print("check PASS: constant-N 2N/i and n=2 Tavare identities; "
          + "; ".join("%s %.1e" % kv for kv in out.items()))


def cmd_pin(args):
    extra = {}
    if args.extra:
        with open(args.extra) as fh:
            extra = json.load(fh)
    make_pins(args.out, seed=args.seed, map_path=args.map, chrom=args.chrom,
              r_per_bp=args.r_per_bp, classes=extra.get("classes"),
              arms=extra.get("arms"))


def cmd_truth(args):
    pins = load_pins(args.pins)
    run_truth(args.cls, pins, args.reps, args.out)


def cmd_arm(args):
    pins = load_pins(args.pins)
    run_arm(args.arm, pins, args.out, args.out, n_dip=args.n_dip, L=args.length)


def cmd_all(args):
    pins = load_pins(args.pins)
    classes = list(CLASSES) + [c for c in pins.get("classes", {})
                               if c not in CLASSES]
    arms = list(ARMS) + [a for a in pins.get("arms", {}) if a not in ARMS]
    for cls in classes:
        run_truth(cls, pins, args.reps, args.out)
    for arm in arms:
        run_arm(arm, pins, args.out, args.out)


def cmd_selftest(args):
    rec = selftest(args.out, load_pins(args.pins) if args.pins else None)
    z = rec["zcheck_exact_vs_mc"]
    c = rec["subscale_arm_TK1_rON_GC"]["counts"]
    print("selftest %s: zcheck max|z| %.4f (< 4), subscale arm kept_snv %d, "
          "trees %d (%.1f s)%s"
          % ("PASS" if rec["pass"] else "FAIL", z["max_abs_z"],
             c["kept_biallelic_snv"], c["num_trees"], rec["wall_s_total"],
             "" if rec["out_dir"] is None else "; files in " + rec["out_dir"]))
    if not rec["pass"]:
        sys.exit(1)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="planted truths for SFS / two-site estimators: closed-form "
                    "(float) expected spectra and seeded msprime generators")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pins_help = "pins.json to run under (default: the shipped reference pins)"
    p = sub.add_parser("list", help="list truth classes and arms")
    p.add_argument("--pins", help=pins_help)
    p = sub.add_parser("check", help="analytic checks vs the reference (no msprime)")
    p.add_argument("--pins", help=pins_help)
    p = sub.add_parser("pin", help="write DIR/pins.json (base seed, map rate, fitted sizes)")
    p.add_argument("--out", required=True, metavar="DIR")
    p.add_argument("--seed", type=int, help="base seed (default: random)")
    p.add_argument("--map", help="bedGraph genetic map (chrom start end cM/Mb), plain or .gz")
    p.add_argument("--chrom", default="chr21", help="chromosome to average in --map")
    p.add_argument("--r-per-bp", type=float, dest="r_per_bp",
                   help="map-average crossover rate per bp per generation "
                        "(default without --map: %g)" % FALLBACK_R_MAP)
    p.add_argument("--extra", help="JSON file {'classes': {...}, 'arms': {...}} to add")
    p = sub.add_parser("truth", help="planted class truth -> DIR/truth_<CLS>.json")
    p.add_argument("--cls", required=True)
    p.add_argument("--out", required=True, metavar="DIR")
    p.add_argument("--pins", help=pins_help)
    p.add_argument("--reps", type=int, default=TRUTH_REPS_DEFAULT)
    p = sub.add_parser("arm", help="one arm's genotype + polarization tables -> DIR")
    p.add_argument("--arm", required=True)
    p.add_argument("--out", required=True, metavar="DIR")
    p.add_argument("--pins", help=pins_help)
    p.add_argument("--n-dip", type=int, dest="n_dip", default=None,
                   help="diploid individuals (default %d)" % N_DIP)
    p.add_argument("--length", type=float, default=None,
                   help="region length in bp (default: the arm's)")
    p = sub.add_parser("all", help="every class truth and every arm -> DIR")
    p.add_argument("--out", required=True, metavar="DIR")
    p.add_argument("--pins", help=pins_help)
    p.add_argument("--reps", type=int, default=TRUTH_REPS_DEFAULT)
    p = sub.add_parser("selftest", help="reduced-scale end to end (needs msprime)")
    p.add_argument("--out", default=None, metavar="DIR",
                   help="keep the generated files here (default: temporary)")
    p.add_argument("--pins", help=pins_help)
    args = ap.parse_args(argv)
    {"list": cmd_list, "check": cmd_check, "pin": cmd_pin, "truth": cmd_truth,
     "arm": cmd_arm, "all": cmd_all, "selftest": cmd_selftest}[args.cmd](args)


if __name__ == "__main__":
    main()

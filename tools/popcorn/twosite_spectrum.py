#!/usr/bin/env python3
"""twosite_spectrum — two-site frequency spectrum of a diploid genotype
panel at a fixed projected sample size, from genotype tables.

WHAT IT COMPUTES
  Given per-individual genotypes of biallelic SNVs along one chromosome
  (VCF, or a 'pos ref alt gts' TSV with gts a string of 0/1/2 alt dosages),
  an optional polarization table 'pos ref alt anc conf' (the format
  popcorn.ancestral / epo_join.py writes) and an optional genetic map
  (bedGraph 'chrom start end rate', rate = local cM/Mb), estimate() scans
  every pair of retained sites within max_bp of each other and accumulates
  the EXACT projected pair spectrum of twosite_projection.project_grid at
  each requested n = 2m (expectation over all C(N, m) subsets of m of the N
  individuals; every individual must be typed at every retained site):

    unfolded  phi[i-1][j-1] = sum over pairs of P(derived count i at
              site 1, j at site 2), i, j = 1..n-1, over pairs whose two
              sites both orient (ancestral base called; see orient_site);
    folded    the per-site minor-allele fold, (n/2) x (n/2), over all
              pairs (no polarization needed);

  binned by the distance between the two sites — genetic distance in cM
  from the map when one is given (cumulative cM by integrating the local
  rate over covered intervals; a site in an interval the map does not
  cover has no genetic position and its pairs are excluded and counted,
  not imputed), else physical distance in bp — with edges[b] <= d <
  edges[b+1]. Per bin it reports phi as exact fraction strings, the
  normalized spectrum q_hat = phi / sum(phi) (exact strings and floats),
  the number of pairs, and a delete-one-block jackknife standard error and
  95% interval for q_hat over blocks of the pair midpoint ((pos1+pos2)/2 in
  windows of block_mb megabases). Pair classes, tested in this order and
  mutually exclusive:
    mnv_control       bp distance <= mnv_bp (default 2): candidate
                      multi-nucleotide mutation pairs (Schrider et al.
                      2011), set aside from the signal and reported as
                      their own spectrum;
    na_gap_excluded   (map given) either site not covered by the map:
                      excluded, counted;
    hotspot_control   (map given) the maximum local rate over the map
                      intervals from site 1's to site 2's exceeds
                      hotspot_cm_per_mb: set aside, reported separately;
    out_of_bin_range  distance outside [edges[0], edges[-1]): counted;
    signal            binned and jackknifed.
  The two control classes are additionally binned by distance under the
  same edges and block rule (unbinnable control pairs are counted under
  bins_unbinned). The scan aggregates pairs by (class, bin, block, 3x3
  table, orientation flips) and adds each distinct group's projected
  matrices once, scaled by the exact integer multiplicity; since every
  accumulator is a sum of per-pair matrices that are pure functions of the
  table and Fraction arithmetic is exact, this is identical to per-pair
  accumulation, only faster.

  orient_site(pos, ref, alt, anc_map, high_only) decides the polarization
  of one site from the 'pos ref alt anc conf' table: +1 (REF ancestral,
  ALT derived), -1 (ALT ancestral: the site's axis is flipped, i -> n-i) or
  None with a reason (no_anc_row, refalt_mismatch, unpolarizable = no
  ancestral call, low_conf_skipped when high_only, anc_neither = the
  ancestral base is a third allele). The ancestral base is compared
  case-insensitively, so a low-confidence lowercase EPO call orients like
  its uppercase form and `conf` alone carries the confidence; pass
  high_only=True to use high-confidence calls only.

  brute_check_file(path) runs twosite_projection.brute_check on a small
  genotype file (DP operator vs explicit enumeration, exact ==).

REGISTER: DATA-UTILITY built on an EXACT core. The projected spectra phi
and q_hat are exact rationals (twosite_projection, fractions.Fraction end
to end; emitted as 'p/q' strings alongside floats); the jackknife standard
errors and intervals are plain double-precision summaries of the exact
block sums; distances in cM are floats from the map. Nothing here is a
fitted model or a certified bound. Standard library only.

INPUTS
  geno  VCF/VCF.gz (biallelic SNVs kept; a site with any missing GT is
        dropped so every retained site is typed in all N individuals; the
        dosage is the number of '1' alleles in GT) or TSV/TSV.gz with header
        'pos ref alt gts'. One chromosome per run: with a VCF pass chrom
        ('chrN' and 'N' both match) unless the file holds one chromosome.
  anc   'pos ref alt anc conf' TSV(.gz), conf in {high, low, none}; rows
        must match the genotype REF/ALT. Optional: without it only the
        folded spectra are filled.
  map   bedGraph(.gz) 'chrom start end rate' with rate in cM/Mb, intervals
        sorted and disjoint per chromosome, positions 0-based half-open;
        '#', 'track' and 'browser' lines are skipped. Optional: without it
        distances are in bp and the na_gap/hotspot classes never occur.
  n     even projected chromosome counts (m = n/2 individuals <= N).
  bins  ascending distance edges in the distance unit (cM with a map, bp
        without); defaults DEFAULT_BINS_CM / DEFAULT_BINS_BP.

OUTPUT (estimate() dict; the CLI writes it as JSON)
  inputs (paths as given, chrom, format, row counts), n_individuals,
  n_sites, site_counts, dist_unit ('cM'|'bp'), params,
  site_orientation_counts, pair_class_counts,
  signal_pairs_unfolded_eligible, mnv_control_pairs (positions, capped),
  aggregation (pairs grouped, distinct groups, cache counters), and
  results['n4'|'n6'|...][ 'unfolded' | 'folded' ] =
    {"bins": {"<b>": {edges, n_pairs, n_blocks, empty, phi_exact, qhat,
                      qhat_exact, jackknife{n_blocks_used, se, ci95_lo,
                      ci95_hi | note}}},
     "controls": {"<class>": {phi_exact, n_pairs, qhat, bins{...},
                              bins_unbinned{reason: count}}}}
  Matrices are row = site-1 count, column = site-2 count (site 1 is the
  smaller position); unfolded indices run 1..n-1, folded 1..n/2.

COST: the pair scan is O(pairs within max_bp) dictionary work; the exact
work is one project_grid per DISTINCT (3x3 table, n), so panels whose
variants are mostly rare (few distinct tables) run much faster than the
worst case. Measured worst case, random synthetic genotypes where almost
every pair has its own table, n = (4, 6), one core: about 3e4 pairs/s at
N = 12 and about 3e3 pairs/s at N = 50-100. Memory is the group dictionary
(one small tuple per distinct (class, bin, block, table, flips)) plus the
grid cache (peak RSS 1.6 GB in that test, at up to 3e5 distinct
(table, n) entries); cap it with cache_max_entries on large N (values
unchanged).

REFERENCES
  Hudson (2001), Genetics 159:1805-1817 (two-locus sampling distributions).
  Marth, Czabarka, Murvai & Sherry (2004), Genetics 166:351-372; Gutenkunst
    et al. (2009), PLoS Genet. 5:e1000695 (hypergeometric projection of
    frequency spectra).
  Schrider, Hourmozdi & Hahn (2011), Curr. Biol. 21:1051-1054 (pervasive
    multinucleotide mutational events; why pairs a few bp apart are set
    aside).
  Efron & Tibshirani (1993), An Introduction to the Bootstrap, Chapman &
    Hall (delete-one jackknife; here over genomic blocks).

USAGE
  from twosite_spectrum import estimate, write_json
  doc = estimate("calls.chr21.vcf.gz", anc="anc_chr21.tsv.gz",
                 gmap="genetic_map.bedGraph.gz", chrom="chr21", n=(4, 6))
  doc = estimate("panel.tsv", n=(8,), bins=[1, 100, 1000, 10000, 100000])
  write_json("spectrum.json", doc)
  # CLI
  python3 twosite_spectrum.py --geno calls.vcf.gz --anc anc.tsv.gz \\
          --map genetic_map.bedGraph.gz --chrom chr21 --n 4,6 --out spectrum.json
  python3 twosite_spectrum.py --geno panel.tsv --n 8 --bins 1,100,1000,10000,100000 --out s.json
  python3 twosite_spectrum.py --brute-check --geno tiny_geno.tsv --n 4,6,8
Standard library only.
"""
import argparse
import gzip
import json
import math
import os
import re
import sys

from twosite_projection import (GridCache, add_into, add_scaled,
                                brute_check, cells_key, nine_cell, zeros)

BASES = frozenset("ACGT")
DEFAULT_BINS_CM = "0,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1"
DEFAULT_BINS_BP = "1,10,100,1000,10000,100000"
PAIR_CLASSES = ("mnv_control", "na_gap_excluded", "hotspot_control",
                "out_of_bin_range", "signal")

__all__ = ["BASES", "DEFAULT_BINS_CM", "DEFAULT_BINS_BP", "PAIR_CLASSES",
           "load_geno_vcf", "load_geno_tsv", "load_geno", "load_anc",
           "load_map", "annotate_cm", "span_max_rate", "orient_site",
           "summarize_bin", "estimate", "brute_check_file", "write_json",
           "main"]


# -------------------------------------------------------------------- input
def _opener(path):
    path = os.fspath(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def _chrom_forms(chrom):
    """The accepted spellings of a chromosome name: with and without 'chr'."""
    c = str(chrom)
    bare = c[3:] if c.lower().startswith("chr") else c
    return {c, bare, "chr" + bare}


def load_geno_vcf(path, chrom):
    """Biallelic SNV genotypes from a VCF(.gz).

    Returns (sites, n_ind, counts) with sites = [(pos, ref, alt,
    dosages)], dosages a tuple of per-individual alt-allele counts read
    from the GT field ('/' or '|' separated; count of '1' alleles). A
    record with ',' in ALT, a non-ACGT or multi-base REF/ALT, or ANY
    missing GT ('.') is skipped and counted; with chrom given ('chrN' and
    'N' both match) records of other chromosomes are skipped and counted."""
    keep = _chrom_forms(chrom) if chrom is not None else None
    sites, counts = [], {"nonbiallelic_or_multiallelic": 0, "non_snv": 0,
                         "missing_gt": 0, "other_chrom": 0, "retained": 0}
    n_ind = None
    with _opener(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                n_ind = len(line.rstrip("\n").split("\t")) - 9
                continue
            t = line.rstrip("\n").split("\t")
            if keep is not None and t[0] not in keep:
                counts["other_chrom"] += 1
                continue
            ref, alt = t[3], t[4]
            if "," in alt:
                counts["nonbiallelic_or_multiallelic"] += 1
                continue
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES \
                    or alt not in BASES:
                counts["non_snv"] += 1
                continue
            gt_i = t[8].split(":").index("GT")
            dosages, ok = [], True
            for s in t[9:]:
                alleles = re.split(r"[/|]", s.split(":")[gt_i])
                if "." in alleles:
                    ok = False
                    break
                dosages.append(sum(a == "1" for a in alleles))
            if not ok:
                counts["missing_gt"] += 1
                continue
            sites.append((int(t[1]), ref, alt, tuple(dosages)))
            counts["retained"] += 1
    if n_ind is None and sites:
        n_ind = len(sites[0][3])
    return sites, n_ind, counts


def load_geno_tsv(path):
    """Genotypes from a TSV(.gz) with header 'pos ref alt gts', gts a string
    of 0/1/2 alt dosages, one character per individual. Non-ACGT REF/ALT and
    rows with any other dosage character are skipped and counted.
    Returns (sites, n_ind, counts)."""
    sites = []
    counts = {"non_snv": 0, "bad_dosage": 0, "retained": 0}
    with _opener(path) as f:
        header = f.readline().split()
        if header[:4] != ["pos", "ref", "alt", "gts"]:
            raise ValueError(f"TSV geno header must be 'pos ref alt gts', "
                             f"got {header}")
        for line in f:
            pos, ref, alt, gts = line.split()
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES \
                    or alt not in BASES:
                counts["non_snv"] += 1
                continue
            if not set(gts) <= set("012"):
                counts["bad_dosage"] += 1
                continue
            sites.append((int(pos), ref, alt, tuple(int(c) for c in gts)))
            counts["retained"] += 1
    n_ind = len(sites[0][3]) if sites else 0
    return sites, n_ind, counts


def load_geno(path, fmt="auto", chrom=None):
    """Dispatch to the VCF or TSV reader ('auto': VCF iff '.vcf' in the file
    name), sort sites by position and check a constant individual count.
    Returns (sites, n_ind, counts, fmt_used)."""
    path = os.fspath(path)
    if fmt == "auto":
        fmt = "vcf" if ".vcf" in os.path.basename(path) else "tsv"
    if fmt == "vcf":
        sites, n_ind, counts = load_geno_vcf(path, chrom)
    else:
        sites, n_ind, counts = load_geno_tsv(path)
    sites.sort(key=lambda s: s[0])
    for s in sites:
        if len(s[3]) != n_ind:
            raise ValueError("inconsistent individual count across sites")
    return sites, n_ind, counts, fmt


def load_anc(path):
    """Polarization table 'pos ref alt anc conf' (TSV, .gz ok; '#' lines
    skipped) -> {pos: (ref, alt, anc, conf)}."""
    anc = {}
    with _opener(path) as f:
        header = f.readline().split()
        if header != ["pos", "ref", "alt", "anc", "conf"]:
            raise ValueError(f"anc TSV header must be "
                             f"'pos ref alt anc conf', got {header}")
        for line in f:
            if line.startswith("#"):
                continue
            pos, ref, alt, a, conf = line.split()
            anc[int(pos)] = (ref, alt, a, conf)
    return anc


# ---------------------------------------------------------------------- map
def load_map(path, chrom):
    """Genetic map in bedGraph form (chrom start end rate, rate = local
    cM/Mb; '#', 'track', 'browser' lines skipped; 'chrN' and 'N' both match
    chrom) -> sorted intervals + cumulative cM at interval starts.

    Cumulative cM by integration: cM(x) = sum over covered intervals below
    x of (end-start)*rate/1e6; uncovered gaps contribute nothing, and a
    SITE inside a gap has no genetic position (annotate_cm gives None: the
    site is excluded, not imputed). Intervals must be sorted and disjoint."""
    keep = _chrom_forms(chrom)
    starts, ends, rates = [], [], []
    with _opener(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e, r = line.split()[:4]
            if c not in keep:
                continue
            starts.append(int(s))
            ends.append(int(e))
            rates.append(float(r))
    cum, acc = [], 0.0
    for s, e, r in zip(starts, ends, rates):
        cum.append(acc)
        acc += (e - s) * r / 1e6
    for k in range(1, len(starts)):
        if starts[k] < ends[k - 1]:
            raise ValueError("map intervals not sorted/disjoint")
    return {"starts": starts, "ends": ends, "rates": rates, "cum": cum,
            "n_intervals": len(starts)}


def annotate_cm(sites, m):
    """Sorted-merge local lookup: one pass over sorted sites x sorted
    intervals. Returns per-site (cM, interval_idx) with None for NA gaps."""
    out, k, K = [], 0, m["n_intervals"]
    for pos, _, _, _ in sites:
        while k < K and m["ends"][k] <= pos:
            k += 1
        if k < K and m["starts"][k] <= pos < m["ends"][k]:
            out.append((m["cum"][k]
                        + m["rates"][k] * (pos - m["starts"][k]) / 1e6, k))
        else:
            out.append((None, None))
    return out


def span_max_rate(m, k1, k2):
    """Maximum local rate over map intervals k1..k2 inclusive."""
    return max(m["rates"][k1:k2 + 1])


# ------------------------------------------------------------------ orient
def orient_site(pos, ref, alt, anc_map, high_only):
    """Polarization of one site from the 'pos ref alt anc conf' table.

    Returns (+1, 'ok') when the ancestral base equals REF (ALT derived),
    (-1, 'ok') when it equals ALT (REF derived; flip this site's axis), or
    (None, reason) with reason in no_anc_row | refalt_mismatch |
    unpolarizable (conf 'none' or anc not a base) | low_conf_skipped
    (high_only and conf 'low') | anc_neither (third allele). The ancestral
    base is upper-cased for the comparison: a lowercase (low-confidence)
    EPO call orients like its uppercase form; `conf` carries the
    confidence and high_only filters on it."""
    rec = anc_map.get(pos)
    if rec is None:
        return None, "no_anc_row"
    aref, aalt, anc, conf = rec
    if (aref, aalt) != (ref, alt):
        return None, "refalt_mismatch"
    anc = anc.upper()
    if conf == "none" or anc not in BASES:
        return None, "unpolarizable"
    if high_only and conf == "low":
        return None, "low_conf_skipped"
    if anc == ref:
        return +1, "ok"
    if anc == alt:
        return -1, "ok"
    return None, "anc_neither"


# --------------------------------------------------------------- jackknife
def summarize_bin(blocks, dim):
    """blocks: {block_id: (matrix of Fractions, n_pairs)}. Exact totals phi,
    q^hat = phi / sum(phi) (exact strings + floats), and a delete-one-block
    jackknife on q^hat: se^2 = (G-1)/G * sum_g (q_(-g) - mean)^2 over the G
    blocks whose removal leaves a nonzero total; ci95 = q^hat +- 1.96 se.
    Fewer than 2 usable blocks -> no CI (a note instead)."""
    tot = zeros(dim)
    n_pairs = 0
    for mat, np_ in blocks.values():
        add_into(tot, mat)
        n_pairs += np_
    S = sum(sum(row) for row in tot)
    if S == 0:
        return {"n_pairs": n_pairs, "n_blocks": len(blocks), "empty": True}
    qhat = [[float(tot[r][c] / S) for c in range(dim)] for r in range(dim)]
    qhat_exact = [[str(tot[r][c] / S) for c in range(dim)]
                  for r in range(dim)]
    phi_exact = [[str(tot[r][c]) for c in range(dim)] for r in range(dim)]
    out = {"n_pairs": n_pairs, "n_blocks": len(blocks), "empty": False,
           "phi_exact": phi_exact, "qhat": qhat, "qhat_exact": qhat_exact}
    reps = []
    for b, (mat, _) in sorted(blocks.items()):
        Sg = sum(sum(row) for row in mat)
        if S - Sg == 0:
            continue
        reps.append([[float((tot[r][c] - mat[r][c]) / (S - Sg))
                      for c in range(dim)] for r in range(dim)])
    G = len(reps)
    if G >= 2:
        se = [[0.0] * dim for _ in range(dim)]
        for r in range(dim):
            for c in range(dim):
                vals = [rep[r][c] for rep in reps]
                mu = sum(vals) / G
                var = (G - 1) / G * sum((v - mu) ** 2 for v in vals)
                se[r][c] = math.sqrt(var)
        out["jackknife"] = {
            "n_blocks_used": G,
            "se": se,
            "ci95_lo": [[qhat[r][c] - 1.96 * se[r][c] for c in range(dim)]
                        for r in range(dim)],
            "ci95_hi": [[qhat[r][c] + 1.96 * se[r][c] for c in range(dim)]
                        for r in range(dim)],
        }
    else:
        out["jackknife"] = {"n_blocks_used": G,
                            "note": "fewer than 2 nonempty blocks; no CI"}
    return out


# ----------------------------------------------------------------- estimate
def _edges(bins, unit):
    if bins is None:
        bins = DEFAULT_BINS_CM if unit == "cM" else DEFAULT_BINS_BP
    if isinstance(bins, str):
        edges = [float(x) for x in bins.split(",")]
    else:
        edges = [float(x) for x in bins]
    if edges != sorted(edges) or len(edges) < 2:
        raise ValueError("bin edges must be ascending, >=2 values")
    return edges


def estimate(geno, anc=None, gmap=None, chrom=None, n=(4, 6), bins=None,
             max_bp=100_000, mnv_bp=2, hotspot_cm_per_mb=10.0, block_mb=1.0,
             high_conf_only=False, fmt="auto", mnv_pairs_cap=200_000,
             cache_max_entries=0):
    """Two-site spectra of a genotype panel at projected sizes n (see the
    module docstring for the classes, binning, registers and output).

    geno   path (VCF/TSV, see load_geno) or a pre-loaded
           (sites, n_ind[, counts]) tuple with sites = [(pos, ref, alt,
           dosages)].
    anc    path of a 'pos ref alt anc conf' table, a {pos: (ref, alt, anc,
           conf)} dict, or None (folded spectra only).
    gmap   path of a bedGraph genetic map, a load_map() dict, or None
           (distances in bp).
    chrom  chromosome for the VCF and map filters ('chrN'/'N' both match).
    n      iterable of even projected chromosome counts, or 'a,b' string.
    bins   distance edges (list or 'a,b,...' string; cM with a map, bp
           without); None -> DEFAULT_BINS_CM / DEFAULT_BINS_BP.
    Returns the output dict (JSON-serializable)."""
    if isinstance(geno, (str, os.PathLike)):
        geno_label = os.fspath(geno)
        sites, n_ind, geno_counts, fmt_used = load_geno(geno, fmt, chrom)
    else:
        sites, n_ind = list(geno[0]), int(geno[1])
        geno_counts = dict(geno[2]) if len(geno) > 2 else {"retained":
                                                            len(sites)}
        fmt_used, geno_label = "sites", None
        sites.sort(key=lambda s: s[0])
        for s in sites:
            if len(s[3]) != n_ind:
                raise ValueError("inconsistent individual count across sites")
    if anc is None:
        anc_map, anc_label = {}, None
    elif isinstance(anc, dict):
        anc_map, anc_label = anc, None
    else:
        anc_map, anc_label = load_anc(anc), os.fspath(anc)
    use_map = gmap is not None
    if use_map:
        if isinstance(gmap, dict):
            m_map, map_label = gmap, None
        else:
            m_map, map_label = load_map(gmap, chrom), os.fspath(gmap)
        cm_ann = annotate_cm(sites, m_map)
        unit = "cM"
    else:
        m_map, map_label, cm_ann, unit = None, None, None, "bp"
    if isinstance(n, str):
        n_list = sorted(int(x) for x in n.split(","))
    else:
        n_list = sorted(int(x) for x in n)
    for nn in n_list:
        if nn % 2 or nn < 2:
            raise ValueError(f"n must be even (n=2m individuals-projected), "
                             f"got {nn}")
        if nn // 2 > n_ind:
            raise ValueError(f"n={nn} needs m={nn//2} <= N={n_ind} "
                             f"individuals")
    edges = _edges(bins, unit)
    orient = []
    orient_counts = {}
    for pos, ref, alt, _ in sites:
        o, why = orient_site(pos, ref, alt, anc_map, high_conf_only)
        orient.append(o)
        orient_counts[why] = orient_counts.get(why, 0) + 1

    block_bp = int(block_mb * 1e6)
    nbins = len(edges) - 1
    # acc[n][register]: 'bins'->{bin:{block:[mat,n_pairs]}}, 'ctl'->{cls:mat}
    acc = {nn: {"unfolded": {"bins": {}, "ctl": {}},
                "folded": {"bins": {}, "ctl": {}}} for nn in n_list}
    cls_counts = {c: 0 for c in PAIR_CLASSES}
    ctl_unfolded_pairs = {nn: {} for nn in n_list}
    # per-bin control-class accumulators (same binning rule as signal pairs)
    ctlbins = {nn: {"unfolded": {}, "folded": {}} for nn in n_list}
    ctl_unbinned = {"unfolded": {}, "folded": {}}  # reg -> cls -> why -> ct
    signal_unfolded_pairs = 0
    mnv_pairs = []
    # Aggregation: insertion-ordered multiplicity dict
    #   groups[(cls, bslot, block, ckey, flips)] = exact int pair count
    # bslot = signal bin (cls == "signal") | control bin or None (controls)
    # flips = (row_flip, col_flip) when both sites orient, else None.
    # ckey tuples are interned so equal tables share one tuple object
    # (memory only; equality/hashing unchanged). n is NOT in the key: every
    # n sees the same pair set, so one multiplicity serves all n.
    groups = {}
    _key_intern = {}
    S = len(sites)
    for i in range(S):
        pos1, _, _, d1 = sites[i]
        for j in range(i + 1, S):
            pos2, _, _, d2 = sites[j]
            d_bp = pos2 - pos1
            if d_bp > max_bp:
                break
            if use_map:
                cm1, k1 = cm_ann[i]
                cm2, k2 = cm_ann[j]
            if d_bp <= mnv_bp:
                cls, bi = "mnv_control", None
                if len(mnv_pairs) < mnv_pairs_cap:
                    mnv_pairs.append([pos1, pos2])
            elif use_map and (cm1 is None or cm2 is None):
                cls_counts["na_gap_excluded"] += 1
                continue  # gap sites EXCLUDED, not imputed
            elif use_map and span_max_rate(m_map, k1, k2) > hotspot_cm_per_mb:
                cls, bi = "hotspot_control", None
            else:
                d = (cm2 - cm1) if use_map else d_bp
                if d < edges[0] or d >= edges[-1]:
                    cls_counts["out_of_bin_range"] += 1
                    continue
                bi = max(b for b in range(nbins) if edges[b] <= d)
                cls = "signal"
            cls_counts[cls] += 1
            both_orient = orient[i] is not None and orient[j] is not None
            if cls == "signal" and both_orient:
                signal_unfolded_pairs += 1
            ckey = cells_key(nine_cell(d1, d2))
            ckey = _key_intern.setdefault(ckey, ckey)
            block = ((pos1 + pos2) // 2) // block_bp
            if cls != "signal":
                # bin the control pair per the IDENTICAL signal rule (same
                # edges, same block rule); unbinnable pairs counted, never
                # binned. Class totals are accumulated regardless.
                if use_map and (cm1 is None or cm2 is None):
                    bi_c, why_c = None, "na_gap"
                else:
                    d_c = (cm2 - cm1) if use_map else d_bp
                    if d_c < edges[0] or d_c >= edges[-1]:
                        bi_c, why_c = None, "out_of_bin_range"
                    else:
                        bi_c = max(b for b in range(nbins)
                                   if edges[b] <= d_c)
                if bi_c is None:
                    dct = ctl_unbinned["folded"].setdefault(cls, {})
                    dct[why_c] = dct.get(why_c, 0) + 1
                    if both_orient:
                        dct = ctl_unbinned["unfolded"].setdefault(cls, {})
                        dct[why_c] = dct.get(why_c, 0) + 1
                bslot = bi_c
            else:
                bslot = bi
            flips = ((orient[i] < 0, orient[j] < 0)
                     if both_orient else None)
            gk = (cls, bslot, block, ckey, flips)
            groups[gk] = groups.get(gk, 0) + 1

    # Replay: one pass over distinct groups in insertion order (= first-
    # occurrence order of the pair scan). Per group and per n this performs
    # exactly the additions a per-pair accumulation would perform for each
    # of its `mult` pairs; add_scaled multiplies by the exact multiplicity.
    # `cent` is the READ-ONLY cache slot, `aent` the mutable accumulator
    # entry: never pass a cached object as the first argument of
    # add_scaled/add_into.
    cache = GridCache(cache_max_entries)
    for (cls, bslot, block, ckey, flips), mult in groups.items():
        both_orient = flips is not None
        for nn in n_list:
            cent = cache.entry(ckey, n_ind, nn // 2)
            fmat = cache.fold(cent)
            if both_orient:
                umat = cache.umat(cent, flips[0], flips[1])
            if cls == "signal":
                dim = nn - 1
                blocks = acc[nn]["folded"]["bins"].setdefault(bslot, {})
                aent = blocks.setdefault(block, [zeros(nn // 2), 0])
                add_scaled(aent[0], fmat, mult)
                aent[1] += mult
                if both_orient:
                    blocks = acc[nn]["unfolded"]["bins"].setdefault(bslot, {})
                    aent = blocks.setdefault(block, [zeros(dim), 0])
                    add_scaled(aent[0], umat, mult)
                    aent[1] += mult
            else:
                ctl = acc[nn]["folded"]["ctl"].setdefault(
                    cls, zeros(nn // 2))
                add_scaled(ctl, fmat, mult)
                if bslot is not None:  # binned control (folded)
                    blocks = ctlbins[nn]["folded"].setdefault(
                        cls, {}).setdefault(bslot, {})
                    aent = blocks.setdefault(block, [zeros(nn // 2), 0])
                    add_scaled(aent[0], fmat, mult)
                    aent[1] += mult
                if both_orient:
                    ctl = acc[nn]["unfolded"]["ctl"].setdefault(
                        cls, zeros(nn - 1))
                    add_scaled(ctl, umat, mult)
                    ctl_unfolded_pairs[nn][cls] = \
                        ctl_unfolded_pairs[nn].get(cls, 0) + mult
                    if bslot is not None:  # binned control (unfolded)
                        blocks = ctlbins[nn]["unfolded"].setdefault(
                            cls, {}).setdefault(bslot, {})
                        aent = blocks.setdefault(block, [zeros(nn - 1), 0])
                        add_scaled(aent[0], umat, mult)
                        aent[1] += mult

    def ctl_report(mat, npairs=None):
        Ssum = sum(sum(r) for r in mat)
        rep = {"phi_exact": [[str(v) for v in row] for row in mat]}
        if npairs is not None:
            rep["n_pairs"] = npairs
        if Ssum:
            rep["qhat"] = [[float(v / Ssum) for v in row] for row in mat]
        return rep

    results = {}
    for nn in n_list:
        rn = {}
        for reg in ("unfolded", "folded"):
            dim = (nn - 1) if reg == "unfolded" else nn // 2
            bins_out = {}
            for bi, blocks in sorted(acc[nn][reg]["bins"].items()):
                bd = {b: (m_[0], m_[1]) for b, m_ in blocks.items()}
                bins_out[str(bi)] = {
                    "edges": [edges[bi], edges[bi + 1]],
                    **summarize_bin(bd, dim)}
            ctl_out = {cls: ctl_report(mat,
                       ctl_unfolded_pairs[nn].get(cls)
                       if reg == "unfolded" else cls_counts[cls])
                       for cls, mat in acc[nn][reg]["ctl"].items()}
            for cls, rep in ctl_out.items():  # per-bin control spectra
                bl = ctlbins[nn][reg].get(cls, {})
                rep["bins"] = {str(bi): {
                    "edges": [edges[bi], edges[bi + 1]],
                    **summarize_bin({b: (m_[0], m_[1])
                                     for b, m_ in blocks.items()}, dim)}
                    for bi, blocks in sorted(bl.items())}
                rep["bins_unbinned"] = ctl_unbinned[reg].get(cls, {})
            rn[reg] = {"bins": bins_out, "controls": ctl_out}
        results[f"n{nn}"] = rn

    agg = {"pairs_grouped": sum(groups.values()),
           "distinct_groups": len(groups), "cache": cache.report()}
    return {
        "step": "twosite_spectrum",
        "inputs": {"geno": geno_label, "anc": anc_label, "map": map_label,
                   "chrom": chrom, "format": fmt_used,
                   "n_anc_rows": len(anc_map),
                   "n_map_intervals": (m_map["n_intervals"]
                                       if use_map else None)},
        "n_individuals": n_ind,
        "n_sites": len(sites),
        "site_counts": geno_counts,
        "dist_unit": unit,
        "params": {
            "n_projected": n_list, "max_bp": max_bp, "mnv_bp": mnv_bp,
            "hotspot_cm_per_mb": hotspot_cm_per_mb if use_map else None,
            "block_mb": block_mb, "bin_edges": edges,
            "high_conf_only": high_conf_only,
        },
        "site_orientation_counts": orient_counts,
        "pair_class_counts": cls_counts,
        "mnv_control_pairs": mnv_pairs,
        "mnv_control_pairs_cap": mnv_pairs_cap,
        "signal_pairs_unfolded_eligible": signal_unfolded_pairs,
        "projection": "exact rational (fractions), expectation over all "
                      "C(N,m) subsamples of individuals; phi sums exact end "
                      "to end; qhat also emitted as exact fraction strings; "
                      "jackknife se/ci are floats",
        "aggregation": agg,
        "results": results,
    }


def brute_check_file(path, n_list=(4, 6, 8), fmt="auto", chrom=None):
    """Read a small genotype file and run twosite_projection.brute_check on
    its dosage vectors (DP operator vs explicit enumeration, exact ==)."""
    sites, n_ind, _counts, fmt_used = load_geno(path, fmt, chrom)
    r = brute_check([s[3] for s in sites], n_ind, n_list)
    r["geno"] = os.fspath(path)
    r["format"] = fmt_used
    return r


def write_json(path, obj):
    """JSON with indent 1 and a trailing newline; '-' writes to stdout."""
    if path == "-":
        json.dump(obj, sys.stdout, indent=1)
        sys.stdout.write("\n")
        return
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")


# ---------------------------------------------------------------------- CLI
def main(argv=None):
    p = argparse.ArgumentParser(
        description="Two-site frequency spectrum of a diploid genotype "
                    "panel at projected sample sizes n = 2m (exact "
                    "hypergeometric projection over individuals), binned "
                    "by genetic (cM, with --map) or physical (bp) distance, "
                    "with block-jackknife intervals.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--geno", required=True,
                   help="genotype VCF(.gz) or TSV(.gz) ('pos ref alt gts')")
    p.add_argument("--anc", default=None,
                   help="polarization TSV 'pos ref alt anc conf' "
                        "(without it: folded spectra only)")
    p.add_argument("--map", default=None,
                   help="genetic map bedGraph(.gz) 'chrom start end cM/Mb' "
                        "(without it: distances in bp)")
    p.add_argument("--chrom", default=None,
                   help="chromosome (VCF + map filter; 'chrN'/'N' match)")
    p.add_argument("--n", default="4,6",
                   help="projected chromosome counts n=2m, comma-separated")
    p.add_argument("--max-bp", type=int, default=100_000,
                   help="pair retention: bp distance <= this")
    p.add_argument("--mnv-bp", type=int, default=2,
                   help="pairs at bp distance <= this go to mnv_control")
    p.add_argument("--hotspot-cm-per-mb", type=float, default=10.0,
                   help="with --map: pairs spanning a local rate above "
                        "this go to hotspot_control")
    p.add_argument("--bins", default=None,
                   help="distance bin edges, comma-separated (cM with "
                        "--map, default " + DEFAULT_BINS_CM + "; bp "
                        "without, default " + DEFAULT_BINS_BP + ")")
    p.add_argument("--block-mb", type=float, default=1.0,
                   help="jackknife block size in Mb (pair midpoint)")
    p.add_argument("--high-conf-only", action="store_true",
                   help="orient only sites whose conf is 'high'")
    p.add_argument("--format", choices=("auto", "vcf", "tsv"),
                   default="auto")
    p.add_argument("--out", default="-",
                   help="output JSON path ('-' = stdout)")
    p.add_argument("--cache-max-entries", type=int, default=0,
                   help="cap on memoized distinct tables (0 = unbounded; "
                        "past the cap new tables are recomputed, never "
                        "evicted; values are identical for any cap)")
    p.add_argument("--brute-check", action="store_true",
                   help="instead of estimating: check the projection "
                        "operator against explicit enumeration of ALL "
                        "subsamples on --geno (<=30 individuals, <=200 "
                        "sites); exact equality required, rc=1 otherwise")
    args = p.parse_args(argv)
    n_list = [int(x) for x in args.n.split(",")]
    log = sys.stderr if args.out == "-" else sys.stdout
    try:
        if args.brute_check:
            r = brute_check_file(args.geno, n_list, args.format, args.chrom)
            for k, v in r["per_n"].items():
                print(f"{k}: m={v['m']} subsamples={v['n_subsamples']} "
                      f"pairs={v['n_pairs']} cells={v['cells_checked']} "
                      f"mismatches={v['mismatches']} "
                      f"sum_to_one_failures={v['sum_to_one_failures']}",
                      file=log)
            print(f"[brute-check] EXACT_EQUAL_ALL = {r['EXACT_EQUAL_ALL']}",
                  file=log, flush=True)
            write_json(args.out, r)
            return 0 if r["EXACT_EQUAL_ALL"] else 1
        doc = estimate(args.geno, anc=args.anc, gmap=args.map,
                       chrom=args.chrom, n=n_list, bins=args.bins,
                       max_bp=args.max_bp, mnv_bp=args.mnv_bp,
                       hotspot_cm_per_mb=args.hotspot_cm_per_mb,
                       block_mb=args.block_mb,
                       high_conf_only=args.high_conf_only, fmt=args.format,
                       cache_max_entries=args.cache_max_entries)
    except ValueError as e:
        p.error(str(e))
    write_json(args.out, doc)
    print(f"[twosite_spectrum] N={doc['n_individuals']} "
          f"sites={doc['n_sites']} unit={doc['dist_unit']} "
          f"pairs: {doc['pair_class_counts']} "
          f"groups: {doc['aggregation']['distinct_groups']}"
          + (f" -> {args.out}" if args.out != "-" else ""),
          file=log, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

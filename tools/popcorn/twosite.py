"""popcorn.twosite — exact hypergeometric projection of two-site genotype
tables and the two-site frequency spectrum of a diploid panel.

Aliases, in place and by identity (see popcorn._loader):
  twosite_projection — the EXACT core: nine_cell (3x3 genotype table of a
                       pair of sites from two 0/1/2 dosage vectors),
                       project_grid(cells, N, m) = P(alt counts (i, j) in a
                       uniformly random subset of m of the N individuals),
                       i, j = 0..2m, as Fractions summing to exactly 1
                       (integer DP over the nine cells, one exact
                       division); brute_grid (explicit enumeration, the
                       reference route); flip_grid / seg_matrix /
                       fold_matrix (orientation, unfolded segregating
                       block, minor-allele fold); add_into / add_scaled /
                       zeros / cells_key; GridCache (memo by distinct
                       table); brute_check (DP vs enumeration, exact ==).
  twosite_spectrum   — the DATA-UTILITY estimator: read genotypes (VCF or
                       'pos ref alt gts' TSV), an optional 'pos ref alt anc
                       conf' polarization table (the table popcorn.ancestral
                       writes) and an optional bedGraph genetic map; scan
                       site pairs within max_bp; classify (mnv_control /
                       na_gap_excluded / hotspot_control / out_of_bin_range
                       / signal); accumulate the exact projected pair
                       spectra per distance bin (cM with a map, bp without),
                       unfolded over oriented pairs and folded over all;
                       block-jackknife intervals; CLI with --brute-check.

Registers: twosite_projection EXACT (fractions.Fraction end to end);
twosite_spectrum DATA-UTILITY on that exact core (phi and q_hat exact, the
jackknife se/ci plain floats, cM distances floats from the map). Standard
library only. Reference data: reference/twosite/ (a synthetic 12-individual,
24-site panel with its polarization table and genetic map, the expected
brute-check counts and the expected spectra).
"""
from ._loader import import_in_place

twosite_projection = import_in_place("twosite_projection")
twosite_spectrum = import_in_place("twosite_spectrum")

# Convenience re-exports: exact core.
nine_cell = twosite_projection.nine_cell
project_grid = twosite_projection.project_grid
brute_grid = twosite_projection.brute_grid
flip_grid = twosite_projection.flip_grid
seg_matrix = twosite_projection.seg_matrix
fold_matrix = twosite_projection.fold_matrix
add_into = twosite_projection.add_into
add_scaled = twosite_projection.add_scaled
zeros = twosite_projection.zeros
cells_key = twosite_projection.cells_key
GridCache = twosite_projection.GridCache
brute_check = twosite_projection.brute_check

# Convenience re-exports: estimator.
estimate = twosite_spectrum.estimate
orient_site = twosite_spectrum.orient_site
summarize_bin = twosite_spectrum.summarize_bin
load_geno = twosite_spectrum.load_geno
load_geno_vcf = twosite_spectrum.load_geno_vcf
load_geno_tsv = twosite_spectrum.load_geno_tsv
load_anc = twosite_spectrum.load_anc
load_map = twosite_spectrum.load_map
annotate_cm = twosite_spectrum.annotate_cm
brute_check_file = twosite_spectrum.brute_check_file
write_json = twosite_spectrum.write_json
PAIR_CLASSES = twosite_spectrum.PAIR_CLASSES

__all__ = ["twosite_projection", "twosite_spectrum", "nine_cell",
           "project_grid", "brute_grid", "flip_grid", "seg_matrix",
           "fold_matrix", "add_into", "add_scaled", "zeros", "cells_key",
           "GridCache", "brute_check", "estimate", "orient_site",
           "summarize_bin", "load_geno", "load_geno_vcf", "load_geno_tsv",
           "load_anc", "load_map", "annotate_cm", "brute_check_file",
           "write_json", "PAIR_CLASSES"]

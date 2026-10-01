"""popcorn.ancestral — ancestral-state join against an EPO ancestral FASTA.

Aliases, in place and by identity (see popcorn._loader):
  epo_join — stream a sites table (chrom pos ref alt / pos ref alt) against
             a per-chromosome Ensembl-EPO-style ancestral-allele FASTA;
             classify each site by EPO confidence (uppercase = high,
             lowercase = low, 'N'/'-'/'.' = no call), orient REF/ALT to
             ancestral/derived over biallelic SNVs, write the polarized
             'pos ref alt anc conf' table and return the counters. Standard
             library only; also a CLI (python3 epo_join.py --fasta --sites
             --out [--chrom] [--summary]).

Register: DATA-UTILITY (integer counters and string columns; no estimator,
nothing certified). Reference fixture: reference/ancestral/ (a synthetic
12-bp ancestor and ten sites exercising every branch).
"""
from ._loader import import_in_place

epo_join = import_in_place("epo_join")

# Convenience re-exports.
join_sites = epo_join.join_sites
polarization_rates = epo_join.polarization_rates
load_fasta = epo_join.load_fasta
ancestral_base = epo_join.ancestral_base
classify = epo_join.classify
orient = epo_join.orient
iter_sites = epo_join.iter_sites
HIGH = epo_join.HIGH
LOW = epo_join.LOW

__all__ = ["epo_join", "join_sites", "polarization_rates", "load_fasta",
           "ancestral_base", "classify", "orient", "iter_sites", "HIGH",
           "LOW"]

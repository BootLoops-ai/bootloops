"""popcorn.planted — planted truths for SFS / two-site estimators.

Aliases, in place and by identity (see popcorn._loader):
  synth_truths — (i) closed-form (float) expected branch-length spectrum
                 E[L_i] under a piecewise-constant population size (Tavare's
                 ancestral process, double precision; ploidy-1 units:
                 epochs [(t_start_gen, N_chrom), ...], pair rate 1/N per
                 generation) and the epoch-size fit that builds an
                 SFS-matched pair of histories; (ii) seeded msprime
                 generators for named truth classes (TK1/TK2 multi-epoch
                 Kingman, B13/B17 Beta-coalescent, D005/D020 Dirac, or
                 classes declared in the pins file) that plant q_ij =
                 E[L_i L_j]/E[L_tot^2] and xi_i = E[L_i]/E[L_tot] at
                 n = 4, 6 with jackknife SEs (run_truth) and write
                 two-site-ready synthetic data per arm — a VCF-like unphased
                 genotype table plus the true 'pos ref alt anc conf'
                 polarization table — under recombination ON/OFF and one
                 gene-conversion arm, with per-arm calibrated mu and r
                 (run_arm); one pins file (make_pins) fixes the base seed
                 every stream derives from. numpy always; scipy for the
                 fit; msprime + tskit for the generators (imported lazily:
                 the analytic side works without them).

Register: SIMULATION / DATA-UTILITY (Monte Carlo expectations known to their
quoted SE; the closed form is float, not certified). Reference set under
reference/planted/ (pins, six class truths, the gene-conversion arm's record,
the reduced selftest's expected numbers); all outputs go to a caller-given
directory, never into the package tree.

CLI: synth_truths.py {list,check,pin,truth,arm,all,selftest} (see --help).
"""
from ._loader import import_in_place

synth_truths = import_in_place("synth_truths")

# Convenience re-exports.
expected_branch_sfs = synth_truths.expected_branch_sfs
tavare_coefs = synth_truths.tavare_coefs
analytic_selfchecks = synth_truths.analytic_selfchecks
match_epoch_sizes = synth_truths.match_epoch_sizes
reference_checks = synth_truths.reference_checks
class_params = synth_truths.class_params
arm_spec = synth_truths.arm_spec
sim_kwargs = synth_truths.sim_kwargs
branch_class_lengths = synth_truths.branch_class_lengths
derive_seed = synth_truths.derive_seed
make_pins = synth_truths.make_pins
load_pins = synth_truths.load_pins
measure_map_mean = synth_truths.measure_map_mean
run_truth = synth_truths.run_truth
run_arm = synth_truths.run_arm
selftest = synth_truths.selftest
load_reference_truth = synth_truths.load_reference_truth
load_reference_arm = synth_truths.load_reference_arm
gz_payload_sha256 = synth_truths.gz_payload_sha256
versions = synth_truths.versions
CLASSES = synth_truths.CLASSES
ARMS = synth_truths.ARMS
TRUTH_NS = synth_truths.TRUTH_NS
TK1_EPOCHS = synth_truths.TK1_EPOCHS
TK2_TIMES = synth_truths.TK2_TIMES
REFERENCE_DIR = synth_truths.REFERENCE_DIR
REGISTER = synth_truths.REGISTER


def normalized_sfs(epochs, n):
    """Normalized closed-form (float) expected SFS xi_i = E[L_i] / sum_j E[L_j],
    i = 1..n-1, under the piecewise-constant history `epochs` (list of floats,
    sums to 1)."""
    xi = expected_branch_sfs(epochs, n)
    return (xi / xi.sum()).tolist()


__all__ = ["synth_truths", "expected_branch_sfs", "tavare_coefs",
           "analytic_selfchecks", "match_epoch_sizes", "reference_checks",
           "class_params", "arm_spec", "sim_kwargs", "branch_class_lengths",
           "derive_seed", "make_pins", "load_pins", "measure_map_mean",
           "run_truth", "run_arm", "selftest", "load_reference_truth",
           "load_reference_arm", "gz_payload_sha256", "versions",
           "normalized_sfs", "CLASSES", "ARMS", "TRUTH_NS", "TK1_EPOCHS",
           "TK2_TIMES", "REFERENCE_DIR", "REGISTER"]

"""clinch.adapters.certpass — certificate-pass adapter package for
externally fitted point-process optima (coded-Poisson multiplier fits,
conditional-spatial multinomial fits, incomplete-data ETAS reference
fits). Consumes clinch (engine/cert/jets) and
baller.certify.block_krawczyk BY IDENTITY.

Some runner/assembler scripts in this package additionally read the
originating study's fit artifacts (params, fitted-multiplier JSONs,
feature cubes) and helper modules (fit_multipliers, playoff), which are
NOT shipped here: those scripts resolve that tree through the
CERTPASS_REFERENCE_DIR environment variable and refuse loudly when it is
unset. The oracles and gates (oracle_*.py, gates_*.py) are
self-contained.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RECEIPTS_DIR = os.path.join(HERE, "receipts")
DATA_DIR = os.path.join(HERE, "data")

# clinch is the great-grandparent of this file; baller lives beside it in
# the tools tree (env-overridable), matching the etienne adapter pattern.
CLINCH_DIR = os.path.dirname(os.path.dirname(HERE))
BALLER_DIR = os.environ.get("CLINCH_BALLER_DIR") or os.path.join(
    os.path.dirname(CLINCH_DIR), "baller")
for p in (CLINCH_DIR, BALLER_DIR):
    if p not in sys.path:
        sys.path.insert(0, p)


def reference_root():
    """Root of the originating study's (un-shipped) artifact tree.

    Required by the run_*/assemble_* scripts that read fit artifacts or
    import the study's helper modules; refuses loudly when unset.
    """
    root = (os.environ.get("CERTPASS_REFERENCE_DIR")
            or os.environ.get("CERTPASS_LANE_ROOT", ""))
    if not root or not os.path.isdir(root):
        raise SystemExit(
            "certpass: CERTPASS_REFERENCE_DIR is not set (or not a directory). "
            "This script reads un-shipped study artifacts; point "
            "CERTPASS_REFERENCE_DIR at the study tree that holds round4/, "
            "stage2/, and stack/ before running it.")
    return root

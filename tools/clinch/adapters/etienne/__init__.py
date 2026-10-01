"""clinch.adapters.etienne — CLINCH oracle adapter for an ecology study's
published Etienne-model fitted optima (retro-certification adapter).

CONSUMPTION CONTRACT (read-only, by identity):
  - reference likelihood/derivative code: the study's etienne reference tree
    ({ball_engine,phase2_certify,multisample_ball,phase2_multisample_eqI}.py,
    env CLINCH_ETIENNE_DIR) imported in place (sys.dont_write_bytecode set BEFORE import — the
    etienne tree is READ-ONLY, no .pyc lands there), sha256-pinned in
    ident.PINS and verified fail-closed at first import.
  - data: the study's own data files, paths + sha256 logged per run.
  - fits of record: the study's receipted JSONs (paths in runs).
Nothing in the etienne tree is modified; certificates and the report land
in the receipts dir (env CLINCH_RECEIPTS_DIR).

ORACLES (all present the block-arrow contract of
baller.certify.block_krawczyk; NLL = -lnP so certificates are strict local
MINIMA of the negative log-likelihood):
  oracle_ss     single-sample Etienne (lntheta border, lnI block) — the
                phase2_certify.lnP_derivs2 exact term-wise gradients/
                Hessians, consumed by identity.
  oracle_eqi    equal-I multi-sample (same 2-dim structure; kernel Mtilde
                built once, I-independent) — phase2_multisample_eqI.
  oracle_corner the m=1 (Ewens) boundary fits in (lntheta, q=1/I)
                coordinates: q=0 is a declared lower bound, the corner
                (KKT) register; closed-form K_J=1, K_{J-1}=sum_{n_i>=2}
                n_i/2 at the boundary (gated against K_ball).
  oracle_diffi  the Panama-65 different-I fit: 66 parameters
                (lntheta border + ONE dense 65-dim lnI block). NOTE: the
                I_p x I_q Hessian couplings are NONZERO (shared species +
                the (theta)_A coupling), so the sound presentation is one
                dense latent block, NOT 65 independent dim-1 blocks;
                derivative machinery is the
                product-tree context scheme documented in oracle_diffi.
"""
import sys as _sys

_sys.dont_write_bytecode = True  # law: the etienne tree is READ-ONLY

import os as _os
# Dirs resolve relative to this file / from the environment; the study's
# etienne tree and the receipts store are env-configured and not shipped —
# consumers fail loudly on absent paths.
ADAPTER_DIR = _os.path.dirname(_os.path.abspath(__file__))
CLINCH_DIR = _os.path.dirname(_os.path.dirname(ADAPTER_DIR))
ETIENNE_DIR = _os.environ.get("CLINCH_ETIENNE_DIR", "")
PILOT_DIR = ETIENNE_DIR + "/pilot"
DATA_DIR = ETIENNE_DIR + "/data"
SCIENCE_DIR = ETIENNE_DIR + "/science"
RECEIPTS_DIR = _os.environ.get("CLINCH_RECEIPTS_DIR", "")

if CLINCH_DIR not in _sys.path:
    _sys.path.insert(0, CLINCH_DIR)

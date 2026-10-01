r"""membound.spec — family/region registry for 2SF memory-region boundaries:
the WORKED EXAMPLE of the general omega-core engine (membound.core). A user
core needs no registry entry: describe it in JSON and run
`python3 -m membound --spec core.json` (GUIDE.md INPUTS).

Registry contract:
  input  = family + region tag + sector (+ prescription routing)
  output = boundary vector in the |omega|^{-k eps} layer basis with per-layer
           constants (Gamma/hypergeometric closed forms where they exist,
           high-precision numerics + PSLQ closure slots where they do not).

REGISTRY STATUS (honesty flags):
  P/memory/t40-static  : VALIDATED path — the c_M trivial 1-slot case;
                         omega-cores I1M/I3M below reproduce the reference
                         c_M = 1 (arXiv:2601.16256 memory sector).
  PX2/memory/*         : SLOT-DEFINED only. Same planar RPPR/PRRP omega-core
                         family (crossed graviton reroutes which D_ij carry the
                         retarded factors); constants not provided in
                         this release (needs the crossed-family master list).
  NP/RR/*              : SLOT-DEFINED only. Three active gravitons at the S3
                         vertex -> 3 K-lines with independent frequencies. The
                         n_freq=3 core engine it needs (sign octants,
                         w4 = w1+w2+w3) IS built — OmegaCoreSpec(..., n_freq=3);
                         constants not provided in this release (needs
                         the NP master list).
"""
from .core import OmegaCoreSpec

# --- omega-cores for the P-family memory sector (gamma-3 prescription) -------
# I1^(M) core: (0+ - i w1)^{1-2e} (0+ - i w2)^{1-2e} prod_{i=1..3} |w_i|^{-e} K_e
I1M_CORE = OmegaCoreSpec(
    name="I1M",
    ret=[(1, 1, 2, 'ret'), (2, 1, 2, 'ret')],
    klines=[1, 2, 3],
    poly=None,
    proj='re',
)

# I3^(M) core: w1^3 w2^3 (0+ - i w1)^{-1-2e} (0+ - i w2)^{-1-2e} prod |w_i|^{-e} K_e
I3M_CORE = OmegaCoreSpec(
    name="I3M",
    ret=[(1, -1, 2, 'ret'), (2, -1, 2, 'ret')],
    klines=[1, 2, 3],
    poly={1: 3, 2: 3},
    proj='re',
)

REGISTRY = {
    ("P", "memory", "t40-static"): {
        "status": "VALIDATED (c_M reproduction path)",
        "cores": {"I1M": I1M_CORE, "I3M": I3M_CORE},
        "prescription": {"D31": "ret->vertex", "D24": "ret->vertex"},
        "layer_k": I1M_CORE.layer_k,  # = 7: matches the 2^{7eps} prefactors
        "masters": ["I1M", "I2M(IBP)", "I3M"],
    },
    ("PX2", "memory", None): {
        "status": "SLOT-DEFINED; constants not provided in this release",
        "cores": {},
        "notes": "planar RPPR/PRRP with one crossing; same 2-freq engine",
    },
    ("NP", "RR", None): {
        "status": "SLOT-DEFINED; n_freq=3 engine available - "
                  "constants not provided in this release",
        "cores": {},
        "n_freq": 3,
        "notes": "S3-symmetric 3-graviton vertex; sign octants, w4=w1+w2+w3; "
                 "build cores with OmegaCoreSpec(..., n_freq=3)",
    },
}

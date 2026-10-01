"""Seed derivation for the control sampler.

RNG seeds are fixed in advance in a seeds JSON file ({"seeds":
{"control_sampler_master": "<hex>"}}; not included in the package -- the
tests read its path from TERRIER_SEEDS_JSON); per-bin stream =
PCG64(seed XOR sha256(geometry|bin)[0:16]).
PINNED READING:
  master = int(seeds["control_sampler_master"], 16)            # 256-bit
  tag    = int(sha256("{geometry}|{bin}".encode()).hexdigest()[0:16], 16)
                                                               # 64-bit
  stream seed = master XOR tag; bit generator = numpy PCG64(stream seed).
Draws used by the sampler are raw 64-bit integers mapped to exact dyadic
rationals k/2^64 (fractions.Fraction) -- no float ever enters the sampler.
NO census data is touched anywhere in this module.
"""
import hashlib
import json
from fractions import Fraction

import numpy as np
from numpy.random import PCG64, Generator

TWO64 = 2 ** 64


def load_master(seeds_json_path, key="control_sampler_master"):
    """Load the master seed (hex string) from the seeds JSON as an int."""
    with open(seeds_json_path) as fh:
        seeds = json.load(fh)["seeds"]
    return int(seeds[key], 16)


def stream_seed(master_int, geometry, bin_label):
    """Pinned derivation: master XOR first-16-hex of sha256(geometry|bin)."""
    tag_hex = hashlib.sha256(
        "{}|{}".format(geometry, bin_label).encode("utf-8")
    ).hexdigest()[0:16]
    return master_int ^ int(tag_hex, 16)


def bin_generator(master_int, geometry, bin_label):
    """The fixed per-(geometry, bin) PCG64 stream."""
    return Generator(PCG64(stream_seed(master_int, geometry, bin_label)))


def raw64(gen):
    """One raw 64-bit draw from the stream, as a python int."""
    return int(gen.integers(0, TWO64, dtype=np.uint64))


def dyadic_uniform(gen):
    """Exact uniform dyadic rational in [0, 1): raw64/2^64 as a Fraction."""
    return Fraction(raw64(gen), TWO64)

"""Kira SYSTEM_*.gz weight decode (empirical, validated on production artifacts).

Kira 3.1, integral_ordering 5, WEIGHTBITS (A,B,C,D)=(5,4,17,13) as recorded in
results/kira.db of the validation run (a production 2-loop family).

Empirically established layout of the 64-bit term weight w (direct probes
against kira's own tables):

    w >> 39          : sector-block ordinal (monotone in t; NOT needed — the
                       SYSTEM term line carries the sector mask explicitly)
    (w >> 34) & 0x1F : dots  (= r - t; 5 bits = WEIGHTBITS.A)
    (w >> 30) & 0xF  : s     (sum |negative powers|; 4 bits = WEIGHTBITS.B)
    w & (2^30 - 1)   : intra-class enumeration index (C+D = 17+13 = 30 bits)

    w < 2^30         : MASTER-class ordinal (masters/preferred get small
                       weights; dots/s NOT encoded — resolve via the explicit
                       masters file of the run).

Validation:
  - all 54 masters' explicit (sector, dots, s) appear among decoded weights;
  - decoded s <= s_box + 1 and dots <= dots_box(t) + 1 for every one of the
    ~120k distinct (w, sector) pairs — exactly the one-step IBP width, and a
    bit-shifted decode breaks both bounds;
  - 161 distinct sector masks == kira's nontrivial reduce tree.
"""

MASTER_WEIGHT_CEILING = 1 << 30


def popcount(x: int) -> int:
    return bin(x).count("1")


def decode(w: int):
    """Return (dots, s) for a packed weight, or None if master-class."""
    if w < MASTER_WEIGHT_CEILING:
        return None
    return (w >> 34) & 0x1F, (w >> 30) & 0xF

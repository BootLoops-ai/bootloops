"""Kira SYSTEM_*.gz weight decode (empirical, version-pinned).

The packing decoded here is Kira's own term-weight convention (Kira:
P. Maierhöfer, J. Usovitsch, P. Uwer; J. Klappert, F. Lange; Z. Wu —
[Kira1, Kira2, Kira3] in REFERENCES.md; gitlab.com/kira-pyred/kira,
GPL-3.0-or-later). The layout below was established empirically from run
artifacts (SYSTEM_*.gz dumps and kira.db metadata), as the validation notes
say; this file contains no Kira source code.

Kira 3.1, integral_ordering 5, WEIGHTBITS (A,B,C,D)=(5,4,17,13) as
recorded in the results/kira.db of the pinned artifact trees this decode
was validated on.

Empirically established layout of the 64-bit term weight w:

    w >> 39          : sector-block ordinal (monotone in t; NOT needed — the
                       SYSTEM term line carries the sector mask explicitly)
    (w >> 34) & 0x1F : dots  (= r - t; 5 bits = WEIGHTBITS.A)
    (w >> 30) & 0xF  : s     (sum |negative powers|; 4 bits = WEIGHTBITS.B)
    w & (2^30 - 1)   : intra-class enumeration index (C+D = 17+13 = 30 bits)

    w < 2^30         : MASTER-class ordinal (masters/preferred get small
                       weights; dots/s NOT encoded — resolve via the explicit
                       masters file of the run).

Validation (see M-P runlog):
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

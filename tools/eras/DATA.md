# ERAS — package data

One data file comes with the package. It is plain JSON produced for this package; it contains no code.

| file | size | sha256[:16] | what it is |
|---|---|---|---|
| `QUARTET_COLLAPSE.json` | 1 kB | `7ff6be4be9e3f519` | the 318 site-pattern counts (quartet-collapse summary), the fixed input of the `eras.py` worked example and of `verify_eras_adversarial.py` |

## Source and license

The counts are derived from the DravLex v1.0 Dravidian lexical database: Kolipakam, Jordan,
Dunn, Greenhill, Bouckaert, Gray & Verkerk, R. Soc. Open Sci. 5 (2018) 171504; dataset
lexibank/dravlex CLDF v1.0, doi:10.5281/zenodo.5121580; distributed under CC BY 4.0 and
used here with attribution (see also THIRD_PARTY.md at the repository root). Only the
derived quartet-collapse summary ships; the database itself is not redistributed. The same
file ships, sha-pinned, as `tools/baller/vendor/ling_onesided/QUARTET_COLLAPSE.json`;
it is reference data and is not meant to be edited.

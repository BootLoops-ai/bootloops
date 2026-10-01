"""popcorn._pins — SHA-256 pins of the shipped certified files.

The certified engine and certificate modules carry the package's correctness
claims, so their shipped bytes are pinned and popcorn.verify() fails closed on
any drift: a tampered or half-edited copy is a named refusal, never a silent
fallback. Regenerate pins (only when deliberately changing shipped bytes)
with:  python3 -c "import popcorn._pins as p; p.regen()"
"""
import hashlib
from pathlib import Path

PKG_DIR = Path(__file__).resolve().parent

# Engine files (certified SFS/DFE/dominance machinery; exact Lambda-coalescent;
# float-validated two-locus and transient engines; EPO join utility; exact
# linked 2-SFS engine, hull projection and class certificates; exact folded
# two-site deciders; exact two-site projection + spectrum estimator;
# certified transient-SFS enclosure engine; planted-truth generator and
# two-window simulation harness (simulation register)).
ENGINE_SHA256 = {
    "sfs_engine.py":        "cff99bdfd1f530a9ca9e4bd0328ca176a47adafb5465810d7998e09bf195a0df",
    "arb_route.py":         "8ecd6a94cd2d898259229e85607533b871a99d623544075b7f962915517cea92",
    "certquad.py":          "a8aa92e113d6a58c53939e0e5954d70741bc2f12f90dd46e4f05a0608f73bbb2",
    "dfe_layer.py":         "ee35819b363eed18f0a9553697301de115b4f1e82eb0583c4c6597b1683396d5",
    "dominance_oracle.py":  "9b4fe845b8887ee935eb052994c277c20b761f921bbf19f35b920a644888271b",
    "dominance_qseries.py": "db7961cf79c18458089863ff2f66ad7e5780c45387da415961d012d2a3c6d1fb",
    "gate_phase1.py":       "e9d4e991f5c818163de1057d6722b5d8a4fe9f318fcae767f1cafa163a95a114",
    "lambda_exact.py":      "493d22396b631e100158635543a9573756967b3285f052193a55438e80c25d34",
    "twolocus_engine.py":   "5bb9dc13c70e3fb53a891e53c2095da5a80c9c47f43dc136f9da9961027e12a7",
    "transient_sfs_engine.py": "f669a1152fc92ab26e096ac342e537ee8025ddc1ffe4ea2356a514937e67b40b",
    "epo_join.py":          "e5cbecdf871b5d18d486dcaad4877365a251c0762035aacc7986fd7811ccfd6c",
    "twosfs_engine.py":        "26da059b0c71db87496236e49ff6c59743ddaed2248cbf1caf6ac53b04156296",
    "twosfs_hull.py":          "b26128b08fb8a49e0e021e9690f04818c6766795ee4b9ce2990e19263dc03c20",
    "twosfs_certificates.py":  "ec24cc1a9e1ba61e681473bd09e3ebe05438250041a869a6b72b154084f16144",
    "foldgate_engine.py":      "a4b537617a7479237e53af96950a3c5a5a0844df47697b84fafe19e885ba217d",
    "twosite_projection.py":   "aa72890298fbbcf1a36598175baeb5a31d2d608788ce3d5fa59968a657b67c26",
    "twosite_spectrum.py":     "6d16ef086f3f8a91ac3dd69b166f3236b6355f7717e83f2e5efd81ddbf7c0a3c",
    "transient_enclosure_engine.py":  "60a7c023ce57ef6c71638f3d9ce73c9b1b20f97033d7dd4db4fa7e3e8122b13c",
    "synth_truths.py":                "6a2e8bc301beec10e8ebdb55d44895b9a3a7101fb4f52841df32db7d2b4cd55e",
    "twowindow_sim.py":               "a5f0efddc5ceef856a81584d9265094cb07445f0ae657fa1a6c6078cd57e4acf",
}

# Certificate-instrument files.
CERT_SHA256 = {
    "positivity.py":           "30c65025008ef443815009ca254162b57352d0bf0473dea0e820c2709e75792c",
    "exact_lp.py":             "cbef74287483ee9e04abad6b102e6387be2ff6c56e4070cfc2073e3a4f8a0d4c",
    "fast_lp.py":              "93b3f2903c5c74a9206498d88915f32b15843e8bd77a8ac78a98c37ee39d02b1",
    "cone_lp_b.py":            "0297708871ad5450c62864e22599a72deac4cd78954f77a5eae2575da0acbc79",
    "cert_lp.py":              "2e0c48f13ace6d7e63dab5c31d08cc7bff2a252a851eae05a736e31a3e16f5fb",
    "region_certificates.py":  "2df2d14a6e7ae395556adf6b211a7f66568d064d343df6b15a68098d9d34f1ad",
}

# Reference data shipped with the package.
PKG_SHA256 = {
    "REFERENCE_VALUES.json": "22c019a601a6efafbebad0df0b036bcf272690df7a38a7615d929a6b499794fc",
    "reference/ancestral/epo_mini.fa": "8cfe0ba0af07fcd94a36f5eb6f20f2b3942c9fbfc2896a1ee7898a80a12ba698",
    "reference/ancestral/epo_mini_expected.json": "b0565620926b96b1645b860dfb998c1690ab649a52a2e87565ad8e320a063248",
    "reference/ancestral/epo_mini_sites.tsv": "76b9e6ae9d0dd21cbf6ac1a555cfa35266a28644c7cfc85d1a6d82f1e233cd3d",
    "reference/foldgate/folded_reference_n456.json":    "c5d77ec726a006dbea728e0630d93c7ae02547eb3dc623f0d9128506aaf7df81",
    "reference/foldgate/kingman_class_polys_n456.json": "e352f6be7134810127a2cced9d3a49a6735dbaf1ec0b94972bfed8fa8baec2df",
    "reference/lambda_coalescent/kingman_class_functional_n20.json": "1646ececdb74ad5bab6a9efe255c4b2aa842cd175a75932baf772b2301da977e",
    "reference/lambda_coalescent/spectra_n20.json": "aaab32601e123d7a1357cb0561a4d90f7fc1154a58cad87813067160af714df0",
    "reference/transient/certref_n1000_S-1.json": "83ae72115bdb636691f0df83f8ae97feb530457e5c07677c7e9596897593b539",
    "reference/transient/certref_n1000_S-100.json": "00adbba5182fe2525ff68437f889bcb1a192dba49d497be03a63fb3b87afc3d4",
    "reference/transient/certref_n1000_S-20.json": "f6a85fe3c9ed78543192e9ee7660d769de3184fdbecbe340faa3259df9a750ad",
    "reference/transient/certref_n1000_S-5.json": "68b0873a78dbce11b99d98dc39d10240a707a9f531ed50fb8228378882ac57e0",
    "reference/transient/certref_n1000_S0.json": "025522f81aec4d2aceb79957c0c826d5f994770022be0768e96ff0dbf97106dc",
    "reference/transient/eq_valid_n1000.json": "6732d8faee65b5d26b041236a2286880baa82f9d6b1a049dc8aa25552fabc218",
    "reference/transient/selfconv_n1000.json": "ea8b5c1ed7b32080a37af3faff5bea51bd13b8a2d449c06aa78071e65b7700cc",
    "reference/twolocus/kingman_n8_exact.json": "74a708d300cf10f717fae095a56214b20314bb5157a87deca11bef4c9ef66548",
    "reference/twosfs/class_certificates.json":         "0e0e0524ae58bfa4a01876191a39012434cdf7977af5c24ccef65738529ace02",
    "reference/twosfs/engine_checks.json":              "f16f22e97c89810f8ebc1bdf386e77762440418b1d7b60d2afeace37adac2a72",
    "reference/twosfs/msprime_moments.json":            "70474e5baccadc89672e2c688a0e2b05a253173df5895f2218adddc9efbe2681",
    "reference/twosfs/support_hull.json":               "0e5a0f6fd316200fa4926b081ca366ff0f56021b63fb407530cee8ca79f7041b",
    "reference/twosite/brutecheck_expected.json":       "e6a215e0228643b79caf61e3d6fbacd967aab1f1b6fe3f54f808e0e0b7f667ec",
    "reference/twosite/spectrum_expected.json":         "b6d92b7aad4191c54e321be3426a6e163fa79bda7672e924719956a6e6f3280d",
    "reference/twosite/tiny_anc.tsv":                   "465230ca9e2306979652836adc870f7c51065f1b608976e31bbbf669fff357b2",
    "reference/twosite/tiny_geno.tsv":                  "b124a85aaabc44a169e84ca1765423552eb3f0e253ac6fe2cd13451f78be32e9",
    "reference/twosite/tiny_map.bedgraph":              "600c79b5b2980f0418bdea7145ac64c30cb225895db10b31f65342194972c403",
    "reference/enclosure/cell_M200_S0_n20.json":      "0fb0c457c262be8da84a9ddf1f67abdbd28becc46d7d8eaf1295727fe1042c8b",
    "reference/enclosure/cell_M200_S-5_n20.json":     "7d48cfd9edf70c847a8326e38d83d31f35d74e7a352a4f169f83912d6987422d",
    "reference/enclosure/cell_M200_S-50_n20.json":    "f0f478d3fb18c3cf2302abcd140fe2d765d5b7aa1f67cc7447e71cf7a851b968",
    "reference/enclosure/cell_M500_S0_n20.json":      "65597c5f89306958cb66ff8170021144409e36dfa2b586a1f5dd8b67bc4ca996",
    "reference/enclosure/cell_M500_S-5_n20.json":     "56ee6675e41467624a6351438b38a146986f2a1085b778c573e26aa6a5d0d421",
    "reference/enclosure/cell_M500_S-50_n20.json":    "0219f6fd8db5969495e88a416faa351355bfc1a31cd486ea2da686e2cabbf441",
    "reference/enclosure/xcheck_M40_S-5.json":        "d66a3cc1c390832e38303963ad9cef614ea41f6f42734ae4a7284894af447298",
    "reference/planted/pins.json":                    "a3e7d22bcdf6d496a10e2098ad1737e9467cf756198f9d01632c0840bfefe3f3",
    "reference/planted/truth_TK1.json":               "681442af6ca2e23a8944a2c5d57df5a143a580b4d63af49af2dde166c802a75c",
    "reference/planted/truth_TK2.json":               "e1281c4e6ddfc1892c452069b8d779a515339a928202ad43b4b2121ea25f45a5",
    "reference/planted/truth_B13.json":               "0959a3c789331a0e659f55945bf57054edb6c4b37a77c5ac5730b9095280ac7a",
    "reference/planted/truth_B17.json":               "6929a4f843d0c16acd229bcf18a86514736128d4ceb19b5126e66cc84adcac23",
    "reference/planted/truth_D005.json":              "8877b4a6e539b66665978f5b62cebb9936e92af6ae6cf17f2b03ed433883ea79",
    "reference/planted/truth_D020.json":              "c3dd2fe7fd53e299e57a912fb9757f9aa2a564a852cd1e7cd86854e504cd368c",
    "reference/planted/arm_TK1_rON_GC.json":          "1f9cf5df444d8fd01652148d2434df7da5a641ec219a67afc97126cd73eecde0",
    "reference/planted/selftest_expected.json":       "c0c69ded29de1ec800ff867168ce70c4babf9d2848a2e4d62286dfe26ef76146",
    "reference/twowindow/linked_kingman_n8.json":     "05e106c75fc8da5d09d9488b672c3b32c2a141365543e00d48697fc919ae2e34",
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(verbose=False):
    """Check every pinned file against its SHA-256. Returns [] on clean, else a
    list of 'name: got != pinned' strings. Callers should treat ANY mismatch
    as fatal — the pinned bytes back the certification claims."""
    bad = []
    for pins in (ENGINE_SHA256, CERT_SHA256, PKG_SHA256):
        for name, want in pins.items():
            p = PKG_DIR / name
            if not p.is_file():
                bad.append(f"{name}: MISSING at {p}")
                continue
            got = sha256_file(p)
            if got != want:
                bad.append(f"{name}: sha256 {got[:16]}... != pinned {want[:16]}...")
            elif verbose:
                print(f"PASS sha256 {name}")
    return bad


def regen():
    """Print a fresh pin block for the current shipped bytes (maintainer aid;
    paste the values back into this file AND certsfs.PINNED_SHA256)."""
    for pins in (ENGINE_SHA256, CERT_SHA256, PKG_SHA256):
        for name in pins:
            print(f'    "{name}": "{sha256_file(PKG_DIR / name)}",')

"""trust._pins — sha256 pins for the vendored engine set (fail-closed)."""

PINS = {
    "lpsyz/bessel_oracle2.py":
        "befa2fea68d37584ca0eb9add5512fd4b7e6889e3c94f174a8c62c16a48a9213",
    "lpsyz/lp_syz.py": 
        "2abf0192c1bc60f23be2dd6ecbda54bd6cdc05a7573985520770b6675e34f0ce",
    "lpsyz/lp_syz_431.py": 
        "f2ef826821d5f7350ffe9b1468d512d34c5c2e90e60592d4c806c732efc59892",
    "lpsyz/lp_syz_prod.py": 
        "56b86d89dffa761ec3904ad3d69e332bf29c1ee148238d228ef02aa3295b2d7a",
}

# Identity-linked LIVE tools (fronted in place, never shadowed). Path-pinned;
# shas are ADVISORY (live tools advance independently; drift is
# reported loudly by trust.verify(), never fatal).
# strata and receipt resolve package-relative (both are MEMBER trees beside
# this package: <tool>/strata and <tool>/receipt); TRUST_STRATA_ROOT /
# TRUST_RECEIPT_ROOT override, TRUST_IBPLAPPER_ROOT comes from the
# environment. Unset entries are skipped by the advisory linked-drift check
# (they refuse fail-closed at import time through _core.IMPORT_ROOTS instead).
import os as _os
_TOOLPKG = _os.path.abspath(_os.path.join(
    _os.path.dirname(_os.path.abspath(__file__)), ".."))
LINKED = {k: v for k, v in {
    "strata": (_os.path.join(_os.environ["TRUST_STRATA_ROOT"], "strata")
               if _os.environ.get("TRUST_STRATA_ROOT")
               else _os.path.join(_TOOLPKG, "strata")),
    "ibplapper": _os.environ.get("TRUST_IBPLAPPER_ROOT", ""),
    "receipt": _os.environ.get("TRUST_RECEIPT_ROOT") or _os.path.join(_TOOLPKG, "receipt"),
}.items() if v}

# ADVISORY baseline shas of the LINKED live-tool files trust depends on
# (drift is definable only against a recorded baseline). Relative to
# LINKED[name]. verify() REPORTS drift loudly in report["linked_drift"] and
# the advisory print, NEVER raises (the advisory pattern: live tools
# advance independently). Re-baseline with a recorded reason only.
# Re-baseline record: receipt/adapters/strata.py and WITNESS_FORMAT.md (both
# copies, still byte-identical) re-baselined when the receipt package became
# the receipt MEMBER of this package (tools/receipt -> tools/trust/receipt):
# the adapter's default strata path and one path mention in the format spec
# moved with it; no code bytes of core.py changed.
LINKED_SHAS = {
    "strata": {  # loader/CoeffEvaluator + eliminator + certify: the fronted files
        "__init__.py": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        "loader.py": "049c369343dd906a5f09ff03d5d616ba05874aa9666dfed8f1476fc4a70f4c44",
        "fp_eliminate.py": "1d3dc417449347c93c4fc605307fc8081b76ec55cf6ca369bd647bee69a94846",
        "certify.py": "47a505dfc74dbff5f673f02106cd811e087364cbf40fcd5a13606f6cb0aa7f1d",
    },
    "ibplapper": {  # winnow's vendored receipt core (the site's winnow.receipt)
        "ibplapper/receipt/core.py": "96b0a223365cfb0be0a513a52880b34ced31269130c111ee63db689fa322beb0",
        "ibplapper/receipt/WITNESS_FORMAT.md": "59562276a88e11a25f482705bf241938a67340476a4e315c5ddb583d96725052",
    },
    "receipt": {  # solver-agnostic standalone core + the adapter check-2 rides
        "core.py": "96b0a223365cfb0be0a513a52880b34ced31269130c111ee63db689fa322beb0",
        "adapters/strata.py": "86896eb8f743289cbbc2fb033bd5b12df76e95180ca47980c16f2a520ae1778a",
        "WITNESS_FORMAT.md": "59562276a88e11a25f482705bf241938a67340476a4e315c5ddb583d96725052",
    },
}

"""ibplapper.receipt — the RECEIPT verifier/emitter component, VENDORED.

RECEIPT is the library's native verifier/emitter; the solver-agnostic
standalone CLI lives at tools/trust/receipt in this repository (the
receipt member of the trust package; companion, same contract).

Vendoring rule: core.py / emitter.py / detector.py / WITNESS_FORMAT.md are
BYTE-IDENTICAL copies of the canonical tools/trust/receipt files (import,
never fork). VENDOR_MANIFEST.json records the source sha256s
and tests/test_vendor_integrity.py fails if either side drifts.

Convenience API (glue only — no verification logic lives here):

    from ibplapper import receipt
    ok, detail = receipt.verify(rows, witness)      # dict or path
    rc, report = receipt.verify_many(rows, wits)
    fpr = receipt.system_fingerprint(rows, p)
"""
from . import core, detector, emitter  # noqa: F401
from .core import (MalformedInput, RC_FAIL, RC_MALFORMED, RC_OK,  # noqa: F401
                   WITNESS_VERSION, load_system_jsonl, load_witness,
                   system_fingerprint, verify_many, verify_row,
                   write_witness)
from .detector import detect  # noqa: F401


def verify(rows, witness, table_row=None, check_fingerprint=None):
    """Verify one witness (loaded dict OR path to a v1 witness JSON) against
    the caller's OWN parse of the system rows. Returns (passed, detail).

    The independence contract (the wrong-table lesson): `rows` should come
    from the caller's independent parse of the generating system, not from
    the solver that produced the table being checked.
    """
    if isinstance(witness, str):
        witness = load_witness(witness)
    return verify_row(rows, witness, table_row=table_row,
                      check_fingerprint=check_fingerprint)

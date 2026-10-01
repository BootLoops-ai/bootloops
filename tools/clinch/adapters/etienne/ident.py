"""Identity consumption of the etienne pilot modules (read-only, sha-pinned)
+ the read-audit log every certificate run appends to.

The reference derivative code is IMPORTED IN PLACE (never copied); the
sha256 pins below are verified fail-closed before the first import. A pin
mismatch is a typed refusal — certifying against drifted reference code
would be untestable.
"""
import hashlib
import json
import os
import sys
import time

sys.dont_write_bytecode = True

from . import PILOT_DIR, RECEIPTS_DIR


class AdapterIdentityError(RuntimeError):
    pass


# sha256 of the consumed reference files (disk-verified pins)
PINS = {
    "ball_engine.py":
        "54c06b31cb1e010d70525c499d4fb97f348853594fd460cdd6e88706e6c18c6f",
    "phase2_certify.py":
        "0b33aeb7ec65d5898ed638373b38ae4e04232a9036bee3c6c5d0f83d5e3855cb",
    "multisample_ball.py":
        "4e69ae92ca3fcf3dd9917de997847c38a59de4ef33b27c211306f364d3a8d521",
    "phase2_multisample_eqI.py":
        "3aa5917add556c0bfa1d30d7987b7be2072d195801793fa118b1dd3b6e06f1e1",
}

_READ_LOG = []


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log_read(path, why):
    """Audit-log every etienne-tree path this adapter reads."""
    _READ_LOG.append(dict(path=os.path.abspath(path), why=why,
                          sha256=(sha256(path) if os.path.isfile(path)
                                  else None),
                          stamp=time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                              time.gmtime())))


def flush_read_log(tag):
    os.makedirs(RECEIPTS_DIR, exist_ok=True)
    p = os.path.join(RECEIPTS_DIR, f"READ_LOG_{tag}.json")
    json.dump(_READ_LOG, open(p, "w"), indent=1)
    return p


_verified = False
_mods = {}


def verify():
    global _verified
    if _verified:
        return
    for name, pin in PINS.items():
        path = os.path.join(PILOT_DIR, name)
        got = sha256(path)
        if got != pin:
            raise AdapterIdentityError(
                f"reference pin mismatch on {path}: got {got}, pinned "
                f"{pin} — the etienne reference moved; re-pin only after "
                f"re-gating")
        log_read(path, "identity-consumed reference module (sha verified)")
    _verified = True


def pilot():
    """Import the pinned pilot modules in place; returns a namespace dict."""
    verify()
    if _mods:
        return _mods
    if PILOT_DIR not in sys.path:
        sys.path.insert(0, PILOT_DIR)
    assert sys.dont_write_bytecode, "read-only law: bytecode writes off"
    import ball_engine
    import phase2_certify
    import multisample_ball
    import phase2_multisample_eqI
    for m in (ball_engine, phase2_certify, multisample_ball,
              phase2_multisample_eqI):
        assert m.__file__.startswith(PILOT_DIR), m.__file__
    _mods.update(ball_engine=ball_engine, phase2_certify=phase2_certify,
                 multisample_ball=multisample_ball,
                 phase2_multisample_eqI=phase2_multisample_eqI)
    return _mods

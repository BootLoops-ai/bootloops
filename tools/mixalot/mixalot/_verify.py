"""Fail-closed vendor integrity + advisory-loud live-source drift check."""
import hashlib
import os

from ._pins import PINS

_VENDOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor")
_verified = {"ok": False}


class VendorTamperError(RuntimeError):
    pass


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(quiet=False):
    """Re-hash the vendored engine set against the pins. FAIL-CLOSED on any
    vendor mismatch (raises VendorTamperError). Upstream sources are
    checked ADVISORY-LOUD: drift is reported in the returned dict (and
    printed unless quiet), never raises. Returns the integrity report."""
    report = {"vendor_ok": [], "vendor_tampered": [],
              "source_drift": [], "source_missing": []}
    for rel, (source, pinned) in PINS.items():
        vp = os.path.join(_VENDOR, rel)
        if not os.path.exists(vp):
            raise VendorTamperError(
                f"vendored file MISSING: {rel} (incomplete copy or deleted "
                "file — re-obtain the full MIXALOT bundle)")
        h = _sha(vp)
        if pinned is None or h != pinned:
            report["vendor_tampered"].append((rel, h, pinned))
        else:
            report["vendor_ok"].append(rel)
        if os.path.exists(source):
            if _sha(source) != pinned:
                report["source_drift"].append((rel, source))
        else:
            report["source_missing"].append((rel, source))
    # ANY unpinned file under vendor/ is a tamper class (planted .pyc/.so/.pth/
    # case-variants), and SYMLINKS anywhere under vendor/ break self-containment
    # (a dir symlink hides payloads from the walk). Kept self-contained by
    # design — each package verifies itself without importing another.
    for root, dirs, files in os.walk(_VENDOR):
        for d in list(dirs):
            p = os.path.join(root, d)
            if os.path.islink(p):
                report["vendor_tampered"].append(
                    (os.path.relpath(p, _VENDOR) + "/", "SYMLINKED-DIR", None))
                dirs.remove(d)
            elif d == "__pycache__":
                import shutil as _sh
                _sh.rmtree(p, ignore_errors=True)  # never read, never trusted
                dirs.remove(d)
        for f in files:
            p = os.path.join(root, f)
            rel2 = os.path.relpath(p, _VENDOR)
            if os.path.islink(p):
                report["vendor_tampered"].append((rel2, "SYMLINK", None))
            elif rel2 not in PINS:
                report["vendor_tampered"].append((rel2, "UNPINNED", None))
    if report["vendor_tampered"]:
        _verified["ok"] = False
        raise VendorTamperError(
            "MIXALOT vendor integrity FAILED (tool refuses to compute): "
            + "; ".join(f"{r} sha {h[:16]} != pin "
                        f"{(p or 'UNSTAMPED')[:16]}"
                        for r, h, p in report["vendor_tampered"]))
    _verified["ok"] = True
    if not quiet:
        shown_missing = [s for _, s in report["source_missing"] if not s.startswith("origin:")]
        if report["source_drift"] or shown_missing:
            print("MIXALOT ADVISORY: upstream sources differ from the "
                  "vendored snapshot (tool remains valid; the upstream tree "
                  "is the record of origin):")
            for rel, src in report["source_drift"]:
                print(f"  drift: {rel} vs {src}")
            for rel, src in report["source_missing"]:
                if not src.startswith("origin:"):   # 'origin:<name>' ids name the unshipped tree of origin; absence is expected
                    print(f"  missing source: {src}")
    return report


def require_verified():
    """Compute paths call this: runs verify(quiet=True) once per process."""
    if not _verified["ok"]:
        verify(quiet=True)

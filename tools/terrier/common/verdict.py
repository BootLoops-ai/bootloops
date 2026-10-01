#!/usr/bin/env python3
"""
verdict.py (common chassis) — two-dps verdict chassis + gate primitives:
asserting print-scoreboard gates (ball_gate, route_gate, two_dps_gate),
the ball primitives parse_ball / matched_digits, and the sha-stamped atomic
receipt.  Wing modules (periods/verdicts.py) delegate here.
Battery: selftest_verdict.py.
"""
import hashlib, json, os

def verdict(tag, ok, msg=""):
    print(f"[{tag}] {'PASS' if ok else 'FAIL'}  {msg}", flush=True)
    assert ok, f"{tag} FAILED: {msg}"
    return True

def parse_ball(sst):
    """'[m +/- r]' | '[m]' | 'm' -> (mpf mid, mpf rad) at dps 250."""
    import mpmath as mp
    mp.mp.dps = 250
    sst = sst.strip()
    if sst.startswith("["):
        body = sst[1:-1]
        if "+/-" in body:
            m, r = body.split("+/-")
            return mp.mpf(m.strip()), mp.mpf(r.strip())
        return mp.mpf(body), mp.mpf(0)
    return mp.mpf(sst), mp.mpf(0)

def matched_digits(a, b):
    import mpmath as mp
    if a == b:
        return 999            # exact midpoint match sentinel
    return int(mp.floor(-mp.log10(abs(a - b) / abs(b))))

def ball_gate(tag, mine_s, other_s, bar_digits, max_radius_ratio=16.0):
    """LIKE-FOR-LIKE regression: overlap + matched digits + comparable radius.
    NEVER apply across routes (the radius bar is meaningless between two
    different certified routes); use route_gate for cross-route checks."""
    m1, r1 = parse_ball(mine_s); m2, r2 = parse_ball(other_s)
    ov = bool(abs(m1 - m2) <= r1 + r2)
    dm = matched_digits(m1, m2)
    rat = float(r1 / r2) if r2 > 0 else 1.0
    verdict(tag, ov and dm >= bar_digits and rat <= max_radius_ratio,
            f"overlap={ov}, matched digits={dm} (bar {bar_digits}), "
            f"radius ratio={rat:.2f}")
    return {"overlap": ov, "matched_digits": dm, "radius_ratio": rat}

def route_gate(tag, a_s, b_s, bar_digits):
    """CROSS-ROUTE agreement: overlap + digits only (no radius bar)."""
    m1, r1 = parse_ball(a_s); m2, r2 = parse_ball(b_s)
    ov = bool(abs(m1 - m2) <= r1 + r2)
    dm = matched_digits(m1, m2)
    verdict(tag, ov and dm >= bar_digits,
            f"overlap={ov}, matched digits={dm} (bar {bar_digits})")
    return {"overlap": ov, "matched_digits": dm}

def code_sha(*paths):
    """sha256 stamp of the exact code files a receipt certifies.

    Delegates to common/receipts.py code_sha (SORTED concatenation), so the
    same file set always stamps the same sha regardless of argument order."""
    from receipts import code_sha as _canonical
    return _canonical(list(paths))

def two_dps_gate(tag, lo_s, hi_s, bar_digits):
    """two-dps cert: lo-dps midpoint must match hi-dps midpoint beyond bar."""
    d = matched_digits(parse_ball(lo_s)[0], parse_ball(hi_s)[0])
    verdict(tag, d >= bar_digits, f"{d} matched digits (bar {bar_digits})")
    return d

def receipt(payload, out_path, code_paths=()):
    """Receipt schema + atomic write: {payload, code_sha256, schema}.
    Atomic tmp+fsync+rename, the same receipts pattern as
    common/receipts.py."""
    rec = {"schema": "terrier-receipt-v1", "payload": payload,
           "code_sha256": code_sha(*code_paths) if code_paths else None}
    tmp = out_path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(rec, f, indent=1, default=str)
        f.flush(); os.fsync(f.fileno())
    os.replace(tmp, out_path)
    return rec

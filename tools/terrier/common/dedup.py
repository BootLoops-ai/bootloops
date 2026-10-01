#!/usr/bin/env python3
"""dedup.py (common chassis) — known-point pullback-orbit dedup for
attractor-point escalation chains.  Motivating case: a candidate at
AESZ118 z=-1/32 that turned out to be a card-documented pullback of the
known point T3.  Hence the law: dedup runs BEFORE any escalation compute,
and a newform-ID match is a REDISCOVERY, never a NEW point.

orbit_match(op, spec): EXACT match of a candidate (op, z* | minpoly)
against the packaged KNOWN_ORBIT (identity rows for the published points
+ their closure under documented exact coordinate/pullback relations)
-> "REDISCOVERY-OF-<Tk>" or None.

newform_backstop(label): identification-time backstop — the verified
newform labels of each known point (incl. 14.4.a.b (T4 wt-4) and the
wt-2 forms 34.2.b.a / 11.2.a.a).

Orbit data is packaged at common/data/KNOWN_ORBIT.jsonl (env override
TERRIER_KNOWN_ORBIT); emit_rediscovery takes an EXPLICIT receipt path
(no built-in default — the caller owns its ledger). Battery:
selftest_dedup.py (5 must-fire + 3 must-not-fire orbit controls + 5
newform controls)."""
import json, os, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ORBIT = os.environ.get("TERRIER_KNOWN_ORBIT",
                       os.path.join(HERE, "data", "KNOWN_ORBIT.jsonl"))

# verified newform labels per known point ("(twist)" = the
# documented 825 twist used in T3's rank-1 L-identities).
KNOWN_NEWFORMS = {
    "14.2.a.a": ["T1", "T4"], "14.4.a.a": ["T1"],
    "34.2.b.a": ["T2"], "34.4.b.a": ["T2"],
    "11.2.a.a": ["T3"], "33.4.a.b": ["T3"], "825.4.a.f": ["T3 (twist)"],
    "14.4.a.b": ["T4"],
}


def _load_orbit():
    return [json.loads(l) for l in open(ORBIT)]


def _isqrt_exact(n):
    import math                            # exact integer sqrt, no floats
    s = math.isqrt(n)
    return s if s * s == n else None


def _cand_keys(spec):
    """Candidate spec -> exact match keys. 'z=p/q' -> rat key.
    'minpoly=a,b,c' (a x^2+b x+c, quadratic points) -> alg key, normalized
    ascending [c,b,a], primitive, positive leading coeff; if disc is a
    perfect square the roots are rational -> rat keys instead."""
    if spec.startswith("z="):
        return ["rat:" + str(Fr(spec[2:]))]
    a, b, c = (int(x) for x in spec.split("=", 1)[1].split(","))
    disc = b * b - 4 * a * c
    s = _isqrt_exact(disc) if disc >= 0 else None
    if s is not None:
        return ["rat:" + str(Fr(-b + s, 2 * a)),
                "rat:" + str(Fr(-b - s, 2 * a))]
    import math
    g = math.gcd(math.gcd(abs(a), abs(b)), abs(c)) or 1
    cf = [c // g, b // g, a // g]
    if cf[-1] < 0:
        cf = [-v for v in cf]
    return ["alg:" + ",".join(str(v) for v in cf)]


def orbit_match(op, spec):
    """-> (verdict, matched_rows) — verdict 'REDISCOVERY-OF-<Tk>' or
    None. Run BEFORE any escalation compute (front-door law)."""
    keys = set(_cand_keys(spec))
    hits = [r for r in _load_orbit()
            if r["operator_key"] == op and r["match_key"] in keys]
    if not hits:
        return None, []
    pts = sorted({r["source_point"].rstrip("pm") for r in hits})
    return "REDISCOVERY-OF-" + "/".join(pts), hits


def newform_backstop(label):
    """-> 'REDISCOVERY-OF-<Tk>' or None. Call at newform-ID time; a match
    means the candidate is a known point: label REDISCOVERY, never NEW."""
    pts = KNOWN_NEWFORMS.get(str(label).strip())
    if not pts:
        return None
    return "REDISCOVERY-OF-" + "/".join(sorted(
        {p.split()[0].rstrip("pm") for p in pts}))


def emit_rediscovery(op, spec, verdict, hits, receipt, status=None):
    """Append the dedup verdict to the CALLER'S receipt ledger (JSONL;
    explicit path — no built-in default) and optionally a status log."""
    row = {"op": op, "spec": spec, "ts": time.strftime("%F %T"),
           "verdict": verdict, "stage": "dedup-pre-escalation",
           "matched_orbit_rows": hits,
           "note": "known-point pullback dedup: no transport, no fresh "
                   "gate, no escalation compute (documented orbit "
                   "of a known point)"}
    with open(receipt, "a") as f:
        f.write(json.dumps(row) + "\n")
    if status:
        with open(status, "a") as f:
            f.write(f"**DEDUP {time.strftime('%H:%M')}** op={op} {spec} "
                    f"-> {verdict} (dedup pre-step; escalation SKIPPED — "
                    f"documented orbit of a known point)\n")
    return row


def check(op, spec, receipt=None, status=None):
    """Front door: orbit_match, emitting to the caller's ledger iff a
    receipt path is given. -> verdict or None."""
    v, hits = orbit_match(op, spec)
    if v and receipt:
        emit_rediscovery(op, spec, v, hits, receipt, status)
    return v


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("usage: dedup.py newform <label> | dedup.py check <op> <spec> "
              "[--receipt PATH]")
        sys.exit(2)
    if sys.argv[1] == "newform":            # dedup.py newform <label>
        v = newform_backstop(sys.argv[2])
        print(v or "NO-MATCH")
        sys.exit(0 if v is None else 10)
    # dedup.py check <op> <spec> [--receipt PATH]: rc 10 = REDISCOVERY
    # (skip escalation), rc 0 = no match (proceed).
    op, spec = sys.argv[2], sys.argv[3]
    rp = sys.argv[sys.argv.index("--receipt") + 1] \
        if "--receipt" in sys.argv else None
    v = check(op, spec, receipt=rp)
    print(v or "NO-MATCH")
    sys.exit(0 if v is None else 10)

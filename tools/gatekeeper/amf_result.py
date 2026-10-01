"""
Canonical parser for raw AMFlow-cpp `solve_integrals` JSON result entries.

CONVENTION (verified against amflow-cpp/src/pipeline/solve_integrals.cpp and
1-loop-box reference outputs):

    result[i] = {
      "integral": {"family":..., "indices":[...]},
      "leading_order": <int>,             # = min(c['order'] for c in coefficients)
      "coefficients": [
         {"order": <int>,                 # ABSOLUTE eps-power (NOT relative!)
          "value": {"re":"[X +/- E]", "im":"[Y +/- E]"}},
         ...
      ]
    }

So `c['order']` IS the absolute Laurent index; `leading_order` is REDUNDANT
metadata (= the smallest order present). Adding them double-counts.

Footgun (caught live): a parser that did
`{lo+c['order']: ...}` shifted a 3L Laurent by lo and gated at 1d instead of
139d. The bug is silent when lo==0 (e.g. most 2L integrals).
"""
import re as _re

import mpmath as mp

__all__ = ["amf_laurent", "amf_arb", "amf_ball"]

# Arb's printer (arb_get_str(x, d, 0), which amflow-cpp's acb_to_json calls)
# emits exactly these shapes; the separator is a literal " +/- ".
_BALL_RE = _re.compile(r"^\[?\s*(?P<mid>[^\[\],]*?)\s*\+/-\s*(?P<rad>[^\[\]]*?)\s*\]?$")
_INTERVAL_RE = _re.compile(r"^\[\s*(?P<a>[^,\]]+?)\s*,\s*(?P<b>[^,\]]+?)\s*\]$")


def amf_ball(v):
    """Split one printed value into (midpoint, radius) DECIMAL STRINGS.

    The ONE place the ball form is parsed (oracle_integrity reuses this; do
    not inline a second parser).  Accepted forms:

        "0" / "0.25" / "-1.0e-5"    exact value            -> (s, None)
        "[X +/- E]"                  Arb ball               -> ("X", "E")
        "X +/- E"                    same, no brackets      -> ("X", "E")
        "[+/- E]"                    Arb zero-midpoint ball -> ("0", "E")
        "[a, b]"                     interval (defensive)   -> ((a+b)/2, |b-a|/2)

    The interval form is never emitted by amflow-cpp; it is converted with
    mpmath at a working precision set from the string length and returned
    as decimal strings.  Everything else is a pure string split: no
    arithmetic, no precision loss, the caller decides the dps.
    """
    s = str(v).strip()
    m = _BALL_RE.match(s)
    if m and "+/-" in s:
        mid = m.group("mid").strip() or "0"
        rad = m.group("rad").strip() or "0"
        return mid, rad
    m = _INTERVAL_RE.match(s)
    if m:
        a_s, b_s = m.group("a"), m.group("b")
        ndig = sum(c.isdigit() for c in a_s + b_s)
        with mp.workdps(max(30, ndig + 10)):
            a, b = mp.mpf(a_s), mp.mpf(b_s)
            mid, rad = (a + b) / 2, abs(b - a) / 2
            n = max(20, ndig + 5)
            return mp.nstr(mid, n), mp.nstr(rad, n)
    return s.lstrip("[").rstrip("]"), None


def amf_arb(v):
    """Parse one Arb-ball string '[X +/- E]' or {'re':...,'im':...} -> mpc
    (the MIDPOINT; the radius is dropped — use amf_ball to keep it).
    Reads at the CALLER's mp.mp.dps (set it >= the Arb mantissa width before
    calling, or you silently truncate)."""
    if isinstance(v, dict):
        return mp.mpc(amf_arb(v.get("re", "0")), amf_arb(v.get("im", "0")))
    mid, _rad = amf_ball(v)
    return mp.mpf(mid)


def amf_laurent(result_entry):
    """result_entry -> {abs_eps_order: mpc}.
    `result_entry` is one element of out_json['result']."""
    out = {}
    for c in result_entry["coefficients"]:
        out[int(c["order"])] = amf_arb(c["value"])
    # sanity: leading_order should equal min key (warn, don't fail)
    lo = result_entry.get("leading_order")
    if lo is not None and out and lo != min(out):
        import warnings
        warnings.warn(f"amf_laurent: leading_order={lo} != min(order)={min(out)}; "
                      f"using c['order'] as absolute (correct).")
    return out

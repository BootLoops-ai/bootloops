#!/usr/bin/env python3
"""
gate.py — digit gates against vendored reference strings.

Gate rule (the two-precision discipline):
  For each vendored key with a stored decimal string of n significant
  digits, a computed pair (value at dps_pair[0], value at dps_pair[1])
  PASSes iff
    (a) matched significant digits vs the FULL stored string
        >= n - 1   (the whole vendored string must be reproduced, not a
        prefix — no partial credit), AND
    (b) the two computed precisions agree to >= n digits (a value that
        moves with working precision is noise, not a result).
  Overall PASS iff every vendored key passes and none is missing.

Degenerate vendored strings (the relative-digit metric is
meaningless for them, and the old code passed them vacuously):
  * A vendored string that parses to ZERO (e.g. '0', '0.0' — the common
    vanishing-imaginary-part reference entry) is gated ABSOLUTELY: PASS iff
    BOTH computed pair members satisfy |computed| <= scale_ref *
    10^-zero_digits, where scale_ref defaults to the max |value| over the
    NONZERO vendored sibling keys (override per key or globally with
    zero_scale=) and zero_digits defaults to the max stored-digit claim over
    those siblings (override with zero_digits=).  With no nonzero sibling
    and no explicit scale the key HARD-FAILS — an un-scaled zero cannot be
    gated.
  * A NONZERO vendored string with only 1 significant digit HARD-FAILS:
    matched >= n_stored - 1 = 0 is vacuously true for ANY computed value, so
    a 1-digit reference gates nothing. Store more digits.

dps_pair must be strictly increasing (dlo < dhi; suggest dhi >= dlo + 10) —
equal precisions silently convert the two-precision integrity gate into a
one-precision gate, so they raise ValueError.  Byte-identical computed pair
strings are flagged per key ('pair_identical_strings': True): they defeat
the two-precision check by construction (a caller that computed once and
passed (v, v)).

Provenance of the design rules:
* Two-precision stability + full-stored-string matching: the solved-bar
  audit standard — a symbolic form must be checkable by an independent
  arbitrary-precision evaluation, and a value that moves with working
  precision is noise, not a result.  Proven on production 79-entry
  two-precision gate runs.
* Matched-digit metric agree_d = -log10(|a-b|/scale): the standard
  measured-agreement metric; digit counts are always measured, never
  invented.
* gate_strings WRITES NOTHING (pure function); write_gate_report is the
  only writer and only touches the path it is given.
* Precision discipline: vendored/computed strings are parsed inside
  mp.workdps at a working precision covering their full length (the mpf
  import-dps footgun — parsing at ambient mp.dps silently
  truncates).  Global mp.dps is never touched.

Vendored strings are real decimal strings (optionally signed, optional
exponent).  Complex quantities gate as two
keys (re/im).  Computed entries may be str, mpf, or mpc.
"""

import json

import mpmath as mp

__all__ = ["gate_strings", "write_gate_report", "feed_gate_table"]


def _sig_digits(s):
    """Count significant digits of a decimal string.

    Sign, decimal point and exponent are ignored; leading zeros are not
    significant ('0.000123' -> 3); trailing zeros ARE ('1.0' -> 2, they are
    a stored-precision claim).
    """
    t = str(s).strip()
    if not t:
        raise ValueError("gate: empty vendored string")
    for e in ("e", "E"):
        if e in t:
            t = t.split(e)[0]
            break
    t = t.lstrip("+-").replace(".", "")
    if not t.isdigit():
        raise ValueError("gate: vendored string %r is not a plain decimal" % (s,))
    stripped = t.lstrip("0")
    if not stripped:            # "0", "0.0", ... — a bare zero claims 1 digit
        return 1
    return len(stripped)


def _parse(x):
    """Parse str/mpf/mpc AT THE CURRENT working precision (call inside
    workdps — import-dps footgun)."""
    if isinstance(x, (mp.mpf, mp.mpc)):
        return +mp.mpc(x)
    return mp.mpc(mp.mpmathify(str(x).strip()))


def _agree_digits(a, b, cap):
    """Matched significant digits between a and b, capped at `cap`.
    Call inside workdps.  Both zero -> cap (perfect agreement)."""
    diff = abs(a - b)
    scale = max(abs(a), abs(b))
    if diff == 0:
        return float(cap)
    if scale == 0:
        return float(cap)
    d = float(-mp.log10(diff / scale))
    return max(0.0, min(d, float(cap)))


def gate_strings(computed, vendored, dps_pair, zero_scale=None,
                 zero_digits=None):
    """Gate computed two-precision pairs against vendored digit strings.

    Args:
      computed: {key: (value_at_dps_pair[0], value_at_dps_pair[1])}.
                Values may be str (preferred — full digits preserved),
                mpf, or mpc.
      vendored: {key: stored decimal string} — the FULL reference string;
                its significant-digit count sets the bar.
      dps_pair: (dps_lo, dps_hi) — the two precisions the pair was
                computed at, STRICTLY increasing (dlo < dhi; suggest
                dhi >= dlo + 10).  Order must match the computed tuples.
      zero_scale: scale reference for ZERO-vendored keys: a scalar
                (str/mpf) applied to every zero key, or {key: scalar}.
                Default: max |value| over the NONZERO vendored siblings.
      zero_digits: absolute digit bar for zero-vendored keys (int).
                Default: max stored-digit claim over the nonzero siblings.

    Returns a report dict (WRITES NOTHING):
      {'overall': 'PASS'|'FAIL',
       'rule': ...,
       'dps_pair': [lo, hi],
       'n_keys': int, 'n_pass': int,
       'per_key': {key: {'stored_digits': int,
                         'matched_digits': float,     # vs full stored string
                         'pair_agree_digits': float,  # two-precision check
                         'pair_identical_strings': bool,  # (v, v) flag
                         'PASS': bool, 'reason': str-or-None,
                         # zero-vendored keys additionally carry:
                         'zero_bar': str, 'abs_computed': [str, str]}},
       'missing_keys': [...], 'extra_keys': [...]}

    matched_digits compares the vendored value against the HIGHER-precision
    member of the pair.  PASS iff matched_digits >= stored_digits - 1 AND
    pair_agree_digits >= stored_digits.  Zero-vendored keys gate absolutely
    (|both members| <= zero_bar); nonzero 1-digit vendored strings hard-fail
    (see module docstring).

    Raises ValueError on a malformed dps_pair, including dhi <= dlo (equal
    precisions would silently degrade the gate to one-precision).
    """
    if len(dps_pair) != 2:
        raise ValueError("gate_strings: dps_pair must be a (lo, hi) tuple")
    dlo, dhi = int(dps_pair[0]), int(dps_pair[1])
    if dhi <= dlo:
        raise ValueError(
            "gate_strings: dps_pair must be strictly increasing, got "
            "(%d, %d) — equal/reversed precisions silently convert the "
            "two-precision integrity gate into a one-precision gate; use "
            "dhi >= dlo + 10" % (dlo, dhi))

    # pre-pass: parse every vendored string once to find the zero keys and
    # the sibling scale/digit references for them.
    all_len = max([len(str(v).strip()) for v in vendored.values()],
                  default=0)
    wp_pre = max(dlo, dhi, all_len) + 20
    vend_zero = {}
    with mp.workdps(wp_pre):
        sib_scale = mp.mpf(0)
        sib_digits = 0
        for key in vendored:
            sv = str(vendored[key]).strip()
            vz = (abs(_parse(sv)) == 0)
            vend_zero[key] = vz
            if not vz:
                sib_scale = max(sib_scale, abs(_parse(sv)))
                sib_digits = max(sib_digits, _sig_digits(sv))
        sib_scale_str = mp.nstr(sib_scale, 25) if sib_scale > 0 else None

    per_key = {}
    missing = []
    n_pass = 0
    for key in vendored:
        sv = str(vendored[key]).strip()
        n_stored = _sig_digits(sv)
        if key not in computed:
            missing.append(key)
            per_key[key] = {"stored_digits": n_stored,
                            "matched_digits": None,
                            "pair_agree_digits": None,
                            "PASS": False,
                            "reason": "missing from computed"}
            continue
        pair = computed[key]
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            per_key[key] = {"stored_digits": n_stored,
                            "matched_digits": None,
                            "pair_agree_digits": None,
                            "PASS": False,
                            "reason": "computed[%r] is not a 2-tuple "
                                      "(need values at both dps)" % (key,)}
            continue
        pair_identical = str(pair[0]).strip() == str(pair[1]).strip()

        if vend_zero[key]:
            # ---- zero-vendored key: ABSOLUTE gate against a scale ref ----
            zs = None
            if isinstance(zero_scale, dict):
                zs = zero_scale.get(key)
            elif zero_scale is not None:
                zs = zero_scale
            if zs is None:
                zs = sib_scale_str
            nd = zero_digits if zero_digits is not None else \
                (sib_digits if sib_digits > 0 else None)
            str_len = max(len(sv), max(len(str(p)) for p in pair))
            wp = max(dlo, dhi, str_len, nd or 0) + 20
            with mp.workdps(wp):
                a = _parse(pair[0])
                b = _parse(pair[1])
                if zs is None or nd is None:
                    per_key[key] = {
                        "stored_digits": n_stored,
                        "matched_digits": None,
                        "pair_agree_digits": None,
                        "pair_identical_strings": pair_identical,
                        "PASS": False,
                        "reason": ("zero-vendored key has no scale "
                                   "reference (no nonzero sibling key and "
                                   "no explicit zero_scale=/zero_digits=) — "
                                   "an un-scaled zero cannot be gated")}
                    continue
                bar = abs(_parse(zs)) * mp.mpf(10) ** (-int(nd))
                ok = (abs(a) <= bar) and (abs(b) <= bar)
                rec = {"stored_digits": n_stored,
                       "matched_digits": None,
                       "pair_agree_digits": None,
                       "pair_identical_strings": pair_identical,
                       "zero_bar": mp.nstr(bar, 8),
                       "abs_computed": [mp.nstr(abs(a), 8),
                                        mp.nstr(abs(b), 8)],
                       "PASS": ok,
                       "reason": None if ok else
                       ("zero-vendored: |computed| = (%s, %s) exceeds the "
                        "absolute zero bar %s = scale_ref * 10^-%d"
                        % (mp.nstr(abs(a), 4), mp.nstr(abs(b), 4),
                           mp.nstr(bar, 4), int(nd)))}
            if ok:
                n_pass += 1
            per_key[key] = rec
            continue

        if n_stored < 2:
            # ---- nonzero 1-digit vendored string: un-gateable ----
            per_key[key] = {"stored_digits": n_stored,
                            "matched_digits": None,
                            "pair_agree_digits": None,
                            "pair_identical_strings": pair_identical,
                            "PASS": False,
                            "reason": ("vendored string %r has only 1 "
                                       "significant digit — matched >= "
                                       "n_stored-1 = 0 is vacuous, the key "
                                       "gates nothing; store >= 2 digits"
                                       % (sv,))}
            continue

        # working precision must cover the full stored string, both computed
        # strings, and both declared precisions — nothing may truncate.
        str_len = max(len(sv), max(len(str(p)) for p in pair))
        wp = max(n_stored, dlo, dhi, str_len) + 20
        with mp.workdps(wp):
            v = _parse(sv)
            a = _parse(pair[0])
            b = _parse(pair[1])
            pair_agree = _agree_digits(a, b, wp)
            hi_val = b
            matched = _agree_digits(hi_val, v, wp)
        ok_match = matched >= n_stored - 1
        ok_pair = pair_agree >= n_stored
        ok = ok_match and ok_pair
        reason = None
        if not ok:
            bits = []
            if not ok_match:
                bits.append("matched %.2f < stored-1 = %d"
                            % (matched, n_stored - 1))
            if not ok_pair:
                bits.append("two-precision agreement %.2f < stored = %d "
                            "(unstable — noise, not a result)"
                            % (pair_agree, n_stored))
            reason = "; ".join(bits)
        if ok:
            n_pass += 1
        per_key[key] = {"stored_digits": n_stored,
                        "matched_digits": round(matched, 2),
                        "pair_agree_digits": round(pair_agree, 2),
                        "pair_identical_strings": pair_identical,
                        "PASS": ok,
                        "reason": reason}

    extra = sorted(set(computed) - set(vendored))
    overall = "PASS" if (n_pass == len(vendored) and not missing) else "FAIL"
    return {
        "overall": overall,
        "rule": ("PASS iff matched >= stored_digits-1 vs FULL stored string "
                 "AND two computed precisions agree to >= stored_digits; "
                 "zero-vendored keys gate ABSOLUTELY (both pair members <= "
                 "scale_ref*10^-zero_digits); nonzero 1-digit vendored "
                 "strings hard-fail (un-gateable)"),
        "dps_pair": [dlo, dhi],
        "n_keys": len(vendored),
        "n_pass": n_pass,
        "per_key": per_key,
        "missing_keys": missing,
        "extra_keys": extra,
    }


def feed_gate_table(entries, digit_bar=30):
    """Per-literal FEED gate table.

    Pure formatting + thresholds: NOTHING is computed from values here —
    feed it MEASURED numbers (e.g. gate_strings per_key entries, or the
    run's own recompute log).  One row per literal:

        literal | stored digits | recomputed-vs-stored matched |
        two-precision agreement | verdict

    Verdict thresholds (the charter bar):
      PASS iff  stored_digits    >= digit_bar          (>= 30 genuine
                                                        digits, charter)
           AND  matched_digits   >= stored_digits - 1  (the FULL stored
                                                        string reproduced —
                                                        no partial credit)
           AND  pair_agree_digits >= stored_digits     (two-precision
                                                        stable — a value
                                                        that moves with
                                                        working precision
                                                        is noise)
    A literal with a missing measurement (None matched/pair) FAILs — an
    unmeasured literal is not a gated literal.

    Args:
      entries: iterable of dicts, each with keys
               'literal' (str), 'stored_digits' (int),
               'matched_digits' (float or None),
               'pair_agree_digits' (float or None),
               and optionally 'note' (appended to the row).
      digit_bar: minimum stored-digit claim (default 30 — the charter
               ">= 30 genuine held-out digits" bar).

    Returns (WRITES NOTHING):
      {'overall': 'PASS'|'FAIL', 'n_literals': int, 'n_pass': int,
       'digit_bar': int, 'rule': str,
       'per_literal': [{'literal','stored_digits','matched_digits',
                        'pair_agree_digits','PASS','reason'}, ...],
       'table': str}   # the formatted table, ready to print/log
    """
    digit_bar = int(digit_bar)
    per = []
    for e in entries:
        try:
            name = str(e["literal"])
            stored = int(e["stored_digits"])
            matched = e["matched_digits"]
            pair = e["pair_agree_digits"]
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("feed_gate_table: malformed entry %r (%s)"
                             % (e, exc))
        note = str(e.get("note", "") or "")
        reasons = []
        if matched is None or pair is None:
            reasons.append("unmeasured (matched/pair is None) — an "
                           "unmeasured literal is not a gated literal")
            matched_f = None if matched is None else float(matched)
            pair_f = None if pair is None else float(pair)
        else:
            matched_f = float(matched)
            pair_f = float(pair)
            if stored < digit_bar:
                reasons.append("stored %d < charter bar %d"
                               % (stored, digit_bar))
            if matched_f < stored - 1:
                reasons.append("matched %.2f < stored-1 = %d"
                               % (matched_f, stored - 1))
            if pair_f < stored:
                reasons.append("two-precision agreement %.2f < stored = %d "
                               "(noise, not a result)" % (pair_f, stored))
        ok = not reasons
        per.append({"literal": name,
                    "stored_digits": stored,
                    "matched_digits": matched_f,
                    "pair_agree_digits": pair_f,
                    "PASS": ok,
                    "reason": None if ok else "; ".join(reasons),
                    "note": note or None})
    n_pass = sum(1 for r in per if r["PASS"])
    overall = "PASS" if (per and n_pass == len(per)) else "FAIL"

    wname = max([len(r["literal"]) for r in per] + [len("literal")])
    fmt = "%-{w}s  %6s  %8s  %10s  %-7s%s".replace("{w}", str(wname))

    def _num(v):
        return "-" if v is None else "%.1f" % v

    lines = ["FEED gate (charter bar: stored >= %dd, matched >= stored-1, "
             "two-precision agree >= stored)" % digit_bar,
             fmt % ("literal", "stored", "matched", "pair-agree",
                    "verdict", "")]
    for r in per:
        tail = ""
        if r["note"]:
            tail = "  " + r["note"]
        if not r["PASS"]:
            tail += "  [%s]" % r["reason"]
        lines.append(fmt % (r["literal"], r["stored_digits"],
                            _num(r["matched_digits"]),
                            _num(r["pair_agree_digits"]),
                            "PASS" if r["PASS"] else "FAIL", tail))
    lines.append("FEED GATE: %s (%d/%d literals)"
                 % (overall, n_pass, len(per)))

    return {"overall": overall,
            "n_literals": len(per),
            "n_pass": n_pass,
            "digit_bar": digit_bar,
            "rule": ("PASS iff stored >= %d AND matched >= stored-1 AND "
                     "pair_agree >= stored; unmeasured literals FAIL"
                     % digit_bar),
            "per_literal": per,
            "table": "\n".join(lines)}


def write_gate_report(report, out_json):
    """Write a gate_strings report to out_json (the ONLY writer in this
    module).  mpf/mpc leftovers are stringified, never rounded silently."""
    with open(out_json, "w") as fh:
        json.dump(report, fh, indent=1, default=str)
        fh.write("\n")

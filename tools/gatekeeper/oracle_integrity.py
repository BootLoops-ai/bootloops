#!/usr/bin/env python3
"""
oracle_integrity.py — hash-based dedup auditor for oracle per-master-per-point
JSON dumps.

The bug this guards against (seen in practice): two ISP/numerator masters
mis-collected with the wrong propagator indices, so their oracle files were
*byte-identical* to two other masters'. A naive value-fit then "certified"
more masters than were genuinely distinct — the fit certifies its own
training data.

Detection strategy: parse each JSON, DROP the identity fields (tag/k/indices),
re-serialize the remaining numeric payload (s, t, leading, coeffs) with sorted
keys, SHA256 it. Any hash group containing >1 distinct (tag, indices) is a
DUPLICATE — same numerics under two different master labels.

Also flags:
  * empty / all-zero coeffs   (the eps_order footgun: amflow returned nothing
                               but the file was still written),
  * unreadable / non-JSON files.
Two dump forms are read — the tagged per-file form and the AMFlow out.json
form (one record per result entry); see "Dump forms" below.

Value strings.  Coefficient entries are parsed through amf_result.amf_ball,
so every form amflow-cpp's acb_to_json (Arb's arb_get_str) emits is read:
plain decimals ("0", "0.25", "-1.0e-5"), balls "[mid +/- rad]" (also
"mid +/- rad" without brackets and Arb's zero-midpoint "[+/- rad]"), and
defensively "[a, b]" intervals.  The |z| < tol zero test is applied to the
MIDPOINT with mpmath at a working precision set from the string length
(never float() on a 100-digit string).  Before this, a ball string failed
float(), the fallback stripped "+-0." and saw "[" -> the entry counted as
ZERO and a valid nonzero output was quarantined.

Consistent-with-zero.  A ball whose radius exceeds |mid| (rad > |mid| > tol)
does not certify a nonzero value, but it is not the eps_order footgun
either.  It is reported in its OWN field, `consistent_with_zero`
(one record per such entry: file/tag/order/part/mid/rad), never folded
into `suspicious_zero`.  Exit-code semantics: REPORT-ONLY by default —
listed and printed, verdict and exit code unchanged; `--strict-zero-balls`
(scan(..., strict_zero_balls=True)) promotes those files to the quarantine
set and to QUARANTINE_NEEDED.  Chosen this way because quarantine is
destructive (files move) and a wide radius is a precision statement about
the run, not evidence the file is a mis-write; the strict mode exists for
fits that must not see an uncertified zero.

Dump forms.  Two file shapes are read (see _records_of):
  A. the tagged per-file form — top-level identity keys (tag/k/indices, ...)
     beside "coeffs" (or "coefficients") = {"<order>": {"re": s, "im": s}};
     one record per file, behaviour unchanged;
  B. the AMFlow solve_integrals out.json form — no top-level coefficient key
     and "result" = [{"integral": {"family", "indices"}, "leading_order",
     "coefficients": [{"order": n, "value": {"re": s, "im": s}}, ...]}, ...];
     EVERY result entry is one record: tag "<file-tag>#<indices joined by
     commas>" (file-tag = the file's tag/k/master/name, else the entry's
     integral.family, else the file's basename stem), indices = the entry's
     indices (entry.indices, else entry.integral.indices), coeffs =
     {str(order): value} re-keyed from entry.coefficients — the shape the
     classifier reads — and the duplicate hash over that record (the entry
     minus its identity, coefficients re-keyed the same way).
  A file with neither shape, or an out.json whose result list is empty, is one
  record with empty coeffs: "empty coeffs", exactly as before.  n_files counts
  files; n_records counts records (form A: one per readable file; form B: one
  per result entry).  Before this, an out.json dump had no top-level key, read
  as {} and was quarantined "empty coeffs" while its balls were finite.

FAST: pure hashing; the only arbitrary-precision arithmetic is the
midpoint/radius comparison of the zero test.

CLI:
    python oracle_integrity.py DIR [DIR ...] \
        [--report integrity_report.json] [--quarantine] [--zero-tol 1e-80] \
        [--strict-zero-balls]
"""
from __future__ import annotations
import argparse
import glob
import hashlib
import json
import os
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict

import mpmath as mp

try:                                  # `from gatekeeper import ...` / -m
    from .amf_result import amf_ball
except ImportError:                   # `python oracle_integrity.py`
    from amf_result import amf_ball


# ---------------------------------------------------------------- data model
IDENTITY_KEYS = ("tag", "k", "indices", "master", "name")


@dataclass
class IntegrityReport:
    dirs: list
    n_files: int = 0
    n_records: int = 0   # classifier records: 1 per tagged file, 1 per out.json result entry
    ok_count: int = 0
    duplicates: list = field(default_factory=list)      # [{"hash","files","tags","indices"}]
    suspicious_zero: list = field(default_factory=list) # [{"file","tag","reason"}]
    unreadable: list = field(default_factory=list)      # [{"file","error"}]
    # balls with rad > |mid|: [{"file","tag","order","part","mid","rad"}]
    consistent_with_zero: list = field(default_factory=list)
    strict_zero_balls: bool = False   # True => consistent_with_zero quarantines
    verdict: str = "CLEAN"

    def quarantine_set(self):
        """All files that should be EXCLUDED from a fit (duplicates beyond
        the first per hash, plus suspicious-zero, plus unreadable; plus the
        consistent-with-zero files ONLY when strict_zero_balls is set)."""
        bad = set()
        for grp in self.duplicates:
            for f in grp["files"][1:]:       # keep the first
                bad.add(f)
        for s in self.suspicious_zero:
            bad.add(s["file"])
        for u in self.unreadable:
            bad.add(u["file"])
        if self.strict_zero_balls:
            for c in self.consistent_with_zero:
                bad.add(c["file"])
        return bad

    def to_json(self):
        return asdict(self)


# ---------------------------------------------------------------- helpers
def _tag_of(d):
    """Best-effort master identifier from a record."""
    for k in ("tag", "k", "master", "name"):
        if k in d:
            return str(d[k])
    return None


def _normalized_payload(d):
    """Return a canonical JSON string of the *numeric* content only
    (identity fields stripped)."""
    payload = {k: v for k, v in d.items() if k not in IDENTITY_KEYS}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _parse_value(s):
    """One printed value string -> (mid, rad, dps): mpf midpoint, mpf radius
    (None if exact), and the working precision they were read at.
    Ball/interval forms are split by amf_result.amf_ball; the numbers are
    read with mpmath at a working precision set from the string length, so
    a 110-digit midpoint is never truncated through float().  Callers that
    compare them do so inside `mp.workdps(dps)` (mpmath rounds abs() to the
    context precision)."""
    mid_s, rad_s = amf_ball(s)
    ndig = sum(c.isdigit() for c in mid_s) + (sum(c.isdigit() for c in rad_s)
                                              if rad_s else 0)
    dps = max(30, ndig + 10)
    with mp.workdps(dps):
        mid = mp.mpf(mid_s)
        rad = mp.mpf(rad_s) if rad_s is not None else None
    return mid, rad, dps


def _classify_coeffs(coeffs, tol):
    """coeffs: {"-3":{"re":"...","im":"..."}, ...}  ->  (all_zero, cwz)

    all_zero : True if every entry's MIDPOINT has |mid| <= tol (or is
               missing) — the eps_order footgun class (suspicious_zero).
    cwz      : [(order, part, mid_str, rad_str)] for every ball entry with
               |mid| > tol but rad > |mid| — consistent with zero, reported
               in its own category, never counted as zero here.
    Plain-float entries behave exactly as before (|z| > tol => nonzero)."""
    if not coeffs:
        return True, []
    all_zero, cwz = True, []
    tol = mp.mpf(tol)
    for order, v in coeffs.items():
        for part in ("re", "im"):
            s = str(v.get(part, "0")) if isinstance(v, dict) else "0"
            try:
                mid, rad, dps = _parse_value(s)
            except (ValueError, TypeError):
                # unparseable string with nonzero leading digits -> nonzero
                # (legacy fallback; kept for non-Arb producers)
                stripped = s.lstrip("+-0.[ ")
                if stripped and stripped[0] in "123456789":
                    all_zero = False
                continue
            with mp.workdps(dps):
                if abs(mid) > tol:
                    all_zero = False
                    if rad is not None and rad > abs(mid):
                        cwz.append((str(order), part, mp.nstr(mid, 20),
                                    mp.nstr(rad, 8)))
    return all_zero, cwz


def _is_effectively_zero(coeffs, tol):
    """True if every entry's midpoint parses to |z| <= tol (or is missing).
    Kept as the compatibility name; see _classify_coeffs."""
    return _classify_coeffs(coeffs, tol)[0]


def _iter_json_files(dirs):
    for d in dirs:
        if os.path.isfile(d) and d.endswith(".json"):
            yield d
            continue
        for f in sorted(glob.glob(os.path.join(d, "**", "*.json"),
                                  recursive=True)):
            yield f


# ------------------------------------------------ dump forms (see docstring)
ENTRY_IDENTITY_KEYS = ("integral", "coefficients")   # out.json per-entry identity + payload key


def _entry_indices(entry):
    """indices of one out.json result entry: entry['indices'] if present,
    else entry['integral']['indices'] (amflow-cpp's solve_integrals form)."""
    idx = entry.get("indices")
    if idx is None and isinstance(entry.get("integral"), dict):
        idx = entry["integral"].get("indices")
    return idx


def _entry_coeffs(entry):
    """result[i].coefficients ([{order, value}, ...]) -> {str(order): value},
    the shape _classify_coeffs reads; value kept as emitted ({re, im})."""
    out = {}
    for c in entry.get("coefficients") or []:
        if isinstance(c, dict) and "order" in c:
            out[str(c["order"])] = c.get("value")
    return out


def _entry_tag(file_tag, entry, idx, i, f):
    """'<file-tag>#<indices joined by commas>': file-tag = the file's own
    tag/k/master/name, else the entry's integral.family, else the file's
    basename stem; indices absent -> 'entry<i>'."""
    if file_tag is None and isinstance(entry.get("integral"), dict):
        file_tag = entry["integral"].get("family")
    if file_tag is None:
        file_tag = os.path.splitext(os.path.basename(f))[0]
    if isinstance(idx, (list, tuple)):
        suffix = ",".join(str(x) for x in idx)
    elif idx is not None:
        suffix = str(idx)
    else:
        suffix = "entry%d" % i
    return "%s#%s" % (file_tag, suffix)


def _records_of(d, f):
    """Yield (tag, indices_tuple, coeffs, payload_dict) — the classifier
    records of one parsed dump d (file f).

    Form B (AMFlow out.json): no top-level 'coeffs'/'coefficients' key and
    d['result'] a non-empty list whose dict entries carry 'coefficients'
    (an entry without them is a record with empty coeffs).  One record per
    entry; payload = the entry minus 'integral' and IDENTITY_KEYS with
    'coefficients' re-keyed to {str(order): value} under 'coeffs'.

    Otherwise form A (or an unreadable shape): ONE record, tag = _tag_of(d),
    indices = d['indices'], coeffs = the top-level key (or {}), payload = d
    — byte-for-byte the pre-out.json behaviour."""
    result = d.get("result")
    if ("coeffs" not in d and "coefficients" not in d
            and isinstance(result, list) and result
            and any(isinstance(e, dict) and "coefficients" in e
                    for e in result)):
        file_tag = _tag_of(d)
        for i, e in enumerate(result):
            e = e if isinstance(e, dict) else {}
            idx = _entry_indices(e)
            idx_t = tuple(idx) if isinstance(idx, list) else idx
            coeffs = _entry_coeffs(e)
            payload = {k: v for k, v in e.items()
                       if k not in IDENTITY_KEYS
                       and k not in ENTRY_IDENTITY_KEYS}
            payload["coeffs"] = coeffs
            yield _entry_tag(file_tag, e, idx, i, f), idx_t, coeffs, payload
        return
    idx = d.get("indices")
    idx_t = tuple(idx) if isinstance(idx, list) else idx
    coeffs = d.get("coeffs") or d.get("coefficients") or {}
    yield _tag_of(d), idx_t, coeffs, d


# ---------------------------------------------------------------- core
def scan(dirs, zero_tol: float = 1e-80,
         strict_zero_balls: bool = False) -> IntegrityReport:
    """Walk DIRS, hash normalized numeric payloads, group, flag duplicates
    and suspicious-zero records; list consistent-with-zero balls (rad >
    |mid|) in their own field — report-only unless strict_zero_balls, which
    promotes them to the quarantine set / QUARANTINE_NEEDED.
    Two dump forms are read (module docstring, _records_of): the tagged
    per-file form (one record per file) and the AMFlow out.json form (one
    record per result entry, tagged '<file-tag>#<indices>').  n_files counts
    files, n_records the records classified.
    Returns an IntegrityReport."""
    if isinstance(dirs, (str, os.PathLike)):
        dirs = [dirs]
    rep = IntegrityReport(dirs=[str(d) for d in dirs],
                          strict_zero_balls=bool(strict_zero_balls))

    groups = defaultdict(list)   # hash -> [(file, tag, indices_tuple)]

    for f in _iter_json_files(dirs):
        rep.n_files += 1
        try:
            with open(f) as fh:
                d = json.load(fh)
        except Exception as e:                       # noqa: BLE001
            rep.unreadable.append({"file": f, "error": str(e)})
            continue

        tag = _tag_of(d)

        # tag/indices consistency (soft check: tag embeds an int that should
        # match k if both present)
        if "k" in d and tag is not None:
            try:
                if str(d["k"]) != tag and not tag.endswith(str(d["k"])):
                    # only warn — different conventions exist
                    pass
            except Exception:
                pass

        # one record per tagged file (form A); one per result entry of an
        # AMFlow out.json (form B) — see _records_of
        for tag, idx_t, coeffs, payload in _records_of(d, f):
            rep.n_records += 1

            # suspicious-zero / empty / consistent-with-zero balls
            if not coeffs:
                rep.suspicious_zero.append(
                    {"file": f, "tag": tag, "reason": "empty coeffs"})
            else:
                all_zero, cwz = _classify_coeffs(coeffs, zero_tol)
                if all_zero:
                    rep.suspicious_zero.append(
                        {"file": f, "tag": tag,
                         "reason": "all coeffs |.|<%g (eps_order footgun?)" % zero_tol})
                for order, part, mid_s, rad_s in cwz:
                    rep.consistent_with_zero.append(
                        {"file": f, "tag": tag, "order": order, "part": part,
                         "mid": mid_s, "rad": rad_s})

            h = hashlib.sha256(_normalized_payload(payload).encode()).hexdigest()
            groups[h].append((f, tag, idx_t))

    # classify groups
    for h, members in groups.items():
        identities = {(m[1], m[2]) for m in members}
        if len(identities) > 1:
            rep.duplicates.append({
                "hash": h,
                "files": [m[0] for m in members],
                "tags": sorted({str(m[1]) for m in members}),
                "indices": sorted({str(m[2]) for m in members}),
                "n": len(members),
            })
            continue
        # same identity appearing multiple times is harmless (re-runs) but
        # still note it if user passed the same dir twice — count once as ok
        rep.ok_count += len(members)

    rep.verdict = ("CLEAN"
                   if not rep.duplicates
                   and not rep.suspicious_zero
                   and not rep.unreadable
                   and not (rep.strict_zero_balls and rep.consistent_with_zero)
                   else "QUARANTINE_NEEDED")
    return rep


def quarantine(rep: IntegrityReport, verbose=True):
    """Move every file in rep.quarantine_set() into a sibling
    `<parent>_QUARANTINE/` directory (created if absent)."""
    moved = []
    for f in sorted(rep.quarantine_set()):
        if not os.path.exists(f):
            continue
        parent = os.path.dirname(os.path.abspath(f))
        qdir = parent.rstrip("/").rstrip(os.sep) + "_QUARANTINE"
        os.makedirs(qdir, exist_ok=True)
        dst = os.path.join(qdir, os.path.basename(f))
        shutil.move(f, dst)
        moved.append((f, dst))
        if verbose:
            print(f"  quarantined {f} -> {dst}")
    return moved


# ---------------------------------------------------------------- CLI
def _main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("dirs", nargs="+", help="oracle data directories")
    p.add_argument("--report", default="integrity_report.json",
                   help="output JSON path")
    p.add_argument("--quarantine", action="store_true",
                   help="move duplicate/zero files to <dir>_QUARANTINE/")
    p.add_argument("--zero-tol", type=float, default=1e-80,
                   help="threshold below which a coeff counts as zero")
    p.add_argument("--strict-zero-balls", action="store_true",
                   help="treat balls with radius > |midpoint| "
                        "(consistent_with_zero) as quarantine-grade; "
                        "default lists them without changing the verdict")
    a = p.parse_args(argv)

    rep = scan(a.dirs, zero_tol=a.zero_tol,
               strict_zero_balls=a.strict_zero_balls)

    print(f"[oracle_integrity] scanned {rep.n_files} files in {len(a.dirs)} dir(s)")
    print(f"  ok           : {rep.ok_count}")
    print(f"  duplicates   : {len(rep.duplicates)} group(s)  "
          f"({sum(g['n'] for g in rep.duplicates)} files)")
    print(f"  suspicious-0 : {len(rep.suspicious_zero)}")
    print(f"  consistent-0 : {len(rep.consistent_with_zero)} ball(s) with rad>|mid|  "
          + ("(STRICT: quarantine-grade)" if rep.strict_zero_balls
             else "(report-only; --strict-zero-balls to quarantine)"))
    print(f"  unreadable   : {len(rep.unreadable)}")
    print(f"  VERDICT      : {rep.verdict}")
    for c in rep.consistent_with_zero[:12]:
        print(f"    CWZ {c['tag']} order={c['order']} {c['part']}: "
              f"mid={c['mid']} rad={c['rad']}  {c['file']}")
    if len(rep.consistent_with_zero) > 12:
        print(f"    ... (+{len(rep.consistent_with_zero)-12} more)")
    for g in rep.duplicates:
        print(f"    DUP {g['hash'][:12]}  tags={g['tags']}  "
              f"indices={g['indices']}  n={g['n']}")
        for f in g["files"][:6]:
            print(f"        {f}")
        if g["n"] > 6:
            print(f"        ... (+{g['n']-6} more)")

    with open(a.report, "w") as fh:
        json.dump(rep.to_json(), fh, indent=1)
    print(f"  wrote {a.report}")

    if a.quarantine and rep.verdict != "CLEAN":
        print("[oracle_integrity] --quarantine: moving offending files")
        quarantine(rep)

    return 0 if rep.verdict == "CLEAN" else 2


if __name__ == "__main__":
    sys.exit(_main())

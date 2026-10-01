"""emitall.engine — load a claims spec, emit every value, compare, report.

The engine's contract (the emit-all law): every quoted headline is RE-EMITTED
from its banked receipt through the spec's field expression and compared
against the quoted string. Any discrepancy is a FINDING — listed, never
smoothed — and the process exits nonzero on any finding. Outward displays are
EMITTED, never transcribed by hand (the sequencing rule: "display values
are computed outward by the emitting script, never transcribed from a
nearest-rounding recount").

Statuses per claim row:
    MATCH            emitted value equals the quote under the claim's mode
    FINDING          it does not (also: outward-band violations, eval errors)
    EMITTED          no quote given — value emitted for the record only
    MISSING-RECEIPT  the claim touched a receipt whose file is absent

Receipt-level findings (independent of claims):
    missing-receipt  the file is absent (one finding per receipt; every claim
                     that touches it is additionally marked MISSING-RECEIPT)
    stale-receipt    mtime check failed: the receipt is OLDER than one of its
                     declared `newer_than` inputs, or older than
                     `max_age_days`. Claims still evaluate (the data exists);
                     the staleness itself is the finding.

Comparison modes:
    printed       (default) compare at the quote's own printed precision —
                  byte-compatible with the emitting scripts it audits:
                  numbers are formatted with exactly the quote's decimal
                  count (str(int(round(v))) when the quote has no decimal
                  point); strings are compared after stripping "," and a
                  trailing "%" from both sides. An explicit "decimals" key
                  overrides the inferred precision. NOTE Python bools count
                  as numeric here (deliberately); wrap booleans in
                  str(...) in the spec.
    exact         str(emitted) == str(quoted), no normalization at all.
    outward-band  the printed interval [quoted_lo, quoted_hi] must be the
                  OUTWARD rounding of the computed interval at the stated
                  granularity: quoted_lo == floor(lo/g)*g and
                  quoted_hi == ceil(hi/g)*g, exact Fraction arithmetic.
                  Violations are graded: "inward-edge" when a printed edge
                  cuts into the computed interval (the register violation —
                  the 3-of-6 inward-edge class), "outward-loose" when the
                  band encloses but is wider than floor/ceil.

Spec format (JSON; YAML accepted when PyYAML is installed):
    {
      "campaign":     "...",
      "register":     "free-text register statement",
      "quote_source": "which document the quoted strings come from",
      "base":         "/abs/path prefix for relative receipt paths",
      "receipts": {
        "alias": {"path": "runs/x.json",
                  "format": "json"|"csv",         # default: by extension
                  "newer_than": ["other/path"],   # mtime ordering to enforce
                  "max_age_days": 30}             # optional absolute bound
      },
      "defaults": {"mode": "printed"},
      "claims": [
        {"section": "S", "label": "L",
         "lets": {"name": "expr", ...},           # evaluated in order
         "expr": "f('alias', '/ptr')",            # single value
           — or —
         "format": "{a}/{b}", "exprs": {"a": "...", "b": "..."},
           — or (outward-band) —
         "mode": "outward-band",
         "lo_expr": "...", "hi_expr": "...",
         "quoted_lo": 66, "quoted_hi": 100,
         "granularity": 1, "scale": "1e-6",
         "quoted": "3182",                        # omit => EMITTED-only row
         "decimals": 1,                           # optional printed override
         "note": "static rider text",
         "note_format": "...{x}...", "note_exprs": {"x": "..."}}
      ]
    }
"""
import csv
import json
import math
import os
import subprocess
from fractions import Fraction

from . import expr as E


class SpecError(ValueError):
    """The spec file itself is malformed (not a data finding — fix the spec)."""


# ---------------------------------------------------------------- loading

def load_spec(path):
    with open(path) as fh:
        text = fh.read()
    if path.endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:
            raise SpecError("YAML spec given but PyYAML is not installed")
        spec = yaml.safe_load(text)
    else:
        spec = json.loads(text)
    if not isinstance(spec, dict) or "claims" not in spec or "receipts" not in spec:
        raise SpecError("spec must be a dict with 'receipts' and 'claims'")
    return spec


def _receipt_path(base, path):
    return path if os.path.isabs(path) else os.path.join(base, path)


def load_receipts(spec, spec_dir):
    """Returns (env, receipt_report_rows, receipt_findings)."""
    base = spec.get("base", spec_dir)
    loaded, rep, findings = {}, [], []
    for alias, cfg in spec["receipts"].items():
        if isinstance(cfg, str):
            cfg = {"path": cfg}
        p = _receipt_path(base, cfg["path"])
        row = {"alias": alias, "path": p}
        if not os.path.exists(p):
            loaded[alias] = E.Env.MISSING
            row["status"] = "MISSING"
            findings.append({"kind": "missing-receipt", "receipt": alias,
                             "path": p, "note": "receipt file absent"})
            rep.append(row)
            continue
        mt = os.path.getmtime(p)
        row["mtime_utc"] = _iso_utc(mt)
        fmt = cfg.get("format") or ("csv" if p.endswith(".csv") else "json")
        if fmt == "json":
            with open(p) as fh:
                loaded[alias] = json.load(fh)
        elif fmt == "csv":
            with open(p) as fh:
                loaded[alias] = list(csv.DictReader(fh))
        else:
            raise SpecError(f"receipt {alias!r}: unknown format {fmt!r}")
        row["status"] = "loaded"
        for dep in cfg.get("newer_than", []):
            dp = _receipt_path(base, dep)
            if not os.path.exists(dp):
                findings.append({"kind": "stale-receipt", "receipt": alias,
                                 "path": p,
                                 "note": f"declared input {dp} absent — "
                                         "ordering unverifiable"})
                row["status"] = "STALE?"
            elif mt < os.path.getmtime(dp):
                findings.append({"kind": "stale-receipt", "receipt": alias,
                                 "path": p,
                                 "note": f"older than declared input {dp} "
                                         f"({row['mtime_utc']} < "
                                         f"{_iso_utc(os.path.getmtime(dp))})"})
                row["status"] = "STALE"
        max_age = cfg.get("max_age_days")
        if max_age is not None:
            import time
            if (time.time() - mt) > float(max_age) * 86400:
                findings.append({"kind": "stale-receipt", "receipt": alias,
                                 "path": p,
                                 "note": f"older than max_age_days={max_age}"})
                row["status"] = "STALE"
        rep.append(row)
    return E.Env(loaded), rep, findings


def _iso_utc(t):
    import datetime
    return datetime.datetime.fromtimestamp(
        t, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ------------------------------------------------------------- comparison

def fmt_like(quote, value, decimals=None):
    """Format value with the same decimal places as the quoted string.

    Byte-compatible with the fmt_like semantics of the emitting scripts it
    audits; an explicit `decimals` overrides the inference.
    """
    if decimals is not None:
        if decimals > 0:
            return f"{value:.{decimals}f}"
        return str(int(round(value)))
    q = str(quote).replace(",", "").rstrip("%")
    if "." in q:
        dec = len(q.split(".")[1])
        return f"{value:.{dec}f}"
    return str(int(round(value)))


def compare_printed(emitted, quoted, decimals=None):
    qs = str(quoted).replace(",", "").rstrip("%")
    if isinstance(emitted, Fraction):
        # frac() values print at the quote's register like any number;
        # exactness lives in `exact` mode (str(Fraction) == "p/q")
        es = fmt_like(quoted, float(emitted), decimals)
    elif isinstance(emitted, (int, float)):
        es = fmt_like(quoted, emitted, decimals)
    else:
        es = str(emitted).replace(",", "").rstrip("%")
    return es == qs, es


def _frac(x, what):
    """Exact Fraction from a JSON number or decimal/scientific string."""
    try:
        if isinstance(x, float):
            # JSON floats like 1e-6 arrive as binary doubles; str() recovers
            # the shortest decimal, which is what the spec author wrote.
            return Fraction(str(x))
        return Fraction(str(x))
    except (ValueError, ZeroDivisionError) as e:
        raise SpecError(f"outward-band {what} not a number: {x!r} ({e})") from None


def _gran_decimals(g):
    """Decimal places needed to print multiples of granularity g (cap 12)."""
    for d in range(0, 13):
        if (g * 10 ** d).denominator == 1:
            return d
    return 12


def _fdisp(x, dp):
    return f"{float(x):.{dp}f}" if dp > 0 else str(int(x))


def compare_outward_band(lo, hi, quoted_lo, quoted_hi, granularity=1, scale=1):
    """Enforce the outward rounding law: floor lo / ceil hi at granularity.

    Returns (ok, emitted_display, violation_note_or_None). Exact Fraction
    arithmetic throughout — no float participates in the accept/reject.
    """
    g = _frac(granularity, "granularity")
    if g <= 0:
        raise SpecError("outward-band granularity must be positive")
    s = _frac(scale, "scale")
    lo_v = _frac(lo, "computed lo") * s
    hi_v = _frac(hi, "computed hi") * s
    if lo_v > hi_v:
        raise SpecError(f"outward-band computed lo > hi ({lo_v} > {hi_v})")
    q_lo = _frac(quoted_lo, "quoted_lo")
    q_hi = _frac(quoted_hi, "quoted_hi")
    flo = math.floor(lo_v / g) * g
    chi = math.ceil(hi_v / g) * g
    dp = _gran_decimals(g)
    emitted = (f"outward [{_fdisp(flo, dp)}-{_fdisp(chi, dp)}] at g={g} "
               f"(computed lo={float(lo_v)!r}, hi={float(hi_v)!r})")
    if q_lo == flo and q_hi == chi:
        return True, emitted, None
    inward = []
    if q_lo > lo_v:
        inward.append(f"lo edge {_fdisp(q_lo, dp)} > computed lo {float(lo_v)!r}")
    if q_hi < hi_v:
        inward.append(f"hi edge {_fdisp(q_hi, dp)} < computed hi {float(hi_v)!r}")
    if inward:
        return False, emitted, ("inward-edge: " + "; ".join(inward)
                                + " — printed band cuts into the computed "
                                  "interval (rounding-law violation)")
    return False, emitted, ("outward-loose: band encloses the interval but is "
                            f"not the outward rounding "
                            f"[{_fdisp(flo, dp)}-{_fdisp(chi, dp)}] at g={g}")


# ------------------------------------------------------------ claim eval

def _eval_lets(claim, env):
    variables = {}
    for name, ex in (claim.get("lets") or {}).items():
        variables[name] = E.evaluate(ex, env, variables)
    return variables


def _eval_emitted(claim, env, variables):
    if "expr" in claim:
        return E.evaluate(claim["expr"], env, variables)
    if "format" in claim:
        parts = {k: E.evaluate(ex, env, variables)
                 for k, ex in (claim.get("exprs") or {}).items()}
        try:
            return claim["format"].format(**parts)
        except (KeyError, ValueError) as e:
            raise SpecError(f"format template failed: {e}") from None
    raise SpecError(f"claim {claim.get('label')!r} has neither expr nor format")


def _eval_note(claim, env, variables):
    if "note_format" in claim:
        parts = {k: E.evaluate(ex, env, variables)
                 for k, ex in (claim.get("note_exprs") or {}).items()}
        return claim["note_format"].format(**parts)
    return claim.get("note", "")


def run_spec(spec_path, out_path=None, quiet=False):
    """Evaluate every claim; write the report; return (report, exit_code)."""
    spec = load_spec(spec_path)
    spec_dir = os.path.dirname(os.path.abspath(spec_path))
    env, receipt_rows, findings = load_receipts(spec, spec_dir)
    default_mode = (spec.get("defaults") or {}).get("mode", "printed")

    rows = []
    for claim in spec["claims"]:
        section = claim.get("section", "")
        label = claim.get("label")
        if not label:
            raise SpecError("every claim needs a label")
        mode = claim.get("mode", default_mode)
        quoted = claim.get("quoted")
        note = claim.get("note", "")
        row = {"section": section, "label": label, "mode": mode,
               "quoted": quoted, "note": note}
        try:
            variables = _eval_lets(claim, env)
            if mode == "outward-band":
                lo = E.evaluate(claim["lo_expr"], env, variables)
                hi = E.evaluate(claim["hi_expr"], env, variables)
                ok, emitted, violation = compare_outward_band(
                    lo, hi, claim["quoted_lo"], claim["quoted_hi"],
                    claim.get("granularity", 1), claim.get("scale", 1))
                row["quoted"] = quoted = claim.get(
                    "quoted", f"[{claim['quoted_lo']}-{claim['quoted_hi']}]")
                row["emitted"] = emitted
                row["note"] = _eval_note(claim, env, variables)
                if ok:
                    row["status"] = "MATCH"
                else:
                    row["status"] = "FINDING"
                    findings.append({"kind": "outward-band", "section": section,
                                     "label": label, "emitted": emitted,
                                     "quoted": quoted, "note": violation})
                    row["note"] = (violation + (" | " + row["note"]
                                                if row["note"] else ""))
                rows.append(row)
                continue

            emitted = _eval_emitted(claim, env, variables)
            row["emitted"] = emitted
            row["note"] = _eval_note(claim, env, variables)
            if quoted is None:
                row["status"] = "EMITTED"
            else:
                if mode == "printed":
                    ok, printed = compare_printed(emitted, quoted,
                                                  claim.get("decimals"))
                    row["emitted_printed"] = printed
                elif mode == "exact":
                    ok = str(emitted) == str(quoted)
                else:
                    raise SpecError(f"unknown mode {mode!r}")
                row["status"] = "MATCH" if ok else "FINDING"
                if not ok:
                    findings.append({"kind": "mismatch", "section": section,
                                     "label": label, "emitted": emitted,
                                     "quoted": quoted, "note": row["note"]})
        except E.ReceiptMissing as e:
            row["status"] = "MISSING-RECEIPT"
            row["emitted"] = None
            row["note"] = str(e)
        except (E.ExprError, SpecError, TypeError, ValueError,
                ZeroDivisionError) as e:
            row["status"] = "FINDING"
            row["emitted"] = None
            row["note"] = f"eval-error: {e}"
            findings.append({"kind": "eval-error", "section": section,
                             "label": label, "quoted": quoted,
                             "note": str(e)})
        rows.append(row)

    n_match = sum(1 for r in rows if r["status"] == "MATCH")
    report = {
        "run": "emitall",
        "spec": os.path.abspath(spec_path),
        "campaign": spec.get("campaign", ""),
        "register": spec.get("register", ""),
        "quote_source": spec.get("quote_source", ""),
        "date_utc": subprocess.run(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"],
                                   capture_output=True, text=True).stdout.strip(),
        "receipts": receipt_rows,
        "headlines_emitted": len(rows),
        "matches": n_match,
        "findings_count": len(findings),
        "FINDINGS": findings,
        "table": rows,
    }
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(report, fh, indent=1, default=str)
    if not quiet:
        print_report(report, out_path)
    return report, (0 if not findings else 1)


def print_report(report, out_path=None):
    for r in report["table"]:
        q = "" if r.get("quoted") is None else f"  [QUOTED: {r['quoted']}]"
        print(f"{r['status']:16s} {r['section']:12s} {r['label']}: "
              f"{r.get('emitted')}{q}")
        if r.get("note"):
            print(f"                 note: {r['note']}")
    print(f"\n{report['headlines_emitted']} headlines emitted; "
          f"{report['matches']} MATCH; {report['findings_count']} FINDINGS")
    for i, f_ in enumerate(report["FINDINGS"], 1):
        loc = f_.get("label") or f_.get("receipt") or "?"
        print(f"FINDING {i} [{f_['kind']}] {loc}: "
              f"emitted {f_.get('emitted')!r} vs quoted {f_.get('quoted')!r}"
              f" — {f_.get('note', '')}")
    if out_path:
        print("report ->", out_path)

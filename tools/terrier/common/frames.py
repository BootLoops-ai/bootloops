#!/usr/bin/env python3
"""frames.py — terrier common chassis: convention-registry template.

A frozen convention dictionary as data: named rows with evidence classes
[SRC]/[PIL]/[REF], executed-toggle tables, and tamper-evident freeze/load,
so a project records its sign/normalization conventions in one pinned table
instead of in comments. stdlib only.

Freeze semantics:
- a PINNED row means the recomputation stack MUST adopt the frozen choice;
  deviation is a bug on OUR side;
- a float match at published precision is a convention LOCK, not a
  confirmation of any published value;
- pin by number, not by symbol; toggles record the executed one-code-both-
  numbers evidence.
"""
import hashlib, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from receipts import atomic_write_json, utc_now  # noqa: E402

STATUSES = ("PINNED-V", "PINNED-N", "PINNED-V+N", "OURS", "PIN-REQ")
EVIDENCE_CLASSES = ("SRC", "PIL", "REF")   # source-verbatim / pilot / referee


def _rows_sha(rows):
    blob = json.dumps(rows, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


class Registry:
    """Named convention rows with evidence class + toggles; freezable."""

    def __init__(self, name):
        self.name, self.rows, self.frozen_sha = name, {}, None

    def add(self, row_id, statement, status, group="", evidence=(), toggles=()):
        assert self.frozen_sha is None, "registry is FROZEN — no new rows"
        assert status in STATUSES, "status %r not in %s" % (status, STATUSES)
        for ev in evidence:
            assert ev.get("class") in EVIDENCE_CLASSES, "bad evidence class"
            assert ev.get("anchor"), "evidence needs a verbatim anchor"
        self.rows[row_id] = {"statement": statement, "status": status,
                             "group": group, "evidence": list(evidence),
                             "toggles": list(toggles)}

    def add_toggle(self, row_id, setting, observed, matches):
        """Executed toggle: one code, each setting, what it matches."""
        assert self.frozen_sha is None, "registry is FROZEN"
        self.rows[row_id]["toggles"].append(
            {"setting": setting, "observed": str(observed), "matches": matches})

    def unpinned(self):
        return sorted(i for i, r in self.rows.items() if r["status"] == "PIN-REQ")

    def require_pinned(self, row_ids):
        """Gate: every row a computation relies on must exist and not be
        PIN-REQ. Raises AssertionError naming the offenders (fail-closed)."""
        bad = [i for i in row_ids if i not in self.rows]
        bad += [i for i in row_ids
                if i in self.rows and self.rows[i]["status"] == "PIN-REQ"]
        assert not bad, "unpinned/missing convention rows: %s" % sorted(set(bad))

    def freeze(self, path):
        """Write the frozen registry (atomic) with a tamper-evident rows-sha.
        PIN-REQ rows freeze as OPEN items (listed in the header, §4 pattern)."""
        self.frozen_sha = _rows_sha(self.rows)
        atomic_write_json(path, {
            "schema": "terrier-convention-registry-v1", "name": self.name,
            "frozen_utc": utc_now(), "sha256_rows": self.frozen_sha,
            "n_rows": len(self.rows), "open_pin_req": self.unpinned(),
            "rows": self.rows})
        return self.frozen_sha

    @classmethod
    def load(cls, path):
        """Load + verify the tamper-evident sha; returns a FROZEN registry."""
        d = json.load(open(path))
        assert d.get("schema") == "terrier-convention-registry-v1", "bad schema"
        got = _rows_sha(d["rows"])
        assert got == d["sha256_rows"], (
            "registry TAMPERED: rows sha %s != frozen %s" % (got[:12],
                                                             d["sha256_rows"][:12]))
        reg = cls(d["name"])
        reg.rows, reg.frozen_sha = d["rows"], d["sha256_rows"]
        return reg

    def to_markdown(self):
        """Markdown convention table — the next project's starting page."""
        out = ["# %s — convention registry (%d rows, %d PIN-REQ open)" %
               (self.name, len(self.rows), len(self.unpinned())),
               "", "| row | group | frozen convention | status | evidence |",
               "|---|---|---|---|---|"]
        for i in sorted(self.rows):
            r = self.rows[i]
            ev = "; ".join("[%s] %s" % (e["class"], e["anchor"])
                           for e in r["evidence"]) or "OPEN"
            out.append("| %s | %s | %s | %s | %s |" %
                       (i, r["group"], r["statement"], r["status"], ev))
        for i in sorted(self.rows):
            for t in self.rows[i]["toggles"]:
                out.append("- toggle %s[%s]: observed %s -> matches %s" %
                           (i, t["setting"], t["observed"], t["matches"]))
        return "\n".join(out) + "\n"

#!/usr/bin/env python3
"""bank.py — the interface table T (mmap, append-only, demand-scoped,
content-hash keyed).

Banks reduced interface rows per (family, cell, node) so repeat grid points
are lookups. Binary layout, append-only:

  data file  <base>.bin : records, each
      u64 lead_col | u32 n_terms | n_terms * (u64 col, u64 val)
  index file <base>.idx.json : {"meta": {...}, "keys": {key: {"off": int,
      "n_rows": int, "len": int, "sha256": hex}}}   (rewritten per append —
      v0 scale; the .bin is never rewritten, only appended)

Key = "family|cell|node". Demand scoping is the CALLER's contract: only rows
for demanded integrals are appended (STRATA Pass-1 census decides demand).
Every deliverable read back is hash-verified against the index (mmap read).
"""
import hashlib
import json
import mmap
import os
import struct

_REC_HDR = struct.Struct("<QI")
_TERM = struct.Struct("<QQ")


class InterfaceTable:
    def __init__(self, base, meta=None):
        self.bin_path = base + ".bin"
        self.idx_path = base + ".idx.json"
        if os.path.exists(self.idx_path):
            with open(self.idx_path) as fh:
                self.idx = json.load(fh)
            if meta:
                assert self.idx["meta"] == meta, "meta mismatch on reopen"
        else:
            self.idx = {"meta": meta or {}, "keys": {}}
            open(self.bin_path, "ab").close()
            self._flush_idx()

    # ------------------------------------------------------------- internals
    def _flush_idx(self):
        tmp = self.idx_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(self.idx, fh)
        os.replace(tmp, self.idx_path)

    @staticmethod
    def _encode(rows):
        parts = []
        for lead, terms in rows:
            parts.append(_REC_HDR.pack(int(lead), len(terms)))
            for c, v in terms:
                parts.append(_TERM.pack(int(c), int(v)))
        return b"".join(parts)

    # ------------------------------------------------------------------- API
    @staticmethod
    def key(family, cell, node):
        return f"{family}|{cell}|{node}"

    def append_rows(self, cell_key, node_key, rows, family=None):
        """rows: [(lead_col, [(col, val), ...]), ...]; returns the key.
        Append-only: appending an existing key is REFUSED (immutability —
        a changed reduction result must be a new node/cell key)."""
        k = self.key(family or self.idx["meta"].get("family", "fam"),
                     cell_key, node_key)
        assert k not in self.idx["keys"], f"key exists (append-only): {k}"
        blob = self._encode(rows)
        with open(self.bin_path, "ab") as fh:
            off = fh.tell()
            fh.write(blob)
        self.idx["keys"][k] = {"off": off, "len": len(blob), "n_rows": len(rows),
                               "sha256": hashlib.sha256(blob).hexdigest()}
        self._flush_idx()
        return k

    def get(self, cell_key, node_key, family=None, verify=True):
        k = self.key(family or self.idx["meta"].get("family", "fam"),
                     cell_key, node_key)
        ent = self.idx["keys"][k]
        with open(self.bin_path, "rb") as fh:
            mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
            blob = mm[ent["off"]:ent["off"] + ent["len"]]
            mm.close()
        if verify:
            assert hashlib.sha256(blob).hexdigest() == ent["sha256"], \
                f"T corruption at {k}"
        rows, pos = [], 0
        for _ in range(ent["n_rows"]):
            lead, n = _REC_HDR.unpack_from(blob, pos)
            pos += _REC_HDR.size
            terms = []
            for _ in range(n):
                c, v = _TERM.unpack_from(blob, pos)
                pos += _TERM.size
                terms.append((c, v))
            rows.append((lead, terms))
        assert pos == ent["len"]
        return rows

    def has(self, cell_key, node_key, family=None):
        return self.key(family or self.idx["meta"].get("family", "fam"),
                        cell_key, node_key) in self.idx["keys"]

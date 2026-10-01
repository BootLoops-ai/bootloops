#!/usr/bin/env python3
"""sweepchassis.py — generic receipts / blind-controls / void-on-miss / halt / resume
sweep chassis (stdlib only), for ANY (entry stream, predicate chain). The hash
laws below are deterministic and byte-stable; selftest gate SC7 checks them
against a reference shard receipt.

Deterministic laws:
  CTRL-ID-V1   anonymous check_ids (the evaluator cannot tell control from witness)
  CTRL-POS-V1  hash-determined blind control placement within a batch
  receipts     atomic tmp+fsync+rename JSON; a shard is done iff its receipt parses
               COMPLETE with matching code_sha256 + stream_sha256 (resume law)
  void-on-miss any missed control -> VOID_CONTROL_FAIL, dispatch halts, rc=8;
               any per-entry integrity mismatch -> VOID_INTEGRITY, rc=8
  HALT law     survivor candidate -> evidence frozen FIRST (candidate JSON), HALT
               marker + a line appended to attention.jsonl under the sweep root,
               shard quarantined (never a COMPLETE receipt), all dispatch
               refused (rc=9) until clear_halt(); the quarantined shard re-runs
               on resume
  replay       re-run a completed shard -> byte-identical receipt under a pinned clock
  single-writer lock with stale-pid takeover; NO RNG anywhere.

Interface (complete toy example in selftest_sweepchassis.py):
  evaluate(req) -> result dict          # predicate-chain evaluator (subprocess or fn)
  adjudicate(evaluate, sid, i, entry) -> {"verdict": str, optional "integrity":
      [msgs], "survivor": bool, "record": {...}}
  controls: [{"control_id","kind","req","expected",("tier"),("match_fn")}]; default
      expectation law = every expected k == result[k] (tier C).
Run law: ulimit -v 32505856, nice >= 5.
"""
import hashlib, json, os, tempfile, time


def sha_file(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def atomic_write(path, obj):
    """Atomic JSON write: tmp + fsync + rename."""
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True)
        f.flush()
        os.fsync(f.fileno())
    os.rename(tmp, path)


# DELEGATION LAW: the CTRL-ID-V1 / CTRL-POS-V1 hash primitives are OWNED by
# ../../common/controls.py (selftest SC7 pins them against a reference shard
# receipt). Delegate, never duplicate.
# hslot (per-control batch slot) has no common/ owner: kept here.
import importlib.util as _ilu
_cspec = _ilu.spec_from_file_location("_svc_controls", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "common", "controls.py"))
_svc = _ilu.module_from_spec(_cspec)
_cspec.loader.exec_module(_svc)
hid, hpos = _svc.hid, _svc.hpos


def hslot(shard_id, control_id, nbatch):
    """CTRL-POS-V1 per-control batch slot."""
    return int(hashlib.sha256(("CTRL-POS-V1:%s:%s" % (shard_id, control_id)).encode()
                              ).hexdigest(), 16) % (nbatch + 1)


class HaltSurvivor(Exception):
    def __init__(self, candidate_path, sid):
        self.candidate_path, self.sid = candidate_path, sid


def match_control(c, r):
    """Default tier-C expectation law: every expected key equals the result."""
    got = {k: r.get(k) for k in c["expected"]}
    return all(r.get(k) == v for k, v in c["expected"].items()), got


class Chassis:
    """One sweep instance rooted at `root` (receipts/HALT/lock/out live under it).
    code_files pin code_sha256; `clock` is injectable so a replayed shard is
    BYTE-identical (selftest SC3) — production uses the wall clock."""

    def __init__(self, root, code_files, label="", scope_line="", semantics="",
                 schema="sweepchassis-shard-receipt-v1", clock=time.time):
        self.root, self.label, self.schema = root, label, schema
        self.scope_line, self.semantics, self.clock = scope_line, semantics, clock
        self.code_files = list(code_files)
        for d in ("shards", os.path.join("shards", "quarantine"), "out"):
            os.makedirs(os.path.join(root, d), exist_ok=True)

    def code_sha(self):
        h = hashlib.sha256()
        for f in sorted(self.code_files):
            h.update(open(f, "rb").read())
        return h.hexdigest()

    def receipt_path(self, sid):
        return os.path.join(self.root, "shards", "receipt_%s.json" % sid)

    def receipt_ok(self, sid, stream_sha):
        """Resume law: COMPLETE + code sha + stream sha."""
        try:
            r = json.load(open(self.receipt_path(sid)))
            return (r.get("status") == "COMPLETE"
                    and r["pins"]["code_sha256"] == self.code_sha()
                    and r["pins"]["stream_sha256"] == stream_sha)
        except Exception:
            return False

    def halt_file(self):
        return os.path.join(self.root, "HALT.json")

    def dispatch_allowed(self):
        if os.path.exists(self.halt_file()):
            h = json.load(open(self.halt_file()))
            return False, ("HALT active since %s (reason=%s, shard=%s, candidate=%s)"
                           " — adjudicate + clear_halt before any dispatch"
                           % (h.get("utc"), h.get("reason"), h.get("shard_id"),
                              h.get("candidate_file")))
        return True, ""

    def clear_halt(self):
        """Adjudication acknowledged: archive candidate under out/adjudicated/ and
        remove the HALT marker."""
        if not os.path.exists(self.halt_file()):
            return False
        h = json.load(open(self.halt_file()))
        src = os.path.join(self.root, h["candidate_file"])
        dstdir = os.path.join(self.root, "out", "adjudicated")
        os.makedirs(dstdir, exist_ok=True)
        if os.path.exists(src):
            os.rename(src, os.path.join(dstdir, h["candidate_file"]))
        os.remove(self.halt_file())
        return True

    def emit_candidate(self, sid, payload, rehearsal=False):
        """HALT law: freeze evidence FIRST, then HALT marker + attention line
        (rehearsal artifacts are marked loudly synthetic)."""
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime(self.clock()))
        tag = "REHEARSAL_" if rehearsal else ""
        path = os.path.join(self.root, "CANDIDATE_%s%s_%s.json" % (tag, sid, stamp))
        atomic_write(path, dict(payload,
            schema=self.schema + "-candidate", REHEARSAL_SYNTHETIC=bool(rehearsal),
            label=("REHEARSAL-SYNTHETIC drill artifact — NOT a discovery" if rehearsal
                   else "!! SURVIVOR CANDIDATE — triple-verify before any claim"),
            shard_id=sid, utc=stamp, pins={"code_sha256": self.code_sha()}))
        atomic_write(self.halt_file(), {"reason": "SURVIVOR_CANDIDATE",
            "rehearsal": rehearsal, "candidate_file": os.path.basename(path),
            "shard_id": sid, "utc": stamp})
        with open(os.path.join(self.root, "attention.jsonl"), "a") as f:
            f.write(json.dumps({"utc": stamp, "kind": "SURVIVOR_CANDIDATE",
                                "rehearsal": rehearsal, "shard": sid,
                                "file": os.path.basename(path)}, sort_keys=True) + "\n")
            f.flush()
            os.fsync(f.fileno())
        return path

    def run_shard(self, evaluate, adjudicate, sid, stream_tag, stream_sha, batch,
                  controls, pins_extra=None, rehearsal=False):
        """One shard over `batch` with blind control injection (CTRL-POS-V1 slots,
        CTRL-ID-V1 ids). Raises HaltSurvivor on a survivor (evidence frozen +
        quarantined FIRST); otherwise writes + returns the atomic receipt."""
        t0 = self.clock()
        inj = sorted(controls, key=lambda c: hid(sid, c["control_id"]))
        slots = {}
        for c in inj:
            slots.setdefault(hslot(sid, c["control_id"], len(batch)), []).append(c)
        ctrl_res, miss, integ, records, verdicts = [], [], [], [], {}
        for i in range(len(batch) + 1):
            for c in slots.get(i, []):
                cid = hid(sid, c["control_id"])
                r = evaluate(dict(c["req"], check_id=cid))
                mf = c.get("match_fn")
                ok, got = mf(c, r) if mf else match_control(c, r)
                ctrl_res.append({"control_id": c["control_id"], "check_id": cid,
                                 "pass": ok, "tier": c.get("tier", "C"), "got": got})
                if not ok:
                    miss.append(c["control_id"])
            if i >= len(batch):
                break
            a = adjudicate(evaluate, sid, i, batch[i])
            verdicts[a["verdict"]] = verdicts.get(a["verdict"], 0) + 1
            integ += a.get("integrity", [])
            records.append(a.get("record", {"i": i, "verdict": a["verdict"]}))
            if a.get("survivor"):
                cand = self.emit_candidate(sid, {"stream": stream_tag, "index": i,
                    "entry": batch[i], "adjudication": a.get("record")}, rehearsal)
                atomic_write(os.path.join(self.root, "shards", "quarantine",
                                          "receipt_%s.json" % sid),
                             {"schema": self.schema, "shard_id": sid,
                              "status": "QUARANTINED_SURVIVOR"
                                        + ("_REHEARSAL" if rehearsal else ""),
                              "REHEARSAL_SYNTHETIC": rehearsal,
                              "candidate_file": os.path.basename(cand),
                              "pins": {"code_sha256": self.code_sha()}})
                raise HaltSurvivor(cand, sid)
        status = "COMPLETE"
        if miss:
            status = "VOID_CONTROL_FAIL"
        if integ:
            status = "VOID_INTEGRITY"
        receipt = {"schema": self.schema, "label": self.label,
                   "scope_line": self.scope_line, "survivor_semantics": self.semantics,
                   "shard_id": sid, "stream": stream_tag, "n_entries": len(batch),
                   "pins": dict(pins_extra or {}, code_sha256=self.code_sha(),
                                stream_sha256=stream_sha),
                   "verdicts": verdicts, "records": records,
                   "integrity": {"mismatches": integ, "n": len(integ)},
                   "controls": {"tier": ctrl_res, "pass": not miss, "missed": miss,
                                "injection_pos_hash": hpos(sid, len(inj))},
                   "rng": "none", "status": status,
                   "timing": {"wall_s": round(self.clock() - t0, 2),
                              "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                           time.gmtime(t0))}}
        atomic_write(self.receipt_path(sid), receipt)
        return receipt

    def run(self, evaluate, adjudicate, streams, controls, batch_size=1000, **kw):
        """streams = [(tag, entries, stream_sha)]. Returns (rc, n_run, n_skipped):
        rc 0 clean, 8 void-on-miss/integrity (dispatch halted), 9 survivor-halt or
        halt-marker active (dispatch refused)."""
        ok, why = self.dispatch_allowed()
        if not ok:
            print("!! DISPATCH REFUSED: %s" % why)
            return 9, 0, 0
        done = skipped = 0
        for tag, entries, ssha in streams:
            nb = (len(entries) + batch_size - 1) // batch_size
            for ib in range(nb):
                sid = "%sB%03d" % (tag, ib)
                if self.receipt_ok(sid, ssha):
                    skipped += 1
                    continue
                batch = entries[ib * batch_size:(ib + 1) * batch_size]
                try:
                    r = self.run_shard(evaluate, adjudicate, sid, tag, ssha, batch,
                                       controls, **kw)
                except HaltSurvivor as h:
                    print("!! SURVIVOR CANDIDATE %s — halted; evidence frozen at %s"
                          % (h.sid, os.path.basename(h.candidate_path)))
                    return 9, done, skipped
                done += 1
                if r["status"] != "COMPLETE":
                    print("!! VOID (%s) shard %s — halting dispatch (rc=8)"
                          % (r["status"], sid))
                    return 8, done, skipped
        return 0, done, skipped

    def replay(self, evaluate, adjudicate, sid, stream_tag, stream_sha, batch,
               controls, **kw):
        """Determinism audit: re-run a completed shard, diff all
        but timing. Under a pinned clock the receipt FILE is byte-identical (SC3)."""
        old = json.load(open(self.receipt_path(sid)))
        new = self.run_shard(evaluate, adjudicate, sid, stream_tag, stream_sha,
                             batch, controls, **kw)
        a = {k: v for k, v in old.items() if k != "timing"}
        b = {k: v for k, v in new.items() if k != "timing"}
        return json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)

    def lock(self, name="chassis.lock"):
        return _Lock(os.path.join(self.root, name))


class _Lock:
    """Single-writer lock, stale-pid takeover, ownership-checked removal
    (never delete another process's lock)."""

    def __init__(self, path):
        self.path = path

    def _take(self):
        fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, str(os.getpid()).encode())
        os.close(fd)

    def __enter__(self):
        try:
            self._take()
        except FileExistsError:
            pid = open(self.path).read().strip()
            try:
                os.kill(int(pid), 0)
                raise RuntimeError("lock held by live pid %s — refusing double launch" % pid)
            except (ProcessLookupError, ValueError):
                os.remove(self.path)
                self._take()
        return self

    def __exit__(self, *exc):
        try:
            if open(self.path).read().strip() == str(os.getpid()):
                os.remove(self.path)
        except FileNotFoundError:
            pass
        return False

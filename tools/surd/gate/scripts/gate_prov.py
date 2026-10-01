"""Producer block for the run records ("receipts") written by the gate/scripts command-line tools: script path + sha256, launch line,
inputs with sha256, host, module shas (taken at import and again at write; a difference means a file was edited during the run),
package versions and a UTC stamp.  Also the output root shared with the slice scripts and a write-once JSON writer."""
import hashlib, os, socket, subprocess, sys, json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MODULES = ["oracle_cells.py", "oracle_core.py", "oracle.py", "gate_prov.py"]
# output root (same rule as slice/scripts/slice_prov.py): $SURD_WORK, default ./surd_work; nothing is created at import
WORK = os.path.abspath(os.environ.get("SURD_WORK") or os.path.join(os.getcwd(), "surd_work"))


def wdir(*parts):
    """directory WORK/<parts>, created on first use; refuses a WORK inside the package tree"""
    w = os.path.realpath(WORK); r = os.path.realpath(ROOT)
    if w == r or w.startswith(r + os.sep):
        sys.exit("surd: the output root %s is inside the package tree; set SURD_WORK to a scratch directory" % WORK)
    d = os.path.join(WORK, *parts); os.makedirs(d, exist_ok=True); return d


def sha(path):
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception as e:
        return "ERR:%s" % e


def stamp():
    try:
        return subprocess.check_output(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip()
    except Exception:
        import time; return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# module hashes are taken AT IMPORT (when the producing process starts and loads its code), not when the record is written;
# the producer block reports both when they differ (an in-place edit during a run).
IMPORT_STAMP = None
try:
    IMPORT_SHAS = {m: sha(os.path.join(HERE, m)) for m in sorted(os.listdir(HERE)) if m.endswith(".py") or m.endswith(".sh")}
    IMPORT_STAMP = stamp()
except Exception:
    IMPORT_SHAS = {}
_SCRIPT_SHA_AT_IMPORT = sha(os.path.abspath(sys.argv[0])) if sys.argv and sys.argv[0] and os.path.exists(sys.argv[0]) else None


def result_hash(response_path):
    """reproducible pin of a HyperFLINT response = sha256 of the canonical JSON of its 'result' array (the file itself embeds
    timing fields and is not bit-reproducible across re-runs)."""
    try:
        r = json.loads(open(response_path).read().strip().splitlines()[-1])
        return hashlib.sha256(json.dumps(r.get("result"), sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    except Exception as e:
        return "ERR:%s" % e


def dump_json(obj, path, indent=1):
    """atomic JSON write (tmp + rename); returns the sha256 of the written text"""
    data = json.dumps(obj, indent=indent, default=str)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp%d" % os.getpid()
    with open(tmp, "w") as fh: fh.write(data)
    os.replace(tmp, path)
    return hashlib.sha256(data.encode()).hexdigest()


def producer(inputs=(), extra=None, modules=None):
    script = os.path.abspath(sys.argv[0])
    mods = list(modules or MODULES)
    at_import = {m: IMPORT_SHAS.get(m, "not-hashed-at-import") for m in mods}
    at_write = {m: sha(os.path.join(HERE, m)) for m in mods}
    changed = {m: {"at_import": at_import[m], "at_write": at_write[m]} for m in mods if at_import[m] != at_write[m]}
    ssha_now = sha(script)
    blk = {"script": script, "script_sha256": _SCRIPT_SHA_AT_IMPORT or ssha_now, "launch": " ".join(sys.argv), "cwd": os.getcwd(),
           "host": socket.gethostname(), "python": sys.version.split()[0],
           "stamp_utc": stamp(), "import_stamp_utc": IMPORT_STAMP,
           "modules": at_import, "modules_hashed_at": "import",
           "inputs": {p: sha(p) for p in inputs},
           "inputs_result_hash": {p: result_hash(p) for p in inputs if str(p).endswith(".response.json")}}
    if changed or (_SCRIPT_SHA_AT_IMPORT and _SCRIPT_SHA_AT_IMPORT != ssha_now):
        blk["WARNING_modules_changed_during_run"] = changed
        if _SCRIPT_SHA_AT_IMPORT and _SCRIPT_SHA_AT_IMPORT != ssha_now: blk["WARNING_script_changed_during_run"] = {"at_import": _SCRIPT_SHA_AT_IMPORT, "at_write": ssha_now}
    for mod, key in (("flint", "python_flint"), ("sympy", "sympy"), ("mpmath", "mpmath")):
        try:
            blk[key] = __import__(mod).__version__
        except Exception:
            pass
    if extra:
        blk.update(extra)
    return blk

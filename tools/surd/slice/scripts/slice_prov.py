"""Run-record ("receipt") producer blocks, output/data roots and the pid log for the slice pipeline.  Engine modules are imported from
symbolic/scripts and gate/scripts; their sha256 are recorded in every receipt.

Roots (SLICE / ROOT locate CODE and package data only: scripts, gate/cells, gate/lr, data/):
  WORK = $SURD_WORK (default ./surd_work under the current directory): EVERY write of the slice scripts --
         receipts SLICE_*.json (+ superseded/), ckpt/, logs/, out/ -- and every read of this pipeline's own intermediate products
         (ckpt/<tag>/stage*.pkl, per-group receipts, SLICE_CLASSES.json, S0/eval caches).  Nothing is written into the package tree.
  DATA = $SURD_DATA (default: WORK): a directory holding assembled files out/*.json.gz and the receipts that the structure
         tools and tests/run_tests.py READ (same meaning as run_tests.py --data).  Reads of single receipts (SLICE_S5_LETTERS.json,
         SLICE_S0_ORACLE_t*.json, ckpt/cmsyz_cubics_slice.json) go through rpath()/rglob(): the WORK copy if this work tree produced one,
         else the DATA copy.
No directory is created at import (directories are made at first write)."""
import hashlib, os, socket, subprocess, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); SLICE = os.path.dirname(HERE); ROOT = os.path.dirname(SLICE)
WORK = os.path.abspath(os.environ.get("SURD_WORK") or os.path.join(os.getcwd(), "surd_work"))
DATA = os.path.abspath(os.environ.get("SURD_DATA") or WORK)
def _inside_package(p):
    p = os.path.realpath(p); r = os.path.realpath(ROOT); return p == r or p.startswith(r + os.sep)
def wdir(*parts):
    """directory WORK/<parts>, created on first use (write time), never at import; refuses a WORK inside the package tree"""
    if _inside_package(WORK): sys.exit("surd: the output root %s is inside the package tree; set SURD_WORK to a scratch directory" % WORK)
    d = os.path.join(WORK, *parts); os.makedirs(d, exist_ok=True); return d
def rpath(*parts):
    """READ path of a receipt or assembled file: the WORK copy if it exists (produced in this work tree), else the DATA copy"""
    p = os.path.join(WORK, *parts)
    return p if (WORK == DATA or os.path.exists(p)) else os.path.join(DATA, *parts)
def rglob(pattern):
    """sorted glob of receipts: under WORK if any match there, else under DATA"""
    r = sorted(glob.glob(os.path.join(WORK, pattern)))
    return r if (r or WORK == DATA) else sorted(glob.glob(os.path.join(DATA, pattern)))
SYMB_S = os.path.join(ROOT, "symbolic", "scripts"); GATE_S = os.path.join(ROOT, "gate", "scripts")
for p in (HERE, SYMB_S, GATE_S):
    if p not in sys.path: sys.path.insert(0, p)
sys.path.insert(0, HERE)
def sha(path):
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""): h.update(chunk)
        return h.hexdigest()
    except Exception as e: return "ERR:%s" % e
def stamp():
    try: return subprocess.check_output(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip()
    except Exception:
        import time; return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
ENGINE = {}
for d, mods in ((SYMB_S, ("fibr.py", "alphabet.py")),
                (GATE_S, ("oracle.py", "oracle_cells.py", "oracle_core.py", "hf_piece.py", "onefold_arb.py", "hlog_eval.py", "piece_select.py", "gate_prov.py"))):
    for m in mods: ENGINE[os.path.relpath(os.path.join(d, m), ROOT)] = sha(os.path.join(d, m))
def own_shas(): return {m: sha(os.path.join(HERE, m)) for m in sorted(os.listdir(HERE)) if m.endswith(".py") or m.endswith(".sh")}
IMPORT_SHAS = own_shas(); IMPORT_STAMP = stamp()
def producer(inputs=(), extra=None):
    import flint, sympy, mpmath
    now = own_shas()
    p = {"script": os.path.abspath(sys.argv[0]) if sys.argv and sys.argv[0] else None, "script_sha256": sha(os.path.abspath(sys.argv[0])) if sys.argv and sys.argv[0] and os.path.exists(sys.argv[0]) else None,
         "launch": " ".join(sys.argv), "cwd": os.getcwd(), "host": socket.gethostname(), "pid": os.getpid(),
         "inputs": [{"path": os.path.abspath(x), "sha256": sha(x)} for x in inputs], "engine_shas": ENGINE,
         "module_shas_at_import(slice/scripts)": IMPORT_SHAS, "modules_changed_during_run": sorted(k for k in now if IMPORT_SHAS.get(k) != now[k]), "import_stamp_utc": IMPORT_STAMP,
         "versions": {"python": sys.version.split()[0], "python-flint": flint.__version__, "sympy": sympy.__version__, "mpmath": mpmath.__version__}, "stamp_utc": stamp()}
    if extra: p.update(extra)
    return p
def write_receipt(path, rec, inputs=(), extra=None):
    """write-once: an existing receipt with different content is archived to <dir>/superseded/<name>_<sha12>.json first; atomic rename"""
    rec = dict(rec); rec["producer"] = producer(inputs, extra)
    if _inside_package(os.path.dirname(os.path.abspath(path))): sys.exit("surd: refusing to write %s inside the package tree; set SURD_WORK to a scratch directory" % path)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)          # lazy: the WORK root may not exist yet
    if os.path.exists(path):
        old = sha(path)[:12]; sd = os.path.join(os.path.dirname(path), "superseded"); os.makedirs(sd, exist_ok=True)
        os.replace(path, os.path.join(sd, os.path.basename(path)[:-5] + "_" + old + ".json"))
        with open(os.path.join(sd, "LEDGER.txt"), "a") as f: f.write("%s %s superseded sha12=%s by pid %d (%s)\n" % (stamp(), os.path.basename(path), old, os.getpid(), " ".join(sys.argv)[:200]))
    tmp = path + ".tmp%d" % os.getpid()
    with open(tmp, "w") as f: json.dump(rec, f, indent=1, default=str)
    os.replace(tmp, path); return path
def log_pid(tag):
    with open(os.path.join(wdir("logs"), "pids.log"), "a") as f:          # WORK/logs, created on first use
        f.write("%s pid=%d tag=%s argv=%s\n" % (stamp(), os.getpid(), tag, " ".join(sys.argv)))

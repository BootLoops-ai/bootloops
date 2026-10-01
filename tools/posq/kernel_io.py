# kernel_io.py -- python driver <-> posq_kernel C file protocol.
# Serialization is EXACT: arb balls as (mid mantissa, mid exp, rad mantissa,
# rad exp) integer 4-tuples; the kernel reconstructs via arf_set_fmpz_2exp +
# arf_get_mag (radius rounds UP: sound).
import json, os, subprocess
from fractions import Fraction as Fr
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
KERNEL = f"{HERE}/posq_kernel"          # a local build beside the source, if present (see GUIDE)
KERNEL_SRC = f"{HERE}/posq_kernel.c"


class KernelUnavailable(RuntimeError):
    """POSQ-KERNEL-UNAVAILABLE: the C kernel neither runs nor builds on this
    machine. Named, fail-closed -- kernel entry points refuse, never guess."""


_RESOLVED = {}                          # one probe/build per process


def _runs(path):
    """True iff `path` executes on this machine: with no job file the kernel
    prints its usage line and exits 2. A binary linked against a FLINT
    shared-library version this machine lacks (or built for another CPU
    architecture) fails the exec/loader step instead and never gets there."""
    try:
        r = subprocess.run([path], capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return r.returncode == 2 and "usage" in (r.stderr + r.stdout)


def _build_dir():
    d = os.environ.get("POSQ_KERNEL_CACHE")
    if d:
        return d
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(
        os.path.expanduser("~"), ".cache")
    return os.path.join(base, "posq")


def _try_build(log):
    """Compile posq_kernel.c (the GUIDE build line) into the cache dir, keyed
    by source sha so a source change rebuilds. Returns the path or None."""
    import hashlib
    try:
        tag = hashlib.sha256(open(KERNEL_SRC, "rb").read()).hexdigest()[:16]
    except OSError as e:
        log.append(f"kernel source unreadable: {e}")
        return None
    out = os.path.join(_build_dir(), f"posq_kernel-{tag}")
    if _runs(out):
        return out                      # cached earlier build
    try:
        os.makedirs(_build_dir(), exist_ok=True)
    except OSError as e:
        log.append(f"build cache dir unwritable: {e}")
        return None
    ccs = [os.environ["CC"]] if os.environ.get("CC") else ["cc", "gcc"]
    tmp = f"{out}.tmp.{os.getpid()}"
    for cc in ccs:
        cmd = [cc, "-O2", "-o", tmp, KERNEL_SRC,
               "-lflint", "-lmpfr", "-lgmp", "-lm"]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=300)
        except OSError:
            log.append(f"{cc}: compiler not found")
            continue
        except subprocess.TimeoutExpired:
            log.append(f"{cc}: build timed out")
            continue
        if r.returncode == 0:
            os.replace(tmp, out)        # atomic vs a concurrent builder
            if _runs(out):
                return out
            log.append(f"{cc}: built but does not execute")
        else:
            log.append(f"{cc}: rc={r.returncode} {r.stderr.strip()[:200]}"
                       .rstrip())
    if os.path.exists(tmp):
        os.remove(tmp)
    return None


def kernel_path():
    """Resolve the kernel binary, once per process: $POSQ_KERNEL override,
    then a posq_kernel binary built beside the source if it executes here
    (it is linked against the FLINT it was built with), then an automatic
    build from posq_kernel.c into the cache directory. Raises
    KernelUnavailable (named) when none of these run."""
    if "path" in _RESOLVED:
        return _RESOLVED["path"]
    override = os.environ.get("POSQ_KERNEL")
    if override:
        if not _runs(override):
            raise KernelUnavailable(
                f"POSQ-KERNEL-UNAVAILABLE: $POSQ_KERNEL={override!r} does not"
                f" execute here -- unset it or point it at a working build")
        _RESOLVED["path"] = override
        return override
    log = []
    if _runs(KERNEL):
        _RESOLVED["path"] = KERNEL
        return KERNEL
    log.append("no runnable posq_kernel binary beside the source (not built "
               "yet, or FLINT shared-library / CPU architecture mismatch)")
    built = _try_build(log)
    if built:
        _RESOLVED["path"] = built
        return built
    raise KernelUnavailable(
        "POSQ-KERNEL-UNAVAILABLE: " + "; ".join(log)
        + " -- build it yourself (needs FLINT development headers): "
          "cc -O2 -o posq_kernel posq_kernel.c -lflint -lmpfr -lgmp -lm, "
          "then set POSQ_KERNEL to the result")

def man_exp4(x):
    m = x.mid(); r = x.rad()
    mm, me = m.man_exp(); rm, re = r.man_exp()
    return (int(mm), int(me), int(rm), int(re))

def from_man_exp4(q):
    mm, me, rm, re = (int(t) for t in q)
    prec_save = ctx.prec
    ctx.prec = max(prec_save, 2048 + 64)     # exact reconstruction
    out = arb(mm) * arb(2) ** me + (arb(rm) * arb(2) ** re) * arb(0, 1)
    ctx.prec = prec_save
    return out

def write_rules(path, rules):
    """rules: dict u1/u2/v3 -> {nodes: [arb], weights: [arb]}."""
    with open(path, "w") as f:
        for tag in ("u1", "u2", "v3"):
            r = rules[tag]
            f.write(f"{len(r['nodes'])}\n")
            for x in r["nodes"]:
                f.write(" ".join(str(t) for t in man_exp4(x)) + "\n")
            for w in r["weights"]:
                f.write(" ".join(str(t) for t in man_exp4(w)) + "\n")

def write_tables(path, tabs_by_value):
    """tabs_by_value: dict value_index y (0..15) -> {(i,j,l): Fraction}.
    Missing values/entries -> 0. Layout fixed: y, i 0..4, j 0..3, l 0..1."""
    with open(path, "w") as f:
        for y in range(16):
            D = tabs_by_value.get(y, {})
            for i in range(5):
                for j in range(4):
                    for l in range(2):
                        v = D.get((i, j, l), Fr(0))
                        f.write(f"{v.numerator} {v.denominator}\n")

def write_job(path, mode, prec, pnum, pden, i0, i1, rulepath, outpath,
              vecs, tabpath="-"):
    with open(path, "w") as f:
        f.write(f"PREC {prec}\nMODE {mode}\nPNUM {pnum}\nPDEN {pden}\n")
        f.write(f"I0 {i0}\nI1 {i1}\nRULES {rulepath}\nTABLES {tabpath}\n")
        f.write(f"OUT {outpath}\nNVEC {len(vecs)}\n")
        for v in vecs:
            assert len(v) == 16
            f.write(" ".join(str(int(t)) for t in v) + "\n")

def run_kernel(jobpath, timeout=1800):
    """kernel cap: 1800 s per sweep."""
    r = subprocess.run([kernel_path(), jobpath], capture_output=True,
                       text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"posq_kernel failed: {r.returncode} {r.stderr[:500]}")
    return r.stdout

def run_chunked(tagpath, mode, prec, pnum, pden, n1, rulepath, vecs,
                tabpath="-", nchunk=2, timeout=1800):
    """run the kernel in nchunk i-range chunks (kernel cap 1800 s/sweep law);
    returns (summed balls per vec, total kernel cpu_s). Sum of partial sums
    == the streamed total (same register; rung-1 chunked assembly law)."""
    from flint import arb, ctx
    parts = []
    cpu = 0.0
    for ci in range(nchunk):
        i0 = ci * n1 // nchunk
        i1 = (ci + 1) * n1 // nchunk
        jp = f"{tagpath}.job{ci}"
        op = f"{tagpath}.out{ci}"
        write_job(jp, mode, prec, pnum, pden, i0, i1, rulepath, op, vecs,
                  tabpath=tabpath)
        run_kernel(jp, timeout=timeout)
        hdr, balls = read_out(op)
        cpu += float(hdr["CPU_S"])
        parts.append(balls)
        os.remove(jp)
    prec_save = ctx.prec
    ctx.prec = 256
    tots = []
    for s in range(len(vecs)):
        t = parts[0][s]
        for ci in range(1, nchunk):
            t = t + parts[ci][s]
        tots.append(t)
    ctx.prec = prec_save
    return tots, cpu

def read_out(path):
    """returns (header dict, [arb ball per vec])."""
    lines = [l.split() for l in open(path) if l.strip()]
    h = lines[0]
    hdr = {h[i]: h[i + 1] for i in range(0, len(h), 2)}
    tots = {}
    for l in lines[1:]:
        assert l[0] == "TOT"
        tots[int(l[1])] = from_man_exp4(l[2:6])
    n = len(tots)
    return hdr, [tots[s] for s in range(n)]

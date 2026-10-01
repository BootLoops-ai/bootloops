#!/usr/bin/env python3
"""amflow_kit.keypred — offline ibp-cache key predictor, yaml-byte canonicalized.
Member of amflow-kit, the house kit around the AMFlow.cpp fork.

Predicts the exact key amflow-cpp's ibp_cache will compute for a
kira input directory, so cache injections land where the live run looks.

Source of truth (the convention is DERIVED from, not guessed):
  amflow-cpp src/ibp/ibp_cache.cpp (ships in the AMFlow.cpp fork,
  the sibling repository amflow-cpp)
    - kInputFiles order + "\\0@" name "\\0=" framing + FNV-1a 64
    - canonicalize_jobs_yaml(): jobs.yaml is hashed AFTER (1) dropping
      n_thread:/n_threads:/threads: lines and (2) stripping ", d: <int>"
      from seed-spec lines (topologies+sectors) BEFORE select_integrals.
      The d: inside select_mandatory_recursively (Masters mode, after
      select_integrals) is RESULT-determining and is KEPT.
  amflow-cpp src/ibp/kira_yaml.cpp (same tree)
    - the writer whose output is the ONE pinned byte convention: LF line
      endings, trailing final newline, fixed indentation.  AMFLOW_REDUCE_DCAP=1
      appends ", d: N" to the Reduce seed line AT WRITE TIME; the hash-time
      canonicalizer strips it again — so post-restore work-dir bytes (with
      ", d: 1") and key-time conventions (without) key IDENTICALLY under the
      canonical form, and ONLY under the canonical form.
  amflow-cpp src/ibp/reduce.cpp + src/ibp/numd_subst.cpp (same tree)
    - numeric-d Reduce legs (AMFLOW_DIFFEQ_NUMD set, FireFly mode) SALT
      the key with the string "AMFLOW_DIFFEQ_NUMD=<value>": the d ->
      P/Q substitution rewrites tmp/SYSTEM files that are NOT hashed
      inputs, so the unsalted key of a numeric-d run is byte-identical
      to the symbolic run's and the two must never share a cache entry.
      The salt is folded AFTER the input files under its own
      "\\0@salt\\0=" frame; empty salt folds NOTHING, keeping unsalted
      keys byte-identical.  Masters-mode legs stay UNSALTED (their
      result is d-independent) so numeric-d and symbolic runs share
      preheat entries.  <value> is "P" or "P/Q", each within signed
      64-bit, Q > 0 (validate_numd).

Threat class this closes: offline key predictions that hash RAW jobs.yaml
bytes of one d-cap convention miss live keys computed from the other.
The canonical form predicts the same key from EITHER byte convention —
checked by --selftest on synthetic fixture pairs in ../tests/fixtures/keypred/
(one toy
input dir per byte convention, planted so the equalities/inequalities
between their keys are known in advance).

HONEST SCOPE: canonicalization removes the BYTE-CONVENTION miss class only.
A mirror whose semantic shape differs from what amflow writes (sectors/r
box, run_firefly vs triangular mode, target set) still misses; that class
stays covered by the run-read re-point law.  --explain prints the
canonical bytes so shape drift is inspectable.

CLI: (run from tools/amflow-kit as `python3 -m amflow_kit.keypred ARGS`
      or `python3 amflow_kit/keypred.py ARGS`)
  keypred.py key DIR [--raw] [--explain] [--salt=S | --numd[=P/Q]]
                                             # canonical (default) or raw key
  keypred.py predict DIR [--explain] [--salt=S | --numd[=P/Q]]
                                             # mirror door: pin layer (CRLF/
                                             # final-newline repair) + canonical
  keypred.py --selftest                      # fixture + mutation + salt battery

  --numd=P/Q  salts with the exact live string "AMFLOW_DIFFEQ_NUMD=P/Q"
              (a numeric-d Reduce leg's key); bare --numd reads the value
              from the AMFLOW_DIFFEQ_NUMD env var, as the live run does.
  --salt=S    folds S verbatim (the generic extra-salt door).

Exit codes: 0 ok, 1 selftest/validation failure, 2 usage/missing dir.
"""
import os
import sys

# Canonical kira input list — MUST mirror ibp_cache.cpp kInputFiles.
KINPUT = [
    "config/integralfamilies.yaml",
    "config/kinematics.yaml",
    "jobs.yaml",
    "jobs_kira2math.yaml",
    "preferred",
    "target",
]

FNV_OFFSET = 0xcbf29ce484222325
FNV_PRIME = 0x100000001b3
MASK = (1 << 64) - 1

# Extra-salt framing — MUST mirror ibp_cache.cpp ("\0@salt\0=", 8 bytes).
SALT_FRAME = b"\x00@salt\x00="
NUMD_ENV = "AMFLOW_DIFFEQ_NUMD"
INT64_MIN = -(1 << 63)
INT64_MAX = (1 << 63) - 1


def _fnv1a(h, data: bytes):
    for b in data:
        h ^= b
        h = (h * FNV_PRIME) & MASK
    return h


def canonicalize_jobs_yaml(raw: str) -> str:
    """Exact port of ibp_cache.cpp canonicalize_jobs_yaml (line-surgery,
    byte-preserving elsewhere).  Do NOT add normalizations here — this
    function must track the C++ byte-for-byte; convention repair for
    hand-built mirrors lives in pin_jobs_yaml()."""
    out = []
    past_select = False
    i = 0
    n = len(raw)
    while i < n:
        eol = raw.find("\n", i)
        end = n if eol == -1 else eol
        line = raw[i:end]

        if "select_integrals:" in line:
            past_select = True

        stripped = line.lstrip(" \t")
        if (stripped.startswith("n_thread:") or stripped.startswith("n_threads:")
                or stripped.startswith("threads:")):
            i = n if eol == -1 else eol + 1
            continue  # drop line entirely (incl. newline)

        if (not past_select and "topologies:" in line
                and "sectors:" in line):
            dp = line.find(", d: ")
            if dp != -1:
                q = dp + 5
                while q < len(line) and (line[q] == "-" or line[q].isdigit()):
                    q += 1
                line = line[:dp] + line[q:]

        out.append(line)
        if eol != -1:
            out.append("\n")
        i = n if eol == -1 else eol + 1
    return "".join(out)


def pin_jobs_yaml(raw: str):
    """Pin hand-built mirror bytes onto the amflow writer convention
    BEFORE canonicalization: kira_yaml.cpp emits LF-only lines and a
    trailing final newline; a mirror edited on another convention would
    hash differently in amflow (which reads raw bytes) so the prediction
    must be made over the repaired form.  Returns (pinned_text, repairs)."""
    repairs = []
    if "\r" in raw:
        raw = raw.replace("\r\n", "\n").replace("\r", "\n")
        repairs.append("crlf->lf")
    if raw and not raw.endswith("\n"):
        raw += "\n"
        repairs.append("final-newline-added")
    return raw, repairs


def validate_numd(v: str) -> str:
    """Exact port of numd_subst.cpp validate_numd: the value is "P" or
    "P/Q" with P a signed and Q an unsigned run of ASCII digits, each
    within signed 64-bit (the std::stol guard Kira's coefficient parser
    imposes), Q > 0.  Returns the value unchanged; raises ValueError on
    anything else."""
    slash = v.find("/")
    p_str = v if slash == -1 else v[:slash]
    q_str = "1" if slash == -1 else v[slash + 1:]

    def is_int(s: str, allow_sign: bool) -> bool:
        if not s:
            return False
        i = 1 if allow_sign and s[0] in "+-" else 0
        digits = s[i:]
        return bool(digits) and all(c in "0123456789" for c in digits)

    if not is_int(p_str, allow_sign=True) or not is_int(q_str, allow_sign=False):
        raise ValueError(f"{NUMD_ENV}: value must be an integer or rational "
                         f"\"P/Q\" (got '{v}')")
    p, q = int(p_str), int(q_str)
    if not (INT64_MIN <= p <= INT64_MAX) or not (0 < q <= INT64_MAX):
        raise ValueError(f"{NUMD_ENV}: P and Q must each fit in a signed "
                         f"64-bit integer with Q > 0 (got '{v}')")
    return v


def numd_salt(value: str) -> str:
    """The exact salt string a numeric-d Reduce leg folds into its key —
    the env assignment as reduce.cpp builds it ("AMFLOW_DIFFEQ_NUMD=" +
    validated value)."""
    return NUMD_ENV + "=" + validate_numd(value)


def ibp_cache_key(dirpath: str, raw: bool = False, pin: bool = False,
                  explain=None, salt: str = "") -> str:
    """Replica of ibp_cache.cpp ibp_cache_key().  raw=True hashes jobs.yaml
    byte-exact (forensics only — reproduces the OLD defective predictions);
    pin=True runs pin_jobs_yaml first (mirror door).  salt="" (the
    default) folds nothing — unsalted keys are byte-identical to the
    salt-free convention; a non-empty salt is folded AFTER the input
    files under the SALT_FRAME framing (numeric-d Reduce legs pass
    numd_salt(value))."""
    h = FNV_OFFSET
    for rel in KINPUT:
        p = os.path.join(dirpath, rel)
        if not os.path.exists(p):
            continue
        h = _fnv1a(h, b"\x00@")
        h = _fnv1a(h, rel.encode())
        h = _fnv1a(h, b"\x00=")
        with open(p, "rb") as f:
            data = f.read()
        if rel == "jobs.yaml" and not raw:
            text = data.decode()
            reps = []
            if pin:
                text, reps = pin_jobs_yaml(text)
            text = canonicalize_jobs_yaml(text)
            data = text.encode()
            if explain is not None:
                explain[rel] = {"canonical_bytes": text, "pin_repairs": reps}
        h = _fnv1a(h, data)
        if explain is not None and rel != "jobs.yaml":
            explain[rel] = {"bytes": len(data)}
    if salt:
        h = _fnv1a(h, SALT_FRAME)
        h = _fnv1a(h, salt.encode())
        if explain is not None:
            explain["salt"] = salt
    return format(h, "016x")


# ---------------------------------------------------------------------------
# Selftest: synthetic fixtures (../tests/fixtures/keypred/) — a toy one-loop
# family ("toybox") written in
# the kira input format, one directory per byte convention, PLANTED so the
# key relations are known in advance:
#   dcap_a / dcap_b / threads share ONE canonical key (the core claim: the
#   ", d: 1" seed-cap byte and n_thread lines do not change the canonical
#   key) while the raw keys of dcap_a and threads differ from it (the miss
#   mechanism this tool exists to kill); masters_d1 vs masters_d2 differ in
#   the RESULT-determining select_mandatory_recursively d:, which the
#   canonicalizer must KEEP, so their canonical keys must differ.
# The hex constants are regression pins generated from these fixtures by
# this replica at fixture-creation time (the convention itself derives from
# ibp_cache.cpp / kira_yaml.cpp, which ship in the AMFlow.cpp fork,
# the sibling repository amflow-cpp);
# the cross-fixture equalities/inequalities asserted below are the
# load-bearing checks, and the FNV-1a core is anchored to its published
# public test vector.
# ---------------------------------------------------------------------------

# the kit's test data: tools/amflow-kit/tests/fixtures/keypred/ (one dir up
# from this package)
FIXDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "tests", "fixtures", "keypred")

BASE_KEY = "44061221bf77d053"    # shared canonical key of the dcap/threads trio
SALTED_KEY = "433bdd9b7654c7cc"  # the same trio under numd_salt("-3/7")

# (fixture, expected canonical key, expected raw key or None)
CONTROLS = [
    # d-cap pair: a carries ", d: 1", b does not; canonical keys identical,
    # raw-hashing a gives the (planted) miss key.
    ("dcap_a", BASE_KEY, "04eba56dbe20b0d6"),
    ("dcap_b", BASE_KEY, BASE_KEY),
    # n_thread line present: dropped by the canonicalizer, kept by raw.
    ("threads", BASE_KEY, "5674ecfbce1f5ecd"),
    # Masters mode (no target/jobs_kira2math files: exercises the
    # skip-missing-file path); the select_mandatory_recursively d: is
    # result-determining and KEPT — d:1 vs d:2 key differently.
    ("masters_d1", "9ca1d74c87666a98", "9ca1d74c87666a98"),
    ("masters_d2", "4bfcc30a829be7fb", "4bfcc30a829be7fb"),
]


def selftest() -> int:
    import shutil
    import tempfile
    fails = []
    ok = []

    # hash-core anchor: published FNV-1a 64 test vector (fail-only check)
    if format(_fnv1a(FNV_OFFSET, b"foobar"), "016x") != "85944171f73967e8":
        fails.append("FNV-1a core: public test vector 'foobar' mismatch")

    keys = {}
    for fix, want_canon, want_raw in CONTROLS:
        d = os.path.join(FIXDIR, fix)
        if not os.path.isdir(d):
            fails.append(f"{fix}: fixture dir missing")
            continue
        got_c = ibp_cache_key(d)
        got_r = ibp_cache_key(d, raw=True)
        keys[fix] = (got_c, got_r)
        if got_c != want_canon:
            fails.append(f"{fix}: canonical {got_c} != pinned {want_canon}")
        else:
            ok.append(f"{fix}: canonical == pinned {want_canon}")
        if want_raw is not None and got_r != want_raw:
            fails.append(f"{fix}: raw {got_r} != pinned {want_raw}")

    # planted cross-fixture relations (fail-only; the load-bearing content)
    if len(keys) == len(CONTROLS):
        if not (keys["dcap_a"][0] == keys["dcap_b"][0] == keys["threads"][0]):
            fails.append("planted equality broken: dcap_a/dcap_b/threads "
                         "canonical keys differ")
        if keys["dcap_a"][1] == keys["dcap_a"][0]:
            fails.append("planted inequality broken: raw(dcap_a) == canonical "
                         "(the d-cap byte no longer matters?)")
        if keys["threads"][1] == keys["threads"][0]:
            fails.append("planted inequality broken: raw(threads) == canonical "
                         "(the n_thread line no longer matters?)")
        if keys["masters_d1"][0] == keys["masters_d2"][0]:
            fails.append("planted inequality broken: masters d:1 == d:2 "
                         "(result-determining d: must be KEPT)")

    # Mutation battery on the dcap_a fixture: byte-convention drifts a mirror
    # can carry; `predict` (pin+canonical) must recover the base key.
    src = os.path.join(FIXDIR, "dcap_a")
    live = BASE_KEY
    muts = {
        "crlf": lambda t: t.replace("\n", "\r\n"),
        "no-final-newline": lambda t: t.rstrip("\n"),
        "dcap-stripped": lambda t: t.replace(", d: 1}", "}"),
        "thread-line": lambda t: t.replace(
            "    preferred_masters:",
            "    n_thread: 16\n    preferred_masters:", 1),
    }
    if os.path.isdir(src):
        for name, f in muts.items():
            tmp = tempfile.mkdtemp(prefix="keypred_mut_")
            try:
                shutil.copytree(src, os.path.join(tmp, "d"))
                jp = os.path.join(tmp, "d", "jobs.yaml")
                with open(jp) as fh:
                    t = fh.read()
                with open(jp, "w", newline="") as fh:
                    fh.write(f(t))
                got = ibp_cache_key(os.path.join(tmp, "d"), pin=True)
                if got != live:
                    fails.append(f"mutation {name}: predict {got} != base {live}")
                else:
                    ok.append(f"mutation {name}: predict == base")
                # negative control: raw hashing of the mutated bytes must
                # NOT equal the base key (else the mutation is vacuous)
                if name in ("crlf", "thread-line"):
                    if ibp_cache_key(os.path.join(tmp, "d"), raw=True) == live:
                        fails.append(f"mutation {name}: raw unexpectedly == base"
                                     " (vacuous control)")
            finally:
                shutil.rmtree(tmp, ignore_errors=True)

    # Salt battery — the extra-salt fold (numeric-d Reduce legs).  Planted
    # relations: empty salt reproduces the legacy key byte-identically; a
    # non-empty salt moves the key (the numeric-d/symbolic collision class
    # this fold exists to kill); distinct numd values key distinctly; the
    # salt composes with canonicalization (the dcap/threads trio still
    # shares ONE salted key) and with the pin layer (salted predict on a
    # CRLF mirror recovers the pinned salted key).
    s_numd = numd_salt("-3/7")
    if s_numd != "AMFLOW_DIFFEQ_NUMD=-3/7":
        fails.append(f"numd_salt byte drift: {s_numd!r}")
    d_a = os.path.join(FIXDIR, "dcap_a")
    if os.path.isdir(d_a):
        got_s = ibp_cache_key(d_a, salt=s_numd)
        if got_s != SALTED_KEY:
            fails.append(f"salt: dcap_a salted {got_s} != pinned {SALTED_KEY}")
        else:
            ok.append(f"salt: dcap_a numd salted == pinned {SALTED_KEY}")
        if ibp_cache_key(d_a, salt="") != ibp_cache_key(d_a):
            fails.append("salt: empty salt no longer byte-identical to the "
                         "unsalted key")
        else:
            ok.append("salt: empty salt == unsalted key")
        if got_s == ibp_cache_key(d_a):
            fails.append("planted inequality broken: salted == unsalted "
                         "(numeric-d and symbolic entries would collide)")
        else:
            ok.append("salt: salted != unsalted")
        if ibp_cache_key(d_a, salt=numd_salt("1/2")) == got_s:
            fails.append("planted inequality broken: numd -3/7 and 1/2 "
                         "share a key")
        else:
            ok.append("salt: distinct numd values key distinctly")
        trio = {ibp_cache_key(os.path.join(FIXDIR, f), salt=s_numd)
                for f in ("dcap_a", "dcap_b", "threads")
                if os.path.isdir(os.path.join(FIXDIR, f))}
        if trio != {SALTED_KEY}:
            fails.append("planted equality broken: dcap_a/dcap_b/threads "
                         "salted canonical keys differ")
        else:
            ok.append("salt: dcap/threads trio shares one salted key")
        tmp = tempfile.mkdtemp(prefix="keypred_salt_")
        try:
            shutil.copytree(d_a, os.path.join(tmp, "d"))
            jp = os.path.join(tmp, "d", "jobs.yaml")
            with open(jp) as fh:
                t = fh.read()
            with open(jp, "w", newline="") as fh:
                fh.write(t.replace("\n", "\r\n"))
            if ibp_cache_key(os.path.join(tmp, "d"), pin=True,
                             salt=s_numd) != SALTED_KEY:
                fails.append("salt: salted predict on crlf mirror missed the "
                             "pinned salted key")
            else:
                ok.append("salt: salted predict recovers pinned key on "
                          "crlf mirror")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    else:
        fails.append("salt: dcap_a fixture dir missing")

    # numd format gate — exact port of the C++ validator; refusals are
    # fail-only, the aggregate gets one pass line.
    numd_good = ["301", "-3/7", "+5", "1/2", "9223372036854775807"]
    numd_bad = ["", "1/-2", "1/0", "a", "1.5", "3/", "/7", "--3",
                "9223372036854775808", "1/9223372036854775808"]
    numd_fails_before = len(fails)
    for v in numd_good:
        try:
            validate_numd(v)
        except ValueError as e:
            fails.append(f"numd gate: rejected valid '{v}': {e}")
    for v in numd_bad:
        try:
            validate_numd(v)
            fails.append(f"numd gate: accepted malformed '{v}'")
        except ValueError:
            pass
    if len(fails) == numd_fails_before:
        ok.append(f"numd format gate: {len(numd_good)} accepted, "
                  f"{len(numd_bad)} refused")

    for line in ok:
        print(f"  PASS {line}", flush=True)
    for line in fails:
        print(f"  FAIL {line}", flush=True)
    print(f"keypred selftest: {len(ok)} pass, {len(fails)} fail", flush=True)
    return 1 if fails else 0


def main(argv):
    if "--selftest" in argv:
        return selftest()
    args = [a for a in argv if not a.startswith("--")]
    flags = {a for a in argv if a.startswith("--")}
    if len(args) != 2 or args[0] not in ("key", "predict"):
        sys.stderr.write(__doc__.split("CLI:")[1])
        return 2
    verb, d = args
    if not os.path.isdir(d):
        sys.stderr.write(f"keypred: not a directory: {d}\n")
        return 2
    salt_flags = [a for a in flags
                  if a == "--numd" or a.startswith(("--salt=", "--numd="))]
    if len(salt_flags) > 1:
        sys.stderr.write("keypred: give at most one of --salt= / --numd\n")
        return 2
    salt = ""
    if salt_flags:
        opt = salt_flags[0]
        try:
            if opt.startswith("--salt="):
                salt = opt[len("--salt="):]
                if not salt:
                    raise ValueError("--salt= needs a non-empty value")
            elif opt.startswith("--numd="):
                salt = numd_salt(opt[len("--numd="):])
            else:  # bare --numd: read the env var the live run reads
                v = os.environ.get(NUMD_ENV, "")
                if not v:
                    raise ValueError(
                        f"--numd without a value needs {NUMD_ENV} set")
                salt = numd_salt(v)
        except ValueError as e:
            sys.stderr.write(f"keypred: {e}\n")
            return 2
    explain = {} if "--explain" in flags else None
    key = ibp_cache_key(d, raw=("--raw" in flags and verb == "key"),
                        pin=(verb == "predict"), explain=explain, salt=salt)
    print(key, flush=True)
    if explain is not None:
        j = explain.get("jobs.yaml", {})
        if j.get("pin_repairs"):
            sys.stderr.write(f"# pin repairs: {j['pin_repairs']}\n")
        if "salt" in explain:
            sys.stderr.write(f"# salt: {explain['salt']}\n")
        if "canonical_bytes" in j:
            sys.stderr.write("# canonical jobs.yaml bytes as hashed:\n")
            sys.stderr.write(j["canonical_bytes"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

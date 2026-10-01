#!/usr/bin/env bash
# amflow_smoke.sh — post-rebuild smoke test for an amflow_cli binary.
#
# Usage:  amflow_smoke.sh BINARY_PATH [JOB_TIMEOUT_SECONDS]
#
# Runs 2 known-answer jobs through BINARY_PATH and validates the outputs with
# tools/amflow-kit/amflow_output_lint.py plus a >=20-significant-digit numeric compare
# against reference outputs:
#   job vac : in_vac2Bprobe.json   (2L vacuum, 5 dotted ints, leads -2..-1)
#             reference = fixture out_vac2Bprobe.json (vendor-confirmed).
#             KNOWN BUG PROBE: a broken build (double-precision
#             Frobenius-resonance underflow) returns all-zero here.
#   job tri : 1-loop triangle, a copy of in_m3soft.json cut to its first 2
#             integrals; reference generated ONCE with the vendor binary
#             into $REF_DIR and reused thereafter.
#
# Originals are never modified: job specs are COPIED + rewritten (work_dir),
# each run gets a fresh run dir under $AMFLOW_SMOKE_ROOT, and each binary
# gets its OWN IBP-cache subdir (a shared cache would let a broken binary
# coast on cache HITs from a healthy one).
#
# Exit 0 iff BOTH jobs pass (binary rc + lint + numeric compare).
set -uo pipefail

# needs bash >= 4 (associative arrays); macOS ships bash 3.2
[ "${BASH_VERSINFO[0]}" -ge 4 ] || { echo "FATAL: amflow_smoke.sh needs bash >= 4 (macOS stock bash is 3.2 — brew install bash)"; exit 2; }

BIN="${1:?usage: amflow_smoke.sh BINARY_PATH [JOB_TIMEOUT_SECONDS]}"
JOB_TIMEOUT="${2:-1800}"

# Wall clamp: GNU timeout(1) is absent on stock macOS (coreutils installs it
# as gtimeout). Absent both, run unclamped — loudly, since a hung job then
# holds the smoke open.
if command -v timeout >/dev/null 2>&1; then TIMEOUT_BIN="timeout"
elif command -v gtimeout >/dev/null 2>&1; then TIMEOUT_BIN="gtimeout"
else
  TIMEOUT_BIN=""
  echo "WARN: no timeout(1)/gtimeout(1) on PATH — jobs run WITHOUT the ${JOB_TIMEOUT}s wall clamp (macOS: brew install coreutils)"
fi

TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LINT="$TOOLS_DIR/amflow_output_lint.py"
# Paths are env-overridable; defaults assume the repo layout (fixtures beside this
# kit) and a scratch root under $TMPDIR. AMFLOW_SMOKE_WRAP may name a launch
# wrapper (e.g. a jemalloc wrapper); unset = run the binary directly.
WRAP="${AMFLOW_SMOKE_WRAP:-}"
FIX="${AMFLOW_SMOKE_FIX:-$TOOLS_DIR/../fixtures/amflow_smoke}"    # read-only fixtures (see README.md there)
VENDOR_BIN="${AMFLOW_SMOKE_VENDOR_BIN:-}"   # reference binary; only needed the first time the triangle reference is generated
SMOKE_ROOT="${AMFLOW_SMOKE_ROOT:-${TMPDIR:-/tmp}/amflow_smoke}"
REF_DIR="${AMFLOW_SMOKE_REF_DIR:-$SMOKE_ROOT/ref}"
CACHE_ROOT="${AMFLOW_SMOKE_CACHE_ROOT:-$SMOKE_ROOT/cache}"
RUN_ROOT="${AMFLOW_SMOKE_RUN_ROOT:-$SMOKE_ROOT/runs}"
DIGITS=20

# FERMATPATH must point at a fer64 binary for amflow_cli; inherited from the
# caller's environment (not set here).

[ -x "$BIN" ]  || { echo "FATAL: binary not executable: $BIN"; exit 2; }
[ -f "$LINT" ] || { echo "FATAL: lint script missing: $LINT"; exit 2; }
[ -d "$FIX" ]  || { echo "FATAL: fixtures dir missing: $FIX (set AMFLOW_SMOKE_FIX)"; exit 2; }
if [ -n "$WRAP" ] && [ ! -x "$WRAP" ]; then echo "FATAL: AMFLOW_SMOKE_WRAP not executable: $WRAP"; exit 2; fi

BTAG="$(echo -n "$BIN" | sha1sum | cut -c1-8)"
export AMFLOW_IBP_CACHE="$CACHE_ROOT/$BTAG"      # per-binary cache isolation
mkdir -p "$AMFLOW_IBP_CACHE" "$REF_DIR"
RUN="$RUN_ROOT/run_$(date +%Y%m%dT%H%M%S)_${BTAG}_$$"
mkdir -p "$RUN"
echo "binary    : $BIN"
echo "run dir   : $RUN"
echo "ibp cache : $AMFLOW_IBP_CACHE"

# ---- job specs (copy + rewrite work_dir; cut tri to first 2 integrals) ----
python3 - "$FIX" "$RUN" <<'PYEOF'
import json, sys
fix, run = sys.argv[1], sys.argv[2]
j = json.load(open(f"{fix}/in_vac2Bprobe.json"))
j["work_dir"] = f"{run}/work_vac"
json.dump(j, open(f"{run}/in_vac.json", "w"), indent=1)
j = json.load(open(f"{fix}/in_m3soft.json"))
j["integrals"] = j["integrals"][:2]
j["work_dir"] = f"{run}/work_tri"
json.dump(j, open(f"{run}/in_tri.json", "w"), indent=1)
PYEOF

# ---- triangle reference: generate ONCE with the vendor binary ----
TRI_REF="$REF_DIR/out_tri2_ref.json"
if [ ! -s "$TRI_REF" ]; then
  echo "[ref] generating triangle reference with vendor binary (one-time)..."
  [ -x "$VENDOR_BIN" ] || { echo "FATAL: vendor binary missing: '$VENDOR_BIN' (set AMFLOW_SMOKE_VENDOR_BIN)"; exit 2; }
  mkdir -p "$RUN/refwork"
  python3 -c "import json,sys; j=json.load(open('$RUN/in_tri.json')); j['work_dir']='$RUN/refwork/work_tri_ref'; json.dump(j,open('$RUN/in_tri_ref.json','w'),indent=1)"
  ( export AMFLOW_IBP_CACHE="$CACHE_ROOT/vendor_ref"; mkdir -p "$AMFLOW_IBP_CACHE"
    ${TIMEOUT_BIN:+"$TIMEOUT_BIN" "$JOB_TIMEOUT"} ${WRAP:+"$WRAP"} "$VENDOR_BIN" "$RUN/in_tri_ref.json" "$TRI_REF" \
      > "$RUN/log_tri_ref.log" 2>&1 )
  rc=$?
  if [ $rc -ne 0 ] || [ ! -s "$TRI_REF" ]; then
    echo "FATAL: vendor triangle-reference generation failed (rc=$rc), log: $RUN/log_tri_ref.log"
    rm -f "$TRI_REF"; exit 2
  fi
  python3 "$LINT" -q "$TRI_REF" || { echo "FATAL: vendor triangle reference fails lint"; rm -f "$TRI_REF"; exit 2; }
  echo "[ref] triangle reference stored: $TRI_REF"
fi

# ---- high-precision numeric compare (Decimal prec=120; floats cap at ~16d) ----
numcompare () {  # numcompare OUT REF DIGITS
python3 - "$1" "$2" "$3" <<'PYEOF'
import json, re, sys
from decimal import Decimal, getcontext
getcontext().prec = 120
BALL = re.compile(r"^\[(?P<mid>[^\s\]]*)\s*\+/-\s*[^\]]*\]$")
def mid(s):
    s = str(s).strip()
    m = BALL.match(s)
    if m:
        s = m.group("mid").strip()
        if s in ("", "+", "-"): return Decimal(0)
    if "/" in s:
        a, b = s.split("/", 1); return Decimal(a) / Decimal(b)
    return Decimal(s)
def coefs(path):
    r = json.load(open(path)).get("result", [])
    out = []
    for integ in r:
        d = {}
        for c in integ.get("coefficients", []):
            v = c["value"]
            re_s, im_s = (v.get("re","0"), v.get("im","0")) if isinstance(v, dict) else (v[0], v[1])
            d[c["order"]] = (mid(re_s), mid(im_s))
        out.append(d)
    return out
out_p, ref_p, need = sys.argv[1], sys.argv[2], int(sys.argv[3])
out, ref = coefs(out_p), coefs(ref_p)
if len(out) != len(ref):
    print(f"COMPARE FAIL: integral count {len(out)} != ref {len(ref)}"); sys.exit(1)
worst, worst_at, fails, ncmp = Decimal("Infinity"), None, [], 0
for k, (do, dr) in enumerate(zip(out, ref)):
    scale = max([abs(a) + abs(b) for a, b in dr.values()] + [Decimal(0)])
    for order, (rr, ri) in dr.items():
        if rr == 0 and ri == 0: continue
        ncmp += 1
        if order not in do:
            fails.append(f"int{k} eps^{order}: missing in output"); continue
        oo = do[order]
        for part, r, o in (("re", rr, oo[0]), ("im", ri, oo[1])):
            err = abs(o - r)
            den = abs(r) if r != 0 else scale
            rel = err / den if den != 0 else err
            digs = Decimal(999) if rel == 0 else -rel.log10()
            if digs < worst: worst, worst_at = digs, f"int{k} eps^{order} {part}"
            if digs < need: fails.append(f"int{k} eps^{order} {part}: only {digs:.1f} digits")
print(f"compared {ncmp} nonzero ref coefficients; worst agreement "
      f"{'>120' if worst > 500 else f'{worst:.1f}'} digits at {worst_at}")
if fails:
    print("COMPARE FAIL (need >= %d digits):" % need)
    for f in fails[:10]: print("  " + f)
    sys.exit(1)
print(f"COMPARE PASS (>= {need} digits)")
sys.exit(0)
PYEOF
}

# ---- run one job: run_job NAME INSPEC REF LINTARGS... ----
declare -A STATUS WALL
run_job () {
  local name="$1" inspec="$2" ref="$3"; shift 3
  local out="$RUN/out_${name}.json" log="$RUN/log_${name}.log"
  echo; echo "=== job $name : $(basename "$inspec") ==="
  local t0=$SECONDS
  ${TIMEOUT_BIN:+"$TIMEOUT_BIN" "$JOB_TIMEOUT"} ${WRAP:+"$WRAP"} "$BIN" "$inspec" "$out" > "$log" 2>&1
  local rc=$? ; WALL[$name]=$((SECONDS - t0))
  echo "binary rc=$rc, wall=${WALL[$name]}s, log=$log"
  if [ $rc -ne 0 ] || [ ! -s "$out" ]; then
    [ $rc -eq 124 ] && echo "TIMED OUT after ${JOB_TIMEOUT}s"
    tail -3 "$log" | sed 's/^/  | /'
    STATUS[$name]="FAIL(run rc=$rc)"; return
  fi
  if ! python3 "$LINT" "$@" "$out"; then STATUS[$name]="FAIL(lint)"; return; fi
  if ! numcompare "$out" "$ref" "$DIGITS"; then STATUS[$name]="FAIL(compare)"; return; fi
  STATUS[$name]="PASS"
}

run_job vac "$RUN/in_vac.json" "$FIX/out_vac2Bprobe.json" --class vacuum --in "$RUN/in_vac.json"
run_job tri "$RUN/in_tri.json" "$TRI_REF" --class kinematic

echo; echo "===== SMOKE SUMMARY for $BIN ====="
fail=0
for name in vac tri; do
  printf "  %-4s %-18s wall=%ss\n" "$name" "${STATUS[$name]:-FAIL(notrun)}" "${WALL[$name]:--}"
  [ "${STATUS[$name]:-FAIL}" = "PASS" ] || fail=1
done
[ $fail -eq 0 ] && echo "RESULT: PASS" || echo "RESULT: FAIL"
exit $fail

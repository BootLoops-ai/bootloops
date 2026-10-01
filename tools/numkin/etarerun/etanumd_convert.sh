#!/bin/bash
[ "$(uname -s)" = Linux ] || { echo "FATAL: $(basename "$0") requires Linux — its kira-liveness checks read /proc" >&2; exit 2; }
# etanumd_convert.sh SRC_TARGET_REDUCE DST EPS_PQ [PARALLEL] [MEM_GB] [--nolaunch]
#
# Convert an amflow η-DE target_reduce/ shard (2-var: d, eta; Fermat back-sub
# stalled) into a numeric-d 1-var (eta) FireFly rerun, reusing the already-
# generated SYSTEM files.  d is fixed to 4 - 2*eps, eps = P/Q rational.
#
# Operates on the amflow target_reduce/ dir shape (config/, jobs.yaml, preferred, target,
# sectormappings/, results/{preferredMasters,<fam>/masters},
# tmp/<fam>/{SYSTEM_*.gz,SYSTEMconfig,masters,VER_*}).
#
#  - copies config, preferred, target, sectormappings verbatim
#  - copies tmp/<fam>/{SYSTEMconfig,masters}; substitutes d -> (NUM/DEN) inside
#    every SYSTEM_*.gz (token-safe: alphabet must be exactly {Eq,d,eta})
#  - does NOT copy tmp/<fam>/VER_* (Fermat forward-solve intermediates for the
#    symbolic-d run — inconsistent with substituted system)
#  - copies results/{preferredMasters,<fam>/masters} so kira sees master list;
#    does NOT copy results/kira.db (fresh FireFly)
#  - rewrites jobs.yaml: run_initiate: true, run_firefly: true (drop
#    run_triangular / run_back_substitution — FireFly path replaces them)
#  - launches kira detached (setsid+nohup) with jemalloc, FERMATPATH, prlimit
#
# Footguns applied:
#   * kira --set_value only touches Fermat path; FireFly reads SYSTEM_*.gz
#     directly => sed-substitute is the only working route
#   * std::stol: NUM,DEN must each be < 9.2e18  (checked)
#   * FERMATPATH must be set even for pure-FF runs (kira preflight)
#   * setsid nohup </dev/null so kira survives parent-process exit
#
set -euo pipefail
SRC=${1:?need SRC target_reduce dir}
DST=${2:?need DST dir}
EPS_PQ=${3:?need eps as P/Q e.g. 1/10000}
PAR=${4:-4}
MEM_GB=${5:-20}
NOLAUNCH=0; [ "${6:-}" = "--nolaunch" ] && NOLAUNCH=1

# --- sanity on source shape -------------------------------------------------
[ -d "$SRC/tmp" ]            || { echo "[etanumd] ERR: no $SRC/tmp"; exit 1; }
[ -d "$SRC/config" ]         || { echo "[etanumd] ERR: no $SRC/config"; exit 1; }
[ -d "$SRC/sectormappings" ] || { echo "[etanumd] ERR: no $SRC/sectormappings"; exit 1; }
[ -f "$SRC/jobs.yaml" ]      || { echo "[etanumd] ERR: no $SRC/jobs.yaml"; exit 1; }

# --- compute d = 4 - 2*eps as reduced NUM/DEN -------------------------------
P=${EPS_PQ%%/*}; Q=${EPS_PQ##*/}
[[ "$P" =~ ^-?[0-9]+$ && "$Q" =~ ^[0-9]+$ ]] || { echo "[etanumd] ERR: eps must be integer P/Q, got '$EPS_PQ'"; exit 1; }
read NUM DEN < <(python3 -c "from fractions import Fraction as F; d=4-2*F($P,$Q); print(d.numerator, d.denominator)")
# std::stol guard (kira parses coefficient tokens with stol)
python3 -c "import sys; n,d=$NUM,$DEN; sys.exit(0 if abs(n)<9_200_000_000_000_000_000 and d<9_200_000_000_000_000_000 else 1)" \
  || { echo "[etanumd] ERR: d=$NUM/$DEN overflows std::stol"; exit 1; }
DVAL="($NUM\\/$DEN)"   # sed-safe (escape / for s/// replacement)
echo "[etanumd] eps=$EPS_PQ  =>  d=$NUM/$DEN"

# --- build DST --------------------------------------------------------------
rm -rf "$DST"; mkdir -p "$DST"
for f in config sectormappings; do cp -a "$SRC/$f" "$DST/"; done
for f in preferred preferred_masters target de_targets; do
  [ -e "$SRC/$f" ] && cp -a "$SRC/$f" "$DST/"
done
# results/: masters + preferredMasters only (NEVER kira.db / kira2*)
mkdir -p "$DST/results"
[ -f "$SRC/results/preferredMasters" ] && cp -a "$SRC/results/preferredMasters" "$DST/results/"
for fam in "$SRC"/results/*/; do
  [ -d "$fam" ] || continue
  famname=$(basename "$fam")
  [ -f "$fam/masters" ] && { mkdir -p "$DST/results/$famname"; cp -a "$fam/masters" "$DST/results/$famname/"; }
done

# --- token-alphabet verify + substitute d -> (NUM/DEN) ----------------------
NSYS_TOTAL=0
for fam in "$SRC"/tmp/*/; do
  famname=$(basename "$fam"); mkdir -p "$DST/tmp/$famname"
  for aux in SYSTEMconfig masters; do
    [ -f "$fam$aux" ] && cp -a "$fam$aux" "$DST/tmp/$famname/"
  done
  # alphabet check on largest SYSTEM file
  BIG=$(ls -S "$fam"SYSTEM_*.gz 2>/dev/null | head -1) || true
  [ -n "$BIG" ] || { echo "[etanumd] ERR: no SYSTEM_*.gz in $fam"; exit 2; }
  ALPHA=$(zcat "$BIG" | grep -ao '[a-zA-Z][a-zA-Z0-9_]*' | sort -u | tr '\n' ' ')
  if [ "$ALPHA" != "Eq d eta " ]; then
    echo "[etanumd] ERR: token alphabet '$ALPHA' != 'Eq d eta' in $famname"; exit 3
  fi
  n=0
  for g in "$fam"SYSTEM_*.gz; do
    zcat "$g" | sed -E "s/\bd\b/$DVAL/g" | gzip -1 > "$DST/tmp/$famname/$(basename "$g")"
    n=$((n+1))
  done
  echo "[etanumd] $famname: $n SYSTEM files substituted d -> $DVAL"
  # residual standalone-d check
  res=$(for g in "$DST/tmp/$famname/"SYSTEM_*.gz; do zcat "$g" | grep -c '\bd\b' || true; done | awk '{s+=$1} END{print s+0}')
  [ "$res" = 0 ] || { echo "[etanumd] ERR: $res residual 'd' tokens after sed"; exit 4; }
  # file-count check
  a=$(ls "$fam"SYSTEM_*.gz | wc -l); b=$(ls "$DST/tmp/$famname/"SYSTEM_*.gz | wc -l)
  [ "$a" = "$b" ] || { echo "[etanumd] ERR: SYSTEM count $a != $b"; exit 5; }
  NSYS_TOTAL=$((NSYS_TOTAL+n))
done
echo "[etanumd] total SYSTEM files: $NSYS_TOTAL"

# --- jobs.yaml: add run_firefly, drop run_triangular/run_back_substitution --
python3 - "$SRC/jobs.yaml" "$DST/jobs.yaml" <<'PYEOF'
import sys, re
src, dst = sys.argv[1], sys.argv[2]
lines = open(src).read().splitlines()
out, saw_ff = [], False
for ln in lines:
    s = ln.strip()
    if s.startswith('run_triangular') or s.startswith('run_back_substitution'):
        continue
    if s.startswith('run_firefly'):
        m = re.match(r'^(\s*)', ln); indent = m.group(1)
        out.append(f'{indent}run_firefly: true'); saw_ff = True; continue
    out.append(ln)
    if s.startswith('run_initiate') and not saw_ff:
        m = re.match(r'^(\s*)', ln); indent = m.group(1)
        out.append(f'{indent}run_firefly: true'); saw_ff = True
# ensure kira2math export step present (a SRC jobs.yaml can lose it during
# a FireFly restage)
txt = '\n'.join(out)
if 'kira2math' not in txt:
    fam = None
    for ln in out:
        m = re.search(r'topologies:\s*\[(\w+)\]', ln)
        if m: fam = m.group(1); break
    if fam:
        out.append(' - kira2math:')
        out.append('    target:')
        out.append(f'     - [{fam}, target]')
open(dst, 'w').write('\n'.join(out) + '\n')
if not saw_ff:
    sys.exit("[etanumd] ERR: could not place run_firefly in jobs.yaml")
PYEOF
echo "[etanumd] jobs.yaml written (run_firefly: true, kira2math ensured)"

# --- provenance stamp -------------------------------------------------------
cat > "$DST/ETANUMD_SLICE.json" <<EOF
{"src":"$SRC","eps":"$EPS_PQ","d_num":$NUM,"d_den":$DEN,"nsys":$NSYS_TOTAL,"par":$PAR,"mem_gb":$MEM_GB,"built_utc":"$(date -u +%FT%TZ)"}
EOF

# --- launch -----------------------------------------------------------------
if [ "$NOLAUNCH" = 1 ]; then
  echo "[etanumd] --nolaunch: DST prepared at $DST (not launching kira)"
  exit 0
fi
JEMALLOC=${JEMALLOC:-/usr/lib/x86_64-linux-gnu/libjemalloc.so.2}  # env override; Debian/Ubuntu path is the fallback
[ -f "$JEMALLOC" ] || JEMALLOC=/lib/x86_64-linux-gnu/libjemalloc.so.2
FERMAT=${FERMATPATH:-fer64}  # fer64 (Fermat) on PATH or via FERMATPATH
[ -x "$FERMAT" ] || FERMAT=/usr/share/Ferl7/fer64  # fallback install location
MEM_KB=$((MEM_GB*1024*1024))
cd "$DST"
env LD_PRELOAD="$JEMALLOC" FERMATPATH="$FERMAT" OPENBLAS_NUM_THREADS=2 \
  setsid nohup nice -n 5 prlimit --as=$((MEM_KB*1024)) \
  kira jobs.yaml --parallel="$PAR" > kira_numd.log 2>&1 < /dev/null &
FORKPID=$!
sleep 2
# find real kira child by cwd (guard against set -e spurious exit)
REAL=$( (for p in $(pgrep -x kira); do [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$(readlink -f "$DST")" ] && echo $p; done | head -1) || true )
echo "[etanumd] launched kira in $DST (fork pid=$FORKPID, real pid=${REAL:-?}, parallel=$PAR, mem_cap=${MEM_GB}GB)"
echo "${REAL:-$FORKPID}" > "$DST/kira.pid"
exit 0

#!/bin/bash
# Surd smoke test: ONE support group of the QCD quark-jet channel q_qbpqpgq (sigma = 1234, first pair 34, chart x_2 = 1) end to
# end from the included cell pickles, writing ONLY under the target directory  $1  (default ${TMPDIR:-/tmp}/surd_smoke.<stamp>;
# a target inside the package tree, or an existing one, is refused):
#   S1  slice_run.py     two Brown-linear fibration steps, exact in Q[x,t,j]/(j^2+1)          -> ckpt/<tag>/stage{1,2}.pkl, SLICE_S12_<tag>.json
#   G1  slice_gate12.py  stage 1-2 check: symbolic one-fold object vs the independent residue-route fiber function (Arb) at 2 t x 4 x  -> SLICE_G12_<tag>.json
#   S3  stage3.py        function-level third integration, exact; value at t = 1/2, 4/5 vs the two-chart Arb piece evaluator -> SLICE_S3_<tag>.json
#   A   slice_classes.py + assemble.py + evaluate.py: the one-group (incomplete) channel file in the documented format, re-evaluated by the
#       file evaluator (canonical keys, polyroots, hpath) and compared with S3's own evaluator (independent numeric route on the same exact object)
# and checks the results against frozen 30-digit values.  Prints "[smoke] RESULT: PASS" and exits 0 on success.
# Single thread; wall ~3-5 min; peak RSS ~0.2 GB.  Needs python3 with python-flint, sympy, mpmath (see GUIDE.md REQUIREMENTS).
set -euo pipefail
TESTS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); ROOT=$(dirname "${TESTS:?}"); S="${ROOT:?}/slice/scripts"
[ -f "${S}/slice_prov.py" ] || { echo "smoke_test.sh: cannot locate slice/scripts under ${ROOT}" >&2; exit 2; }
STAMP=$(date -u +%Y%m%dT%H%M%SZ); TARGET=${1:-${TMPDIR:-/tmp}/surd_smoke.${STAMP}}
[ ! -e "${TARGET}" ] || { echo "smoke_test.sh: target ${TARGET} exists; give a fresh directory" >&2; exit 2; }
mkdir -p "${TARGET:?}"; TARGET=$(cd "${TARGET:?}" && pwd)
case "${TARGET}/" in "$(cd "${ROOT}" && pwd)"/*) rmdir "${TARGET}"; echo "smoke_test.sh: refusing a target inside the package tree (${ROOT}); give a scratch directory as \$1" >&2; exit 2;; esac
TAG=q_qbpqpgq_s1234_p34_g2; CH=q_qbpqpgq; SG=1234; PR=34
# frozen expectations: closed form / ext(t) of this group at t = 1/2 and 4/5 (32 digits shown; the checks certify >= 30)
EXP_T12="0.12095802771852580968992232548419"; EXP_T45="0.01955527772608710742851706675049"; MIN_DIGITS=30
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export SURD_WORK="${TARGET:?}" SURD_DATA="${TARGET:?}"          # every write of the pipeline goes under the target
ulimit -v "${SURD_VMEM_KB:-8388608}" 2>/dev/null || true
mkdir -p "${TARGET:?}/logs"; L="${TARGET:?}/logs"; cd "${S:?}"
echo "[smoke] package ${ROOT}; target ${TARGET}; group ${TAG}; $(date -u +%Y-%m-%dT%H:%M:%SZ)"
step() { local name=$1; shift; local t0; t0=$(date +%s); echo "[smoke] ${name} ..."; nice -n 19 "$@" > "${L}/${name}.out" 2>&1 || { echo "[smoke] ${name} FAILED (rc $?); log ${L}/${name}.out" >&2; tail -5 "${L}/${name}.out" >&2; exit 1; }; echo "[smoke] ${name} done ($(( $(date +%s) - t0 )) s)"; }
step S1_slice_run   python3 slice_run.py --channel "$CH" --sigma "$SG" --pair "$PR"
step G1_slice_gate12 python3 slice_gate12.py --tag "$TAG"
step S3_stage3      python3 stage3.py --tag "$TAG" --t 1/2,4/5
step A1_slice_classes python3 slice_classes.py --channels "$CH"
step A2_assemble    python3 assemble.py --channels "$CH"
step A3_evaluate    python3 evaluate.py --file "${TARGET}/out/E4C_LO_QCD_dipole_slice_${CH}.json.gz" --t 1/2 --dps 60 --digits 30
python3 - "$TARGET" "$TAG" "$EXP_T12" "$EXP_T45" "$MIN_DIGITS" "$ROOT" "${L}/A3_evaluate.out" <<'PY'
import json, sys, os
from fractions import Fraction as Fr
import mpmath as mp
root, tag, e12, e45, mind, pkg, evalout = sys.argv[1:8]; ok = True; mp.mp.dps = 60
g = json.load(open(os.path.join(root, "SLICE_G12_%s.json" % tag))); s3 = json.load(open(os.path.join(root, "SLICE_S3_%s.json" % tag)))
print("[smoke] stage 1-2 check (symbolic one-fold object vs residue-route fiber function): min digits %.1f PASS %s" % (g["min_agree_digits"], g["PASS(>=30)"])); ok &= bool(g["PASS(>=30)"]) and g["min_agree_digits"] >= float(mind)
print("[smoke] stage-3 check (closed form vs two-chart Arb piece evaluator): min digits %s PASS %s" % (s3["min_agree_digits"], s3["PASS(>=30)"])); ok &= bool(s3["PASS(>=30)"])
v12 = None
for row, exp in zip(s3["gate_rows"], (e12, e45)):
    val = row["closed_form_value/ext"].strip("()").split(" ")[0]; same = val.startswith(exp)
    if row["t"] == "1/2": v12 = mp.mpf(val)
    print("[smoke] t=%s closed form %s... frozen %s... match %s" % (row["t"], val[:24], exp, same)); ok &= same
ins = os.path.realpath(s3["producer"]["script"]).startswith(os.path.realpath(pkg) + os.sep)
print("[smoke] producer script recorded in the stage-3 run record is inside the package:", ins); ok &= ins
# assembled one-group file re-evaluated by the file evaluator: F_file(1/2) = mult * w_flavor(5)*S_ch * (closed form/ext)(1/2)
cl = json.load(open(os.path.join(root, "SLICE_CLASSES.json")))["channels"]["q_qbpqpgq"]; mult = [c["mult"] for c in cl["classes"] if c["rep"] == "1234"][0]
ev = json.load(open(evalout)); F = mp.mpf(ev["value"]); pred = mult * 4 * mp.mpf(e12)          # w_flavor(nf=5) = nf - 1 = 4, S_ch = 1 for q_qbpqpgq
rel = abs(F / pred - 1); dg = 99.0 if rel == 0 else float(-mp.log10(rel))
print("[smoke] assembled file at t=1/2 by evaluate.py: %s (mode %s, %s digits est.); / (mult %d x 4) vs frozen: %.1f digits" % (mp.nstr(F, 32), ev["mode"], ev["digits_est"], mult, dg)); ok &= dg >= float(mind) - 1
print("[smoke] RESULT:", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
PY

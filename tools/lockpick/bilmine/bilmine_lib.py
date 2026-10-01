#!/usr/bin/env python3
# lockpick bilmine member — miner library: prereg-sealed bilinear-relation lattice miner on certified period data (homogeneous LLL reduction + HNF lattice identity, capacity-lawful height ladder, windows/holdout law, script-emitted receipts).
"""bilmine_lib.py — bilinear-lattice miner library.

extends: lockpick (mplll_lattice.reduce_rows backend chain + pslq_gate.canonicalize
and the capacity/two-precision/controls discipline of the package guide).
Exact stage: sympy Rational/HNF + python-flint integer linear algebra.

Laws wired here (sealed by the member's own prereg ceremony, seal_prereg.py):
- ambient-dps input law: every decimal string parsed inside workdps(len+50);
- trim-order law: C entries rounded to the point's certified cap before any fit;
- capacity law printed per rung: (n_unknowns+1)*log10(H) + 20 <= d_leg;
- floors: |residual| <= row_scale * 10^-(d_leg-30), full-precision residual arbiter;
- canonicalization: pslq_gate.canonicalize (gcd 1, first nonzero > 0);
- producer blocks (+ optional producer-lint check) on every emitted JSON.

Driver data contract (no default paths are shipped; drivers refuse loudly):
  BILMINE_LEG    — the leg working directory (receipts land in $BILMINE_LEG/work)
  BILMINE_T5     — the certified connection-matrix input dir (conn_*.json, linkb_*.json)
  BILMINE_PFRAME — the exact generator-frame JSON (G_integer, invariant_form.P)
  BILMINE_C1_PI / BILMINE_C1_N — the C1 control's certified period 5-vector
                   JSON and its integer certificate JSON
  BILMINE_PRODUCER_LINT — OPTIONAL path to a producer-lint script; when unset
                   the lint step is skipped
"""
import hashlib, json, os, subprocess, sys
from fractions import Fraction
from math import gcd

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))  # package parent (tools/)
from lockpick.mplll_lattice import reduce_rows          # noqa: E402
from lockpick.pslq_gate import canonicalize             # noqa: E402


# ---------- data-contract plumbing (no defaults shipped) ----------

def require_env(name, purpose):
    v = os.environ.get(name)
    if not v:
        raise SystemExit(f"bilmine: environment variable {name} is required "
                         f"({purpose}); no default path is shipped")
    return v


def leg_dir():
    return require_env("BILMINE_LEG", "the leg working directory — receipts "
                       "land in $BILMINE_LEG/work")


def t5_dir():
    return require_env("BILMINE_T5", "the certified connection-matrix input "
                       "directory (conn_*.json, linkb_*.json)")


# ---------- producer/receipt plumbing ----------

def sha_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def emit(obj, out_path, script_path):
    stamp = subprocess.check_output(["date", "-u"]).decode().strip()
    obj["producer"] = {
        "leg": "BILMINE",
        "script": os.path.basename(script_path),
        "script_sha256": sha_file(script_path),
        "stamp_utc": stamp,
    }
    json.dump(obj, open(out_path, "w"), indent=1, default=str)
    lint = os.environ.get("BILMINE_PRODUCER_LINT")
    if not lint:
        return None
    r = subprocess.run([sys.executable, lint, "--check", out_path],
                       capture_output=True, text=True)
    line = (r.stdout or r.stderr).strip()
    print(line)
    if r.returncode != 0:
        raise RuntimeError(f"producer_lint failed on {out_path}: {line}")
    return line


# ---------- parsing (ambient-dps law + trim-order law) ----------

def parse_real(s, cap=None):
    with mp.workdps(len(s) + 50):
        x = mp.mpf(s)
        if cap is not None:
            x = mp.mpf(mp.nstr(x, cap, strip_zeros=False))
        return x


def parse_c_matrix(path, cap):
    """Stored conn_*.json -> 6x6 list of mpc (trimmed to cap sig digits), plus meta."""
    d = json.load(open(path))
    C = d["C_rows_local_cols_x0"]
    M = [[None] * 6 for _ in range(6)]
    im_max = mp.mpf(0)
    with mp.workdps(cap + 80):
        for i in range(6):
            for j in range(6):
                re_s, im_s = C[i][j]
                re = parse_real(re_s, cap)
                im = parse_real(im_s, cap)
                M[i][j] = mp.mpc(re, im)
                sc = max(abs(re), mp.mpf(1))
                im_max = max(im_max, abs(im) / sc)
    return M, {"tier": d["tier"], "point": d["point"], "dps_land": d["dps_land"],
               "row_labels": d["row_labels"], "col_labels_x0": d["col_labels_x0"],
               "max_rel_im": mp.nstr(im_max, 6), "sha256": sha_file(path)}


# ---------- window construction ----------

def sym_pairs():
    return [(k, l) for k in range(6) for l in range(k, 6)]      # 21


def anti_pairs():
    return [(k, l) for k in range(6) for l in range(k + 1, 6)]  # 15


def entry_list(sector, holdout_index=0):
    pairs = sym_pairs() if sector == "sym" else anti_pairs()
    fit = [(i, j) for (i, j) in pairs if i != holdout_index and j != holdout_index]
    hold = [(i, j) for (i, j) in pairs if i == holdout_index or j == holdout_index]
    return fit, hold


def build_rows(C, sector, entries, use_im, dps):
    """Rows of C^T S_loc C - S_mum = 0 on the given (i,j) FIT entries.

    Unknowns: [S_loc coords over sector pairs] + [S_mum coords over the FIT
    entries ONLY] (held-out S_mum entries are never unknowns — they are
    determined post-hoc and gated by the near-integer holdout check; the C3
    plant caught the free-coordinate defect of the wider parameterization).
    All products computed inside workdps(dps) — never at ambient dps.
    """
    loc = sym_pairs() if sector == "sym" else anti_pairs()
    rows, tags = [], []
    with mp.workdps(dps):
        zero = mp.mpc(0, 0)
        for (i, j) in entries:
            coeffs = []
            for (k, l) in loc:
                if sector == "sym":
                    v = C[k][i] * C[l][j] + (C[l][i] * C[k][j] if k != l else zero)
                else:
                    v = C[k][i] * C[l][j] - C[l][i] * C[k][j]
                coeffs.append(v)
            for (a, b) in entries:
                coeffs.append(mp.mpc(-1, 0) if (a, b) == (i, j) else zero)
            rows.append([c.real for c in coeffs]); tags.append((i, j, "re"))
            if use_im:
                rows.append([c.imag for c in coeffs]); tags.append((i, j, "im"))
    return rows, tags


def holdout_eval(C, sector, hold_entries, S_loc_vec, dps):
    """For each held-out entry (i,j): v = (C^T S_loc C)[i,j] at dps.

    Returns [{entry, re, im, nearest_int, dist_to_int_rel, im_rel}] with
    row-scale-relative distances (scale = max |product coeff| of that entry).
    """
    loc = sym_pairs() if sector == "sym" else anti_pairs()
    out = []
    with mp.workdps(dps):
        for (i, j) in hold_entries:
            acc = mp.mpc(0, 0)
            scale = mp.mpf(0)
            for t, (k, l) in enumerate(loc):
                if sector == "sym":
                    co = C[k][i] * C[l][j] + (C[l][i] * C[k][j] if k != l else 0)
                else:
                    co = C[k][i] * C[l][j] - C[l][i] * C[k][j]
                acc += co * S_loc_vec[t]
                scale = max(scale, abs(co))
            scale = scale or mp.mpf(1)
            ni = mp.nint(acc.real)
            out.append({"entry": (i, j),
                        "nearest_int": int(ni),
                        "dist_to_int_rel": abs(acc.real - ni) / scale,
                        "im_rel": abs(acc.imag) / scale})
    return out


# ---------- the homogeneous simultaneous miner ----------

def capacity_line(n_unknowns, H, d_leg):
    need = float((n_unknowns + 1) * mp.log10(mp.mpf(H)) + 20)
    return {"H": H, "need_digits": round(need, 1), "d_leg": d_leg,
            "lawful": need <= d_leg}


def mine_homogeneous(rows, n_unknowns, d_scale, hmax, bkz_beta=20):
    """LLL/BKZ mine for integer c (0 < max|c| <= hmax) with rows.c ~ 0."""
    m = len(rows)
    with mp.workdps(d_scale + 60):
        S = mp.mpf(10) ** d_scale
        scales = [max(abs(x) for x in r) or mp.mpf(1) for r in rows]
        lat = []
        for i in range(n_unknowns):
            right = [int(mp.nint(S * rows[r][i] / scales[r])) for r in range(m)]
            lat.append([1 if j == i else 0 for j in range(n_unknowns)] + right)
    red, backend, wall = reduce_rows(lat, None, bkz_beta)
    cands, seen = [], set()
    for r in red:
        c = r[:n_unknowns]
        h = max(abs(x) for x in c)
        if h == 0 or h > hmax:
            continue
        cc = canonicalize(c)
        if cc and cc not in seen:
            seen.add(cc)
            cands.append(list(cc))
    return cands, backend, wall


def residuals(rows, c, d_arbiter):
    out = []
    with mp.workdps(d_arbiter):
        for r in rows:
            s = max(abs(x) for x in r) or mp.mpf(1)
            v = mp.fsum(r[k] * c[k] for k in range(len(c)))
            out.append(abs(v) / s)
    return out


def passes_floor(res, d_leg, margin=30):
    thr = mp.mpf(10) ** (-(d_leg - margin))
    worst = max(res) if res else mp.mpf(0)
    return all(x <= thr for x in res), (float(mp.log10(worst)) if worst > 0 else -9999.0)


# ---------- basket degeneracy audit (G5) ----------

def audit_rows(rows, d_leg):
    zero_rows = [i for i, r in enumerate(rows) if max(abs(x) for x in r) == 0]
    n = len(rows[0]) if rows else 0
    zero_cols = [j for j in range(n)
                 if all(abs(r[j]) == 0 for r in rows)]
    dup = []
    with mp.workdps(80):
        keys = {}
        for i, r in enumerate(rows):
            s = max(abs(x) for x in r) or mp.mpf(1)
            key = tuple(mp.nstr(x / s, 30) for x in r)
            if key in keys:
                dup.append((keys[key], i))
            keys[key] = i
    return {"n_rows": len(rows), "zero_rows": zero_rows, "zero_cols": zero_cols,
            "duplicate_row_pairs": dup,
            "audit_pass": not zero_rows and not dup and not zero_cols,
            "floor_log10": float(-(d_leg - 30))}


# ---------- exact integer lattice assembly ----------

def hnf_basis(vectors):
    """Primitive canonical basis of the Z-span of integer row vectors."""
    if not vectors:
        return []
    import flint
    M = flint.fmpz_mat([[int(x) for x in v] for v in vectors])
    H = M.hnf()
    out = []
    for i in range(H.nrows()):
        row = [int(H[i, j]) for j in range(H.ncols())]
        if any(row):
            cv = canonicalize(row)
            out.append(list(cv))
    return out


def vec_to_mats(c, sector):
    loc = sym_pairs() if sector == "sym" else anti_pairs()
    n = len(loc)
    Sl = [[0] * 6 for _ in range(6)]
    Sm = [[0] * 6 for _ in range(6)]
    for t, (k, l) in enumerate(loc):
        Sl[k][l] = int(c[t])
        Sl[l][k] = int(c[t]) if sector == "sym" else -int(c[t])
    for t, (a, b) in enumerate(loc):
        Sm[a][b] = int(c[n + t])
        Sm[b][a] = int(c[n + t]) if sector == "sym" else -int(c[n + t])
    return Sl, Sm


def mats_to_vec(Sl, Sm, sector):
    loc = sym_pairs() if sector == "sym" else anti_pairs()
    return [Sl[k][l] for (k, l) in loc] + [Sm[a][b] for (a, b) in loc]


# ---------- exact invariant spaces (stage E) ----------

def _frac(x):
    return Fraction(int(x.p), int(x.q)) if hasattr(x, "p") else Fraction(str(x))


def _mat6(vals):
    return [[int(vals[i][j]) for j in range(6)] for i in range(6)]


def _mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(6)) for j in range(6)]
            for i in range(6)]


def _tr(A):
    return [[A[j][i] for j in range(6)] for i in range(6)]


def invariant_space(gens, sector, orientation="GtSG"):
    """Exact primitive HNF basis of {S sector-symmetric: g-invariant for all 10}.
    Pure-int construction + flint fmpz_mat nullspace (exact over Q)."""
    import flint
    loc = sym_pairs() if sector == "sym" else anti_pairs()
    n = len(loc)
    rows = []
    for g in gens:
        G = _mat6(g)
        Gt = _tr(G)
        img = []
        for t in range(n):
            k, l = loc[t]
            S = [[0] * 6 for _ in range(6)]
            S[k][l] += 1
            if k != l:
                S[l][k] += 1 if sector == "sym" else -1
            T = _mm(_mm(Gt, S), G) if orientation == "GtSG" else _mm(_mm(G, S), Gt)
            img.append([T[a][b] - S[a][b] for (a, b) in loc])
        for r in range(n):
            rows.append([img[t][r] for t in range(n)])
    A = flint.fmpz_mat(rows)
    X, nullity = A.nullspace()
    basis = []
    for c in range(nullity):
        iv = [int(X[r, c]) for r in range(n)]
        cv = canonicalize(iv)
        if cv and any(cv):
            basis.append(list(cv))
    return hnf_basis(basis)


def commutant_dim(gens):
    import flint
    rows = []
    for g in gens:
        G = _mat6(g)
        img = []
        for t in range(36):
            X = [[0] * 6 for _ in range(6)]
            X[t // 6][t % 6] = 1
            T = [[sum(G[i][k] * X[k][j] - X[i][k] * G[k][j] for k in range(6))
                  for j in range(6)] for i in range(6)]
            img.append([T[r // 6][r % 6] for r in range(36)])
        for r in range(36):
            rows.append([img[t][r] for t in range(36)])
    A = flint.fmpz_mat(rows)
    _, nullity = A.nullspace()
    return nullity


def exact_invariance_check(S_int, gens, orientation="GtSG"):
    """Exact integer check g^T S g == S (or g S g^T == S) for all gens."""
    import flint
    S = flint.fmpz_mat(S_int)
    for g in gens:
        G = flint.fmpz_mat(g)
        T = (G.transpose() * S * G) if orientation == "GtSG" else (G * S * G.transpose())
        if T != S:
            return False
    return True


def in_span_exact(vec, basis):
    """Is integer vec in the Q-span of integer basis rows? exact."""
    if not basis:
        return not any(vec)
    from sympy import Matrix
    A = Matrix(basis).T
    b = Matrix(vec)
    sol, params = A.gauss_jordan_solve(b)
    try:
        _ = sol
        return True
    except Exception:
        return False


# ---------- rational detection (gate G4a) ----------

def detect_rational(x, den_cap=10 ** 6, tol_log10=-40, near_zero_log10=None):
    """x (mpf) -> (Fraction or None, abs_frac_err_log10).

    FRACTIONAL-PART ABSOLUTE criterion (an instrument fix, caught by the G4
    smoke test): a q<=den_cap rational approximation of a large number is
    meaningless under a relative test (any 1e60-scale value passes at 1e-40
    relative). Here: split x = n + f (n integer, f in [0,1)), best-rational
    the FRACTIONAL part at denominator cap, accept iff |f - p/q| <= 10^tol
    ABSOLUTE. Zero detection (near_zero_log10, absolute vs caller scale)
    short-circuits."""
    if near_zero_log10 is not None and (x == 0 or abs(x) <= mp.mpf(10) ** near_zero_log10):
        return Fraction(0), -9999.0
    with mp.workdps(mp.mp.dps + 20):
        n = int(mp.floor(x))
        f = x - n
        fs = mp.nstr(f, 100, strip_zeros=False)
    fr = Fraction(fs).limit_denominator(den_cap)
    with mp.workdps(160):
        err = abs(mp.mpf(fr.numerator) / mp.mpf(fr.denominator) - f)
        ok = err <= mp.mpf(10) ** tol_log10
    out = Fraction(n) + fr if ok else None
    return out, (float(mp.log10(err)) if err > 0 else -9999.0)

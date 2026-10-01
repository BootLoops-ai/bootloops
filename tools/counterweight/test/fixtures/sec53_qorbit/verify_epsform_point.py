#!/usr/bin/env python3
"""verify_epsform_point.py -- INDEPENDENT exact check of a counterweight block result. Inputs: the block JSON (A in (d, x), the
converter's form) and the engine's out JSON (T, Atilde as Nemo strings in (eps, x) with `//`; nested lists COLUMN-major:
json[j][i] = M[i,j]). At each rational (eps, y) point the identity  T A T^-1 + (dT/dy) T^-1 == eps * Atilde  is checked
in exact Fractions, dT/dy by exact dual-number evaluation of the strings (no numerics, no sympy); Atilde is checked eps-free by
re-evaluating it at a second eps. Default points: (eps, y) = (1/223, 29/211) and (5/6, 17/139) (d = 4 - 2 eps = 890/223 and 7/3). --no-transpose runs the row-major reading as a CONTROL (must FAIL). Main-guarded; rc 0 PASS / 3 FAIL."""
import argparse, ast, hashlib, json, os, subprocess, sys
from fractions import Fraction as F


class Dual:
    __slots__ = ('a', 'b')

    def __init__(self, a, b=F(0)):
        self.a = F(a); self.b = F(b)

    def __add__(s, o): o = _d(o); return Dual(s.a + o.a, s.b + o.b)
    __radd__ = __add__

    def __sub__(s, o): o = _d(o); return Dual(s.a - o.a, s.b - o.b)

    def __rsub__(s, o): o = _d(o); return Dual(o.a - s.a, o.b - s.b)

    def __neg__(s): return Dual(-s.a, -s.b)

    def __mul__(s, o): o = _d(o); return Dual(s.a * o.a, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__

    def __truediv__(s, o): o = _d(o); return Dual(s.a / o.a, (s.b * o.a - s.a * o.b) / (o.a * o.a))

    def __rtruediv__(s, o): return _d(o) / s

    def __pow__(s, k):
        k = int(k); r = Dual(1)
        for _ in range(k): r = r * s
        return r


def _d(o):
    return o if isinstance(o, Dual) else Dual(o)


def eval_dual(s, env):
    tree = ast.parse(s.replace('^', '**'), mode='eval')

    def w(n):
        if isinstance(n, ast.Expression): return w(n.body)
        if isinstance(n, ast.Constant) and isinstance(n.value, int) and not isinstance(n.value, bool): return Dual(n.value)
        if isinstance(n, ast.Name) and n.id in env: return env[n.id]
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            v = w(n.operand); return -v if isinstance(n.op, ast.USub) else v
        if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.FloorDiv, ast.Pow)):
            a, b = w(n.left), w(n.right)
            if isinstance(n.op, ast.Add): return a + b
            if isinstance(n.op, ast.Sub): return a - b
            if isinstance(n.op, ast.Mult): return a * b
            if isinstance(n.op, (ast.Div, ast.FloorDiv)): return a / b
            assert b.b == 0 and b.a.denominator == 1 and b.a >= 0, 'power'
            return a ** int(b.a)
        raise ValueError(f'unsupported node {type(n).__name__} in {s[:80]!r}')
    return w(tree)


def mat_from_json(rows, transpose):
    n = len(rows)
    return [[str(rows[j][i]) if transpose else str(rows[i][j]) for j in range(n)] for i in range(n)]


def eval_mat(M, env):
    return [[eval_dual(M[i][j], env) if M[i][j] != '0' else Dual(0) for j in range(len(M))] for i in range(len(M))]


def val(M): return [[x.a for x in r] for r in M]


def der(M): return [[x.b for x in r] for r in M]


def matmul(A, B):
    n = len(A); return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def matinv(A):
    n = len(A); M = [list(r) + [F(int(i == j)) for j in range(n)] for i, r in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        piv = M[c][c]; M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [r[n:] for r in M]


def check(block, out, eps, y, transpose):
    d = 4 - 2 * eps
    envA = {'d': Dual(d), 'x': Dual(y, 1)}; envT = {'eps': Dual(eps), 'x': Dual(y, 1)}; envT2 = {'eps': Dual(eps + F(1, 97)), 'x': Dual(y, 1)}
    A = val(eval_mat(block['A'], envA))
    Tm = eval_mat(mat_from_json(out['T'], transpose), envT); T = val(Tm); dT = der(Tm)
    At = val(eval_mat(mat_from_json(out['Atilde'], transpose), envT)); At2 = val(eval_mat(mat_from_json(out['Atilde'], transpose), envT2))
    Ti = matinv(T); lhs = matmul(matmul(T, A), Ti); dTTi = matmul(dT, Ti)
    n = len(A); ok_id = all(lhs[i][j] + dTTi[i][j] == eps * At[i][j] for i in range(n) for j in range(n))
    ok_free = At == At2
    return ok_id, ok_free


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('block'); ap.add_argument('out'); ap.add_argument('--receipt', default=None)
    ap.add_argument('--points', default='1/223,29/211;5/6,17/139'); ap.add_argument('--no-transpose', action='store_true')
    a = ap.parse_args(); block = json.load(open(a.block)); out = json.load(open(a.out))
    if not out.get('ok') or 'T' not in out:
        print(f"verify: {os.path.basename(a.out)} carries no rotation (ok={out.get('ok')}; stop={str(out.get('stop'))[:120]}) -- nothing to verify"); sys.exit(2)
    res = []
    for p in a.points.split(';'):
        e, y = (F(t) for t in p.split(',')); ok_id, ok_free = check(block, out, e, y, not a.no_transpose)
        res.append({'eps': str(e), 'd': str(4 - 2 * e), 'y': str(y), 'identity_TAT_plus_dTT_eq_eps_Atilde': ok_id, 'Atilde_eps_free': ok_free})
    allok = all(r['identity_TAT_plus_dTT_eq_eps_Atilde'] and r['Atilde_eps_free'] for r in res)
    tag = 'row-major CONTROL' if a.no_transpose else 'column-major (file convention)'
    print(f"verify {os.path.basename(a.block)} n={block['n']} rows={block['rows']} [{tag}]: " + '; '.join(f"(eps={r['eps']}, y={r['y']}): identity {r['identity_TAT_plus_dTT_eq_eps_Atilde']} eps-free {r['Atilde_eps_free']}" for r in res) + f" -> {'PASS' if allok else 'FAIL'}")
    if a.receipt:
        json.dump({'stamp_utc': subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'], capture_output=True, text=True).stdout.strip(), 'block': a.block, 'block_sha256': hashlib.sha256(open(a.block, 'rb').read()).hexdigest(),
                   'out': a.out, 'out_sha256': hashlib.sha256(open(a.out, 'rb').read()).hexdigest(), 'transpose': not a.no_transpose, 'points': res, 'verdict': 'PASS' if allok else 'FAIL',
                   'method': 'exact Fraction arithmetic; dT/dy by dual numbers on the strings; T A T^-1 + T\' T^-1 == eps*Atilde entrywise; Atilde re-evaluated at eps + 1/97'}, open(a.receipt, 'w'), indent=1)
    sys.exit(0 if allok else 3)


if __name__ == '__main__':
    main()

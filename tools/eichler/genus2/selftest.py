#!/usr/bin/env python3
"""eichler.genus2 smoke battery — run: python3 genus2/selftest.py (from
tools/eichler/), or python3 selftest.py in this directory. The package-level
`tools/eichler/selftest.py` runs these legs together with the other members.

Exercises the public-runnable members and verifies the documented fail-closed
refusals of the reference-only members (see GUIDE.md BATTERY):

  S1  mestre_port: Clebsch/Igusa-Clebsch invariants vs the Sage doctest
      vectors (exact sympy equality) + Mestre conic for (I2,I4,I6,I10) =
      (1,2,3,4) proportional to the Sage Mestre_conic coefficient vector.
  S2  conic_fast: the built-in worked example runs in a scratch directory;
      the congruence diagonalization is re-checked exactly (P L P^T diagonal,
      diagonal == d) and, when a K-point is found, its back-substitution into
      the conic must be exactly zero.
  R1  invariant_harness without G2KIT_TRUE_JSON (or its alias
      EICHLER_GENUS2_TRUE_JSON): must refuse loudly, naming the env var and
      the required oracle-json schema (never a silent run).
  R2  pslq_law_template without `invariants_pipeline`: must refuse loudly as
      a read-and-adapt template (never a partial import).

Scratch: S2 writes its worked-example json into a temporary directory under
$EICHLER_WORK when that is set (else the system temp dir); nothing is written
into the package tree.

Exit 0 all PASS; assertion failure (nonzero exit) otherwise.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

ORACLE_ENV = ('G2KIT_TRUE_JSON', 'EICHLER_GENUS2_TRUE_JSON')


def _scratch_parent():
    w = os.environ.get('EICHLER_WORK')
    if w:
        os.makedirs(w, exist_ok=True)
        return w
    return None


def s1_mestre_port():
    import sympy as sp
    import mestre_port as mp_
    f1 = mp_.X**6 + mp_.Y**6
    assert mp_.clebsch_ABCD(f1) == (2, sp.Rational(2, 3), sp.Rational(-2, 9), 0)
    assert mp_.igusa_clebsch(f1) == (-240, 1620, -119880, -46656)
    t = sp.Symbol('t')
    p = t**6 + t**5 + t**4 + t**2 + 2
    ph = sp.expand(sum(sp.Poly(p, t).coeff_monomial(t**i) * mp_.X**i * mp_.Y**(6 - i)
                       for i in range(7)))
    assert mp_.igusa_clebsch(ph) == (-496, 6220, -955932, -1111784)
    # Mestre conic for (1,2,3,4) vs Sage Mestre_conic([1,2,3,4]) coefficients
    x_, y_, z_ = mp_.mestre_xyz(*map(sp.Integer, (1, 2, 3, 4)))
    L = mp_.mestre_conic_matrix(x_, y_, z_)
    got = (L[0, 0], 2 * L[0, 1], L[1, 1], 2 * L[0, 2], 2 * L[1, 2], L[2, 2])
    want = (-2572155000, -317736000, 1250755459200,
            2501510918400, 39276887040, 2736219686912)
    ratio = sp.nsimplify(got[0] / want[0])
    assert all(sp.simplify(g - ratio * w) == 0 for g, w in zip(got, want)), \
        f'Mestre conic not proportional to the Sage vector: {got}'
    print('S1 PASS  mestre_port invariants + conic vs Sage doctest vectors (exact)')


def s2_conic_fast():
    old = os.getcwd()
    with tempfile.TemporaryDirectory(dir=_scratch_parent()) as td:
        os.chdir(td)          # the worked example banks a json into its cwd
        try:
            import conic_fast as cf
        finally:
            os.chdir(old)
    # re-check the congruence diagonalization exactly: P L P^T == diag(d)
    A = cf.matmul(cf.matmul(cf.P, cf.L), cf.transpose(cf.P))
    for i in range(3):
        for j in range(3):
            if i == j:
                assert (A[i][j] - cf.d[i]).iszero(), f'diagonal mismatch at {i}'
            else:
                assert A[i][j].iszero(), f'off-diagonal residue at ({i},{j})'
    if cf.found is not None:
        assert cf.chk.iszero(), 'K-point back-check nonzero'
        print('S2 PASS  conic_fast: diagonalization exact, K-point on the conic')
    else:
        print('S2 PASS  conic_fast: diagonalization exact '
              '(no K-point in box -> obstruction candidate recorded, as documented)')


def r1_harness_refusal():
    env = {k: v for k, v in os.environ.items() if k not in ORACLE_ENV}
    r = subprocess.run([sys.executable, os.path.join(HERE, 'invariant_harness.py')],
                       capture_output=True, text=True, env=env,
                       cwd=_scratch_parent() or HERE)
    msg = r.stdout + r.stderr
    assert r.returncode != 0, 'invariant_harness ran without an oracle json'
    assert 'G2KIT_TRUE_JSON' in msg and 'rosenhain' in msg, \
        f'refusal did not name the requirement: {msg[-300:]}'
    print('R1 PASS  invariant_harness refuses loudly without G2KIT_TRUE_JSON')


def r2_template_refusal():
    r = subprocess.run([sys.executable, os.path.join(HERE, 'pslq_law_template.py')],
                       capture_output=True, text=True,
                       cwd=_scratch_parent() or HERE)
    msg = r.stdout + r.stderr
    assert r.returncode != 0, 'pslq_law_template ran without its pipeline'
    assert 'invariants_pipeline' in msg, \
        f'refusal did not name invariants_pipeline: {msg[-300:]}'
    print('R2 PASS  pslq_law_template refuses loudly without invariants_pipeline')


LEGS = (s1_mestre_port, s2_conic_fast, r1_harness_refusal, r2_template_refusal)


def run_all():
    """Run every leg; returns the number of legs (raises on any failure)."""
    for leg in LEGS:
        leg()
    return len(LEGS)


if __name__ == '__main__':
    n = run_all()
    print(f'eichler.genus2 selftest: ALL PASS ({n}/{n})')

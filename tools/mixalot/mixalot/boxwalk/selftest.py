"""mixalot.boxwalk selftest — must PASS before any production use (< 10 min).

 A. 2-variable synthetic with a sympy closed form: full produce() chain
    (plan -> refill+slab walk -> gates G1-G3 -> CRT -> exact Fraction)
    must equal sympy's exact integral.  Exercises the refill machinery.
 B. JC-quartet toy (lysozyme-shaped N=20, xxxx:15 then xxyy:5): walk window
    == dense oracle at EVERY step, BOTH primes (pilot p8 discipline).
 C. Mutation control: corrupted program/polynomial -> loud disagreement.
 D. Fiber recurrence: fit/lift/verify chain + recurrence extension mod p;
    shift-operator calculus (GCRD/right-factor recovery on closed-form
    operators, LCLM data semantics mod p, mutation must-fail, fitted-
    operator lift-and-reduce against the walk fiber).
 E. Singular syz smoke (skipped gracefully if Singular missing).
 F. Manifest verifier: independent replay of the A manifest must PASS;
    corrupted residue / corrupted value / corrupted program hash must
    FAIL at the named check; over-budget dense replay must REFUSE loudly.

Run:  python3 -m mixalot.boxwalk selftest      (from tools/mixalot/)
  or  python3 selftest.py                       (from this directory)
BOXWALK_SELFTEST_OUT=/path/file.json redirects the record off the tree.
"""
import os, sys, time, json
import numpy as np

if __package__ in (None, ''):
    # run as a plain script: tools/mixalot (three levels up) holds `mixalot`
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))))
from mixalot import boxwalk
from mixalot.boxwalk import core, driver, planner, fiber, syz
from mixalot.boxwalk.relations import RelationSystem

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = {}


def check(name, ok, detail=''):
    RESULTS[name] = bool(ok)
    print(f"[selftest] {name}: {'PASS' if ok else 'FAIL'} {detail}",
          flush=True)


def t_A_sympy():
    import sympy as sp
    from fractions import Fraction
    spec = boxwalk.load_spec({
        'name': 'selftest2var', 'variables': ['x', 'y'],
        'polynomials': {'a': {'1 0': 1, '0 1': 1},
                        'b': {'0 0': 1, '1 1': 1}},
        'target_u': {'a': 3, 'b': 2}, 'order': 'auto'})
    pl = planner.plan(spec, syz_timeout=60, log=lambda *a: None)
    man = driver.produce(spec, plan=pl, nprimes=6, batch=3, procs=1,
                         log=lambda *a: None)
    x, y = sp.symbols('x y')
    Zs = sp.nsimplify(sp.integrate(sp.integrate(
        (x + y) ** 3 * (1 + x * y) ** 2, (x, 0, 1)), (y, 0, 1)))
    ref = Fraction(int(Zs.p), int(Zs.q))
    mechs = [(s['class'], s['mech']) for s in pl['segments']]
    ok = man['ok'] and man['Z_fraction'] == ref
    check('A_sympy_exact', ok, f"Z={man.get('Z',{}).get('num')}/"
          f"{man.get('Z',{}).get('den')} ref={ref} mechs={mechs}")
    check('A_gates', all(g.get('pass') for g in man['gates'].values()),
          str({k: v.get('pass') for k, v in man['gates'].items()}))
    return spec, pl, man


def t_B_jc_toy():
    d = json.load(open(f'{HERE}/examples/jc_quartet.json'))
    d['target_u'] = {'xxxx': 15, 'xxyy': 5}
    d['order'] = ['xxxx', 'xxyy']
    spec = boxwalk.load_spec(d)
    pl = planner.plan(spec, syz_timeout=240, log=lambda *a: None)
    primes = core.primes31(2)
    bad = [0]

    def cb(u, M, valid, pvec):
        v = min(valid)
        for k, p in enumerate(pvec):
            ref = core.moment_window(spec, u, v - 1, int(p))
            sub = M[k][tuple(slice(0, v) for _ in range(spec.n))]
            if not np.array_equal(sub % int(p), ref % int(p)):
                bad[0] += 1
    res = driver.walk(spec, pl, primes, step_cb=cb)
    ok = (bad[0] == 0 and res['u'] == spec.target_u)
    check('B_jc_toy_walk_vs_oracle', ok,
          f"20 steps x 2 primes, mismatches={bad[0]}, "
          f"side={pl['side']}, mechs="
          f"{[(s['class'], s['mech']) for s in pl['segments']]}")
    return spec, pl


def t_C_mutation(spec2, pl2, spec_jc, pl_jc):
    g_refill = driver.gate_mutation(spec2, pl2, log=lambda *a: None)
    g_slab = driver.gate_mutation(spec_jc, pl_jc, log=lambda *a: None)
    check('C_mutation_refill', g_refill['pass'], f"mode={g_refill['mode']}")
    check('C_mutation_slab', g_slab['pass'], f"mode={g_slab['mode']}")


def t_D_fiber(spec2):
    rec = fiber.emit_fiber_recurrence(spec2, ray_class='a', depth=80,
                                      rmax=12, dmax=8,
                                      log=lambda *a: None)
    ok = rec['ok'] and rec['verified']
    check('D_fiber_fit_lift_verify', ok,
          f"(r,d)=({rec.get('r')},{rec.get('d')}) "
          f"bits={rec.get('max_coeff_bits')}")
    if not ok:
        return None
    # extension: recompute prefix mod a fresh prime, extend, compare oracle
    p = core.primes31(5)[-1]
    r, d = rec['r'], rec['d']
    Zs = fiber.fiber_mod(spec2, rec['u_base'], 'a', r + 5, p)
    ext = fiber.extend_mod(Zs, rec['coeffs'], r, d, p, 40)
    ref = fiber.fiber_mod(spec2, rec['u_base'], 'a', 40, p)
    check('D_fiber_extend', ext == ref, "extend to n=40 == oracle")
    return rec


def t_D_gcrd(spec2, rec):
    # closed-form shift operators: G = S-(n+1) annihilates n!,
    # A = S-2 annihilates 2^n, B = S-3 annihilates 3^n
    G = {'coeffs': {'0,0': '-1', '0,1': '-1', '1,0': '1'}, 'r': 1, 'd': 1}
    A = {'coeffs': {'0,0': '-2', '1,0': '1'}, 'r': 1, 'd': 0}
    B = {'coeffs': {'0,0': '-3', '1,0': '1'}, 'r': 1, 'd': 0}
    L1, L2 = fiber.op_mul(A, G), fiber.op_mul(B, G)
    got = fiber.gcrd_reduce(L1['coeffs'], L1['r'], L1['d'], other=L2)
    refG = fiber.gcrd_reduce(G['coeffs'], G['r'], G['d'])
    ok = (got['ok'] and got['reduced']
          and got['coeffs'] == refG['coeffs']
          and fiber.right_divides(L1, got) and fiber.right_divides(L2, got)
          and not fiber.right_divides(got, L1))
    check('D_gcrd_recover', ok, f"(r,d)=({got.get('r')},{got.get('d')})")
    # data semantics mod p: L1 = A.G kills n!; lclm(A,G) kills 2^n + n!
    # while A alone must not
    p = core.primes31(4)[-1]
    N = 30
    fact = [1] * (N + 1)
    for i in range(1, N + 1):
        fact[i] = fact[i - 1] * i % p
    mix = [(fact[i] + pow(2, i, p)) % p for i in range(N + 1)]
    M = fiber.lclm(A, G)
    ok2 = (fiber.annihilates_mod(L1, fact, p) and M['r'] == 2
           and fiber.annihilates_mod(M, mix, p)
           and not fiber.annihilates_mod(A, mix, p))
    check('D_gcrd_lclm', ok2, f"lclm (r,d)=({M['r']},{M['d']})")
    # mutation: one corrupted coefficient -> right division must fail and
    # the GCRD must degrade away from the planted factor
    bad = dict(L1['coeffs'])
    bad['0,0'] = str(int(bad['0,0']) + 1)
    badop = {'coeffs': bad, 'r': L1['r'], 'd': L1['d']}
    gbad = fiber.gcrd_reduce(bad, L1['r'], L1['d'], other=L2)
    ok3 = (not fiber.right_divides(badop, got)
           and gbad['coeffs'] != refG['coeffs'])
    check('D_gcrd_mutation', ok3,
          f"corrupted gcrd (r,d)=({gbad.get('r')},{gbad.get('d')})")
    # fitted operator: lift the emitted recurrence by a left factor,
    # reduce back, and check the result annihilates the walk fiber mod p
    if not (rec and rec.get('ok')):
        check('D_gcrd_fiber', False, 'no fitted recurrence available')
        return
    L = fiber.op_mul(A, rec)
    red = fiber.gcrd_reduce(L['coeffs'], L['r'], L['d'], other=rec)
    refR = fiber.gcrd_reduce(rec['coeffs'], rec['r'], rec['d'])
    Zs = fiber.fiber_mod(spec2, rec['u_base'], rec['ray_class'], 40, p)
    ok4 = (red['ok'] and red['reduced']
           and red['coeffs'] == refR['coeffs']
           and fiber.annihilates_mod(red, Zs, p))
    check('D_gcrd_fiber', ok4, f"(r,d)=({red.get('r')},{red.get('d')})")


def t_E_singular(spec_jc):
    import shutil
    if not shutil.which('Singular'):
        check('E_singular_smoke', True, 'SKIPPED (Singular not found)')
        return
    system = RelationSystem(spec_jc)
    fams, meta = syz.certified_families(spec_jc, system, ['xxxx'], 'xxyy',
                                        timeout=240, log=lambda *a: None)
    ok = (meta['n_certified'] > 0 and meta['n_certified'] == meta['n_raw']
          and meta.get('mutation_control_fails', False))
    check('E_singular_smoke', ok,
          f"source={meta['source']} certified={meta['n_certified']}/"
          f"{meta['n_raw']} mutation_fails={meta.get('mutation_control_fails')}")


def t_F_verify(spec2, pl2, man):
    from mixalot.boxwalk import verify
    quiet = lambda *a: None                                    # noqa: E731
    rep = verify.verify_manifest(spec2, man, nfresh=2, plan=pl2, log=quiet)
    check('F_verify_pass', rep['ok'],
          str({k: v['pass'] for k, v in rep['checks'].items()}))
    base = {k: v for k, v in man.items() if k != 'Z_fraction'}
    # corrupted residue -> V2 and V3 must fail
    m1 = json.loads(json.dumps(base))
    p0 = next(iter(m1['residues']))
    m1['residues'][p0] = (m1['residues'][p0] + 1) % int(p0)
    r1 = verify.verify_manifest(spec2, m1, nfresh=2, log=quiet)
    # corrupted banked value -> V2 and V4 must fail
    m2 = json.loads(json.dumps(base))
    m2['Z']['num'] = str(int(m2['Z']['num']) + 1)
    r2 = verify.verify_manifest(spec2, m2, nfresh=2, log=quiet)
    # over-budget dense replay -> named refusal, verdict False
    r3 = verify.verify_manifest(spec2, base, nfresh=1, max_cells=1,
                                log=quiet)
    # corrupted program hash -> V6 must fail (plan supplied)
    m4 = json.loads(json.dumps(base))
    m4['program_sha'] = '0' * 16
    r4 = verify.verify_manifest(spec2, m4, nfresh=2, plan=pl2, log=quiet)
    ok = (not r1['ok'] and not r1['checks']['V2_residues']['pass']
          and not r1['checks']['V3_reconstruction']['pass']
          and not r2['ok'] and not r2['checks']['V2_residues']['pass']
          and not r2['checks']['V4_fresh_replay']['pass']
          and not r3['ok']
          and r3['checks']['V4_fresh_replay'].get('refused')
          and not r4['ok'] and not r4['checks']['V6_program_sha']['pass'])
    check('F_verify_mutation', ok,
          'residue/value/budget/program corruptions all failed loudly')


def main():
    t0 = time.time()
    spec2, pl2, man2 = t_A_sympy()
    spec_jc, pl_jc = t_B_jc_toy()
    t_C_mutation(spec2, pl2, spec_jc, pl_jc)
    rec = t_D_fiber(spec2)
    t_D_gcrd(spec2, rec)
    t_E_singular(spec_jc)
    t_F_verify(spec2, pl2, man2)
    ok = all(RESULTS.values())
    out = {'PASS': ok, 'results': RESULTS,
           'wall_s': round(time.time() - t0, 1),
           'date': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    # BOXWALK_SELFTEST_OUT redirects the record off the tracked tree (e.g. a
    # mktemp file); unset, the record is written beside the package.
    out_path = os.environ.get('BOXWALK_SELFTEST_OUT') or f'{HERE}/SELFTEST.json'
    json.dump(out, open(out_path, 'w'), indent=1)
    print(f"[selftest] TOTAL: {'PASS' if ok else 'FAIL'} "
          f"({out['wall_s']}s) -> {out_path}", flush=True)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

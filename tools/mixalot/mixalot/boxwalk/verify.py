"""boxwalk.verify — independent replay verifier for banked manifests.

Re-checks a produce() MANIFEST.json against its spec WITHOUT the recorded
pivot program: the banked exact value is replayed mod fresh primes through
the dense oracle (direct convolution + separable moment contraction — an
algorithmically independent path from the recorded-program walk that made
the manifest), the CRT + Wang rational reconstruction is re-derived here
from the banked residues with local integer arithmetic, and every banked
residue is checked directly against the banked fraction.  Checks:

  V1_spec_sha        recomputed spec hash == banked spec_sha
  V2_residues        banked Z reduces to EVERY banked per-prime residue
  V3_reconstruction  local CRT + Wang ratrec over the banked residues
                     reproduces the banked Z exactly
  V4_fresh_replay    dense-oracle Z mod nfresh primes DISJOINT from the
                     banked set == banked Z mod those primes
  V5_gates           manifest ok=True, every recorded gate pass=True,
                     banked target_u matches the spec
  V6_program_sha     (only when a plan is supplied) recomputed program
                     hash == banked program_sha

verify_manifest() is read-only.  The dense replay's cost is estimated
first (coefficient-array cells) and the check REFUSES loudly — named
reason, overall verdict False — rather than run past max_cells.
"""
import json, math
from fractions import Fraction

from . import core


# ------------------------------------------------- local reconstruction --
# Deliberately local re-implementations (not shared with the producing
# path): the point of V3 is an independent derivation of the banked value.

def _crt(residues, moduli):
    r, M = 0, 1
    for a, p in zip(residues, moduli):
        t = (a - r) * pow(M % p, -1, p) % p
        r, M = r + M * t, M * p
    return r % M, M


def _ratrec(a, m):
    a %= m
    u0, u1, v0, v1 = m, a, 0, 1
    bound = math.isqrt(m // 2)
    while u1 > bound:
        q = u0 // u1
        u0, u1 = u1, u0 - q * u1
        v0, v1 = v1, v0 - q * v1
    if v1 == 0 or abs(v1) > bound or math.gcd(u1, abs(v1)) != 1:
        return None
    return Fraction(u1, v1) if v1 > 0 else Fraction(-u1, -v1)


def fresh_primes(k, exclude):
    """k primes descending from 2^31, skipping the excluded set."""
    from sympy import prevprime
    excl = {int(x) for x in exclude}
    out, p = [], (1 << 31)
    while len(out) < k:
        p = prevprime(p)
        if int(p) not in excl:
            out.append(int(p))
    return out


def dense_cells(spec, u):
    """Coefficient-array size of the dense replay at exponent vector u."""
    cells = 1
    for e in range(spec.n):
        cells *= 1 + sum(int(u[lab]) * spec.degs(lab)[e] for lab in u)
    return cells


def _as_spec(spec):
    if isinstance(spec, core.Spec):
        return spec
    if isinstance(spec, dict):
        return core.Spec(spec)
    return core.Spec.load(spec)


def _as_manifest(manifest):
    if isinstance(manifest, str):
        return json.load(open(manifest))
    return {k: v for k, v in manifest.items() if k != 'Z_fraction'}


def _print(msg):
    print(msg, flush=True)


def verify_manifest(spec, manifest, nfresh=2, plan=None,
                    max_cells=int(5e7), log=_print):
    """Verify a banked manifest (dict or MANIFEST.json path) against its
    spec (Spec, dict, or path).  Returns {'ok': bool, 'checks': {...}};
    ok is True only when every check ran and passed."""
    assert nfresh >= 1, "nfresh >= 1: a verdict needs at least one fresh prime"
    spec = _as_spec(spec)
    man = _as_manifest(manifest)
    checks = {}
    missing = [k for k in ('spec_sha', 'target_u', 'residues', 'Z', 'gates')
               if k not in man]
    if missing:
        checks['V0_shape'] = {'pass': False,
                              'detail': f'manifest missing {missing}'}
        return {'ok': False, 'name': man.get('name'), 'checks': checks}

    # V1 spec hash
    sha = spec.sha()
    checks['V1_spec_sha'] = {'pass': sha == man['spec_sha'],
                             'recomputed': sha, 'banked': man['spec_sha']}

    num = int(man['Z']['num'])
    den = int(man['Z']['den'])
    residues = {int(p): int(z) % int(p) for p, z in man['residues'].items()}

    # V2 banked fraction reduces to every banked residue
    bad = []
    for p, z in residues.items():
        if den % p == 0 or num % p * pow(den % p, p - 2, p) % p != z:
            bad.append(p)
    checks['V2_residues'] = {'pass': not bad, 'n_residues': len(residues),
                             'bad_primes': bad}

    # V3 independent CRT + ratrec over the banked residues
    ps = sorted(residues)
    r, M = _crt([residues[p] for p in ps], ps)
    fr = _ratrec(r, M)
    checks['V3_reconstruction'] = {
        'pass': fr is not None and fr == Fraction(num, den),
        'reconstructed': str(fr)}

    # V4 dense-oracle replay at fresh primes (cost-guarded)
    u = {k: int(v) for k, v in man['target_u'].items()}
    cells = dense_cells(spec, u)
    log(f"[verify] dense replay estimate: {cells} cells x {nfresh} primes "
        f"(max_cells={max_cells})")
    if cells > max_cells:
        checks['V4_fresh_replay'] = {
            'pass': False, 'refused': True, 'cells': cells,
            'detail': f'dense replay needs {cells} cells > '
                      f'max_cells={max_cells}; raise max_cells to force'}
    else:
        exclude = set(residues)
        fp = []
        while len(fp) < nfresh:
            cands = fresh_primes(nfresh - len(fp), exclude)
            exclude |= set(cands)
            fp += [p for p in cands if den % p]
        bad4 = []
        for p in fp:
            zo = core.Z_mod(spec, u, p)
            if num % p * pow(den % p, p - 2, p) % p != zo:
                bad4.append(p)
        checks['V4_fresh_replay'] = {'pass': not bad4, 'primes': fp,
                                     'bad_primes': bad4}

    # V5 recorded gates + target consistency
    gates = man.get('gates') or {}
    checks['V5_gates'] = {
        'pass': (bool(man.get('ok')) and bool(gates)
                 and all(g.get('pass') for g in gates.values())
                 and u == spec.target_u),
        'gates': {k: bool(g.get('pass')) for k, g in gates.items()},
        'target_u_matches_spec': u == spec.target_u}

    # V6 program hash (optional, when the plan is available)
    if plan is not None:
        from .planner import program_sha
        psha = program_sha(plan)
        checks['V6_program_sha'] = {
            'pass': (psha == man.get('program_sha')
                     and plan.get('spec_sha') == man['spec_sha']),
            'recomputed': psha, 'banked': man.get('program_sha')}

    ok = all(c['pass'] for c in checks.values())
    log(f"[verify] {man.get('name')}: {'PASS' if ok else 'FAIL'} "
        f"{ {k: c['pass'] for k, c in checks.items()} }")
    return {'ok': bool(ok), 'name': man.get('name'), 'checks': checks}

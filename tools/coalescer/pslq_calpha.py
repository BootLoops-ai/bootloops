#!/usr/bin/env python3
r"""pslq_calpha.py — PSLQ closure of the coalescence coefficients c_alpha through lockpick.pslq_gate,
in the sealrun (t6_close.py) form, for any family / mass tuple / ring.

The scan engine is lockpick.pslq_gate, imported by PATH from the sibling package (tools/lockpick/pslq_gate.py, or
--engine PATH).  Two gates before any number is read: (1) the engine FILE's sha256 must be listed in KNOWN_ENGINE_FILES
(an unlisted file is REFUSED rc 3 by name unless --allow-unlisted-engine is passed, which prints a loud UNLISTED
warning line); (2) its PROTOCOL units (DEFAULT, the three error classes, canonicalize, ndigits, check_digits,
capacity, _pools_from_strings, two_prec_stable, members_audit, graded_basket, reverify, _headroom, controls) are
hashed at every run (token hash, docstrings and comments stripped) and compared with the pin below; a mismatch is
REFUSED (rc 3) whether or not the file is listed.  The token pin sees only the 15 units: a statement placed after
them (a monkeypatch through globals()) is invisible to it, which is why the listed-file gate is the default.  An
unparseable engine file is REFUSED rc 3 by name.  The engine file's own sha256 is printed beside the pin.

Protocol (sealrun PROVENANCE / t6_close.py, re-instantiated for the coalescer's ring baskets):
  1. the basket is DECLARED from the ring spec (calpha_rings: family -> natural | extended | a json) and its
     member strings are built at run time; member names and values are written to the output;
  2. members_audit at scan dps BEFORE every scan — a nonempty audit is DegeneratePoolError by name (rc 3) and
     no NULL is reported;
  3. controls INSIDE the scan basket: the engine's synthetic (3*m0-7*m1)/5, a planted target' = (3*m_a - 7*m_b
     + 2*m_c)/5 (height 7, >= 10x below every ladder rung), and a matched-magnitude negative control
     (sha512-derived) that must NULL at the top rung — controls fail => rc 2, NO verdict of any kind;
  4. two-precision canonicalized scans (two_prec_stable) over the height ladder 1e4 .. maxcoeff, each rung
     capacity()-checked at the low leg (over-capacity rungs REFUSED and printed; no lawful rung at all =>
     CapacityError by name, rc 4); the dps pair is the sealrun plan (fit legs <= target_digits-40, >= 30 held-out
     digits) or --dps-pair LO,HI (the held-out law is still asserted);
  5. HIT => reverify at (target_digits-5) dps against digits no fit leg saw, pigeonhole line, canonical vector
     recorded; the closed-form EXPECTATION (calpha_rings: the log-form vector is DERIVED from the recorded structure,
     the additive-form vector is the scan of the closed form's own value) is compared: MATCH or MISMATCH by name;
  6. NULL through the ladder => exclusion floor with the capacity arithmetic (NULL-with-floor).
  7. the K3 (1,1,1,9) c_{3/2} control target runs in the SAME log basket (--k3-routeB / --k3-c) and must be
     CLOSED on the vector derived from -sqrt(3)/(36 pi).

Targets: --c ALPHA=STRING (repeatable), --from-routeB FILE (the Route B JSON of record; the midpoint string is
trimmed to floor(acc_bits*log10 2) digits, the trim-order law), --from-gate GATE.json (gate_compare's
"targets" block).  Forms: log (target log|c| over {1, log pi, log p, log Gamma(a)}) and additive (target
c*pi^k over {1, sqrt p, Gamma(a)^2, ...}); --form log,additive.

Controls on the command line (each must fail BY NAME, never a spurious relation):
  --planted-wrong NAME=EXPR   substitute a wrong constant for one ring member (e.g. 'logG(1/4)=log(gamma(1/3))')
  --add-member NAME=EXPR      add a member (a Q-dependent one, e.g. 'logG(3/4)=log(gamma(3/4))', is caught by
                              members_audit -> DegeneratePoolError, rc 3)
  --closed FILE.json          a wrong closed form -> EXPECTATION MISMATCH by name (rc 1)
  --engine PATH               an engine whose file sha is not listed, or whose protocol units differ from the pin,
                              is REFUSED rc 3 by name; --allow-unlisted-engine names the exception (the pin still gates)

Exit: 0 every target CLOSED (reverified; expectation matched where one exists)
      1 a NULL-with-floor, an expectation MISMATCH, or a HIT that failed reverify (all named in the output)
      2 controls failed (no verdict) / usage
      3 REFUSED by name: unlisted or unparseable engine file, engine pin mismatch, degenerate pool, missing input,
        closed form outside the ring
      4 DigitsError / CapacityError by name
"""
import argparse, hashlib, importlib.util, io, json, os, re, subprocess, sys, time, tokenize
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import calpha_rings as R   # noqa: E402

# ---- the pin: token hash of the engine's protocol units (unit_token_hashes below), re-derived at every run
PSLQ_GATE_PROTOCOL_SHA256 = 'ede58bdc5aed46be395f99155d1bd97a15249e611ae318d3f0d5aa4638b1b219'
PROTOCOL_UNITS = ['DEFAULT', 'DigitsError', 'CapacityError', 'DegeneratePoolError', 'canonicalize', 'ndigits',
                  'check_digits', 'capacity', '_pools_from_strings', 'two_prec_stable', 'members_audit',
                  'graded_basket', 'reverify', '_headroom', 'controls']
KNOWN_ENGINE_FILES = {   # file sha256 -> label: the FIRST gate (an unlisted file is refused unless --allow-unlisted-engine)
    '24f79bed90c9c166d922bfb23d0f54cb4164ed9d82c91540580338a869c0c864': 'tools/lockpick/pslq_gate.py (earlier pin)',
    '952b5a0dc75c80490b2d9d9a29b99d4e8140c970ca3576d73e1dbc35962c392b': 'tools/lockpick/pslq_gate.py (earlier pin)',
    '09c9afcf338265103935d9015b74745e0cf4c33cd19760637248601f6fc82f0b': 'tools/lockpick/pslq_gate.py (this repository)',
}
LADDER_DEFAULT = '1e4,1e5,1e6,1e8,1e10,1e12'
SELF_SHA256 = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
CALPHA_RINGS_SHA256 = hashlib.sha256(open(R.__file__, 'rb').read()).hexdigest()


class Refused(Exception):
    pass


class ControlsFailed(Exception):
    pass


# ---------------------------------------------------------------- engine binding + pin
def unit_token_hashes(path):
    import ast
    src = open(path).read(); tree = ast.parse(src); out = {}
    for n in tree.body:
        name = None
        if isinstance(n, (ast.FunctionDef, ast.ClassDef)): name = n.name
        elif isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name): name = n.targets[0].id
        if name not in PROTOCOL_UNITS: continue
        seg = ast.get_source_segment(src, n)
        has_doc = (isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.body and isinstance(n.body[0], ast.Expr)
                   and isinstance(n.body[0].value, ast.Constant) and isinstance(n.body[0].value.value, str))
        toks, skipped = [], False
        for t in tokenize.generate_tokens(io.StringIO(seg).readline):
            if t.type in (tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT,
                          tokenize.ENCODING, tokenize.ENDMARKER): continue
            if has_doc and t.type == tokenize.STRING and not skipped:
                skipped = True; continue
            toks.append(t.string)
        out[name] = hashlib.sha256(' '.join(toks).encode()).hexdigest()
    return out


def protocol_sha(path):
    h = unit_token_hashes(path)
    missing = [u for u in PROTOCOL_UNITS if u not in h]
    if missing: raise Refused(f'engine {path}: protocol units missing {missing}')
    return hashlib.sha256('\n'.join(f'{u}:{h[u]}' for u in PROTOCOL_UNITS).encode()).hexdigest()


def bind_engine(spec, allow_unlisted=False):
    if spec == 'lockpick': path = os.path.join(HERE, '..', 'lockpick', 'pslq_gate.py')
    else: path = spec
    path = os.path.abspath(path)
    if not os.path.exists(path): raise Refused(f'engine {spec!r}: no file at {path}')
    try:
        fsha = hashlib.sha256(open(path, 'rb').read()).hexdigest()
        psha = protocol_sha(path)
    except (SyntaxError, ValueError, OSError) as e:   # an unparseable / unreadable engine file: refused by name, never a traceback
        raise Refused(f'engine {path}: unparseable engine file ({type(e).__name__}: {e})')
    label = KNOWN_ENGINE_FILES.get(fsha, 'UNLISTED')
    pin_ok = psha == PSLQ_GATE_PROTOCOL_SHA256
    line = (f'ENGINE {path} file_sha256={fsha} ({label}) protocol_sha256={psha} '
            f'pin={PSLQ_GATE_PROTOCOL_SHA256[:16]}')
    if label == 'UNLISTED' and not allow_unlisted:     # gate 1: the file sha must be listed (the token pin sees only the 15 units)
        print(line + ' UNLISTED FILE REFUSED')
        raise Refused(f'engine file sha256 {fsha} not in KNOWN_ENGINE_FILES; protocol units '
                      f'{"match" if pin_ok else "DIFFER FROM"} the pin; pass --allow-unlisted-engine to run an unlisted engine '
                      f'(the protocol pin still gates; a statement after the 15 units is invisible to it)')
    if not pin_ok:                                       # gate 2: the protocol pin, listed or not
        print(line + ' PIN MISMATCH')
        raise Refused(f'pslq_gate protocol pin mismatch: engine {path} protocol_sha256={psha} != pinned '
                      f'{PSLQ_GATE_PROTOCOL_SHA256} (file sha256 {fsha}, {label})')
    print(line + ' PIN OK')
    if label == 'UNLISTED':
        print(f'WARNING: UNLISTED ENGINE FILE {path} (sha256 {fsha}) RUNS UNDER --allow-unlisted-engine: only its 15 protocol '
              f'units are pinned; any statement after them is unchecked. Do not bank a verdict from this run.')
    ms = importlib.util.spec_from_file_location('pslq_gate_engine', path)
    eng = importlib.util.module_from_spec(ms); ms.loader.exec_module(eng)
    for v in ([9, -4, 0], [-9, 4, 0], [2, 0, 4, 6], [0, 0]):
        if eng.canonicalize(v) != R.canonicalize(v):
            raise Refused(f'engine canonicalize({v}) = {eng.canonicalize(v)} != ring module {R.canonicalize(v)}')
    return eng, {'path': path, 'file_sha256': fsha, 'label': label, 'protocol_sha256': psha,
                 'pinned_protocol_sha256': PSLQ_GATE_PROTOCOL_SHA256, 'unlisted_allowed': bool(allow_unlisted)}


# ---------------------------------------------------------------- targets
def load_B_str(path, alpha):
    """Route B JSON (tolerant parser): (midpoint string, acc_bits) for alpha ('5/4' -> key '5//4')."""
    if not os.path.exists(path): raise Refused(f'Route B input missing: {path}')
    txt = open(path).read()
    akey = alpha.replace('/', '//')
    blk_m = re.search(r'"' + re.escape(akey) + r'":\s*\{([^}]+)\}', txt)
    if not blk_m: raise Refused(f'{path}: no block for alpha {alpha} (key {akey})')
    blk = blk_m.group(1)
    mid_m = re.search(r'"c_re_mid":\s*"\[([+\-0-9.eE]+)\s*(?:\+/-|±)', blk)
    acc_m = re.search(r'"c_acc_bits":\s*(\d+)', blk)
    if not (mid_m and acc_m): raise Refused(f'{path}: alpha {alpha}: c_re_mid / c_acc_bits not found')
    return mid_m.group(1), int(acc_m.group(1))


def trim_to_quotable(mid_str, acc_bits):
    """The trim-order law: the midpoint string cut to floor(acc_bits*log10 2) significant digits."""
    nd = int(acc_bits * 0.30102999566398120)
    with mp.workdps(nd + 20):
        return mp.nstr(mp.mpf(mid_str), nd), nd


def stamp():
    return subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'], capture_output=True, text=True).stdout.strip()


# ---------------------------------------------------------------- the sealed scan
def scan_target(eng, label, ts, names, members, lo, hi, rev, ladder, maxsteps, out_rungs):
    """Ladder scan; returns (hit_height, vector) or (None, None). Raises DegeneratePoolError (engine's) / Refused."""
    n = len(names)
    for h in ladder:
        cap = eng.capacity(lo, n, h)
        rung = {'height': h, 'capacity': {k: (round(v, 1) if isinstance(v, float) else v) for k, v in cap.items()}}
        if not cap['ok']:
            rung['skipped'] = 'over capacity at fit legs'
            out_rungs.append(rung)
            print(f'    rung {h:g}: REFUSED over capacity ({cap["required"]:.0f}d > {cap["available"]}d)')
            break
        v = eng.two_prec_stable(ts, names, members, dps_pair=(lo, hi), maxcoeff=int(h), maxsteps=maxsteps)
        rung['result'] = list(v) if v else None
        out_rungs.append(rung)
        print(f'    rung {h:g}: {"HIT " + str(tuple(v)) if v else "null"}')
        if v: return h, v
    return None, None


def run(a):
    t0 = time.time()
    eng, engine_rec = bind_engine(a.engine, a.allow_unlisted_engine)
    if a.family:
        spec = R.family_spec(a.family)
        if a.alphas: spec = R.spec_from_masses(spec['msq'], a.alphas.split(','), spec['pi_power'], spec['closed'], a.family)
    else:
        if not (a.masses and a.alphas): raise Refused('--masses needs --alphas')
        spec = R.spec_from_masses(a.masses.split(','), a.alphas.split(','), pi_power=a.pi_power or '1')
    if a.pi_power: spec['pi_power'] = R.frac_str(a.pi_power)
    fam = spec['family'] or 'masses'
    forms = [f.strip() for f in a.form.split(',') if f.strip()]
    for f in forms:
        if f not in ('log', 'additive'): raise Refused(f'--form {f!r}: log | additive')

    # ---- targets (alpha -> string), trim-order law on Route B inputs
    targets, sources = {}, {}
    for spec_c in a.c or []:
        al, s = spec_c.split('=', 1)
        targets[R.frac_str(al)] = s.strip(); sources[R.frac_str(al)] = 'command line'
    if a.from_routeB:
        for al in spec['alphas']:
            s, bits = load_B_str(a.from_routeB, al)
            ts, nd = trim_to_quotable(s, bits)
            targets[al] = ts; sources[al] = f'{os.path.basename(a.from_routeB)} (acc_bits {bits} -> {nd} d)'
    if a.from_gate:
        if not os.path.exists(a.from_gate): raise Refused(f'--from-gate input missing: {a.from_gate}')
        g = json.load(open(a.from_gate))
        for al, t in g.get('targets', {}).items():
            targets[al] = t['value']; sources[al] = f'{os.path.basename(a.from_gate)} targets[{al}] ({t.get("digits")} d)'
    if not targets: raise Refused('no target: give --c ALPHA=STRING, --from-routeB FILE or --from-gate GATE.json')
    controls_targets = {}
    if a.k3_routeB:
        s, bits = load_B_str(a.k3_routeB, '3/2')
        ts, nd = trim_to_quotable(s, bits)
        controls_targets['K3_control'] = ts; sources['K3_control'] = f'{os.path.basename(a.k3_routeB)} (acc_bits {bits} -> {nd} d)'
    if a.k3_c:
        controls_targets['K3_control'] = a.k3_c.strip(); sources['K3_control'] = 'command line'
    k3spec = R.family_spec('K3')
    k3cf = R.closed_form(k3spec, '3/2')

    # ---- closed-form expectations
    closed = {}
    for al in targets:
        cf = R.closed_form(spec, al, a.closed)     # 'auto' -> the recorded form or None; 'none'; or a json path
        if cf is not None: closed[al] = cf

    # ---- rings
    ring_json = R.ring_from_json(a.ring) if a.ring not in ('natural', 'extended') else None
    ladder = sorted({float(x) for x in a.ladder.split(',')} | {float(a.maxcoeff)})
    ladder = [int(h) for h in ladder if h <= float(a.maxcoeff)]
    ladder_add = sorted({float(x) for x in a.ladder.split(',')} | {float(a.maxcoeff_additive)})
    ladder_add = [int(h) for h in ladder_add if h <= float(a.maxcoeff_additive)]
    add_exprs = {}
    for spec_m in a.add_member or []:
        nm, ex = spec_m.split('=', 1); add_exprs[nm.strip()] = ex.strip()
    wrong = {}
    for spec_w in a.planted_wrong or []:
        nm, ex = spec_w.split('=', 1); wrong[nm.strip()] = ex.strip()

    out = {'producer': {'script': os.path.basename(__file__), 'script_sha256': SELF_SHA256,
                        'calpha_rings_sha256': CALPHA_RINGS_SHA256, 'engine': engine_rec, 'stamp_utc': stamp(),
                        'argv': sys.argv[1:], 'python': sys.version.split()[0], 'mpmath': mp.__version__},
           'family': fam, 'msq': list(spec['msq']), 'alphas': list(spec['alphas']), 'pi_power': spec['pi_power'],
           'ring': a.ring, 'forms': forms, 'targets': {}, 'sources': sources, 'scans': [], 'controls': {},
           'verdicts': {}}
    for al, ts in list(targets.items()) + list(controls_targets.items()):
        out['targets'][al] = {'digits': eng.ndigits(ts), 'head40': ts[:42]}
    print(f'PSLQ c_alpha closure: family {fam} msq {list(spec["msq"])} alphas {list(spec["alphas"])} ring {a.ring} forms {forms}')
    print(f'  targets: ' + ', '.join(f'{al} ({eng.ndigits(ts)} d, {sources[al]})' for al, ts in list(targets.items()) + list(controls_targets.items())))
    if wrong: print(f'  PLANTED-WRONG members: {wrong}')
    if add_exprs: print(f'  ADDED members: {add_exprs}')

    fails, refusals = [], []
    for form in forms:
        # basket
        if ring_json is not None:
            if form not in ring_json: raise Refused(f'--ring {a.ring}: no {form!r} ring')
            names = list(ring_json[form]) if form == 'log' else list(ring_json[form]['members'])
            pi_power = spec['pi_power'] if form == 'log' else ring_json[form]['pi_power']
        else:
            names = R.ring_names(spec, form, a.ring)
            pi_power = spec['pi_power']
        for nm in add_exprs:
            if nm not in names: names.append(nm)
        # per-form target strings
        form_targets = {}
        for al, ts in list(targets.items()) + ([] if form != 'log' else list(controls_targets.items())):
            nd = eng.ndigits(ts)
            if nd < 150: raise eng.DigitsError(f'target {al} carries only {nd} digits; refusing (need >= 150)')
            with mp.workdps(nd + 20):
                v = mp.mpf(ts)
                if form == 'log': w = mp.log(abs(v))
                else: w = v * mp.pi ** R._mpfrac(pi_power)
                form_targets[al] = (mp.nstr(w, nd, strip_zeros=False), nd)
        max_nd = max(nd for _, nd in form_targets.values())
        min_nd = min(nd for _, nd in form_targets.values())
        member_d = max_nd + 10
        members = {}
        with mp.workdps(member_d + 20):
            for nm in names:
                if nm in wrong: val = R.parse_expr(wrong[nm])
                elif nm in add_exprs: val = R.parse_expr(add_exprs[nm])
                else: val = R.member_value(nm)
                members[nm] = mp.nstr(val, member_d, strip_zeros=False)
        # dps plan (sealrun form) PER TARGET: fit legs <= target_digits-40, >= 30 held-out digits; an explicit
        # --dps-pair applies wherever the held-out law allows it, else that target's own plan is used and named
        def plan_for(nd):
            rev = min(nd - 5, member_d - 5)
            auto = (min(nd - 40, rev - 30) - 30, min(nd - 40, rev - 30), rev)
            if a.dps_pair == 'auto': return auto, 'sealrun plan'
            lo, hi = (int(x) for x in a.dps_pair.split(','))
            if lo < 20 or hi <= lo: raise Refused(f'--dps-pair {a.dps_pair}: not a usable pair')
            if rev - hi < 30:
                return auto, f'--dps-pair {a.dps_pair} violates the held-out law at {nd} d (reverify {rev} - hi {hi} < 30); sealrun plan used'
            return (lo, hi, rev), f'--dps-pair {a.dps_pair}'
        plans = {al: plan_for(nd) for al, (ts, nd) in form_targets.items()}
        top = ladder[-1] if form == 'log' else ladder_add[-1]
        rungs_here = ladder if form == 'log' else ladder_add
        print(f'\n=== form {form}: basket ({len(names)}) {names}; member strings {member_d} d; ladder {[f"{h:g}" for h in rungs_here]} ===')
        for al, ((lo, hi, rev), note) in plans.items():
            print(f'  plan {al}: fit legs ({lo},{hi}) d; reverify {rev} d; held-out {rev - hi} d [{note}]')
        # 2. members_audit at the highest scan leg
        hi_max = max(p[0][1] for p in plans.values())
        aud = eng.members_audit(members, names, dps=min(300, hi_max), h=1e5)
        print(f'  members_audit @ {min(300, hi_max)}d h=1e5: {"clean" if not aud else aud}')
        if aud:
            rel = aud[0]['relation']
            err = eng.DegeneratePoolError(rel, [n for n in names if rel.get(n)], residual=aud[0]['residual'], dps=aud[0]['dps'])
            print(f'  REFUSED: {type(err).__name__}: {err}')
            out['scans'].append({'form': form, 'names': names, 'members_audit': aud, 'refused': f'{type(err).__name__}: {err}'})
            refusals.append(f'{form}: {type(err).__name__} — no NULL reported')
            continue
        # 3. controls inside the scan basket, once per distinct dps plan
        out['controls'][form] = {}
        plan_groups = {}
        for (lo, hi, rev), _ in plans.values():
            plan_groups[(lo, hi)] = min(rev, plan_groups.get((lo, hi), rev))
        for (lo, hi), rev_dps in sorted(plan_groups.items()):
            proto = {**eng.DEFAULT, 'dps_pair': (lo, hi), 'maxcoeff': int(top)}
            known = []
            if len(names) >= 3:
                ia, ib = 1, 2
                ic = len(names) - 1 if len(names) - 1 not in (ia, ib) else 0   # three DISTINCT members (a 3-member basket takes m_0)
                with mp.workdps(rev_dps + 40):
                    synth = (3 * mp.mpf(members[names[ia]]) - 7 * mp.mpf(members[names[ib]]) + 2 * mp.mpf(members[names[ic]])) / 5
                    sstr = mp.nstr(synth, rev_dps + 30, strip_zeros=False)
                want = [5] + [0] * len(names); want[1 + ia] = -3; want[1 + ib] = 7; want[1 + ic] = -2
                known.append({'name': f'planted (3*{names[ia]} - 7*{names[ib]} + 2*{names[ic]})/5', 'target': sstr,
                              'members': names, 'vector': want})
            ctl = eng.controls(members, proto, known_hits=known)
            with mp.workdps(rev_dps + 40):
                seed = hashlib.sha512(form_targets[list(form_targets)[0]][0].encode()).hexdigest()
                neg = mp.mpf('0.' + seed.translate(str.maketrans('abcdef', '123456')) * 8)
                negs = mp.nstr(neg * mp.mpf(members[names[-1]]), rev_dps + 30, strip_zeros=False)
            nv = eng.two_prec_stable(negs, names, members, dps_pair=(lo, hi), maxcoeff=int(top), maxsteps=a.maxsteps)
            ctl['negative'] = {'control': f'sha512-derived x {names[-1]} at height {top:g}', 'got': nv, 'ok': nv is None}
            ctl['ok'] = ctl['ok'] and nv is None
            ctl['dps_pair'] = [lo, hi]
            print(f'  controls at fit legs ({lo},{hi}) d:')
            for d in ctl['detail'] + [ctl['negative']]:
                print(f'    control {d["control"]}: {"PASS" if d["ok"] else "FAIL"} got={d.get("got")}'
                      + (f' want={d["want"]}' if 'want' in d else '') + (f' [{d["headroom_warn"]}]' if d.get('headroom_warn') else ''))
            out['controls'][form][f'{lo},{hi}'] = ctl
            if not ctl['ok']:
                print('  CONTROLS FAILED — refusing any verdict for this form')
                raise ControlsFailed(f'{form} at fit legs ({lo},{hi}): controls failed')
        # 4./5./6. scans
        for al, (ts, nd) in form_targets.items():
            (lo, hi, rev_dps), plan_note = plans[al]
            label = f'{form}:{al}'
            print(f'  -- target {label} ({nd} d; fit legs ({lo},{hi}), reverify {rev_dps})')
            rec = {'form': form, 'alpha': al, 'target_digits': nd, 'names': names,
                   'dps_plan': {'fit_pair': [lo, hi], 'reverify_dps': rev_dps, 'held_out_digits': rev_dps - hi, 'note': plan_note},
                   'rungs': []}
            try:
                h, v = scan_target(eng, label, ts, names, members, lo, hi, rev_dps, rungs_here, a.maxsteps, rec['rungs'])
            except eng.DegeneratePoolError as e:
                print(f'  REFUSED: {type(e).__name__}: {e}')
                rec['refused'] = f'{type(e).__name__}: {e}'
                out['scans'].append(rec); refusals.append(f'{label}: DegeneratePoolError'); continue
            # expectation
            exp_vec, exp_note = None, None
            cf = k3cf if al == 'K3_control' else closed.get(al)
            if cf is not None:
                if form == 'log':
                    try:
                        exp_vec = cf.log_vector(names); exp_note = f'derived from the closed form {cf.display}'
                    except R.NotRingMember as e:
                        raise Refused(f'closed form for {al} outside the ring: {e}')
                else:
                    with mp.workdps(nd + 20):
                        cstr = mp.nstr(cf.value() * mp.pi ** R._mpfrac(pi_power), nd, strip_zeros=False)
                    ev = None
                    for hh in rungs_here:
                        if not eng.capacity(lo, len(names), hh)['ok']: break
                        ev = eng.two_prec_stable(cstr, names, members, dps_pair=(lo, hi), maxcoeff=int(hh), maxsteps=a.maxsteps)
                        if ev: break
                    exp_vec = tuple(ev) if ev else None
                    exp_note = f'scan of the closed form value {cf.display} x pi^{pi_power}' + ('' if ev else ' -> NULL (closed form not in this additive ring at these heights)')
            if v:
                rv = eng.reverify({'target': ts, 'members': names, 'vector': list(v)}, members, dps=rev_dps, tol=f'1e-{rev_dps - 60}')
                maxc = max(abs(int(c)) for c in v)
                with mp.workdps(50):
                    pfloor = -(len(names) * mp.log10(mp.mpf(maxc))) if maxc > 1 else mp.mpf(0)
                rec['HIT'] = {'height_rung': h, 'vector': list(v), 'vector_named': dict(zip(['T'] + names, v)),
                              'height': maxc, 'reverify': rv, 'pigeonhole_floor_log10': float(pfloor)}
                verdict = 'CLOSED' if rv['ok'] else 'HIT-UNVERIFIED (reverify failed — treat as NULL-class)'
                print(f'    HIT at rung {h:g}: {dict(zip(["T"] + names, v))} height {maxc}; reverify @{rev_dps}d relresid {rv["relresid"]} {"ok" if rv["ok"] else "FAILED"}; pigeonhole floor 1e{float(pfloor):.0f}')
            else:
                hmax = max((r['height'] for r in rec['rungs'] if 'result' in r), default=0)
                if hmax == 0:
                    raise eng.CapacityError(f'{label}: no lawful rung at fit legs ({lo},{hi}) d for {len(names)} members '
                                            f'(rung {rungs_here[0]:g} needs {eng.capacity(lo, len(names), rungs_here[0])["required"]:.0f} d)')
                rec['NULL'] = {'exclusion_height': hmax,
                               'statement': (f'no integer relation with |coeff| <= {hmax:g} links {label} to the {len(names)}-member '
                                             f'basket at fit legs ({lo},{hi})d (controls green, audit clean); capacity-checked per rung')}
                verdict = 'NULL-with-floor'
                print(f'    NULL-with-floor: {rec["NULL"]["statement"]}')
            if cf is not None:
                rec['expectation'] = {'vector': list(exp_vec) if exp_vec else None, 'note': exp_note}
                if exp_vec is not None and v and tuple(v) == tuple(exp_vec):
                    rec['expectation']['match'] = True
                    print(f'    EXPECTATION MATCH: {exp_note} -> {tuple(exp_vec)}')
                else:
                    rec['expectation']['match'] = False
                    verdict = f'EXPECTATION MISMATCH ({verdict})'
                    print(f'    EXPECTATION MISMATCH: expected {tuple(exp_vec) if exp_vec else None} ({exp_note}), got {tuple(v) if v else None}')
            rec['verdict'] = verdict
            out['verdicts'][label] = verdict
            out['scans'].append(rec)
            if verdict != 'CLOSED': fails.append(f'{label}: {verdict}')

    out['wall_s'] = round(time.time() - t0, 2)
    out['fails'] = fails; out['refusals'] = refusals
    if refusals: out['status'] = 'REFUSED'
    elif fails: out['status'] = 'FAIL'
    else: out['status'] = 'CLOSED'
    json.dump(out, open(a.out, 'w'), indent=1, default=str)
    print(f'\nwrote {a.out}')
    print(f'PSLQ_VERDICT family={fam} ring={a.ring} forms={",".join(forms)} status={out["status"]} '
          f'verdicts={json.dumps(out["verdicts"])}' + (f' refusals={refusals}' if refusals else '')
          + (f' fails={fails}' if fails else '') + f' wall_s={out["wall_s"]}')
    print(f'PRODUCER pslq_calpha.py sha256={SELF_SHA256[:16]} calpha_rings.py sha256={CALPHA_RINGS_SHA256[:16]} '
          f'engine={engine_rec["file_sha256"][:16]} protocol={engine_rec["protocol_sha256"][:16]}')
    if refusals: return 3
    return 0 if not fails else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--family', help=f'one of {R.families()} (or --masses + --alphas)')
    ap.add_argument('--masses'); ap.add_argument('--alphas')
    ap.add_argument('--c', action='append', help='ALPHA=STRING target (repeatable)')
    ap.add_argument('--from-routeB', help='Route B JSON of record (midpoint trimmed to the certified digits)')
    ap.add_argument('--from-gate', help="GATE.json written by gate_compare (its 'targets' block)")
    ap.add_argument('--k3-routeB', help='Route B JSON of the K3 control (scanned in the same log basket)')
    ap.add_argument('--k3-c', help='K3 control target string')
    ap.add_argument('--ring', default='natural', help='natural | extended | FILE.json')
    ap.add_argument('--form', default='log,additive')
    ap.add_argument('--pi-power', help='additive-form pi power override (family default: K3 1, CY3 5/2, CY4 2)')
    ap.add_argument('--closed', default='auto', help='auto | none | FILE.json (expectations)')
    ap.add_argument('--dps-pair', default='auto', help='auto (sealrun plan) | LO,HI')
    ap.add_argument('--maxcoeff', default='1e6', help='ladder top for the log form')
    ap.add_argument('--maxcoeff-additive', default='1e8', help='ladder top for the additive form')
    ap.add_argument('--ladder', default=LADDER_DEFAULT)
    ap.add_argument('--maxsteps', type=int, default=120000)
    ap.add_argument('--engine', default='lockpick', help='lockpick | PATH to a pslq_gate.py')
    ap.add_argument('--allow-unlisted-engine', action='store_true',
                    help='run an engine file whose sha256 is not in KNOWN_ENGINE_FILES (refused by default; the protocol pin still gates; a loud warning is printed)')
    ap.add_argument('--planted-wrong', action='append', help='NAME=EXPR control: wrong constant for a member')
    ap.add_argument('--add-member', action='append', help='NAME=EXPR: add a member (degenerate-basket control)')
    ap.add_argument('--out', default='PSLQ_CALPHA.json')
    a = ap.parse_args(argv)
    if not a.family and not a.masses: ap.error('--family or --masses required')
    try:
        return run(a)
    except Refused as e:
        print(f'REFUSED: {e}'); return 3
    except (R.NotRingMember, R.RingSpecError) as e:
        print(f'REFUSED: {type(e).__name__}: {e}'); return 3
    except ControlsFailed as e:
        print(f'CONTROLS FAILED: {e} — no verdict'); return 2
    except ValueError as e:      # the engine's DigitsError / CapacityError (and only those by name)
        if type(e).__name__ in ('DigitsError', 'CapacityError'):
            print(f'{type(e).__name__}: {e}'); return 4
        raise


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
r"""
shift_opt.py — per-sector loop-momentum-shift optimizer for AMFlow η-DE walls.

WHY. AMFlow's η-DE reduction size is dominated by how the loop-momentum
routing couples the loops in the *present* propagators of a sector. A
unimodular affine shift  (l,k1,k2) → M·(l,k1,k2) + Σcᵢpᵢ  leaves the
physical zero-ISP integral value invariant (Jacobian=1) but can decouple
loops in the sector's η-deformed IBP system by orders of magnitude.

Measured on a 3-loop 15-propagator reference family:
  k1sh (k1→l+p1−k1) shrinks sec1018 η-DE 815k→86k (9.4×).
  k2sh (k2→l+p1−k2) shrinks sec510 537k→99k (5.4×), sec1014 3.0×, sec1018 2.4× further.
  Direct 2-var Fermat: k2sh sec510 items 25k→29 in 487s vs orig ETA 60+h (>200×).
  All shifts gated ≥168d vs independent oracle at sec763 control.

HOW.
  1. Shift D1..D_{n_prop} by direct sympy substitution.
  2. Re-choose ISPs (D_{n_prop+1}..D15) as the simplest quadratic forms
     completing the rank-15 loop-bilinear span (never shift ISPs directly:
     the shifted form is complex and irrelevant for zero-ISP integrals).
  3. Emit an amflow input JSON per (shift,sector).
  4. [--probe] Launch under the Fermat-mode env, poll η-masters + Fermat
     items₀, write SHIFT_PROBE.json. 8-min cap per probe.

USAGE:
  # build shifted-family JSONs from a template amflow input
  shift_opt.py build TEMPLATE.json --shift 'k2sh:k2=l+p1-k2' \
      --sec 510:0,2,1,1,1,1,1,1,1,0 --out DIR

  # enumerate the standard shift catalog for a family and probe
  shift_opt.py probe TEMPLATE.json --sec 1022:0,1,2,1,1,1,1,1,1,1 \
      --catalog --out DIR --amf /path/to/amflow_cli

  # standard catalog (auto-generated for a 3-loop family with loops L=[l0,l1,l2]):
  #   {li}sh:    li → l0+P−li          (P = external momentum in D1)
  #   {li}{lj}sh:  both
  #   lshp:      l0 → l0−P
  #   composites lshp∘{li}sh
  # A general GL(3,ℤ) enumeration is NOT attempted: the affine 1-shifts above
  # cover the "decouple one loop from another in present props" mechanism,
  # which is empirically the load-bearing structure.

CAVEATS:
  - Zero-ISP integrals ONLY are shift-invariant. If your target has ISP
    numerators, the value changes and you must re-express in the shifted
    basis (not implemented here).
  - The rank check is over the 15 loop-bilinear scalar products; ext-ext
    products are kinematic constants and dropped.
  - Fermat-mode env (DCAP=1 RUN_FIREFLY=0 PREHEAT_DOT=0) is the FAST route
    for shifted families; FF 2-var (d,η) reconstruction is ~10× slower here.
"""
import json, os, sys, argparse, subprocess, time, re, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "amflow-kit")); from amflow_kit import memfence  # noqa: E402 — tools/amflow-kit/amflow_kit/memfence.py, same hop in the repo and a mirror tree laid out as tools/<pkg>/
import sympy as sp


def _pid_alive(pid):
    """True if pid is a live process. Signal-0 probe — portable (/proc is
    Linux-only); PermissionError means alive but foreign-owned."""
    try:
        os.kill(int(pid), 0)
    except (ProcessLookupError, TypeError, ValueError):
        return False
    except PermissionError:
        pass
    return True

# ---------------------------------------------------------------------------
def _syms(loops, exts):
    return {n: sp.Symbol(n) for n in loops+exts}

def _sp_basis(loops, exts):
    B=[]
    for i,a in enumerate(loops):
        for b in loops[i:]: B.append((a,b))
    for a in loops:
        for b in exts: B.append((a,b))
    return B

def _mom_vec(expr_str, syms, keys):
    e = sp.sympify(expr_str, locals=syms)
    poly = sp.Poly(e, *[syms[n] for n in keys], domain='QQ')
    return {n: poly.coeff_monomial(syms[n]) for n in keys}

def _sq_row(prop_str, syms, loops, exts, spb):
    s = prop_str.strip()
    for tail in ('- msq','-msq','- m2','-m2'):
        if s.endswith(tail): s = s[:-len(tail)].strip(); break
    if not (s.startswith('(') and s.endswith(')^2')):
        raise ValueError(f"bad prop: {prop_str!r}")
    vec = _mom_vec(s[1:-3], syms, loops+exts)
    row = {ab:0 for ab in spb}
    for a in loops+exts:
        for b in loops+exts:
            c = vec[a]*vec[b]
            if c==0 or (a in exts and b in exts): continue
            k = (a,b) if (a,b) in row else (b,a)
            row[k] += c
    return [sp.Rational(row[ab]) for ab in spb]

def apply_shift(prop_str, shift, syms):
    s = prop_str.strip(); tail=''
    for t in ('- msq','-msq','- m2','-m2'):
        if s.endswith(t): s,tail = s[:-len(t)].strip(),' - '+t.lstrip('- '); break
    mom = sp.sympify(s[1:-3], locals=syms)
    subs = {syms[k]: sp.sympify(v, locals=syms) for k,v in shift.items()}
    return f"({sp.expand(mom.subs(subs, simultaneous=True))})^2{tail}"

def complete_isps(props_phys, loops, exts, syms, spb):
    M = sp.Matrix([_sq_row(p, syms, loops, exts, spb) for p in props_phys])
    cands = ([f"({a} + {b})^2" for a in loops for b in exts] +
             [f"({a})^2" for a in loops] +
             [f"({loops[i]} + {loops[j]})^2" for i in range(len(loops)) for j in range(i+1,len(loops))] +
             [f"({a} + {b1} + {b2})^2" for a in loops for i,b1 in enumerate(exts) for b2 in exts[i+1:]])
    isps=[]
    for c in cands:
        M2 = M.col_join(sp.Matrix([_sq_row(c, syms, loops, exts, spb)]))
        if M2.rank() > M.rank():
            M = M2; isps.append(c)
        if M.rank()==len(spb): break
    if M.rank()!=len(spb):
        raise RuntimeError(f"rank {M.rank()}/{len(spb)}; ISP completion failed")
    return isps

def build_family(tmpl_fam, shift, tag, n_phys=None):
    loops, legs = tmpl_fam['loops'], tmpl_fam['legs']
    # eliminate conserved leg
    exts = [p for p in legs if p not in tmpl_fam.get('conservation',{})]
    syms = _syms(loops, exts)
    spb  = _sp_basis(loops, exts)
    props = tmpl_fam['propagators']
    n_phys = n_phys or (len(props) - (len(spb)-len(props)) if len(props)<len(spb) else len(spb))
    # heuristic: physical props are those before ISPs; caller should pass n_phys
    if n_phys is None or n_phys<=0 or n_phys>len(props): n_phys = len(props)
    phys = [apply_shift(p, shift, syms) for p in props[:n_phys]]
    isps = complete_isps(phys, loops, exts, syms, spb)
    fam = dict(tmpl_fam)
    fam['name'] = f"{tmpl_fam['name']}_{tag}"
    fam['propagators'] = phys + isps
    return fam

def std_catalog(loops, P='p1'):
    """Standard affine 1-shift catalog for a 3-loop family."""
    l0 = loops[0]; sub = loops[1:]
    cat = {}
    for k in sub:
        cat[f'{k}sh'] = {k: f'{l0} + {P} - {k}'}
    cat['lshp'] = {l0: f'{l0} - {P}'}
    for k in sub:
        cat[f'lsh{k[0].upper()}'] = {l0: f'{l0} + {k} - {P}'}
        cat[f'{k}shp'] = {l0: f'{l0} - {P}', k: f'{l0} - {k}'}
    if len(sub)>=2:
        cat['k12sh']  = {sub[0]: f'{l0} + {P} - {sub[0]}', sub[1]: f'{l0} + {P} - {sub[1]}'}
        cat['k12shp'] = {l0: f'{l0} - {P}', sub[0]: f'{l0} - {sub[0]}', sub[1]: f'{l0} - {sub[1]}'}
    return cat

# ---------------------------------------------------------------------------
FERMAT_ENV = dict(
    AMFLOW_REDUCE_DCAP='1', AMFLOW_RUN_FIREFLY='0', AMFLOW_NO_AUTO_FIREFLY='1',
    AMFLOW_PREHEAT_DOT='0', AMFLOW_KIRA_MEM_CAP_GB='15', OPENBLAS_NUM_THREADS='2')

def probe(inp, out_json, log, work_dir, amf, wrap, cache, extra_env=None):
    env = dict(os.environ, **FERMAT_ENV, AMFLOW_IBP_CACHE=cache)
    if 'FERMATPATH' not in env:
        env.setdefault('FERMATPATH', 'fer64')  # fer64 on PATH or pre-set FERMATPATH
    if extra_env: env.update(extra_env)
    # memfence: FERMAT_ENV's cap is the memory BUDGET (15 GiB). Inside a
    # finite-memory.max cgroup leaf the helper exports '0' (no fuse; the leaf
    # is the fence); otherwise the 2.5 x budget ADDRESS-SPACE fuse (37.5 GiB)
    # the binary applies as RLIMIT_AS to its Kira child — never the 1x budget.
    memfence.apply_policy_to_env(env, memfence.fuse_policy(env.get(memfence.CAP_ENV)))
    env[memfence.LAUNCH_ID_ENV] = memfence.new_launch_id()  # this probe's identity
    if os.path.isdir(work_dir):
        import shutil; shutil.rmtree(work_dir)
    os.makedirs(work_dir, exist_ok=True)
    t0 = time.time()
    with open(log,'w') as lf:
        p = subprocess.Popen(['setsid', wrap, amf, inp, out_json],
                             env=env, stdout=lf, stderr=subprocess.STDOUT,
                             stdin=subprocess.DEVNULL, start_new_session=True)
    # env readback from /proc: setsid(1) forks and exits at once (p.pid is the
    # fork parent) and the pre-exec stages (setsid, the wrapper) still carry
    # the launcher's env, so resolve the EXEC'D amflow: argv tail
    # [amf, inp, out_json] + start time + THIS probe's nonce in its environ
    # (argv + time alone is ambiguous on a re-run), skipping cmdlines that
    # still name the wrapper or setsid. The readback is written beside the
    # log as <log>.memfence.json (shift_opt has no per-probe receipt); a
    # mismatch aborts BY THAT PID, never a pattern, and only if the pid
    # still carries the nonce.
    nonce = {memfence.LAUNCH_ID_ENV: env[memfence.LAUNCH_ID_ENV]}
    intended = memfence.intended_from_env(env, nonce)
    rb = memfence.env_readback(
        memfence.find_spawned_pid([amf, inp, out_json], t0, window=15.0,
                                  require_env=nonce, exclude=(wrap, 'setsid')),
        intended)
    with open(log + '.memfence.json', 'w') as fh:
        json.dump({'argv': [wrap, amf, inp, out_json], 'direct_pid': p.pid,
                   'fence_mode': env.get(memfence.MODE_ENV),
                   'cap_export': env.get(memfence.CAP_ENV), 'readback': rb},
                  fh, indent=1, default=str)
    if rb['diff']:
        ab = memfence.abort_by_name(rb['pid'], [amf, inp, out_json], reason='env readback mismatch',
                                    require_env=nonce)
        raise RuntimeError(f"memfence: env readback mismatch for pid {rb['pid']}: {rb['diff']} "
                           f"(killed by pid: gone={ab['gone']}, refused={ab['refused']}); see {log}.memfence.json")
    return p.pid

def harvest(work_dir):
    """Read η-masters, seed, mandatory, Fermat items from a probe workdir."""
    d={}
    ph = glob.glob(f'{work_dir}/part_0/amf/system_0_diffeq/masters_preheat/results/*/masters')
    if ph: d['n_eta_masters'] = sum(1 for _ in open(ph[0]))
    klog = f'{work_dir}/part_0/amf/system_0_diffeq/target_reduce/kira.log'
    if os.path.exists(klog):
        txt = open(klog).read()
        m = re.search(r'rmax: (\d+) smax: (\d+) dmax: (\d+)', txt)
        if m: d['seed']=f"r{m.group(1)}s{m.group(2)}d{m.group(3)}"
        m = re.search(r'mandatory list: (\d+)', txt)
        if m: d['mandatory']=int(m.group(1))
        m = re.search(r'Loading (\d+) equations', txt)
        if m: d['n_eqns_ff']=int(m.group(1))
        its = re.findall(r'iteration: (\d+); items left: (\d+); time: \( ([\d.]+) s \)', txt)
        if its:
            d['fermat_items0']=int(its[0][1]); d['fermat_items_last']=int(its[-1][1])
            d['fermat_t_last']=float(its[-1][2]); d['fermat_iters']=int(its[-1][0])
    return d

# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('verb', choices=['build','probe','harvest'])
    ap.add_argument('template', help='template amflow input JSON (family + options)')
    ap.add_argument('--shift', action='append', default=[],
                    help='tag:loop=expr[,loop=expr] (repeatable)')
    ap.add_argument('--catalog', action='store_true', help='use std_catalog shifts')
    ap.add_argument('--sec', action='append', default=[],
                    help='sec:idx1,idx2,...  (D1..Dn indices, ISP zeros appended)')
    ap.add_argument('--n-phys', type=int, help='# physical (non-ISP) propagators')
    ap.add_argument('--out', default='.', help='output dir')
    ap.add_argument('--amf', default=os.environ.get('AMFLOW_CLI', 'amflow_cli'))
    ap.add_argument('--wrap', default=os.environ.get('AMFLOW_WRAP', 'run_amflow_jemalloc.sh'))
    ap.add_argument('--cap', type=int, default=480, help='per-probe wall cap (s)')
    args = ap.parse_args()

    tmpl = json.load(open(args.template))
    fam0 = tmpl['family']
    loops = fam0['loops']
    n_isp = len(_sp_basis(loops, [p for p in fam0['legs'] if p not in fam0.get('conservation',{})])) - (args.n_phys or 10)

    shifts = {}
    if args.catalog:
        shifts.update(std_catalog(loops))
    for s in args.shift:
        tag, spec = s.split(':',1)
        shifts[tag] = dict(kv.split('=') for kv in spec.split(','))

    secs = {}
    for s in args.sec:
        sid, idx = s.split(':',1)
        v = [int(x) for x in idx.split(',')]
        v += [0]*(len(fam0['propagators'])-len(v))
        secs[int(sid)] = v

    os.makedirs(args.out, exist_ok=True)
    written=[]; pids={}
    for tag, sh in shifts.items():
        fam = build_family(fam0, sh, tag, n_phys=args.n_phys)
        for sec, idx in secs.items():
            job = dict(tmpl)
            job['family'] = fam
            job['integrals'] = [{'indices': idx}]
            job['work_dir'] = os.path.join(args.out, f'work_{tag}_sec{sec}')
            fn = os.path.join(args.out, f'in_{tag}_sec{sec}.json')
            json.dump(job, open(fn,'w'), indent=2)
            written.append(fn)
            print(f'[build] {tag} sec{sec} -> {fn}')
            if args.verb=='probe':
                pid = probe(fn, os.path.join(args.out,f'out_{tag}_sec{sec}.json'),
                            os.path.join(args.out,f'log_{tag}_sec{sec}.log'),
                            job['work_dir'], args.amf, args.wrap,
                            os.path.join(args.out,f'cache_{tag}'))
                pids[(tag,sec)] = pid
                print(f'[probe] {tag} sec{sec} pid={pid}')

    if args.verb=='probe':
        t0=time.time(); res={}
        while pids and time.time()-t0 < args.cap:
            time.sleep(20)
            for (tag,sec),pid in list(pids.items()):
                wd = os.path.join(args.out, f'work_{tag}_sec{sec}')
                d = harvest(wd)
                if 'fermat_items0' in d or not _pid_alive(pid):
                    res[f'{tag}_sec{sec}']=d
                    print(f'[probe] {tag} sec{sec}: {d}')
        # cap: kill + harvest remaining
        for (tag,sec),pid in pids.items():
            try:
                pgid = os.getpgid(pid)
                os.killpg(pgid, 15); time.sleep(0.5); os.killpg(pgid, 9)
            except Exception: pass
            res.setdefault(f'{tag}_sec{sec}', harvest(os.path.join(args.out,f'work_{tag}_sec{sec}')))
        json.dump(res, open(os.path.join(args.out,'SHIFT_PROBE.json'),'w'), indent=2)
        print(f'[probe] wrote SHIFT_PROBE.json ({len(res)} entries)')

    if args.verb=='harvest':
        res={}
        for tag in shifts:
            for sec in secs:
                res[f'{tag}_sec{sec}'] = harvest(os.path.join(args.out,f'work_{tag}_sec{sec}'))
        json.dump(res, open(os.path.join(args.out,'SHIFT_PROBE.json'),'w'), indent=2)
        print(json.dumps(res, indent=2))

if __name__=='__main__':
    main()

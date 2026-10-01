#!/usr/bin/env python3
r"""quad_probe.py - per-chain quadrature convergence prober (probe BEFORE farming).

A multi-D fixed-node quadrature is a product of per-chain rules, each with
its own level knob (env-var levels like LEVEL_OUT / LEVEL_K).
Varying ONE chain's level with
the others frozen cancels every other chain's error EXACTLY in consecutive
differences, so the diff sequence measures that chain's convergence in
isolation.  Run at LOW dps (25-30) on 1-2 sample points BEFORE any farm; all
level configs are evaluated in PARALLEL.  A chain gaining < 5 digits/level
is flagged METHOD-MISMATCHED; standing advice:
an integrand factor that is exactly a classical weight (1/sqrt at the
endpoints = Chebyshev/Jacobi) wants the matching Gauss rule, not tanh-sinh
-- tools/gatekeeper/probes/ell_subst.py selects it.

Worked example (a measured production case):
a farm's inner K-chain was tanh-sinh between the two roots of a
concave quadratic P2 with a 1/sqrt(P0 P1 P2 P6) integrand; P0,P1 vanish just
outside the interval (pinching pair throttles ts).  Measured: L_K 4->5 rel
diff 3.0e-14, 5->6 rel diff 3.7e-17  =>  ~3 digits/level  =>
METHOD-MISMATCHED (a 16-17d farm ceiling; an 80-pt farm + a 24-shard re-farm
were wasted before a 3-minute probe exposed it).  1/sqrt(P2) IS the
Chebyshev weight.
Call pattern:
  python3 tools/gatekeeper/probes/quad_probe.py --workdir <dir with your farm module> \
    --module my_farm --func my_moments --args="-3,-5,0" \
    --base DPS=26,LEVEL_OUT=1,LEVEL_K=4 --chains LEVEL_K,LEVEL_OUT \
    --span 2 --extract "r[0][0]" --jobs 8
"""
import os, sys, json, argparse, subprocess
import mpmath as mp
ADVICE = ("  ==> METHOD-MISMATCHED (<5 digits/level): change METHOD, not level."
          "\n      Classical-weight factor (1/sqrt at endpoints = Chebyshev/"
          "Jacobi) -> matching\n      Gauss rule: "
          "tools/gatekeeper/probes/ell_subst.py.")


def probe(evaluator, base_levels, chains=None, span=2, noise_dps=None):
    """evaluator(levels_dict) -> list of mpf (or has .many(list_of_dicts)).
    Vary each chain over base..base+span (others at base); returns
    {chain: {levels, diffs, rates(digits/level), at_noise_floor, mismatched}}."""
    chains = list(chains or base_levels)
    cfgs, idx = [], {}
    for ch in chains:
        for lv in range(int(base_levels[ch]), int(base_levels[ch]) + span + 1):
            cfg = dict(base_levels); cfg[ch] = lv
            k = tuple(sorted(cfg.items()))
            if k not in idx:
                idx[k] = len(cfgs); cfgs.append(cfg)
    vals = (evaluator.many(cfgs) if hasattr(evaluator, 'many')
            else [evaluator(c) for c in cfgs])
    floor = mp.mpf(10)**(-(noise_dps or mp.mp.dps) + 2)
    rep = {}
    for ch in chains:
        lvls = list(range(int(base_levels[ch]), int(base_levels[ch]) + span + 1))
        vv = []
        for lv in lvls:
            cfg = dict(base_levels); cfg[ch] = lv
            vv.append(vals[idx[tuple(sorted(cfg.items()))]])
        sc = [max(abs(x), mp.mpf('1e-300')) for x in vv[-1]]
        diffs = [max(abs(x - y)/s for x, y, s in zip(va, vb, sc))
                 for va, vb in zip(vv, vv[1:])]
        rates = [float(mp.log10(d0/d1)) if d1 > 0 else float('inf')
                 for d0, d1 in zip(diffs, diffs[1:])]
        meas = [r for d, r in zip(diffs[1:], rates) if d > floor]
        rep[ch] = dict(levels=lvls, value=mp.nstr(vv[-1][0], 20),
                       diffs=[mp.nstr(d, 3) for d in diffs],
                       rates=[round(r, 2) for r in rates],
                       at_noise_floor=not meas,
                       mismatched=bool(meas) and min(meas) < 5.0)
    return rep


class EnvEvaluator:
    """module.func(*args) in a fresh subprocess per config; chain levels are
    passed as ENV VARS (covers the nodes-built-at-import pattern).  .many()
    runs configs in parallel (<= jobs)."""
    def __init__(self, workdir, module, func, args='', extract='r', jobs=8):
        self.wd, self.mod, self.func = workdir, module, func
        self.args, self.extract, self.jobs = args, extract, jobs

    def _code(self):
        return ("import os, sys, mpmath as mp\n"
                f"sys.path.insert(0, {self.wd!r})\n"
                f"import {self.mod} as M\n"
                f"r = M.{self.func}({self.args})\n"
                f"v = ({self.extract})\n"
                "def fl(x):\n"
                "    try: return [z for e in x for z in fl(e)]\n"
                "    except TypeError: return [x]\n"
                "print('QPROBE: ' + ' '.join(mp.nstr(z, mp.mp.dps)"
                " for z in fl(v)), flush=True)\n")

    def many(self, cfgs):
        procs = []
        for i, cfg in enumerate(cfgs):
            env = dict(os.environ, **{k: str(v) for k, v in cfg.items()})
            while sum(p.poll() is None for _, _, p in procs) >= self.jobs:
                next(p for _, _, p in procs if p.poll() is None).wait()
            procs.append((i, cfg, subprocess.Popen(
                [sys.executable, '-c', self._code()], env=env, cwd=self.wd,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)))
        out = [None]*len(cfgs)
        for i, cfg, p in procs:
            so, se = p.communicate()
            lines = [l for l in so.splitlines() if l.startswith('QPROBE: ')]
            if p.returncode or not lines:
                raise RuntimeError(f"eval failed for {cfg}:\n{se[-800:]}")
            out[i] = [mp.mpf(t) for t in lines[-1][8:].split()]
        return out

    def __call__(self, cfg):
        return self.many([cfg])[0]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for o in ('--workdir', '--module', '--func'):
        ap.add_argument(o, required=True)
    ap.add_argument('--args', default='')
    ap.add_argument('--base', required=True, help='NAME=int,... incl. DPS')
    ap.add_argument('--chains', required=True, help='subset of base to vary')
    ap.add_argument('--span', type=int, default=2)
    ap.add_argument('--extract', default='r', help='expr in r, e.g. r[0][0]')
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--json', default=None)
    a = ap.parse_args(argv)
    base = {k: int(v) for k, v in (kv.split('=') for kv in a.base.split(','))}
    mp.mp.dps = base.get('DPS', 28) + 10
    ev = EnvEvaluator(os.path.abspath(a.workdir), a.module, a.func, a.args,
                      a.extract, a.jobs)
    rep = probe(ev, base, a.chains.split(','), a.span,
                noise_dps=base.get('DPS', 28))
    bad = 0
    for ch, r in rep.items():
        print(f"[quad_probe] chain {ch}: levels {r['levels']} "
              f"(others at base)\n  consecutive max-rel diffs: {r['diffs']}\n"
              f"  digits/level: {r['rates']}"
              + ('  (diffs at dps noise floor - converged, rate unmeasurable;'
                 ' healthy)' if r['at_noise_floor'] else ''))
        if r['mismatched']:
            bad += 1
            print(ADVICE)
    if a.json:
        json.dump(rep, open(a.json, 'w'), indent=1)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())

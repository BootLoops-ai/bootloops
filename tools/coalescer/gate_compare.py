#!/usr/bin/env python3
"""Two-method gate + closed-form verification for the coalescence coefficients c_alpha.

Reads Route A (mpmath/s-coord) and Route B (Arb/z-coord) outputs for a family
(or a mass tuple + exponent list), compares every sample against the closed
forms from calpha_rings (evaluated at run time), computes s_dec-independence
per route and the cross-route gate (min agreeing digits over all A x B
pairs), runs the K3 positive control on its own Route A/B files when given,
and writes GATE.json.  The printed lines keep the form of the CY3 record
(the gate log of record, tests/fixtures/gate.log): on the CY3 inputs of record the
output is byte-identical up to the appended TARGET / GATE_VERDICT / PRODUCER lines.

Usage:
    gate_compare.py --family CY3 --routeA ROUTEA_CY3.json \
        --routeB ROUTEB_CY3_sdec025.json,ROUTEB_CY3_sdec0175.json \
        --k3-routeA ROUTEA_K3.json --k3-routeB ROUTEB_K3_ctrl.json --out GATE.json
    gate_compare.py --masses 1,1,1,1,1,25 --alphas 3/2 --closed none --routeA ... --routeB ...
    --closed auto (default) = the family's recorded closed forms (K3 c_{3/2}, CY3 c_{5/4} and c_{7/4},
    CY4 c_{3/2} = -sqrt(5)/(40 pi^2); a family without one -> digits vs closed omitted, cross-route
    gate + s_dec-independence only); --closed FILE.json = caller's structures {alpha: {sign, rational,
    pi, radicals, gamma}} (a form outside the ring is REFUSED); --closed none.  --bar D: the gate bar
    (default 30, the mandatory two-method floor).

Targets (GATE.json "targets", the strings pslq_calpha.py --from-gate reads): per alpha the best Route B
midpoint (largest c_acc_bits) trimmed to `digits` significant digits, where digits =
int(min(floor(acc_bits*log10 2), B s_dec-independence)) whenever the Route B samples carry two or more
distinct s_dec values, and floor(acc_bits*log10 2) alone (the ball's quotable digits) with a printed CAVEAT
when only one s_dec is present.  The ball's digit count bounds the string only when the series truncation
(NTHR, NSER at that s_dec) sits below the ball radius: at CY4 the s_dec = 1/4 midpoints at 896 and 1280 bits
carried the same 3.4e-251 offset from the closed form while the s_dec = 7/40 midpoint agreed to 257.9 d, so
the 1280-bit ball's 261 digits certified 248 = int(min(261, 248.22)).  One TARGET line per alpha is printed
after the "wrote" line (the K3-form mask drops it with GATE_VERDICT / PRODUCER).

Exit: 0 PASS (every cross-route gate >= bar, every sample vs closed >= bar, K3 control >= bar)
      1 FAIL (a gate below the bar, or no cross-route gate because a route is absent)
      2 usage
      3 REFUSED by name (a named input missing / unparseable, a --closed form outside the ring,
        an alpha absent from every input).
"""
import argparse, hashlib, json, os, re, subprocess, sys
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import calpha_rings as R   # noqa: E402

CALPHA_RINGS_SHA256 = hashlib.sha256(open(R.__file__, 'rb').read()).hexdigest()
SELF_SHA256 = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()


class Refused(Exception):
    pass


def digits(x, y):
    x = mp.mpc(x); y = mp.mpc(y)
    if x == y: return float(mp.mp.dps)
    return float(-mp.log10(abs(x - y) / max(abs(y), mp.mpf(10)**(-mp.mp.dps))))


def parse_arb_mid(s):
    """Parse Arblib midpoint string like '[-0.123456 +/- 1.2e-40]' -> mpf."""
    m = re.match(r'\[?\s*([+\-0-9.eE]+)\s*(?:\+/-|±)', s)
    if m: return mp.mpf(m.group(1))
    return mp.mpf(s)


def load_routeA(path, alpha):
    """Route A JSON: {tag: {per_sdec: {sdec: {c_alpha: {alpha: {c_re: str}}}}}} -> {tag@sdec: mpf}."""
    if not os.path.exists(path): raise Refused(f'Route A input missing: {path}')
    d = json.load(open(path))
    if 'per_sdec' in d:      # a single-run dict (routeA_monoproj.run output written bare, e.g. a probe file)
        d = {d.get('tag', os.path.basename(path)): d}
    out = {}
    for tag, run in d.items():
        for sdec, sd in run["per_sdec"].items():
            c = sd["c_alpha"].get(alpha)
            if c:
                out[f"{tag}@{sdec}"] = mp.mpf(c["c_re"])
    return out


def load_routeB(path, alpha):
    """Tolerant parser for the (slightly malformed) Julia JSON output.
    Returns {B@sdec: mpf midpoint, _acc_bits: int, _mid_str: str} for alpha ('5/4' -> key '5//4')."""
    if not os.path.exists(path): raise Refused(f'Route B input missing: {path}')
    txt = open(path).read()
    sdec_m = re.search(r'"s_dec":\s*"([^"]+)"', txt)
    sdec = sdec_m.group(1) if sdec_m else "?"
    akey = alpha.replace("/", "//")
    blk_m = re.search(r'"' + re.escape(akey) + r'":\s*\{([^}]+)\}', txt)
    if not blk_m: return {}
    blk = blk_m.group(1)
    mid_m = re.search(r'"c_re_mid":\s*"\[([+\-0-9.eE]+)\s*(?:\+/-|±)', blk)
    acc_m = re.search(r'"c_acc_bits":\s*(\d+)', blk)
    if not mid_m: return {}
    return {f"B@{sdec}": mp.mpf(mid_m.group(1)),
            "_acc_bits": int(acc_m.group(1)) if acc_m else 0, "_mid_str": mid_m.group(1)}


def quotable_digits(acc_bits):
    """Certified decimal digits of an Arb ball with acc_bits accurate bits: floor(acc_bits*log10(2))."""
    return int(acc_bits * 0.30102999566398120)


def stamp():
    return subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'], capture_output=True, text=True).stdout.strip()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--family', help=f'one of {R.families()} (or --masses + --alphas)')
    ap.add_argument('--masses', help='squared masses, e.g. 1,1,1,1,1,25')
    ap.add_argument('--alphas', help='fractional exponents, e.g. 3/2 or 5/4,7/4 (default: the family table)')
    ap.add_argument('--closed', default='auto', help='auto | none | FILE.json of closed-form structures')
    ap.add_argument('--routeA', default='', help='Route A JSON path(s), comma-separated')
    ap.add_argument('--routeB', default='', help='Route B JSON path(s), comma-separated')
    ap.add_argument('--k3-routeA', default='', help='Route A JSON of the K3 (1,1,1,9) positive control')
    ap.add_argument('--k3-routeB', default='', help='Route B JSON of the K3 positive control')
    ap.add_argument('--dps', type=int, default=300)
    ap.add_argument('--bar', type=float, default=30.0, help='two-method gate bar in digits (mandatory floor 30)')
    ap.add_argument('--out', default='GATE.json')
    a = ap.parse_args(argv)
    if not a.family and not a.masses:
        ap.error('--family or --masses required')
    mp.mp.dps = a.dps   # set BEFORE any mp.mpf(string) — import-dps footgun
    try:
        return run(a)
    except Refused as e:
        print(f'REFUSED: {e}')
        return 3
    except R.NotRingMember as e:
        print(f'REFUSED: NotRingMember: {e}')
        return 3
    except R.RingSpecError as e:
        print(f'REFUSED: RingSpecError: {e}')
        return 3


def run(a):
    if a.family:
        spec = R.family_spec(a.family)
        if a.alphas: spec = R.spec_from_masses(spec['msq'], a.alphas.split(','), spec['pi_power'], spec['closed'], a.family)
    else:
        if not a.alphas: raise Refused('--alphas required with --masses')
        spec = R.spec_from_masses(a.masses.split(','), a.alphas.split(','))
    fam = spec['family'] or 'masses'
    alphas = list(spec['alphas'])
    if not a.routeA and not a.routeB: raise Refused('no --routeA / --routeB input named')
    closed = {}
    for alpha in alphas:
        cf = R.closed_form(spec, alpha, a.closed)
        if cf is not None:
            cf.log_vector(R.ring_names(spec, 'log', 'extended'))   # ring-membership check -> NotRingMember
            closed[alpha] = cf
    A_files = [f for f in a.routeA.split(',') if f]
    B_files = [f for f in a.routeB.split(',') if f]
    for f in A_files + B_files + [a.k3_routeA, a.k3_routeB]:
        if f and not os.path.exists(f): raise Refused(f'named input missing: {f}')

    print("="*78)
    print(f"TWO-METHOD GATE: {', '.join('c_{'+al+'}' for al in alphas)} for {fam} ({','.join(str(m) for m in spec['msq'])}) banana coalescence")
    print("="*78)

    result = {"closed_form": {al: closed[al].display for al in alphas if al in closed}}
    gates, targets, fails, target_lines = [], {}, [], []
    for alpha in alphas:
        print(f"\n--- c_{{{alpha}}} ---")
        cf = closed.get(alpha)
        if cf is not None:
            cval = cf.value()
            print(f"  closed form: {cf.display}")
            print(f"  closed value: {mp.nstr(cval, 200)}")
        else:
            cval = None
            print("  closed form: NONE RECORDED (digits vs closed omitted; cross-route gate + s_dec-independence only)")
        A = {}
        for f in A_files:
            A.update(load_routeA(f, alpha))
        B, best = {}, None
        for f in B_files:
            b = load_routeB(f, alpha)
            B.update({k: v for k, v in b.items() if not k.startswith('_')})
            if b and (best is None or b['_acc_bits'] > best[0]):
                best = (b['_acc_bits'], b['_mid_str'], f)
        if not A and not B: raise Refused(f'alpha {alpha} absent from every named input')
        print(f"  Route A samples: {list(A.keys())}")
        print(f"  Route B samples: {list(B.keys())}")
        rep = {}
        if cval is not None:
            for src, vals in [("A", A), ("B", B)]:
                for k, v in vals.items():
                    d = digits(v, cval)
                    rep[f"{src}:{k}"] = round(d, 2)
                    print(f"    {src}:{k} = {mp.nstr(v,60)}  vs closed: {d:.2f} d")
                    if d < a.bar: fails.append(f'c_{{{alpha}}} {src}:{k} vs closed {d:.2f} d < bar {a.bar}')
        # cross-route gate (best A vs best B)
        if A and B:
            Avals = list(A.values()); Bvals = list(B.values())
            gate = min(digits(x, y) for x in Avals for y in Bvals)
            print(f"  CROSS-ROUTE GATE (min A vs B): {gate:.2f} d")
            rep["cross_route_gate_d"] = round(gate, 2)
            gates.append((alpha, gate))
            if gate < a.bar: fails.append(f'c_{{{alpha}}} cross-route gate {gate:.2f} d < bar {a.bar}')
        else:
            fails.append(f'c_{{{alpha}}}: one route only ({"A" if A else "B"}) — no cross-route gate')
        # s_dec independence
        for src, vals in [("A", A), ("B", B)]:
            keys = list(vals.keys())
            if len(keys) >= 2:
                ind = min(digits(vals[keys[i]], vals[keys[j]])
                          for i in range(len(keys)) for j in range(i+1, len(keys)))
                print(f"  {src} s_dec-independence: {ind:.2f} d")
                rep[f"{src}_sdec_ind_d"] = round(ind, 2)
        entry = dict(per_sample_vs_closed=rep)
        if cval is not None:
            entry = dict(closed_value=mp.nstr(cval, 250), per_sample_vs_closed=rep)
        result[f"c_{alpha}"] = entry
        if best is not None:
            nd_ball = quotable_digits(best[0])
            b_ind = rep.get("B_sdec_ind_d")     # present only when the Route B samples carry >= 2 distinct s_dec
            if b_ind is not None:
                nd = int(min(nd_ball, b_ind))
                rule = (f"int(min(ball digits floor(acc_bits*log10 2) = {nd_ball}, B s_dec-independence = {b_ind} d)) = {nd}; "
                        f"the best Route B midpoint trimmed to {nd} digits")
            else:
                nd = nd_ball
                rule = (f"ball digits floor(acc_bits*log10 2) = {nd_ball}; the best Route B midpoint trimmed to {nd} digits; "
                        f"CAVEAT single s_dec: the ball alone bounds this string; no s_dec-independence member -- run a second s_dec")
            with mp.workdps(nd + 20):
                tstr = mp.nstr(mp.mpf(best[1]), nd)
            targets[alpha] = {"value": tstr, "digits": nd, "acc_bits": best[0], "source": os.path.basename(best[2]),
                              "ball_digits": nd_ball, "B_sdec_ind_d": b_ind, "rule": rule}
            target_lines.append(f"TARGET c_{{{alpha}}}: {nd} d [{os.path.basename(best[2])}, acc_bits {best[0]}] rule: {rule}")

    # product check (two closed forms whose product is rational x pi^k)
    if len(alphas) == 2 and all(al in closed for al in alphas):
        prod = R.product_identity(closed[alphas[0]], closed[alphas[1]])
        if prod is not None:
            a0, a1 = alphas
            print(f"\n--- product c_{{{a0}}}*c_{{{a1}}} = {prod.render_rational_pi('log')} ---")
            prod_closed = prod.value()
            print(f"  closed: {mp.nstr(prod_closed, 60)}")
            print(f"  c{a0.replace('/', '')}*c{a1.replace('/', '')}: {mp.nstr(closed[a0].value()*closed[a1].value(), 60)}")
            print(f"  agree: {digits(closed[a0].value()*closed[a1].value(), prod_closed):.1f} d (algebraic identity)")
            result["product_identity"] = f"c_{{{a0}}}*c_{{{a1}}} = {prod.render_rational_pi('json')}"

    # K3 positive control
    if a.k3_routeA or a.k3_routeB:
        k3 = R.family_spec('K3')
        cf3 = R.closed_form(k3, '3/2')
        CLOSED_K3 = cf3.value()
        print(f"\n--- K3 positive control ---")
        A_K3 = load_routeA(a.k3_routeA, "3/2") if a.k3_routeA else {}
        B_K3 = load_routeB(a.k3_routeB, "3/2") if a.k3_routeB else {}
        ctrl = {}
        for src, vals in [("A", A_K3), ("B", {k: v for k, v in B_K3.items() if not k.startswith('_')})]:
            for k, v in vals.items():
                d = digits(v, CLOSED_K3)
                print(f"    {src}:{k} = {mp.nstr(v,60)}  vs {cf3.display}: {d:.2f} d")
                ctrl[f"{src}:{k}"] = round(d, 2)
                if d < a.bar: fails.append(f'K3 control {src}:{k} {d:.2f} d < bar {a.bar}')
        if not ctrl: fails.append('K3 control named but empty')
        result["K3_control"] = ctrl

    result["targets"] = targets
    result["family"] = fam; result["msq"] = list(spec['msq']); result["alphas"] = alphas
    result["gate"] = {"bar_d": a.bar, "cross_route_gate_d": {al: round(g, 2) for al, g in gates},
                      "verdict": "PASS" if not fails else "FAIL", "fails": fails}
    result["producer"] = {"script": os.path.basename(__file__), "script_sha256": SELF_SHA256,
                          "calpha_rings_sha256": CALPHA_RINGS_SHA256, "dps": a.dps, "stamp_utc": stamp(),
                          "inputs": {"routeA": A_files, "routeB": B_files, "k3_routeA": a.k3_routeA, "k3_routeB": a.k3_routeB,
                                     "closed": a.closed}}
    json.dump(result, open(a.out, "w"), indent=1)
    print(f"\nwrote {a.out}")
    for line in target_lines: print(line)
    gmin = min((g for _, g in gates), default=None)
    print(f"GATE_VERDICT family={fam} alphas={','.join(alphas)} cross_route_min_d={gmin if gmin is None else f'{gmin:.2f}'} "
          f"bar_d={a.bar} closed={'yes' if closed else 'none'} k3_control={'yes' if a.k3_routeA or a.k3_routeB else 'no'} "
          f"{'PASS' if not fails else 'FAIL'}" + ('' if not fails else ' :: ' + ' | '.join(fails)))
    print(f"PRODUCER gate_compare.py sha256={SELF_SHA256[:16]} calpha_rings.py sha256={CALPHA_RINGS_SHA256[:16]} dps={a.dps}")
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())

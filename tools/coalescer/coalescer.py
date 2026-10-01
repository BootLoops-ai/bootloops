#!/usr/bin/env python3
"""coalescer.py — spectral-projector coalescence-coefficient extractor.
Extracts c_α (coefficient of ϖ₀ onto non-unipotent threshold branch Φ_α) via
P_λ = (M₀−I)^n · ∏_{μ≠λ}(M₀−μI) / [(λ−1)^n · ∏(λ−μ)].
Generic over PF-order/mass-tuple/eigenvalue-set.
Two doors onto one kernel (routeA_monoproj.run_operator):
  --op FILE.json   GENERIC: your Fuchsian operator L = Σ_j P_j(z) θ_z^j (exact polynomial
                   coefficients) + the exact Taylor coefficients of its holomorphic solution at
                   z = 0 (the seed) [+ alpha, label, seed_radius, knobs]; schema: GUIDE.md INPUTS,
                   and dump_pf_data.py writes exactly this layout.
  --masses a,b,..  the banana worked example: the L-loop banana PF operator for squared masses
                   (a,b,..) and its BFKNS series seed are built for you.
Usage: coalescer.py --masses 1,1,1,9 --alpha 3/2 --dps 200
       coalescer.py --op my_operator.json --alpha 3/2 --dps 200
Example --op file (and template): examples/gauss_2f1.json (Gauss 2F1(1/4,3/4;1;z)).
Self-test: coalescer.py --test  (K3 c_{3/2}=−√3/(36π) probe-mode through BOTH doors: the banana
door in-process and the same operator re-exported to JSON through --op; plus the Gauss 2F1
example through --op against Gauss's connection formula; full 195d in GATE.json)
Anchors: --sb S [--sb2 S2] [--lift L] override the defaults 2.2*thr / 4.2*thr / thr/4 (--sb alone
sets sb2 = 2*sb); the anchor law is checked from the operator's own singularities BEFORE any
transport and a violating value is REFUSED by name (rc 3); --check-only prints the check and exits;
--anchor-law seed|roots selects the floor (see --help).  Exit codes: 0 ok, 1 gate fail, 2 usage,
3 refused by the anchor law."""
import sys, os, json, time, argparse, subprocess, shutil, tempfile
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from routeA_monoproj import run as runA, run_operator, load_operator_json
from pf_engine import find_minimal_op, z as _z
import math
import sympy as sp

SELFTEST_BAR_D = 20


def selftest():
    """K3(1,1,1,9) probe: c_{3/2} vs −√3/(36π) through BOTH doors, plus one operator that is
    not a banana.  Leg 1 (banana door, in-process): routeA_monoproj.run builds the operator and
    the BFKNS seed.  Leg 2 (generic door, round trip): the same operator and seed are
    re-exported to JSON by dump_pf_data.dump and fed to this CLI's --op in a child process; its
    value must match the closed form to the same bar and is compared with leg 1's.  Leg 3
    (generic door, independent closed form): the shipped examples/gauss_2f1.json (Gauss
    2F1(1/4,3/4;1;z), hand-written operator + Pochhammer seed) through --op in a second child;
    c_a must equal e^{i pi a} Gamma(b-a)/(Gamma(b)Gamma(1-a)) (Gauss's connection formula) for
    a = 1/4 and 3/4 at both matching points.  The children start first so the three legs
    overlap.  Full 195d/223d gate: GATE.json.  Scratch in a mkdtemp() dir (honours TMPDIR),
    removed on exit."""
    t0 = time.time()
    mp.mp.dps = 60
    ref = -mp.sqrt(3) / (36 * mp.pi)
    tmp = tempfile.mkdtemp(prefix="coalescer_selftest_")
    try:
        # ---- leg 3 setup: the shipped non-banana example through --op (child, started first) ----
        gauss_json = os.path.join(HERE, "examples", "gauss_2f1.json")
        cmd3 = [sys.executable, os.path.abspath(__file__), "--op", gauss_json, "--dps", "30",
                "--out", "gauss_out.json"]
        so3 = open(os.path.join(tmp, "gauss_stdout.txt"), "w")
        child3 = subprocess.Popen(cmd3, stdout=so3, stderr=subprocess.STDOUT, cwd=tmp)
        # ---- leg 2 setup: K3 operator + seed -> JSON (the --op layout), then the --op child ----
        from dump_pf_data import dump
        _, op_json, summ = dump((1, 1, 1, 9), nthr=40, nser=500, out=os.path.join(tmp, "pf_k3_op.jl"),
                                verbose=False)
        print(f"[selftest] leg 2: K3 operator re-exported to {os.path.basename(op_json)} "
              f"(order={summ['order']} degz={summ['degz']} residual={summ['residual']} "
              f"frac_exps={summ['frac_exps']} seed terms={len(summ['bfkns_a']) - 1} "
              f"seed_radius={summ['seed_radius']})  [{time.time()-t0:.1f}s]")
        op_out = os.path.join(tmp, "op_out.json")
        cmd = [sys.executable, os.path.abspath(__file__), "--op", os.path.basename(op_json), "--dps", "30",
               "--sb", "80", "--sb2", "150", "--lift", "10", "--sdec", "0.15", "--nthr", "40",
               "--out", os.path.basename(op_out)]
        with open(os.path.join(tmp, "op_stdout.txt"), "w") as so:
            child = subprocess.Popen(cmd, stdout=so, stderr=subprocess.STDOUT, cwd=tmp)
            # ---- leg 1: banana door, in-process (unchanged) ----
            out = runA("K3_selftest", (1, 1, 1, 9), dps=30, sb=80, sb2=150, lift=10,
                       sdec_list=[mp.mpf('0.15')], Nser=500, Nthr=40, frac=0.3, nloop=32)
            rc2 = child.wait()
            rc3 = child3.wait(); so3.close()
        ca = out['per_sdec']['0.15']['c_alpha']['3/2']
        c = mp.mpc(mp.mpf(ca['c_re']), mp.mpf(ca['c_im']))
        d = -mp.log10(abs((c - ref) / ref))
        print(f"\n[selftest] c_3/2 = {mp.nstr(c, 25)}")
        print(f"[selftest] ref   = {mp.nstr(ref, 25)}")
        print(f"[selftest] match = {float(d):.2f} d   wall = {time.time()-t0:.1f}s")
        assert d >= SELFTEST_BAR_D, f"K3 probe self-test FAILED (banana door): {float(d):.2f}d < {SELFTEST_BAR_D}d"
        # ---- leg 2 verdict ----
        child_out = open(os.path.join(tmp, "op_stdout.txt"), encoding="utf-8").read()
        print("\n[selftest] leg 2: --op child output follows (coalescer.py " + " ".join(cmd[2:]) + ")")
        print("\n".join("  | " + ln for ln in child_out.rstrip("\n").splitlines()))
        assert rc2 == 0, f"--op child exited rc={rc2}"
        J = json.load(open(op_out, encoding="utf-8"))
        ca2 = J['per_sdec']['0.15']['c_alpha']['3/2']
        c2 = mp.mpc(mp.mpf(ca2['c_re']), mp.mpf(ca2['c_im']))
        d2 = -mp.log10(abs((c2 - ref) / ref))
        d12 = -mp.log10(abs((c2 - c) / ref)) if c2 != c else mp.mpf(30)
        print(f"\n[selftest] --op door c_3/2 = {mp.nstr(c2, 25)}")
        print(f"[selftest] --op door match = {float(d2):.2f} d vs closed form; "
              f"{float(d12):.1f} d vs the banana door (same kernel; seeds differ only in "
              f"exact-int vs mp-float construction)   wall = {time.time()-t0:.1f}s")
        assert d2 >= SELFTEST_BAR_D, f"K3 probe self-test FAILED (--op door): {float(d2):.2f}d < {SELFTEST_BAR_D}d"
        # ---- leg 3 verdict: Gauss 2F1(1/4,3/4;1;z), c_a = e^{i pi a} Gamma(b-a)/(Gamma(b)Gamma(1-a)) ----
        child3_out = open(os.path.join(tmp, "gauss_stdout.txt"), encoding="utf-8").read()
        print("\n[selftest] leg 3: --op child output follows (coalescer.py --op examples/gauss_2f1.json "
              + " ".join(cmd3[4:]) + ")")
        print("\n".join("  | " + ln for ln in child3_out.rstrip("\n").splitlines()))
        assert rc3 == 0, f"--op child (gauss) exited rc={rc3}"
        G = json.load(open(os.path.join(tmp, "gauss_out.json"), encoding="utf-8"))
        qa, qb = mp.mpf(1) / 4, mp.mpf(3) / 4
        refs = {"1/4": mp.expj(mp.pi * qa) * mp.gamma(qb - qa) / (mp.gamma(qb) * mp.gamma(1 - qa)),
                "3/4": mp.expj(mp.pi * qb) * mp.gamma(qa - qb) / (mp.gamma(qa) * mp.gamma(1 - qb))}
        d3 = []
        for al, rf in refs.items():
            for sd, rec in G['per_sdec'].items():
                c3 = mp.mpc(mp.mpf(rec['c_alpha'][al]['c_re']), mp.mpf(rec['c_alpha'][al]['c_im']))
                d3.append(float(-mp.log10(abs((c3 - rf) / rf))))
            print(f"[selftest] gauss c_{{{al}}}: ref = {mp.nstr(rf, 25)}   match = "
                  f"{min(d3[-len(G['per_sdec']):]):.2f} d over s_dec in {list(G['per_sdec'])}")
        print(f"[selftest] --op door (Gauss 2F1, not a banana) match = {min(d3):.2f} d vs Gauss's "
              f"connection formula   wall = {time.time()-t0:.1f}s")
        assert min(d3) >= SELFTEST_BAR_D, f"self-test FAILED (--op door, Gauss 2F1): {min(d3):.2f}d < {SELFTEST_BAR_D}d"
        print("[selftest] PASS (probe mode, both doors, K3 round trip + Gauss 2F1; full 195d Route-B gate in GATE.json)")
        return float(d), float(d2), min(d3)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


EXIT_OK, EXIT_GATE_FAIL, EXIT_USAGE, EXIT_REFUSED = 0, 1, 2, 3

ANCHOR_LAW_HELP = """\
Anchor law (README): the anchor s_b must satisfy s_b > 1/min|z_sing| (the MUM
convergence radius of the seed series), NOT 1.5*max|z_sing|.  Units: z = -1/s
(pf_engine); the BFKNS seed is sum_n a_n z_b^n at z_b = -1/s_b, so the disk
|z_b| < R reads s_b > 1/R.  Two readings of min|z_sing|, both computed from the
operator this run builds and both printed on the [anchor-law] line:
  seed   (default) R = 1/thr with thr = (sum_i sqrt(m_i^2))^2: the leading
         threshold z = -1/thr is the nearest true singularity of the seed
         series (exact for the multinomial series: thr^n/(n+1)^(L-1) <= |a_n|
         <= thr^n), so the floor is s_b > thr.  Leading-polynomial roots
         inside |z| < 1/thr are apparent for the seed.
  roots  R = min|z| over ALL roots of the operator's leading polynomial (the
         list Route B's auto-anchor ceil(1.6/min|z_sing|) uses); floor
         s_b > 1/min|z_sing|.  Conservative: it counts apparent singularities.
Example, --masses 1,1,1,1,1,25 (thr = 100): the leading polynomial has
14 roots with min|z_sing| = 0.0013249, so --anchor-law roots demands
s_b > 754.80 and REFUSES the defaults 220/420 (rc 3); the seed law
accepts them (floor 100; the seed agrees to 40 of 40 digits at s_b = 220
between 600 and 400 terms).  The recorded probe of that tuple ran
--sb 1210 --sb2 2420; Route B's auto-anchor for it is 1208.
Both sb and sb2 are checked; a value at or below the floor in force is refused
by name before anything is transported.
Exit codes: 0 ok; 1 gate fail (--test below its bar, or an error); 2 usage;
3 REFUSED by the anchor law.
"""


def leading_zsings(Pj, order):
    """Roots of the leading polynomial P_order(z): the same computation as pf_engine.build_ode_s
    (sp.nroots at 50 digits, roots with |z| > 1e-40 kept), sorted by |z|."""
    roots = []
    for r in sp.nroots(sp.Poly(Pj[order], _z), n=50):
        rr = complex(r)
        if abs(rr) > 1e-40:
            roots.append(rr)
    return sorted(roots, key=abs)


def anchor_law(msq, Pj, order, thr, sb, sb2, law, what=None, seed_radius=None):
    """Both floors from the operator's own singularities, and the verdict on (sb, sb2) under `law`.
    Banana door: thr = (sum sqrt m)^2 and the seed disk is |z| < 1/thr (z = -1/thr a leading root).
    Generic door (--op): msq is None; seed_radius R (the JSON's, optional) gives the seed floor
    s_b > 1/R, honoured only if some leading root has |z| = R; without it only the roots law exists."""
    roots = leading_zsings(Pj, order)
    min_abs = min(abs(r) for r in roots)
    note = None
    law_in_force = law
    if msq is not None:
        thr_f = float(thr)
        z_thr = -1.0 / thr_f
        thr_is_root = any(abs(r - z_thr) <= 1e-9 * abs(z_thr) for r in roots)
        R = 1.0 / thr_f
        if law == "seed" and not thr_is_root:
            law_in_force = "roots"
            note = "z=-1/thr is not a root of the leading polynomial: the roots law is in force"
    elif seed_radius is not None:
        R = float(seed_radius)
        thr_f = 1.0 / R
        thr_is_root = any(abs(abs(r) - R) <= 1e-9 * R for r in roots)
        if law == "seed" and not thr_is_root:
            law_in_force = "roots"
            note = "no leading root has |z| = seed_radius: the roots law is in force"
    else:
        R = None; thr_f = None; thr_is_root = False
        if law == "seed":
            law_in_force = "roots"
            note = "the operator JSON declares no seed_radius: the roots law is in force"
    inside = 0 if R is None else sum(1 for r in roots if abs(r) < R * (1 - 1e-9))
    floor_roots = 1.0 / min_abs
    floor_seed = thr_f
    floor = floor_seed if law_in_force == "seed" else floor_roots
    if what is None:
        what = f"masses={_fmt_masses(msq)}" if msq is not None else "op"
    return dict(msq=(tuple(msq) if msq is not None else None), what=what, thr=thr_f, seed_radius=R,
                order=order, n_roots=len(roots), min_abs=min_abs,
                floor_roots=floor_roots, floor_seed=floor_seed, thr_is_root=thr_is_root,
                roots_inside_seed_disk=inside, law=law, law_in_force=law_in_force, floor=floor,
                sb=sb, sb2=sb2, sb_ok=sb > floor, sb2_ok=sb2 > floor,
                sb_auto_routeB=int(math.ceil(1.6 / min_abs)), note=note)


def _fmt_masses(msq):
    return "(" + ",".join(str(m) for m in msq) + ")"


def anchor_law_line(r, src, route):
    if r['msq'] is not None:
        head = (f"[anchor-law] {r['what']} thr={r['thr']!r} order={r['order']} "
                f"leading_roots={r['n_roots']} min|z_sing|={r['min_abs']!r} floor_roots={r['floor_roots']!r} "
                f"floor_seed={r['floor_seed']!r} (z=-1/thr is a leading root: {'yes' if r['thr_is_root'] else 'no'}; "
                f"leading roots inside the seed disk |z|<1/thr: {r['roots_inside_seed_disk']}, apparent for the seed) ")
    else:
        head = (f"[anchor-law] {r['what']} seed_radius={r['seed_radius']!r} order={r['order']} "
                f"leading_roots={r['n_roots']} min|z_sing|={r['min_abs']!r} floor_roots={r['floor_roots']!r} "
                f"floor_seed={r['floor_seed']!r} (1/seed_radius; a leading root on |z|=seed_radius: "
                f"{'yes' if r['thr_is_root'] else 'no'}; leading roots inside the seed disk: "
                f"{r['roots_inside_seed_disk']}, apparent for the seed) ")
    s = (head +
         f"law={r['law_in_force']} floor={r['floor']!r} sb={r['sb']!r} ({src['sb']}) sb2={r['sb2']!r} ({src['sb2']}) "
         f"lift={src['lift_val']!r} ({src['lift']}) "
         f"verdict={'OK' if (r['sb_ok'] and r['sb2_ok']) else 'REFUSED'}")
    if route == "B":
        s += (f" route=B sb_auto(routeB)={r['sb_auto_routeB']} (ceil(1.6/min|z_sing|); routeB_monoproj.jl recomputes it "
              f"from its own singularity list at run time; env SB overrides it)")
    if r['note']:
        s += f" note: {r['note']}"
    return s


def anchor_law_refusal(r, src):
    bad = []
    if not r['sb_ok']:
        bad.append(f"sb = {r['sb']!r} ({src['sb']})")
    if not r['sb2_ok']:
        bad.append(f"sb2 = {r['sb2']!r} ({src['sb2']})")
    if r['law_in_force'] == "seed":
        rule = (f"law=seed floor s_b > thr = {r['floor']!r} (the BFKNS seed series converges for |z| < 1/thr; "
                f"z_b = -1/s_b)") if r['msq'] is not None else \
               (f"law=seed floor s_b > 1/seed_radius = {r['floor']!r} (the seed series converges for "
                f"|z| < seed_radius = {r['seed_radius']!r}; z_b = -1/s_b)")
        alt = ""
    else:
        rule = (f"law=roots floor s_b > 1/min|z_sing| = {r['floor']!r} (min|z_sing| = {r['min_abs']!r} over the "
                f"{r['n_roots']} roots of the order-{r['order']} operator's leading polynomial)")
        alt = ((f", or --anchor-law seed (floor thr = {r['floor_seed']!r})" if r['msq'] is not None else
                f", or --anchor-law seed (floor 1/seed_radius = {r['floor_seed']!r})")
               if (r['law'] != "seed" and r['floor_seed'] is not None) else "")
    return (f"[anchor-law] REFUSED (rc {EXIT_REFUSED}): {r['what']}: {rule}; "
            + " and ".join(bad) + (" is" if len(bad) == 1 else " are")
            + f" not above the floor. Nothing transported. Override: --sb S --sb2 S2 with S, S2 > {r['floor']!r}{alt}.")


OP_SCHEMA_HELP = """\
--op FILE.json (the generic door; dump_pf_data.py writes this layout, key aliases in brackets):
  "theta_form_Pj" ["theta_coeffs", "Pj"]: {"0": [c00, c01, ...], ..., "r": [...]} or a list of
      r+1 lists -- P_j(z) exact coefficients (ints, "p/q" strings or [p, q] pairs), LOW -> HIGH
      power of z, for the operator L = sum_j P_j(z) theta_z^j, theta_z = z d/dz.  z = 0 must be
      the point where your seed converges; the branches are read at z = infinity (t = 1/z).
  "bfkns_a" ["seed"]: [a0, a1, ...] -- exact Taylor coefficients of the holomorphic solution
      f(z) = sum_n a_n z^n at z = 0 (ints / "p/q" / [p, q]); all are used unless --nser N.
  "order": r (optional, checked), "alpha": "p/q" (optional: the exponent to headline; every
      fractional exponent at z = infinity is extracted), "label": "...",
  "seed_radius": R (optional: the seed's convergence radius in z; enables --anchor-law seed
      with floor s_b > 1/R), "sb", "sb2", "lift", "sdec": [...], "nthr", "nser" (optional run
      knobs; command-line flags override them).  Defaults without them: sb = 2.2*floor,
      sb2 = 4.2*floor, lift = floor/4 with floor = the anchor-law floor in force.
"""


def _parse_sdec(s, dps):
    """Comma-separated string or a JSON list -> mp numbers made at the run's precision."""
    items = [x.strip() for x in s.split(",")] if isinstance(s, str) else [str(x) for x in s]
    with mp.workdps(max(int(dps), 15)):
        return [mp.mpf(x) for x in items if x]


def main():
    ap = argparse.ArgumentParser(description=__doc__, epilog=OP_SCHEMA_HELP + "\n" + ANCHOR_LAW_HELP,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    door = ap.add_mutually_exclusive_group()
    door.add_argument("--masses", help="banana door: comma-separated m² tuple, e.g. 1,1,1,9")
    door.add_argument("--op", metavar="FILE.json",
                      help="generic door: operator + seed JSON (schema below; dump_pf_data.py writes it)")
    ap.add_argument("--alpha", help="fractional exponent p/q (optional; all extracted if omitted)")
    ap.add_argument("--dps", type=int, default=80)
    ap.add_argument("--route", choices=["A", "B"], default="A")
    ap.add_argument("--out", default=None)
    ap.add_argument("--test", action="store_true", help="run the K3 c_{3/2} probe self-test (both doors)")
    ap.add_argument("--sb", type=float, default=None,
                    help="anchor s_b (default 2.2*thr, thr=(sum sqrt m)^2 [--op: 2.2*floor]); given alone, sb2 = 2*sb")
    ap.add_argument("--sb2", type=float, default=None, help="second anchor s_b2 (default 4.2*thr [4.2*floor], or 2*sb with --sb)")
    ap.add_argument("--lift", type=float, default=None, help="imaginary lift of the transport path (default thr/4 [floor/4])")
    ap.add_argument("--sdec", default=None, help="comma-separated Euclidean matching points s_dec (default 0.15,0.25)")
    ap.add_argument("--nthr", type=int, default=None, help="threshold Frobenius terms (default max(60, 2*dps))")
    ap.add_argument("--nser", type=int, default=None,
                    help="seed series terms (default max(600, 15*dps) [--op: every coefficient in the file])")
    ap.add_argument("--anchor-law", choices=["seed", "roots"], default=None,
                    help="floor in force: seed = thr (exact for the BFKNS seed) [--op: 1/seed_radius]; roots = "
                         "1/min|z_sing| over all leading-polynomial roots (conservative) "
                         "(default: seed for --masses, roots for --op)")
    ap.add_argument("--check-only", action="store_true",
                    help="build/load the operator, print the [anchor-law] line and exit 0 (ok) / 3 (refused); no transport")
    args = ap.parse_args()
    if args.test:
        selftest(); return
    if not args.masses and not args.op:
        ap.error("--masses required (or --op FILE.json, or --test)")
    src = {}
    t0 = time.time()
    if args.op:
        # ---------------- generic door ----------------
        if args.route == "B":
            ap.error("--op runs Route A (pure python); Route B reads a pf_*_data.jl (dump_pf_data.py) via env PFDATA")
        try:
            OP = load_operator_json(args.op)
        except (OSError, ValueError, TypeError, KeyError) as e:
            print(f"[op] cannot load {args.op}: {e}", file=sys.stderr)
            sys.exit(EXIT_USAGE)
        knobs = OP['knobs']
        Pj, order, degz = OP['Pj'], OP['order'], OP['degz']
        what = f"op={os.path.basename(args.op)} ({OP['label']})"
        law = args.anchor_law or "roots"
        # floor first (defaults hang off it), then the anchors
        r0 = anchor_law(None, Pj, order, None, 1.0, 1.0, law, what=what, seed_radius=OP['seed_radius'])
        floor = r0['floor']
        def pick(name, dflt, dsrc):
            v = getattr(args, name)
            if v is not None:
                return float(v), "override"
            if name in knobs:
                return float(knobs[name]), "from the JSON"
            return float(dflt), dsrc
        sb, src['sb'] = pick("sb", 2.2 * floor, "default 2.2*floor")
        if args.sb2 is None and args.sb is not None:
            sb2, src['sb2'] = 2.0 * sb, "derived 2*sb"
        else:
            sb2, src['sb2'] = pick("sb2", 4.2 * floor, "default 4.2*floor")
        lift, src['lift'] = pick("lift", floor / 4, "default floor/4")
        src['lift_val'] = lift
        r = anchor_law(None, Pj, order, None, sb, sb2, law, what=what, seed_radius=OP['seed_radius'])
        print(anchor_law_line(r, src, "A") + f"  [operator order={order} degz={degz} from {os.path.basename(args.op)}, "
              f"seed terms={len(OP['seed']) - 1}, loaded in {time.time()-t0:.1f}s]", flush=True)
        if not (r['sb_ok'] and r['sb2_ok']):
            print(anchor_law_refusal(r, src), flush=True)
            sys.exit(EXIT_REFUSED)
        if args.check_only:
            print("[anchor-law] check only: nothing transported", flush=True)
            sys.exit(EXIT_OK)
        sdec_list = _parse_sdec(args.sdec if args.sdec else knobs.get('sdec', "0.15,0.25"), args.dps)
        Nthr = args.nthr if args.nthr is not None else int(knobs.get('nthr', max(60, 2 * args.dps)))
        Nser = args.nser if args.nser is not None else (int(knobs['nser']) if 'nser' in knobs else None)
        out = run_operator(f"coalescer_op_{os.path.splitext(os.path.basename(args.op))[0]}", Pj, order,
                           OP['seed'], dps=args.dps, sb=sb, sb2=sb2, lift=lift, sdec_list=sdec_list,
                           Nthr=Nthr, frac=0.3, nloop=32, Nser=Nser, label=OP['label'])
        out['op_file'] = os.path.abspath(args.op)
        alpha = args.alpha or OP['alpha']
        if alpha:
            alpha = str(sp.Rational(alpha))
            print(f"\nc_{{{alpha}}} = {out['c_alpha_best'].get(alpha)}")
        if args.out:
            json.dump(out, open(args.out, "w"), indent=1)
            print(f"WROTE {args.out}")
        return
    # ---------------- banana door (the worked example) ----------------
    law = args.anchor_law or "seed"
    msq = tuple(int(x) for x in args.masses.split(","))
    thr = sum(mp.sqrt(m) for m in msq) ** 2
    # anchors: the written defaults, or the overrides (--sb alone derives sb2 = 2*sb)
    if args.sb is None:
        sb = float(2.2 * thr); src['sb'] = "default 2.2*thr"
    else:
        sb = float(args.sb); src['sb'] = "override"
    if args.sb2 is None:
        if args.sb is None:
            sb2 = float(4.2 * thr); src['sb2'] = "default 4.2*thr"
        else:
            sb2 = 2.0 * sb; src['sb2'] = "derived 2*sb"
    else:
        sb2 = float(args.sb2); src['sb2'] = "override"
    if args.lift is None:
        lift = float(thr / 4); src['lift'] = "default thr/4"
    else:
        lift = float(args.lift); src['lift'] = "override"
    src['lift_val'] = lift
    if args.route == "B" and args.sb is not None and not float(args.sb).is_integer():
        ap.error(f"--route B takes an integer --sb (routeB_monoproj.jl parses SB as Int): got {args.sb!r}")
    # the anchor-law check, from the operator this run would use, before anything is transported
    found = find_minimal_op(list(msq), ord_max=10, degz_max=18)
    if found is None:
        print(f"[anchor-law] no PF operator found for masses={_fmt_masses(msq)} (ord_max=10, degz_max=18)")
        sys.exit(EXIT_GATE_FAIL)
    order, degz, Pj, resid, _ = found
    r = anchor_law(msq, Pj, order, thr, sb, sb2, law)
    print(anchor_law_line(r, src, args.route) + f"  [operator order={order} degz={degz} residual={resid}, built in {time.time()-t0:.1f}s]", flush=True)
    if not (r['sb_ok'] and r['sb2_ok']):
        print(anchor_law_refusal(r, src), flush=True)
        sys.exit(EXIT_REFUSED)
    if args.check_only:
        print("[anchor-law] check only: nothing transported", flush=True)
        sys.exit(EXIT_OK)
    if args.route == "B":
        # Route B (Arb/z-coord): needs pf_*_data.jl (regenerate via dump_pf_data.py if absent);
        # PFDATA (env) selects the data file, which must be this tuple's; --sb reaches it as env SB;
        # --sb2/--lift are not used by Route B (single anchor; its lift is fixed in routeB_monoproj.jl)
        env = dict(os.environ)
        if args.sb is not None:
            env["SB"] = str(int(args.sb))
        print(f"[anchor-law] route B: PFDATA={env.get('PFDATA', 'pf_cy3_data.jl')} SB={env.get('SB', '(auto)')}", flush=True)
        cmd = ["julia", os.path.join(HERE, "routeB_monoproj.jl"), "--dps", str(args.dps)]
        sys.exit(subprocess.call(cmd, cwd=HERE, env=env))
    sdec_list = _parse_sdec(args.sdec or "0.15,0.25", args.dps)
    out = runA(f"coalescer_{args.masses}", msq, dps=args.dps,
               sb=sb, sb2=sb2, lift=lift,
               sdec_list=sdec_list,
               Nser=(args.nser if args.nser is not None else max(600, 15 * args.dps)),
               Nthr=(args.nthr if args.nthr is not None else max(60, 2 * args.dps)),
               frac=0.3, nloop=32)
    if args.alpha:
        print(f"\nc_{{{args.alpha}}} = {out['c_alpha_best'].get(args.alpha)}")
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
        print(f"WROTE {args.out}")


if __name__ == "__main__":
    main()

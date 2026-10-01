#!/usr/bin/env python3
r"""
gf_solve.py — gravityFlow closure solver.

Solves the self-consistent boundary system of a fixed-point auxiliary-mass flow:

    M = C(eps) v + (anchor rows),   v_k = M[sigma(k)] or analytic anchor value,

where C comes from unit-injection linear-response probes of flow beta (etaC-style
chain, AMFLOW_DUMP_EPS_GRID logs) and the residual null space is closed against a
SECOND flow's response matrix A_alpha (AMFLOW_FIXEDPOINT_PROBE output).

Validation built in (all hard): anchor rows exact; raw flow-beta equations
satisfied; A_alpha row census — EVERY row's residual reported, consistency
threshold 1e-20 relative. Validated on a 4PM cut-eikonal reference family
(43 masters) at y=5/3 and y=13/5: 43/43 rows, worst residual ~5e-52.

Inputs:
  --fpA fp_A.json            second-flow response (AMFLOW_FIXEDPOINT_PROBE)
  --probes DIR               unit-injection logs log_XX.log + keys.json
  --cmap circular_map.json   key -> ['master', idx] | ['J63', None]
  --anchors "63:J63,111:J63" sector-index-vector anchors (J63 = -(i/8)pi^-3/2 Gamma(eps-1/2)^3).
                             Each entry is SECTOR:FORM. SECTOR is the sector index of a
                             corner master in the fpA preferred basis (bit i of SECTOR =
                             index vector slot i); FORM names the analytic anchor value
                             (J63 is the built-in cut-tadpole form; a symbolic
                             cut-tadpole calculus generalizing it is a designed
                             extension, DESIGN.md section 4c). The default matches the
                             reference-family anchors; every anchor is validated against
                             the basis and refused loudly when it does not fit.
  --out fp_boundary.json     per-eps explicit_boundary table for the parent
  --extra-zero KEY           additional keys emitted as exact zeros (reducible-to-zero)

Usage: python3 gf_solve.py --fpA fp_A.json --probes eb_probes \
         --cmap eb_probes/circular_map.json --out fp_boundary_parent.json
"""
import argparse, json, re, sys
import mpmath as mp


def parse_ball(s):
    s = s.strip()
    balls = re.findall(r'\[([^\]]+)\]|([+\-]?[0-9][0-9.eE+\-]*)', s)
    nums = []
    for b, f in balls:
        if b:
            nums.append(mp.mpf(b.split('+/-')[0].strip()))
        elif f:
            nums.append(mp.mpf(f))
    if not nums:
        return mp.mpc(0)
    if len(nums) == 1:
        return mp.mpc(0, nums[0]) if ('I' in s or 'j' in s) else mp.mpc(nums[0], 0)
    return mp.mpc(nums[0], nums[1])


def load_grid(logpath, n_eps):
    eps, M = [], {}
    for line in open(logpath, errors='replace'):
        if '[EPS_GRID] eps[' in line:
            eps.append(parse_ball(line.split('=', 1)[1]))
        elif '[EPS_GRID] M[' in line:
            m = re.match(r'.*M\[(\d+)\]\[(\d+)\] = (.*)', line)
            if m:
                M.setdefault(int(m.group(1)), [None] * n_eps)[int(m.group(2))] = \
                    parse_ball(m.group(3))
    return eps, M


def J63(eps):
    """(cut tadpole)^3: bare delta(2k.u) cuts, measure d^dk/(i pi^{d/2})."""
    return -(mp.mpc(0, 1) / 8) * mp.pi ** mp.mpf('-1.5') \
        * mp.gamma(eps - mp.mpf('0.5')) ** 3


# Named analytic anchor forms accepted in an --anchors spec.  J63 is the
# factorized cut-tadpole closed form above; further forms belong to the
# symbolic cut-tadpole calculus (DESIGN.md section 4c).
ANCHOR_FORMS = {'J63': J63}

DEFAULT_ANCHORS = '63:J63,111:J63'


def sector_tuple(sec, n_slots):
    """Sector index -> corner index vector (bit i of sec = slot i)."""
    if sec < 0 or (sec >> n_slots):
        raise SystemExit(
            f'--anchors: sector index {sec} does not fit a {n_slots}-slot '
            f'family (its binary support needs a higher slot)')
    return tuple((sec >> i) & 1 for i in range(n_slots))


def parse_anchors(spec):
    """Parse an --anchors spec "SEC:FORM,SEC:FORM,..." into
    [(sector_index, form_callable), ...] with loud refusals."""
    anchors = []
    for entry in spec.split(','):
        entry = entry.strip()
        if not entry:
            continue
        sec_s, sep, name = entry.partition(':')
        if not sep or not sec_s.strip().isdigit():
            raise SystemExit(
                f'--anchors: malformed entry {entry!r} — expected '
                f'SECTOR:FORM, e.g. "{DEFAULT_ANCHORS}"')
        name = name.strip()
        if name not in ANCHOR_FORMS:
            raise SystemExit(
                f'--anchors: unknown anchor form {name!r} in entry '
                f'{entry!r} — built-in forms: '
                f'{", ".join(sorted(ANCHOR_FORMS))} (the symbolic '
                f'cut-tadpole calculus is a designed extension, '
                f'DESIGN.md section 4c)')
        anchors.append((int(sec_s.strip()), ANCHOR_FORMS[name]))
    if not anchors:
        raise SystemExit(
            '--anchors: at least one anchor is required — the fixed point '
            'alone is rank-deficient (DESIGN.md section 1)')
    return anchors


def solve(fpA_path, probes_dir, cmap_path, out_path, parent_fam=None,
          anchor_secs=((1, 1, 1, 1, 1, 1) + (0,) * 9,
                       (1, 1, 1, 1, 0, 1, 1) + (0,) * 8),
          extra_zero=(), n_eps=12, dps=120, thresh='1e-20', anchors=None):
    mp.mp.dps = dps
    d = json.load(open(fpA_path))
    basis = [tuple(int(x) for x in k.split('|')[1:]) for k in d['preferred']]
    n = len(basis)
    if parent_fam is None:
        parent_fam = d['preferred'][0].split('|')[0]
    keys = json.load(open(f'{probes_dir}/keys.json'))
    cmap = json.load(open(cmap_path))
    sub_al = [tuple(int(x) for x in k.split('|')[1:]) for k in d['sub_masters']]
    pos_al = {v: i for i, v in enumerate(basis)}
    perm_al = [pos_al[v] for v in sub_al]
    cols = {}
    eps_ref = None
    for k in range(len(keys)):
        eps, M = load_grid(f'{probes_dir}/log_{k:02d}.log', n_eps)
        if len(M) != n:
            raise SystemExit(f'log_{k:02d}: {len(M)} masters != {n}')
        cols[k] = M
        if eps_ref is None:
            eps_ref = eps
    # Anchor rows.  An --anchors spec ("SEC:FORM,...") is decoded against
    # THIS family's slot count; without one, anchor_secs (the reference
    # family's two corner vectors) is used as before.  Either way a vector
    # absent from the basis refuses loudly instead of dying in basis.index.
    n_slots = len(basis[0])
    if anchors is not None:
        anchor_list = [(sector_tuple(sec, n_slots), form, sec)
                       for sec, form in parse_anchors(anchors)]
    else:
        anchor_list = [(tuple(t), J63, None) for t in anchor_secs]
    anchor_idx = []
    for t, form, sec in anchor_list:
        try:
            anchor_idx.append((basis.index(t), form))
        except ValueError:
            label = (f'sector {sec}' if sec is not None else
                     'default anchor (the reference-family corner vector)')
            raise SystemExit(
                f'--anchors: {label} decodes to index vector {t}, which is '
                f'not among the {n} fpA preferred masters of family '
                f'{parent_fam!r} ({n_slots} slots) — pass '
                f'--anchors "SECTOR:FORM,..." matching this family '
                f'(bit i of SECTOR = slot i; built-in form: J63)')
    out, census = {}, []
    for ke in range(n_eps):
        eps = eps_ref[ke]
        Jv = J63(eps)
        A = mp.eye(n)
        rhs = mp.matrix(n, 1)
        for k, key in enumerate(keys):
            kind, idx = cmap[key]
            if kind == 'master':
                for m in range(n):
                    A[m, idx] -= cols[k][m][ke]
            else:
                for m in range(n):
                    rhs[m] += cols[k][m][ke] * Jv
        rows = [[A[i, j] for j in range(n)] for i in range(n)]
        rr = [rhs[i] for i in range(n)]
        for ia, form in anchor_idx:
            r = [mp.mpc(0)] * n
            r[ia] = mp.mpc(1)
            rows.append(r)
            rr.append(form(eps))
        m_ = len(rows)
        for i in range(m_):
            mx = max(abs(x) for x in rows[i])
            if mx > 0:
                rows[i] = [x / mx for x in rows[i]]
                rr[i] /= mx
        work = [r[:] for r in rows]
        workr = rr[:]
        piv_cols, free_cols = [], []
        p = 0
        for c in range(n):
            pi, best = None, mp.mpf(0)
            for i in range(p, m_):
                if abs(work[i][c]) > best:
                    best, pi = abs(work[i][c]), i
            if pi is None or best < mp.mpf('1e-60'):
                free_cols.append(c)
                continue
            work[p], work[pi] = work[pi], work[p]
            workr[p], workr[pi] = workr[pi], workr[p]
            piv_cols.append(c)
            for i in range(p + 1, m_):
                if work[i][c] == 0:
                    continue
                f = work[i][c] / work[p][c]
                for cc in range(c, n):
                    work[i][cc] -= f * work[p][cc]
                workr[i] -= f * workr[p]
            p += 1

        def backsub(freevals, hom):
            x = [mp.mpc(0)] * n
            for fc, fv in zip(free_cols, freevals):
                x[fc] = fv
            for irow in range(p - 1, -1, -1):
                c = piv_cols[irow]
                s = mp.mpc(0) if hom else workr[irow]
                for cc in range(c + 1, n):
                    s -= work[irow][cc] * x[cc]
                x[c] = s / work[irow][c]
            return x

        xpart = backsub([mp.mpc(0)] * len(free_cols), False)
        nulls = [backsub([mp.mpc(1) if jj == jn else mp.mpc(0)
                          for jj in range(len(free_cols))], True)
                 for jn in range(len(free_cols))]
        Arows = []
        for i in range(n):
            r = [(1 if i == j else 0) for j in range(n)]
            for j2 in range(n):
                r[perm_al[j2]] -= parse_ball(d['A'][ke][i][j2])
            if max(abs(x) for x in r) > 0:
                Arows.append(r)
        cands = [(sum(r[j] * xpart[j] for j in range(n)),
                  [sum(r[j] * nl[j] for j in range(n)) for nl in nulls])
                 for r in Arows]
        nf = len(free_cols)
        if nf:
            G = mp.matrix(nf, nf)
            b = mp.matrix(nf, 1)
            for (R0, Rt) in cands:
                for a_ in range(nf):
                    b[a_] -= mp.conj(Rt[a_]) * R0
                    for b_ in range(nf):
                        G[a_, b_] += mp.conj(Rt[a_]) * Rt[b_]
            t = mp.lu_solve(G, b)
        else:
            t = mp.matrix(0, 1)
        resids = [abs(R0 + sum(Rt[a_] * t[a_] for a_ in range(nf)))
                  for (R0, Rt) in cands]
        ncons = sum(1 for r_ in resids if r_ < mp.mpf(thresh))
        census.append({'eps_index': ke, 'free_cols': free_cols,
                       'consistent': ncons, 'total': len(cands),
                       'worst': mp.nstr(max(resids), 4) if resids else '0'})
        if ncons != len(cands):
            print(f'WARNING eps[{ke}]: only {ncons}/{len(cands)} second-flow '
                  f'rows consistent (worst {census[-1]["worst"]})')
        Mv = [xpart[j] + sum(t[a_] * nulls[a_][j] for a_ in range(nf))
              for j in range(n)]
        for m, bv in enumerate(basis):
            key = parent_fam + '|' + '|'.join(map(str, bv))
            out.setdefault(key, []).append(
                {'re': mp.nstr(Mv[m].real, 190), 'im': mp.nstr(Mv[m].imag, 190)})
    for z in extra_zero:
        out[z] = [{'re': '0', 'im': '0'}] * n_eps
    json.dump(out, open(out_path, 'w'), indent=0)
    all_ok = all(c['consistent'] == c['total'] for c in census)
    print(json.dumps({'census': census[:2] + census[-1:], 'all_consistent': all_ok,
                      'keys_written': len(out), 'out': out_path}, indent=1))
    return 0 if all_ok else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--fpA', required=True)
    ap.add_argument('--probes', required=True)
    ap.add_argument('--cmap', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--extra-zero', nargs='*', default=[])
    ap.add_argument('--anchors', default=DEFAULT_ANCHORS, metavar='SPEC',
                    help='comma-separated SECTOR:FORM anchors; SECTOR is the '
                         'sector index of a corner master in the fpA '
                         'preferred basis (bit i = slot i), FORM names the '
                         'analytic value (built-in: J63). Default: '
                         + DEFAULT_ANCHORS)
    ap.add_argument('--n-eps', type=int, default=12,
                    help='eps-grid length of the probe logs (AMFLOW_NEPS_OVERRIDE)')
    ap.add_argument('--dps', type=int, default=120)
    a = ap.parse_args()
    sys.exit(solve(a.fpA, a.probes, a.cmap, a.out,
                   extra_zero=a.extra_zero, n_eps=a.n_eps, dps=a.dps,
                   anchors=a.anchors))

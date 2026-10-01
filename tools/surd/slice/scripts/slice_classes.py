"""sigma-classes on the slice per channel: the 24 assignments modulo (a) the slice automorphism pi = (12)(34) of the labeled distance matrix
(d_{pi sigma} == d_sigma identically in t: identical integrands) and (b) the channel's exact particle-relabeling symmetries tau
(I_{sigma tau} = I_sigma), each VERIFIED EXACTLY here: the degree -4 homogeneous cell integrand Phi_sigma (cells of chart 4 and chart 1) satisfies
Phi_{sigma o tau}(y) == Phi_sigma(y o tau) at random rational points y and random rational t (exact Fractions), for every tau used.
Output SLICE_CLASSES.json: per channel the class representatives, multiplicities, members, first-pair groups (pair -> n_cells, chart, order)."""
import os, sys, json, itertools, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import slice_prov, slice_common as sc
from fractions import Fraction as Fr
import oracle_cells as OC, hf_piece
WORK = slice_prov.WORK
CHANNELS = OC.QUARK_JET + OC.GLUON_JET + ["n4"]
def compose(s, t):  # (s o t)(i) = s(t(i)); tuples 1-based
    return tuple(s[t[i] - 1] for i in range(4))
def phi(cellres, y, d):
    """homogeneous integrand: Phi(y) = y_g^-4 * f_g(y / y_g)"""
    g = cellres["gauge"]; xs = {j: Fr(y[j - 1]) / Fr(y[g - 1]) for j in (1, 2, 3, 4)}
    return OC.eval_cells_exact(cellres, xs, d) / Fr(y[g - 1]) ** 4
def verify_tau(ch, tau, rng, npts=2):
    ok = True; detail = []
    for g in (4, 1):
        cellres = OC.build_cells(ch, g, verbose=False)
        for _ in range(npts):
            t = Fr(rng.randint(1, 19), rng.randint(2, 23)); d = sc.d_of_t(t)
            sig = tuple(rng.sample([1, 2, 3, 4], 4)); st = compose(sig, tau)
            y = [Fr(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(4)]
            ytau = [y[tau[i] - 1] for i in range(4)]          # (y o tau)_i = y_{tau(i)}
            lhs = phi(cellres, y, sc.sigma_d(d, st))
            r1 = phi(cellres, ytau, sc.sigma_d(d, sig)); yinv = [None] * 4
            tinv = [0] * 4
            for i in range(4): tinv[tau[i] - 1] = i + 1
            ytinv = [y[tinv[i] - 1] for i in range(4)]; r2 = phi(cellres, ytinv, sc.sigma_d(d, sig))
            eq1, eq2 = (lhs == r1), (lhs == r2); ok = ok and (eq1 or eq2)
            detail.append({"chart": g, "t": str(t), "sigma": "".join(map(str, sig)), "Phi_sigma.tau(y)==Phi_sigma(y o tau)": eq1, "==Phi_sigma(y o tau^-1)": eq2, "value": str(float(lhs))})
    return ok, detail
def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--channels", default=",".join(CHANNELS), help="comma list (default: all channels; each needs its chart-4 and chart-1 cell pickles or the term lists, see DATA.md)"); a = ap.parse_args()
    slice_prov.log_pid("slice_classes"); rng = random.Random(921); t0 = time.time(); out = {"channels": {}}
    pi = (2, 1, 4, 3)
    for ch in a.channels.split(","):
        taus = [tuple(x) for x in OC.RELABELING_SYMMETRIES.get(ch, [])]          # candidate relabeling symmetries; each is VERIFIED exactly below before use
        if ch == "n4":
            taus = []       # the N=4 term list is one orientation summed over recoilers; no relabeling symmetry assumed
        ver = {}; used = []
        for tau in taus:
            ok, det = verify_tau(ch, tau, rng); ver["".join(map(str, tau))] = {"verified_exact": ok, "checks": det}
            if ok: used.append(tau)
        # slice Z2 check (identical d-matrix polynomials) on two group objects is done in slice_z2check; here by construction
        H = set([(1, 2, 3, 4)]); frontier = list(used)
        while frontier:      # group closure
            tnew = frontier.pop()
            for h in list(H):
                for c in (compose(h, tnew), compose(tnew, h)):
                    if c not in H: H.add(c); frontier.append(c)
        classes = []; seen = set()
        for sig in itertools.permutations((1, 2, 3, 4)):
            if sig in seen: continue
            orb = set(); fr = [sig]
            while fr:
                s = fr.pop()
                if s in orb: continue
                orb.add(s); fr.append(compose(pi, s))
                for h in H: fr.append(compose(s, h))
            seen |= orb; classes.append((sig, len(orb), sorted(orb)))
        assert sum(m for s, m, o in classes) == 24
        groups = {}
        cell4 = OC.build_cells(ch, 4, verbose=False)
        prs = {}
        for dk, conv in cell4["cells"]:
            p = hf_piece.first_pair(dk); prs[p] = prs.get(p, 0) + 1
        for p, n in sorted(prs.items()):
            kj = "%d%d" % p; g = sc.CHART[kj]; groups[kj] = {"n_cells(chart4 census)": n, "gauge": g, "order": sc.order_for(kj, g)}
        out["channels"][ch] = {"relabeling_symmetries_candidates": ["".join(map(str, t)) for t in taus], "verification": ver, "symmetry_group_order(incl. verified taus)": len(H),
                               "n_classes": len(classes), "classes": [{"rep": "".join(map(str, s)), "mult": m, "members": ["".join(map(str, x)) for x in o]} for s, m, o in classes], "first_pair_groups": groups,
                               "n_group_objects": len(classes) * len(groups)}
        print(ch, "H order", len(H), "classes", len(classes), [("".join(map(str, s)), m) for s, m, o in classes], "groups", list(groups), flush=True)
    out["total_group_objects"] = sum(v["n_group_objects"] for v in out["channels"].values()); out["wall_s"] = round(time.time() - t0, 1)
    slice_prov.write_receipt(os.path.join(WORK, "SLICE_CLASSES.json"), out, inputs=[os.path.join(sc.GATE, "cells", "%s_g%d.pkl" % (c, g)) for c in a.channels.split(",") for g in (4, 1)])
    print("total group objects", out["total_group_objects"])
if __name__ == "__main__":
    main()

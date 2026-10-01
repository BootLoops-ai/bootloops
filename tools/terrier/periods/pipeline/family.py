#!/usr/bin/env python3
"""
family.py — model card for an h_eff-parameter effective family.

A card carries everything the generalized restriction + certification pipeline
needs about one effective family / vacuum:
  * charge_vectors: h Mori vectors (entry 0..n_num-1 = degree entries) — these
    determine the effective Gamma-series A(nu) = prod_num Gamma(1+v.nu) /
    prod_den Gamma(1+v.nu) (CY condition per vector asserted);
  * kappa (triple intersections), a_mat, c2D (b = c2D/24), chi_A (xi = zeta3
    chi_A/(2 v^3); the sign of xi is pinned by the stage1 xi-flip mutation
    gate of run_pipe.py), pinned GV leaders + window;
  * PFFV flux data M, K, Q_D3 (curve_from_flux layer);
  * pins: sigma/eps sign frame, s_star, tau_pin (exact rational strings);
  * transport route parameters (s0/rho0/waypoint hints per route);
  * the effective PF ideal + module basis (vacuum layer), sympy strings.
All rational data is exact (Fraction); no floats anywhere in a card.
"""
from fractions import Fraction as Fr
from collections import Counter
import json, os


class Card:
    def __init__(self, d):
        self.raw = d
        self.name = d["name"]
        self.h = int(d["h"])
        self.charge = [list(map(int, v)) for v in d["charge_vectors"]]
        assert len(self.charge) == self.h
        self.n_num = int(d.get("n_num", 1))
        ncol = len(self.charge[0])
        assert all(len(v) == ncol for v in self.charge)
        for a in range(self.h):
            assert sum(self.charge[a]) == 0, f"CY condition fails, vector {a}"
        num, den = Counter(), Counter()
        for i in range(ncol):
            v = tuple(self.charge[a][i] for a in range(self.h))
            if i < self.n_num:
                num[tuple(-x for x in v)] += 1
            else:
                den[v] += 1
        # factors: (v, multiplicity, is_numerator);  gamma-cancellation:
        # sum_num m*v - sum_den m*v = 0 is exactly the CY condition above.
        self.factors = [(v, m, True) for v, m in sorted(num.items())] + \
                       [(v, m, False) for v, m in sorted(den.items())]
        chk = [0] * self.h
        for v, m, isn in self.factors:
            for a in range(self.h):
                chk[a] += (m if isn else -m) * v[a]
        assert all(x == 0 for x in chk), "gamma-cancellation bookkeeping"
        self.kappa = {tuple(sorted(map(int, k.split(",")))): Fr(vv)
                      for k, vv in d["kappa"].items()}
        self.a_mat = [[Fr(x) for x in row] for row in d["a_mat"]]
        self.c2D = [int(x) for x in d["c2D"]]
        self.b_vec = [Fr(x, 24) for x in self.c2D]
        self.chi_A = int(d["chi_A"])
        self.xi_coef = Fr(self.chi_A, 2)          # xi = xi_coef * z3 * v^-3
        self.gv_pinned = {tuple(map(int, k.split(","))): int(v)
                          for k, v in d["gv_pinned"].items()}
        self.gv_window = tuple(int(x) for x in d["gv_window"])
        # non-simplicial effective-cone route (geff_series.eff_route):
        # eff_rays = extreme rays of the CY-effective cone in the card basis;
        # gv_table = provenance-gated GV classes (spec table) in the card
        # basis — REQUIRED when the cone is non-simplicial (gv_extract
        # fail-closes there).  Both optional; absent = legacy simplicial card.
        self.eff_rays = [tuple(int(x) for x in r)
                         for r in d["eff_rays"]] if d.get("eff_rays") else None
        self.gv_table = {tuple(map(int, k.split(","))): Fr(v)
                         for k, v in d["gv_table"].items()} \
            if d.get("gv_table") else None
        self.M = [Fr(x) for x in d["M_flux"]]
        self.K = [Fr(x) for x in d["K_flux"]]
        self.Q_D3 = int(d["Q_D3"])
        # CONDITIONAL-C2 mode (half-integral tadpole, see curve_from_flux):
        # opt-in card flag; REQUIRES a pin of the paper's TRUE -M.K/2 value
        # with a provenance string (the true value recomputed from the
        # paper's flux vectors — NOT the printed value when they differ,
        # cf. arXiv:2107.09064 eq. (6.26) where 83/2 is a typo).  Absent =
        # stock integer-tadpole gate, byte-identical behavior.
        aht = d.get("allow_half_tadpole")
        if aht is not None:
            assert isinstance(aht, dict) and aht.get("paper_true_tadpole") \
                and aht.get("provenance"), \
                "allow_half_tadpole needs paper_true_tadpole + provenance"
            self.allow_half_tadpole = {
                "pin": Fr(aht["paper_true_tadpole"]),
                "provenance": str(aht["provenance"])}
        else:
            self.allow_half_tadpole = None
        self.sigma = tuple(int(x) for x in d["sigma"]) if d.get("sigma") else None
        self.eps = int(d["eps"]) if d.get("eps") else None
        self.s_star = Fr(d["s_star"]) if d.get("s_star") else None
        self.tau_pin = Fr(d["tau_pin"]) if d.get("tau_pin") else None
        self.racetrack = {int(k): Fr(v)
                          for k, v in d.get("racetrack_Nm", {}).items()}
        self.routes = d.get("routes", {})          # per-route transport params
        self.ideal = d.get("ideal")                # sympy strings (vacuum layer)
        self.module_basis = [tuple(map(int, b)) for b in d["module_basis"]] \
            if d.get("module_basis") else None
        self.op_search = d.get("op_search", {})    # rmax/smax/primes/fit_window

    def kap(self, i, j, k):
        return self.kappa.get(tuple(sorted((i, j, k))), Fr(0))


def load_card(path):
    with open(path) as f:
        return Card(json.load(f))

#!/usr/bin/env python3
"""Unit tests for tailcut.py (WP2 certified exponential tail cut).

Controls:
  1. bound VALIDITY: certified bound >= true numeric tail (K^4 and K^2,
     both gammainc branches), and not absurdly loose (< 10^6 x true).
  2. cutoff selection: bound < tol at accepted X; X grows with dps;
     legacy seed kept when already sufficient.
  3. fail-closed: validity-domain raise, cap raise.
  4. eps-fan worst-node semantics.
  5. numeric envelope check at the cutoff (verify=True).
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from mpmath import mp, mpf, mpc

from wayfinder.tailcut import (
    besselk_env_C, besselk_tail_bound, choose_bessel_cutoff,
    TailDomainError, TailNonConvergence, TailEnvelopeError)


def true_tail(nu, sigma, X, m, dps=40):
    with mp.workdps(dps):
        f = lambda x: x ** sigma * mp.besselk(nu, x) ** m
        return abs(mp.quad(f, [X, X + 5, X + 20, mp.inf]))


class TestBound(unittest.TestCase):
    def test_bound_covers_true_tail_K4(self):
        # specimen shape: sigma = 2*eps-1 near -1, nu = 1-eps near 1, m=4
        with mp.workdps(40):
            for X in (5, 10, 20):
                b = besselk_tail_bound(1, -1, X, m=4)
                t = true_tail(1, -1, X, 4)
                self.assertGreater(b, t)
                self.assertLess(b, t * 1e6)

    def test_bound_covers_true_tail_K2_gammainc_branch(self):
        # a = Re(sigma) - m/2 = 3 - 1 = 2 > 0 -> incomplete-Gamma branch
        with mp.workdps(40):
            # nu=1/2: envelope is EXACT (K_{1/2}=sqrt(pi/2x)e^-x, C=0) ->
            # bound == true tail (equality is the sharpest possible case)
            b = besselk_tail_bound(mpf(1) / 2, 3, 12, m=2)
            t = true_tail(mpf(1) / 2, 3, 12, 2)
            self.assertGreaterEqual(b, t * (1 - mpf(10) ** -30))
            self.assertLess(b, t * 1e6)
            # nu=1: strict cover on the same branch
            b1 = besselk_tail_bound(1, 3, 12, m=2)
            t1 = true_tail(1, 3, 12, 2)
            self.assertGreater(b1, t1)
            self.assertLess(b1, t1 * 1e6)

    def test_complex_order_and_sigma(self):
        # complex nu (eps-fan node) and complex sigma: bound must bound the
        # modulus of the complex integral
        with mp.workdps(40):
            eps = mpc(0.01, 0.03)
            nu = 1 - eps
            sigma = 2 * eps - 1
            b = besselk_tail_bound(nu, sigma, 15, m=4)
            f = lambda x: x ** sigma * mp.besselk(nu, x) ** 4
            t = abs(mp.quad(f, [15, 25, 50, mp.inf]))
            self.assertGreater(b, t)

    def test_env_C_explicit(self):
        self.assertEqual(besselk_env_C(mpf(1) / 2), 0)      # |4/4-1|/8
        self.assertEqual(besselk_env_C(1), mpf(3) / 8)

    def test_domain_raise(self):
        with self.assertRaises(TailDomainError):
            besselk_tail_bound(3, -1, 5, m=4)   # X=5 < |nu|^2=9


class TestChoose(unittest.TestCase):
    def test_bound_below_tol_and_dps_scaling(self):
        nodes = [(1, -1)]
        c60 = choose_bessel_cutoff(60, nodes, 50, m=4)
        c120 = choose_bessel_cutoff(120, nodes, 50, m=4)
        self.assertLess(c60.bound, c60.tol)
        self.assertLess(c120.bound, c120.tol)
        self.assertGreater(c120.X, c60.X)      # X = f(dps)

    def test_seed_kept_when_sufficient(self):
        # at dps 30 the legacy X=50 is far more than enough: keep it
        c = choose_bessel_cutoff(30, [(1, -1)], 50, m=4)
        self.assertEqual(c.X, mpf(50))

    def test_seed_lifted_to_validity(self):
        # |nu|^2 = 9 > seed 2: start lifted, no raise
        c = choose_bessel_cutoff(30, [(3, -1)], 2, m=4)
        self.assertGreaterEqual(c.X, 9)
        self.assertLess(c.bound, c.tol)

    def test_cap_raise(self):
        with self.assertRaises(TailNonConvergence):
            choose_bessel_cutoff(200, [(1, -1)], 50, m=4, x_cap=60)

    def test_worst_node_fan(self):
        # node with larger Re(sigma) has the larger tail -> worst
        nodes = [(1, -1), (1, mpc(3, 1))]
        c = choose_bessel_cutoff(50, nodes, 50, m=4)
        self.assertEqual(c.node_index, 1)
        # every node individually satisfies the accepted X
        with mp.workdps(80):
            for nu, sigma in nodes:
                self.assertLess(besselk_tail_bound(nu, sigma, c.X, m=4),
                                c.tol)

    def test_verify_envelope(self):
        c = choose_bessel_cutoff(40, [(1, -1)], 50, m=4, verify=True)
        self.assertTrue(c.verified)

    def test_bound_line_prints(self):
        c = choose_bessel_cutoff(40, [(1, -1)], 50, m=4)
        s = c.bound_line("B tail")
        self.assertIn("[certified]", s)
        self.assertIn("X=", s)


if __name__ == "__main__":
    unittest.main(verbosity=2)

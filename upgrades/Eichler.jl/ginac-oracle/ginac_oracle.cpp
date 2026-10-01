// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by
// Claude (Anthropic) under his supervision. This program links GiNaC/CLN
// (GPL-2.0-or-later); a compiled ginac_oracle binary is a GPL-covered combined work
// and is not redistributed here. GiNaC is not bundled.
//
// ginac_oracle.cpp
//
// Independent cross-check oracle for Eichler.jl, using GiNaC's
// implementation of iterated integrals of modular forms
// (Walden-Weinzierl, arXiv:2010.05271; GiNaC >= 1.8).
//
// Everything is evaluated with Digits = 140 (CLN arbitrary precision);
// truncation cross-checks below guarantee well over 100 correct digits.
//
// Objects evaluated (conventions of AW arXiv:1704.08895 / docs/conventions.md):
//
//  Part 1:  E1(tau; chi0, chi1)      = Eisenstein_kernel(k=1, N=6, a=1, b=-3, K=1)
//           E1(2tau; chi0, chi1)     = Eisenstein_kernel(k=1, N=6, a=1, b=-3, K=2)
//           at qbar = exp(2 pi i tau), tau = 1/10 + i.
//
//  Part 2:  f3 = -3 sqrt(3) [ E3(tau; chi1, chi0) - 8 E3(2tau; chi1, chi0) ]
//              = modular_form_kernel(3, -3*sqrt(3)*(Eisenstein_kernel(3,6,-3,1,1)
//                                                   - 8*Eisenstein_kernel(3,6,-3,1,2)))
//           f4 = (3 sqrt(2) E1(tau; chi0, chi1))^4
//              = modular_form_kernel(4, 324*Eisenstein_kernel(1,6,1,-3,1)^4)
//           "1" = basic_log_kernel()  (the one-form dqbar/qbar)
//
//           I(f3; qbar), I(1, f3; qbar), I(1, f4, 1, f3; qbar)
//           via GiNaC iterated_integral(lst{...}, qbar, N_trunc).
//           List order convention (GiNaC manual 5.12.4 / WW paper):
//           first list element = OUTERMOST integration, i.e.
//           iterated_integral(lst{f1,...,fn}, y) = I(f1,...,fn; y)
//             = int_0^y dy1/y1 f1(y1) I(f2,...,fn; y1),
//           lower limit 0 with the AW tangential-base-point prescription.
//
//  Evaluation points:
//    qbar_t  = -0.079691981555040530318632444762388080332923188146499300401626486412595114957926
//              (exact rational: -79691981555040530318632444762388080332923188146499300401626486412595114957926/10^78,
//               nome at t = -1)
//    qbar_p  = 1/20  = +0.05
//    qbar_c  = 1/50 + (3/100) I = 0.02 + 0.03 i
//
// Compile:  g++ -O2 -o ginac_oracle ginac_oracle.cpp -lginac -lcln

#include <ginac/ginac.h>
#include <ginac/version.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

using namespace GiNaC;

static std::string num2str(const ex& e)
{
	std::ostringstream s;
	s << e;
	return s.str();
}

// pretty-print a (complex) number: Re and Im on separate lines
static void print_value(std::ostream& out, const std::string& label, const ex& val)
{
	ex v = val.evalf();
	if (!is_a<numeric>(v)) {
		out << label << " = (NOT NUMERIC) " << v << "\n";
		return;
	}
	numeric n = ex_to<numeric>(v);
	out << label << "\n";
	out << "  Re = " << n.real() << "\n";
	out << "  Im = " << n.imag() << "\n";
}

static void print_abs(std::ostream& out, const std::string& label, const ex& a, const ex& b)
{
	ex d = (a - b).evalf();
	if (is_a<numeric>(d))
		out << label << " = " << abs(ex_to<numeric>(d)) << "\n";
	else
		out << label << " = (NOT NUMERIC) " << d << "\n";
}

// turn a q_expansion_modular_form result (pseries or polynomial) into a plain
// polynomial truncated at degree N.  NOTE: the topmost coefficient(s) returned
// by q_expansion_modular_form can be incomplete for K>1 kernels (boundary
// truncation artifact), so callers always request a few orders MORE than N.
static ex poly_through(const ex& qexp, const symbol& q, int N)
{
	ex p = qexp;
	if (is_a<pseries>(p))
		p = series_to_poly(p);
	p = p.expand();
	ex r = 0;
	for (int n = 0; n <= N; ++n)
		r += p.coeff(q, n) * pow(q, n);
	return r;
}

// degree-ordered string form of a truncated q-expansion
static std::string poly_ordered_str(const ex& qexp, const symbol& q, int N)
{
	ex p = qexp;
	if (is_a<pseries>(p))
		p = series_to_poly(p);
	p = p.expand();
	std::ostringstream s;
	bool first = true;
	for (int n = 0; n <= N; ++n) {
		ex c = p.coeff(q, n);
		if (c.is_zero())
			continue;
		if (!first) s << " + ";
		s << "(" << c << ")";
		if (n >= 1) s << "*q^" << n;
		first = false;
	}
	if (first) s << "0";
	s << " + O(q^" << (N + 1) << ")";
	return s.str();
}

// extract numeric coefficients 0..N of an (expanded) polynomial in q
static std::vector<numeric> coeff_vector(const ex& poly, const symbol& q, int N)
{
	ex p = poly;
	if (is_a<pseries>(p))
		p = series_to_poly(p);
	p = p.expand();
	std::vector<numeric> v;
	v.reserve(N + 1);
	for (int n = 0; n <= N; ++n) {
		ex c = p.coeff(q, n).evalf();
		if (!is_a<numeric>(c))
			throw std::runtime_error("coefficient is not numeric: " + num2str(c));
		v.push_back(ex_to<numeric>(c));
	}
	return v;
}

// truncated Cauchy product
static std::vector<numeric> conv(const std::vector<numeric>& a, const std::vector<numeric>& b, int N)
{
	std::vector<numeric> c(N + 1, numeric(0));
	for (int i = 0; i <= N; ++i) {
		if (a[i].is_zero()) continue;
		for (int j = 0; j + i <= N; ++j)
			c[i + j] += a[i] * b[j];
	}
	return c;
}

// Horner evaluation of sum_{n=0}^{N} c_n q^n
static numeric horner(const std::vector<numeric>& c, const numeric& q)
{
	numeric s(0);
	for (int n = (int)c.size() - 1; n >= 0; --n)
		s = s * q + c[n];
	return s;
}

int main(int argc, char** argv)
{
	const std::string outpath = (argc > 1) ? argv[1] : "results.txt";
	std::ofstream out(outpath);
	if (!out) { std::cerr << "cannot open " << outpath << "\n"; return 1; }

	Digits = 140;

	out << "============================================================================\n";
	out << "GiNaC cross-check oracle for Eichler.jl  (independent of the Julia code)\n";
	out << "GiNaC version: " << GINACLIB_VERSION << "   (Walden-Weinzierl arXiv:2010.05271 kernels)\n";
	out << "CLN working precision: Digits = " << Digits << " decimal digits\n";
	out << "All printed values are floating point at that precision; the truncation\n";
	out << "cross-checks reported below bound the error far below 1e-100.\n";
	out << "============================================================================\n\n";

	symbol q("q");

	// ------------------------------------------------------------------
	// kernels
	// ------------------------------------------------------------------
	ex E1K1 = Eisenstein_kernel(1, 6, 1, -3, 1);   // E1(tau;   chi0, chi1)
	ex E1K2 = Eisenstein_kernel(1, 6, 1, -3, 2);   // E1(2 tau; chi0, chi1)
	ex E3K1 = Eisenstein_kernel(3, 6, -3, 1, 1);   // E3(tau;   chi1, chi0)
	ex E3K2 = Eisenstein_kernel(3, 6, -3, 1, 2);   // E3(2 tau; chi1, chi0)

	ex P3 = -3 * sqrt(ex(3)) * (E3K1 - 8 * E3K2);
	ex f3 = modular_form_kernel(3, P3);

	ex P4 = pow(3 * sqrt(ex(2)) * E1K1, 4);        // = 324 * E1K1^4
	ex f4 = modular_form_kernel(4, P4);

	ex one = basic_log_kernel();                   // dqbar/qbar

	out << "Kernel definitions (GiNaC objects):\n";
	out << "  E1K1 = " << E1K1 << "\n";
	out << "  E1K2 = " << E1K2 << "\n";
	out << "  E3K1 = " << E3K1 << "\n";
	out << "  E3K2 = " << E3K2 << "\n";
	out << "  f3   = " << f3 << "\n";
	out << "  f4   = " << f4 << "\n";
	out << "  1    = " << one << "  (basic log kernel dq/q)\n\n";

	out << "q-expansions (exact, truncated to the shown order; obtained from\n";
	out << "q_expansion_modular_form with a safety margin, since its topmost\n";
	out << "coefficient can be incomplete for K>1 kernels):\n";
	out << "  E1K1: " << poly_ordered_str(ex_to<Eisenstein_kernel>(E1K1).q_expansion_modular_form(q, 17), q, 12) << "\n";
	out << "  E1K2: " << poly_ordered_str(ex_to<Eisenstein_kernel>(E1K2).q_expansion_modular_form(q, 17), q, 12) << "\n";
	out << "  E3K1: " << poly_ordered_str(ex_to<Eisenstein_kernel>(E3K1).q_expansion_modular_form(q, 17), q, 12) << "\n";
	out << "  E3K2: " << poly_ordered_str(ex_to<Eisenstein_kernel>(E3K2).q_expansion_modular_form(q, 17), q, 12) << "\n";
	out << "  f3  : " << poly_ordered_str(ex_to<modular_form_kernel>(f3).q_expansion_modular_form(q, 17), q, 12) << "\n";
	out << "  f4  : " << poly_ordered_str(ex_to<modular_form_kernel>(f4).q_expansion_modular_form(q, 17), q, 12) << "\n\n";

	// consistency check: f3 in the Sebbar-generator form of docs/conventions.md
	//   f3 = 36 sqrt(3) (e1^3 - e1^2 e2 - 4 e1 e2^2 + 4 e2^3),  e1 = E1K1, e2 = E1K2
	{
		ex P3b = 36 * sqrt(ex(3)) * (pow(E1K1, 3) - pow(E1K1, 2) * E1K2
		                             - 4 * E1K1 * pow(E1K2, 2) + 4 * pow(E1K2, 3));
		ex f3b = modular_form_kernel(3, P3b);
		ex pa = poly_through(ex_to<modular_form_kernel>(f3).q_expansion_modular_form(q, 41), q, 38);
		ex pb = poly_through(ex_to<modular_form_kernel>(f3b).q_expansion_modular_form(q, 41), q, 38);
		out << "Identity check  -3*sqrt(3)*(E3K1-8*E3K2) == 36*sqrt(3)*(e1^3-e1^2*e2-4*e1*e2^2+4*e2^3)\n";
		out << "  difference of q-expansions through q^38: " << (pa - pb).expand() << "\n\n";
	}

	// ==================================================================
	// Part 1: Eisenstein series values at tau = 1/10 + i
	// ==================================================================
	out << "============================================================================\n";
	out << "PART 1: Eisenstein series at tau = 1/10 + i,  qbar = exp(2 pi i tau)\n";
	out << "============================================================================\n";
	ex tau    = numeric(1, 10) + I;
	ex qbar_tau = exp(2 * Pi * I * tau).evalf();
	print_value(out, "qbar = exp(2 pi i (1/10 + i))", qbar_tau);
	out << "\n";

	// values via the kernel's own numerical evaluator (adaptive truncation)
	ex v_E1K1 = ex_to<Eisenstein_kernel>(E1K1).get_numerical_value(qbar_tau);
	ex v_E1K2 = ex_to<Eisenstein_kernel>(E1K2).get_numerical_value(qbar_tau);
	print_value(out, "[1a] E1(tau; chi0, chi1)   = Eisenstein_kernel(1,6,1,-3,1) at qbar", v_E1K1);
	print_value(out, "[1b] E1(2 tau; chi0, chi1) = Eisenstein_kernel(1,6,1,-3,2) at qbar", v_E1K2);

	// independent in-program check: truncated exact q-expansion, order 70
	// (|qbar| = e^{-2 pi} ~ 1.87e-3, so 70 terms give ~ 190 digits)
	{
		ex p1 = poly_through(ex_to<Eisenstein_kernel>(E1K1).q_expansion_modular_form(q, 75), q, 70);
		ex p2 = poly_through(ex_to<Eisenstein_kernel>(E1K2).q_expansion_modular_form(q, 75), q, 70);
		ex s1 = p1.subs(q == qbar_tau).evalf();
		ex s2 = p2.subs(q == qbar_tau).evalf();
		print_abs(out, "  check |[1a] - (q-expansion through q^70)|", v_E1K1, s1);
		print_abs(out, "  check |[1b] - (q-expansion through q^70)|", v_E1K2, s2);
	}
	// also the values of f3 and f4 at this point, for free
	ex v_f3tau = ex_to<modular_form_kernel>(f3).get_numerical_value(qbar_tau);
	ex v_f4tau = ex_to<modular_form_kernel>(f4).get_numerical_value(qbar_tau);
	print_value(out, "[1c] f3(tau) at qbar (same point, bonus)", v_f3tau);
	print_value(out, "[1d] f4(tau) at qbar (same point, bonus)", v_f4tau);
	out << "\n";

	// ==================================================================
	// Part 2: iterated integrals
	// ==================================================================
	out << "============================================================================\n";
	out << "PART 2: iterated integrals I(f3; q), I(1, f3; q), I(1, f4, 1, f3; q)\n";
	out << "  iterated_integral(lst{k1,...,kn}, qbar, N_trunc); first element outermost;\n";
	out << "  base point 0 with tangential base point regularisation (AW prescription).\n";
	out << "  All three integrals here are convergent (f3 has no constant term), so no\n";
	out << "  log(qbar) / branch ambiguity arises at negative real qbar; this is verified\n";
	out << "  numerically below by the real-coefficient series reconstruction.\n";
	out << "============================================================================\n\n";

	// the three evaluation points (exact)
	ex qbar_t = numeric("-79691981555040530318632444762388080332923188146499300401626486412595114957926")
	            / pow(numeric(10), 78);
	ex qbar_p = numeric(1, 20);
	ex qbar_c = numeric(1, 50) + numeric(3, 100) * I;

	struct Pt { const char* name; const char* desc; ex val; };
	std::vector<Pt> points = {
		{"qbar_t", "exact rational -79691981555040530318632444762388080332923188146499300401626486412595114957926/10^78 (nome at t=-1)", qbar_t},
		{"qbar_p", "exact rational 1/20 = +0.05", qbar_p},
		{"qbar_c", "exact 1/50 + (3/100)*I = 0.02 + 0.03 i", qbar_c},
	};

	lst l1 = {f3};
	lst l2 = {one, f3};
	lst l3 = {one, f4, one, f3};

	const int N0 = 96, N1 = 256, N2 = 320;   // truncation orders for stability checks

	// ---- independent in-program reconstruction by direct series algebra ----
	// coefficients a3[n] of f3 and b4[n] of f4 to order N2, built from the
	// exact Eisenstein q-expansions (NOT from iterated_integral):
	const int NC = N2;
	std::vector<numeric> a3, b4, e31, e32, e1v;
	{
		ex p31 = ex_to<Eisenstein_kernel>(E3K1).q_expansion_modular_form(q, NC + 5);
		ex p32 = ex_to<Eisenstein_kernel>(E3K2).q_expansion_modular_form(q, NC + 5);
		ex p1  = ex_to<Eisenstein_kernel>(E1K1).q_expansion_modular_form(q, NC + 5);
		e31 = coeff_vector(p31, q, NC);
		e32 = coeff_vector(p32, q, NC);
		e1v = coeff_vector(p1,  q, NC);
		numeric m3 = ex_to<numeric>((-3 * sqrt(ex(3))).evalf());
		a3.resize(NC + 1);
		for (int n = 0; n <= NC; ++n)
			a3[n] = m3 * (e31[n] - numeric(8) * e32[n]);
		std::vector<numeric> e2 = conv(e1v, e1v, NC);
		std::vector<numeric> e4 = conv(e2, e2, NC);
		b4.resize(NC + 1);
		for (int n = 0; n <= NC; ++n)
			b4[n] = numeric(324) * e4[n];
	}
	// I(f3; q)          = sum c[n] q^n,  c[n] = a3[n]/n
	// I(1, f3; q)       = sum g[n] q^n,  g[n] = c[n]/n
	// I(f4, 1, f3; q)   = sum h[n] q^n,  h[n] = (1/n) sum_{m=1}^{n} b4[n-m] g[m]
	// I(1, f4, 1, f3;q) = sum j[n] q^n,  j[n] = h[n]/n
	std::vector<numeric> cc(NC + 1, numeric(0)), gg(NC + 1, numeric(0)),
	                     hh(NC + 1, numeric(0)), jj(NC + 1, numeric(0));
	for (int n = 1; n <= NC; ++n) {
		cc[n] = a3[n] / numeric(n);
		gg[n] = cc[n] / numeric(n);
	}
	for (int n = 1; n <= NC; ++n) {
		numeric s(0);
		for (int m = 1; m <= n; ++m)
			s += b4[n - m] * gg[m];
		hh[n] = s / numeric(n);
		jj[n] = hh[n] / numeric(n);
	}

	for (const auto& P : points) {
		out << "----------------------------------------------------------------------------\n";
		out << "Evaluation point " << P.name << " = " << P.desc << "\n";
		print_value(out, std::string("  ") + P.name + " (numerical)", P.val);
		out << "\n";

		numeric qn = ex_to<numeric>(ex(P.val).evalf());

		// depth 1
		ex I1z = iterated_integral(l1, P.val, N0).evalf();
		ex I1a = iterated_integral(l1, P.val, N1).evalf();
		ex I1b = iterated_integral(l1, P.val, N2).evalf();
		ex I1d = iterated_integral(l1, P.val).evalf();          // adaptive default
		print_value(out, std::string("[2a] I(f3; ") + P.name + ")          [N_trunc=320]", I1b);
		print_abs(out, "     truncation check |N=320 - N=256|", I1b, I1a);
		print_abs(out, "     sensitivity      |N=320 - N=96| (proves N_trunc honored)", I1b, I1z);
		print_abs(out, "     adaptive default |N=320 - default|", I1b, I1d);
		print_abs(out, "     series reconstruction |GiNaC - sum a3[n]/n q^n|", I1b, horner(cc, qn));

		// depth 2
		ex I2z = iterated_integral(l2, P.val, N0).evalf();
		ex I2a = iterated_integral(l2, P.val, N1).evalf();
		ex I2b = iterated_integral(l2, P.val, N2).evalf();
		ex I2d = iterated_integral(l2, P.val).evalf();
		print_value(out, std::string("[2b] I(1, f3; ") + P.name + ")       [N_trunc=320]", I2b);
		print_abs(out, "     truncation check |N=320 - N=256|", I2b, I2a);
		print_abs(out, "     sensitivity      |N=320 - N=96|", I2b, I2z);
		print_abs(out, "     adaptive default |N=320 - default|", I2b, I2d);
		print_abs(out, "     series reconstruction |GiNaC - sum a3[n]/n^2 q^n|", I2b, horner(gg, qn));

		// depth 4
		ex I3z = iterated_integral(l3, P.val, N0).evalf();
		ex I3a = iterated_integral(l3, P.val, N1).evalf();
		ex I3b = iterated_integral(l3, P.val, N2).evalf();
		ex I3d = iterated_integral(l3, P.val).evalf();
		print_value(out, std::string("[2c] I(1, f4, 1, f3; ") + P.name + ") [N_trunc=320]", I3b);
		print_abs(out, "     truncation check |N=320 - N=256|", I3b, I3a);
		print_abs(out, "     sensitivity      |N=320 - N=96|", I3b, I3z);
		print_abs(out, "     adaptive default |N=320 - default|", I3b, I3d);
		print_abs(out, "     series reconstruction |GiNaC - convolution series|", I3b, horner(jj, qn));
		out << "\n";
	}

	out << "============================================================================\n";
	out << "Notes:\n";
	out << " * 'truncation check' compares iterated_integral with N_trunc=256 vs 320;\n";
	out << "   agreement at the <=1e-130 level certifies > 100 correct digits (a printed\n";
	out << "   0.0 means the two evaluations are bitwise identical at 140-digit precision).\n";
	out << " * 'sensitivity' compares N_trunc=96 vs 320: a small nonzero difference\n";
	out << "   demonstrates that N_trunc is honored and the series has converged.\n";
	out << " * 'adaptive default' compares against iterated_integral(l, qbar) without\n";
	out << "   explicit truncation (GiNaC's internal convergence strategy).\n";
	out << " * 'series reconstruction' recomputes each integral from the exact\n";
	out << "   q-expansion coefficients by elementary series algebra inside this\n";
	out << "   program (a3[n]/n, a3[n]/n^2, and the b4*g convolution), independent\n";
	out << "   of GiNaC's iterated_integral code path.\n";
	out << " * All three iterated integrals are given by convergent power series with\n";
	out << "   REAL coefficients (times the printed prefactors); no log(qbar) terms\n";
	out << "   occur, hence the negative real point qbar_t is unambiguous.\n";
	out << "============================================================================\n";

	out.close();
	std::cout << "wrote " << outpath << "\n";
	return 0;
}

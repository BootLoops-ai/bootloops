"""Arb (python-flint acb) integrator for the LAST integration of a numeric-kinematics one-fold object
      E(x) = sum_i R_i(x) * prod_words ZIP_reg[word_i(x)]        (HyperFLINT nint=2 output; words of weight <= 2)
over x in (0, inf), closing a PIECE against its two-chart piece-oracle; and the engine-independent residue-route fiber function
fibre_residue() used by slice/scripts/slice_gate12.py as the semantics gate of the symbolic one-fold objects.

Why a plain fixed-precision real-axis quadrature fails here (it returns e.g. 1.7e55 / 8.3e41 instead of 3.6e-2):
  (1) exp-sinh nodes on the REAL axis reach x ~ 1e+-41 where the termwise magnitudes |R_i ZIP| ~ x^(1..2) log^2 x exceed the
      sum (~1/x^2) by 120+ digits -> at dps 50 the node values are pure cancellation garbage;
  (2) on the real axis the representation is termwise SINGULAR: poles of R_i of order up to 13 at ~12 positive rational points
      and fiber letters that lie ON the inner contour (positive letters over whole x-intervals), so a real-axis evaluation needs
      a prescription the intermediate object does not carry.
Method:
  (a) the fiber function F(x) = int int (piece integrand) dx_k dx_j is HOLOMORPHIC on the cut plane C minus (-inf,0] (every
      denominator is multilinear with positive coefficients: x*(pos)+(pos) != 0 for |arg x| < pi), so the last contour is
      ROTATED to the ray x = t e^{i theta}: no termwise singularities, letters off the inner contour, and the result must be
      theta-independent and REAL (two built-in checks);
  (b) every node is evaluated in ball arithmetic (acb) with per-node precision escalation until rad(w*E) < 10^-(digits+4)
      (cancellation is measured, not guessed); hyperlog words by the recursive-Hoelder Taylor evaluator ported to acb
      (rounding rigorous; series-truncation allowance 2^-prec added as radius, documented as the non-rigorous component);
  (c) exp-sinh in t with two step sizes (h, h/2) and truncation |s| <= S (tails ~ e^{-(pi/2) sinh S} from the 1/x^2 decay,
      recorded);
  (d) SEMANTICS GATE before integrating: E(x) at complex x is compared with an INDEPENDENT evaluation of the fiber function
      (residue route from the exact cells: analytic x_a-integration by -sum Res[R log(-x_a)] in acb + exp-sinh over x_i),
      which shares no code and no input with HyperFLINT.
Usage: onefold_arb.py --piece PIECE.json --chart g2 [--digits 32 --theta 0.45 --theta2 0.8 --levels 1/10,1/20 --S 4.7] [--out FILE]
PIECE.json (see MANUAL.md, 'piece file'): {channel, sigma, pair_group, d_sigma {"ab": "p/q"}, onefold {g<gauge>: {gauge, order, response}},
optional oracle {route: {value_mid}}}; 'response' is the path of a HyperFLINT eval-json response with nint=2 (last variable left symbolic)."""
import os, sys, json, re, time, math, argparse, itertools, resource
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import sympy
from flint import acb, arb, acb_poly, acb_series, fmpq, fmpq_poly, ctx
from fibrate_dep import zip_num
import oracle_cells as OC, oracle_core as CO, hf_piece, gate_prov
GATE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ----------------------------------------------------------------------------- parsing
def qpoly(expr, x):
    P = sympy.Poly(sympy.expand(expr), x); cs = [sympy.Rational(c) for c in reversed(P.all_coeffs())]
    return fmpq_poly([fmpq(int(c.p), int(c.q)) for c in cs])

def load_onefold(path):
    r = json.loads(open(path).read().strip().splitlines()[-1])
    assert "result" in r and not any(k in r for k in ("failed", "divergent", "error")), [k for k in r if k != "result"]
    var = r["vars"][-1]; x = sympy.Symbol(var)
    terms = []; letters = {}
    for t in r["result"]:
        R = sympy.together(sympy.sympify(t["coef"].replace("^", "**"), locals={var: x})); N, D = sympy.fraction(R)
        words = []
        for w in t["key"]:
            toks = []
            for l in w:
                if l not in letters:
                    e = sympy.together(sympy.sympify(l.replace("^", "**"), locals={var: x}))
                    if e.free_symbols:
                        n_, d_ = sympy.fraction(e); letters[l] = ("V", qpoly(n_, x), qpoly(d_, x))
                    else:
                        q = sympy.Rational(e); letters[l] = ("E", Fr(int(q.p), int(q.q)))
                toks.append(l)
            words.append(tuple(toks))
        terms.append((qpoly(N, x), qpoly(D, x), words))
    return var, terms, letters

# ----------------------------------------------------------------------------- acb hyperlog evaluator (weight-generic)
def A(q):  # Fraction -> acb
    return acb(fmpq(q.numerator, q.denominator))

class PrecisionError(RuntimeError):
    """raised when the current working precision cannot decide a branch side / separate letters; caller doubles ctx.prec."""

class HArb:
    """G(w;1) and ZIP_reg(word) with letters = exact Fractions or acb values (tokens); zip_num semantics."""
    def __init__(self, mode="closed"):
        self.mode = mode          # "closed": weight-1/2 words by log/Li2 closed forms (branch-validity checked per call), "series": Hoelder-Taylor only
        self.vals = {}; self.gcache = {}; self.zcache = {}; self.depth = 0; self.nfallback = 0; self.uncertain = False; self.ncoincident = 0
    def set_node(self, vals):
        self.vals = vals; self.gcache = {}; self.zcache = {}; self.uncertain = False; self.depth = 0
    # tokens: ("E", Fraction) exact letter | ("V", key) acb value self.vals[key] (key hashable; derived letters get derived keys)
    def num(self, tk):
        return A(tk[1]) if tk[0] == "E" else self.vals[tk[1]]
    @staticmethod
    def is_exact(tk, v): return tk[0] == "E" and tk[1] == v
    def mag(self, tk):
        return abs(float(tk[1])) if tk[0] == "E" else float(abs(self.num(tk)).mid())
    def lin(self, tk, al, be):   # al*a + be  (al, be Fractions) as token
        if tk[0] == "E": return ("E", al * tk[1] + be)
        key = ("lin", tk[1], al, be)
        if key not in self.vals: self.vals[key] = A(al) * self.vals[tk[1]] + A(be)
        return ("V", key)
    def mob(self, tk):           # a/(1+a)
        if tk[0] == "E": return ("E", tk[1] / (1 + tk[1]))
        key = ("mob", tk[1])
        if key not in self.vals: a = self.vals[tk[1]]; self.vals[key] = a / (1 + a)
        return ("V", key)
    def key(self, w):
        return tuple(w)

    def series(self, w):
        amin = min(self.mag(tk) for tk in w if not self.is_exact(tk, 0)); depth = len(w)
        NT = int((ctx.prec + 24 + 3 * depth) * math.log(2) / math.log(amin)) + 20 + 5 * depth
        b = self.num(w[-1]); ib = 1 / b; p = acb(1); c = [acb(0)] * (NT + 1)
        for n in range(1, NT + 1):
            p *= ib; c[n] = -p / n
        for tk in reversed(w[:-1]):
            if self.is_exact(tk, 0):
                c = [acb(0)] + [c[n] / n for n in range(1, NT + 1)]
            else:
                av = self.num(tk); ia = 1 / av; cp = [acb(0)] * (NT + 1); s = acb(0); ap = acb(1); iap = acb(1)
                for j in range(1, NT + 1):
                    s += c[j - 1] * ap; iap *= ia; cp[j] = -s * iap / j; ap *= av
                c = cp
        tot = acb(0)
        for v in c[1:]: tot += v
        eps = arb(2) ** (-(ctx.prec - 4))                      # series-truncation allowance (non-rigorous component)
        return tot + acb(arb(0, eps), arb(0, eps)) * (abs(tot) + 1)

    @staticmethod
    def _sgn(x):      # certain sign of an arb: +1 / -1 / 0 (uncertain or zero)
        if x.lower() > 0: return 1
        if x.upper() < 0: return -1
        return 0
    def _closed_ok(self, a1, a2):
        """validity of the dilog closed form along the straight path z in [0,1]:
        the segment v(z) = (z-a1)/(a2-a1) must avoid the closed negative real axis."""
        d = a2 - a1
        if d.contains(0): return False
        v0 = -a1 / d; v1 = (1 - a1) / d
        s0, s1 = self._sgn(v0.imag), self._sgn(v1.imag)
        if s0 != 0 and s0 == s1: return True
        if s0 == 0 or s1 == 0:
            ex0, ex1 = v0.imag.is_zero(), v1.imag.is_zero()          # exactly real endpoints (exact rational letters)
            if (s0 == 0 and not ex0) or (s1 == 0 and not ex1):
                self.uncertain = True; return False               # sign not decidable at this precision
            for vv, exact in ((v0, ex0), (v1, ex1)):
                if exact and not (vv.real.lower() > 0): return False
            return True
        # opposite certain signs: crossing point
        sstar = v0.imag / (v0.imag - v1.imag); R = v0.real + sstar * (v1.real - v0.real)
        return bool(R.lower() > 0)
    def G2_closed(self, a1, a2):
        """G(a1,a2;1) = Li2(u1) - Li2(u0) + log(1-1/a2) log(v1),  u = 1 - v, v(z) = (z-a1)/(a2-a1); None if not valid."""
        if not self._closed_ok(a1, a2): return None
        d = a2 - a1; v0 = -a1 / d; v1 = (1 - a1) / d
        return (1 - v1).polylog(2) - (1 - v0).polylog(2) + (1 - 1 / a2).log() * v1.log()

    def G1(self, w):
        w = tuple(w)
        if not w: return acb(1)
        k = self.key(w)
        if k in self.gcache: return self.gcache[k]
        assert not self.is_exact(w[0], 1) and not self.is_exact(w[-1], 0), w
        for tk in w:
            if tk[0] == "E": assert not (0 < tk[1] < 1), ("exact letter on (0,1)", tk)
            else:
                z = self.num(tk)
                if self._sgn(z.imag) == 0 and not (z.real.upper() < 0 or z.real.lower() > 1):
                    raise PrecisionError("acb letter on the path (0,1) with uncertain side: %s" % z.str(10))
        val = None
        if self.mode == "closed" and len(w) == 1:
            a = self.num(w[0]); val = (1 - 1 / a).log()
        elif self.mode == "closed" and len(w) == 2:
            if self.is_exact(w[1], 1):                       # endpoint letter: G(a1,1;1) = -Li2(1/(1-a1))  (a1 = 0: -zeta(2))
                a1 = self.num(w[0]); val = -acb(2).zeta() if self.is_exact(w[0], 0) else -(1 / (1 - a1)).polylog(2)
            elif self.is_exact(w[0], 0):                     # G(0,a2;1) = -Li2(1/a2)
                a2 = self.num(w[1]); val = -(1 / a2).polylog(2)
            elif self.key(w[:1]) == self.key(w[1:]):        # repeated letter: G(a,a;1) = log(1-1/a)^2 / 2
                a1 = self.num(w[0]); L = (1 - 1 / a1).log(); val = L * L / 2
            else:
                a1, a2 = self.num(w[0]), self.num(w[1])
                if (a2 - a1).contains(0):                    # numerically coincident distinct letters: symmetric approximant (error O(|a1-a2|) <= ball scale)
                    self.ncoincident += 1; val = (1 - 1 / a1).log() * (1 - 1 / a2).log() / 2
                else:
                    val = self.G2_closed(a1, a2)
                if val is None and not any(self.is_exact(t, 1) for t in w):      # swap via shuffle: G(a1,a2) = G(a1)G(a2) - G(a2,a1)
                    alt = self.G2_closed(a2, a1)
                    if alt is not None: val = (1 - 1 / a1).log() * (1 - 1 / a2).log() - alt
                    elif self.uncertain: raise PrecisionError("closed-form branch side undecidable at %d bits" % ctx.prec)
                    else: self.nfallback += 1
        if val is None:
            nz = [self.mag(tk) for tk in w if not self.is_exact(tk, 0)]
            if min(nz) >= 1.98:
                val = self.series(w)
            else:
                self.depth += 1
                if self.depth > 40: self.depth = 0; raise PrecisionError("Hoelder recursion depth > 40 (letter too close to the path/endpoints): %s" % [self.num(t).str(8) if t[0] != "E" else str(t[1]) for t in w])
                n = len(w); tot = acb(0)
                for kk in range(n + 1):
                    left = tuple(self.lin(tk, Fr(-2), Fr(2)) for tk in reversed(w[:kk]))    # 2(1-a)
                    right = tuple(self.lin(tk, Fr(2), Fr(0)) for tk in w[kk:])              # 2a
                    tot += (-1) ** kk * self.G1(left) * self.G1(right)
                val = tot; self.depth -= 1
        self.gcache[k] = val
        return val

    # ---- zip_num.G_one semantics (shuffle-regularize leading exact 1's, strip trailing exact 0's)
    def reg_leading_ones(self, w):
        k = 0
        while k < len(w) and self.is_exact(w[k], 1): k += 1
        if k == 0: return {w: Fr(1)}
        v = w[k:]
        if not v: return {}
        sh = zip_num.shuffle(tuple(("E", Fr(1)) for _ in range(k)), v); acc = {}
        for u, c in sh.items():
            if u == w: continue
            for uu, cc in self.reg_leading_ones(u).items(): acc[uu] = acc.get(uu, Fr(0)) - Fr(c) * cc
        return {u: c for u, c in acc.items() if c}
    def strip_trailing_zeros(self, w):
        if not w or not self.is_exact(w[-1], 0): return {w: Fr(1)}
        if all(self.is_exact(tk, 0) for tk in w): return {}
        u = w[:-1]; ins = {}
        for i in range(len(u) + 1):
            cand = u[:i] + (("E", Fr(0)),) + u[i:]; ins[cand] = ins.get(cand, 0) + 1
        t = ins.pop(w); out = {}
        for cand, m in ins.items():
            for v, c in self.strip_trailing_zeros(cand).items(): out[v] = out.get(v, Fr(0)) - Fr(m) * c / t
        return out
    def G_one(self, w):
        tot = acb(0)
        for w1, c1 in self.strip_trailing_zeros(tuple(w)).items():
            for w2, c2 in self.reg_leading_ones(w1).items():
                tot += A(c1 * c2) * self.G1(w2)
        return tot
    def zip_value(self, word):
        """ZIP_reg[a1..an] = int_0^inf reg: x=t/(1-t): each letter a -> G-letters {1 (coef -1)} + {a/(1+a) (coef +1) unless a == -1}."""
        word = tuple(word); k = self.key(word)
        if k in self.zcache: return self.zcache[k]
        exps = [((), Fr(1))]
        for tk in word:
            opts = [(("E", Fr(1)), Fr(-1))]
            if not self.is_exact(tk, -1):
                opts.append((self.mob(tk), Fr(1)))
            exps = [(ww + (l,), c * s) for (ww, c) in exps for (l, s) in opts]
        tot = acb(0)
        for ww, c in exps: tot += A(c) * self.G_one(ww)
        self.zcache[k] = tot
        return tot

class OneFold:
    def __init__(self, path):
        self.path = path; self.var, self.terms, self.letters = load_onefold(path); self.H = HArb(); self._pc = {}
    def _polys(self):
        if ctx.prec not in self._pc:
            cv = lambda p: acb_poly([acb(c) for c in p.coeffs()])
            self._pc[ctx.prec] = ([(cv(N), cv(D), words) for N, D, words in self.terms], {l: (v[0],) + ((cv(v[1]), cv(v[2])) if v[0] == "V" else (v[1],)) for l, v in self.letters.items()})
        return self._pc[ctx.prec]
    def E(self, x):
        terms, lets = self._polys()
        vals = {}; toks = {}
        for l, v in lets.items():
            if v[0] == "V": vals[l] = v[1](x) / v[2](x); toks[l] = ("V", l)
            else: toks[l] = ("E", v[1])
        self.H.set_node(vals)
        tot = acb(0); mx = arb(0)
        for N, D, words in terms:
            v = N(x) / D(x)
            for w in words: v *= self.H.zip_value(tuple(toks[l] for l in w))
            tot += v; av = abs(v)
            if av.upper() > mx.upper(): mx = av.upper()
        return tot, mx      # value ball, max |term| (cancellation monitor)

# ----------------------------------------------------------------------------- independent residue-route fiber function (acb)
class RouteEvalC(CO.RouteEval):
    """acb port of RouteEval.outer_setup/node: complex x_o allowed (poles in x_a complex, principal log(-x_a))."""
    def outer_setup_c(self, xo):
        Pi = [acb_poly([acb_poly([acb(c) for c in p.coeffs()])(xo) if p.degree() >= 0 else acb(0) for p in row]) for row in self.T]
        ats = []
        for Aa in self.atoms:
            def topoly(dct):
                if not dct: return None
                deg = max(ei for (eo, ei) in dct); cs = [acb(0)] * (deg + 1)
                for (eo, ei), c in dct.items(): cs[ei] = cs[ei] + acb(c) * xo ** eo
                return acb_poly(cs)
            ats.append((topoly(Aa["alpha"]), topoly(Aa["beta"])))
        return Pi, ats
    def node_c(self, setup, xi):
        Pi, ats = setup; ncoef = [p(xi) for p in Pi]; C = acb(1); pole_ab = []
        for Aa, (alp, bep) in zip(self.atoms, ats):
            if Aa["afree"]: C *= bep(xi) ** Aa["M"]
            elif Aa["is_xa"]: pass
            else: pole_ab.append((alp(xi), bep(xi), Aa["M"]))
        M0 = self.M0
        for k in range(min(M0, len(ncoef))):
            if not ncoef[k].contains(0): raise RuntimeError("x_a^-%d pole does not cancel" % (M0 - k))
        nt = ncoef[M0:]; top = self.degD_rest - 2
        for k in range(top + 1, len(nt)):
            if not nt[k].contains(0): raise RuntimeError("degree excess at x_a^inf")
        nt = nt[:top + 1]; Npoly = acb_poly(nt) if nt else acb_poly([acb(0)]); total = acb(0); npole = len(pole_ab)
        for k in range(npole):
            al, be, M = pole_ab[k]; p = -be / al
            sh = Npoly(acb_poly([p, acb(1)])); cs = [sh[i] for i in range(sh.degree() + 1)] if sh.degree() >= 0 else [acb(0)]
            cs = list(cs)[:M] + [acb(0)] * max(0, M - len(cs)); ser = acb_series(cs, M)
            for j in range(npole):
                if j == k: continue
                al2, be2, M2 = pole_ab[j]; fac = acb_series([al2 * p + be2, al2], M); ser = ser * (fac.inv() ** M2 if M2 > 1 else fac.inv())
            if M0: ser = ser * (acb_series([p, acb(1)], M).inv() ** M0)
            L = acb_series([-p, acb(-1)], M).log(); ser = ser * L
            total += ser[M - 1] / al ** M
        return -total / C

_CELLCACHE = {}
def fibre_residue(piece, chart_gauge, last, xo, h, S, digits, prec0=256, precmax=8192, swap_roles=False):
    """F(x_o) = int_0^inf dx_i g(x_o,x_i) at complex x_o by the residue route (independent of HyperFLINT)."""
    ch = piece["channel"]; sig = tuple(int(c) for c in piece["sigma"]); kj = piece["pair_group"]
    pair = "all" if kj == "all" else tuple(sorted(int(c) for c in kj))
    dsig = {(int(k[0]), int(k[1])): Fr(v) for k, v in piece["d_sigma"].items()}
    key = (ch, chart_gauge)
    if key not in _CELLCACHE: _CELLCACHE[key] = OC.build_cells(ch, chart_gauge, verbose=False)
    cellres = _CELLCACHE[key]
    import piece_select
    sub = piece_select.select_cells(cellres, pair, piece.get("support_atoms"))
    rf = CO.InstanceRF(sub, dsig)
    others = [k for k in rf.xvars if k != last]; i_, a_ = (others[1], others[0]) if swap_roles else (others[0], others[1])
    rev = RouteEvalC(rf, last, i_, a_)
    c = arb.pi() / 2; K = int(math.ceil(S / float(h))); hh = arb(fmpq(h.numerator, h.denominator))
    tot = acb(0); stats = {"nodes": 0, "max_prec": prec0}
    xc = abs(xo).sqrt()                   # inner map centered at sqrt|x_o| (oracle_core convention)
    for k in range(-K, K + 1):
        prec = prec0
        while True:
            ctx.prec = prec
            xo_p = acb(xo)
            setup = rev.outer_setup_c(xo_p)
            s = hh * (arb(k) + arb(1) / 2); xi = abs(xo_p).sqrt().mid() * (c * s.sinh()).exp(); w = hh * c * s.cosh() * xi
            try:
                g = rev.node_c(setup, acb(xi)) * w
                ok = g.rad() < 10.0 ** (-(digits + 4))
            except RuntimeError:
                ok = False
            if ok or prec >= precmax: break
            prec *= 2
        stats["max_prec"] = max(stats["max_prec"], prec); stats["nodes"] += 1
        tot += g
    ctx.prec = prec0
    return tot, stats, (i_, a_)

# ----------------------------------------------------------------------------- ray quadrature of E
def ray_integral(OFo, theta, h, S, digits, prec0=256, precmax=16384, log=None):
    c = math.pi / 2; K = int(math.ceil(S / float(h))); tot = acb(0)
    st = {"nodes": 0, "max_prec": prec0, "escalations": 0, "max_cancel_digits": 0.0, "prec_hist": {}, "precmax_hits": 0}
    prec = prec0
    for idx, k in enumerate(range(-K, K + 1)):
        if idx % 6 == 0 and prec > prec0: prec //= 2          # sticky precision, relaxed every 6th node
        while True:
            ctx.prec = prec
            hh = arb(fmpq(h.numerator, h.denominator)); s = hh * (arb(k) + arb(1) / 3)
            t = (arb(c) * s.sinh()).exp(); w = hh * arb(c) * s.cosh() * t
            eth = acb(0, 1) * theta; eth = eth.exp()
            x = eth * t
            try:
                val, mx = OFo.E(x); g = val * w * eth
                if g.rad() < 10.0 ** (-(digits + 4)): break
                if prec >= precmax: st["precmax_hits"] += 1; break
            except PrecisionError:
                if prec >= precmax: raise
            prec *= 2; st["escalations"] += 1
        st["prec_hist"][prec] = st["prec_hist"].get(prec, 0) + 1
        canc = float((mx / (abs(val).mid() if abs(val).mid() != 0 else arb(1))).log().mid()) / math.log(10) if mx != 0 else 0.0
        st["max_cancel_digits"] = max(st["max_cancel_digits"], canc); st["max_prec"] = max(st["max_prec"], prec); st["nodes"] += 1
        tot += g
        if log and idx % 20 == 0: log("   node %d/%d t=%s prec=%d cancel=%.1f dig  partial=%s" % (idx, 2 * K + 1, t.mid().str(6), prec, canc, tot.real.mid().str(20)))
    ctx.prec = prec0
    st["hlog_closed_form_fallbacks"] = OFo.H.nfallback; st["hlog_coincident_letter_approximants"] = OFo.H.ncoincident
    return tot, st

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--piece", required=True); ap.add_argument("--chart", default="g2", help="which one-fold object of the piece json (g<gauge>)")
    ap.add_argument("--digits", type=int, default=32); ap.add_argument("--theta", type=float, default=0.45); ap.add_argument("--theta2", type=float, default=0.8)
    ap.add_argument("--levels", default="1/10,1/20"); ap.add_argument("--S", type=float, default=4.7); ap.add_argument("--prec0", type=int, default=256)
    ap.add_argument("--gate-points", dest="gpts", default="0.6+0.4j,2.2+1.7j"); ap.add_argument("--out", default=None); ap.add_argument("--skip-theta2", action="store_true")
    a = ap.parse_args()
    T0 = time.time(); logf = lambda m: print("[%7.1fs] %s" % (time.time() - T0, m), flush=True)
    piece = json.load(open(a.piece)); of = piece["onefold"][a.chart]; gauge = of["gauge"]; last = of["order"][-1]
    OFo = OneFold(of["response"])
    res = {"piece_json": os.path.abspath(a.piece), "channel": piece["channel"], "sigma": piece["sigma"], "pair_group": piece["pair_group"], "shape": piece.get("shape"), "shape_str": piece.get("shape_str"),
           "onefold_response": of["response"], "onefold_sha256": gate_prov.sha(of["response"]), "chart_gauge": gauge, "last_var": "x%d" % last, "n_terms": len(OFo.terms), "n_letters": len(OFo.letters),
           "digits_target": a.digits, "prec0": a.prec0}
    # (d) semantics gate: E(x) vs independent residue-route fiber function at complex points
    gate = []
    for ps in a.gpts.split(","):
        z = complex(ps.replace(" ", "")); x0 = acb(z.real, z.imag)
        t0 = time.time(); prec = a.prec0
        while True:
            ctx.prec = prec
            try:
                e, mx = OFo.E(acb(z.real, z.imag))
                if e.rad() < 10.0 ** (-(a.digits + 3)) or prec >= 8192: break
            except PrecisionError:
                if prec >= 8192: raise
            prec *= 2
        we = time.time() - t0; eprec = prec
        t0 = time.time(); fvals = {}
        for lv in ("1/16", "1/32"):
            fv, fst, roles = fibre_residue(piece, gauge, last, acb(z.real, z.imag), Fr(lv), 4.7, a.digits, prec0=a.prec0); fvals[lv] = (fv, fst)
        wf = time.time() - t0
        f1, f2 = fvals["1/16"][0], fvals["1/32"][0]
        ctx.prec = 512
        rel = abs(e - f2) / abs(f2); lvd = abs(f1 - f2) / abs(f2)
        row = {"x": ps, "E_onefold": e.str(40, radius=True), "E_prec_bits": eprec, "E_wall_s": round(we, 1), "max_term_over_sum_log10": round(float((mx / abs(e).mid()).log().mid()) / math.log(10), 2),
               "F_residue_route(h=1/32)": f2.str(40, radius=True), "F_residue_level_diff(1/16 vs 1/32)": lvd.mid().str(5, radius=False), "F_wall_s": round(wf, 1), "F_roles_i_a": list(roles),
               "rel_diff_E_vs_F": rel.mid().str(5, radius=False), "agree_digits": (int(math.floor(-math.log10(float(rel.mid().str(10, radius=False))))) if rel.mid() != 0 else 99)}
        gate.append(row); logf("semantics gate x=%s: E=%s F=%s rel=%s (%s digits)" % (ps, e.str(25), f2.str(25), row["rel_diff_E_vs_F"], row["agree_digits"]))
    res["semantics_gate(onefold_vs_residue_route_at_complex_x)"] = gate
    # (c) ray quadrature at theta, two levels; then theta2 at the finer level
    runs = {}
    levels = [Fr(x) for x in a.levels.split(",")]
    for th, lvs in ((a.theta, levels), (a.theta2, levels[-1:] if not a.skip_theta2 else [])):
        for lv in lvs:
            t0 = time.time(); OFo.H = HArb()
            tot, st = ray_integral(OFo, th, lv, a.S, a.digits, prec0=a.prec0, log=logf); w = time.time() - t0
            ctx.prec = 512
            runs["theta=%g,h=%s" % (th, lv)] = {"theta": th, "h": str(lv), "S": a.S, "Re": tot.real.str(42, radius=True), "Im": tot.imag.str(8, radius=True), "Re_mid": tot.real.mid().str(50, radius=False), "abs_Im_mid": abs(tot.imag).mid().str(5, radius=False), "wall_s": round(w, 1), **st}
            logf("ray theta=%g h=%s: Re=%s Im=%s wall=%.0fs maxprec=%d cancel<=%.0f dig" % (th, lv, tot.real.str(36), tot.imag.str(5), w, st["max_prec"], st["max_cancel_digits"]))
    res["ray_integrals"] = runs
    # comparisons
    ctx.prec = 512
    keys = list(runs)
    fin = runs["theta=%g,h=%s" % (a.theta, levels[-1])]; v_fin = arb(fin["Re_mid"])
    cmp_ = {}
    if len(levels) > 1:
        v_coarse = arb(runs["theta=%g,h=%s" % (a.theta, levels[0])]["Re_mid"]); cmp_["level_diff_rel(h vs h/2, theta)"] = (abs(v_coarse - v_fin) / abs(v_fin)).mid().str(5, radius=False)
    if not a.skip_theta2:
        v_t2 = arb(runs["theta=%g,h=%s" % (a.theta2, levels[-1])]["Re_mid"]); cmp_["theta_independence_rel(theta vs theta2)"] = (abs(v_t2 - v_fin) / abs(v_fin)).mid().str(5, radius=False)
    cmp_["abs_Im_over_Re(theta, finest)"] = (arb(fin["abs_Im_mid"]) / abs(v_fin)).mid().str(5, radius=False)
    if piece.get("oracle"):
        orc = {r: arb(piece["oracle"][r]["value_mid"]) for r in piece["oracle"]}
        v_ball = arb(fin["Re"])
        for r, ov in orc.items():
            rel = (abs(v_ball - ov) / abs(ov)).upper(); cmp_["rel_diff_vs_piece_oracle_route_%s(upper bound incl. ball radius)" % r] = rel.str(5, radius=False)
        relB = float((abs(v_ball - orc.get("B", list(orc.values())[0])) / abs(list(orc.values())[0])).upper().str(10, radius=False))
        cmp_["agree_digits_vs_piece_oracle"] = int(math.floor(-math.log10(relB))) if relB > 0 else 99
        cmp_["piece_oracle_A_vs_B_rel"] = piece.get("A_vs_B_rel")
    if len(levels) > 2:
        v_mid = arb(runs["theta=%g,h=%s" % (a.theta, levels[-2])]["Re_mid"]); cmp_["level_diff_rel(finest two levels, theta)"] = (abs(v_mid - v_fin) / abs(v_fin)).mid().str(5, radius=False)
    res["value_onefold_arb"] = fin["Re"]; res["comparisons"] = cmp_
    res["method"] = {"contour": "ray x = t exp(i theta), t in (0,inf) by exp-sinh t = exp((pi/2) sinh((k+1/3)h)), |k h| <= S; fiber function holomorphic on C minus (-inf,0]",
                     "hyperlogs": "weight<=2 ZIP_reg words (zip_num semantics: x=t/(1-t), leading-1/trailing-0 shuffle regularization) -> G(a;1)=log(1-1/a), G(a1,a2;1)=Li2(1-v1)-Li2(1-v0)+log(1-1/a2)log(v1) with v(z)=(z-a1)/(a2-a1), used only when the segment v([0,1]) avoids (-inf,0] (else the shuffle-swapped ordering; proof: one ordering is always valid off the path); endpoint/zero/repeated letters special-cased; Hoelder-Taylor series evaluator kept as cross-check mode (--mode series)",
                     "balls": "acb throughout; node accepted when rad(w E) < 10^-(digits+4), precision doubled otherwise (sticky, relaxed every 6th node); series-truncation allowance 2^-prec is the only non-rigorous radius component (closed-form mode has none)",
                     "why_not_real_axis_fixed_precision": ["real-axis exp-sinh nodes reach x~1e+-41 where termwise |R_i ZIP| exceed the O(1/x^2) sum by up to ~156 digits (measured here: max_cancel_digits) -> dps 50 returned cancellation garbage (1.7e55, 8.3e41)", "on the real axis the one-fold representation is termwise singular (poles of R_i up to order 13 at 12 positive rational points; fiber letters ON the inner contour over whole x-intervals), so real-axis evaluation needs a prescription the intermediate object does not carry; on the rotated ray every letter is off the contour with a certain side"]}
    if a.skip_theta2:
        res["PASS(>=30 digits vs piece-oracle and theta-independent and Im~0)"] = None; res["PASS_note"] = "theta2 skipped: theta-independence not tested, no PASS label issued"
    else:
        res["PASS(>=30 digits vs piece-oracle and theta-independent and Im~0)"] = bool(cmp_.get("agree_digits_vs_piece_oracle", 0) >= 30 and float(cmp_.get("theta_independence_rel(theta vs theta2)", "0")) < 1e-30 and float(cmp_["abs_Im_over_Re(theta, finest)"]) < 1e-30)
    res["certified_digits_note"] = "agree_digits use the UPPER bound of |closure ball - oracle mid|/|oracle| (ball radius included)"
    res["wall_total_s"] = round(time.time() - T0, 1); res["self_maxrss_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    res["producer"] = gate_prov.producer(inputs=[of["response"], os.path.abspath(a.piece)], modules=["onefold_arb.py", "oracle_cells.py", "oracle_core.py", "hf_piece.py", "gate_prov.py"])
    out = a.out or os.path.join(gate_prov.wdir("onefold"), "ONEFOLD_ARB_%s_s%s_p%s_%s.json" % (piece["channel"], piece["sigma"], piece["pair_group"], a.chart))
    json.dump(res, open(out, "w"), indent=1); logf("wrote %s" % out)
    print(json.dumps(res["comparisons"], indent=1))

if __name__ == "__main__":
    main()

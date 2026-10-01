#!/usr/bin/env python3
"""posq.py -- posq (positive-quadrature evidence) packaged tool surface. [GATED]

FLAG: GATED -- packaged from a validated build ladder plus a 15-topology
production table (the table itself is not shipped here). Promotion =
reviewed merge. See README.md beside this file for the object class and the
BINDING CAVEATS.

WHAT THIS FILE IS: the packaging surface for the posq engine. The engine
itself ships VERBATIM beside this file:

  posq_kernel.c / posq_kernel   C flint/arb multi-topology certified positive
                                Gauss-Jacobi tensor sweep kernel (cat/bal/ms
                                modes; per-node VALUE evaluation P_y = A_y +
                                B_y*v, never coefficient expansion; positive-
                                weighted ball sum).
  kernel_io.py                  exact ball serialization python <-> kernel.
  bench_rung1_gauss.py          certified Gauss-Jacobi rule builder (interval-
                                Newton at prec 2048, receipt gates), python
                                reference sweep, exact-rational N=8 surrogate.
  topos.py exact_polys.py exact_polys_bal.py bal_fiber.py   exact structure
                                (15-topology enumeration, pattern polynomials,
                                balanced min/diff decomposition).
  derive_posq_sentences.py      certified-sentence derivation law (ROUND_FLOOR
                                / ROUND_CEILING display discipline).
The production driver (run_quartet.py), its panel specs, result tables, and
the production pattern-count data file (QUARTET_COLLAPSE.json) are NOT
included in this distribution; modules that would read the data file refuse
by name unless POSQ_STAGE0 points at a directory containing it.

THIS MODULE provides, for consumers and for the battery rerun law:
  * SurrogateExact  -- the N=8 exact-rational surrogate object Z(p, lambda)
    (the pinned rung-1 cross-check object: counts 0011:6, 1110:1, 1010:1,
    rules (17,13,5), degrees (32,24,16), y3-parity exact) as an exact
    closed form: symbolic-in-p expansion of F over (u1,u2,v), exact Fraction
    evaluation, and exact-derivative jets (dZ/dp, dZ/dlambda) evaluated in
    ball arithmetic. Gate: at (p0, lambda) = (709/2048, 10) the Fraction
    value reproduces the pinned RUNG1 Z_exact_str EXACTLY (== on Fraction;
    sha256 of the "num/den" string pinned below).
  * PosqKernelAdaptation -- the eras-battery adaptation (consumer protocol of
    tools/eras/verify_eras_adversarial.py) with the C kernel as the
    point-value engine:
      point_eval(p, prec) = posq_kernel sweep at the point's exact dyadic
        (p, lambda) mids -- rules (17,13,5) rebuilt + receipt-gated at each
        exact rational c = 2*lambda*p*(1-p) (the production register's
        rebuild-per-p law), PREC = prec.
      enclose(p0, radii, K) = rigorous CLOSED-BOX enclosure of Z over
        p0 +/- radii, vector (p, lambda), via the MONOTONE c-BRACKET x
        ball-q register: I(c,q) = int u1^(c-1) u2^(c-1) v^(c/2-1) F(u;q) du
        is completely monotone (decreasing) in c because ln u <= 0 and
        F >= 0, so with c_lo = c(p_lo,l_lo) <= c(p,l) <= c_hi = c(p_hi,l_hi)
        (valid on p-boxes inside (0, 1/2), lambda > 0),
          Z in [ (c_lo^3/2) * lb I(c_hi, Q),  (c_hi^3/2) * ub I(c_lo, Q) ],
        where I(c, Q) is the exact positive Gauss sum swept in python ball
        arithmetic (same ab_fiber chain as the kernel) with the pi1-channel
        p replaced by the ball Q covering [p_lo, p_hi]. Every K returns this
        one form (protocol-permitted cap). ADAPTATION SCAFFOLDING NOTE: this
        enclosure form exists to feed the battery; it is NOT a shipped
        production p-layer (production p-integration is the arm-C ellipse
        register; pointwise production-N p-enclosures were KILLed on price
        at rung 3 -- see README caveats).
      deriv(p, dim) = exact closed-form partial (dZ/dp or dZ/dlambda) of the
        surrogate rational function, evaluated in ball arithmetic at prec
        1280 -- an INDEPENDENT derivative path (exact calculus, no
        quadrature), so the battery's FD arm cross-checks the kernel against
        exact calculus at the 1.1e-12 gate.
      degraded() = the same adaptation with the rule builder forced to
        53-bit. Measured behavior: interval-Newton certification FAILS at
        node 0 (contraction impossible at 53-bit arithmetic noise), so the
        degraded enclosure is NON-FINITE -- the battery counts that as
        caught (no-information can never pass). posq has no import-time
        float-constant pathway (all inputs are exact rationals); insufficient
        build precision is stopped by the rule receipt gates, and this probe
        demonstrates the battery sees it.

WHY THE SURROGATE OBJECT: the battery needs ~700 point evaluations per run;
a production 318-count sweep is 48.7M nodes = 0.2 CPU-h per point (measured
at rung 1). The N=8 surrogate runs the SAME rule builder at the same spec
(prec 2048, same receipt gates), the SAME C kernel in the same mode, the
SAME evaluation chain -- only the pattern counts and rule sizes shrink, and
its exact rational value is a pinned cross-check target. Production-scale
validation of the kernel is separate and already recorded: exact-rational
surrogate containment (cat + bal), ball overlap vs the python engine on
subranges, full-sweep containment of the recorded rung-1 ball at 1e-13 nats,
and 3 planted count corruptions firing (recorded with the production table,
which is not shipped here).

Usage: python3 verify_posq_adaptation.py (beside this file) runs the posq
gates + the eras battery against this adaptation and prints the verdict.
Quick check: python3 posq.py --selftest runs the exact-surrogate pin and one
live kernel-containment leg (seconds).
"""
import hashlib
import os
import sys
import tempfile
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
# guarded sys.path insert — only if absent, and REMOVED
# after our own imports complete (below), so generic names here (kernel_io)
# can no longer shadow user modules imported after `import posq`/
# baller.quad.
_PATH_INSERTED = HERE not in sys.path
if _PATH_INSERTED:
    sys.path.insert(0, HERE)

import bench_rung1_gauss as bg          # noqa: E402
import exact_polys as ep                # noqa: E402
import kernel_io as kio                 # noqa: E402
from flint import arb, ctx              # noqa: E402
if _PATH_INSERTED and HERE in sys.path:
    sys.path.remove(HERE)  # leave no residue

# ---- pinned surrogate object (rung-1 cross-check object) ----
SURR_COUNTS = [((0, 0, 1, 1), 6), ((1, 1, 1, 0), 1), ((1, 0, 1, 0), 1)]
SURR_VEC16 = [0] * 16
for _x, _n in SURR_COUNTS:
    SURR_VEC16[(_x[0] << 3) | (_x[1] << 2) | (_x[2] << 1) | _x[3]] = _n
P0 = Fr(709, 2048)                       # pinned dyadic p0
LAM0 = Fr(10)                            # pinned Exp rate
# sha256 of the rung-1 surrogate exact value as "num/den",
# pinned at package time:
Z_EXACT_STR_SHA256_PIN = \
    "e6fc044cdd0dabfbb3137bc6494bd3cf86eb29c24b988568a1a154660f0b7bb1"
LN_Z_EXACT_PIN = -25.97566273247789      # pinned rung-1 surrogate ln Z


def frs(s):
    """exact rational from 'a/b' or decimal-float string / float."""
    if isinstance(s, Fr):
        return s
    if isinstance(s, float):
        return Fr(s)                     # exact dyadic
    t = str(s)
    if "/" in t:
        n, d = t.split("/")
        return Fr(int(n), int(d))
    return Fr(float(t))                  # exact dyadic of the parsed float


def c_of(p, lam):
    return 2 * lam * p * (1 - p)


def arb_of_fr(q, prec=None):
    old = ctx.prec
    if prec:
        ctx.prec = prec
    try:
        return arb(q.numerator) / arb(q.denominator)
    finally:
        ctx.prec = old


def fr_of_mid(x):
    """exact Fraction of an arb midpoint (dyadic)."""
    m = x.mid()
    mm, me = m.man_exp()
    mm, me = int(mm), int(me)
    return Fr(mm) * Fr(2) ** me


class SurrogateExact:
    """N=8 surrogate Z(p, lambda) as an exact closed form.
    Z = (c^3/2) * sum_{(i,j,l)} cf_{ijl}(p) / ((c+i)(c+j)(c/2+l)),
    c = 2*lambda*p*(1-p); cf from the symbolic-in-p expansion of
    F = prod pat_poly_uv(x)^n over the surrogate counts (exact_polys
    machinery, v = u3^2 parity-reduced). Built once (~2 s)."""

    def __init__(self):
        F = ep.pconst({0: Fr(1)})
        for x, n in SURR_COUNTS:
            F = ep.pmul(F, self._ppow(ep.pat_poly_uv(x), n))
        self.F = F
        degs = (max(m[0] for m in F), max(m[1] for m in F),
                max(m[2] for m in F))
        assert degs == (32, 24, 8), f"surrogate v-degrees {degs}"

    @staticmethod
    def _ppow(poly, n):
        r = ep.pconst({0: Fr(1)})
        b = poly
        while n:
            if n & 1:
                r = ep.pmul(r, b)
            b = ep.pmul(b, b) if n > 1 else b
            n >>= 1
        return r

    def Z_fr(self, p, lam):
        """exact Fraction value (use for gates; big-denominator p is slow)."""
        c = c_of(p, lam)
        s = Fr(0)
        for (i, j, l), cp in self.F.items():
            cf = sum(v * p ** k for k, v in cp.items())
            s += cf / ((c + i) * (c + j) * (c / 2 + l))
        return c ** 3 * s / 2

    def jet_arb(self, p, lam, prec=1280):
        """(Z, dZ/dp, dZ/dlam) balls of the exact closed form at exact
        rational (p, lam); rigorous ball evaluation at prec."""
        old = ctx.prec
        ctx.prec = prec
        try:
            c = c_of(p, lam)
            cA = arb_of_fr(c)
            cp_A = arb_of_fr(2 * lam * (1 - 2 * p))     # dc/dp
            cl_A = arb_of_fr(2 * p * (1 - p))           # dc/dlam
            pA = arb_of_fr(p)
            # p powers
            pdeg = max(max(cp) for cp in self.F.values())
            pw = [arb(1)]
            for _ in range(pdeg):
                pw.append(pw[-1] * pA)
            # denominator inverse tables
            di = max(m[0] for m in self.F)
            dj = max(m[1] for m in self.F)
            dl = max(m[2] for m in self.F)
            invA = [1 / (cA + i) for i in range(di + 1)]
            invB = [1 / (cA + j) for j in range(dj + 1)]
            invC = [1 / (cA / 2 + l) for l in range(dl + 1)]
            S = arb(0)      # sum cf * D^-1
            Sp = arb(0)     # sum cf' * D^-1
            Sc = arb(0)     # sum cf * dD^-1/dc
            for (i, j, l), cp in self.F.items():
                Dm1 = invA[i] * invB[j] * invC[l]
                cf = arb(0)
                cfp = arb(0)
                for k, v in cp.items():
                    vA = arb_of_fr(v)
                    cf += vA * pw[k]
                    if k:
                        cfp += vA * k * pw[k - 1]
                S += cf * Dm1
                Sp += cfp * Dm1
                Sc += cf * (-(invA[i] + invB[j] + invC[l] / 2)) * Dm1
            pre = cA ** 3 / 2
            dpre = 3 * cA ** 2 / 2
            Z = pre * S
            dZdp = dpre * cp_A * S + pre * (Sp + Sc * cp_A)
            dZdl = dpre * cl_A * S + pre * Sc * cl_A
            return Z, dZdp, dZdl
        finally:
            ctx.prec = old


class PosqKernelAdaptation:
    """eras-battery adaptation: posq C kernel as point engine, monotone
    c-bracket x ball-q sweep as the closed-box enclosure, exact closed-form
    jets as the derivative path. Vector parameters (p, lambda). See module
    docstring; the enclosure form is battery scaffolding, not a shipped
    production p-layer."""
    label = "posq-kernel-surrogate@2048/192"
    dim = 2

    def __init__(self, exact=None, prec_sweep=192, degraded_build=None):
        self.exact = exact or SurrogateExact()
        self.prec_sweep = prec_sweep
        self.degraded_build = degraded_build   # None | 53 (probe variant)
        self._enc_cache = {}
        self._rule_cache = {}
        self._jet_cache = {}
        self._tmp = tempfile.mkdtemp(prefix="posq_adapt_",
                                     dir=os.environ.get("TMPDIR", "/tmp"))
        self._ctr = 0

    # ---- certified rules at exact rational c (receipt-gated) ----
    def _rules_at(self, c):
        key = str(c)
        if key in self._rule_cache:
            return self._rule_cache[key]
        old_prec = ctx.prec
        old_build = bg.PREC_BUILD
        try:
            if self.degraded_build:
                bg.PREC_BUILD = self.degraded_build
            spec = [(17, c - 1, "a_u1"), (13, c - 1, "a_u2"),
                    (5, c / 2 - 1, "a_v3")]
            rules = {}
            for n, al, tag in spec:
                r = bg.build_rule(n, al, tag)
                rec = r["receipts"]
                ok = (rec["sum_w_contains_moment"]
                      and rec["all_weights_positive"]
                      and all(v["contains_exact"]
                              for v in rec["monomial_exactness"].values()))
                if not ok:
                    raise RuntimeError(f"rule receipts FAILED {tag} c={c}")
                rules[tag.split("_")[1]] = r
            self._rule_cache[key] = rules
            return rules
        finally:
            bg.PREC_BUILD = old_build
            ctx.prec = old_prec

    # ---- kernel point sweep at exact rational (p, lam) ----
    def _kernel_Z(self, p, lam, prec):
        c = c_of(p, lam)
        rules = self._rules_at(c)
        self._ctr += 1
        base = os.path.join(self._tmp, f"pt{self._ctr}")
        rp, jp, op = base + ".rules", base + ".job", base + ".out"
        old = ctx.prec
        try:
            ctx.prec = bg.PREC_BUILD + 64
            kio.write_rules(rp, {"u1": rules["u1"], "u2": rules["u2"],
                                 "v3": rules["v3"]})
            kio.write_job(jp, "cat", prec, p.numerator, p.denominator,
                          0, 17, rp, op, [SURR_VEC16])
            kio.run_kernel(jp, timeout=600)
            _, balls = kio.read_out(op)
            ctx.prec = prec + 64
            Z = arb_of_fr(c) ** 3 * balls[0] / 2
            return Z
        finally:
            ctx.prec = old
            for f in (rp, jp, op):
                if os.path.exists(f):
                    os.remove(f)

    # ---- python ball-q sweep of I(c, Q) (enclosure inner engine) ----
    def _sweep_ballq(self, rules, pball, prec):
        old = ctx.prec
        ctx.prec = prec
        try:
            pi1 = pball
            pi0 = 1 - pball
            x1, w1 = rules["u1"]["nodes"], rules["u1"]["weights"]
            x2, w2 = rules["u2"]["nodes"], rules["u2"]["weights"]
            x3, w3 = rules["v3"]["nodes"], rules["v3"]["weights"]
            pats = [x for x, _ in SURR_COUNTS]
            exps = [n for _, n in SURR_COUNTS]
            tot = arb(0)
            for i in range(len(x1)):
                row = arb(0)
                for j in range(len(x2)):
                    A, Bv = bg.ab_fiber(x1[i], x2[j], pi0, pi1, pats)
                    fib = arb(0)
                    for k in range(len(x3)):
                        Fv = None
                        for t in range(len(pats)):
                            Pt = A[t] + Bv[t] * x3[k]
                            Pn = Pt if exps[t] == 1 else Pt ** exps[t]
                            Fv = Pn if Fv is None else Fv * Pn
                        fib += w3[k] * Fv
                    row += w2[j] * fib
                tot += w1[i] * row
            return tot
        finally:
            ctx.prec = old

    # ---------------- consumer protocol ----------------
    def enclose(self, p0, radii, K):
        key = (tuple(str(t) for t in p0), tuple(str(t) for t in radii), 0)
        if key in self._enc_cache:
            return self._enc_cache[key]
        try:
            pc, lc = frs(p0[0]), frs(p0[1])
            rp, rl = frs(radii[0]), frs(radii[1])
            pl, ph = pc - rp, pc + rp
            ll, lh = lc - rl, lc + rl
            if not (0 < pl and ph < Fr(1, 2) and ll > 0):
                # outside the monotone-bracket domain: no claim
                return (arb("nan"), arb("nan"))
            clo, chi = c_of(pl, ll), c_of(ph, lh)
            Rlo = self._rules_at(clo)
            Rhi = self._rules_at(chi)
            old = ctx.prec
            try:
                ctx.prec = 256
                pball = arb_of_fr(pc) + arb(0, 1) * arb_of_fr(rp)
                Ihi = self._sweep_ballq(Rhi, pball, self.prec_sweep)
                Ilo = self._sweep_ballq(Rlo, pball, self.prec_sweep)
                ctx.prec = 256
                lo = bg.lb_of(arb_of_fr(clo) ** 3 * Ihi / 2)
                hi = bg.ub_of(arb_of_fr(chi) ** 3 * Ilo / 2)
                out = (lo, hi)
            finally:
                ctx.prec = old
        except RuntimeError:
            # rule build/receipts failed (e.g. degraded 53-bit variant):
            # no-information, never a silent pass
            out = (arb("nan"), arb("nan"))
        self._enc_cache[key] = out
        return out

    def point_eval(self, p, prec=None):
        pr = prec or 256
        pf, lf = fr_of_mid(p[0]), fr_of_mid(p[1])
        return self._kernel_Z(pf, lf, pr)

    def deriv(self, p, dim):
        pf, lf = fr_of_mid(p[0]), fr_of_mid(p[1])
        key = (str(pf), str(lf))
        if key not in self._jet_cache:
            self._jet_cache[key] = self.exact.jet_arb(pf, lf, prec=1280)
        _, dp, dl = self._jet_cache[key]
        return dp if dim == 0 else dl

    def degraded(self):
        """same adaptation, rule builder forced to 53-bit. Measured: the
        interval-Newton certification fails at node 0, so the degraded
        enclosure is non-finite (caught by the probe as no-information)."""
        return PosqKernelAdaptation(exact=self.exact,
                                    prec_sweep=self.prec_sweep,
                                    degraded_build=53)


def z_exact_sha256(exact=None):
    """sha256 of the exact surrogate value at (p0, lam0) as 'num/den' --
    must equal the sha256 of the pinned RUNG1 surrogate Z_exact_str."""
    ex = exact or SurrogateExact()
    Z = ex.Z_fr(P0, LAM0)
    s = f"{Z.numerator}/{Z.denominator}"
    return hashlib.sha256(s.encode()).hexdigest(), Z


def _selftest():
    """Quick live battery entry (python3 posq.py --selftest, seconds):

      S1  the exact closed-form surrogate value at (p0, lam0) reproduces the
          sha-pinned "num/den" string EXACTLY (pure exact arithmetic), and
          its ln matches the pinned float to 1e-9.
      S2  the C kernel sweeps the same point live and its certified ball
          CONTAINS the exact rational, beating the 2^-100 relwidth line.
          When no kernel binary runs and none can be built, S2 SKIPs with
          the named reason (kernel-absent is the designed cold-clone
          behavior; the full adversarial battery is verify_posq_adaptation).

    Exit 0 = S1 passed and S2 passed or named-SKIPped; 1 otherwise.
    """
    from flint import arb, ctx
    fails = []
    sha, Zx = z_exact_sha256()
    ok1 = sha == Z_EXACT_STR_SHA256_PIN
    old = ctx.prec
    ctx.prec = 320
    Zex = arb(Zx.numerator) / arb(Zx.denominator)
    lnZ = float(Zex.log().mid())
    ctx.prec = old
    ok1b = abs(lnZ - LN_Z_EXACT_PIN) < 1e-9
    print(f"[{'PASS' if ok1 else 'FAIL'}] S1 exact surrogate sha256 pin "
          f"({sha[:16]}..)")
    print(f"[{'PASS' if ok1b else 'FAIL'}] S1b ln Z = {lnZ:.14f} vs pinned "
          f"{LN_Z_EXACT_PIN}")
    if not (ok1 and ok1b):
        fails.append("S1")
    try:
        kio.kernel_path()
    except kio.KernelUnavailable as e:
        print(f"[SKIP] S2 kernel containment: {e}")
    else:
        adapt = PosqKernelAdaptation()
        Zk = adapt._kernel_Z(P0, LAM0, 256)
        ctx.prec = 320
        contains = bool(Zk.contains(Zex))
        relw2 = float((Zk.rad() / Zk.mid()).log() / arb(2).log())
        ctx.prec = old
        ok2 = contains and relw2 < -100
        print(f"[{'PASS' if ok2 else 'FAIL'}] S2 live kernel ball contains "
              f"exact rational (relwidth 2^{relw2:.1f}, line 2^-100)")
        if not ok2:
            fails.append("S2")
    print("posq selftest: " + ("PASS" if not fails else f"FAIL {fails}"))
    return 1 if fails else 0


if __name__ == "__main__" and "--selftest" in sys.argv:
    sys.exit(_selftest())

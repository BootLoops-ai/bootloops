#!/usr/bin/env python3
r"""calpha_rings.py — family / mass-tuple / ring table for the coalescer closure legs.

Shared by gate_compare.py (closed-form comparison, two-route gate) and
pslq_calpha.py (PSLQ closure through lockpick.pslq_gate).  Nothing numeric
is pasted here: every closed form is a STRUCTURE (sign, rational, pi power,
radical exponents, Gamma exponents) evaluated at run time at the caller's
dps, and every expected integer relation is DERIVED from that structure.

Family data (msq = squared masses, alpha = fractional threshold exponent):

    family  msq                 alphas     threshold (sum sqrt m)^2   sqrt(m_L)   Gamma directions
    K3      (1,1,1,9)           3/2        36 = 2^2 3^2               3           {}   (1/2 -> sqrt(pi), absorbed)
    CY3     (1,1,1,1,16)        5/4, 7/4   64 = 2^6                   4           {1/4}
    CY4     (1,1,1,1,1,25)      3/2        100 = 2^2 5^2              5           {}   (3/2 - {1,2,3} = half-integers)

Gamma directions = the fractional parts of alpha minus the integer threshold
exponents, reduced by the reflection formula to arguments in (0, 1/2):
Gamma(1-a) = pi / (sin(pi a) Gamma(a)), and Gamma(1/2) = sqrt(pi) lives in the
{pi} direction already.  This is the local Barnes-type statement: the
connection coefficient between the fractional branch
and the integer tower is built from Gamma of the indicial differences.

Rings (member NAMES follow the record's pslq_production.py so the CY3 fixture
reproduces PSLQ_PRODUCTION.json name for name):

    log (multiplicative form, target log|c|):
        natural  = ['1', 'log_pi', 'log2', 'log3', 'log5', 'log7'] + ['logG(a)' for a in Gamma directions]
        extended = natural + the control directions logG(1/4), logG(1/3) not already present
                   (a control direction must carry coefficient 0 in a genuine hit)
    additive (target c * pi^k, k = the family's pi_power):
        {1} x {sqrt(p) : p in radicals} x {G(a)^2, G(1-a)^2 : a in Gamma directions},
        radicals = squarefree parts of sqrt(m_L), of the threshold, and of the reflection
        sines 1/sin(pi a)^2 (a = 1/4 -> 2, a = 1/3 -> 3).

Recorded closed forms (the record's expressions; the display strings are
labels, the values come from the structures, and load-time selftest evaluates
the display expressions and asserts they agree with the structures):

    K3  c_{3/2} = -sqrt(3)/(36 pi)
    CY3 c_{5/4} = -1/(sqrt(2 pi) Gamma(1/4)^2) = -Gamma(3/4)^2/(2 sqrt(2) pi^{5/2})
    CY3 c_{7/4} = -5 Gamma(1/4)^2/(384 sqrt(2) pi^{5/2});  product 5/(768 pi^3)
    CY4 c_{3/2} = -sqrt(5)/(40 pi^2)   (recorded: 248 certified digits, PSLQ height 6 with 35 held-out
        digits, the K3 control in the same basket; |c|^2 = 1/(320 pi^4), so the additive pi power is 2 -- the
        K3-class prediction |c|^2 in Q pi^{-2}, i.e. pi power 1, did not hold)

API:
    family_spec(name) / spec_from_masses(msq, alphas)   -> dict (msq, alphas, threshold, primes, gamma_args, radicals, pi_power)
    ring_names(spec, form='log'|'additive', ring='natural'|'extended')  -> [names]
    member_value(name)                       -> mpf at the current mp.dps (name grammar below)
    members(names, digits)                   -> {name: decimal string of `digits` digits}
    closed_form(spec, alpha, source='auto'|dict) -> ClosedForm or None
    ClosedForm.value()                       -> mpf at the current dps
    ClosedForm.log_vector(names)             -> canonical integer tuple over [T]+names, or raises NotRingMember
    ClosedForm.display / .display_json
    product_identity(cf_a, cf_b)             -> ClosedForm (rational x pi^k) or None
    parse_expr(text)                         -> mpf (restricted expression grammar: numbers, pi, sqrt, gamma/Gamma, log, exp, + - * / ^ ( ))
    selftest()                               -> 0/1 (structures vs display expressions, derived vectors, ring audits)

Name grammar (member_value): '1'; 'pi'; 'log_pi'; 'log<p>' (p an integer); 'logG(a)' with a = p/q;
'sqrt<n>'; 'G(a)^k'; products joined by '*' of the multiplicative atoms; 'pi^k', '1/pi^k'.
"""
import hashlib, io, json, math, os, re, sys, tokenize
from fractions import Fraction
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
SMALL_PRIMES = (2, 3, 5, 7)             # the K3-class log ring (record: pslq_production.py log_basis_names)
CONTROL_GAMMA = (Fraction(1, 4), Fraction(1, 3))   # extended-ring control directions
SIN2_RADICAL = {Fraction(1, 2): 1, Fraction(1, 4): 2, Fraction(1, 3): 3, Fraction(1, 6): 1}
# 1/sin(pi a)^2 = 2 (a=1/4), 4/3 (a=1/3), 4 (a=1/6), 1 (a=1/2): squarefree radical part.

# 'display' = the record's printed label (gate.log / pslq_production.log); 'expr' = a parseable form of the
# same expression (both are evaluated by selftest() and asserted equal to the structure's value).
_FAMILIES = {
    'K3':  {'msq': (1, 1, 1, 9), 'alphas': ('3/2',), 'pi_power': '1',
            'closed': {'3/2': {'sign': -1, 'rational': '1/36', 'pi': '-1', 'radicals': {'3': '1/2'}, 'gamma': {},
                               'display': '-sqrt(3)/(36pi)', 'expr': '-sqrt(3)/(36*pi)'}}},
    'CY3': {'msq': (1, 1, 1, 1, 16), 'alphas': ('5/4', '7/4'), 'pi_power': '5/2',
            'closed': {'5/4': {'sign': -1, 'rational': '1', 'pi': '-1/2', 'radicals': {'2': '-1/2'}, 'gamma': {'1/4': -2},
                               'display': '-1/(sqrt(2*pi)*Gamma(1/4)^2) = -Gamma(3/4)^2/(2*sqrt(2)*pi^(5/2))',
                               'expr': '-1/(sqrt(2*pi)*Gamma(1/4)^2)'},
                       '7/4': {'sign': -1, 'rational': '5/384', 'pi': '-5/2', 'radicals': {'2': '-1/2'}, 'gamma': {'1/4': 2},
                               'display': '-5*Gamma(1/4)^2/(384*sqrt(2)*pi^(5/2))',
                               'expr': '-5*Gamma(1/4)^2/(384*sqrt(2)*pi^(5/2))'}}},
    'CY4': {'msq': (1, 1, 1, 1, 1, 25), 'alphas': ('3/2',), 'pi_power': '2',   # the additive-form power the record found (the K3-class prediction was 1)
            'closed': {'3/2': {'sign': -1, 'rational': '1/40', 'pi': '-2', 'radicals': {'5': '1/2'}, 'gamma': {},
                               'display': '-sqrt(5)/(40pi^2)', 'expr': '-sqrt(5)/(40*pi^2)'}}},
}


class NotRingMember(ValueError):
    """A closed form (or a closed-form spec) needs a direction the ring does not carry."""


class RingSpecError(ValueError):
    """Malformed family / masses / alpha / ring specification."""


# ---------------------------------------------------------------- arithmetic helpers
def factor_int(n):
    """Prime factorisation of a positive integer -> {p: e}."""
    n = int(n)
    if n <= 0: raise RingSpecError(f'factor_int({n})')
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out


def isqrt_exact(n):
    r = math.isqrt(int(n))
    return r if r * r == n else None


def squarefree_part(n):
    """n = s^2 * f with f squarefree -> f."""
    f = 1
    for p, e in factor_int(n).items():
        if e % 2: f *= p
    return f


def frac(x):
    if isinstance(x, Fraction): return x
    if isinstance(x, int): return Fraction(x)
    return Fraction(str(x).strip())


def frac_str(f):
    f = frac(f)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def reflect_gamma_arg(a):
    """Reduce a Gamma argument to (0, 1/2]: returns (a_reduced, flipped) — Gamma(1-a) = pi/(sin(pi a) Gamma(a))."""
    a = frac(a) % 1
    if a == 0: raise RingSpecError('Gamma at an integer argument')
    if a > Fraction(1, 2): return 1 - a, True
    return a, False


# ---------------------------------------------------------------- family / spec
def spec_from_masses(msq, alphas, pi_power=None, closed=None, family=None):
    """Build the ring spec from the squared-mass tuple and the fractional exponents (strings 'p/q')."""
    msq = tuple(int(m) for m in msq)
    roots = [isqrt_exact(m) for m in msq]
    if any(r is None for r in roots):
        raise RingSpecError(f'msq {msq}: every squared mass must be a perfect square for the (sum sqrt m)^2 threshold')
    alphas = tuple(frac_str(a) for a in alphas)
    for a in alphas:
        if frac(a).denominator == 1: raise RingSpecError(f'alpha {a} is an integer, not a fractional threshold exponent')
    thr = sum(roots) ** 2
    mL = roots[-1]
    # Gamma directions: alpha - (integer exponents 1..order) mod 1, reflection-reduced, Gamma(1/2) dropped
    gam = set()
    for a in alphas:
        for k in range(1, len(msq) + 1):
            d = (frac(a) - k) % 1
            if d == 0: continue
            r, _ = reflect_gamma_arg(d)
            if r == Fraction(1, 2): continue
            gam.add(r)
        for b in alphas:            # differences between fractional exponents (5/4 vs 7/4 -> 1/2)
            d = (frac(a) - frac(b)) % 1
            if d:
                r, _ = reflect_gamma_arg(d)
                if r != Fraction(1, 2): gam.add(r)
    gam = tuple(sorted(gam))
    rad = set()
    for n in (mL, thr):
        f = squarefree_part(n)
        if f > 1: rad.add(f)
    for a in gam:
        f = SIN2_RADICAL.get(a)
        if f is None: raise RingSpecError(f'Gamma direction {a}: reflection radical not tabulated (add to SIN2_RADICAL)')
        if f > 1: rad.add(f)
    rad = tuple(sorted(rad))
    if pi_power is None:
        pi_power = '1'
    return {'family': family, 'msq': msq, 'alphas': alphas, 'threshold': thr, 'sqrt_mL': mL,
            'threshold_primes': sorted(factor_int(thr)), 'gamma_args': gam, 'radicals': rad,
            'pi_power': frac_str(pi_power), 'closed': dict(closed or {})}


def family_spec(name):
    if name not in _FAMILIES:
        raise RingSpecError(f'unknown family {name!r}; known: {sorted(_FAMILIES)} (or pass --masses + --alphas)')
    f = _FAMILIES[name]
    return spec_from_masses(f['msq'], f['alphas'], pi_power=f['pi_power'], closed=f['closed'], family=name)


def families():
    return sorted(_FAMILIES)


# ---------------------------------------------------------------- rings
def ring_names(spec, form='log', ring='natural'):
    if form == 'log':
        names = ['1', 'log_pi'] + [f'log{p}' for p in SMALL_PRIMES] + [f'logG({frac_str(a)})' for a in spec['gamma_args']]
        if ring == 'extended':
            for a in CONTROL_GAMMA:
                n = f'logG({frac_str(a)})'
                if n not in names: names.append(n)
        elif ring != 'natural':
            raise RingSpecError(f'ring {ring!r}: natural | extended | <json path>')
        return names
    if form == 'additive':
        rads = [''] + [f'sqrt{p}' for p in spec['radicals']]
        gams = ['']
        for a in spec['gamma_args']:
            gams += [f'G({frac_str(a)})^2', f'G({frac_str(1 - a)})^2']
        names = []
        for g in gams:
            for r in rads:
                atoms = [x for x in (r, g) if x]
                names.append('*'.join(atoms) if atoms else '1')
        # record order: 1, sqrt2, G(1/4)^2, G(3/4)^2, sqrt2*G(1/4)^2, sqrt2*G(3/4)^2
        base = [n for n in names if n == '1' or n.startswith('sqrt') and '*' not in n]
        pure_g = [n for n in names if n.startswith('G(')]
        mixed = [n for n in names if n.startswith('sqrt') and '*' in n]
        out = base + pure_g + mixed
        if ring == 'extended':
            for p in SMALL_PRIMES:
                n = f'sqrt{p}'
                if n not in out and squarefree_part(p) == p: out.append(n)
        elif ring != 'natural':
            raise RingSpecError(f'ring {ring!r}: natural | extended | <json path>')
        return out
    raise RingSpecError(f'form {form!r}: log | additive')


def ring_from_json(path):
    """{'log': [names], 'additive': {'pi_power': 'k', 'members': [names]}} — names in the grammar above."""
    d = json.load(open(path))
    out = {}
    if 'log' in d: out['log'] = list(d['log'])
    if 'additive' in d: out['additive'] = {'pi_power': frac_str(d['additive'].get('pi_power', '1')),
                                           'members': list(d['additive']['members'])}
    if not out: raise RingSpecError(f'{path}: no "log" or "additive" ring')
    for form, names in out.items():
        for n in (names if form == 'log' else names['members']):
            member_value(n)   # grammar check at load
    return out


_ATOM = re.compile(r'^(1|pi|pi\^(-?\d+)|1/pi\^(\d+)|1/pi|log_pi|log(\d+)|logG\(([0-9/]+)\)|sqrt(\d+)|G\(([0-9/]+)\)(?:\^(-?\d+))?)$')


def member_value(name):
    """mpf value of a ring member NAME at the current mp.dps."""
    v = mp.mpf(1)
    for atom in name.split('*'):
        m = _ATOM.match(atom.strip())
        if not m: raise RingSpecError(f'member name {name!r}: atom {atom!r} not in the grammar')
        a = atom.strip()
        if a == '1': continue
        if a == 'pi': v *= mp.pi; continue
        if a == '1/pi': v /= mp.pi; continue
        if a == 'log_pi': v *= mp.log(mp.pi); continue
        if m.group(2): v *= mp.pi ** int(m.group(2)); continue
        if m.group(3): v /= mp.pi ** int(m.group(3)); continue
        if m.group(4): v *= mp.log(int(m.group(4))); continue
        if m.group(5): v *= mp.log(mp.gamma(_mpfrac(m.group(5)))); continue
        if m.group(6): v *= mp.sqrt(int(m.group(6))); continue
        if m.group(7):
            k = int(m.group(8)) if m.group(8) else 1
            v *= mp.gamma(_mpfrac(m.group(7))) ** k; continue
    return v


def _mpfrac(s):
    f = frac(s)
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)


def members(names, digits):
    """{name: decimal string carrying `digits` significant digits}, built at dps digits+20."""
    out = {}
    with mp.workdps(digits + 20):
        for n in names:
            out[n] = mp.nstr(member_value(n), digits, strip_zeros=False)
    return out


# ---------------------------------------------------------------- closed forms
class ClosedForm:
    """c = sign * rational * pi^e_pi * prod p^e_p (radicals, rational exponents) * prod Gamma(a)^g_a."""

    def __init__(self, sign, rational, pi, radicals, gamma, display=None, expr=None, label=None):
        self.sign = int(sign)
        if self.sign not in (-1, 1): raise RingSpecError('sign must be +1 or -1')
        self.rational = frac(rational)
        if self.rational <= 0: raise RingSpecError('rational must be positive (sign carries the sign)')
        self.pi = frac(pi)
        self.radicals = {int(p): frac(e) for p, e in (radicals or {}).items() if frac(e) != 0}
        self.gamma = {frac(a): int(g) for a, g in (gamma or {}).items() if int(g) != 0}
        self.display = display or self.render()
        self.expr = expr
        self.label = label

    # exact prime-exponent map of |c| without the Gamma part: {p: Fraction}, plus pi
    def prime_exponents(self):
        ex = {}
        for p, e in factor_int(self.rational.numerator).items(): ex[p] = ex.get(p, 0) + e
        for p, e in factor_int(self.rational.denominator).items(): ex[p] = ex.get(p, 0) - e
        for n, e in self.radicals.items():
            for p, k in factor_int(n).items(): ex[p] = ex.get(p, 0) + k * e
        return {p: frac(e) for p, e in ex.items() if e != 0}

    def gamma_reduced(self):
        """Gamma exponents reduced to arguments in (0,1/2]: returns ({a: g}, extra_pi, extra_primes)."""
        gam, extra_pi, extra_p = {}, Fraction(0), {}
        for a, g in self.gamma.items():
            r, flipped = reflect_gamma_arg(a)
            if not flipped:
                gam[r] = gam.get(r, 0) + g
            else:   # Gamma(1-r)^g = pi^g sin(pi r)^-g Gamma(r)^-g ; sin(pi r)^-2 = SIN2_RADICAL ratio
                gam[r] = gam.get(r, 0) - g
                extra_pi += g
                # sin(pi r)^{-g}: r=1/4 -> (sqrt2)^g = 2^{g/2}; r=1/3 -> (2/sqrt3)^g = 2^g 3^{-g/2}; r=1/6 -> 2^g; r=1/2 -> 1
                if r == Fraction(1, 4): extra_p[2] = extra_p.get(2, 0) + Fraction(g, 2)
                elif r == Fraction(1, 3): extra_p[2] = extra_p.get(2, 0) + g; extra_p[3] = extra_p.get(3, 0) - Fraction(g, 2)
                elif r == Fraction(1, 6): extra_p[2] = extra_p.get(2, 0) + g
                elif r == Fraction(1, 2): pass
                else: raise NotRingMember(f'Gamma({a}) reflection not tabulated')
        gam = {a: g for a, g in gam.items() if g}
        if Fraction(1, 2) in gam:      # Gamma(1/2) = sqrt(pi)
            extra_pi += Fraction(gam.pop(Fraction(1, 2)), 2)
        return gam, extra_pi, {p: e for p, e in extra_p.items() if e}

    def value(self):
        v = mp.mpf(self.rational.numerator) / mp.mpf(self.rational.denominator)
        v *= mp.pi ** _mpfrac(self.pi)
        for p, e in self.radicals.items(): v *= mp.mpf(p) ** _mpfrac(e)
        for a, g in self.gamma.items(): v *= mp.gamma(_mpfrac(a)) ** g
        return self.sign * v

    def log_vector(self, names):
        """Canonical integer relation over [log|c|] + names (names in the log grammar): D*T - D*sum = 0."""
        ex = self.prime_exponents()
        gam, extra_pi, extra_p = self.gamma_reduced()
        for p, e in extra_p.items(): ex[p] = ex.get(p, 0) + e
        ex = {p: e for p, e in ex.items() if e}
        e_pi = self.pi + extra_pi
        coeffs = {}
        have = set(names)
        for p, e in ex.items():
            n = f'log{p}'
            if n not in have: raise NotRingMember(f'closed form needs log{p} (exponent {e}); ring carries {names}')
            coeffs[n] = e
        if e_pi:
            if 'log_pi' not in have: raise NotRingMember('closed form needs log_pi; ring lacks it')
            coeffs['log_pi'] = e_pi
        for a, g in gam.items():
            n = f'logG({frac_str(a)})'
            if n not in have: raise NotRingMember(f'closed form needs {n} (exponent {g}); ring carries {names}')
            coeffs[n] = frac(g)
        D = 1
        for e in coeffs.values(): D = D * e.denominator // math.gcd(D, e.denominator)
        vec = [D] + [int(-D * coeffs.get(n, 0)) for n in names]
        return canonicalize(vec)

    def render(self):
        """Parseable expression of the structure (sign*rational*pi^e*p^e*Gamma(a)^g), for forms without a label."""
        parts = []
        if self.rational != 1: parts.append(f'({frac_str(self.rational)})')
        if self.pi: parts.append('pi' if self.pi == 1 else f'pi^({frac_str(self.pi)})')
        for p, e in sorted(self.radicals.items()): parts.append(f'{p}^({frac_str(e)})')
        for a, g in sorted(self.gamma.items()): parts.append(f'Gamma({frac_str(a)})' + ('' if g == 1 else f'^({g})'))
        s = '*'.join(parts) if parts else '1'
        return ('-' if self.sign < 0 else '') + s

    def render_rational_pi(self, style='json'):
        """Rational x pi^k rendering for the product identity ('json': 5/(768*pi^3); 'log': 5/(768 pi^3))."""
        if self.gamma or any(e.denominator != 1 for e in self.radicals.values()):
            raise RingSpecError('render_rational_pi() is for rational x pi^k forms only')
        num = self.rational.numerator
        den = self.rational.denominator
        for p, e in self.radicals.items():
            if e > 0: num *= p ** int(e)
            else: den *= p ** int(-e)
        s = '-' if self.sign < 0 else ''
        k = self.pi
        if k.denominator != 1: raise RingSpecError('render(): integer pi power only')
        k = int(k)
        sep = '*' if style == 'json' else ' '
        if k < 0:
            pip = f'pi^{-k}' if -k != 1 else 'pi'
            return f'{s}{num}/({den}{sep}{pip})' if den != 1 else f'{s}{num}/{pip}'
        if k > 0:
            pip = f'pi^{k}' if k != 1 else 'pi'
            return f'{s}{num}{sep}{pip}/{den}' if den != 1 else f'{s}{num}{sep}{pip}'
        return f'{s}{num}/{den}' if den != 1 else f'{s}{num}'

    def to_dict(self):
        return {'sign': self.sign, 'rational': frac_str(self.rational), 'pi': frac_str(self.pi),
                'radicals': {str(p): frac_str(e) for p, e in self.radicals.items()},
                'gamma': {frac_str(a): g for a, g in self.gamma.items()}, 'display': self.display}

    @classmethod
    def from_dict(cls, d, label=None):
        return cls(d.get('sign', 1), d.get('rational', '1'), d.get('pi', '0'), d.get('radicals', {}),
                   d.get('gamma', {}), display=d.get('display'), expr=d.get('expr'), label=label)


def closed_form(spec, alpha, source='auto'):
    """ClosedForm for alpha: 'auto' = the family's recorded form (None if unrecorded); 'none' = None;
    a dict {alpha: structure} or a path to such a JSON = the caller's form."""
    alpha = frac_str(alpha)
    if source == 'none' or source is None: return None
    if source == 'auto':
        d = spec['closed'].get(alpha)
        return ClosedForm.from_dict(d, label=alpha) if d else None
    if isinstance(source, str):
        source = json.load(open(source))
    d = source.get(alpha)
    return ClosedForm.from_dict(d, label=alpha) if d else None


def product_identity(cf_a, cf_b):
    """Product of two closed forms if it is rational x pi^k (Gamma and radicals cancel) — else None."""
    if cf_a is None or cf_b is None: return None
    gam = {}
    for cf in (cf_a, cf_b):
        for a, g in cf.gamma.items(): gam[a] = gam.get(a, 0) + g
    rad = {}
    for cf in (cf_a, cf_b):
        for p, e in cf.radicals.items(): rad[p] = rad.get(p, 0) + e
    # reflection pairs Gamma(a)Gamma(1-a) = pi/sin(pi a) collapse to pi x radical
    extra_pi, extra_rad = Fraction(0), {}
    for a in list(gam):
        b = 1 - a
        if a < b and b in gam and gam[a] == gam[b]:
            g = gam.pop(a); gam.pop(b)
            extra_pi += g
            f = SIN2_RADICAL.get(a)
            if a == Fraction(1, 4): extra_rad[2] = extra_rad.get(2, 0) + Fraction(g, 2)
            elif a == Fraction(1, 3): extra_rad[2] = extra_rad.get(2, 0) + g; extra_rad[3] = extra_rad.get(3, 0) - Fraction(g, 2)
            elif a == Fraction(1, 6): extra_rad[2] = extra_rad.get(2, 0) + g
    gam = {a: g for a, g in gam.items() if g}
    if gam: return None
    for p, e in extra_rad.items(): rad[p] = rad.get(p, 0) + e
    rad = {p: e for p, e in rad.items() if e}
    if any(frac(e).denominator != 1 for e in rad.values()): return None
    rational = cf_a.rational * cf_b.rational
    for p, e in rad.items():
        rational *= Fraction(p) ** int(e)
    return ClosedForm(cf_a.sign * cf_b.sign, rational, cf_a.pi + cf_b.pi + extra_pi, {}, {})


# ---------------------------------------------------------------- expression grammar (display strings, controls)
_ALLOWED = re.compile(r'^[0-9a-zA-Z_+\-*/^().,\s]+$')
_NAMES = {'pi': lambda: mp.pi, 'e': lambda: mp.e}
_FUNCS = {'sqrt': mp.sqrt, 'gamma': mp.gamma, 'Gamma': mp.gamma, 'log': mp.log, 'ln': mp.log, 'exp': mp.exp}


def parse_expr(text):
    """Evaluate a restricted closed-form expression at the current mp.dps (integers become exact mpf;
    1/3 is mpf(1)/mpf(3), never a float)."""
    if not _ALLOWED.match(text): raise RingSpecError(f'expression {text!r}: characters outside the grammar')
    out = []
    for t in tokenize.generate_tokens(io.StringIO(text.replace('^', '**')).readline):
        if t.type == tokenize.NUMBER: out.append(f"mpf('{t.string}')")
        elif t.type == tokenize.NAME:
            if t.string in _NAMES or t.string in _FUNCS or t.string == 'mpf': out.append(t.string)
            else: raise RingSpecError(f'expression {text!r}: name {t.string!r} not allowed')
        elif t.type == tokenize.OP: out.append(t.string)
        elif t.type in (tokenize.NEWLINE, tokenize.NL, tokenize.ENDMARKER, tokenize.INDENT, tokenize.DEDENT): pass
        else: raise RingSpecError(f'expression {text!r}: token {t.string!r} not allowed')
    ns = {'__builtins__': {}, 'mpf': mp.mpf}
    ns.update({k: v() for k, v in _NAMES.items()}); ns.update(_FUNCS)
    try:
        return mp.mpf(eval(' '.join(out), ns))
    except (SyntaxError, TypeError, NameError) as e:
        raise RingSpecError(f'expression {text!r}: {type(e).__name__}: {e}')


def canonicalize(vec):
    """gcd-divide + first nonzero positive (the lockpick.pslq_gate rule, re-stated here so the ring module
    needs no engine import; pslq_calpha asserts equality with the engine's canonicalize on every run)."""
    if not vec: return None
    g = 0
    for v in vec: g = math.gcd(g, abs(int(v)))
    if g == 0: return None
    out = [int(v) // g for v in vec]
    for v in out:
        if v:
            if v < 0: out = [-q for q in out]
            break
    return tuple(out)


def sha16(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()[:16]


# ---------------------------------------------------------------- selftest
def selftest(dps=120, verbose=True):
    fails, total = 0, 0

    def rep(name, ok):
        nonlocal fails, total
        total += 1; fails += (not ok)
        if verbose: print(('PASS' if ok else 'FAIL'), name, flush=True)
    with mp.workdps(dps):
        for fam in families():
            spec = family_spec(fam)
            for alpha, d in spec['closed'].items():
                cf = closed_form(spec, alpha)
                v = cf.value()
                forms = [d['expr']] + [s.strip() for s in d['display'].split(' = ') if s.strip() != d['expr']]
                for e in forms:
                    try:
                        w = parse_expr(e)
                    except RingSpecError:
                        if e == d['expr']:
                            rep(f'{fam} c_{{{alpha}}} expr {e!r} parses', False)
                        continue      # a display label that is not in the grammar (e.g. "36pi") is a label only
                    ok = abs(w - v) <= abs(v) * mp.mpf(10) ** (-(dps - 10))
                    rep(f'{fam} c_{{{alpha}}} structure == expression {e!r} to {dps-10} d', ok)
                # the derived log relation must annihilate log|c| numerically over the natural log ring
                names = ring_names(spec, 'log', 'natural')
                vec = cf.log_vector(names)
                x = [mp.log(abs(v))] + [member_value(n) for n in names]
                res = abs(mp.fsum(c * y for c, y in zip(vec, x)))
                rep(f'{fam} c_{{{alpha}}} derived log vector {vec} over {names} annihilates log|c| ({mp.nstr(res, 3)})',
                    res < mp.mpf(10) ** (-(dps - 10)))
        # ring facts (the table in the docstring), computed
        k3, cy3, cy4 = family_spec('K3'), family_spec('CY3'), family_spec('CY4')
        rep('K3 threshold 36, sqrt_mL 3, no Gamma direction, radicals (3,)',
            (k3['threshold'], k3['sqrt_mL'], k3['gamma_args'], k3['radicals']) == (36, 3, (), (3,)))
        rep('CY3 threshold 64, sqrt_mL 4, Gamma direction (1/4,), radicals (2,)',
            (cy3['threshold'], cy3['sqrt_mL'], cy3['gamma_args'], cy3['radicals']) == (64, 4, (Fraction(1, 4),), (2,)))
        rep('CY4 threshold 100, sqrt_mL 5, no Gamma direction, radicals (5,)',
            (cy4['threshold'], cy4['sqrt_mL'], cy4['gamma_args'], cy4['radicals']) == (100, 5, (), (5,)))
        rep("CY3 natural log ring = ['1','log_pi','log2','log3','log5','log7','logG(1/4)']",
            ring_names(cy3, 'log') == ['1', 'log_pi', 'log2', 'log3', 'log5', 'log7', 'logG(1/4)'])
        rep("CY3 natural additive ring = ['1','sqrt2','G(1/4)^2','G(3/4)^2','sqrt2*G(1/4)^2','sqrt2*G(3/4)^2']",
            ring_names(cy3, 'additive') == ['1', 'sqrt2', 'G(1/4)^2', 'G(3/4)^2', 'sqrt2*G(1/4)^2', 'sqrt2*G(3/4)^2'])
        rep("K3 natural additive ring = ['1','sqrt3']; CY4 = ['1','sqrt5']",
            ring_names(k3, 'additive') == ['1', 'sqrt3'] and ring_names(cy4, 'additive') == ['1', 'sqrt5'])
        rep('CY4 extended log ring adds the logG(1/4) and logG(1/3) control directions',
            ring_names(cy4, 'log', 'extended') == ring_names(cy4, 'log') + ['logG(1/4)', 'logG(1/3)'])
        prod = product_identity(closed_form(cy3, '5/4'), closed_form(cy3, '7/4'))
        rep("CY3 product identity derived: 5/(768*pi^3)", prod is not None and prod.render_rational_pi('json') == '5/(768*pi^3)'
            and abs(prod.value() - closed_form(cy3, '5/4').value() * closed_form(cy3, '7/4').value())
            < mp.mpf(10) ** (-(dps - 10)))
        # a closed form outside the ring is refused by name
        try:
            ClosedForm(-1, '1/11', '-1', {}, {}).log_vector(ring_names(cy3, 'log'))
            rep('closed form with prime 11 outside the log ring -> NotRingMember', False)
        except NotRingMember:
            rep('closed form with prime 11 outside the log ring -> NotRingMember', True)
        try:
            ClosedForm(-1, '1', '-1', {}, {'1/5': 2}).log_vector(ring_names(cy3, 'log'))
            rep('closed form with Gamma(1/5) outside the ring -> NotRingMember', False)
        except NotRingMember:
            rep('closed form with Gamma(1/5) outside the ring -> NotRingMember', True)
        # reflection: the second display form of c_{5/4} as a structure gives the SAME log vector
        alt = ClosedForm(-1, '1/2', '-5/2', {'2': '-1/2'}, {'3/4': 2})
        rep('reflection: -Gamma(3/4)^2/(2 sqrt2 pi^{5/2}) reduces to the same log vector as -1/(sqrt(2pi)Gamma(1/4)^2)',
            alt.log_vector(ring_names(cy3, 'log')) == closed_form(cy3, '5/4').log_vector(ring_names(cy3, 'log')))
        rep('parse_expr keeps rationals exact (gamma(1/3) at dps: matches mp.gamma(mpf(1)/3))',
            abs(parse_expr('gamma(1/3)') - mp.gamma(mp.mpf(1) / 3)) < mp.mpf(10) ** (-(dps - 5)))
        for fam in families():
            spec = family_spec(fam)
            for alpha in spec['closed']:
                cf = closed_form(spec, alpha)
                rep(f'{fam} c_{{{alpha}}} render() {cf.render()!r} parses back to the structure value',
                    abs(parse_expr(cf.render()) - cf.value()) <= abs(cf.value()) * mp.mpf(10) ** (-(dps - 10)))
    print('CALPHA_RINGS SELFTEST', 'FAIL' if fails else 'PASS', f'({total - fails}/{total}, 0 skipped)', flush=True)
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(selftest())

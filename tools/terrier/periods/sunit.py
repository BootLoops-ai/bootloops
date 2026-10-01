#!/usr/bin/env python3
"""sunit.py (periods wing) — refined S-integrality candidate
generator (census version; generator code in
pipeline/sunit/{gen_sunit,gen_refined}.py). Facade only routes +
packages; no law duplicated.

LAW (refined-conductor-boxS): the fibre at z
has small conductor <=> z is S-integral on P1 minus the singular divisor
— for ONE S inside the FROZEN box {2,3,5,7,11,13}, ALL of lead(f_z),
const(f_z) and N(m_i(z)) for EVERY singular factor m_i are S-supported
(distances-to-singular-values law). ACCEPTANCE is family-free and
box-exact: S_min(z) = union of the supports, accept iff S_min <= box.
Generation nets are BOUNDED (rational exponent/height/support caps;
quadratic-field Q(sqrt D) group enumeration with pool caps) — every run
emits the exponent-bound truncation statement as MACHINE-READABLE
output (truncation_statement(); a bounded net, NOT a finiteness proof).

Surfaces:
  configure_net / truncation_statement   net caps as data, in + out
  load_cards / sing_factors              card atlas + singular divisor
  smin_rational / smin_quadratic         family-free acceptance -> S_min
  fund_unit / quad_cell / q_* helpers    Q(sqrt D) enumeration machinery
  generate_op                            one op -> candidates JSONL +
                                         stats (known-orbit strips via
                                         packaged common/data orbit)
  known_point_gate                       falsifiable known-point gate +
                                         SIZES receipt (its default gate
                                         ops reference an AESZ operator
                                         atlas not included in the
                                         package — point
                                         TERRIER_SUNIT_ATLAS at your own
                                         atlas to use it)
Env overrides: TERRIER_SUNIT_ATLAS, TERRIER_KNOWN_ORBIT,
TERRIER_SUNIT_OUT (set per-call by generate_op/known_point_gate).
Battery: selftest_sunit.py — planted-point recovery with hand-derived
S_min on the shipped synthetic control card, plus must-fail controls."""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUNIT_DIR = os.path.join(HERE, "pipeline", "sunit")
if SUNIT_DIR not in sys.path:
    sys.path.insert(0, SUNIT_DIR)

CENSUS_NET = dict(A_MAX=12, H_RAT=40.0, RAT_SUPP_MAX=4, H_POOL=20.0,
                  POOL_CAP=600, H_QUAD=40.0, SQFREE_D="squarefree 2..100")


def _mods():
    import gen_sunit as GS
    import gen_refined as GR
    return GS, GR


def truncation_statement():
    """The exponent-bound truncation statement of the CURRENT net, as
    machine-readable data (census SIZES2 meta law + completeness text)."""
    GS, GR = _mods()
    return {
        "law": "refined-conductor-boxS (family-free box-exact acceptance)",
        "box": sorted(GR.BOX),
        "net": dict(A_MAX=GR.A_MAX, H_RAT=GR.H_RAT,
                    RAT_SUPP_MAX=GR.RAT_SUPP_MAX, H_POOL=GR.H_POOL,
                    POOL_CAP=GR.POOL_CAP, H_QUAD=GR.H_QUAD,
                    SQFREE_D=list(GS.SQFREE_D),
                    quad_gen_families="S in box, |S|<=3 (net only)"),
        "completeness": (
            "rational tier complete for z = +-prod p^a_p, supp(z) <= "
            f"{GR.RAT_SUPP_MAX} box primes, 1 <= |a_p| <= {GR.A_MAX}, "
            f"Weil height <= {GR.H_RAT}; quadratic tier = bounded "
            "group-generated net per Q(sqrt D). BOUNDED NET, NOT A "
            "FINITENESS PROOF; acceptance itself is exact."),
        "blowup_alarm_threshold": GR.BLOWUP}


def configure_net(**kw):
    """Override net caps (A_MAX, H_RAT, RAT_SUPP_MAX, H_POOL, POOL_CAP,
    H_QUAD on gen_refined; SQFREE_D list on gen_sunit). Acceptance is
    untouched — only GENERATION nets move. Returns the machine-readable
    truncation statement of the resulting net."""
    GS, GR = _mods()
    for k, v in kw.items():
        if k == "SQFREE_D":
            GS.SQFREE_D = list(v)
        elif hasattr(GR, k):
            setattr(GR, k, v)
        else:
            raise KeyError(f"unknown net parameter {k!r}")
    return truncation_statement()


def load_cards():
    """{op: card} + alias map from the packaged atlas (env-overridable
    via TERRIER_SUNIT_ATLAS)."""
    GS, _GR = _mods()
    return GS.load_cards()


def sing_factors(card):
    """Primitive ascending integer factors of the singular divisor
    (carded singular_locus UNION leadpoly factors; origin dropped)."""
    _GS, GR = _mods()
    return GR.sing_factors(card)


def smin_rational(factors, num, den, sgn=1):
    """Family-free acceptance for z = sgn*num/den -> S_min set or None
    (None = out of box, or z IS a singular point)."""
    _GS, GR = _mods()
    return GR.check_rational(num, den, sgn, factors)


def smin_quadratic(triple, D, factors):
    """Family-free acceptance for quadratic z = (a + b sqrt D)/c ->
    (S_min, height, is_unit, minpoly_key) or None."""
    _GS, GR = _mods()
    return GR.check_quad(triple, D, math.sqrt(D), factors)


def fund_unit(D):
    GS, _GR = _mods()
    return GS.fund_unit(D)


def quad_cell(S, D, cache):
    """Bounded S-unit-group enumeration cell for (S, Q(sqrt D)) —
    gen_sunit best-first law; caps read from gen_sunit globals."""
    GS, _GR = _mods()
    return GS.quad_cell(S, D, cache)


def generate_op(op, outdir, cache=None):
    """Run the refined generator for one op into outdir/<op>.jsonl
    (rows carry S_min, height, unit flag, KNOWN-ORBIT-REDISCOVERY
    dedup tags from the packaged orbit). Returns (stats, cache);
    stats includes blowup_alarm (per-op list alarm law)."""
    GS, GR = _mods()
    os.makedirs(outdir, exist_ok=True)
    os.environ["TERRIER_SUNIT_OUT"] = outdir
    cards, _al = GS.load_cards()
    if cache is None:
        cache = {"known": GS.load_known()}
    return GR.assemble_op(op, cards[op], cache), cache


def known_point_gate(outdir, ops=None):
    """The falsifiable known-point gate: generate all gate ops, then
    check the 8 GATE_POINTS (writes the gate rows and a SIZES receipt
    with the truncation statement into outdir). Returns
    (all_pass, gate_rows, sizes)."""
    _GS, GR = _mods()
    os.makedirs(outdir, exist_ok=True)
    sizes, cache = {}, None
    for op in (ops or GR.GATE_OPS):
        sizes[op], cache = generate_op(op, outdir, cache)
    os.environ["TERRIER_SUNIT_OUT"] = outdir
    ok = GR.run_gate()
    rows = [json.loads(l) for l in
            open(os.path.join(outdir, "GATE_REFINED.jsonl"))]
    receipt = {"law": "refined-conductor-boxS",
               "truncation_statement": truncation_statement(),
               "gate": {r["gate"]: r["status"] for r in rows},
               "ops": sizes}
    json.dump(receipt, open(os.path.join(outdir, "SIZES_gate.json"),
                            "w"), indent=1)
    return ok, rows, sizes

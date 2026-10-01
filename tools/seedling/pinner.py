"""Pinned jobs.yaml emission.

Emits the exact kira jobs.yaml shape used by the deployed pin recipe
(golden-file validated byte-exact against four stored production staging
cells).

s-FLOOR RULE (measured exhibit): support-only pinning of s is NOT safe. A
stage-2 pin at (r:12, s:0, d:5) produced a system that SOLVED CLEANLY but gave
wrong values (dropped super-sector-chain equations => false-master
termination), caught only by the two-slice value gate. The reference s:1 pin
on the same family passed 442/442. Rule: the pinned s is floored at the super-sector chain
demand, derived conservatively from the target set:
  - any target with ISP degree >= 1            -> chain floor >= 1
  - any target strictly below a top sector     -> chain floor >= 1
    (its reduction traverses super-sector chains, which carry ISP-bearing seeds)
  - only if ALL targets are top-sector and s=0 -> floor 0 permitted
Callers may force s_floor=0 only explicitly (UNSAFE; must be logged); the runner's
two-slice differential gate remains mandatory regardless.

Margins: (r,s,d) offsets added on top of the floored support. s-margins should only
ever come from the escalation ladder, last (s is the explosive dimension: x5.78
family bulk per unit at the reference pin point).

Cache note: any change to jobs.yaml changes kira ibp-cache keys; always emit into a
fresh run directory.
"""

_TEMPLATE = """\
jobs:
 - reduce_sectors:
    reduce:
     - {{topologies: [{fam}], sectors: [{sectors}], r: {r}, s: {s}, d: {d}}}
{trunc}    select_integrals:
      select_mandatory_list:
        - [{fam}, {target_name}]
{preferred_line}    integral_ordering: {ordering}
    run_initiate: true
    run_firefly: {firefly}
"""

_TRUNC = """\
    truncate_sp:
     - {{topologies: [{fam}], l: {l}}}
"""


def chain_s_floor(targets, top_sectors):
    """Minimum safe s given super-sector chain demand (conservative).

    targets: iterable of support.Target; top_sectors: iterable of sector ids.
    Returns 0 only when every target sits IN a top sector with no ISP powers.
    """
    targets = list(targets)
    tops = set(top_sectors)
    floor = 0
    for tg in targets:
        if tg.s >= 1:
            floor = max(floor, tg.s, 1)
        elif tg.sector not in tops:
            floor = max(floor, 1)
    return floor


def pinned_cell(support_g, targets, top_sectors, margins=(0, 0, 0), s_floor=None):
    """Compute the pinned (r,s,d). s_floor=None -> auto (chain_s_floor);
    explicit s_floor=0 is honored but UNSAFE (caller must log it)."""
    if s_floor is None:
        s_floor = chain_s_floor(targets, top_sectors)
    r = support_g["r"] + margins[0]
    s = max(support_g["s"], s_floor) + margins[1]
    d = support_g["d"] + margins[2]
    return r, s, d


def render_jobs_yaml(family, top_sectors, r, s, d, *, truncate_sp_l=None,
                     target_name="target", preferred_name="preferred",
                     integral_ordering=5, run_firefly=False):
    trunc = _TRUNC.format(fam=family, l=truncate_sp_l) if truncate_sp_l else ""
    return _TEMPLATE.format(
        fam=family,
        sectors=", ".join(str(x) for x in top_sectors),
        r=r, s=s, d=d, trunc=trunc,
        target_name=target_name,
        preferred_line=(f"    preferred_masters: {preferred_name}\n"
                        if preferred_name else ""),
        ordering=integral_ordering,
        firefly="true" if run_firefly else "false",
    )


def emit_jobs_yaml(out_path, family, top_sectors, support_g, targets=None,
                   margins=(0, 0, 0), s_floor=None, **kw):
    """support_g: global_support() dict; targets: list for the s-floor derivation
    (required unless s_floor given explicitly)."""
    if s_floor is None and targets is None:
        raise ValueError("need targets for chain_s_floor, or an explicit s_floor")
    r, s, d = pinned_cell(support_g, targets or [], top_sectors, margins, s_floor)
    text = render_jobs_yaml(family, top_sectors, r, s, d, **kw)
    with open(out_path, "w") as fh:
        fh.write(text)
    return text


# --- per-sector schedules ------------------------------------------------------
# Only 21% of the global pin's bulk sits in target sectors on the reference
# family; a schedule gives target sectors their own support and closure sectors a
# chain cell. Expressible in stock kira: multiple reduce entries. GATING (chain-
# floor lesson) is mandatory before production use: closure demand is set by what the
# reduction maps down, not by the target list.

_SCHED_TEMPLATE = """\
jobs:
 - reduce_sectors:
    reduce:
{entries}{trunc}    select_integrals:
      select_mandatory_list:
        - [{fam}, {target_name}]
{preferred_line}    integral_ordering: {ordering}
    run_initiate: true
    run_firefly: {firefly}
"""


def build_schedule(sector_support, census, closure_cell, s_floor=1):
    """Build a per-sector schedule.

    sector_support: support.support_table() output (target sectors).
    census: dict sector_id -> t (ALL nontrivial sectors to seed).
    closure_cell: callable t -> (r, s, d) for non-target sectors.
    s_floor: chain floor applied to target sectors (default 1).
    Returns dict cell(r,s,d) -> sorted [sector ids].
    """
    groups = {}
    for sid, t in census.items():
        if sid in sector_support:
            e = sector_support[sid]
            cell = (max(e["r"], t), max(e["s"], s_floor), e["d"])
        else:
            cell = tuple(closure_cell(t))
        groups.setdefault(cell, []).append(sid)
    return {c: sorted(v) for c, v in groups.items()}


def render_schedule_jobs_yaml(family, groups, *, truncate_sp_l=None,
                              target_name="target", preferred_name="preferred",
                              integral_ordering=5, run_firefly=False):
    """groups: dict (r,s,d) -> [sector ids] (build_schedule output)."""
    lines = []
    for (r, s, d) in sorted(groups, reverse=True):
        secs = ", ".join(str(x) for x in groups[(r, s, d)])
        lines.append(f"     - {{topologies: [{family}], sectors: [{secs}], "
                     f"r: {r}, s: {s}, d: {d}}}\n")
    trunc = _TRUNC.format(fam=family, l=truncate_sp_l) if truncate_sp_l else ""
    return _SCHED_TEMPLATE.format(
        fam=family, entries="".join(lines), trunc=trunc,
        target_name=target_name,
        preferred_line=(f"    preferred_masters: {preferred_name}\n"
                        if preferred_name else ""),
        ordering=integral_ordering,
        firefly="true" if run_firefly else "false",
    )

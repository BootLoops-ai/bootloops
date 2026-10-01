"""canonical_form — UT-basis / ε-factorized-DE builder.

The minimum-viable surface (no external engine needed):

    from canonical_form import Family, leading_singularities, ut_rotation
    rot = ut_rotation(family_spec, masters, point={"s":-3,"t":-7})
    heldout_certify(values, kin, ..., ut_rotation=rot)

Full ε-factorization (wraps the Counterweight engine one directory up;
override with env COUNTERWEIGHT_ROOT):

    from canonical_form import factor_epsilon, fuchsify, to_dlog
"""
from .canonical_form import (                                # noqa: F401
    Family, LSResult,
    leading_singularities, ut_rotation, ut_rotation_matrix,
    elliptic_period, verify_unit_pole,
    fuchsify, factor_epsilon, to_dlog, wrap_epsfactor,
    eps,
)

__all__ = [
    "Family", "LSResult",
    "leading_singularities", "ut_rotation", "ut_rotation_matrix",
    "elliptic_period", "verify_unit_pole",
    "fuchsify", "factor_epsilon", "to_dlog", "wrap_epsfactor",
    "eps",
]

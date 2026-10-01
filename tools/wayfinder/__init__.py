#!/usr/bin/env python3
"""
wayfinder — saved-DE loading, fixed-eps endpoint transport, Frobenius
landing (matrix, with integer-gap shearing; scalar resonant towers in
frob_scalar; var=infinity branch classification via the boundary_branches
delegator; two-branch apparent-pole series in two_sector), eps-fan Laurent
extraction, AMFlow-convention boundary seeds, digit gates, input manifests,
the refine-until-gate adaptive quadrature (quad.quad_refine),
and the conic third-kind reduction core
(conic_thirdkind — exact extended
t-system builder over third-kind atoms J=INT x^k/(f*Y) on a CONIC
Y^2 = q(x,t) [deg-2 scope NAMED, fail-closed], certified Hermite/partial-
fraction engines, ExtSys duck-typed DESystem for the moving-pole march),
the engine-agnostic eps->0 limit extractor (epslimit — Laurent/Richardson
extrapolation of an EXTERNAL engine's (eps_i, value_i) grid; the sibling of
epsfan, which owns its evaluator and picks its own nodes) and the FLINT
exact rational-function exporter (flintexport — writes the monomial-json
connections de_load reads; needs python-flint, lazy).
One shared evaluator contract:

    DESystem: .n, .var, .A(x, eps, dps) -> n x n list[list[mpc]], .meta
    optional fast-path protocol (duck-typed by transport_fixed_eps):
        .A_series(x, eps, dps, M), .singular_points
    — attach via DESystem.enable_fast_path(eps, wp) (ratfun RF layer)
    or provide your own.

See README.md for the API summary
and the design rules (no Laurent window,
R=1e-3 Cauchy default, real-eps Frobenius matching only, explicit dps
everywhere — global mp.dps is never touched by any module).

Import styles supported:
    from wayfinder import transport_fixed_eps            # package
    sys.path.insert(0, <this dir>); import transport       # flat (tests)
"""

from .de_load import (DESystem, load_monomial_json, load_amatrix_json,
                      load_kira_targets)
from .transport import transport_fixed_eps, HAVE_FLINT
from .frobenius import frobenius_basis, land
from . import frob_scalar          # scalar resonant-tower transport (module)
from . import boundary_branches    # var=inf delegator (lazy, import-guarded)
from .epsfan import (cauchy_laurent, vandermonde_laurent, anchor_roundtrip,
                     cauchy_laurent_certified, vandermonde_laurent_certified,
                     eps_analyticity_radius, desystem_eps_polys,
                     EpsFanCertifyError)
from .boundary import (tadpole, vacuum_known, singlemass_vacuum,
                       single_mass_prefactor, vacuum_ending_seed)
from .gate import gate_strings, write_gate_report, feed_gate_table
from .manifest import sha256_file, write_manifest, check_manifest
from .ratfun import RF
from .two_sector import two_sector_series, eval_two_sector
from .quad import quad_refine, QuadCert, QuadNonConvergence
from .acb_fast_attach import attach_fast_path, attach_acb_fast

__version__ = "1.0.0"


def __getattr__(name):
    # conic_thirdkind is LAZY (frob_scalar discipline: package import stays
    # sympy-free; the member imports sympy at ITS module load). NB: must use
    # importlib here — `from . import X` inside __getattr__ re-enters this
    # hook via _handle_fromlist's hasattr probe = infinite recursion.
    # epslimit (pure mpmath) and flintexport (python-flint + sympy) follow the
    # same lazy discipline: package import never pulls python-flint or sympy.
    if name in ("conic_thirdkind", "epslimit", "flintexport"):
        import importlib
        return importlib.import_module("." + name, __name__)
    raise AttributeError(f"module 'wayfinder' has no attribute {name!r}")

__all__ = [
    # de_load — saved DE systems -> one DESystem contract
    "DESystem",
    "load_monomial_json",
    "load_amatrix_json",
    "load_kira_targets",
    # transport — fixed-eps Taylor endpoint transport
    "transport_fixed_eps",
    # frobenius — regular-singular local basis + landing (REAL eps0 only;
    # nonzero-integer indicial gaps handled by Moser/Turrittin shearing)
    "frobenius_basis",
    "land",
    # frob_scalar — SCALAR operators / P-recurrences with resonant integer
    # towers, value-at-singularity transport (use frob_scalar.transport_value)
    "frob_scalar",
    # boundary_branches — thin import-guarded delegator to
    # tools/frobenius_boundary (branch classification at var=infinity;
    # eig(D1) caveat flagged in its docstring)
    "boundary_branches",
    # epsfan — Cauchy Laurent fan over eps (R=1e-3, 32 nodes, measured diag;
    # schwarz_real mirror is OPT-IN, Schwarz-real rows only), real-grid
    # Vandermonde mode (production Laurent-gate node design), anchor round-trip
    "cauchy_laurent",
    "vandermonde_laurent",
    "anchor_roundtrip",
    # epsfan certified layer — runtime R0 from exact eps-singularity data (REFUSES without
    # it), control-circle Cauchy alias bound with (r,M) escalation to caps,
    # cross-radius belt gate, two-grid + cond-monitored Vandermonde,
    # conditioning-derived working precision; fail-closed
    # EpsFanCertifyError naming (bound, tol, r, M, R0)
    "cauchy_laurent_certified",
    "vandermonde_laurent_certified",
    "eps_analyticity_radius",
    "desystem_eps_polys",
    "EpsFanCertifyError",
    # boundary — AMFlow-convention vacuum/tadpole boundary values
    "tadpole",
    "vacuum_known",
    "singlemass_vacuum",
    "single_mass_prefactor",
    "vacuum_ending_seed",
    # gate — two-precision digit gates against vendored strings + the
    # per-literal FEED gate table
    "gate_strings",
    "write_gate_report",
    "feed_gate_table",
    # manifest — sha256 input-artifact manifests
    "sha256_file",
    "write_manifest",
    "check_manifest",
    # ratfun — fixed-eps rational-function layer
    # (backs DESystem.enable_fast_path's A_series/singular_points protocol)
    "RF",
    # two_sector — two-branch (s^n, s^{n-eps}) Frobenius series with
    # apparent s^-k poles + slot scheduling (single-pin autopin)
    "two_sector_series",
    "eval_two_sector",
    # quad — adaptive-quadrature refine-until-gate (double-refinement
    # agreement gated at 10^-(dps+guard), depth seeds scale with dps,
    # fail-closed QuadNonConvergence with named diagnostics)
    "quad_refine",
    "QuadCert",
    "QuadNonConvergence",
    # acb_fast_attach — attacher pair for scalar-operator companion systems
    # with EXACT rational coefficients (the scalar-companion route de_load
    # does not build): duck-typed mpmath A_series fast path + flint acb
    # step-kernel table enabling backend="acb"; fixed at eps=0
    "attach_fast_path",
    "attach_acb_fast",
    # conic_thirdkind — the conic third-kind reduction core (promotion,
    # holding-class consolidation item 3): certified conic Hermite +
    # squarefree partial fractions + third-kind-atom connection builder +
    # ExtSys DESystem wrapper. SCOPE NAMED fail-closed: conic (deg-2 in the
    # fiber var) only — ConicScopeError on deg!=2; T-den fiber factors
    # squarefree only. Lazy module (sympy loads on first touch).
    "conic_thirdkind",
    # epslimit — engine-agnostic eps->0 Laurent extractor over an external
    # engine's (eps_i, value_i) grid (auto pole order, log(eps) resonance
    # detection by leave-one-out, cross-method achieved-digit estimate) +
    # richardson_boundary (series summed AT its radius of convergence with
    # known tail exponents). Lazy module; NOT epsfan (see docstring above).
    "epslimit",
    # flintexport — FLINT fmpq_mpoly fraction walker: bulk exact rational-
    # function algebra for DE-connection assembly; the exporter behind the
    # monomial-json format load_monomial_json reads. Lazy module; requires
    # python-flint (ImportError on first touch otherwise).
    "flintexport",
    # capability flag (python-flint OPTIONAL accel; pure mpmath otherwise)
    "HAVE_FLINT",
]

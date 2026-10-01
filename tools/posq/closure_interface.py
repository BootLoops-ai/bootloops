"""posq shared core — closure interface, READ-ONLY evaluation tier.

An evaluation boundary for plug-in log-likelihood closures: a consumer
writes a closure to the harness spec and gets register-labeled evidence
evaluations back; the closure objects served here are the harness spec,
spec-compliance checking, and MC-register evidence evaluation.

THE BOUNDARY. This module is the sanctioned surface: consumers import THIS
module and nothing behind it. The backing engine (the pinned R3
certified-evidence build) is imported privately, sha-pinned at load, and
never re-exported. Everything returned is a plain JSON-able dict. All calls
are read-only evaluation — nothing here mutates any stored artifact, writes
any file, or reaches production drivers.

WHAT IS STAGED vs WHAT IS MISSING — see coverage() for the mechanical answer.
Objects that do not exist in the pinned build are not improvised: their calls
raise ClosureObjectMissing carrying the precise reason and the smallest
unlock (no silent fallbacks).

REGISTER. Every evidence number leaving this surface is
register='R3-certified-calibrated': float-class, graded BY exact/enclosure
anchors, never feeding a certified claim. The shipped calibration was
measured on JC anchors; IT DOES NOT TRANSFER to plug-in closures — every
result carries this caveat verbatim (calibration_transfer field). Re-grade
against exact/enclosure anchors of the consumer's own object before leaning
on plug-in error bars.
"""

import hashlib
import os
import sys

__all__ = [
    "HARNESS_SPEC", "MC_REGISTER", "REGISTER_LAW",
    "ClosureInterfaceError", "ClosureObjectMissing",
    "verify_pins", "check_closure", "evaluate_evidence",
    "infer_hazard_field", "posterior_draws", "coverage",
]

# --------------------------------------------------------------------------
# The sanctioned build (the pinned R3 certified-evidence package) — private,
# sha-pinned.
# --------------------------------------------------------------------------

# backing-engine root comes from the environment (the loader
# below fails LOUD when it is unset or the pinned files are absent)
_R3_ROOT = os.environ.get("POSQ_R3_ROOT", "")

# sha256 pins of every upstream file the evaluation path can execute, taken
# against a build already exercised end-to-end (import/run SUCCESS,
# PlugInLoglik wired with no code changes). Any drift fails LOUD.
_PINNED = {
    "certified_evidence/__init__.py":
        "728818785d56e584dcfda27209c11ada6b90d7231203d5b96f2d337f09916d46",
    "certified_evidence/calibration.py":
        "39b3ac6ac470fef6ebe9b1d94d1e5d12b04170402c12959a61453dd6dd769f31",
    "certified_evidence/dataset.py":
        "b589e050286fa20194eb4cfed9c54d1de03055cfe66d6756bb06b96bc6aa3231",
    "certified_evidence/engine.py":
        "221563bd4d322afd989815dea1eba106c81dd43980f25442ee99d2c90a8a3d83",
    "certified_evidence/estimator.py":
        "03de8ae1b5e24701b8745eca9f65e3ab9fa0ad07e3dcfa944969bb0c562ac80d",
    "certified_evidence/fastlik.py":
        "c2f5fd60316018bf88e1e55d5dc6751e5ed0506c99140f3aa5cd47a3344cc39e",
    "certified_evidence/fastlik_k2p.py":
        "7f3a16c0ca668d599f3fa15b93ed0214aae046a23785bca719273ae9c53fb1ea",
    "certified_evidence/gate.py":
        "213345842afe157794b73f20df1954f7da636fb7a4b40d2c824c29f1360e9837",
    "certified_evidence/k2p.py":
        "466c9da486e8a5b504287217b169fb7278f5d84923e991adbef158858237f2e0",
    "certified_evidence/models.py":
        "52cdba060e208dd29d5e4d28ff7102d43551079cf69aa549dbb7eebfa5233352",
    "certified_evidence/prior.py":
        "a06bca54426f48c93b5e5776eba988be1acc7bd793f8d8595e8ceb9ac51b4b97",
    "certified_evidence/registers.py":
        "83cba08e17d3e72855aab378b0201c93c5d6a89f262f362cea591ea33aa2893c",
    "certified_evidence/survey.py":
        "0dee9d4244c4d327e205054f7fcbc92b6fe4d6dee03c66491a16eab664d31fe1",
    "certified_evidence/data/calibration_v1/calibration.json":
        "10f1258b880c89f5cb3c00d23442cf631efba7a8ba4f6ff16831ec33952c7038",
}

MC_REGISTER = "R3-certified-calibrated"

REGISTER_LAW = (
    "Float-class law (R3, applied to every number from this surface): "
    "(1) an R3 estimate never feeds a certified (R0/R1/R2) claim; "
    "(2) R3 is graded BY exact/enclosure anchors, never the reverse; "
    "(3) every estimate carries register='R3-certified-calibrated' and a "
    "versioned calibration_ref; "
    "(4) the shipped calibration was measured on JC anchors and does NOT "
    "transfer to plug-in closures — re-grade before leaning on plug-in "
    "error bars."
)

_CALIBRATION_TRANSFER_CAVEAT = (
    "Plug-in closures get no free calibration transfer: the shipped "
    "calibration (JC anchors) labels the estimator, not YOUR closure's "
    "error. Re-grade against exact/enclosure anchors of your own object "
    "before leaning on error bars (the R3 register law)."
)

HARNESS_SPEC = """loglik-closure harness spec (the contract plug-in closures are written to)
===========================================================================
A closure is ONE Python callable:

    loglik(U) -> L

  U : numpy float array, shape (n, d), points in the OPEN unit cube
      [0,1]^d (the harness may probe any point of the cube).
  L : numpy float array, shape (n,), log-likelihood values. -inf is
      allowed (hard zeros); NaN is a spec violation. The callable must
      be vectorized, deterministic, and side-effect-free.

The closure OWNS the map from the cube to its parameter space and folds
that map's Jacobian into the returned values, so that

    Z = E_{u ~ U[0,1]^d} [ exp(loglik(u)) ]

is exactly the evidence the closure author means. Covariates, offsets,
and hierarchy (e.g. a hazard field's per-slot covariates, per-bin exposure
offsets, gamma slot heterogeneity — marginalized or
raised to explicit cube dimensions) live INSIDE the closure; the harness
never sees them. Validate with check_closure(); evaluate with
evaluate_evidence(). Register of every output: R3-certified-calibrated
(float-class; see REGISTER_LAW and the calibration-transfer caveat)."""


class ClosureInterfaceError(Exception):
    """Interface-level failure (pin drift, spec violation, bad arguments)."""


class ClosureObjectMissing(ClosureInterfaceError):
    """A requested closure object does not exist in the pinned build.

    Carries the precise ``reason`` and the smallest ``unlock``.
    Never improvised around.
    """

    def __init__(self, obj, reason, unlock):
        super().__init__(f"{obj}: {reason} [unlock: {unlock}]")
        self.obj = obj
        self.reason = reason
        self.unlock = unlock


_ce = None  # cached private handle to the pinned upstream package


def verify_pins():
    """sha256-verify every pinned upstream file. Returns the pin dict.

    Raises ClosureInterfaceError on any mismatch or missing file — the
    sanctioned build is the ONLY build this surface will evaluate on.
    """
    if not _R3_ROOT:
        raise ClosureInterfaceError(
            "POSQ_R3_ROOT unset — path to the pinned R3 build is "
            "required (fail-closed; no default)")
    for rel, want in _PINNED.items():
        path = os.path.join(_R3_ROOT, rel)
        if not os.path.isfile(path):
            raise ClosureInterfaceError(
                f"pinned upstream file missing: {path} — the pinned "
                f"R3 build is not present; do not substitute another build")
        with open(path, "rb") as fh:
            got = hashlib.sha256(fh.read()).hexdigest()
        if got != want:
            raise ClosureInterfaceError(
                f"upstream drift: {rel} sha256 {got} != pinned {want}. "
                f"The interface refuses to run on an unpinned build; "
                f"re-pin only after re-validating against the producing build")
    return dict(_PINNED)


def _load():
    """Private: verify pins, then import the upstream package (cached)."""
    global _ce
    if _ce is None:
        verify_pins()
        if _R3_ROOT not in sys.path:
            sys.path.insert(0, _R3_ROOT)
        import certified_evidence as ce  # noqa: private, never re-exported
        _ce = ce
    return _ce


# --------------------------------------------------------------------------
# STAGED OBJECT 1 — the harness spec + spec-compliance check
# --------------------------------------------------------------------------

def check_closure(loglik, dim, seed=0, n=64):
    """Read-only spec-compliance check of a loglik closure (HARNESS_SPEC).

    Probes the closure at ``n`` deterministic cube points (plus the n=1
    edge case), checks: callable, (n,)-shaped float output, no NaN
    (-inf allowed), vectorization consistency (batch == pointwise), and
    determinism (identical repeated call). Raises ClosureInterfaceError
    naming the first failed check; returns a JSON-able receipt on pass.
    """
    import numpy as np
    if not callable(loglik):
        raise ClosureInterfaceError("spec: loglik must be callable")
    d = int(dim)
    if d < 1:
        raise ClosureInterfaceError("spec: dim must be >= 1")
    rng = np.random.default_rng(int(seed))
    # interior points only (open cube); boundary behavior is the closure's
    # own business as long as values are float/-inf, never NaN
    U = 0.5 * (rng.random((int(n), d)) + 0.25)
    L1 = np.asarray(loglik(U))
    if L1.shape != (int(n),):
        raise ClosureInterfaceError(
            f"spec: loglik((n,{d})) must return shape ({int(n)},), "
            f"got {L1.shape}")
    if not np.issubdtype(L1.dtype, np.floating):
        raise ClosureInterfaceError(
            f"spec: loglik must return float values, got dtype {L1.dtype}")
    if np.isnan(L1).any():
        raise ClosureInterfaceError(
            "spec: loglik returned NaN (hard zeros must be -inf, never NaN)")
    L2 = np.asarray(loglik(U))
    if not np.array_equal(L1, L2, equal_nan=False):
        raise ClosureInterfaceError(
            "spec: loglik is not deterministic (identical input, "
            "different output)")
    Lp = np.asarray(loglik(U[:1]))
    if Lp.shape != (1,) or not (
            np.isneginf(L1[0]) and np.isneginf(Lp[0])
            or np.isclose(Lp[0], L1[0], rtol=1e-12, atol=0.0)):
        raise ClosureInterfaceError(
            "spec: batch and pointwise evaluation disagree "
            f"({Lp!r} vs {L1[0]!r}) — closure is not consistently vectorized")
    return {
        "ok": True, "dim": d, "n_probes": int(n), "seed": int(seed),
        "finite_fraction": float(np.isfinite(L1).mean()),
        "checks": ["callable", "shape", "float-dtype", "no-NaN",
                   "deterministic", "vectorization-consistent"],
    }


# --------------------------------------------------------------------------
# STAGED OBJECT 2 — MC-register evidence evaluation of a closure
# --------------------------------------------------------------------------

def evaluate_evidence(loglik, dim, budget=100000, seed=0, name="closure",
                      mode="v2"):
    """Read-only MC-register evidence evaluation of a loglik closure.

    Runs check_closure first (small probe), then the pinned R3
    estimator (pilot -> product-Beta/copula proposal -> scrambled-Sobol
    RQMC-IS) on the closure. Returns a plain dict:

      logZ                 float estimate of ln E_u[exp(loglik(u))]
      register             'R3-certified-calibrated' (float-class)
      register_law         the float-class law, verbatim
      calibration_version  the versioned calibration the ESTIMATOR carries
      calibration_transfer the non-transfer caveat (binding for closures)
      certificate          the estimator's machine-checkable claim record
      + mode/budget/n_evals/seed/name/dim bookkeeping.

    Nothing beyond this dict crosses the boundary; no upstream object is
    returned. Read-only: no files written, no state kept.
    """
    ce = _load()
    spec = check_closure(loglik, dim, seed=0, n=32)
    if mode not in ("v2", "v1"):
        raise ClosureInterfaceError(
            "mode must be 'v2' (default) or 'v1' (the frozen earlier spec); "
            "the survey mode is JC-fold-native and not part of this surface")
    model = ce.PlugInLoglik(loglik=loglik, dim=int(dim), name=str(name))
    res = ce.Estimator(model=model, mode=mode).evidence(
        None, budget=int(budget), seed=int(seed))
    return {
        "logZ": float(res.logZ),
        "register": res.register,
        "register_law": REGISTER_LAW,
        "calibration_version": res.calibration_ref.version,
        "calibration_transfer": _CALIBRATION_TRANSFER_CAVEAT,
        "certificate": res.certificate(),
        "spec_check": spec,
        "mode": res.mode, "budget": int(res.budget),
        "n_evals": int(res.n_evals), "seed": int(res.seed),
        "name": str(name), "dim": int(dim),
    }


# --------------------------------------------------------------------------
# MISSING OBJECTS — mechanical refusals, never improvisation
# --------------------------------------------------------------------------

_MISSING_POSTERIOR = dict(
    reason=(
        "no sanctioned posterior-sampling surface exists behind this "
        "boundary: the validated packages ship EVIDENCE values (logZ) at "
        "the R3 register, not posterior draws; the pilot MH behind the "
        "boundary is proposal scaffolding — uncalibrated, ungated, and it "
        "stays behind this boundary"),
    unlock=(
        "a posterior-draw surface must be built, calibrated "
        "against exact anchors, and gated — then this interface grows the "
        "call; until then a hazard-field posterior over {h_b} "
        "cannot be served from this surface"),
)


def infer_hazard_field(*_args, **_kwargs):
    """Posterior over a hazard field {h_b} — NOT PROVIDED. Raises, precisely.

    Also note: the hazard-field LOGLIK itself is by design the
    consumer's contribution — write it to
    HARNESS_SPEC and it evaluates here via evaluate_evidence today.
    """
    raise ClosureObjectMissing("hazard-field posterior inference",
                               **_MISSING_POSTERIOR)


def posterior_draws(*_args, **_kwargs):
    """Posterior draws for any closure — NOT PROVIDED. Raises, precisely."""
    raise ClosureObjectMissing("posterior draws", **_MISSING_POSTERIOR)


def coverage():
    """The closure objects, object by object: STAGED / MISSING."""
    return {
        "loglik-closure harness spec": dict(
            status="STAGED",
            call="HARNESS_SPEC + check_closure()",
            note="the contract a plug-in closure is written to; "
                 "covariates/offsets/heterogeneity live inside the closure"),
        "MC-register evaluation": dict(
            status="STAGED",
            call="evaluate_evidence() -> logZ @ R3-certified-calibrated",
            note="evidence values only; calibration-transfer caveat binds"),
        "register law / labeling": dict(
            status="STAGED",
            call="MC_REGISTER, REGISTER_LAW; every result self-labels"),
        "hazard-field loglik closure": dict(
            status="MISSING-BY-DESIGN",
            call="none",
            note="the closure is the consumer's contribution; a "
                 "hazard-shaped toy fixture "
                 "(closure_fixtures F3) demonstrates the wiring"),
        "posterior over {h_b} / rate re-inference posteriors":
            dict(status="MISSING",
                 call="infer_hazard_field()/posterior_draws() raise "
                      "ClosureObjectMissing",
                 note=_MISSING_POSTERIOR["reason"]),
        "posterior-draw MC runner": dict(
            status="MISSING",
            call="none",
            note="no sanctioned posterior-draw runner exists behind this "
                 "boundary"),
    }

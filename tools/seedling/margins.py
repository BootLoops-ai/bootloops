"""Registered-M2 margin prediction.

WRAPS a registered frozen model — never reimplements it: the ridge
fit/predict and the feature definitions are IMPORTED from a model directory
you supply (SEEDLING_ML_DIR; it must contain train_eval.py + features.py).
ANY change to those files or to the fit recipe below requires a NEW review
gate; this module only routes data into them.

FROZEN MODEL (the only fit this module ever performs): M2 ridge (lam=1.0)
on the registered per-t-class training labels (SEEDLING_LABELS JSON:
labels_meff = per-t-class minimal sufficient uniform closure margin m in the
(t+m, 1, m) staging). Registered reference verdicts the wrapped model
reproduces exactly: held-out class m2_raw = 4.047, pred ceil = 5 == truth;
family-transfer margins {3:4, 4:4, 5:5, 6:5} on the transfer family.

FAMILY TRANSFER: per-class features for a NEW family are computed by the
executed transfer recipe (generic n_sp), identical to
features.class_features on the training family. Registered feature set only;
NO posterior quantities.

The margin model + labels are external registered artifacts, NOT part of this
package: `--margins m2` requires SEEDLING_ML_DIR and SEEDLING_LABELS to be
set (loud error otherwise). `--margins const` (the registered baseline arm)
is fully self-contained.

HARD FLOORS (chain-floor rule; the model may only ADD margin above floors):
  - meff floor: margin = max(1, ceil(raw - 1e-9))  (train_eval line, frozen)
  - s chain floor: emitted from pinner.chain_s_floor (REUSED, not
    duplicated): s >= 1 whenever any target has ISP powers or sits below a
    top sector. Consumers (the run subcommand / build_schedule) must apply
    it; an unsafe s=0 override must be explicit and logged.

--margins=const (B-const baseline, the registered comparison arm): a single
registered constant margin (B_CONST = 6, the max over the frozen train
labels), same floors. The const arm is fully self-contained: it loads
neither the model directory nor the labels.

Every emission carries provenance: model_hash (sha256 over the fitted
w/mu/sd + labels sha + feature keys), label provenance (path + sha256),
floor-applied flags.
"""
import hashlib
import json
import math
import os
import sys

from . import pinner

_ML_DIR = os.environ.get("SEEDLING_ML_DIR", "")
DEFAULT_LABELS = os.environ.get("SEEDLING_LABELS", "")
HOLD_OUT_T = 5          # frozen train/hold-out split (registered; never re-rolled)
LAM = 1.0               # registered ridge lambda (train_eval.m2_ridge default)
B_CONST = 6             # registered constant margin (B-const baseline arm):
                        # max over the frozen train labels; the model-gated
                        # battery leg cross-asserts it against
                        # fit_frozen()["b_const"] whenever the labels are
                        # present, so drift cannot pass silently.

_CACHE = {}


def _ml():
    """Import the registered modules (they are the model; never copied)."""
    if not _ML_DIR:
        raise RuntimeError(
            "margins m2 mode needs the registered frozen model: set "
            "SEEDLING_ML_DIR to the directory holding train_eval.py + "
            "features.py (and SEEDLING_LABELS to the labels JSON). "
            "The const arm (--margins const) is fully self-contained "
            "and never loads them.")
    if _ML_DIR not in sys.path:
        sys.path.insert(0, _ML_DIR)
    import features as _features            # registered features (frozen)
    import train_eval as _train_eval        # registered M2 ridge (frozen)
    return _features, _train_eval


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def fit_frozen(labels_path=DEFAULT_LABELS):
    """Fit (or return cached) the FROZEN registered M2 model. Returns dict:
    model (m2_ridge output), train_labels, b_const, feature_keys,
    provenance."""
    if not labels_path:
        raise RuntimeError(
            "margins m2 mode needs the registered training labels: set "
            "SEEDLING_LABELS to the labels JSON (labels_meff), or pass "
            "--labels. The const arm (--margins const) is fully "
            "self-contained and never loads them.")
    if labels_path in _CACHE:
        return _CACHE[labels_path]
    features, train_eval = _ml()
    data = json.load(open(labels_path))
    raw = data.get("labels_meff", data)
    labels = {int(k): v for k, v in raw.items() if v is not None}
    train_t = [t for t in sorted(labels) if t != HOLD_OUT_T]
    train_labels = {t: labels[t] for t in train_t}
    _tg, tab, census = features.load_family()
    cls = features.class_features(census, tab)
    X = [[cls[t][k] for k in train_eval.FEATURE_KEYS] for t in train_t]
    y = [train_labels[t] for t in train_t]
    m2 = train_eval.m2_ridge(X, y, lam=LAM)
    labels_sha = _sha256_file(labels_path)
    model_hash = hashlib.sha256(json.dumps({
        "w": [repr(float(v)) for v in m2["w"]],
        "mu": [repr(float(v)) for v in m2["mu"]],
        "sd": [repr(float(v)) for v in m2["sd"]],
        "labels_sha256": labels_sha,
        "feature_keys": train_eval.FEATURE_KEYS,
        "lam": LAM, "holdout_t": HOLD_OUT_T,
    }, sort_keys=True).encode()).hexdigest()
    out = {
        "model": m2,
        "train_labels": train_labels,
        "b_const": max(train_labels.values()),
        "feature_keys": list(train_eval.FEATURE_KEYS),
        "train_class_features": cls,
        "provenance": {
            "model": "M2-ridge (registered; SEEDLING_ML_DIR/train_eval.py)",
            "model_hash": model_hash,
            "labels_path": labels_path,
            "labels_sha256": labels_sha,
            "labels_provenance": ("registered labels_meff JSON "
                                  "(grid + perturbations)"),
            "train_classes": train_t,
            "holdout_t": HOLD_OUT_T,
            "lam": LAM,
            "gate": "registered training gate (frozen)",
        },
    }
    _CACHE[labels_path] = out
    return out


def family_class_features(census, tab, n_sp):
    """Per-class features for ANY family — the executed family-transfer
    recipe, generic n_sp. Registered feature set only. Identical to
    features.class_features on n_sp=15 families (asserted at registration)."""
    nontrivial = set(census)
    per = {}
    for sid, t in census.items():
        children = sum(1 for b in range(n_sp)
                       if (sid >> b) & 1 and (sid & ~(1 << b)) in nontrivial)
        parents = sum(1 for b in range(n_sp)
                      if not (sid >> b) & 1 and (sid | (1 << b)) in nontrivial)
        dists = [census[ts] - t for ts in tab
                 if ts != sid and (ts & sid) == sid and ts in census]
        cdist = min(dists) if dists else (0 if sid in tab else n_sp)
        per[sid] = dict(t=t, m=n_sp - t, closure_dist=cdist,
                        n_children=children, n_parents=parents)
    classes = {}
    for sid, f in per.items():
        if sid in tab:
            continue                 # closure classes only (registered unit)
        classes.setdefault(f["t"], []).append(f)
    out = {}
    for t, rows in sorted(classes.items()):
        n = len(rows)
        out[t] = dict(t=t, m=n_sp - t, n_sectors=n,
                      mean_closure_dist=sum(r["closure_dist"]
                                            for r in rows) / n,
                      mean_children=sum(r["n_children"] for r in rows) / n,
                      mean_parents=sum(r["n_parents"] for r in rows) / n)
    return out


def predict_margins(census, tab, targets, top_sectors, mode="m2", n_sp=15,
                    labels_path=DEFAULT_LABELS):
    """Per-closure-class margin prediction for a family.

    census: dict sector -> t (nontrivial sectors); tab: support.support_table
    output; targets: list[support.Target]; top_sectors: iterable.
    mode: "m2" (frozen registered model) | "const" (B-const baseline arm).

    Returns dict with per_class {t: {raw, margin, floor_applied}}, s_floor
    (pinner.chain_s_floor — REUSED chain-floor rule), b_const, provenance."""
    if mode not in ("m2", "const"):
        raise ValueError(f"margins mode {mode!r} not in ('m2','const')")
    cls = family_class_features(census, tab, n_sp)
    s_floor = pinner.chain_s_floor(targets, top_sectors)
    common_prov = {
        "mode": mode,
        "n_sp": n_sp,
        "features_extrapolated": (n_sp != 15),
        "s_floor_rule": "pinner.chain_s_floor (chain-floor rule; "
                        "reused, not duplicated)",
        "s_floor": s_floor,
        "meff_floor": 1,
    }
    if mode == "const":
        # B-const baseline arm: fully self-contained — the registered
        # constant margin only, never the frozen model or labels (those are
        # external artifacts; the const arm must work without them).
        per = {t: {"raw": None, "margin": B_CONST, "floor_applied": False}
               for t in cls}
        prov = {
            "model": "B-const baseline (registered constant margin arm; "
                     "self-contained)",
            "b_const": B_CONST,
            "b_const_provenance": ("registered constant = max over the "
                                   "frozen train labels; cross-asserted "
                                   "against fit_frozen()['b_const'] by the "
                                   "model-gated battery leg"),
        }
        prov.update(common_prov)
        return {"mode": mode, "per_class": per, "s_floor": s_floor,
                "b_const": B_CONST, "provenance": prov}
    frozen = fit_frozen(labels_path)
    _features, train_eval = _ml()
    per = {}
    for t, f in cls.items():
        raw = train_eval.m2_predict(frozen["model"],
                                    [f[k] for k in frozen["feature_keys"]])
        if raw is None or not math.isfinite(raw):
            # Degenerate-input guard: a non-finite prediction means
            # degenerate features reached the frozen model — fail LOUDLY,
            # never emit a NaN-derived margin. (Input guard only; the
            # registered model itself is untouched.)
            raise ValueError(
                f"margins: non-finite m2 raw prediction {raw!r} for "
                f"class t={t} (features {f}) — degenerate input")
        margin = max(1, int(math.ceil(raw - 1e-9)))   # frozen meff floor
        floor_applied = margin > math.ceil(raw - 1e-9)
        per[t] = {"raw": raw, "margin": margin,
                  "floor_applied": bool(floor_applied)}
    prov = dict(frozen["provenance"])
    prov.update(common_prov)
    return {"mode": mode, "per_class": per, "s_floor": s_floor,
            "b_const": frozen["b_const"], "provenance": prov}

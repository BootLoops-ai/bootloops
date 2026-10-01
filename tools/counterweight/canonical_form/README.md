# canonical_form — UT-basis / ε-factorized-DE builder

Sits between the IBP cache (Kira masters + DE matrix) and the value-fit
(`gatekeeper.heldout_certify`). Implements the principle: **fit in the LS-defined canonical
UT basis, not the raw Kira basis**.

## What's implemented vs wrapped vs stubbed

| piece | status | engine |
|---|---|---|
| `Family(spec)` — parse amflow `family` block → Baikov rep | **implemented** | sympy |
| `leading_singularities()` — maximal-cut LS per master | **implemented** | sympy (iterated discriminant) |
| `ut_rotation()` → `{ut_tag:{raw_tag:coeff}}` for `heldout_certify` | **implemented** | — |
| `elliptic_period()` — ϖ₀ of y²=quartic via ellipk | **implemented** | mpmath |
| `verify_unit_pole()` — ε⁻ᴸ-is-constant gate on oracle data | **implemented** | reuses `heldout_cv._apply_ut_rotation` |
| `to_dlog()` — Ã = Σ A_a ∂log a → constant matrices | **implemented** | sympy partial fractions |
| `fuchsify()` / `factor_epsilon()` — Lee's algorithm | **wrapped** | `Counterweight` (Julia/Nemo) via `counterweight_bridge.jl` |
| Moser reduction (Barkatou–Moser shearing) of higher-order poles | **implemented in the Julia engine** | `Fuchsia.moser_reduce`, invoked by `Normalize.try_epsfactor` on non-Fuchsian input (`test/test_moser.jl` T1–T4). Genuinely irregular (Moser-irreducible) points remain a designed STOP. The Python `fuchsify()` wrapper keys on the pre-reduction certificate and still refuses non-Fuchsian input by name (`NotImplementedError` + recipe) |

## Minimum-viable usage (no engine needed)

```python
from canonical_form import Family, ut_rotation
from gatekeeper import heldout_certify, load_oracle_values

fam = Family(target["family"])
rot = ut_rotation(fam, target["masters"], point={"s": -3, "t": -7})
# rot = {'ut_box': {'box': '21/4'}, ...}   ← exact heldout_certify format

values, kin, _ = load_oracle_values([oracle_dir])
res = heldout_certify(values, kin, weight=2, basis_fn=..., basis_names=...,
                      ut_rotation=rot)
```

The rotation is **diagonal**: `J_tag = (1/LS_tag(point)) · I_tag`. Off-diagonal
D⁻ subsector subtraction (the full UT completion) is not shipped with this tool; a downstream
UT-completion step consumes the `LSResult` this tool
produces as its `LS` field.

## Per-master `maxcut` override

Automatic Baikov is exact for the **top sector**. For subsectors the
top-family Gram over-counts external momenta (residual-ISP discriminant chain
picks up spurious factors), and for γ₀=0 cases (e.g. equal-mass sunrise,
L=2,E=1) the standard rep is degenerate. In both cases supply the curve
directly — same schema as `tools/geotriage/graph_specs/*["maxcut"]`:

```python
masters = [{"tag": "sun", "indices": [1,1,1],
            "maxcut": {"poly": "z*(z-4)*((tt-z)**2 - 4*z)", "vars": ["z"]}}]
```

## Full ε-factorization (wrapped engine)

```python
from canonical_form import factor_epsilon, to_dlog
Atilde, T, report = factor_epsilon(A, x)       # → Julia Counterweight (Lee 1411.0911)
letters = to_dlog(Atilde, x, alphabet=[x, 1+x])
```

The bridge serializes A as strings, runs
`julia --project=$COUNTERWEIGHT_ROOT counterweight_bridge.jl`, and parses the JSON
back. The engine is validated against Henn's massless planar double box
(`$COUNTERWEIGHT_ROOT/test/validate_henn.jl`). If the engine STOPs (non-Fuchsian
input, non-semisimple obstruction) you get a `NotImplementedError` with the
exact recipe pointer.

## CLI

```
python canonical_form.py box1l.json --point '{"s":-3,"t":-7}' \
    --out out/box1l/ut_rotation.json
```

## Tests

```
python test_canonical_form.py
```

T1 box LS=1/(st) · T2 to_dlog exact · T3 sunrise elliptic ϖ₀ · T4 oracle wire
· T5 Counterweight bridge certify. Set `CANONICAL_FORM_SKIP_JULIA=1` to skip T5.

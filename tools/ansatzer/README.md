# ansatzer

Predicts, **before** you farm high-precision oracle points, the residual
symbol-ansatz dimension per weight (and per master, if a connection is
supplied) after the standard structural cuts. The residual is the number of
unknowns the value-fit must determine, hence — times a 1.3 safety factor — the
number of oracle points the fit needs. Run this first to know whether a
value-fit is affordable at all: a 7-letter alphabet can collapse to a handful
of unknowns per master while a 19-letter one leaves ~80 points to farm.

Pure Python; depends only on `numpy` and `sympy`.

## Method

At weight w the raw tensor symbol space over an n-letter alphabet has
dimension n^w. Six linear constraint families cut it down, in order:

1. **Integrability** (dA = A∧A): for every adjacent slot and environment, the
   wedge Σ c · dlog L_a ∧ dlog L_b must vanish as a 2-form — an exact-ℚ (or
   sampled mod-p) nullspace. This is the engine.
2. **First-entry**: leading letter restricted to physical thresholds
   (Mandelstams, masses).
3. **Last-entry**: trailing letter restricted by the connection/discontinuity
   structure.
4. **Steinmann**: double discontinuities in overlapping channels vanish, so
   two *adjacent* entries may not carry disjoint channel-tag sets. Letters
   sharing a tag are same-channel (iterated discontinuities in one channel
   are allowed); untagged letters are unconstrained. Active whenever the
   alphabet supplies `channel` tags (disable with `--no-steinmann`).
5. **Parity/reality**: if a per-letter ℤ₂ parity is supplied (radicals,
   σ-symmetry), project to the target sector.
6. **Coproduct / Δ^{(w−1,1)}** (sequential discontinuity): if the canonical DE
   is dF = ε (Σ_a A_a dlog ℓ_a) · F, then Disc_a F^(w) ∝ A_a · F^(w−1); for
   master m the last letter a is allowed only if row m of A_a is non-zero
   *and* the (w−1)-prefix lies in the already-cut weight-(w−1) space.
   Recursive; only the support of the connection (which entries are non-zero)
   is used.

The report gives the surviving dimension after each cut, the final residual,
the implied point count, and a verdict (`COLLAPSES` or `FULL`).

## Package contents

| file | role |
|---|---|
| `ansatzer.py` | the predictor (CLI + importable functions) |
| `alphabets/` | input fixtures, validated by the test battery |
| `test_ansatzer.py` | test battery |
| `p4_J25_collapse.json` | sample output report |

## Input

`alphabet.json` (native schema; the output of
[`../landau-alphabet/landau_alphabet.py`](../landau-alphabet/) is
auto-adapted):

```json
{
  "name": "...", "vars": ["s","t"],
  "letters": {"L_s":"s", "L_t":"t", ...},
  "first_entry": [...], "last_entry": [...],
  "parity": {"O_A3":[1]}, "parity_target": [1],
  "channel": {"s":["S"], ...}
}
```

Algebraic letters must be supplied on a **rationalized chart**, so that all
dlogs are rational in `vars` — letters left on the raw chart give wrong
counts. Example: `O_A3 = ((1+y)/(1-y))²` on the chart `t = 4/(1+y²)`
reproduces the radical-aware counts exactly (see `alphabets/p4_J25.json`).

`channel` (optional) maps a letter to its channel tags — a list, or a single
bare tag string. Supplying it activates the Steinmann cut; omit the block (or
pass `--no-steinmann`) to skip it. A letter listed under several tags is
compatible with a neighbor sharing any one of them.

`conn.json` (optional): `{"masters":[...], "A":{letter:[[i,j,"q"],...]}}`.
Only the support (which `(i,j)` are non-zero) is used.

## Usage

```bash
python3 ansatzer.py --alphabet A.json [--connection C.json] \
        --wmax 4 [--method auto|exact_Q|modp] [--timeout 120] \
        [--compare B.json] --out report.json
```

- `--timeout` caps each weight's integrability computation (seconds).
- `--compare` prints a side-by-side report for a second alphabet.

A per-orbit mode replaces the global collapse table with per-orbit residuals
K_w, computed by a recursive integrable-symbol builder with an adaptive mod-p
convergence check (6 up to `--npts` sample points):

```bash
python3 ansatzer.py --per-orbit --alphabet A.json \
        [--wmax 4] [--npts 96] [--out r.json]
```

Output: `{"Kw_LE": ..., "Kw_full_noLE": ..., "converged": ..., "wallclock_s": ...}`.

The per-orbit path loads its input through the same adapter as the global
path (collapse schema or `landau_alphabet.py` output; `first_entry`/
`last_entry` default to all letters). Channel tags are not applied on this
path — the Steinmann cut lives in the global table only.

## Output

Table columns:

```
weight | raw | integrability | first | last | steinmann | parity | coproduct | RESIDUAL | pts_needed
```

- **exact_Q** mode: each column is the subspace *dimension* after that cut.
- **modp** mode: `first/last/steinmann/parity/coproduct` are word-prefilter
  *counts*; the `integrability` = `RESIDUAL` column is the converged mod-p
  nullspace dimension on the surviving words.

## Honesty flags

The report states exactly what it did and did not compute:

- `coproduct: SKIPPED (no connection supplied)` — integrability + entry
  conditions alone.
- `steinmann: applied (…) | vacuous (no channel tags) | SKIPPED (disabled)` —
  whether the adjacent-channel veto ran.
- `rank_method: exact_Q | modp | trivial` — per weight.
- `bound: ≥N` — integrability timed out at that weight, so the residual is an
  **upper bound only**. Do not budget a point farm off a bounded residual;
  re-run with a longer timeout, exact mode, or `--per-orbit` first.

## Validated fixtures (`alphabets/`)

Each fixture's expected counts are checked by the battery. Fixture names are
physics-family tags.

| file | checks |
|---|---|
| `trivial_2var.json` | hand-checkable 2-letter case: w2 4→3, w3 8→4 |
| `p4_J23.json` | 3-letter rational: integrable tower (3,7,15,31), first-entry residual (2,4,8,16) |
| `p4_J25.json` | 4-letter, parity-odd: integrable (12,32,80) at w2–4, residual (3,8,20) — exact match on the rationalized chart |
| `p4_J23_conn.json` | toy connection — demonstrates the per-master coproduct cut |
| `synth_7letter.json` | synthetic 7-letter alphabet (invented letters, formal parity grading, three channel tags): w1 cuts hand-checkable (first 5, residual 3), w2 Steinmann word count hand-checkable (9), w2/w3 exact-Q dims pinned with and without the channel block, strong collapse |
| `5pt_2mass_nonplanar.json` | 19-letter two-mass five-point family: little collapse, large residual |

Further family alphabets in `alphabets/` follow the same schema.

## Tests

```bash
python3 test_ansatzer.py
```

T1–T4, T6 (Steinmann cut: hand-checkable 2-letter case, exact-Q/mod-p
agreement, honesty flags) and T7 (per-orbit input adapter) are
self-contained (the fixtures above). T5 exercises the
`landau_alphabet.py`-output adapter against an engine-produced
`alphabet.json` (point `COLLAPSE_T5_ALPHABET` at one) and skips cleanly
when none is supplied.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

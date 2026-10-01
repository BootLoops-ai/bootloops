# RECEIPT witness format — v1.0

One witness file certifies ONE reduction-table row at ONE (kinematic point,
prime). The witness is the contract between any solver and any verifier.

## Design invariant

**The verify path is strictly solver-agnostic.** A witness makes no
reference to — and a verifier may make no assumption about — how the table
was produced (kira, B2FT, Wiedemann, cohomology pairing, a neural guesser,
anything). The only objects in the contract are:

- the **system**: an ordered list of sparse F_p rows `R_i`
  (`dict col_id -> value mod p`), evaluated at the witness's (point, prime);
- the **witness**: target column `t`, coefficient dict `c`, sparse
  multiplier vector `lam`;
- the **identity**:

```
sum_i lam[i] * R_i  ==  e_t  -  sum_m c[m] * e_m      (mod p)
```

If the identity holds, the row `I_t = sum_m c[m] I_m` is an exact linear
consequence of the system at this (point, prime), whoever produced it. This
is what makes any future reduction algorithm a drop-in replacement: trust
attaches to the receipt, not to the engine (see README — the
wrong-table exhibit is the proof-by-demonstration).

Column ids are opaque integers. Solver-specific meaning (kira weights,
index tuples, sector masks) lives in OPTIONAL label fields and in adapters;
the core checker (`core.py`) never reads them and never imports adapter
code.

## Schema (JSON, one object per file)

```json
{
  "receipt_version": "1.0",
  "kind": "lambda-witness",
  "identity": "sum_i lam[i]*R_i == e_target - sum_m c[m]*e_m (mod p)",
  "p": 2147483647,
  "point": {"d": 1234577, "eta": 87654321},
  "family": "myfam",
  "target": {"col": 141321603907593, "label": [1,0,1,0,1,1,3,0,0]},
  "c": {"<col_id>": 307413002, "...": "..."},
  "labels": {"<col_id>": [0,0,0,0,0,1,1,0,0]},
  "lam": {"n_rows": 7135, "idx": [3, 17, "..."], "val": [9021, 411, "..."]},
  "system": {"n_rows": 7135,
             "fingerprint": "<sha256, optional>",
             "source": "<free-text provenance, optional>"}
}
```

Field rules:

| field | required | meaning |
|---|---|---|
| `receipt_version` | yes | major version `1` accepted by this tool |
| `p` | yes | prime modulus (> 1) |
| `point` | recommended | free-form parameter dict (e.g. `{d, eta}`); needed by adapters to re-evaluate the system; the core treats it as opaque |
| `target.col` | yes | target column id (integer) |
| `target.label` | no | human-readable integral indices |
| `c` | yes | `col_id -> value mod p`, convention `I_t = sum_m c[m] I_m`; zero values are dropped |
| `labels` | no | `col_id -> indices` for the `c` keys |
| `lam` | yes | sparse multiplier: `n_rows` = system row count, parallel `idx`/`val` arrays; indices unique, in `[0, n_rows)`; **`idx` indexes the system's row ORDER**, which is part of the system definition |
| `system.fingerprint` | no | sha256 of the evaluated system (`core.system_fingerprint`: row order as given, columns sorted per row, `p` prefixed). When present and the verifier computed its own, a mismatch is a FAIL |
| `family`, `kind`, `identity`, `system.source` | no | provenance decoration only |

Size class (measured): lambda density mean 3.3% at the 7,135-row
reference-family scale, flat certificate <= 32 KB/row/prime.

## Verification semantics (normative)

1. `len(rows) == lam.n_rows`, else FAIL (`reason: lam n_rows mismatch`).
2. If both fingerprints exist and differ: FAIL.
3. Residual `r = sum_i lam[i]*R_i mod p` (sparse accumulation, zero entries
   dropped).
4. `expect = {t: 1} + {m: -c[m] mod p for nonzero c[m]}`.
5. PASS iff `r == expect` exactly. Any deviation — one coefficient, one
   lambda entry, one extra residual column — is a FAIL. There are no
   tolerances.
6. Optional table mode: given a claimed row `c_claimed`, additionally
   require `c_claimed == c` after mod-p normalization and zero-dropping.
   A witness that HOLDS while `c_claimed != c` means the claimed table row
   is wrong (the wrong-table exhibit class) — reported as FAIL with reason
   `claimed table row != certified row`.

## rc discipline (CLI)

- `0` — every row verified PASS / no detector alarms
- `1` — at least one row FAILED / detector alarm fired
- `2` — malformed input (unreadable JSON, missing required fields,
  duplicate or out-of-range `lam.idx`, unsupported version, no witnesses)

Malformed = "cannot be interpreted"; a well-formed witness that does not
match the system is a FAIL (rc 1), never rc 2.

## Scope

- Certificates are **per-(kinematic point, prime)**. PASS at two primes at
  one point is exactly that; there is **no cross-point (Schwartz-Zippel)
  claim** and no symbolic-table claim.
- The witness certifies membership of the row in the ROW SPAN of the given
  system. It does not certify that the system itself is the right one or
  complete — that is the `detect` layer (false-master audit)
  plus staging provenance gates.
- 2-loop measurements do not derisk 3-loop claims (pre-registered scope).

## Legacy format (read-only)

Earlier-format pairs `c_<fam>_<p>_<w>.json` + `lam_<fam>_<p>_<w>.npy`
(dense int64 npy, tuple-keyed `c`) are readable through
`adapters/strata.read_legacy_pair` in `tools/trust/receipt` (tuple -> col-id
resolution via the accompanying masters map). New emissions always write v1.

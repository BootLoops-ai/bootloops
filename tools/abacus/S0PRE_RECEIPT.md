# S0PRE — the Q3/Q6 input-contract pins (binding)

The S0-PRE row of the M3T spec is satisfied by the two design pins below;
S0 runs under them.
Binding on abcount: Q3 additions (c) actual-input invariants on every output
receipt + declared isogeny data; (d) no inferred equivalence to the pinned
lattice — identification-grade claims need the pinned lattice as input or a
separately-receipted stage-level isometry. Q6: GLOBAL pin = P;
P(-1) internal-only with mandatory sign_convention field +
single-site receipted boundary conversion.

---

**Q3 — PERMISSIVE INPUT, STRICT LABELING.** The tool MAY accept any member of the isogeny
class (finite-index sub/overlattice) PROVIDED: (a) the induced polarization
type rides with it, and (b) if not principal, a declared isogeny to a principally
polarized member. TWO ADDITIONS, binding: (c) every output receipt carries the ACTUAL
input lattice's own invariants — (A,q), det, signature — plus the declared isogeny data,
so no downstream reader can mistake which object was counted; (d) the tool NEVER infers
equivalence to the pinned lattice: any identification-grade claim against P requires
either the pinned lattice itself as input or a separately-receipted isometry, checked
at STAGE level (a proper overlattice keeps the signature but changes
(A,q) — a mistaken identification stays dead by labeling, not by input prohibition).

**Q6 — GLOBAL PIN: P, the convention of record.** Rationale: every
receipt that crosses stages ultimately feeds one downstream record, and that record
is stated in the P convention — receipts of record match the record. P(-1) is PERMITTED as
engine-native internal representation under two conditions:
the mandatory sign_convention field on every artifact (as contracted, N5-exercised), and
boundary conversion at a SINGLE site — one conversion function, receipted, not per-stage
ad hoc flips. This pin is GLOBAL and inherited verbatim by every consuming stage.

Every stage proceeds under both pins.

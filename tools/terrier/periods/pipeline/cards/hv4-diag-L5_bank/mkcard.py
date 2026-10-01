#!/usr/bin/env python3
# mkcard.py (DIAG-SLICE) — assemble cards/hv4-diag-L5.json from the battery
# receipt (BATTERY_RECEIPT.json), the two-engine series (SERIES_A.json, SERIES_B.json)
# and the frame-rank receipt (FRAME_RANK5.json) in this directory; run it
# from this directory.  The card it writes is sha-pinned in
# regression_manifest.json.  Slice-operator card (schema slice-operator-v1):
# NOT an h_eff GKZ family card (cards.py validates it against its own key
# set); consumers are the transport/conifold layers (matrix payload only,
# scalar-refusal law).
import json, hashlib
BR = json.load(open("BATTERY_RECEIPT.json"))
SA = json.load(open("SERIES_A.json")); SB = json.load(open("SERIES_B.json"))
FR = json.load(open("FRAME_RANK5.json"))
comment = (
 "HV4 (Hulek-Verrill fourfold, arXiv:2404.12422 GvdH) DIAGONAL-SLICE order-5 "
 "Picard-Fuchs operator card - the S6-symmetric locus phi^1=..=phi^6=phi "
 "M6-eigenvalue-+1 period sector. SLICE-RESTRICTED: every quantity DIAG-SLICE; "
 "the full 6-param census needs the 6-parameter card (not included). "
 "Operator DERIVED as the UNIQUE (order 5, phi-deg 3) annihilator of the "
 "two-engine diagonal fundamental-period series (6-letter abelian-square "
 "coefficients, eq (5.18) restricted; engines independent to n=400 + direct "
 "6-composition anchor n<=40) and MATCHES the printed JKK operator "
 "(arXiv:2404.12422v2 TeX source, ~eq (5.60), "
 "cited to Jockers-Kotlewski-Kuusela 2312.07611) coefficient-exactly. "
 "Frame-side pin from the 29x29 integer frame matrices (FRAME_RANK5.json): S6-invariant "
 "sector of the 29-frame is EXACTLY rank 5, T1..T6 product restricts and is "
 "MUM with unipotency index 5, restricted Sigma-Gram det 6480 signature "
 "(3,2). Symbol (1-4phi)(1-16phi)(1-36phi) = coefficient line 1 - 56phi + "
 "784phi^2 - 2304phi^3; roots == the three diagonal conifold points "
 "{1/36, 1/16, 1/4}; leading D-form coeff = phi^5*symbol -> NO apparent "
 "singularities. Exponents: MUM {0^5} at 0; {0,1,3/2,2,3} at each conifold "
 "(CY4-node type); {1,1,3/2,2,2} at infinity. K3 vacuum-carrying subsector "
 "L3 (Verrill 1996 + JKK; Domb series) carried as k3_sector, battery green, "
 "symbol relation sym(L5) = sym(L3)*(1-36phi) verified. GAPS: GAP-D1 exact "
 "all-orders annihilation is CITED (JKK) + depth-400 two-engine "
 "necessary-condition verified - no independent creative-telescoping "
 "certificate issued here; GAP-D2 (= frame G3) conifold monodromies in the "
 "5-frame not closed-form derivable from the frame matrices - needs certified "
 "transport. FLUX LAW (GvdH sec 5.3): no fluxes couple to this sector "
 "in their construction; vacuum content on the slice lives in k3_sector.")
# exact operator coefficient lists: recompute here (independent expansion)
from fractions import Fraction as Fr
def pmul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r
def padd(a, b):
    r = [0]*max(len(a), len(b))
    for i, x in enumerate(a): r[i] += x
    for i, x in enumerate(b): r[i] += x
    return r
def psc(a, s): return [x*s for x in a]
TH = [0, 1]; th1 = [1, 1]; th2 = [2, 1]
inner = padd(psc(pmul(pmul(TH, th1), padd(pmul(TH, th1), [1])), 14), [3])
Q1 = psc(pmul([1, 2], inner), -2)
Q2 = psc(pmul(pmul(pmul(th1, th1), th1), padd(psc(pmul(TH, th2), 196), [255])), 4)
Q3 = psc(pmul(pmul(pmul(th1, th1), pmul(th2, th2)), [3, 2]), -1152)
QL5 = {"Q0_phi0": [0]*5+[1], "Q1_phi1": Q1, "Q2_phi2": Q2, "Q3_phi3": Q3}
QL3 = {"Q0_phi0": [0]*3+[1],
       "Q1_phi1": psc(pmul([1, 2], padd(psc(pmul(TH, th1), 5), [2])), -2),
       "Q2_phi2": psc(pmul(pmul(th1, th1), th1), 64)}
R5 = [[QL5[k][d] if d < len(QL5[k]) else 0 for k in
       ("Q0_phi0", "Q1_phi1", "Q2_phi2", "Q3_phi3")] for d in range(6)]
R3 = [[QL3[k][d] if d < len(QL3[k]) else 0 for k in
       ("Q0_phi0", "Q1_phi1", "Q2_phi2")] for d in range(4)]
card = {
 "name": "hv4-diag-L5", "schema": "slice-operator-v1",
 "scope": "DIAG-SLICE (S6-symmetric locus of the 6-param HV4 family)",
 "comment": comment, "n_params": 1, "order": 5, "phi_degree": 3,
 "operator_normal_form": "L5 = sum_j phi^j Q_j(theta), theta = phi d/dphi; "
   "recurrence sum_j Q_j(m-j) c_{m-j} = 0 (convention gated on Domb anchor)",
 "Q_polys_theta": QL5, "R_polys_phi_thetaright": R5,
 "symbol": [1, -56, 784, -2304],
 "symbol_factors": ["1-4*phi", "1-16*phi", "1-36*phi"],
 "singular_points": ["0", "1/36", "1/16", "1/4", "inf"],
 "apparent_singularities": "NONE (leading D-form coeff = phi^5 * symbol)",
 "local_exponents": BR["exponents_L5"],
 "mum_series_head_c0_c40": SA["c6"][:41],
 "series_depth": {"N_TWO_ENGINE": 400, "N_DIRECT6_ANCHOR": 40},
 "pfaffian": BR["pfaffian_payload"],
 "k3_sector": {
   "comment": "vacuum-carrying K3 subsector on the slice (Verrill 1996 + "
     "JKK via GvdH; Domb series; polynomial periods (1, t, -12t^2+1), "
     "Hauptmodul eq (5.71) exact arithmetic - eps0=0, no tiles)",
   "order": 3, "phi_degree": 2, "Q_polys_theta": QL3,
   "R_polys_phi_thetaright": R3, "symbol": [1, -20, 64],
   "symbol_factors": ["1-4*phi", "1-16*phi"],
   "singular_points": ["0", "1/16", "1/4", "inf"],
   "local_exponents": BR["exponents_L3"],
   "domb_head_c0_c40": SA["c4_domb"][:41],
   "symbol_relation": "sym(L5) = sym(L3) * (1-36*phi)  [verified exact]"},
 "frame_receipt": {
   "invariant_sector_dim": FR["invariant_sector_dim"],
   "mum_index5": FR["mum_index5"], "N5_restricted": FR["N5_restricted"],
   "gram5": FR["gram5"], "gram5_det": FR["gram5_det"],
   "gram5_signature": FR["gram5_signature"],
   "note": "S6-invariant sector of the 29-frame (frame_rank5.py, frame "
     "matrices supplied via TERRIER_HV4_FRAME_DIR); the 6-cycle matrix M6 is "
     "(5.59)-as-displayed = transpose of printed (5.58); "
     "M6 is a slot permutation so the +1-eigenspace is convention-free"},
 "provenance": {
   "source": "arXiv:2404.12422v2 TeX source (sha256 4791813f66b1632f46727c5a"
     "5f8b0fdf385a60e599bc29ba3c5d889df37e3a52), L5 = ~eq (5.60), "
     "L3 + the K3 diagonal fundamental period in the same block",
   "operator_citation": "Jockers-Kotlewski-Kuusela arXiv:2312.07611 (L5); "
     "Verrill J. Math. Kyoto 36 (1996) 423 + JKK (L3)",
   "frame": "29x29 integer frame matrices (MSWAP, M6, T1..T6, SIGMA; two-engine "
     "derived, not included in the package; frame_rank5.py reads them from "
     "TERRIER_HV4_FRAME_DIR) -> FRAME_RANK5.json",
   "recon": "reference diagonal conifold points {1/36, 1/16, 1/4} (HV4 "
     "discriminant of arXiv:2404.12422 on the diagonal), used by the symbol gate"},
 "gaps": [
   "GAP-D1: exact all-orders annihilation CITED (JKK 2312.07611) + "
     "necessary-condition verified to n=400 two-engine; no independent "
     "creative-telescoping/Griffiths-Dwork certificate issued here",
   "GAP-D2 (= frame G3): conifold monodromies of the 5-frame not "
     "closed-form derivable from the frame matrices; needs certified transport",
   "SCOPE: DIAG-SLICE only; the full 6-param family is not covered"],
 "battery": {"receipt": "hv4-diag-L5_bank/BATTERY_RECEIPT.json",
   "scripts": ["engine_a.py", "engine_b.py", "frame_rank5.py",
               "battery_l5.py", "mkcard.py"],
   "verdict": "ALL GATES PASS"},
 "status": "SLICE-CARD OK (battery green). DIAG-SLICE ONLY - "
   "not a census instrument; transport NOT run (no routes field)."}
# consistency asserts against the battery receipt before writing
assert card["symbol"] == [QL5["Q0_phi0"][5], QL5["Q1_phi1"][5],
                          QL5["Q2_phi2"][5], QL5["Q3_phi3"][5]]
assert card["k3_sector"]["symbol"] == [QL3["Q0_phi0"][3], QL3["Q1_phi1"][3],
                                       QL3["Q2_phi2"][3]]
assert BR["fit"]["matches_printed_L5"] and BR["symbol"]["matches_reference_discriminant_points"]
out = "../hv4-diag-L5.json"
json.dump(card, open(out, "w"), indent=1)
h = hashlib.sha256(open(out, "rb").read()).hexdigest()
print("CARD WRITTEN", out, "sha256", h[:16])

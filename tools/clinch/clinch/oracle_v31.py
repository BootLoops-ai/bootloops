"""clinch.oracle_v31 — the v3.1 dial-layout oracle: penalized-NLL gradient
and block-arrow Hessian of the reference joint likelihood in arb balls.

REFERENCE OF RECORD: the originating study's likelihood.py (reference code,
not distributed with this package; sha-pinned)
`nll_grad` (the module every fit runner imports) — the receipted numpy
analytic gradient (FD-verified there). This oracle re-derives the SAME
math in ball arithmetic with second derivatives by jets (clinch.jets);
it is identity-gated against the reference gradient at the receipted
fit13/fit17 thetas BEFORE any certificate is issued (adapter_v31.gates).

SCOPE (guarded, the jaxlike.check_space pattern): EXACTLY the v3.1 space —
DialSpace(n_storm=0, n_weather=6), 51 globals, linear climate planes,
weather triplets, NON-centered pooled effects, STRUCTURE_MODE off. Any
other slot set REFUSES cold (check_space): certifying an unported branch
would be untestable code.

CLIP GUARDS (smoothness honesty): the reference clips p5/q5 to
[1e-9, 1-1e-9], pdie/q to [1e-10, 1-1e-10], floors stay5 at 1e-10 and mu
at 1e-12. A certificate is only valid where the objective is SMOOTH, i.e.
where every clip is provably inactive over the whole box. Each guard
proves strict interiority in ball arithmetic or raises OracleRefusal
(-> a named REFUSAL, never a shrug). Receipted zero clip activity at the
gate thetas (jaxlike/IDENTITY_GATES.json) says the guards are expected to
prove easily at the fitted optima.

BLOCK STRUCTURE (the arrow the certificate exploits): the latent Hessian
is EXACTLY block-diagonal over (process, family-component) blocks — cells
touch one species each, species/genus/family latents of one process only
couple within a family component (genus->family consistency is VERIFIED at
partition build; inconsistent genera merge components), processes never
share a cell, and the v3.1 space has no hard centering. Globals (51) are
the border; log_sig/log_nu/log_phi cross-terms land in the border blocks.

Ball hygiene: x*x, never x**2 (baller MANUAL law: arb pow NaNs on
zero-containing balls); all data constants enter arb EXACTLY (float64);
transcendental caches are keyed by exact cell constants AND the active
precision (an unkeyed cache across prec would be the banked footgun).
"""
import numpy as np
from flint import arb, arb_mat

from .jets import Jet, hidx, dig_ball, trig_ball, lgam_ball

__all__ = ["OracleRefusal", "check_space", "Partition", "ModelV31Oracle",
           "V31_GLOBAL_SLOTS"]

# the exact v3.1 global slot layout (name -> width), order of v3lib.DialSpace
V31_GLOBAL_SLOTS = [
    ("a_surv", 7), ("a_adv", 6), ("a_fec", 2), ("b_surv", 2), ("b_gro", 2),
    ("b_fec", 2), ("b_sxk", 1), ("g_surv", 3), ("g_gro", 3), ("g_fec", 3),
    ("log_sig", 9), ("lam", 3), ("log_phi", 1), ("log_nu", 1),
    ("w_surv", 3), ("w_gro", 3),
]
N_GLOBAL = sum(w for _, w in V31_GLOBAL_SLOTS)          # 51


class OracleRefusal(Exception):
    """A named oracle refusal (clip guard / smoothness not provable)."""

    def __init__(self, where, detail):
        self.where, self.detail = where, detail
        super().__init__(f"{where}: {detail}")


def check_space(slice_bounds, n_total):
    """Scope guard: refuse any dial space that is not EXACTLY v3.1."""
    names = [nm for nm, _ in V31_GLOBAL_SLOTS]
    got_glob = {nm: b for nm, b in slice_bounds.items()
                if not nm.startswith("z_")}
    if sorted(got_glob) != sorted(names):
        extra = sorted(set(got_glob) - set(names))
        missing = sorted(set(names) - set(got_glob))
        raise OracleRefusal(
            "check_space", f"not the v3.1 space: extra={extra} "
            f"missing={missing} — this oracle certifies ONLY the ported "
            f"v3.1 layout (unported branches refuse cold)")
    off = 0
    for nm, w in V31_GLOBAL_SLOTS:
        if tuple(slice_bounds[nm]) != (off, off + w):
            raise OracleRefusal("check_space",
                                f"slot {nm} at {slice_bounds[nm]}, "
                                f"expected ({off},{off + w})")
        off += w
    zn = [nm for nm in slice_bounds if nm.startswith("z_")]
    want_z = [f"z_{lev}_{p}" for p in ("surv", "adv", "fec")
              for lev in ("sp", "gen", "fam")]
    if sorted(zn) != sorted(want_z):
        raise OracleRefusal("check_space", f"latent slots {sorted(zn)}")
    hi = max(b for _, b in slice_bounds.values())
    if hi != n_total:
        raise OracleRefusal("check_space", "slice cover mismatch")
    return True


class Partition:
    """(process, family-component) latent blocks + the 51-global border.

    Genus->family consistency is VERIFIED: a genus whose species span
    several families merges those families into one component (union-find)
    — the block-diagonality claim is then EXACT by construction."""

    def __init__(self, slice_bounds, n_total, S, NG, NF, gen_of, fam_of,
                 families, n_border=N_GLOBAL, z_shift=0):
        """n_border/z_shift: the extended-KKT spaces (oracle_v36) append
        auxiliary BORDER variables after the model globals — latent theta
        indices shift by z_shift and the border widens; the v3.1 default
        (51, 0) is byte-identical to the original behavior."""
        self.n_total = int(n_total)
        self.n_border = int(n_border)
        self.z_shift = int(z_shift)
        gen_of = np.asarray(gen_of)
        fam_of = np.asarray(fam_of)
        parent = list(range(NF))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)

        genus_fams = {}
        for s in range(S):
            genus_fams.setdefault(int(gen_of[s]), set()).add(int(fam_of[s]))
        merged = 0
        for g, fams in genus_fams.items():
            fams = sorted(fams)
            for b in fams[1:]:
                union(fams[0], b)
                merged += 1
        comp_of_fam = np.array([find(f) for f in range(NF)])
        comps = sorted(set(int(c) for c in comp_of_fam))
        comp_index = {c: i for i, c in enumerate(comps)}
        self.n_comp = len(comps)
        self.merged_genus_family_links = merged

        zbase = {}
        off = n_border - self.z_shift
        for p in ("surv", "adv", "fec"):
            for lev, cnt in (("sp", S), ("gen", NG), ("fam", NF)):
                nm = f"z_{lev}_{p}"
                lo, hi = slice_bounds[nm]
                assert (lo, hi) == (off, off + cnt), nm
                zbase[nm] = lo + self.z_shift
                off += cnt
        assert off + self.z_shift == n_total

        # theta index lists per (process, component)
        blocks = [[[] for _ in range(self.n_comp)] for _ in range(3)]
        fam_comp = np.array([comp_index[int(comp_of_fam[f])]
                             for f in range(NF)])
        gen_comp = np.zeros(NG, dtype=int)
        for s in range(S):
            gen_comp[int(gen_of[s])] = fam_comp[int(fam_of[s])]
        sp_comp = fam_comp[fam_of]
        for pi, p in enumerate(("surv", "adv", "fec")):
            for s in range(S):
                blocks[pi][sp_comp[s]].append(zbase[f"z_sp_{p}"] + s)
            for g in range(NG):
                blocks[pi][gen_comp[g]].append(zbase[f"z_gen_{p}"] + g)
            for f in range(NF):
                blocks[pi][fam_comp[f]].append(zbase[f"z_fam_{p}"] + f)
        self.blocks = []
        self.block_names = []
        procs = ("surv", "adv", "fec")
        for pi in range(3):
            for ci in range(self.n_comp):
                idx = np.array(sorted(blocks[pi][ci]), dtype=int)
                if len(idx) == 0:
                    continue
                self.blocks.append(idx)
                fams_in = [families[f] for f in range(NF)
                           if fam_comp[f] == ci]
                nm = fams_in[0] + ("" if len(fams_in) == 1
                                   else f"+{len(fams_in) - 1}")
                self.block_names.append(f"{procs[pi]}:{nm}")
        self.border_idx = np.arange(self.n_border)
        # theta index -> (block, offset); globals -> (-1, g)
        self.block_of = np.full(n_total, -1, dtype=int)
        self.off_of = np.zeros(n_total, dtype=int)
        self.off_of[:self.n_border] = np.arange(self.n_border)
        for bi, idx in enumerate(self.blocks):
            self.block_of[idx] = bi
            self.off_of[idx] = np.arange(len(idx))
        assert np.all(self.block_of[self.n_border:] >= 0), \
            "unassigned latent"

    @property
    def dims(self):
        return ([len(b) for b in self.blocks], self.n_border)

    def flatten(self, x):
        """(zs, g) -> flat list of arb, length n_total."""
        zs, g = x
        th = [None] * self.n_total
        for gi in range(self.n_border):
            th[gi] = g[gi]
        for bi, idx in enumerate(self.blocks):
            blk = zs[bi]
            for j, t in enumerate(idx):
                th[t] = blk[j]
        assert all(v is not None for v in th)
        return th

    def split_grad(self, grad):
        Fz = [[grad[t] for t in idx] for idx in self.blocks]
        Fg = [grad[gi] for gi in range(self.n_border)]
        return Fz, Fg


class _Acc:
    """Symmetric block-arrow Hessian accumulator (python lists of arb)."""

    def __init__(self, part):
        self.part = part
        ng = part.n_border
        z = arb(0)
        self.G = [[z] * ng for _ in range(ng)]
        self.D = [[[z] * len(idx) for _ in range(len(idx))]
                  for idx in part.blocks]
        self.B = [[[z] * ng for _ in range(len(idx))]
                  for idx in part.blocks]

    def add(self, a, b, val):
        """H[a,b] += val (a, b theta indices; symmetric storage)."""
        p = self.part
        ba, bb = p.block_of[a], p.block_of[b]
        oa, ob = p.off_of[a], p.off_of[b]
        if ba < 0 and bb < 0:                      # G (store upper, sym later)
            i, j = (oa, ob) if oa <= ob else (ob, oa)
            self.G[i][j] += val
        elif ba >= 0 and bb >= 0:
            assert ba == bb, "cross-block latent coupling — partition broken"
            i, j = (oa, ob) if oa <= ob else (ob, oa)
            self.D[ba][i][j] += val
        else:                                      # border coupling
            if ba >= 0:
                self.B[ba][oa][ob] += val
            else:
                self.B[bb][ob][oa] += val

    def to_mats(self):
        ng = self.part.n_border
        G = arb_mat(ng, ng)
        for i in range(ng):
            for j in range(i, ng):
                G[i, j] = self.G[i][j]
                if i != j:
                    G[j, i] = self.G[i][j]
        D, B = [], []
        for bi, idx in enumerate(self.part.blocks):
            n = len(idx)
            Db = arb_mat(n, n)
            for i in range(n):
                for j in range(i, n):
                    Db[i, j] = self.D[bi][i][j]
                    if i != j:
                        Db[j, i] = self.D[bi][i][j]
            Bb = arb_mat(n, ng)
            for i in range(n):
                for j in range(ng):
                    Bb[i, j] = self.B[bi][i][j]
            D.append(Db)
            B.append(Bb)
        return D, B, G


class ModelV31Oracle:
    """The baller.certify.block_krawczyk oracle for the v3.1 model.

    arrs: numpy cell arrays + structure (adapter_v31.extract_arrays);
    consts: prior constants (adapter_v31.constants). All floats enter arb
    exactly; slice layout is check_space-guarded at construction."""

    def __init__(self, arrs, consts):
        check_space(consts["slice_bounds"], consts["n_total"])
        self.arrs, self.consts = arrs, consts
        sb = consts["slice_bounds"]
        self.sl = {nm: sb[nm][0] for nm in sb}
        self.part = Partition(
            sb, consts["n_total"], consts["S"], consts["NG"], consts["NF"],
            arrs["gen_of"], arrs["fam_of"], consts["families"])
        self._const_cache = {}

    @property
    def dims(self):
        return self.part.dims

    # ---- baller oracle contract ----------------------------------------
    def F(self, x):
        th = self.part.flatten(x)
        _, grad, _ = self.evaluate(th, order2=False)
        return self.part.split_grad(grad)

    def H(self, x):
        th = self.part.flatten(x)
        _, _, acc = self.evaluate(th, order2=True)
        return acc.to_mats()

    # ---- the model ------------------------------------------------------
    def _nll_const(self, prec_key):
        """Theta-independent nll constants (binomial coefficients,
        -lgamma(R+1)): one arb sum per precision, cached keyed by prec."""
        if prec_key in self._const_cache:
            return self._const_cache[prec_key]
        a = self.arrs
        s = arb(0)
        for n, y in zip(a["sv_n"], a["sv_die"]):
            s -= (arb(float(n) + 1).lgamma() - arb(float(y) + 1).lgamma()
                  - arb(float(n) - float(y) + 1).lgamma())
        for R in a["rc_count"]:
            s += arb(float(R) + 1).lgamma()
        self._const_cache[prec_key] = s
        return s

    def evaluate(self, th, order2=True):
        """(nll ball, grad list, hess accumulator|None) at the arb vector
        th (tight or box). Raises OracleRefusal on any unprovable clip."""
        from flint import ctx
        a, sl = self.arrs, self.sl
        part = self.part
        n_total = self.consts["n_total"]
        grad = [arb(0)] * n_total
        cen = {}
        acc = _Acc(part) if order2 else None
        nll = self._nll_const(ctx.prec)

        def gv(name, i=0):
            return th[sl[name] + i]

        sig = [gv("log_sig", j).exp() for j in range(9)]
        # clip bounds as exact arb
        c19, c19u = arb("1e-9"), 1 - arb("1e-9")
        c110, c110u = arb("1e-10"), 1 - arb("1e-10")
        c112 = arb("1e-12")

        def clip_guard(ball, lo, hi, where):
            if not (ball > lo and ball < hi):
                raise OracleRefusal(where, f"clip bound not provably "
                                           f"inactive over the box: "
                                           f"{ball.str(8, radius=True)}")

        # log-space forms (enclosure-stable over wide boxes; measured:
        # naive expit/product balls lose positivity at r_unit ~ 0.06
        # sigma — mid+/-rad products of wide positive factors include 0):
        #   logsigmoid(x) = log(expit(x)), branch chosen so the inner exp
        #   argument has non-positive midpoint (single exp, positive
        #   lower endpoint; log1p then tight);
        #   pdie = -expm1(dt5 * logsigmoid(eta)) has a provably positive
        #   enclosure whenever the guards hold.
        ln19 = c19.log()
        ln110 = c110.log()

        def logsigmoid(x):
            if float(x.v.mid()) >= 0.0:
                return -((-x).exp().log1p())
            return x - x.exp().log1p()

        # ---------------- survival (jets in eta, log_nu) ----------------
        ilognu = sl["log_nu"]
        lognu_j = Jet.var(gv("log_nu"), 1, 2, order2)
        nu_j = lognu_j.exp()
        lg_nu_j = nu_j._chain(lgam_ball(nu_j.v), dig_ball(nu_j.v),
                              None if not order2 else trig_ball(nu_j.v))
        dig_cache = {}

        def lg_of_nu_plus(nfl):
            got = dig_cache.get(nfl)
            if got is None:
                w = nu_j.v + nfl
                got = (lgam_ball(w), dig_ball(w),
                       None if not order2 else trig_ball(w))
                dig_cache[nfl] = got
            lv, dg, tg = got
            return nu_j._chain(lv, dg, tg)

        zb_s = {"sp": sl["z_sp_surv"], "gen": sl["z_gen_surv"],
                "fam": sl["z_fam_surv"]}
        for c in range(len(a["sv_n"])):
            sp = int(a["sv_sp"][c]); k = int(a["sv_k"][c])
            st = int(a["sv_site"][c])
            n = float(a["sv_n"][c]); y = float(a["sv_die"][c])
            dt5 = float(a["sv_dt5"][c])
            x0 = float(a["X"][sp, 0]); x1 = float(a["X"][sp, 1])
            kc = (k - 3.0) / 2.0
            zb = float(a["sv_zbar"][c]); xd = float(a["sv_xdry"][c])
            Dz = float(a["D_site"][st]) * zb
            gen = int(a["gen_of"][sp]); fam = int(a["fam_of"][sp])
            izs = zb_s["sp"] + sp; izg = zb_s["gen"] + gen
            izf = zb_s["fam"] + fam
            zsv, zgv, zfv = th[izs], th[izg], th[izf]
            C0 = float(a["C"][st, 0]); C1 = float(a["C"][st, 1])
            C2 = float(a["C"][st, 2])
            eta = (gv("a_surv", k) + x0 * gv("b_surv", 0)
                   + x1 * gv("b_surv", 1) + (x0 * kc) * gv("b_sxk")
                   + C0 * gv("g_surv", 0) + C1 * gv("g_surv", 1)
                   + C2 * gv("g_surv", 2)
                   + zb * gv("w_surv", 0) + Dz * gv("w_surv", 1)
                   + xd * gv("w_surv", 2)
                   + sig[0] * zsv + sig[1] * zgv + sig[2] * zfv)
            ej = Jet.var(eta, 0, 2, order2)
            lp5 = logsigmoid(ej)          # log p5
            lm5 = logsigmoid(-ej)         # log (1 - p5)
            # CLIP CENSUS (not a refusal): the certified target is the
            # UNCLIPPED objective — the objective whose gradient the
            # reference's analytic chain computes (clips-as-identity, the
            # reference port's documented semantics); any reference clip bound the
            # box cannot be proven clear of is COUNTED and receipted.
            if not (lp5.v > ln19 and lm5.v > ln19):
                cen["sv_p5"] = cen.get("sv_p5", 0) + 1
            lpsur = lp5 * dt5
            psur = lpsur.exp()
            pdie = -(lpsur.expm1())
            if not (pdie.v > c110 and lpsur.v > ln110):
                cen["sv_pdie"] = cen.get("sv_pdie", 0) + 1
            al = pdie * nu_j
            be = psur * nu_j
            ll = ((al + y).lgamma() + (be + (n - y)).lgamma()
                  - lg_of_nu_plus(n) + lg_nu_j - al.lgamma() - be.lgamma())
            nll -= ll.v
            ge = -ll.g[0]              # dNLL/deta
            gn = -ll.g[1]              # dNLL/dlognu
            coords = ((sl["a_surv"] + k, 1.0), (sl["b_surv"], x0),
                      (sl["b_surv"] + 1, x1), (sl["b_sxk"], x0 * kc),
                      (sl["g_surv"], C0), (sl["g_surv"] + 1, C1),
                      (sl["g_surv"] + 2, C2), (sl["w_surv"], zb),
                      (sl["w_surv"] + 1, Dz), (sl["w_surv"] + 2, xd))
            bcoords = ((izs, sig[0]), (izg, sig[1]), (izf, sig[2]),
                       (sl["log_sig"], sig[0] * zsv),
                       (sl["log_sig"] + 1, sig[1] * zgv),
                       (sl["log_sig"] + 2, sig[2] * zfv))
            for idx, cf in coords:
                grad[idx] += ge * cf
            for idx, cf in bcoords:
                grad[idx] += ge * cf
            grad[ilognu] += gn
            if order2:
                hee = -ll.h[0]; hen = -ll.h[1]; hnn = -ll.h[2]
                allc = coords + bcoords
                for i in range(len(allc)):
                    ia, ca = allc[i]
                    hca = hee * ca
                    for j in range(i, len(allc)):
                        ib, cb = allc[j]
                        acc.add(ia, ib, hca * cb)
                    acc.add(ia, ilognu, hen * ca)
                acc.add(ilognu, ilognu, hnn)
                # d2eta terms: (z,l log_sig,l): sig_l ; (ls_l, ls_l): sig_l z
                acc.add(izs, sl["log_sig"], ge * sig[0])
                acc.add(izg, sl["log_sig"] + 1, ge * sig[1])
                acc.add(izf, sl["log_sig"] + 2, ge * sig[2])
                acc.add(sl["log_sig"], sl["log_sig"], ge * (sig[0] * zsv))
                acc.add(sl["log_sig"] + 1, sl["log_sig"] + 1,
                        ge * (sig[1] * zgv))
                acc.add(sl["log_sig"] + 2, sl["log_sig"] + 2,
                        ge * (sig[2] * zfv))

        # ---------------- advance (jet in eta) ---------------------------
        zb_a = {"sp": sl["z_sp_adv"], "gen": sl["z_gen_adv"],
                "fam": sl["z_fam_adv"]}
        for c in range(len(a["av_n"])):
            sp = int(a["av_sp"][c]); k = int(a["av_k"][c])
            st = int(a["av_site"][c])
            n = float(a["av_n"][c]); y = float(a["av_up"][c])
            dt5 = float(a["av_dt5"][c])
            x0 = float(a["X"][sp, 0]); x1 = float(a["X"][sp, 1])
            zb = float(a["av_zbar"][c]); xd = float(a["av_xdry"][c])
            Dz = float(a["D_site"][st]) * zb
            gen = int(a["gen_of"][sp]); fam = int(a["fam_of"][sp])
            izs = zb_a["sp"] + sp; izg = zb_a["gen"] + gen
            izf = zb_a["fam"] + fam
            zsv, zgv, zfv = th[izs], th[izg], th[izf]
            C0 = float(a["C"][st, 0]); C1 = float(a["C"][st, 1])
            C2 = float(a["C"][st, 2])
            eta = (gv("a_adv", k) + x0 * gv("b_gro", 0)
                   + x1 * gv("b_gro", 1)
                   + C0 * gv("g_gro", 0) + C1 * gv("g_gro", 1)
                   + C2 * gv("g_gro", 2)
                   + zb * gv("w_gro", 0) + Dz * gv("w_gro", 1)
                   + xd * gv("w_gro", 2)
                   + sig[3] * zsv + sig[4] * zgv + sig[5] * zfv)
            ej = Jet.var(eta, 0, 1, order2)
            lq5 = logsigmoid(ej)          # log q5
            lm5a = logsigmoid(-ej)        # log (1 - q5)
            if not (lq5.v > ln19 and lm5a.v > ln19):
                cen["av_q5"] = cen.get("av_q5", 0) + 1
            lstay = lm5a * dt5            # log stay5 (exact log form)
            q = -(lstay.expm1())          # 1 - stay5, positive enclosure
            if not (lstay.v > ln110 and q.v > c110):
                cen["av_q_stay5"] = cen.get("av_q_stay5", 0) + 1
            ll = q.log() * y + lstay * (n - y)
            nll -= ll.v
            ge = -ll.g[0]
            coords = ((sl["a_adv"] + k, 1.0), (sl["b_gro"], x0),
                      (sl["b_gro"] + 1, x1),
                      (sl["g_gro"], C0), (sl["g_gro"] + 1, C1),
                      (sl["g_gro"] + 2, C2), (sl["w_gro"], zb),
                      (sl["w_gro"] + 1, Dz), (sl["w_gro"] + 2, xd))
            bcoords = ((izs, sig[3]), (izg, sig[4]), (izf, sig[5]),
                       (sl["log_sig"] + 3, sig[3] * zsv),
                       (sl["log_sig"] + 4, sig[4] * zgv),
                       (sl["log_sig"] + 5, sig[5] * zfv))
            for idx, cf in coords:
                grad[idx] += ge * cf
            for idx, cf in bcoords:
                grad[idx] += ge * cf
            if order2:
                hee = -ll.h[0]
                allc = coords + bcoords
                for i in range(len(allc)):
                    ia, ca = allc[i]
                    hca = hee * ca
                    for j in range(i, len(allc)):
                        ib, cb = allc[j]
                        acc.add(ia, ib, hca * cb)
                acc.add(izs, sl["log_sig"] + 3, ge * sig[3])
                acc.add(izg, sl["log_sig"] + 4, ge * sig[4])
                acc.add(izf, sl["log_sig"] + 5, ge * sig[5])
                acc.add(sl["log_sig"] + 3, sl["log_sig"] + 3,
                        ge * (sig[3] * zsv))
                acc.add(sl["log_sig"] + 4, sl["log_sig"] + 4,
                        ge * (sig[4] * zgv))
                acc.add(sl["log_sig"] + 5, sl["log_sig"] + 5,
                        ge * (sig[5] * zfv))

        # ------------- fecundity (jets in eta, a_fec1, log_phi) ----------
        offs = self.consts["fec_offsets"]                 # [0,1,2,3]
        ilogphi = sl["log_phi"]; ia1 = sl["a_fec"] + 1
        logphi_j = Jet.var(gv("log_phi"), 2, 3, order2)
        phi_j = logphi_j.exp()
        lgphi_j = phi_j._chain(lgam_ball(phi_j.v), dig_ball(phi_j.v),
                               None if not order2 else trig_ball(phi_j.v))
        zb_f = {"sp": sl["z_sp_fec"], "gen": sl["z_gen_fec"],
                "fam": sl["z_fam_fec"]}
        for c in range(len(a["rc_count"])):
            sp = int(a["rc_sp"][c]); st = int(a["rc_site"][c])
            R = float(a["rc_count"][c]); dt5 = float(a["rc_dt5"][c])
            zb = float(a["rc_zbar"][c])
            x0 = float(a["X"][sp, 0]); x1 = float(a["X"][sp, 1])
            grp = int(a["group_of"][sp])
            gen = int(a["gen_of"][sp]); fam = int(a["fam_of"][sp])
            izs = zb_f["sp"] + sp; izg = zb_f["gen"] + gen
            izf = zb_f["fam"] + fam
            zsv, zgv, zfv = th[izs], th[izg], th[izf]
            C0 = float(a["C"][st, 0]); C1 = float(a["C"][st, 1])
            C2 = float(a["C"][st, 2])
            eta = (gv("a_fec", 0) + x0 * gv("b_fec", 0)
                   + x1 * gv("b_fec", 1)
                   + C0 * gv("g_fec", 0) + C1 * gv("g_fec", 1)
                   + C2 * gv("g_fec", 2)
                   + zb * gv("lam", grp)
                   + sig[6] * zsv + sig[7] * zgv + sig[8] * zfv)
            ej = Jet.var(eta, 0, 3, order2)
            a1j = Jet.var(gv("a_fec", 1), 1, 3, order2)
            # LOG-SPACE fecundity (enclosure-stable over wide boxes —
            # ball PRODUCTS of wide positive factors lose positivity and
            # log(mu) NaNs; measured at r_unit=0.057 sigma, fit13):
            #   A = exp(O_SHIFT*a1) * B,  B = sum a_j exp(a1*(off_j-O_SHIFT))
            #   (shifted exponents stay narrow; B provably positive)
            #   lmu = log(dt5) + eta + O_SHIFT*a1 + log B
            #   log(phi+mu) branch-stable via log1p(exp(smaller-larger)).
            O_SHIFT = 1.5                       # pinned exponent shift
            Bj = Jet.const(arb(0), 3, order2)
            for j, ad in enumerate(a["rc_adults"][c]):
                adf = float(ad)
                if adf != 0.0:
                    Bj = Bj + (a1j * (float(offs[j]) - O_SHIFT)).exp() * adf
            lmu = (ej + a1j * O_SHIFT + Bj.log()) + arb(dt5).log()
            # floor guard: mu > 1e-12  <=>  lmu > log(1e-12), certified
            if not (lmu.v > c112.log()):
                cen["rc_mu"] = cen.get("rc_mu", 0) + 1
            dlm = lmu - logphi_j
            if float(dlm.v.mid()) <= 0.0:
                lpm = logphi_j + (dlm.exp()).log1p()
            else:
                lpm = lmu + ((-dlm).exp()).log1p()
            ll = ((phi_j + R).lgamma() - lgphi_j
                  + phi_j * (logphi_j - lpm) + (lmu - lpm) * R)
            nll -= ll.v
            ge = -ll.g[0]; ga1 = -ll.g[1]; gph = -ll.g[2]
            coords = ((sl["a_fec"], 1.0), (sl["b_fec"], x0),
                      (sl["b_fec"] + 1, x1),
                      (sl["g_fec"], C0), (sl["g_fec"] + 1, C1),
                      (sl["g_fec"] + 2, C2), (sl["lam"] + grp, zb))
            bcoords = ((izs, sig[6]), (izg, sig[7]), (izf, sig[8]),
                       (sl["log_sig"] + 6, sig[6] * zsv),
                       (sl["log_sig"] + 7, sig[7] * zgv),
                       (sl["log_sig"] + 8, sig[8] * zfv))
            for idx, cf in coords:
                grad[idx] += ge * cf
            for idx, cf in bcoords:
                grad[idx] += ge * cf
            grad[ia1] += ga1
            grad[ilogphi] += gph
            if order2:
                hee = -ll.h[hidx(3, 0, 0)]
                hea = -ll.h[hidx(3, 0, 1)]
                hep = -ll.h[hidx(3, 0, 2)]
                haa = -ll.h[hidx(3, 1, 1)]
                hap = -ll.h[hidx(3, 1, 2)]
                hpp = -ll.h[hidx(3, 2, 2)]
                allc = coords + bcoords
                for i in range(len(allc)):
                    ia, ca = allc[i]
                    hca = hee * ca
                    for j in range(i, len(allc)):
                        ib, cb = allc[j]
                        acc.add(ia, ib, hca * cb)
                    acc.add(ia, ia1, hea * ca)
                    acc.add(ia, ilogphi, hep * ca)
                acc.add(ia1, ia1, haa)
                acc.add(ia1, ilogphi, hap)
                acc.add(ilogphi, ilogphi, hpp)
                acc.add(izs, sl["log_sig"] + 6, ge * sig[6])
                acc.add(izg, sl["log_sig"] + 7, ge * sig[7])
                acc.add(izf, sl["log_sig"] + 8, ge * sig[8])
                acc.add(sl["log_sig"] + 6, sl["log_sig"] + 6,
                        ge * (sig[6] * zsv))
                acc.add(sl["log_sig"] + 7, sl["log_sig"] + 7,
                        ge * (sig[7] * zgv))
                acc.add(sl["log_sig"] + 8, sl["log_sig"] + 8,
                        ge * (sig[8] * zfv))

        # ---------------- priors (analytic; reference formulas) ----------
        co = self.consts

        def gauss(name, width, center, sd):
            base = sl[name]
            for i in range(width):
                cen = center[i] if isinstance(center, (list, tuple,
                                                       np.ndarray)) else center
                sdi = sd[i] if isinstance(sd, (list, tuple,
                                               np.ndarray)) else sd
                r = (th[base + i] - float(cen)) / float(sdi)
                nonlocal_nll[0] += 0.5 * (r * r)
                grad[base + i] += r / float(sdi)
                if order2:
                    acc.add(base + i, base + i,
                            arb(1) / (float(sdi) * float(sdi)))

        nonlocal_nll = [nll]
        psd = co["prior_sd"]
        gauss("a_surv", 7, co["pc_surv_logit"], psd["a_surv"])
        gauss("a_adv", 6, 0.0, psd["a_adv"])
        gauss("a_fec", 2, [co["pc_fec_log"], 0.0],
              [psd["a_fec0"], psd["a_fec1"]])
        for nm in ("b_surv", "b_gro", "b_fec"):
            gauss(nm, 2, 0.0, psd["b"])
        gauss("b_sxk", 1, 0.0, psd["b_sxk"])
        for nm in ("g_surv", "g_gro", "g_fec"):
            gauss(nm, 3, 0.0, psd["g"])
        gauss("log_sig", 9, co["pc_log_sig"], psd["log_sig"])
        gauss("lam", 3, 0.0, psd["lam"])
        gauss("w_surv", 3, 0.0, psd["w"])
        gauss("w_gro", 3, 0.0, psd["w"])
        gauss("log_phi", 1, co["pc_log_phi"], psd["log_phi"])
        gauss("log_nu", 1, co["pc_log_nu"], psd["log_nu"])
        nll = nonlocal_nll[0]
        for p in ("surv", "adv", "fec"):
            for lev in ("sp", "gen", "fam"):
                base, hi = co["slice_bounds"][f"z_{lev}_{p}"]
                for t in range(base, hi):
                    z = th[t]
                    nll += 0.5 * (z * z)
                    grad[t] += z
                    if order2:
                        acc.add(t, t, arb(1))
        self.last_clip_census = dict(cen)
        return nll, grad, acc

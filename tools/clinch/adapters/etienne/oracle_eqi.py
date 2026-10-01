"""Equal-I multi-sample Etienne oracle (2-dim): u = (lntheta, lnI).

For equal I the multi-sample kernel factorizes M(D,A,I) = I^A Mtilde(A)
with Mtilde parameter-independent (built ONCE with I=1) — single-sample
structure with the prefactor I-terms summed over samples. Derivatives:
phase2_multisample_eqI.lnP_derivs_ms by identity. Same layout as
oracle_ss: border=[lntheta], one latent block [lnI]. NLL sign convention.
"""
from . import ident
from .common import Part, cur_prec, mat1

from clinch.oracle_v31 import OracleRefusal


class EtienneEqIOracle:
    def __init__(self, Dmat, kernel_prec=208, tag=""):
        p = ident.pilot()
        self._derivs = p["phase2_multisample_eqI"].lnP_derivs_ms
        self.Dmat = [tuple(int(x) for x in row) for row in Dmat]
        self.N = len(self.Dmat[0])
        self.S = len(self.Dmat)
        self.Js = [sum(r[i] for r in self.Dmat) for i in range(self.N)]
        self.tag = tag
        self.kernel_prec = int(kernel_prec)
        self.MT = p["phase2_multisample_eqI"].build_Mtilde(
            self.Dmat, self.kernel_prec)
        self.part = Part(blocks=[[1]], block_names=["local:lnI"],
                         border_idx=[0])
        self._centers = None

    @property
    def dims(self):
        return ([1], 1)

    def _eval(self, x):
        zs, g = x
        th, I = g[0].exp(), zs[0][0].exp()
        if not (th.is_finite() and I.is_finite()):
            raise OracleRefusal("domain", "exp(u) not finite over the box")
        return self._derivs(self.Dmat, self.MT, th, I, cur_prec(),
                            centers=self._centers)

    def nll(self, x):
        lnP, _, _ = self._eval(x)
        return -lnP

    def F(self, x):
        _, gr, _ = self._eval(x)
        return [[-gr[1]]], [-gr[0]]

    def H(self, x):
        _, _, Hm = self._eval(x)
        return [mat1(-Hm[1][1])], [mat1(-Hm[1][0])], mat1(-Hm[0][0])

#!/usr/bin/env python3
"""fault_harness.py — DECLARED PLANTED FAULT for the deform_transport battery.
Re-creates a real historical kit bug exactly as recorded:

  the kit's local recurrence evaluated the W_k polynomials at the
  TARGET index m where the engine letter evaluates at the SOURCE index
  m-k — i.e. the kit was marching the FORMAL-DUAL operator.

This wrapper converts the kit's source-index call sites into target-index
evaluations WITHOUT touching the kit's bytes: the kit calls
W_jets(k, m-k, J); the wrapper forwards to W_jets(k, (m-k)+k, J) =
W_jets(k, m, J) for k >= 1 (k = 0 agrees in both conventions).  theta_form is
dualized CONSISTENTLY (the historical bug built the theta form to the same
wrong convention, which is why gate GD did not catch it): the target-index
recurrence sum_k W_k(m) c_{m-k} = 0 is the operator sum_k x^k W_k(theta + k),
so the dual theta rows are c_i(x) = sum_k x^k [n^i] W_k(n+k).

The battery uses this harness for the planted-fault row ONLY; it is not part
of the instrument.  The instrument law it exercises: the target-index
convention computes the FORMAL-DUAL operator — invisible on the self-dual
rank-2 Legendre control, CAUGHT by the rank-6 Sym^5 control (the control-first
design's whole point)."""
import math


class TargetIndexFault:
    """Wraps a copied-kit OpData; delegates everything except the recurrence
    index convention (and the theta form, dualized to match)."""

    def __init__(self, op):
        self._op = op

    def __getattr__(self, attr):
        return getattr(self._op, attr)

    def W_jets(self, k, m, J):
        # kit call sites pass the SOURCE index m-k as `m`; forwarding with a
        # +k shift evaluates the jets at the TARGET index — the planted fault.
        if k:
            return self._op.W_jets(k, m + k, J)
        return self._op.W_jets(0, m, J)

    def theta_form(self):
        """Dual theta rows: c_i(x) = sum_k x^k [n^i] W_k(n+k) (the operator
        L~ = sum_k x^k W_k(theta+k) the faulted recurrence actually marches)."""
        op = self._op
        ci = [[0] * (op.dmax + 1) for _ in range(op.order + 1)]
        for k in range(op.dmax + 1):
            poly = op.Wk_poly[k]           # coeffs of W_k(n) in n
            # shift n -> n+k: [n^i] W_k(n+k) = sum_{t>=i} C(t,i) k^{t-i} poly[t]
            for i in range(min(len(poly), op.order + 1)):
                acc = 0
                for t in range(i, len(poly)):
                    acc += poly[t] * math.comb(t, i) * (k ** (t - i))
                ci[i][k] = acc
        return ci

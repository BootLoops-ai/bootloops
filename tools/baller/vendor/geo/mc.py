# baller engine — vectorized Monte Carlo with Clopper-Pearson certification.
"""Vectorized large-sample Monte Carlo for CSEP consistency-test quantiles,
with Clopper-Pearson certification. Reproduces the sampling semantics of
pyCSEP _poisson_likelihood_test (fixed-N multinomial via cumsum+searchsorted,
Poisson joint LL with log(n!) penalty) but batched for 10^6+ simulations.
"""
import numpy
from scipy.special import gammaln
from scipy.stats import beta as _beta


def sim_ll_fixed_n(log_rates, cum_weights, n_events, const_term, nsims, rng, batch=200000):
    # batch<=0 would spin an INFINITE loop — refused
    if batch <= 0:
        raise ValueError(f"batch must be > 0 (got {batch})")
    """Simulate joint LL for catalogs of fixed size n_events.

    log_rates: (nbins,) natural log of (possibly scaled) bin rates
    cum_weights: (nbins,) cumsum(rates)/sum(rates)  (same array pyCSEP uses)
    const_term: expected_forecast_count subtracted in the LL
    Returns array (nsims,) of LL values (float64).
    """
    out = numpy.empty(nsims)
    done = 0
    while done < nsims:
        b = min(batch, nsims - done)
        r = rng.random((b, n_events))
        pnts = numpy.searchsorted(cum_weights, r, side='right')
        # sum of log rates over events
        s = log_rates[pnts].sum(axis=1)
        # multiplicity penalty: sum over bins ln(m!) ; via sorted run positions
        ps = numpy.sort(pnts, axis=1)
        run = numpy.ones(b)
        pen = numpy.zeros(b)
        for j in range(1, n_events):
            eq = ps[:, j] == ps[:, j - 1]
            run = numpy.where(eq, run + 1, 1.0)
            pen += numpy.where(eq, numpy.log(run), 0.0)
        out[done:done + b] = s - pen - const_term
        done += b
    return out


def sim_ll_poisson_n(log_rates, cum_weights, lam, const_term, nsims, rng, batch=100000):
    # batch<=0 would spin an INFINITE loop — refused
    if batch <= 0:
        raise ValueError(f"batch must be > 0 (got {batch})")
    """Simulate joint LL with N ~ Poisson(lam) (pyCSEP unconditional L-test)."""
    out = numpy.empty(nsims)
    done = 0
    while done < nsims:
        b = min(batch, nsims - done)
        ns = rng.poisson(lam, size=b)
        res = numpy.empty(b)
        for n_val in numpy.unique(ns):
            sel = numpy.where(ns == n_val)[0]
            if n_val == 0:
                res[sel] = -const_term
                continue
            r = rng.random((len(sel), int(n_val)))
            pnts = numpy.searchsorted(cum_weights, r, side='right')
            s = log_rates[pnts].sum(axis=1)
            ps = numpy.sort(pnts, axis=1)
            run = numpy.ones(len(sel))
            pen = numpy.zeros(len(sel))
            for j in range(1, int(n_val)):
                eq = ps[:, j] == ps[:, j - 1]
                run = numpy.where(eq, run + 1, 1.0)
                pen += numpy.where(eq, numpy.log(run), 0.0)
            res[sel] = s - pen - const_term
        out[done:done + b] = res
        done += b
    return out


def quantile_ci(k, n, conf=0.99):
    # out-of-domain k would return silent (nan, nan) — refused
    if not (0 <= k <= n and n > 0):
        raise ValueError(f"quantile_ci requires 0 <= k <= n, n > 0 (got k={k}, n={n})")
    """Clopper-Pearson CI for a binomial proportion k/n."""
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else _beta.ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else _beta.ppf(1 - a, k + 1, n - k)
    return lo, hi


def obs_ll(log_rates, counts_idx, counts_val, const_term):
    s = (log_rates[counts_idx] * counts_val).sum()
    pen = gammaln(numpy.asarray(counts_val, dtype=float) + 1).sum()
    return s - pen - const_term

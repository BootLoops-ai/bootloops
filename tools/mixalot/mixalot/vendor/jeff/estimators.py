#!/usr/bin/env python3
"""estimators.py — MC marginal-likelihood estimator suite (practitioner defaults).

Model family (LSX conventions, matching pilot/closed_form_1var.py, verified by exact
rational arithmetic at build time):
  - family m1: one binary variable, cells x in {0,1}, per-observation pmf of a
    k-mixture of Bernoulli(theta_j).
  - family m4: four-flip binomial sufficient statistic, cells x in {0..4},
    per-observation pmf of a k-mixture of Binomial(4, theta_j)  (binomial
    coefficient C(4,x) included in the cell probability).
  Likelihood of the count vector u is the BARE product over cells
      L(params) = prod_x p_x(params)^{u_x}
  (iid-sequence likelihood; NO multinomial coefficient over orderings — this is
  the convention under which pilot/closed_form_1var.py gives the exact evidence).
  Priors: uniform on (0,1) for every parameter.  k=1: params=(theta,);
  k=2: params=(w, theta1, theta2), all iid U(0,1).
  Evidence: Z = integral of L over the unit cube.

Estimators (all at literature-standard practitioner defaults, documented in the
emitted settings dict; no tuning-to-win):
  hm      harmonic mean (Newton & Raftery 1994)         — expected to be terrible
  bridge  iterative bridge sampling (Meng & Wong 1996), moment-matched Gaussian
          proposal in logit space, no warp (bridgesampling package default)
  ss      stepping-stone, K=20 Beta(0.3,1)-quantile rungs (Xie et al. 2011)
  chib    Chib (1995) from data-augmentation Gibbs; RAW value emitted as
          logZ_hat (the practitioner-default failure mode), permutation-
          corrected value as second field.  Label-switching caveat: Neal (1999).
  nested  dynesty static NestedSampler, pure defaults (nlive=500)

Backend: self-contained JAX HMC / numpy Gibbs ("--backend jax", the only one
implemented here).  "--backend pymc" is a stub that will be slotted in on a larger host
without interface change.

CLI:
  python estimators.py --family m4|m1 --k 1|2 --counts u0,u1,... \
         --estimator hm|bridge|ss|chib|nested --seed S [--backend jax]
Output: single-line JSON
  {"logZ_hat": float, "err_est": float|null, "diagnostics": {...},
   "settings": {...}, "wall_s": float}
(chib additionally emits "logZ_perm_corrected").
Determinism: same seed => same answer (all RNGs seeded from --seed).
"""

import argparse
import json
import math
import os
import sys
import time

# Cap thread pools BEFORE numpy/jax import: hosts with many cores run probe-class
# python runs live under ulimit -v 31GB; a full-width XLA/Eigen/LLVM thread pool
# (96 threads x stacks + glibc arenas) flakily exhausts the address-space cap at
# JIT-compile time (observed: "LLVM ERROR: pthread_create failed"). 4 threads is
# plenty for these tiny models and keeps single runs polite on the shared box.
os.environ.setdefault("XLA_FLAGS",
                      "--xla_cpu_multi_thread_eigen=false "
                      "intra_op_parallelism_threads=4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MALLOC_ARENA_MAX", "2")

import numpy as np
from scipy.special import expit, logsumexp
from scipy.stats import beta as beta_dist
from scipy.optimize import minimize

# ----------------------------------------------------------------------------- model


class Model:
    """k-mixture of Binomial(n_flips, theta) on (n_flips+1)-cell counts."""

    def __init__(self, family, k, counts):
        if family == "m4":
            self.n_flips = 4
        elif family == "m1":
            self.n_flips = 1
        else:
            raise ValueError(f"unknown family {family!r}")
        self.family = family
        self.cells = self.n_flips + 1
        if len(counts) != self.cells:
            raise ValueError(f"family {family} needs {self.cells} counts, got {len(counts)}")
        if k not in (1, 2):
            raise ValueError("k must be 1 or 2")
        self.k = k
        self.dim = 1 if k == 1 else 3  # (theta,) or (w, theta1, theta2)
        self.counts = np.asarray(counts, dtype=np.float64)
        self.x = np.arange(self.cells, dtype=np.float64)
        self.logC = np.array([math.lgamma(self.n_flips + 1) - math.lgamma(i + 1)
                              - math.lgamma(self.n_flips - i + 1) for i in range(self.cells)])

    # ---- numpy log-likelihood, vectorized over rows of theta (m, dim) -> (m,)
    def loglik_np(self, params):
        params = np.atleast_2d(np.asarray(params, dtype=np.float64))
        t_all = np.clip(params, 1e-300, 1.0 - 1e-16)

        def comp_logpmf(t):  # t (m,1) -> (m, cells)
            with np.errstate(divide="ignore"):
                return (self.logC[None, :] + self.x[None, :] * np.log(t)
                        + (self.n_flips - self.x)[None, :] * np.log1p(-t))

        if self.k == 1:
            logp = comp_logpmf(t_all[:, 0:1])
        else:
            w = t_all[:, 0:1]
            with np.errstate(divide="ignore"):
                la = np.log(w) + comp_logpmf(t_all[:, 1:2])
                lb = np.log1p(-w) + comp_logpmf(t_all[:, 2:3])
            logp = np.logaddexp(la, lb)
        contrib = np.where(self.counts[None, :] > 0, self.counts[None, :] * logp, 0.0)
        return contrib.sum(axis=1)

    # ---- jax scalar log-likelihood (theta vector (dim,)) built lazily
    def loglik_jax(self):
        import jax.numpy as jnp
        x = jnp.asarray(self.x)
        logC = jnp.asarray(self.logC)
        u = jnp.asarray(self.counts)
        nf = float(self.n_flips)
        k = self.k

        def comp(t):
            return logC + x * jnp.log(t) + (nf - x) * jnp.log1p(-t)

        def ll(params):
            if k == 1:
                logp = comp(params[0])
            else:
                la = jnp.log(params[0]) + comp(params[1])
                lb = jnp.log1p(-params[0]) + comp(params[2])
                logp = jnp.logaddexp(la, lb)
            return jnp.sum(jnp.where(u > 0, u * logp, 0.0))

        return ll


def _logmeanexp(a, axis=None):
    a = np.asarray(a)
    n = a.size if axis is None else a.shape[axis]
    return logsumexp(a, axis=axis) - math.log(n)


def _ess(x):
    """Effective sample size of per-chain series x (chains, draws); initial
    positive-sequence truncation of the autocorrelation, summed over chains."""
    total = 0.0
    for c in np.atleast_2d(x):
        c = c - c.mean()
        n = len(c)
        v = c.var()
        if v <= 0:
            total += float(n)
            continue
        acf = np.correlate(c, c, "full")[n - 1:] / (v * n)
        s = 0.0
        for t in range(1, n):
            if acf[t] < 0:
                break
            s += acf[t]
        total += n / (1.0 + 2.0 * s)
    return float(total)


# ----------------------------------------------------------------------------- JAX HMC

HMC_DEFAULTS = dict(chains=4, leapfrog_steps=10, target_accept=0.8,
                    adaptation="dual averaging (Hoffman & Gelman 2014 Alg.5: "
                               "gamma=0.05, t0=10, kappa=0.75, eps0=0.1)")


def run_hmc(model, beta, seed, warmup, draws, chains=4, leapfrog=10, target=0.8):
    """Hand-rolled HMC in logit (unconstrained) space on the power posterior
    prior(theta) * L(theta)^beta.  Returns dict with theta draws, loglik values,
    z draws, acceptance rate.  Deterministic in seed (JAX CPU)."""
    import jax
    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp

    d = model.dim
    ll = model.loglik_jax()

    def logpost(z, b):
        t = jax.nn.sigmoid(z)
        # uniform prior on (0,1)^d => density 1; Jacobian of sigmoid:
        lj = -jnp.sum(jax.nn.softplus(z) + jax.nn.softplus(-z))
        return b * ll(t) + lj

    grad = jax.grad(logpost, argnums=0)
    T = warmup + draws
    mu, gamma, t0, kappa = math.log(10 * 0.1), 0.05, 10.0, 0.75

    def one_step(carry, t):
        z, key, leps, leps_bar, hbar, b = carry
        key, k1, k2 = jax.random.split(key, 3)
        eps = jnp.exp(jnp.where(t < warmup, leps, leps_bar))
        p0 = jax.random.normal(k1, (d,))
        lp0 = logpost(z, b)
        h0 = -lp0 + 0.5 * jnp.dot(p0, p0)

        def leap(i, st):
            zz, pp = st
            pp = pp + 0.5 * eps * grad(zz, b)
            zz = zz + eps * pp
            pp = pp + 0.5 * eps * grad(zz, b)
            return zz, pp

        z1, p1 = jax.lax.fori_loop(0, leapfrog, leap, (z, p0))
        lp1 = logpost(z1, b)
        h1 = -lp1 + 0.5 * jnp.dot(p1, p1)
        log_alpha = jnp.minimum(0.0, h0 - h1)
        log_alpha = jnp.where(jnp.isfinite(log_alpha), log_alpha, -jnp.inf)
        alpha = jnp.exp(log_alpha)
        accept = jnp.log(jax.random.uniform(k2)) < log_alpha
        z_new = jnp.where(accept, z1, z)
        # dual averaging (warmup only)
        m = t + 1.0
        hbar_new = (1.0 - 1.0 / (m + t0)) * hbar + (target - alpha) / (m + t0)
        leps_new = mu - jnp.sqrt(m) / gamma * hbar_new
        eta = m ** (-kappa)
        leps_bar_new = eta * leps_new + (1.0 - eta) * leps_bar
        upd = t < warmup
        leps = jnp.where(upd, leps_new, leps)
        leps_bar = jnp.where(upd, leps_bar_new, leps_bar)
        hbar = jnp.where(upd, hbar_new, hbar)
        return (z_new, key, leps, leps_bar, hbar, b), (z_new, alpha)

    def run_chain(key, z0, b):
        carry0 = (z0, key, jnp.log(0.1), jnp.log(0.1), 0.0, b)
        _, (zs, alphas) = jax.lax.scan(one_step, carry0, jnp.arange(T, dtype=jnp.float64))
        return zs, alphas

    runner = jax.jit(jax.vmap(run_chain, in_axes=(0, 0, None)))
    key = jax.random.PRNGKey(seed)
    kinit, krun = jax.random.split(key)
    theta0 = jax.random.uniform(kinit, (chains, d), minval=0.05, maxval=0.95)
    z0 = jnp.log(theta0) - jnp.log1p(-theta0)
    chain_keys = jax.random.split(krun, chains)
    zs, alphas = runner(chain_keys, z0, float(beta))
    zs = np.asarray(zs)[:, warmup:, :]                   # (chains, draws, d)
    accept = float(np.asarray(alphas)[:, warmup:].mean())
    theta = expit(zs)
    llk = model.loglik_np(theta.reshape(-1, d)).reshape(chains, draws)
    return dict(z=zs, theta=theta, llk=llk, accept=accept)


# ----------------------------------------------------------------------------- estimators


def est_hm(model, seed):
    s = dict(HMC_DEFAULTS)
    s.update(warmup=1000, draws=1000,
             estimator="harmonic mean of the likelihood over posterior draws",
             citation="Newton & Raftery (1994), JRSS-B 56:3-48 (eq. 15). Known "
                      "pathology: the estimator has infinite variance whenever the "
                      "prior is heavier-tailed than the posterior; expected to be "
                      "biased upward. Emitted as-is: that is data, not a bug.",
             err_note="err_est = std over the 4 per-chain harmonic-mean estimates "
                      "/ sqrt(4); unreliable for the same reason the estimator is.")
    r = run_hmc(model, 1.0, seed, 1000, 1000)
    llk = r["llk"]
    logZ = float(-_logmeanexp(-llk.ravel()))
    per_chain = -_logmeanexp(-llk, axis=1)
    err = float(np.std(per_chain, ddof=1) / math.sqrt(llk.shape[0]))
    diag = dict(acceptance=r["accept"], ess_loglik=_ess(llk),
                n_draws=int(llk.size), per_chain_logZ=[float(v) for v in per_chain])
    return logZ, err, diag, s


def _mvn_logpdf(z, mean, cov_chol):
    from scipy.linalg import solve_triangular
    d = z.shape[1]
    y = solve_triangular(cov_chol, (z - mean).T, lower=True)
    return (-0.5 * np.sum(y * y, axis=0) - 0.5 * d * math.log(2 * math.pi)
            - np.sum(np.log(np.diag(cov_chol))))


def est_bridge(model, seed):
    s = dict(HMC_DEFAULTS)
    s.update(warmup=1000, draws=1000,
             estimator="iterative (optimal) bridge sampling",
             proposal="moment-matched multivariate Gaussian in logit space, warp: "
                      "none (method='normal' default of the R bridgesampling package)",
             split="first half of each chain fits the proposal, second half enters "
                   "the bridge (bridgesampling default)",
             tol=1e-8, max_iter=1000,
             citation="Meng & Wong (1996), Statist. Sinica 6:831-860; Gronau et al. "
                      "(2017), J. Math. Psych. 81:80-97 (bridgesampling defaults).",
             err_note="err_est = sqrt(RE^2) of Fruhwirth-Schnatter (2004) as in "
                      "bridgesampling, with the autocorrelation term estimated by "
                      "batch means over the posterior half-chains.")
    r = run_hmc(model, 1.0, seed, 1000, 1000)
    z = r["z"]                                  # (chains, 1000, d)
    chains, draws, d = z.shape
    half = draws // 2
    z_fit = z[:, :half, :].reshape(-1, d)
    z1 = z[:, half:, :]                         # bridge samples, keep chain structure
    n1 = chains * half
    mean = z_fit.mean(axis=0)
    cov = np.cov(z_fit.T).reshape(d, d)
    chol = np.linalg.cholesky(cov + 1e-12 * np.eye(d))
    rng = np.random.default_rng(seed + 10**6)
    n2 = n1
    z2 = mean + (chol @ rng.standard_normal((d, n2))).T

    def log_q1(zz):  # unnormalized posterior in z-space
        t = expit(zz)
        lj = -np.sum(np.logaddexp(0, zz) + np.logaddexp(0, -zz), axis=1)
        return model.loglik_np(t) + lj

    l1 = (log_q1(z1.reshape(-1, d)) - _mvn_logpdf(z1.reshape(-1, d), mean, chol))
    l2 = (log_q1(z2) - _mvn_logpdf(z2, mean, chol))
    lstar = float(np.median(l1))
    s1, s2 = n1 / (n1 + n2), n2 / (n1 + n2)
    r_t = 1.0  # in units of exp(lstar)
    e1 = np.exp(l1 - lstar)
    e2 = np.exp(l2 - lstar)
    n_iter = 0
    for n_iter in range(1, s["max_iter"] + 1):
        num = np.mean(e2 / (s1 * e2 + s2 * r_t))
        den = np.mean(1.0 / (s1 * e1 + s2 * r_t))
        r_new = num / den
        if abs(r_new - r_t) / r_t < s["tol"]:
            r_t = r_new
            break
        r_t = r_new
    logZ = float(math.log(r_t) + lstar)
    # RE^2 error estimate
    f1 = 1.0 / (s1 * np.exp(l1 - logZ) + s2)            # at posterior samples
    f2 = np.exp(l2 - logZ) / (s1 * np.exp(l2 - logZ) + s2)
    f1c = f1.reshape(chains, half)
    nb = 20
    bs = half // nb
    bm = f1c[:, :nb * bs].reshape(chains, nb, bs).mean(axis=2).ravel()
    var_mean_f1 = np.var(bm, ddof=1) / bm.size
    re2 = var_mean_f1 / np.mean(f1) ** 2 + np.var(f2, ddof=1) / (n2 * np.mean(f2) ** 2)
    err = float(math.sqrt(max(re2, 0.0)))
    diag = dict(acceptance=r["accept"], ess_loglik=_ess(r["llk"]),
                iterations=n_iter, n1=n1, n2=n2, converged=bool(n_iter < s["max_iter"]))
    return logZ, err, diag, s


def est_ss(model, seed):
    K = 20
    s = dict(HMC_DEFAULTS)
    s.update(estimator="stepping-stone sampling",
             K=K, beta_spacing="beta_k = (k/K)^(1/0.3), k=0..K, i.e. quantiles of "
                               "Beta(0.3,1) (Xie et al. default, as in MrBayes/Phycas)",
             per_rung="4 chains x (500 warmup + 500 draws) HMC at each of the K "
                      "sampling rungs beta_0..beta_{K-1}; beta_0=0 rung sampled iid "
                      "from the prior (exact)",
             citation="Xie, Lewis, Fan, Kuo & Chen (2011), Syst. Biol. 60:150-160.",
             err_note="err_est = sqrt(sum over rungs of var(per-chain log ratio "
                      "estimates)/4); simple between-chain variance, no "
                      "autocorrelation correction beyond chain independence.")
    betas = (np.arange(K + 1) / K) ** (1.0 / 0.3)
    chains, warm, draws = 4, 500, 500
    logZ = 0.0
    var_sum = 0.0
    acc, ess_min = [], np.inf
    rng = np.random.default_rng(seed + 2 * 10**6)
    for k in range(K):
        dlt = betas[k + 1] - betas[k]
        if betas[k] == 0.0:
            theta = rng.uniform(size=(chains, draws, model.dim))
            llk = model.loglik_np(theta.reshape(-1, model.dim)).reshape(chains, draws)
        else:
            r = run_hmc(model, betas[k], seed + 1000 + k, warm, draws,
                        chains=chains)
            llk = r["llk"]
            acc.append(r["accept"])
            ess_min = min(ess_min, _ess(llk))
        logr_pooled = float(_logmeanexp(dlt * llk.ravel()))
        logr_chain = _logmeanexp(dlt * llk, axis=1)
        logZ += logr_pooled
        var_sum += float(np.var(logr_chain, ddof=1)) / chains
    err = float(math.sqrt(var_sum))
    diag = dict(n_rungs=K, betas_first3=[float(b) for b in betas[:3]],
                acceptance_mean=float(np.mean(acc)), ess_loglik_min=float(ess_min),
                draws_per_rung=chains * draws)
    return float(logZ), err, diag, s


def _gibbs_chain(model, seed, warmup, draws):
    """Data-augmentation Gibbs for the k=2 mixture on cell counts.
    State: allocations a_x (count of cell-x observations in component 1), w, t1, t2.
    Returns per-draw conditional Beta parameters (for the Chib ordinate) and states."""
    rng = np.random.default_rng(seed)
    u = model.counts.astype(np.int64)
    x = model.x
    nf = model.n_flips
    w, t1, t2 = rng.uniform(), rng.uniform(), rng.uniform()
    N = int(u.sum())
    S = float((x * u).sum())
    out = []
    for it in range(warmup + draws):
        f1 = np.exp(model.logC + x * np.log(t1 + 1e-300) + (nf - x) * np.log1p(-min(t1, 1 - 1e-16)))
        f2 = np.exp(model.logC + x * np.log(t2 + 1e-300) + (nf - x) * np.log1p(-min(t2, 1 - 1e-16)))
        pi = w * f1 / (w * f1 + (1 - w) * f2 + 1e-300)
        a = rng.binomial(u, pi)
        n1 = int(a.sum()); n2 = N - n1
        S1 = float((x * a).sum()); S2 = S - S1
        aw, bw = 1 + n1, 1 + n2
        a1, b1 = 1 + S1, 1 + nf * n1 - S1
        a2, b2 = 1 + S2, 1 + nf * n2 - S2
        w = rng.beta(aw, bw)
        t1 = rng.beta(a1, b1)
        t2 = rng.beta(a2, b2)
        if it >= warmup:
            out.append((aw, bw, a1, b1, a2, b2, w, t1, t2))
    return np.array(out)


def est_chib(model, seed):
    s = dict(estimator="Chib (1995) marginal likelihood from the Gibbs output",
             gibbs="4 independent chains x (1000 warmup + 1000 draws) data-"
                   "augmentation Gibbs (allocation counts per cell; conjugate "
                   "Beta full conditionals, uniform=Beta(1,1) priors)",
             eval_point="posterior mode (flat prior => MLE); Nelder-Mead in logit "
                        "space started from the best Gibbs draw",
             citation="Chib (1995), JASA 90:1313-1321.",
             label_switching_caveat="Neal (1999, 'Erroneous results in Marginal "
                 "likelihood from the Gibbs output'): for mixtures, the Rao-"
                 "Blackwellized posterior-ordinate estimate is biased whenever the "
                 "Gibbs chain fails to visit all k! labelings of the components "
                 "(typically overstating the ordinate by ~k!, i.e. logZ low by "
                 "~log k!). logZ_hat is the RAW UNCORRECTED value — the "
                 "practitioner-default failure mode being graded. "
                 "logZ_perm_corrected averages the conditional ordinate over all "
                 "k! label permutations (Berkhof, van Mechelen & Gelman 2003) and "
                 "is emitted as a second field.",
             err_note="err_est = std of per-chain raw logZ estimates / sqrt(4) "
                      "(k=2); exact conjugate case k=1 has err_est=0.")
    if model.k == 1:
        # conjugate: posterior is exactly Beta(1+S, 1+nf*N-S); Chib identity is exact
        u = model.counts
        N, S = float(u.sum()), float((model.x * u).sum())
        a, b = 1 + S, 1 + model.n_flips * N - S
        res = minimize(lambda z: -float(model.loglik_np(expit(z[None, :]))[0]),
                       x0=np.array([math.log((S + 0.5) / (model.n_flips * N - S + 0.5))]),
                       method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-12))
        tstar = expit(res.x)[None, :]
        llk_star = float(model.loglik_np(tstar)[0])
        ordinate = float(beta_dist.logpdf(tstar[0, 0], a, b))
        logZ = llk_star - ordinate  # log prior = 0
        diag = dict(theta_star=[float(tstar[0, 0])], exact_conjugate=True,
                    loglik_at_mode=llk_star)
        return logZ, 0.0, diag, s, logZ

    chains, warm, draws = 4, 1000, 1000
    runs = [_gibbs_chain(model, seed + 3 * 10**6 + c, warm, draws) for c in range(chains)]
    allp = np.concatenate(runs, axis=0)   # columns: aw,bw,a1,b1,a2,b2,w,t1,t2
    states = allp[:, 6:9]
    # evaluation point: posterior mode (flat prior => MLE), started from best draw
    llk_states = model.loglik_np(states)
    z0 = np.clip(states[int(np.argmax(llk_states))], 1e-6, 1 - 1e-6)
    z0 = np.log(z0) - np.log1p(-z0)
    res = minimize(lambda z: -float(model.loglik_np(expit(z[None, :]))[0]), x0=z0,
                   method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-12, maxiter=5000))
    tstar = expit(res.x)
    llk_star = float(model.loglik_np(tstar[None, :])[0])
    wst, t1st, t2st = tstar

    def ordinates(p):
        aw, bw, a1, b1, a2, b2 = (p[:, i] for i in range(6))
        lo_id = (beta_dist.logpdf(wst, aw, bw) + beta_dist.logpdf(t1st, a1, b1)
                 + beta_dist.logpdf(t2st, a2, b2))
        lo_sw = (beta_dist.logpdf(wst, bw, aw) + beta_dist.logpdf(t1st, a2, b2)
                 + beta_dist.logpdf(t2st, a1, b1))
        return lo_id, np.logaddexp(lo_id, lo_sw) - math.log(2.0)

    lo_id, lo_corr = ordinates(allp)
    logZ_raw = llk_star - float(_logmeanexp(lo_id))
    logZ_corr = llk_star - float(_logmeanexp(lo_corr))
    per_chain = [llk_star - float(_logmeanexp(ordinates(r)[0])) for r in runs]
    err = float(np.std(per_chain, ddof=1) / math.sqrt(chains))
    diag = dict(theta_star=[float(v) for v in tstar], loglik_at_mode=llk_star,
                per_chain_logZ_raw=per_chain, n_draws=int(allp.shape[0]),
                mean_alloc_frac=float((allp[:, 0] - 1).mean() / model.counts.sum()))
    return float(logZ_raw), err, diag, s, float(logZ_corr)


def est_nested(model, seed):
    import dynesty
    s = dict(estimator="static nested sampling (dynesty defaults)",
             sampler="dynesty.NestedSampler, nlive=500, all other constructor and "
                     "run_nested arguments at dynesty 3.0.0 defaults (default "
                     "bound='multi', default dlogz with add_live)",
             nlive=500,
             citation="Skilling (2006), Bayesian Anal. 1:833-859; Speagle (2020), "
                      "MNRAS 493:3132-3158 (dynesty).",
             err_note="err_est = dynesty's own logzerr (information-based).")

    def loglike(theta):
        return float(model.loglik_np(theta[None, :])[0])

    def ptform(u):
        return u

    rstate = np.random.default_rng(seed)
    sampler = dynesty.NestedSampler(loglike, ptform, model.dim, nlive=500,
                                    rstate=rstate)
    sampler.run_nested(print_progress=False)
    res = sampler.results
    logZ = float(res.logz[-1])
    err = float(res.logzerr[-1])
    diag = dict(niter=int(res.niter), ncall=int(np.sum(res.ncall)),
                eff=float(res.eff))
    return logZ, err, diag, s


# ----------------------------------------------------------------------------- main


def run(family, k, counts, estimator, seed, backend="jax"):
    if backend != "jax":
        raise NotImplementedError(
            f"backend {backend!r} not available on this host; only the self-"
            f"contained 'jax' backend is implemented (pymc backend slots in on "
            f"a larger host with the same interface).")
    model = Model(family, k, counts)
    t0 = time.perf_counter()
    extra = {}
    if estimator == "hm":
        logZ, err, diag, s = est_hm(model, seed)
    elif estimator == "bridge":
        logZ, err, diag, s = est_bridge(model, seed)
    elif estimator == "ss":
        logZ, err, diag, s = est_ss(model, seed)
    elif estimator == "chib":
        logZ, err, diag, s, logZ_corr = est_chib(model, seed)
        extra["logZ_perm_corrected"] = logZ_corr
    elif estimator == "nested":
        logZ, err, diag, s = est_nested(model, seed)
    else:
        raise ValueError(f"unknown estimator {estimator!r}")
    wall = time.perf_counter() - t0
    s.update(family=family, k=k, counts=[int(c) for c in counts], seed=seed,
             backend=backend, prior="uniform U(0,1) on every parameter",
             likelihood="bare product over cells prod_x p_x^{u_x} (iid-sequence "
                        "likelihood, no multinomial coefficient; matches "
                        "pilot/closed_form_1var.py conventions)")
    out = dict(logZ_hat=logZ, err_est=err, diagnostics=diag, settings=s,
               wall_s=round(wall, 3))
    out.update(extra)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--family", required=True, choices=["m4", "m1"])
    ap.add_argument("--k", required=True, type=int, choices=[1, 2])
    ap.add_argument("--counts", required=True,
                    help="comma-separated cell counts u0,...  (5 cells for m4, 2 for m1)")
    ap.add_argument("--estimator", required=True,
                    choices=["hm", "bridge", "ss", "chib", "nested"])
    seedgrp = ap.add_mutually_exclusive_group(required=True)
    seedgrp.add_argument("--seed", type=int)
    seedgrp.add_argument("--seeds",
                         help="comma-separated seed list; runs all in ONE process "
                              "(amortizes import+jit) and emits a batch JSON")
    ap.add_argument("--backend", default="jax",
                    help="sampling backend; 'jax' (implemented) or 'pymc' (stub)")
    args = ap.parse_args()
    counts = [int(v) for v in args.counts.split(",")]
    if args.seeds:
        seeds = [int(s) for s in args.seeds.split(",")]
        results = [run(args.family, args.k, counts, args.estimator, s, args.backend)
                   for s in seeds]
        print(json.dumps({"batch_seeds": seeds, "results": results}))
    else:
        out = run(args.family, args.k, counts, args.estimator, args.seed, args.backend)
        print(json.dumps(out))


if __name__ == "__main__":
    main()

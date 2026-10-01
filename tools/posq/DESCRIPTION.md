# posq — tool description

```latex
\subsection{\texttt{posq}: certified evidence integrals by exact positive
quadrature}
\texttt{posq} computes a Bayesian evidence integral as a certified positive
value-sum --- a degree-matched exact Gauss sum on the raw likelihood product,
plus log-space patch enclosures for what survives --- never expanding
coefficients, never subtracting, so two-sided machine-width certificates come
from positivity alone. The object class is likelihoods of continuous-time Markov
chains on small graphs with clock priors: after the substitution
$u_k=e^{-\beta s_k}$, exponential prior increments marginalize exactly into
weights $u^{c-1}$ with $c$ rational in the remaining parameters, and the
integrand $F$ is the raw product of pattern probabilities --- positive on the
domain and polynomial of known multidegree. A Gauss--Jacobi tensor rule built at
that degree (interval-Newton node certification at 2048 bits; receipts require
strictly positive weight balls, exact-moment reproduction, and monomial
exactness) integrates $F$ \emph{exactly}: the clock dimensions contribute zero
quadrature error, and because positive weights multiply positive values the
streamed ball sum is cancellation-free, giving certified lower and upper
endpoints at machine width (measured relative width $2^{-179}$ at 192-bit
working precision on the production 48.7-million-node object). The single
remaining dimension is integrated by a non-adaptive Cauchy--ellipse register
(certified Gauss--Legendre nodes per panel, an exactly integrable
Bernstein-basis absolute majorant for the remainder, and exact endpoint bands).
On the production four-taxon data set this yielded, first, one complete
two-sided evidence row of width $0.063$ nats at $27.4$ CPU-h, and then a full
fifteen-topology certified table (all rows two-sided, widths
$0.0095$--$1.16$ nats, $93$ CPU-h), from which certified decision sentences
follow by the winner-lower-minus-rival-upper rule: a maximum-evidence topology
with $\ln\mathrm{BF}\ge 0.47$ over the runner-up and $\ge 5.2$ over the losing
split classes, agreeing with an independently computed model's ranking.
Four usage rules are binding. The exact clock-marginalization is
quartet-native: four-leaf shapes only, with balanced shapes through a separately
gated min/diff decomposition (exactness gates on degrees, parity, and per-panel
pole caps must pass before any sweep). The one engineered error dimension is the
$p$-layer, whose register menu is priced by measurement, not by convergence
constants --- the shipped ellipse register meets a $0.1$-nat width target at
$27$--$38$ CPU-h per topology and its panel counts are set by the absolute
majorant, which pays a measured factor-of-several penalty over signed
sensitivity. The assembly must remain sign-positive: the only mixed-sign step is
the per-node value $P_y=A_y+B_yv$, every $P$ and every assembled $Z$ ball must
certify strictly positive, per-sweep relative width worse than $2^{-100}$ is a
hard failure, and boundary-clustered nodes are screened explicitly. Coefficient
expansion is prohibited --- the two failure mechanisms this tool exists to avoid
are expansion cost and expansion cancellation. The packaged engine was
additionally run, as a consumer adaptation, through the shipped \texttt{eras}
adversarial battery (closed-box containment with exact corner, face, and
near-face band sampling; planted narrowing, shifting, and dropped-dimension
corruptions; a precision probe; finite-difference cross-checks of an
exact-calculus derivative path), passing at three seeds with all plants firing
and zero violations.
Positioning against the nearest lines of work --- auto-validating interval
methods in phylogenetics, guaranteed inference for probabilistic programs,
validated quadrature (used here as substrate, not claimed), and the
robust-Bayes lineage for bounds on evidence ratios --- is given with full
references on the tool page, https://www.bootloops.ai/tools/posq.html.
```

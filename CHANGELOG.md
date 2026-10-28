# Changelog

## Unreleased (based on 2.7-6)

### Exact binary latent draws

- In `src/jomo1C.c`, draw observed binary latent values from their full
  one-sided conditional normal distribution. Category 1 requires a positive
  value; category 2 requires a negative value. The former rejection sampler
  also imposed bounds of -3 and 4 and stopped after 10,000 attempts, changing
  the probit augmentation target and potentially retaining a stale value.
- Use R's normal CDF and quantile functions with log-tail probabilities, so
  even observations far into the conditional tail can be sampled without
  probability underflow. The conditional mean and variance are unchanged.
- Add a one-sweep regression test of signs, conditional distributions, draws
  beyond the former bounds, and conditional means 40 standard deviations
  from the threshold.

### Reject invalid covariance proposals without retrying

- In the categorical covariance branch of `src/jomo1C.c`, make one normal
  proposal per free covariance entry. If the candidate matrix is not positive
  definite, restore the entry and skip its likelihood evaluation. Previously,
  proposals were redrawn up to 100 times, and an invalid final proposal could
  reach factorization and inversion.
- For binary correlations the unit diagonal makes the proposal scale constant,
  so the single normal proposal is symmetric. Retrying until feasible instead
  conditions it on a feasible interval `(l, h)`. Its normalizing probability
  `A(r) = Phi((h-r)/s) - Phi((l-r)/s)` depends on the current correlation;
  the former likelihood-only Metropolis ratio omitted `A(r) / A(r')`.
  Immediate rejection avoids that conditioning and preserves the intended
  target without an additional proposal-ratio term for binary correlations.
- Add a regression comparison with the analytic correlation posterior for
  two concordant binary observations, together with its reflected counterpart
  after relabeling one response. Check posterior means, boundary-tail mass,
  symmetry, fixed unit diagonals, and positive definiteness.
- The categorical covariance acceptance ratio uses the augmented normal
  likelihood. `l1cov.prior` enters only the continuous-only covariance branch;
  it does not impose an inverse-Wishart or Jeffreys prior on the binary
  correlations. The binary model's intended priors remain flat in the
  regression coefficients and the free correlations (LKJ(1)).

### Scope

These changes originate in the binary draw-sampler audit in rbmi's
`misc/design_binary_extension.qmd`. They live in the shared single-level
`jomo1C` backend used by `jomo1cat`, `jomo1mix`, and their `.MCMCchain`
variants (including dispatch through `jomo1`/`jomo`). Binary latent updates
also apply to binary columns in mixed-outcome models. The latent sampler for
outcomes with more than two categories, the continuous-only covariance
branch, and the other C backends are unchanged. This is not a general repair
or validation of all mixed-outcome or multilevel samplers. Public arguments
and intended priors are unchanged; random-number sequences will change.
The single-attempt covariance rule applies throughout `jomo1C`'s categorical
branch, including mixed outcomes. It does not resolve other issues such as
state-dependent proposal scales for free continuous variances in that branch.

### Validation

The regression scripts in `tests/` use base R and run during `R CMD check`.
After installing this source, individual scripts can also be run with
`Rscript tests/<script>.R`.

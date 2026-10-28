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

### Validation

The regression scripts in `tests/` use base R and run during `R CMD check`.
After installing this source, individual scripts can also be run with
`Rscript tests/<script>.R`.

library(jomo)

# One sweep leaves the observed latent values at the draws made using the
# supplied beta/covariance/start.imp, before beta and omega are updated.
# The first column's conditional distribution is therefore known exactly.
set.seed(20261028)
n <- 4000L
category <- rep(1:2, length.out = n)
y <- data.frame(a = factor(category), b = factor(rev(category)))
rho <- 0.7
for (mu in c(0, -40, 40)) {
  start <- cbind(rep(0, n), rep(c(-1, 1), each = n / 2))
  # The C backend can modify start.imp in place; keep the expected
  # conditional parameters separately before entering the sampler.
  conditional_mean <- mu + rho * start[, 2]
  conditional_sd <- sqrt(1 - rho^2)
  fit <- jomo1cat.MCMCchain(
    y, Y.numcat = c(2L, 2L), beta.start = matrix(c(mu, 0), 1L),
    l1cov.start = matrix(c(1, rho, rho, 1), 2L), start.imp = start,
    nburn = 1L, output = 0
  )
  z <- fit$finimp.latnorm
  stopifnot(all(is.finite(z)), all(z[category == 1, 1] > 0),
            all(z[category == 2, 1] < 0),
            all(z[rev(category) == 1, 2] > 0),
            all(z[rev(category) == 2, 2] < 0))
  for (cat in 1:2) {
    rows <- category == cat
    lower <- cat == 2L
    # Conditional tail probability ratios must be uniform, including when
    # ordinary (non-log) normal tail probabilities would underflow.
    u <- exp(pnorm(z[rows, 1], conditional_mean[rows], conditional_sd,
                   lower.tail = lower, log.p = TRUE) -
             pnorm(0, conditional_mean[rows], conditional_sd,
                   lower.tail = lower, log.p = TRUE))
    stopifnot(all(is.finite(u)), all(u > 0 & u < 1),
              abs(mean(u) - 0.5) < 0.025,
              abs(mean(u^2) - 1/3) < 0.025)
  }
  if (mu == -40) stopifnot(any(z[, 1] < -3))
  if (mu == 40) stopifnot(any(z[, 1] > 4))
}

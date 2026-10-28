library(jomo)

# Two concordant observations (1,1) and (2,2) give a known marginal posterior
# for rho under flat beta and rho priors. Integrating out both intercepts
# gives E[max(D1,0) * max(D2,0)], where D ~ N(0, 2 * Omega). Thus the
# unnormalized density is
#   (sqrt(1-rho^2) + rho * (pi/2 + asin(rho))) / pi,
# its integral is 3/4, and E[rho] = 4/9. This target puts appreciable mass
# near the positive-definite boundary, where retry-until-feasible is biased.
density <- function(r) {
  (sqrt(1 - r^2) + r * (pi/2 + asin(r))) / pi
}
tail_probability <- integrate(density, 0.9, 1)$value / (3/4)

for (direction in c(1, -1)) {
  set.seed(20261028)
  y <- data.frame(a = factor(1:2),
                  b = factor(if (direction == 1) 1:2 else 2:1))
  fit <- jomo1cat.MCMCchain(y, Y.numcat = c(2L, 2L),
                           nburn = 1000000L, output = 0)
  omega <- fit$collectomega
  stopifnot(all(is.finite(fit$collectbeta)), all(is.finite(omega)),
            all(omega[1, 1, ] == 1), all(omega[2, 2, ] == 1),
            all(omega[1, 2, ] == omega[2, 1, ]),
            all(abs(omega[1, 2, ]) < 1))
  # Relabeling the second response reflects the posterior around zero.
  r <- direction * omega[1, 2, -seq_len(10000L)]
  stopifnot(abs(mean(r) - 4/9) < 0.02,
            abs(mean(r > 0.9) - tail_probability) < 0.015)
}

testthat::test_that("gm works", {
  n <- 300  # Number of observations
  p <- 1000  # Number of predictors included in model
  q <- 60  # Number of true predictors
  amplitude = 20 # signal amplitude (for noise level = 1)

  # Generate the variables from a multivariate normal distribution
  times <- 1:p
  sigma <- 1
  rho = 0.5
  H <- abs(outer(times, times, "-"))
  V <- sigma * rho^H
  set.seed(0)
  x = MASS::mvrnorm(n, mu = rep(0,p), Sigma=V)

  # Generate the response.
  beta_q = rnorm(q, 0, amplitude)/sqrt(n)
  true_index= sample(p,q)
  beta = numeric(p)
  beta[true_index]=beta_q
  mu = x %*% beta
  y = mu + rnorm(n)

  # Fit the Gaussian mirror model
  library(parallel)
  fit = GM::gm(y, x, ncores=15, do_simultaneous = TRUE)
  plot(fit$gm_statistics, beta, main = "do_simultaneous = TRUE")
  abline(v=0)
  fit_fast = GM::gm(y, x, ncores=15, do_simultaneous = FALSE)
  plot(fit_fast$gm_statistics, beta, main = "do_simultaneous = FALSE")
  abline(v=0)
  plot(fit$gm_statistics, fit_fast$gm_statistics, col = ifelse(beta!=0, "red", "black"), main = "Both")
})






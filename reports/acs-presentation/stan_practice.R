library(rstan)
library(bayesplot)
library(truncnorm)

data = list(
  N_obs = 100,
  N_mis = 5,
  y_obs = rtruncnorm(100, mean = 2, sd = 5, a = 1.5) + rnorm(100, mean = 0, sd = 5)
)


fit <- rstan::stan(file = "8schools.stan", data = data, init = 1)

data

//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
// saved as 8schools.stan

data {
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  real y_obs[N_obs];
}
parameters {
  real mu;
  real<lower=0> sigma_err;
  real<lower=0> sigma_true;
  real y_mis[N_mis];
  real<lower=1.5> y_mis_true[N_mis];
  real<lower=1.5> y_true[N_obs];
}
model {
  sigma_err ~ exponential_lpdf(1);
  sigma_true ~ exponential_lpdf(1);
  mu ~ normal(0, 1);
  for (n in 1:N_obs) 
  y_obs[n] ~ normal(y_true, sigma_err);
  for (n in 1:N_obs)
    y_true[n] ~ normal(mu, sigma_true) T[1.5,];
  for (n in 1:N_mis)
    y_mis[n] ~ normal(y_mis_true, sigma_err);
  for (n in 1:N_mis)
    y_mis_true[n] ~ normal(mu, sigma_true) T[1.5,];
}

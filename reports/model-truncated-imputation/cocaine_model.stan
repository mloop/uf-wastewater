// generated with brms 2.9.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  real lloq;
  real lb[N];  // lower truncation bounds;
  int<lower=0> Nmi;  // number of missings
  int<lower=1> Jmi[Nmi];  // positions of missings
  // data for group-level effects of ID 1
  int<lower=1> N_1;
  int<lower=1> M_1;
  int<lower=1> J_1[N];
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;
  int<lower=1> M_2;
  int<lower=1> J_2[N];
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;
  int<lower=1> M_3;
  int<lower=1> J_3[N];
  vector[N] Z_3_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;
  int<lower=1> M_4;
  int<lower=1> J_4[N];
  vector[N] Z_4_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector<lower=0, upper = lloq>[Nmi] Ymi;  // estimated missings
  real temp_Intercept;  // temporary intercept
  real<lower=0> sigma;  // residual SD
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // unscaled group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // unscaled group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // unscaled group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // unscaled group-level effects
}
transformed parameters {
  // group-level effects
  vector[N_1] r_1_1 = (sd_1[1] * (z_1[1]));
  // group-level effects
  vector[N_2] r_2_1 = (sd_2[1] * (z_2[1]));
  // group-level effects
  vector[N_3] r_3_1 = (sd_3[1] * (z_3[1]));
  // group-level effects
  vector[N_4] r_4_1 = (sd_4[1] * (z_4[1]));
}
model {
  vector[N] Yl = Y;
  vector[N] mu = temp_Intercept + rep_vector(0, N);
  Yl[Jmi] = Ymi;
  for (n in 1:N) {
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_4_1[J_4[n]] * Z_4_1[n];
  }
  // priors including all constants
  target += student_t_lpdf(temp_Intercept | 3, 0, 10);
  target += student_t_lpdf(sigma | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_1[1] | 0, 1);
  target += student_t_lpdf(sd_2 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_2[1] | 0, 1);
  target += student_t_lpdf(sd_3 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_3[1] | 0, 1);
  target += student_t_lpdf(sd_4 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_4[1] | 0, 1);
  // likelihood including all constants
  if (!prior_only) {
    for (n in 1:N) {
      target += lognormal_lpdf(Yl[n] | mu[n], sigma) -
        lognormal_lccdf(lb[n] | mu[n], sigma);
    }
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = temp_Intercept;
}

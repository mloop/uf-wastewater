# Purpose: the goal of this program is to fit a model for the concentrations of metabolites in the following nested groups: location, time, run, and metabolite. This particular analysis will attempt to fit an interaction between the time and metabolite groups.


library(tidyverse)
library(rethinking)
library(ggridges)


water <- read_tsv("../../data/water_cleaned.txt") %>%
  mutate(run = factor(run),
         has_value = if_else(is.na(value) == TRUE, 0, 1),
         time_pretty = as.character(time_pretty)) %>%
  group_by(metabolite_name) %>%
  mutate(total_values = sum(has_value)) %>%
  filter(value < 100000 | is.na(value) == TRUE,
         total_values > 6) 


# Fit Bayesian model

## Priors
options(mc.cores = parallel::detectCores())
set.seed(123)

# Use rethinking package to fit models
f <- alist(
  # Model
  y_obs ~ dnorm(y_est, sigma_error),# T[y_obs,],
  sigma_error <- sigma_error_metabolite[metabolite_name],
  y_est ~ dnorm( mu , sigma ),# T[0,],
  mu <- b + bj[location] + bk[time_pretty] + bl[run] + bm[metabolite_name],
  bk[time_pretty] ~ dnorm(0, sigma_time_pretty),
  sigma_error_metabolite[metabolite_name] ~ dcauchy( 0 , 1),
  b ~ dnorm(0, 5),
  bj[location] ~ dnorm(0, sigma_location),
  bl[run] ~ dnorm(0, sigma_location),
  bm[metabolite_name] ~ dnorm(0, sigma_metabolite_name),
   # Hyper priors
  mu_metabolite[metabolite_name] ~ dnorm(0, 1),
  epsilon ~ dcauchy(0, 1),
  sigma ~ dcauchy(0, 1),
  sigma_time_pretty ~ dcauchy(0, 1),
  sigma_location ~ dcauchy(0, 1),
  sigma_run ~ dcauchy(0, 1),
  sigma_metabolite_name ~ dcauchy(0, 1),
  sigma_error_metabolite ~ dcauchy(0, 1))

# Simpler one
# I'm still not sure this model is identifiable
f <- alist(
  # Model
  y_obs ~ dnorm(y_est, sigma_error),# T[y_obs,],
  #sigma_error <- sigma_error_metabolite[metabolite_name],
  y_est ~ dnorm( mu , sigma ),# T[0,],
  mu <- b + bj[location] + bk[time_pretty] + bl[run] + bm[metabolite_name],
  bk[time_pretty] ~ dnorm(0, sigma_time_pretty),
  #sigma_error_metabolite[metabolite_name] ~ dcauchy( 0 , 1),
  b ~ dnorm(0, 5),
  bj[location] ~ dnorm(0, sigma_location),
  bl[run] ~ dnorm(0, sigma_location),
  bm[metabolite_name] ~ dnorm(0, sigma_metabolite_name),
  # Hyper priors
  mu_metabolite[metabolite_name] ~ dnorm(0, 1),
  epsilon ~ dcauchy(0, 1),
  sigma ~ dcauchy(0, 1),
  sigma_time_pretty ~ dcauchy(0, 1),
  sigma_location ~ dcauchy(0, 1),
  sigma_run ~ dcauchy(0, 1),
  sigma_metabolite_name ~ dcauchy(0, 1),
  #sigma_error_metabolite ~ dcauchy(0, 1)
  sigma_error ~ dcauchy(0, 1)
)
  
fit <- map2stan(
  f,
  data = list(y_obs = water$value,
              location = water$location,
              time_pretty = water$time_pretty,
              run = water$run,
              metabolite_name = water$metabolite_name),
  constraints = list(mu = "real<lower = 0>"),
  iter = 3000,
  warmup = 2500,
  chains = 2,
  cores = 2
)
  ## Had some issues with convergence on default priors. I will use slightly more informative priors and more iterations in order to try and speed convergence.

prior <- c(prior(normal(0, 5), class = "Intercept"),
           prior(exponential(5), class = "sd"))

fit <- brm(value | trunc(lb = 0) + mi() ~ (1 | location) + (1 | time_pretty) + (1 | run) + (1 + time_pretty | metabolite_name),
           family = lognormal(link = "identity", link_sigma = "log"),
           prior = prior,
           data = water,
           chains = 2,
           iter = 5000)

# Save the model
save(fit, file = "bayesian-model-concentration-combine-metabolites-inxn-time-metabolite.RData")


# Summarize the model
fit


# Check convergence
plot(fit)

# Posterior-predictive checks
## So many missing values that complete case posterior predictive checks may not be super helpful


pred_data <- modelr::data_grid(water, location, time_pretty, run, metabolite_name)

posterior_pred_time_metabolite <- posterior_predict(fit, newdata = pred_data, re_formula = ~ (1 | time_pretty) + (1 + time_pretty | metabolite_name)) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(pred_data, .) %>%
  group_by(time_pretty, metabolite_name) %>%
  slice(1) %>%
  select(-run, -location) %>%
  gather(draw, prediction, starts_with("V"))

# Time lapse and metabolite
p <- posterior_pred_time_metabolite %>%
  ggplot(aes(x = prediction, y = factor(time_pretty))) +
  geom_density_ridges(color = "gray") +
  facet_wrap(~ metabolite_name) +
  labs(
    title = "Predicted densities of metabolite concentrations",
    y = "Time of collection",
    x = "Concentration of metabolite (ng/mL)"
  ) +
  ggthemes::theme_tufte()

ggsave(filename = "posterior-predictions-time-metabolites-density-ixn-time-metabolite.png", p, dpi = 600, width = 12, height = 8, units = "in")

p <- posterior_pred_time_metabolite %>%
  summarise(
    mean = quantile(prediction, probs = 0.5),
    q_25 = quantile(prediction, probs = 0.25),
    q_75 = quantile(prediction, probs = 0.75)
  ) %>%
  ggplot(aes(x = time_pretty, y = mean)) +
  geom_pointrange(aes(ymin = q_25, ymax = q_75), size = 0.1) +
  facet_wrap(~ metabolite_name) +
  ggthemes::theme_tufte() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Median, 25th, and 75th percentiles of concentration of each metabolite",
    x = "Time of collection",
    y = "Concentration of metabolite (ng/mL)"
  )

ggsave(filename = "posterior-predictions-time-metabolites-pointrange-inxn-time-metabolite.png", p, dpi = 600)

# Posterior means
posterior_mean_time_metabolite <- posterior_linpred(fit, newdata = pred_data, re_formula = ~ (1 | time_pretty) + (1 + time_pretty | metabolite_name)) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(pred_data, .) %>%
  group_by(time_pretty, metabolite_name) %>%
  slice(1) %>%
  select(-run, -location) %>%
  gather(draw, prediction, starts_with("V"))

posterior_mean_time_metabolite %>%
  summarise(
    mean = quantile(prediction, probs = 0.5),
    q_25 = quantile(prediction, probs = 0.25),
    q_75 = quantile(prediction, probs = 0.75)
  ) %>%
  ggplot(aes(x = time_pretty, y = mean)) +
  geom_pointrange(aes(ymin = q_25, ymax = q_75), size = 0.1) +
  facet_wrap(~ metabolite_name) +
  ggthemes::theme_tufte() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Median, 25th, and 75th percentiles of concentration of each metabolite",
    x = "Time of collection",
    y = "Concentration of metabolite (ng/mL)"
  )

# Site by drug

posterior_pred_time_metabolite <- posterior_predict(fit, newdata = pred_data, re_formula = ~ (1 | location) + (1 + time_pretty | metabolite_name)) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(pred_data, .) %>%
  group_by(location, time_pretty, metabolite_name) %>%
  slice(1) %>%
  select(-run) %>%
  gather(draw, prediction, starts_with("V"))

posterior_pred_time_metabolite %>%
  summarise(
    mean = quantile(prediction, probs = 0.5),
    q_25 = quantile(prediction, probs = 0.25),
    q_75 = quantile(prediction, probs = 0.75)
  ) %>%
  ggplot(aes(x = time_pretty, y = mean, color = factor(location))) +
  geom_pointrange(aes(ymin = q_25, ymax = q_75), size = 0.1) +
  facet_wrap(~ metabolite_name, scales = "free_y") +
  ggthemes::theme_tufte() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Median, 25th, and 75th percentiles of concentration of each metabolite",
    x = "Time of collection",
    y = "Concentration of metabolite (ng/mL)"
  )

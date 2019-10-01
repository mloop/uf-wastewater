# Purpose: the goal of this program is to fit a model for the concentrations of metabolites in the following nested groups: location, time, run, and metabolite. This particular analysis will attempt to fit an interaction between the time and metabolite groups.


library(tidyverse)
library(brms)
library(ggridges)
library(truncnorm)

water <- read_tsv("../../data/water_cleaned.txt") %>%
  mutate(extraction = factor(extraction),
         machine = factor(machine),
         has_value = if_else(is.na(value) == TRUE, 0, 1),
         time_pretty = as.character(time_pretty)) %>%
  group_by(metabolite) %>%
  mutate(total_values = sum(has_value)) %>%
  filter(total_values > 6) 


# Fit Bayesian model

## Priors
options(mc.cores = parallel::detectCores())
set.seed(123)
## Had some issues with convergence on default priors. I will use slightly more informative priors and more iterations in order to try and speed convergence.

prior <- c(prior(normal(0, 5), class = "Intercept"),
           prior(cauchy(0, 1), class = "sd"),
           prior(cauchy(0, 1), class = "sigma"))

fit_cocaine <- brm(value | trunc(lb = 0, ub = 20) + mi() ~ (1 | location) + (1 | time_pretty) + (1 | extraction) + (1 | machine),
           family = lognormal(link = "identity", link_sigma = "log"),
           prior = prior,
           data = filter(water, metabolite == "Cocaine"),
           chains = 4,
           iter = 10000,
           sample_prior = TRUE, 
           control = list(adapt_delta = 0.9))

fit_benzo <- brm(value | trunc(lb = 0, ub = 20) + mi() ~ (1 | location) + (1 | time_pretty) + (1 | extraction) + (1 | machine),
                   family = lognormal(link = "identity", link_sigma = "log"),
                   prior = prior,
                   data = filter(water, metabolite == "Benzoylecgonine"),
                   chains = 4,
                   iter = 10000,
                   sample_prior = TRUE, 
                 control = list(adapt_delta = 0.9))

# Save the model
save(fit_cocaine, file = "bayesian-model-concentration-cocaine.RData")
save(fit_benzo, file = "bayesian-model-concentration-benzo.RData")


# Summarize the model
fit_cocaine
priors_cocaine <- prior_samples(fit_cocaine)

fit_benzo


# Check convergence
plot(fit_cocaine)

# Posterior-predictive checks
## So many missing values that complete case posterior predictive checks may not be super helpful
pp_check(fit_cocaine, newdata = filter(water, metabolite == "Cocaine") %>% drop_na()) +
  geom_vline(xintercept = c(0, 20))

pp_check(fit_benzo, newdata = filter(water, metabolite == "Benzoylecgonine") %>% drop_na()) +
  geom_vline(xintercept = c(0, 20))

# Check predicted values that were missing
as_tibble(fit_cocaine) %>%
  select(contains("Ymi")) %>%
  gather(obs, prediction) %>%
  ggplot(aes(x = prediction, y = obs)) +
  ggridges::geom_density_ridges()


pred_data <- modelr::data_grid(filter(water, metabolite == "Cocaine"), location, time_pretty, extraction, machine)

pred_data_benzo <- modelr::data_grid(filter(water, metabolite == "Benzoylecgonine"), location, time_pretty, extraction, machine)

posterior_mean_cocaine <- posterior_linpred(fit_cocaine, newdata = pred_data, re_formula = ~ (1 | time_pretty)) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(pred_data, .) %>%
  group_by(time_pretty, metabolite) %>%
  slice(1) %>%
  select(-extraction, -machine, -location) %>%
  gather(draw, prediction, starts_with("V")) %>%
  nest() %>%
  mutate(
    sigma_draws = map(data, ~as_tibble(fit_cocaine) %>%
                        select(sigma))
  ) %>%
  unnest() %>%
  mutate(
    post_mean = etruncnorm(a = 0, b = 20, mean = prediction, sd = sigma)
  )

posterior_mean_benzo <- posterior_linpred(fit_benzo, newdata = pred_data_benzo, re_formula = ~ (1 | time_pretty)) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(pred_data, .) %>%
  group_by(time_pretty, metabolite) %>%
  slice(1) %>%
  select(-extraction, -machine, -location) %>%
  gather(draw, prediction, starts_with("V")) %>%
  nest() %>%
  mutate(
    sigma_draws = map(data, ~as_tibble(fit_benzo) %>%
                        select(sigma))
  ) %>%
  unnest() %>%
  mutate(
    post_mean = etruncnorm(a = 0, b = 20, mean = prediction, sd = sigma)
  )

posterior_pred_time_cocaine <- posterior_predict(fit_cocaine, newdata = pred_data, re_formula = ~ (1 | time_pretty)) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(pred_data, .) %>%
  group_by(time_pretty, metabolite) %>%
  slice(1) %>%
  select(-extraction, -machine, -location) %>%
  gather(draw, prediction, starts_with("V"))

posterior_pred_time_benzo <- posterior_predict(fit_benzo, newdata = pred_data_benzo, re_formula = ~ (1 | time_pretty)) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(pred_data, .) %>%
  group_by(time_pretty, metabolite) %>%
  slice(1) %>%
  select(-extraction, -machine, -location) %>%
  gather(draw, prediction, starts_with("V"))

p <- posterior_pred_time_benzo %>%
  group_by(time_pretty) %>%
  summarise(
    mean = quantile(prediction, probs = 0.5),
    q_25 = quantile(prediction, probs = 0.025),
    q_75 = quantile(prediction, probs = 0.975)
  ) %>%
  mutate(metabolite = "Benzoylecgonine",
         lloq = 0.05) %>%
  bind_rows(., posterior_pred_time_cocaine %>%
              group_by(time_pretty) %>%
              summarise(
                mean = quantile(prediction, probs = 0.5),
                q_25 = quantile(prediction, probs = 0.025),
                q_75 = quantile(prediction, probs = 0.975)
              ) %>%
              mutate(metabolite = "Cocaine",
                     lloq = 0.15)) %>%
  ggplot(aes(x = time_pretty, y = mean)) +
  geom_pointrange(aes(ymin = q_25, ymax = q_75)) +
  geom_hline(aes(yintercept = lloq), linetype = "dashed") +
  facet_wrap(~ metabolite) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Median, 2.5th, and 97.5th percentiles of concentration of metabolite across all sites and runs",
    x = "Time of collection",
    y = "Concentration of metabolite (ng/mL)"
  )

ggsave(filename = "posterior-predictions-time-cocaine-benzo.png", p, dpi = 600)

p <- posterior_mean_benzo %>%
  ungroup() %>%
  group_by(time_pretty) %>%
  summarise(
    mean = mean(post_mean),
    q_025 = quantile(post_mean, probs = 0.025),
    q_975 = quantile(post_mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Benzoylecgonine",
         lloq = 0.05) %>%
  bind_rows(., posterior_mean_cocaine %>%
              group_by(time_pretty) %>%
              summarise(
                mean = mean(post_mean),
                q_025 = quantile(post_mean, probs = 0.025),
                q_975 = quantile(post_mean, probs = 0.975)
              ) %>%
              mutate(metabolite = "Cocaine",
                     lloq = 0.15)) %>%
  ggplot(aes(x = time_pretty, y = mean)) +
  geom_pointrange(aes(ymin = q_025, ymax = q_975)) +
  geom_hline(aes(yintercept = lloq), linetype = "dashed") +
  facet_wrap(~ metabolite) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Posterior mean (95% credible intervals) concentration of metabolite across all sites and runs",
    x = "Time of collection",
    y = "Concentration of metabolite (ng/mL)"
  )

ggsave(filename = "posterior-means-time-cocaine-benzo.png", p, dpi = 600)

# Standard deviation plots

sd_cocaine <- as_tibble(fit_cocaine) %>%
  select(contains("sd"), -contains("prior")) %>%
  mutate(iteration = seq(1, n())) %>%
  gather(source, posterior_sd, -iteration) %>%
  group_by(source) %>%
  summarise(median_sd = median(posterior_sd)) %>%
  mutate(metabolite = "Cocaine")

sd_benzo <- as_tibble(fit_benzo) %>%
  select(contains("sd"), -contains("prior")) %>%
  mutate(iteration = seq(1, n())) %>%
  gather(source, posterior_sd, -iteration) %>%
  group_by(source) %>%
  summarise(median_sd = median(posterior_sd)) %>%
  mutate(metabolite = "Benzoylecgonine")

p <- sd_cocaine %>%
  bind_rows(sd_benzo) %>%
  mutate(source = factor(source) %>% fct_recode("extraction" = "sd_extraction__Intercept",
                             "location" = "sd_location__Intercept",
                             "machine" = "sd_machine__Intercept",
                             "time" = "sd_time_pretty__Intercept")) %>%
  ggplot(aes(x = median_sd, y = fct_reorder(source, median_sd), color = metabolite)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    y = "Source of variation",
    x = "Median posterior SD",
    title = "Sources of variation in metabolite concentration"
  ) +
  scale_color_viridis_d()

ggsave(filename = "posterior-sd-by-source.png", p, dpi = 600, height = 5, units = "in")

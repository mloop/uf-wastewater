# Purpose: summarize the results of the fitted models

# EXCLUDE PSEUDOEPHEDRINE RIGHT NOW. THE MODEL RESULTS ARE NOT ANALYZABLE. RHAT > 1.1 IS COMMON.

library(tidybayes)
library(bayesplot)
library(tidyverse)
library(tidyselect)
library(rstan)

load("cocaine_fit.RData")
load("amphetamine_fit.RData")
load("benzoylecgonine_fit.RData")
load("codeine_fit.RData")
load("tramadol_fit.RData")
load("phentermine_fit.RData")
load("noroxycodone_fit.RData")
load("norhydrocodone_fit.RData")
load("oxycodone_fit.RData")

cocaine_intervals <- cocaine_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Cocaine")

amphetamine_intervals <- amphetamine_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Amphetamine")

benzoylecgonine_intervals <- benzoylecgonine_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Benzoylecgonine")

codeine_intervals <- codeine_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Codeine")

tramadol_intervals <- tramadol_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Tramadol")

phentermine_intervals <- phentermine_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Phentermine")

noroxycodone_intervals <- noroxycodone_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Noroxycodone")

norhydrocodone_intervals <- norhydrocodone_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Norhydrocodone")

oxycodone_intervals <- oxycodone_fit %>%
  rstan::extract(pars = c("b_Intercept", vars_select(names(.), contains("r_")), "sigma")) %>%
  as_tibble() %>%
  mutate(draw = seq(1, n())) %>%
  group_by(draw, b_Intercept, sigma) %>%
  select(contains("r_4")) %>%
  gather(time_point, effect, -draw, -b_Intercept, -sigma) %>%
  mutate(mu = b_Intercept + effect,
         mean = truncnorm::etruncnorm(a = 0, mean = mu, sd = sigma)) %>%
  ungroup() %>%
  group_by(time_point) %>%
  summarise(
    median_posterior = quantile(mean, probs = 0.5),
    low = quantile(mean, probs = 0.025),
    high = quantile(mean, probs = 0.975)
  ) %>%
  mutate(metabolite = "Oxycodone")

metabolite_intervals <- bind_rows(cocaine_intervals, amphetamine_intervals, codeine_intervals, benzoylecgonine_intervals, tramadol_intervals, phentermine_intervals, noroxycodone_intervals, norhydrocodone_intervals, oxycodone_intervals) %>%
  mutate(time_point = time_point %>% factor() %>%
           fct_recode("6:30" = "r_4_1[1]",
                      "7:00" = "r_4_1[2]",
                      "7:30" = "r_4_1[3]",
                      "8:00" = "r_4_1[4]",
                      "8:30" = "r_4_1[5]",
                      "9:00" = "r_4_1[6]",
                      "9:30" = "r_4_1[7]",
                      "10:00" = "r_4_1[8]",
                      "10:30" = "r_4_1[9]",
                      "11:00" = "r_4_1[10]",
                      "11:30" = "r_4_1[11]") %>%
           fct_relevel("6:30", "7:00", "7:30", "8:00", "8:30", "9:00", "9:30", "10:00", "10:30", "11:00", "11:30"))

metabolite_intervals %>%
  arrange(metabolite, time_point)

metabolite_intervals %>%
  ggplot(aes(x = time_point, y = median_posterior)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  facet_wrap(~metabolite, scales = "free_y") +
  labs(
    title = "Distribution of posterior mean concentration"
  )

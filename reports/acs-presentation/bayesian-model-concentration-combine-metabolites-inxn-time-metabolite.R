# Purpose: the goal of this program is to fit a model for the concentrations of metabolites in the following nested groups: location, time, run, and metabolite. This particular analysis will attempt to fit an interaction between the time and metabolite groups.


library(tidyverse)
library(brms)
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
  facet_wrap(~ metabolite_name, scales = "free_x") +
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

ggsave(filename = "posterior-predictions-time-metabolites-pointrange-inxn-time-metabolite.png", p, dpi = 600)

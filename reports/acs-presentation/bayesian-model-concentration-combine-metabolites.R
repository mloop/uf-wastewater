# Purpose: the goal of this program is to fit a model for the concentrations of metabolites in the following nested groups: location, time, run, and metabolite. Because we have such a small effective sample size, we will attempt to partially pool information across all metabolites. This procedure essentially assumes that if we know something about the concentration of one metabolite, we know about about the concentrations of other metabolites.


library(tidyverse)
library(brms)
library(ggridges)

water <- read_tsv("../../data/water_cleaned.txt") %>%
  mutate(run = factor(run),
         has_value = if_else(is.na(value) == TRUE, 0, 1)) %>%
  group_by(metabolite_name) %>%
  mutate(total_values = sum(has_value)) %>%
  filter(value < 100000 | is.na(value) == TRUE,
         total_values > 6) 


# Fit Bayesian model

## Priors
## Had some issues with convergence on default priors. I will use slightly more informative priors and more iterations in order to try and speed convergence.

prior <- c(prior(normal(0, 5), class = "Intercept"),
           prior(exponential(5), class = "sd"))

options(mc.cores = parallel::detectCores())
set.seed(123)

fit <- brm(value | trunc(lb = 0) + mi() ~ (1 | location) + (1 | time_lapse_hours) + (1 | run) + (1 | metabolite_name),
           family = lognormal(link = "identity", link_sigma = "log"),
           prior = prior,
           data = water,
           chains = 2,
           iter = 5000)

# Save the model
save(fit, file = "bayesian-model-concentration-combine-metabolites.RData")


# Summarize the model
fit


# Check convergence
plot(fit)

# Posterior-predictive checks
## So many missing values that complete case posterior predictive checks may not be super helpful


pred_data <- modelr::data_grid(water, location, time_lapse_hours, run, metabolite_name)

posterior_mean_time_metabolite <- posterior_linpred(fit, newdata = pred_data, re_formula = ~ (1 | time_lapse_hours) + (1 | metabolite_name)) %>%
  t() %>%
  as_tibble() %>%
  bind_cols(pred_data, .) %>%
  group_by(time_lapse_hours, metabolite_name) %>%
  slice(1) %>%
  select(-run, -location) %>%
  gather(draw, prediction, starts_with("V"))

# Look into effects of each grouping variable

b = data.frame(fit) %>% as_tibble()
b %>% select(contains("metabolite_name"), -contains("sd")) %>% summarise_all(mean) %>% gather(metabolite, mean_concentration) %>% arrange(mean_concentration)
b %>% select(contains("time_lapse"), -contains("sd")) %>% summarise_all(mean) %>% gather(time, mean_concentration) %>% arrange(mean_concentration) %>%
  ggplot(aes(x = time, y = mean_concentration)) +
  geom_point() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Time lapse and metabolite
p <- posterior_mean_time_metabolite %>%
  ggplot(aes(x = prediction, y = factor(time_lapse_hours))) +
  geom_density_ridges() +
  facet_wrap(~ metabolite_name)

ggsave(filename = "posterior-mean-time-metabolites.png", p, dpi = 600)

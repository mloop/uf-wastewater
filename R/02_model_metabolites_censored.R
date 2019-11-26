# Purpose: model each of 10 most common metabolites individually, in their own models

# Author: Matthew Shane Loop, PhD

# Packages
library(tidyverse)
library(brms)

# Options
options(mc.cores = parallel::detectCores())
set.seed(789301784)

# Data
water <- read_tsv("../data/water_cleaned.txt") %>% mutate_if(is.character, ~na_if(., "")) %>%
  mutate(time_pretty = as.character(time_pretty)) %>%
  group_by(metabolite) %>%
  mutate(non_missing = if_else(is.na(value) == FALSE, 1, 0),
         total_non_missing = sum(non_missing)) %>%
  filter(total_non_missing > 50) %>%
  mutate(censored_value = if_else(is.na(value) == TRUE, lloq,
                                  if_else(value < lloq, lloq,
                                          if_else(value > uloq, uloq, value)
                                  )
  ),
  log_value = log(censored_value),
  censored = if_else(censored_value == lloq, "left",
                     if_else(censored_value == uloq, "right", "none"))
  )

# Create model fitting functions
fit_metabolites <- function(metabolite_to_model, prior, adapt_delta = 0.99){
  
  fit <- brm(log_value | cens(censored) ~ time_pretty + (1 | machine) + (1 | extraction) + (1 | location),
             family = gaussian(),
             data = filter(water, metabolite == metabolite_to_model),
             iter = 2000,
             chains = 4,
             sample_prior = TRUE,
             prior = prior,
             control = list(adapt_delta = adapt_delta))
  saveRDS(fit, file = paste0("../output/02_model_metabolites_censored_", metabolite_to_model, ".rds", sep = ""))
}

# Fit the models

## Priors
prior_most_metabolites <- c(prior(normal(0, 1), class = "Intercept"),  # These priors were chosen from iteratively doing prior predictive checks
           prior(cauchy(0, 0.5), class = "sd"),
           prior(cauchy(0, 0.5), class = "sigma"),
           prior(normal(0, 1), class = "b"))

prior_noroxycodone <- c(prior(normal(0, 0.05), class = "Intercept"),  # These priors were chosen from iteratively doing prior predictive checks
                           prior(cauchy(0, 0.25), class = "sd"),
                           prior(cauchy(0, 0.25), class = "sigma"),
                           prior(normal(0, 0.05), class = "b"))

prior_pseudoephedrine <- c(prior(normal(0, 5), class = "Intercept"),  # These priors were chosen from iteratively doing prior predictive checks
                           prior(cauchy(0, 0.5), class = "sd"),
                           prior(cauchy(0, 0.5), class = "sigma"),
                           prior(normal(0, 1), class = "b"))


fit_metabolites("Amphetamine", prior_most_metabolites)
fit_metabolites("Benzoylecgonine", prior_most_metabolites)
fit_metabolites("Cocaine", prior_most_metabolites)
fit_metabolites("Hydrocodone", prior_most_metabolites)
fit_metabolites("Norhydrocodone", prior_most_metabolites)
fit_metabolites("Noroxycodone", prior_noroxycodone)
fit_metabolites("Oxycodone", prior_most_metabolites)
fit_metabolites("Phentermine", prior_most_metabolites, adapt_delta = 0.999)
fit_metabolites("Pseudoephedrine", prior_pseudoephedrine)
fit_metabolites("Tramadol", prior_most_metabolites)
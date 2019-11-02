# Purpose: create a Bayesian hierarchical model for predicting metabolite concentration by time, location, extraction, and machine. Average over censored values, making sure that you include the lower limit of detection for each metabolite appropriately. Use only the metabolites with 7 or more non-missing values

# We will actually just do the metabolites with the 10 highest number of observed values. Some of the low values are just too hard.

# Load packages
library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Read in data

water <- read_tsv("../../data/water_cleaned.txt") %>%
  mutate(time_pretty = as.character(time_pretty)) %>%
  mutate(censored = if_else(is.na(value) == TRUE, "left", "none"),
         log_value = log(value),
         log_lloq = log(lloq))

water_common <- water %>%
  group_by(metabolite) %>%
  mutate(measured = if_else(is.na(value) == FALSE, 1, 0),
         total_measurements = sum(measured)) %>%
  filter(total_measurements >= 7) %>%
  ungroup() %>%
  mutate(value = na_if(value, 0),
         censored = if_else(is.na(value) == TRUE, "left", "none"))
  

# Cocaine model

model <- make_stancode(
  value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
  data = water_common %>%
    filter(metabolite == "Cocaine"),
  family = lognormal(link = "identity", link_sigma = "log")
)

cocaine_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                           data = water_common %>%
                             filter(metabolite == "Cocaine"),
                           family  =lognormal())
cocaine_data$lloq <- 0.15
cocaine_fit <- stan(file = "cocaine_model.stan", 
                    data = cocaine_data,
                    chains = 4,
                    iter = 4000,
                    control = list(max_treedepth = 12))
save(cocaine_fit, file = "cocaine_fit.RData")

# Amphetamine
amphetamine_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                                  data = water_common %>%
                                    filter(metabolite == "Amphetamine"),
                                  family  =lognormal())
amphetamine_data$lloq <- 0.05
amphetamine_fit <- stan(file = "cocaine_model.stan", 
                        data = amphetamine_data,
                        chains = 4,
                        iter = 4000)
save(amphetamine_fit, file = "amphetamine_fit.RData")

# Benzoylecgonine
benzoylecgonine_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                                  data = water_common %>%
                                    filter(metabolite == "Benzoylecgonine"),
                                  family  =lognormal())
benzoylecgonine_data$lloq <- 0.05
benzoylecgonine_fit <- stan(file = "cocaine_model.stan", 
                        data = benzoylecgonine_data,
                        chains = 4,
                        iter = 4000)
save(benzoylecgonine_fit, file = "benzoylecgonine_fit.RData")

# Codeine
codeine_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                                      data = water_common %>%
                                        filter(metabolite == "Codeine"),
                                      family  =lognormal())
codeine_data$lloq <- 0.05
codeine_fit <- stan(file = "cocaine_model.stan", 
                            data = codeine_data,
                            chains = 4,
                            iter = 4000)
save(codeine_fit, file = "codeine_fit.RData")

# Tramadol
tramadol_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                              data = water_common %>%
                                filter(metabolite == "Tramadol"),
                              family  =lognormal())
tramadol_data$lloq <- 0.05
tramadol_fit <- stan(file = "cocaine_model.stan", 
                    data = tramadol_data,
                    chains = 4,
                    iter = 4000)
save(tramadol_fit, file = "tramadol_fit.RData")

# Hydrocodone
hydrocodone_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                               data = water_common %>%
                                 filter(metabolite == "Hydrocodone"),
                               family  =lognormal())
hydrocodone_data$lloq <- 0.05
hydrocodone_fit <- stan(file = "cocaine_model.stan", 
                     data = hydrocodone_data,
                     chains = 4,
                     iter = 4000)
save(hydrocodone_fit, file = "hydrocodone_fit.RData")

# Pseudoephedrine

## The pseudoephedrine data is *really* hard to fit right now, just to capture that one extreme value (I think).
pseudoephedrine_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                               data = water_common %>%
                                 filter(metabolite == "Pseudoephedrine"),
                               family  =lognormal())
pseudoephedrine_data$lloq <- 0.05
pseudoephedrine_fit <- stan(file = "cocaine_model.stan", 
                     data = pseudoephedrine_data,
                     chains = 4,
                     iter = 8000,
                     warmup = 5000)
save(pseudoephedrine_fit, file = "pseudoephedrine_fit.RData")

# Phentermine
phentermine_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                                      data = water_common %>%
                                        filter(metabolite == "Phentermine"),
                                      family  =lognormal())
phentermine_data$lloq <- 0.05
phentermine_fit <- stan(file = "cocaine_model.stan", 
                            data = phentermine_data,
                            chains = 4,
                            iter = 4000,
                            warmup = 2000)
save(phentermine_fit, file = "phentermine_fit.RData")

# Noroxycodone
noroxycodone_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                                  data = water_common %>%
                                    filter(metabolite == "Noroxycodone"),
                                  family  =lognormal())
noroxycodone_data$lloq <- 0.05
noroxycodone_fit <- stan(file = "cocaine_model.stan", 
                        data = noroxycodone_data,
                        chains = 4,
                        iter = 4000,
                        warmup = 2000)
save(noroxycodone_fit, file = "noroxycodone_fit.RData")

# Norhydrocodone
norhydrocodone_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                                   data = water_common %>%
                                     filter(metabolite == "Norhydrocodone"),
                                   family  =lognormal())
norhydrocodone_data$lloq <- 0.05
norhydrocodone_fit <- stan(file = "cocaine_model.stan", 
                         data = norhydrocodone_data,
                         chains = 4,
                         iter = 6000,
                         warmup = 3000,
                         control = list(adapt_delta = 0.9))
save(norhydrocodone_fit, file = "norhydrocodone_fit.RData")

# Oxycodone
oxycodone_data <- make_standata(value | mi() + trunc(0, ) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
                                     data = water_common %>%
                                       filter(metabolite == "Oxycodone"),
                                     family  =lognormal())
oxycodone_data$lloq <- 0.05
oxycodone_fit <- stan(file = "cocaine_model.stan", 
                           data = oxycodone_data,
                           chains = 4,
                           iter = 8000,
                           warmup = 4000,
                           control = list(adapt_delta = 0.9))
save(oxycodone_fit, file = "oxycodone_fit.RData")

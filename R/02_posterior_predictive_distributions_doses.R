# Create posterior predictive distributions of doses per 1,000

# Packages
library(tidyverse)
library(brms)
library(cowplot)

# Options
set.seed(98589534)

# Dataset
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

# Read in models
fit_amphetamine <- readRDS(file = "../output/02_model_metabolites_censored_Amphetamine.rds")
fit_benzo <- readRDS(file = "../output/02_model_metabolites_censored_Benzoylecgonine.rds")
fit_cocaine <- readRDS(file = "../output/02_model_metabolites_censored_Cocaine.rds")
fit_hydrocodone <- readRDS(file = "../output/02_model_metabolites_censored_Hydrocodone.rds")
fit_norhydrocodone <- readRDS(file = "../output/02_model_metabolites_censored_Norhydrocodone.rds")
fit_noroxycodone <- readRDS(file = "../output/02_model_metabolites_censored_Noroxycodone.rds")
fit_oxycodone <- readRDS(file = "../output/02_model_metabolites_censored_Oxycodone.rds")
fit_phentermine <- readRDS(file = "../output/02_model_metabolites_censored_Phentermine.rds")
fit_pseudoephedrine <- readRDS(file = "../output/02_model_metabolites_censored_Pseudoephedrine.rds")
fit_tramadol <- readRDS(file = "../output/02_model_metabolites_censored_Tramadol.rds")

results <- tibble(models = list(fit_amphetamine, fit_benzo, fit_cocaine, fit_hydrocodone, fit_norhydrocodone, fit_noroxycodone, fit_oxycodone, fit_phentermine, fit_pseudoephedrine, fit_tramadol)) %>%
  mutate(
    metabolite = c("Amphetamine", "Benzoylecgonine", "Cocaine", "Hydrocodone", "Norhydrocodone", "Noroxycodone", "Oxycodone", "Phentermine", "Pseudoephedrine", "Tramadol")
  )

# Read in ancillary data from stadium needed to estimate concentrations
metabolite_info <- read_tsv("../data/03_metabolism_data.txt")

flow <- read_tsv("../data/flow_rate.txt") %>%
  mutate(time_pretty = as.character(time_pretty),
         liters_previous_30_minutes = gallons_previous_30_minutes * 3.78541)

stadium_info <- read_tsv("../data/stadium_seating.txt")

# Posterior predictive distribution, accounting for all random effects
predicted_consumption_locations <- results %>%
  group_by(metabolite) %>%
  mutate(
    posterior_predictions = map2(models, metabolite, ~posterior_predict(.x, re_formula = ~NULL) %>%
                                   t() %>%
                                   as_tibble() %>%
                                   bind_cols(filter(water, metabolite == .y), .) %>%
                                   ungroup() %>%
                                   dplyr::select(metabolite, time_pretty, location, machine, extraction, V1:V4000) %>%
                                   gather(iteration, mean_log_conc, -time_pretty, -location, -machine, -extraction, -metabolite) %>%
                                   mutate(
                                     pred_concentration = exp(mean_log_conc) * 1000,  # original data reported as ng/mL; convert to liters for mass load calculations
                                     metabolite = .y
                                   ) %>%
                                   left_join(., metabolite_info, by = "metabolite") %>%
                                   left_join(., flow, by = "time_pretty") %>%  # flow rate isn't specific to location, which is a limitation
                                   left_join(., stadium_info, by = "location") %>%
                                   mutate(location_weight = 1 / proportion_of_stadium) %>% 
                                   mutate(
                                     mass_load = pred_concentration * (liters_previous_30_minutes / location_weight) * (100 / (100 + stability)) * 1e-6,  # water flow is for the whole stadium, so need to split the flow appropriately across locations. That's why (liters / location_weight) is used
                                     consumption_per_1000 = mass_load * 100 * (1 / excretion) * mwpar_mwmet * 1000 / (80651 * proportion_of_stadium) / typical_dose_mg, # unit is doses per 1000,
                                     consumption_missing = if_else(is.na(consumption_per_1000) == TRUE, 1, 0)
                                     )
                                   )
  ) %>%
  unnest(posterior_predictions) %>%
  ungroup()

predicted_consumption_stadium <- predicted_consumption_locations %>%
  group_by(metabolite, time_pretty, machine, extraction, iteration) %>%
  summarise(mass_load_stadium = sum(mass_load),
            consumption_per_1000_stadium = sum(consumption_per_1000)) %>%
  mutate(consumption_missing_stadium = if_else(is.na(consumption_per_1000_stadium) == TRUE, 1, 0)) %>%
  ungroup()

saveRDS(predicted_consumption_locations, "../output/02_posterior_predictive_doses_locations.rds")

saveRDS(predicted_consumption_stadium, "../output/02_posterior_predictive_doses_stadium.rds")

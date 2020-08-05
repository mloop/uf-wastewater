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

## This calculation is a bit complicated and there have been a lot of false starts. Therefore, here is a description of the calculation in words, before doing it in code. Steps 1 - 2 are done in code below, but the other steps are performed in individual figure or table files. But they are located here so that the entire thought process is documented in one place.

## 1. Predict the concentration of the substance at a given location, time, machine, and extraction

## 2. Calculate mass load at each of these groups, based upon water flow in that part of the stadium over the previous 30 minutes. Flow is provided for the whole stadium, so multiply that flow by the proportion of the stadium covered by each location. And concentrations are provided in ng/mL, so you have to multiply by 1,000 to get the per L to work with the flow data (in L). Also, to get to mg instead of ng, multiply by 1e-6.

## 3. To calculate mass load over all locations for each sample time, simply sum over the locations within a given time, machine, and extraction group.

## 4. To get mass load over the entire game and whole stadium, sum the mass load over all locations and times, within a given machine and extraction group.

## 5. To calculate consumption, for the whole stadium over the entire course of the game (only calculate consumption at this stage, per Bikram Subedi), use the following formula:

    ## doses = mass_load_in_mg * 100 * (1 / excretion) * (MW_parent / MW_metabolite) / typical_dose_in_mg


predicted_mass_load <- results %>%
  group_by(metabolite) %>%
  mutate(
    posterior_predictions = map2(models, metabolite, ~posterior_predict(.x, re_formula = ~NULL) %>%
                                   as_tibble(.name_repair = "unique") %>%  # Not including the .name_repair argument eventually creates an error. It's an odd error, because it only apprears if you have a clean R session.
                                   mutate(iteration = seq(1, n())) %>% 
                                   gather(observation, mean_log_conc, -iteration) %>% 
                                   spread(iteration, mean_log_conc) %>% 
                                   select(-observation) %>%
                                   bind_cols(filter(water, metabolite == .y), .) %>%
                                   ungroup() %>%
                                   dplyr::select(metabolite, time_pretty, location, machine, extraction, `1`:`10000`) %>%
                                   gather(iteration, mean_log_conc, -time_pretty, -location, -machine, -extraction, -metabolite) %>%
                                   mutate(
                                     pred_concentration = exp(mean_log_conc) * 1000,  # original data reported as ng/mL; convert to liters for mass load calculations
                                     metabolite = .y
                                   ))
  ) %>%
  ungroup() %>%
  dplyr::select(-metabolite) %>%
  unnest(posterior_predictions) %>%
  dplyr::select(-models) %>%
  left_join(., metabolite_info, by = "metabolite") %>%
 left_join(., flow, by = "time_pretty") %>%  # flow rate isn't specific to location, which is a limitation
 left_join(., stadium_info, by = "location") %>%
 mutate(
   mass_load = pred_concentration * liters_previous_30_minutes * proportion_of_stadium * (100 / (100 + stability)) * 1e-6,  # 
   mass_load_missing = if_else(is.na(mass_load) == TRUE, 1, 0)
   )

saveRDS(predicted_mass_load, "../output/02_posterior_predictive_mass_load.rds")

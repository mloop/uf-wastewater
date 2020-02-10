# Purpose: model the 10 most common metabolites multivariately

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



# The difficulty right now is figuring out how to structure the dataset so that we have both the value and the censoring for each variable

amphetamine <- water %>%
  filter(metabolite == "Amphetamine") %>%
  rename(amphetamine_log_value = log_value,
         amphetamine_censored = censored)%>%
  ungroup() %>%
  select(time_pretty, location, extraction, machine, contains("_log_value"), contains("_censored"))

benzo <- water %>%
  filter(metabolite == "Benzoylecgonine") %>%
  rename(benzo_log_value = log_value,
         benzo_censored = censored) %>%
  ungroup() %>%
  select(time_pretty, location, extraction, machine, contains("_log_value"), contains("_censored"))

water_mv_data <- inner_join(amphetamine, benzo, by = c("time_pretty", "location", "extraction", "machine"))

water_mv_model <- water %>%
  select(time_pretty, location, extraction, machine, metabolite, value) %>%
  spread(metabolite, value)

# You can't do censored multivariate models. One element of the vector may belong to the CDF (left censored), while another element may be part of the PDF (no censoring). You'd have to have the "counterfactual data point" of what that value would have been had it been censored or observed. Weird.

bf_amphetamine <- bf(amphetamine_log_value | cens(amphetamine_censored) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine))
bf_benzo <- bf(benzo_log_value | cens(benzo_censored) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine))

fit_mv <- brm(
  mvbind(Amphetamine + Benzoylecgonine + Cocaine) | mi() ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
  data = water_mv_model,
  chains = 4,
  cores = 4
)

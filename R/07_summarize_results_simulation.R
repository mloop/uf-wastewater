library(tidyverse)
library(arm)
library(truncnorm)

sim_fits <- read_rds("../output/06_simulated_fits.rds")

sim_fits %>%
  dplyr::select(prop_missing, iteration, int_low, int_high) %>%
  mutate(
    true_mean = etruncnorm(mean = 2.5, sd = sqrt(3.5 ^ 2 + 1 + 1), a = 0, b = Inf),
    covers = if_else(int_low < true_mean & int_high > true_mean, 1, 0)
  ) %>%
  ungroup() %>%
  group_by(prop_missing) %>%
  summarise(
    coverage = mean(covers)
  )

# So right now the coverage is 1, so Amanda appears to have been correct. I haven't implemented the missing data part yet, so that's a problem.
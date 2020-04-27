library(tidyverse)
library(arm)
library(truncnorm)

sim_data <- read_rds("../output/05_simulated_data.rds")

sim_fits <- sim_data %>%
  slice(1:10) %>%
  mutate(
    fits = map(data, ~lmer(y ~ (1 | location) + (1 | time), data = .,
                           control = lmerControl(optimizer = "bobyqa"))),
    confidence_int = map(fits, ~confint(.)),
    int_low = map_dbl(confidence_int, ~.[4, 1]),
    int_high = map_dbl(confidence_int, ~.[4, 2])
  )

write_rds(sim_fits, "../output/06_simulated_fits.rds")
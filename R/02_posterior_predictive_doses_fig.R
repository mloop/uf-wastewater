library(tidyverse)
library(brms)
library(cowplot)

# Read in models
readRDS("../output/02_posterior_predictive_doses_stadium.rds") %>%
  ungroup() %>%
  filter(consumption_missing_stadium == 0) %>%
  group_by(metabolite, iteration, extraction, machine) %>%
  summarise(consumption_per_1000_stadium_over_game = sum(consumption_per_1000_stadium)) %>%
  ungroup() %>%
  group_by(metabolite) %>%
  summarise(median_consumption = quantile(consumption_per_1000_stadium_over_game, probs = 0.5, na.rm = FALSE),
            low_consumption = quantile(consumption_per_1000_stadium_over_game, probs = 0.25, na.rm = FALSE),
            high_consumption = quantile(consumption_per_1000_stadium_over_game, probs = 0.75, na.rm = FALSE),
  ) %>%
  ggplot(aes(y = factor(metabolite) %>% fct_reorder(median_consumption), x = median_consumption * 80.651)) +
  geom_point() +
  ggpubr::theme_pubr() +
  labs(
    y = "Substance",
    x = "Doses",
    title = "Median predicted number of doses over entire game"
  ) -> p

ggsave(file = "../figs/02_posterior_predictive_doses.png", p, width = 7, height = 5, units = "in")

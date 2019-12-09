library(tidyverse)
library(brms)
library(cowplot)

# Read in models
readRDS("../output/02_posterior_predictive_doses.rds") %>%
  ungroup() %>%
  group_by(metabolite, iteration, extraction, machine, location) %>%
  mutate(consumption_per_1000_over_locations_over_game = sum(consumption_per_1000_over_locations)) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(metabolite) %>%
  filter(consumption_missing == 0) %>%
  summarise(median_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000_over_locations_over_game, weights = location_weight, probs = 0.5, na.rm = FALSE),
            low_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000_over_locations_over_game, weights = location_weight, probs = 0.25, na.rm = FALSE),
            high_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000_over_locations_over_game, weights = location_weight, probs = 0.75, na.rm = FALSE),
  ) %>%
  ggplot(aes(x = metabolite, y = median_mean_consumption * 80.651)) +
  geom_pointrange(aes(ymin = low_mean_consumption * 80.651, ymax = high_mean_consumption * 80.651)) +
  theme_bw() +
  ggrepel::geom_label_repel(aes(x = metabolite, y = median_mean_consumption * 80.651, label = prettyNum(round(median_mean_consumption * 80.651, digits = 0), big.mark = ","))) +
  labs(
    y = "Doses in stadium over entire game",
    x = "Compound",
    title = "Median and interquartile range of estimated doses for each compound over\nentire game"
  ) -> p

ggsave(file = "../figs/02_posterior_predictive_doses.png", p, width = 7, height = 4, units = "in")

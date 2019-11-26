library(tidyverse)
library(brms)
library(cowplot)

# Read in models
readRDS("../output/02_posterior_predictive_doses.rds") %>%
  ungroup() %>%
  group_by(metabolite) %>%
  filter(consumption_missing == 0) %>%
  summarise(median_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000, weights = location_weight, probs = 0.5, na.rm = FALSE),
            low_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000, weights = location_weight, probs = 0.25, na.rm = FALSE),
            high_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000, weights = location_weight, probs = 0.75, na.rm = FALSE),
  ) %>%
  ggplot(aes(x = metabolite, y = median_mean_consumption * 80.651)) +
  geom_pointrange(aes(ymin = low_mean_consumption * 80.651, ymax = high_mean_consumption * 80.651)) +
  ggrepel::geom_label_repel(aes(x = metabolite, y = median_mean_consumption * 80.651, label = prettyNum(round(median_mean_consumption * 80.651, digits = 0), big.mark = ","))) +
  labs(
    y = "Doses in stadium",
    x = "Compound",
    title = "Median and interquartile range of estimated distribution of doses for each compound"
  ) -> p

ggsave(file = "../figs/02_posterior_predictive_doses.png", p)

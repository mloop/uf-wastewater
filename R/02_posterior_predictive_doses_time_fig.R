library(tidyverse)
library(brms)
library(cowplot)

readRDS("../output/02_posterior_predictive_doses.rds") %>%
  ungroup() %>%
  group_by(metabolite, time_pretty) %>%
  filter(consumption_missing == 0) %>%
  summarise(median_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000, weights = location_weight, probs = 0.5, na.rm = FALSE),
            low_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000, weights = location_weight, probs = 0.25, na.rm = FALSE),
            high_mean_consumption = Hmisc::wtd.quantile(consumption_per_1000, weights = location_weight, probs = 0.75, na.rm = FALSE),
  ) %>%  
ggplot(aes(x = time_pretty, y = median_mean_consumption * 80.651)) +
  geom_pointrange(aes(ymin = low_mean_consumption * 80.651, ymax = high_mean_consumption * 80.651)) +
  facet_wrap(~ metabolite) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    y = "Doses in stadium",
    x = "Compound",
    title = "Median and interquartile range of estimated distribution of doses\nfor each compound, by time of collection"
  ) -> p

ggsave(file = "../figs/02_posterior_predictive_doses_time.png", p)

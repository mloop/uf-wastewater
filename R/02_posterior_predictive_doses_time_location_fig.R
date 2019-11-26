library(tidyverse)
library(brms)
library(cowplot)

readRDS("../output/02_posterior_predictive_doses.rds") %>%
  ungroup() %>%
  group_by(metabolite, time_pretty, location) %>%
  filter(consumption_missing == 0) %>%
  summarise(median_mean_consumption = quantile(consumption_per_1000, probs = 0.5, na.rm = FALSE),
            low_mean_consumption = quantile(consumption_per_1000, probs = 0.25, na.rm = FALSE),
            high_mean_consumption = quantile(consumption_per_1000, probs = 0.75, na.rm = FALSE),
  ) %>%  
  ggplot(aes(x = time_pretty, y = median_mean_consumption * 80.651, color = factor(location))) +
  geom_pointrange(aes(ymin = low_mean_consumption * 80.651, ymax = high_mean_consumption * 80.651), position = position_dodge(0.8)) +
  facet_wrap(~ metabolite) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_brewer(name = "Location") +
  theme(legend.position = "bottom") +
  labs(
    y = "Doses in stadium",
    x = "Compound",
    title = "Median and interquartile range of estimated distribution of doses\nfor each compound, by time of collection and location"
  ) -> p

ggsave(file = "../figs/02_posterior_predictive_doses_time_location.png", p)

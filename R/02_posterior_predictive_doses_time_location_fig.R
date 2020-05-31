library(tidyverse)
library(brms)
library(cowplot)

predicted_consumption <- readRDS("../output/02_posterior_predictive_doses_locations.rds") %>%
  ungroup() %>%
  filter(consumption_missing == 0) %>%
  group_by(metabolite, time_pretty, location) %>%
  summarise(median_mass_load = quantile(pred_concentration, probs = 0.5, na.rm = FALSE),
            low_mass_load = quantile(pred_concentration, probs = 0.25, na.rm = FALSE),
            high_mass_load = quantile(pred_concentration, probs = 0.75, na.rm = FALSE),
  )

predicted_consumption %>%
  ggplot(aes(x = time_pretty, y = median_mass_load, color = factor(location))) +
  geom_pointrange(aes(ymin = low_mass_load, ymax = high_mass_load), position = position_dodge(0.5), size = 0.3) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_color_grey(name = "Location") +
  labs(
    y = "Predicted concentration (ng/mL)",
    x = "Time of collection",
    title = "Median and interquartile range of estimated concentration (ng/mL) for each substance that passed\nthrough system over previous 30 minutes, by location"
  ) -> p

ggsave(file = "../figs/02_posterior_predictive_doses_time_location.png", p, width = 9, height = 5, units = "in")

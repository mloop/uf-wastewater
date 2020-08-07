library(tidyverse)
library(brms)
library(cowplot)

predicted_consumption <- readRDS("../output/02_posterior_predictive_mass_load.rds") %>%
  filter(mass_load_missing == 0) %>%
  ungroup() %>%
  group_by(metabolite, time_pretty, machine, extraction, iteration) %>%
  summarise(mass_load_stadium = sum(mass_load)) %>%
  ungroup() %>%
  group_by(metabolite, time_pretty) %>%
  summarise(median_mass_load = quantile(mass_load_stadium, probs = 0.5, na.rm = FALSE),
            low_mass_load = quantile(mass_load_stadium, probs = 0.25, na.rm = FALSE),
            high_mass_load = quantile(mass_load_stadium, probs = 0.75, na.rm = FALSE),
  )

predicted_consumption %>%
ggplot(aes(x = time_pretty, y = median_mass_load)) +
  geom_point(size = 0.3) +
  geom_path(aes(group = metabolite), size = 0.3) +
  facet_wrap(~ metabolite, ncol = 2) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8)
  ) +
  labs(
    y = "Estimated mass load (mg)",
    x = "Time of collection",
    title = "Median estimated mass load (mg) for each\nsubstance through system over previous\n30 minutes"
  ) -> p

ggsave(file = "../figs/02_posterior_predictive_mass_load_time.png", p, width = 5, height = 5, units = "in")
library(tidyverse)
library(brms)
library(cowplot)

predicted_consumption <- readRDS("../output/02_posterior_predictive_mass_load.rds") %>%
  ungroup() %>%
  filter(mass_load_missing == 0) %>%
  group_by(metabolite, time_pretty, location) %>%
  summarise(median_mass_load = quantile(mass_load, probs = 0.5, na.rm = FALSE),
            low_mass_load = quantile(mass_load, probs = 0.25, na.rm = FALSE),
            high_mass_load = quantile(mass_load, probs = 0.75, na.rm = FALSE),
  )

# The following line can adjust the order of analytes appearing on the figure.
predicted_consumption$metabolite <- factor(predicted_consumption$metabolite, levels=c("Oxycodone", "Norhydrocodone", "Hydrocodone", "Noroxycodone", 
                                                                        "Phentermine", " ", "Cocaine", "Benzoylecgonine", "Amphetamine", "  ",
                                                                        "Tramadol", "   ", "Pseudoephedrine", "    "))


predicted_consumption %>%
  ggplot(aes(x = time_pretty, y = median_mass_load, color = factor(location))) +
  geom_pointrange(aes(ymin = low_mass_load, ymax = high_mass_load), position = position_dodge(0.5), size = 0.15) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 2) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  ) +
  scale_color_grey(name = "Location") +
  labs(
    y = "Estimated mass load (mg)",
    x = "Time of collection" #,
#    title = "Median and interquartile range of estimated mass load (mg) for each substance that passed\nthrough system over previous 30 minutes, by location"
  ) -> p

ggsave(file = "../figs/02_posterior_predictive_mass_load_time_location.png", p, width = 9, height = 9, units = "in")


library(tidyverse)

water <- read_tsv("../../data/water_cleaned.txt")

p <- water %>%
  mutate(has_value = if_else(is.na(value) == TRUE, 0, 1)) %>%
  group_by(metabolite) %>%
  mutate(total_values = sum(has_value)) %>%
  filter(value < 100000 | is.na(value) == TRUE,
         total_values > 6) %>%
  ggplot(aes(x = time_pretty, y = value, color = metabolite, linetype = extraction, group = interaction(extraction, metabolite))) +
  geom_path() +
  facet_grid(machine ~location) +
  theme(
    legend.position = "bottom"
  ) +
  scale_color_viridis_d(name = "") +
  labs(
    title = "Observed concentrations of each metabolite over time, by run and location",
    subtitle = "Values are limited to those less than 100,000. Metabolites limited to those with at least 7 non-missing values.",
    y = "Concentration (ng/mL)"
  )

ggsave("observed-metabolite-time-series.png", p, dpi = 600,
       width = 10,
       height = 5,
       units = "in")

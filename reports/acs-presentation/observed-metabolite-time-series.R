library(tidyverse)

water <- read_tsv("../../data/water_cleaned.txt") %>%
  mutate(has_value = if_else(is.na(value) == TRUE, 0, 1),
         time_pretty = as.character(time_pretty)) %>%
  group_by(metabolite) %>%
  mutate(total_values = sum(has_value)) %>%
  filter(total_values > 6)

p <- water %>%
<<<<<<< Updated upstream
  mutate(has_value = if_else(is.na(value) == TRUE, 0, 1)) %>%
  group_by(metabolite) %>%
  mutate(total_values = sum(has_value)) %>%
  filter(value < 100000 | is.na(value) == TRUE,
         total_values > 6) %>%
  ggplot(aes(x = time_pretty, y = value, color = metabolite, linetype = extraction, group = interaction(extraction, metabolite))) +
  geom_path() +
  facet_grid(machine ~location) +
=======
  ungroup() %>%
  group_by(time_pretty, metabolite, location) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  ggplot(aes(x = time_pretty, y = mean_value, color = metabolite, group = metabolite)) +
  geom_path() +
  facet_wrap(~location) +
  theme_bw() +
>>>>>>> Stashed changes
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_viridis_d(name = "") +
  labs(
    title = "Observed concentrations of each metabolite over time, by location",
    subtitle = "Values are averaged over extractions and machine runs. Metabolites limited to those with at least 7 non-missing values.",
    y = "Concentration (ng/mL)",
    x = "Time"
  )

ggsave("observed-metabolite-time-series.png", p, dpi = 600)

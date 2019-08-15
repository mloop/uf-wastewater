library(tidyverse)
library(naniar)

water <- read_tsv("../../data/water_cleaned.txt") %>%
  mutate(run = factor(run),
         has_value = if_else(is.na(value) == TRUE, 0, 1),
         time_pretty = as.character(time_pretty)) %>%
  group_by(metabolite_name) %>%
  mutate(total_values = sum(has_value)) %>%
  filter(value < 100000 | is.na(value) == TRUE,
         total_values > 6) 

water %>%
  ggplot(aes(x = time_pretty, y = value)) +
  geom_miss_point() +
  facet_grid(metabolite_name ~ location)


water %>%
  spread(metabolite_name, value) %>%
  gg_miss_upset()

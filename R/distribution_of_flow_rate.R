library(skimr)
library(tidyverse)
library(haven)


water <- read_tsv("../data/water_cleaned.txt")

water 
library(Hmisc)

water %>%
  filter(date == "2018-09-08") %>%
  group_by(location, timestamp) %>%
  slice(1) %>%
  ggplot(aes(x = timestamp, y = flowrate)) +
  geom_point() +
  geom_smooth() +
  labs(
    title = "Flow rate across sites at each time point",
    subtitle = "Perhaps a time trend with breaks in game?"
  ) +
  facet_wrap(~ location) -> p

p
ggsave("flow_rate_distribution.png", p)

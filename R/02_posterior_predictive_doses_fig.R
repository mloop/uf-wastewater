library(tidyverse)
library(brms)
library(cowplot)

# Read in models
predicted_mass_loads <- readRDS("../output/02_posterior_predictive_mass_load.rds")

predicted_mass_loads %>%
  ungroup() %>%
  mutate(
    doses = mass_load * 100 / excretion * mwpar_mwmet / typical_dose_mg
  ) %>%
  filter(mass_load_missing == 0, is.na(doses) == FALSE) %>%
  group_by(metabolite, iteration, extraction, machine) %>%
  summarise(doses_whole_stadium_entire_game = sum(doses)) %>%
  ungroup() %>%
  group_by(metabolite) %>%
  summarise(median_doses = quantile(doses_whole_stadium_entire_game, probs = 0.5, na.rm = FALSE)
  ) %>%
  ggplot(aes(x = factor(metabolite) %>% fct_reorder(median_doses), y = median_doses)) +
  geom_bar(stat = "identity") +
  ggrepel::geom_text_repel(aes(label = round(median_doses, digits = 0) %>% prettyNum(., big.mark = ",")), nudge_y = 100) +
  ggpubr::theme_pubr() +
  labs(
    y = "Doses",
    x = "Substance",
    title = "Estimated median number of doses throughout stadium over entire game"
  ) +
  theme(
    text = element_text(size = 18)
  ) +
  coord_flip() -> p

ggsave(file = "../figs/02_posterior_predictive_doses.png", p, width = 15, height = 10, units = "in")

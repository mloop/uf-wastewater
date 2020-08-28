library(tidyverse)
library(brms)
library(cowplot)

# Read in models
predicted_mass_loads <- readRDS("../output/02_posterior_predictive_mass_load.rds")

predicted_mass_loads %>%
  ungroup() %>%
  mutate(
    doses = mass_load * 100 / excretion * mwpar_mwmet / typical_dose_parent_mg
  ) %>%
  filter(!(metabolite %in% c("Cocaine", "Oxycodone", "Hydrocodone"))) %>%
  group_by(metabolite, iteration, extraction, machine) %>%
  summarise(doses_whole_stadium_entire_game = sum(doses)) %>%
  ungroup() %>%
  group_by(metabolite) %>%
  summarise(median_doses = quantile(doses_whole_stadium_entire_game, probs = 0.5, na.rm = FALSE)
  ) %>%
  ggplot(aes(x = factor(metabolite) %>% fct_recode("Cocaine*" = "Benzoylecgonine", "Oxycodone*" = "Noroxycodone", "Hydrocodone*" = "Norhydrocodone") %>% fct_reorder(median_doses), y = median_doses)) +
  geom_bar(stat = "identity") +
  ggrepel::geom_text_repel(aes(label = round(median_doses, digits = 0) %>% prettyNum(., big.mark = ",")), nudge_y = 100) +
  ggpubr::theme_pubr() +
  labs(
    y = "Doses",
    x = "Substance",
    title = "Estimated median number of doses throughout stadium over entire game",
    caption = "*Number of doses are estimated from concentrations of the metabolites, excretion rate, ratio of molecular weights of\nparent compound and metabolite, and typical dose of the parent compound."
  ) +
  theme(
    text = element_text(size = 18),
    plot.caption = element_text(size = 10, hjust = 1)
  ) +
  coord_flip() -> p

ggsave(file = "../figs/02_posterior_predictive_doses.png", p, width = 15, height = 9, units = "in")

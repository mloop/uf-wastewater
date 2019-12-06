
library(Hmisc)
library(tidyverse)
library(GGally)
library(brms)
library(cowplot)



water <- read_tsv("../data/water_cleaned.txt") %>% mutate_if(is.character, funs(na_if(., ""))) %>%
  mutate(time_pretty = as.character(time_pretty),
         extraction = factor(extraction) %>% fct_recode("A" = "zach", "B" = "austin")) %>%
  group_by(metabolite) %>%
  mutate(non_missing = if_else(is.na(value) == FALSE, 1, 0),
         total_non_missing = sum(non_missing)) %>%
  mutate(censored_value = if_else(is.na(value) == TRUE, lloq, 
                                  if_else(value < lloq, lloq,
                                          if_else(value > uloq, uloq, value)
                                  )
  ),
  log_value = log(censored_value),
  censored = if_else(censored_value == lloq, "left",
                     if_else(censored_value == uloq, "right", "none"))
  )

water_analysis <- water %>%
  mutate(time_pretty = as.character(time_pretty)) %>%
  group_by(metabolite) %>%
  mutate(non_missing = if_else(is.na(value) == FALSE, 1, 0),
         total_non_missing = sum(non_missing)) %>%
  filter(total_non_missing > 50) %>%
  mutate(censored_value = if_else(is.na(value) == TRUE, lloq, 
                                  if_else(value < lloq, lloq,
                                          if_else(value > uloq, uloq, value)
                                  )
  ),
  log_value = log(censored_value),
  censored = if_else(censored_value == lloq, "left",
                     if_else(censored_value == uloq, "right", "none"))
  )

# Summary statistics for results

## What drugs were tested for?
water %>%
  distinct(metabolite)

## Drugs with no observed values?
water %>%
  mutate(prop_nonmissing = total_non_missing / n()) %>%
  filter(prop_nonmissing < 0.5, total_non_missing > 0) %>%
  distinct(metabolite) %>%
  nrow()

## Drugs observed in less than 50% of samples?
water %>%
  filter(total_non_missing == 0) %>%
  distinct(metabolite) %>%
  nrow()

## Drugs detected in every sample?
water %>%
  filter(total_non_missing == 98) %>%
  distinct(metabolite) %>%
  nrow()

# Figures

water %>%
  ggplot(aes(x = time_pretty, y = value, color = metabolite, linetype = factor(extraction) %>% fct_relevel("A"), group = interaction(extraction, metabolite))) +
  geom_path() +
  facet_grid(machine ~location) +
  scale_linetype(name = "Extraction") +
  scale_color_discrete(name = "Metabolite") +
  labs(
    x = "Time of collection",
    y = "Concentration (ng/mL)",
    title = "Observed concentration of metabolites with at least 1 nonmissing value, by location and\nmass spectrometer platform"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16)
  ) -> p
ggsave(filename = "../figs/01_longitudinal_all_metabolites.png", p)

water_analysis %>%
  ggplot(aes(x = censored_value)) +
  geom_histogram(bins = 100) +
  geom_vline(aes(xintercept = lloq), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = uloq), linetype = "dashed", color = "red") +
  facet_wrap(~ metabolite, scales = "free") +
  labs(x = "Measured concentration (ng/mL)") -> p
ggsave(filename = "../figs/01_histogram_top_metabolites.png", p)

water %>%
  select(metabolite, time_pretty, location, extraction, machine, value) %>%
  spread(metabolite, value) %>%
  select(-time_pretty, -location, -extraction, -machine) %>%
  visdat::vis_miss() +
  theme(axis.text.x.top = element_text(angle = 90)) -> p
ggsave(filename = "../figs/01_missing_data_pattern.png", p)

water %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100) +
  geom_vline(aes(xintercept = lloq), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = uloq), linetype = "dashed", color = "red") +
  facet_wrap(~ metabolite, scales = "free_x") +
  cowplot::theme_cowplot() +
  labs(x = "Observed concentration (ng/mL)") -> p

ggsave(filename = "prez-pics/histogram_all_metabolites.png", p, width = 15, height = 10, units = "in")



water %>%
  ggplot(aes(x = time_pretty, y = value, color = metabolite, linetype = extraction, group = interaction(extraction, metabolite))) +
  geom_path() +
  facet_grid(machine ~location) +
  cowplot::theme_cowplot() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) -> p
ggsave(filename = "prez-pics/longitudinal_all_metabolites.png", p, width = 15, height = 10, units = "in")

water %>%
  group_by(metabolite) %>%
  mutate(non_missing = if_else(is.na(value) == FALSE, 1, 0),
         total_non_missing = sum(non_missing),
         percent_non_missing = total_non_missing / n()) %>%
  select(metabolite, percent_non_missing) %>%
  filter(percent_non_missing != 0)



water %>%
  select(metabolite, time_pretty, location, extraction, machine, value) %>%
  spread(metabolite, value) %>%
  select(-time_pretty, -location, -extraction, -machine) %>%
  visdat::vis_miss() -> p
ggsave(filename = "prez-pics/missing_data.png", p, width = 10, height = 10, units = "in")



water_nomiss <- water %>%
  group_by(metabolite) %>%
  mutate(non_missing = if_else(is.na(value) == FALSE, 1, 0),
         total_non_missing = sum(non_missing)) %>%
  filter(total_non_missing > 50)

water_nomiss %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100) +
  geom_vline(aes(xintercept = lloq), linetype = "dashed", color = "red") +
  geom_vline(aes(xintercept = uloq), linetype = "dashed", color = "red") +
  facet_wrap(~ metabolite, scales = "free_x") +
  cowplot::theme_cowplot()+
  labs(x = "Observed concentration (ng/mL)") -> p
ggsave(filename = "prez-pics/histogram_common_metabolites.png", p, width=12, heigh=10, units = "in")



water_nomiss %>%
  select(time_pretty, location, extraction, machine, metabolite, value) %>%
  spread(metabolite, value) %>%
  mutate(Noroxycodone = if_else(Noroxycodone == 0, 0.01, Noroxycodone)) %>%
  select(5:14) %>%
  mutate_all(~log(.)) %>%
  ggscatmat(corMethod = "spearman") -> p

ggsave(filename = "prez-pics/matrix_plot.png", p, width  =12, height = 12, units = "in")

library(Hmisc)
library(tidyverse)
library(GGally)
library(brms)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

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

water_metabolites_grouped <- read_tsv("../data/water_cleaned_grouped.txt") %>% mutate_if(is.character, funs(na_if(., ""))) %>%
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

water_metabolites_group_sorted <- read_tsv("../data/water_cleaned_group_sorted.txt") %>% mutate_if(is.character, funs(na_if(., ""))) %>%
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

water_analysis <- water_metabolites_group_sorted %>%
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

## Drugs with at least 1 observed value
water %>%
  filter(total_non_missing > 0) %>%
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

# The following line can adjust the order of analytes appearing on the figure.
water_analysis$metabolite <- factor(water_analysis$metabolite, levels=c("Oxycodone", "Norhydrocodone", "Hydrocodone", "Noroxycodone", 
                                                                        "Phentermine", " ", "Cocaine", "Benzoylecgonine", "Amphetamine", "  ",
                                                                        "Tramadol", "   ", "Pseudoephedrine", "    "))

water_analysis %>%
  ggplot(aes(x = censored_value)) +
  geom_histogram(bins = 100) +
  geom_vline(aes(xintercept = lloq, color = "lower"), linetype = "dashed") +
  geom_vline(aes(xintercept = uloq, color = "higher"), linetype = "dashed") +
  scale_colour_manual(name = "limits", values = c(lower="green", higher="red")) +
  theme_bw() +
  facet_wrap(~ metabolite, scales = "free", ncol = 2, drop=FALSE) +
  labs(x = "Measured concentration (ng/mL)") -> p
ggsave(filename = "../figs/01_histogram_top_metabolites.png", p, width = 10, height = 10, units = "in")

water %>%
  select(metabolite, time_pretty, location, extraction, machine, value) %>%
  spread(metabolite, value) %>%
  select(-time_pretty, -location, -extraction, -machine) %>%
  visdat::vis_miss() +
  theme(axis.text.x.top = element_text(angle = 90)) -> p
ggsave(filename = "../figs/01_missing_data_pattern.png", p)

water_metabolites_grouped %>%
  select(metabolite, time_pretty, location, extraction, machine, value) %>%
  spread(metabolite, value) %>%
  select(-time_pretty, -location, -extraction, -machine) %>%
  visdat::vis_miss() +
  theme(axis.text.x.top = element_text(angle = 90, size = 10)) -> p
ggsave(filename = "../figs/01_missing_data_pattern_grouped.png", p)

###################################
# Advanced missing pattern figure #
###################################
col_order <- water_metabolites_group_sorted$metabolite[0:56]
# We have to choose Benzoylecgonine to match and resort rows, since others have lots of duplicated "NA"
row_order_Benzoylecgonine <- water_metabolites_group_sorted[water_metabolites_group_sorted$metabolite=="Benzoylecgonine",]$value
  
water_metabolites_group_sorted_spread <- water_metabolites_group_sorted %>%
  select(metabolite, time_pretty, location, extraction, machine, value) %>%
  spread(metabolite, value) %>%
  select(-time_pretty, -location, -extraction, -machine) 

water_metabolites_group_sorted_spread <- water_metabolites_group_sorted_spread[,col_order]
water_metabolites_group_sorted_spread <- water_metabolites_group_sorted_spread[match(row_order_Benzoylecgonine, water_metabolites_group_sorted_spread$Benzoylecgonine),]

water_metabolites_group_sorted_spread_rounded <- water_metabolites_group_sorted_spread
water_metabolites_group_sorted_spread_rounded[!is.na(water_metabolites_group_sorted_spread_rounded)] <- 1
water_metabolites_group_sorted_spread_rounded[is.na(water_metabolites_group_sorted_spread_rounded)] <- 0
heatmap_matrix <- data.matrix(water_metabolites_group_sorted_spread_rounded)
split_group_col <- water_metabolites_group_sorted$metabolite_group[0:56]
split_group_row <- c(rep(c("Shimadzu"), times = 33), rep(c("Waters_SPE1"), times = 33), rep(c("Waters_SPE2"), times = 33))
p <- Heatmap(heatmap_matrix, cluster_columns = FALSE, cluster_rows = FALSE, name = "missing pattern", 
             col = colorRamp2(c(0, 1), c("gray", "black")), column_split = split_group_col, row_split = split_group_row,
             column_gap = unit(1.5, "mm"))
png(filename="../figs/01_missing_data_pattern_splitted.png", width = unit(1000, "mm"), height = unit(500, "mm"))
p
dev.off()

# In the following dataframe, 0 is missing and 1 is detected
water_metabolites_group_sorted_spread_rounded$machine <- c(rep(c("Shimadzu"), times = 33), rep(c("Waters_SPE1"), times = 33), rep(c("Waters_SPE2"), times = 33))
write.csv(water_metabolites_group_sorted_spread_rounded, "../data/metabolites_missing_pattern.csv", row.names = FALSE)
######################## end #######################

#####################
# Advanced bar plot #
#####################



######################## end #######################

water %>%
  select(metabolite, time_pretty, location, extraction, machine, value) %>%
  spread(metabolite, value) %>%
  select(-time_pretty, -location, -extraction, -machine) %>%
  mutate(machine=rep(c("shimadzu", "waters_1", "waters_2"), times=33)) %>%
  subset(machine=="shimadzu") %>%
  select(-machine) %>%
  visdat::vis_miss() +
  theme(axis.text.x.top = element_text(angle = 90)) -> p
ggsave(filename = "../figs/01_missing_data_pattern_shimadzu.png", p)

water %>%
  select(metabolite, time_pretty, location, extraction, machine, value) %>%
  spread(metabolite, value) %>%
  select(-time_pretty, -location, -extraction, -machine) %>%
  mutate(machine=rep(c("shimadzu", "waters_1", "waters_2"), times=33)) %>%
  subset(machine=="waters_1") %>%
  select(-machine) %>%
  visdat::vis_miss() +
  theme(axis.text.x.top = element_text(angle = 90)) -> p
ggsave(filename = "../figs/01_missing_data_pattern_waters_spe1.png", p)

water %>%
  select(metabolite, time_pretty, location, extraction, machine, value) %>%
  spread(metabolite, value) %>%
  select(-time_pretty, -location, -extraction, -machine) %>%
  mutate(machine=rep(c("shimadzu", "waters_1", "waters_2"), times=33)) %>%
  subset(machine=="waters_2") %>%
  select(-machine) %>%
  visdat::vis_miss() +
  theme(axis.text.x.top = element_text(angle = 90)) -> p
ggsave(filename = "../figs/01_missing_data_pattern_waters_spe2.png", p)

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
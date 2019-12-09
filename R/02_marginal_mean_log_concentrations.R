
library(Hmisc)
library(tidyverse)
library(GGally)
library(brms)
library(cowplot)

set.seed(73984)

# Read in models
fit_amphetamine <- readRDS(file = "../output/02_model_metabolites_censored_Amphetamine.rds")
fit_benzo <- readRDS(file = "../output/02_model_metabolites_censored_Benzoylecgonine.rds")
fit_cocaine <- readRDS(file = "../output/02_model_metabolites_censored_Cocaine.rds")
fit_hydrocodone <- readRDS(file = "../output/02_model_metabolites_censored_Hydrocodone.rds")
fit_norhydrocodone <- readRDS(file = "../output/02_model_metabolites_censored_Norhydrocodone.rds")
fit_noroxycodone <- readRDS(file = "../output/02_model_metabolites_censored_Noroxycodone.rds")
fit_oxycodone <- readRDS(file = "../output/02_model_metabolites_censored_Oxycodone.rds")
fit_phentermine <- readRDS(file = "../output/02_model_metabolites_censored_Phentermine.rds")
fit_pseudoephedrine <- readRDS(file = "../output/02_model_metabolites_censored_Pseudoephedrine.rds")
fit_tramadol <- readRDS(file = "../output/02_model_metabolites_censored_Tramadol.rds")


results <- tibble(models = list(fit_amphetamine, fit_benzo, fit_cocaine, fit_hydrocodone, fit_norhydrocodone, fit_noroxycodone, fit_oxycodone, fit_phentermine, fit_pseudoephedrine, fit_tramadol)) %>%
  mutate(
    metabolite = c("Amphetamine", "Benzoylecgonine", "Cocaine", "Hydrocodone", "Norhydrocodone", "Noroxycodone", "Oxycodone", "Phentermine", "Pseudoephedrine", "Tramadol")
  )

marginal_plots <- results %>%
  group_by(metabolite) %>%
  mutate(
    means = map(models, ~marginal_effects(.x, method = "fitted", re_formula = ~NULL)),
    plots = map2(means, metabolite, ~plot(.x, plot = FALSE)[[1]] +
                   labs(
                     title = .y,
                     x = "Time of day",
                     y = "Mean natural log concentration [log(ng/mL)]"
                   ) +
                   geom_hline(color = "red", yintercept = log(0.05)) +
                   geom_hline(color = "red", yintercept = log(20)) +
                   theme_bw() +
                   theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
  )

plot_grid(plotlist = marginal_plots$plots) -> p

ggsave(filename = "../figs/02_marginal_means.png", p, width = 16, height = 12, units = "in")
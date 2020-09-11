library(Hmisc)
library(tidyverse)
library(GGally)
library(brms)
library(cowplot)

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

# Prior predictive checks

## It would take some serious thought to figure out prior predictive checks within the brms framework. Posterior predictive checks look fine.

# Posterior predictive checks

set.seed(8276343)

ppc_p <- results %>%
  group_by(metabolite) %>%
  mutate(
    pp = map2(models, metabolite, ~pp_check(.x, type = "hist", nsamples = 5) + labs(title = .y) + ggthemes::theme_tufte() + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 20), legend.position = "none", strip.background = element_blank()) + labs(title = .y))
  )

# Define order of subplots
p_list <- list(ppc_p$pp[[7]], ppc_p$pp[[6]], ppc_p$pp[[4]], ppc_p$pp[[5]],
            ppc_p$pp[[8]], NULL, ppc_p$pp[[3]], ppc_p$pp[[2]],
            ppc_p$pp[[1]], NULL, ppc_p$pp[[10]], NULL, ppc_p$pp[[9]], NULL)
# Plot
h_lines <- (1:6)/7
plots <- plot_grid(plotlist = p_list, ncol = 2, label_size = 10)
cowplot::plot_grid(plots, ncol = 1, rel_heights = c(1, 9)) + geom_hline(yintercept = h_lines) -> p

p
ggsave(filename = "../figs/02_posterior_predictive_checks.png", p, width = 15, height = 18)

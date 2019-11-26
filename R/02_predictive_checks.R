library(Hmisc)
library(tidyverse)
library(GGally)
library(brms)
library(cowplot)

# Read in models
load(file = "../output/02_model_metabolites_censored_Amphetamine.RData")
load(file = "../output/02_model_metabolites_censored_Benzoylecgonine.RData")
load(file = "../output/02_model_metabolites_censored_Cocaine.RData")
load(file = "../output/02_model_metabolites_censored_Hydrocodone.RData")
load(file = "../output/02_model_metabolites_censored_Norhydrocodone.RData")
load(file = "../output/02_model_metabolites_censored_Noroxycodone.RData")
load(file = "../output/02_model_metabolites_censored_Oxycodone.RData")
load(file = "../output/02_model_metabolites_censored_Phentermine.RData")
load(file = "../output/02_model_metabolites_censored_Pseudoephedrine.RData")
load(file = "../output/02_model_metabolites_censored_Tramadol.RData")

results <- tibble(models = list(fit_amphetamine, fit_benzo, fit_cocaine, fit_hydrocodone, fit_norhydrocodone, fit_noroxycodone, fit_oxycodone, fit_phentermine, fit_pseudoephedrine, fit_tramadol)) %>%
  mutate(
    metabolite = c("Amphetamine", "Benzoylecgonine", "Cocaine", "Hydrocodone", "Norhydrocodone", "Noroxycodone", "Oxycodone", "Phentermine", "Pseudoephedrine", "Tramadol")
  )

# Prior predictive checks

## It would take some serious thought to figure out prior predictive checks within the brms framework. Posterior predictive checks look fine.

# Posterior predictive checks

ppc_p <- results %>%
  group_by(metabolite) %>%
  mutate(
    pp = map2(models, metabolite, ~pp_check(.x, type = "hist", nsamples = 5) + labs(title = .y) + ggthemes::theme_tufte() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", strip.background = element_blank()) + labs(title = .y))
  )

title <- ggdraw() + 
  draw_label(
    "Posterior predictive distributions",
    fontface = 'bold',
    x = 0.25,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

plots <- plot_grid(plotlist = ppc_p$pp, ncol = 5)
cowplot::plot_grid(title, plots, ncol = 1, rel_heights = c(1, 9)) -> p

p
ggsave(filename = "../figs/02_posterior_predictive_checks.png", p, width = 10)
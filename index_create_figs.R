
library(Hmisc)
library(tidyverse)
library(GGally)
library(brms)
library(cowplot)



water <- read_tsv("data/water_cleaned.txt") %>% mutate_if(is.character, funs(na_if(., ""))) %>%
  mutate(time_pretty = as.character(time_pretty),
         extraction = factor(extraction) %>% fct_recode("A" = "zach", "B" = "austin"))

water_analysis <- read_tsv("data/water_cleaned.txt") %>% mutate_if(is.character, ~na_if(., "")) %>%
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



options(mc.cores = parallel::detectCores())
prior_amphetamine <- c(prior(normal(0, 10), class = "Intercept"),  # These priors were chosen from iteratively doing prior predictive checks
                       prior(cauchy(0, 5), class = "sd"),
                       prior(cauchy(0, 5), class = "sigma"),
                       prior(normal(0, 5), class = "b"))

fit_amphetamine_priors <- brm(log_value | cens(censored) ~ time_pretty + (1 | machine) + (1 | extraction) + (1 | location),
                              family = gaussian(),
                              data = filter(water_analysis, metabolite == "Amphetamine"),
                              iter = 2000,
                              chains = 1, 
                              sample_prior = "only",
                              prior = prior_amphetamine)

pp_check(fit_amphetamine_priors, type = "hist", binwidth = 0.1) -> p

ggsave(filename = "prez-pics/prior_predictive_checks_bad.png", p)



prior_amphetamine <- c(prior(normal(0, 1), class = "Intercept"),  # These priors were chosen from iteratively doing prior predictive checks
                       prior(cauchy(0, 0.5), class = "sd"),
                       prior(cauchy(0, 0.5), class = "sigma"),
                       prior(normal(0, 1), class = "b"))
fit_amphetamine_priors <- brm(log_value | cens(censored) ~ time_pretty + (1 | machine) + (1 | extraction) + (1 | location),
                              family = gaussian(),
                              data = filter(water_analysis, metabolite == "Amphetamine"),
                              iter = 2000,
                              chains = 1, 
                              sample_prior = "only",
                              prior = prior_amphetamine)



pp_check(fit_amphetamine_priors, type = "hist", binwidth = 0.1) -> p
ggsave(filename = "prez-pics/prior_predictive_checks_good.png", p)


load(file = "reports/02_model_metabolites_censored_amphetamine.RData")
load(file = "reports/02_model_metabolites_censored_benzoylecgonine.RData")
load(file = "reports/02_model_metabolites_censored_cocaine.RData")
load(file = "reports/02_model_metabolites_censored_hydrocodone.RData")
load(file = "reports/02_model_metabolites_censored_norhydrocodone.RData")
load(file = "reports/02_model_metabolites_censored_noroxycodone.RData")
load(file = "reports/02_model_metabolites_censored_oxycodone.RData")
load(file = "reports/02_model_metabolites_censored_phentermine.RData")
load(file = "reports/02_model_metabolites_censored_pseudoephedrine.RData")
load(file = "reports/02_model_metabolites_censored_tramadol.RData")

results <- tibble(models = list(fit_amphetamine, fit_benzo, fit_cocaine, fit_hydrocodone, fit_norhydrocodone, fit_noroxycodone, fit_oxycodone, fit_phentermine, fit_pseudoephedrine, fit_tramadol)) %>%
  mutate(
    metabolite = c("Amphetamine", "Benzoylecgonine", "Cocaine", "Hydrocodone", "Norhydrocodone", "Noroxycodone", "Oxycodone", "Phentermine", "Pseudoephedrine", "Tramadol")
  )

pp_plots <- results %>%
  group_by(metabolite) %>%
  mutate(
    pp = map2(models, metabolite, ~pp_check(.x, type = "hist", nsamples = 3) + labs(title = .y) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  )

plot_grid(plotlist = pp_plots$pp) -> p
ggsave(filename = "prez-pics/posterior_predictive_checks.png", p, width = 10, height = 8, units = "in")

marginal_plots <- results %>%
  group_by(metabolite) %>%
  mutate(
    means = map(models, ~marginal_effects(.x, method = "fitted", re_formula = ~NULL)),
    plots = map2(means, metabolite, ~plot(.x, plot = FALSE)[[1]] +
                   labs(
                     title = .y,
                     x = "Time of day",
                     y = "Mean log concentration [log(ng/mL)]"
                   ) +
                   geom_hline(color = "red", yintercept = log(0.05)) +
                   geom_hline(color = "red", yintercept = log(20)) +
                   cowplot::theme_cowplot() +
                   theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
  )

plot_grid(plotlist = marginal_plots$plots) -> p
ggsave(filename = "prez-pics/marginal_means.png", p, width = 16, height = 12, units = "in")

sigma_plots <- results %>%
  group_by(metabolite) %>%
  mutate(
    sd_plots = map2(models, metabolite, ~posterior_samples(.x) %>%
                      select(contains("sd"), sigma) %>%
                      mutate(iteration = seq(1, n())) %>%
                      gather(standard_deviation, draw, -iteration) %>%
                      ggplot(aes(x = draw, fill = standard_deviation)) +
                      geom_density() +
                      scale_fill_viridis_d() +
                      labs(title = .y)
    )
  )

plot_grid(plotlist = sigma_plots$sd_plots) -> p
ggsave(filename = "prez-pics/sd_plots.png", p, width = 16, height = 8, units = "in")

metabolite_info <- read_tsv("reports/03_metabolism_data.txt")
flow <- read_tsv("data/flow_rate.txt") %>%
  mutate(time_pretty = as.character(time_pretty),
         flow_rate_gallons_per_day = flow_rate * 1e6,
         flow_rate_liters_per_day = flow_rate_gallons_per_day * 3.78541)

predicted_consumption_overall <- results %>%
  group_by(metabolite) %>%
  mutate(
    posterior_predictions = map2(models, metabolite, ~posterior_predict(.x, re_formula = NA) %>%
                                   t() %>%
                                   as_tibble() %>%
                                   bind_cols(filter(water, metabolite == .y), .) %>%
                                   group_by(time_pretty) %>%
                                   slice(1) %>%
                                   ungroup() %>%
                                   dplyr::select(time_pretty, V1:V12000) %>%
                                   gather(iteration, mean_log_conc, -time_pretty) %>%
                                   mutate(
                                     mean_concentration = exp(mean_log_conc) * 1000,  # original data reported as ng/mL; convert to liters for mass load calculations
                                     metabolite = .y
                                   ) %>%
                                   left_join(., metabolite_info, by = "metabolite") %>%
                                   left_join(., flow, by = "time_pretty") %>%
                                   mutate(
                                     mass_load = mean_concentration * flow_rate_liters_per_day * (100 / (100 + stability)) * 1e-6,
                                     consumption_per_1000 = mass_load * 100 * (1 / excretion) * mwpar_mwmet * 1000 / 80651 / typical_dose_mg, # unit is doses per 1000
                                     consumption_missing = if_else(is.na(consumption_per_1000) == TRUE, 1, 0))
    ),
    quantiles_predictions = map(posterior_predictions, ~group_by(., time_pretty) %>%
                                  filter(consumption_missing == 0) %>%
                                  summarise(
                                    median_mean_consumption = median(consumption_per_1000),
                                    low_mean_consumption = quantile(consumption_per_1000, probs = 0.25),
                                    high_mean_consumption = quantile(consumption_per_1000, probs = 0.75),
                                  )
    ),
    plots = map2(quantiles_predictions, metabolite, ~ggplot(data = .x, aes(x = time_pretty, y = median_mean_consumption)) +
                   geom_pointrange(aes(ymin = low_mean_consumption, ymax = high_mean_consumption), position = position_dodge(0.8)) +
                   theme(
                     axis.text.x = element_text(angle = 45, hjust = 1)
                   ) +
                   labs(
                     y = "Doses per 1,000 attendees",
                     x = "Time of wastewater collection",
                     title = .y
                   ))
  )

predicted_consumption_overall %>%
  unnest(posterior_predictions) %>%
  ungroup() %>%
  group_by(metabolite) %>%
  filter(consumption_missing == 0) %>%
  summarise(median_mean_consumption = median(consumption_per_1000) * 80.651,
            low_mean_consumption = quantile(consumption_per_1000, probs = 0.25) * 80.651,
            high_mean_consumption = quantile(consumption_per_1000, probs = 0.75) * 80.651,
  ) %>%
  ggplot(aes(x = factor(metabolite) %>% fct_reorder(median_mean_consumption) %>% fct_rev(), y = median_mean_consumption)) +
  geom_pointrange(aes(ymin = low_mean_consumption, ymax = high_mean_consumption)) +
  theme_cowplot()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  labs(
    y = "Doses for 80,651 attendees",
    x = "Metabolite",
    title = "Estimated 25th, 50th, and 75th percentiles of number of\ndoses per 1,000 attendees"
  ) -> p

ggsave(filename = "prez-pics/dosage_grand.png", p, dpi = 300, width = 7, height = 7, units = "in")


predicted_consumption_overall %>%
  unnest(posterior_predictions) %>%
  ungroup() %>%
  group_by(metabolite) %>%
  filter(consumption_missing == 0, consumption_per_1000 < 100) %>%
  ggplot(aes(x = consumption_per_1000)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ metabolite) +
  theme_cowplot()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  labs(
    y = "Doses per 1,000 attendees",
    x = "Time of wastewater collection",
    title = "Estimated 25th, 50th, and 75th percentiles of number of\ndoses per 1,000 attendees"
  ) -> p

ggsave(filename = "prez-pics/dosage_grand_hist.png", p, width = 9, height = 5, units = "in")


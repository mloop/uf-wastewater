water_mv_model <- water %>%
  select(time_pretty, location, extraction, machine, metabolite, log_value) %>%
  spread(metabolite, log_value)

fit_mv <- brm(
  mvbind(Amphetamine, Benzoylecgonine, Cocaine, Hydrocodone, Norhydrocodone, Noroxycodone, Oxycodone, Phentermine, Pseudoephedrine, Tramadol) ~ (1 | time_pretty) + (1 | location) + (1 | extraction) + (1 | machine),
  data = water_mv_model,
  chains = 2,
  cores = 4
)

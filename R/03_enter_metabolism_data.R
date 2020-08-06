library(tidyverse)

metabolite_performance <- tribble(
  ~metabolite, ~mwpar_mwmet, ~stability, ~excretion, ~typical_dose_parent_mg,
  "Benzoylecgonine", 1.05, 5.5, 34.6, 50,
  "Norhydrocodone", 1.05, 0, 5, 10,  # Assume perfect stability; we don't actually know,
  "Noroxycodone", 1.05, 0, 23, 10  # Assume perfect stability; we don't actually know
)

write_tsv(metabolite_performance, "../data/03_metabolism_data.txt")
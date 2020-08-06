library(tidyverse)

metabolite_performance <- tribble(
  ~metabolite, ~mwpar_mwmet, ~stability, ~excretion, ~typical_dose_parent_mg,
  "Amphetamine", 1, 46.8, 33.2, 30,
  "Benzoylecgonine", 1.05, 5.5, 34.6, 50,
  "Cocaine", 1, -7.7, 5, 50,
  "Hydrocodone", 1, 9.6, 19.6, 10,	
  "Norhydrocodone", 1.05, 0, 5, 10,  # Assume perfect stability; we don't actually know,	
  "Noroxycodone", 1.05, 0, 23, 10,  # Assume perfect stability; we don't actually know	
  "Oxycodone", 1, 9.6, 11.5, 10,  # from https://www.drugs.com/dosage/oxycodone.html	
  "Phentermine", 1, 0, 48, 22.5, # Assume perfect stability; we don't actually know	
  "Pseudoephedrine", 1, -40.4, 88, 240,	
  "Tramadol", 1, -11, 29, 100
)

write_tsv(metabolite_performance, "../data/03_metabolism_data.txt")
library(skimr)
library(tidyverse)
library(haven)
library(readxl)
library(lubridate)

set.seed(123)

water <- read_xlsx("../data/20190508_UFWW_Results_Comparison.xlsx",
                   sheet = "analysis_data") %>% 
  mutate_if(is.character, ~na_if(., ""))

skim(water)

water_cleaned <- water %>%
  rename_all(funs(tolower(.))) %>%
  select(-time) %>%
  mutate(date = as_date(date),
         outofrange = if_else(outofrange == "+", "yes", outofrange),
         outofrange = if_else(is.na(concentration) == FALSE & is.na(outofrange), "no", outofrange)) %>%
  mutate(observation = seq(1, nrow(water))) %>%
  group_by(observation) %>%
  mutate(
    concentration_simulated = if_else(is.na(concentration), runif(1, 0, cutoff), concentration)
  ) %>%
  ungroup()

#water_cleaned %>%
 # write_tsv(path = "../data/water_cleaned.txt")
water %>%
   write_tsv(path = "../data/water_cleaned.txt")
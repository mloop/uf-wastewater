library(skimr)
library(tidyverse)
library(haven)
library(readxl)
library(lubridate)


water <- read_xlsx("../data/20181031 UF Waste Water Results Cleaned.xlsx") %>% mutate_if(is.character, funs(na_if(., "")))

skim(water)

water_cleaned <- water %>%
  rename(date_time = `Date&Time`) %>%
  rename_all(funs(tolower(.))) %>%
  select(-time) %>%
  mutate(date = as_date(date),
         outofrange = if_else(outofrange == "+", "yes", outofrange),
         outofrange = if_else(is.na(concentration) == FALSE & is.na(outofrange), "no", outofrange))

water_cleaned %>%
  write_tsv(path = "../data/water_cleaned.txt")

# Purpose: To import the provided datasets and provide a clean version for future analysis
# Authors: Matthew Loop and Dominick Lemas

library(skimr)
library(tidyverse)
library(lubridate)
library(readxl)
library(Hmisc)

water <- read_excel("../data/20190508_UFWW_Results_Comparison.xlsx") %>%
  mutate(hr = hour(time) %>% as.character(),
         min = minute(time) %>% as.character(),
         min = if_else(min == "0", "00", min),
         value = as.numeric(value)) %>%
  unite(time_pretty, c("hr", "min"), sep = ":") %>%
  rename_all(~tolower(.)) %>%
  mutate(lloq_ng_ml = na_if(lloq_ng_ml, "N/A"),
         uloq_ng_ml = na_if(lloq_ng_ml, "N/A"),
         barcode = as.character(barcode),
         location = factor(location),
         metabolite_number = factor(metabolite_number),
         run = factor(run),
         lloq_ng_ml = as.numeric(lloq_ng_ml),
         uloq_ng_ml = as.numeric(uloq_ng_ml)
  ) %>%
  group_by(location, metabolite_name) %>%
  arrange(time) %>%
  mutate(time_lapse_hours = (time - first(time)) %>% as.numeric(.) / 3600) %>%
  ungroup()

water %>%
   write_tsv(path = "../data/water_cleaned.txt")

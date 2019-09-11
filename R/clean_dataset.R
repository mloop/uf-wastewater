# Purpose: To import the provided datasets and provide a clean version for future analysis
# Authors: Matthew Loop and Dominick Lemas

library(tidyverse)
library(lubridate)
library(readxl)
library(Hmisc)

water <- read_excel("../data/UFWW Results Comparison-8-14-2019(1).xlsx",
                    sheet = "Multiple analyses") 

intermediate <- water %>% 
  slice(3:59) %>%
  rename(
    metabolite = Legend,
    uloq = `...2`,
    lloq = `Sample ID`
  ) %>%
  na_if("-") %>%
  na_if("N/A") %>%
  na_if("ND") %>%
  select(-starts_with("UF"))

times <- intermediate %>% slice(1) %>%
  gather(run, time, -metabolite, -uloq, -lloq) %>%
  select(run, time)

stacked <- intermediate %>%
  slice(2:59) %>%
  group_by(metabolite) %>%
  gather(run, value, -lloq, -uloq, -metabolite) %>%
  mutate(
    value = if_else(str_detect(value, "OQ"), 
                    str_sub(value, start = 13, end = -2),
                    value),
    value = if_else(str_detect(value, "OQ"), 
                    str_sub(value, start = 14, end = -1),
                    value),
    value = na_if(value, ""),
    value = if_else(str_detect(value, "Ratio"), 
                    str_sub(value, start = 15, end = -2),
                    value),
    value = na_if(value, "/Below LLO"),
    value = as.numeric(value)
  ) %>%
  left_join(., times, by = "run") %>%
  separate(time, into = c("date", "time_pretty", "location"), sep = " ") %>%
  mutate(
    location = str_sub(location, -1, -1),
    below_lloq = if_else(value < lloq, 1, 0),
    above_uloq = if_else(value > uloq, 1, 0),
    run = if_else(str_detect(run, "Shimadzu") == TRUE, "shimadzu",
                  if_else(str_detect(run, "Austin") == TRUE, "austin", "prior")),
    time_pretty = factor(time_pretty) %>% fct_relevel("6:30", "7:00", "7:30", "8:00", "8:30", "9:00", "9:30", "10:00", "10:30", "11:00"),
    extraction = if_else(run == "prior", "zach", "austin"),
    machine = if_else(run == "shimadzu", "shimadzu", "waters")
  ) %>%
  select(-run)

stacked %>%
   write_tsv(path = "../data/water_cleaned.txt")

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

## Now do flow rate
water <- read_tsv("../data/water_cleaned.txt") %>% mutate_if(is.character, funs(na_if(., ""))) %>%
  mutate(time_pretty = as.character(time_pretty),
         extraction = factor(extraction) %>% fct_recode("A" = "zach", "B" = "austin"))

flow <- read_excel("../data/Copy of WRF_Inf-flows.xls", sheet = "20180908WRF_Inf", skip = 2) %>%
  rename(millions_gallons_per_day = `09/08/18 MGD`) %>%
  select(Date, millions_gallons_per_day) %>%
  mutate(
    Date = factor(Date)
  ) %>%
  separate(Date, into = c("date", "time_pretty"), sep = " ") %>%
  mutate(
    gallons_per_minute = millions_gallons_per_day / 24 / 60 * 1e6,
    time = lubridate::hms(time_pretty) %>% as.numeric(),
    time_interval_minutes = seq(0, 5 * (n() - 1), by = 5),
  ) %>%
  filter(time_interval_minutes <= 330) %>%
  mutate(ahead = lead(gallons_per_minute),
         area = if_else(gallons_per_minute < ahead, 5 * gallons_per_minute + 0.5 * 5 * (ahead - gallons_per_minute),
                        5 * ahead + 0.5 * 5 * (gallons_per_minute - ahead))) %>%
  mutate(
    group = c(rep(1:11, each = 6), 12)
  ) %>%
  group_by(group) %>%
  mutate(gallons_previous_30_minutes = sum(area)) %>%
  ungroup() %>%
  mutate(prev = lag(gallons_previous_30_minutes)) %>%
  select(time_pretty, prev) %>%
  mutate(time_pretty = str_replace(time_pretty, "18", "06"),
         time_pretty = str_replace(time_pretty, "19", "07"),
         time_pretty = str_replace(time_pretty, "^20", "08"),
         time_pretty = str_replace(time_pretty, "21", "09"),
         time_pretty = str_replace(time_pretty, "22", "10"),
         time_pretty = str_replace(time_pretty, "23", "11")) %>%
  right_join(., water %>% ungroup() %>% distinct(time_pretty), by = "time_pretty") %>%
  rename(gallons_previous_30_minutes = prev)

flow %>%
  write_tsv("../data/flow_rate.txt")


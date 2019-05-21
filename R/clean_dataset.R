##-------------- 
# **************************************************************************** #
# ***************                Project Overview              *************** #
# **************************************************************************** #

# Author:      Matt Loop & Dominick Lemas 
# Date:        May 21, 2019 
# IRB:
# Description: identify missing data for UF waste water
# Data: '20190508 UFWW Results Comparison.xlsx'

# **************************************************************************** #
# ***************                Directory Variables           *************** #
# **************************************************************************** #

# directory variables
work.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\");work.dir
data.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\");data.dir
out.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\results\\");out.dir

# set working directory
setwd(work.dir)
list.files()

# **************************************************************************** #
# ***************                Library                       *************** #
# **************************************************************************** #

library(skimr)
library(tidyverse)
library(haven)
library(readxl)
library(lubridate)

# **************************************************************************** #
# ***************                Read Data                     *************** #
# **************************************************************************** #

set.seed(123)

# read data
data.file.name="20190508_UFWW_Results_Comparison.xlsx";data.file.name
data.file.path=paste0(data.dir,"",data.file.name);data.file.path
water <- read_xlsx(data.file.path, sheet="analysis_data") %>% 
  mutate_if(is.character, funs(na_if(., "")))

# **************************************************************************** #
# ***************               format data                    *************** #
# **************************************************************************** #

names(water)
str(water)

# factors
water$location=as.factor(water$location)
water$metabolite_number=as.factor(water$metabolite_number)

# format time
water$Location=as.factor(waste.dat$Location)


# Start here
skim(water)

water_cleaned <- water %>%
  rename(date_time = `Date&Time`) %>%
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

water_cleaned %>%
  write_tsv(path = "../data/water_cleaned.txt")

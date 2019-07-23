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
library(openxlsx)
library(xlsx)
library(lubridate)

# **************************************************************************** #
# ***************                Read Data                     *************** #
# **************************************************************************** #

set.seed(123)

# read data
data.file.name="20190508_UFWW_Results_Comparison.xlsx";data.file.name
data.file.path=paste0(data.dir,"",data.file.name);data.file.path

water=read.xlsx2(data.file.path, sheetIndex=1, header=TRUE, colClasses="character")%>%
  mutate_if(is.character, funs(na_if(., "")))


########################
# Matthew code to read in data
########################
water <- read_excel("20190508_UFWW_Results_Comparison.xlsx") %>%
  mutate(hr = hour(time) %>% as.character(),
         min = minute(time) %>% as.character(),
         min = if_else(min == "0", "00", min),
         value = as.numeric(value)) %>%
  unite(time_pretty, c("hr", "min"), sep = ":") %>%
  rename_all(~tolower(.)) %>%
  mutate(lloq_ng_ml = na_if(lloq_ng_ml, "N/A"),
         uloq_ng_ml = na_if(lloq_ng_ml, "N/A"))
water

# **************************************************************************** #
# ***************               format data : ugly!!                   ******* #
# **************************************************************************** #

names(water)
head(water)
str(water)

# character
water$metabolite_number=as.character(water$metabolite_number)
water$metabolite_name=as.character(water$metabolite_name)
water$ULOQ_ng_mL=as.character(water$ULOQ_ng_mL)
water$LLOQ_ng_mL=as.character(water$LLOQ_ng_mL)
water$sample_id=as.character(water$sample_id)
water$sample_id2=as.character(water$sample_id2)
water$barcode=as.character(water$barcode)
water$value=as.character(water$value)
water$FlowRate=as.character(water$FlowRate)

# numeric
water$ULOQ_ng_mL=as.numeric(water$ULOQ_ng_mL)
water$LLOQ_ng_mL=as.numeric(water$LLOQ_ng_mL)
water$run=as.numeric(water$run)
water$value=as.numeric(water$value)
water$FlowRate=as.numeric(water$FlowRate)

# check data
head(water)
str(water)

# factors
water$time=as.factor(water$time)
unique(water$time)

water.df=water%>%
  rename_all(funs(tolower(.)))%>%
  mutate(time = factor(time, levels = c("6_30_00 AM", "7_00_00 AM",
                                            "7_30_00 AM", "8_00_00 AM",
                                            "8_30_00 AM", "9_00_00 AM",
                                            "9_30_00 AM", "10_00_00 AM",
                                            "10_30_00 AM", "11_00_00 AM",
                                            "11_30_00 AM")),
         time=recode(time,
                      "6_30_00 AM"="6:30",
                      "7_00_00 AM"="7:00",
                      "7_30_00 AM"="7:30",
                      "8_00_00 AM"="8:00",
                      "8_30_00 AM"="8:30",
                      "9_00_00 AM"="9:00",
                      "9_30_00 AM"="9:30",
                      "10_00_00 AM"="10:00",
                      "10_30_00 AM"="10:30",
                      "11_00_00 AM"="11:00",
                      "11_30_00 AM"="11:30"))
# check levels
levels(water.df$time)
names(water.df)

# export
# save(water.df, file=paste0(data.dir,"water.RData"))

# **************************************************************************** #
# ***************               cleaned                   ******* #
# **************************************************************************** #

# Dom: cant quite figure out what you got going below.

# Matthew: Below was just indicating whether a value was out of range, and then simulating a value it was below the cutoff. This was just for the grant analysis, and we won't exactly do this in the real analysis.

# Start here

water_cleaned <- water %>%
  rename_all(funs(tolower(.))) %>%
  mutate(
    value_simulated = if_else(is.na(value), runif(1, 0, as.numeric(lloq_ng_ml)), value)
    ) %>%
  ungroup()

#water_cleaned %>%
 # write_tsv(path = "../data/water_cleaned.txt")
water_cleaned %>%
   write_tsv(path = "../data/water_cleaned.txt")

##-------------- 
# **************************************************************************** #
# ***************                Project Overview              *************** #
# **************************************************************************** #

# Author:      Dominick Lemas 
# Date:        May 22, 2019 
# IRB:
# Description: exploratory data analysis 

# **************************************************************************** #
# ***************                Directory Variables           *************** #
# **************************************************************************** #


# Directory Locations
work.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\");work.dir
data.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\");data.dir
out.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\results\\");out.dir

# Set Working Directory
setwd(work.dir)
list.files()

# **************************************************************************** #
# ***************                Library                       *************** #
# **************************************************************************** #

# library(readxl)
# install.packages(data.table)
# library(data.table)
library(tidyr)
library(dplyr)
library(tidyverse)

# **************************************************************************** #
# ***************  UF_WasteWater_ResultsCleaned.csv                                              
# **************************************************************************** # 

#Read Data
data.file.name="UF_WasteWater_ResultsCleaned.csv";data.file.name
data.file.path=paste0(data.dir,"",data.file.name);data.file.path
waste.dat<- read.csv(data.file.path);

# look at data
head(waste.dat); str(waste.dat); names(waste.dat)

# format time
waste.dat$Location=as.factor(waste.dat$Location)

# format time
unique(waste.dat$Time)
levels(waste.dat$Time)

# 11 time points
waste.dat$Time = factor(waste.dat$Time,
           levels = c("6:30:00 PM", 
                      "7:00:00 PM", 
                      "7:30:00 PM",
                      "8:00:00 PM",
                      "8:30:00 PM",
                      "9:00:00 PM",
                      "9:30:00 PM",
                      "10:00:00 PM",
                      "10:30:00 PM",
                      "11:00:00 PM",
                      "11:30:00 PM"))

# take a look at the data
time.by.location.count=waste.dat %>%
  group_by(Time) %>%
  count(Location)%>%
  write_tsv(path = paste0(out.dir,"time_by_location.txt"))

# missing 7:30 at location 1 (24 observations): Added 3/6/2019 by DJL
# missing 8:00 at location 3 (24 observations): Added 3/6/2019 by DJL
# missing 10:30 at location 3 (24 observations): Still missing as of 3/6/2019
# 
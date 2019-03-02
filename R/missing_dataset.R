##-------------- 
# **************************************************************************** #
# ***************                Project Overview              *************** #
# **************************************************************************** #

# Author:      Dominick Lemas 
# Date:        Fedb 21, 2019 
# IRB:
# Description: identify missing data for UF waste water
# Data: C:\Users\djlemas\Dropbox (UFL)\02_Projects\BEACH_STUDY\RedCap

# **************************************************************************** #
# ***************                Directory Variables           *************** #
# **************************************************************************** #

# Computer
location="djlemas";location

# Directory Locations

work.dir=paste("C:\\Users\\",location,"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\",sep="");work.dir
data.dir=paste("C:\\Users\\",location,"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\",sep="");data.dir
out.dir=paste("C:\\Users\\",location,"\\Dropbox (UFL)\\02_Projects\\WASTE_WATER\\2018_KY_Sep18\\",sep="");out.dir

# Set Working Directory
setwd(work.dir)
list.files()

# **************************************************************************** #
# ***************                Library                       *************** #
# **************************************************************************** #

# library(readxl)
# install.packages(data.table)
library(data.table)
library(tidyr)
library(dplyr)
# library(reshape2)

# **************************************************************************** #
# ***************  UF_WasteWater_ResultsCleaned.csv                                              
# **************************************************************************** # 

#Read Data
data.file.name="UF_WasteWater_ResultsCleaned.csv";data.file.name
data.file.path=paste0(data.dir,"",data.file.name);data.file.path
waste.dat<- read.csv(data.file.path);

# look at data
head(waste.dat); str(waste.dat); names(waste.dat)

# take a look at he data
waste.dat %>%
  group_by(Location) %>%
  count(Time)

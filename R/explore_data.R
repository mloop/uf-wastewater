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
library(ggplot2)

# **************************************************************************** #
# ***************                  load.RData                                              
# **************************************************************************** # 

#Read Data
data.file.name="water.RData";data.file.name
data.file.path=paste0(data.dir,"",data.file.name);data.file.path
load(data.file.path);

# look at data
water=water.df
head(water); str(water); names(water)

# count by location
time.by.location.count=water %>%
  group_by(time) %>%
  count(location)

# run 1 vs run 2
df=water%>%
  filter(run==1)%>
  
  # Interleaved histograms
ggplot(water, aes(x=value , color=run)) +
  geom_histogram(fill="white", position="dodge")+
  theme(legend.position="top")
# Add mean lines
p<-ggplot(df, aes(x=weight, color=sex)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")+
  theme(legend.position="top")
p
  

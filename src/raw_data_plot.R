#Date: 06/12/2018
#author: Aminath Shausan
#This program reads in raw data from "daily viremia levels.csv" and cleans it to be used 
#in myModel_fit.R
#-------------------------------------------------------
#clear history 
rm(list=ls())

#save the packages in this library
.libPaths("C:/Users/shausan/r-libraries") #use this on qut desktop
#.libPaths("C:/Users/shau/r-libraries") #use this on laptop
#-----------------------------------

#install packages 
install.packages("tidyverse")
library(tidyverse)
library(dplyr)
library(tidyr)

#-------------------------------------------------------
#set working directory
#first set wd as the main "dengue" project and then use the following
setwd(paste0(getwd(), "/src")) #set working directory to src folder
load('../data/DENV1_data.RData') #loads the saved data 
#plot viral load for each subject (with subject name)
ggplot(data = DENV1_data)+ 
  geom_point(mapping = aes(x=FeverDay, y=log10(viremia) , shape = viremia_sign) , show.legend = FALSE)+
  facet_wrap(~ StudyNo, ncol=9) 

#plot viral load for each subject (without subject name)
jpeg('../results/DENV1/DENV1_before_fit.jpeg')
postscript('../results/DENV1/DENV1_before_fit.eps', horizontal = FALSE, onefile = FALSE, paper = "special", height = 10, width = 10)
qplot(
  FeverDay,
  log10(viremia),
  data = DENV1_data, shape = viremia_sign,  color = Serology,
  group = StudyNo,
  geom = c('point'),
  xlab = 'Fever Day', xlim = c(-7,4),
  ylab = 'Viremia (log10 - copies/ml)',
    show.legend = FALSE
) + 
  facet_wrap(~ StudyNo) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
)
 # axis(1, at = -7:4)
#axis(1, -7:4)
dev.off()
###############################
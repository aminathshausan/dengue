#Date: 06/12/2018
#author: Aminath Shausan
#This program reads in raw data from "daily viremia levels.csv" and cleans it to be used 
#in myModel_fit.R
#-------------------------------------------------------
#clear history 
rm(list=ls())

#save the packages in this library
.libPaths("C:/Users/shausan/r-libraries") #use this on qut desktop
.libPaths("C:/Users/shau/r-libraries") #use this on laptop
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
#read dataset in "data" folder
dataset <- read.csv("../data/daily viremia levels.csv")

#filter DENV1 data only 
DENV1_data <- filter(dataset, Serotype == "DENV1")

#plot viral load for each subject without facet_wrap text
qplot(
  DOI,
  log10(viremia),
  data = DENV1_data, color = viremia_sign,
  group = StudyNo,
  geom = c('point'),
  xlab = 'Day of illness',
  ylab = 'Viremia (log10 - copies/ml)'
) + 
  facet_wrap(~ StudyNo) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

#Extract each patient data from "DENV1_data" set
subjects <- c(603, 604, 605, 607, 610, 615, 616, 
              618, 622, 623, 624, 627, 631, 632, 
              633, 635, 637, 639, 646, 647, 650, 
              651, 653, 654, 656, 658, 659, 660, 
              662, 663, 667, 668, 669, 673, 676,
              686, 687, 690, 695, 702, 706, 708, 
              710, 711, 714, 716, 721, 724, 806, 
              812, 813, 817, 818, 823, 826, 829,
              833, 834, 836, 837, 840, 845, 846, 
              849, 850, 851, 852, 855, 857, 859,
              865, 872, 874, 875, 877, 879, 883,
              888) 

for (i in 1:length(subjects)) {
  object_name <- paste0("p", i, "_data")  #paste0("subject", 1) results "subject1"
  assign(object_name, filter(DENV1_data, StudyNo == subjects[i])) #assign(x, value) assigns value to x
}

#save extracted patient data
save.image(file = "patient_data.RData")
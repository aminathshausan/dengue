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

#plot viral load for each subject (with subject name)
ggplot(data = DENV1_data)+ 
  geom_point(mapping = aes(x=DOI, y=log10(viremia) , shape = viremia_sign) , show.legend = FALSE)+
  facet_wrap(~ StudyNo, ncol=9) 

#plot viral load for each subject (without subject name)
qplot(
  DOI,
  log10(viremia),
  data = DENV1_data, shape = viremia_sign, 
  group = StudyNo,
  geom = c('point'),
  xlab = 'Day of illness',
  ylab = 'Viremia (log10 - copies/ml)',
  show.legend = FALSE
) + 
  facet_wrap(~ StudyNo) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

#Extract each patient data from "DENV1_data" set (use this for individual fits)
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
save.image(file = '../data/patient_data.RData')
#########################################
#extract viremia and DOI measurements from each subject (three subjects with 4 observations are removed from analysis)
#cluster subjects into two clusters 
#cl1 = linear profile. (=48 subjects)
subjects_cl1 <- c(603, 604,  607, 610, 615, 616, 
                  623, 627,  632, 637, 651,  654,
                  656, 658,   662, 663, 667, 668, 
                  673, 676,  686,  690, 695,  706, 
                  708, 711,  724, 806,  812, 
                  813,  818, 823, 826, 829,  833,  
                  836,   845, 846, 852,   859,865, 
                  872, 874, 875, 877, 879, 883, 888) 
viremia_cl1 <- matrix(data =NA, ncol =length(subjects_cl1), nrow = 5);
times_cl1 <- matrix(data =NA, ncol =length(subjects_cl1), nrow = 5);
is_censored <- matrix(data =NA, ncol =length(subjects_cl1), nrow = 5)
#length(subjects_cl1)
for (i in 1:length(subjects_cl1) ) {
  viremia_cl1[,i] <- filter(DENV1_data, StudyNo == subjects_cl1[i]) %>% select(viremia)  %>% unlist
  times_cl1[,i] <- filter(DENV1_data, StudyNo == subjects_cl1[i]) %>% select(DOI)  %>% unlist
  is_censored[,i] <- filter(DENV1_data, StudyNo == subjects_cl1[i]) %>% select(viremia_sign)  %>% unlist
  is_censored[,i] <- ifelse(is_censored[,i] %in% c('1'), 0, 1)
}

#change values in is_censored matrix so that 1 =0 and 2 =1
#ifelse(condition, result if TRUE, result if FALSE)
#is_censored[,1][5] <-2
#is_censored[,2][1] <-2

#for (i in 1:2){
#is_censored[,i] <- ifelse(is_censored[,i] %in% c('1'), 0, 1)
#}

save.image(file = '../data/cluster1_data.RData')

#cl2 = nonlinear profile (=27 subjects)
subjects_cl2 <- c(605, 618, 622,  624,  631, 
                  633,  639, 646, 647, 650, 
                  653,  659, 660, 669, 687, 702,   
                  710,  716, 817, 834,  837, 
                  840,  849, 850, 851,  855, 857) 



viremia_cl2 <- matrix(data =NA, ncol =length(subjects_cl2), nrow = 5);
times_cl2 <- matrix(data =NA, ncol =length(subjects_cl2), nrow = 5);
is_censored <- matrix(data =NA, ncol =length(subjects_cl2), nrow = 5)

for (i in 1:length(subjects_cl2) ) {
  viremia_cl2[,i] <- filter(DENV1_data, StudyNo == subjects_cl2[i]) %>% select(viremia)  %>% unlist
  times_cl2[,i] <- filter(DENV1_data, StudyNo == subjects_cl2[i]) %>% select(DOI)  %>% unlist
  is_censored[,i] <- filter(DENV1_data, StudyNo == subjects_cl2[i]) %>% select(viremia_sign)  %>% unlist
  is_censored[,i] <- ifelse(is_censored[,i] %in% c('1'), 0, 1)
}
save.image(file = '../data/cluster2_data.RData')

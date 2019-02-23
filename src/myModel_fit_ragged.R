#' ---
#' title: "my proposed model"
#' author: "Aminath Shausan"
#' date: "February, 12, 2019"
#' This R code is used to estimate parameters upto beta using "myModel-ragged.stan"
#'  parameters estimated are delta, std, gamma, kappa, eta, omega, beta
#'  ---------------
#clear history 
rm(list=ls())

#save the packages in this library
.libPaths("C:/Users/shausan/r-libraries") #use this on desktop
.libPaths("C:/Users/shau/r-libraries") #use this on laptop
#-----------------------------------
# install required packages 
library("rstan")
#options(width = 90)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
#--------------------------------------------------
#set working directory to src folder
#first set wd as the main "dengue" project and then use the following
setwd(paste0(getwd(), "/src")) #set working directory to src folder

#--------load cluster 1 ragged data and prepare it to fit model------------------------------
load('../data/cluster1_ragged_data.RData') #loads the saved data

#extract viremia and measurement times 
#J <- nrow(num_obs)                #number of patients
J <-2
N <- 10
#N <- nrow(cl1_data_ragged)        #total number of observations 
y_obs <- cl1_data_ragged$viremia[1:10]       #observed viremia measurements as a vector
sample_time <- cl1_data_ragged$FeverDay[1:10]   #measurement times 
idx <- rep(1:2,each=5) #index of patients (all secondary cases have 5 observations)

#gp_sizes <-  num_obs$n
#t_pred <- t_pred[, 1:J]
#is_censored <- is_censored[,1:J]

#create prediction time(s) matrix by adding IP of 6 days to each observed time(s) vector
t_pred <- matrix(data =NA, ncol =48, nrow = 11);
for (i in 1:48) {
  t_pred[,i] <- (times_cl1[,i][1]-6):times_cl1[,i][5]
}

#---------Estimating delta, std, gamma, kappa, eta, omega, beta------------
stan_data <- list(J = J,                     #number of patients
                  N = N,                 #number of total observed measurements
                  ts = sample_time,           #time points where observed data are collected
                   y_hat = y_obs,         #observed data
                   idx = idx         #index of patients
                  )
               
# Test / debug the model:                                    
test <- stan("myModel_ragged.stan",
             data = stan_data,
             chains = 1, iter = 1, 
             algorithm = "Fixed_param") 
fit_model = stan(fit = test,
                 data = stan_data,
                 chains = 4,     #number of Markov Chains
                 warmup = 1000, #number of warmup iterations per chain
                 iter = 2500,   #total number of iterations per chain
                 refresh = 100,
                 control = list(adapt_delta = 0.8, max_treedepth = 10))

save.image(file = "p1_upto_beta.RData")
#-------------------------------------------------------


load("patient1_upto_beta.RData") #loads the saved data

#' ---
#' title: "my proposed model"
#' author: "Aminath Shausan"
#' date: "December 20, 2018"
#' This R code is used to estimate parameters upto beta using "myModel-hierarchical.stan"
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

#--------load cluster 1 data and prepare it to fit model------------------------------
load('../data/cluster1_data.RData') #loads the saved data

#create prediction time(s) matrix by adding IP of 6 days to each observed time(s) vector
t_pred <- matrix(data =NA, ncol =48, nrow = 11);
for (i in 1:48) {
  t_pred[,i] <- (times_cl1[,i][1]-6):times_cl1[,i][5]
}

J <- 2                            #number of patients
y_obs <- viremia_cl1[, 1:J]       #observed viremia measurements 
sample_time <- times_cl1[, 1:J]   #measurement times 
t_pred <- t_pred[, 1:J]
is_censored <- is_censored[,1:J]


#---------Estimating delta, std, gamma, kappa, eta, omega, beta------------
#(date: 21/12/18)
stan_data <- list(J = J,                     #number of patients
                  n_obs = 5,                 #number of observed measurements
                  n_pred = 11,                #number of predicted measurements (from generated quantities block)
                  y0 = c(1e8, 0, 0, 0, 1, 1),  #initial state of ode 
                  t0 = -4,                     #starting point of ode 
                  ts = sample_time,           #time points where observed data are collected
                  t_pred =t_pred,           # time point where predictions are made from ode (in the generated quantities block)
                  A= 1.4e6,  
                  sigma = 0.5,   
                  y_hat = y_obs,               #observed data
                  is_censored = is_censored,    #indicator matrix 
                  L = 357)
#   run_estimation=1)       #switch on likelihood 

# Test / debug the model:                                    
test <- stan("myModel-hierarchical.stan",
             data = stan_data,
             chains = 1, iter = 100, 
             control = list(adapt_delta = 0.8, max_treedepth = 10)) 
fit_model = stan(fit = test,
                 data = stan_data,
                 chains = 4,     #number of Markov Chains
                 warmup = 1000, #number of warmup iterations per chain
                 iter = 2500,   #total number of iterations per chain
                 refresh = 100,
                 control = list(adapt_delta = 0.8, max_treedepth = 10))

save.image(file = '../results/p1_to_p2_noPool.RData')
#-------------------------------------------------------


load("patient1_upto_beta.RData") #loads the saved data

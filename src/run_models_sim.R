# Author: Aminath Shausan 
# Date: May 2022
# This program runs either  'model a' or 'model b' simulation using 
#'model_a_sim.stan' program or 'model_b_sim.stan' program, respectively
#------------------------------------------------------------
#clear history
rm(list=ls())

#-------------------------------------
#required pakages
library(tidyverse)
library(tidyr)
library(dplyr)
library(deSolve)
library("ggplot2")
library('loo')
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#================================================

#1. load pre-processed data and extract what is required for the Stan program
#DENV1_all data is used for simulation in both models and b
load('../data/DENV1_all.RData') #loads the preprocessed DENV1 data from Nguyen et al (2013)

#Clapham_DENV1 primary and secondary data used in the simulation
#loads the preprocessed DENV1 secondary/primary  data from Clapham et al (2016)
#load('../data/Clapham_DENV1_secondary_data.RData')
#load('../data/Clapham_DENV1_primary_data.RData')

#the 'DENV1_all.RData' contains the following matrices: 
#'viremia' = DENV1 meaurements for all patients (5 by 67 matrix)
#'times' = measured time points (5 by 67 matrix)
#'is_censored' = indicator matrix for each observation (0=above LOD, 1=below LOD) (5 by 67 matrix)
#'serology' = indicator matrix for each patients' type (1 = primary, 2 = secondary) (1 by 67 matrix)
n_patients <- ncol(viremia)
#t_pred <- matrix(data =NA, ncol =n_patients, nrow = 11); #predicted times
t_pred <- matrix(data =NA, ncol =n_patients, nrow = 19); #uncomment for Clapham data
t0 <- matrix(data = NA,  ncol = n_patients, nrow =1);   #initial times of ODE
for (i in 1:n_patients) {
 #t_pred[,i] <- (times[,i][1]-6):times[,i][5]   # IP=6days
  t_pred[,i] <- (times[,i][1]-6):times[,i][13] #uncomment for Clapham data 
  t0[i] <- times[,i][1]-9                       #initial time points of ODE 
}

#2. prepare data for either 'model_a_sim.stan' or 'model_b_sim.stan' program
J<- n_patients                          #number of patients
y_obs <- viremia[, 1:J]       #observed viremia measurements 
sample_time <- times[, 1:J]   #measurement times 
t_pred <- t_pred[, 1:J]
t0 <- t0[1:J]
is_censored <- is_censored[,1:J]
#serology <- serology[1:J]   #uncomment for model b simulation

sim_data <- list(J = J, #numb. of subjects
                 A = 1.4e6, omega = 1e4, sigma =0.5, eta = 4e4, #these are fixed values
                 n_obs = length(y_obs[,1]), #5,   #number of observed measurements per subject
                 n_pred = length(t_pred[,1]),    #number of predicted measurements (from generated quantities block)
                 t0 = t0,                     #starting point of ode 
                 ts = sample_time,           #time points where observed data are collected
                 t_pred =t_pred,          # time point where predictions are made from ode (in the generated quantities block)
                 y_hat = y_obs,
                 is_censored = is_censored,    #indicator matrix 
                 L = 357,                      #LOD: DENV1=357, 
                 #K = 2,    #number of groups (uncomment for model b simulation)
                 #serology = serology,       # uncomment for model b simulations
                 run_estimation=1)   #switch on the likelihood

#3. fit model a to DENV1 data
fit_model <- stan("model_a_sim.stan",   #change to 'model_b_sim.stan' for model b simulation         
                  data = sim_data,
                  chains = 4, iter = 2500,   warmup = 1000,  
                  control = list(adapt_delta = 0.8, max_treedepth = 10)
                  # seed=12345
)
save(fit_model,  sim_data,  file = "../results/DENV1_all_model_a_fit.RData", envir = .GlobalEnv)
#uncomment for the corresponding model simulation
#save(fit_model,  sim_data,  file = "../results/DENV1_all_model_b_fit.RData", envir = .GlobalEnv)
#save(fit_model,  sim_data,  file = "../results/Clapham_primary_fit.RData", envir = .GlobalEnv)
#save(fit_model,  sim_data,  file = "../results/Clapham_secondary_fit.RData", envir = .GlobalEnv)

#4. summarize results and check fit-
print(fit_model, pars = c("delta", "gamma", "kappa", "beta", "mu", 'std','lp__'), prob = c(0.025, 0.5, 0.975),digits=3)
#print(fit_model, pars = c("V0"), prob = c(0.025, 0.5, 0.975),digits=3)
traceplot(fit_model, pars = c("delta",  "gamma", "kappa", "beta","mu", 'std'), inc_warmup =TRUE, nrow=2)
traceplot(fit_model, pars = c("V0[1]",  "V0[2]", "V0[3]", "V0[4]","V0[5]", 'V0[6]'), inc_warmup =TRUE, nrow=2)
pairs(fit_model, pars = c("delta","gamma","kappa", "beta","mu", 'V0[1]', 'V0[2]'), las = 1)

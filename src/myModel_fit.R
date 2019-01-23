#' ---
#' title: "my proposed model"
#' author: "Aminath Shausan"
#' date: "December 20, 2018"
#' This R code is used to estimate parameters upto beta using "myModel-fit_upto_beta.stan"
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

#--------load patient data and extract observed viremia and meaured time------------------------------
#Use a function to load all patient data and get required patient data
Load_Patient_Data <- function(RData, env = new.env()){
  load(RData, env)
  p_data =  get("p1_data",pos=env) #change p1 to which patient data to select
  return(p_data)
}
p_data <- Load_Patient_Data('../data/patient_data.RData') #returns required patient data
rm(Load_Patient_Data) # remove the temporary environment to free up memory

sample_time  =  p_data$DOI  # viremia measurement times 
y_obs =  p_data$viremia     #viremia measurements


#---------Estimating delta, std, gamma, kappa, eta, omega, beta------------
#(date: 21/12/18)
stan_data_beta <- list(n_obs = 5,                 #number of observed measurements
                       n_pred = 11,                #number of predicted measurements (from generated quantities block)
                       y0 = c(1e8, 0, 0, 0, 1, 1),  #initial state of ode 
                       t0 = -4,                     #starting point of ode 
                       ts = sample_time,           #time points where observed data are collected
                       t_pred =-3:7,           # time point where predictions are made from ode (in the generated quantities block)
                       A= 1.4e6,  
                       sigma = 0.5,   
                       y_hat = y_obs,           #observed data
                       run_estimation=1)       #switch on likelihood 

# Test / debug the model:                                    
test <- stan("myModel-fit_upto_beta.stan",
             data = stan_data_beta,
             chains = 1, iter = 100, 
             control = list(adapt_delta = 0.8, max_treedepth = 10)) 
fit_model = stan(fit = test,
                 data = stan_data_beta,
                 chains = 4,     #number of Markov Chains
                 warmup = 1000, #number of warmup iterations per chain
                 iter = 2500,   #total number of iterations per chain
                 refresh = 1000,
                 control = list(adapt_delta = 0.8, max_treedepth = 10))
save.image(file = "p1_upto_beta.RData")
#-------------------------------------------------------


load("patient1_upto_beta.RData") #loads the saved data

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

#---------------use this for individual patient fit----------------------------------------------
#Use a function to load all patient data and get required patient data
Load_Patient_Data <- function(RData, env = new.env()){
  load(RData, env)
  p_data =  get("p1_data", pos=env) #change p1 to which patient data to select
  return(p_data)
}
p1_data <- Load_Patient_Data('../data/patient_data.RData') #returns required patient data
rm(Load_Patient_Data) # remove the temporary environment to free up memory

#create a matrix of viremia measurements 
J <- 2
y_obs <- matrix(data = c(p1_data$viremia, p2_data$viremia), 5, 2)
sample_time <- matrix(data = c(p1_data$DOI, p2_data$DOI), 5, 2)
t_pred <- matrix(data =c(-3:7, -3:7), 11,2)
#sample_time  =  as.matrix(p_data$DOI)  # viremia measurement times 
#y_obs =  as.matrix(p_data$viremia)    #viremia measurements


#---------Estimating delta, std, gamma, kappa, eta, omega, beta------------
#(date: 21/12/18)
stan_data <- list(J = J,                     #number of patients
                       n_obs = 5,                 #number of observed measurements
                       n_pred = 11,                #number of predicted measurements (from generated quantities block)
                       y0 = c(1e8, 0, 0, 0, 1, 1),  #initial state of ode 
                       t0 = -4,                     #starting point of ode 
                       ts = sample_time,           #time points where observed data are collected
                       t_pred =t_pred,    # time point where predictions are made from ode (in the generated quantities block)
                       A= 1.4e6,  
                       sigma = 0.5,   
                       y_hat = y_obs,           #observed data
                       run_estimation=1)       #switch on likelihood 

# Test / debug the model:                                    
test <- stan("myModel-fit_upto_beta.stan",
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

save.image(file = '../results/p1_to_p2_exponential.RData')
#-------------------------------------------------------


load("patient1_upto_beta.RData") #loads the saved data

# Author: Aminath Shausan 
# Date: May 2022
# This program simulates Experiment 1 or Experiment 2 studies for treating with TIPs
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
#------TIPs therapy (First and second simulation studies)------------------------
#first and second simulation study: administer TIPs on a particular day with specific effectiveness

#define the ode
ode_func <- function(t, y, params) {
  with(as.list(c(params, y)), {
    
    dS = A- beta * y[5] * y[1] - beta * y[6] * y[1];                  #Susceptible cells ;
    
    dI_V = beta * y[1] * y[5] - beta * y[6] * y[2]- delta * y[2] -mu*y[2];      #virus only infected cells 
    
    dI_T = beta * y[6] * y[1] - beta * y[5] * y[3] ;                   #TIP only infected  cells 
    
    dI_B = beta * y[6] * y[2] +  beta * y[5] * y[3] - gamma * y[4]+mu*y[2];     #Co-infected cells
    
    dV = omega * y[2] + omega * (1-sigma) * y[4] - kappa * y[5];          #free Virus
    
    dT = eta * sigma* y[4] - kappa * y[6];                              #free TIPs
    
    res1 <- c(dS,dI_V ,dI_T,dI_B,dV,dT)
    list(res1)
  })
}

#these are median posterior estimates for DENV1 from model a
params_median <- list(A =1.4e6, omega = 1e4, sigma = 0.5, eta = 4e4, 
                      beta = 1.283e-10, delta = 8.844,  gamma = 6.585,  kappa = 8.855,  mu = 5.13e-4)

#initial conditions with no TIPs and with TIPs
inits_noTIPs<- c(S0=1e8 , I_V0=0, I_T0=0, I_B0=0, V0=55.29, T0=0)

times_noTIPs = -10:5


# Run the integration and Store the output in a data frame:
#1. no TIPS
out_noTIPs <- ode(inits_noTIPs, times_noTIPs, ode_func, params_median, method="ode45")
out_noTIPs <- data.frame(out_noTIPs)
colnames(out_noTIPs) <- c("time", "S", "I_V", "I_T", "I_B", "V", "T")
out_noTIPs <- as.data.frame(out_noTIPs)

#change sigma = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
params_sigma <- list(A =1.4e6, omega = 1e4, sigma = 1.0, eta = 4e4, 
                     beta = 1.283e-10, delta = 8.844,  gamma = 6.585,  kappa = 8.855,  mu = 5.13e-4)

#change value of i corresponding to the specific Fever Day
i = 7; # (3=FeverDay -8, 5=FeverDay -6, 7=FeverDay -4, 9 = FeverDay -2, 11 =FeverDay 0)
inits_1e8TIPs <- c(S0=out_noTIPs[i,2], I_V0=out_noTIPs[i,3], I_T0=out_noTIPs[i,4], I_B0=out_noTIPs[i,5], V0=out_noTIPs[i,6], T0=1e8+out_noTIPs[i,7])
inits_1e10TIPs <- c(S0=out_noTIPs[i,2], I_V0=out_noTIPs[i,3], I_T0=out_noTIPs[i,4], I_B0=out_noTIPs[i,5], V0=out_noTIPs[i,6], T0=1e10+out_noTIPs[i,7])
inits_1e12TIPs <- c(S0=out_noTIPs[i,2], I_V0=out_noTIPs[i,3], I_T0=out_noTIPs[i,4], I_B0=out_noTIPs[i,5], V0=out_noTIPs[i,6], T0=1e12+out_noTIPs[i,7])

times_withTIPs <- -4:5 #change this to -8:5, -6:5, -4:5, -2:5, 0:5
#2. 1e8 TIPs
out_1e8TIPs <- ode(inits_1e8TIPs, times_withTIPs, ode_func, params_sigma, method="ode45")
out_1e8TIPs <- data.frame(out_1e8TIPs)
colnames(out_1e8TIPs) <- c("time", "S", "I_V", "I_T", "I_B", "V", "T")
out_1e8TIPs <- as.data.frame(out_1e8TIPs)


#3. 1e10 TIPs
out_1e10TIPs <- ode(inits_1e10TIPs, times_withTIPs, ode_func, params_sigma, method="ode45")
out_1e10TIPs <- data.frame(out_1e10TIPs)
colnames(out_1e10TIPs) <- c("time", "S", "I_V", "I_T", "I_B", "V", "T")
out_1e10TIPs <- as.data.frame(out_1e10TIPs)

#4. 1e12 TIPs
out_1e12TIPs <- ode(inits_1e12TIPs, times_withTIPs, ode_func, params_sigma, method="ode45")
out_1e12TIPs <- data.frame(out_1e12TIPs)
colnames(out_1e12TIPs) <- c("time", "S", "I_V", "I_T", "I_B", "V", "T")
out_1e12TIPs <- as.data.frame(out_1e12TIPs)

#plot subplots in Figs 5 to 7 and figs S3 to S4
ggplot(out_noTIPs, aes(x=time, y=log10(V))) + ylim(-1, 13)+
  geom_line(col="red") +  geom_hline(yintercept=log10(357), color = "black")+ 
  geom_line(data = out_1e8TIPs, aes(x=time, y=log10(V)), color = "magenta")+
  geom_line(data = out_1e8TIPs, aes(x=time, y=log10(T)), color = "magenta", linetype=2)+
  geom_line(data = out_1e10TIPs, aes(x=time, y=log10(V)), color = "orange") +
  geom_line(data = out_1e10TIPs, aes(x=time, y=log10(T)), color = "orange", linetype= 2) +
  geom_line(data = out_1e12TIPs, aes(x=time, y=log10(V)), color = "green") +
  geom_line(data = out_1e12TIPs, aes(x=time, y=log10(T)), color = "green", linetype= 2)+
  geom_line(data = out_noTIPs, aes(x=time, y=log10(V)), color = "red")+
  labs(x = "Fever Day", y = "log10(titre)")

#For first simulation study: change fever_day-number and effec_number according to the feverday and sigma value used and 
#save all results. plots can be recovered from the .RData
save.image(file = '../results/first_sim_study/fever_day-4_effec_100.RData')

#For second simulation study: change fever_day-number and effec_number according to the feverday and sigma value used and 
#save all results. plots can be recovered from the .RData
save.image(file = '../results/second_sim_study/fever_day-4_effec_100.RData')
#=============================================

#-------Individual patient therapy  (S6 Fig to S16 Fig) --------------------------------------------
#use TIPs therapy for each patient, for a particular value of sigma and for a given fever day 
#use doses of: 1e8, 1e10 and 1e12
#use model a fit

rm(list=ls())  #clear workspace

load('../results/DENV1_all_model_a_fit.RData')  #load model a fit

y_pred <- extract(fit_model, pars= c('y_pred_S','y_pred_IV', 'y_pred_IT', 'y_pred_IB','y_pred_viremia', 'y_pred_T'))

J <- 67

ode_func <- function(t, y, params) {
  with(as.list(c(params, y)), {
    
    dS = A- beta * y[5] * y[1] - beta * y[6] * y[1];                  #Susceptible cells ;
    
    dI_V = beta * y[1] * y[5] - beta * y[6] * y[2]- delta * y[2] -mu*y[2];      #virus only infected cells 
    
    dI_T = beta * y[6] * y[1] - beta * y[5] * y[3] ;                   #TIP only infected  cells 
    
    dI_B = beta * y[6] * y[2] +  beta * y[5] * y[3] - gamma * y[4]+mu*y[2];     #Co-infected cells
    
    dV = omega * y[2] + omega * (1-sigma) * y[4] - kappa * y[5];          #free Virus
    
    dT = eta * sigma* y[4] - kappa * y[6];                              #free TIPs
    
    res1 <- c(dS,dI_V ,dI_T,dI_B,dV,dT)
    list(res1)
  })
}

#t_pred_onset <- matrix(data =NA, ncol =J, nrow = 5);    #predicted times,   observed times                   
#t_pred_prior <- matrix(data =NA, ncol =J, nrow = 9); #predicted times, 4 days prior to 1st day of onset 
t_pred_post <- matrix(data =NA, ncol =J, nrow = 3); #predicted times, 2 days after the 1st day of onset 

for (j in 1:J) {
#  t_pred_onset[,j] <- (sim_data$ts[,j][1]):sim_data$ts[,j][5]   #onset
#  t_pred_prior[,j] <- (sim_data$ts[,j][1]-4):sim_data$ts[,j][5]   #add 4 days prior to onset  
  t_pred_post[,j] <- (sim_data$ts[,j][1]+2):sim_data$ts[,j][5]   #deduce 2 days post  onset 
}

#change name of 'ode_output_sigma_time_dose' according to sigma time and TIPs dose values 
#ode_output_0.8_onset_1e8 <-  matrix(data =NA, ncol =J, nrow = length(t_pred_onset[,1])) 
#ode_output_0.8_prior_4days_1e8 <-  matrix(data =NA, ncol =J, nrow = length(t_pred_prior[,1])) 
ode_output_0.5_post_2days_1e12 <-  matrix(data =NA, ncol =J, nrow = length(t_pred_post[,1])) 

for (j in 1:J){
  #change sigma value as required (0.5,0.6,0.7,0.8,1.0)
  params <- list(A =1.4e6, omega = 1e4, sigma = 0.5, eta = 4e4, 
                 beta = 1.283e-10, delta = 8.844,  gamma = 6.585,  kappa = 8.855,  mu = 5.13e-4)
  
  
  cred_S = apply(y_pred$y_pred_S[ , ,j], 2, median)
  cred_IV = apply(y_pred$y_pred_IV[ , ,j], 2, median)
  cred_IT = apply(y_pred$y_pred_IT[ , ,j], 2, median)
  cred_IB = apply(y_pred$y_pred_IB[ , ,j], 2, median)
  cred_viremia = apply(y_pred$y_pred_viremia[ , ,j], 2, median)
  cred_T = apply(y_pred$y_pred_T[ , ,j], 2, median)
  
  #change TIPs dose at T0 from 1e8, 1e10, 1e12
 # inits_onset <- c(S0=cred_S[7] , I_V0=cred_IV[7], I_T0=cred_IT[7],   #use these inits for  symptoms onset
#             I_B0=cred_IB[7], V0=cred_viremia[7], T0=cred_T[7]+1e8)
 
 #  inits_prior <- c(S0=cred_S[3] , I_V0=cred_IV[3], I_T0=cred_IT[3],   #use these inits for 4days prior onset
#             I_B0=cred_IB[3], V0=cred_viremia[3], T0=cred_T[3]+1e8)
  
  inits_post <- c(S0=cred_S[9] , I_V0=cred_IV[9], I_T0=cred_IT[9],   #use these inits for 2days post onset
            I_B0=cred_IB[9], V0=cred_viremia[9], T0=cred_T[9]+1e12)
  
  #out_onset <- ode(inits_onset, t_pred_onset[,j], ode_func, params, method="ode45")
  #out_prior <- ode(inits_prior, t_pred_prior[,j], ode_func, params, method="ode45")
  out_post <- ode(inits_post, t_pred_post[,j], ode_func, params, method="ode45")
  #ode_output_0.8_onset_1e8[,j] <- out_onset[,6] 
  #ode_output_0.8_prior_4days_1e8[,j] <- out_prior[,6] 
  ode_output_0.5_post_2days_1e12[,j] <- out_post[,6] 
 
}
#change name 'ode_output_sigma_onset_dose' according to sigma and dose used
#save(ode_output_0.8_onset_1e8,  t_pred_onset, fit_model, sim_data, file = "../results/individual_therapy/ode_output_0.8_onset_1e8.RData", envir = .GlobalEnv)
#save(ode_output_0.8_prior_4days_1e8,  t_pred_prior, fit_model, sim_data, file = "../results/individual_therapy/ode_output_0.8_prior_4days_1e8.RData", envir = .GlobalEnv)
save(ode_output_0.5_post_2days_1e12,  t_pred_post, fit_model, sim_data, file = "../results/individual_therapy/ode_output_0.5_post_2days_1e12.RData", envir = .GlobalEnv)


################################################
#compute trajectory for each patient with TIPs therapy 
#change name of ode_output as required for the plot
rm(list=ls())
e1 <- new.env(parent = emptyenv())

#load("../results/individual_therapy/ode_output_1.0_onset_1e8.RData", envir = e1)
#load("../results/individual_therapy/ode_output_1.0_onset_1e10.RData", envir = e1)
#load("../results/individual_therapy/ode_output_1.0_onset_1e12.RData", envir = e1)

#load("../results/individual_therapy/ode_output_1.0_prior_4days_1e8.RData", envir = e1)
#load("../results/individual_therapy/ode_output_1.0_prior_4days_1e10.RData", envir = e1)
#load("../results/individual_therapy/ode_output_1.0_prior_4days_1e12.RData", envir = e1)

load("../results/individual_therapy/ode_output_0.5_post_2days_1e8.RData", envir = e1)
load("../results/individual_therapy/ode_output_0.5_post_2days_1e10.RData", envir = e1)
load("../results/individual_therapy/ode_output_0.5_post_2days_1e12.RData", envir = e1)


#ode_output_1e8 <- e1$ode_output_0.7_onset_1e8
#ode_output_1e10 <- e1$ode_output_0.7_onset_1e10
#ode_output_1e12 <- e1$ode_output_0.7_onset_1e12

#ode_output_1e8 <- e1$ode_output_1.0_prior_4days_1e8
#ode_output_1e10 <- e1$ode_output_1.0_prior_4days_1e10
#ode_output_1e12 <- e1$ode_output_1.0_prior_4days_1e12

ode_output_1e8 <- e1$ode_output_0.5_post_2days_1e8
ode_output_1e10 <- e1$ode_output_0.5_post_2days_1e10
ode_output_1e12 <- e1$ode_output_0.5_post_2days_1e12

#t_pred <- e1$t_pred_onset
#t_pred <- e1$t_pred_prior
t_pred <- e1$t_pred_post
sim_data <-e1$sim_data


#params2 <- #e1$y_pred  #no such file
y_pred <-  e1$fit_model #extract model a fit
y_pred <- extract(y_pred, pars= c('y_pred_viremia')) #viremia measurements from fit

#change value of j corresponding to the 1st 5 primary and 1st 5 secondary patients
#j = 3, 5, 8, 9, 12 are primary patients
#j= 1,2, 4, 6, 7, 12 are secondary patients
prim_pat <- list(3, 5, 8, 9, 12 )

for (j in prim_pat) {
  
  df_fit = data.frame(cred_median_viremia = apply(y_pred$y_pred_viremia[ , ,j], 2, median),  Time = sim_data$t_pred[,j])
  df_1e8 <- data.frame( cred_median_1e8 = ode_output_1e8[,  j], Time =  t_pred[,j])     #
  df_1e10 <- data.frame(cred_median_1e10 =ode_output_1e10[,  j], Time =  t_pred[,j])    #
  df_1e12 <- data.frame(cred_median_1e12 =ode_output_1e12[,  j], Time =  t_pred[,j])         #
  
   pl = ggplot(df_fit, aes(x=Time, y=log10(cred_median_viremia))) + 
           geom_hline(yintercept=log10(357), color = "black")+  ylim(-1, 10)+
           geom_line(data = df_fit, aes(x=Time, y=log10(cred_median_viremia)), color = "blue") +
           geom_line(data = df_1e8, aes(x=Time, y=log10(cred_median_1e8)), color = "magenta") +
           geom_line(data = df_1e10, aes(x=Time, y=log10(cred_median_1e10)), color = "orange") +
           geom_line(data = df_1e12, aes(x=Time, y=log10(cred_median_1e12)), color = "green") +
           labs(x = "Fever Day", y = "log10(titre)") 
   #change sigma value as required
   ggsave(pl, file= paste0("../images/individual_therapy/post/0.5_primary_p", j,".jpeg"), width = 4.5, height = 3.75)
   print(pl)
}

#clear all plots
dev.off(dev.list()["RStudioGD"]) 

sec_pat <- list( 1,2, 4, 6, 7, 12)
for (j in sec_pat) {#j =8;
  
  df_fit = data.frame(cred_median_viremia = apply(y_pred$y_pred_viremia[ , ,j], 2, median),  Time = sim_data$t_pred[,j])
  df_1e8 <- data.frame( cred_median_1e8 = ode_output_1e8[,  j], Time =  t_pred[,j])     #
  df_1e10 <- data.frame(cred_median_1e10 =ode_output_1e10[,  j], Time =  t_pred[,j])    #
  df_1e12 <- data.frame(cred_median_1e12 =ode_output_1e12[,  j], Time =  t_pred[,j])         #
  
  pl = ggplot(df_fit, aes(x=Time, y=log10(cred_median_viremia))) + 
    geom_hline(yintercept=log10(357), color = "black")+  ylim(-1, 10)+
    geom_line(data = df_fit, aes(x=Time, y=log10(cred_median_viremia)), color = "blue") +
    geom_line(data = df_1e8, aes(x=Time, y=log10(cred_median_1e8)), color = "magenta") +
    geom_line(data = df_1e10, aes(x=Time, y=log10(cred_median_1e10)), color = "orange") +
    geom_line(data = df_1e12, aes(x=Time, y=log10(cred_median_1e12)), color = "green") +
    labs(x = "Fever Day", y = "log10(titre)") 
  ggsave(pl, file= paste0("../images/individual_therapy/post/0.5_secondary_p", j,".jpeg"), width = 4.5, height = 3.75)
  print(pl)
}


# Author: Aminath Shausan 
# Date: May 2022
# This program analyzes results from either 'model a' or 'model b' fit 
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

#===(Table 3, Table 4, S1 Table) summarize parameters for model  fits====
load('../results/DENV1_all_model_a_fit.RData')  #load model a fit
#load('../results/DENV1_all_model_b_fit.RData')  #uncomment for model b fit
#load('../results/Clapham_primary_fit.RData')  #uncomment for Clampham_DENV1_primary fit
#load('../results/Clapham_secondary_fit.RData')  #uncomment for Clampham_DENV1_secondary fit

#print(fit_model, pars=c("delta", "gamma")) #also works

#summarize fixed effect parameters
fixed_effect_params = as.data.frame(fit_model, pars = c("delta", "gamma", "kappa", "beta", "mu", "std"))
fixed_eff_median = apply(fixed_effect_params, 2, median)
#fixed_eff_mean = apply(fixed_effect_params, 2, mean)
fixed_effect_low = apply(fixed_effect_params, 2,quantile, probs = c(0.025)) #lower value of Credible Interval
fixed_effect_high = apply(fixed_effect_params, 2,quantile, probs = c(0.975))#upper value of Credible Interval

#quantile(fixed_effect_params$mu, 0.975) # also works

#model a fit: summarize V0: (V0 is patient specific, so take median of the medians) 
V0_post = as.data.frame(fit_model, pars = c("V0"))
V0_medians = apply(V0_post, 2, median) #this is a named numbered list
summary(V0_medians) #give min, max, IQR and median of delta_medians

#model b fit: summarise V0 for primary and secondary separately patients 
#V0 is patient specific (take median of the medians)
V0_post = as.data.frame(fit_model, pars = c("V0"))
idx_prim <- which(serology == '1')
idx_sec <- which(serology == '2')
V0_prim <- V0_post[, idx_prim]
V0_sec <- V0_post[, idx_sec]  
V0_medians_prim = apply(V0_prim, 2, median) #this is a named numbered list
V0_medians_sec = apply(V0_sec, 2, median)
summary(V0_medians_prim) #give min, max, IQR and median of delta_medians
summary(V0_medians_sec)
##########################################
#  plot 95% posterior credible intervals for each patient
#(Fig 3, S1 Fig. S2 Fig, S5 Fig. subplots are from DENV1_all_model_a_fit.RData)
#(Fig 4 subplots are from DENV1_all_model_b_fit.RData)

post <- extract(fit_model, pars= c('y_pred_S',  'y_pred_viremia', 'y_pred_T'))

#change value of i = 1:67
# i = { 3  5  8  9 12 13 14 16 17 20 22 23 24 28 29 31 32 40 42 43 46 51 56 59 62} are primary patients
# i ={1  2  4  6  7 10 11 15 18 19 21 25 26 27 30 33 34 35 36 37 38 39 41 44 45 47 48 49 50 52 
  # 53 54 55 57 58 60 61 63 64 65 66 67} are secondary patients
{i=1;
  
  cred_median_S = apply(post$y_pred_S[, ,i], 2, median)               #susceptible cells
  cred_low_S=  apply(post$y_pred_S[, ,i], 2, quantile, probs = c(0.025))
  cred_high_S =  apply(post$y_pred_S[, ,i], 2, quantile, probs = c(0.975))
  
  cred_median_viremia = apply(post$y_pred_viremia[, ,i], 2, median)     #free wT-virus
  cred_low_viremia=  apply(post$y_pred_viremia[, ,i], 2, quantile, probs = c(0.025))
  cred_high_viremia =  apply(post$y_pred_viremia[, ,i], 2, quantile, probs = c(0.975))
  
  cred_median_T = apply(post$y_pred_T[, ,i], 2, median)               #TIPs
  cred_low_T=  apply(post$y_pred_T[, ,i], 2, quantile, probs = c(0.025))
  cred_high_T =  apply(post$y_pred_T[, ,i], 2, quantile, probs = c(0.975))
  
  df_fit = data.frame(cred_low_S, cred_median_S, cred_high_S,
                      cred_low_viremia, cred_median_viremia, cred_high_viremia,
                      cred_low_T, cred_median_T, cred_high_T,
                      Time = sim_data$t_pred[,i]) 
  
  df_sample =  data.frame(Time = sim_data$ts[,i],  viremia = sim_data$y_hat[,i]) #observed data
  
  print(ggplot(df_sample, aes(x=Time, y=log10(viremia))) +
          geom_point(col="black", shape = 19, size = 1.5)+
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_median_viremia)), color = "blue") +
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_high_viremia)), color = "blue", linetype=3) +
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_low_viremia)), color = "blue", linetype=3) +
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_median_T)), color = "green") +
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_high_T)), color = "green", linetype=3) +
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_low_T)), color = "green", linetype=3) +
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_median_S)), color = "purple") +
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_high_S)), color = "purple", linetype=3) +
          geom_line(data = df_fit, aes(x=Time, y=log10(cred_low_S)), color = "purple", linetype=3)+ 
          # Aesthetics
          labs(x = "Fever Day", y = "log10(Viremia)") )
  
}

#=====prepare the S-by-J matrix of loglikelihoods to be used in  waic() function============
# ---- Table 2 results ------------

#S= (total number of MCMC iterations)
#J = (number of patients)
#this structure of loglikelihood used with loo() implies we are leaving out one patient
#use the following code with either model a or  b fits 
post2 <-  extract(fit_model, pars =c("std", 'y_new'))

n_iter <- length(post2$std)#6000#10000     #number of MCMC post-warmup iterations 
J=sim_data$J
n_obs <-5
log_like_ind <- matrix(data =NA, ncol =J, nrow = 5); #a 5-by-J matrix of log_likelihood for each iteration
log_lik_a <- matrix(data =NA, ncol =J, nrow = n_iter);#log_likelihood matrix of size n_iter x J
#uncomment for model b fit
#log_lik_b <- matrix(data =NA, ncol =J, nrow = n_iter);#log_likelihood matrix of size n_iter x J
for(i in 1:n_iter){
  for(j in 1:J){
    for(n in 1:n_obs){
      if(is_censored[n,j]==0){
        log_like_ind[n,j]<-    dlnorm(sim_data$y_hat[n,j], meanlog = log(post2$y_new[i,,j][n]), sdlog = post2$std[i], log = TRUE)
      }else{
        log_like_ind[n,j]<-  plnorm(sim_data$y_hat[n,j], meanlog = log(post2$y_new[i,,j][n]), sdlog = post2$std[i], lower.tail=TRUE, log.p = TRUE)
      }
    }
  } 
   log_lik_a[i,] <- colSums(log_like_ind)  
  #log_lik_b[i,] <- colSums(log_like_ind)  #uncomment for model b fit
}

#save the log_like matrices so that they can be used to compute waic.
# log_lik_a is the log-likelihood of DENV1 model a
# log_lik_b is the log-likelihood of DENV1 model b

save(log_lik_a,   file = "../results/logLike_a.RData", envir = .GlobalEnv)
#save(log_lik_b,   file = "../results/logLike_b.RData", envir = .GlobalEnv) #unomment for model b fit

####### Compute WAIC (Table 2 result ) ########
waic_model_a <- waic(log_lik_a)
#waic_model_b <- waic(log_lik_b)  #uncomment for model b fit
waic_diff <- compare(waic_model_a, waic_model_b) 

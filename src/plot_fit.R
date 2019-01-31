#' ---
#' title: "Extract Stan object and plot relavent statistics"
#' author: "Aminath Shausan"
#' date: "December 20, 2018"
#' ---
#clear history 
rm(list=ls())

#==================================================================
#This program loads each patient's Stan fit and extracts posterior draws for the estimated parameters 
#==================================================================
#save the packages in this library
.libPaths("C:/Users/shausan/r-libraries") #use this on qut desktop
.libPaths("C:/Users/shau/r-libraries") #use this on laptop
#---------------------------------------------------------------
#load required packages 
library("bayesplot")
library("ggplot2")
library(dplyr)
install.packages("shinystan")
library("shinystan")
#=======================================================
#load simulated data, which are computed from prior distribution
p_fit <- load('p4_upto_beta.RData')

#=====================================================================
#Inference for fit using real data

print(fit_model, pars = c("delta", "std","gamma", "kappa","eta", "omega", "beta", "lp__"), digits=3)
traceplot(fit_model, pars = c("delta", "std", "gamma", "kappa","eta", "omega", "beta"), inc_warmup =TRUE, nrow=4)
pairs(fit_model, pars = c("delta","std","gamma", "eta","kappa", "eta", "omega", "beta"), las = 1)

#compute credible interval and plot them for "fake_fit"
params <- extract(fit_model)

cred_median = apply(params$y_ppc[, ,2], 2, median)
cred_low =  apply(params$y_ppc[, ,2], 2, quantile, probs = c(0.025))
cred_high =  apply(params$y_ppc[, ,2], 2, quantile, probs = c(0.975))

cred_median_ode = apply(params$y_pred[, , 5], 2, median)
cred_low_ode=  apply(params$y_pred[, , 5], 2, quantile, probs = c(0.025))
cred_high_ode =  apply(params$y_pred[, , 5], 2, quantile, probs = c(0.975))


df_sample =  data.frame(sample_time[,2],  y_obs[,2]) #observed data
colnames(df_sample) = c("Time", "Viremia")
#df_fit = data.frame(cred_low, cred_median, cred_high, cred_low_ode, cred_median_ode, cred_high_ode, Time = stan_data_beta$t_pred) #predicted credible interval
df_fit = data.frame(cred_low, cred_median, cred_high,  Time = stan_data_beta$t_pred[,2]) #predicted credible interval


#plot  posterior predictive interval for fake data
ggplot(df_sample, aes(x=Time, y=log10(Viremia))) +
  geom_point(col="black", shape = 19, size = 1.5) +
  # Error in integration:
  geom_line(data = df_fit, aes(x=Time, y=log10(cred_median)), color = "red") +
  geom_line(data = df_fit, aes(x=Time, y=log10(cred_high)), color = "red", linetype=3) +
  geom_line(data = df_fit, aes(x=Time, y=log10(cred_low)), color = "red", linetype=3) +
 # geom_line(data = df_fit, aes(x=Time, y=log10(cred_median_ode)), color = "blue") +
#  geom_line(data = df_fit, aes(x=Time, y=log10(cred_high_ode)), color = "blue", linetype=3) +
#  geom_line(data = df_fit, aes(x=Time, y=log10(cred_low_ode)), color = "blue", linetype=3) +
  # Aesthetics
  labs(x = "Day of illness", y = "log10(Viremia)") +
  scale_x_continuous(limits=c(-3, 8)) #+
#scale_y_continuous(limits=c(0,20)) 
#===========shiny stan plots==============
launch_shinystan(fit_model)

#=======================================
#plot posterior and prior density for each parameter for each patient 
#1. delta patient 1 
x <- seq(0, to = 30, by = 1)
delta_prior <- dexp(x, rate = 1/5)
plot(density(params$delta), xlab = "x", main= "delta")
lines(x, delta_prior, type = "l", col ="red")

#2. std patient 1 (plot is okay)
std_prior <- dexp(x, rate = 1)
plot(density(params$std), xlab = "x", main= "std")
lines(x, std_prior, type = "l", col ="red")

#3. gamma patient 1 (plot is okay)
gamma_prior <- dexp(x, rate = 1/5)
plot(density(params$gamma), xlab = "x", main= "gamma")
lines(x, gamma_prior, type = "l", col ="red")

#4. kappa  patient 1 (plot is okay)
kappa_prior <- dexp(x, rate = 1/5)
plot(density(params$kappa), xlab = "x", main= "kappa")
lines(x, kappa_prior, type = "l", col ="red")

#5. eta  patient 1 (plot is okay)
x_eta <- seq(0, 5e5, by = 100)
eta_prior <- dlnorm(x_eta, mean=9.9, sd = 1)
plot(density(params$eta), xlab = "x", main= "eta")
lines(x_eta, eta_prior, type = "l", col ="red")

#6. omega  patient 1 (plot is okay)
x_omega <- seq(0, 5e5, by = 100)
omega_prior <- dlnorm(x_omega, mean=9.2, sd = 1)
plot(density(params$omega), xlab = "x", main= "omega",  ylim=c(0,7e-05))
lines(x_omega, omega_prior, type = "l", col ="red")

#7. beta  patient 1 (plot is okay)
x_beta <- seq(0, 6e-10, by = 1e-11)
beta_prior <- dlnorm(x_beta, mean= -24, sd = 1)
plot(density(params$beta), xlab = "x", main= "beta",  ylim=c(0,1.8e+10))
lines(x_beta, beta_prior, type = "l", col ="red")
dev.off() #clears current figure

x_beta <- seq(0, 5e-10, by = 1e-11)
beta_raw_prior <- rnorm(x_beta, mean= 0, sd = 1)
beta_prior <-exp(-24 +1*beta_raw_prior)
plot(density(p1_post$beta), xlab = "x", main= NA)
lines(x_beta, beta_prior, type = "l", col ="red")
dev.off() #clears current figure
#--------------------------------------------
#check if prior is informative; condition is [post sd (for any parameter) > 0.1 * prior sd , then prior is informative]
#change patient number accordingly
sd(p10_post$delta) > 0.1* 5
sd(p10_post$std) > 0.1* 1
sd(p10_post$gamma) > 0.1* 5
sd(p10_post$kappa) > 0.1* 5
sd(log(p10_post$eta)) > 0.1* 1 #this comparision is on the log scale
sd(log(p10_post$omega)) > 0.1* 1 #this comparision is on the log scale


#------------------------------------------------
ggplot(data  = p1_post, aes(x = delta)) + #patient 1
  geom_density() + 
  stat_function(fun=dexp,geom = "line",size=1,col="blue",args = (mean=5)) #this is not good plot
# Aesthetics
labs(x = "delta", y = "density")

#=========plot post density for all patients ==============
#Use a function to extract posterior draws for each parameter
LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  p_fit =  get("fit_model",pos=env)
  p_params = as.data.frame(p_fit, pars = c("delta", "std", "gamma", "kappa", "eta", "omega", "beta"))
  return(p_params)
  # return(env) 
}
p1_post <- LoadToEnvironment('p1_upto_beta.RData')
p2_post <- LoadToEnvironment('p2_upto_beta.RData')
p3_post <- LoadToEnvironment('p3_upto_beta.RData')
p4_post <- LoadToEnvironment('p4_upto_beta.RData')
p5_post <- LoadToEnvironment('p5_upto_beta.RData')
p6_post <- LoadToEnvironment('p6_upto_beta.RData')
p7_post <- LoadToEnvironment('p7_upto_beta.RData')
p8_post <- LoadToEnvironment('p8_upto_beta.RData')
p9_post <- LoadToEnvironment('p9_upto_beta.RData')
p10_post <- LoadToEnvironment('p10_upto_beta.RData')

rm(LoadToEnvironment) # remove the temporary environment to free up memory


#plot posterior density of delta for all patients
ggplot(data  = p1_post, aes(x = delta)) + #patient 1
  geom_density() + 
  geom_density(data = p2_post, color = "purple") + # p2
  geom_density(data = p3_post, color = "orange") + # p3
  geom_density(data = p4_post, color = "blue") + # p4
  geom_density(data = p5_post, color = "red") + # p5
  geom_density(data = p6_post, color = "green") + # p6
  geom_density(data = p7_post, color = "maroon") + # p7
  geom_density(data = p8_post, color = "brown") + #p8
  geom_density(data = p9_post, color = "azure") + # p9
  geom_density(data = p10_post, color = "coral") + # p10
  # Aesthetics
  labs(x = "delta", y = "density")
#posterior delta with prior 
x <- seq(0, to = 30, by = 1)
delta_prior <- dexp(x, rate = 1/5)
plot(density(p1_post$delta), xlab = "x", main= "delta", ylim =c(0,0.3))
lines(density(p2_post$delta), col ="purple")
lines(density(p3_post$delta), col ="orange")
lines(density(p4_post$delta), col ="blue")
lines(density(p5_post$delta), col ="red")
lines(density(p6_post$delta), col ="green")
lines(density(p7_post$delta), col ="maroon")
lines(density(p8_post$delta), col ="brown")
lines(density(p9_post$delta), col ="azure")
lines(density(p10_post$delta), col ="coral")
lines(x, delta_prior, type = "o", col ="red")

#plot posterior density of std for all patients
ggplot(data  = p1_post, aes(x = std)) + #patient 1
  geom_density() + 
  geom_density(data = p2_post, color = "purple") + # p2
  geom_density(data = p3_post, color = "orange") + # p3
  geom_density(data = p4_post, color = "blue") + # p4
  geom_density(data = p5_post, color = "red") + # p5
  geom_density(data = p6_post, color = "green") + # p6
  geom_density(data = p7_post, color = "maroon") + # p7
  geom_density(data = p8_post, color = "brown") + #p8
  geom_density(data = p9_post, color = "azure") + # p9
  geom_density(data = p10_post, color = "coral") + # p10
  # Aesthetics
  labs(x = "std", y = "density") 

x <- seq(0, to = 30, by = 1)
std_prior <- dexp(x, rate = 1)
plot(density(p1_post$std), xlab = "x", main= "std", ylim =c(0,1.5))
lines(density(p2_post$std), col ="purple")
lines(density(p3_post$std), col ="orange")
lines(density(p4_post$std), col ="blue")
lines(density(p5_post$std), col ="red")
lines(density(p6_post$std), col ="green")
lines(density(p7_post$std), col ="maroon")
lines(density(p8_post$std), col ="brown")
lines(density(p9_post$std), col ="azure")
lines(density(p10_post$std), col ="coral")
lines(x, std_prior, type = "o", col ="red")

#plot posterior density of gamma for all patients
ggplot(data  = p1_post, aes(x = gamma)) + #patient 1
  geom_density() + 
  geom_density(data = p2_post, color = "purple") + # p2
  geom_density(data = p3_post, color = "orange") + # p3
  geom_density(data = p4_post, color = "blue") + # p4
  geom_density(data = p5_post, color = "red") + # p5
  geom_density(data = p6_post, color = "green") + # p6
  geom_density(data = p7_post, color = "maroon") + # p7
  geom_density(data = p8_post, color = "brown") + #p8
  geom_density(data = p9_post, color = "azure") + # p9
  geom_density(data = p10_post, color = "coral") + # p10
  # Aesthetics
  labs(x = "gamma", y = "density") 

x <- seq(0, to = 30, by = 1)
gamma_prior <- dexp(x, rate = 1/5)
plot(density(p1_post$gamma), xlab = "x", main= "gamma", ylim =c(0,0.2))
lines(density(p2_post$gamma), col ="purple")
lines(density(p3_post$gamma), col ="orange")
lines(density(p4_post$gamma), col ="blue")
lines(density(p5_post$gamma), col ="red")
lines(density(p6_post$gamma), col ="green")
lines(density(p7_post$gamma), col ="maroon")
lines(density(p8_post$gamma), col ="brown")
lines(density(p9_post$gamma), col ="azure")
lines(density(p10_post$gamma), col ="coral")
lines(x, gamma_prior, type = "o", col ="red")

#plot posterior density of kappa for all patients
ggplot(data  = p1_post, aes(x = kappa)) + #patient 1
  geom_density() + 
  geom_density(data = p2_post, color = "purple") + # p2
  geom_density(data = p3_post, color = "orange") + # p3
  geom_density(data = p4_post, color = "blue") + # p4
  geom_density(data = p5_post, color = "red") + # p5
  geom_density(data = p6_post, color = "green") + # p6
  geom_density(data = p7_post, color = "maroon") + # p7
  geom_density(data = p8_post, color = "brown") + #p8
  geom_density(data = p9_post, color = "azure") + # p9
  geom_density(data = p10_post, color = "coral") + # p10
  # Aesthetics
  labs(x = "kappa", y = "density") 

x <- seq(0, to = 30, by = 1)
kappa_prior <- dexp(x, rate = 1/5)
plot(density(p1_post$kappa), xlab = "x", main= "kappa", ylim =c(0,0.31))
lines(density(p2_post$kappa), col ="purple")
lines(density(p3_post$kappa), col ="orange")
lines(density(p4_post$kappa), col ="blue")
lines(density(p5_post$kappa), col ="red")
lines(density(p6_post$kappa), col ="green")
lines(density(p7_post$kappa), col ="maroon")
lines(density(p8_post$kappa), col ="brown")
lines(density(p9_post$kappa), col ="azure")
lines(density(p10_post$kappa), col ="coral")
lines(x, kappa_prior, type = "o", col ="red")


#plot posterior density of eta for all patients
ggplot(data  = p1_post, aes(x = eta)) + #patient 1
  geom_density() + 
  geom_density(data = p2_post, color = "purple") + # p2
  geom_density(data = p3_post, color = "orange") + # p3
  geom_density(data = p4_post, color = "blue") + # p4
  geom_density(data = p5_post, color = "red") + # p5
  geom_density(data = p6_post, color = "green") + # p6
  geom_density(data = p7_post, color = "maroon") + # p7
  geom_density(data = p8_post, color = "brown") + #p8
  geom_density(data = p9_post, color = "azure") + # p9
  geom_density(data = p10_post, color = "coral") + # p10
  # Aesthetics
  labs(x = "eta", y = "density") 

x_eta <- seq(0, 6e5, by = 10000)
eta_prior <- dlnorm(x_eta, mean=9.9, sd = 1)
plot(density(p1_post$eta), xlab = "x", main= "eta")
lines(density(p2_post$eta), col ="purple")
lines(density(p3_post$eta), col ="orange")
lines(density(p4_post$eta), col ="blue")
lines(density(p5_post$eta), col ="red")
lines(density(p6_post$eta), col ="green")
lines(density(p7_post$eta), col ="maroon")
lines(density(p8_post$eta), col ="brown")
lines(density(p9_post$eta), col ="azure")
lines(density(p10_post$eta), col ="coral")
lines(x_eta, eta_prior, type = "o", col ="red")


#plot posterior density of omega for all patients 
plot_omega_all = ggplot(data  = p1_post, aes(x = omega)) + #patient 1
  geom_density() + 
  geom_density(data = p2_post, color = "purple") + # p2
  geom_density(data = p3_post, color = "orange") + # p3
  geom_density(data = p4_post, color = "blue") + # p4
  geom_density(data = p5_post, color = "red") + # p5
  geom_density(data = p6_post, color = "green") + # p6
  geom_density(data = p7_post, color = "maroon") + # p7
  geom_density(data = p8_post, color = "brown") + #p8
  geom_density(data = p9_post, color = "azure") + # p9
  geom_density(data = p10_post, color = "coral") + # p10
  # Aesthetics
  labs(x = "omega", y = "density") 
plot_omega_all 

x_omega <- seq(0, 5e5, by = 10000)
omega_prior <- dlnorm(x_omega, mean=9.2, sd = 1)
plot(density(p1_post$omega), xlab = "x", main= "omega",  ylim=c(0,7e-05))
lines(density(p2_post$omega), col ="purple")
lines(density(p3_post$omega), col ="orange")
lines(density(p4_post$omega), col ="blue")
lines(density(p5_post$omega), col ="red")
lines(density(p6_post$omega), col ="green")
lines(density(p7_post$omega), col ="maroon")
lines(density(p8_post$omega), col ="brown")
lines(density(p9_post$omega), col ="azure")
lines(density(p10_post$omega), col ="coral")
lines(x_omega, omega_prior, type = "o", col ="red")

#plot beta
x_beta <- seq(0, 6e-10, by = 1e-11)
beta_prior <- dlnorm(x_beta, mean= -24, sd = 1)
plot(density(p1_post$beta), xlab = "x", main= "beta",  ylim=c(0,1.8e+10))
lines(density(p2_post$beta), col ="purple")
lines(density(p3_post$beta), col ="orange")
lines(density(p4_post$beta), col ="blue")
lines(density(p5_post$beta), col ="red")
lines(density(p6_post$beta), col ="green")
lines(density(p7_post$beta), col ="maroon")
lines(density(p8_post$beta), col ="brown")
lines(density(p9_post$beta), col ="azure")
lines(density(p10_post$beta), col ="coral")
lines(x_beta, beta_prior, type = "", col ="red")
#=======================================================================
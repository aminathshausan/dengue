#' ---
#' title: "Plot raw data and sensitivity of parameters"
#' author: "Aminath Shausan"
#' date: "May , 2022"
#' This program plots DENV1 primary and secondary data and sensitivity of ode parameters
#'  ---------------
#clear history 
rm(list=ls())
#-----------------------------------
# install required packages 
#install.packages("deSolve","dplyr","ggplot2")
#====================
library(deSolve)
library(dplyr)
library(ggplot2)
#==========================================================
# (Fig. 1 ): plot observed DENV1 data 

load('../data/DENV1_data.RData') #loads the saved pre-processed DENV1 data
#Figure 1: plot DENV1 data
ggplot(data = DENV1_data)+  ylim(0, 12)+
  geom_line(mapping = aes(x=FeverDay, y=log10(viremia) , group = StudyNo , color = Serology) , show.legend = FALSE)+ 
  geom_point(mapping = aes(x=FeverDay, y=log10(viremia) , group = StudyNo, color = Serology, shape = viremia_sign)
             , show.legend = FALSE)+ 
labs(x = "Fever day", y = "log10(Viremia)-copies/ml")

ggsave(filename = "../images/raw_data_and_sensitivity_plots/DENV1_raw_data.jpeg", width = 4.5, height = 3.75)
##############################################

# --- plot sensitivity of parameters ---------

#1. Local (Visual) sensitivity Analysis:vary one parameter at a time, while fixing all others at a baseline values and asses the effect of varying this 
#parameter on the output 
# For the baseline trajectory, use the simulation of ode for the values: 
params <- list(A = 1.4e6, beta = 6.65757e-11, delta = 5.927266,
               gamma = 5.078851 , omega = 11332.68, sigma = 0.5, kappa = 3.711207, eta = 42759.61 , mu =0.0001, std =0.01197159 )


inits_base <- c(S=1e8 , I_V=0, I_T=0, I_B=0, V0= 57.122, T=0)
times = -10:5

ode_func <- function(t, y, params) {
  with(as.list(c(params, y)), {
    
    dS = A- beta * y[5] * y[1] - beta * y[6] * y[1] #- d*y[1];                  #Susceptible cells ;
    
    dI_V = beta * y[1] * y[5] - beta * y[6] * y[2]- delta * y[2] -mu*y[2];      #virus only infected cells 
    
    dI_T = beta * y[6] * y[1] - beta * y[5] * y[3] ;                   #TIP only infected  cells 
    
    dI_B = beta * y[6] * y[2] +  beta * y[5] * y[3] - gamma * y[4]+mu*y[2];     #Co-infected cells
    
    dV = omega * y[2] + omega * (1-sigma) * y[4] - kappa * y[5];          #free Virus
    
    dT = eta * sigma* y[4] - kappa * y[6];                              #free TIPs
    
    res1 <- c(dS,dI_V ,dI_T,dI_B,dV,dT)
    list(res1)
  })
}

#1. baseline profile
out_base <- ode(inits_base, times, ode_func, params, method="ode45")
out_base <- data.frame(out_base) #rlnorm(1,meanlog = log(out_base[1:16,6]), sdlog = params$std)
colnames(out_base) <- c("time", "S", "I_V", "I_T", "I_B", "V", "T")
out_base <- as.data.frame(out_base)

ggplot(data = DENV1_data)+  ylim(0, 12)+
  geom_line(mapping = aes(x=FeverDay, y=log10(viremia) , group = StudyNo, color = Serology) , show.legend = FALSE)+ 
  geom_line(data =  out_base, aes(x=time, y=log10(V)), color = "black")+
  geom_hline(yintercept=log10(357), color = "black")+
  labs(x = "Fever Day", y = "log10(Viremia)-copies/ml")

ggsave(filename = "../images/raw_data_and_sensitivity_plots/base_line.jpeg", width = 4.5, height = 3.75)

#=======plot sensitivity of fixed parameters (Figure 2, sub plots)=================

#2. plot sensitivity of A. 
#Note: vary parameter A in prams list above using 1e2 to 1e9. 
#label  each ode output as 'out_A_value' where value = A value
#Save all resulting ode outputs as 'vary_A.RData'. 
#We have saved this simulation prior and saved results

load('../results/vary_A.RData') #loads the saved ODE output by varying A  

ggplot(out_base, aes(x=time, y=log10(V))) + ylim(-0.7, 12)+   #baseline measurement
  geom_line(col="black") +    
  geom_line(data = out_A_1e9 , aes(x=time, y=log10(V)), color = "blue") +
  geom_line(data = out_A_1e8  , aes(x=time, y=log10(V)), color = "green") +
  geom_line(data = out_A_1e7 , aes(x=time, y=log10(V)), color = "pink") +
  geom_line(data = out_A_1e6 , aes(x=time, y=log10(V)), color = "red") +
  geom_line(data = out_A_1e5 , aes(x=time, y=log10(V)), color = "purple")+
  geom_line(data = out_A_1e4 , aes(x=time, y=log10(V)), color = "magenta")+ 
  geom_line(data = out_A_1e3 , aes(x=time, y=log10(V)), color = "brown")+
  geom_line(data = out_A_1e2 ,aes(x=time, y=log10(V)), color = "orange")+
  geom_line(data = out_base, aes(x=time, y=log10(V)), color = "black")+ #baseline
  labs(x = "Fever Day", y = "log10(Viremia)-copies/ml")

ggsave(filename = "../images/raw_data_and_sensitivity_plots/sens_A.jpeg", width = 4.5, height = 3.75)

#3. plot sensitivity of omega

#Note: vary parameter omega in prams list above using values:
#  1e5, 8e5, 6e4, 3e4, 2e4, 1.5e4, 1e4 and 800. 
#label  each ode output as 'out_omega_value' where value = omega value
#Save all resulting ode outputs as 'vary_omega.RData'
# We have saved this simulation prior and saved results

load('../results/vary_omega.RData') #loads the saved ODE output by varying omega  

ggplot(out_base, aes(x=time, y=log10(V))) + ylim(-0.7, 12)+   #baseline measurement
  geom_line(col="black") +    
  geom_line(data = out_omega_1e5 , aes(x=time, y=log10(V)), color = "blue") +
  geom_line(data = out_omega_8e4  , aes(x=time, y=log10(V)), color = "green") +
  geom_line(data = out_omega_6e4 , aes(x=time, y=log10(V)), color = "pink") +
  geom_line(data = out_omega_4e4 , aes(x=time, y=log10(V)), color = "red") +
  geom_line(data = out_omega_2e4 , aes(x=time, y=log10(V)), color = "purple")+
  geom_line(data = out_omega_1.5e4 , aes(x=time, y=log10(V)), color = "magenta")+ 
  geom_line(data = out_omega_1e4 , aes(x=time, y=log10(V)), color = "brown")+
  geom_line(data = out_omega_8e3 ,aes(x=time, y=log10(V)), color = "orange")+
  geom_line(data = out_base, aes(x=time, y=log10(V)), color = "black")+ #baseline
  labs(x = "Fever Day", y = "log10(Viremia)-copies/ml")

ggsave(filename = "../images/raw_data_and_sensitivity_plots/sens_omega.jpeg", width = 4.5, height = 3.75)

#4. plot sensitivity of eta
#Note: vary parameter eta in prams list above using values:
#  1e6, 1e5, 8e4, 5e4, 1e4, 8e3, and 800. 
#label  each ode output as 'out_eta_value' where value = eta value
#Save all resulting ode outputs as 'vary_eta.RData'
# We have saved this simulation prior and saved results

load('../results/vary_eta.RData') #loads the saved ODE output by varying eta 

ggplot(out_base, aes(x=time, y=log10(V))) + ylim(-0.7, 12)+   #baseline measurement
  geom_line(col="black") +    
  geom_line(data = out_eta_1e6 , aes(x=time, y=log10(V)), color = "blue") +
  geom_line(data = out_eta_1e5   , aes(x=time, y=log10(V)), color = "green") +
  geom_line(data = out_eta_8e4 , aes(x=time, y=log10(V)), color = "pink") +
  geom_line(data = out_eta_5e4  , aes(x=time, y=log10(V)), color = "red") +
  geom_line(data = out_eta_1e4  , aes(x=time, y=log10(V)), color = "purple")+
  geom_line(data = out_eta_8e3  , aes(x=time, y=log10(V)), color = "magenta")+ 
  geom_line(data = out_eta_8e2  , aes(x=time, y=log10(V)), color = "brown")+
  geom_line(data = out_base, aes(x=time, y=log10(V)), color = "black")+ #baseline
  labs(x = "Fever Day", y = "log10(Viremia)-copies/ml")

ggsave(filename = "../images/raw_data_and_sensitivity_plots/sens_eta.jpeg", width = 4.5, height = 3.75)

#4. plot sensitivity of sigma
#Note: vary parameter sigma in prams list above using values:
#  1, 0.8, 0.6, 0.4, 0.2, 0.5, and 0. 
#label  each ode output as 'out_sigma_value' where value = sigma value
#Save all resulting ode outputs as 'vary_sigma.RData'
# We have saved this simulation prior and saved results

load('../results/vary_sigma.RData') #loads the saved ODE output by varying sigma 

ggplot(out_base, aes(x=time, y=log10(V))) + ylim(-0.7, 12)+   #baseline measurement
  geom_line(col="black") +    
  geom_line(data = out_sigma_1 , aes(x=time, y=log10(V)), color = "blue") +
  geom_line(data = out_sigma_0.8   , aes(x=time, y=log10(V)), color = "green") +
  geom_line(data = out_sigma_0.6 , aes(x=time, y=log10(V)), color = "pink") +
  geom_line(data = out_sigma_0.4  , aes(x=time, y=log10(V)), color = "red") +
  geom_line(data = out_sigma_0.2  , aes(x=time, y=log10(V)), color = "purple")+
  geom_line(data = out_sigma_0  , aes(x=time, y=log10(V)), color = "magenta")+ 
  geom_line(data = out_base, aes(x=time, y=log10(V)), color = "black")+ #baseline
  labs(x = "Fever Day", y = "log10(Viremia)-copies/ml")

ggsave(filename = "../images/raw_data_and_sensitivity_plots/sens_sigma.jpeg", width = 4.5, height = 3.75)

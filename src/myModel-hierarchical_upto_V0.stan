
//This stan code can either generate fake data from priors (when likelihood is switched off) and compute  prior predictions 
// or fit model to observed data and compute posterior predictions (when likelihood is switched on)
//Author: Aminath Shausan 
//Date created: 22/2/19
///////////////////////////
//Some priors are reparametrised using _raw as there were warnings about divergent transitions 
//In this code, only V0 is patient specific
functions {
     real[] ode(real t,
     real[] y,
     real[] params,
     real[] x_r, // x_r:real data, x_i: integer data
     int[] x_i) {

     real dydt[6];
    
     real A;
     real beta;
     real delta;
     real gamma;
     real omega;
     real sigma;
     real kappa;
     real eta;
     delta = params[1]; // death rate of cells infected with virus(this is estimated first), and std is estimated 2nd
     gamma = params[2]; //death rate of co-infected cells. this is estimated 3rd  
     kappa = params[3]; //virus clearance rate. estimated 4th
     eta   = params[4]; //TIP production rate from co-infected cells. estimated 5th
     omega = params[5]; //virus production rate . estimated 6th
     beta  = params[6]; // rate of infection by virus and TIP particles . estimated 7th
    
     A = x_r[1];
     sigma = x_r[2];
    
     
     dydt[1] =  A - beta * y[5] * y[1] - beta * y[6] * y[1];                 //Susceptible cells ;
     dydt[2] = beta * y[1] * y[5] - beta * y[6] * y[2]- delta * y[2];      //virus only infected cells 
     dydt[3] = beta * y[6] * y[1] - beta * y[5] * y[3] ;                   //TIP only infected  cells 
     dydt[4]=  beta * y[6] * y[2] +  beta * y[5] * y[3] - gamma * y[4];        //Co-infected cells
     dydt[5] = omega * y[2] + omega * (1-sigma) * y[4] - kappa * y[5];          //free Virus
     dydt[6] = eta * sigma* y[4] - kappa * y[6];                          //free TIPs
    
  return dydt;
 }
}

data {
    int<lower=1> J;                  //number of patients
    int<lower=1> n_obs;               //number of observed measurements per patient
    // real y0[6];                      //initial condition of ode 
     real t0[J];                         //initial time 
     real ts[n_obs, J];               //time points of observed measurements
     real A;
     real sigma;
     real<lower =0> y_hat[n_obs,J];        //observed viremia measurements only
     int<lower =1> n_pred;               // number of predicted values from posterior, computed in the generated quantities block 
     real  t_pred[n_pred,J];                //time span at which posterior predictions are made
     int <lower =0, upper =1> is_censored[n_obs, J];  // Indicator matrix: 0 if viremia > LOD, 1 otherwise 
     real  L;                        // censoring point 
    //int<lower=0, upper = 1> run_estimation; // a switch to evaluate the likelihood (use this to simulate data from priors)
}

transformed data {
      real x_r[2];
      int x_i[0];
     
      x_r[1] = A;
      x_r[2]= sigma;
}


parameters{ //the following parameters are to be estimated
        real<lower =0>  delta;         
        real<lower = 0> std;           
        real<lower = 0> gamma;
        real<lower = 0> kappa;        
        real eta_raw;
        real omega_raw;
        real beta_raw; 
        //real<lower =0> mu_delta;         //population parameter of delta
        //real<lower =0> mu_kappa;          //population parameter of kappa
        //real<lower =0> mu_std;          //population parameter of std
        real<lower = 0> V0[J];             //initial viral load is patient specific 
        real<lower =0> mu_V0;            //population mean of V0
        real<lower =0> sigma_V0;         //population standard deviation of V0
}

transformed parameters {
              real y[n_obs,6];       //output from ODE solver
              real y_new[n_obs, J]; //vector to hold machine precision solutions from ODE (only for viremia measurements) 
              real y0[J,6];          //initial condition of ODE (depends on incubation periods)
           //   real<lower = 0> V0[J];           // initial concentration of virus (patient specific)
              real<lower = 0> beta;
              real<lower = 0> eta;
              real<lower = 0> omega;
         
              //reparametrise parameters
              eta   = exp(9.9 + 1*eta_raw);     //this computes log(eta) ~ normal(9.9, 1) hence eta ~ lognormal(9.9,1)
              omega = exp(9.2 + 1*omega_raw);
              beta  = exp(-24 + 1*beta_raw);
              
              
                for (j in 1:J){
                   real params[6];
                  params[1] = delta;
                  params[2] = gamma;
                  params[3] = kappa;
                  params[4] = eta;
                  params[5] = omega;
                  params[6] = beta;
                 y0[j,1] = 1e8;    // initial number of susceptible (uninfected) cells 
                 y0[j,2] = 0;    // initial number of cells infected with only virus
                 y0[j,3] = 0;    // initial number of cells infected with only TIPs
                 y0[j,4] = 0;    // initial number of co-infected cells
                 y0[j,5] = V0[j];    // initial concentration of virus particles
                 y0[j,6] = 1;    // initial concentration of TIPs 
              
               //compute solutions only at the observed time points (ts)
                 y = integrate_ode_rk45(ode, y0[j,], t0[j], ts[,j], params, x_r, x_i, 1e-5, 1e-6, 1e4);
              
            //extract viremia vector and add machine precision
               for (t in 1:n_obs) {
                    y_new[t,j]= y[, 5][t] + 10*1e-6; //add 10*Absolute Tolerance
                 } 
              }
}

model {
      //population parameters (hyper priors)
     //  mu_delta ~ exponential(0.2);       //mean of delta
      // mu_kappa ~ exponential(0.2);      //mean of kappa
     //  mu_std ~ exponential(1);           //mean of std
       mu_V0 ~ normal(0, 100);           //mean of V0
       sigma_V0 ~ exponential(0.1);     // sd of V0
       //sigma_V0 ~ normal(0, 10);        //sd of V0
        
       //priors
        delta ~ exponential(0.2);       
        std ~ exponential(1);          
        gamma ~ exponential(0.2);
        kappa ~ exponential(0.2);      //  original kappa ~ exponential(0.2);
        eta_raw ~ normal(0,1);               //use this to reparametrise eta
        omega_raw ~ normal(0, 1);            //use this to reparametrise omega
        beta_raw ~ normal(0, 1);             //use this to reparametrise beta
  //      V0_raw ~ normal(0,1);
        V0 ~ normal(mu_V0, sigma_V0);         //original V0 ~ normal(0,100)
        
        //for each patient, compute likelihood 
          for (j in 1:J){
             for (n in 1:n_obs){
                if (is_censored[n,j] == 0){
                     y_hat[n,j] ~ lognormal(log(y_new[n,j]), std);                //likelihood for observed data
                } else {
                target += lognormal_lcdf(L | log(y_new[n,j]), std);          //likelihood for censored data
                       }
                }
        }
}

generated quantities {
   
       real y_pred[n_pred, 6];              //predicted ODE solution
       real y_pred_viremia[n_pred, J];      //vector to hold predicted viremia trajectory for each patient
       real y_pred_TIPs[n_pred, J];         //vector to hold predicted TIPs  trajectory for each patient   
       real y_pred_coInf[n_pred, J];        //vector to hold predicted number of co-infected cells for each patient
       real y_pred_Sus[n_pred, J];         //vector to hold predicted number of susceptible cells
       real y_ppc[n_pred,J];                //vector to hold predicted  viremia levels for each patient with noise  
  
       for (j in 1:J){
         real params[6];
         params[1] = delta;
         params[2] = gamma;
         params[3] = kappa;
         params[4] = eta;
         params[5] = omega;
         params[6] = beta;
   
       //compute predictions including incubation period (t_pred), without measurement variablity
         y_pred = integrate_ode_rk45(ode, y0[j,], t0[j], t_pred[,j], params, x_r, x_i, 1e-5, 1e-6, 1e4);
         y_pred_viremia[,j] = y_pred[, 5];
         y_pred_TIPs[,j] = y_pred[,6]; 
         y_pred_Sus[,j] = y_pred[,1];
         y_pred_coInf[, j] = y_pred[,4];
         
         //compute predictions using measurement varaiability (std)
         for (t in 1:n_pred) {
         y_ppc[t,j] = lognormal_rng(log(y_pred[, 5][t]+ 10*1e-6), std);
        
      }
     
       }
} 

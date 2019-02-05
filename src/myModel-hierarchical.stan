
//This stan code can either generate fake data from priors (when likelihood is switched off) and compute  prior predictions 
// or fit model to observed data and compute posterior predictions (when likelihood is switched on)
//Author: Aminath Shausan 
//Date created: 20/12/18
///////////////////////////

//Some priors are reparametrised using _raw as there were warnings about divergent transitions 

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
     real y0[6];                      //initial condition of ode
     real t0;                         //initial time 
     real ts[n_obs, J];               //time points of observed measurements
     real A;
     real sigma;
    real<lower =0> y_hat[n_obs,J];        //observed viremia measurements only
    int<lower =1> n_pred;               // number of predicted values from posterior, computed in the generated quantities block 
    real  t_pred[n_pred,J];                //time span at which posterior predictions are made
   // int<lower = 0> n_cens;                  //number of censered observations
    int <lower =0, upper =1> is_censored[n_obs, J];  // Indicator matrix: 0 if viremia > LOD, 1 otherwise 
     real <upper = min(y_hat[,J])> L;                        // censoring point 
    //int<lower=0, upper = 1> run_estimation; // a switch to evaluate the likelihood (use this to simulate data from priors)
}

transformed data {
      real x_r[2];
      int x_i[0];
     
      x_r[1] = A;
      x_r[2]= sigma;
}


parameters{ //the following parameters are to be estimated
        real<lower =0>  delta[J];      //delta is patient specific
       // real<lower =0>  delta_raw[J];  //hyper prior for each delta
        real<lower = 0> std;           // making std as patient specific made R session to abort
        real<lower = 0> gamma;
        real<lower = 0> kappa[J];     //patient specific
       // real<lower = 0> kappa_raw[J]; //hyper prior for each kappa
        real eta_raw;
        real omega_raw;
        real beta_raw; 
        real<lower =0> mu_delta; //population parameter of delta
        real<lower =0> sigma_delta; //population parameter of delta
        real<lower =0> mu_kappa; //population parameter of kappa
        real<lower =0> sigma_kappa; //population parameter of delta
      //  real<lower =0> lambda_std; //population parameter of std
}

transformed parameters {
              real y[n_obs,6];       //output from ODE solver
              real y_new[n_obs, J]; //vector to hold machine precision solutions from ODE (only for viremia measurements) 
              
              real<lower = 0> beta;
              real<lower = 0> eta;
              real<lower = 0> omega;
              
              //reparametrise parameters
              eta   = exp(9.9 + 1*eta_raw);     //this computes log(eta) ~ normal(9.9, 1) hence eta ~ lognormal(9.9,1)    
              omega = exp(9.2 + 1*omega_raw);
              beta  = exp(-24 + 1*beta_raw);
              
              for (j in 1:J){
              
                  real params[6];
                  params[1] = delta[j];
                  params[2] = gamma;
                  params[3] = kappa[j];
                  params[4] = eta;
                  params[5] = omega;
                  params[6] = beta;
              
              //compute solutions only at the observed time points (ts)
                 y = integrate_ode_rk45(ode, y0, t0, ts[,j], params, x_r, x_i, 1e-5, 1e-6, 1e4);
              
            //extract viremia vector and add machine precision
               for (t in 1:n_obs) {
            //      for (z in 1:6){
                      y_new[t,j]= y[, 5][t] + 10*1e-6; //add 10*Absolute Tolerance
                 } 
           //  }
              }
}

model {
      //population parameters (hyper priors)
      //with normal hyper priors for delta and kappa, it took about 43 mins with 2 patients data and many div. transitions
       mu_delta ~ normal(0, 10);        //mean of delta
       sigma_delta ~ exponential(0.2); //standard deviation of delta
       mu_kappa ~ normal(0, 10);        //mean of delta
       sigma_kappa ~ exponential(0.2); //standard deviation of delta
     // lambda_std ~ exponential(1);   //this also made R session to abort 
      //lambda_std ~ normal(0, 1); //95%CI for lambda_kappa in [0, 0.6] (with this, R session aborted)
      
     // lambda_delta ~ exponential(0.2); //implies mean of lambda_delta = 1/0.2 =5
     // lambda_kappa ~ exponential(0.2); 
      
      //hyper priors 
 //     delta_raw ~ exponential(1/lambda_delta); 
  //    kappa_raw ~ exponential(1/lambda_kappa);
      
      //priors
      //  delta ~ exponential(1/lambda_delta); //  original delta ~ exponential(0.2);
        delta ~ normal(mu_delta, sigma_delta);
          
        std ~ exponential(1);
       // std ~ exponential(lambda_std);     //  original std ~ exponential(1);
        gamma ~ exponential(0.2);
    //    kappa ~ exponential(1/lambda_kappa); //  original kappa ~ exponential(0.2);
        kappa ~ normal(mu_kappa, sigma_kappa);
        eta_raw ~ normal(0,1);
        omega_raw ~ normal(0, 1);
        beta_raw ~ normal(0, 1); //original: beta_raw ~normal(0,1)
        
        //evaluate likelihood conditionally. run_estimation =0 switch off likelihood. run_estimation =1 switch on likelihood
       // if(run_estimation==1){
          //for each patient, compute likelihood 
          for (j in 1:J){
             for (n in 1:n_obs){
                if (is_censored[n,j] == 0){
                     y_hat[n,j] ~ lognormal(log(y_new[n,j]), std);                //likelihood for observed data
                } else {
                target += lognormal_lcdf(L | log(y_new[n,j]), std);          //likelihood of censored data
                       }
                }
        }
}

generated quantities {
   
       real y_pred[n_pred, 6];
        real y_ppc[n_pred,J];
       
       
       for (j in 1:J){
         real params[6];
         params[1] = delta[j];
         params[2] = gamma;
         params[3] = kappa[j];
         params[4] = eta;
         params[5] = omega;
         params[6] = beta;
      
       //compute predictions including incubation period (t_pred), without measurement variablity
         y_pred = integrate_ode_rk45(ode, y0, t0, t_pred[,j], params, x_r, x_i, 1e-5, 1e-6, 1e4);
         
         //compute predictions using measurement varaiability (std)
         for (t in 1:n_pred) {
         y_ppc[t,j] = lognormal_rng(log(y_pred[, 5][t]+ 10*1e-6), std);
        
      }
     
       }
} 


//This stan code can either generate fake data from priors (when likelihood is switched off) and compute  prior predictions 
// or fit model to observed data and compute posterior predictions (when likelihood is switched on)
//Author: Aminath Shausan 
//Date created: 27/03/19
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
     real mu;
     delta = params[1]; // death rate of cells infected with virus(this is estimated first), and std is estimated 2nd
     gamma = params[2]; //death rate of co-infected cells. this is estimated 3rd  
     kappa = params[3]; //DENV clearance rate. estimated 4th
     beta  = params[4]; // rate of infection by virus and TIP particles . estimated 7th
     mu    = params[5]; // rate virus mutate to form TIps within cells infected with only TIPs
     A = x_r[1];
     omega = x_r[2];
     sigma = x_r[3];
     eta = x_r[4];
     
     dydt[1] =  A - beta * y[5] * y[1] - beta * y[6] * y[1];                 //Susceptible cells ;
     dydt[2] = beta * y[1] * y[5] - beta * y[6] * y[2]- delta * y[2]- mu*y[2];      //virus only infected cells 
     dydt[3] = beta * y[6] * y[1] - beta * y[5] * y[3] ;                   //TIP only infected  cells 
     dydt[4]=  beta * y[6] * y[2] +  beta * y[5] * y[3] - gamma * y[4]+ mu*y[2];        //Co-infected cells
     dydt[5] = omega * y[2] + omega * (1-sigma) * y[4] - kappa * y[5];          //free Virus
     dydt[6] = eta * sigma* y[4] - kappa * y[6];                          //free TIPs
    
  return dydt;
 }
}

 data {
       int<lower=1> J;                  //number of patients
       real<lower =0> A;
       real<lower =0> omega;
       real<lower =0>sigma;
       real<lower =0> eta;
       int<lower=1> n_obs;          //number of observed measurements 
       real t0[J];                     //initial time of ODE
       real ts[n_obs,J];
       real<lower =0> y_hat[n_obs,J];
       int<lower =1> n_pred;            // number of predicted values from posterior, computed in the generated quantities block 
       real  t_pred[n_pred,J];
       int <lower =0, upper =1> is_censored[n_obs,J];  // Indicator matrix: 0 if viremia > LOD, 1 otherwise 
       real  L;                                      // censoring point 
       int<lower=0, upper = 1> run_estimation; // a switch to evaluate the likelihood
       int<lower =1> K;                        //number of groups
       int<lower=1, upper=2> serology[J]; // 1=primary, 2=secondary
}

 transformed data {
       real x_r[4];
       int x_i[0];
       x_r[1] = A;
       x_r[2]= omega;
       x_r[3]= sigma;
       x_r[4]= eta;
 }

 parameters{ //the following parameters are to be estimated
     
     real <lower =0, upper =20> delta[K];
      real <lower =0, upper =20> gamma[K];
     real <lower =0, upper =20> kappa[K];
      real <lower =1e-11, upper =5e-10> beta[K];
      real <lower =1e-5, upper =1e-2> mu[K];
      real<lower =0> V0[J];
      real<lower = 0> std;
}

transformed parameters {
                real y[n_obs,6];  //output from ODE solver
               real y_new[n_obs, J]; //vector to hold machine precision solutions from ODE (for only viremia measurements)
               real y0[J,6];        //initial condition of ODE
               real params[K,5];
               
               
               for(k in 1:K){ //parameters for each group
                        params[k,1] = delta[k];
                        params[k,2] = gamma[k];
                        params[k,3] = kappa[k];
                        params[k,4] = beta[k];
                        params[k,5] = mu[k];
               }
            
            for(j in 1:J){
                        y0[j,1] = 1e8;    // initial number of susceptible (uninfected) cells
                        y0[j,2] = 0;    // initial number of cells infected with only virus
                        y0[j,3] = 0;    // initial number of cells infected with only TIPs
                        y0[j,4] = 0;    // initial number of co-infected cells
                        y0[j,5] = V0[j];    // initial concentration of virus particles
                        y0[j,6] = 0;    // initial concentration of TIPs 
            
               //compute solutions only at the observed time points (ts)
                 y = integrate_ode_rk45(ode, y0[j,], t0[j], ts[,j], params[serology[j],], x_r, x_i, 1e-5, 1e-6, 1e4); //use this for primary/secondary grouping
               
                //add machine precision to vector y (only for viremia trajectory)
                  for (t in 1:n_obs) {
                     y_new[t,j]= y[, 5][t] + 10*1e-6; //add 10*Absolute Tolerance
                    } 
              }
 }

model {
          std ~ normal(1, 2);
          V0 ~ normal(0, 100);
       
        //evaluate likelihood conditionally. run_estimation =0 switch off likelihood. run_estimation =1 switch on likelihood
         if(run_estimation==1){
           for (j in 1:J){
              for (n in 1:n_obs){
                 if (is_censored[n,j] == 0){
                    y_hat[n,j] ~ lognormal(log(y_new[n,j]), std);            //likelihood for uncensored observed data
               } else {
              
                target += lognormal_lcdf(L |log(y_new[n,j]), std);          //likelihood of censoreddata
               }
            }
         }
       }
     
}
   generated quantities {
               real y_pred[n_pred, 6];
           real y_pred_S[n_pred, J];
           real y_pred_IV[n_pred, J];
           real y_pred_IT[n_pred, J];
           real y_pred_IB[n_pred, J];
          real y_pred_viremia[n_pred, J];
          real y_pred_T[n_pred, J];
          real y_ppc[n_pred, J];
     
      //compute predictions including incubation period (t_pred), without measurement variablity
          for (j in 1:J){
           y_pred = integrate_ode_rk45(ode, y0[j,], t0[j], t_pred[,j], params[serology[j],], x_r, x_i, 1e-5, 1e-6, 1e4);
         
            y_pred_S[,j] = y_pred[, 1];         //target cells
            y_pred_IV[,j] = y_pred[, 2];        // cells infected with only viruses
            y_pred_IT[,j] = y_pred[, 3];        //cells infected with only TIPs
            y_pred_IB[,j] = y_pred[, 4];       //co-infected cells
            y_pred_viremia[,j] = y_pred[, 5];       //free virus
            y_pred_T[,j] = y_pred[, 6];       //free TIPs
         //compute predictions using measurement varaiability (std) only for viremia measurements
           for (t in 1:n_pred) {
              y_ppc[t,j] = lognormal_rng(log(y_pred[, 5][t]+ 10*1e-6), std);
         }
       }
}

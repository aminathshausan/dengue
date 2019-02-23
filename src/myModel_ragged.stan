//This stan code is used to check if ragged and censored data structure works appropriately
//Author: Aminath Shausan 
//Date created: 12/02/19
///////////////////////////

//define ragged data structure
data {
    int<lower=1> J;                  //number of patients
    int<lower=1> N;                 //total number of observations in the data set
  //  real<lower =0> y_hat[N];        //observed viremia measurements.all observed measurements pushed together in a flat array 
     real ts[N];            //time points of observed measurements
     vector[N] y_hat;
   //  vector[N] ts;
    int<lower =1, upper =J>idx[N];  //index to identify each patient
   // int gp_sizes[J];          //number of observations for each patient
    //later try changing vector to arrays
}


parameters{ 
  
}

model {
      
}

generated quantities {
         real viremia[N];
        real obs_times[N];
       // vector[N] viremia;
       // vector[N] obs_times;
         
      //check if the data extracts correct measurements
      for (n in 1:N){
        viremia[n] = y_hat[n];
        obs_times[n] = ts[n];
      }
       // int count;
      //  count = 1;
      //  for (j in 1:J){//this is nor working yet 
      //    viremia[count:gp_sizes[j]] = segment(y_hat, count, gp_sizes[j]);
      //    obs_times[count:gp_sizes[j]] = segment(ts, count, gp_sizes[j]);
        //  viremia[count:gp_sizes[j]] = segment(y_hat, count, gp_sizes[j]);
         // obs_times[count:gp_sizes[j]] = segment(ts, count, gp_sizes[j]);
      //    count = count + gp_sizes[j];
     //   }
} 

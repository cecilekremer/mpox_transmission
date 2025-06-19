data{
  int<lower=1> N; // number of observations: individuals x visits x possible contacts
  int<lower=1> M; // number of individuals x visits
  
  // Data
  array[N] real<lower=0> t_s; // symptom onset
  array[N] real<lower=0> t_e_r; // upper bound exposure
  array[N] real<lower=0> t_e_l; // lower bound exposure
  array[N] int<lower=0,upper=1> censored; // single (0) or multiple (1) exposure
  
  array[M] int<lower=1> n_contacts; // number of contacts for each case
  
  // prior on duration of lesions; assuming infectious period corresponds to lesion period
  real<lower=0> log_mu_lesions; 
  real<lower=0> sigma_lesions;
  
  // Auxiliary data
  array[N] int<lower=0> n_steps_prior; // number of days over which to integrate
}
transformed data{
  vector[N] log_prior_contacts; // prior on who is the infecting contact
  array[N] vector[max(n_steps_prior)] log_prior_exposures; // prior on when exposure occurred
  
  for(m in 1:M){
    // for(i in 1:n_contacts[m]){
      log_prior_contacts[m] = log(1./n_contacts[m]); // uniform prior over contacts
      // }
  }
  
  for(i in 1:N){
    int n = n_steps_prior[i];
    array[n] real ts;
    ts = linspaced_array(n, t_e_r[i] - n, t_e_r[i]);
    
    if(n > 0){
      for(j in 1:n){
        log_prior_exposures[i][j] = lognormal_lccdf(t_e_r[i] - ts[j] + 1 | log_mu_lesions, sigma_lesions);
        if (is_nan(log_prior_exposures[i][j] )) {
          print("ll prior, i", i, " j", j, " ll", log_prior_exposures[i][j]);
        }
      }
      // normalize
      log_prior_exposures[i][1:n] = log(softmax(log_prior_exposures[i][1:n]));
      // if(i == 1){
      //   print(log_prior_exposures[i]);
      // }
    }
  }
}
parameters{
  real<lower=0> par[2]; // distributional parametsrs
  // vector<lower=0, upper=1>[N] uE;	// Uniform value for sampling between start and end exposure
}
transformed parameters{
  // vector[M] tE; 	// infection moment
  vector[M] log_lik; // log-likelihood vector
  
  for(m in 1:M){
    int n = n_contacts[m];
    int cnt = 1;
    vector[n] ll; // placeholder for logliks for each contact
    
    // integrate over all contacts
    for(i in 1:n){
      
      real ll_c;
      
      if(censored[i] == 1){ // multiple occasions
      
        int k = n_steps_prior[i];
        vector[k] ll_exp; // loglik of exposures prior to most recent contact
        vector[k] ll_exp2 = rep_vector(-log(k), k); // ll denominator to account for truncation
        array[k] real ts; // times of possible exposure to integrate over
        
        ts = linspaced_array(k, t_e_r[i] - k, t_e_r[i]);
        
        for(j in 1:k){
          
          real tE = t_s[i] - ts[j];
          // real bound = t_e_r[i] - ts[j] + 21;
          real this_ll;
          
          this_ll = %s_lpdf(tE | par[1], par[2]);
          
          ll_exp[j] = log_prior_exposures[i][j] + this_ll;
          
        }
        
        ll_c = log_sum_exp(ll_exp);
        
      }else{ // single occasion
        
        real tE = t_s[i] - t_e_l[i]; // exposure time = infection time
        ll_c = %s_lpdf(tE | par[1], par[2]);
        
      }
      
      // this contact's loglik
      ll[cnt] = log_prior_contacts[i] + ll_c;
      cnt += 1;
      
    }
    
    // sum over contacts
    log_lik[m] = log_sum_exp(ll);
  }
  
}
model{
  // Contribution to likelihood of incubation period
  // this is the log-likelihood
  target += sum(log_lik);
  // target += %s_lpdf(tSymptomOnset -  tE  | par[1], par[2]);
  
}
// generated quantities {
//   // likelihood for calculation of looIC
//   vector[M] log_lik;
//   for (i in 1:M) {
//     log_lik[i] = %s_lpdf(tSymptomOnset[i] -  tE[i]  | par[1], par[2]);
//   }
// }
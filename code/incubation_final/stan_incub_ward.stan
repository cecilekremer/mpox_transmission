// distribution: weibull
// truncation:   right
// doubly interval censored

data {
  int<lower = 0> N;                     // number of records
  vector<lower = 0>[N] a_minus;         // lower limit of event A = start exposure
  vector<lower = 0>[N] a_plus;          // upper limit of event A = end exposure
  vector<lower = 0>[N] b_minus;         // lower limit of event B = symptom onset - 0.5
  vector<lower = 0>[N] b_plus;          // upper limit of event B = symptom onset + 0.5
  // int<lower = 0, upper = 1> incubation; // inclusion of time from O to A (incubation period)
  real<lower = 0> upper_bound;          // the latest time of observation
  int<lower=1, upper =2> incl_cov;   // indicator whether covariate is included
  int<lower = 0, upper = 1> sexual[N]; // indicator for sexual transmission
}

transformed data {
  int X_i[0];                           // empty array
}

parameters {
  real<lower = 0> mean_[incl_cov];                    // the distribution mean
  real<lower = 0> alpha[incl_cov];                    // the alpha parameter
  vector<lower = 0, upper = 1>[N] a_window; // where time a lies in the event A window
  vector<lower = 0, upper = 1>[N] b_window; // where time b lies in the event B window
  // vector<lower = 0, upper = 1>[N] a2_window; // where time a2 lies in the event A window
  // vector<lower = 0>[N] t0;                  // time from O to A
}

transformed parameters {
  real<lower = 0> beta[incl_cov]; 
  
  for(i in 1:incl_cov){
    beta[i] = mean_[i]/tgamma(1.0 + 1.0/alpha[i]);
  }
  
  vector<lower = min(a_minus), upper = max(a_plus)>[N] a;
  // vector<lower = min(a_minus), upper = max(a_plus)>[N] a2;
  vector<lower = min(b_minus), upper = max(b_plus)>[N] b;
  vector[N] ub;
  vector<lower = 0>[N] tstar;
  
  b = b_minus + (b_plus - b_minus) .* b_window; // sample exact symptom onset time
  
  
  for (n in 1:N){
    ub[n] = min([a_plus[n], b[n]]'); // take minimum of upper exposure and sampled symptom onset
  }
  a = a_minus + (ub - a_minus) .* a_window; // sample time of infection from exposure window
  // a2 = a_minus + (ub - a_minus) .* a2_window;
  
  for (n in 1:N){
    ub[n] = min([upper_bound,a[n]+11]'); // ?? take minimum of last observation date and sampled time of infection+11
  }
  tstar = upper_bound - a;    
  
  
}

model {
  for(i in 1:incl_cov){
    mean_[i] ~ normal(5.0, 10.0);
    alpha[i] ~ exponential(0.0001);    
  }
  
  // if (incubation)
  // t0 ~ lognormal(1.63, 0.5);
  // else
  // t0 ~ normal(0, 1e-10);
  
  if(incl_cov == 2){
    for(n in 1:N){
      if(sexual[n] == 0){
        target += weibull_lpdf((b[n] - a[n]) | alpha[1], beta[1]) - weibull_lcdf(upper_bound - a[n] | alpha[1], beta[1]);
      }else if(sexual[n] == 1){
        target += weibull_lpdf((b[n] - a[n]) | alpha[2], beta[2]) - weibull_lcdf(upper_bound - a[n] | alpha[2], beta[2]);
      }
    }
  }else{
    target += weibull_lpdf((b - a) | alpha[1], beta[1]) - weibull_lcdf(upper_bound - a | alpha[1], beta[1]);
  }
  
}

generated quantities {
  real sd_[incl_cov];
  real limit_val_[incl_cov];
  real median_[incl_cov];
  
  for(i in 1:incl_cov){
    sd_[i] = beta[i]*sqrt(tgamma(1.0+2.0/alpha[i])-(tgamma(1.0+1.0/alpha[i]))^2);
    limit_val_[i] = beta[i]*(-log(1-0.95))^(1/alpha[i]);
    median_[i] = beta[i]*(-log(1-0.5))^(1/alpha[i]);    
  }
  
  vector[N] log_likelihood;
  for (n in 1:N){
    if(incl_cov == 2){
      if(sexual[n] == 0){
        log_likelihood[n] = weibull_lpdf(b[n] - a[n] | alpha[1], beta[1]) - log(weibull_cdf(upper_bound - a[n] | alpha[1], beta[1]) -
        weibull_cdf(1 | alpha[1], beta[1]) // conditioning on incubation period being > 1 day: https://www.medrxiv.org/content/10.1101/2024.07.22.24310801v1.full
        );
      }else if(sexual[n] == 1){
        log_likelihood[n] = weibull_lpdf(b[n] - a[n] | alpha[2], beta[2]) - log(weibull_cdf(upper_bound - a[n] | alpha[2], beta[2]) -
        weibull_cdf(1 | alpha[2], beta[2]) // conditioning on incubation period being > 1 day: https://www.medrxiv.org/content/10.1101/2024.07.22.24310801v1.full
        );
      }
    }else{
      log_likelihood[n] = weibull_lpdf(b[n] - a[n] | alpha[1], beta[1]) - log(weibull_cdf(upper_bound - a[n] | alpha[1], beta[1]) -
        weibull_cdf(1 | alpha[1], beta[1]) // conditioning on incubation period being > 1 day: https://www.medrxiv.org/content/10.1101/2024.07.22.24310801v1.full
        );
    }
  }
  
}

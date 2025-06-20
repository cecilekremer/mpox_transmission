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
  // int<lower=1, upper =2> incl_cov;   // indicator whether covariate is included
  int<lower = 0, upper = 1> sexual[N]; // indicator for sexual transmission
  // real<lower = 0> age[N];
  int<lower = 0, upper = 1> age[N]; // age group (1 = adult)
  int<lower = 0, upper = 1> kamituga[N];
}

transformed data {
  int X_i[0];                           // empty array
}

parameters {
  real theta[3]; // regression parameters
  real<lower = 0> alpha;                    // the alpha parameter (Weibull shape)
  vector<lower = 0, upper = 1>[N] a_window; // where time a lies in the event A window
  vector<lower = 0, upper = 1>[N] b_window; // where time b lies in the event B window
  // vector<lower = 0, upper = 1>[N] a2_window; // where time a2 lies in the event A window
  // vector<lower = 0>[N] t0;                  // time from O to A
}

transformed parameters {
  real<lower = 0> beta[N]; // Weibull scale
  real<lower = 0> mean_[N];
  
  for(i in 1:N){
    mean_[i] = theta[1] + theta[2]*sexual[i] + theta[3]*age[i]; //+ 
                // theta[4]*kamituga[i] +
                //theta[4]*age[i]*sexual[i];
    beta[i] = mean_[i]/tgamma(1.0 + 1.0/alpha);
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
  
  alpha ~ exponential(0.0001);
  for(k in 1:3){
    theta[k] ~ normal(0, 10);
  }
  
  // if (incubation)
  // t0 ~ lognormal(1.63, 0.5);
  // else
  // t0 ~ normal(0, 1e-10);
  
  for(i in 1:N){
    target += weibull_lpdf((b[i] - a[i]) | alpha, beta[i]) - weibull_lcdf(upper_bound - a[i] | alpha, beta[i]);
  }
  
}
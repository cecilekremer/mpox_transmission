data{
  int N; // number of observations
  int K; // number of factor levels (sexual vs non-sexual)
  // real X[N]; // incubation period values
  int<lower=1, upper=K> sex[N]; // factor (sexual transmission)
  real tStartExposure[N];
  real tEndExposure[N];
  real tSymptomOnset[N];
}
parameters{
  real<lower=0> par[2*K];
  vector<lower=0, upper=1>[N] uE;	// Uniform value for sampling between start and end exposure
}
transformed parameters{
  vector[N] tE; 	// infection moment
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
}
model{
  // Contribution to likelihood of incubation period
  target += %s_lpdf(tSymptomOnset -  tE  | par[1], par[2]);
}
generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = %s_lpdf(tSymptomOnset[i] -  tE[i]  | par[1], par[2]);
  }
}
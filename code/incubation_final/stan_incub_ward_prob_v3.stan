transformed parameters {
  real<lower = 0> beta[incl_cov];
  for(i in 1:incl_cov) {
    beta[i] = mean_[i] / tgamma(1.0 + 1.0 / alpha[i]);
  }
}

model {
  int num_steps = 10; // Number of integration points for numerical approximation

  for (i in 1:incl_cov) {
    mean_[i] ~ normal(5.0, 10.0);
    alpha[i] ~ exponential(0.0001);
  }

  real delta; // Step size for integration
  vector[N] exposure_likelihood;

  for (n in 1:N) {
    real likelihood_sum = 0;
    
    delta = (a_plus[n] - a_minus[n]) / num_steps; // Step size
    
    for (i in 0:(num_steps - 1)) {
      real candidate_a = a_minus[n] + i * delta; // Discretized infection time
      real weight = lognormal_lpdf(candidate_a | log(14.1), 0.5); // Prior probability of infection time

      if (incl_cov == 2) {
        if (sexual[n] == 0) {
          likelihood_sum += exp(weibull_lpdf(b[n] - candidate_a | alpha[1], beta[1]) + weight) * delta;
        } else if (sexual[n] == 1) {
          likelihood_sum += exp(weibull_lpdf(b[n] - candidate_a | alpha[2], beta[2]) + weight) * delta;
        }
      } else {
        likelihood_sum += exp(weibull_lpdf(b[n] - candidate_a | alpha[1], beta[1]) + weight) * delta;
      }
    }

    exposure_likelihood[n] = log(likelihood_sum);
  }

  target += sum(exposure_likelihood);
}

data {
  int<lower=0> N_unc;
  int<lower=0> N_cen;
  int L;  // minimum duration
  real<lower=0> dur_unc[N_unc];
  real<lower=0> dur_cen[N_cen];
}
parameters {
  real<lower=0> sigma;
  real<lower=0> alpha;
}
model {
  sigma ~ gamma(1, 0.0001);
  alpha ~ gamma(1, 0.0001);

  for(i in 1:N_unc) {
    dur_unc[i] ~ weibull(alpha, sigma) T[L,];
  }
  increment_log_prob(weibull_ccdf_log(dur_cen, alpha, sigma));
}


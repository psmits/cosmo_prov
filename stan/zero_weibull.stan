data {
  int<lower=0> N_unc;
  int<lower=0> N_cen;
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

  dur_unc ~ weibull(alpha, sigma);
  increment_log_prob(weibull_ccdf_log(dur_cen, alpha, sigma));
}


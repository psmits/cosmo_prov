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
  sigma ~ cauchy(0, 2.5);
  alpha ~ cauchy(0, 2.5);

  for(i in 1:N_unc) {
    if(dur_unc[i] == L) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha, sigma));
    } else {
      dur_unc[i] ~ weibull(alpha, sigma);
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha, sigma));
  }
}


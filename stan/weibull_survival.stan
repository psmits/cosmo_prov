data {
  int<lower=0> N;
  int<lower=0> N_unc;
  int<lower=0> N_cen;
  real<lower=0> dur_unc[N_unc];
  real<lower=0> dur_cen[N_cen];
  vector[N_unc] size_unc;
  vector[N_cen] size_cen;
  vector[N_unc] occ_unc;
  vector[N_cen] occ_cen;
  vector[N_unc] herb_unc;
  vector[N_cen] herb_cen;
  vector[N_unc] insect_unc;
  vector[N_cen] insect_cen;
  vector[N_unc] omni_unc;
  vector[N_cen] omni_cen;
  vector[N_unc] ground_unc;
  vector[N_cen] ground_cen;
  vector[N_unc] scans_unc;
  vector[N_cen] scans_cen;
}
parameters {
  vector[8] beta;
  real<lower=0> alpha;
}
model {
  alpha ~ gamma(1, 0.0001);

  dur_unc ~ weibull(alpha, exp(-(beta[1] + beta[2] * size_unc + beta[3] * occ_unc +
  beta[4] * herb_unc + beta[5] * insect_unc + beta[6] * omni_unc + beta[7] * ground_unc + beta[8] * scans_unc) / alpha));
  increment_log_prob(weibull_ccdf_log(dur_cen, alpha, 
        exp(-(beta[1] + beta[2] * size_cen + beta[3] * occ_cen + beta[4] * herb_cen + beta[5] * insect_cen + beta[6] * omni_cen + beta[7] * ground_cen + beta[8] * scans_cen) / alpha)));
}

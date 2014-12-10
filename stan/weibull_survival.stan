data {
  int<lower=0> N;
  int<lower=0> N_unc;
  int<lower=0> N_cen;
  int D;  // number of diet categories
  int M;  // number of diet categories
  int L;  // minimum duration
  int C;  // number of cohorts
  real<lower=0> dur_unc[N_unc];
  real<lower=0> dur_cen[N_cen];
  vector[N_unc] size_unc;
  vector[N_cen] size_cen;
  vector[N_unc] occ_unc;
  vector[N_cen] occ_cen;
  matrix[N_unc, D] diet_unc;
  matrix[N_cen, D] diet_cen;
  matrix[N_unc, M] move_unc;
  matrix[N_cen, M] move_cen;
  int coh_unc[N_unc];
  int coh_cen[N_cen];
  int samp_unc[N_unc];
  int samp_cen[N_cen];
  matrix[N, N] vcv;
}
transformed data {
 matrix[N, N] vcv_inv;
 vector[N] mumu;

 vcv_inv <- inverse(vcv);

 for(i in 1:N) mumu[i] <- 0;
}
parameters {
  real beta_inter;
  real beta_occ;
  real beta_size;
  real<lower=0> alpha;
  vector[M] beta_move;
  vector[D] beta_diet;
  real<lower=0> fv;
  real rando[C];
  real<lower=0> sigma_phy;
  vector[N] phy;
}
model {
  beta_inter ~ normal(0, 10);
  beta_occ ~ normal(0, 10);
  beta_size ~ normal(0, 10);
  for(i in 1:M) {
    beta_move[i] ~ normal(0, 10);
  }
  for(i in 1:D) {
    beta_diet[i] ~ normal(0, 10);
  }

  fv ~ cauchy(0, 2.5);
  for(i in 1:C) {
    rando[i] ~ normal(0, fv);
  }

  sigma_phy ~ cauchy(0, 2.5);
  // non-constant part of log(det(sigma_phy * vcv) ^ -0.5
  lp__ <- lp__ - 0.5 * N * log(sigma_phy);  
  // log of kernal of mulinorm
  lp__ <- lp__ - (phy' * vcv_inv * phy) / (2 * sigma_phy);  

  alpha ~ cauchy(0, 2.5);

  for(i in 1:N_unc) {
    if(dur_unc[i] == L) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha,
            exp(-(beta_inter +
                beta_occ * occ_unc[i] + 
                beta_size * size_unc[i] + 
                diet_unc[i] * beta_diet +
                move_unc[i] * beta_move +
                rando[coh_unc[i]] + phy[samp_unc[i]]) / alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha,
            exp(-(beta_inter +
                beta_occ * occ_unc[i] + 
                beta_size * size_unc[i] + 
                diet_unc[i] * beta_diet +
                move_unc[i] * beta_move +
                rando[coh_unc[i]] + phy[samp_unc[i]]) / alpha)));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(beta_inter +
              beta_occ * occ_cen[i] + 
              beta_size * size_cen[i] + 
              diet_cen[i] * beta_diet +
              move_cen[i] * beta_move +
              rando[coh_cen[i]] + phy[samp_cen[i]]) / alpha)));
  }
}

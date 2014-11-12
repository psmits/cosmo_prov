data {
  int<lower=0> N;
  int<lower=0> N_unc;
  int<lower=0> N_cen;
  int D;  // number of diet categories
  int M;  // number of diet categories
  int R;  // number of interactions
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
  matrix[N_unc, R] rac_unc;
  matrix[N_cen, R] rac_cen;
  int coh_unc[N_unc];
  int coh_cen[N_cen];
}
parameters {
  real beta_inter;
  real beta_occ;
  real beta_size;
  real<lower=0> alpha;
  vector[M] beta_move;
  vector[D] beta_diet;
  vector[R] beta_rac;
}
model {
  beta_inter ~ normal(0, 100);
  beta_occ ~ normal(0, 10);
  beta_size ~ normal(0, 10);
  for(i in 1:M) {
    beta_move[i] ~ normal(0, 10);
  }
  for(i in 1:D) {
    beta_diet[i] ~ normal(0, 10);
  }
  for(i in 1:R) {
    beta_rac[i] ~ normal(0, 10);
  }

  alpha ~ cauchy(0, 2.5);

  for(i in 1:N_unc) {
    if(dur_unc[i] == L) {
      increment_log_prob(weibull_cdf_log(dur_unc[i], alpha,
            exp(-(beta_inter +
                beta_occ * occ_unc[i] + 
                beta_size * size_unc[i] + 
                diet_unc[i] * beta_diet +
                move_unc[i] * beta_move + 
                rac_unc[i] * beta_rac) / alpha)));
    } else {
      increment_log_prob(weibull_log(dur_unc[i], alpha,
            exp(-(beta_inter +
                beta_occ * occ_unc[i] + 
                beta_size * size_unc[i] + 
                diet_unc[i] * beta_diet +
                move_unc[i] * beta_move + 
                rac_unc[i] * beta_rac) / alpha)));
    }
  }
  for(i in 1:N_cen) {
    increment_log_prob(weibull_ccdf_log(dur_cen[i], alpha,
          exp(-(beta_inter +
              beta_occ * occ_cen[i] + 
              beta_size * size_cen[i] + 
              diet_cen[i] * beta_diet +
              move_cen[i] * beta_move + 
              rac_cen[i] * beta_rac) / alpha)));
  }
}

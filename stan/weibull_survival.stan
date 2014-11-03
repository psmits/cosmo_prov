data {
  int<lower=0> N;
  int<lower=0> N_unc;
  int<lower=0> N_cen;
  int D;  // number of diet categories
  int M;  // number of diet categories
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
}
parameters {
  real beta_inter;
  vector[2] beta;
  vector[D] beta_diet;
  vector[M] beta_move;
  real<lower=0> alpha;
}
model {
  beta_inter ~ normal(0, 100);
  beta ~ normal(0, 100);
  beta_diet ~ normal(0, 100);
  beta_move ~ normal(0, 100);

  alpha ~ gamma(1, 0.0001);

  dur_unc ~ weibull(alpha, exp(-(beta_inter +
          beta[1] * occ_unc + 
          beta[2] * size_unc + 
          diet_unc * beta_diet + 
          move_unc * beta_move) / alpha));
  increment_log_prob(weibull_ccdf_log(dur_cen, alpha, 
        exp(-(beta_inter +
            beta[1] * occ_cen + 
            beta[2] * size_cen + 
            diet_cen * beta_diet + 
            move_cen * beta_move) / alpha)));
}

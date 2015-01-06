data {
  int<lower=0> N;  // sample size
  int D;  // number of diet categories
  int M;  // number of move categories
  int<lower=0> degree[N];  // # i co-occurs with
  vector[N] mass;
  matrix[N, D] diet;
  matrix[N, M] move;
  
  matrix[N, N] vcv;  // phylogenetic covariance
}
transformed data {
  matrix[N, N] vcv_inv;

  vcv_inv <- inverse(vcv);
}
parameters {
  real beta_inter;
  real beta_mass;
  vector[M] beta_move;
  vector[D] beta_diet;

  real<lower=0> sigma_phy;
  vector[N] phy;

  real<lower = 0> omega;
}
transformed parameters {
  // make a variance
  real<lower=0> sq_sigma;
  real<lower=0> phi;

  sq_sigma <- sigma_phy^2;  

  phi <- 1 / omega;
}
model {
  vector[N] mu;
  
  beta_inter ~ normal(0, 10);
  beta_mass ~ normal(0, 10);
  for(i in 1:M) {
    beta_move[i] ~ normal(0, 10);
  }
  for(i in 1:D) {
    beta_diet[i] ~ normal(0, 10);
  }

  omega ~ cauchy(0, 2.5);

  // phylogenetic effect
  sigma_phy ~ cauchy(0, 2.5);
  // non-constant part of log(det(sigma_phy * vcv) ^ -0.5
  increment_log_prob(-0.5 * N * log(sq_sigma));
  // log of kernal of mulinorm
  increment_log_prob(-(transpose(phy) * vcv_inv * phy) / (2 * sq_sigma));
  
  mu <- (beta_inter + beta_mass * mass + 
        diet * beta_diet + move * beta_move +
        phy);

  degree ~ neg_binomial_2_log(mu, phi);
}




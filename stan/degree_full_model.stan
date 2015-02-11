data {
  int<lower=0> N;  // sample size
  int D;  // number of diet categories
  int M;  // number of move categories
  vector[N] off;
  int<lower=0> degree[N];  // # i co-occurs with
  vector[N] mass;
  matrix[N, D] diet;
  matrix[N, M] move;
  
  matrix[N, N] vcv;  // phylogenetic covariance
  matrix[N, N] adj;  // adjacency matrix
}
transformed data {
//  vector[N] zeroes;
  matrix[N, N] vcv_inv;
//  matrix[N, N] DS;

//  for(i in 1:N) zeroes[i] <- 0;

  vcv_inv <- inverse(vcv);
//  for(i in 1:N)
//    for(j in 1:N)
//      DS[i, j] <- if_else(i==j, sum(row(adj, i)), 0.0);
}
parameters {
  real beta_inter;
  real beta_mass;
  vector[M] beta_move;
  vector[D] beta_diet;

  real<lower=0> sigma_phy;
  vector[N] phy;

//  real<lower=0> sigma_spt;  // prec of spatial
//  real<lower=0,upper=1> p;  // stength of spatial
//  vector[N] spatial;
}
transformed parameters {
  real<lower=0> sig_phy_sq;
//  real<lower=0> tau;

  sig_phy_sq <- sigma_phy * sigma_phy;
//  tau <- 1 / sigma_spt;
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

  // phylogenetic effect
  sigma_phy ~ cauchy(0, 2.5);
  // non-constant part of log(det(sigma_phy * vcv) ^ -0.5
  increment_log_prob(-0.5 * N * log(sig_phy_sq));
  // log of kernal of mulinorm
  increment_log_prob(-(transpose(phy) * vcv_inv * phy) / (2 * sig_phy_sq));
  
  // spatial effect
//  sigma_spt ~ cauchy(0, 2.5);
//  p ~ uniform(0, 1);
//  spatial ~ multi_normal_prec(zeroes, (tau * tau) * (DS - p * adj));

  mu <- (beta_inter + beta_mass * mass + 
        diet * beta_diet + move * beta_move +
        phy + log(off));

  degree ~ poisson_log(mu);
}
generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  mu <- (beta_inter + beta_mass * mass + 
         diet * beta_diet + move * beta_move + 
         phy + log(off));

  for(i in 1:N) {
    log_lik[i] <- poisson_log_log(degree[i], mu[i]);
  }
}

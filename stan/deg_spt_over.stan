data {
  int<lower=0> N;  // sample size
  int D;  // number of diet categories
  int M;  // number of move categories
  int<lower=0> degree[N];  // # i co-occurs with
  vector[N] mass;
  matrix[N, D] diet;
  matrix[N, M] move;
  
  matrix[N, N] adj;  // adjacency matrix
}
transformed data {
  matrix[N, N] DS;
  vector[N] zeros;

  for(i in 1:N)
    for(j in 1:N)
      DS[i, j] <- if_else(i==j, sum(row(adj, i)), 0.0);
  
  for(i in 1:N)
    zeros[i] <- 1;
}
parameters {
  real beta_inter;
  real beta_mass;
  vector[M] beta_move;
  vector[D] beta_diet;

  real<lower=0> tau;  // prec of spatial
  real<lower=0,upper=1> p;  // stength of spatial
  vector[N] spatial;

  real<lower = 0> omega;
}
transformed parameters {
  real<lower=0> sigma;
  real<lower=0> phi;
  
  sigma <- 1 / tau;

  phi <- 1 / omega;
}
model {
  vector[N] mu;
  real spatial_mean;
  vector[N] spatial_std;
  
  beta_inter ~ normal(0, 10);
  beta_mass ~ normal(0, 10);
  for(i in 1:M) {
    beta_move[i] ~ normal(0, 10);
  }
  for(i in 1:D) {
    beta_diet[i] ~ normal(0, 10);
  }

  omega ~ cauchy(0, 2.5);

  // spatial effect
  sigma ~ cauchy(0, 2.5);
  spatial ~ multi_norm_prec(0 vector, (tau * tau) * (DS - p * adj));

  spatial_mean <- mean(spatial);
  spatial_std <- spatial - spatial_mean;

  mu <- (beta_inter + beta_mass * mass + 
        diet * beta_diet + move * beta_move +
        spatial_std);

  degree ~ neg_binomial_2_log(mu, phi);
}

data {
  int<lower=0> N;  // sample size
  int D;  // number of diet categories
  int M;  // number of move categories
  vector[N] off;
  int<lower=0> degree[N];  // # i co-occurs with
  vector[N] mass;
  matrix[N, D] diet;
  matrix[N, M] move;
}
parameters {
  real beta_inter;
  real beta_mass;
  vector[M] beta_move;
  vector[D] beta_diet;

  real<lower = 0> phi;
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

  phi ~ cauchy(0, 2.5);

  mu <- (beta_inter + beta_mass * mass + 
        diet * beta_diet + move * beta_move +
        log(off));

  degree ~ neg_binomial_2_log(mu, phi);
}
generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  mu <- (beta_inter + beta_mass * mass + 
         diet * beta_diet + move * beta_move + 
         log(off));

  for(i in 1:N) {
    log_lik[i] <- neg_binomial_2_log_log(degree[i], mu[i], phi);
  }
}


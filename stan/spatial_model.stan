data {
  int<lower = 1> N;
  int diversity[N];
  matrix[N,N] adj;
  vector[N] neighbors;
}
transformed data {
  vector[N] init;
  vector[N] eigen;
  vector[N] neighbor_trans;
  matrix[N,N] mat;
  
  for(n in 1:N) {
    init[n] <- 0;
    neighbor_trans[n] <- neighbors[n] ^ (- 0.5);
  }

  mat <- diag_matrix(neighbor_trans) * adj * diag_matrix(neighbor_trans);
  eigen <- eigenvalues_sym(mat);
}
parameters {
  real intercept;

  vector[N] spatial;  // storage
  real<lower=0> sigma;  // spatial variance
  real<lower=min(eigen),upper=max(eigen)> p;  // spatial "strength"
}
transformed parameters {
  real<lower=0> tau;  // spatial precision

  tau <- 1 / sigma;
}
model {
  // CAR prior
  spatial ~ multi_normal_prec(init, tau * (diag_matrix(neighbors) - p * adj));
  sigma ~ cauchy(0, 1);
  intercept ~ normal(0, 5);

  for(n in 1:N) {
    diversity[n] ~ poisson_log(intercept + spatial[n]);
  }
}
generated quantities {
  vector[N] log_lik;

  for(n in 1:N) {
    log_lik[n] <- poisson_log_log(diversity[n], intercept + spatial[n]);
  }
}

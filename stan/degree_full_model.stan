data {
  int<lower=0> N;  // sample size
  int D;  // number of diet categories
  int M;  // number of move categories
  int<lower=0> degree[N];  // # i co-occurs with
  vector[N] mass;
  matrix[N, D] diet;
  matrix[N, M] move;
  
  matrix[N, N] vcv;  // phylogenetic covariance
  matrix[N, N] adj;  // adjacency matrix
}
transformed data {
  matrix[N, N] vcv_inv;
  matrix[N, N] DS;

  vcv_inv <- inverse(vcv);
  for(i in 1:N)
    for(j in 1:N)
      DS[i, j] <- if_else(i==j, sum(row(adj, i)), 0.0);
}
parameters {
  real beta_inter;
  real beta_mass;
  vector[M] beta_move;
  vector[D] beta_diet;

  real<lower=0> sigma_phy;
  vector[N] phy;

  real<lower=0> tau;  // var of spatial
  real<lower=0,upper=1> p;  // stength of spatial
  vector[N] spatial;
}
transformed parameters {
  // make a variance
  real<lower=0> sq_sigma;
  real<lower=0> tau_sq;
  sq_sigma <- sigma_phy^2;  
  tau_sq <- tau^2;
}
model {
  real spatial_mean;
  vector[N] spatial_std;
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
  increment_log_prob(-0.5 * N * log(sq_sigma));
  // log of kernal of mulinorm
  increment_log_prob(-(transpose(phy) * vcv_inv * phy) / (2 * sq_sigma));
  
  // spatial effect
  tau ~ cauchy(0, 2.5);
  increment_log_prob(-0.5 / (tau_sq) * (transpose(spatial) * 
                     DS * spatial - p * (transpose(spatial) * 
                     adj * spatial)));
  increment_log_prob(-0.5 * N * log(tau_sq) + 0.5 * log(determinant(DS - p * adj)));
  
  spatial_mean <- mean(spatial);  // sum to zero constraint
  spatial_std <- spatial - spatial_mean;


  mu <- (beta_inter + beta_mass * mass + 
        diet * beta_diet + move * beta_move +
        phy + spatial);

  degree ~ poisson_log(mu);
}



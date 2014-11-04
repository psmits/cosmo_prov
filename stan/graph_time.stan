data {
  int T;  // number of time bins
  vector[T] N;  // number of taxa per bin
  vector[T] L;  // number of localities per bin
  vector[T] C;  // code length per bin
  int np;  // number of predictor categories
  matrix[T, np] preds; 
}
parameters {
//  real beta_node;
//  real beta_loc;
//  real beta_edge;
  real intercept;
  vector[np] beta_trait;
  real<lower=0> sigma;
}
model {
//  beta_node ~ student_t(4, 0, 100);
//  beta_loc ~ student_t(4, 0, 100);
// beta_edge ~ student_t(4, 0, 100);
  beta_trait ~ normal(0, 100);
  intercept ~ normal(0, 100);

  C ~ normal(intercept + preds * beta_trait, sigma);
}

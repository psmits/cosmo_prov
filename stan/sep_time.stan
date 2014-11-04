data {
  int T;  // number of time bins
  vector[T] N;  // number of taxa per bin
  vector[T] L;  // number of localities per bin
  vector[T] C;  // code length per bin
  int G;  // number of groups
  int M[T];  // continental membership
  int np;  // number of predictor categories
  row_vector[np] preds[T]; 
}
parameters {
  real intercept[G];
  vector[np] beta_trait[G];
  real<lower=0> sigma;
}
model {
  for (d in 1:np) {
    for (g in 1:G) {
      intercept[g] ~ normal(0, 100);
      beta_trait[g,d] ~ normal(0, 100);
    }
  }
  
  for(n in 1:T) {
    C[n] ~ normal(intercept[M[n]] + preds[n] * beta_trait[M[n]], sigma);
  }
}

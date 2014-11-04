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
  real mu[np];
  real<lower=0> hsig[np];
  real i_m;
  real<lower=0> i_v;
  
  vector[np] beta_trait[G];
  real intercept[G];
  
  real<lower=0> sigma;
}
model {
  intercept ~ normal(i_m, i_v);

  for (d in 1:np) {
    mu[d] ~ normal(0, 100);
    hsig[d] ~ cauchy(0, 2.5);
    for (g in 1:G) {
      beta_trait[g,d] ~ normal(mu[d], hsig[d]);
    }
  }
  
  for(n in 1:T) {
    C[n] ~ normal(intercept[M[n]] + preds[n] * beta_trait[M[n]], sigma);
  }
}

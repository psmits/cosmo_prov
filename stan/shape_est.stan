data {
  int N;
  real y[N];
}
parameters {
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  alpha ~ cauchy(0, 2.5);
  sigma ~ cauchy(0, 2.5);

  y ~ weibull(alpha, sigma);
}

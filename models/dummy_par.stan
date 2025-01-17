functions {
  real partial_sum_lpdf(array[] real slice_y, int start, int end, vector mu, real sigma) {
    return normal_lupdf(slice_y | mu[start:end], sigma);
  }
}
data {
  int<lower=0> N ;
  int<lower=0> P ;
  vector[N] y ;
  vector[N] year ;
  array[N] int<lower=0, upper=P> plot;
}
parameters {
  vector[N] alpha_p ;
  vector[N] beta_p ;
  real alpha ;
  real beta ;
  real<lower=0> sigma ;
  real<lower=0> sigma_a ;
  real<lower=0> sigma_b ;
}
transformed parameters {
  vector[N] mu = alpha_p[plot] + beta_p[plot] .* year ;
}
model {
  int grainsize = 1;
  target += reduce_sum(partial_sum_lpdf, y, grainsize, mu, sigma) ;
  alpha_p ~ normal(alpha, sigma_a) ;
  beta_p ~ normal(beta, sigma_b) ;
}
generated quantities {
  vector[N] log_lik ;
  for(n in 1:N)
    log_lik[n] = normal_cdf(y[n] | mu[n], sigma) ;
}

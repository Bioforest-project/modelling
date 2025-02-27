data {
  int<lower=0> n_obs;
  int<lower=0> n_site;
  int<lower=0> n_plot;
  vector[n_obs] y;
  vector[n_obs] time;
  array[n_obs] int<lower=0, upper=n_site> site;
  array[n_obs] int<lower=0, upper=n_plot> plot;
  array[n_plot] int<lower=0, upper=n_site> site_plot;
  array[2] real thetaInf_bounds;
  vector[n_site] mu_thetaInf_s;
  vector[n_site] sigma_thetaInf_s;
  array[2] real lambda_bounds;
  array[2] real dist_bounds;
  array[2] real delta_bounds;
  array[2] real tau_bounds;
  int<lower=1> time_max;
}
transformed data {
  vector[time_max] time_pred;
  for(t in 1:time_max)
    time_pred[t] = t;
}
parameters {
  real<lower=dist_bounds[1], upper=dist_bounds[2]> mu_dist; // starting point
  real<lower=0> sigma_dist;
  vector<lower=dist_bounds[1], upper=dist_bounds[2]>[n_plot] dist_p;
  real<lower=lambda_bounds[1], upper=lambda_bounds[2]> mu_lambda; // recovery rate
  real<lower=0> sigma_lambda;
  vector<lower=lambda_bounds[1], upper=lambda_bounds[2]>[n_plot] lambda_p;
  vector<lower=thetaInf_bounds[1], upper=thetaInf_bounds[2]>[n_site] thetaInf_s;
  real<lower=0> sigma;
  real<lower=delta_bounds[1], upper=delta_bounds[2]> mu_delta; // str var
  vector<lower=-dist_p, upper=delta_bounds[2]>[n_plot] delta_p;
  real<lower=0> sigma_delta;
  real<lower=tau_bounds[1]-3, upper=tau_bounds[2]-3> mu_tau_0; //str time
  vector<lower=tau_bounds[1]-3, upper=tau_bounds[2]-3>[n_site] tau_0_s;
  real<lower=0> sigma_tau;
}
transformed parameters {
  vector[n_obs] ltp = 1 - exp(-lambda_p[plot] .* time);
  vector[n_obs] stp = delta_p[plot] .* 
                      (time ./ tau_0_s[site] .* exp( 1 - time ./ tau_0_s[site] )) .*
                      (time ./ tau_0_s[site] .* exp( 1 - time ./ tau_0_s[site] ));
  vector[n_obs] mu = thetaInf_s[site] .* (dist_p[plot] + ltp .* (1 - dist_p[plot]) + stp);
}
model {
  log(y) ~ normal(log(mu), sigma);
  dist_p ~ normal(mu_dist, sigma_dist);
  for(s in 1:n_site)
    thetaInf_s[s] ~ normal(mu_thetaInf_s[s], sigma_thetaInf_s[s]);
  lambda_p ~ normal(mu_lambda, sigma_lambda);
  delta_p ~ normal(mu_delta, sigma_delta);
  tau_0_s ~ normal(mu_tau_0, sigma_tau);
  sigma ~ std_normal();
  sigma_dist ~ std_normal();
  sigma_thetaInf_s ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
}
generated quantities {
  matrix[time_max, n_plot] y_pred;
  vector[n_plot] y15_rel_p;
  real mu_tau = mu_tau_0 + 3;
  vector[n_site] tau_s = tau_0_s + 3;
  real mu_t90 = log(10) / mu_lambda + 3;
  vector[n_plot] t90_p = log(10) / lambda_p + 3;
  real mu_delta_pct = mu_delta*100;
  vector[n_plot] delta_pct_p = delta_p*100;
  real mu_dist_pct = mu_dist*100;
  vector[n_plot] dist_pct_p = (1-dist_p)*100;
  for(p in 1:n_plot)
    y_pred[,p] = thetaInf_s[site_plot[p]] .* 
                  (dist_p[p] + 
                  (1 - exp(-lambda_p[p] * time_pred)) 
                  .* (1 - dist_p[p]) + 
                  (delta_p[p] * 
                 (time_pred / tau_0_s[site_plot[p]] .* 
                 exp( 1 - time_pred ./ tau_0_s[site_plot[p]] )) .*
                 (time_pred / tau_0_s[site_plot[p]] .* 
                 exp( 1 - time_pred ./ tau_0_s[site_plot[p]] ))));
   y15_rel_p = to_vector(y_pred[(15+3),]) ./ thetaInf_s[site_plot];
}

data {
  int<lower=0> n_rec; // obs reco
  int<lower=0> n_old; // obs old
  int<lower=0> n_pre; // obs prelog
  int<lower=0> n_site;
  int<lower=0> n_plot_rec;
  
  int<lower=0> n_plot_old;
  
  vector[n_rec] y_rec;
  vector[n_old] y_old;
  vector[n_pre] y_pre;
  vector[n_rec] time;
  array[n_rec] int<lower=0, upper=n_site> site_rec;
  array[n_old] int<lower=0, upper=n_site> site_old;
  array[n_pre] int<lower=0, upper=n_site> site_pre;
  array[n_rec] int<lower=0, upper=n_plot_rec> plot_rec;
  
  array[n_old] int<lower=0, upper=n_plot_old> plot_old;
  array[n_pre] int<lower=0, upper=n_plot_rec> plot_pre;
  
  array[n_plot_rec] int<lower=0, upper=n_site> site_plot_rec;
  
  array[n_plot_old] int<lower=0, upper=n_site> site_plot_old;
  
  array[2] real thetaInf_bounds;
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
  // real<lower=dist_bounds[1], upper=dist_bounds[2]> mu_dist; // starting point
  // real<lower=0> sigma_dist;
  vector<lower=dist_bounds[1], upper=dist_bounds[2]>[n_plot_rec] dist_p;
  // real<lower=lambda_bounds[1], upper=lambda_bounds[2]> mu_lambda; // recovery rate
  // real<lower=0> sigma_lambda;
  vector<lower=lambda_bounds[1], upper=lambda_bounds[2]>[n_plot_rec] lambda_p;
  // real<lower=thetaInf_bounds[1], upper=thetaInf_bounds[2]> mu_thetaInf; // ending point
  // real<lower=0> sigma_thetaInf_s;
  vector<lower=thetaInf_bounds[1], upper=thetaInf_bounds[2]>[n_site] thetaInf_s;
  
  real<lower=0> sigma_noise;
  vector[n_plot_old] noiseold_p;
  vector[n_plot_rec] noiserec_p;
  
  real<lower=0> sigma_old;
  real<lower=0> sigma_pre;
  real<lower=0> sigma_rec;
  real<lower=delta_bounds[1], upper=delta_bounds[2]> mu_delta; // str var
  vector<lower=delta_bounds[1], upper=delta_bounds[2]>[n_plot_rec] delta_p;
  real<lower=0> sigma_delta;
  real<lower=tau_bounds[1]-3, upper=tau_bounds[2]-3> mu_tau_0; //str time
  vector<lower=tau_bounds[1]-3, upper=tau_bounds[2]-3>[n_site] tau_0_s;
  real<lower=0> sigma_tau;
}
transformed parameters {
  vector[n_plot_rec] theta0_p = thetaInf_s[site_plot_rec] .* dist_p;
  
  vector[n_old] mu_old = exp(log(thetaInf_s[site_old]) + noiseold_p[plot_old]);
  
  vector[n_pre] mu_pre = exp(log(thetaInf_s[site_pre]) + noiserec_p[plot_pre]);
  
  vector[n_rec] ltp_rec = 1 - exp(-lambda_p[plot_rec] .* time);
  vector[n_rec] stp_rec = delta_p[plot_rec] .* 
                           (time ./ tau_0_s[site_rec] .* exp( 1 - time ./ tau_0_s[site_rec] )) .*
                           (time ./ tau_0_s[site_rec] .* exp( 1 - time ./ tau_0_s[site_rec] ));
  vector[n_rec] mu_rec = exp(log(theta0_p[plot_rec] + 
                         ltp_rec .* (thetaInf_s[site_rec] - theta0_p[plot_rec]) + 
                         stp_rec .* thetaInf_s[site_rec]) + noiserec_p[plot_rec]);
}
model {
  log(y_old) ~ normal(log(mu_old), sigma_old);
  log(y_pre) ~ normal(log(mu_pre), sigma_pre);
  log(y_rec) ~ normal(log(mu_rec), sigma_rec);
  // dist_p ~ cauchy(mu_dist, sigma_dist);
  // thetaInf_s ~ cauchy(mu_thetaInf, sigma_thetaInf_s);
  
  noiseold_p ~ normal(0, sigma_noise);
  noiserec_p ~ normal(0, sigma_noise);
  
  // lambda_p ~ cauchy(mu_lambda, sigma_lambda);
  delta_p ~ cauchy(mu_delta, sigma_delta);
  tau_0_s ~ cauchy(mu_tau_0, sigma_tau);
  sigma_old ~ std_normal();
  sigma_pre ~ std_normal();
  sigma_rec ~ std_normal();
  // sigma_dist ~ std_normal();
  // sigma_thetaInf_s ~ std_normal();
  
  sigma_noise ~ std_normal();
  
  // sigma_lambda ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
}
generated quantities {
  matrix[time_max, n_plot_rec] y_pred;
  vector[n_plot_rec] y15_rel_p;
  real mu_tau = mu_tau_0 + 3;
  vector[n_site] tau_s = tau_0_s + 3;
  // real mu_t90 = log(10) / mu_lambda + 3;
  vector[n_plot_rec] t90_p = log(10) / lambda_p + 3;
  // real mu_delta_pct = mu_delta*100;
  vector[n_plot_rec] delta_pct_p = delta_p*100;
  // real mu_dist_pct = mu_dist*100;
  vector[n_plot_rec] dist_pct_p = (1-dist_p)*100;
  for(p in 1:n_plot_rec)
    y_pred[,p] = theta0_p[p] + 
                 (thetaInf_s[site_plot_rec[p]] - theta0_p[p]) * 
                 (1 - exp(-lambda_p[p] * time_pred)) +
                 thetaInf_s[site_plot_rec[p]] * 
                 (delta_p[p] * 
                 (time_pred / tau_0_s[site_plot_rec[p]] .* 
                 exp( 1 - time_pred ./ tau_0_s[site_plot_rec[p]] )) .*
                 (time_pred / tau_0_s[site_plot_rec[p]] .* 
                 exp( 1 - time_pred ./ tau_0_s[site_plot_rec[p]] )));
   y15_rel_p = to_vector(y_pred[(15+3),]) ./ thetaInf_s[site_plot_rec];
}

data {
  int<lower=0> n_rec; // obs reco
  int<lower=0> n_old; // obs old
  int<lower=0> n_pre; // obs prelog
  int<lower=0> n_site;
  int<lower=0> n_plot_rec;
  vector[n_rec] y_rec;
  vector[n_old] y_old;
  vector[n_pre] y_pre;
  vector[n_rec] time;
  array[n_rec] int<lower=0, upper=n_site> site_rec;
  array[n_old] int<lower=0, upper=n_site> site_old;
  array[n_pre] int<lower=0, upper=n_site> site_pre;
  array[n_rec] int<lower=0, upper=n_plot_rec> plot_rec;
  array[n_plot_rec] int<lower=0, upper=n_site> site_plot_rec;
  array[2] real mu_thetaInf_bounds;
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
  real<lower=dist_bounds[1], upper=dist_bounds[2]> mu_dist; // starting point
  real<lower=0> sigma_dist;
  vector<lower=dist_bounds[1], upper=dist_bounds[2]>[n_plot_rec] dist_p;
  real<lower=lambda_bounds[1], upper=lambda_bounds[2]> mu_lambda; // recovery rate
  real<lower=0> sigma_lambda;
  vector<lower=lambda_bounds[1], upper=lambda_bounds[2]>[n_plot_rec] lambda_p;
  real<lower=mu_thetaInf_bounds[1], upper=mu_thetaInf_bounds[2]> mu_thetaInf; // ending point
  real<lower=0> sigma_thetaInf;
  vector<lower=thetaInf_bounds[1], upper=thetaInf_bounds[2]>[n_site] thetaInf_s;
  real<lower=0> sigma_old;
  real<lower=0> sigma_pre;
  real<lower=0> sigma_rec;
}
transformed parameters {
  vector[n_plot_rec] theta0_p = thetaInf_s[site_plot_rec] .* dist_p;
  vector[n_old] mu_old = thetaInf_s[site_old];
  vector[n_pre] mu_pre = thetaInf_s[site_pre];
  vector[n_rec] ltp_rec = 1 - exp(-lambda_p[plot_rec] .* time);
  vector[n_rec] mu_rec = theta0_p[plot_rec] + 
                         (thetaInf_s[site_rec] - theta0_p[plot_rec]) .* 
                         (ltp_rec);
}
model {
  log(y_old) ~ normal(log(mu_old), sigma_old);
  log(y_pre) ~ normal(log(mu_pre), sigma_pre);
  log(y_rec) ~ normal(log(mu_rec), sigma_rec);
  dist_p ~ cauchy(mu_dist, sigma_dist);
  thetaInf_s ~ cauchy(mu_thetaInf, sigma_thetaInf);
  lambda_p ~ cauchy(mu_lambda, sigma_lambda);
  sigma_old ~ std_normal();
  sigma_pre ~ std_normal();
  sigma_rec ~ std_normal();
  sigma_dist ~ std_normal();
  sigma_thetaInf ~ std_normal();
  sigma_lambda ~ std_normal();
}
generated quantities {
  matrix[time_max, n_plot_rec] y_pred;
  for(p in 1:n_plot_rec)
    y_pred[,p] = theta0_p[p] + 
                 (thetaInf_s[site_plot_rec[p]] - theta0_p[p]) * 
                 (1 - exp(-lambda_p[p] * time_pred));
}

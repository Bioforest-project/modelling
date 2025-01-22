data {
  int<lower=0> n_rec; // obs reco
  int<lower=0> n_old; // obs old
  int<lower=0> n_pre; // obs prelog
  int<lower=0> n_site;
  int<lower=0> n_plot_rec;
  vector[n_rec] stem_rec;
  vector[n_old] stem_old;
  vector[n_pre] stem_pre;
  vector[n_rec] time;
  array[n_rec] int<lower=0, upper=n_site> site_rec;
  array[n_old] int<lower=0, upper=n_site> site_old;
  array[n_pre] int<lower=0, upper=n_site> site_pre;
  array[n_rec] int<lower=0, upper=n_plot_rec> plot_rec;
  array[n_plot_rec] int<lower=0, upper=n_site> site_plot_rec;
}
parameters {
  real<lower=0.5, upper=1.49> mu_dist ; // starting point
  real<lower=0> sigma_dist ;
  vector<lower=0.5, upper=1.49>[n_plot_rec] dist_p ;
  real<lower=0, upper=0.5> mu_lambda ; // recovery rate
  real<lower=0> sigma_lambda ;
  vector<lower=0, upper=0.5>[n_plot_rec] lambda_p ;
  real<lower=1, upper=2000> mu_thetaInf ; // ending point
  real<lower=0> sigma_thetaInf ;
  vector<lower=1, upper=2000>[n_site] thetaInf_s ;
  real<lower=0> sigma_old ;
  real<lower=0> sigma_pre ;
  real<lower=0> sigma_rec ;
  real<lower=0, upper=2> mu_delta ; // str var
  vector<lower=0, upper=2>[n_plot_rec] delta_p ;
  real<lower=0> sigma_delta ;
  real<lower=5, upper=30> mu_tau ; //str time
  vector<lower=5, upper=30>[n_site] tau_s ;
  real<lower=0> sigma_tau ;
}
transformed parameters {
  vector[n_plot_rec] theta0_p = thetaInf_s[site_plot_rec] .* dist_p;
  vector[n_old] mu_old = thetaInf_s[site_old];
  vector[n_pre] mu_pre = thetaInf_s[site_pre];
  vector[n_rec] ltp_rec = 1 - exp(-lambda_p[plot_rec] .* time);
  vector[n_rec] stp_rec = delta_p[plot_rec] .* 
                           (time ./ tau_s[site_rec] .* exp( 1 - time ./ tau_s[site_rec] )) .*
                           (time ./ tau_s[site_rec] .* exp( 1 - time ./ tau_s[site_rec] ));
  vector[n_rec] mu_rec = theta0_p[plot_rec] + 
                         (thetaInf_s[site_rec] - theta0_p[plot_rec]) .* 
                         (ltp_rec+stp_rec) ;
}
model {
  log(stem_old) ~ normal(log(mu_old), sigma_old) ;
  log(stem_pre) ~ normal(log(mu_pre), sigma_pre) ;
  log(stem_rec) ~ normal(log(mu_rec), sigma_rec) ;
  dist_p ~ lognormal(log(mu_dist), sigma_dist) ;
  thetaInf_s ~ lognormal(log(mu_thetaInf), sigma_thetaInf) ;
  lambda_p ~ lognormal(log(mu_lambda), sigma_lambda) ;
  delta_p ~ lognormal(log(mu_delta), sigma_delta) ;
  tau_s ~ lognormal(log(mu_tau), sigma_tau) ;
  sigma_old ~ std_normal();
  sigma_pre ~ std_normal();
  sigma_rec ~ std_normal();
  sigma_dist ~ std_normal();
  sigma_thetaInf ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
}
generated quantities {
  vector[n_rec] log_lik ;
  for(n in 1:n_rec)
    log_lik[n] = normal_cdf(log(stem_rec[n]) | mu_rec[n], sigma_rec) ;
}

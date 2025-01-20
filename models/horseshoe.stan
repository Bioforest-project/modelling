// horseshoe example https://github.com/to-mi/stan-survival-shrinkage/blob/master/wei_hs.stan
// https://statmodeling.stat.columbia.edu/2015/02/17/bayesian-survival-analysis-horseshoe-priors/
// https://github.com/mcol/hsstan
// https://github.com/mcol/hsstan/blob/master/inst/stan/hs.stan
// https://github.com/mcol/hsstan/blob/master/inst/stan/chunks/hs.fun
// https://avehtari.github.io/modelselection/bodyfat.html
functions {
  vector hs_prior_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real nu) {
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);
    return (r1_global * sqrt(r2_global)) * r1_local .* sqrt_vec(r2_local);
  }
}
data {
  int<lower=0> n_rec; // obs reco
  int<lower=0> n_old; // obs old
  int<lower=0> n_pre; // obs prelog
  int<lower=0> n_site;
  int<lower=0> n_plot;
  vector[n_rec] stem_rec;
  vector[n_old] stem_old;
  vector[n_pre] stem_pre;
  vector[n_rec] time;
  array[n_rec] int<lower=0, upper=n_site> site_rec;
  array[n_old] int<lower=0, upper=n_site> site_old;
  array[n_pre] int<lower=0, upper=n_site> site_pre;
  array[n_rec] int<lower=0, upper=n_plot> plot_rec;
  array[n_old] int<lower=0, upper=n_plot> plot_old;
  array[n_pre] int<lower=0, upper=n_plot> plot_pre;
}
parameters {
  real<lower=0, upper=2000> mu_theta0 ; // starting point
  real<lower=0> sigma_theta0 ;
  vector<lower=0, upper=2000>[n_plot] theta0_p ;
  real<lower=0, upper=0.5> mu_lambda ; // recovery rate
  real<lower=0> sigma_lambda ;
  vector<lower=0, upper=0.5>[n_plot] lambda_p ;
  real<lower=0, upper=2000> mu_thetaInf ; // ending point
  real<lower=0> sigma_thetaInf ;
  vector<lower=0, upper=2000>[n_site] thetaInf_s ;
  real<lower=0> sigma_old ;
  real<lower=0> sigma_pre ;
  real<lower=0> sigma_rec ;
  real<lower=0, upper=1> mu_delta ; // str var
  vector<lower=0, upper=1>[n_plot] delta_p ;
  real<lower=0> sigma_delta ;
  real<lower=0, upper=30> mu_tau ; //str time
  vector<lower=0, upper=30>[n_site] tau_s ;
  real<lower=0> sigma_tau ;
}
transformed parameters {
  vector[n_old] mu_old = thetaInf_s[site_old];
  vector[n_pre] mu_pre = thetaInf_s[site_pre];
  vector[n_rec] mu_rec;
  for(k in 1:n_rec)
    mu_rec[k] = theta0_p[plot_rec[k]] + 
                (thetaInf_s[site_rec[k]] - theta0_p[plot_rec[k]]) *
                (1 - exp(-lambda_p[plot_rec[k]] * time[k])) *
                (delta_p[plot_rec[k]] * 
                (time[k]/tau_s[site_rec[k]] * exp( 1 - time[k]/tau_s[site_rec[k]] )) *
                (time[k]/tau_s[site_rec[k]] * exp( 1 - time[k]/tau_s[site_rec[k]] ))
                );
}
model {
  stem_old ~ lognormal(log(mu_old), sigma_old) ;
  stem_pre ~ lognormal(log(mu_pre), sigma_pre) ;
  stem_rec ~ lognormal(log(mu_rec), sigma_rec) ;
  theta0_p ~ lognormal(log(mu_theta0), sigma_theta0) ;
  thetaInf_s ~ lognormal(log(mu_thetaInf), sigma_thetaInf) ;
  lambda_p ~ lognormal(log(mu_lambda), sigma_lambda) ;
}
generated quantities {
  vector[n_rec] log_lik ;
  for(n in 1:n_rec)
    log_lik[n] = lognormal_cdf(stem_rec[n] | log(mu_rec[n]), sigma_rec) ;
}

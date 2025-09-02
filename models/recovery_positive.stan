data {
  int<lower=0> n; // # obs
  int<lower=0> s; // # sites
  int<lower=0> p; // # plots
  vector[n] y; // forest attribute
  vector[n] t; // time
  array[n] int<lower=0, upper=n> site; // site index
  array[n] int<lower=0, upper=n> plot; // plot index
  array[p] int<lower=0, upper=n> site_plot; // site index in plots
  vector[s] mu_theta_s; // mean site asymptotic value in controls
  vector[s] sigma_theta_s; // sd site asymptotic value in controls
}
transformed data {
  vector[55] t_pred;
  for(i in 1:55)
    t_pred[i] = i;
}
parameters {
  vector<lower=0.1, upper=1>[p] phi_p; // disturbance intensity
  vector<lower=0.001, upper=0.5>[p] lambda_p; // recovery rate
  real<lower=0.001, upper=0.5> mu_lambda;
  real<lower=0> sigma_lambda;
  vector<lower=0, upper=2>[p] delta_p; // bump height
  real<lower=0, upper=2> mu_delta;
  real<lower=0> sigma_delta;
  vector<lower=6-3, upper=40-3>[s] tau0_s; // bump time
  real<lower=6-3, upper=40-3> mu_tau0;
  real<lower=0> sigma_tau;
  vector<lower=(mu_theta_s-sigma_theta_s)*0.1, upper=(mu_theta_s+sigma_theta_s)*10>[s] theta_s; // asymptotic value
  real<lower=0> sigma; // residual variation
}
transformed parameters {
  vector[n] ltp = 1 - exp(-lambda_p[plot] .* t);
  vector[n] stp = delta_p[plot] .*
                  (t ./ tau0_s[site] .* exp(1 - t ./ tau0_s[site])) .*
                  (t ./ tau0_s[site] .* exp(1 - t ./ tau0_s[site]));
  vector[n] mu = theta_s[site] .* (phi_p[plot] + ltp.*(1 - phi_p[plot]) + stp);
}
model {
  log(y) ~ normal(log(mu), sigma);
  lambda_p ~ normal(mu_lambda, sigma_lambda);
  delta_p ~ lognormal(mu_delta, sigma_delta);
  tau0_s ~ normal(mu_tau0, sigma_tau);
  for(i in 1:s)
    theta_s[i] ~ normal(mu_theta_s[i], sigma_theta_s[i]);
  sigma ~ std_normal();
  sigma_lambda ~ std_normal();
  sigma_delta ~ std_normal();
  sigma_tau ~ std_normal();
}
generated quantities {
  matrix[55, p] y_pred;
  vector[p] y15_p;
  real mu_tau = mu_tau0 + 3; // tau in real time
  vector[s] tau_s = tau0_s + 3;
  real mu_t90 = log(10) / mu_lambda + 3; // 90pct recovery time
  vector[p] t90_p = log(10) / lambda_p + 3;
  real mu_deltapct = mu_delta*100; // bump intensity in pct
  vector[p] deltapct_p = delta_p*100;
  vector[p] phipct_p = (1-phi_p)*100;
  for(i in 1:p)
    y_pred[,i] = theta_s[site_plot[i]] .* 
                  (phi_p[i] + 
                  (1 - exp(-lambda_p[i] * t_pred)) 
                  .* (1 - phi_p[i]) + 
                  (delta_p[i] * 
                 (t_pred / tau0_s[site_plot[i]] .* 
                 exp( 1 - t_pred ./ tau0_s[site_plot[i]] )) .*
                 (t_pred / tau0_s[site_plot[i]] .* 
                 exp( 1 - t_pred ./ tau0_s[site_plot[i]] ))));
   y15_p = to_vector(y_pred[(15+3),]) ./ theta_s[site_plot];
}

library(data.table)
library(dplyr)
library(rstan) 
library(xts)
library(lubridate)
library(ggplot2)
library(shinystan)
library(loo)

number = c(8.0, 9.0, 9.1, 10.2, 10.4, 10.0, 10.3, 12.2, 12.6, 13.9)
day = rep(1:5, 2)
type = factor(rep(c("A", "B"), each = 5))
bacteria = data.frame(day, number, type)

####  Model 1: A one-way ANOVA using type = “A” as the reference category.

# Create a dummy variable
type_b = c(type == 'B')

# Create model
stanmod = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  vector[n] type_b; //indicator for type b
  real<lower=0> v_sample; // sample variance of y
  real mu_sample; // sample mean of y

}

parameters {
  real mu_0;
  real alpha_b;
  real<lower=0> sigmasq;

}

transformed parameters {
  vector[n] mu;        // mean of observations
  mu = mu_0 + alpha_b*type_b;

}
model {

  // priors
  mu_0 ~ normal(mu_sample, 10 * sqrt(v_sample));
  alpha_b ~ normal(0.0, mu_sample); 
  sigmasq ~ inv_gamma(0.01, 0.01);

  // data distribution
  for(i in 1:n){
    y[i] ~ normal(mu[i], sqrt(sigmasq));
  }
}
generated quantities {
  real Rbsq;          // goodness-of-fit
  vector[n] log_lik;  // log likelihood of data
  Rbsq = 1 - sigmasq/v_sample;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i]|mu[i], sqrt(sigmasq));
}
"

stan_dat = list(n = length(number), y = number,
                type_b = type_b, v_sample = var(number), mu_sample = mean(number))

stan_fit = stan(model_code = stanmod, data = stan_dat, iter = 10^4)
#sso <- launch_shinystan(stan_fit)

# Mean, 0.025 quantile, and 0.975 quantile
summary(stan_fit)$summary[,c("mean", "2.5%", "97.5%")]

summary(stan_fit)$summary

# Gelman-Rubin statistic
summary(stan_fit)$summary[,"Rhat"]

# check convergence with trace plots
stan_trace(stan_fit, c("mu_0", "alpha_b", "sigmasq"))

# Check the ACF of draws
stan_ac(stan_fit, c("mu_0", "alpha_b", "sigmasq"))


# WAIC and LOOC
log_lik = extract_log_lik(stan_fit, merge_chains = FALSE)
r_eff = exp(relative_eff(log_lik))

waic(log_lik)
loo(log_lik, r_eff = r_eff)

####  Model 1: A one-way ANOVA using type = “A” as the reference category.

# Create a dummy variable
type_b = c(type == 'B')

# Create model
stanmod = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  vector[n] type_b; //indicator for type b
  real<lower=0> v_sample; // sample variance of y
  real mu_sample; // sample mean of y

}

parameters {
  real mu_0;
  real alpha_b;
  real<lower=0> sigmasq;

}

transformed parameters {
  vector[n] mu;        // mean of observations
  mu = mu_0 + alpha_b*type_b;

}
model {

  // priors
  mu_0 ~ normal(mu_sample, 10 * sqrt(v_sample));
  alpha_b ~ normal(0.0, mu_sample); 
  sigmasq ~ inv_gamma(0.01, 0.01);

  // data distribution
  for(i in 1:n){
    y[i] ~ normal(mu[i], sqrt(sigmasq));
  }
}
generated quantities {
  real Rbsq;          // goodness-of-fit
  vector[n] log_lik;  // log likelihood of data
  Rbsq = 1 - sigmasq/v_sample;
  for (i in 1:n) log_lik[i] = normal_lpdf(y[i]|mu[i], sqrt(sigmasq));
}
"

stan_dat = list(n = length(number), y = number,
                type_b = type_b, v_sample = var(number), mu_sample = mean(number))

stan_fit = stan(model_code = stanmod, data = stan_dat, iter = 10^4)
#sso <- launch_shinystan(stan_fit)

# Mean, 0.025 quantile, and 0.975 quantile
summary(stan_fit)$summary[,c("mean", "2.5%", "97.5%")]

summary(stan_fit)$summary

# Gelman-Rubin statistic
summary(stan_fit)$summary[,"Rhat"]

# check convergence with trace plots
stan_trace(stan_fit, c("mu_0", "alpha_b", "sigmasq"))

# Check the ACF of draws
stan_ac(stan_fit, c("mu_0", "alpha_b", "sigmasq"))


# WAIC and LOOC
log_lik = extract_log_lik(stan_fit, merge_chains = FALSE)
r_eff = exp(relative_eff(log_lik))

waic(log_lik)
loo(log_lik, r_eff = r_eff)


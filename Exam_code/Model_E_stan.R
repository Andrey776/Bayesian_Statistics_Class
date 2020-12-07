library(data.table)
library(truncnorm)
library(coda)
devtools::install_github("stan-dev/loo")
library(loo)
library(rstan) 
library(shinystan)
library(invgamma)
library(bayesplot)
library(posterior)

devtools::install_github("jfrench/bayesutils")

df_covid = bayesutils::covid_dec4

stan_mod_E = "
data {
  int T;
  int deaths[T]; // COVID-19 deaths
  vector[T] log_pop;    // log of U.S. states population,
  vector[T] income; // median income (USD)
  vector[T] bs;     // percentage of the population with bachelor’s degrees.
}

parameters {
  real beta0;
  real beta1;
  real beta2;
}
transformed parameters {
  real log_lambda[T];
  for (i in 1:T) {
   log_lambda[i] = log_pop[i] + beta0 + beta1 * income[i] + beta2 * bs[i];
  }

}

model {
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
 
  for (i in 1:T) {
    deaths[i] ~ poisson_log(log_lambda[i]);
  }
}
generated quantities {
  vector[T] log_lik;

  for (i in 1:T) {
    log_lik[i] = poisson_lpmf(deaths[i] | exp(log_lambda[i]));
  }  
  
}

"

T <- length(df_covid$bs)
set.seed(94)

model_name = "Model_E_"
stan_dat_E = list(T = T, deaths = df_covid$deaths, log_pop = log(df_covid$population), 
                  income = df_covid$income, bs = df_covid$bs)
r_fit_E = stan(model_code = stan_mod_E, data = stan_dat_E, 
               iter = 50000, chains = 2, 
               control = list(max_treedepth = 21))

summary(r_fit_E, pars = c("beta0", "beta1", "beta2"), prob = c(0.025, 0.975))$summary

sso <- launch_shinystan(r_fit_E)
file_name = paste(model_name, as.character(Sys.Date()), ".rda", sep="")
setwd('/Users/AM/Documents/_CU Masters/2020 fall Bayesian_7393/code/Bayesian_Statistics_Class_Code/Exam_code')
#readRDS(r_fit, file = "") 
saveRDS(r_fit_E, file = file_name, compress = "xz") 
waic_loo = function (fit) {
  log_lik = extract_log_lik(fit, merge_chains = FALSE)
  r_eff = exp(relative_eff(log_lik))
  wl = c(waic(log_lik)$waic, loo(log_lik, r_eff = r_eff)$looic)
  return(wl)
}

cat("STAN Model E WAIC", waic_loo(r_fit_E)[1])
cat("STAN Model E LOOIC", waic_loo(r_fit_E)[2])


samples = extract(r_fit_E)
ncycles = length(samples[[1]])
T <- length(df_covid$bs)

# each row of yrep is a sample from the pp distribution
yrep = matrix(0, ncol = T, nrow = ncycles)
for (i in seq_len(T)) {
  log_l = log(df_covid$population[i]) +
    samples$beta0 + 
    samples$beta1 * df_covid$income[i] + 
    samples$beta2 * df_covid$bs[i]
  yrep[, i] = rpois(ncycles, lambda = exp(log_l))
}

ppc_intervals(y = df_covid$deaths, yrep = yrep) + ggtitle("Model A: Posterior predictive intervals to the observed data values") + xlab("State") + ylab("Covid deaths")


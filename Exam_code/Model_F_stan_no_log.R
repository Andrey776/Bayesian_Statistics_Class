stan_mod_F = "
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
  real phi;
}
transformed parameters {
  real eta[T];
  for (i in 1:T) {
   eta[i] = exp(log_pop[i] + beta0 + beta1 * income[i] + beta2 * bs[i]);
  }

}

model {
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  phi ~ gamma(0.01, 0.01);

  for (i in 1:T) {
    deaths[i] ~ neg_binomial_2(exp(eta[i]), phi);
  }
}
generated quantities {
  vector[T] log_lik;

  for (i in 1:T) {
    log_lik[i] = neg_binomial_2_lpmf(deaths[i] | exp(eta[i]),  phi);
  }  
  
}

"

T <- length(df_covid$bs)
set.seed(95)

model_name = "Model_F_"
stan_dat_F = list(T = T, deaths = df_covid$deaths, log_pop = log(df_covid$population), 
                  income = df_covid$income, bs = df_covid$bs)
r_fit_F = stan(model_code = stan_mod_F, data = stan_dat_F, 
               iter = 5000, chains = 2, 
               control = list(max_treedepth = 21))

summary(r_fit_F, pars = c("beta0", "beta1", "beta2", "phi"), prob = c(0.025, 0.975))$summary

#sso <- launch_shinystan(r_fit_F)
file_name = paste(model_name, as.character(Sys.Date()), ".rda", sep="")
setwd('/Users/AM/Documents/_CU Masters/2020 fall Bayesian_7393/code/Bayesian_Statistics_Class_Code/Exam_code')
#readRDS(r_fit, file = "") 
saveRDS(r_fit_F, file = file_name, compress = "xz") 


cat("STAN Model F WAIC", waic_loo(r_fit_F)[1])
cat("STAN Model F LOOIC", waic_loo(r_fit_F)[2])

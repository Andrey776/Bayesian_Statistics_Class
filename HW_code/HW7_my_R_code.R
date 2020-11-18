setwd('/Users/AM/Documents/_CU Masters/2020 fall Bayesian_7393/HW')
load(file = "bank_salary.rda")
options(digits=4)
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mvtnorm)
library(bayesplot)

# obtain the number of observations
n = length(bank_salary$bsalary)
# obtain the sample variance of the data
v = var(bank_salary$bsalary)
N=10^5 # number of iterations
X = cbind(1, bank_salary$male, bank_salary$education, bank_salary$experience, bank_salary$time)

# Create model. 

## Model_1
# Prior distributions:
# beta_j ~ N(0, 10^4) (sigmasq = 10^4)
# tau ~ Gamma(0.01,0.01).
stanmod_1 = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  vector[n] male; //covariate
  vector[n] edu; //covariate
  vector[n] exper; //covariate
  vector[n] t; //covariate
  real<lower=0> v; // sample variance of y  
}
parameters {
  real<lower=0> prec; // tau
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;  
}

transformed parameters{
  real<lower=0> sigma; //get sigma from the precision
  sigma = sqrt(1/prec);
}
model {
  //specify priors
  beta0 ~ normal(0.0, 100);
  beta1 ~ normal(0.0, 100);
  beta2 ~ normal(0.0, 100);
  beta3 ~ normal(0.0, 100);
  beta4 ~ normal(0.0, 100);
  prec ~ gamma(0.01, 0.01);

  // data distribution
  for(i in 1:n){
    y[i] ~ normal(beta0 + beta1*male[i] + beta2*edu[i] + beta3*exper[i] + beta4*t[i], sigma);
  }
}
generated quantities {
  real<lower=0> sigmasq; //get sigmasq from the precision
  real Rbsq;
  sigmasq = 1/prec;
  Rbsq = 1 - sigmasq / v;
}
"
## Model_2
# Prior distributions:
# beta ~ N(0, sigma^2 * c^2 * (X'X)^(-1)) with c^2 = n
# tau ~ Gamma(0.01,0.01).
stanmod_2 = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  matrix[n, 5] X; //covariates
  vector[5] mu0; //prior mean for beta
  cov_matrix[5] V; //V part of Zellner's g-prior
  ////cov_matrix[n] I; //nxn identity matrix
  real<lower=0> csq; //constant for Zellner's g-prior
  real<lower=0> v; //sample variance
}
parameters {
  real<lower=0> prec;
  vector[5] beta;
}
transformed parameters{
  real<lower=0> sigmasq; //get sigmasq from the precision
  ////vector[n] mu;  // mean of responses
  sigmasq = 1/prec;
  ////mu = X * beta;
}
model {
  vector[n] mu;  // mean of responses.  Temporary
  mu = X * beta;
  // specify priors
  beta ~ multi_normal(mu0, sigmasq*csq*V);
  prec ~ gamma(0.01, 0.01);

  // specify data distribution
  // y ~ multi_normal(mu, sigmasq * I);
  // The statement below is equivalent to what is above,
  // but faster since we don't sample from a multivariate
  // normal distribution internally.
  for(i in 1:n) {
    y[i] ~ normal(mu[i], sqrt(sigmasq));
  }
}
generated quantities {
  real Rbsq;
  Rbsq = 1 - sigmasq / v;
}
"

## Model_3
# Create model.  Notice the quotes
stanmod_3 = "
data {
  int<lower=1> n; // number of observations
  vector[n] y; // data
  matrix[n, 5] X; //covariates
  vector[5] mu0; //prior mean for beta
  cov_matrix[5] V; //V part of Zellner's g-prior
  cov_matrix[n] I; //nxn identity matrix
  real<lower=0> csq; //constant for Zellner's g-prior
  real<lower=0> v; //sample variance
}
parameters {
  real<lower=0> prec;
  vector[5] beta;
}
transformed parameters{
  real<lower=0> sigmasq; //get sigmasq from the precision
  sigmasq = 1/prec;
}
model {
  vector[n] mu;  // mean of responses.  Temporary
  mu = X * beta;
  //prior distributions
  beta ~ multi_normal(mu0, sigmasq*csq*V);
  prec ~ gamma(0.01, 0.01);
  // data distribution
  //y ~ multi_normal(X*beta, sigmasq*I);
  for(i in 1:n) {
     y[i] ~ normal(mu[i], sqrt(sigmasq));
  }
}
generated quantities {
  real Rbsq;
  Rbsq = 1 - sigmasq / v;
}
"

# fit models using stan with 4 chains
stan_dat_1 = list(n = n, y = bank_salary$bsalary, male = bank_salary$male,
                edu = bank_salary$education, exper = bank_salary$experience,
                t = bank_salary$time, v = v)
salary_fit_1 = stan(model_code = stanmod_1, data = stan_dat_1, iter = N, chains = 4)

stan_dat_2 = list(n = n, y = bank_salary$bsalary,
                X = X, mu0 = c(0, 0, 0, 0, 0), V = solve(crossprod(X)),
                I = diag(n), csq = n, v = v)
salary_fit_2 = stan(model_code = stanmod_2, data = stan_dat_2, iter = N, chains = 4)

stan_dat_3 = list(n = n, y = bank_salary$bsalary,
                X = X, mu0 = c(0, 0, 0, 0, 0), V = solve(crossprod(X)),
                I = diag(n), csq = 10^3, v = v)
salary_fit_3 = stan(model_code = stanmod_3, data = stan_dat_3, iter = N, chains = 4)

dim(salary_fit_3)# i.e. a warmup of 50,000 iterations and 50,000 saved values  

# calculate the Gelman-Rubin statistic
summary(salary_fit_1)$summary[,"Rhat"]
summary(salary_fit_2)$summary[,"Rhat"]
summary(salary_fit_3)$summary[,"Rhat"]

# check the ACF of draws from the MCMC chain
stan_ac(salary_fit_1)
stan_ac(salary_fit_2)
stan_ac(salary_fit_3)

#effectiv sample size
summary(salary_fit_1)$summary[, "n_eff"]
summary(salary_fit_2)$summary[, "n_eff"]
summary(salary_fit_3)$summary[, "n_eff"]

# traceplots
#traceplot(salary_fit_1)
#traceplot(salary_fit_2)
#traceplot(salary_fit_3)

# output
a = array(0, dim = c(3,6,4), dimnames = list(c("Model_1", "Model_2", "Model_3"),
                                             c("beta_0", "beta_1","beta_2", "beta_3", "beta_4", "sigmasq"),
                                             c("Mean", "Median","2.5%", "97.5%")))
a[1,,] = summary(salary_fit_1)$summary[,c("mean", "50%","2.5%", "97.5%")][c(2,3,4,5,6,8),]
a[2,,] = summary(salary_fit_2)$summary[,c("mean", "50%","2.5%", "97.5%")][c(2,3,4,5,6,7),]
a[3,,] = summary(salary_fit_3)$summary[,c("mean", "50%","2.5%", "97.5%")][c(2,3,4,5,6,7),]

for (i in dimnames(a)[[2]]) {
  cat('\n',i,'\n')
  print(a[,i,])
}
###_2_distribution of Rb^2
fitss = dplyr::bind_rows(cbind(as.data.frame(salary_fit_1), model = "Model_1"),
                          cbind(as.data.frame(salary_fit_2), model = "Model_2"),
                          cbind(as.data.frame(salary_fit_3), model = "Model_3"))
fitss_df = as.data.frame(fitss)
# comparison of Rbsq
ggplot(fitss, aes(x = Rbsq)) + theme_bw() +
  geom_density(aes(fill = model), alpha = 0.4)

###_3_Provide a point estimate
x_names = names(bank_salary)
for (m in 1:3){
  cat("bsalary_hut = ", a[m,1,1], "+ ")
  for (i in 2:4) {
    cat(a[m,i,1], "*", x_names[i], "+ ")
  }
  cat(a[m,5,1], "*", x_names[5])
  cat("\n")
}
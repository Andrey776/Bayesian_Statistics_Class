
stanmod_3a = "
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

stan_dat_3a = list(n = n, y = bank_salary$bsalary,
                  X = X, mu0 = c(0, 0, 0, 0, 0), V = solve(crossprod(X)),
                  I = diag(n), csq = 10^6, v = v)
salary_fit_3a = stan(model_code = stanmod_3a, data = stan_dat_3a, iter = N, chains = 4)

summary(salary_fit_3a)$summary[,c("mean", "50%", "2.5%", "97.5%")][c(2,3,4,5,6,7),]

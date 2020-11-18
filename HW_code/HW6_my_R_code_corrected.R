# source("http://math.ucdenver.edu/~jfrench/data/bayesdists.R")
# devtools::install_github("jfrench/bayesutils")
library(mvtnorm)
library(bayesplot)
library(invgamma)
library(coda)
library(rstan)

set.seed(77)
y = bayesutils::rinvgamma(100, 2, 3)
range(y)
mean(y)
sd(y)
### 1b
log_q_alpha = function (y, a_i, b_i) {
  100*a_i*log(b_i) +
    (a_i+1)*sum(log(1/y)) -
    100 * log(factorial(a_i-1))
}
alpha = seq(1,5, length=50)
b=7

alpha_map = function(b_i) {
  log_q_alpha_v = sapply(alpha, log_q_alpha, y=y, b_i=b_i)
  alpha[which.max(log_q_alpha_v)]
}

a_map = alpha_map(b_i=b)
log_q_alpha_scaled = function(y, a_i, b_i) {log_q_alpha(y, a_i, b_i) -
    log_q_alpha(y, alpha_map(b_i=b_i), b_i)}

proposal = dnorm(alpha, a_map, sqrt(0.04))

log_q_alpha_scaled_v = sapply(alpha, log_q_alpha_scaled, y=y, b_i=b)
plot(alpha, exp(log_q_alpha_scaled_v), ylab = "density",
     type = "l", col = "red",
     ylim= c(0, max(proposal)*1.05))
lines(alpha, proposal, col = "blue")
title("full conditional distribution of alpha; beta = 7")

### 2c
mh_gibbs = function(B, theta_start, jump_parm, alpha_left, alpha_right) { # USE sigma^2 here for normal
  sigma = sqrt(jump_parm[2])
  theta = array(0, c((B+1), 2))
  sum_inv = sum(1/y)
  theta[1,1] = theta_start[1]
  theta[1,2] = theta_start[2]
  mu_current = alpha_map(theta[1,2])

  for (i in 2:dim(theta)[1]) {
    alpha_star = theta[(i-1),1]
    beta_star = rgamma(1, shape = (100*alpha_star+4), rate = (4 + sum_inv))
    theta[i,2] = beta_star
    mu_star = alpha_map(beta_star) # new proposal mean
    alpha_star = rnorm(1, mu_star, sigma) # generate alpha within left-right borders
    while ((alpha_star < alpha_left)|(alpha_star > alpha_right)) {
      alpha_star = rnorm(1, mu_star, sigma)
    }
    num_logr = log_q_alpha_scaled(y, a_i=alpha_star, b_i=beta_star) -
      dnorm(alpha_star, mu_star, sigma, log = TRUE)
    den_logr = (log_q_alpha_scaled(y, a_i=theta[i-1, 1], b_i=beta_star) -
      dnorm(theta[(i-1),1], mu_star, sigma, log = TRUE))
    logr = num_logr - den_logr
    if (log(runif(1)) <= min(logr, 0)) {
      theta[i,1] = alpha_star
      mu_current = mu_star
    } else {
      theta[i,1] = theta[(i - 1), 1]
    }
  }
  return(theta)
}
set.seed(77)
B = 100000
keep = (B/2 + 1):(B + 1)
jpar = c(a_map, 0.022)
chain1 = mh_gibbs(B, theta_start = c(2.3,4), jump_parm = jpar,1,5)
chain2 = mh_gibbs(B, theta_start = c(2.8,5), jump_parm = jpar,1,5)
chain3 = mh_gibbs(B, theta_start = c(3.3,6), jump_parm = jpar,1,5)
chain4 = mh_gibbs(B, theta_start = c(3.8,7), jump_parm = jpar,1,5)
chain5 = mh_gibbs(B, theta_start = c(3.9,8), jump_parm = jpar,1,5)

chain = rbind(chain1[keep,], chain2[keep,], chain3[keep,], chain4[keep,], chain5[keep,])
mc_a = mcmc.list(mcmc(chain1[keep,1]), mcmc(chain2[keep,1]),
                 mcmc(chain3[keep,1]), mcmc(chain4[keep,1]),
                 mcmc(chain5[keep,1]))
summary (mc_a)
mc_b= mcmc.list(mcmc(chain1[keep,2]), mcmc(chain2[keep,2]),
                 mcmc(chain3[keep,2]), mcmc(chain4[keep,2]),
                 mcmc(chain5[keep,2]))
summary (mc_b)
plot(density(chain[,1]), main = "", xlab = "alpha", xlim = c(1, 5))
plot(density(chain[,2]), main = "", xlab = "beta")

### 2a
stanmod = "
data {
  int<lower = 1> n;
  vector<lower = 0>[n] y;
  real<lower = 0> beta_gamma_a;
  real<lower = 0> beta_gamma_b;
  real alpha_unif_left;
  real alpha_unif_right;
}
parameters {
  real<lower = 0> alpha;
  real<lower = 0> beta;
}
model {
  y ~ inv_gamma(alpha, beta);
  alpha ~ uniform(alpha_unif_left, alpha_unif_right);
  beta ~ gamma(beta_gamma_a, beta_gamma_b);
}"
stan_dat = list(n = length(y), y = y,
                alpha_unif_left = 1, alpha_unif_right = 5,
                beta_gamma_a = 4, beta_gamma_b = 4)
fit = stan(model_code = stanmod, data = stan_dat, iter = 100000, chains = 5)
codasamples = As.mcmc.list(fit)
summary(codasamples, quantiles = c(.025,.25,.5,.75, .975))
fit_array = as.array(fit)
dim(fit_array)

### 2b
# summary of fit object
summary(fit)

# rstan plotting functions
# posterior intervals and point estimates
# traceplot of chains
stan_trace(fit)
# density plot of theta
stan_dens(fit, "alpha")
stan_dens(fit, "beta")

# ACF plot of chain
stan_ac(fit, "alpha")
stan_ac(fit, "beta")

# coda plots
# convert samples to coda object
codasamples = As.mcmc.list(fit)

# trace plots of posterior samples
coda::traceplot(codasamples)

# density plots of samples
densplot(codasamples)

# summarize codasamples
summary(codasamples, quantiles = c(.025,.25,.5,.75, .975))

# bayesplot plots
posterior = as.array(fit) # convert to array for plotting
# central posterior interval (median shown as point, by default)
# shows 50% and 95% posterior intervals, by default
mcmc_intervals(posterior, pars = "theta")
# area plot of each parameter (shows 50% credible interval by default)
# shows median as vertical line
mcmc_areas(posterior, pars = "theta")
mcmc_areas(posterior, pars = "ytilde")
# histogram
mcmc_hist(posterior, pars = c("theta", "ytilde"))
# histogram by chain
mcmc_hist_by_chain(posterior, pars = c("theta", "ytilde"))
# density plot
mcmc_dens(posterior, pars = c("theta", "ytilde"))
# overlay density for each chain
mcmc_dens_overlay(posterior, pars = c("theta", "ytilde"))
# violin plots
mcmc_violin(posterior, pars = c("theta", "ytilde"))
# trace plots
mcmc_trace(posterior, pars = "theta")
# highlight a single chain
mcmc_trace_highlight(posterior, pars = "theta", highlight = 1)


# source("http://math.ucdenver.edu/~jfrench/data/bayesdists.R")
# devtools::install_github("jfrench/bayesutils")
library(mvtnorm)
library(bayesplot)
library(invgamma)
library(coda)

### 1a
set.seed(77)
y = bayesutils::rinvgamma(100, 2, 3)
range(y)
mean(y)
sd(y)
### 1b option A
log_b_posterior_1 = function(y, alpha, beta_i) { 
  sum(dinvgamma(y, alpha, beta_i, log = TRUE)) +
    dgamma(beta_i, 4, 4, log = TRUE)
}# beta_i - just one point
beta = seq(3, 6, length = 1000)

log_b_posterior = sapply(beta, log_b_posterior_1, y = y, alpha = 3)
log_b_posterior = log_b_posterior - max(log_b_posterior)
b_posterior = exp(log_b_posterior)
plot(beta, b_posterior, ylab = "density_option_A", 
     type = "l", col = "orange") 
### 1b option B
log_b_posterior_2 = function (y, beta_i) {
  303*log(beta_i) - 
    beta_i*(4+sum(sapply(y, function(x) {1/x})))
}
log_b_posterior = sapply(beta, log_b_posterior_2, y = y)
log_b_posterior = log_b_posterior - max(log_b_posterior)
b_posterior = exp(log_b_posterior)
plot(beta, b_posterior, ylab = "density_option_B", 
     type = "l", col = "blue") 
### 1c
beta_map = beta[which.max(b_posterior)]
m=beta_map
var = (beta_map^2)/300
plot(beta, b_posterior, ylab = "density", 
     type = "l", col = "orange", 
     ylim=c(0,1.6*max(b_posterior))) 
b_clt_approxim = dnorm(beta, m, sqrt(var))
lines(beta, b_clt_approxim, col = "blue")
legend("topright", legend = c("unnorm posterior", "proposal"),
       col = c("orange", "blue"), lty = 1)
title("proposal distribution for MPH algorithm")
# NB: 
# norm_const is omitted as 
# pnorm(0,4.30131, sqrt( 0.06164194))=1.534768e-67

### 1d
mh = function(B, start, jump_parm) { # USE sigma^2 here for normal
  mu = jump_parm[1]
  sigma = sqrt(jump_parm[2])
  theta = numeric(B + 1)
  theta[1] = start
  
  for (i in 2:length(theta)) {
    theta_star = rnorm(1, mu, sigma)
    while (theta_star < 0) {
      theta_star = rnorm(1, mu, sigma)
    }
    num_logr = (log_b_posterior_2(y, theta_star) - 
                  max(log_b_posterior)) -
      dnorm(theta_star, mu, sigma, log = TRUE)
    den_logr = (log_b_posterior_2(y, theta[i-1]) - 
                  max(log_b_posterior)) -
      dnorm(theta[i-1], mu, sigma, log = TRUE)
    logr = num_logr - den_logr
    
    if (log(runif(1)) <= min(logr, 0)) {
      theta[i] = theta_star
    } else {
      theta[i] = theta[i - 1]
    }
  }
  return(theta)
}
set.seed(77)
B = 100000
keep = (B/2 + 1):(B + 1)
jpar = c(m, sqrt(var))
chain1 = mh(B, start = 3.3, jump_parm = jpar)
chain2 = mh(B, start = 3.8, jump_parm = jpar)
chain3 = mh(B, start = 4.3, jump_parm = jpar)
chain4 = mh(B, start = 4.8, jump_parm = jpar)
chain5 = mh(B, start = 5.2, jump_parm = jpar) 

mc = mcmc.list(mcmc(chain1[keep]),
               mcmc(chain2[keep]),
               mcmc(chain3[keep]),
               mcmc(chain4[keep]),
               mcmc(chain5[keep]))
# trace plot for the 5 chains
traceplot(mc)
summary (mc)

# acf cplot
autocorr.plot(mc)

# effective sample size
effectiveSize(mc)

# Calculate scale reduction factor
gelman.diag(mc, autoburnin = FALSE)

# Plot Gelman statistic
gelman.plot(mc, autoburnin = FALSE)

# Heidelberg & Welch
# If the halfwidth test fails, extend the chain(s)
heidel.diag(mc)

# Raftery-Lewis
# Want to the dependence factor less than 5.  If not the case,
# we'll need more samples
raftery.diag(mc)

# Geweke
# z scores for a test of equality (should be between -2 and 2)
# if converged
geweke.diag(mc)

### 2a option A
log_q_alpha_1 = function (y, alpha_i) {
  sum(dinvgamma(y, alpha_i, 4, log = TRUE)) + log(1/4)
}
alpha = seq(1, 5, length = 1000)
log_q_alpha_v_a = sapply(alpha, log_q_alpha_1, y=y)
plot(alpha, exp(log_q_alpha_v_a), 
     ylab = "density_option_A", 
     type = "l", col = "orange")

### 2a option B
log_q_alpha_2 = function (y, a_i) {
  100*a_i*log(4) + 
    (a_i+1)*sum(log(1/y)) -
    100 * log(factorial(a_i-1)) - 
    log(4) - 
    4 * sum(1/y)
}
log_q_alpha_v_b = sapply(alpha, log_q_alpha_2, y=y)
plot(alpha, exp(log_q_alpha_v_b), ylab = 
       "density_option_B", 
     type = "l", col = "blue")

### 2b
alpha = seq(2, 3.2, length = 1000)
log_q_alpha_v = sapply(alpha, log_q_alpha_1, y=y)
m_a=alpha[which.max(log_q_alpha_v)]
var_a = 0.022
proposal = dnorm(alpha, m_a, sqrt(var_a))
scale_const_a = 2.0e-76
plot(alpha, exp(log_q_alpha_v)/(scale_const_a), 
     ylab = "unnorm_post_dens_alpha_scaled", 
     type = "l", col = "green",
     ylim = c(0, max(proposal)*1.05))
 
lines(alpha, proposal, col = "blue")
abline(v=m_a, col = "black", lty=1)
legend("topright", 
       legend = c("unnorm_post", "proposal", "max"),
       col = c("green", "blue", "black"), lty = 1)
title("proposal distribution for the MPH algorithm")


### 2c
mh_a = function(B, start, jump_parm) { # USE sigma^2 here for normal
  mu = jump_parm[1]
  sigma = sqrt(jump_parm[2])
  theta = numeric(B + 1)
  theta[1] = start
  log_scale_const_a = log(scale_const_a)
  for (i in 2:length(theta)) {
    theta_star = rnorm(1, mu, sigma)
    while ((theta_star < 1)|(theta_star > 5)) {
      theta_star = rnorm(1, mu, sigma)
    }
    num_logr = log_q_alpha_1(y, theta_star) - 
                  log_scale_const_a -
      dnorm(theta_star, mu, sigma, log = TRUE)
    den_logr = (log_q_alpha_1(y, theta[i-1]) - 
                  log_scale_const_a) -
      dnorm(theta[i-1], mu, sigma, log = TRUE)
    logr = num_logr - den_logr
    
    if (log(runif(1)) <= min(logr, 0)) {
      theta[i] = theta_star
    } else {
      theta[i] = theta[i - 1]
    }
  }
  return(theta)
}
set.seed(77)
B = 100000
keep = (B/2 + 1):(B + 1)
jpar = c(m_a, sqrt(var_a))
chain1 = mh_a(B, start = 2.3, jump_parm = jpar)
chain2 = mh_a(B, start = 2.8, jump_parm = jpar)
chain3 = mh_a(B, start = 3.3, jump_parm = jpar)
chain4 = mh_a(B, start = 3.8, jump_parm = jpar)
chain5 = mh_a(B, start = 4.2, jump_parm = jpar) 

mc_a = mcmc.list(mcmc(chain1[keep]),
               mcmc(chain2[keep]),
               mcmc(chain3[keep]),
               mcmc(chain4[keep]),
               mcmc(chain5[keep]))
# trace plot for the 5 chains
traceplot(mc_a)
summary (mc_a)
# acf cplot
autocorr.plot(mc_a)
# effective sample size
effectiveSize(mc_a)
# Calculate scale reduction factor
gelman.diag(mc_a, autoburnin = FALSE)
# Plot Gelman statistic
gelman.plot(mc_a, autoburnin = FALSE)
# Heidelberg & Welch
# If the halfwidth test fails, extend the chain(s)
heidel.diag(mc_a)
# Raftery-Lewis
# Want to the dependence factor less than 5.  If not the case,
# we'll need more samples
raftery.diag(mc_a)
# Geweke
# z scores for a test of equality (should be between -2 and 2)
# if converged
geweke.diag(mc_a)


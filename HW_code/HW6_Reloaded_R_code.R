#source("http://math.ucdenver.edu/~jfrench/data/bayesdists.R")
#devtools::install_github("jfrench/bayesutils")
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
sum_inv = sum(1/y)
sigma = 0.3
beta_a = 1.01
beta_b = 1.01
### 1b
log_q_alpha = function (a, b) {
  100*a*log(b) + 
    (a+1)*sum(log(1/y)) -
    100 * log(factorial(a-1))
}
alpha_seq = seq(1,5, length=500)
beta_point = 6
plot(alpha_seq, log_q_alpha(a = alpha_seq, b=beta_point), ylab = "log_q_alpha", 
     type = "l", col = "blue")
title("beta_point = 6")

plot(alpha_seq, dbeta((alpha_seq-1)/4, 1.01,1.01), ylab = "density", 
     type = "l", col = "blue")
title("Beta(1.01, 1.01)*4 + 1")


### 1b
alpha_point = 2
beta_seq = seq(1, 7, length = 1000)

d_beta = function(a, b) {
  dgamma(b, shape = (100*a + 4), rate = (4 + sum_inv))
}
plot(beta_seq, d_beta(a = alpha_point, b=beta_seq), ylab = "d_beta", 
     type = "l", col = "orange")
legend("topright", 
       legend = c("d_beta"), 
       col = c("orange"), lty = 1, bty = "n")
title("alpha_point = 4")



### 2c
mh_gibbs = function(B, theta_start) {
  theta = array(0, c((B+1), 2))
  theta[1,1] = theta_start[1]
  theta[1,2] = theta_start[2]

  for (i in 2:dim(theta)[1]) {
    alpha_star = theta[(i-1),1]
    beta_star = rgamma(1, shape = (100*alpha_star+4), rate = (4 + sum_inv))
    theta[i,2] = beta_star
    alpha_star = rbeta(1, beta_a, beta_b)*4+1
    #alpha_star = runif(1, 1, 5)
    #alpha_star = rnorm(1, theta[(i-1),1], sigma) 
    #while ((alpha_star < 1)|(alpha_star > 5)) {
    #  alpha_star = rnorm(1, theta[(i-1),1], sigma) 
    #}
    ### why does rnorm work better than runif??

    num_logr = log_q_alpha(a=alpha_star, b=beta_star) -
      dbeta((alpha_star-1)/4,beta_a,beta_b)
      #dnorm(alpha_star, theta[i-1, 1], sigma, log = TRUE) 
      #- 
      #dunif(alpha_star, 1, 5, log = TRUE) 
    ### think it is possible  not to correct for the proposal/ jumping (both for unif and norm) ###
    ### as it is symmetric and thus cancels in den_logr ###
    den_logr = log_q_alpha(a=theta[i-1, 1], b=beta_star) -
      dbeta((theta[i-1, 1]-1)/4,beta_a,beta_b)
    
      #dnorm(theta[i-1, 1], alpha_star, sigma, log = TRUE) 
    #-
      #dunif(alpha_star, 1, 5, log = TRUE)
    ### think it is possible not to correct for the proposal/ jumping (both for unif and norm) ###
    ### as it is symmetric and thus cancels in den_logr ###
    logr = num_logr - den_logr
    if (log(runif(1)) <= min(logr, 0)) {
      theta[i,1] = alpha_star
    } else {
      theta[i,1] = theta[(i - 1), 1]
    }
  }
  return(theta)
}
B = 10^5
keep = (B/2 + 1):(B + 1)
chain1 = mh_gibbs(B, theta_start = c(2.3,4))
chain2 = mh_gibbs(B, theta_start = c(2.8,5))
chain3 = mh_gibbs(B, theta_start = c(3.3,6))
chain4 = mh_gibbs(B, theta_start = c(3.9,8))
chain5 = mh_gibbs(B, theta_start = c(4.2,10)) 

chain = rbind(chain1[keep,], chain2[keep,], chain3[keep,], chain4[keep,], chain5[keep,])
mc_a = mcmc.list(mcmc(chain1[keep,1]), mcmc(chain2[keep,1]),
                 mcmc(chain3[keep,1]), mcmc(chain4[keep,1]),
                 mcmc(chain5[keep,1]))
summary (mc_a)
mc_b= mcmc.list(mcmc(chain1[keep,2]), mcmc(chain2[keep,2]),
                 mcmc(chain3[keep,2]), mcmc(chain4[keep,2]),
                 mcmc(chain5[keep,2]))
summary (mc_b)



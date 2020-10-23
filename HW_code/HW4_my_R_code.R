### HW4 ###
#library(MCMCvis)
#library(mcmcplots)
library(bayesplot)
library(ggplot2)
library(rstan)

### Problem 1 ###
### Pr_1_a_b ###
theta = seq(10^-12, 0.9999999, length = 1000)
qtheta = function(theta, y) {
  dbeta(theta, 1.1, 1.5) * dnbinom(y, 200, theta)
  }
y=303

pr = dnorm(theta, 0.4, 0.04)
qt = qtheta(theta, y)

m=max(qt)/ max(pr)
while (sum(m * pr - qt < 0)>0) {
  m=m*1.01
}
Mg = function(theta) {
  m * dnorm(theta, mean = 0.4, sd = 0.04)
}
 
# plot qtheta
plot(theta, qt, ylab = "density", type = "l", col = "orange", 
     ylim=c(0,max(1.1*max(m*pr), 1.1*max(qt))))
# plot Mg(theta)
lines(theta, m * pr, col = "blue")
legend("topright", legend = c("qtheta", "Mg"),
       col = c("orange", "blue"), lty = 1, bty = "n")
title("qtheta and Mg vs theta")

### Pr_1_c ###
B = 10^6
mytheta = numeric(B)

i = 0 # the samples accepted
while (i < B) {
  x = rnorm(1, 0.4, 0.04) # sample from g distribution
  
  # accept x with probability q(x)/Mg(x)
  if (runif(1) <= qtheta(x, y)/Mg(x)) {
    i = i + 1
    mytheta[i] = x
  }
}
theta = seq(0.35, 0.45, length = 1000)

dmytheta = density(mytheta)
dtruth = qtheta(theta, y)/0.0009461865

plot(theta, dtruth, col = "orange", type = "l")
lines(dmytheta)
legend("topleft", legend = c("approximation", "truth"), 
       col = c("black", "orange"), lwd = c(1, 1), bty = "n")

### Problem 2 ###
### Pr_2_c ###

gibbs = function(theta, B) {
  #create matrix to store samples
  theta_sims = matrix(0, nrow = B, ncol = 2)
  # run gibbs sampler for B cycles
  for (i in 1:B) {
    # simulate from full conditional distribution for theta1
    theta[1] = rgamma(1,shape=3, rate=theta[2]^2+4)
    # simulate from full conditional distribution for theta2
    theta[2] = rnorm(1, mean = 1/(theta[1]+1), sd = sqrt(0.5/(theta[1]+1)))
    # save sample
    theta_sims[i, ] = theta
  }
  return(theta_sims)
}
theta = c(2,-1)
b=10^5
theta_gibbs = gibbs(theta, B=b)
plot(theta_gibbs, pch = ".", xlab = expression(theta[1]),
     ylab = expression(theta[2]))
title("Samples from Gibbs sampler")
plot(density(theta_gibbs[,1]), main = "", xlab = "theta_1", xlim = c(0, 3))
plot(density(theta_gibbs[,2]), main = "", xlab = "theta_2", xlim = c(-3, 4))

### Pr_2_d ###

# mean, median, var and 95% posterior interval for theta_1 and thrta_2
# we take only the second part of the chain as it is closer to the convergence
summary(theta_gibbs[(b/2):b,1])
var(theta_gibbs[(b/2):b,1])
quantile(theta_gibbs[(b/2):b,1], c(.025, .975))

summary(theta_gibbs[(b/2):b,2])
var(theta_gibbs[(b/2):b,2])
quantile(theta_gibbs[(b/2):b,2], c(.025, .975))

### Pr_2_e ###

#plot path of gibbs sampler for 100 iterations
plot.mcmc.path(theta_gibbs[1:100,], x0 = theta,
               xlab = expression(theta_1), ylab = expression(theta_2))

#plot path of cycles for 100 iterations
plot(rbind(theta, theta_gibbs[1:100,]),
     xlab = expression(theta_1), ylab = expression(theta_2), type = "l")



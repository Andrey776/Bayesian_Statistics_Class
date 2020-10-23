# range of theta values
theta = seq(10^-12, 0.9999999, length = 1000)
# plot prior
prior = dbeta(theta, 1.1, 1.5)
plot(theta, prior , ylab = "density", type = "l", col = "orange")
# plot likelihood
likelihood = dnbinom(303, 200, theta)
k = max(prior)/ max(likelihood)
lines(theta, likelihood * k, col = "blue")
legend("topright", legend = c("prior", "likelihood"),
       col = c("orange", "blue"), lty = 1, bty = "n")
title("Likelihood (scaled) and prior density vs theta")

### Problem 2 ###
ltheta = function(theta, y) {
  dnbinom(y, 200, theta)
}
qtheta = function(theta, y) {
  dbeta(theta, 1.1, 1.5) * dnbinom(y, 200, theta)
  }
y=303
theta = seq(0.39, 0.405, length = 1000)
qt = qtheta(theta, y)
lt = ltheta(theta, y)
plot(theta, qt, ylab = "density", type = "l", col = "orange", ylim=c(0,1.25*max(qt)))
lines(theta, lt, col = "blue")
legend("topleft", legend = c("q_posterior", "likelihood"),
       col = c("orange", "blue"), lty = 1, bty = "n")
title("Likelihood and posterior density vs theta")

map = optimize(qtheta, interval = c(0.2, 0.6), y = y, maximum = TRUE)
map$maximum

mle = optimize(ltheta, interval = c(0.2, 0.6), y = y, maximum = TRUE)
mle$maximum

abline(v=map$maximum, col = "orange", lty = 2)
abline(v=mle$maximum, col = "blue", lty = 2)
legend("topright", legend = c("map", "mle"),
       col = c("orange", "blue"), lty = 2, bty = "n")

### problem 3a ###
theta = seq(10^-12, 0.9999999, length = 1000)
y = 303
poster_scale_robust = 10

qtheta_robust = function(theta, y) {
  log_pprior = dbeta(theta, 1.1, 1.5, log = TRUE) 
  log_likelihood = sum(dnbinom(y, 200, theta, log = TRUE))
  exp(log_pprior + log_likelihood + poster_scale_robust)
  }
vqtheta_robust = Vectorize(qtheta_robust, vectorize.args = "theta")
nconst = integrate(f = vqtheta_robust, lower = 0, upper = 1, y = y)

dpost = function(theta, y) {
  qtheta_robust(theta, y)/nconst$value
}

p_data = nconst$value / exp(poster_scale_robust)
p_data

vdpost = Vectorize(dpost, vectorize.args = "theta")
nconst_double_check = integrate(f = vdpost, lower = 0, upper = 1, y = y)
nconst_double_check
  
### Problem 3b ###
# define function to integrate over to find the posterior mean
e = function(theta) {
  theta * vdpost(theta, y = y)
}

# Mean of posterior
e_theta = integrate(e, 0, 1)
e_theta$value

### Problem 3c ###
# define function to integrate over to find the posterior VAR
e = function(theta) {
  (theta^2) * vdpost(theta, y = y)
}
# VAR of posterior
e_theta_sqr = integrate(e, 0, 1)
var_theta = e_theta_sqr$value - e_theta$value^2
var_theta

### Problem 4 ###
# normal approximation posterior
observed_inf = function(theta) {
  200.1/theta^2 + 303.5/(1-theta)^2 
}

dpapprox_norm = function(theta, y, thetahat) {
  dnorm(theta, mean = thetahat, sd = sqrt(1/observed_inf(thetahat)))
}
# vectorize over theta
vdapprox_norm = Vectorize(dpapprox_norm, vectorize.args = "theta")
theta = seq(10^-12, 0.9999999, length = 1000)
e = function(theta) {
  vdapprox_norm(theta, y, thetahat = map$maximum)
}
integrate(e,0,1) # to double-check
# plot true density
theta1 = seq(0.392, 0.404, length = 1000)

plot(theta1, vdpost(theta1, y), ylab = "density", type = "l", col = "orange")
# plot normal approximation
lines(theta1, vdapprox_norm(theta1, y, thetahat = map$maximum), col = "blue")
legend("topright", legend = c("true", "approximation"),
       col = c("orange", "blue"), lty = 1)
title("Normal approximation to posterior")

# define function to integrate over to find the posterior mean
e = function(theta) {
  theta * vdapprox_norm(theta, y, thetahat = map$maximum)
}
# Mean of posterior
e_theta_approx_norm = integrate(e, 0, 1)
e_theta_approx_norm$value
# define function to integrate over to find the posterior VAR
e = function(theta) {
  (theta^2) * vdapprox_norm(theta, y = y, thetahat = map$maximum)
}
# VAR of posterior
e_theta_approx_norm_sqr = integrate(e, 0, 1)
var_theta_approx_norm = e_theta_approx_norm_sqr$value - e_theta_approx_norm$value^2
var_theta_approx_norm

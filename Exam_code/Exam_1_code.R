### problem 6b ###
beta = seq(10^-1, 20, length = 1000)
l = function(beta) {
  303*log(beta) - 70*beta
}
plot(beta, l(beta), ylab = "L(beta)", 
     type = "l", col = "orange")
legend("bottomright", 
       legend = c("303*log(beta) - 70*beta"), 
       col = c("orange"), lty = 1, bty = "n")

### problem 6c ###
beta = seq(10^-1, 20, length = 1000)
l = function(beta) {
  303*log(beta) - 70*beta
}
  
map = optimize(l, interval = c(1, 10), maximum = TRUE)
map$maximum
plot(beta, l(beta), ylab = "L(beta)", 
     type = "l", col = "orange")
abline(v=map$maximum, col = "blue")
legend("bottomright", 
       legend = c("303*log(beta) - 70*beta", "map"), 
       col = c("orange", "blue"), lty = 1, bty = "n")

### Problem 7c ###
# normal approximation posterior
m = 4.328577
var =  0.06245526
a = 304
b = 70

dpapprox_norm = function(beta, m, var) {
  dnorm(beta, mean = m, sd = sqrt(var))
}

beta_0 = seq(2, 7, length = 1000)
beta_1 = seq(4.0, 4.6, length = 1000)

plot(beta_0, dgamma(beta_0, a, b), ylab = "density", 
     type = "l", col = "orange")

# plot normal approximation
lines(beta_0, dpapprox_norm(beta_0, m, var), col = "blue")
legend("topright", legend = c("true", "approximation"),
       col = c("orange", "blue"), lty = 1)
title("Normal approximation to posterior")

plot(beta_1, dgamma(beta_1, a, b), ylab = "density", 
     type = "l", col = "orange")
# plot normal approximation
lines(beta_1, dpapprox_norm(beta_1, m, var), col = "blue")
legend("bottomright", legend = c("true", "approximation"),
       col = c("orange", "blue"), lty = 1)
title("Normal approximation to posterior")






### Problem 8c ###
m = 4.328577
var =  0.06245526
a = 304
b = 70
precision=10^-5
q = qgamma(c(precision, 1-precision), a, b)
q
beta_0 = seq(q[1], q[2], length = 1000)
qp_beta = function(beta) {
  (beta^(a-1))*exp(-b*beta)
}
dpapprox_norm = function(beta) {
  dnorm(beta, mean = m, sd = sqrt(var))
}
beta_appr = dpapprox_norm(beta_0)
q_beta = qp_beta(beta_0)

M=max(q_beta)/ max(beta_appr)
while (sum(M * beta_appr - q_beta < 0)>0) {
  M=M*1.01
}
M

plot(beta_0, q_beta, ylab = "density", 
     type = "l", col = "orange", 
     ylim=c(0,1.05*max(M*beta_appr)))
lines(beta_0, M*beta_appr, col = "blue")
legend("topright", legend = c("q(b|y)", "Mg(b)"),
       col = c("orange", "blue"), lty = 1)
title("Unnormalized posterior and bounding function")

### Problem_8d ###
m = 4.328577
var =  0.06245526
a = 304
b = 70
sd = sqrt(var)
size = 10^5
beta = seq(q[1], q[2], length = 1000) #q and M are from #8c

qp_beta_unnorm = function(beta) {
  (beta^(a-1))*exp(-b*beta)
}

qp_beta_unnorm_rej = numeric(size)
n=1
while (n < size) {
  x = rnorm(1, mean=m, sd=sd)
  if (runif(1,0,1) < qp_beta_unnorm(x)/(M*dnorm(x, mean=m, sd=sd))) {
    qp_beta_unnorm_rej[n] = x
    n=n+1
  }
}
dtruth = dgamma(beta, a, b)
plot(beta, dtruth, col = "black", type = "l")
lines(density(qp_beta_unnorm_rej), col = "red")
legend("topright", legend = c("truth", "simulation"),
       lwd = c(1, 1), col = c("black", "red"))

### Problem_9b ###
p_alpha_unnorm = function(a) {
  (((4^a)/factorial(a-1))^100) * exp((a+1)*(-64))
}
alpha = seq(3, 5, length = 10^3)
density = p_alpha_unnorm(alpha)
plot(alpha, density, col = "blue", type = "l")
legend("topright", legend = c("unnorm. full conditional density for alpha"),
       lwd = c(1, 1), col = c("blue"))

optimize(p_alpha_unnorm,interval = c(1,5), maximum = TRUE)


### Problem_9c ###
a=3
p_beta_unnorm = function(b, a) {
  (b^(100*a+3)) * exp((b+4)*(-256))
}

p_beta_unnorm(4.32,3)

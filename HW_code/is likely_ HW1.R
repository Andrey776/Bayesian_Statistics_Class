t <- seq(0.00001,0.99999,.01)
dens1 <- dbeta(t,0.5,1)# * (dbinom(0,10,t) + dbinom(1,10,t) + dbinom(2,10,t))
#dens2 <- t^3 * (1-t)^13 + 10* t^4 * (1-t)^12 + 45*t^5 * (1-t)^11

plot (t, dens1, ylim=c(0,1.1*max(dens1)), type="l", cex=2)
#plot (t, dens2, ylim=c(0,1.1*max(dens2)), type="l", cex=2)

y1 = 0.5
d = 0.3*dnorm(y1,5,2) + 0.7*dchisq(y1,3)
n = 0.3*dnorm(y1,5,2)
n/d

pbinom(5,10,.5) - pbinom(4,10,.5)
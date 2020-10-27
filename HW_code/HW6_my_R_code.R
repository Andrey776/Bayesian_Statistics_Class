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

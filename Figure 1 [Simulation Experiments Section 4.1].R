################
### Figure 1 ###
################

simulation <- function(n, d, sigma) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  treatment_propensity <- pmax(0.1, pmin(sin(pi*X[,1]*X[,2]), 0.9))
  main_effect <- sin(pi*X[,1]*X[,2]) + 2*(X[,3] - 0.5)^2 + X[,4] + 0.5*X[,5]
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  zeta1 <- 1 + 1/(1 + exp(-20*(X[,1] - (1/3))))
  zeta2 <- 1 + 1/(1 + exp(-20*(X[,2] - (1/3))))
  Tau <- matrix(zeta1 * zeta2, nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, treatment_propensity = treatment_propensity, W = W, Y = Y, Tau = Tau)
  return(output)
}

s <- 1
n <- 10000
d = 25

dat <- simulation(n = n, d = d, sigma = s)
treatment_propensity.for.plot <- dat$treatment_propensity

minp = min(treatment_propensity.for.plot)
maxp = max(treatment_propensity.for.plot)
rngp = maxp-minp
ncol = 100
true.scl = pmax(ceiling(ncol * (treatment_propensity.for.plot - minp) / rngp), 1)
hc = heat.colors(ncol)
plot(X.for.plot[,1], X.for.plot[,2], pch = 16, col = hc[true.scl], xlab = "X1", ylab = "X2", main = "")


simulation <- function(n, d, sigma) {
  X <- matrix(runif(n*d, min = 0, max = 1), nrow = n, ncol = d)
  treatment_propensity <- pmax(0.1, pmin(sin(pi*X[,1]*X[,2]), 0.9))
  main_effect <- sin(pi*X[,1]*X[,2]) + 2*(X[,3] - 0.5)^2 + X[,4] + 0.5*X[,5]
  noise <- sigma*rnorm(n=n)
  W <- matrix(rbinom(n = n, size = 1, prob = treatment_propensity), nrow = n, ncol = 1)
  zeta1 <- 1 + 1/(1 + exp(-20*(X[,1] - (1/3))))
  zeta2 <- 1 + 1/(1 + exp(-20*(X[,2] - (1/3))))
  Tau <- matrix(zeta1 * zeta2, nrow = n, ncol = 1)
  Y <- matrix(main_effect + (W - 0.5) * Tau + noise, nrow = n, ncol = 1)
  output <- list(X = X, treatment_propensity = treatment_propensity, W = W, Y = Y, Tau = Tau)
  return(output)
}
s <- 1
n <- 10000
d = 25

tau.for.plot <- dat$Tau

minp = min(tau.for.plot)
maxp = max(tau.for.plot)
rngp = maxp-minp
ncol = 100
true.scl = pmax(ceiling(ncol * (tau.for.plot - minp) / rngp), 1)
hc = heat.colors(ncol)
plot(X.for.plot[,1], X.for.plot[,2], pch = 16, col = hc[true.scl], xlab = "X1", ylab = "X2", main = "")

